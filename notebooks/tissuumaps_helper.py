import warnings
warnings.filterwarnings("ignore")

import os
import json
import numpy as np
import pandas as pd
import pyvips
import geobuf
from typing import Dict
from skimage.measure import approximate_polygon
from spatialdata.models import get_table_keys
from spatialdata._core.operations._utils import transform_to_data_extent
from spatialdata import SpatialData, get_centroids, get_extent, to_circles
from spatialdata_io.experimental import to_legacy_anndata
from tissuumaps import read_h5ad
from tissuumaps_schema import CURRENT_SCHEMA_MODULE
from tissuumaps_schema.base import RootSchemaBaseModel
from anndata import AnnData

# Constants
VIPS_JPEG_COMPRESSION = 100
VIPS_EXCLUDE_MIN_INTENSITY = False
VIPS_FORCE_RESCALE = True

def add_image(sdata: SpatialData, image_element: str, output_folder: str, tmap_state: Dict) -> None:
    """
    Adds an image element from a SpatialData object to a specified output folder in a tiled and rescaled format. 
    Updates TissUUmaps `tmap_state` dictionary with information about the newly added image.

    Parameters:
        sdata (SpatialData): The SpatialData object containing image data and metadata. 
                             The image data is accessed through `sdata[image_element].data`.
        image_element (str): The name of the image element to be processed and saved.
        output_folder (str): Path to the folder where the output image and its associated tiled representation will be stored.
        tmap_state (Dict): A TissUUmaps project state dictionary. 
                           The function appends the image information to the "layers" key of this dictionary.

    Functionality:
        1. Extracts the image data corresponding to `image_element` from the `sdata` object.
        2. Saves the image in the output folder.
        3. Applies 8 bit rescaling to the image data.
           - Handles cases where the image has negative values or values outside this range.
           - Applies an exclusion for minimum intensity if configured.
        4. Saves the image as a tiled, JPEG-compressed TIFF file using `pyvips`.
        5. Adds an entry to `tmap_state["layers"]`, linking the image element to its tiled source.

    Returns:
        None
    """
    print ("  - Adding image", image_element)
    coordinate_system = sdata.coordinate_systems[0]

    output_path = f"{output_folder}/{coordinate_system}_coord/images/{image_element}.tif"
    if not os.path.exists(output_path):
        os.makedirs(f"{output_folder}/{coordinate_system}_coord/images", exist_ok=True)
        image_data = sdata[image_element].data.compute()[0, :, :].astype(np.float32)

        img_vips = pyvips.Image.new_from_memory(
            image_data.tobytes(),
            image_data.shape[1],
            image_data.shape[0],
            bands=1,
            format="float",
        )

        min_val = img_vips.percent(0)
        max_val = img_vips.percent(100)

        if VIPS_EXCLUDE_MIN_INTENSITY:
            absolute_min_val = img_vips.min()
            img_vips_tmp = (img_vips == absolute_min_val).bandand().ifthenelse(max_val + 1, img_vips)
            min_val = img_vips_tmp.percent(0)

        if min_val == max_val:
            min_val, max_val = 0, 255

        if VIPS_FORCE_RESCALE or img_vips.min() < 0 or img_vips.max() > 255:
            img_vips = (255.0 * (img_vips - min_val)) / (max_val - min_val)
            img_vips = (img_vips < 0).ifthenelse(0, img_vips)
            img_vips = (img_vips > 255).ifthenelse(255, img_vips)

        img_vips = img_vips.scaleimage()
        img_vips.tiffsave(
            output_path,
            pyramid=True,
            tile=True,
            tile_width=256,
            tile_height=256,
            compression="jpeg",
            Q=VIPS_JPEG_COMPRESSION,
            properties=True,
        )

    tmap_state["layers"].append({
        "name": image_element,
        "tileSource": f"{coordinate_system}_coord/images/{image_element}.tif.dzi",
    })


def add_shape(sdata: SpatialData, shape_element: str, output_folder: str, tmap_state: Dict) -> None:
    """
    Adds a shape element from a SpatialData object to a specified output folder in an approximate format. 
    Updates TissUUmaps `tmap_state` dictionary with information about the newly added shape.

    Parameters:
        sdata (SpatialData): The SpatialData object containing shape data and metadata. 
                             The shape data is accessed through `sdata[shape_element]`.
        shape_element (str): The name of the shape element to be processed and saved.
                             The special case "cells" is ignored and skipped by the function.
        output_folder (str): Path to the folder where the output shape file will be stored.
        tmap_state (Dict): A TissUUmaps project state dictionary. 
                           The function appends the shape information to the "regionFiles" key of this dictionary.

    Functionality:
        1. Converts the shape data to GeoJSON format and approximates the polygons in the shape data with a tolerance of 0.75.
        2. Encodes the approximated GeoJSON data into a Geobuf format with a precision of 3 and saves it as a `.pbf` file.
        3. Adds an entry to `tmap_state["regionFiles"]`, linking the shape element to its `.pbf` source.

    Returns:
        None
    """

    print ("  - Adding shape", shape_element)
    coordinate_system = sdata.coordinate_systems[0]

    os.makedirs(f"{output_folder}/{coordinate_system}_coord/shapes", exist_ok=True)

    output_path = f"{output_folder}/{coordinate_system}_coord/shapes/{shape_element}.pbf"
    if not os.path.exists(output_path):
        geojson_data = json.loads(sdata[shape_element].to_json())

        geojson_data_approx = {
            "type": "FeatureCollection",
            "features": [
                {
                    **f,
                    "properties": {**f["properties"], "collectionIndex": 0},
                    "geometry": {
                        **f["geometry"],
                        "coordinates": [
                            approximate_polygon(np.array(c), tolerance=0.75).tolist()
                            for c in f["geometry"]["coordinates"]
                        ],
                    },
                }
                for f in geojson_data["features"]
            ],
        }

        with open(output_path, "wb") as f:
            pbf = geobuf.encode(geojson_data_approx, 3)
            data = geobuf.geobuf_pb2.Data()
            data.ParseFromString(pbf)
            data.precision = 3
            f.write(data.SerializeToString())

    tmap_state["regionFiles"].append({
        "title": shape_element,
        "path": f"{coordinate_system}_coord/shapes/{shape_element}.pbf",
    })


def add_table(sdata: SpatialData, table_element: str, output_folder: str, tmap_state: Dict) -> None:
    """
    Adds a table element from a SpatialData object to a specified output folder in a `.h5ad` format. 
    Updates TissUUmaps `tmap_state` dictionary with information about the newly added table.

    Parameters:
        sdata (SpatialData): The SpatialData object containing table data and metadata. 
                             The table data is accessed through `sdata[table_element]`.
        table_element (str): The name of the table element to be processed and saved.
        output_folder (str): Path to the folder where the output table file will be stored.
        tmap_state (Dict): A TissUUmaps project state dictionary. 
                           The function appends the table information to the "markerFiles" key of this dictionary.

    Functionality:
        1. Retrieves the table data (`table_adata`) and determines the associated spatial region.
        2. Computes the centroids of the spatial region and stores them in `table_adata.obsm["spatial"]`.
        3. Saves the table as an `.h5ad` file in the specified output folder.
        4. Appends the processed marker file entry to `tmap_state["markerFiles"]`.

    Returns:
        None
    """

    print ("  - Adding table", table_element)
    coordinate_system = sdata.coordinate_systems[0]
    os.makedirs(f"{output_folder}/{coordinate_system}_coord/tables", exist_ok=True)
    table_adata = to_legacy_anndata(
        sdata=sdata,
        coordinate_system=coordinate_system,
        table_name=table_element,
        include_images=False
    )

    output_path = f"{output_folder}/{coordinate_system}_coord/tables/{table_element}.h5ad"
    table_adata.write(output_path)

    tmap = read_h5ad.h5ad_to_tmap("", output_path, library_id=None)
    single_markerfile = tmap["markerFiles"][0]
    single_markerfile["dropdownOptions"] = [
        option for markerfile in tmap["markerFiles"]
        for option in markerfile["dropdownOptions"]
    ]
    single_markerfile.update({
        "path": f"{coordinate_system}_coord/tables/{single_markerfile['path']}",
        "name": table_element,
        "title": table_element,
        "uid": table_element.replace("_",""),
    })
    single_markerfile["expectedHeader"]["cb_cmap"] = "interpolateViridis"
    extent = get_extent(sdata, coordinate_system=coordinate_system, exact=False)
    
    region, _, instance_key  = get_table_keys(table_adata)
    if isinstance(region, list):
        region = region[0]
    region_circles = to_circles(sdata[region])
    radius = region_circles["radius"].iloc[0]
    tmap_factor = 0.0037244394752190928
    scale = radius * 2 / (extent["x"][1] - extent["x"][0])
    scale = scale / tmap_factor
    single_markerfile["expectedHeader"]["scale_factor"] = scale

    tmap_state["markerFiles"].append(single_markerfile)


def add_points(sdata: SpatialData, points_element: str, output_folder: str, tmap_state: Dict) -> None:
    """
    Adds a points element from a SpatialData object to a specified output folder in a `.h5ad` format.
    Updates TissUUmaps `tmap_state` dictionary with information about the newly added points.

    Parameters:
        sdata (SpatialData): The SpatialData object containing points data and metadata. 
                             The points data is accessed through `sdata[points_element].compute()`.
        points_element (str): The name of the points element to be processed and saved.
        output_folder (str): Path to the folder where the output points file will be stored.
        tmap_state (Dict): A TissUUmaps project state dictionary. 
                           The function appends the points information to the "markerFiles" key of this dictionary.

    Functionality:
        1. Retrieves the points data and computes the centroids of the spatial region.
        2. Saves the points as an `.h5ad` file in the specified output folder.
        3. Appends the processed marker file entry to `tmap_state["markerFiles"]`.

    Returns:
        None
    """
    print ("  - Adding points", points_element)
    coordinate_system = sdata.coordinate_systems[0]
    os.makedirs(f"{output_folder}/{coordinate_system}_coord/points", exist_ok=True)
    points = sdata[points_element].compute()
    adata = AnnData(obs=points)

    centroids = get_centroids(sdata[points_element], coordinate_system=sdata.coordinate_systems[0]).compute().to_numpy()
    adata.obsm["spatial"] = centroids

    for column in adata.obs.columns:
        if not pd.api.types.is_numeric_dtype(adata.obs[column]):
            adata.obs[column] = adata.obs[column].astype(str).astype("category")

    output_path = f"{output_folder}/{coordinate_system}_coord/points/{points_element}.h5ad"
    adata.write(output_path)

    tmap = read_h5ad.h5ad_to_tmap("", output_path, library_id=None)
    single_markerfile = tmap["markerFiles"][-1].copy()
    single_markerfile["dropdownOptions"] = [
        option for markerfile in tmap["markerFiles"] for option in markerfile["dropdownOptions"]
    ]

    if not single_markerfile["dropdownOptions"]:
        return

    single_markerfile.update({
        "path": f"{coordinate_system}_coord/points/{single_markerfile['path']}",
        "name": points_element,
        "title": points_element,
        "uid": points_element.replace("_",""),
    })
    single_markerfile["expectedHeader"]["cb_cmap"] = "interpolateViridis"

    tmap_state["markerFiles"].append(single_markerfile)


def get_tmap_state() -> Dict:
    model_type = getattr(CURRENT_SCHEMA_MODULE, "Project", None)
    if model_type is None or not issubclass(model_type, RootSchemaBaseModel):
        raise ValueError("Invalid model type in the current schema version.")

    json_data = model_type.construct(by_alias=True)
    return json_data.dict(by_alias=True, exclude_none=True)


def spatialdata_to_tissuumaps(
    sdata: SpatialData, 
    elements: Dict[str, str],
    output_folder: str, 
    coordinate_system: str
) -> None:
    os.makedirs(output_folder, exist_ok=True)
    if len(sdata.coordinate_systems) > 1:
        sdata = transform_to_data_extent(
            sdata, coordinate_system=coordinate_system, target_unit_to_pixels=1, maintain_positioning=False
        )
    tmap_state = get_tmap_state()

    for image_element in elements["images"]:# or sdata.images:
        try:
            add_image(sdata, image_element, output_folder, tmap_state)
        except Exception as e:
            print(f"Failed to add image {image_element}: {e}")

    for shape_element in elements["shapes"]:# or sdata.shapes:
        try:
            add_shape(sdata, shape_element, output_folder, tmap_state)
        except Exception as e:
            print(f"Failed to add shape {shape_element}: {e}")

    for table_element in elements["tables"]:# or sdata.tables:
        try:
            add_table(sdata, table_element, output_folder, tmap_state)
        except Exception as e:
            print(f"Failed to add table {table_element}: {e}")

    for points_element in elements["points"]:# or sdata.points:
        try:
            add_points(sdata, points_element, output_folder, tmap_state)
        except Exception as e:
            print(f"Failed to add points {points_element}: {e}")
    
    # Save tmap file
    with open(f"{output_folder}/tissuumaps_{coordinate_system}.tmap", "w") as f:
        json.dump(tmap_state, f, indent=2)
    
    return sdata, f"{output_folder}/tissuumaps_{coordinate_system}.tmap"

# Add main function with parameters for input file and output folder
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Convert SpatialData to TissUUmaps")
    parser.add_argument("input_file", type=str, help="Path to the input SpatialData file")
    parser.add_argument("output_folder", type=str, help="Path to the output folder")
    parser.add_argument("--coordinate_systems", type=str, nargs="+", help="List of coordinate systems to convert to")

    args = parser.parse_args()

    sdata = SpatialData.read(args.input_file)

    if args.coordinate_systems is None:
        coordinate_systems = sdata.coordinate_systems
    else:
        coordinate_systems = args.coordinate_systems

    for cs in coordinate_systems:
        print(f"Converting to TissUUmaps in coordinate system: {cs}")
        sdata_extended, tmap_file = spatialdata_to_tissuumaps(
            sdata, 
            elements = {
                "images": sdata.images,
                "shapes": sdata.shapes,
                "tables": sdata.tables,
                "points": sdata.points
            },
            output_folder = args.output_folder,
            coordinate_system = cs
        )