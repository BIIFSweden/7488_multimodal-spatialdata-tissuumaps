from collections import defaultdict
from pathlib import Path
from typing import Any, Dict

import anndata
import geojson
import numpy as np
import os
import pandas as pd
import tifffile
from geopandas import GeoDataFrame
from roifile import ImagejRoi, ROI_TYPE
from shapely.geometry import Polygon
from spatialdata.models import Image2DModel, PointsModel, ShapesModel, TableModel
from spatialdata import concatenate
from tqdm import tqdm
import json

def load_table(
    sdata: Any, 
    input_file: str, 
    element_name: str,
    transformations: Dict[str, Any]
) -> None:
    """
    Load cell types data into an anndata table and add it to SpatialData.

    Parameters:
        sdata (Any): SpatialData object.
        input_file (str): Path to the CSV file containing cell type data.
        transformations (Dict[str, Any]): Transformation dictionary for spatial coordinates.
    """
    table_name = element_name + "_cells"
    cell_types = pd.read_csv(input_file)

    obsm = {"spatial": cell_types[["X", "Y"]].to_numpy()}
    # Keep columns in ["roi", "Cell.type", "Number", "Area"] if they exist
    obs = cell_types[[col for col in ["roi", "Cell.type", "Number", "Area"] if col in cell_types.columns]]
    obs["region_key"] = table_name
    obs["Index"] = obs.index
    
    #X = cell_types.drop(columns=["X", "Y", "roi", "Cell.type", "Number", "Area"])
    # Check existence of other columns
    exclude_cols = ["X", "Y", "roi", "Cell.type", "Number", "Area", "Unnamed: 0"]
    X = cell_types.drop(columns=[col for col in exclude_cols if col in cell_types.columns])

    adata = anndata.AnnData(X=X, obs=obs, obsm=obsm)
    
    sdata.tables[element_name] = TableModel.parse(
        adata, region=table_name, region_key="region_key", instance_key="Index"
    )

    circles = ShapesModel.parse(
        adata.obsm["spatial"],
        geometry=0,
        radius=40,
        index=adata.obs["Index"],
        transformations=transformations,
    )
    sdata.shapes[table_name] = circles

def load_image(
    sdata: Any, tif_path: str, element_name: str, transformations: Dict[str, Any]
) -> None:
    """
    Load a TIFF image and add it to SpatialData.

    Parameters:
        sdata (Any): SpatialData object.
        tif_path (str): Path to the TIFF image file.
        element_name (str): Name of the image element in SpatialData.
        transformations (Dict[str, Any]): Transformation dictionary for spatial coordinates.
    """
    image_data = tifffile.imread(tif_path)

    if len(image_data.shape) == 2:
        image_data = np.expand_dims(image_data, axis=0)

    elif image_data.ndim == 3:
        y, x = image_data.shape[:2]
        c_last = image_data.shape[2]

        # Detect if input is (y, x, c) → transpose to (c, y, x)
        if c_last <= 4 and image_data.shape[0] > 4:
            # (y, x, c) format
            image_data = np.transpose(image_data, (2, 0, 1))
        elif image_data.shape[0] <= 4:
            # Likely already (c, y, x)
            pass
        else:
            raise ValueError(
                f"Cannot determine whether TIFF is cyx or xyc. Shape: {image_data.shape}"
            )

    else:
        raise ValueError(
            f"Unsupported number of dimensions: {image_data.ndim}. "
            "Expected 2D or 3D TIFF."
        )

    metadata = {
        "scale": {"x": 1, "y": 1},
        "translation": {"x": 0, "y": 0},
        "dims": ("c", "y", "x"),
    }
    if image_data.shape[0] == 3:
        metadata["c_coords"] = ["r", "g", "b"]
    elif image_data.shape[0] == 4:
        metadata["c_coords"] = ["r", "g", "b", "a"]
        
    image = Image2DModel.parse(
        data=image_data, transformations=transformations, **metadata
    )
    sdata.images[element_name] = image

def load_rois(
    sdata: Any,
    input_folder: str,
    element_name: str,
    output_file: str,
    class_name: str,
    class_color: str,
    transformations: Dict[str, Any],
) -> None:
    """
    Convert ImageJ ROI files to GeoJSON and add them to SpatialData.

    Parameters:
        sdata (Any): SpatialData object.
        input_folder (str): Path to the folder containing ROI files.
        output_file (str): Path to the output GeoJSON file.
        class_name (str): Classification name for the ROIs.
        class_color (str): Classification color for the ROIs.
        transformations (Dict[str, Any]): Transformation dictionary for spatial coordinates.
    """
    
    if not os.path.exists(output_file):
        features = []
        print (input_folder, output_file, os.listdir(input_folder))
        for filename in tqdm(os.listdir(input_folder)):
            
            if filename.endswith(".roi"):

                roi_path = os.path.join(input_folder, filename)
                try:
                    
                    roi = ImagejRoi.fromfile(roi_path)
                    
                    if roi.roitype in (ROI_TYPE.TRACED, ROI_TYPE.POLYGON, ROI_TYPE.FREEHAND):
                        print(os.path.splitext(filename)[0])
                        coords = roi.coordinates().tolist()
                        coords.append(coords[0])
                        feature = geojson.Feature(
                            geometry=geojson.Polygon([coords]),
                            properties={
                                "name": os.path.splitext(filename)[0],
                                "classification": {
                                    "name": class_name,
                                    "color": class_color,
                                },
                            },
                        )
                        features.append(feature)
                except Exception as e:
                    
                    print(f"Error processing {filename}: {e}")
                             
        feature_collection = geojson.FeatureCollection(features)
        
        with open(output_file, "w") as f:
            geojson.dump(feature_collection, f, indent=2)

    shapes_model = ShapesModel.parse(output_file, transformations=transformations)
    
    # Manually assign the 'name' property from the GeoJSON features
    with open(output_file, "r") as f:
        geojson_data = geojson.load(f)

    names = [feature["properties"]["name"] for feature in geojson_data["features"]]
    shapes_model["name"] = names
    sdata.shapes[element_name] = shapes_model


def load_points(
    sdata: Any, input_folder: str, element_name: str, transformations: Dict[str, Any]
) -> None:
    """
    Load points from CSV files and add them to SpatialData.

    Parameters:
        sdata (Any): SpatialData object.
        input_folder (str): Path to the folder containing CSV files.
        element_name (str): Name of the points element in SpatialData.
        transformations (Dict[str, Any]): Transformation dictionary for spatial coordinates.
    """
    dfs = []

    for filename in os.listdir(input_folder):
        if filename.endswith(".csv"):
            df = pd.read_csv(os.path.join(input_folder, filename))
            df["target"] = Path(filename).stem
            dfs.append(df)

    spots = pd.concat(dfs, ignore_index=True)

    points = PointsModel.parse(
        spots, coordinates={"x": "X", "y": "Y"}, transformations=transformations
    )
    sdata.points[element_name] = points

def load_numpy_grid(
    sdata: Any, obs_array: np.ndarray, obs_names: list, element_name: str, transformations: Dict[str, Any]
) -> None:
    """
    Load 3D numpy array containing channels and add it to SpatialData as a table and image.

    Parameters:
        sdata (Any): SpatialData object.
        obs_array (np.ndarray): Array containing 3D numpy array, where the first dimension is the number of channels.
        obs_names (list): List of channel names.
        element_name (str): Name of the element in SpatialData.
        transformations (Dict[str, Any]): Transformation dictionary for spatial coordinates.
    """
    grid = []
    x_cell_size, y_cell_size = 1, 1
    for x0 in np.arange(0, obs_array.shape[-1], x_cell_size):
        for y0 in np.arange(0, obs_array.shape[-2], y_cell_size):
            cell = Polygon(
                [
                    (x0, y0),
                    (x0 + x_cell_size, y0),
                    (x0 + x_cell_size, y0 + y_cell_size),
                    (x0, y0 + y_cell_size),
                ]
            )
            grid.append(cell)

    grid_gdf = GeoDataFrame(grid, columns=["geometry"])
    grid_spatial = np.array([list(i.coords[0]) for i in grid_gdf.centroid])

    shapes_model = ShapesModel.parse(
        grid_spatial,
        geometry=0,
        radius=x_cell_size / 2,
        transformations=transformations,
    )
    sdata.shapes[f"{element_name}_grid"] = shapes_model
    adata_X = pd.DataFrame()
    
    for obs_name, obs_channel in zip(obs_names, obs_array):
        adata_X[obs_name] = obs_channel.T.flatten()

    adata = anndata.AnnData(X=adata_X)
    adata.obs["region_key"] = f"{element_name}_grid"
    adata.obs["Index"] = adata.obs.index.astype(int)
    adata.obsm["spatial"] = grid_spatial

    sdata.tables[f"{element_name}_obs"] = TableModel.parse(
        adata, region=f"{element_name}_grid", region_key="region_key", instance_key="Index"
    )

    # Perform PCA on channels to convert N channels to RGB channels
    image_data = obs_array
    
    from sklearn.decomposition import PCA
    pca = PCA(n_components=3)
    # (n_channels, y, x) -> (n_channels, y*x)
    image_data_pca = pca.fit_transform(image_data.reshape(image_data.shape[0], -1).T).T
    image_data_pca = image_data_pca.reshape(3, image_data.shape[1], image_data.shape[2])
    # Convert to RGB uint8
    image_data_pca = (image_data_pca - image_data_pca.min()) / (image_data_pca.max() - image_data_pca.min()) * 255
    image_data_pca = image_data_pca.astype(np.uint8)
    metadata = {
        "scale": {"x": 1, "y": 1},
        "translation": {"x": 0, "y": 0},
        "dims": ("c", "y", "x"),
        "c_coords": ["r", "g", "b"],
    }

    image = Image2DModel.parse(
        image_data_pca, transformations=transformations, **metadata
    )

    sdata.images[f"{element_name}_image"] = image

def load_DESI(
    sdata: Any, input_folder: str, element_name: str, transformations: Dict[str, Any]
) -> None:
    """
    Load DESI data and add it to SpatialData.

    Parameters:
        sdata (Any): SpatialData object.
        input_folder (str): Path to the folder containing DESI text files.
        element_name (str): Name of the DESI element in SpatialData.
        transformations (Dict[str, Any]): Transformation dictionary for spatial coordinates.
    """
    image_data = []
    obs_names = []
    for filename in os.listdir(input_folder):
        if filename.endswith(".txt"):
            obs_name = filename.replace(input_folder, "").replace(".txt", "")
            obs_name = "".join(e if e.isalnum() or e == "-" else "_" for e in obs_name)
            image_data.append(
                pd.read_csv(os.path.join(input_folder, filename), sep="\t", header=None).to_numpy()
            )
            obs_names.append(obs_name)
    image_data = np.array(image_data)
    #image_data = np.where(image_data == 0, np.nan, image_data) # Set 0 to NA
    load_numpy_grid(sdata, image_data, obs_names, element_name, transformations)

def merge(
    sdata_list
):
    sdata = concatenate(
        sdata_list
    )
    # Loop through all _DESI_obs tables, and _cell_types tables
    # Merge them into a single anndata object

    # Known coordinate system IDs
    ids = sorted(sdata.coordinate_systems, key=len, reverse=True)  # longest first

    suffix_groups = defaultdict(list)

    for name, adata in sdata.tables.items():
        matched_id = None
        
        for id_ in ids:
            prefix = id_ + "_"
            if name.startswith(prefix):
                matched_id = id_
                suffix = name[len(prefix):]
                suffix_groups[suffix].append(adata)
                break

        if matched_id is None:
            # Table that matches none of the IDs
            # Could be logged or ignored
            pass

    for suffix, adatas in suffix_groups.items():
        print (suffix, len(adatas))
        if len(adatas) == 1:
            continue
        merged_tables = anndata.concat(adatas, join="outer", label="batch")
        regions = merged_tables.obs["region_key"].unique().tolist()
        sdata.tables[suffix] = TableModel.parse(merged_tables, region=regions, region_key="region_key", instance_key="Index")

        for id_ in ids:
            table_name = f"{id_}_{suffix}"
            if table_name in sdata.tables:
                del sdata.tables[table_name]
    return sdata