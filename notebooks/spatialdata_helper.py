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
from tqdm import tqdm
import json

def load_table(
    sdata: Any, 
    input_file: str, 
    transformations: Dict[str, Any]
) -> None:
    """
    Load cell types data into an anndata table and add it to SpatialData.

    Parameters:
        sdata (Any): SpatialData object.
        input_file (str): Path to the CSV file containing cell type data.
        transformations (Dict[str, Any]): Transformation dictionary for spatial coordinates.
    """
    cell_types = pd.read_csv(input_file)

    obsm = {"spatial": cell_types[["X", "Y"]].to_numpy()}
    obs = cell_types[["roi", "Cell.type", "Number", "Area"]]
    obs["region_key"] = "cells"
    obs["Index"] = obs.index

    X = cell_types.drop(columns=["X", "Y", "roi", "Cell.type", "Number", "Area"])

    adata = anndata.AnnData(X=X, obs=obs, obsm=obsm)
    
    sdata.tables["cell_types"] = TableModel.parse(
        adata, region="cells", region_key="region_key", instance_key="Index"
    )

    circles = ShapesModel.parse(
        adata.obsm["spatial"],
        geometry=0,
        radius=40,
        index=adata.obs["Index"],
        transformations=transformations,
    )
    sdata.shapes["cells"] = circles

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

    metadata = {
        "scale": {"x": 1, "y": 1},
        "translation": {"x": 0, "y": 0},
        "dims": ("c", "y", "x"),
    }

    image = Image2DModel.parse(
        data=image_data, transformations=transformations, **metadata
    )
    sdata.images[element_name] = image

def load_rois(
    sdata: Any,
    input_folder: str,
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

        for filename in tqdm(os.listdir(input_folder)):
            if filename.endswith(".roi"):
                roi_path = os.path.join(input_folder, filename)
                try:
                    roi = ImagejRoi.fromfile(roi_path)
                    if roi.roitype in (ROI_TYPE.TRACED, ROI_TYPE.POLYGON):
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
    sdata.shapes[class_name] = shapes_model

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
    DESI_X = pd.DataFrame()
    DESI_spatial = None

    if "DESI" in sdata.points:
        del sdata.points["DESI"]
    image_data = []
    for filename in os.listdir(input_folder):
        if filename.endswith(".txt"):
            obs_name = filename[:-4].split("13Oct_P3_lung_DESI_SCRINSHOT_")[1]
            obs_name = "".join(e if e.isalnum() or e == "-" else "_" for e in obs_name)
            image_data.append(
                pd.read_csv(os.path.join(input_folder, filename), sep="\t", header=None).to_numpy()
            )

            if DESI_spatial is None:
                grid = []
                x_cell_size, y_cell_size = 1, 1
                for x0 in np.arange(0, image_data[-1].shape[1], x_cell_size):
                    for y0 in np.arange(0, image_data[-1].shape[0], y_cell_size):
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
                DESI_spatial = np.array([list(i.coords[0]) for i in grid_gdf.centroid])

                shapes_model = ShapesModel.parse(
                    DESI_spatial,
                    geometry=0,
                    radius=x_cell_size / 2,
                    transformations=transformations,
                )
                sdata.shapes["DESI"] = shapes_model

            intensity = image_data[-1].T.flatten()
            DESI_X[obs_name] = intensity

    adata = anndata.AnnData(X=DESI_X)
    adata.obs["region_key"] = "DESI"
    adata.obs["Index"] = adata.obs.index.astype(int)
    adata.obsm["spatial"] = DESI_spatial

    #adata.write_h5ad("../data/output/DESI_obs.h5ad")

    sdata.tables[f"{element_name}_obs"] = TableModel.parse(
        adata, region="DESI", region_key="region_key", instance_key="Index"
    )

    # Perform PCA on DESI data to convert N channels to 1 channel
    image_data = np.array(image_data)
    
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
