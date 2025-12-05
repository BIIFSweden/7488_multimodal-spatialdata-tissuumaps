"""
build_sdata_from_yaml.py

Utility that reads the YAML specification produced for the
“DESI + SCRINSHOT” project and builds a ready-to-use SpatialData
object.

Dependencies
------------
pip install pyyaml spatialdata napari-spatialdata  # plus your helpers
"""

from __future__ import annotations
import os
from pathlib import Path
import re
import yaml
import pandas as pd
import numpy as np
import json

# SpatialData & transformations
from spatialdata import SpatialData, get_extent
from spatialdata.transformations import (
    Affine,
    Identity,
    Scale,
    get_transformation,
    set_transformation,
    Sequence,
)

# Helper modules that already exist in the repository
import spatialdata_helper as sd_h
import tissuumaps_helper as tm_h          # not used here, but handy later


# ---------------------------------------------------------------------
# internal helpers
# ---------------------------------------------------------------------
def _resolve(path_like: str | Path, root: Path) -> Path:
    """
    Expand "${root}" and return an absolute Path.
    """
    path_like = str(path_like).replace("${root}", str(root))
    path = Path(path_like)
    return path if path.is_absolute() else (root / path).resolve()


# ---------------------------  TRANSFORMS  -----------------------------
def _build_affine_from_csv(csv_file: Path) -> Affine:
    """
    Implements exactly the R-style parsing snippet you provided.
    """
    import numpy as np
    import pandas as pd

    df = pd.read_csv(csv_file, header=None, dtype=str)
    df = df.applymap(lambda x: re.sub(r"Matrix:|\[|\]", "", x))
    mat = df.astype(float).to_numpy()
    mat = mat.reshape((2, 3), order="F")
    mat = np.vstack([mat, [0, 0, 1]])
    return Affine(mat, input_axes=("x", "y"), output_axes=("x", "y"))

# ---------------------------------------------------------------------
# 1-bis.  *Extremely* simple transform parser
# ---------------------------------------------------------------------
def _parse_all_transforms(t_yaml: dict, root: Path):
    """
    The YAML contains only one primitive transform declared as

        transformations:
          trans_DESI_to_SCRINSHOT:
            type: affine_csv
            file: ".../something.csv"

    Therefore we do **not** distinguish between “primitive transforms”
    and “transformation-sets” any more.  We simply read the single entry
    and wrap it in a tiny dictionary whose only key is the literal
    string "default".  This keeps the public API identical to the old
    version, i.e. every loader still receives a *dictionary* that maps
    coordinate-system names to the actual SpatialData transform object.

    Returns
    -------
    transform_sets : dict[str, dict[str, BaseTransform]]
        e.g.   { "trans_DESI_to_SCRINSHOT": { "default": Affine(...) } }
    """
    transform_sets: dict[str, dict[str, object]] = {}

    for name, spec in t_yaml.items():
        ttype = spec["type"].lower()

        if ttype == "identity":
            tr = Identity()

        elif ttype == "affine_csv":
            csv_file = _resolve(spec["file"], root)
            tr = _build_affine_from_csv(csv_file)

        elif ttype == "inversed_affine_csv":
            csv_file = _resolve(spec["file"], root)
            tr = _build_affine_from_csv(csv_file).inverse()

        elif ttype == "affine":
            mat = np.asarray(spec["matrix"], dtype=float)
            tr = Affine(mat, input_axes=("x", "y"), output_axes=("x", "y"))

        else:
            raise ValueError(f"Unsupported transformation type: {ttype}")

        # wrap the single transform in a dict so that the downstream code
        # still gets a mapping  { coord_system : transform }
        transform_sets[name] = {"global": tr}

    return transform_sets

# ---------------------------------------------------------------------
# main builder
# ---------------------------------------------------------------------
def build_spatialdata_from_yaml(
        yaml_file: str | Path,
        output_dir: str | Path = None
    ) -> tuple[str, SpatialData]:
    """
    Parameters
    ----------
    yaml_file
        Path to the YAML specification.

    Returns
    -------
    sdata : SpatialData
        Fully populated object ready for visualisation or export.
    """
    yaml_file = Path(yaml_file).resolve()
    cfg = yaml.safe_load(open(yaml_file, "r"))
    sample_id = cfg.get("sample_id", ".")          # string as given in YAML

    # ------------------------------------------------------------------
    # 0.  basic paths & folders
    # ------------------------------------------------------------------
    root_spec = cfg.get("root", ".")          # string as given in YAML
    root_path = Path(root_spec)

    # If the path is NOT absolute, prepend the YAML folder
    if not root_path.is_absolute():
        root_path = (yaml_file.parent / root_path)

    root = root_path.resolve()
    print("root:", root)
    output_dir = Path(output_dir) if output_dir is not None else (yaml_file.parent / "output").resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # 1.  build transformation objects
    # ------------------------------------------------------------------
    transform_sets = _parse_all_transforms(cfg["transformations"], root)
    print("Defined transformation-sets:", transform_sets)
    # helper to fetch transformation-set by name
    def get_tset(name: str):
        """
        Retrieve the transformation-set requested by an image / shape / ….
        If the name does not exist in `transform_sets`, fall back to an
        Identity() for every coordinate system we encounter later.
        """
        if name in transform_sets:
            return transform_sets[name]

        # -------  fallback  ------------------------------------------------
        # We return a dictionary that maps to Identity for any CS key.
        # Using defaultdict lets us stay flexible – any missing key that
        # napari-spatialdata asks for will still return Identity().
        from collections import defaultdict
        return defaultdict(lambda: Identity())

    # ------------------------------------------------------------------
    # 2.  create empty SpatialData object
    # ------------------------------------------------------------------
    sdata = SpatialData()

    # ------------------------------------------------------------------
    # 3.  IMAGES
    # ------------------------------------------------------------------
    for img in cfg.get("images", []):
        t_set = get_tset(img["transformations"])
        path = _resolve(img["file"], root)
        sd_h.load_image(sdata, path, sample_id + "_" + img["id"], transformations=t_set)

    # ------------------------------------------------------------------
    # 4.  SHAPES
    # ------------------------------------------------------------------
    for shp in cfg.get("shapes", []):
        folder = _resolve(shp["folder"], root)
        t_set = get_tset(shp["transformations"])
        out_file = output_dir / f"{sample_id}_{shp['id']}.geojson"
        sd_h.load_rois(
            sdata,
            folder,
            sample_id + "_" + shp["id"],
            output_file=str(out_file),
            class_name=shp.get("class_name", shp["id"]),
            class_color=shp.get("class_color", [255, 255, 255]),
            transformations=t_set,
        )

    # ------------------------------------------------------------------
    # 5.  POINTS
    # ------------------------------------------------------------------
    for pts in cfg.get("points", []):
        folder = _resolve(pts["folder"], root)
        t_set = get_tset(pts["transformations"])
        sd_h.load_points(
            sdata, 
            folder, 
            sample_id + "_" + pts["id"], 
            transformations=t_set
        )

    # ------------------------------------------------------------------
    # 6.  TABLES
    # ------------------------------------------------------------------
    for tbl in cfg.get("tables", []):
        if "file" in tbl:
            csv_file = _resolve(tbl["file"], root)
            t_set = get_tset(tbl["transformations"])
            sd_h.load_table(
                sdata, 
                csv_file, 
                element_name=sample_id + "_" + tbl["id"],
                transformations=t_set
            )
        # tables with "generated_by" are created automatically by load_DESI

    # ------------------------------------------------------------------
    # 7.  DESI block (mass-spec specific)
    # ------------------------------------------------------------------
    if "DESI" in cfg:
        ds = cfg["DESI"]
        t_set = get_tset(ds["transformations"])
        print ("DESI transform-set:", t_set)
        in_path = _resolve(ds["input_path"], root)
        sd_h.load_DESI(
            sdata,
            str(in_path),
            element_name=sample_id + "_" + ds.get("id", "DESI"),
            transformations=t_set,
        )

    # Get full sample width from sdata
    sdata_extent = get_extent(sdata, coordinate_system="global")
    full_width = sdata_extent['x'][1] - sdata_extent['x'][0]
    print (f"Sample '{sample_id}' full width: {full_width:.1f} units")
    if ("min_width" in cfg) and (full_width < cfg["min_width"]):
        scale_factor = cfg["min_width"] / full_width
        print (f"Scaling up sample '{sample_id}' by factor {scale_factor:.3f} to reach min_width={cfg['min_width']}")
        scale_transform = Scale(
            scale=[scale_factor, scale_factor],
            axes=("x", "y")
        )
        for key in (
            list(sdata.images)
            + list(sdata.labels)
            + list(sdata.shapes)
            + list(sdata.points)
        ):
            old = get_transformation(sdata[key])
            set_transformation(sdata[key], Sequence([old, scale_transform]), to_coordinate_system="global")
        sdata_extent = get_extent(sdata, coordinate_system="global")
        full_width = sdata_extent['x'][1] - sdata_extent['x'][0]
        print (f"Sample '{sample_id}' full width after rescaling: {full_width:.1f} units")
        
    sdata.rename_coordinate_systems({'global': sample_id})
    return sample_id, sdata

# ---------------------------------------------------------------------
# CLI helper (updated)
# ---------------------------------------------------------------------
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Build SpatialData objects from one or more YAML descriptions."
    )
    parser.add_argument(
        "--yaml",
        type=str,
        nargs="+",
        required=True,
        help="One or more YAML configuration files.",
    )
    parser.add_argument(
        "--out",
        type=str,
        default=None,
        help="Output path for the resulting SpatialData object "
             "(if --merge is used, output is for the merged object). "
             "If --merge is NOT used, each sdata is written next to its YAML file.",
    )
    parser.add_argument(
        "--merge",
        action="store_true",
        help="Merge all created SpatialData objects into a single one.",
    )
    parser.add_argument(
        "--col",
        type=int,
        default=None,
        help="Number of columns for arranging coordinate systems in TissUUmaps.",
    )

    args = parser.parse_args()

    # --------------------------------------------------------------
    # 1. Build an sdata for each YAML file
    # --------------------------------------------------------------
    sdata_list = []
    ids = []
    for yml in args.yaml:
        print(f"\nProcessing {yml} ...")
        id, sdata = build_spatialdata_from_yaml(yml, output_dir=args.out)
        print (sdata)
        sdata_list.append(sdata)
        ids.append(id)
    print (ids)
    # --------------------------------------------------------------
    # 2. If merge requested → merge all sdata objects
    # --------------------------------------------------------------
    print("\nMerging all SpatialData objects...")
    merged = sd_h.merge(sdata_list)
    print (merged)
    os.makedirs(args.out, exist_ok=True)
    
    merged.write(args.out + "/merged_spatialdata.zarr", overwrite=True)
    print(f"Merged SpatialData written to: {args.out}")
    
    sdata_extended, tmap_file = tm_h.spatialdata_to_tissuumaps(
        merged, 
        elements = {
            "images": merged.images,
            "shapes": merged.shapes,
            "tables": merged.tables,
            "points": merged.points,
        },
        output_folder = args.out,
        coordinate_systems = ids,
        columns = args.col if hasattr(args, "col") else None,
    )
    print (tmap_file)
    # Post-process the TissUUmaps file:
    with open(tmap_file, "r") as f:
        tmap_state = json.load(f)

    #tmap_state["markerFiles"][0]["expectedHeader"]["scale_factor"] = 2
    #tmap_state["markerFiles"][1]["expectedHeader"]["scale_factor"] = 0.4
    #tmap_state["markerFiles"][2]["expectedHeader"]["scale_factor"] = 0.1
    for markerFile in tmap_state["markerFiles"]:
        for dropdownOption in markerFile["dropdownOptions"]:
            dropdownOption["name"] = dropdownOption["name"].replace("Gene expression: ", "")

    # Make the opacity column the same as the color column for the first marker file
    for markerFile in tmap_state["markerFiles"]:
        if not "DESI" in markerFile["name"]:
            continue
        for dropdownOption in markerFile["dropdownOptions"]:
            if not "expectedHeader.cb_col" in dropdownOption:
                continue
            if dropdownOption["expectedHeader.cb_col"][:2] == "/X":
                dropdownOption["expectedHeader.opacity_col"] = dropdownOption["expectedHeader.cb_col"]
                dropdownOption["expectedRadios.opacity_check"] = True
            else:
                dropdownOption["expectedHeader.opacity_col"] = ""
                dropdownOption["expectedRadios.opacity_check"] = False

    tmap_state["plugins"] = ["addLabels"]
    
    with open(tmap_file, "w") as f:
        json.dump(tmap_state, f, indent=2)

    print("\nDone.")
