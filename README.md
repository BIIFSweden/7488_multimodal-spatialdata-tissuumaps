# Spatial mapping of oxylipins in lung tissue

Combine multiple imaging modalities of the same tissue section — DESI mass spectrometry imaging, microscopy images (e.g. H&E), SCRINSHOT data and ImageJ ROI annotations — into a single [SpatialData](https://spatialdata.scverse.org/) object, applying affine transformations to align the modalities in a shared coordinate system. The resulting object is then exported as a [TissUUmaps](https://tissuumaps.github.io/) project for interactive visualization.

Inputs are described in a YAML configuration file (see [notebooks/params_examples/](notebooks/params_examples/)), and multiple samples can be merged into a single SpatialData object for joint visualization.

[Download code](../../archive/refs/heads/main.zip)

## Installation

1. Clone this repository locally.
2. Install [pixi](https://pixi.sh/latest/#installation).

> **Note:** The pixi environment is currently defined for Windows (win-64) only.
> Reading whole-slide images additionally requires [OpenSlide](https://openslide.org/download/) installed in `C:\openslide`.

## Usage

### Batch processing (recommended)

1. Create a YAML configuration file describing your input data (images, ROIs, transformations, DESI data). Examples are provided in [notebooks/params_examples/](notebooks/params_examples/).
2. In a terminal or powershell, run:
    ```
    pixi run create_spatialdata
    ```

### Notebooks

1. In a terminal or powershell, run:
    ```
    pixi run jupyter lab
    ```
    A Jupyter-lab tab should open in your browser. In Jupyter-lab, navigate to the notebook [notebooks/create_spatialdata.ipynb](notebooks/create_spatialdata.ipynb).
    Change input and output folders depending on the data location on your disk.
2. Run the notebook. Use [notebooks/merge_sdata.ipynb](notebooks/merge_sdata.ipynb) to merge multiple SpatialData objects into one.
3. Test TissUUmaps output in [TissUUmaps](https://tissuumaps.github.io/installation/).

## License

[MIT](LICENSE)

## Contact

[SciLifeLab BioImage Informatics Facility (BIIF)](https://www.scilifelab.se/units/bioimage-informatics/)

Developed by [Christophe Avenel](mailto:christophe.avenel@scilifelab.se)
