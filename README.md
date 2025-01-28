# Spatial mapping of oxylipins in lung tissue

Co-registration of images from (i) mass spectrometry imaging and (ii) microscope images of the same tissue section to identify regions of interest in the microscope images to annotate the mass spectrometry image.

[Download code](../../archive/refs/heads/main.zip)

## Installation

1. Clone this repository locally.
2. Install uv package from [https://docs.astral.sh/uv/getting-started/installation/](https://docs.astral.sh/uv/getting-started/installation/).

## Usage

1. In a terminal or powershell, run:
    ```
    uv run jupyter-lab
    ```
    A Jupyter-lab tab should open in your browser. Navigate to the notebook `notebooks/create_spatialdata.ipynb`.
    Change input folder depending on the data location.
2. Run the notebook.
3. Test TissUUmaps output in [TissUUmaps](https://tissuumaps.github.io/installation/).

## License

[MIT](LICENSE)

## Contact

[SciLifeLab BioImage Informatics Facility (BIIF)](https://www.scilifelab.se/units/bioimage-informatics/)

Developed by [Christophe Avenel](mailto:christophe.avenel@scilifelab.se)
