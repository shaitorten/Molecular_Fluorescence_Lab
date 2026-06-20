# Fluorescence Spectroscopy and Image Analysis

This repository contains a Python-based physical chemistry/optics data analysis tool designed to process and analyze the emission spectra of organic fluorophores (Fluorescein, Rhodamine 6G, and Rhodamine B). 

The project is divided into two primary analytical components:
1. **Spectrometer Data Analysis:** Validating the Beer-Lambert law at the low-concentration linear regime and quantitatively identifying the onset of non-radiative decay mechanisms (fluorescence quenching) at higher concentrations.
2. **Photographic Image Analysis:** Extracting and quantifying fluorescence intensity directly from smartphone images.

---

## Key Features

### Part 1: Spectrometer Data Analysis
* **Automated Data Processing:** Batch reads and parses spectrometer data from `.ods` (OpenDocument Spreadsheet) files.
* **Numerical Integration:** Calculates the total integrated fluorescence intensity across a specified wavelength range using the trapezoidal rule.
* **Rigorous Error Propagation:** Computes absolute measurement uncertainties by combining statistical Poisson noise (shot noise) with instrumental wavelength resolution errors ($\Delta\lambda$).
* **Weighted Linear Fitting:** Performs Weighted Least Squares (WLS) regression on the low-concentration regime to extract the photophysical proportionality constant ($k\epsilon_{\lambda}$), taking into account the calculated error bars.
* **Data Visualization:** Generates publication-ready dual plots showing both the full concentration range and the isolated linear fit for the first few concentration data points.

### Part 2: Photographic Image Analysis
* **Smartphone Data Extraction:** Processes photographs (e.g., captured via a smartphone) to extract quantitative fluorescence intensity values directly from pixel data.
* **Spatial Intensity Profiling:** Maps pixel intensity across the specified sample area (e.g., along the cuvette) to evaluate fluorescence distribution.
* **Linear Fitting:** Conducts a best-fit linear model to Intensity Vs Optical Length.

---

## Requirements

To run the scripts, ensure you have the following Python libraries installed:

* `pandas` (with the `odfpy` engine for reading `.ods` files)
* `numpy`
* `scipy`
* `matplotlib`
* `opencv-python` or `Pillow` (for the image processing module)

---

## Usage

### Part 1: Spectrometer Analysis
1. Organize your raw spectrometer data (`.ods` format) into the working directory. Ensure the filenames contain the concentration values (e.g., `6G_0.05mm.ods`).
2. Update the `materials_files` dictionary in the script to match your local file names.
3. Configure the physical parameters:
   * `I0`: Incident LED intensity (derived from the reference integration).
   * `l_path`: Optical path length of the cuvette.
   * `delta_lambda`: Wavelength accuracy/resolution of your spectrometer.
4. Run the script to output the calculated linear slopes, $k\epsilon_{\lambda}$ constants, and analytical plots.

### Part 2: Image Analysis
1. Place your raw image files into the designated image directory.
2. Define the Region of Interest (ROI) coordinates in the script to crop the image to the relevant fluorescent area (e.g., the illuminated section of the cuvette).
3. Run the image processing module to output the integrated pixel intensities and generate spatial intensity profile plots.

```
