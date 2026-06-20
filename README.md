# Molecular_Fluorescence_Lab
Python-based data analysis tool for molecular fluorescence spectroscopy

# Fluorescence Spectroscopy Data Analysis

This repository contains a Python-based physical chemistry/optics data analysis tool designed to process and analyze the emission spectra of organic fluorophores (Fluorescein, Rhodamine 6G, and Rhodamine B). 

The primary goal of this project is to validate the Beer-Lambert law at the low-concentration linear regime and quantitatively identify the onset of non-radiative decay mechanisms, such as fluorescence quenching, at higher concentrations.

## Key Features

* **Automated Data Processing:** Batch reads and parses spectrometer data from `.ods` (OpenDocument Spreadsheet) files.
* **Numerical Integration:** Calculates the total integrated fluorescence intensity across a specified wavelength range using the trapezoidal rule.
* **Rigorous Error Propagation:** Computes absolute measurement uncertainties by combining statistical Poisson noise (shot noise) with instrumental wavelength resolution errors ($\Delta\lambda$).
* **Weighted Linear Fitting:** Performs Weighted Least Squares (WLS) regression on the low-concentration regime to extract the photophysical proportionality constant ($k\epsilon_{\lambda}$), taking into account the calculated error bars.
* **Data Visualization:** Generates publication-ready dual plots showing both the full concentration range (highlighting non-linear quenching effects) and the isolated linear fit for the first few concentration data points.

## Requirements

To run the script, ensure you have the following Python libraries installed:

* `pandas` (with `odfpy` engine for reading `.ods` files)
* `numpy`
* `scipy`
* `matplotlib`

## Usage

1. Organize your raw spectrometer data (`.ods` format) into the working directory. Ensure the filenames contain the concentration values (e.g., `6G_0.05mm.ods`).
2. Update the `materials_files` dictionary in the script to match your local file names.
3. Configure the physical parameters:
   * `I0`: Incident LED intensity (derived from the reference integration).
   * `l_path`: Optical path length of the cuvette (typically 1.0 cm).
   * `delta_lambda`: Wavelength accuracy/resolution of your spectrometer.
4. Run the script. It will output the calculated linear slopes, the $k\epsilon_{\lambda}$ products with their respective uncertainties, and render the analytical plots.
