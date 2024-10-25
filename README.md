# Readme

Datasets and `R` scripts to reproduce `Sharma et al., 2024` (under revision). Hyperspectral Sensing for High-Throughput Chloride Detection in Grapevines. The Plant Phenome Journal.
A preprint of the submitted manuscript is available at: https://www.biorxiv.org/content/10.1101/2024.06.20.599906v1.abstract.

This repository contains three sub-folders: 

1. `00_rawdata`
2. `01_codes`
3. `02_processeddata`.

**00_rawdata** contains the spectral data collected during the experiment, organized by trial.

 - `trial01` contains spectral data from two devices, the **SVC HR-1024i** and the **NIR-S-G1**. 
 
 `230606_svc_reflectance.csv` was collected with the SVC HR-1024i and is the dataset used in the provided codes. innospectra_reflectance.csv contains spectral data collected with the **NIR-S-G1**.

- `trial02` contains data collected only with the SVC instrument (`230718_svc_reflectance.csv`).
	
- Each trial folder contains a file named `chloridometer_readings.csv`, which corresponds to the ground trait data recorded with the chloridometer.
	
**01_codes** contains the `R` scripts used to process the spectral data and perform the analyses presented in the manuscript.

**02_processeddata** contains clean spectral data processed with the waves package in `R`.