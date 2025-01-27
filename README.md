# polarization-splitting on-chip spectrometer
| Folder | Description |
| :---: | :--- |
| common | functions used by scripts in two or more other folders |
| data | raw data underlying plots in main paper or supplement |
| experiment code | code used to control equipment and generate data |
| figures | code used to generate plots from raw data for main paper figures |
| supplement-figures | code used to generate plots from raw data for supplement figures |
| utilities | code not used to generate plots for paper, but useful for visualizing raw data |

Note: in order for any MATLAB code to run, this entire repository should be added to the MATLAB path. Additionally, the repository at https://github.com/ucsd-chip-scale-photonics-lab/matlab-equipment-code should be downloaded and added to the path in order for equipment control code to work. Finally, all code should be run with this repository root as the working directory for hard-coded relative paths to data files to work in figure-generating code.
