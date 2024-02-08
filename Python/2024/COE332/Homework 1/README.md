Unit Testing and Error Logging Statistic Functions

Description:
This folder contins four scripts of code- two of which outline functions that derrive summary statistics of a dataset, specifically for meteorite landings on Earth. The objective of this project is to develop functions capable of detailing to the reader relavent overview data, as well as getting the author accustomed to logging, unit testing and error handling of code.

The ml_data_analysis is the main file where summary statistics are outputed. This file calls the gcd_agorithm file and imports a function calculating the great circle from gcd_algorithm. The ml_data_analysis script then reads the dataset and, using said functions, compiles a summary of the set. 

The latter two files act as unit testing scripts for the former files. Specifically, the test_gcd_algorithm tests the main script and the test_gcd_algorithm runs a unit test for the gcd_algorithm file. 

Installation:
The dataset this code was intended for can be found in the following link: https://data.nasa.gov/Space-Science/Meteorite-Landings/gh4g-9sfh/about_data. To use the files, be sure to download the following libraries: csv, logging, matplotlib.pyplot, numpy, math, and norm from scipy.states. The main file will import and call these libraries when run. 

The meteorite landing data file must be a CSV imported to the same folder the other code files are house. The file must have the name 'Meteorite_Landings' for the script to open the file appropriately. 

Usage:
Make sure to import all libraries as depcited in the main script (ml_data_analysis). After downloading all files into the same folder and aquiring the dataset file in the same folder as well, run the main script to see the analysis. Run the seperate unit testing files to verify functions are working as should. 
