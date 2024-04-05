Unit Testing and Error Logging Statistic Functions

Description:
This folder contains four scripts of code- two of which outline functions that derive summary statistics of a dataset, 
specifically for meteorite landings on Earth. The objective of this project is to develop functions capable of detailing 
to the reader relevant overview data, as well as getting the author accustomed to logging, unit testing, and error handling of code. 

The ml_data_analysis is the main file where summary statistics are outputted. This file calls the gcd_agorithm file and imports a 
function calculating the great circle from gcd_algorithm. The ml_data_analysis script then reads the dataset and, using said functions, 
compiles a summary of the set. 

The latter two files act as unit testing scripts for the former files. Specifically, the test_gcd_algorithm tests the main script and 
the test_gcd_algorithm runs a unit test for the gcd_algorithm file. 

Installation:
Please clone this repository to obtain everything needed for creating your docker image of the container. 

Upon installation, ensure the dependencies inside the Dockerfile include downloading the: csv library, logging library, matplotlib.pyplot 
function, numpy library, math library, and the norm function from the scipy.states module. 

Furthermore, ensure the dataset is a CSV.

![Alt text](https://github.com/AaronPandian/coe323-homeworks/blob/main/homework03/diagram.png)
Image Details: 
The software diagram above represents the sequential process for the creation of your container, housing the code in this repository. As evident, the user connects to a virtual machine that then pulls information from the web for the container docker file and dataset. These files are then set up in the virtual machine terminal; the separate scripts are then run for two purposes: to conduct data analysis on the dataset and to test the scripts doing said data analysis. 

Usage:
Initialize
Clone the GitHub repository to your virtual machine. Then move the source code to a serperate folder "code/" inside a root folder. 

How to build the image from a Dockerfile
Implement "docker build" to construct an image of the container using the path of the Dockerfile from the moved source code. Then use "docker 
tag" to set a tag for the image. Then run "docker images" to ensure the instance was created successfully with the corresponding tag.  

How to get the input data from the web
The dataset this code was intended for can be found in the following website: https://data.nasa.gov/Space-Science/Meteorite-Landings/gh4g-9sfh/about_data. 
For the above dataset, use "wget https://raw.githubusercontent.com/tacc/coe-332-sp24/main/docs/unit02/sample-data/Meteorite_Landings.csv" and move the imported "Meteorite_Landings.csv" to a folder "data/" also within the root folder. 

How to mount the data inside the container at run time
Now using the following command, we can mount our data into the container and run the main scipt ml_data_analysis.py concurrently. 
"docker run --rm -v $PWD/Meteorite_Landings.csv:/data/Meteorite_Landings.csv \
       username/ml_data_analysis:1.0 \
       ml_data_analysis.py /data/Meteorite_Landings.csv"
The "username/" above represents the image and tag name on local machine. The command above uses the "-v" command to enable data mounting. 

How to run the containerized code
The scipt above manages to run the ml_data_analysis.py and gcd_algorithm.py scripts with the provided data to see the dataset analysis. 

How to run the containerized unit tests
To run the unit test files, no input data is necessary. As a result, use the "docker run command to run the test_ml_data_analysis.py and 
test_gcd_algorithm.py scripts without mounting the dataset. Run the separate unit testing files to verify functions are working as should. 
