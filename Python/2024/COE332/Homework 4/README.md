# Capture and Analysis of ISS Trajectory Data 

## High-Level Description
This folder contains a Dockerfile, allowing the user to create and run an image of the ISS tracking and analysis code. The code consists of a main script "iss_tracker" printing current summaries of the ISS data taken online, and a unit test script "test_iss_tracker" that tests the main script functions. 

## How to Access the Data
The ISS tracking data this code requests can be found on the NASA website:
"https://spotthestation.nasa.gov/trajectory_data.cfm." This ephemeris data, compiled by the NASA Johnson Space Center, contains a header section and a primary data section. The header primarily, concerning the utility of this code, contains the ISS mass in kg, drag area in m^2, and drag coefficient used in generating the subsequent data. The subsequent section contains data from the last 15-day interval. The timesteps vary from 4 minutes to 2 seconds and each house state vectors containing the time in UTC; position X, Y, and Z in km; and velocity X, Y, and Z in km/s.

## How to Build the Container
Clone the GitHub repository to your machine, then log in to docker from your machine. 

Implement "docker build" to construct an image of the container using the path of the Dockerfile from the source code. It will look something like `docker build -t <dockerhubusername>/<code>:<version> .` to run the Dockerfile. 

Subsequently, use `docker tag` to set a tag for the image. Next run "docker images" to ensure the instance was created successfully with the corresponding tag.

## How to Run the Containerized Scripts
Now using the following command, we can run the main script iss_tracker.py as `docker run username/iss_tracker:1.0 iss_tracker.py`. The `username/` above represents the image and tag name on a local machine.

To run the unit test, use the same `docker run` command replacing the filename with "test_iss_tracker.py" to test the main script functions. 

## Output and What to Expect
In running the main script from an image, the user should receive a summary printed out to the terminal. This summary includes logging statements relaying the code is running at various sections of the script. In between these milestone logs, the function will output three output statements- one for each function present in the main script. The first produces the time range in ISO format of the dataset at runtime, the second prints the final recorded state values, and the final function states the average speed across the dataset as well as the instantaneous speed at the final epoch. 
