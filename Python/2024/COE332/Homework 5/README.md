# REST API Analysis of ISS Trajectory Data 

## High-Level Description
This folder contains a Dockerfile, allowing the user to create and run an image of a Flask app analyzing ISS tracking data. The code consists of a main script "iss_tracker" hosting the web application functions returning current pieces of the ISS data taken online, and a unit test script "test_iss_tracker" that tests the main script functions. 

## Software Diagram
![Alt text](https://github.com/AaronPandian/coe323-homeworks/blob/main/homework05/diagram.png)

## How to Access the Data
The ISS tracking data this code requests can be found on the NASA website:
"https://spotthestation.nasa.gov/trajectory_data.cfm." This ephemeris data, compiled by the NASA Johnson Space Center, contains a header section and a primary data section. The header primarily, concerning the utility of this code, contains the ISS mass in kg, drag area in m^2, and drag coefficient used in generating the subsequent data. The subsequent section contains data from the last 15-day interval. The timesteps vary from 4 minutes to 2 seconds and each house state vectors containing the time in UTC; position X, Y, and Z in km; and velocity X, Y, and Z in km/s.

## How to Build the Container
Clone the GitHub repository to your machine, then log in to docker from your machine. 

Implement "docker build" to construct an image of the container using the path of the Dockerfile from the source code. It will look something like `docker build -t <dockerhubusername>/<code>:<version> .` to build the Dockerfile. The `<dockerhubusername>` above represents the image and tag name on a local machine. The `<code>` represents the filename, for example "iss_tracker.py".

Subsequently, use `docker tag` to set a tag for the image. Next run `docker images` to ensure the instance was created successfully with the corresponding tag.

## How to Deploy Containerized Code as a Flask App
Now using the following command, we can run the main script iss_tracker.py as `docker run --name "iss-tracker-app" -d -p 5000:5000 <dockerhubusername>/<code>:<version>` and obtain an image instance. 

To run the unit test, use the same `docker run <dockerhubusername>/<code>:<version>` command replacing the code file with "test_iss_tracker.py" to test the main script functions. 

## Accessing Routes
Once the image is running, the terminal will be waiting for requests to be made using specific URL routes. Using the HTTPS URL displayed in the terminal, open another window in your command prompt and paste the URL. Then append the following routes at the end of the URL to obtain the desired functions. 

* `/epochs` returns the whole dataset.
* `/epochs?limit=int&offset=int` returns a modified dataset, given the integer query parameters. The offset parameter sets the start and the limit sets the final index for which the data is returned, thus the limit must be greater than the offset to query successfully.
* `/epochs/<epoch>` returns the state vector for the specific epoch at index `<epoch>`.
* `/epochs/<epoch>/speed` returns the instantaneous speed of the ISS at the specified epoch index `<epoch>`.
* `/now` returns the state vector and instantaneous speed for the most current epoch.

## Output and What to Expect
In running the main script from an image, once running the routes above, the user should receive the respective information printed out to the terminal. 
