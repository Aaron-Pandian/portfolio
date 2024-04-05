# Flask ISS Trajectory Analysis Service

### High-Level Description
This repository contains the source code and instructions to build, run, and test a developed containerized flask application. The service relays current information on the International Space Station (ISS) trajectory to the user. The methods of engaging with this data can be found in the various flask routes detailed below. The ephemeris data the service utilizes is from NASA's public ISS trajectory dataset [[1]](#citations). 

### Table of Contents
1. [Software Diagram](#software-diagram)
3. [Data Description](#data-description)
4. [Repository Description](#repository-description)
5. [Build and Deploy](#build-and-deploy)
   1. [How to Build the Container](#how-to-build-the-container)
   2. [How to Deploy Containerized Code as a Flask App](#how-to-deploy-containerized-code-as-a-flask-app)
   3. [How to Run Unit Tests](#how-to-run-unit-tests)
6. [Service Functionality](#service-functionality)
   1. [Accessing Routes](#accessing-routes)
   2. [What Outputs to Expect](#what-outputs-to-expect)
7. [Citations](#citations)

### Software Diagram
![Alt text](https://github.com/AaronPandian/iss-tracker-service/blob/main/diagram.png)

### Data Description
The ISS tracking data this app requests can be found on the NASA website [[1]](#citations). This ephemeris dataset, compiled by the NASA Johnson Space Center, contains a header section and a primary data section. The header contains the ISS mass in kg, drag area in m<sup>2</sup>, and drag coefficient used in generating the subsequent data. The primary data section contains information from the last 15-day interval. The timesteps vary from 4 minutes to 2 seconds and timestep notes state vectors detailing the time in UTC ISO date format; position X, Y, and Z in km; and velocity X, Y, and Z in km/s.

### Repository Description
In this repository, there are six critical files. A short description of each is bulleted in the following list. 
* _iss_tracker.py_ is the main service file running the flask application that enables a user to request certain compiled information from the ISS trajectory dataset.
* _docker-compose.yml_ is the YAML file used to state rules created to **replace** complicated docker run commands the user must input to create an image instance. 
* _Dockerfile_ is the recipe for our application installation process used by _docker-compose.yml_. 
* _requirements.txt_ is a text file managing package dependencies for the application used by _Dockerfile_. 
* _test_iss_tracker.py_ is the testing script that runs unit tests on the routes and functions developed within _iss_tracker.py_.

### Build and Deploy
First, ensure the environment you are using has Docker installed. Second, you should be conducting the following within the root folder you imported the source code into before.
#### How to Build the Container
Create a folder in your directory to input the code, for example, _iss_app_, by running `mkdir iss_app`.

Run `cd iss_app` to enter the created folder then run the `wget <linktofile>` command to import all the files from this repository into your directory. 

Once all the files are gathered, double-check with `ls`. To build the image, run the command `docker build -t <dockerhubusername>/iss_tracker:1.0 .` and check if the build was successful with `docker images`.

By implementing the command above, you should see the image created with the tag `<dockerhubusername>/iss_tracker:1.0`.

#### How to Deploy Containerized Code as a Flask App
After building the image, you can run the instance as a container. To do so, enter the command `docker-compose up -d`. This runs the Docker Compose file which deploys the image in the background.  

At this point, your container is running the main iss_tracker.py script in the background of your terminal. Use the following section to interact with the application. 

Once you are done running the Flask app, to clean up your interface, remove the image using the container ID found when running `docker images`. Once the ID is found, run `docker stop <containerID>` to stop the application from running in the background, and then `docker rm <containerID>` to remove the instance from your list of images.

#### How to Run Unit Tests
Just to note, you can run a unit test script to ensure the main iss_tracker.py script is running as it should. After the image is built and while the main script is not running, use the `docker run <dockerhubusername>/iss_tracker:1.0 test_iss_tracker.py` command to run the test. If no output is seen, then the main service script is working as it should be.  

### Service Functionality
#### Accessing Routes
After running the `docker-compose up -d` command, a background terminal will be waiting for requests to be made using specific URL routes. Using the HTTPS URL displayed in your main terminal, type `curl <URL>`, then append the following routes at the end of the URL to induce the desired dataset analysis. 

1. `/epochs` returns the whole dataset.
2. `/epochs?limit=int&offset=int` returns a modified dataset given the integer query parameters. Note that to run this route, specifically, the command will include quotation marks like `curl '<URL>/epochs?limit=int&offset=int'`.
     * The offset parameter denotes the number of epochs you want to return and the limit sets the final epoch index. If the offset value is greater than the limit, _all_ preceding epoch data points will be returned. 
3. `/epochs/<epoch>` returns the state vector information for the specific epoch index `<epoch>`.
     * The `<epoch>` entered into the route will throw an error if this value is not an integer.
4. `/epochs/<epoch>/speed` returns the instantaneous speed of the ISS at the specified epoch at index `<epoch>`.
     * The `<epoch>` entered into the route will throw an error if this value is not an integer.
5. `/epochs/<epoch>/location` returns the latitude, longitude, altitude, and geoposition- in order, at the specified epoch index `<epoch>`.
     * The `<epoch>` entered into the route will throw an error if this value is not an integer.
6. `/now` returns the instantaneous speed for the latest epoch alongside the current latitude, longitude, altitude, and geoposition- in order.
7. `/comment` returns the comments from the ISS ephemeris dataset.
8. `/header` returns the header information, as detailed prior, from the ISS ephemeris dataset.
9. `/metadata` returns the metadata from the ISS ephemeris dataset.
     * This information includes the ISS object name, ID, center name, data reference frame, time system, start time, and end time. 

#### What Outputs to Expect
In running the main script from an image, once running the routes above, the user should receive the respective information printed out to the terminal. A few example commands and their output can be seen below. Note that the dataset is updating constantly so _example outputs may not be exactly what you see_, and host HTTPS URL links vary. 

After running the `docker-compose up -d` command, the link provided was _http://10.146.109.153:5000_. Using this, the following example commands were run. 

##### `curl http://10.146.109.153:5000/metadata`
```
{
  "CENTER_NAME": [
    "EARTH"
  ],
  "OBJECT_ID": [
    "1998-067-A"
  ],
  "OBJECT_NAME": [
    "ISS"
  ],
  "REF_FRAME": [
    "EME2000"
  ],
  "START_TIME": [
    "2024-064T12:00:00.000Z"
  ],
  "STOP_TIME": [
    "2024-079T12:00:00.000Z"
  ],
  "TIME_SYSTEM": [
    "UTC"
  ]
}
```

##### `curl http://10.146.109.153:5000/epochs/10/speed`
```
[
  7.664049347556523
]
```

##### `curl http://10.146.109.153:5000/epochs/20/location`
```
[
  [
    -21.416201983874277,
    -49.391407404888184,
    429.995041495049
  ],
  "Vale Formoso, Novo Horizonte, Regio Imediata de So Jos do Rio Preto, Regio Geogrfica Intermediria de So Jos do Rio Preto, San Pablo, Regin Sudeste, Brasil"
]
```

### Citations
<a id="1">[1]</a>
“Spot the Station.” ISS Trajectory Data, NASA, https://spotthestation.nasa.gov/trajectory_data.cfm. Accessed 6 Mar. 2024. 
