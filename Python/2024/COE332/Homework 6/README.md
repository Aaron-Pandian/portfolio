# REST API Analysis of Human Genome Data using Redis 

## High-Level Description
### Project
The HUGO Gene Nomenclature Committee (HGNC) selects a unique name for each gene sequence. In this project, HGNC data is populated into a Redis database through a Flask interface to enable efficient data analysis for a user to conduct.  

### Files
This folder contains a **Dockerfile** and **requirements.txt** file, which holds library dependencies of the code, allowing the user to create and run an API image. Furthermore, the **docker-compose.yaml** file provides a swift method to build said image. The code consists of a main script **genome_api.py** hosting the web application functions- returning analytical information from the HGNC dataset online. 

## How to Access the Data
The genome data this code requests can be found on the [HGNC website](https://spotthestation.nasa.gov/trajectory_data.cfm). The construction of the data can be found on the website. The specific JSON genome data referenced can be found [here](https://g-a8b222.dd271.03c0.data.globus.org/pub/databases/genenames/hgnc/json/hgnc_complete_set.json). The format of said data starts with a basic header including the number of genomes listed. After this, each entry in the JSON dictionary represents a new genome indicated with a name. Each genome lists different characteristics from the parameters indicated on the website, like the **hgnc_id**, **symbol**, and **locus_group**.

## How to Build the Container
First, ensure the environment you are using has Docker installed and ensure the working directory has given Redis permissions.

### How to Build the Container
Create a folder in your directory to input the code found in this folder by running `mkdir genome_app`. In this directory run `mkdir data`.

Run `cd genome_app` to enter the created folder then run the `wget <linktofile>` command to import all the files from this repository into your directory. 

Once all the files are gathered, double-check with `ls`. To build the image, run the command `docker build -t <dockerhubusername>/genome_api:1.0 .` and check if the build was successful with `docker images` or `docker ps -a`.

By implementing the command above, you should see the image created with the tag `<dockerhubusername>/genome_api:1.0`.

### How to Deploy Containerized Code as a Flask App
After building the image, you can run the instance as a container. To do so, enter the command `docker-compose up -d --build`. This runs the Docker Compose file which deploys the image in the background. 

At this point, your container is running the main **genome_api.py** script in the terminal background. Use the following section to interact with the application. 

Once you are done running the Flask app remove the image using the container ID found when running `docker images`. Once the ID is found, run `docker stop <containerID>` to stop the application from running in the background, and then `docker rm <containerID>` to remove the instance from your list of images.

## Accessing Routes
Once the image is running, the terminal will wait for requests to be made using specific URL routes in the background. Using the HTTPS URL displayed in the terminal (you can use `localhost:5000` if you executed in the background), paste the URL and append the following routes at the end of the URL to obtain the desired functions. 

* A POST request to `/data` loads the HGNC data to a Redis database.
*   The command will look like `curl -X POST <URL>/data`.
* A GET request to `/data` should return all populated data from the Redis database as a JSON list.
*   The command will look like `curl -X GET <URL>/data`.
* A DELETE request to `/data` should delete all data from the Redis database.
*   The command will look like `curl -X DELETE <URL>/data`.
* `/genes` returns a json-formatted list of all hgnc_id fields.
* `/genes/<hgnc_id>` route should return all data associated with a given `<hgnc_id>`. 

## Output and What to Expect
In running the main script from an image, the user should receive the respective information printed out to the terminal once running the routes above. Some example commands are shown below.
##### `curl -X POST localhost:5000/data`
```
["The POST request is completed"]
```
##### `curl localhost:5000/genes`
Here, the output is a large list which was shortened for this example. 
```
[ . . . HGNC:13093, HGNC:13104, HGNC:4929,HGNC:21193, HGNC:21961, HGNC:12978, HGNC:26673, HGNC:33517, HGNC:14097, HGNC:20812, HGNC:16155, HGNC:30990, HGNC:16157, HGNC:25704, HGNC:29299, HGNC:43768, HGNC:43769, HGNC:43770, HGNC:29316, HGNC:26993, HGNC:23528, HGNC:45103, HGNC:34495, HGNC:21224, HGNC:13194, HGNC:25468, HGNC:13195, HGNC:13198, HGNC:13199, HGNC:28160, HGNC:32058, HGNC:38032, HGNC:25820, HGNC:13200, HGNC:51695, HGNC:29027, HGNC:24523 ]
```
##### `curl localhost:5000/genes/HGNC:37133`
```
{'gene_group': ['Antisense RNAs'], 'locus_type': 'RNA, long non-coding', 'locus_group': 'non-coding RNA', 'ucsc_id': 'uc002qse.3', 'date_name_changed': '2012-08-15', 'date_modified': '2013-06-27', 'ensembl_gene_id': 'ENSG00000268895', 'vega_id': 'OTTHUMG00000183508', 'lncipedia': 'A1BG-AS1', 'prev_symbol': ['NCRNA00181', 'A1BGAS', 'A1BG-AS'], 'refseq_accession': ['NR_015380'], 'hgnc_id': 'HGNC:37133', 'symbol': 'A1BG-AS1', 'status': 'Approved', 'agr': 'HGNC:37133', 'location_sortable': '19q13.43', 'rna_central_id': ['URS00007E4F6E'], 'ena': ['BC040926'], 'location': '19q13.43', 'gene_group_id': [1987], 'name': 'A1BG antisense RNA 1', 'date_symbol_changed': '2010-11-25', 'date_approved_reserved': '2009-07-20', 'entrez_id': '503538', 'prev_name': ['non-protein coding RNA 181', 'A1BG antisense RNA (non-protein coding)', 'A1BG antisense RNA 1 (non-protein coding)'], 'alias_symbol': ['FLJ23569'], '_version_': 1793942586732314624, 'uuid': '493a38e2-a822-4a93-9228-ddd05504fb4b'}
```
