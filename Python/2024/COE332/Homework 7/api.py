#!/usr/bin/env python3

# Imports
import requests
from flask import Flask, request
import numpy as np
import redis
import json
import csv
from jobs import add_job, get_job_by_id, get_job_ids

# Global variables/constants
app = Flask(__name__)
rd = redis.Redis(host='redis-db', port=6379, db=0)

# Function definitions
@app.route('/data', methods=['GET', 'POST', 'DELETE'])
def handle_data():
    if request.method == 'POST':
        # Get the genome data from the website, use the CSV file with times
        response = requests.get("https://data.austintexas.gov/api/views/dx9v-zd7x/rows.csv?accessType=DOWNLOAD")
        # Create CSV dictionary as JSON format
        data = {}
        data['gene_data'] = [] 
        open('data.csv', 'wb').write(response.content)
        with open('data.csv', 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                data['gene_data'].append(dict(row))
        # Iterate over the data to store in redis
        for item in data['gene_data']:
            rd.set(item['Traffic Report ID'], json.dumps(item)) # Item at index one of datapoint is ID number
        # Return response
        os.remove("data.csv")
        return "The POST request is completed\n"

    elif request.method == 'GET':
        return_value = []
        # Iterate over keys in redis
        for item in rd.keys():
            return_value.append(json.loads(rd.get(item)))
        # Return everything as a JSON object
        return return_value

    elif request.method == 'DELETE':
        # Delete everything in redis
        for item in rd.keys():
            rd.delete(item)
        # Return response
        return "The DELETE request is completed\n"

@app.route('/ids', methods=['GET'])
def get_ids():
    return_value = []
    # Iterate over keys in redis
    for item in rd.keys():
    # Append the keys
        value = return_value.append(item.decode('utf8'))
    # Return everything as a JSON object
    return return_value

@app.route('/ids/<desired_id>', methods=['GET'])
def get_id_data(desired_id):
    return_value = []
    # Iterate over keys in redis
    for item in rd.keys():
        # Check if specified key
        if desired_id == item.decode('utf8'):
            return_value.append(json.loads(rd.get(item)))
    # Return response once found
    return return_value

@app.route('/jobs', methods=['POST', 'GET'])
def submit_job():
    if request.method == 'POST':
        try:
            data = request.get_json()
            # Set parameters to be the start and end dates
            job_dict = add_job(data['start'], data['end'])
            return job_dict
        except TypeError:
            return 'There is a problem with the job input, make sure to follow\
 the correct date format: "XX/XX/XXXX" for each start and end date.\n'

@app.route('/jobs/<jobid>', methods=['GET'])
def get_job(jobid):
    return get_job_by_id(jobid)

# Main function definition

# Run Flask
if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=5000)
