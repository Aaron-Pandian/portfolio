#!/usr/bin/env python3

# Imports
import requests
from flask import Flask, request
import numpy as np
import redis
import json
import logging

# Global variables / constants
app = Flask(__name__)

rd = redis.Redis(host='redis-db', port=6379, db=0)

# Class definitions

# Function definitions
@app.route('/data', methods=['GET', 'POST', 'DELETE'])
def handle_data():
    if request.method == 'POST':
        # Get the genome data from the website
        response = requests.get("https://g-a8b222.dd271.03c0.data.globus.org/pub/databases/genenames/hgnc/json/hgnc_complete_set.json")
        data = response.json()
        # Iterate over the data to store in redis
        for item in data['response']['docs']:
            rd.set(item['hgnc_id'], json.dumps(item))
        # Return response
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

@app.route('/genes', methods=['GET'])
def get_ids():
    return_value = []
    # Iterate over keys in redis
    for item in rd.keys():
    # Append the keys
        value = return_value.append(item.decode('utf8'))
    # Return everything as a JSON object
    return return_value

@app.route('/genes/<desired_id>', methods=['GET'])
def get_id_data(desired_id):
    return_value = []
    # Iterate over keys in redis
    for item in rd.keys():
        # Check if specified key
        if desired_id == item.decode('utf8'):
            return_value.append(json.loads(rd.get(item)))
    # Return response once found
    return return_value

# Main function definition

# Run Flask
if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=5000)
