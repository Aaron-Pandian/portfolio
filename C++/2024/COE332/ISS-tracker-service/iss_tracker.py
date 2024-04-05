#!/usr/bin/env python3

# Imports
import xmltodict
import logging
import statistics
from statistics import mean
import requests
import math
from math import sqrt
from flask import Flask, request
from geopy.geocoders import Nominatim
import numpy as np
from functools import partial

# Global variables / constants
app = Flask(__name__)
geolocator = Nominatim(user_agent = 'agent')

# Class definitions

# Function definitions
def get_dataset(url: str):
    """
    Ingests an url for a website with an xml dataset. Then gets the dataset using the requests library and splits the dataset into two list-dictionaries (a list of dictionaries)- one for summary information and on for data. 

    Args:
        url (str): The website url accessing the xml dataset

    Returns:
        states (list): A list of iterable python dictionaries for the states of the spacecraft at each timestamp. 

        summary (list): A list of iterable python dictionaries for initial comments in the dataset. 

        items (int): The integer number of timestamp recordings of spacecraft state data. 
    """
    try:
        response = requests.get(url)
        response.status_code
        response.content
        open('ISS.xml', 'wb').write(response.content)
    except TypeError:
        logging.warning('The input value is not a valid string')

    # logging state vectors of [UTC time, position, and velocity] in order
    states = {}
    states['newtime'] = []
    summary = {}
    summary['comment'] = []
    items = 0

    try:
        with open('ISS.xml', 'r') as f:
            reader = xmltodict.parse(f.read())
            for row in reader['ndm']['oem']['body']['segment']['data']['stateVector']:
                states['newtime'].append(row)
                items += 1
            for row in reader['ndm']['oem']['body']['segment']['data']['COMMENT']:
                summary['comment'].append(row)
    except KeyError:
        logging.warning('The input dataset is not what this function is intended for.')

    return states, summary, items

@app.route('/epochs', methods=['GET']) # Fix this, it only takes in the first query parameter defined in the route. 
def return_iss_dataset():
    """
    Finds ISS tracking dataset from xml website dataset. Then aquires the dataset using the requests library and outputs the list-dictionaries (a list of dictionaries). The optional query parameters allow for an indexing of the whole dataset. Using the limit or end index, an offset value is provided to denote the number of datapoints. The final output is a dictionary that ends at index 'limit' and starts at 'offset' points before the 'limit' index. 

    Args:
        None
        
    Returns:
        dataset (dict): A list of iterable python dictionaries that make up the ISS tracking dataset.  
    """
    url = 'https://nasa-public-data.s3.amazonaws.com/iss-coords/current/ISS_OEM/ISS.OEM_J2K_EPH.xml'
    states, summary, items = get_dataset(url)

    try:
        offset = int(request.args.get(key='offset', default=0))
        limit = int(request.args.get(key='limit', default=items))
    except TypeError:
        logging.warning("Invalid offset or limit parameter; both must be an integer.")

    if limit >= offset:
        start = limit - offset
    else:
        start = 0

    aug_states = {}
    aug_states['newtime'] = []
    index = 0

    for row in states['newtime']:
        if index >= start and index < limit:
            aug_states['newtime'].append(row)
        index += 1

    return aug_states

def time_range(time1: str, time2: str):
    """
    Calculates the time range of the dataset using a given initial and final epoch timestamp. 

    Args:
        time1 (str): The timestamp of the first epoch in the dataset as a string. 

        time2 (str): The timestamp of the last epoch in the dataset as a string.

    Returns:
        None
    """
    print("\nThe range of data spans a time from {} to {}\n".format(str(time1), str(time2)))

def full_epoch(timeStamp:str, stateVector: list, velocities: list):
    """
    Prints the epoch charecteristics closest to the time of capturing the dataset from web server. 

    Args:
        timeStamp (str): The timestamp of the last epoch in the dataset as a string. 

        stateVector (list): The X, Y, and Z position values of the spacecraft in its last recorded timestamp as a string list.

        stateVector (list): The X, Y, and Z velcoity values of the spacecraft in its last recorded timestamp as a string list.

    Returns:
        None
    """
    try:
        [x, y, z] = stateVector 
        [vx, vy, vz] = velocities
        print("\nThe most recent state of the spacecraft was recorded at {} with a coordinate position of [{}, {}, {}] km and velocity components of [{}, {}, {}] km/s both in X, Y, Z order\n".format(timeStamp, x, y, z, vx, vy, vz))
    except TypeError:
        logging.warning('Either the state or the velocity vector input was not a list')

def calculate_speed(states: list, items: int):
    """
    Calculates the average and final instantaneous speed of the spacecraft from the dataset. 

    Args:
        states (list): A list of dictionaries containing all state values at each timestamp across the dataset as strings.  

        items (int): The integer number of timestamp recordings of spacecraft state data. 

    Returns:
        avgSpeed (float): The average speed of the spacecraft across the time range of the dataset. 

        instSpeed (float): The instantaneous speed of the spacecraft closest to the time of capturing the dataset from web server.
    """
    speeds = []
    instSpeed = 0

    last_index = items - 1
    counter = 0

    try:
        for row in states['newtime']:
            # Finding velocity elements in provided list dictionary
            vx = float(row['X_DOT']['#text'])
            vy = float(row['Y_DOT']['#text'])
            vz = float(row['Z_DOT']['#text'])

            try:
                vmag = sqrt((vx**2)+(vy**2)+(vz**2))
            except OverflowError:
                logging.warning('The computation is too large to be represented')

            if counter == last_index:
                instSpeed = vmag
            speeds.append(vmag)
            counter += 1
    except ValueError:
        logging.warning('Please input the correct ISS tracking dataset.')

    speed = mean(speeds)
    return speed, instSpeed

@app.route('/now', methods=['GET']) # update using location function
def return_iss_now():
    """
    Finds ISS tracking dataset from xml website dataset. Then analyzes the state at the final recorded epoch using the requests library.

    Args:
        None
        
    Returns:
        instSpeed (float): A value representative of the instantaneous speed of the ISS in its final recorded epoch. 

        posVec (list): A list of of the latitude, longitude, and altitude of the ISS in its final recorded epoch.

        location.address (str): A string of the nearest address of the ISS based on longitude and latitude in its final recorded epoch. 
    """
    url = 'https://nasa-public-data.s3.amazonaws.com/iss-coords/current/ISS_OEM/ISS.OEM_J2K_EPH.xml'
    states, summary, items = get_dataset(url)

    final_index = items - 1
    end_state = states['newtime'][final_index]

    # Printing final epoch state variables
    stateVector = [end_state['X']['#text'], end_state['Y']['#text'], end_state['Z']['#text'], end_state['X_DOT']['#text'], end_state['Y_DOT']['#text'], end_state['Z_DOT']['#text']]

    # Calculating speeds
    x = float(stateVector[0]) 
    y = float(stateVector[1])
    z = float(stateVector[2])
    vx = float(stateVector[3])
    vy = float(stateVector[4])
    vz = float(stateVector[5])

    try:
        instSpeed = sqrt((vx**2)+(vy**2)+(vz**2))
    except OverflowError:
        logging.warning('The computation is too large to be represented')

    # Because all epoch times are same format and length, slice strings to find the minutes and hours
    isoTime = end_state['EPOCH']
    hrs = int(isoTime[slice(9, 11)])
    mins = int(isoTime[slice(12, 14)])

    # Determine lat, lon, alt, geo pos and return all in list. 
    R = 6371 # km 
    lat = math.degrees(math.atan2(z, math.sqrt(x**2 + y**2))) # degrees
    lon = math.degrees(math.atan2(y, x)) - ((hrs-12)+(mins/60))*(360/24) + 19

    if lon > 180: 
        lon = -180 + (lon - 180)
    if lon < -180: 
        lon = 180 + (lon + 180)

    alt = math.sqrt(x**2 + y**2 + z**2) - R

    # Set posVec
    posVec = [lat, lon, alt]

    # Finding geolocation
    coordstr = str(lat) + ', ' + str(lon)
    reverse = partial(geolocator.reverse, language="es")
    location = reverse(coordstr)
    geolocUnicode = str(location)
    geolocEncoded = geolocUnicode.encode("ascii", "ignore")
    geoloc = geolocEncoded.decode()

    # If no address is found geopy returns None, this is the case when the ISS is over an oceon. 
    if geoloc == 'None':
        locationString = 'The ISS is over the ocean'
        return [instSpeed, lat, lon, alt, locationString]
    else:
        return [instSpeed, posVec, geoloc]

@app.route('/epochs/<epoch>', methods=['GET'])
def return_iss_state(epoch):
    """
    Finds ISS tracking dataset from xml website dataset. Then aquires the full state information from a specific epoch using the requests library.

    Args:
        None
        
    Returns:
        stateVector (list): A list of of the state position and velocity values for a specific epoch in the dataset.
    """
    url = 'https://nasa-public-data.s3.amazonaws.com/iss-coords/current/ISS_OEM/ISS.OEM_J2K_EPH.xml'
    states, summary, items = get_dataset(url)

    try:
        epoch = int(epoch)
    except TypeError:
        logging.warning('Please enter a integer for the specific epoch value.')

    state = states['newtime'][epoch]
    stateVector = [state['X']['#text'], state['Y']['#text'], state['Z']['#text'], state['X_DOT']['#text'], state['Y_DOT']['#text'], state['Z_DOT']['#text']]

    return stateVector

@app.route('/epochs/<epoch>/speed', methods=['GET'])
def return_iss_speed(epoch):
    """
    Finds ISS tracking dataset from xml website dataset. Then aquires the instantaneous speed from a specific epoch using the requests library.

    Args:
        None
        
    Returns:
        vmag (list): A value representative of the instantaneous speed of the ISS for the specific epoch defined.  
    """
    url = 'https://nasa-public-data.s3.amazonaws.com/iss-coords/current/ISS_OEM/ISS.OEM_J2K_EPH.xml'
    states, summary, items = get_dataset(url)

    try:
        epoch = int(epoch)
    except TypeError:
        logging.warning('Please enter a integer for the specific epoch value.')

    state = states['newtime'][epoch]
    vx = float(state['X_DOT']['#text'])
    vy = float(state['Y_DOT']['#text'])
    vz = float(state['Z_DOT']['#text'])

    vmag = sqrt((vx**2)+(vy**2)+(vz**2))
    return [vmag]

@app.route('/comment', methods=['GET'])
def return_iss_comment():
    """
    Finds ISS tracking dataset from xml website dataset. Then aquires summary commentary using the requests library.

    Args:
        None
        
    Returns:
        comments (list): The values denoted in the 'comment' key of the ISS dataset.  
    """
    url = 'https://nasa-public-data.s3.amazonaws.com/iss-coords/current/ISS_OEM/ISS.OEM_J2K_EPH.xml'
    states, summary, items = get_dataset(url)

    return summary['comment']

@app.route('/header', methods=['GET'])
def return_iss_header():
    """
    Finds ISS tracking dataset from xml website dataset. Then aquires header data using the requests library.

    Args:
        None
        
    Returns:
        header (dict): The values denoted in the 'header' key of the ISS dataset.  
    """
    url = 'https://nasa-public-data.s3.amazonaws.com/iss-coords/current/ISS_OEM/ISS.OEM_J2K_EPH.xml'

    response = requests.get(url)
    response.status_code
    response.content
    open('ISS.xml', 'wb').write(response.content)

    header = {}
    header['CREATION_DATE'] = []
    header['ORIGINATOR'] = []

    with open('ISS.xml', 'r') as f:
        reader = xmltodict.parse(f.read())
        header['CREATION_DATE'].append(reader['ndm']['oem']['header']['CREATION_DATE'])
        header['ORIGINATOR'].append(reader['ndm']['oem']['header']['ORIGINATOR'])
        
    return header

@app.route('/metadata', methods=['GET'])
def return_iss_metadata():
    """
    Finds ISS tracking dataset from xml website dataset. Then aquires the metadata using the requests library.

    Args:
        None
        
    Returns:
        metadata (dict): The values denoted in the 'metadata' key of the ISS dataset.  
    """
    url = 'https://nasa-public-data.s3.amazonaws.com/iss-coords/current/ISS_OEM/ISS.OEM_J2K_EPH.xml'

    response = requests.get(url)
    response.status_code
    response.content
    open('ISS.xml', 'wb').write(response.content)

    metadata = {}
    metadata['OBJECT_NAME'] = []
    metadata['OBJECT_ID'] = []
    metadata['CENTER_NAME'] = []
    metadata['REF_FRAME'] = []
    metadata['TIME_SYSTEM'] = []
    metadata['START_TIME'] = []
    metadata['STOP_TIME'] = []

    with open('ISS.xml', 'r') as f:
        reader = xmltodict.parse(f.read())
        index = reader['ndm']['oem']['body']['segment']['metadata']
        
        # Getting all metadata using index
        metadata['OBJECT_NAME'].append(index['OBJECT_NAME'])
        metadata['OBJECT_ID'].append(index['OBJECT_ID'])
        metadata['CENTER_NAME'].append(index['CENTER_NAME'])
        metadata['REF_FRAME'].append(index['REF_FRAME'])
        metadata['TIME_SYSTEM'].append(index['TIME_SYSTEM'])
        metadata['START_TIME'].append(index['START_TIME'])
        metadata['STOP_TIME'].append(index['STOP_TIME'])
        
    return metadata

@app.route('/epochs/<epoch>/location', methods=['GET']) 
def return_iss_location(epoch):
    """
    Finds ISS tracking dataset from xml website dataset. Then analyzes the state at a specific epoch using the requests library.

    Args:
        None
        
    Returns:
        instSpeed (float): A value representative of the instantaneous speed of the ISS in its final recorded epoch. 

        posVec (list): A list of of the latitude, longitude, and altitude of the ISS in its final recorded epoch.

        location.address (str): A string of the nearest address of the ISS based on longitude and latitude in its final recorded epoch. 
    """
    url = 'https://nasa-public-data.s3.amazonaws.com/iss-coords/current/ISS_OEM/ISS.OEM_J2K_EPH.xml'
    states, summary, items = get_dataset(url)

    try:
        epoch = int(epoch)
    except TypeError:
        logging.warning('Please enter a integer for the specific epoch value.')

    state = states['newtime'][epoch]
    x = float(state['X']['#text']) 
    y = float(state['Y']['#text'])
    z = float(state['Z']['#text'])

    # Because all epoch times are same format and length, slice strings to find the minutes and hours
    isoTime = state['EPOCH']
    hrs = int(isoTime[slice(9, 11)])
    mins = int(isoTime[slice(12, 14)])

    # Determine lat, lon, alt, geo pos and return all in list. 
    R = 6371 # km 
    lat = math.degrees(math.atan2(z, math.sqrt(x**2 + y**2))) # degrees
    lon = math.degrees(math.atan2(y, x)) - ((hrs-12)+(mins/60))*(360/24) + 19

    # Correction if calculated outside the bounds
    if lon > 180: 
        lon = -180 + (lon - 180)
    if lon < -180: 
        lon = 180 + (lon + 180)

    alt = math.sqrt(x**2 + y**2 + z**2) - R

    # Set posVec
    posVec = [lat, lon, alt]

    # Finding geolocation and setting it to appropriate type
    coordstr = str(lat) + ', ' + str(lon)
    reverse = partial(geolocator.reverse, language="es")
    location = reverse(coordstr, zoom=20)
    geolocUnicode = str(location)
    geolocEncoded = geolocUnicode.encode("ascii", "ignore")
    geoloc = geolocEncoded.decode()

    # If no address is found geopy returns None, this is the case when the ISS is over an oceon. 
    if geoloc == 'None':
        locationString = 'The ISS is over the ocean'
        return [posVec, locationString]
    else:
        return [posVec, geoloc]

# Main function definition
def main():

    logging.basicConfig(level='DEBUG')
    logging.debug('Starting main script')

    # Obtaining url from website link
    url = 'https://nasa-public-data.s3.amazonaws.com/iss-coords/current/ISS_OEM/ISS.OEM_J2K_EPH.xml'
    states, summary, items = get_dataset(url)

    logging.error('Recieved and interpreted the dataset successfully')

    # Printing time range
    final_index = items - 1 
    start_state = states['newtime'][0]
    end_state = states['newtime'][final_index]
    time = time_range(str(start_state['EPOCH']), str(end_state['EPOCH']))

    logging.debug('Function time_range ran successfully')

    # Printing final epoch state variables
    final_time = str(end_state['EPOCH'])
    stateVector = [end_state['X']['#text'], end_state['Y']['#text'], end_state['Z']['#text']]
    velocities = [end_state['X_DOT']['#text'], end_state['Y_DOT']['#text'], end_state['Z_DOT']['#text']]
    full_epoch(final_time, stateVector, velocities)

    logging.debug('Function full_epoch ran successfully')

    # Calculating speeds
    speed, instSpeed = calculate_speed(states, items)
    print("\nThe average speed of the ISS is {} km/s and the instantaneous speed at the final recording of data is {} km/s\n".format(speed, instSpeed))
    
    logging.debug('Function calculate_speed ran successfully')
    
    logging.info('Ending main script')

# Run Flask
if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0')
