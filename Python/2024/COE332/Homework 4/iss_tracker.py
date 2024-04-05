#!/usr/bin/env python3

# Imports
import xmltodict
import logging
import statistics
from statistics import mean
import requests
import datetime
from datetime import datetime, timedelta
import math
from math import sqrt

# Global variables / constants

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
        print('The input dataset is not what this function is intended for.')

    return states, summary, items

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
        print('Either the state or the velocity vector input was not a list')

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
    for row in states['newtime']:
        # Finding velocity elements in provided list dictionary
        vx = float(row['X_DOT']['#text'])
        vy = float(row['Y_DOT']['#text'])
        vz = float(row['Z_DOT']['#text'])

        try:
            vmag = sqrt((vx**2)+(vy**2)+(vz**2))
        except OverflowError:
            print('The computation is too large to be represented')

        if counter == last_index:
            instSpeed = vmag
        speeds.append(vmag)
        counter += 1

    speed = mean(speeds)
    return speed, instSpeed

# Main function definition
def main():

    logging.basicConfig(level='DEBUG')
    logging.debug('Starting main script')

    # obtaining url from website link
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

# Call to main function
if __name__ == '__main__':
    main()