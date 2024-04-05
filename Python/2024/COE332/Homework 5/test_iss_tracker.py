#!/usr/bin/env python3

# Imports
import xmltodict
import xml
import logging
import statistics
from statistics import mean
import requests
import datetime
from datetime import datetime, timedelta
import math
from math import sqrt
from iss_tracker import get_dataset, full_epoch, time_range, calculate_speed, return_iss_dataset, return_iss_now, return_iss_speed, return_iss_state, app
import pytest
from flask import Flask, request

# Global variables / constants
app.testing = True
client = app.test_client()
# Class definitions

# Function definitions
def test_get_dataset_exceptions():
    """
    Testing how the get_dataset function handles errors. 

    Args:
        None

    Returns:
        None
    """
    with pytest.raises(ValueError):
        get_dataset(int('10')) 
    with pytest.raises(TypeError):
        get_dataset('too','many','inputs')

def test_time_range_exceptions():
    """
    Testing how the time_range function handles errors. 
 
    Args:
        None

    Returns:
        None
    """
    with pytest.raises(TypeError):
        time_range('too','many','inputs')

def test_full_epoch_exceptions():
    """
    Testing how the full_epoch function handles errors. 

    Args:
        None

    Returns:
        None 
    """
    with pytest.raises(ValueError):
        full_epoch('12th', [], ['3', '4', '5'])     

def test_calculate_speed():
    """
    Testing truths to validate the calculate_speed funciton.

    Args:
        None

    Returns:
        None
    """
    test_states = {'newtime': [{'EPOCH': '2024-045T12:00:00.000Z', 'X': {'@units': 'km', '#text': '-3'}, 'Y': {'@units': 'km', '#text': '-4'}, 'Z': {'@units': 'km', '#text': '5'}, 'X_DOT': {'@units': 'km/s', '#text': '6'}, 'Y_DOT': {'@units': 'km/s', '#text': '-7'}, 'Z_DOT': {'@units': 'km/s', '#text': '8'}} ]}
    test_items = 2
    avg, inst = calculate_speed(test_states,test_items)
    assert round(avg) == 12
    assert round(inst) == 0 # because no second timestep called
    assert isinstance(avg, float) == True

def test_calculate_speed_exceptions():
    """
    Testing how the calculate_speed function handles errors. 

    Args:
        None

    Returns:
        None
    """
    with pytest.raises(TypeError):
        calculate_speed([], 0)     

def test_return_iss_dataset(client):
    """
    Testing truths to validate the return_iss_dataset funciton.

    Args:
        None

    Returns:
        None
    """
    response = client.get('/epochs')
    datasetBytes = response.data
    dataset = datasetBytes.decode('utf8')
    assert isinstance(dataset, str) == True

def test_return_iss_speed(client):
    """
    Testing truths to validate the return_iss_speed funciton.

    Args:
        None

    Returns:
        None
    """
    response = client.get('/epochs/3/speed')
    finalSpeedBytes = response.data
    finalSpeed = finalSpeedBytes.decode('utf8')
    assert isinstance(finalSpeed, str) == True

def test_return_iss_now(client):
    """
    Testing truths to validate the return_iss_now funciton.

    Args:
        None

    Returns:
        None
    """
    response = client.get('/now')
    nowStateBytes = response.data
    nowState = nowStateBytes.decode('utf8')
    assert isinstance(nowState, str) == True

def test_return_iss_state(client):
    """
    Testing truths to validate the return_iss_state funciton.

    Args:
        None

    Returns:
        None
    """
    response = client.get('/epoch/10')
    stateBytes = response.data
    stateVector = stateBytes.decode('utf8')
    assert isinstance(stateVector, str) == True

def unit_test():
    """
    Run all unit test functions.

    Args:
        None

    Returns:
        None
    """
    test_get_dataset_exceptions()
    test_time_range_exceptions()
    test_full_epoch_exceptions()
    test_calculate_speed()
    test_calculate_speed_exceptions()
    test_return_iss_dataset(client=client)
    test_return_iss_now(client=client)
    test_return_iss_speed(client=client)
    test_return_iss_state(client=client)

if __name__ == '__main__':
    unit_test()