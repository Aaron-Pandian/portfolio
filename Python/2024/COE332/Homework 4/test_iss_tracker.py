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
from iss_tracker import get_dataset, full_epoch, time_range, calculate_speed
import pytest

# Global variables / constants

# Class definitions

# Function definitions
def test_get_dataset_exceptions(url: str):
    """
    Testing how the get_dataset function handles errors. 

    Args:
        None

    Returns:
        None
    """
    with pytest.raises(ValueError):
        get_dataset(int('10')) 
    with pytest.raises(RuntimeError):
        get_dataset('too','many','inputs')

def test_time_range_exceptions(time1: str, time2: str):
    """
    Testing how the time_range function handles errors. 
 
    Args:
        None

    Returns:
        None
    """
    with pytest.raises(ValueError):
        time_range(int('10'), 4.432) 
    with pytest.raises(RuntimeError):
        time_range('too','many','inputs')

def test_full_epoch_exceptions(timeStamp:str, stateVector: list, velocities: list):
    """
    Testing how the full_epoch function handles errors. 

    Args:
        None

    Returns:
        None 
    """
    with pytest.raises(ValueError):
        full_epoch(int('10'), 4.432, [2, 3, 4]) 
    with pytest.raises(ZeroDivisionError):
        full_epoch('12th', [], ['3', '4', '5'])     
    with pytest.raises(RuntimeError):
        full_epoch('too','many','inputs')

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
    with pytest.raises(ZeroDivisionError):
        calculate_speed([], 0)     
    with pytest.raises(RuntimeError):
        calculate_speed('too','many','inputs')
    with pytest.raises(KeyError):
        compute_average_mass([{'a': 1}, {'b': 1}], 'a')