#!/usr/bin/env python3

# Imports
import csv
import logging
import gcd_algorithm as gcd
import matplotlib.pyplot as plt 
import numpy as np
import math
from scipy.stats import norm 

# Global variables / constants

# Class definitions

# Function definitions


def compute_statistics(a_list_of_dicts: dict, a_key_string: str):
    """
    Iterates through a list of dictionaries, pulling out values associated with a given key. Returns the mean and standard deviation of those values as well as a list of all samples. 

    Args:
        a_list_of_dicts (list): A list of dictionaries, each dict should have the same set of keys.

        a_key_string (string): A key that appears in each dictionary associated
        with the desired value (will enforce float type).

    Returns:
        result (list): A list containing the average value and population standard deviation as floats, and a list of all the masses. 
    """
    samples = []
    total_mass = 0.
    for i in a_list_of_dicts:
        try:
            total_mass += float(i[a_key_string])
            samples.append(float(i[a_key_string]))
        except TypeError:
            logging.warning(f'encountered non-float value {i[a_key_string]} in compute_statistics')
        except KeyError:
            logging.warning(f'could not find key {a_key_string} in the {a_list_of_dicts} dictionary for compute_statistics')

    avg = (total_mass / len(a_list_of_dicts))

    total_deviation = 0
    for i in a_list_of_dicts:
        total_deviation += (float(i[a_key_string]) - avg)**2
    sd = (total_deviation/(len(a_list_of_dicts)))**(1/2)

    stats = [samples, avg, sd]

    return stats


def most_frequent(list: list):
    """
    Parses through given list and return the most frequent 
    value within that list.  

    Args:
        list (list): Contains values of the same type, must be one- dimentional. (Will not enforce type).

    Returns:
        mode (N/A): The value of the most common variable within
        the input list. The most common variable depends on what type
        of variables are held within the list and the frequency of 
        said variable. 
    """
    counter = 0
    mode = []
    try:
        for i in list:
            length = list.count(i)
            if length > counter:
                counter = length
                mode = i
    except TypeError:
        logging.warning(f'encountered non-list value input in most_frequent')

    return mode


# Main function definition
def main():

    logging.basicConfig(level='DEBUG')
    logging.debug('Starting main script')

    ml_data = {}
    ml_data['meteorite_landings'] = []

    with open('Meteorite_Landings.csv', 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            ml_data['meteorite_landings'].append(dict(row))
        logging.warning('The file has been read')

    print("\n")
    print("Summary Statistics:")

    # Function 1
    meteor_class = []
    for row in ml_data['meteorite_landings']:
        meteor_class.append(row['recclass'])

    print("Most meterorites are of class {}\n".format(most_frequent(meteor_class)))
    logging.error('Function most_frequent has executed successfully\n')

    # Function 2
    max_distance = ['','',0]
    for meteor1 in ml_data['meteorite_landings']:
        for meteor2 in ml_data['meteorite_landings']:
            if meteor1['name'] == meteor2['name']:
                pass
            else:
                distance = gcd.great_circle(
                    float(meteor1['reclat']), float(meteor1['reclong']), float(meteor2['reclat']), float(meteor2['reclong']))
                if distance > max_distance[2]:
                    max_distance = [meteor1['name'],meteor2['name'], distance]
    
    print("The maximum distance is between meteoroid {} and {}, with a distance of {}\n".format(max_distance[0], max_distance[1], max_distance[2]))
    logging.error('Function great_circle has executed successfully')

    # Function 3
    logging.getLogger('matplotlib').propagate=False
    x = compute_statistics(ml_data['meteorite_landings'], 'mass (g)')
    x_axis = np.arange(min(x[0]), max(x[0]), 10)
    plt.plot(x_axis, norm.pdf(x_axis, x[1], x[2]))
    plt.xlabel('Probabilities')
    plt.ylabel('Masses')
    plt.title("Normal Distribution of Meteor Landings' Mass")
    plt.show()

    logging.info('Ending main script')


# Call to main function
if __name__ == '__main__':
    main()
