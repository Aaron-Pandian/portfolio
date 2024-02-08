#!/usr/bin/env python3

# Imports
import csv
import logging
import math
from math import acos, asin, radians, cos, sin

# Global variables / constants

# Class definitions

# Function definitions


def great_circle(latitude1: float, longitude1: float, latititude2: float, longitude2: float):
    """
    Calculates the geometric great circle, or circular intersection of a sphere and a plane, created by two points, passing through the sphere center. This function is used in the context of Earth, so the two points are coordinate points by longitude and latitude. Requires the math library.

    Args:
        latitude1 (float): The latitudinal coordinate of the first point.

        longitude1 (float): The longitudinal coordinate of the first point.

        latitude2 (float): The latitudinal coordinate of the second point.

        longitude2 (float): The longitudinal coordinate of the second point.

    Returns:
        result (float): The distance of the great circle generated from the given two points, where a point is one longitude and latitude set.
    """
    a = radians(latitude1)
    b = radians(latititude2)

    if longitude1 > longitude2:
        x = radians(longitude1)
        y = radians(longitude2)
    else:
        y = radians(longitude1)
        x = radians(longitude2)

    r = 6378  # Earth radius in km
    try:
        distance = r*acos((cos(a)*cos(b)*cos(x-y))+(sin(a)*sin(b)))
    except TypeError:
        logging.warning(f'encountered non-float value input in great_circle')

    return distance


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

    # Implement function here
    max_distance = ['','',0]
    for meteor1 in ml_data['meteorite_landings']:
        for meteor2 in ml_data['meteorite_landings']:
            if meteor1['name'] == meteor2['name']:
                pass
            else:
                distance = great_circle(
                    float(meteor1['reclat']), float(meteor1['reclong']), float(meteor2['reclat']), float(meteor2['reclong']))
                if distance > max_distance[2]:
                    max_distance = [meteor1['name'], meteor2['name'], distance]

    print("The biggest circle distance is {}".format(max_distance[2]))

    logging.error('Function great_circle has executed successfully')

    logging.info('Ending main script')


# Call to main function
if __name__ == '__main__':
    main()
