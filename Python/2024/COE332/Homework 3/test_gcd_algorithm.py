#!/usr/bin/env python3

# Imports
from gcd_algorithm import great_circle
import pytest

# Global variables / constants

# Class definitions

# Function definitions
def test_great_circle():
    """
    Tests the most_frequent function from ml_data_analysis script.

    Args:
        None

    Returns:
        None
    """
    # value not a float
    assert round(great_circle(36,85,19,83)) == 1903
    assert round(great_circle(36.3,85.1,19.9,83.3)) == 1834
    assert round(great_circle(-36.3,85.1,19.9,-83.3)) == 17889
    assert isinstance(great_circle(-36.3,85.1,19.9,-83.3), float) == True



def test_great_circle_exceptions():
    """
    Tests the most_frequent function from ml_data_analysis script for exception errors.

    Args:
        None

    Returns:
        None
    """
    # value not a float
    with pytest.raises(TypeError):
        great_circle(1,'a',2,3)

# Main function definition
def main():

    # Running tests
    test_great_circle()
    test_great_circle_exceptions()

# Call to main function
if __name__ == '__main__':
    main()
