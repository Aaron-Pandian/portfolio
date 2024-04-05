#!/usr/bin/env python3

# Imports
from ml_data_anaysis import compute_statistics
from ml_data_anaysis import most_frequent 
import pytest

# Global variables / constants

# Class definitions

# Function definitions
def test_compute_statistics():
    """
    Tests the compute_statistics function from ml_data_analysis script.

    Args:
        None

    Returns:
        None
    """
    # check for average
    assert compute_statistics([{'a': 1}], 'a')[1] == 1
    assert compute_statistics([{'a': 1}, {'a': 2}], 'a')[1] == 1.5
    assert compute_statistics([{'a': 1}, {'a': 2}, {'a': 3}], 'a')[1] == 2
    assert compute_statistics([{'a': 10}, {'a': 1}, {'a': 1}], 'a')[1] == 4
    assert isinstance(compute_statistics([{'a': 1}, {'a': 2}], 'a')[1], float) == True
    
    # check for standard deviation
    assert compute_statistics([{'a': 1}], 'a')[2] == 0
    assert compute_statistics([{'a': 1}, {'a': 2}], 'a')[2] == 0.5
    assert round(compute_statistics([{'a': 1}, {'a': 2}, {'a': 3}], 'a')[2],2) == 0.82
    assert isinstance(compute_statistics([{'a': 1}, {'a': 2}], 'a')[2], float) == True

def test_compute_statistics_exceptions():
    """
    Tests the compute_statistics function from ml_data_analysis script for exception errors.

    Args:
        None

    Returns:
        None
    """
    # sent an empty list
    with pytest.raises(ZeroDivisionError):
        compute_statistics([], 'a')
    # dictionaries not uniform
    with pytest.raises(KeyError):
        compute_statistics([{'a': 1}, {'b': 1}], 'a')
    # value not a float
    with pytest.raises(ValueError):
        compute_statistics([{'a': 1}, {'a': 'x'}], 'a')
    # key not in dicts
    with pytest.raises(KeyError):
        compute_statistics([{'a': 1}, {'a': 2}], 'b')


def test_most_frequent():
    """
    Tests the most_frequent function from ml_data_analysis script.  

    Args:
        None

    Returns:
        None 
    """
    assert most_frequent([1,2,3,4,5,3,2,3,5,6,4,2,2]) == 2
    assert most_frequent(["one","two","two","one","three","four","one"]) == "one"
    assert most_frequent(["one","two","two","one","three","four","one",1,2,5,6,7,8]) == "one"
    assert isinstance(most_frequent([1,2,3,4,5,3,2,3,5,6,4,2,2]), int) == True
    
def test_most_frequent_exceptions():
    """
    Tests the most_frequent function from ml_data_analysis script for exception errors.  

    Args:
        None

    Returns:
        None 
    """
    with pytest.raises(ZeroDivisionError):
        ml.most_frequent([])

# Main function definition


def main():

    # Running tests
    test_compute_statistics()
    test_most_frequent()

# Call to main function
if __name__ == '__main__':
    main()
