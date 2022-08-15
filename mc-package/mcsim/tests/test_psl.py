"""
Tests for Python Standard Library code.
"""

import mcsim.monte_carlo

def test_calculate_distance():

    # Define test data
    point1 = [0, 0, 0]
    point2 = [1, 0, 0]

    # Define expected behavior
    expected_distance = 1

    # Use function to get observed behavior.
    observed_distance = mcsim.monte_carlo.calculate_distance(point1, point2)

    # Assert that our expected behavior matches our observed behavior.
    assert observed_distance == expected_distance

def test_calculate_distance_periodic():

    # Define test data
    point1 = [0, 0, 0]
    point2 = [8, 0, 0]
    box_length = 10

    # Define expected behavior
    expected_distance = 3

    # Use function to get observed behavior.
    observed_distance = mcsim.monte_carlo.calculate_distance(point1, point2, box_length=10)

    # Assert that our expected behavior matches our observed behavior.
    assert observed_distance == expected_distance, "My test failed."