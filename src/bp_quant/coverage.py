"""
Library for custom methods.
"""

def perc_true(values) -> float:
    """Calculate the percentage of positive non-zero values in list.

    Args:
        values (list): Input list with values

    Returns:
        float: Percentage of positive non-zero values in list
    """

    num_true = sum(1 for val in values if val > 0)
    return float(num_true) / len(values)
