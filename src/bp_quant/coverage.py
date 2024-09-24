"""
Library for custom methods.
"""

def perc_true(lst) -> float:
    """Calculate the percentage of positive non-zero values in list.

    Args:
        lst (list): Input list with values

    Returns:
        float: Percentage of positive non-zero values in list
    """
    n = len(lst)
    num_true = sum(1 for val in lst if val > 0)
    return float(num_true)/n
