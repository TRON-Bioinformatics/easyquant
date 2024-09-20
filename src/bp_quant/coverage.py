
def perc_true(lst):
    n = len(lst)
    num_true = sum(1 for val in lst if val > 0)
    return float(num_true)/n
