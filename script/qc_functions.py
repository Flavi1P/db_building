import pandas as pd
import numpy as np
import statistics


def range_test(chla):
    bad_index = [i for i, x in enumerate(chla) if x >= 90]
    qc = np.repeat(3, len(chla))
    qc[bad_index] = 4
    return(qc)

def stuck_value(chla, qc_vector, window = 10):
    chla = np.array(chla)
    ii = 0
    new_qc = qc_vector
    while ii + window < len(chla):
        if np.std(chla[ii: ii + window]) == 0:
            new_qc = np.repeat(4, len(chla))
        ii += 1
    return(new_qc)

def dark_correction(chla, pressure):
    if max(pressure) < 100:
        new_qc = np.repeat(3, len(chla))


if __name__ == "__main__":
    chla = np.array([-0.0146, -0.0146, -0.0073, -0.0146, -0.0073, -0.0219, -0.0292, -0.0219, -0.0146, -0.0073,  0.0000, -0.0219, -0.0292, 100, -0.0219, -0.0146, -0.0146, -0.0146, 0.0219, 0.2044,
     0.3796, 0.6351, 0.6643, 0.9490, 1.0731, 1.0074, 1.0512, 1.2045, 1.2921, 1.2118, 1.2483, 1.5330, 1.5914, 1.2775, 1.2191, 1.1607, 1.2848, 1.1680, 1.2410])
    qc = range_test(chla)
    qc = stuck_value(chla, qc)
    print(qc)

