import math
import sys
import numpy as np
import sdf

sys.path.insert(0,'../lare3d-tools')
from print_energy import Energy

def delta(connectivities, j):
    return ((connectivities[j] - connectivities[j-1]) != 0).astype(int)\
#          + ((connectivities[j] - connectivities[j-1]) < 0).astype(int)

def fix_current_restart(data):
    minimum = data[1:].argmin()
    if data[minimum+1] < 0.001:
        data[minimum+1] = 0.5*(data[minimum] + data[minimum+2])
        
def fix_ohmic_heating_restart(data):
    for i in range(2):
        diff = data[1:] - data[:-1]
        if np.any(np.abs(diff) > 0.01):
            max_point = np.argmax(np.abs(diff))
#             print(max_point, data[max_point-1:max_point+1], diff[max_point-1:max_point+1])
            if diff[max_point] > 0:
                data[max_point+1:] += diff[max_point]
            else:
                data[max_point+1:] -= diff[max_point]

def get_maxes(filenames, index, min_time=-1, max_time=-1):
    energies = [Energy(f) for f in filenames]
    return [get_max(en.data, index, min_time, max_time) for en in energies]

def get_max(data, index, min_time=-1, max_time=-1):
    if max_time > 0:
        max_index = find_nearest(data[:,0], max_time)
    else:
        max_index = len(data[:,0])+1
        
    if min_time > 0:
        min_index = find_nearest(data[:,0], min_time)
    else:
        min_index = 0
    
    return np.max(data[min_index:max_index, index])

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

