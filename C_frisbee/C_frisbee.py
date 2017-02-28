"""
This is a python wrapper for the c_frisbee.so shared library, which can be compiled in this directory from the C code that is present.
"""

import numpy as np
import ctypes, os, sys, inspect
from ctypes import c_double, c_int, POINTER, cdll

#This is to allow this to work from anywhere, in case we implement a setup.py
library_path = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))+"/c_frisbee.so"
cflib = cdll.LoadLibrary(library_path)

def get_trajectory(initial_positions,flight_time,N_times,params):
    """Get the trajectory from the driver.

    Args:
        initial_positions (array_like): All initial positions and velocities of the disc.
        flight_time (float): Total flight time; seconds.
        N_times (float): The number of times to evaluate the positions at.
        params (array_like): List of force/torque parameters.

    Returns:
        times (array_like): List of all times that the positions are evaluated at.
        trajectory (array_like): 2d array of dimensions 12 by N_times that has all coordinates of the disc over the course of the trajectory.

    """
    #Outline the argument types
    driver = cflib.driver
    driver.argtypes = [POINTER(c_double),POINTER(c_double),c_double,c_int,POINTER(c_double)]

    #Create an array for the trajectory+time
    #Note: we add an extra dimension, to hold the times.
    #This is the only way to send an array through ctypes.
    N = len(initial_positions)
    all_positions = np.zeros(((N+1)*N_times))

    #Create pointers for the input arrays
    initial_positions_in = initial_positions.ctypes.data_as(POINTER(c_double))
    params_in = params.ctypes.data_as(POINTER(c_double))
    all_positions_in = all_positions.ctypes.data_as(POINTER(c_double))

    #Call the driver
    driver(initial_positions_in,params_in,flight_time,N_times,all_positions_in)
    all_positions = all_positions.reshape((N_times,N+1))
    times = np.copy(all_positions[:,0])
    all_positions = np.copy(all_positions[:,1:])

    #Return an array of all of the positions at all the times
    return times,all_positions
