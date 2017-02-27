C_frisbee
=========
This has a C implentation of the frisbee simulation code above, but with a custom RK4 ODE integrator instead of a pre-packaged integrator. This is being built for faster calculations.

Because C is not object oriented, it doesn't make sense to represent a frisbee in the same way as in the python code. Instead, the user is presented with a function call that calculates a trajectory based on a set of parameters and initial conditions. This can then be called in the get_trajectory() function in frisbee.py, as an alternative to running the odeint() integrator that comes from scipy.