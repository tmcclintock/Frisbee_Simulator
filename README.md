Frisbee_Simulator
=================

A python code that simulates frisbee throws. Much of the code is 
from Elizabeth Hannah's repository called FrisbeeResearch.

This repository contains code that can:
- calculate trajectories of throws given some initial conditions
- plot trajectories for visualization purposes
- perform an analysis of a throw to find the optimal set of aerodynamic parameters that describe that disc (i.e. the force and torque coupling constants)

A very simple example of how to use the code can be found in example/example_throw.py.

It's authors can be contacted at thmsmcclintock@gmail.com (Tom McClintock) and ehannah@email.arizona.edu (Elizabeth Hannah).

Requirements
============

The Python dependencies of this repository can be installed with
```
pip install -r requirements.txt
```

Additionally, in order to compile the fast back-end, you need a C compiler. This project has been tested with gcc 4.4.7, 4.4.8, and clang.