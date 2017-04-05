Here we have some files that perform some tests on our code.

For instance, the compare_codes.py file comares the same throw calculated
using the python routine and the C code. These aren't expected to be exactly
the same because scipy.odeint() chooses it's own integrator, while the C
code uses RK4.