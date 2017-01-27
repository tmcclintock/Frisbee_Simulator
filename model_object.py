"""
This file creates a 'model object' which is initalized with 
a Frisbee's parameter values. This file calculates the
amplitudes of the lift force, drag force, and
torques about each principle axis.

For more information see, e.g. Hummel 2003.
"""

import numpy as np
alpha_0 = 4.0*np.pi/180. #radians (4 degrees), constant value reported in Hummel 2003

class Model(object):
	def __init__(self, PL0, PLa, PD0, PDa, PTya, PTywy, PTy0, PTxwx, PTxwz, PTzwz):
                """
                Constructor

                PL0: lift parameter at 0 angle of attack (alpha)
                PLa: linear lift parameter that multiplies angle of attack
                PD0: drag parameter at alpha_0, or the angle of attack at minimum drag
                PDa: quadratic drag parameter, multiplies square angle of attack
                PTy0: y-axis torque parameter (pitch) at alpha = 0
                PTya: y-axis torque parameter linear in alpha
                PTywy: y-axis torque parameter linear in y-axis angular velocity (wy)
                PTxwx: x-axis torque parameter linear in x-axis angular velocity (wx)
                PTxwz: x-axis torque parameter linear in z-axis angular velocity (wz)
                PTzwz: z-axis torque parameter linear in z-axis angular velocity (wz)
                """
		self.PL0=PL0 
		self.PLa=PLa 
		self.PD0=PD0 
		self.PDa=PDa 
		self.PTy0=PTy0 
		self.PTya=PTya 
		self.PTywy=PTywy 
		self.PTxwx=PTxwx 
		self.PTxwz=PTxwz 
		self.PTzwz=PTzwz 

        def __str__(self):
                return "Model: PL0=%.2e\t PLa=%.2e\t PD0=%.2e\t PDa=%.2e\t PTya=%.2e\t PTywy=%.2e\t PTy0=%.2e\t PTxwx=%.2e\t PTxwz=%.2e\t PTzwz%.2e"%(self.PL0, self.PLa, self.PD0, self.PDa, self.PTya, self.PTywy, self.PTy0, self.PTxwx, self.PTxwz, self.PTzwz)

	def C_lift(self, alpha):
		return self.PL0 + self.PLa*alpha

	def C_drag(self,alpha):
		return self.PD0 + self.PDa*(alpha-alpha_0)*(alpha-alpha_0)

	def C_y(self, alpha, wy):
		return self.PTy0 + self.PTywy*wy + self.PTya*alpha

	def C_x(self, wx, wz):
		return self.PTxwz*wz + self.PTxwx*wx

	def C_z(self, wz):
		return self.PTzwz*wz

