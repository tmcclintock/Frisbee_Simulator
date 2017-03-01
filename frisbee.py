"""
The frisbee object. This contains the code to
initialize a frisbee object along with
its model for the forces that act on it,
and an equations of motion function for use
in an integrator.
"""
import numpy as np
import coefficient_model
from scipy.integrate import odeint
import sys, os, inspect
sys.path.insert(0,os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))+"/C_frisbee/")
c_frisbee_available = True
try: import C_frisbee
except OSError: 
  print "You must compile the C code if you want to use it. See the README."
  c_frisbee_available = False
"""
Constants:
PI
rho: density of air; kg/m^3
area: area of a disc; m^2
diameter: diameter of a disc; m
mass: mass of a disc; kg
g: gravitational acceleration; m/s^2
Izz: z-axis moment of inertial; kg*m^2
Ixx,Iyy,Ixy; x(y,cross)-axis moment of inertia; kg*m^2
F_gravity: gravitational force on the disc; kg*m/s^2
"""
PI=3.141592653589793
rho=1.225
area=0.057
mass=0.175
g=9.81
Izz=0.002352 
Ixx=Iyy=Ixy=0.001219 
diameter=2*(area/PI)**0.5
F_gravity=mass*g*np.array([0.,0.,-1.])


class Frisbee(object):

  def __init__(self,x,y,z,vx,vy,vz,phi,theta,gamma,
               phidot,thetadot,gammadot,use_C=False,debug=False):
    """
    Constructor

    x,y,z: euclidean positions
    vx,vy,vz: velocities
    phi,theta,gamma: euler angles
    phidot,thetadot,gammadot: euler angle time derivatives
    debug: default False; for printing

    Note: all kinematic variables are collectively
    referred to as 'coordinates' since they are the
    variables that change in the equations of motion (EOM).

    Calls update_data_fields to calculate quantities 
    that depend only on positions.
    """
    self.x=x 
    self.y=y
    self.z=z
    self.vx=vx
    self.vy=vy
    self.vz=vz
    self.phi=phi
    self.theta=theta
    self.gamma=gamma
    self.phidot=phidot
    self.thetadot=thetadot
    self.gammadot=gammadot
    self.use_C=use_C
    self.debug=debug
    self.update_data_fields()
    

  def __str__(self):
    outstr = "Frisbee:\n"+\
             "Position: (%f,%f,%f)\n"%(self.x,self.y,self.z)+\
             "Velocity: (%f,%f,%f)\n"%(self.vx,self.vy,self.vz)+\
             "Angles:   (%f,%f,%f)\n"%(self.phi,self.theta,self.gamma)+\
             "AngVelos: (%f,%f,%f)\n"%(self.phidot,self.thetadot,self.gammadot)
    return outstr


  def initialize_model(self,*args):#PL0,PLa,PD0,PDa,PTya,PTywy,PTy0,PTxwx,PTxwz,PTzwz):
    """
    Used to create a model for the forces and
    torques on this frisbee.
    """
    if len(args) == 1 and isinstance(args[0],np.ndarray):
      self.coefficients=coefficient_model.Model(args[0])
    elif len(args) == 10:
      PL0,PLa,PD0,PDa,PTxwx,PTxwz,PTy0,PTya,PTywy,PTzwz = args
      self.coefficients=coefficient_model.Model(PL0,PLa,PD0,PDa,PTxwx,PTxwz,PTy0,PTya,PTywy,PTzwz)
    else: raise Exception("Usage error: Model initialized incorrectly.")
    return

  def update_coordinates(self,coordinates):
    """
    Given some new coordinates, update the coordinates
    """
    x,y,z,vx,vy,vz,phi,theta,gamma,phidot,thetadot,gammadot = coordinates
    self.x=x 
    self.y=y
    self.z=z
    self.vx=vx
    self.vy=vy
    self.vz=vz
    self.phi=phi
    self.theta=theta
    self.gamma=gamma
    self.phidot=phidot
    self.thetadot=thetadot
    self.gammadot=gammadot
    self.update_data_fields()
    return


  def get_coordinates(self):
    """
    Return the current coordinates of the frisbee.
    """
    return np.array([self.x,self.y,self.z,self.vx,self.vy,self.vz,
                     self.phi,self.theta,self.gamma,
                     self.phidot,self.thetadot,self.gammadot])


  def update_data_fields(self):
    """
    Update the data fields in the frisbee.
    """
    self.calc_trig_functions()
    self.velocity = np.array([self.vx,self.vy,self.vz])
    self.angle_dots = np.array([self.phidot,self.thetadot,self.gammadot])
    self.vhat = self.velocity/np.linalg.norm(self.velocity)
    self.v2 = np.dot(self.velocity,self.velocity)
    self.rotation_matrix = self.calc_rotation_matrix()
    self.xbhat,self.ybhat,self.zbhat = self.calc_body_hat_vectors()
    self.angle_of_attack = self.calc_angle_of_attack()
    self.angular_velocity_frisframe = self.calc_angular_velocity_frisframe()
    self.angular_velocity_labframe = np.dot(self.angular_velocity_frisframe,self.rotation_matrix)
    self.wxb,self.wyb,self.wzb = self.calc_angular_velocity()
    return


  def calc_trig_functions(self):
    """
    Calculates the trig functions
    of the euler angles of the frisbee.
    """
    self.sin_phi = np.sin(self.phi)
    self.cos_phi = np.cos(self.phi)
    self.sin_theta = np.sin(self.theta)
    self.cos_theta = np.cos(self.theta)
    return


  def calc_rotation_matrix(self):
    """
    Calculates the euler rotation matrix.
    R(phi,theta) = Ry(theta)Rx(phi)
    
    See https://en.wikipedia.org/wiki/Davenport_chained_rotations. 
    """
    sp,cp = self.sin_phi,self.cos_phi
    st,ct = self.sin_theta,self.cos_theta
    return np.array([[ct,sp*st,-st*cp],
                     [0,cp,sp],
                     [st,-sp*ct,cp*ct]])


  def calc_body_hat_vectors(self):
    """
    Calculates the unit (hat) vectors
    fixed to the disc (excluding spin)
    in terms of the lab frame.
    """
    v = self.velocity
    zbhat = self.rotation_matrix[2]
    v_dot_zbhat = np.dot(v,zbhat)
    v_in_plane = v - (zbhat*v_dot_zbhat)
    xbhat = v_in_plane/np.linalg.norm(v_in_plane)
    ybhat = np.cross(zbhat,xbhat)
    return [xbhat,ybhat,zbhat]


  def calc_angle_of_attack(self):
    """
    Calculates angle of attack (AOA).

    AOA is defined as between plane of the disc
    and the velocity vector.
    """
    v = self.velocity
    zbhat = self.zbhat
    v_dot_zbhat = np.dot(v,zbhat)
    v_in_plane = v - zbhat*v_dot_zbhat
    return -np.arctan(v_dot_zbhat/np.linalg.norm(v_in_plane))


  def calc_angular_velocity_frisframe(self):
    """
    Calculates the angular velocity as seen in the
    frisbee frame. This is \vec{w}, not \vec{w}_F.
    See Hummel 2003 page 34.

    Note: \vec{w} \dot R gives the angular
    velocity in the lab frame.
    """
    st,ct = self.sin_theta,self.cos_theta
    return np.array([self.phidot*ct,
                     self.thetadot,
                     self.phidot*st + self.gammadot])


  def calc_angular_velocity(self):
    """
    Calculates the angular velocity along the
    body unit vectors as expressed in the lab frame.
    """
    av_labframe = self.angular_velocity_labframe
    xbhat,ybhat,zbhat = self.xbhat,self.ybhat,self.zbhat
    wxb = np.dot(av_labframe,xbhat)
    wyb = np.dot(av_labframe,ybhat)
    wzb = np.dot(av_labframe,zbhat)
    return [wxb,wyb,wzb]


  def get_acceleration(self):
    """
    Calculate acceleration of the positions.
    """
    alpha,v2 = self.angle_of_attack,self.v2
    vhat,ybhat = self.vhat,self.ybhat
    force_amplitude = 0.5*rho*area*v2
    C_lift =  self.coefficients.C_lift(alpha)
    C_drag =  self.coefficients.C_drag(alpha)
    F_lift = C_lift*force_amplitude*np.cross(vhat,ybhat)
    F_drag = C_drag*force_amplitude*(-vhat)
    total_force=F_lift+F_drag+F_gravity
    if self.debug:
      print "In get_acceleration:"
      print "\tC_lift:",C_lift
      print "\tC_drag:",C_drag
      print "\tAmplitude:",force_amplitude
      print "\tF_lift:",F_lift
      print "\tF_drag:",F_drag
      print "\tF_grav:",F_gravity
    return total_force/mass


  def get_torque(self):
    """
    Calculates the torques (moments) on the frisbee.
    """
    alpha=self.angle_of_attack
    v2=self.v2
    torque_amplitude = 0.5*rho*diameter*area*v2
    wxb,wyb,wzb = self.wxb,self.wyb,self.wzb
    C_x = self.coefficients.C_x(wxb,wzb)
    C_y = self.coefficients.C_y(alpha,wyb)
    C_z = self.coefficients.C_z(wzb)
    torque_x = C_x*torque_amplitude*self.xbhat
    torque_y = C_y*torque_amplitude*self.ybhat
    torque_z = C_z*torque_amplitude*np.array([0,0,1.])
    #Rotate into the lab frame.
    total_torque=np.dot(self.rotation_matrix,torque_x+torque_y)+torque_z
    #Optional: turn off rotations
    #total_torque *= 0
    if self.debug:
      print "In get_torque"
      print "\tAmplitude:",torque_amplitude
      xbhat,ybhat,zbhat = self.xbhat,self.ybhat,self.zbhat
      avf = self.angular_velocity_frisframe
      avl = self.angular_velocity_labframe
      R = self.rotation_matrix
      print "\tavf: ",avf
      print "\tC1:",R[0]
      print "\tC2:",R[1]
      print "\tC3:",R[2]
      print "\tavl: ",avl
      print "\txbhat = ",xbhat
      print "\tybhat = ",ybhat
      print "\tzbhat = ",zbhat
      print "\twx = %.2e\twy = %.2e\twz = %.2e"%(wxb,wyb,wzb)
      print "\tRoll amp:",C_x*torque_amplitude
      print "\tPitch amp:",C_y*torque_amplitude
      print "\tSpin amp:",C_z*torque_amplitude
      print "\ttotal_torque:",total_torque
    return total_torque


  def ang_acceleration(self):
    """
    Calculate angular accelerations in radians/s^2.
    
    See, e.g. Hummel 2003.
    """
    total_torque = self.get_torque()
    st,ct = self.sin_theta,self.cos_theta
    phidot,thetadot,gammadot = self.phidot,self.thetadot,self.gammadot

    phi_dd = (total_torque[0]+2*Ixy*thetadot*phidot*st-Izz*thetadot*(phidot*st+gammadot))/(ct*Ixy)
    theta_dd = (total_torque[1]+Izz*phidot*ct*(phidot*st+gammadot)-Ixy*phidot*phidot*ct*st)/Ixy
    gamma_dd = (total_torque[2]-Izz*(phidot*thetadot*ct+phi_dd*st))/Izz

    if self.debug:
      print "In ang_acceleration:"
      print "\tphi_dd:",2*Ixy*thetadot*phidot*st,Izz*thetadot*(phidot*st+gammadot)
      print "\ttheta_dd:",Izz*phidot*ct*(phidot*st+gammadot),Ixy*phidot**2*ct*st
      print "\tgamma_dd:",-Izz*(phidot*thetadot*ct+phi_dd*st)
    return np.array([phi_dd,theta_dd,gamma_dd])


  def derivatives_array(self):
    """
    Compute the derivatives of all coordinates.
    """
    derivatives = np.zeros(12)
    derivatives[0:3] = self.velocity
    derivatives[3:6] = self.get_acceleration()
    derivatives[6:9] = self.angle_dots
    derivatives[9:12]= self.ang_acceleration()
    if self.debug:
      print "In derivatives_array:"
      print "\tvelocities: ",derivatives[0:3]
      print "\tforces/m: ",derivatives[3:6]
      print "\tangle dots: ",derivatives[6:9]
      print "\tang_accs: ",derivatives[9:12]
    return derivatives


  def equations_of_motion(self,coordinates,t):
    """
    Return the equations of motion.
    For use with scipy integrators.
    """
    self.update_coordinates(coordinates)
    if self.z <= 0.0: return np.zeros_like(coordinates)
    return  self.derivatives_array()

  def get_trajectory(self,time_initial,time_final,dt=0.01):
    """
    Get a frisbee's trajectory give an initial
    and final time. The timestep size can be specified.
    This requires that the frisbee hass been properly
    initialized with a model.
    """
    coordinates = self.get_coordinates()
    flight_time = time_final-time_initial
    N_times = int(flight_time/dt)
    if self.use_C:
      if c_frisbee_available:
        return C_frisbee.get_trajectory(coordinates, flight_time, N_times, self.coefficients.get_model_as_array())
      else: raise Exception("Attempted to use C_frisbee without compiling.")
    times = np.linspace(time_initial,time_final,N_times)

    return times,odeint(self.equations_of_motion,coordinates,times)
    

#An example of initializing, printing, and calling a function of the frisbee
if __name__ == "__main__":
  #Using a single frisbee
  x,y,z = 0.0, 0.0, 1.0
  vx,vy,vz = 10.0,0.0,0.0
  phi,theta,gamma = 0.0,0.0,0.0
  phidot,thetadot,gammadot = 0.0,0.0,50.0
  test_frisbee = Frisbee(x,y,z,
                         vx,vy,vz,
                         phi,theta,gamma,
                         phidot,thetadot,gammadot,debug=True)
  test_frisbee.initialize_model(0.33,1.9,0.18,0.69,-1.3e-2,-1.7e-3,-8.2e-2,0.43,-1.4e-2,-3.4e-5)
  print test_frisbee
  test_frisbee.derivatives_array()
  print "\n"

  #Integrating a frisbee's equations of motion
  from scipy.integrate import odeint
  test_frisbee = Frisbee(x,y,z,
                         vx,vy,vz,
                         phi,theta,gamma,
                         phidot,thetadot,gammadot,use_C=False,debug=False)
  model = np.array([0.33,1.9,0.18,0.69,-1.3e-2,-1.7e-3,-8.2e-2,0.43,-1.4e-2,-3.4e-5])
  test_frisbee.initialize_model(model)
  times,trajectory = test_frisbee.get_trajectory(0.0,3.0,dt=0.1)
  print "Integrating the equations of motion at 30 time steps gives"
  print "all 12 kinematic variables at 30 times:",trajectory.shape
  
