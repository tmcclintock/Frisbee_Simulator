#include <stdio.h>
#include <math.h>

#include "equations_of_motion.h"
#include "coefficient_model.h"

#define pi M_PI //pi
#define m 0.175 //Frisbee mass; kg
#define g 9.81  //acceleration due to gravity; m/s/s
#define Area 0.057 //planform area of a frisbee (top down); m^2
#define d 0.269396834 //diameter of frisbee; m
#define rho 1.225 //average density of air; kg/m^3
#define Izz 0.002352 //moment of inertia about ZZ; kg*m^2
#define Ixy 0.001219 //moment of inertia about either XX or YY; kg*m^2

/*Units convention (unless otherwise noted in the code):
  - all Euclidean positions are in meters
  - all angles are in radians
  - times are all in seconds
  - masses are all in kilograms
*/

/*Naming convention:
  - any variable with a 'Dot' at the end means it is a time derivative
  - 'positions' refers to x, y, z as well as the euler angles in addition
  to the time derivatives of all of them
  - 'derivs' are the time derivatives of all of the positions
  - 't' has the current time, which is not needed in this calculation,
  but could in principle account for a time-varying wind
  - 'params' contains all of the parameters we are fitting for. That is,
  all of the flight parameters that are unknown
  - 'tau' refers to a torque (aka moment) along one of the inertial axes
*/

double dot(double*A,double*B){
  /*Take the dot product of two vectors.

    Args:
        A: a vector
	B: a vector

    Returns:
        A dot B: the dot product
   */
  return A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
}

void cross(double*A,double*B,double*C){
  /*Given two vectors A and B, calculate the cross product and store it in C.

    Note: Uses a temp array so that cross product can be done in place

    Args:
        A: a vector
	B: a vector
	C: a vector that contains the cross product
   */
  double D[3];
  D[0] = A[1]*B[2]-A[2]*B[1];
  D[1] = A[2]*B[0]-A[0]*B[2];
  D[2] = A[0]*B[1]-A[1]*B[0];
  C[0]=D[0],C[1]=D[1],C[2]=D[2];
  return;
}

void add(double*A,double *B,double *C){
  /*Given two vectors A and B, add them and store it in C.
    
    Note: Uses a temp array so that cross produces can be done in place

    Args:
        A: a vector
	B: a vector
	C: a vector that contains the sum
   */
  double D[3];
  D[0] = A[0]+B[0];
  D[1] = A[1]+B[1];
  D[2] = A[2]+B[2];
  C[0]=D[0],C[1]=D[1],C[2]=D[2];
  return;
}

void print_vector(double*A){
  printf("A[0]=%.3e\tA[1]=%.3e\tA[2]=%.3e\n",A[0],A[1],A[2]);
  //printf("A[0]=%e\tA[1]=%e\tA[2]=%e\tlen=%e\n",A[0],A[1],A[2],sqrt(dot(A,A)));
}

void equations_of_motion(double*positions,double*derivs,double t,double*params){
  /*The equations of motion of a frisbee.

    Args:
        positions: current value of all coordinates and velocities
	derivs: array used to contain the derivatives of the positions
	t: current time
	params: all of the parameters that control the scaling of the forces and torques on the disc. These are the quantities of interest.
   */
  //Declare some temporary variables that may be used
  double temp;
  int i,j,k;

  /*The six coordinates and velocities in the lab frame
    arrive in the arrays 'positions'.
   */
  double x = positions[0], y = positions[1], z = positions[2];
  double vx= positions[3], vy= positions[4], vz= positions[5];
  double phi     = positions[6], theta    = positions[7];
  double gamma = positions[8];
  double phiDot  = positions[9], thetaDot = positions[10];
  double gammaDot= positions[11];

  //Take in the various force and moment coefficients
  double PL0 = params[0], PLa = params[1];
  double PD0 = params[2], PDa = params[3];
  double PTxwx = params[4], PTxwz = params[5];
  double PTy0 = params[6], PTya = params[7], PTywy = params[8];
  double PTzwz = params[9];

  //This is the angle of attack that with minimal drag
  double alpha0 = 4*M_PI/180.;

  //Check if the disc has hit the ground
  if (z <= 0){
    for(i=0; i<12; i++){
      derivs[i] = 0;
    }
    return;
  }

  //TODO: include a wind. Probably as a static set of variables.
  /*
    vx = vx + v_wind_x;
    vy = vy + v_wind_y;
    vz = vz + v_wind_z;
  */

  //Calculate the sine and cosine of theta and phi to be used later
  //This is a convenient short hand
  double s_phi = sin(phi),  c_phi = cos(phi);
  double s_tht = sin(theta),c_tht = cos(theta);

  //Construct the Euler rotation matrix to go from 
  //the inertial frame (N) to the body frame (C)
  //NOTE: Tnc takes you from N to C
  //while {Tnc}^T = Tcn takes you from C to N
  double Tnc[3][3];
  Tnc[0][0] = c_tht, Tnc[0][1] = s_phi*s_tht,  Tnc[0][2] = -c_phi*s_tht;
  Tnc[1][0] = 0,     Tnc[1][1] = c_phi,        Tnc[1][2] = s_phi;
  Tnc[2][0] = s_tht, Tnc[2][1] = -s_phi*c_tht, Tnc[2][2] = c_phi*c_tht;

  //Define the rotation matrix rows
  double C1[] = {Tnc[0][0],Tnc[0][1],Tnc[0][2]};
  double C2[] = {Tnc[1][0],Tnc[1][1],Tnc[1][2]};
  double C3[] = {Tnc[2][0],Tnc[2][1],Tnc[2][2]};

  //Find the velocity components in the body frame
  double v[]={vx,vy,vz};
  double v_dot_C3 = dot(v,C3);
  double v_plane[] = {v[0]-C3[0]*v_dot_C3, v[1]-C3[1]*v_dot_C3, v[2]-C3[2]*v_dot_C3};

  //Find the length of these various velocities
  double norm_v_plane = sqrt(dot(v_plane,v_plane));
  double norm_v = sqrt(dot(v,v));
  //double norm_vp      = sqrt(dot(v_plane,v_plane));

  //Find the unit vectors for various directions
  double v_hat[] = {v[0]/norm_v,v[1]/norm_v,v[2]/norm_v};
  double vp_hat[] = {v_plane[0]/norm_v_plane,v_plane[1]/norm_v_plane,v_plane[2]/norm_v_plane};
  double ulat[3];
  cross(C3,vp_hat,ulat);

  //Calculate the angle of attack, alpha
  double alpha = -atan(v_dot_C3/norm_v_plane);

  //Make copies of some of the arrays, mostly so that it is easier to see what is happening
  double x_C_hat[] = {vp_hat[0],vp_hat[1],vp_hat[2]};
  double y_C_hat[] = {ulat[0],ulat[1],ulat[2]};
  double z_C_hat[] = {C3[0],C3[1],C3[2]};

  //Calculate the lift force
  double aerodynamic_Force_amp = 0.5*rho*Area*dot(v,v);
  double Flift_amp = C_lift(alpha,PL0,PLa)*aerodynamic_Force_amp;
  double Flift[3];
  cross(v_hat,y_C_hat,Flift);
  Flift[0] = Flift[0]*Flift_amp; Flift[1] = Flift[1]*Flift_amp; Flift[2] = Flift[2]*Flift_amp;

  //Calculate the drag force
  double Fdrag_amp = C_drag(alpha,PD0,PDa,alpha0)*aerodynamic_Force_amp;
  double Fdrag[]={-Fdrag_amp*v_hat[0],-Fdrag_amp*v_hat[1],-Fdrag_amp*v_hat[2]};

  //Calculate the gravitational force
  double Fgrav[]={0,0,-m*g};

  //Sum the forces to get the total force
  double Ftot[3];
  add(Flift,Fdrag,Ftot);
  add(Ftot,Fgrav,Ftot);

  //Write the angular velocities in the frisbee frame
  double w_in_C[] = {phiDot*c_tht, thetaDot, phiDot*s_tht + gammaDot};

  //Calculate the angular velocity in the lab frame
  double w_in_N[]={0,0,0};
  for(i=0;i<3;i++){
    for(j=0;j<3;j++){// w_in_C dot Tnc
      w_in_N[i] +=  w_in_C[j] * Tnc[j][i];
    }
  }

  //Find the components of the angular velocity
  double w_x = dot(w_in_N,x_C_hat);
  double w_y = dot(w_in_N,y_C_hat);
  double w_z = dot(w_in_N,z_C_hat);

  //Calculate the spin parameter, which could be used later
  //to incorporate the Robins-Magnus force or if we
  //have more accurate data about drag coefficients
  double spin_parameter = w_z*d/(2*norm_v);

  //Calculate each torque (aka moment)
  double tau_amp = 0.5*rho*d*Area*dot(v,v);
  double tau_x_amp = C_x(w_x,w_z,PTxwx,PTxwz)*tau_amp;
  double tau_y_amp = C_y(alpha,w_y,PTy0,PTywy,PTya)*tau_amp;
  double tau_z_amp = C_z(w_z,PTzwz)*tau_amp;

  double tau_x[] = {x_C_hat[0]*tau_x_amp, x_C_hat[1]*tau_x_amp, x_C_hat[2]*tau_x_amp};
  double tau_y[] = {y_C_hat[0]*tau_y_amp, y_C_hat[1]*tau_y_amp, y_C_hat[2]*tau_y_amp};
  double tau_z_N[] = {0, 0, tau_z_amp}; //Already in the correct frame

  //Calculate the total torque
  double tau_total[]={0,0,0};
  for(i=0;i<3;i++)
    {
      for(j=0;j<3;j++)
	{//Tnc dot tau_x+y
	  tau_total[i] += Tnc[i][j] * (tau_x[j] + tau_y[j]);
	}
    }
  add(tau_total, tau_z_N, tau_total);

  //If we want, set the torques all to 0
  //tau_total[0] = 0, tau_total[1] = 0, tau_total[2] = 0;

  //TODO: Remove the wind from the velocities
  /*
    vx = vx + v_wind_x;
    vy = vy + v_wind_y;
    vz = vz + v_wind_z;
  */

  //Calculate all of the derivatives
  derivs[0] = vx;
  derivs[1] = vy;
  derivs[2] = vz;
  derivs[3] = Ftot[0]/m;
  derivs[4] = Ftot[1]/m;
  derivs[5] = Ftot[2]/m;
  derivs[6] = phiDot;
  derivs[7] = thetaDot;
  derivs[8] = gammaDot;
  derivs[9] = (tau_total[0]+2*Ixy*thetaDot*phiDot*s_tht-Izz*thetaDot*(phiDot*s_tht+gammaDot))/(Ixy*c_tht);
  derivs[10] = (tau_total[1]+Izz*phiDot*c_tht*(phiDot*s_tht+gammaDot)-Ixy*phiDot*phiDot*c_tht*s_tht)/Ixy;
  derivs[11]= (tau_total[2]-Izz*(phiDot*thetaDot*c_tht + derivs[9]*s_tht))/Izz;

  if (t<0.0){
    printf("\nCCCCC\n\tC_lift = %f\n",C_lift(alpha,PL0,PLa));
    printf("\tC_drag = %f\n",C_drag(alpha,PD0,PDa,alpha0));
    printf("\tAero Amplitude = %f\n",aerodynamic_Force_amp);
    printf("\tF_lift");print_vector(Flift);
    printf("\tF_drag");print_vector(Fdrag);
    printf("\tF_grav");print_vector(Fgrav);
    printf("\tTorque Amplitude = %f\n",tau_amp);
    printf("\tavf:");print_vector(w_in_C);
    printf("\tC1:");print_vector(C1);
    printf("\tC2:");print_vector(C2);
    printf("\tC3:");print_vector(C3);
    printf("\tavl:");print_vector(w_in_N);
    //printf("\txbhat");print_vector(x_C_hat);
    //printf("\tybhat");print_vector(y_C_hat);
    //printf("\tzbhat");print_vector(z_C_hat);
    printf("\twx = %.3e\twy = %.3e\twz = %.3e\n",w_x,w_y,w_z);
    printf("\tRoll amp = %f\n",tau_x_amp);
    printf("\tPitch amp = %f\n",tau_y_amp);
    printf("\tSpin amp = %f\n",tau_z_amp);
    printf("\ttotal torque");print_vector(tau_total);
    printf("\tphi_dd: %e\n",2*Ixy*thetaDot*phiDot*s_tht-Izz*thetaDot*(phiDot*s_tht+gammaDot));
    printf("\ttheta_dd: %e\n",Izz*phiDot*c_tht*(phiDot*s_tht+gammaDot)-Ixy*phiDot*phiDot*c_tht*s_tht);
    printf("\tgamma_dd: %e\n",-Izz*(phiDot*thetaDot*c_tht + derivs[9]*s_tht));

    for(i=0;i<12;i++)
      printf("cderiv[%d] = %f\n",i,derivs[i]);
    //for(i=0;i<10;i++)
    //  printf("cparams[%d] = %f\n",i,params[i]);

    //printf("t= %f\td6,d7,d11: %f, %f, %f\n",t,derivs[6],derivs[7],derivs[11]);
    //printf("t=%f\td8,d9,d10: %f, %f, %f\n",t,derivs[8],derivs[9],derivs[10]);
    //printf("gd = %f\n",gammaDot);
  }
  //End the function
  return;
}
