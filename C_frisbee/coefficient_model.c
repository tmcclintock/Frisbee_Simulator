/*This contains the coefficient functions.
 */
#include "coefficient_model.h"

double C_lift(double alpha, double PL0, double PLa){
  /*Calculate the lift coefficient.

    Args:
        alpha: angle of attack
	PL0: C_lift value at no angle of attack
	PLa: linear dependencs of alpha

    Returns:
        C_lift: lift coefficient
   */
  return PL0 + PLa*alpha;
}

double C_drag(double alpha, double PD0, double PDa, double alpha0){
  /*Calculate the drag coefficient.

    Args:
        alpha: angle of attack
	PD0: minimum value of C_drag
	PDa: quadratic dependence of alpha
	alpha0: angle of minimum drag

    Returns:
        C_drag: drag coefficient
   */
  return PD0 + PDa*(alpha-alpha0)*(alpha-alpha0);
}

double C_x(double wx, double wz, double PTxwx, double PTxwz){
  /*Calculate the x-torque coefficient.

    Args:
        wx: x-angular velocity
	wz: z-angular velocity
	PTxwx: linear dependence of wx
	PTxwz: linear dependence of wz

    Returns:
        C_x: x-torque coefficient
   */
  return wx*PTxwx + wz*PTxwz;
}

double C_y(double alpha, double wy, double PTy0, double PTywy, double PTya){
  /*Calculate the y-torque coefficient.

    Args:
        alpha: angle of attack
        wy: y-angular velocity
	PTy0: C_y value at no angle of attack and no y-angular velocity
	PTywy: linear dependence of wy
	PTya: linear dependence of alpha

    Returns:
        C_y: y-torque coefficient
   */
  return PTy0 + alpha*PTya + wy*PTywy;
}

double C_z(double wz, double PTzwz){
  /*Calculate the z-torque coefficient.

    Args:
        wz: z-angular velocity
	PTzwz: linear dependence of wz

    Returns:
        C_z: z-torque coefficient
   */
  return wz*PTzwz;
}
