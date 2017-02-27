#include <stdio.h>
#include "coefficient_model.h"

int main(){
  printf("Testing the connection to each c file.\n");
  double alpha = 3.1;
  double PL0 = 1.0;
  double PLa = 2.0;
  printf("\tC_lift(%.1f,%.1f,%.1f) = %.2f\n",alpha,PL0,PLa,C_lift(alpha,PL0,PLa));
}
