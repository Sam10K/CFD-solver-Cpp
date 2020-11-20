#ifndef Solver_params_H
#define Solver_params_H

#include "headers.hpp"


class Solver_params
{
public:
  char steady,write_files;
  double mu,rho,T,CFL,write_start,write_end;
  int write_freq; // In terms of number of iterations/time_steps

  void Read_Solver_setup();
};

#endif
