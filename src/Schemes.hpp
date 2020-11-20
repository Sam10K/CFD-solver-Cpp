#ifndef Schemes_H
#define Schemes_H

#include "headers.hpp"

struct matrix_solver
{
  string solver;
  int max_iter,nsweeps;
  double tol;
};

class Schemes
{

public:
  string div_scheme,grad_scheme;

  matrix_solver u_matrix_solver,v_matrix_solver,p_matrix_solver;

  void Read_Schemes();


};

#endif
