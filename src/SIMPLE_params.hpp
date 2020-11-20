#ifndef SIMPLE_params_H
#define SIMPLE_params_H

#include "headers.hpp"

class SIMPLE_params
{
public:
  double alpha_p,alpha_u,u_tol,v_tol,p_tol;
  int max_iter;

  void Read_SIMPLE_params();
};

#endif
