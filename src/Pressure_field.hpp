#ifndef Pressure_field_H
#define Pressure_field_H

#include "headers.hpp"
#include "Schemes.hpp"
#include "Solver_setup.hpp"
#include "SIMPLE_params.hpp"
#include "Mesh.hpp"


class Pressure
{
  MatrixXd p_c,p_f;

  Mesh mesh;
  Solver_params solver_setup;
  Schemes schemes;
  SIMPLE_params SIMPLE;
public:
  MatrixXd grad_p_c,grad_p_f;

  Pressure()
  { } //Dummy constructor

  void set_zero() // To set everything to zero. Used for Pdash initialization
  {
    p_c.setZero();
    p_f.setZero();
    grad_p_c.setZero();
    grad_p_f.setZero();
  }

  Pressure(Solver_params,Schemes,SIMPLE_params,Mesh);

  void initialize();

  Mesh get_mesh()
  {
    return mesh;
  }

  Solver_params get_solver_params()
  {
    return solver_setup;
  }

  SIMPLE_params get_simple_params()
  {
    return SIMPLE;
  }

  Schemes get_schemes()
  {
    return schemes;
  }

  VectorXd get_pressure()
  {
    return p_c;
  }

  void set_c(VectorXd p);

  double get_c(int cell_id)
  {
    return p_c(cell_id);
  }

  double get_n(int node_id)
  {
    return calc_nodal_pressure(node_id);
  }

  double get_f(int face_id,int cell_id);

  Matrix<double,1,2> grad_c(int cell_id)
  {
    return grad_p_c.row(cell_id);
  }

  Matrix<double,1,2> grad_f(int face_id,int cell_id);

// Calculate facial pressure values
  double calc_p_f(int face_id,int cell_id);

// Calculate centroid pressure gradient values
  Matrix<double,1,2> calc_grad_p_c(int cell_id);

// Calculate facial pressure gradient values
  Matrix<double,1,2> calc_grad_p_f(int face_id,int cell_id);

// Calculate nodal pressure values
  double calc_nodal_pressure(int node_id);

};

#endif
