#ifndef Velocity_field_H
#define Velocity_field_H

#include "headers.hpp"
#include "Schemes.hpp"
#include "Solver_setup.hpp"
#include "SIMPLE_params.hpp"
#include "Mesh.hpp"
#include "Pressure_field.hpp"



class Velocity
{
  MatrixXd velocity,grad_vel_c; //size is n_cells
  MatrixXd grad_vel_f,velocity_f; //size is n_faces. Contains values from both the straddling cells of the faces
  Mesh mesh;
  Solver_params solver_setup;
  Schemes schemes;
  SIMPLE_params SIMPLE;

public:
  MatrixXd Ap;
  Pressure P;

  Velocity() //Dummy constructor
  {

  }

  Velocity(Pressure p);

  void initialize();

  RowVectorXd face_interp_vel(int face_id,int cell_id); //To calculate the facail values without Rhie-Chow interpolation

  void set_c(MatrixXd vel);

  MatrixXd get_vel()
  {
    return velocity;
  }

  RowVectorXd get_c(int cell_id)
  {
    return velocity.row(cell_id);
  }

  RowVectorXd get_n(int node_id)
  {
    return calc_nodal_velocity(node_id);
  }

  RowVectorXd get_f(int face_id,int cell_id);

  MatrixXd grad_c(int cell_id);

  MatrixXd grad_f(int face_id,int cell_id);

  double get_continuity();


// Calculate facial velocity values
  RowVectorXd calc_vel_f(int face_id,int cell_id);

// Calculate centroid pressure gradient values
  MatrixXd calc_grad_vel_c(int cell_id);

// Calculate facial velocity gradient values
  MatrixXd calc_grad_vel_f(int face_id,int cell_id);

// Calculate continuity for every cell
  double calc_cell_continuity(int cell_id);

// Calculatre nodal velocity for every node
  RowVectorXd calc_nodal_velocity(int node_id);

};

#endif
