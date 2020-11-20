#ifndef Post_processing_H
#define Post_processing_H

#include "headers.hpp"
#include "Schemes.hpp"
#include "Solver_setup.hpp"
#include "SIMPLE_params.hpp"
#include "Mesh.hpp"
#include "Pressure_field.hpp"
#include "Velocity_field.hpp"



class Post_process
{
  Mesh mesh;
  Solver_params solver_setup;
  Schemes schemes;
  SIMPLE_params SIMPLE;
  Pressure P;
  Velocity U;

public:

  Post_process(Pressure p,Velocity u);

  void writeVTK(int folder); //Create folder with iter number

};

#endif
