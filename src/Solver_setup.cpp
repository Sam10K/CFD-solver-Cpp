#include "Solver_setup.hpp"

void Solver_params::Read_Solver_setup()
{
  ifstream fin;
  string line;
  fin.open(file_solver_setup);
  cout<<"Reading Solver_setup file......."<<endl;
  if(!fin)
  {
    cout<<"Solver_setup file not found. Exiting .........\n";
    exit(0);
  }

  while(fin)
  {
    getline(fin,line);
    char buffer[line.length()+1];
    strcpy(buffer,line.c_str());

    char c;
    double x;
    int y;

    if(sscanf(buffer, "steady = %c;", &c) && strlen(buffer)!=0)
    {
      steady = c;
    }

    if(sscanf(buffer, "mu = %lf;", &x) && strlen(buffer)!=0)
    {
      mu = x;
    }

    if(sscanf(buffer, "rho = %lf;", &x) && strlen(buffer)!=0)
    {
      rho = x;
    }

    if(sscanf(buffer, "T = %lf;", &x) && strlen(buffer)!=0)
    {
      T = x;
    }

    if(sscanf(buffer, "CFL = %lf;", &x) && strlen(buffer)!=0)
    {
      CFL = x;
    }

    if(sscanf(buffer, "write_files = %c;", &c) && strlen(buffer)!=0)
    {
      write_files = c;
    }

    if(sscanf(buffer, "write_start = %lf;", &x) && strlen(buffer)!=0)
    {
      write_start = x;
    }

    if(sscanf(buffer, "write_end = %lf;", &x) && strlen(buffer)!=0)
    {
      write_end = x;
    }

    if(sscanf(buffer, "write_freq = %d;", &y) && strlen(buffer)!=0)
    {
      write_freq = y;
    }


  }

  fin.close();
  cout<<"Done"<<endl<<endl;

}
