#include "Schemes.hpp"

void Schemes::Read_Schemes()
{
  ifstream fin;
  string line;
  fin.open(file_schemes);
  cout<<"Reading Schemes file ..............."<<endl;
  if(!fin)
  {
    cout<<"Schemes file not found. Exiting .........\n";
    exit(0);
  }

  //Schemes Schemes;

  while(fin)
  {
    getline(fin,line);
    char buffer[line.length()+1];
    strcpy(buffer,line.c_str());

    char c[10];
    double x;

    if(sscanf(buffer, "div_scheme = %[^;]", c) && strlen(buffer)!=0)
    {
      string str = c;
      div_scheme = str;
    }

    if(sscanf(buffer, "grad_scheme = %[^;]", c) && strlen(buffer)!=0)
    {
      string str = c;
      grad_scheme = str;
    }

    if(sscanf(buffer, "u_matrix_solver = %[^;]", c) && strlen(buffer)!=0)
    {
      string str = c;
      u_matrix_solver.solver = str;
      while(fin)
      {

        getline(fin,line);
        char buffer[line.length()+1];
        strcpy(buffer,line.c_str());

        int iter,nsweep;
        double tol;

        if(sscanf(buffer, "max_iter = %d;", &iter) && strlen(buffer)!=0)
        {
          u_matrix_solver.max_iter = iter;
        }

        if(sscanf(buffer, "tol = %lf;", &tol) && strlen(buffer)!=0)
        {
          u_matrix_solver.tol = tol;
        }

        if(sscanf(buffer, "nsweeps = %d;", &nsweep) && strlen(buffer)!=0)
        {
          u_matrix_solver.nsweeps = nsweep;
        }

        if(strcmp(buffer,"};")==0)
        {
          break;
        }

      }


    }

    if(sscanf(buffer, "v_matrix_solver = %[^;]", c) && strlen(buffer)!=0)
    {
      string str = c;
      v_matrix_solver.solver = str;
      while(fin)
      {

        getline(fin,line);
        char buffer[line.length()+1];
        strcpy(buffer,line.c_str());

        int iter,nsweep;
        double tol;

        if(sscanf(buffer, "max_iter = %d;", &iter) && strlen(buffer)!=0)
        {
          v_matrix_solver.max_iter = iter;
        }

        if(sscanf(buffer, "tol = %lf;", &tol) && strlen(buffer)!=0)
        {
          v_matrix_solver.tol = tol;
        }

        if(sscanf(buffer, "nsweeps = %d;", &nsweep) && strlen(buffer)!=0)
        {
          v_matrix_solver.nsweeps = nsweep;
        }

        if(strcmp(buffer,"};")==0)
        {
          break;
        }

      }


    }

    if(sscanf(buffer, "p_matrix_solver = %[^;]", c) && strlen(buffer)!=0)
    {
      string str = c;
      p_matrix_solver.solver = str;
      while(fin)
      {

        getline(fin,line);
        char buffer[line.length()+1];
        strcpy(buffer,line.c_str());

        int iter,nsweep;
        double tol;

        if(sscanf(buffer, "max_iter = %d;", &iter) && strlen(buffer)!=0)
        {
          p_matrix_solver.max_iter = iter;
        }

        if(sscanf(buffer, "tol = %lf;", &tol) && strlen(buffer)!=0)
        {
          p_matrix_solver.tol = tol;
        }

        if(sscanf(buffer, "nsweeps = %d;", &nsweep) && strlen(buffer)!=0)
        {
          p_matrix_solver.nsweeps = nsweep;
        }

        if(strcmp(buffer,"};")==0)
        {
          break;
        }

      }


    }

  }


  fin.close();
  cout<<"Done"<<endl<<endl;

}
