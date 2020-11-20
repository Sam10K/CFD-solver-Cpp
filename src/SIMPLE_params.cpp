#include "SIMPLE_params.hpp"


void SIMPLE_params::Read_SIMPLE_params()
{
  ifstream fin;
  string line;
  fin.open(file_SIMPLE_params);
  cout<<"Reading SIMPLE_params file......."<<endl;
  if(!fin)
  {
    cout<<"SIMPLE_params file not found. Exiting .........\n";
    exit(0);
  }


  while(fin)
  {
    getline(fin,line);
    char buffer[line.length()+1];
    strcpy(buffer,line.c_str());

    double x;int xi;

    if(sscanf(buffer, "alpha_p = %lf;", &x) && strlen(buffer)!=0)
    {
      alpha_p = x;
    }

    if(sscanf(buffer, "alpha_u = %lf;", &x) && strlen(buffer)!=0)
    {
      alpha_u = x;
    }

    if(sscanf(buffer, "u_tol = %lf;", &x) && strlen(buffer)!=0)
    {
      u_tol = x;
    }

    if(sscanf(buffer, "v_tol = %lf;", &x) && strlen(buffer)!=0)
    {
      v_tol = x;
    }

    if(sscanf(buffer, "p_tol = %lf;", &x) && strlen(buffer)!=0)
    {
      p_tol = x;
    }

    if(sscanf(buffer, "Max_iter = %d;", &xi) && strlen(buffer)!=0)
    {
      max_iter = xi;
    }

  }

  fin.close();
  cout<<"Done"<<endl<<endl;

}
