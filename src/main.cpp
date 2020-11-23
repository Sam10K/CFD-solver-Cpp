#include "headers.hpp"
#include "Schemes.hpp"
#include "Solver_setup.hpp"
#include "SIMPLE_params.hpp"
#include "Mesh.hpp"
#include "Pressure_field.hpp"
#include "Velocity_field.hpp"
#include "Momentum.hpp"
#include "Pressure_poisson.hpp"
#include "Post_processing.hpp"


string file_solver_setup = "./Solver_setup/solver_setup";
string file_schemes = "./Solver_setup/Schemes";
string file_mesh = "./Mesh_file/mesh.msh";
string file_bc_p = "./Boundary_conditions/P.bc";
string file_bc_u = "./Boundary_conditions/U.bc";
string file_SIMPLE_params = "./Solver_setup/SIMPLE_params";


int main()
{
  Solver_params Solver_params;
  Schemes Scheme_params;
  SIMPLE_params SIMPLE_params;
  Mesh mesh;

  mesh.Read_Mesh();
  mesh.Read_P_BC();
  mesh.Read_U_BC();

  Scheme_params.Read_Schemes();
  Solver_params.Read_Solver_setup();
  SIMPLE_params.Read_SIMPLE_params();

  Pressure P(Solver_params,Scheme_params,SIMPLE_params,mesh);
  Velocity U(P);

  Pressure Pdash=P;

  VectorXd zeros(mesh.n_cells);zeros.setZero();

  double ures,vres,pres;
  double utol = SIMPLE_params.u_tol;
  double vtol = SIMPLE_params.v_tol;
  double ptol = SIMPLE_params.p_tol;

  cout<<"Solver started.............."<<endl<<endl;

  for(int i=0;i<Solver_params.T;i++)
  {
    cout<<"Iteration: "<<i+1<<endl;
    clock_t t;t = clock();
    Momentum_eq momentum(P,U);

    U = momentum.matrix_solve();

    Pdash.set_zero();//Pdash.set_c(zeros);
    Pressure_poisson_eq pressure_poisson(Pdash,U);
    Pdash = pressure_poisson.matrix_solve();

    MatrixXd vol(mesh.n_cells,2);vol<<mesh.volume,mesh.volume;
    MatrixXd vel_dash = -Pdash.grad_p_c.cwiseProduct(vol.cwiseQuotient(U.Ap));
    MatrixXd new_vel = U.get_vel()+vel_dash;
    VectorXd new_p = P.get_pressure()+SIMPLE_params.alpha_p*Pdash.get_pressure();
    P.set_c(new_p);
    U.P = P;
    U.set_c(new_vel);


    ures = vel_dash.col(0).norm()/new_vel.col(0).norm();
    vres = vel_dash.col(1).norm()/new_vel.col(1).norm();
    pres = Pdash.get_pressure().norm()/new_p.norm();

    t = clock() - t;cout<<"Time elapsed for one complete iteration : "<<((double)t)/CLOCKS_PER_SEC<<" s"<<endl;

    cout<<"SIMPLE_Ures: "<<ures<<" SIMPLE_Vres: "<<vres<<" SIMPLE_Pres: "<<pres<<endl;

    cout<<"Volume averaged continuity: "<<U.get_continuity()<<endl<<endl;

    if((Solver_params.write_files=='T')&&(i+1>=Solver_params.write_start)&&(i+1<=Solver_params.write_end))
    {
      if(fmod(i+1-Solver_params.write_start, Solver_params.write_freq) == 0)
      {
        cout<<"Writing VTK data........."<<endl;
        Post_process post(P,U);
        post.writeVTK(i+1);
        cout<<"Done"<<endl<<endl;
      }
    }


    if((ures<utol)&&(vres<vtol)&&(pres<ptol))
    {
      cout<<"Solution converged !!!"<<endl<<endl;
      cout<<"Writing VTK data........."<<endl;
      Post_process post(P,U);
      post.writeVTK(i+1);
      cout<<"Done"<<endl<<endl;

      break;
    }

    if(isnan(ures)||isnan(vres)||isnan(pres))
    {
      cout<<"Divergence detected. Exiting ........"<<endl;
      exit(0);
    }

    if(i==Solver_params.T-1)
    {
      cout<<"Maximum number of SIMPLE iterations reached"<<endl<<endl;
      cout<<"Writing VTK data........."<<endl;
      Post_process post(P,U);
      post.writeVTK(i+1);
      cout<<"Done"<<endl<<endl;

    }

  }


  return 0;
}
