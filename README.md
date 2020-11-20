# CFD-solver-MATLAB
* A 2D Navier-Stokes solver for solving steady, laminar, incompressible flows using finite-volume method and collocated grid arrangement coded in C++
* Unsteady solver to be implemented soon
* Pressure-velocity coupling implemented using SIMPLE algorithm
* Spatial discretization for divergence schemes - Available choices include Upwind, Central differencing, Second order upwind, QUICK and FROMM schemes
* Cell centred gradient algorithms : Least squares gradient scheme
* Matrix solvers available : Gauss Siedel, BICGSTAB
* Accepts all-tri mesh as well as all-quad mesh in 2D ASCII Ansys-Fluent mesh file format (.msh)
* Files are output in .vtk ASCII format for visualizing in Paraview

## Instructions:

* The solver is coded on Ubuntu using g++ 7.5.0 and uses Eigen matrix library (https://eigen.tuxfamily.org/dox/group__QuickRefPage.html) to perform matrix operations
* To compile the solver on your PC, navigate to the downloaded folder (containing src, tutorials folders and Makefile) from a terminal window and run make. The make file is written assuming your .bashrc file is in your home directory
* If successfull, you should have a binary file called run within the same folder
* Now you can run the solver from anywhere on your PC using the command "Sam10K-solver". You have to run this command within the folder containing the solver parameters, just like OpenFOAM. An example of such a older can be seen in any of the ones available within the tutorial folder
* Some example mesh files and their boundary condition files are provided within the tutorial folder. Note that each mesh file is named as mesh.msh. Do not change this name. Similarly do not change the format or name of any other solver files
* Set the boundary conditions using the files 'U.bc' and 'P.bc' inside the folder named Boundary_conditions. Check the example boundary condition files. Currently supports fixed value and zero gradient boundary conditions
* You can edit the solver setup parameters using the file 'solver_setup' inside the solver_setup folder
* Schemes and SIMPLE algorithm parameters can alo be edited using the files inside the solver_setup folder  
* The simulation results are output within a folder called Data in .vtk ASCII format. You can use Paraview to visualize them

## References:
* Error Analysis and Estimation for the Finite Volume Method with Applications to Fluid Flows - PhD thesis by Dr. Hrvoje Jasak (https://foam-extend.fsb.hr/wp-content/uploads/2016/12/Jasak_PhD_1996.pdf)

* The Finite Volume Method in Computational Fluid Dynamics - F. Moukalled, L. Mangani, M. Darwish (https://link.springer.com/book/10.1007/978-3-319-16874-6)

* Computational Methods for Fluid Dynamics - Joel H. Ferziger, Milovan Perić (https://link.springer.com/book/10.1007/978-3-642-56026-2)
