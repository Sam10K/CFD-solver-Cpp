#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <cstring>
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include <vector>
#include <numeric>
#include <time.h>
#include <sys/stat.h>

using namespace Eigen;

using namespace std;

extern string file_solver_setup;
extern string file_schemes;
extern string file_mesh;
extern string file_bc_u;
extern string file_bc_p;
extern string file_SIMPLE_params;
