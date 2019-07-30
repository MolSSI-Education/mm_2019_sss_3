#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <math.h>


double LJ(Eigen::VectorXd coord_ij2, 
            double sigma = 1.0, 
            double epsilon = 1.0,
            double cutoff2 = 2.6);

double cutoff_correction(Eigen::Vector3d box_dims,
            int num_particles, 
            double sigma = 1.0, 
            double epsilon = 1.0, 
            double cutoff = 2.6);