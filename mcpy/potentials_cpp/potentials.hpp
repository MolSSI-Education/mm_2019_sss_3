#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <math.h>
#include <thread>


double LJ(std::vector<double> & coord_ij2, 
            double sigma = 1.0, 
            double epsilon = 1.0,
            double cutoff2 = 2.6);

double cutoff_correction(double volume,
            int num_particles, 
            double sigma = 1.0, 
            double epsilon = 1.0, 
            double cutoff = 2.6);