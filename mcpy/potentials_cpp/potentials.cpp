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
            double cutoff2 = 2.6)
{
    const size_t nthreads = std::thread::hardware_concurrency();

    double sigma2 = pow(sigma, 2);
    double energy = 0;
    {
    std::vector<std::thread> threads(nthreads);
    std::mutex critical;
    
    for(int t = 0;t<nthreads;t++)
    {
      threads[t] = std::thread(std::bind(
        [&](const int bi, const int ei, const int t)
        {

            for ( auto rij2 : coord_ij2 )
            {
                if (rij2 < cutoff2) {
                    double sig_by_r6 = pow( sigma2 / rij2, 3 );
                    double sig_by_r12 = pow( sig_by_r6  , 2);
                    energy += 4.0 * epsilon * (sig_by_r12 - sig_by_r6);
                }
            }
            },t*nloop/nthreads,(t+1)==nthreads?nloop:(t+1)*nloop/nthreads,t));
    }
    std::for_each(threads.begin(),threads.end(),[](std::thread& x){x.join();});

    return energy;
}

double cutoff_correction(double volume,
            int num_particles, 
            double sigma,
            double epsilon,
            double cutoff)
{
    //double volume = volume;

    double sig_by_cutoff3 = pow( sigma / cutoff, 3);
    double sig_by_cutoff9 = pow( sig_by_cutoff3, 3);

    double e_correction = sig_by_cutoff9 - ( 3.0 * sig_by_cutoff3 );
    double n_particles_2 = pow(num_particles, 2);

    e_correction *= 8.0 / 9.0 * M_PI * n_particles_2 * epsilon * pow(sigma, 3) / volume;

    return e_correction;

}