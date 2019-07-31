#include "pairwise.hpp"
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <vector>
#include <math.h>


PairwisePotentials::PairwisePotentials(){}

LJ::LJ(double sigma, double epsilon, double cutoff)
        {
            m_sigma = sigma;
            m_epsilon = epsilon;
            m_cutoff = cutoff;
            m_cutoff2 = cutoff * cutoff;
        }

double LJ::potential(std::vector<double> rij2_array)
        {   
            /* Pairwise potential and correction energy by Lennard-Jones potential

            Parameters
            ----------

            sigma : float
                Distance between two particles when interaction potential is zero

            epsilon: float
                Depth of potential well
            */
            m_rij2_array = rij2_array;
            double sigma2 = pow(m_sigma, 2);
            double energy = 0;
            
            for ( auto rij2 : m_rij2_array )
            {
                if (rij2 < m_cutoff2) {
                    double sig_by_r6 = pow( sigma2 / rij2, 3 );
                    double sig_by_r12 = pow( sig_by_r6  , 2 );
                    energy += 4.0 * m_epsilon * (sig_by_r12 - sig_by_r6);
                }
            }
            return energy;
        }

double LJ::cutoff_correction(double volume,
                    int num_particles)
        {
            /*
            The function corrects interaction energy from energy cutoff.

            Parameters
            ----------

            volume : np.array
                Dimensions of the box in reduced units.

            num_particles : float
                Total number of particles in the box

            sigma : float
                Distance between two particles when interaction potential is zero

            epsilon: float
                Depth of potential well

            Return
            ------

            e_correction : float
                Correction energy from truncation
            */
            double sig_by_cutoff3 = pow( m_sigma / m_cutoff, 3);
            double sig_by_cutoff9 = pow( sig_by_cutoff3, 3);
            double e_correction = sig_by_cutoff9 - ( 3.0 * sig_by_cutoff3 );
            double n_particles_2 = pow(num_particles, 2);

            e_correction *= 8.0 / 9.0 * M_PI * n_particles_2 * m_epsilon * pow(m_sigma, 3) / volume;

            return e_correction;
        }


double LJ::operator()(std::vector<double> rij2_array)
        {
            return potential(rij2_array);
        }
