#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <vector>
#include <math.h>


// C++ Potentials Class

class PairwisePotentials
{ 
    protected:
        double m_sigma;
        double m_epsilon;
        double m_cutoff;
        double m_cutoff2 = m_cutoff * m_cutoff;

    public:
        PairwisePotentials();

        virtual double potential(std::vector<double> rij2_array) = 0;
        virtual double operator()(std::vector<double> rij2_array) = 0;
};

class LJ : public PairwisePotentials
{   
    public:
        std::vector<double> m_rij2_array;
        
        LJ(double sigma = 1.0, double epsilon = 1.0, double cutoff = 2.6);

        double potential(std::vector<double> rij2_array);
        double cutoff_correction(double volume,
                    int num_particles);
        double operator()(std::vector<double> rij2_array);

};

