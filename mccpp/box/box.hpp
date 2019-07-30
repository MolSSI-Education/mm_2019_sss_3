
#pragma once

#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <math.h>

class Box
{
  public:
     std::vector<double> box_dims;
   Eigen::MatrixXd coordinates;

  double volume(void);

  Eigen::MatrixXd coord_wrap(void);

  std::vector<double> minimum_image_dist(std::vector<double> i_particle);

};
