#include "box.hpp"
#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <math.h>
// C++ box class

double Box::volume(void)
 /*Calculate the box volume

  Returns
  -------
  box volume : float
      Computed box volume
  */
{
  double tmp = 1.0;
  for(auto it : box_dims)
    tmp *= it;
  return tmp;
}

Eigen::MatrixXd Box::coord_wrap(void)
/*Wraps coordinates within the box dimension

 Returns
 -------
 wrapped : matrix
     Arrays of the wrapped atomic coordinates.
 */    
{
  int nrows = coordinates.rows();
  int ncols = coordinates.cols();
  Eigen::MatrixXd wrapped(nrows, ncols);
  for(int i = 0; i < nrows; i++)
    {
	for(int j = 0; j < ncols; j++)
	  wrapped(i,j) = coordinates(i,j) - box_dims.at(j) * int (coordinates(i,j) / box_dims.at(j));
    }
  return wrapped;
}

std::vector<double> Box::minimum_image_dist(std::vector<double> i_particle)
/*Parameters
 ----------
 i_particle : array
     xyz coordinate of the i-th particle.

 Returns
 -------
 min_dist : array
     Array of the distances between each i-th particle and remaining
     particles
*/
{
  int nrows = coordinates.rows();
  int ncols = coordinates.cols();

  Eigen::MatrixXd coord_dist;
  for(int i = 0; i < nrows; i++)
    {
	for(int j = 0; j < ncols; j++)
	  coord_dist(i,j) = i_particle.at(j) - coordinates(i,j);
    }
  
  for(int i = 0; i < nrows; i++)
    {
	for(int j = 0; j < ncols; j++)
	  coord_dist(i,j) = coord_dist(i,j) - box_dims.at(j) * int (coord_dist(i,j) / box_dims.at(j));
    }

  std::vector<double> min_dist;
  for(int i = 0; i < nrows; i++)
    {
	double tmp = 0;
	for(int j = 0; j < ncols; j++)
	  tmp += coord_dist(i,j);
	min_dist.at(i) = std::sqrt(tmp);
    }
  return min_dist;
}
