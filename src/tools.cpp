#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;

  /* Check the validity of the following inputs
     * The estimation vector size should not be zero
     * The estimation vector size should be equal to the ground truth vector size
  */
  if(estimations.size()!=ground_truth.size() || estimations.size()==0){
    std::cout << "Invalid estimation or ground truth data" << std::endl;
    return rmse;
  }

  // Accumulate the squared error
  for(unsigned int i = 0; i < estimations.size(); ++i){
  	VectorXd residual = estimations[i] - ground_truth[i];
  	// Square the values
  	residual = residual.array() * residual.array();
  	rmse += residual;
  }

  // Calculate the mean of the squared error
  rmse = rmse/estimations.size();

  // Calculate the square root
  rmse = rmse.array().sqrt();

  //return the result
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */

  MatrixXd Hj(3, 4);
  // Get the state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  // Precompute some terms
  float c1 = px*px + py*py;
  float c2 = sqrt(c1);
  float c3 = c1*c2;

  // Check division by zero
  if (fabs(c1) < 0.0001){
    std::cout << "CalculateJacobian () - Error - Division by zero" << std::endl;
    return Hj;
  }
  // Compute the Jacobian
  Hj << (px/c2), (py/c2), 0, 0,
  			-(py/c1), (px/c1), 0, 0,
  			py*(vx*py -vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;

  return Hj;
}