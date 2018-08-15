#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // initially set to false, set to true in first call of ProcessMeasurement
  is_initialized_ = false;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // previous timestamp
  previous_timestamp_ = 0;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;

  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  // State dimension
  n_x_ = 5;

  // Augmented state dimension
  n_aug_ = n_x_ + 2;

  // number of sigma points
  n_sigma_ = 2 * n_aug_ + 1;

  // Sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  // Weights of sigma points
  weights_ = VectorXd(n_sigma_);
	for (int i = 0; i < n_sigma_; ++i) {
		if (i == 0) {
			weights_[i] = lambda_ / (lambda_ + n_aug_);
		} else {
			weights_[i] = 0.5 / (lambda_ + n_aug_);
		}
	}

  // predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_aug_, n_sigma_);
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} measurement_pack The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage measurement_pack) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
      * Initialize the state x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "UKF: " << endl;
    x_.setConstant(1);

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      x_ = h_radar_function_inverse(measurement_pack.raw_measurements_);
      // ToDo: may tune the covariance matrix initialization based on the sensor measurement type
      P_.setIdentity();
    } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      x_ = h_lidar_function_inverse(measurement_pack.raw_measurements_);
      // ToDo: may tune the covariance matrix initialization based on the sensor measurement type
      P_.setIdentity();
    }

    previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
   */

  //compute the time elapsed between the current and previous measurements
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;

  Prediction(dt);

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    UpdateRadar(measurement_pack);
  } else {
    // Laser updates
    UpdateLidar(measurement_pack);
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} dt the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double dt) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} measurement_pack
 */
void UKF::UpdateLidar(MeasurementPackage measurement_pack) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} measurement_pack
 */
void UKF::UpdateRadar(MeasurementPackage measurement_pack) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
}

VectorXd h_radar_function(VectorXd x) {
	double p_x = x(0);
	double p_y = x(1);
	double v = x(2);
	double yaw = x(3);
	double yawd = x(4);

  double rho = sqrt(p_x*p_x + p_y*p_y);
  double theta = atan2(p_y, p_x);
  double rho_dot = v * cos(yaw - theta);

  VectorXd result(3);
  result << rho, theta, rho_dot;
  return result;
}

VectorXd h_radar_function_inverse(VectorXd hx) {
  double rho = hx(0);
  double theta = hx(1);
  double rho_dot = hx(2);

  double p_x = rho * cos(theta);
  double p_y = rho * sin(theta);

  VectorXd result(5);
  result << p_x, p_y, 0, 0, 0;
  return result;
}

VectorXd h_lidar_function(VectorXd x) {
	double p_x = x(0);
	double p_y = x(1);
	double v = x(2);
	double yaw = x(3);
	double yawd = x(4);

  VectorXd result(2);
  result << p_x, p_y;
  return result;
}

VectorXd h_lidar_function_inverse(VectorXd hx) {
  double p_x = hx(0);
  double p_y = hx(1);

  VectorXd result(5);
  result << p_x, p_y, 0, 0, 0;
  return result;
}
