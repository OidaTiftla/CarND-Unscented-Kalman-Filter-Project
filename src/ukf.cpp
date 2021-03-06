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

  // if this is false, lidar measurements will be ignored (except during init)
  use_lidar_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // previous timestamp
  previous_timestamp_ = 0;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 4; // my review recommended values from 1 to 3

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = M_PI_2; // my review recommended values between 0 and 1

  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_lidar_px_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_lidar_py_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radar_rho_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radar_phi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radar_rhod_ = 0.3;
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
  Xsig_pred_ = MatrixXd(n_x_, n_sigma_);

  // NIS (normalized innovation squared) for lidar
  nis_lidar_ = std::vector<double>();

  // 95% NIS (normalized innovation squared) for lidar
  nis_lidar_95_ = 5.991;

  // NIS (normalized innovation squared) for radar
  nis_radar_ = std::vector<double>();

  // 95% NIS (normalized innovation squared) for radar
  nis_radar_95_ = 7.815;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} measurement_pack The latest measurement data of
 * either radar or lidar.
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
  const float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
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

  // print the output
  cout << "x_ = " << x_ << endl;
  cout << "P_ = " << P_ << endl;
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

  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sigma_);

  //create augmented mean state
  x_aug.setZero();
  x_aug.head(n_x_) = x_;

  //create augmented covariance matrix
  P_aug.setZero();
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  // ToDo: reviewer recommended to prepare this part of the matrix (the Q matrix) in the constructor
  P_aug(n_x_ + 0, n_x_ + 0) = std_a_ * std_a_;
  P_aug(n_x_ + 1, n_x_ + 1) = std_yawdd_ * std_yawdd_;

  //create square root matrix
  const MatrixXd P_aug_sqrt = P_aug.llt().matrixL();

  //create augmented sigma points

  // first point
  Xsig_aug.col(0) = x_aug;

  // other points
  const double factor = lambda_ + n_aug_;
  const MatrixXd sigma_points = sqrt(factor) * P_aug_sqrt;

  for (int i = 0; i < n_aug_; ++i) {
    Xsig_aug.col(i + 1) = x_aug + sigma_points.col(i);
    Xsig_aug.col(i + 1 + n_aug_) = x_aug - sigma_points.col(i);
  }

  //predict sigma points
  for (int i = 0; i < n_sigma_; ++i) {
    Xsig_pred_.col(i) = f_process_model_function(Xsig_aug.col(i).head(n_x_), Xsig_aug.col(i).bottomRows(n_aug_ - n_x_), dt);
  }

  //predict state mean
  x_.setZero();
  for (int i = 0; i < n_sigma_; ++i) {
    x_ += weights_[i] * Xsig_pred_.col(i);
  }

  //predict state covariance matrix
  P_.setZero();
  for (int i = 0; i < n_sigma_; ++i) {
    VectorXd Xi_minus_x = Xsig_pred_.col(i) - x_;

    //angle normalization because of difference calculation above
    while (Xi_minus_x(3) > M_PI) {
      Xi_minus_x(3) -= 2. * M_PI;
    }
    while (Xi_minus_x(3) < -M_PI) {
      Xi_minus_x(3) += 2. * M_PI;
    }

    P_ += weights_[i] * (Xi_minus_x * Xi_minus_x.transpose());
  }
}

/**
 * Updates the state and the state covariance matrix using a lidar measurement.
 * @param {MeasurementPackage} measurement_pack
 */
void UKF::UpdateLidar(MeasurementPackage measurement_pack) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

  // ToDo: The reviewer told me, for Lidar I can still use the linear kalman filter
  // equations. This saves computaional burden.

  //set measurement dimension, lidar can measure p_x and p_y
  const int n_z = 2;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, n_sigma_);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);

  //transform sigma points into measurement space
  for (int i = 0; i < n_sigma_; ++i) {
    Zsig.col(i) = h_lidar_function(Xsig_pred_.col(i)); // + 0: mean of measurement noise
  }

  //calculate mean predicted measurement
  z_pred.setZero();
  for (int i = 0; i < n_sigma_; ++i) {
    z_pred += weights_[i] * Zsig.col(i);
  }

  //calculate innovation covariance matrix S
  S.setZero();
  for (int i = 0; i < n_sigma_; ++i) {
    VectorXd Zi_minus_z_pred = Zsig.col(i) - z_pred;

    S += weights_[i] * (Zi_minus_z_pred * Zi_minus_z_pred.transpose());
  }
  //add measurement noise covariance matrix
  // ToDo: the reviewer recommended to prepare this matrix in the constructor, because it does not change
  MatrixXd R = MatrixXd(n_z, n_z);
  R.setZero();
  R(0, 0) = std_lidar_px_ * std_lidar_px_;
  R(1, 1) = std_lidar_py_ * std_lidar_py_;
  S += R;

  //sensor measurement
  const VectorXd z = measurement_pack.raw_measurements_;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  //calculate cross correlation matrix
  Tc.setZero();
  for (int i = 0; i < n_sigma_; ++i) {
    // state difference
    VectorXd Xi_minus_x = Xsig_pred_.col(i) - x_;
    //angle normalization because of difference calculation above
    while (Xi_minus_x(3) > M_PI) {
      Xi_minus_x(3) -= 2. * M_PI;
    }
    while (Xi_minus_x(3) < -M_PI) {
      Xi_minus_x(3) += 2. * M_PI;
    }

    //residual
    const VectorXd Zi_minus_z_pred = Zsig.col(i) - z_pred;

    Tc += weights_[i] * (Xi_minus_x * Zi_minus_z_pred.transpose());
  }

  //calculate Kalman gain K
  const MatrixXd Si = S.inverse();
  const MatrixXd K = Tc * Si;

  //residual
  const VectorXd z_minus_z_pred = z - z_pred;

  //update state mean and covariance matrix
  x_ += K * z_minus_z_pred;
  P_ -= K * S * K.transpose();

  //calculate the NIS (normalized innovation squared)
  const double epsilon = z_minus_z_pred.transpose() * Si * z_minus_z_pred;
  nis_lidar_.push_back(epsilon);
  const int n_below = std::count_if(nis_lidar_.begin(), nis_lidar_.end(), [this](double nis){ return nis < nis_lidar_95_; });
  const double percent_below = n_below * 100. / nis_lidar_.size();
  std::cout << "NIS lidar percent below: " << percent_below << std::endl;
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

  //set measurement dimension, radar can measure rho, phi, and rho_dot
  const int n_z = 3;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, n_sigma_);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);

  //transform sigma points into measurement space
  for (int i = 0; i < n_sigma_; ++i) {
    Zsig.col(i) = h_radar_function(Xsig_pred_.col(i)); // + 0: mean of measurement noise
  }

  //calculate mean predicted measurement
  z_pred.setZero();
  for (int i = 0; i < n_sigma_; ++i) {
    z_pred += weights_[i] * Zsig.col(i);
  }

  //calculate innovation covariance matrix S
  S.setZero();
  for (int i = 0; i < n_sigma_; ++i) {
    VectorXd Zi_minus_z_pred = Zsig.col(i) - z_pred;

    //angle normalization because of difference calculation above
    while (Zi_minus_z_pred(1) > M_PI) {
      Zi_minus_z_pred(1) -= 2. * M_PI;
    }
    while (Zi_minus_z_pred(1) < -M_PI) {
      Zi_minus_z_pred(1) += 2. * M_PI;
    }

    S += weights_[i] * (Zi_minus_z_pred * Zi_minus_z_pred.transpose());
  }
  //add measurement noise covariance matrix
  // ToDo: the reviewer recommended to prepare this matrix in the constructor, because it does not change
  MatrixXd R = MatrixXd(n_z, n_z);
  R.setZero();
  R(0, 0) = std_radar_rho_ * std_radar_rho_;
  R(1, 1) = std_radar_phi_ * std_radar_phi_;
  R(2, 2) = std_radar_rhod_ * std_radar_rhod_;
  S += R;

  //sensor measurement
  const VectorXd z = measurement_pack.raw_measurements_;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  //calculate cross correlation matrix
  Tc.setZero();
  for (int i = 0; i < n_sigma_; ++i) {
    // state difference
    VectorXd Xi_minus_x = Xsig_pred_.col(i) - x_;
    //angle normalization because of difference calculation above
    while (Xi_minus_x(3) > M_PI) {
      Xi_minus_x(3) -= 2. * M_PI;
    }
    while (Xi_minus_x(3) < -M_PI) {
      Xi_minus_x(3) += 2. * M_PI;
    }

    //residual
    VectorXd Zi_minus_z_pred = Zsig.col(i) - z_pred;
    //angle normalization because of difference calculation above
    while (Zi_minus_z_pred(1) > M_PI) {
      Zi_minus_z_pred(1) -= 2. * M_PI;
    }
    while (Zi_minus_z_pred(1) < -M_PI) {
      Zi_minus_z_pred(1) += 2. * M_PI;
    }

    Tc += weights_[i] * (Xi_minus_x * Zi_minus_z_pred.transpose());
  }

  //calculate Kalman gain K
  const MatrixXd Si = S.inverse();
  const MatrixXd K = Tc * Si;

  //residual
  VectorXd z_minus_z_pred = z - z_pred;
  //angle normalization because of difference calculation above
  while (z_minus_z_pred(1) > M_PI) {
    z_minus_z_pred(1) -= 2. * M_PI;
  }
  while (z_minus_z_pred(1) < -M_PI) {
    z_minus_z_pred(1) += 2. * M_PI;
  }

  //update state mean and covariance matrix
  x_ += K * z_minus_z_pred;
  P_ -= K * S * K.transpose();

  //calculate the NIS (normalized innovation squared)
  const double epsilon = z_minus_z_pred.transpose() * Si * z_minus_z_pred;
  nis_radar_.push_back(epsilon);
  const int n_below = std::count_if(nis_radar_.begin(), nis_radar_.end(), [this](double nis){ return nis < nis_radar_95_; });
  const double percent_below = n_below * 100. / nis_radar_.size();
  std::cout << "NIS radar percent below: " << percent_below << std::endl;
}

VectorXd h_radar_function(const VectorXd x) {
  const double p_x = x(0);
  const double p_y = x(1);
  const double v = x(2);
  const double yaw = x(3);
  // const double yawd = x(4);

  const double rho = sqrt(p_x*p_x + p_y*p_y);
  const double theta = atan2(p_y, p_x);
  const double rho_dot = v * cos(yaw - theta); // or: = (p_x * std::cos(yaw) * v + p_y * std::sin(yaw) * v) / sqrt(p_x*p_x + p_y*p_y);

  VectorXd result(3);
  result << rho, theta, rho_dot;
  return result;
}

VectorXd h_radar_function_inverse(const VectorXd hx) {
  const double rho = hx(0);
  const double theta = hx(1);
  // const double rho_dot = hx(2);

  const double p_x = rho * cos(theta);
  const double p_y = rho * sin(theta);

  VectorXd result(5);
  result << p_x, p_y, 0, 0, 0;
  return result;
}

VectorXd h_lidar_function(const VectorXd x) {
  const double p_x = x(0);
  const double p_y = x(1);
  // const double v = x(2);
  // const double yaw = x(3);
  // const double yawd = x(4);

  VectorXd result(2);
  result << p_x, p_y;
  return result;
}

VectorXd h_lidar_function_inverse(const VectorXd hx) {
  const double p_x = hx(0);
  const double p_y = hx(1);

  VectorXd result(5);
  result << p_x, p_y, 0, 0, 0;
  return result;
}

VectorXd f_process_model_function(VectorXd x, const VectorXd nu, const double dt) {
  // const double p_x = x[0];
  // const double p_y = x[1];
  const double v = x[2];
  const double yaw = x[3];
  const double yawd = x[4];

  const double nu_a = nu[0];
  const double nu_yawdd = nu[1];

  //avoid division by zero
  if (std::abs(yawd) < 0.0000001) {
    // yaw rate is zero
    x[0] += v * std::cos(yaw) * dt;
    x[1] += v * std::sin(yaw) * dt;
  } else {
    // yaw rate is not zero
    const double yaw_delta = yawd * dt;
    x[0] += v / yawd * (std::sin(yaw + yaw_delta) - std::sin(yaw));
    x[1] += v / yawd * (-std::cos(yaw + yaw_delta) + std::cos(yaw));
  }
  x[2] += 0;
  x[3] += yawd * dt;
  x[4] += 0;

  // add noise
  const double dt_squared_half = 0.5 * dt * dt;
  x[0] += dt_squared_half * std::cos(yaw) * nu_a;
  x[1] += dt_squared_half * std::sin(yaw) * nu_a;
  x[2] += dt * nu_a;
  x[3] += dt_squared_half * nu_yawdd;
  x[4] += dt * nu_yawdd;
  return x;
}
