#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:

  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  ///* if this is false, lidar measurements will be ignored (except for init)
  bool use_lidar_;

  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;

  ///* state covariance matrix
  MatrixXd P_;

  ///* predicted sigma points matrix
  MatrixXd Xsig_pred_;

  ///* previous timestamp
  long long previous_timestamp_;

  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  ///* Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  ///* Laser measurement noise standard deviation position1 in m
  double std_lidar_px_;

  ///* Laser measurement noise standard deviation position2 in m
  double std_lidar_py_;

  ///* Radar measurement noise standard deviation radius in m
  double std_radar_rho_;

  ///* Radar measurement noise standard deviation angle in rad
  double std_radar_phi_;

  ///* Radar measurement noise standard deviation radius change in m/s
  double std_radar_rhod_;

  ///* Weights of sigma points
  VectorXd weights_;

  ///* State dimension
  int n_x_;

  ///* Augmented state dimension
  int n_aug_;

  ///* number of sigma points
  int n_sigma_;

  ///* Sigma point spreading parameter
  double lambda_;

  ///* NIS (normalized innovation squared) for lidar
  std::vector<double> nis_lidar_;

  ///* 95% NIS (normalized innovation squared) for lidar
  double nis_lidar_95_;

  ///* NIS (normalized innovation squared) for radar
  std::vector<double> nis_radar_;

  ///* 95% NIS (normalized innovation squared) for radar
  double nis_radar_95_;


  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param measurement_pack The latest measurement data of either radar or lidar
   */
  void ProcessMeasurement(MeasurementPackage measurement_pack);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param dt Time between k and k+1 in s
   */
  void Prediction(double dt);

  /**
   * Updates the state and the state covariance matrix using a lidar measurement
   * @param measurement_pack The measurement at k+1
   */
  void UpdateLidar(MeasurementPackage measurement_pack);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param measurement_pack The measurement at k+1
   */
  void UpdateRadar(MeasurementPackage measurement_pack);
};

VectorXd h_radar_function(VectorXd x);
VectorXd h_radar_function_inverse(VectorXd hx);
VectorXd h_lidar_function(VectorXd x);
VectorXd h_lidar_function_inverse(VectorXd hx);
VectorXd f_process_model_function(VectorXd x, VectorXd nu, double dt);

#endif /* UKF_H */
