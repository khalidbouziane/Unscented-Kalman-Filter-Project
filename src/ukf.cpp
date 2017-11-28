#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = false;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.57;

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

  //
  is_initialized_ = false;

  // size of state vector
  n_x_ = 5;

  // size of augmented state
  n_aug_ = 7;

  // lambda component
  lambda_ = 3 - n_aug_;

  // Predicted signma matrix
  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_+1);

  // Weights matrix
  weights_ = VectorXd(2*n_aug_+1);
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage measurement_pack) {
  /**
  TODO:
  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  cout << "UKF:" << endl;

  // check if need to initialize
  if(!is_initialized_){

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // convert radar measurements to cartesian
      float rho = measurement_pack.raw_measurements_(0);
      float phi = measurement_pack.raw_measurements_(1);
      float rho_dot = measurement_pack.raw_measurements_(2);

      float px = rho*cos(phi);
      float py = rho*sin(phi);
      float vx = rho_dot*cos(phi);
      float vy = rho_dot*sin(phi);
      float v = sqrt(vx*vx + vy*vy);

      x_ << px, py, v, 0, 0;
      //cout << "Initial radar measurement." << endl;

    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
        // laser measurements already in cartesian
        float px = measurement_pack.raw_measurements_(0);
        float py = measurement_pack.raw_measurements_(1);
        cout << "px, py: (" << px << ", "<< py << ")" << endl;

        x_ << px, py, 0, 0, 0;
        //cout << "Initial laser measurement." << endl;
    }

    // init covariance matrix
    P_ << 1, 0, 0, 0, 0,
          0, 1, 0, 0, 0,
          0, 0, 1, 0, 0,
          0, 0, 0, 1, 0,
          0, 0, 0, 0, 1;

    // initialize weights_
    weights_(0) = lambda_ / (lambda_+n_aug_);
    for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
      weights_(i) = 0.5 / (n_aug_+lambda_);
    }

    // save timestamp
    time_us_ = measurement_pack.timestamp_;
    is_initialized_ = true;

    return;
  }

  // calculate time delta_t in seconds and save timestamp
  double delta_t = (measurement_pack.timestamp_ - time_us_) / 1000000.0;
  time_us_ = measurement_pack.timestamp_;

  //cout << "UKF::ProcessMeasurement -> "<< "delta_t = " << delta_t << endl;

  // Prediction
  Prediction(delta_t);
  //cout << "Prediction done." << endl;

  // Update based on the type of sensor
  if (measurement_pack.sensor_type_ == MeasurementPackage::LASER){
    //cout << "Updating lidar... " << endl;
    UpdateLidar(measurement_pack);
  }
  else if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR){
    //cout << "Updating radar... " << endl;
    UpdateRadar(measurement_pack);
  }

  // show state
  cout << "x_ = " << x_ << endl;
  cout << "P_ = " << P_ << endl;

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /*
  * Generate sigma points
  */

  // augmented state vector
  VectorXd x_aug = VectorXd(n_aug_);

  // create augmented state civaruace
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);


  // augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  // augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
  }

  //cout << "UKF::Prediction -> "<< "Augmented sigma points created" << endl;

  /*
  * Predict sigma points
  */

  //predict sigma points
  for (int i = 0; i< 2*n_aug_+1; i++)
  {
    //extract values for better readability
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin(yaw+yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    // add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    // save predicted sigma points
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }

  //cout << "UKF::Prediction -> "<< "Predicted sigma points" << endl;


  /****************************************************
  * Calculate mean and covariance
  ****************************************************/

  //predicted state mean
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }

  //cout << "UKF::Prediction -> "<< "Predicted state mean" << endl;

  //predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
  }

  //cout << "UKF::Prediction -> "<< "Predicted state covariance" << endl;
  //cout << "UKF::Prediction -> "<< "Exit" << endl;

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage measurement_pack) {
  // measurement size
  int n_z = 2;

  // convert prediction to measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2*n_aug_+1);
  for(int i = 0; i < 2*n_aug_ + 1; i++){
    Zsig(0, i) = Xsig_pred_(0, i);
    Zsig(1, i) = Xsig_pred_(1, i);
  }

  // calculate measurement mean
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for(int i=0; i < 2*n_aug_+1; i++){
    z_pred = z_pred + weights_(i)*Zsig.col(i);
  }

  //cout << "UKF::UpdateLidar -> " << "calculated measurement mean"<< endl;

  // caclulate measurement covariance
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);
  for(int i=0; i < 2*n_aug_+1; i++){
    VectorXd z_diff = Zsig.col(i) - z_pred;

    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //cout << "UKF::UpdateLidar -> " << "calculated measurement covariance"<< endl;

  // add noise matrix
  MatrixXd R = MatrixXd(n_z, n_z);
  R << std_laspx_*std_laspx_, 0,
      0, std_laspy_*std_laspy_;

  S = S + R;

  //cout << "UKF::UpdateLidar -> " << "added noise"<< endl;

  // calculate cross correlation matrix
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  for(int i =0; i < 2*n_aug_+1; i++){
    VectorXd z_diff = Zsig.col(i) - z_pred;

    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //cout << "UKF::UpdateLidar -> " << "calculated Tc"<< endl;

  // Get current measurement
  VectorXd z = measurement_pack.raw_measurements_;

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();
  NIS_laser_ = z.transpose() * S.inverse() * z;
  //cout << "UKF::UpdateLidar -> " << "Exit"<< endl;
}

  /**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage measurement_pack) {
  // measurement size
  int n_z = 3;

  // convert prediction to measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2*n_aug_+1);
  for(int i = 0; i < 2*n_aug_ + 1; i++){
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0, i) = sqrt(p_x*p_x + p_y*p_y);
    Zsig(1, i) = atan2(p_y, p_x);
    Zsig(2, i) = (p_x*v1 + p_y*v2) / sqrt(p_x*p_x + p_y*p_y);
  }

  // calculate measurement mean
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for(int i=0; i < 2*n_aug_+1; i++){
    z_pred = z_pred + weights_(i)*Zsig.col(i);
  }

  //cout << "UKF::UpdateRadar -> " << "calculated measurement mean"<< endl;

  // calculate measurement covariance
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);
  for(int i=0; i < 2*n_aug_+1; i++){
    VectorXd z_diff = Zsig.col(i) - z_pred;

    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //cout << "UKF::UpdateRadar -> " << "calculated measurement covariance"<< endl;


  // add noise
  MatrixXd R = MatrixXd(n_z, n_z);
  R << std_radr_*std_radr_, 0, 0,
      0, std_radphi_*std_radphi_, 0,
      0, 0, std_radrd_*std_radrd_;

  S = S + R;

  //cout << "UKF::UpdateRadar -> " << "added noise"<< endl;


  // calculate cross correlation
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  for(int i =0; i < 2*n_aug_+1; i++){
    VectorXd z_diff = Zsig.col(i) - z_pred;
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }
  //cout << "UKF::UpdateRadar -> " << "calculated Tc"<< endl;

  // get measurement
  VectorXd z = measurement_pack.raw_measurements_;

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();
  NIS_radar_ = z.transpose() * S.inverse() * z;
  //cout << "UKF::UpdateRadar -> " << "Exit"<< endl;
}

