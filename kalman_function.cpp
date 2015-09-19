#include "kalman_function.h"

Eigen::MatrixXd kalmanFilter::kalmanFunc(Eigen::MatrixXd phi, Eigen::MatrixXd upsilon, Eigen::MatrixXd basis, Eigen::MatrixXd initial, Eigen::MatrixXd initial_cov, int measurements, Eigen::MatrixXd noise){
		Eigen::MatrixXd x(measurements,2), xe(measurements,2), ym(measurements,1), covariance(measurements,2);
		Eigen::MatrixXd gain;
		x.setZero(measurements,2);
		x.row(0) = initial;
		
		
		xe.setZero(measurements,2);
		
		ym.setZero(measurements,1);
		ym.row(0) = (basis*x.row(0).transpose()).transpose()+(noise.cwiseSqrt().row(0))*Eigen::MatrixXd::Random(1,1);

		covariance.setZero(measurements,2);
		covariance.row(0) = initial_cov.diagonal().transpose();
		
		
		
		// Main loop
		
		for(int i=0; i<(measurements-1); i++){
			
			// Truth and Measurements
			x.row(i+1) = (phi*x.row(i).transpose()+upsilon*(noise.cwiseSqrt().row(1))*Eigen::MatrixXd::Random(1,1)).transpose();
			ym.row(i+1) = (basis*x.row(i+1).transpose()).transpose()+(noise.cwiseSqrt().row(0))*Eigen::MatrixXd::Random(1,1);
			

			
			// Update Equations
			
			gain = initial_cov*basis.transpose()*((basis*initial_cov*basis.transpose())+noise.row(0)).cwiseInverse();
			initial_cov = (Eigen::MatrixXd::Identity(2,2)-gain*basis)*initial_cov;
			xe.row(i) = xe.row(i)+(gain*(ym.row(i)-basis*xe.row(i).transpose())).transpose();


			
			// Propagation Equations
			
			xe.row(i+1) = (phi*xe.row(i).transpose()).transpose();
			initial_cov = phi*initial_cov*phi.transpose()+upsilon*noise.row(1)*upsilon.transpose();
			covariance.row(i+1) = initial_cov.diagonal().transpose();
			
		}

		
		Eigen::MatrixXd result(measurements,6);
		result << xe, x, (covariance.cwiseSqrt())*3;
		
		return result;







}