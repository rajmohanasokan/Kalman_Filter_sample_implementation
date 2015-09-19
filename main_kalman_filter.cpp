#include "kalman_function.h"


int main(){
	Eigen::MatrixXd phi(2,2), upsilon(2,1), basis(1,2), initial(1,2), initial_cov(2,2), noise(2,1);
	
	int measurements;
	measurements = 1000; /* Total number of measurement values to be provided from the dynamic model */

	phi<<0.99985, 0.0098510, -0.029553, 0.97030; /* State/System Matrix of the Dynamic Model*/
	std::cout<<"Phi: \n"<<phi<<std::endl;

	upsilon<<4.9502e-5, 9.8510e-3; /*Input Matrix of the Model*/
	std::cout<<"Upsilon: "<<upsilon<<std::endl;

	basis<<1,0; /*Ouput Matrix of the Model*/
	std::cout<<"Basis: "<<basis<<std::endl;

	initial<<1,1; /*Initial state values*/
	initial_cov<<0.444,0,0,0.444; /*Initial Covariance of the state values*/
	std::cout<<"Initial covariance: \n"<<initial_cov<<std::endl;
	noise<<0.01,1; /*Process and Measurement Noise, i.e., noise(0,0) = Process noise; noise(0,1) = Measurement noise*/

	Eigen::MatrixXd output, result(measurements,4);
	kalmanFilter zz;
	output = zz.kalmanFunc(phi, upsilon, basis, initial, initial_cov, measurements, noise);
	
	std::cout<< output.bottomRows<10>()<< std::endl; /*Sample Display of the returned values*/

	result <<  output.col(0)-output.col(2), output.col(1)-output.col(3), output.col(4), output.col(5); /*Determining the state errors between synthetic data and estimated data*/

	std::cout<< result.bottomRows<10>()<< std::endl; /*Sample Display of the returned values*/
	
	Eigen::VectorXd t;
	t.setLinSpaced(1000,0,10);
	
	double xx[1000];
	Eigen::Map<Eigen::MatrixXd>(xx,1000,1) = result.col(0);
	
	double yy[1000];
	Eigen::Map<Eigen::MatrixXd>(yy, 1000,1) = result.col(1);
	
	double tt[1000];
	Eigen::Map<Eigen::MatrixXd>(tt, 1000,1) = t;
	
	double ss[1000];
	Eigen::Map<Eigen::MatrixXd>(ss, 1000,1) = result.col(2);
	

/*
	* Initiating the Matlab engine, and variables
	*/
	Engine *ep;
	mxArray *xax = NULL, *yax=NULL, *time = NULL, *sigma = NULL;

	/*
	 * Call engOpen with a NULL string. This starts a MATLAB process 
     * on the current host using the command "matlab".
	 */
	if (!(ep = engOpen(""))) {
		fprintf(stderr, "\nCan't start MATLAB engine\n");
		return EXIT_FAILURE;
	}

	/*
	 * Plotting the Ellipse in a Matlab plot window
	*/

	/* 
	 * Create a variable for the data
	*/
	xax = mxCreateDoubleMatrix(1000, 1, mxREAL);
	yax = mxCreateDoubleMatrix(1000, 1, mxREAL);
	time = mxCreateDoubleMatrix(1000, 1, mxREAL);
	sigma = mxCreateDoubleMatrix(1000, 1, mxREAL);

	memcpy((void *)mxGetPr(xax), (void *)xx, sizeof(xx));
	memcpy((void *)mxGetPr(yax), (void *)yy, sizeof(yy));
	memcpy((void *)mxGetPr(time), (void *)tt, sizeof(tt));
	memcpy((void *)mxGetPr(sigma), (void *)ss, sizeof(ss));
	
	/*
	 * Place the variables into the MATLAB workspace
	*/
	engPutVariable(ep, "xax", xax);
	engPutVariable(ep, "yax", yax);
	engPutVariable(ep, "time", time);
	engPutVariable(ep, "sigma", sigma);

	/*
	 * Evaluate a MATLAB Stringstring
	 * Plot the result using the MATLAB Plot functions
	*/
	engEvalString(ep, "plot ( time , sigma ,'r' , time , xax , 'b', time ,-sigma , 'r' )");
	engEvalString(ep, "axis([ 0 10 -0.3 0.3])");
	engEvalString(ep,"set(gca,'ytick',[-0.3 -0.2 -0.1 0 0.1 0.2 0.3])");
	engEvalString(ep, "title('Kalman Filter: State Estimates');");
	engEvalString(ep, "xlabel('time(sec)');");
	engEvalString(ep, "ylabel('x_1 error');");



	/*
	 * We're done! Free memory, close MATLAB figure after performing required operations.
	*/
	mxDestroyArray(xax);
	mxDestroyArray(yax);
	mxDestroyArray(time);
	mxDestroyArray(sigma);
	std::system("pause");
	engEvalString(ep, "close;");

	/*
	 * We're done! Free memory, close MATLAB engine and exit.
	 */
	printf("Done!\n");
	engClose(ep);
	
	return EXIT_SUCCESS;
}