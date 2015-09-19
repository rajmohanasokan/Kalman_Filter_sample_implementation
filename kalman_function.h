#pragma once
#ifndef KALMAN_FUNCTION_H
#define KALMAN_FUNCTION_H
#include <iostream>
#include <Dense>
#include <vector>
#include <string>
#include <algorithm>
#include <iterator>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#include "engine.h"

class kalmanFilter{
public:
	Eigen::MatrixXd kalmanFunc(Eigen::MatrixXd phi, Eigen::MatrixXd upsilon, Eigen::MatrixXd basis, Eigen::MatrixXd initial, Eigen::MatrixXd initial_cov, int measurements, Eigen::MatrixXd noise);
};
#endif