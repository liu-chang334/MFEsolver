#ifndef MATH_H
#define MATH_H

#include <iostream>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>

Eigen::VectorXd compute_eigenvalues(const Eigen::MatrixXd& matrix, bool start_from_largest);
Eigen::MatrixXd compute_tensor6x8_principal(const Eigen::MatrixXd& tensor);


#endif