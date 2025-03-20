#ifndef MATH_H
#define MATH_H

#include <iostream>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>

Eigen::VectorXd compute_eigenvalues(const Eigen::MatrixXd& matrix, bool start_from_largest);
Eigen::MatrixXd compute_tensor6x8_principal(const Eigen::MatrixXd& tensor);

std::vector<double> computeEigenvalues3x3(const std::vector<double> &vector);
std::vector<double> computeEigenvalues3x3(const Eigen::VectorXd &vector);

double computerTrace3x3(const std::vector<double> &vector);
double computerTrace3x3(const Eigen::VectorXd &vector);
double computerSecondInvariant3x3(const std::vector<double> &vector);
double computerSecondInvariant3x3(const Eigen::VectorXd &vector);
double computerThirdInvariant3x3(const std::vector<double> &vector);
double computerThirdInvariant3x3(const Eigen::VectorXd &vector);
double computerVolumetric3x3(const Eigen::VectorXd &vector);

Eigen::VectorXd computeDeviatoric(const Eigen::VectorXd &vector);

#endif