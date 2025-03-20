#include "../include/Math.h"

Eigen::VectorXd compute_eigenvalues(const Eigen::MatrixXd& matrix, bool start_from_largest) 
{
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(matrix);
    if (eigensolver.info() != Eigen::Success) {
        std::cerr << "Eigenvalue calculation failed!" << std::endl;
        return Eigen::VectorXd::Zero(matrix.rows());
    }
    if (start_from_largest) {
        return eigensolver.eigenvalues().reverse();
    }else {
        return eigensolver.eigenvalues();
     }
}

Eigen::MatrixXd compute_tensor6x8_principal(const Eigen::MatrixXd& tensor)
{
    Eigen::MatrixXd principal_tensor = Eigen::MatrixXd::Zero(6, 8);
    for(int i = 0; i < tensor.cols(); i++)
    {
        Eigen::VectorXd col = tensor.col(i);
        Eigen::Matrix3d temp;
        temp << col(0), col(3), col(4),
            col(3), col(1), col(5),
            col(4), col(5), col(2);

        Eigen::VectorXd eigenvalues = compute_eigenvalues(temp, true);
        
        principal_tensor.col(i).head(3) = eigenvalues;
    }    
    return principal_tensor;
}

/**
 * @brief Compute the eigenvalues of a 3x3 symmetric matrix from a 6-element vector.
 *
 * @param vector The input vector of size 6.
 * @return std::vector<double> The eigenvalues of the 3x3 symmetric matrix.
 * @note The input vector must be in the following order: [a11, a22, a33, a12, a13, a23].
 * @note The output vector will be from largest to smallest.
 */
std::vector<double> computeEigenvalues3x3(const std::vector<double>& vector)
{
   if (vector.size() != 6) {
       throw std::invalid_argument("Input vector must have size 6.");
   } 
   Eigen::Matrix3d matrix;
   matrix << vector[0], vector[3], vector[4],
             vector[3], vector[1], vector[5],
             vector[4], vector[5], vector[2];

   Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(matrix);
   if (eigensolver.info() != Eigen::Success) {
        throw std::runtime_error("Eigenvalue calculation failed!");
    }

   Eigen::VectorXd eigenvalues = eigensolver.eigenvalues();
   std::vector<double> result = {eigenvalues(2), eigenvalues(1), eigenvalues(0)};
   return result;
}

/** 
 * @brief Compute the eigenvalues of a 3x3 symmetric matrix from a 6-element vector.
 * 
 * @param vector The input vector of size 6.
 * @return std::vector<double> The eigenvalues of the 3x3 symmetric matrix.
 * @note The input vector must be in the following order: [a11, a22, a33, a12, a13, a23].
*/
std::vector<double> computeEigenvalues3x3(const Eigen::VectorXd &vector)
{
   if (vector.size() != 6) {
       throw std::invalid_argument("Input vector must have size 6.");
   } 
   Eigen::Matrix3d matrix;
   matrix << vector(0), vector(3), vector(4),
             vector(3), vector(1), vector(5),
             vector(4), vector(5), vector(2);

   Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(matrix);
   if (eigensolver.info() != Eigen::Success) {
        throw std::runtime_error("Eigenvalue calculation failed!");
    }

   Eigen::VectorXd eigenvalues = eigensolver.eigenvalues();
   std::vector<double> result = {eigenvalues(2), eigenvalues(1), eigenvalues(0)};
   return result;
}

/**
 * @brief Compute the trace of a 3x3 matrix from a 6-element vector.
 * 
 * @param vector The input vector of size 6.
 * @return double The trace of the 3x3 matrix.
 * @note The input vector must be in the following order: [a11, a22, a33, a12, a13, a23].
 */
double computerTrace3x3(const std::vector<double> &vector)
{
    if (vector.size() != 6) {
        throw std::invalid_argument("Input vector must have size 6."); 
    }
    double result = vector[0] + vector[1] + vector[2];
    return result;
}

/**
 * @brief Compute the trace of a 3x3 matrix from a 6-element vector.
 *
 * @param vector The input vector of size 6.
 * @return double The trace of the 3x3 matrix.
 * @note The input vector must be in the following order: [a11, a22, a33, a12, a13, a23].
 */
double computerTrace3x3(const Eigen::VectorXd &vector)
{
    if (vector.size() != 6) {
        throw std::invalid_argument("Input vector must have size 6."); 
    }
    double result = vector(0) + vector(1) + vector(2);
    return result;
}

/**
 * @brief Compute the second invariant of a 3x3 symmetric matrix from a 6-element vector.
 * 
 * @param vector The input vector of size 6.
 * @return double The second invariant of the 3x3 matrix.
 * @note The input vector must be in the following order: [a11, a22, a33, a12, a13, a23].
 */
double computerSecondInvariant3x3(const std::vector<double> &vector)
{
    if (vector.size()!= 6) {
        throw std::invalid_argument("Input vector must have size 6.");
    }
    double result = vector[0] * vector[1] + vector[1] * vector[2] + vector[2] * vector[0] - 
                    vector[3] * vector[3] - vector[4] * vector[4] - vector[5] * vector[5];
    return result;
}

/**
 * @brief Compute the second invariant of a 3x3 symmetric matrix from a 6-element vector.
 * 
 * @param vector The input vector of size 6.
 * @return double The second invariant of the 3x3 matrix.
 * @note The input vector must be in the following order: [a11, a22, a33, a12, a13, a23].
 */
double computerSecondInvariant3x3(const Eigen::VectorXd &vector)
{
    if (vector.size()!= 6) {
        throw std::invalid_argument("Input vector must have size 6.");
    }
    double result = vector(0) * vector(1) + vector(1) * vector(2) + vector(2) * vector(0) -
                    vector(3) * vector(3) - vector(4) * vector(4) - vector(5) * vector(5);
    return result;
}

/**
 * @brief Compute the third invariant of a 3x3 symmetric matrix from a 6-element vector.
 *
 * @param vector The input vector of size 6.
 * @return double The third invariant of the 3x3 matrix.
 * @note The input vector must be in the following order: [a11, a22, a33, a12, a13, a23].
 */
double computerThirdInvariant3x3(const std::vector<double> &vector)
{
    if (vector.size()!= 6) {
        throw std::invalid_argument("Input vector must have size 6.");
    }
    double result = vector[0] * vector[1] * vector[2] + 2 * vector[3] * vector[4] * vector[5] -
                    vector[0] * vector[5] * vector[5] - vector[1] * vector[4] * vector[4] -
                    vector[2] * vector[3] * vector[3];
    return result;
}

/**
 * @brief Compute the third invariant of a 3x3 symmetric matrix from a 6-element vector.
 *
 * @param vector The input vector of size 6.
 * @return double The third invariant of the 3x3 matrix.
 * @note The input vector must be in the following order: [a11, a22, a33, a12, a13, a23].
 */
double computerThirdInvariant3x3(const Eigen::VectorXd &vector)
{
    if (vector.size()!= 6) {
        throw std::invalid_argument("Input vector must have size 6.");
    }
    double result = vector(0) * vector(1) * vector(2) + 2 * vector(3) * vector(4) * vector(5) -
                    vector(0) * vector(5) * vector(5) - vector(1) * vector(4) * vector(4) -
                    vector(2) * vector(3) * vector(3);
    return result;
}

double computerVolumetric3x3(const Eigen::VectorXd &vector)
{
    if (vector.size()!= 6) {
        throw std::invalid_argument("Input vector must have size 6.");
    }
    double result = (vector(0) + vector(1) + vector(2)) / 3.0;
    return result; 
}

/**
 * @brief Compute the deviatoric part of a 6-element vector.
 *
 * @param vector The input vector of size 6.
 * @return Eigen::VectorXd The deviatoric part of the input vector.
 * @note The input vector must be in the following order: [a11, a22, a33, a12, a13, a23].
 */
Eigen::VectorXd computeDeviatoric(const Eigen::VectorXd &vector)
{
    if (vector.size()!= 6) {
        throw std::invalid_argument("Input vector must have size 6.");
    }

    Eigen::VectorXd result = vector;
    double trace = vector(0) + vector(1) + vector(2);
    result(0) -= trace / 3.0;
    result(1) -= trace / 3.0;
    result(2) -= trace / 3.0;

    return result;
}