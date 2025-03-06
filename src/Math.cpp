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