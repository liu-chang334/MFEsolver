#ifndef FINITE_ELEMENT_SOLVER_H
#define FINITE_ELEMENT_SOLVER_H

#include "FiniteElementModel.h"
#include "C3D8.h"
#include <Eigen/Sparse>

/**
 * @class FiniteElementSolver
 * @brief This class represents a finite element solver.
 * - \c feModel: Finite element model --> FiniteElementModel Class
 * - \c K: Stiffness matrix for the system
 * - \c F: Force vector for the system
 * - \c U: Displacement vector for the system
 * - \c AllStrain: Strain tensor at gauss points in all elements
 * - \c AllStress: Stress tensor at gauss points in all elements
 */
class FiniteElementSolver
{
public:
    const FiniteElementModel& feModel;
    std::shared_ptr<Material> material_;

    Eigen::SparseVector<double> F_;
    Eigen::SparseMatrix<double> K_;
    Eigen::SparseVector<double> Q_;
    Eigen::SparseVector<double> R_;
    Eigen::VectorXd U_;
    Eigen::VectorXd DU_;

    std::vector<C3D8> elements;
    std::vector<std::vector<Eigen::MatrixXd>> tmp_matrices;
    std::vector<Eigen::VectorXd> tensor_to_save;

public:
    FiniteElementSolver(const FiniteElementModel& feModel);

    void initializematerial();
    void initializeExternalForce(Eigen::SparseVector<double>& F);
    void applyBoundaryConditions(Eigen::SparseMatrix<double>& K, Eigen::SparseVector<double>& R);
    void updateTangentMatrixAndInternal();
    void solve_linearelastic();
    void solve_adaptive_nonlinear(double& step_size, int& maxIter);
    bool perform_Newton_Raphson(int maxIter, double tol, double scaleFactor = 1.0, double step_size = 0.1);

    Eigen::VectorXd getElementNodes_DU(const int elementID);
    void getSpecifiedElementStress(const int elementID, Eigen::MatrixXd &strain, Eigen::MatrixXd &stress, 
                                    bool extrapolatetoNodes = true); 
    void getSpecifiedElementStress(const int elementID, Eigen::MatrixXd &strain, Eigen::MatrixXd &stress,
                                    Eigen::MatrixXd &principalStrain, Eigen::MatrixXd &principalStress, 
                                    bool extrapolatetoNodes = true); 
    void getAllElementStress(bool is_principal = true);  
    void FiniteElementSolver::avgFieldAtNodes(int matrix_index, const std::string& filename, 
                                                bool is_write = true);

    void FiniteElementSolver::writeToFile();
    
};
#endif