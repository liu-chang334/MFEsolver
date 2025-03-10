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
    FiniteElementModel feModel;

    Eigen::SparseMatrix<double> K;
    Eigen::SparseVector<double> R;
    Eigen::VectorXd U;

    std::vector<C3D8> elements;
    std::vector<std::vector<Eigen::MatrixXd>> tmp_matrices;
    std::vector<Eigen::VectorXd> tensor_to_save;

public:
    FiniteElementSolver(FiniteElementModel feModel);

    void initializeStiffnessMatrix();
    void initializeResidual();
    void applyBoundaryConditions();
    void updateTangentMatrixAndResidual();
    void solve();
    void solveNonlinear();

    Eigen::VectorXd getElementNodesDisplacement(const int elementID);
    void calcuSpecifiedElementStress(const int elementID, Eigen::MatrixXd &strain, Eigen::MatrixXd &stress, bool extrapolatetoNodes = true); 
    void calcuSpecifiedElementStress(const int elementID, Eigen::MatrixXd &strain, Eigen::MatrixXd &stress,
                                    Eigen::MatrixXd &principalStrain, Eigen::MatrixXd &principalStress, 
                                    bool extrapolatetoNodes = true); 
    void getAllElementStress(bool is_principal = true);  
    void FiniteElementSolver::avgFieldAtNodes(int matrix_index, const std::string& filename, bool is_write = true);
    
};
#endif