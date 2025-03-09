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
    std::vector<std::unique_ptr<C3D8>> elements;

    Eigen::SparseMatrix<double> K;
    Eigen::SparseVector<double> F;
    Eigen::SparseVector<double> U;
    std::vector<std::vector<Eigen::MatrixXd>> tmp_matrices;
    std::vector<Eigen::VectorXd> Tensor;

public:
    FiniteElementSolver(FiniteElementModel feModel);
    void assembleStiffnessMatrix();
    void assembleForceVector();
    void applyBoundaryConditions();
    void solve_U();
    void solve();

    Eigen::VectorXd getElementNodesDisplacement(const int elementID);

    Eigen::MatrixXd calcuElementStrain(const int elementID, bool extrapolatetoNodes = true);
    
    std::pair<Eigen::MatrixXd, Eigen::MatrixXd> calcuElementStress(const int elementID, bool extrapolatetoNodes = true); 
    void calcuElementPrincipalStress(const int elementID, Eigen::MatrixXd &strain, Eigen::MatrixXd &stress,
                                    Eigen::MatrixXd &principalStrain, Eigen::MatrixXd &principalStress, bool extrapolatetoNodes = true); 

    void calcuAllElementStrain(bool extrapolatetoNodes = true);
    void calcuAllElementStress(bool extrapolatetoNodes = true, bool is_principal = false);  

    void FiniteElementSolver::avgFieldAtNodes(int matrix_index, const std::string& filename, bool is_write = true);
    
};
#endif