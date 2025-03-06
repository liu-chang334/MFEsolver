#ifndef FINITE_ELEMENT_SOLVER_H
#define FINITE_ELEMENT_SOLVER_H

#include "FiniteElementModel.h"
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
    Eigen::SparseVector<double> F;
    Eigen::SparseVector<double> U;
    std::vector<Eigen::MatrixXd> tmp_matrix1;
    std::vector<Eigen::MatrixXd> tmp_matrix2;
    std::vector<Eigen::MatrixXd> tmp_matrix3;
    std::vector<Eigen::MatrixXd> tmp_matrix4;
    // std::vector<Eigen::MatrixXd> tmp_matrix5;
    // std::vector<Eigen::MatrixXd> tmp_matrix6;
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

    void avgStrainAtNodes(bool is_write = true);
    void avgStressAtNodes(bool is_write = true);
    void avgPrincipalStrainAtNodes(bool is_write = true);
    void avgPrincipalStressAtNodes(bool is_write = true);
    
};
#endif