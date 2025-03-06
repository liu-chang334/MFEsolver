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
    std::vector<Eigen::MatrixXd> AllStrain;
    std::vector<Eigen::MatrixXd> AllStress;
    std::vector<Eigen::VectorXd> Strain_avg_at_node;
    std::vector<Eigen::VectorXd> Stress_avg_at_node;

public:
    FiniteElementSolver(FiniteElementModel feModel);
    void assembleStiffnessMatrix();
    void assembleForceVector();
    void applyBoundaryConditions();
    void solve_U();
    void solve();

    Eigen::VectorXd getElementNodesDisplacement(const int elementID);

    Eigen::MatrixXd calcuElementStrain(const int elementID, bool extrapolatetoNodes = true);
    Eigen::MatrixXd calcuElementStress(const int elementID, bool extrapolatetoNodes = true); 
    void calcuAllElementStrain(bool extrapolatetoNodes = true, bool is_write = false);
    void calcuAllElementStress(bool extrapolatetoNodes = true, bool is_write = false);  
    void avgStrainAtNodes(bool is_write = true);
    void avgStressAtNodes(bool is_write = true);
    
};
#endif