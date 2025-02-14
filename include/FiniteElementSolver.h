#ifndef FINITE_ELEMENT_SOLVER_H
#define FINITE_ELEMENT_SOLVER_H

#include "FiniteElementModel.h"
#include <Eigen/Sparse>

// FEM Solver
class FiniteElementSolver
{
public:
    FiniteElementModel feModel;
    Eigen::SparseMatrix<double> K;
    Eigen::SparseVector<double> F;
    Eigen::SparseVector<double> U;
    std::vector<Eigen::MatrixXd> AllStrain;
    std::vector<Eigen::MatrixXd> AllStress;

public:
    FiniteElementSolver(FiniteElementModel feModel);
    void assembleStiffnessMatrix();
    void assembleForceVector();
    void applyBoundaryConditions();
    void solve();

    // get the displacement at nodes of element
    Eigen::VectorXd getElementNodesDisplacement(const int elementID);
    // calculate strain and stress tensor at gauss points in element
    Eigen::MatrixXd calcuElementStrain(const int elementID, bool interpolatetoNodes = true);
    Eigen::MatrixXd calcuElementStress(const int elementID, bool interpolatetoNodes = true);   
    
};
#endif