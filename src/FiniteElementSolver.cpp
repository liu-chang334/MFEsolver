#include "../include/FiniteElementSolver.h"
#include "../include/C3D8.h"
#include "../include/Tools.h"


// constructor for FiniteElementSolver
FiniteElementSolver::FiniteElementSolver(FiniteElementModel feModel) : feModel(feModel){}

// assemble the Global Stiffness Matrix
void FiniteElementSolver::assembleStiffnessMatrix()
{
    Eigen::MatrixXd Node = feModel.Node;
    Eigen::MatrixXi Element = feModel.Element;
    Eigen::MatrixXd Material = feModel.Material;

    // get the size of the matrix
    int numNodes = static_cast<int>(Node.rows());
    int numElements = static_cast<int>(Element.rows());
    int numMaterials = static_cast<int>(Material.rows());
    
    // initialize the stiffness matrix
    K = Eigen::SparseMatrix<double>(numNodes * 3, numNodes * 3);

    // loop over the elements
    for (int i = 0; i < numElements; i++)
    {
        // get the element type
        std::string elemType = "C3D8"; // not implemented yet

        // get the element nodes coordinates matrix
        Eigen::VectorXi elemnode = Element.row(i);
        Eigen::MatrixXd elemnodeCoor = Eigen::MatrixXd::Zero(8, 3);
        for (int j = 0; j < 8; j++)
        {
            elemnodeCoor.row(j) = Node.row(elemnode(j) - 1);
        }

        // initialize the element with EID i
        C3D8 elem(i + 1);
        elem.setNodes(elemnode, elemnodeCoor);
        elem.setMaterialByYoungPoisson(Material(0, 0), Material(0, 1));
        // get the element stiffness matrix
        Eigen::MatrixXd elemK = elem.calcuStiffnessMatrix();

        // get the element dof in global coordinate system
        Eigen::VectorXi elemnodedof(24);
        for (int j = 0; j < 8; j++)
        {
            int nodeIndex = elemnode(j);
            // for each node
            for (int k = 0; k < 3; k++)
            {
                elemnodedof(j * 3 + k) = nodeIndex * 3 - 3 + k;
            }
        }
        // assemble the element stiffness matrix to the global stiffness matrix
        for (int j = 0; j < 24; j++)
        {
            for (int k = 0; k < 24; k++)
            {
                int globalRow = elemnodedof(j);
                int globalCol = elemnodedof(k);
                K.coeffRef(globalRow, globalCol) += elemK(j, k);
            }
        }
    }
}

// assemble the Global Force Vector
void FiniteElementSolver::assembleForceVector()
{
    Eigen::MatrixXd Force = feModel.Force;
    Eigen::MatrixXd Node = feModel.Node;
    int numNodes = static_cast<int>(Node.rows());
    int numForce = static_cast<int>(Force.rows());
    F = Eigen::SparseVector<double>(numNodes * 3, 1);

    // loop over the force
    for (int i = 0; i < numForce; i++)
    {
        int forcenode = static_cast<int>(Force(i, 0));
        int forcedirection = static_cast<int>(Force(i, 1));
        double forcevalue = static_cast<double>(Force(i, 2));

        F.coeffRef(3 * forcenode - 4 + forcedirection) += forcevalue;  // the F here is not dealed with constraint yet
    }
}

// apply the constraint --method: multiply a large number
void FiniteElementSolver::applyBoundaryConditions()
{
    long double big_num = 1e8;
    Eigen::MatrixXd Constr = feModel.Constraint;
    int numConstr = static_cast<int>(Constr.rows());
    for (int i = 0; i < numConstr; i++)
    {
        int constrnode = static_cast<int>(Constr(i, 0));
        int constrdirection = static_cast<int>(Constr(i, 1));
        int constrglobalindex = 3 * constrnode - 4 + constrdirection;
        double constrvalue = static_cast<double>(Constr(i, 2));

        F.coeffRef(constrglobalindex) = constrvalue * K.coeff(constrglobalindex, constrglobalindex) * big_num;
        K.coeffRef(constrglobalindex, constrglobalindex) = K.coeff(constrglobalindex, constrglobalindex) * big_num;
    }
}
// solve the linear system
void FiniteElementSolver::solve()
{
    // assemble the stiffness matrix
    assembleStiffnessMatrix();
    // assemble the force vector
    assembleForceVector();
    // apply the constraint method
    applyBoundaryConditions();
    
    // solve the linear system
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(K);
    U = solver.solve(F);
}

Eigen::VectorXd FiniteElementSolver::getElementNodesDisplacement(const int elementID)
{
    Eigen::VectorXi elementnodeID;
    feModel.getNodesIDofElement(elementID, elementnodeID);
    Eigen::VectorXd elemU(24);
    for (int i = 0; i < 8; i++) 
    {
        int nodeIndex = elementnodeID(i);
        for (int j = 0; j < 3; j++)
        {
            elemU(i * 3 + j) = U.coeff(3 * nodeIndex - 3 + j);
        } 
    }
    return elemU;
}

Eigen::MatrixXd FiniteElementSolver::calcuElementStrain(const int elementID, bool interpolatetoNodes)
{
    Eigen::VectorXd elemU = getElementNodesDisplacement(elementID);
    Eigen::MatrixXd elemStrain = Eigen::MatrixXd::Zero(6, 8);
    // get the element nodes coordinates matrix
    Eigen::VectorXi elemnode = feModel.Element.row(elementID - 1);
    Eigen::MatrixXd elemnodeCoor = Eigen::MatrixXd::Zero(8, 3);
    for (int j = 0; j < 8; j++)
    {
        elemnodeCoor.row(j) = feModel.Node.row(elemnode(j) - 1);
    }
    C3D8 elem(elementID);
    elem.setNodes(elemnode, elemnodeCoor);
    elemStrain = elem.calcuStrainTensor(elemU);
    if (!interpolatetoNodes)
    {
        return elemStrain; 
    }
    else
    {
        elemStrain = elem.interpolateTensor(elemStrain);
        return elemStrain;
    }
}

Eigen::MatrixXd FiniteElementSolver::calcuElementStress(const int elementID, bool interpolatetoNodes)
{
    Eigen::VectorXd elemU = getElementNodesDisplacement(elementID);
    Eigen::MatrixXd elemStress = Eigen::MatrixXd::Zero(6, 8);
    // get the element nodes coordinates matrix
    Eigen::VectorXi elemnode = feModel.Element.row(elementID - 1);
    Eigen::MatrixXd elemnodeCoor = Eigen::MatrixXd::Zero(8, 3);
    for (int j = 0; j < 8; j++)
    {
        elemnodeCoor.row(j) = feModel.Node.row(elemnode(j) - 1);
    }
    C3D8 elem(elementID);
    elem.setMaterialByYoungPoisson(feModel.Material(0, 0), feModel.Material(0, 1));
    elem.setNodes(elemnode, elemnodeCoor);
    elemStress = elem.calcuStressTensor(elemU);
    if (!interpolatetoNodes)
    {
        return elemStress;
    }
    else
    {
        elemStress = elem.interpolateTensor(elemStress);
        return elemStress;
    }
    
}