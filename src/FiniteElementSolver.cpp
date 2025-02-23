#include "../include/FiniteElementSolver.h"
#include "../include/C3D8.h"
#include "../include/Tools.h"
#include <chrono>

/**
 * @brief Solve the FEA problem
 * @note The element supported is only C3D8 for now, and the Global Stiffness Matrix
 *      is not dealed with constraint yet


/**
 * @brief Construct a new Finite Element Solver object
 *
 * @param feModel Finite element model
 */
FiniteElementSolver::FiniteElementSolver(FiniteElementModel feModel) : feModel(feModel){}

/**
 * @brief Assemble the Global Stiffness Matrix
 * @note The element supported is only C3D8 for now, and the Global Stiffness Matrix 
 *      is not dealed with constraint yet
 * @note Now it takes a long time to assemble the Global Stiffness Matrix, need to be optimized
 */
void FiniteElementSolver::assembleStiffnessMatrix()
{
    Eigen::MatrixXd Node = feModel.Node;
    Eigen::MatrixXi Element = feModel.Element;
    Eigen::MatrixXd Material = feModel.Material;

    int numNodes = static_cast<int>(Node.rows());
    int numElements = static_cast<int>(Element.rows());
    int numMaterials = static_cast<int>(Material.rows());
    
    K = Eigen::SparseMatrix<double>(numNodes * 3, numNodes * 3);
    for (int i = 0; i < numElements; i++)
    {
        std::string elemType = "C3D8"; // FIXME: only C3D8 is supported for now
        Eigen::VectorXi elemnode = Element.row(i);
        Eigen::MatrixXd elemnodeCoor = Eigen::MatrixXd::Zero(8, 3);
        
        for (int j = 0; j < 8; j++)
        {
            elemnodeCoor.row(j) = Node.row(elemnode(j) - 1);
        }

        C3D8 elem(i + 1);
        elem.setNodes(elemnode, elemnodeCoor);
        elem.setMaterialByYoungPoisson(Material(0, 0), Material(0, 1));
        Eigen::MatrixXd elemK = elem.calcuStiffnessMatrix();

        Eigen::VectorXi elemnodedof(24);
        for (int j = 0; j < 8; j++)
        {
            int nodeIndex = elemnode(j);
            for (int k = 0; k < 3; k++)
            {
                elemnodedof(j * 3 + k) = nodeIndex * 3 - 3 + k;
            }
        }

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

/**
 * @brief Assemble the Global Force Vector
 * @note The force is not dealed with constraint yet
 */
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

        // std::cout << "forcenode: " << forcenode << "\n forcedirection: " << forcedirection 
        //     << "\n forcevalue: " << forcevalue << std::endl;
        F.coeffRef(3 * forcenode - 4 + forcedirection) += forcevalue;
    }
}

// apply the constraint --method: multiply a large number
/**
 * @brief Apply the boundary conditions
 * @note The method is multiply a large number to the K matrix and F vector
 */
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

/**
 * @brief Assemble the Global Stiffness Matrix and the Global Force Vector,
 * and apply the boundary conditions, then solve the linear system
 * @note The method is solve the linear system by LU decomposition now
 */
void FiniteElementSolver::solve()
{
    auto t_start = std::chrono::steady_clock::now();
    assembleStiffnessMatrix();
    auto t_end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = t_end - t_start;
    std::cout 
        << "\033[1;33m[Time]\033[0m "  
        << "Global Stiffness Matrix K("
        << K.rows() << ", " << K.cols() << ") assembly time: "
        << "\033[1;32m" << elapsed_seconds.count() << "s\033[0m\n";
    
    t_start = std::chrono::steady_clock::now();
    assembleForceVector();
    t_end = std::chrono::steady_clock::now();
    elapsed_seconds = t_end - t_start;
    std::cout
        << "\033[1;33m[Time]\033[0m "  
        << "Global Force Vector F("
        << F.rows() << ", " << F.cols() << ") assembly time:"
        << "\033[1;32m" << elapsed_seconds.count() << "s\033[0m\n";

    t_start = std::chrono::steady_clock::now();         
    applyBoundaryConditions();
    t_end = std::chrono::steady_clock::now();
    elapsed_seconds = t_end - t_start;
    std::cout 
        << "\033[1;33m[Time]\033[0m "
        << "Boundary conditions application time:"
        << "\033[1;32m" << elapsed_seconds.count() << "s\033[0m\n";
    
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    t_start = std::chrono::steady_clock::now();
    solver.compute(K);
    U = solver.solve(F);
    t_end = std::chrono::steady_clock::now();
    elapsed_seconds = t_end - t_start;
    std::cout
        << "\033[1;33m[Time]\033[0m "
        << "The system solution time:"
        << "\033[1;32m" << elapsed_seconds.count() << "s\033[0m\n";

    std::string current_path = std::filesystem::current_path().string();
    std::string resultpath = current_path + "\\.." + "\\.." + "\\FEoutput";
    if (std::filesystem::exists(resultpath)){
        std::filesystem::remove_all(resultpath);
    }
    std::filesystem::create_directory(resultpath);  
    saveMatrix2TXT(U, resultpath, "U.txt", 8);
}

/**
 * @brief Get the displacement at nodes of specified element from the 
 *      Global Displacement Vector
 * @param elementID Element ID
 * @return Eigen::VectorXd Displacement vector at nodes of specified element
 */
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

/**
 * @brief Calculate the strain tensor at gauss points in specified element
 * @param elementID Element ID
 * @param interpolatetoNodes If true, the strain tensor is interpolated to nodes
 * @return Eigen::MatrixXd Strain tensor at gauss points defaultly at gauss points,
 *      if interpolatetoNodes is true, the strain tensor is interpolated to nodes
 */
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

/**
 * @brief Calculate the stress tensor at gauss points in specified element
 * @param elementID Element ID
 * @param interpolatetoNodes If true, the stress tensor is interpolated to nodes
 * @return Eigen::MatrixXd Stress tensor at gauss points defaultly at gauss points,
 *      if interpolatetoNodes is true, the stress tensor is interpolated to nodes
 */
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