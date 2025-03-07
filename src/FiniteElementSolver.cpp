#include "../include/FiniteElementSolver.h"
#include "../include/C3D8.h"
#include "../include/Tools.h"
#include "../include/Math.h"
#include <chrono>
#include "Spinner.cpp"

/**
 * @brief Solve the FEA problem
 * @note The element supported is only C3D8 for now, and the Global Stiffness Matrix
 *      is not dealed with constraint yet


/**
 * @brief Construct a new Finite Element Solver object
 *
 * @param[in] feModel Finite element model
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

    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(numElements * 24 * 24);
    
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
                // K.coeffRef(globalRow, globalCol) += elemK(j, k);
                tripletList.emplace_back(globalRow, globalCol, elemK(j, k));
            }
        }
    }
    K.setFromTriplets(tripletList.begin(), tripletList.end());
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
void FiniteElementSolver::solve_U()
{
    {
        Spinner spinner(
            "\033[1;32m[Process]\033[0m Assembling Global Stiffness Matrix K: ", 
            "Global Stiffness Matrix K assembly time: ",
            100);
        assembleStiffnessMatrix();
    }
    std::cout << "\033[1;32m[Info]\033[0m "
            << "Global Stiffness Matrix K is "
            << "\033[1;32m" << "[" << K.rows() << ", " << K.cols() << "]\033[0m\n";
    {
        Spinner spinner(
            "\033[1;32m[Process]\033[0m Assembling Global Force Vector F: ",
            "Global Force Vector F assembly time: ",
            1);    
        assembleForceVector();
    }
    std::cout << "\033[1;32m[Info]\033[0m "
            << "Global Force Vector F is "
            << "\033[1;32m" << "[" << F.rows() << ", " << F.cols() << "]\033[0m\n";
    {
        Spinner spinner(
            "\033[1;32m[Process]\033[0m Applying boundary conditions: ",
            "Boundary conditions application time: ",
            1 );
        applyBoundaryConditions();
    }
    std::cout << "\033[1;32m[Info]\033[0m "
            << "The method dealing with BC is \033[1;33m" << "multiply a large number\033[0m\n";
    {
        Spinner spinner(
            "\033[1;32m[Process]\033[0m Solving the system: ",
            "System solution time: ",
            100); 
        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        solver.compute(K);
        U = solver.solve(F);
    }
    std::cout << "\033[1;32m[Info]\033[0m "
            << "The method solving the system is \033[1;33m" << "LU decomposition\033[0m\n";

    std::string current_path = std::filesystem::current_path().string();
    std::string resultpath = current_path + "\\.." + "\\.." + "\\FEoutput";
    if (std::filesystem::exists(resultpath)){
        std::filesystem::remove_all(resultpath);
    }
    std::filesystem::create_directory(resultpath);  
    saveMatrix2TXT(U, resultpath, "Displacement.txt", 8);
}

/**
 * @brief Get the displacement at nodes of specified element from the 
 *      Global Displacement Vector
 * @param[in] elementID Element ID
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
 * @param[in] elementID Element ID
 * @param[in] interpolatetoNodes If true, the strain tensor is interpolated to nodes
 * @return Eigen::MatrixXd Strain tensor at gauss points defaultly at gauss points,
 *      if interpolatetoNodes is true, the strain tensor is interpolated to nodes
 */
Eigen::MatrixXd FiniteElementSolver::calcuElementStrain(const int elementID, bool extrapolatetoNodes)
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
    if (!extrapolatetoNodes)
    {
        return elemStrain; 
    }
    else
    {
        elemStrain = elem.extrapolateTensor(elemStrain);
        return elemStrain;
    }
}

/**
 * @brief Calculate the stress tensor at gauss points in specified element
 * @param[in] elementID Element ID
 * @param[in] interpolatetoNodes If true, the stress tensor is interpolated to nodes
 * @return Eigen::MatrixXd Stress tensor at gauss points defaultly at gauss points,
 *      if interpolatetoNodes is true, the stress tensor is interpolated to nodes
 */
std::pair<Eigen::MatrixXd, Eigen::MatrixXd> FiniteElementSolver::calcuElementStress(const int elementID, bool extrapolatetoNodes)
{
    Eigen::VectorXd elemU = getElementNodesDisplacement(elementID);
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
    // elemStress = elem.calcuStressTensor(elemU);
    auto[elemStrain, elemStress] = elem.calcuStrainStressTensor(elemU);
    if (!extrapolatetoNodes)
    {
        return std::make_pair(elemStrain, elemStress);
    }
    else
    {
        elemStrain = elem.extrapolateTensor(elemStrain);
        elemStress = elem.extrapolateTensor(elemStress);
        return std::make_pair(elemStrain, elemStress);
    }
}

void FiniteElementSolver::calcuElementPrincipalStress(const int elementID, Eigen::MatrixXd &strain, Eigen::MatrixXd &stress,
                                    Eigen::MatrixXd &principalStrain, Eigen::MatrixXd &principalStress, bool extrapolatetoNodes)
{
    Eigen::VectorXd elemU = getElementNodesDisplacement(elementID);
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
    // elemStress = elem.calcuStressTensor(elemU);
    auto[elemStrain, elemStress] = elem.calcuStrainStressTensor(elemU);
    strain = elemStrain;
    stress = elemStress;
    principalStrain = compute_tensor6x8_principal(strain);
    principalStress = compute_tensor6x8_principal(stress);
    if (extrapolatetoNodes)
    {
        strain = elem.extrapolateTensor(strain);
        stress = elem.extrapolateTensor(stress);
        principalStrain = elem.extrapolateTensor(principalStrain);
        principalStress = elem.extrapolateTensor(principalStress);
    }
}

/**
 * @brief Calculate the strain tensor at gauss points in all elements
 * @param[in] extrapolatetoNodes If true, the strain tensor is extrapolated to nodes
 */
void FiniteElementSolver::calcuAllElementStrain(bool extrapolatetoNodes)
{
    std::vector<Eigen::MatrixXd> tmp_matrix1;
    int numElements = static_cast<int>(feModel.Element.rows());
    for (int i = 0; i < numElements; i++)
    {
        Eigen::MatrixXd elemStrain = calcuElementStrain(i + 1, extrapolatetoNodes);
        tmp_matrix1.push_back(elemStrain);
    }
    tmp_matrices.push_back(tmp_matrix1);
}

void FiniteElementSolver::calcuAllElementStress(bool extrapolatetoNodes, bool is_principal)
{
    std::vector<Eigen::MatrixXd> tmp_matrix1;
    std::vector<Eigen::MatrixXd> tmp_matrix2;
    std::vector<Eigen::MatrixXd> tmp_matrix3;
    std::vector<Eigen::MatrixXd> tmp_matrix4;
    int numElements = static_cast<int>(feModel.Element.rows());
    for (int i = 0; i < numElements; i++)
    {
        if(!is_principal){
            auto[elemStrain, elemStress] = calcuElementStress(i + 1, extrapolatetoNodes);
            tmp_matrix1.push_back(elemStrain);
            tmp_matrix2.push_back(elemStress);
        }else{
            // Method-1: calculate principal strain and stress at gauss points and then extrapolate to nodes
            // Eigen::MatrixXd elemStrain, elemStress, principalStrain, principalStress;
            // calcuElementPrincipalStress(i + 1, elemStrain, elemStress, principalStrain, principalStress, extrapolatetoNodes);
            // tmp_matrix1.push_back(elemStrain);
            // tmp_matrix2.push_back(elemStress);
            // tmp_matrix3.push_back(principalStrain);
            // tmp_matrix4.push_back(principalStress);

            // Method-2: extrapolate strain and stress to nodes first, then calculate principal strain and stress at nodes
            auto[elemStrain, elemStress] = calcuElementStress(i + 1, true);
            Eigen::MatrixXd principalStrain, principalStress;
            principalStrain = compute_tensor6x8_principal(elemStrain);
            principalStress = compute_tensor6x8_principal(elemStress);
            tmp_matrix1.push_back(elemStrain);
            tmp_matrix2.push_back(elemStress);
            tmp_matrix3.push_back(principalStrain);
            tmp_matrix4.push_back(principalStress);
        }
    }
    tmp_matrices.push_back(tmp_matrix1);
    tmp_matrices.push_back(tmp_matrix2);
    tmp_matrices.push_back(tmp_matrix3);
    tmp_matrices.push_back(tmp_matrix4);
}

void FiniteElementSolver::avgFieldAtNodes(int matrix_index, const std::string& filename, bool is_write)
{
    int numNodes    = static_cast<int>(feModel.Node.rows());
    int numElements = static_cast<int>(feModel.Element.rows());

    if (matrix_index < 0 || matrix_index >= static_cast<int>(tmp_matrices.size()))
    {
        std::cerr << "Error: matrix_index " << matrix_index << " out of range!" << std::endl;
        return;
    }

    std::vector<Eigen::VectorXd> sumField_at_node(numNodes, Eigen::VectorXd::Zero(6));
    std::vector<int> count_at_node(numNodes, 0);
    
    for (int j = 0; j < numElements; j++)
    {
        Eigen::VectorXi elemNodes = feModel.Element.row(j);
        for (int k = 0; k < 8; k++)
        {
            int nodeId = elemNodes(k) - 1; // Transform to 0-based index
            sumField_at_node[nodeId] += tmp_matrices[matrix_index][j].col(k);
            count_at_node[nodeId]++;
        }
    }    

    Tensor.clear();

    for (int i = 0; i < numNodes; i++)
    {
        if (count_at_node[i] > 0)
        {
            Tensor.push_back(sumField_at_node[i] / double(count_at_node[i]));
        }
        else
        {
            Tensor.push_back(Eigen::VectorXd::Zero(6));
            std::cerr << "Warning: node " << i + 1 << " has no elements!" << std::endl;
        }
    }

    if (is_write)
    {
        std::string current_path = std::filesystem::current_path().string();
        std::string resultpath = current_path + "\\.." + "\\.." + "\\FEoutput";
        if (std::filesystem::exists(resultpath))
        {
            saveMatrix2TXT(Tensor, resultpath, filename);
        }
    }
}

void FiniteElementSolver::solve()
{
    solve_U();

    {
        Spinner spinner(
            "\033[1;32m[Process]\033[0m Calculating strain and stress: ",
            "Strain and stress calculation time: ", 100);

        calcuAllElementStress(true, true);
        std::vector<std::string> filenames = {"Strain.txt", "Stress.txt", "PrincipalStrain.txt", "PrincipalStress.txt"};
        for (int i = 0; i < 4; i++)
        {
            avgFieldAtNodes(i, filenames[i], true);
        }
    }
}