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
FiniteElementSolver::FiniteElementSolver(const FiniteElementModel& feModel) : feModel(feModel){
    const Eigen::MatrixXd& Node = feModel.Node;
    const Eigen::MatrixXi& Element = feModel.Element;
    const Eigen::MatrixXd& material = feModel.Material;
    
    int numNodes = static_cast<int>(Node.rows());
    int numElements = static_cast<int>(Element.rows());

    K = Eigen::SparseMatrix<double>(numNodes * 3, numNodes * 3);
    R = Eigen::SparseVector<double>(numNodes * 3, 1);
    U = Eigen::VectorXd::Zero(numNodes * 3);

    double E = material(0, 0);
    double nu = material(0, 1);
    mat.setproperties(E, nu);

    initializeResidual();
    
    elements.reserve(numElements);

    for (int i = 0; i < numElements; i++)
    {
        std::string elemType = "C3D8"; // FIXME: only C3D8 is supported for now
        Eigen::VectorXi elemnode = Element.row(i);
        Eigen::MatrixXd elemnodeCoor = Eigen::MatrixXd::Zero(8, 3); 
        for (int j = 0; j < 8; j++)
        {
            elemnodeCoor.row(j) = Node.row(elemnode(j) - 1);
        }

        C3D8 elem(i + 1, &mat);
        elem.setNodes(elemnode, elemnodeCoor);
        elem.initMaterialPoints();

        // elem.setMaterial(mat);
        // elem->setMaterialByYoungPoisson(E, nu);
        elements.push_back(elem);
    }   
}

/**
 * @brief Assemble the Global Stiffness Matrix
 * @note The element supported is only C3D8 for now, and the Global Stiffness Matrix 
 *      is not dealed with constraint yet
 * @note Now it takes a long time to assemble the Global Stiffness Matrix, need to be optimized
 */
void FiniteElementSolver::initializeStiffnessMatrix()
{
    const Eigen::MatrixXd& Node = feModel.Node;
    const Eigen::MatrixXi& Element = feModel.Element;
    int numNodes = static_cast<int>(Node.rows());
    int numElements = static_cast<int>(Element.rows());

    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(numElements * 24 * 24);
    
    for (int i = 0; i < numElements; i++)
    {
        C3D8& elem = elements[i];
        Eigen::VectorXd elemQ;
        Eigen::MatrixXd elemK;
        Eigen::VectorXd elemU = getElementNodesDisplacement(i + 1);
        // std::cout << "elemU: " << elemU.transpose() << std::endl;
        elem.calcuTangentAndResidual(elemU, elemQ, elemK);

        Eigen::VectorXi elemNodes = feModel.Element.row(i);
        Eigen::VectorXi elemnodedof(24);
        for (int j = 0; j < 8; j++)
        {
            int nodeIndex = elemNodes(j);
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
void FiniteElementSolver::initializeResidual()
{
    Eigen::MatrixXd Force = feModel.Force;
    Eigen::MatrixXd Node = feModel.Node;
    int numNodes = static_cast<int>(Node.rows());
    int numForce = static_cast<int>(Force.rows());

    // loop over the force
    for (int i = 0; i < numForce; i++)
    {
        int forcenode = static_cast<int>(Force(i, 0));
        int forcedirection = static_cast<int>(Force(i, 1));
        double forcevalue = static_cast<double>(Force(i, 2));

        R.coeffRef(3 * forcenode - 4 + forcedirection) += forcevalue;
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

        R.coeffRef(constrglobalindex) = constrvalue * K.coeff(constrglobalindex, constrglobalindex) * big_num;
        K.coeffRef(constrglobalindex, constrglobalindex) = K.coeff(constrglobalindex, constrglobalindex) * big_num;
    }
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
            elemU(i * 3 + j) = U(3 * nodeIndex - 3 + j);
        } 
    }
    return elemU;
}

void FiniteElementSolver::calcuSpecifiedElementStress(const int elementID, Eigen::MatrixXd &strain, Eigen::MatrixXd &stress, bool extrapolatetoNodes)
{
    C3D8& elem = elements[elementID - 1];
    strain = Eigen::MatrixXd::Zero(6, 8);
    stress = Eigen::MatrixXd::Zero(6, 8);
    for (int j = 0; j < 8; j++)
    {
        strain.col(j) = elem.materialPoints[j].strain;
        stress.col(j) = elem.materialPoints[j].stress;
    }
    if (extrapolatetoNodes)
    {
        strain = elem.extrapolateTensor(strain);
        stress = elem.extrapolateTensor(stress);
    }
}

void FiniteElementSolver::calcuSpecifiedElementStress(const int elementID, Eigen::MatrixXd &strain, Eigen::MatrixXd &stress,
                                    Eigen::MatrixXd &principalStrain, Eigen::MatrixXd &principalStress, bool extrapolatetoNodes)
{
    C3D8& elem = elements[elementID - 1];
    strain = Eigen::MatrixXd::Zero(6, 8);
    stress = Eigen::MatrixXd::Zero(6, 8);
    principalStrain = Eigen::MatrixXd::Zero(6, 8);
    principalStress = Eigen::MatrixXd::Zero(6, 8);
    for (int j = 0; j < 8; j++)
    {
        strain.col(j) = elem.materialPoints[j].strain;
        stress.col(j) = elem.materialPoints[j].stress;
    }
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



void FiniteElementSolver::getAllElementStress(bool is_principal)
{
    std::vector<Eigen::MatrixXd> tmp_strain;
    std::vector<Eigen::MatrixXd> tmp_stress;
    std::vector<Eigen::MatrixXd> tmp_principalstrain;
    std::vector<Eigen::MatrixXd> tmp_principalstress;
    int numElements = static_cast<int>(feModel.Element.rows());
    for (int i = 0; i < numElements; i++)
    {
        C3D8& elem = elements[i];
        Eigen::MatrixXd matPointsStrain(6, 8);
        Eigen::MatrixXd matPointsStress(6, 8);
        Eigen::MatrixXd principalStrain(6, 8);
        Eigen::MatrixXd principalStress(6, 8);

        if (elem.materialPoints.size() != 8)
        {
            std::cerr << "Error: material points size is not 8!" << std::endl;
        }
        for (int j = 0; j < 8; j++)
        {
            if (elem.materialPoints[j].strain.size() !=6){
                std::cerr << "Error: material points strain size is not 6!" << std::endl; 
            }
            auto mp = elem.materialPoints[j];
            // std::cout << "mp.strain: " << mp.strain.transpose() << std::endl;
            matPointsStrain.col(j) = mp.strain;
            matPointsStress.col(j) = mp.stress;
        }

        // std::string folder = "F:\\MFEsolver\\FEoutput";
        // saveMatrix2TXT(matPointsStrain, folder, "strain_at_Gp.txt", 10);
        // saveMatrix2TXT(matPointsStress, folder, "stress_at_Gp.txt", 10);

        tmp_strain.push_back(elem.extrapolateTensor(matPointsStrain));
        tmp_stress.push_back(elem.extrapolateTensor(matPointsStress));
        
            // Method-1: calculate principal strain and stress at gauss points and then extrapolate to nodes
            // Eigen::MatrixXd elemStrain, elemStress, principalStrain, principalStress;
            // calcuElementPrincipalStress(i + 1, elemStrain, elemStress, principalStrain, principalStress, extrapolatetoNodes);
            // tmp_matrix1.push_back(elemStrain);
            // tmp_matrix2.push_back(elemStress);
            // tmp_matrix3.push_back(principalStrain);
            // tmp_matrix4.push_back(principalStress);

            // Method-2: extrapolate strain and stress to nodes first, then calculate principal strain and stress at nodes
        if (is_principal){
            principalStrain = compute_tensor6x8_principal(elem.extrapolateTensor(matPointsStrain));
            principalStress = compute_tensor6x8_principal(elem.extrapolateTensor(matPointsStress));
        }
        tmp_principalstrain.push_back(principalStrain);
        tmp_principalstress.push_back(principalStress);
        
    }
    tmp_matrices.push_back(tmp_strain);
    tmp_matrices.push_back(tmp_stress);
    tmp_matrices.push_back(tmp_principalstrain);
    tmp_matrices.push_back(tmp_principalstress);
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

    tensor_to_save.clear();

    for (int i = 0; i < numNodes; i++)
    {
        if (count_at_node[i] > 0)
        {
            tensor_to_save.push_back(sumField_at_node[i] / double(count_at_node[i]));
        }
        else
        {
            tensor_to_save.push_back(Eigen::VectorXd::Zero(6));
            std::cerr << "Warning: node " << i + 1 << " has no elements!" << std::endl;
        }
    }

    if (is_write)
    {
        std::string current_path = std::filesystem::current_path().string();
        std::string resultpath = current_path + "\\.." + "\\.." + "\\FEoutput";
        if (std::filesystem::exists(resultpath))
        {
            saveMatrix2TXT(tensor_to_save, resultpath, filename);
        }
    }
}

void FiniteElementSolver::solve()
{
    {
        Spinner spinner(
            "\033[1;32m[Process]\033[0m Assembling Global Stiffness Matrix K: ", 
            "Global Stiffness Matrix K assembly time: ",
            100);
        initializeStiffnessMatrix();
    }
    std::cout << "\033[1;32m[Info]\033[0m "
            << "Global Stiffness Matrix K is "
            << "\033[1;32m" << "[" << K.rows() << ", " << K.cols() << "]\033[0m\n";
    // {
    //     Spinner spinner(
    //         "\033[1;32m[Process]\033[0m Assembling Global Force Vector F: ",
    //         "Global Force Vector F assembly time: ",
    //         1);    
    //     initializeResidual();
    // }
    std::cout << "\033[1;32m[Info]\033[0m "
            << "Global Force Vector F is "
            << "\033[1;32m" << "[" << R.rows() << ", " << R.cols() << "]\033[0m\n";
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
        U = solver.solve(R);
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

    {
        Spinner spinner(
            "\033[1;32m[Process]\033[0m Calculating strain and stress: ",
            "Strain and stress calculation time: ", 100);
        
        int numElements = static_cast<int>(feModel.Element.rows());
        for (int i = 0; i < numElements; i++)
        {
            C3D8& elem = elements[i];
            Eigen::VectorXd elemU = getElementNodesDisplacement(i + 1);
            Eigen::VectorXd elemQ;
            Eigen::MatrixXd elemK;
            elem.calcuTangentAndResidual(elemU, elemQ, elemK);
        }

        getAllElementStress(true);
        std::vector<std::string> filenames = {"Strain.txt", "Stress.txt", "PrincipalStrain.txt", "PrincipalStress.txt"};
        for (int i = 0; i < 4; i++)
        {
            avgFieldAtNodes(i, filenames[i], true);
        }
    }
}

void FiniteElementSolver::updateTangentMatrixAndResidual()
{
    const Eigen::MatrixXd& Node = feModel.Node;
    const Eigen::MatrixXi& Element = feModel.Element;
    int numNodes = static_cast<int>(Node.rows());
    int numElements = static_cast<int>(Element.rows());

    std::vector<Eigen::Triplet<double>> tripletList; 
    tripletList.reserve(numElements * 24 * 24);
    Eigen::SparseVector<double> Q = Eigen::SparseVector<double>(numNodes * 3, 1);

    for (int i = 0; i < numElements; i++)
    {
        C3D8& elem = elements[i];
        Eigen::VectorXd elemQ;
        Eigen::MatrixXd elemK;
        Eigen::VectorXd elemU = getElementNodesDisplacement(i + 1);
        elem.calcuTangentAndResidual(elemU, elemQ, elemK);

        Eigen::VectorXi elemNodes = feModel.Element.row(i);
        Eigen::VectorXi elemnodedof(24);

        for (int j = 0; j < 8; j++)
        {
            int nodeIndex = elemNodes(j);
            for (int k = 0; k < 3; k++)
            {
                elemnodedof(j * 3 + k) = nodeIndex * 3 - 3 + k;
            }
        }

        for (int j = 0; j < 24; j++)
        {
            int globalRow = elemnodedof(j);
            Q.coeffRef(globalRow) += elemQ(j);
            for (int k = 0; k < 24; k++)
            {
                int globalCol = elemnodedof(k);
                tripletList.emplace_back(globalRow, globalCol, elemK(j, k));
            }
        }
    }
    K.setFromTriplets(tripletList.begin(), tripletList.end());
    R -= Q;
}
void FiniteElementSolver::solveNonlinear()
{
    const int maxIter = 20;
    const double tol = 1e-6;

    for (int iter = 0; iter < maxIter; iter++)
    {
        updateTangentMatrixAndResidual();
        applyBoundaryConditions();
        double normR = R.norm();
        
        if (normR < tol)
        {
            std::cout << "Converged! Norm of residual: " << normR << std::endl;
            break;
        }

        std::cout << "Iteration: " << iter + 1 << std::endl;
        std::cout << "Norm of residual: " << normR << std::endl;

        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        Eigen::VectorXd du;
        solver.compute(K);
        du = solver.solve(R);
        U += du;
    }

    std::string current_path = std::filesystem::current_path().string();
    std::string resultpath = current_path + "\\.." + "\\.." + "\\FEoutput";
    if (std::filesystem::exists(resultpath)){
        std::filesystem::remove_all(resultpath);
    }
    std::filesystem::create_directory(resultpath);  
    saveMatrix2TXT(U, resultpath, "Displacement.txt", 12);

    getAllElementStress(true);
    std::vector<std::string> filenames = {"Strain.txt", "Stress.txt", "PrincipalStrain.txt", "PrincipalStress.txt"};
    for (int i = 0; i < 4; i++)
    {
        avgFieldAtNodes(i, filenames[i], true);
    }
}