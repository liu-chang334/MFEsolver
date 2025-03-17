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
    const Eigen::MatrixXd& materialData = feModel.Material;
    
    int numNodes = static_cast<int>(Node.rows());
    int numElements = static_cast<int>(Element.rows());

    K_ = Eigen::SparseMatrix<double>(numNodes * 3, numNodes * 3);
    F_ = Eigen::SparseVector<double>(numNodes * 3, 1);
    Q_ = Eigen::SparseVector<double>(numNodes * 3, 1);
    R_ = Eigen::SparseVector<double>(numNodes * 3, 1);
    U_ = Eigen::VectorXd::Zero(numNodes * 3);
    DU_ = Eigen::VectorXd::Zero(numNodes * 3);
    initializeExternalForce(F_);

    initializematerial();
    std::unordered_map<std::string, double> matParams = {
        {"E", materialData(0, 0)},
        {"nu", materialData(0, 1)},
        {"sigma_y", 205.00},
        {"H", 0.00}
    };

    material_->setMaterialParameters(matParams);
    
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

        C3D8 elem(i + 1, material_);
        elem.setNodes(elemnode, elemnodeCoor);
        elem.initMaterialPoints();

        // elem.setMaterial(mat);
        // elem->setMaterialByYoungPoisson(E, nu);
        elements.push_back(elem);
    }   
}

/**
 * @brief Initialize the material
 * @note This function is called automatically when one FiniteElementSolver object is created
 */

void FiniteElementSolver::initializematerial()
{
    const std::string& matType = feModel.MaterialType;

    if (matType == "LinearElastic") {material_ = std::make_shared<LinearElastic>();}
    else if (matType == "IdealElastoplastic") {material_ = std::make_shared<IdealElastoplastic>();}
    else {std::cerr << "Error: Material type not supported" << std::endl; exit(1);}
}


/**
 * @brief Assemble the Residual Vector
 * @note This function is called automatically when one FiniteElementSolver object is created
 */
void FiniteElementSolver::initializeExternalForce(Eigen::SparseVector<double>& F)
{
    const Eigen::MatrixXd& Force = feModel.Force;
    const Eigen::MatrixXd& Node = feModel.Node;
    int numNodes = static_cast<int>(Node.rows());
    int numForce = static_cast<int>(Force.rows());

    // loop over the force
    for (int i = 0; i < numForce; i++)
    {
        int forcenode = static_cast<int>(Force(i, 0));
        int forcedirection = static_cast<int>(Force(i, 1));
        double forcevalue = static_cast<double>(Force(i, 2));

        F.coeffRef(3 * forcenode - 4 + forcedirection) += forcevalue;
    }
}

// apply the constraint --method: multiply a large number
/**
 * @brief Apply the boundary conditions
 * @note The method is multiply a large number to the K matrix and F vector
 */
void FiniteElementSolver::applyBoundaryConditions(Eigen::SparseMatrix<double>& K, Eigen::SparseVector<double>& R)
{
    long double big_num = 1e8;
    const Eigen::MatrixXd& Constr = feModel.Constraint;
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
Eigen::VectorXd FiniteElementSolver::getElementNodes_DU(const int elementID)
{
    Eigen::VectorXi elementnodeID;
    feModel.getNodesIDofElement(elementID, elementnodeID);
    Eigen::VectorXd delta_elemU(24);
    for (int i = 0; i < 8; i++) 
    {
        int nodeIndex = elementnodeID(i);
        for (int j = 0; j < 3; j++)
        {
            delta_elemU(i * 3 + j) = DU_(3 * nodeIndex - 3 + j);
        } 
    }
    return delta_elemU;
}

void FiniteElementSolver::getSpecifiedElementStress(const int elementID, Eigen::MatrixXd &strain, Eigen::MatrixXd &stress, bool extrapolatetoNodes)
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

void FiniteElementSolver::getSpecifiedElementStress(const int elementID, Eigen::MatrixXd &strain, Eigen::MatrixXd &stress,
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

        for (int j = 0; j < 8; j++)
        {
            if (elem.materialPoints[j].strain.size() !=6){
                std::cerr << "Error: material points strain size is not 6!" << std::endl; 
            }
            auto mp = elem.materialPoints[j];
            matPointsStrain.col(j) = mp.strain;
            matPointsStress.col(j) = mp.stress;
        }

        tmp_strain.push_back(elem.extrapolateTensor(matPointsStrain));
        tmp_stress.push_back(elem.extrapolateTensor(matPointsStress));
        
            // Method-1: calculate principal strain and stress at gauss points and then extrapolate to nodes

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


void FiniteElementSolver::updateTangentMatrixAndInternal()
{
    K_.setZero();
    Q_.setZero();

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
        Eigen::VectorXd delta_elemU = getElementNodes_DU(i + 1);
        elem.updateTangentAndInternal(delta_elemU, elemQ, elemK);

        const Eigen::VectorXi& elemNodes = feModel.Element.row(i);
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
            Q_.coeffRef(globalRow) += elemQ(j);
            for (int k = 0; k < 24; k++)
            {
                int globalCol = elemnodedof(k);
                tripletList.emplace_back(globalRow, globalCol, elemK(j, k));
            }
        }
    }
    K_.setFromTriplets(tripletList.begin(), tripletList.end());
}



void FiniteElementSolver::solve_linearelastic()
{
    {
        Spinner spinner(
            "\033[1;32m[Process]\033[0m Assembling Global Stiffness Matrix K: ", 
            "Global Stiffness Matrix K assembly time: ",
            100);
        updateTangentMatrixAndInternal();
    }
    std::cout << "\033[1;32m[Info]\033[0m "
            << "Global Stiffness Matrix K is "
            << "\033[1;32m" << "[" << K_.rows() << ", " << K_.cols() << "]\033[0m\n";

    {
        // Spinner spinner(
        //     "\033[1;32m[Process]\033[0m Applying boundary conditions: ",
        //     "Boundary conditions application time: ",
        //     1 );
        R_ = F_ - Q_;
        applyBoundaryConditions(K_, R_);
    }

    {
        Spinner spinner(
            "\033[1;32m[Process]\033[0m Solving the system: ",
            "System solution time: ",
            100); 
        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        solver.compute(K_);
        U_ = solver.solve(R_);
        DU_ = U_;
    }

    {  
        int numElements = static_cast<int>(feModel.Element.rows());
        for (int i = 0; i < numElements; i++)
        {
            C3D8& elem = elements[i];
            Eigen::VectorXd elemU = getElementNodes_DU(i + 1);
            Eigen::VectorXd elemQ;
            Eigen::MatrixXd elemK;
            elem.updateTangentAndInternal(elemU, elemQ, elemK);
        }
    }
    writeToFile();
}

bool FiniteElementSolver::perform_Newton_Raphson(int maxIter, double tol, double step_end)
{
    Eigen::VectorXd U_backup = U_;
    DU_ = Eigen::VectorXd::Zero(U_.size());

    // Copy the current strain and stress tensor to a temporary vector
    // We will use this vector to recover the strain and stress tensor when the solution is not converged
    std::vector<Eigen::MatrixXd> tmp_elem_strain;
    std::vector<Eigen::MatrixXd> tmp_elem_stress;
    for (auto& elem : elements)
    {
        int i = 0;
        Eigen::MatrixXd tmp_strain(6, 8);
        Eigen::MatrixXd tmp_stress(6, 8);
        for (auto& mp : elem.materialPoints)
        {
            tmp_strain.col(i) = mp.strain;
            tmp_stress.col(i) = mp.stress;
            i++;
        }
        tmp_elem_strain.push_back(tmp_strain);
        tmp_elem_stress.push_back(tmp_stress);
    }
    // To solve the nonlinear problem, we use the Newton-Raphson method
    // The iteration is performed until the residual norm is smaller than the tolerance
    // or the maximum number of iterations is reached
    for (int iter = 0; iter < maxIter; iter++)
    {
        updateTangentMatrixAndInternal(); 
        // Eigen::SparseVector<double> R_scale = (F_ * step_size - Q_) * scaleFactor;
        Eigen::SparseVector<double> R_scale = F_ * step_end - Q_;
        applyBoundaryConditions(K_, R_scale);
        double normR_scale = R_scale.norm();

        std::cout << "Iteration " << iter << ": " << normR_scale << std::endl;

        if (normR_scale < tol)
        {
            std::cout << "Converged! Total iterations: " << iter << "\n"
                      << "Final residual norm: " << normR_scale << std::endl;
            return true;
        }

        Spinner spinner(
            "\033[1;32m[Process]\033[0m Step size = " + std::to_string(step_end) + ". Iteration " + std::to_string(iter+1) + ": ",
            "Step size = " + std::to_string(step_end) + ". Iteration time: ",
            100);
        
        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        Eigen::VectorXd du;
        solver.compute(K_);
        du = solver.solve(R_scale);
        U_ += du;
        DU_ = du;
    }
    // When the maximum number of iterations is reached, the solution is not converged
    // Recover the displacement vector to the previous result
    U_ = U_backup; 
    // When the maximum number of iterations is reached, the solution is not converged
    // Recover the stress and strain tensor to the previous result
    for (int i = 0; i < elements.size(); i++)
    {
        auto& elem = elements[i];
        for (int j = 0; j < 8; j++)
        {
            elem.materialPoints[j].strain = tmp_elem_strain[i].col(j);
            elem.materialPoints[j].stress = tmp_elem_stress[i].col(j);
        }

    }
    return false;  // Convergence failed
}

void FiniteElementSolver::writeToFile()
{
    std::string current_path = std::filesystem::current_path().string();
    std::string resultpath = current_path + "\\.." + "\\.." + "\\FEoutput";
    if (std::filesystem::exists(resultpath)){
        std::filesystem::remove_all(resultpath);
    }
    std::filesystem::create_directory(resultpath);  
    saveMatrix2TXT(U_, resultpath, "Displacement.txt", 12);

    getAllElementStress(true);
    std::vector<std::string> filenames = {"Strain.txt", "Stress.txt", "PrincipalStrain.txt", "PrincipalStress.txt"};
    for (int i = 0; i < 4; i++)
    {
        avgFieldAtNodes(i, filenames[i], true);
    }
}


void FiniteElementSolver::solve_adaptive_nonlinear(double& step_size, int& maxIter)
{
                                            // Maximum number of iterations
    double tol = 1e-2;                      // Tolerance for convergence
    double scaleFactor = 1.0;               // load factor
    const double scaleFactorDecay = 0.5;    // reduction factor
    double minScaleFactor = 1e-3;           // Minimum scale factor
                                            // step size each iteration

    double min_step_size = step_size;
    double step_size_solved = 0.0;

    R_ = F_; // Initial residual force

    while (R_.norm() > tol)
    {
        double step_end = min_step_size + step_size_solved;
        bool converged = perform_Newton_Raphson(maxIter, tol, step_end);

        if (converged)
        {
            R_ = F_ - Q_;   // Calculate the new residual force
            applyBoundaryConditions(K_, R_);
            step_size_solved = step_end;
        } else if (!converged)
        {
            scaleFactor *= scaleFactorDecay;
            if (scaleFactor * step_size < min_step_size)
            {
                min_step_size = scaleFactor * step_size;
            }

            if (scaleFactor < minScaleFactor)
            {
                std::cerr << "Error: Solution did not converge even with minimal scale factor!" << std::endl;
                return;
            }
            std::cout << "Solution not converged, reducing load factor to: " << scaleFactor << std::endl;
        }
    }

    writeToFile(); 
}
