#include "../include/C3D8.h"


/**
 * @brief Constructor for the C3D8 class.
 *
 * @param[in] elemID The ID of the element.
 * @param[in] matID The ID of the material.
 * @note The order of the gaussian points is defined by right-hand rule.
 *       The bottom points are (-1, -1, -1), (1, -1, -1), (1, 1, -1), (-1, 1, -1),
 *       the top points are    (-1, -1,  1), (1, -1,  1), (1, 1,  1), (-1, 1,  1)
 * @note The order of the gaussian weights in Abaqus is defined as follows:
 *       (-1, -1, -1), (1, -1, -1), (-1, 1, -1), (1, 1, -1),
 *       (-1, -1,  1), (1, -1,  1), (-1, 1,  1), (1, 1,  1)
 */
C3D8::C3D8(int elemID, int matID) : SolidElement(elemID, matID), extrapolateMatrix(8,8){
    const double Tpp = (5.0 + 3.0 * sqrt(3.0)) / 4.0;
    const double Tpm = (5.0 - 3.0 * sqrt(3.0)) / 4.0;
    const double Qpm = (-1.0 + sqrt(3.0)) / 4.0;
    const double Qmm = (-1.0 - sqrt(3.0)) / 4.0;
    extrapolateMatrix << 
        /* Row 0 */ Tpp, Qmm, Qpm, Qmm, Qmm, Qpm, Tpm, Qpm,
        /* Row 1 */ Qmm, Tpp, Qmm, Qpm, Qpm, Qmm, Qpm, Tpm,
        /* Row 2 */ Qpm, Qmm, Tpp, Qmm, Tpm, Qpm, Qmm, Qpm,
        /* Row 3 */ Qmm, Qpm, Qmm, Tpp, Qpm, Tpm, Qpm, Qmm,
        /* Row 4 */ Qmm, Qpm, Tpm, Qpm, Tpp, Qmm, Qpm, Qmm,
        /* Row 5 */ Qpm, Qmm, Qpm, Tpm, Qmm, Tpp, Qmm, Qpm,
        /* Row 6 */ Tpm, Qpm, Qmm, Qpm, Qpm, Qmm, Tpp, Qmm,
        /* Row 7 */ Qpm, Tpm, Qpm, Qmm, Qmm, Qpm, Qmm, Tpp;
}

void C3D8::initMaterialPoints()
{
    const double g = 1.0 / sqrt(3.0);
    Eigen::MatrixXd gaussPoints(8,3);
    Eigen::VectorXd gaussWeights(8);
    gaussPoints << -g, -g, -g,
                    g, -g, -g,
                    g,  g, -g,
                   -g,  g, -g,
                   -g, -g,  g,
                    g, -g,  g,
                    g,  g,  g,
                   -g,  g,  g;
    const double w = 1.0;
    gaussWeights << w * w * w, w * w * w, w * w * w, w * w * w, w * w * w, w * w * w, w * w * w, w * w * w;

    materialPoints.resize(8);
    for (int i = 0; i < 8; ++i)
    {
        materialPoints[i].r = gaussPoints(i, 0);
        materialPoints[i].s = gaussPoints(i, 1);
        materialPoints[i].t = gaussPoints(i, 2);
        materialPoints[i].weight = gaussWeights(i);
    }

    for (int i = 0; i < 8; ++i)
    {
        Eigen::MatrixXd Nrst_diff = calcuShapeFunctionDerivatives(materialPoints[i].r, materialPoints[i].s, materialPoints[i].t);
        Eigen::MatrixXd J = calcuJacobian(Nrst_diff);
        materialPoints[i].J_inv = calcuJacobianInverse(J);
        materialPoints[i].detJ = J.determinant();
        materialPoints[i].B = calcuBMatrix(Nrst_diff, materialPoints[i].J_inv);

        materialPoints[i].strain = Eigen::VectorXd::Zero(6);
        materialPoints[i].stress = Eigen::VectorXd::Zero(6);
    }

    calcuBMatrixVolumetricAverage();
    for (int i = 0; i < 8; ++i)
    {
        materialPoints[i].Bbar = calcuBMatrixCorrected(materialPoints[i].B);
    }
}

/**
 * @brief Calculate the matrix of shape functions at a given gaussian point.
 *
 * @param[in] r The r-coordinate of the gaussian point.
 * @param[in] s The s-coordinate of the gaussian point.
 * @param[in] t The t-coordinate of the gaussian point.
 * @return The matrix of shape functions.
 * @note The order of the shape functions is defined by right-hand rule.
 *       The bottom points are (-1, -1, -1), (1, -1, -1), (1, 1, -1), (-1, 1, -1), 
 *       the top points are    (-1, -1,  1), (1, -1,  1), (1, 1,  1), (-1, 1,  1)
 */
Eigen::VectorXd C3D8::calcuShapeFunctions(double r, double s, double t)
{
    Eigen::VectorXd N(8);

    N << 0.125 * (1 - r) * (1 - s) * (1 - t),
         0.125 * (1 + r) * (1 - s) * (1 - t),
         0.125 * (1 + r) * (1 + s) * (1 - t),
         0.125 * (1 - r) * (1 + s) * (1 - t),
         0.125 * (1 - r) * (1 - s) * (1 + t),
         0.125 * (1 + r) * (1 - s) * (1 + t),
         0.125 * (1 + r) * (1 + s) * (1 + t),
         0.125 * (1 - r) * (1 + s) * (1 + t);

    return N; 
}

/**
 * @brief Calculate the matrix of shape function derivatives at a given gaussian point.
 *
 * @param[in] r The r-coordinate of the gaussian point.
 * @param[in] s The s-coordinate of the gaussian point.
 * @param[in] t The t-coordinate of the gaussian point.
 * @return The matrix of shape function derivatives.
 * @note The order of the shape functions is defined by right-hand rule.
 *       The bottom points are (-1, -1, -1), (1, -1, -1), (1, 1, -1), (-1, 1, -1),
 *       the top points are    (-1, -1,  1), (1, -1,  1), (1, 1,  1), (-1, 1,  1)
 */
Eigen::MatrixXd C3D8::calcuShapeFunctionDerivatives(double r, double s, double t) 
{
    // Nrst_diff: 3x8 matrix
    Eigen::MatrixXd Nrst_diff(3, 8);

    Nrst_diff << 
        -0.125 * (1 - s) * (1 - t),  0.125 * (1 - s) * (1 - t),  0.125 * (1 + s) * (1 - t), -0.125 * (1 + s) * (1 - t),
        -0.125 * (1 - s) * (1 + t),  0.125 * (1 - s) * (1 + t),  0.125 * (1 + s) * (1 + t), -0.125 * (1 + s) * (1 + t),
        
        -0.125 * (1 - r) * (1 - t), -0.125 * (1 + r) * (1 - t),  0.125 * (1 + r) * (1 - t),  0.125 * (1 - r) * (1 - t),
        -0.125 * (1 - r) * (1 + t), -0.125 * (1 + r) * (1 + t),  0.125 * (1 + r) * (1 + t),  0.125 * (1 - r) * (1 + t),

        -0.125 * (1 - r) * (1 - s), -0.125 * (1 + r) * (1 - s), -0.125 * (1 + r) * (1 + s), -0.125 * (1 - r) * (1 + s),
         0.125 * (1 - r) * (1 - s),  0.125 * (1 + r) * (1 - s),  0.125 * (1 + r) * (1 + s),  0.125 * (1 - r) * (1 + s);

    return Nrst_diff;
}

/**
 * @brief Calculate the Jacobian matrix at a given gaussian point.
 *
 * @param[in] Nrst_diff The matrix of shape function derivatives.
 * @return The Jacobian matrix 3x3.
 */
Eigen::MatrixXd C3D8::calcuJacobian(const Eigen::MatrixXd& Nrst_diff) 
{
    if (Nrst_diff.rows() != 3 || Nrst_diff.cols() != 8) {
        throw std::invalid_argument("Nrst_diff must be a 3x8 matrix.");
    }
    if (nodeCoorMatrix.rows() != 8 || nodeCoorMatrix.cols() != 3) {
        throw std::invalid_argument("Node coordinates (xyz) must be an 8x3 matrix.");
    }
    Eigen::MatrixXd J = Nrst_diff * nodeCoorMatrix;
    return J; 
}

/**
 * @brief Calculate the inverse of the Jacobian matrix.
 *
 * @param[in] J The Jacobian matrix 3x3.
 * @return The inverse of the Jacobian matrix 3x3.
 */
Eigen::MatrixXd C3D8::calcuJacobianInverse(const Eigen::MatrixXd& J) 
{
    if (J.rows() != 3 || J.cols() != 3) {
        throw std::invalid_argument("Jacobian matrix must be a 3x3 matrix.");
    }

    Eigen::MatrixXd J_inv = J.inverse();

    return J_inv;
}

/**
 * @brief Calculate the B matrix at a given gaussian point.
 *
 * @param[in] r The r-coordinate of the gaussian point.
 * @param[in] s The s-coordinate of the gaussian point.
 * @param[in] t The t-coordinate of the gaussian point.
 * @return The B matrix 6x24.
 * @note The order of the strains is defined as follows:
 *       0: e_11， 1: e_22， 2: e_33， 3: 2*e_12， 4: 2*e_13， 5: 2*e_23
 */
Eigen::MatrixXd C3D8::calcuBMatrix(const Eigen::MatrixXd& Nrst_diff, const Eigen::MatrixXd& J_inv) 
{
    Eigen::MatrixXd B(6, 24);
    B.setZero();

    for (int i = 0; i < 8; i++) {
        B(0, 3 * i) = J_inv(0, 0) * Nrst_diff(0, i) + J_inv(0, 1) * Nrst_diff(1, i) + J_inv(0, 2) * Nrst_diff(2, i);
        B(0, 3 * i + 1) = 0.0;
        B(0, 3 * i + 2) = 0.0;

        B(1, 3 * i) = 0.0;
        B(1, 3 * i + 1) = J_inv(1, 0) * Nrst_diff(0, i) + J_inv(1, 1) * Nrst_diff(1, i) + J_inv(1, 2) * Nrst_diff(2, i);
        B(1, 3 * i + 2) = 0.0;

        B(2, 3 * i) = 0.0;
        B(2, 3 * i + 1) = 0.0;
        B(2, 3 * i + 2) = J_inv(2, 0) * Nrst_diff(0, i) + J_inv(2, 1) * Nrst_diff(1, i) + J_inv(2, 2) * Nrst_diff(2, i);

        B(3, 3 * i) = J_inv(1, 0) * Nrst_diff(0, i) + J_inv(1, 1) * Nrst_diff(1, i) + J_inv(1, 2) * Nrst_diff(2, i);
        B(3, 3 * i + 1) = J_inv(0, 0) * Nrst_diff(0, i) + J_inv(0, 1) * Nrst_diff(1, i) + J_inv(0, 2) * Nrst_diff(2, i);
        B(3, 3 * i + 2) = 0.0;

        B(4, 3 * i) = J_inv(2, 0) * Nrst_diff(0, i) + J_inv(2, 1) * Nrst_diff(1, i) + J_inv(2, 2) * Nrst_diff(2, i);
        B(4, 3 * i + 1) = 0.0;
        B(4, 3 * i + 2) = J_inv(0, 0) * Nrst_diff(0, i) + J_inv(0, 1) * Nrst_diff(1, i) + J_inv(0, 2) * Nrst_diff(2, i);

        B(5, 3 * i) = 0.0;
        B(5, 3 * i + 1) = J_inv(2, 0) * Nrst_diff(0, i) + J_inv(2, 1) * Nrst_diff(1, i) + J_inv(2, 2) * Nrst_diff(2, i);
        B(5, 3 * i + 2) = J_inv(1, 0) * Nrst_diff(0, i) + J_inv(1, 1) * Nrst_diff(1, i) + J_inv(1, 2) * Nrst_diff(2, i);
    }

    return B; 
}

/**
 * @brief Calculate the volumetric B matrix at a given gaussian point.
 *
 * @param[in] r The r-coordinate of the gaussian point.
 * @param[in] s The s-coordinate of the gaussian point.
 * @param[in] t The t-coordinate of the gaussian point.
 * @return The volumetric B matrix 6x24.
 * @note The order of the strains is defined as follows:
 *       0: e_11， 1: e_22， 2: e_33， 3: 2*e_12， 4: 2*e_13， 5: 2*e_23
 */
Eigen::MatrixXd C3D8::calcuBMatrixVolumetric(const Eigen::MatrixXd& B)
{
    Eigen::MatrixXd Bvol = Eigen::MatrixXd::Zero(6, 24);
    for (int i = 0; i < 8; ++i)
    {
        int col_start = 3 *i;
        Eigen::MatrixXd BDil = Eigen::MatrixXd::Zero(6, 3);
    
        for (int j = 0; j < 3; ++j)
        {
            BDil(0, j) = B(j, col_start + j);
            BDil(1, j) = B(j, col_start + j);
            BDil(2, j) = B(j, col_start + j);
            BDil(3, j) = 0.0;
            BDil(4, j) = 0.0;
            BDil(5, j) = 0.0;
        }
        Bvol.block<6, 3>(0, col_start) = BDil / 3.0;
    }

    return Bvol;
}

/**
 * @brief Calculate the average volumetric B matrix.
 *
 * This function calculates the average volumetric B matrix by integrating over the element's
 * gaussian points. The result is stored in the member variable Bvol_average.
 */
void C3D8::calcuBMatrixVolumetricAverage()
{
    Eigen::MatrixXd B_vol_integrate = Eigen::MatrixXd::Zero(6, 24);
    for (int i = 0; i < 8; ++i)
    {
        double r = materialPoints[i].r;
        double s = materialPoints[i].s;
        double t = materialPoints[i].t;
        Eigen::MatrixXd B = materialPoints[i].B;
        double detJ = materialPoints[i].detJ;
        
        Eigen::MatrixXd Bvol = calcuBMatrixVolumetric(B);
        B_vol_integrate += Bvol * detJ * materialPoints[i].weight;
    }

    double elemvol = calcuVolume();
    Eigen::MatrixXd B_vol_avg = B_vol_integrate / elemvol; 
    Bvol_average = B_vol_avg;   
}

/**
 * @brief Calculate the corrected B matrix at a given gaussian point.
 *
 * @param[in] r The r-coordinate of the gaussian point.
 * @param[in] s The s-coordinate of the gaussian point.
 * @param[in] t The t-coordinate of the gaussian point.
 * @return The corrected B matrix 6x24.
 * @note The order of the strains is defined as follows:
 *       0: e_11， 1: e_22， 2: e_33， 3: 2*e_12， 4: 2*e_13， 5: 2*e_23
 */
Eigen::MatrixXd C3D8::calcuBMatrixCorrected(const Eigen::MatrixXd& B)
{
    // Eigen::MatrixXd B_vol_integrate = Eigen::MatrixXd::Zero(6, 24);

    // //gauss integration
    // for (int i = 0; i < 8; ++i)
    // {
    //     double s = gaussPoints(i, 0);
    //     double t = gaussPoints(i, 1);
    //     double r = gaussPoints(i, 2);

    //     Eigen::MatrixXd J = calcuJacobian(calcuShapeFunctionDerivatives(r, s, t));
    //     Eigen::MatrixXd Bvol = calcuBMatrixVolumetric(r, s, t);

    //     double detJ = J.determinant();
    //     B_vol_integrate += Bvol * detJ * gaussWeights(i);
    // }

    // double elemvol = calcuVolume();
    // Eigen::MatrixXd B_vol_avg = B_vol_integrate / elemvol;
    Eigen::MatrixXd Bvol = calcuBMatrixVolumetric(B);
    Eigen::MatrixXd Bbar = B - Bvol + Bvol_average;

    return Bbar;
}

/**
 * @brief Calculate the deviatoric B matrix at a given gaussian point.
 *
 * @param[in] r The r-coordinate of the gaussian point.
 * @param[in] s The s-coordinate of the gaussian point.
 * @param[in] t The t-coordinate of the gaussian point.
 * @return The deviatoric B matrix 6x24.
 * @note The order of the strains is defined as follows:
 *       0: e_11， 1: e_22， 2: e_33， 3: 2*e_12， 4: 2*e_13， 5: 2*e_23
 */
Eigen::MatrixXd C3D8::calcuBMatrixDeviatoric(const Eigen::MatrixXd& Nrst_diff, const Eigen::MatrixXd& J_inv)
{
    Eigen::MatrixXd B = calcuBMatrix(Nrst_diff, J_inv);
    Eigen::MatrixXd Bvol = calcuBMatrixVolumetric(B);
    Eigen::MatrixXd Bdev = B - Bvol;

    return Bdev;
}

/**
 * @brief Calculate the volume of the element.
 *
 * @return The volume of the element.
 */
double C3D8::calcuVolume()
{
    double elemvol = 0.0;
    for (int i = 0; i < 8; i++) {
        double detJ = materialPoints[i].detJ;
        elemvol += detJ * materialPoints[i].weight;
    }
    return elemvol;
}

/**
 * @brief Calculate the constitutive matrix for the element.
 *
 * @return The constitutive matrix --De 6x6.
 * @note Only for isotropic elastic material.
 */
Eigen::MatrixXd C3D8::calcuConstitutiveMatrix()
{
    Eigen::MatrixXd D = Eigen::MatrixXd::Zero(6, 6);
    D << lambda + 2 * mu, lambda, lambda, 0, 0, 0,
        lambda, lambda + 2 * mu, lambda, 0, 0, 0,
        lambda, lambda, lambda + 2 * mu, 0, 0, 0,
        0, 0, 0, mu, 0, 0,
        0, 0, 0, 0, mu, 0,
        0, 0, 0, 0, 0, mu;

    return D;
}

/**
 * @brief Calculate the stiffness matrix for the element.
 *
 * @return The stiffness matrix --Ke 24x24.
 * @note Only for isotropic elastic material.
 */
Eigen::MatrixXd C3D8::calcuStiffnessMatrix()
{   
    /*Something wrong with this method
    // Eigen::MatrixXd D = calcuConstitutiveMatrix();
    // Eigen::MatrixXd Ke = Eigen::MatrixXd::Zero(24, 24);
    // for (int i = 0; i < 8; i++) {
    //     double r = gaussPoints(i, 0);
    //     double s = gaussPoints(i, 1);
    //     double t = gaussPoints(i, 2);
    //     Eigen::MatrixXd Bdev = calcuBMatrixDeviatoric(r, s, t);
    //     Eigen::MatrixXd J = calcuJacobian(calcuShapeFunctionDerivatives(r, s, t));
    //     double detJ = J.determinant();
    //     Ke += Bdev.transpose() * D * Bdev * detJ * gaussWeights(i);
    // }
    // double r = 0.0; double s = 0.0; double t = 0.0; double w = 8.0;
    // Eigen::MatrixXd Bvol = calcuBMatrixVolumetric(r, s, t);
    // Eigen::MatrixXd J = calcuJacobian(calcuShapeFunctionDerivatives(r, s, t));
    // double detJ = J.determinant();
    // Ke += Bvol.transpose() * D * Bvol * detJ * w;
    */
    // calcuBMatrixVolumetricAverage();
    Eigen::MatrixXd Ke = Eigen::MatrixXd::Zero(24, 24);
    Eigen::MatrixXd D = calcuConstitutiveMatrix();
    for (int i = 0; i < 8; i++) {
        Eigen::MatrixXd Bbar = materialPoints[i].Bbar;
        double detJ = materialPoints[i].detJ;
        Ke += Bbar.transpose() * D * Bbar * detJ * materialPoints[i].weight;
    }

    return Ke;
}

/**
 * @brief Calculate the strain tensor at all gauss points of the element.
 *
 * @param[in] u The displacement vector at all nodes of the element.
 * @return The strain tensor 6x8.
 * @note The order of the strains is defined as follows:
 *       0: e_11， 1: e_22， 2: e_33， 3: 2*e_12， 4: 2*e_13， 5: 2*e_23
 */
Eigen::MatrixXd C3D8::calcuStrainTensor(const Eigen::VectorXd& u)
{
    calcuBMatrixVolumetricAverage();
    Eigen::MatrixXd strain = Eigen::MatrixXd::Zero(6, 8);
    for(int i = 0; i < 8; i++)
    {
        Eigen::MatrixXd Bbar = materialPoints[i].Bbar;
        strain.col(i) = Bbar * u;
    }
    return strain;
}

/**
 * @brief Calculate the stress tensor at all gauss points of the element.
 *
 * @param[in] u The displacement vector at all nodes of the element.
 * @return The stress tensor 6x8.
 * @note The order of the strains is defined as follows:
 *       0: s_11， 1: s_22， 2: s_33， 3: s_12， 4: s_13， 5: s_23
 */
Eigen::MatrixXd C3D8::calcuStressTensor(const Eigen::VectorXd& u)
{
    Eigen::MatrixXd D = calcuConstitutiveMatrix();
    Eigen::MatrixXd strain = calcuStrainTensor(u);
    Eigen::MatrixXd stress = D * strain;
    return stress;
}

std::pair<Eigen::MatrixXd, Eigen::MatrixXd> C3D8::calcuStrainStressTensor(const Eigen::VectorXd& u)
{
    Eigen::MatrixXd D = calcuConstitutiveMatrix();
    Eigen::MatrixXd strain = calcuStrainTensor(u);
    Eigen::MatrixXd stress = D * strain;
    return std::make_pair(strain, stress); 
}

void C3D8::calcuTangentAndResidual(const Eigen::VectorXd& u, Eigen::VectorXd& elemQ, Eigen::MatrixXd& elemK)
{
    elemQ = Eigen::VectorXd::Zero(24);
    elemK = Eigen::MatrixXd::Zero(24, 24);    

    for (auto& mp : materialPoints)
    {
        Eigen::MatrixXd Bbar = mp.Bbar;
        double detJ = mp.detJ;
        Eigen::VectorXd strain = Bbar * u;

        Eigen::VectorXd stress;
        Eigen::MatrixXd tangent;
        material_.updateStressAndTangent(strain, stress, tangent);

        elemQ += Bbar.transpose() * stress * detJ * mp.weight;
        elemK += Bbar.transpose() * tangent * Bbar * detJ * mp.weight;

        // std::cout << "elemQ: " << elemQ.transpose() << std::endl;
 
        // update material point
        mp.strain = strain;
        mp.stress = stress;
    }
}


/**
 * @brief Interpolate a tensor from gauss points to nodes.
 *
 * @param[in] tensor_at_Gpoints The tensor at gauss points 6x8.
 * @return The tensor at nodes 6x8.
 */
Eigen::MatrixXd C3D8::extrapolateTensor(const Eigen::MatrixXd& tensor_at_Gpoints)
{
    Eigen::MatrixXd tensor_at_nodes = Eigen::MatrixXd::Zero(6, 8);
    tensor_at_nodes = (extrapolateMatrix * tensor_at_Gpoints.transpose()).transpose();
    return tensor_at_nodes;
}

// Class C3D8 End