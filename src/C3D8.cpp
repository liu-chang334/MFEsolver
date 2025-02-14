#include "../include/C3D8.h"

//****************C3D8 Class
C3D8::C3D8(int elemID, int matID) : SolidElement(elemID, matID), gaussPoints(8,3), gaussWeights(8,1)
{
    double g = 1.0 / sqrt(3.0);
    gaussPoints << -g, -g, -g,
                    g, -g, -g,
                   -g,  g, -g,
                    g,  g, -g,
                   -g, -g,  g,
                    g, -g,  g,
                   -g,  g,  g,
                    g,  g,  g;
                   

    double w = 1.0;
    gaussWeights << w * w * w, w * w * w, w * w * w, w * w * w, w * w * w, w * w * w, w * w * w, w * w * w;
}

// calculate shape functions values
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

// calculate shape functions derivatives
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

// calculate Jacobian matrix
Eigen::MatrixXd C3D8::calcuJacobian(const Eigen::MatrixXd& Nrst_diff) 
{
    if (Nrst_diff.rows() != 3 || Nrst_diff.cols() != 8) {
        throw std::invalid_argument("Nrst_diff must be a 3x8 matrix.");
    }
    if (nodeCoorMatrix.rows() != 8 || nodeCoorMatrix.cols() != 3) {
        throw std::invalid_argument("Node coordinates (xyz) must be an 8x3 matrix.");
    }

    Eigen::MatrixXd J = Nrst_diff * nodeCoorMatrix;

    return J; // 3x3 matrix
}
// calculate Jacobian inverse matrix
Eigen::MatrixXd C3D8::calcuJacobianInverse(const Eigen::MatrixXd& J) 
{
    if (J.rows() != 3 || J.cols() != 3) {
        throw std::invalid_argument("Jacobian matrix must be a 3x3 matrix.");
    }

    Eigen::MatrixXd J_inv = J.inverse();

    return J_inv; // 3x3 matrix
}

// calculate B matrix
Eigen::MatrixXd C3D8::calcuBMatrix(double r, double s, double t) 
{
    Eigen::MatrixXd Nrst_diff = calcuShapeFunctionDerivatives(r, s, t);
    Eigen::MatrixXd J = calcuJacobian(Nrst_diff);
    Eigen::MatrixXd J_inv = calcuJacobianInverse(J);

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

    return B; // 6x24 matrix
}
// calculate volumetric B matrix  -- Bvol
Eigen::MatrixXd C3D8::calcuBMatrixVolumetric(double r, double s, double t)
{
    Eigen::MatrixXd B = calcuBMatrix(r, s, t);
    Eigen::MatrixXd Bvol = Eigen::MatrixXd::Zero(6, 24);
    // Bvol.setZero();
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

// calculate average Bvol matrix -- Bvol_average
void C3D8::calcuBMatrixVolumetricAverage()
{
    Eigen::MatrixXd B_vol_integrate = Eigen::MatrixXd::Zero(6, 24);
    for (int i = 0; i < 8; ++i)
    {
        double s = gaussPoints(i, 0);
        double t = gaussPoints(i, 1);
        double r = gaussPoints(i, 2);

        Eigen::MatrixXd J = calcuJacobian(calcuShapeFunctionDerivatives(r, s, t));
        Eigen::MatrixXd Bvol = calcuBMatrixVolumetric(r, s, t);

        double detJ = J.determinant();
        B_vol_integrate += Bvol * detJ * gaussWeights(i);
    }

    double elemvol = calcuVolume();
    Eigen::MatrixXd B_vol_avg = B_vol_integrate / elemvol; 

    Bvol_average = B_vol_avg;   
}

// calculate corrected B matrix -- Bbar
Eigen::MatrixXd C3D8::calcuBMatrixCorrected(double r, double s, double t)
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
    Eigen::MatrixXd B = calcuBMatrix(r, s, t);
    Eigen::MatrixXd Bvol = calcuBMatrixVolumetric(r, s, t);
    Eigen::MatrixXd Bbar = B - Bvol + Bvol_average;

    return Bbar;
}


// calculate deviatoric B matrix -- Bdev
Eigen::MatrixXd C3D8::calcuBMatrixDeviatoric(double r, double s, double t)
{
    Eigen::MatrixXd B = calcuBMatrix(r, s, t);
    Eigen::MatrixXd Bvol = calcuBMatrixVolumetric(r, s, t);
    Eigen::MatrixXd Bdev = B - Bvol;

    return Bdev;
}


// calculate volume of element
double C3D8::calcuVolume()
{
    double elemvol = 0.0;
    for (int i = 0; i < gaussPoints.rows(); i++) {
        double r = gaussPoints(i, 0);
        double s = gaussPoints(i, 1);
        double t = gaussPoints(i, 2);

        Eigen::MatrixXd Nrst_diff = calcuShapeFunctionDerivatives(r, s, t);
        Eigen::MatrixXd J = calcuJacobian(Nrst_diff);

        double detJ = J.determinant();
        elemvol += detJ * gaussWeights(i);
    }
    return elemvol;
}

// calculate constitutive matrix --De (only for isotropic elastic material)
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

// calculate stiffness matrix --Ke (only for isotropic elastic material)
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
    calcuBMatrixVolumetricAverage();
    Eigen::MatrixXd Ke = Eigen::MatrixXd::Zero(24, 24);
    Eigen::MatrixXd D = calcuConstitutiveMatrix();
    for (int i = 0; i < 8; i++) {
        double r = gaussPoints(i, 0);
        double s = gaussPoints(i, 1);
        double t = gaussPoints(i, 2);
        Eigen::MatrixXd Bbar = calcuBMatrixCorrected(r, s, t);
        Eigen::MatrixXd J = calcuJacobian(calcuShapeFunctionDerivatives(r, s, t));
        double detJ = J.determinant();
        Ke += Bbar.transpose() * D * Bbar * detJ * gaussWeights(i); 
    }

    return Ke;
}

// calculate strain tensor at all gauss points
Eigen::MatrixXd C3D8::calcuStrainTensor(const Eigen::VectorXd& u)
{
    calcuBMatrixVolumetricAverage();
    Eigen::MatrixXd strain = Eigen::MatrixXd::Zero(6, 8);
    for(int i = 0; i < 8; i++)
    {
        double r = gaussPoints(i, 0);
        double s = gaussPoints(i, 1);
        double t = gaussPoints(i, 2);
        Eigen::MatrixXd Bbar = calcuBMatrixCorrected(r, s, t);
        strain.col(i) = Bbar * u;
    }
    return strain;
}

// calculate stress tensor at all gauss points
Eigen::MatrixXd C3D8::calcuStressTensor(const Eigen::VectorXd& u)
{
    Eigen::MatrixXd D = calcuConstitutiveMatrix();
    Eigen::MatrixXd strain = calcuStrainTensor(u);
    Eigen::MatrixXd stress = D * strain;
    return stress;
}

// interpolate tensor to nodes
Eigen::MatrixXd C3D8::interpolateTensor(const Eigen::MatrixXd& tensor_at_Gpoints)
{
    Eigen::MatrixXd tensor_at_nodes = Eigen::MatrixXd::Zero(6, 8);
    Eigen::MatrixXd transferMatrix = Eigen::MatrixXd::Zero(8, 8);
    for (int i = 0; i < 8; i++)
    {
        double r = gaussPoints(i, 0);
        double s = gaussPoints(i, 1);
        double t = gaussPoints(i, 2);
        Eigen::VectorXd N = calcuShapeFunctions(r, s, t);
        transferMatrix.row(i) = N.transpose();
    }
    tensor_at_nodes = (transferMatrix.inverse() * tensor_at_Gpoints.transpose()).transpose();
    return tensor_at_nodes;
}

//*******************C3D8 Class End