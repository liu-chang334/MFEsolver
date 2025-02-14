#ifndef C3D8_H
#define C3D8_H

#include "ElemBase.h"
#include <Eigen/Sparse>

class C3D8 : public SolidElement
{
public:
    Eigen::MatrixXd gaussPoints;
    Eigen::MatrixXd gaussWeights;
    Eigen::MatrixXd Bvol_average;

public:
    C3D8(int elemID, int matID = -1);

    // calculate shape functions
    Eigen::VectorXd calcuShapeFunctions(double r, double s, double t);

    // calculate shape functions derivatives
    Eigen::MatrixXd calcuShapeFunctionDerivatives(double r, double s, double t);

    // calculate Jacobian matrix
    Eigen::MatrixXd calcuJacobian(const Eigen::MatrixXd& Nrst_diff);

    // calculate Jacobian inverse matrix
    Eigen::MatrixXd calcuJacobianInverse(const Eigen::MatrixXd& J);
    
    // calculate B matrix
    Eigen::MatrixXd calcuBMatrix(double r, double s, double t);

    // calculate volumetric B matrix -- Bvol
    Eigen::MatrixXd calcuBMatrixVolumetric(double r, double s, double t);

    // calculate average Bvol matrix -- Bvol_average
    void calcuBMatrixVolumetricAverage();

    // calculate corrected B matrix -- Bbar
    Eigen::MatrixXd C3D8::calcuBMatrixCorrected(double r, double s, double t);

    // calculate deviatoric B matrix -- Bdev
    Eigen::MatrixXd calcuBMatrixDeviatoric(double r, double s, double t);

    // calculate volume of element
    double calcuVolume();

    // calculate constitutive matrix --De (only for isotropic elastic material)
    Eigen::MatrixXd calcuConstitutiveMatrix();

    // calculate stiffness matrix --Ke (only for isotropic elastic material)
    Eigen::MatrixXd calcuStiffnessMatrix();

    // calculate strain tensor at all gauss points
    Eigen::MatrixXd calcuStrainTensor(const Eigen::VectorXd& u);

    // calculate stress tensor at all gauss points
    Eigen::MatrixXd calcuStressTensor(const Eigen::VectorXd& u);

    // interpolate strain or stress at nodes
    Eigen::MatrixXd interpolateTensor(const Eigen::MatrixXd& tensor_at_Gpoints);
};

#endif // C3D8_H