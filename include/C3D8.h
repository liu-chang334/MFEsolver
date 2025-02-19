#ifndef C3D8_H
#define C3D8_H

#include "ElemBase.h"
#include <Eigen/Sparse>

/**
 * @class C3D8
 * @brief This class inherits from the SolidElement class and provides specific implementations
 *        for the C3D8 element type.
 * - \c gaussPoints: Gauss points for integration
 * - \c gaussWeights: Gauss weights for integration
 * - \c Bvol_average: Average Bvol matrix
 */
class C3D8 : public SolidElement
{
public:
    Eigen::MatrixXd gaussPoints;
    Eigen::MatrixXd gaussWeights;
    Eigen::MatrixXd Bvol_average;

public:
    C3D8(int elemID, int matID = -1);

    Eigen::VectorXd calcuShapeFunctions(double r, double s, double t);
    Eigen::MatrixXd calcuShapeFunctionDerivatives(double r, double s, double t);
    Eigen::MatrixXd calcuJacobian(const Eigen::MatrixXd& Nrst_diff);
    Eigen::MatrixXd calcuJacobianInverse(const Eigen::MatrixXd& J);
    Eigen::MatrixXd calcuBMatrix(double r, double s, double t);
    Eigen::MatrixXd calcuBMatrixVolumetric(double r, double s, double t);
    void calcuBMatrixVolumetricAverage();
    Eigen::MatrixXd C3D8::calcuBMatrixCorrected(double r, double s, double t);
    Eigen::MatrixXd calcuBMatrixDeviatoric(double r, double s, double t);
    double calcuVolume();
    Eigen::MatrixXd calcuConstitutiveMatrix();
    Eigen::MatrixXd calcuStiffnessMatrix();
    Eigen::MatrixXd calcuStrainTensor(const Eigen::VectorXd& u);
    Eigen::MatrixXd calcuStressTensor(const Eigen::VectorXd& u);
    Eigen::MatrixXd interpolateTensor(const Eigen::MatrixXd& tensor_at_Gpoints);
};

#endif // C3D8_H