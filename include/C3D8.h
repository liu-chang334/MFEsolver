#ifndef C3D8_H
#define C3D8_H

#include "ElemBase.h"
#include "../material/LinearElasticMaterial.h"
#include <Eigen/Sparse>
#include <vector>


struct MaterialPoint
{
    double r, s, t;
    double weight;

    double detJ;
    Eigen::MatrixXd J_inv;
    Eigen::MatrixXd B;
    Eigen::MatrixXd Bbar;

    Eigen::VectorXd strain;
    Eigen::VectorXd stress;

    MaterialPoint(){
        r = s = t = 0.0;
        weight = 0.0;

        detJ = 0.0;
        B = Eigen::MatrixXd::Zero(6, 24);
        Bbar = Eigen::MatrixXd::Zero(6, 24);

        strain = Eigen::VectorXd::Zero(6);
        stress = Eigen::VectorXd::Zero(6);
    }
};

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

    Eigen::MatrixXd Bvol_average;
    Eigen::MatrixXd extrapolateMatrix;

    std::vector<MaterialPoint> materialPoints;

    std::unique_ptr<Material> material_;

public:
    C3D8(int elemID, int matID = -1);

    void setMaterial(std::unique_ptr<Material> mat) {material_ = std::move(mat);};
    void initMaterialPoints();
    void calcuTangentAndResidual(const Eigen::VectorXd& u, Eigen::VectorXd& elemQ, Eigen::MatrixXd& elemK);

    Eigen::VectorXd calcuShapeFunctions(double r, double s, double t);
    Eigen::MatrixXd calcuShapeFunctionDerivatives(double r, double s, double t);
    Eigen::MatrixXd calcuJacobian(const Eigen::MatrixXd& Nrst_diff);
    Eigen::MatrixXd calcuJacobianInverse(const Eigen::MatrixXd& J);
    Eigen::MatrixXd calcuBMatrix(const Eigen::MatrixXd& Nrst_diff, const Eigen::MatrixXd& J_inv);
    Eigen::MatrixXd calcuBMatrixVolumetric(const Eigen::MatrixXd& B);
    void calcuBMatrixVolumetricAverage();
    Eigen::MatrixXd C3D8::calcuBMatrixCorrected(const Eigen::MatrixXd& B);
    Eigen::MatrixXd calcuBMatrixDeviatoric(const Eigen::MatrixXd& Nrst_diff, const Eigen::MatrixXd& J_inv);
    double calcuVolume();
    Eigen::MatrixXd calcuConstitutiveMatrix();
    Eigen::MatrixXd calcuStiffnessMatrix();
    Eigen::MatrixXd calcuStrainTensor(const Eigen::VectorXd& u);
    Eigen::MatrixXd calcuStressTensor(const Eigen::VectorXd& u);
    std::pair<Eigen::MatrixXd, Eigen::MatrixXd> calcuStrainStressTensor(const Eigen::VectorXd& u);
    Eigen::MatrixXd extrapolateTensor(const Eigen::MatrixXd& tensor_at_Gpoints);
};

#endif // C3D8_H