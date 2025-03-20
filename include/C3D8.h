#ifndef C3D8_H
#define C3D8_H

#include <Eigen/Sparse>
#include <vector>
#include "ElemBase.h"
#include "Math.h"

#include "../material/LinearElastic.h"
#include "../material/IdealElastoplastic.h"


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
    double epstrain;

    MaterialPoint(){
        r = s = t = 0.0;
        weight = 0.0;

        detJ = 0.0;
        B = Eigen::MatrixXd::Zero(6, 24);
        Bbar = Eigen::MatrixXd::Zero(6, 24);

        strain = Eigen::VectorXd::Zero(6);
        stress = Eigen::VectorXd::Zero(6);
        epstrain = 0.0;
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

    std::shared_ptr<Material> material_;

public:
    C3D8(int elemID, std::shared_ptr<Material> material, int matID = -1);

    // void setMaterial(LinearElasticMaterial* material) { material_ = material; };
    void initMaterialPoints();
    void updateTangentAndInternal(const Eigen::VectorXd& u, Eigen::VectorXd& elemQ, Eigen::MatrixXd& elemK);

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
    Eigen::MatrixXd extrapolateTensor(const Eigen::MatrixXd& tensor_at_Gpoints);
};

#endif // C3D8_H