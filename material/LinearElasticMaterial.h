#ifndef LINEAR_ELASTIC_MATERIAL_H
#define LINEAR_ELASTIC_MATERIAL_H

#include <Eigen/Dense>
#include <vector>

class LinearElasticMaterial
{
public:
    LinearElasticMaterial();
    void setproperties(double E, double nu);

    void updateStressAndTangent(const Eigen::VectorXd &strain, Eigen::VectorXd &stress, Eigen::MatrixXd &tangent);

public:
    double E_;
    double nu_;
    Eigen::MatrixXd D_;
};

#endif // LINEAR_ELASTIC_MATERIAL_H