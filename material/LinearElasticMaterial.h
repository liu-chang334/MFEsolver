#ifndef LINEAR_ELASTIC_MATERIAL_H
#define LINEAR_ELASTIC_MATERIAL_H

#include "Material.h"

class LinearElasticMaterial : public Material
{
public:
    LinearElasticMaterial() = default;
    virtual void setMaterialParameters(const std::unordered_map<std::string, double>& parameters) override;
    virtual void updateStressAndTangent(const Eigen::VectorXd &strain, 
                                        Eigen::VectorXd &stress, Eigen::MatrixXd &tangent) override;

public:
    double E_;
    double nu_;
    Eigen::MatrixXd D_;
};

#endif // LINEAR_ELASTIC_MATERIAL_H