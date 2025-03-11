#ifndef MATERIAL_H
#define MATERIAL_H

#include <iostream>
#include <Eigen/Dense>
#include <unordered_map>
#include <string>

class Material
{
public:
    virtual ~Material() = default;
    virtual void setMaterialParameters(const std::unordered_map<std::string, double>& parameters) = 0;
    virtual void updateStressAndTangent(const Eigen::VectorXd &strain, Eigen::VectorXd &stress, Eigen::MatrixXd &tangent) = 0;
protected:
    std::unordered_map<std::string, double> materialParams_;

};


#endif