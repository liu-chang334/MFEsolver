#ifndef MATERIAL_H
#define MATERIAL_H

#include <Eigen/Dense>
#include <vector>


class Material
{
public:
    Material(){}
    virtual ~Material(){}
    virtual void updateStressAndTangent(
        const Eigen::VectorXd &strain,
        Eigen::VectorXd &stress,
        Eigen::MatrixXd &tangent) = 0;
};

#endif // MATERIAL_H