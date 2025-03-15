#ifndef LINEAR_ELASTIC_H
#define LINEAR_ELASTIC_H

#include "Material.h"
#include "../include/C3D8.h"



class LinearElastic : public Material
{
public:
    LinearElastic() = default;

    /**
     * @brief Set the material parameters
     *
     * @param[in] parameters: The material parameters
     * @note The material parameters are stored in a std::unordered_map<std::string, double>
     */
    void setMaterialParameters(const std::unordered_map<std::string, double>& parameters) override;


    /**
     * @brief Update the stress and tangent
     *
     * @param[in] dstrain: The current strain increment
     * @param[out] mp: The material point at the end of the time step
     * @param[out] tangent: The tangent stiffness matrix at the end of the time step
     */    
    void updateStressAndTangent(const Eigen::VectorXd &dstrain, MaterialPoint& mp, Eigen::MatrixXd &tangent) override;
public:
    double E_;
    double nu_;
    Eigen::MatrixXd D_;
};

#endif // LINEAR_ELASTIC_H