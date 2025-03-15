#ifndef MATERIAL_H
#define MATERIAL_H

#include <iostream>
#include <Eigen/Dense>
#include <unordered_map>
#include <string>

struct MaterialPoint; // Forward declaration

class Material
{
public:
    virtual ~Material() = default;

    /**
     * @brief Set the material parameters
     *
     * @param[in] parameters: The material parameters
     * @note The material parameters are stored in a std::unordered_map<std::string, double>
     */
    virtual void setMaterialParameters(const std::unordered_map<std::string, double>& parameters) = 0;

    /**
     * @brief Update the stress and tangent
     *
     * @param[in] dstrain: The current strain increment
     * @param[out] mp: The material point at the end of the time step
     * @param[out] tangent: The tangent stiffness matrix at the end of the time step
     */
    virtual void updateStressAndTangent(const Eigen::VectorXd &dstrain, MaterialPoint& mp, Eigen::MatrixXd &tangent) = 0;
protected:
    std::unordered_map<std::string, double> materialParams_;

};


#endif