#ifndef LINEAR_ELASTIC_MATERIAL_H
#define LINEAR_ELASTIC_MATERIAL_H

#include "Material.h"

class LinearElasticMaterial : public Material
{
public:
    LinearElasticMaterial(double E, double nu)
    {
        E_ = E;
        nu_ = nu;
        // Calculate the Lame parameters
        double lambda_ = E_ * nu_ / ((1.0 + nu_) * (1.0 - 2.0 * nu_));
        double mu_ = E_ / (2.0 * (1.0 + nu_));

        // Calculate the stiffness matrix
        D_ = Eigen::MatrixXd::Zero(6, 6);
        D_ << lambda_ + 2 * mu_, lambda_, lambda_, 0, 0, 0,
        lambda_, lambda_ + 2 * mu_, lambda_, 0, 0, 0,
        lambda_, lambda_, lambda_ + 2 * mu_, 0, 0, 0,
        0, 0, 0, mu_, 0, 0,
        0, 0, 0, 0, mu_, 0,
        0, 0, 0, 0, 0, mu_;
    }

    /**
     * @brief Update the stress and tangent
     *
     * @param[in] strain: The current strain or strain increment
     * @param[out] stress: The current stress or stress increment
     * @param[out] tangent: The current tangent stiffness matrix
     */
    void updateStressAndTangent(
        const Eigen::VectorXd &strain,
        Eigen::VectorXd &stress,
        Eigen::MatrixXd &tangent)
    {
        // Calculate the stress
        stress = D_ * strain;

        // Calculate the tangent
        tangent = D_;
    }
    

private:
    double E_;
    double nu_;

    Eigen::MatrixXd D_;
};

#endif // LINEAR_ELASTIC_MATERIAL_H