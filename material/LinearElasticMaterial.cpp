#include "LinearElasticMaterial.h"

void LinearElasticMaterial::setMaterialParameters(const std::unordered_map<std::string, double>& parameters)
{
    materialParams_ = parameters;

    if (materialParams_.find("E") == materialParams_.end() || materialParams_.find("nu") == materialParams_.end())
    {
        std::cerr << "Error: Material parameters not set" << std::endl;
        return;
    }

	E_ = materialParams_["E"];
	nu_ = materialParams_["nu"];

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
void LinearElasticMaterial::updateStressAndTangent(const Eigen::VectorXd &dstrain, 
                                Eigen::VectorXd &dstress, Eigen::MatrixXd &tangent)
{
    dstress = D_ * dstrain;
    tangent = D_;
}