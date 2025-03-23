#include "IdealElastoplastic.h"

void IdealElastoplastic::setMaterialParameters(const MaterialData& matdata)
{
    E_ = matdata.E;
    nu_ = matdata.nu;
    // H_ = matdata.H;
    // sigma_y_ = matdata.sigma_y;

    double lambda_ = E_ * nu_ / ((1.0 + nu_) * (1.0 - 2.0 * nu_));
    double mu_ = E_ / (2.0 * (1.0 + nu_));

    D_ = Eigen::MatrixXd::Zero(6, 6);
    D_ << lambda_ + 2 * mu_, lambda_, lambda_, 0, 0, 0,
        lambda_, lambda_ + 2 * mu_, lambda_, 0, 0, 0,
        lambda_, lambda_, lambda_ + 2 * mu_, 0, 0, 0,
        0, 0, 0, mu_, 0, 0,
        0, 0, 0, 0, mu_, 0,
        0, 0, 0, 0, 0, mu_;

}


void IdealElastoplastic::updateStressAndTangent(const Eigen::VectorXd &dstrain, MaterialPoint& mp, Eigen::MatrixXd &tangent)
{
    Eigen::VectorXd dstress(6), trial_stress(6);
    dstress = D_ * dstrain;
    trial_stress = mp.stress + dstress;
    
    // Calculate the deviatoric stress
    Eigen::VectorXd dev_stress = trial_stress;
    double hydrostatic_stress = (trial_stress(0) + trial_stress(1) + trial_stress(2)) / 3.0;
    dev_stress(0) -= hydrostatic_stress;
    dev_stress(1) -= hydrostatic_stress;
    dev_stress(2) -= hydrostatic_stress;

    // Calculate the Von Mises stress
    double von_mises_stress = sqrt(1.5 * dev_stress.transpose() * dev_stress);

    if (von_mises_stress <= sigma_y_){
        // Elastic behavior
        mp.stress = trial_stress;
        tangent = D_;
    }else{
        // Plastic behavior
        double plastic_multiplier = (von_mises_stress - sigma_y_) / (3.0 * (E_ / (2.0 * (1 + nu_)) + H_));
        Eigen::VectorXd plastic_correction = (plastic_multiplier / von_mises_stress) * dev_stress;
        mp.stress = trial_stress - plastic_correction * 2.0 * (E_ / (2.0 * (1 + nu_)));

        // Calculate the tangent matrix
        double hardening_factor = (H_ / (H_ + 3.0 * (E_ / (2.0 * (1 + nu_)))));
        tangent = (1.0 - hardening_factor) * D_;
    }

    mp.strain += dstrain;
}