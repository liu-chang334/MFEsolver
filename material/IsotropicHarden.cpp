#include "IsotropicHarden.h"

void IsotropicHarden::setMaterialParameters(const MaterialData& matdata)
{
    E_ = matdata.E;
    nu_ = matdata.nu;
    hardeningCurve_ = matdata.HardeningCurve;

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

void IsotropicHarden::updateStressAndTangent(const Eigen::VectorXd &dstrain, MaterialPoint& mp, Eigen::MatrixXd &tangent)
{
    // 1. compute the trial stress
    Eigen::VectorXd dstress(6), trial_stress(6);
    dstress = D_ * dstrain;
    trial_stress = mp.stress + dstress;

    // 2. compute the elastic factor
    double f = computerSecondInvariant3x3(mp.stress);
    double trial_f = computerSecondInvariant3x3(trial_stress);
    double yield_stress_t = interpolateHardeningStress(mp.epstrain, hardeningCurve_);
    double Ep_t = computePlasticModulus(mp.epstrain, hardeningCurve_);
    double k = yield_stress_t * yield_stress_t / 3.0;
    double m; // 0 <= m < 1

    if (trial_f - k <= 0){ // elastic behavior
        m = 1.0;
        mp.strain += dstrain;
        mp.stress = trial_stress;
        return;
    } else if (trial_f - k > 0 && f - k < 0){ // elastic and elastic-plastic behavior
        Eigen::VectorXd dstress_dev = computeDeviatoric(dstress);
        Eigen::VectorXd stress_dev = computeDeviatoric(mp.stress);
        double a0, a1, a2;
        a0 = f - k;
        a1 = stress_dev.dot(dstress_dev);
        a2 = dstress_dev.dot(dstress_dev);
        m = (-1 * a1 + sqrt(a1 * a1 - 4 * a0 * a2)) / (2 * a2);
    } else if (trial_f - k > 0 && std::abs(f - k) < 1e-9){  // elastic-plastic behavior
        m = 0.0;
    }

    // 3. accumulate elastic part to the material point
    mp.strain += m * dstrain;
    mp.stress += m * dstress;
    Eigen::VectorXd dstrain_ep = (1 - m) * dstrain;

    // 4. sub-increment method: tangential predict and then return along radial direction
    // 4.1 determine the number of sub-increments
    Eigen::VectorXd stress_dev = computeDeviatoric(mp.stress); 
    Eigen::VectorXd df_dsigma = Eigen::VectorXd::Zero(6);
    df_dsigma << stress_dev(0), stress_dev(1), stress_dev(2), 2 * stress_dev(3), 2 * stress_dev(4), 2 * stress_dev(5);
    double dlambda = (df_dsigma * D_ * dstrain_ep.transpose())(0,0) / ((df_dsigma * D_ * df_dsigma.transpose())(0,0) + yield_stress_t * yield_stress_t * Ep_t * 4.0 / 9.0);
    Eigen::VectorXd dstrain_ep_dev = computeDeviatoric(dstrain_ep);
    Eigen::VectorXd dstrain_e_dev = dstrain_e_dev - dlambda * df_dsigma;
    double eq_dstrain_e_dev = sqrt(dstrain_e_dev.dot(dstrain_e_dev) * 2 / 3);
    double M = 0.0002;
    int num_increments = static_cast<int>(std::ceil(1.0 + (eq_dstrain_e_dev / M)));
    Eigen::VectorXd dstrain_ep_i = dstrain_ep / num_increments;

    for (int i = 0; i < num_increments; i++)
    {
        double yield_stress_i = interpolateHardeningStress(mp.epstrain, hardeningCurve_);
        double Ep_i = computePlasticModulus(mp.epstrain, hardeningCurve_);

        // 4.2 compute the elastic-plastic matrix
        Eigen::VectorXd stress_dev_i = computeDeviatoric(mp.stress); 
        Eigen::VectorXd df_dsigma_i = Eigen::VectorXd::Zero(6);
        df_dsigma_i << stress_dev_i(0), stress_dev_i(1), stress_dev_i(2), 2 * stress_dev_i(3), 2 * stress_dev_i(4), 2 * stress_dev_i(5);
        Dp_ = (D_ * df_dsigma_i.transpose() * df_dsigma_i * D_) / ((df_dsigma_i * D_ * df_dsigma_i.transpose())(0,0) + yield_stress_i * yield_stress_i * Ep_i * 4.0 / 9.0);
        Eigen::MatrixXd Dep = D_ - Dp_;

        // 4.3 to predict the stress increment along tangent direction
        Eigen::VectorXd dstress_ep_i = Dep * dstrain_ep_i;
        double dlambda_i = (df_dsigma_i * D_ * dstrain_ep_i.transpose())(0,0) / ((df_dsigma_i * D_ * df_dsigma_i.transpose())(0,0) + yield_stress_i * yield_stress_i * Ep_i * 4.0 / 9.0);
        double depstrain_i = 2 / 3 * dlambda_i * yield_stress_i;
        Eigen::VectorXd stress_pre = mp.stress + dstress_ep_i;
        mp.epstrain += depstrain_i;
        mp.strain += dstrain_ep_i;

        // 4.4 return the stress along radial direction
        Eigen::VectorXd stress_pre_dev = computeDeviatoric(stress_pre);
        yield_stress_i = interpolateHardeningStress(mp.epstrain, hardeningCurve_);
        double r = sqrt((2 / 3) * yield_stress_i * yield_stress_i / (stress_pre_dev.dot(stress_pre_dev)));
        mp.stress = r * stress_pre;
    }

    // 5. update the tangent matrix
    double yield_stress = interpolateHardeningStress(mp.epstrain, hardeningCurve_);
    double Ep = computePlasticModulus(mp.epstrain, hardeningCurve_);
    stress_dev = computeDeviatoric(mp.stress); 
    df_dsigma << stress_dev(0), stress_dev(1), stress_dev(2), 2 * stress_dev(3), 2 * stress_dev(4), 2 * stress_dev(5);
    Dp_ = (D_ * df_dsigma.transpose() * df_dsigma * D_) / ((df_dsigma * D_ * df_dsigma.transpose())(0,0) + yield_stress * yield_stress * Ep * 4.0 / 9.0);
    Eigen::MatrixXd Dep = D_ - Dp_;
    tangent = Dep;
}


double interpolateHardeningStress(double plasticStrain, const std::vector<std::pair<double, double>>& curve)
{
    if (curve.empty())
    {
        std::cerr << "Error: Hardening curve is empty.\n";
        return 0.0;
    }

    // If the strain is below the first entry, return the first stress value
    if (plasticStrain <= curve.front().first)
    {
        return curve.front().second;
    }

    // If the strain is above the last entry, return the last stress value
    if (plasticStrain >= curve.back().first)
    {
        return curve.back().second;
    }

    // Perform linear interpolation
    for (size_t i = 0; i < curve.size() - 1; ++i)
    {
        if (plasticStrain >= curve[i].first && plasticStrain < curve[i + 1].first)
        {
            double x0 = curve[i].first, y0 = curve[i].second;
            double x1 = curve[i + 1].first, y1 = curve[i + 1].second;

            return y0 + (y1 - y0) * (plasticStrain - x0) / (x1 - x0);
        }
    }

    return 0.0; 
}


double computePlasticModulus(double plasticStrain, const std::vector<std::pair<double, double>>& curve)
{
    if (curve.size() < 2)
    {
        std::cerr << "Error: Hardening curve must have at least two points.\n";
        return 0.0;
    }

    if (plasticStrain <= curve.front().first)
    {
        double ep = (curve[1].second - curve[0].second) /
                    (curve[1].first - curve[0].first);
        return ep;
    }

    if (plasticStrain >= curve.back().first)
    {
        size_t last = curve.size() - 1;
        double ep = (curve[last].second - curve[last - 1].second) /
                    (curve[last].first - curve[last - 1].first);
        return ep;
    }

    for (size_t i = 0; i < curve.size() - 1; ++i)
    {
        double e1 = curve[i].first, s1 = curve[i].second;
        double e2 = curve[i + 1].first, s2 = curve[i + 1].second;

        if (plasticStrain >= e1 && plasticStrain < e2)
        {
            double ep = (s2 - s1) / (e2 - e1);
            return ep;
        }
    }

    return 0.0; 
}