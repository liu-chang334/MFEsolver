#include "IsotropicHarden.h"

void IsotropicHarden::setMaterialParameters(const std::unordered_map<std::string, double>& parameters)
{
    materialParams_ = parameters;

    if (materialParams_.find("E") == materialParams_.end() ||
        materialParams_.find("nu") == materialParams_.end())
    {
        std::cerr << "Material parameters not found" << std::endl;
        return; 
    }

    E_ = materialParams_["E"];
    nu_ = materialParams_["nu"];

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
    double yield_stress = interpolateHardeningStress(mp.epstrain, hardeningCurve_);
    double Ep = computePlasticModulus(mp.epstrain, hardeningCurve_);
    double k = yield_stress * yield_stress / 3.0;
    double m; // 0 <= m < 1

    if (trial_f - k <= 0)
    {
        m = 1.0;
        mp.strain += dstrain;
        mp.stress = trial_stress;
        return;

    } else if (trial_f - k > 0 && f - k < 0)
    {
        Eigen::VectorXd dstress_dev = computeDeviatoric(dstress);
        Eigen::VectorXd stress_dev = computeDeviatoric(mp.stress);
        double a0, a1, a2;
        a0 = f - k;
        a1 = stress_dev.dot(dstress_dev);
        a2 = dstress_dev.dot(dstress_dev);
        m = (-1 * a1 + sqrt(a1 * a1 - 4 * a0 * a2)) / (2 * a2);

    } else if (trial_f - k > 0 && std::abs(f - k) < 1e-9){
        m = 0.0;
    }

    // 3. compute the elastic increment of stress and elastic plastic increment of strain
    Eigen::VectorXd dstrain_ep = (1 - m) * dstrain;
    Eigen::VectorXd stress_ep_start = mp.stress + m * dstress;
    Eigen::VectorXd stress_ep_start_dev = computeDeviatoric(stress_ep_start);

    Eigen::VectorXd df_dsigma = Eigen::VectorXd::Zero(6);
    df_dsigma << stress_ep_start_dev(0), stress_ep_start_dev(1), stress_ep_start_dev(2),
        2 * stress_ep_start_dev(3), 2 * stress_ep_start_dev(4), 2 * stress_ep_start_dev(5);

    Dp_ = (D_ * df_dsigma.transpose() * df_dsigma * D_) / 
        ((df_dsigma * D_ * df_dsigma.transpose())(0,0) + yield_stress * yield_stress * Ep * 4.0 / 9.0);

    Eigen::MatrixXd Dep = D_ - Dp_;

    // 4. predict along the tangent direction
    Eigen::VectorXd dstress_ep = Dep * dstrain_ep;

    double dlambda = (df_dsigma * D_ * dstrain.transpose())(0,0) /
        ((df_dsigma * D_ * df_dsigma.transpose())(0,0) + yield_stress * yield_stress * Ep * 4.0 / 9.0);

    double depstrain = 2 / 3 * dlambda * yield_stress;

    Eigen::VectorXd stress_pre = mp.stress + m * dstress + dstress_ep;
    double epstrain = mp.epstrain + depstrain;

    // 5. return along the radial direction
    Eigen::VectorXd stress_pre_dev = computeDeviatoric(stress_pre);
    double yield_stress_end = interpolateHardeningStress(epstrain, hardeningCurve_);
    double r = sqrt((2 / 3) * yield_stress_end * yield_stress_end / (stress_pre_dev.dot(stress_pre_dev)));

    // 6. update the material point
    mp.strain += dstrain;
    mp.stress = r * stress_pre;
    mp.epstrain = epstrain;

    // 7. update the tangent matrix
    double Ep_end = computePlasticModulus(epstrain, hardeningCurve_);
    Eigen::VectorXd stress_end_dev = computeDeviatoric(mp.stress);
    Eigen::VectorXd df_dsigma_end = Eigen::VectorXd::Zero(6);
    df_dsigma_end << stress_end_dev(0), stress_end_dev(1), stress_end_dev(2),
        2 * stress_end_dev(3), 2 * stress_end_dev(4), 2 * stress_end_dev(5);
    Eigen::MatrixXd Dp_end = (D_ * df_dsigma_end.transpose() * df_dsigma_end * D_) /
            ((df_dsigma_end * D_ * df_dsigma_end.transpose())(0,0) + yield_stress_end * yield_stress_end * Ep_end * 4.0 / 9.0);

    tangent = D_ - Dp_end;
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