#include <cons/MaterialBase.h>

MaterialBase::MaterialBase(double E, double nu, double rho, int dim) : E_(E), nu_(nu), rho_(rho), dim_(dim) {
    if (dim_ == 3) {
        De = Eigen::MatrixXd::Zero(6, 6);
        double lambda = E_ * nu_ / ((1 + nu_) * (1 - 2 * nu_));
        double mu = E_ / (2 * (1 + nu_));

        // elastic 
        De << lambda + 2 * mu, lambda, lambda, 0, 0, 0,
            lambda, lambda + 2 * mu, lambda, 0, 0, 0,
            lambda, lambda, lambda + 2 * mu, 0, 0, 0,
            0, 0, 0, mu, 0, 0,
            0, 0, 0, 0, mu, 0,
            0, 0, 0, 0, 0, mu;
    }
}

void MaterialBase::setMaterialProperty(double E, double nu, double rho, int dim) {
    this->E_ = E;
    this->nu_ = nu;
    this->rho_ = rho;
    this->dim_ = dim;
}
