#include "../include/ElemBase.h"

//************************SolidElement Base Class
SolidElement::SolidElement(int elemID, int matID) : elementID(elemID), materialID(matID) {}

void SolidElement::setNodes(const Eigen::VectorXi& ids, const Eigen::MatrixXd& coords) {
    nodeIDs = ids;
    nodeCoorMatrix = coords;
}

void SolidElement::setMaterialByLame(double lambda, double mu) 
{
    this->lambda = lambda;
    this->mu = mu;
    this->E = mu * (3 * lambda + 2 * mu) / (lambda + mu);
    this->nu = lambda / (2 * (lambda + mu));
    this->G = mu;
}

void SolidElement::setMaterialByYoungPoisson(double E, double nu) 
{
    this->E = E;
    this->nu = nu;
    this->lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
    this->mu = E / (2 * (1 + nu));
    this->G = E / (2 * (1 + nu));
}

void SolidElement::setMaterialByYoungShear(double E, double G) 
{
    this->E = E;
    this->G = G;
    this->nu = E / (2 * G) - 1;
    this->lambda = 2 * G * this->nu / (1 - 2 * this->nu);
    this->mu = G;
}

void SolidElement::setDensity(double density) 
{
    this->rho = density;
}

//*************************SolidElement Base Class End