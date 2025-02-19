#include "../include/ElemBase.h"

// Class SolidElement

/**
 * @brief Construct a new SolidElement:: SolidElement object
 *
 * @param elemID Element ID
 * @param matID Material ID
 */
SolidElement::SolidElement(int elemID, int matID) : elementID(elemID), materialID(matID) {}

/**
 * @brief Set the Nodes object
 *
 * @param ids Node IDs
 * @param coords Node coordinates
 */
void SolidElement::setNodes(const Eigen::VectorXi& ids, const Eigen::MatrixXd& coords) {
    nodeIDs = ids;
    nodeCoorMatrix = coords;
}

/**
 * @brief Set the Material By Lame object
 *
 * @param lambda Lame's first parameter
 * @param mu Lame's second parameter
 */
void SolidElement::setMaterialByLame(double lambda, double mu) 
{
    this->lambda = lambda;
    this->mu = mu;
    this->E = mu * (3 * lambda + 2 * mu) / (lambda + mu);
    this->nu = lambda / (2 * (lambda + mu));
    this->G = mu;
}

/**
 * @brief Set the Material By YoungPoisson object
 *
 * @param E Young's Modulus
 * @param nu Poisson's Ratio
 */
void SolidElement::setMaterialByYoungPoisson(double E, double nu) 
{
    this->E = E;
    this->nu = nu;
    this->lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
    this->mu = E / (2 * (1 + nu));
    this->G = E / (2 * (1 + nu));
}

/**
 * @brief Set the Material By YoungShear object
 *
 * @param E Young's Modulus
 * @param G Shear Modulus
 */
void SolidElement::setMaterialByYoungShear(double E, double G) 
{
    this->E = E;
    this->G = G;
    this->nu = E / (2 * G) - 1;
    this->lambda = 2 * G * this->nu / (1 - 2 * this->nu);
    this->mu = G;
}

/**
 * @brief Set the Density object
 *
 * @param density Density
 */
void SolidElement::setDensity(double density) 
{
    this->rho = density;
}

// Class SolidElement End