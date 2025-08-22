#pragma once
#include <Eigen/Dense>

#include <elem/MaterialPoint3D.h>

class MaterialBase {
public:
    double E_;
    double nu_;
    double rho_;
    int dim_;
    Eigen::MatrixXd De;

public:
    MaterialBase(double E, double nu, double rho, int dim);

    void setMaterialProperty(double E, double nu, double rho, int dim);

    virtual void updatemp(Eigen::VectorXd dstrain, MaterialPoint3D& mp) = 0;




};
