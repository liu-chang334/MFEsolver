#pragma once
#include <Eigen/Dense>

#include <cons/MaterialBase.h>
#include <elem/MaterialPoint3D.h>

class IsoElastic : public MaterialBase {
public:
    IsoElastic(double E, double nu, double rho, int dim) : MaterialBase(E, nu, rho, dim) {}

    void updatemp(Eigen::VectorXd dstrain, MaterialPoint3D& mp) override {

        Eigen::VectorXd strain = mp.strain;
        Eigen::VectorXd stress = mp.stress;
        
        Eigen::VectorXd dstress = De * dstrain;
        stress += dstress;

        mp.stress = stress;
        mp.strain = strain + dstrain;

    }
};