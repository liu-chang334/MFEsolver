#ifndef ELEMBASE_H
#define ELEMBASE_H

#include <vector>
#include <Eigen/Dense>
#include <cmath> 


class SolidElement 
{
public:
    int elementID;         // ElemID
    Eigen::VectorXi nodeIDs;   // NodeID
    Eigen::MatrixXd nodeCoorMatrix; // Nodexyz


    int materialID;        // MatID
    double lambda;         // Lame's first parameter
    double mu;             // Lame's second parameter
    double rho;            // density
    double E;              // Young's modulus
    double nu;             // Poisson's ratio
    double G;              // Shear modulus

public:
    SolidElement(int elemID, int matID = -1);

    // set nodes number and coordinates
    void setNodes(const Eigen::VectorXi& ids, const Eigen::MatrixXd& coords);

    // set material properties by Lame's parameters
    void setMaterialByLame(double lambda, double mu);

    // set material properties by Young's modulus and Poisson's ratio
    void setMaterialByYoungPoisson(double E, double nu);

    // set material properties by Young's modulus and shear modulus
    void setMaterialByYoungShear(double E, double G);

    // set material density
    void setDensity(double density);
};

#endif // ELEMEBASE_H