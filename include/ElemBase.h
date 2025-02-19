#ifndef ELEMBASE_H
#define ELEMBASE_H

#include <vector>
#include <Eigen/Dense>
#include <cmath> 

/**
 * @class SolidElement
 * @brief This class represents a solid element in a finite element analysis.
 * - \c elementID: Element ID
 * - \c nodeIDs: Node IDs
 * - \c nodeCoorMatrix: Node coordinates
 * - \c materialID: Material ID
 * - \c lambda: Lame's first parameter
 * - \c mu: Lame's second parameter
 * - \c rho: density
 * - \c E: Young's modulus
 * - \c nu: Poisson's ratio
 * - \c G: Shear modulus
 * @details
 * - The material properties can be set by Lame's parameters, Young's modulus and Poisson's ratio, or Young's modulus and shear modulus.
 * @note This class is an base class for other element types.
 * @see C3D8
 */
class SolidElement 
{
public:
    int elementID;         
    Eigen::VectorXi nodeIDs;  
    Eigen::MatrixXd nodeCoorMatrix; 

    int materialID;        
    double lambda;         // Lame's first parameter
    double mu;             // Lame's second parameter
    double rho;           
    double E;           
    double nu;             
    double G;              

public:
    SolidElement(int elemID, int matID = -1);
    
    void setNodes(const Eigen::VectorXi& ids, const Eigen::MatrixXd& coords);
    void setMaterialByLame(double lambda, double mu);
    void setMaterialByYoungPoisson(double E, double nu);
    void setMaterialByYoungShear(double E, double G);
    void setDensity(double density);
};

#endif // ELEMEBASE_H