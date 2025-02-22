#ifndef FINITEELEMENTMODEL_H
#define FINITEELEMENTMODEL_H

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <fstream>
#include <filesystem>
#include <string>
#include <sstream>
#include <cstdlib>

/**
 * @class FiniteElementModel
 * @brief Represents a finite element model with associated data structures.
 *
 * The following matrices and maps store different pieces of information:
 *  - \c Node:
 *   [x1, y1, z1
 *    x2, y2, z2
 *    ...
 *    xn, yn, zn]
 *
 * - \c Element:
 *   [ele1node1, ele1node2, ..., ele1node8
 *    ele2node1, ele2node2, ..., ele2node8
 *    ...
 *    eleNnode1, eleNnode2, ..., eleNnode8]
 *
 * - \c Material:
 *   [E, nu]
 *
 * - \c Force:
 *   [nodeid, dofid, value
 *    ...      ...   ...]
 *   
 * - \c Constraint:
 *   [nodeid, dofid, value
 *    ...      ...   ...]
 *  
 * - \c Nsets: A mapping from a set name (e.g., "nsetname1") to the corresponding
 *   collection of node indices.
 *   [nsetname1, [nodeid1, nodeid2, ...]
 *    nsetname2, [nodeid1, nodeid2, ...]
 *   ...]
 *
 * - \c Esets**: A mapping from a set name (e.g., "esetname1") to the corresponding
 *   collection of element indices.
 *   [esetname1, [eleid1, eleid2,...]
 *    esetname2, [eleid1, eleid2,...]
 *  ...]
 *
 * @note The class provides methods to add nodes, elements, materials, loads,
 * constraints, node sets, and element sets to the model.
 */
class FiniteElementModel {
public:
    Eigen::MatrixXd Node;
    Eigen::MatrixXi Element;
    Eigen::MatrixXd Material;
    Eigen::MatrixXd Force;
    Eigen::MatrixXd Constraint;
    std::unordered_map<std::string, Eigen::VectorXi> Nsets;
    std::unordered_map<std::string, Eigen::VectorXi> Elsets;

public:
    void addNode(const std::vector<double>& node);
    void addElement(const std::vector<int>& element);
    void addMaterial(const std::vector<double>& material);
    void addLoad(const std::string& name, const int& dofid, const double& value);
    void addConstraint(const std::string& name, const int& dofid, const double& value);
    void addNset(const std::string& nsetName, const std::vector<int>& nodeIDs);
    void addElset(const std::string& esetName, const std::vector<int>& elementIDs);
    void printModelInfo();
    void getNodesIDofElement(int elementID, Eigen::VectorXi &nodeIDs);
};


void ABAQUSFEMReader(const std::filesystem::path& filepath, FiniteElementModel& femModel);

/**
 * @brief Prints information about a given element.
 *
 * This function takes an element as input and prints various properties of the element.
 * The element type is determined by the template parameter ElementType.
 */
template <typename ElementType>
void printElementInfo(ElementType &element){
    std::cout << "Element Information: " << std::endl;
    std::cout << "Element ID: " << element.elementID << std::endl;
    std::cout << "Node IDs: ";
    std::cout << element.nodeIDs << std::endl;
    std::cout << std::endl;
    std::cout << "Node Coordinates:" << std::endl;
    std::cout << element.nodeCoorMatrix << std::endl;

    std::cout << "Material ID: " << element.materialID << std::endl;
    std::cout << "Material Properties: " << std::endl;
    std::cout << "Lame's first parameter: " << element.lambda << std::endl;
    std::cout << "Lame's second parameter: " << element.mu << std::endl;
    std::cout << "Density: " << element.rho << std::endl;
    std::cout << "Young's Modulus: " << element.E << std::endl;
    std::cout << "Poisson's Ratio: " << element.nu << std::endl;
    std::cout << "Shear Modulus: " << element.G << std::endl;


    std::cout << "Internal Gauss Points and Weights:" << std::endl;
    if constexpr (std::is_same<ElementType, C3D8>::value) {
        // Print Gauss points and weights
        std::cout << "Gauss Points:" << std::endl;
        std::cout << element.gaussPoints << std::endl;
        std::cout << "Gauss Weights:" << std::endl;
        std::cout << element.gaussWeights.transpose() << std::endl;
    }
}

#endif // FINITEELEMENTMODEL_H