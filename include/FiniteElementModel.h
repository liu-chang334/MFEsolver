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

// FEA Model
class FiniteElementModel {
public:

    /* Node: [x1, y1, z1
              x2, y2, z2
              ...,...,...
              xn, yn, zn]

    Element: [ele1node1, ele1node2, ..., ele1node8
               ele2node1, ele2node2, ..., ele2node8
               ...,...,...
               eleNnode1, eleNnode2, ..., eleNnode8]

    Material: [E, nu]

    Force: [nodeid, dofid, value
            ...,...,...]
    Constraint: [nodeid, dofid, value
               ...,...,...]
    
    Nsets: {nsetname1: [nodeid1, nodeid2, ...]
             nsetname2: [nodeid1, nodeid2, ...]}

    Esets: {esetname1: [eleid1, eleid2,...]
             esetname2: [eleid1, eleid2,...]}
    */
    Eigen::MatrixXd Node;
    Eigen::MatrixXi Element;
    Eigen::MatrixXd Material;
    Eigen::MatrixXd Force;
    Eigen::MatrixXd Constraint;
    std::unordered_map<std::string, Eigen::VectorXi> Nsets;
    std::unordered_map<std::string, Eigen::VectorXi> Esets;

public:
   /* Todo: add corrdinate of nodes in order to member variable --> Node
      node: [x, y, z]
   */
    void addNode(const std::vector<double>& node);

    /* Todo: add nodeID of element in order to member variable --> Element
       element: [nodeid1, nodeid2,..., nodeid8]
    */
    void addElement(const std::vector<int>& element);

    /* Todo: add material parameter to member variable --> Material
       material: [E, nu]
    */
    void addMaterial(const std::vector<double>& material);

    /* Todo: add load information to member variable --> Force
       name: nset name
       load: [dofid, value]
    */
    void addLoad(const std::string& name, const int& dofid, const double& value);

    /* Todo: add constraint information to member variable --> Constraint
       name: nset name
       constraint: [dofid, value]
    */
    void addConstraint(const std::string& name, const int& dofid, const double& value);

    /* Todo: add nodeID of node set to member variable --> Nsets
       name: name of node set
       nodes: [nodeid1, nodeid2,..., nodeidn]
    */
    void addNset(const std::string& nsetName, const std::vector<int>& nodeIDs);

    /* Todo: add elementID of element set to member variable --> Esets
       name: name of element set
       elements: [eleid1, eleid2,..., eleidn]
    */
    void addEset(const std::string& esetName, const std::vector<int>& elementIDs);

    void printModelInfo();

    /* Todo: get nodes id of element
       elementID: element id
       return: nodes id of element such as [nodeid1, nodeid2,..., nodeid8]
    */
    void getNodesIDofElement(int elementID, Eigen::VectorXi &nodeIDs);
};

// ABAQUS INP reader
/*Todo: read FEM data from inp file
    filepath: inp file path such as "F:/FEsolvercxx/abaqus/dogbone.inp"
    femModel: FEM model
*/
void ABAQUSFEMReader(const std::filesystem::path& filepath, FiniteElementModel& femModel);

// Print element information
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