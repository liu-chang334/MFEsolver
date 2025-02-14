#ifndef FINITEELEMENTMODEL_H
#define FINITEELEMENTMODEL_H

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <fstream>
#include <filesystem>
#include <string>
#include <sstream>

// FEA Model
class FiniteElementModel {
public:
    Eigen::MatrixXd Node;
    Eigen::MatrixXi Element;
    Eigen::MatrixXd Material;
    Eigen::MatrixXd Force;
    Eigen::MatrixXd Constr;
    std::unordered_map<std::string, Eigen::VectorXi> Nsets;
    std::unordered_map<std::string, Eigen::VectorXi> Esets;

public:
    void addNode(const std::vector<double>& node);
    void addElement(const std::vector<int>& element);
    void addMaterial(const std::vector<double>& material);
    void addLoad(const std::vector<double>& load);
    void addConstraint(const std::vector<double>& constraint);
    void addNset(const std::string& nsetName, const std::vector<int>& nodeIDs);
    void addEset(const std::string& esetName, const std::vector<int>& elementIDs);

    void printModelInfo();

    Eigen::VectorXi getNodesIDofElement(int elementID);
};

// ABAQUS INP reader
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