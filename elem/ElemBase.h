// #pragma once
#include <string>

#include <Eigen/Dense>

class ElemBase{
public:
    int eid;
    std::string type;
    Eigen::VectorXi node;
    Eigen::MatrixXi node_coor;
public:
    ElemBase(int eid, std::string type);
    ~ElemBase() = default;
};