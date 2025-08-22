#pragma once
#include <Eigen/Dense>
#include <vector>

// one definition by inline
inline const double g = 1.0 / sqrt(3.0);
inline std::vector<Eigen::Vector3d> gp_rst = {
    {-g, -g, -g},
    { g, -g, -g},
    { g,  g, -g},
    {-g,  g, -g},
    {-g, -g,  g},
    { g, -g,  g},
    { g,  g,  g},
    {-g,  g,  g}
};

class MaterialPoint3D{
public:
    int id;
    
    Eigen::Vector3d rst;
    Eigen::Vector3d xyz;

    Eigen::MatrixXd B;
    Eigen::MatrixXd B_bar;

    Eigen::VectorXd strain;
    Eigen::VectorXd stress;
public:
    MaterialPoint3D(){};
    MaterialPoint3D(int id) : id(id){
        rst = gp_rst[id];

        strain = Eigen::VectorXd::Zero(6);
        stress = Eigen::VectorXd::Zero(6);
        B = Eigen::MatrixXd::Zero(6, 24);
        B_bar = Eigen::MatrixXd::Zero(6, 24);
    };
};
