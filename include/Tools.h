#ifndef TOOLS_H
#define TOOLS_H

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <fstream>
#include <filesystem>
#include <string>
#include <iomanip>
#include <sstream>
#include <string>


void saveMatrix2TXT(const Eigen::MatrixXd& matrix, const std::string& filefolder, 
                    const std::string& filename, int precision = 6);
Eigen::MatrixXd loadMatrixFromTXT(const std::string& filepath);


#endif // TOOLS_H