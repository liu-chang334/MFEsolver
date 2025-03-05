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
void saveMatrix2TXT(std::vector<Eigen::MatrixXd>& matrix, const std::string& filefolder, 
                    const std::string& filename, int precision = 6);
void saveMatrix2TXT(std::vector<Eigen::VectorXd>& matrix, const std::string& filefolder, 
                        const std::string& filename, int precision = 6);
void saveMatrix2TXT(const Eigen::MatrixXi& matrix, const std::string& filefolder, 
                    const std::string& filename);
Eigen::MatrixXd loadMatrixFromTXT(const std::string& filepath);
size_t findCaseInsensitive(const std::string& haystack, const std::string& needle);

#endif // TOOLS_H