#include "../include/Tools.h"

/**
 * @brief Save a matrix to a TXT file with adjustable precision.
 *
 * @param matrix The matrix to be saved.
 * @param filefolder The folder where the file will be saved.
 * @param filename The name of the file to be saved.
 * @param precision The number of decimal places to be saved.
 */
void saveMatrix2TXT(const Eigen::MatrixXd& matrix, const std::string& filefolder, const std::string& filename, int precision)
{
    // check if the filefolder exists, if not, create it
    if (!std::filesystem::exists(filefolder)) {
        std::filesystem::create_directory(filefolder);
    }
    // check if the file exists, if yes, delete it
    std::string filepath = filefolder + '\\' + filename;
    if (std::filesystem::exists(filepath)) {
        std::filesystem::remove(filepath);
    }

    std::ofstream file(filepath);
    if (file.is_open())
    {
        // Set the precision as requested
        file << std::fixed << std::setprecision(precision);
        
        for (int i = 0; i < matrix.rows(); ++i)
        {
            for (int j = 0; j < matrix.cols(); ++j)
            {
                file << std::setw(10 + precision)  // Adjust width to account for increased precision
                     << matrix(i, j);

                if (j < matrix.cols() - 1) {
                    file << " ";
                }
            }
            file << "\n";
        }
        file.close();
        // std::cout << "Matrix saved to " << filepath << " with " << precision << " decimal places." << std::endl;
    }
    else
    {
        std::cerr << "Error: cannot open file " << filepath << std::endl;
    }
}

void saveMatrix2TXT(const Eigen::MatrixXi& matrix, const std::string& filefolder, const std::string& filename)
{
    // check if the filefolder exists, if not, create it
    if (!std::filesystem::exists(filefolder)) {
        std::filesystem::create_directory(filefolder);
    } 
    // check if the file exists, if yes, delete it
    std::string filepath = filefolder + '/' + filename;
    if (std::filesystem::exists(filepath)) {
        std::filesystem::remove(filepath);
    }

    std::ofstream file(filepath);
    if (file.is_open())
    {
        for (int i = 0; i < matrix.rows(); ++i)
        {
            for (int j = 0; j < matrix.cols(); ++j)
            {
                file << matrix(i, j);

                if (j < matrix.cols() - 1) {
                    file << " ";
                }
            } 
            file << "\n";
        } 
    }
    else
    {
        std::cerr << "Error: cannot open file " << filepath << std::endl;
    }
    file.close();
    std::cout << "Matrix saved to " << filepath << std::endl;
}

/**
 * @brief Load a matrix from a TXT file.
 *
 * @param filepath The path of the file to be loaded.
 * @return Eigen::MatrixXd The loaded matrix.
 */
Eigen::MatrixXd loadMatrixFromTXT(const std::string& filepath)
{
    std::ifstream file(filepath);
    if(!file.is_open())
    {
        std::cerr << "Error: cannot open file " << filepath << std::endl;
    }

    std::vector<std::vector<double>> data;
    std::string line;
    size_t cols = 0;
    while(std::getline(file, line))
    {
        std::istringstream lineStream(line);
        std::vector<double> row;
        double value;
        while(lineStream >> value)
        {
            row.push_back(value);
        }
        if(!row.empty())
        {
            if(cols == 0)
            {
                cols = row.size();
            }else if(row.size() != cols)
            {
                std::cerr << "Error: inconsistent number of columns in file " << filepath << std::endl;
            }
            data.push_back(row);
        }
    }
    file.close();

    // Convert the vector of vectors to an Eigen matrix
    Eigen::MatrixXd matrix(data.size(), cols);
    for(size_t i = 0; i < data.size(); ++i)
    {
        for(size_t j = 0; j < cols; ++j)
        {
            matrix(i, j) = data[i][j];
        }
    }
    return matrix;
}

/**
 * @brief Find the index of a substring in a string, ignoring case.
 *
 * @param haystack The string to search in.
 * @param needle The substring to search for.
 */
size_t findCaseInsensitive(const std::string& haystack, const std::string& needle) {
    auto it = std::search(
        haystack.begin(), haystack.end(),
        needle.begin(), needle.end(),
        [](char ch1, char ch2) { return std::tolower(ch1) == std::tolower(ch2); }
    );

    if (it != haystack.end()) {
        return std::distance(haystack.begin(), it);
    } else {
        return std::string::npos;
    }
}