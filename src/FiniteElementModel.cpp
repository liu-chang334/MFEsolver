#include "../include/FiniteElementModel.h"
#include "../include/Tools.h"

/**
 * @brief Expands the Node matrix by one row and assigns the provided [x, y, z] coordinates.
 *
 * @param[in] node A vector of size 3 containing the node coordinates [x, y, z].
 */
void FiniteElementModel::addNode(const std::vector<double>& node) {
    Eigen::Index currentRows = Node.rows();
    Node.conservativeResize(currentRows + 1, 3);
    Node.row(currentRows) = Eigen::Vector3d(node[0], node[1], node[2]);
}

/**
 * @brief Expands the Element matrix by one row and assigns the provided element nodes.
 * 
 * @param[in] element A vector of integers representing the element nodes [node1, node2, node3, ... , nodeN].
 */
void FiniteElementModel::addElement(const std::vector<int>& element) {
    Eigen::Index currentRows = Element.rows();
    Element.conservativeResize(currentRows + 1, element.size());
    for (size_t i = 0; i < element.size(); ++i) {
        Element(currentRows, i) = element[i];
    }
}

void FiniteElementModel::addMaterial(const MaterialData& mat) {
    Materials.push_back(mat);
}

/**
 * @brief Expands the Nsets map by adding a new nset with the provided name and nodes.
 *
 * @param[in] name The name of the nset.
 * @param[in] nodes A vector of integers representing the node IDs in the nset [node1, node2, node3,..., nodeN].
 * @note The nset name must be unique. If a nset with the same name already exists,
 *      the new nodes will be overwritten using the new one.
 */
void FiniteElementModel::addNset(const std::string& name, const std::vector<int>& nodes) {
    Eigen::VectorXi nset(nodes.size());
    for (size_t i = 0; i < nodes.size(); ++i) {
        nset[i] = nodes[i];
    }
    auto it = Nsets.find(name);
    if (it != Nsets.end()){
        Eigen::VectorXi& existingNodes = it->second;
        size_t originalSize = existingNodes.size();
        size_t newSize = originalSize + nset.size();
        existingNodes.conservativeResize(newSize);
        for (int i = 0; i < nset.size(); ++i) {
            existingNodes[originalSize + i] = nset[i];
        }
    } else { Nsets[name] = nset;}
}

/**
 * @brief Expands the Esets map by adding a new eset with the provided name and elements.
 *
 * @param[in] name The name of the elset.
 * @param[in] elements A vector of integers representing the element IDs in the eset [ele1, ele2, ele3,..., eleN].
 * @note The eset name must be unique. If an eset with the same name already exists,
 *      the new elements will be overwritten using the new one.
 */
void FiniteElementModel::addElset(const std::string& name, const std::vector<int>& elements) {
    Eigen::VectorXi elset(elements.size());
    for (size_t i = 0; i < elements.size(); ++i) {
        elset[i] = elements[i];
    }
    auto it = Elsets.find(name);
    if (it!= Elsets.end()){
        Eigen::VectorXi& existingElements = it->second;
        size_t originalSize = existingElements.size();
        size_t newSize = originalSize + elset.size();
        existingElements.conservativeResize(newSize);
        for (int i = 0; i < elset.size(); ++i) {
            existingElements[originalSize + i] = elset[i];
        }
    } else { Elsets[name] = elset;}
}

/**
 * @brief Expands the Force matrix by adding rows for the specified nset and assigns the provided dofid and value.
 *
 * @param[in] name The name of the nset.
 * @param[in] dofid The degree of freedom ID.
 * @param[in] value The value to be assigned to the specified dofid.
 * @note The nset name must exist in the Nsets map. 
*/
void FiniteElementModel::addLoad(const std::string& name, const int& dofid, const double& value) {
    Eigen::VectorXi nodes = Nsets[name];
    Eigen::Index currentRows = Force.rows();
    Force.conservativeResize(currentRows + nodes.size(), 3);
    for (int i = 0; i < nodes.size(); ++i) {
        Force.row(currentRows + i) = Eigen::Vector3d(nodes[i], dofid, value);
    }
}

/**
 * @brief Expands the Constraint matrix by adding rows for the specified nset and assigns the provided dofid and value.
 *
 * @param[in] name The name of the nset.
 * @param[in] dofid The degree of freedom ID.
 * @param[in] value The value to be assigned to the specified dofid.
 * @note The nset name must exist in the Nsets map.
*/
void FiniteElementModel::addConstraint(const std::string& name, const int& dofid, const double& value) {
    Eigen::VectorXi nodes = Nsets[name];
    Eigen::Index currentRows = Constraint.rows();
    Constraint.conservativeResize(currentRows + nodes.size(), 3);
    for (int i = 0; i < nodes.size(); ++i) {
        Constraint.row(currentRows + i) = Eigen::Vector3d(nodes[i], dofid, value);
    }
}


/**
 * @brief Reads the FEM data from an inp file and populates the FiniteElementModel object.
 *
 * @param[in] filepath The path to the inp file.
 * @param[in,out] femModel Returns the populated FiniteElementModel object.
 * @note 
*/
void ABAQUSFEMReader(const std::filesystem::path& filepath, FiniteElementModel& femModel) 
{
    std::ifstream file(filepath);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file " << filepath << std::endl;
        return;
    }

    // Flags control data parsing: *Node, *Element, *Material, *Load, *Boundary, *Nset, *Elset
    bool inSection[7] = {false, false, false, false, false, false, false};
    std::string currentSetName;  // Store the current set name
    std::string line; // Store the current line

    MaterialData tempMat;
    bool readingElastic = false;
    bool readingPlastic = false;
    bool outsideMaterialBlock = false;

    while (std::getline(file, line)) {
        if (line.empty()) continue;// Skip empty lines
        if (line[0] == '*' && line[1] == '*') continue; // Skip comment lines
        if (line.rfind("*Step", 0) == 0) continue; // not implemented yet
        if (line.rfind("*Static", 0) == 0) continue; // not implemented yet

        if (line[0] == '*' && outsideMaterialBlock && !tempMat.name.empty()) {
            femModel.addMaterial(tempMat);
            tempMat = MaterialData();
            outsideMaterialBlock = false;
            readingElastic = false;
            readingPlastic = false;
        }

        // Check section markers
        if (line.rfind("*Node", 0) == 0) {
            std::fill(std::begin(inSection), std::end(inSection), false);
            inSection[0] = true;
        } else if (line.rfind("*Element", 0) == 0) {
            std::fill(std::begin(inSection), std::end(inSection), false);
            inSection[1] = true;
        } else if (line.rfind("*Material", 0) == 0) {
            std::fill(std::begin(inSection), std::end(inSection), false);
            inSection[2] = true;
            readingElastic = false;
            readingPlastic = false;
            outsideMaterialBlock = false;

            size_t namePos = findCaseInsensitive(line, "name=");
            if (namePos != std::string::npos) {
                tempMat.name = line.substr(namePos + 5);
                tempMat.name.erase(remove(tempMat.name.begin(), tempMat.name.end(), ','), tempMat.name.end()); 
            } else {
                std::cerr << "Error: Material name not specified in the line: " << line << std::endl;
            }
        } else if (line.rfind("*Elastic", 0) == 0) {
            readingElastic = true;
            readingPlastic = false; 
        } else if (line.rfind("*Plastic", 0) == 0) {
            readingElastic = false;
            readingPlastic = true; 
        } else if (line.rfind("*Load", 0) == 0) {
            std::fill(std::begin(inSection), std::end(inSection), false);
            inSection[3] = true;
            outsideMaterialBlock = true;
        } else if (line.rfind("*Boundary", 0) == 0) {
            std::fill(std::begin(inSection), std::end(inSection), false);
            inSection[4] = true;
            outsideMaterialBlock = true;
        } else if (line.rfind("*Nset", 0) == 0) {
            std::fill(std::begin(inSection), std::end(inSection), false);
            inSection[5] = true; 
            outsideMaterialBlock = true;
            // size_t namePos = line.find("Nset=");
            size_t namePos = findCaseInsensitive(line, "Nset=");
            if (namePos != std::string::npos) {
                currentSetName = line.substr(namePos + 5);
                // currentSetName.erase(remove(currentSetName.begin(), currentSetName.end(), ','), currentSetName.end());
            }else if (namePos == std::string::npos) {
                std::cerr << "Error in inp file! ";
                std::cerr << "Please check the line: " << line << std::endl;
            }
        } else if (line.rfind("*Elset", 0) == 0) {
            std::fill(std::begin(inSection), std::end(inSection), false);
            inSection[6] = true; 
            outsideMaterialBlock = true;
            // size_t namePos = line.find("Elset=");
            size_t namePos = findCaseInsensitive(line, "Elset=");
            if (namePos!= std::string::npos) {
                currentSetName = line.substr(namePos + 5);
                // currentSetName.erase(remove(currentSetName.begin(), currentSetName.end(), ','), currentSetName.end());
            }else if (namePos == std::string::npos) {
                std::cerr << "Error in inp file! ";
                std::cerr << "Please check the line: " << line << std::endl;
            }
        }
        // Parse data   
        std::istringstream iss(line);
        char comma; // store the comma ","
        if (inSection[0]) { // parse node coordinates data such as: 1, 0.0, 0.0, 0.0
            std::vector<double> node(3);
            int id;
            if (iss >> id >> comma >> node[0] >> comma >> node[1] >> comma >> node[2]) {
                femModel.addNode(node);
            }
        } else if (inSection[1]) { // parse element data such as: 1, 1, 2, 3, 4, 5, 6, 7, 8
            std::vector<int> nodeIDsofelement;
            int id;
            iss >> id >> comma; // skip the element ID and then data such as: 1, 2, 3, 4, 5, 6, 7, 8
            int node_id;
            int count = 0;
            while (count < 8 && iss >> node_id) {
                nodeIDsofelement.push_back(node_id);
                iss >> comma;
                count++;
            }
            if (count == 8) { 
                femModel.addElement(nodeIDsofelement);
            } 
        } else if (inSection[2] && readingElastic) { // parse material data such as: 2.0e11, 0.3
            iss >> tempMat.E >> comma >> tempMat.nu;
            tempMat.type = "LinearElastic";

        } else if (inSection[2] && readingPlastic) {
            double stress, epstrain;
            if (iss >> stress >> comma >> epstrain) {
               tempMat.HardeningCurve.emplace_back(stress, epstrain); 
            }
            tempMat.type = "ElasticPlastic";
        } else if (inSection[3]) { // parse load data such as: Nset-1, 1, 1.0e6
            std::string nodeSet;
            int dofid;
            double value;
            if (iss >> nodeSet) {
                // delete the last character if it's a comma
                if (nodeSet.back() == ',') {
                    nodeSet.pop_back();
                }
                // read dof and value
                if (iss >> dofid >> comma >> value) {
                    femModel.addLoad(nodeSet, dofid, value);  
                }  
            }
        }else if (inSection[4]) { // parse constraint data such as: Nset-1, 1, 0.0
            std::string nodeSet;    
            int dofid;
            double value;
            if (iss >> nodeSet) {
                // delete the last character if it's a comma
                if (nodeSet.back() == ',') {
                    nodeSet.pop_back();
                }
                // read dof and value
                if (iss >> dofid >> comma >> value) {
                    femModel.addConstraint(nodeSet, dofid, value);
                }
            }
        } else if (inSection[5]) { // parse node set data such as: 1, 2, 3, 4, 5, 6, 7, 8, 9, 10
            std::vector<int> nodes;
            int nodeId;
            while (iss >> nodeId) {
                nodes.push_back(nodeId);
                if (iss.peek() == ',' || iss.peek() == ' ') iss.ignore(); // skip comma or space
            }
            if (!nodes.empty() && !currentSetName.empty()) {
                femModel.addNset(currentSetName, nodes); 
            }
        } else if (inSection[6]) { // parse element set data such as: 1, 2, 3, 4, 5, 6, 7, 8, 9, 10
            std::vector<int> elements;
            int elementId;
            while (iss >> elementId) {
                elements.push_back(elementId);
                if (iss.peek() == ',' || iss.peek() ==' ') iss.ignore(); // skip comma or space
            } 
            if (!elements.empty() && !currentSetName.empty()) {
                femModel.addElset(currentSetName, elements); 
            }
        }
    }

    if (outsideMaterialBlock && !tempMat.name.empty()) {
        femModel.addMaterial(tempMat);
    }

    file.close();
}

/**
 * @brief Get the nodes ID of the element.
 *
 * @param[in] elementID The ID of the element.
 * @param[in,out] nodeIDs Return the nodes ID of the element.
 * @note
*/
void FiniteElementModel::getNodesIDofElement(int elementID, Eigen::VectorXi& nodeIDs) const
{
    if (elementID <= 0 || elementID > Element.rows()) {
        std::cerr << "Error: elementID is out of range!" << std::endl;
        std::exit(EXIT_FAILURE);
    }
    const Eigen::VectorXi& elemnode = Element.row(elementID-1);
    nodeIDs = elemnode;
}

void FiniteElementModel::printMaterialInfo() const {
    std::cout << "======== Material Information ========" << std::endl;
    for (const auto& mat : Materials) {
        std::cout << "Material Name: " << mat.name << std::endl;
        std::cout << "Material Type: " << mat.type << std::endl;
        std::cout << "  Young's Modulus (E): " << mat.E << std::endl;
        std::cout << "  Poisson's Ratio (nu): " << mat.nu << std::endl;

        if (!mat.HardeningCurve.empty()) {
            std::cout << "  Plastic Stress-Strain Curve:" << std::endl;
            std::cout << "    sigma       epsilon" << std::endl;
            for (const auto& [sigma, epsilon] : mat.HardeningCurve) {
                std::cout << "    " << sigma << "        " << epsilon << std::endl;
            }
        } else {
            std::cout << "  (Linear Elastic Material)" << std::endl;
        }

        std::cout << "-------------------------------------" << std::endl;
    }
}
