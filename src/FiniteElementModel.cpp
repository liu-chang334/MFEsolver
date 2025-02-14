#include "../include/FiniteElementModel.h"

/* Todo: add corrdinate of nodes in order to member variable --> Node
   node: [x, y, z]
*/
void FiniteElementModel::addNode(const std::vector<double>& node) {
    Eigen::Index currentRows = Node.rows();
    Node.conservativeResize(currentRows + 1, 3);
    Node.row(currentRows) = Eigen::Vector3d(node[0], node[1], node[2]);
}

/* Todo: add nodeID of element in order to member variable --> Element
   element: [nodeid1, nodeid2,..., nodeid8]
*/
void FiniteElementModel::addElement(const std::vector<int>& element) {
    Eigen::Index currentRows = Element.rows();
    Element.conservativeResize(currentRows + 1, element.size());
    for (size_t i = 0; i < element.size(); ++i) {
        Element(currentRows, i) = element[i];
    }
}

/* Todo: add material parameter to member variable --> Material
   material: [E, nu]
*/
void FiniteElementModel::addMaterial(const std::vector<double>& material) {
    Eigen::Index currentRows = Material.rows();
    Material.conservativeResize(currentRows + 1, 2);
    Material.row(currentRows) = Eigen::Vector2d(material[0], material[1]);
}

/* Todo: add nodeID of node set to member variable --> Nsets
   name: name of node set
   nodes: [nodeid1, nodeid2,..., nodeidn]
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

/* Todo: add elementID of element set to member variable --> Esets
   name: name of element set
   elements: [elementid1, elementid2,..., elementidn]
*/
void FiniteElementModel::addEset(const std::string& name, const std::vector<int>& elements) {
    Eigen::VectorXi eset(elements.size());
    for (size_t i = 0; i < elements.size(); ++i) {
        eset[i] = elements[i];
    }
    auto it = Esets.find(name);
    if (it!= Esets.end()){
        Eigen::VectorXi& existingElements = it->second;
        size_t originalSize = existingElements.size();
        size_t newSize = originalSize + eset.size();
        existingElements.conservativeResize(newSize);
        for (int i = 0; i < eset.size(); ++i) {
            existingElements[originalSize + i] = eset[i];
        }
    } else { Esets[name] = eset;}
}

/* Todo: add load information to member variable --> Force
   name: nset name
   dofid: dof id
   value: load value
*/
void FiniteElementModel::addLoad(const std::string& name, const int& dofid, const double& value) {
    Eigen::VectorXi nodes = Nsets[name];
    Eigen::Index currentRows = Force.rows();
    Force.conservativeResize(currentRows + nodes.size(), 3);
    for (int i = 0; i < nodes.size(); ++i) {
        Force.row(currentRows + i) = Eigen::Vector3d(nodes[i], dofid, value);
    }
}

/* Todo: add constraint information to member variable --> Constraint
   name: nset name
   dofid: dof id
   value: constraint value
*/
void FiniteElementModel::addConstraint(const std::string& name, const int& dofid, const double& value) {
    Eigen::VectorXi nodes = Nsets[name];
    Eigen::Index currentRows = Constraint.rows();
    Constraint.conservativeResize(currentRows + nodes.size(), 3);
    for (int i = 0; i < nodes.size(); ++i) {
        Constraint.row(currentRows + i) = Eigen::Vector3d(nodes[i], dofid, value);
    }
}


/* Todo: read inp file and store data to object ---> FiniteElementModel
   filepath: path to inp file such as "D:/fem/beam.inp"
   femModel: FiniteElementModel object to store data, need to be initialized before calling this function
*/
void ABAQUSFEMReader(const std::filesystem::path& filepath, FiniteElementModel& femModel) 
{
    std::ifstream file(filepath);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file " << filepath << std::endl;
        return;
    }

    // Flags control data parsing: *Node, *Element, *Material, *Load, *Constr, *Nset, *Eset
    bool inSection[7] = {false, false, false, false, false, false, false};
    std::string currentSetName;  // Store the current set name
    std::string line; // Store the current line

    while (std::getline(file, line)) // Read each line
    {
        if (line.empty()) continue;// Skip empty lines

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
        } else if (line.rfind("*Load", 0) == 0) {
            std::fill(std::begin(inSection), std::end(inSection), false);
            inSection[3] = true;
        } else if (line.rfind("*Constr", 0) == 0) {
            std::fill(std::begin(inSection), std::end(inSection), false);
            inSection[4] = true;
        } else if (line.rfind("*Nset", 0) == 0) {
            std::fill(std::begin(inSection), std::end(inSection), false);
            inSection[5] = true; 
            size_t namePos = line.find("Nset=");
            if (namePos != std::string::npos) {
                currentSetName = line.substr(namePos + 5);
                currentSetName.erase(remove(currentSetName.begin(), currentSetName.end(), ','), currentSetName.end());
            }else if (namePos == std::string::npos) {
                std::cerr << "Error in inp file! ";
                std::cerr << "Please check the line: " << line << std::endl;
            }
        } else if (line.rfind("*Eset", 0) == 0) {
            std::fill(std::begin(inSection), std::end(inSection), false);
            inSection[6] = true; 
            size_t namePos = line.find("Eset=");
            if (namePos!= std::string::npos) {
                currentSetName = line.substr(namePos + 5);
                currentSetName.erase(remove(currentSetName.begin(), currentSetName.end(), ','), currentSetName.end());
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
        } else if (inSection[2]) { // parse material data such as: 2.0e11, 0.3
            std::vector<double> material(2);
            if (iss >> material[0] >> comma >> material[1]) {
                femModel.addMaterial(material);
            }
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
                femModel.addEset(currentSetName, elements); 
            }
        }
    }
    file.close();
}

/* Todo: get nodes id of element
   elementID: element id
   return: nodeIDs, nodes id of element such as [nodeid1, nodeid2,..., nodeid8]
*/
void FiniteElementModel::getNodesIDofElement(int elementID, Eigen::VectorXi& nodeIDs) {
    if (elementID <= 0 || elementID > Element.rows()) {
        std::cerr << "Error: elementID is out of range!" << std::endl;
        std::exit(EXIT_FAILURE);
    }
    Eigen::VectorXi elemnode = Element.row(elementID-1);
    nodeIDs = elemnode;
}