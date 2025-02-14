#include "../include/FiniteElementModel.h"


void FiniteElementModel::addNode(const std::vector<double>& node) {
    Eigen::Index currentRows = Node.rows();
    Node.conservativeResize(currentRows + 1, 3);
    Node.row(currentRows) = Eigen::Vector3d(node[0], node[1], node[2]);
}

void FiniteElementModel::addElement(const std::vector<int>& element) {
    Eigen::Index currentRows = Element.rows();
    Element.conservativeResize(currentRows + 1, element.size());
    for (size_t i = 0; i < element.size(); ++i) {
        Element(currentRows, i) = element[i];
    }
}

void FiniteElementModel::addMaterial(const std::vector<double>& material) {
    Eigen::Index currentRows = Material.rows();
    Material.conservativeResize(currentRows + 1, 2);
    Material.row(currentRows) = Eigen::Vector2d(material[0], material[1]);
}

void FiniteElementModel::addNset(const std::string& name, const std::vector<int>& nodes) {
    Eigen::VectorXi nset(nodes.size());
    for (size_t i = 0; i < nodes.size(); ++i) {
        nset[i] = nodes[i];
    }
    Nsets[name] = nset;
}

void FiniteElementModel::addEset(const std::string& name, const std::vector<int>& elements) {
    Eigen::VectorXi eset(elements.size());
    for (size_t i = 0; i < elements.size(); ++i) {
        eset[i] = elements[i];
    }
    Esets[name] = eset; 
}

void FiniteElementModel::addLoad(const std::vector<double>& load) {
    Eigen::Index currentRows = Force.rows();
    Force.conservativeResize(currentRows + 1, 3);
    Force.row(currentRows) = Eigen::Vector3d(load[0], load[1], load[2]);
}

void FiniteElementModel::addConstraint(const std::vector<double>& constraint) {
    Eigen::Index currentRows = Constr.rows();
    Constr.conservativeResize(currentRows + 1, 3);
    Constr.row(currentRows) = Eigen::Vector3d(constraint[0], constraint[1], constraint[2]);
}



void ABAQUSFEMReader(const std::filesystem::path& filepath, FiniteElementModel& femModel) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file " << filepath << std::endl;
        return;
    }

    std::string line;
    // Flags control data parsing: *Node, *Element, *Material, *Load, *Constr, *Nset, *Eset
    bool inSection[7] = {false, false, false, false, false, false, false};
    std::string currentSetName;

    while (std::getline(file, line)) {
        if (line.empty()) continue;

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
            }
        } else if (line.rfind("*Eset", 0) == 0) {
            std::fill(std::begin(inSection), std::end(inSection), false);
            inSection[6] = true; 
            size_t namePos = line.find("Eset=");
            if (namePos!= std::string::npos) {
                currentSetName = line.substr(namePos + 5);
                currentSetName.erase(remove(currentSetName.begin(), currentSetName.end(), ','), currentSetName.end());
            }
        }

        // Parse data   
        std::istringstream iss(line);
        char comma;
        if (inSection[0]) {
            std::vector<double> node(3);
            int id;
            if (iss >> id >> comma >> node[0] >> comma >> node[1] >> comma >> node[2]) {
                femModel.addNode(node);
            }
        } else if (inSection[1]) {
            std::vector<int> element;
            int id;
            iss >> id >> comma; // Skip element ID
            int node_id;
            int count = 0;

            while (count < 8 && iss >> node_id) {
                element.push_back(node_id);
                iss >> comma;
                count++;
            }

            if (count == 8) { 
                femModel.addElement(element);
            } 
        } else if (inSection[2]) {
            std::vector<double> material(2);
            if (iss >> material[0] >> comma >> material[1]) {
                femModel.addMaterial(material);
            }
        } else if (inSection[3]) {
            std::vector<double> load(3);
            int nodeId, dof;
            if (iss >> nodeId >> comma >> dof >> comma >> load[2]) {
                load[0] = nodeId;
                load[1] = dof;
                femModel.addLoad(load);
            }
        } else if (inSection[4]) {
            std::vector<double> constraint(3);
            int nodeId, dof;
            if (iss >> nodeId >> comma >> dof >> comma >> constraint[2]) {
                constraint[0] = nodeId;
                constraint[1] = dof;
                femModel.addConstraint(constraint);
            }
        } else if (inSection[5]) {
            std::vector<int> nodes;
            int nodeId;
            while (iss >> nodeId) {
                nodes.push_back(nodeId);
                if (iss.peek() == ',' || iss.peek() == ' ') iss.ignore();
            }
            if (!nodes.empty() && !currentSetName.empty()) {
                femModel.addNset(currentSetName, nodes); 
            }
        } else if (inSection[6]) {
            std::vector<int> elements;
            int elementId;
            while (iss >> elementId) {
                elements.push_back(elementId);
                if (iss.peek() == ',' || iss.peek() ==' ') iss.ignore();
            } 
            if (!elements.empty() && !currentSetName.empty()) {
                femModel.addEset(currentSetName, elements); 
            }
        }
    }
    file.close();
}

Eigen::VectorXi FiniteElementModel::getNodesIDofElement(int elementID) {
    Eigen::VectorXi elemnode = Element.row(elementID-1);
    return elemnode; 
}