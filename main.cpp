#include <iostream> 
#include "include/C3D8.h"
#include "include/FiniteElementModel.h"
#include "include/FiniteElementSolver.h"
#include "include/Tools.h"
#include "include/FEdataModelPost.h"

VTK_MODULE_INIT(vtkRenderingOpenGL2);
VTK_MODULE_INIT(vtkInteractionStyle);

int main()
{
    // read FEM data from inp
    // get the path of the project folder
    std::string filefolder1 = std::filesystem::current_path().string(); 
    std::string filename1 = "dogbone2.inp";
    std::filesystem::path filepath1 = filefolder1 + "\\.." + "\\.." + "\\abaqus" + "\\" + filename1;
    FiniteElementModel feaModel;
    ABAQUSFEMReader(filepath1, feaModel);

    // initialize the FEA solver
    FiniteElementSolver feaSolver(feaModel);

    // std::cout << "Set-2:" << "\n" << feaModel.Nsets["Set-2"] << std::endl;

    // solve the FEA problem
    feaSolver.solve();

    // post-process
    FEDataModelPost feaModelPost(feaModel);
    feaModelPost.FEdataPlotScalar("U", 1);
    // feaModelPost.FEdataPlotScalar("U", 2);
    // feaModelPost.FEdataPlotScalar("U", 3);

    // calculate element strain and stress
    Eigen::MatrixXd elementStrain = feaSolver.calcuElementStrain(50, false);
    std::cout << elementStrain << std::endl;
    Eigen::MatrixXd elementStress = feaSolver.calcuElementStress(50, false);
    std::cout << elementStress << std::endl;

    return 0;
}
