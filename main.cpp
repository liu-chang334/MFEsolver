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
    std::filesystem::path filefolder1 = "D:/liuchang/FEsolvercxx/abaqus";
    std::string filename1 = "dogbone.inp";
    std::filesystem::path filepath1 = filefolder1 / filename1;
    FiniteElementModel feaModel;
    ABAQUSFEMReader(filepath1, feaModel);

    // initialize the FEA solver
    FiniteElementSolver feaSolver(feaModel);

    // solve the FEA problem
    feaSolver.solve();

    // save the results
    std::string filefolder2 = "D:/liuchang/FEsolvercxx/FEoutput";
    saveMatrix2TXT(feaSolver.K, filefolder2, (std::string)"K.txt", 5);
    // saveMatrix2TXT(feaSolver.F, filefolder2, (std::string)"F.txt", 15);
    saveMatrix2TXT(feaSolver.U, filefolder2, (std::string)"U.txt", 15);

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
