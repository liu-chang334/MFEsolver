#include <iostream> 
#include "include/C3D8.h"
#include "include/FiniteElementModel.h"
#include "include/FiniteElementSolver.h"
#include "include/Tools.h"
#include "include/FEdataModelPost.h"
#include "include/Math.h"
#include "src/Spinner.cpp"

VTK_MODULE_INIT(vtkRenderingOpenGL2);
VTK_MODULE_INIT(vtkInteractionStyle);

int main()
{
    // read FEM data from inp
    std::string filefolder1 = std::filesystem::current_path().string(); 
    std::string filename1 = "dogbone2.inp";
    std::filesystem::path filepath1 = filefolder1 + "\\.." + "\\.." + "\\abaqus" + "\\" + filename1;
    FiniteElementModel feaModel;
    ABAQUSFEMReader(filepath1, feaModel);

    // initialize the FEA solver and solve
    FiniteElementSolver feaSolver(feaModel);
    feaSolver.solve_linearelastic();

    // bool is_success = feaSolver.performNewtonIteration(10, 1e-6, 1.0);
    // feaSolver.writeToFile();

    // post-process
    FEDataModelPost feaModelPost(feaModel);
    // feaModelPost.FEdataPlotScalar("U", 1);
    // feaModelPost.FEdataPlotScalar("U", 2);
    // feaModelPost.FEdataPlotScalar("U", 3);
    // feaModelPost.FEdataPlotScalar("E", 11);
    // feaModelPost.FEdataPlotScalar("E", 22);
    // feaModelPost.FEdataPlotScalar("E", 33);
    // feaModelPost.FEdataPlotScalar("E", 12);
    // feaModelPost.FEdataPlotScalar("E", 13);
    // feaModelPost.FEdataPlotScalar("E", 23);
    feaModelPost.FEdataPlotScalar("S", 11);
    feaModelPost.FEdataPlotScalar("S", 22);
    feaModelPost.FEdataPlotScalar("S", 33);
    feaModelPost.FEdataPlotScalar("S", 12);
    feaModelPost.FEdataPlotScalar("S", 13);
    feaModelPost.FEdataPlotScalar("S", 23);
    // feaModelPost.FEdataPlotScalar("S_princ", 1);
    // feaModelPost.FEdataPlotScalar("S_princ", 2);
    // feaModelPost.FEdataPlotScalar("S_princ", 3);
    // feaModelPost.FEdataPlotScalar("E_princ", 1);
    // feaModelPost.FEdataPlotScalar("E_princ", 2);
    // feaModelPost.FEdataPlotScalar("E_princ", 3);

    return 0;
}
