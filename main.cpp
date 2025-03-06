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
    feaSolver.solve();
    
    // post-process
    FEDataModelPost feaModelPost(feaModel);
    // feaModelPost.FEdataPlotScalar("U", 1);
    // feaModelPost.FEdataPlotScalar("S", 11);
    // feaModelPost.FEdataPlotScalar("E", 11);
    feaModelPost.FEdataPlotScalar("S_princ", 1);
    // feaModelPost.FEdataPlotScalar("E_princ", 1);

    return 0;
}
