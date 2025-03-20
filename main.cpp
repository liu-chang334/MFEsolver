#include <iostream> 
#include "include/C3D8.h"
#include "include/FiniteElementModel.h"
#include "include/FiniteElementSolver.h"
#include "include/FEdataModelPost.h"
#include "include/Math.h"


VTK_MODULE_INIT(vtkRenderingOpenGL2);
VTK_MODULE_INIT(vtkInteractionStyle);

int main()
{
    // // read FEM data from inp
    // std::string filefolder1 = std::filesystem::current_path().string(); 
    // std::string filename1 = "dogbone2.inp";
    // std::filesystem::path filepath1 = filefolder1 + "\\.." + "\\.." + "\\abaqus" + "\\" + filename1;
    // FiniteElementModel feaModel;
    // ABAQUSFEMReader(filepath1, feaModel);

    // // initialize the FEA solver and solve
    // FiniteElementSolver feaSolver(feaModel);
    // // feaSolver.solve_linearelastic();
    // double step_size = 0.50;  // (1 / step_size) must be an integer
    // int maxIter = 10;  
    // feaSolver.solve_adaptive_nonlinear(step_size, maxIter);

    // // post-process
    // FEDataModelPost feaModelPost(feaModel);
    // feaModelPost.FEdataPlotScalar("S", 11);
    // feaModelPost.FEdataPlotScalar("S_princ", 1);

    std::vector<double> vector = {1, 2, 3, 4, 5, 6};
    std::vector<double> result = computeEigenvalues3x3(vector);
    std::cout << "Eigenvalues: " << result[0] << " " << result[1] << " " << result[2] << std::endl;

    double trace = computerTrace3x3(vector);
    std::cout << "Trace: " << trace << std::endl;

    double secondInvariant = computerSecondInvariant3x3(vector);
    std::cout << "Second Invariant: " << secondInvariant << std::endl;

    double thirdInvariant = computerThirdInvariant3x3(vector);
    std::cout << "Third Invariant: " << thirdInvariant << std::endl;

    return 0;
}
