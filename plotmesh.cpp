#include <vtkUnstructuredGrid.h>
#include <vtkDataSetMapper.h>
#include <vtkCellType.h>
#include <vtkProperty.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkAutoInit.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkCylinderSource.h>
#include <vtkPolyDataMapper.h>

VTK_MODULE_INIT(vtkRenderingOpenGL2);
VTK_MODULE_INIT(vtkInteractionStyle);

int main()
{
    
    static float x[27][3] = {
        {0,0,0}, {1,0,0}, {2,0,0}, {0,1,0}, {1,1,0}, {2,1,0},
        {0,0,1}, {1,0,1}, {2,0,1}, {0,1,1}, {1,1,1}, {2,1,1},
        {0,1,2}, {1,1,2}, {2,1,2}, {0,1,3}, {1,1,3}, {2,1,3},
        {0,1,4}, {1,1,4}, {2,1,4}, {0,1,5}, {1,1,5}, {2,1,5},
        {0,1,6}, {1,1,6}, {2,1,6}
    };
    static vtkIdType pts[12][8] = {
        {0,1,4,3,6,7,10,9},
        {1,2,5,4,7,8,11,10},
        {6,10,9,12,0,0,0,0},
        {8,11,10,14,0,0,0,0},
        {15,16,17,14,13,12,0,0},
        {18,15,19,16,20,17,0,0},
        {22,23,20,19,0,0,0,0},
        {21,22,18,0,0,0,0,0},
        {22,19,18,0,0,0,0,0},
        {23,26,0,0,0,0,0,0},
        {21,24,0,0,0,0,0,0},
        {25,0,0,0,0,0,0,0}
    };

    vtkNew<vtkPoints> points;
    for (int i = 0; i < 27; i++) points->InsertPoint(i, x[i]);

    vtkNew<vtkUnstructuredGrid> ugrid;
    ugrid->SetPoints(points);

    ugrid->InsertNextCell(VTK_HEXAHEDRON, 8, pts[0]);
    ugrid->InsertNextCell(VTK_HEXAHEDRON, 8, pts[1]);
    ugrid->InsertNextCell(VTK_TETRA, 4, pts[2]);
    ugrid->InsertNextCell(VTK_TETRA, 4, pts[3]);
    ugrid->InsertNextCell(VTK_POLYGON, 6, pts[4]);
    ugrid->InsertNextCell(VTK_TRIANGLE_STRIP, 6, pts[5]);
    ugrid->InsertNextCell(VTK_QUAD, 4, pts[6]);
    ugrid->InsertNextCell(VTK_TRIANGLE, 3, pts[7]);
    ugrid->InsertNextCell(VTK_TRIANGLE, 3, pts[8]);
    ugrid->InsertNextCell(VTK_LINE, 2, pts[9]);
    ugrid->InsertNextCell(VTK_LINE, 2, pts[10]);
    ugrid->InsertNextCell(VTK_VERTEX, 1, pts[11]);

    vtkNew<vtkDataSetMapper> ugridMapper;
    ugridMapper->SetInputData(ugrid);
    vtkNew<vtkActor> ugridActor;
    ugridActor->SetMapper(ugridMapper);
    ugridActor->GetProperty()->SetRepresentationToWireframe();
    ugridActor->GetProperty()->SetColor(0, 0, 0);

    vtkNew<vtkRenderer> renderer;
    renderer->AddActor(ugridActor);
    renderer->SetBackground(1, 1, 1);
    vtkNew<vtkRenderWindow> renWin;
    renWin->AddRenderer(renderer);
    renWin->SetSize(450, 450);
    renWin->Render();

    vtkNew<vtkRenderWindowInteractor> iren;
    iren->SetRenderWindow(renWin);
    vtkNew<vtkInteractorStyleTrackballCamera> style;
    iren->SetInteractorStyle(style);
    iren->Initialize();
    iren->Start();

    return 0;
}
