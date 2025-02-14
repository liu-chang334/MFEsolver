#include <vtkNew.h>
#include <vtkPoints.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCellType.h>

int main()
{
    vtkNew<vtkPoints> points;
    points->InsertNextPoint(1.0, 1.0, 1.0);
    points->InsertNextPoint(1.0, 0.0, 1.0);
    points->InsertNextPoint(1.0, 1.0, 0.0);
    points->InsertNextPoint(1.0, 0.0, 0.0);
    points->InsertNextPoint(0.0, 1.0, 1.0);
    points->InsertNextPoint(0.0, 0.0, 1.0);
    points->InsertNextPoint(0.0, 1.0, 0.0);
    points->InsertNextPoint(0.0, 0.0, 0.0);
    vtkNew<vtkUnstructuredGrid> grid;

    // 1. SetPoints(vtkPoints *): 向vtkUnstructuredGrid中添加点
    grid->SetPoints(points);
    
    // 2. vtkIdType InsertNextCell(int cellType, vtkIdType numPoints, const vtkIdType* pointIds):
    // 向vtkUnstructuredGrid中添加单元
    // cellType: 单元类型的枚举值，如VTK_TRIANGLE, VTK_QUAD, VTK_HEXAHEDRON等
    // numPoints: 单元中点的数量
    // pointIds: 单元中点的索引，从0开始计数  
    vtkIdType pts[8] = {4, 5, 7, 6, 0, 1, 3, 2};
    grid->InsertNextCell(VTK_HEXAHEDRON, 8, pts);
    
    // 3. GetNumberOfPoints(): 获取vtkUnstructuredGrid中点的数量
    int numPoints = grid->GetNumberOfPoints();
    std::cout << "Number of points in grid: " << numPoints << std::endl;

    // 4. GetNumberOfCells(): 获取vtkUnstructuredGrid中单元的数量
    int numCells = grid->GetNumberOfCells();
    std::cout << "Number of cells in grid: " << numCells << std::endl;

    // 5. GetCellType(vtkIdType cellId): 获取vtkUnstructuredGrid中指定单元的类型
    vtkIdType cellType;
    cellType = grid->GetCellType(0);
    const char* cellName = vtkCellTypes::GetClassNameFromTypeId(cellType);
    std::cout << "Cell type of cell 0: " << cellName << std::endl;

    return 0;
}