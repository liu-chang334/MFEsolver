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

    // 1. SetPoints(vtkPoints *): ��vtkUnstructuredGrid����ӵ�
    grid->SetPoints(points);
    
    // 2. vtkIdType InsertNextCell(int cellType, vtkIdType numPoints, const vtkIdType* pointIds):
    // ��vtkUnstructuredGrid����ӵ�Ԫ
    // cellType: ��Ԫ���͵�ö��ֵ����VTK_TRIANGLE, VTK_QUAD, VTK_HEXAHEDRON��
    // numPoints: ��Ԫ�е������
    // pointIds: ��Ԫ�е����������0��ʼ����  
    vtkIdType pts[8] = {4, 5, 7, 6, 0, 1, 3, 2};
    grid->InsertNextCell(VTK_HEXAHEDRON, 8, pts);
    
    // 3. GetNumberOfPoints(): ��ȡvtkUnstructuredGrid�е������
    int numPoints = grid->GetNumberOfPoints();
    std::cout << "Number of points in grid: " << numPoints << std::endl;

    // 4. GetNumberOfCells(): ��ȡvtkUnstructuredGrid�е�Ԫ������
    int numCells = grid->GetNumberOfCells();
    std::cout << "Number of cells in grid: " << numCells << std::endl;

    // 5. GetCellType(vtkIdType cellId): ��ȡvtkUnstructuredGrid��ָ����Ԫ������
    vtkIdType cellType;
    cellType = grid->GetCellType(0);
    const char* cellName = vtkCellTypes::GetClassNameFromTypeId(cellType);
    std::cout << "Cell type of cell 0: " << cellName << std::endl;

    return 0;
}