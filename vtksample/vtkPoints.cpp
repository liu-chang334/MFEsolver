#include <vtkSmartPointer.h>
#include <vtkPoints.h>


int main()
{
    // 1. InsertNextPoint(double x, double y, double z): 向点集添加一个点。
    vtkNew<vtkPoints> points;
    points->InsertNextPoint(1.0, 2.0, 3.0); 
    points->InsertNextPoint(4.0, 5.0, 6.0); 

    // 2. GetPoint(int id, double point[3]): 获取点集中指定 ID 的点坐标。
    double point[3];
    points->GetPoint(1, point); 
    std::cout << "Point 1: " << point[0] << ", " << point[1] << ", " << point[2] << std::endl;

    // 3. SetPoint(int id, double x, double y, double z): 修改点集中指定 ID 的点坐标。
    points->SetPoint(1, 7.0, 8.0, 9.0);  
    points->GetPoint(1, point);
    std::cout << "Modified Point 1: " << point[0] << ", " << point[1] << ", " << point[2] << std::endl;

    // 4. GetNumberOfPoints(): 获取点集中的点数。
    int numPoints = points->GetNumberOfPoints();
    std::cout << "Number of points: " << numPoints << std::endl;

    // 5. InsertPoint(int id, double x, double y, double z): 在点集中指定 ID 的位置插入一个点, 会替换原来的点。
    points->InsertPoint(2, 10.0, 11.0, 12.0); 
    points->GetPoint(1, point);
    std::cout << "Inserted Point 1: " << point[0] << ", " << point[1] << ", " << point[2] << std::endl;
    points->GetNumberOfPoints();
    std::cout << "Number of points: " << numPoints << std::endl;

    // 6. GetBounds(double *bounds): 获取点集中的边界框。
    double bounds[6];
    points->GetBounds(bounds);
    std::cout << "Bounds:" << std::endl;
    std::cout << "Xmin: " << bounds[0] << ", Xmax: " << bounds[1] << std::endl;
    std::cout << "Ymin: " << bounds[2] << ", Ymax: " << bounds[3] << std::endl;
    std::cout << "Zmin: " << bounds[4] << ", Zmax: " << bounds[5] << std::endl;

    // 7. Initialize(): 清空点集。
    points->Initialize();
    numPoints = points->GetNumberOfPoints();
    std::cout << "Number of points after initialization: " << numPoints << std::endl;

    return 0;
}