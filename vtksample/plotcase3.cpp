#include <vtkNew.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCellType.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkAutoInit.h>
#include <vtkProperty.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkNamedColors.h>
#include <vtkDoubleArray.h>
#include <vtkArrayCalculator.h>
#include <vtkScalarBarActor.h>

// 初始化 VTK 模块（仅在静态编译时需要）
VTK_MODULE_INIT(vtkRenderingOpenGL2);
VTK_MODULE_INIT(vtkInteractionStyle);

int main()
{
    /* 1. vtkPoints类: 用于存储几何点的位置
    */
    vtkNew<vtkPoints> points;
    points->InsertNextPoint(1.0, 1.0, 1.0);
    points->InsertNextPoint(1.0, 0.0, 1.0);
    points->InsertNextPoint(1.0, 1.0, 0.0);
    points->InsertNextPoint(1.0, 0.0, 0.0);
    points->InsertNextPoint(0.0, 1.0, 1.0);
    points->InsertNextPoint(0.0, 0.0, 1.0);
    points->InsertNextPoint(0.0, 1.0, 0.0);
    points->InsertNextPoint(0.0, 0.0, 0.0);

    /* 2. vtkDoubleArray类: 用于存储标量或矢量分量数据
    */
    vtkNew<vtkDoubleArray> displacement;
    displacement->SetName("Displacement");
    displacement->SetNumberOfComponents(3);
    displacement->InsertNextTuple3(0.1, 1.5, 0.3);
    displacement->InsertNextTuple3(0.4, 0.5, 0.6);
    displacement->InsertNextTuple3(0.7, 0.8, 0.9);
    displacement->InsertNextTuple3(1.0, 1.1, 1.2);
    displacement->InsertNextTuple3(1.3, 1.4, 1.5);
    displacement->InsertNextTuple3(1.6, 1.5, 1.8);
    displacement->InsertNextTuple3(2.0, 2.1, 2.2);
    displacement->InsertNextTuple3(2.3, 1.4, 2.5);
    

    /* 3. vtkUnstructuredGrid类: 表示非结构化网格数据
    */
    vtkNew<vtkUnstructuredGrid> grid;
    grid->SetPoints(points); 

    vtkIdType pts[8] = {4, 5, 7, 6, 0, 1, 3, 2};
    grid->InsertNextCell(VTK_HEXAHEDRON, 8, pts);
    grid->GetPointData()->SetVectors(displacement); // 设置点的位移(矢量数据)

    /* 4. vtkcalculator类: 用于计算标量或矢量数据
    */
    vtkNew<vtkArrayCalculator> calculator;
    calculator->SetInputData(grid);
    calculator->SetAttributeTypeToPointData();
    calculator->AddVectorArrayName("Displacement");
    calculator->SetFunction("Displacement[1]");
    calculator->SetResultArrayName("DisplacementX");
    calculator->Update();

    double range[2];
    vtkDataSet* dataSet = vtkDataSet::SafeDownCast(calculator->GetOutput());
    vtkDataArray* displacement_ = dataSet->GetPointData()->GetArray("DisplacementX");
    displacement_->GetRange(range);

    /* 5. vtkDataSetMapper类: 将数据集映射为几何图形
    */
    vtkNew<vtkDataSetMapper> mapper;
    mapper->SetInputConnection(calculator->GetOutputPort());
    mapper->SetScalarModeToUsePointFieldData();
    mapper->SelectColorArray("DisplacementX");
    mapper->SetScalarRange(range);

    /* 6. vtkNamedColors类: 用于获取 VTK 中预定义的颜色, 可以通过名称获取颜色值。
    GetColor3d(const vtkStdString &name): 获取指定名称的颜色的 RGB 值。
    GetData(): 获取颜色数组的指针。
    SetColor(const vtkStdString &name, double r, double g, double b): 设置指定名称的颜色的 RGB 值。
    ColorExists(const vtkStdString &name): 检查指定名称的颜色是否存在。
    */
    vtkNew<vtkNamedColors> colors;
    double *peacock;
    peacock = colors->GetColor3d("Peacock").GetData();
    std::cout << "Peacock:" << std::endl;
    std::cout << peacock[0] << " " << peacock[1] << " " << peacock[2] << std::endl;// 0.2 0.631373 0.788235

    colors->SetColor("MyColor", 0.2, 0.631373, 0.788235);

    bool isExist = colors->ColorExists("MyColor");
    if (isExist)
    {
        std::cout << "MyColor exists." << std::endl;
    }
    else
    {
        std::cout << "MyColor does not exist." << std::endl;
    }

    /* 7. vtkActor类: 表示可渲染的图形对象
    */
    vtkNew<vtkActor> actor;
    actor->SetMapper(mapper);
    actor->GetProperty()->SetColor(colors->GetColor3d("MyColor").GetData());
    actor->GetProperty()->SetEdgeVisibility(true);
    actor->GetProperty()->SetEdgeColor(1.0, 0.0, 0.0);
    actor->GetProperty()->SetLineWidth(2.0);

    /* 8. vtkRenderer类: 用于渲染图形
    */
    vtkNew<vtkRenderer> renderer;
    renderer->AddActor(actor); 
    renderer->SetBackground(0.5, 0.5, 0.5); 

    /* 9. vtkRenderWindow类: 用于显示渲染结果
    */
    vtkNew<vtkRenderWindow> renderWindow;
    renderWindow->AddRenderer(renderer); 
    renderWindow->SetSize(640, 480); 
    renderWindow->Render();

    /* 10. vtkScalarBarActor类: 用于显示标量条
    */
    vtkNew<vtkScalarBarActor> scalarBar;
    scalarBar->SetLookupTable(mapper->GetLookupTable());
    scalarBar->SetTitle("DisplacementX");
    scalarBar->SetNumberOfLabels(5);
    renderer->AddActor(scalarBar);


    /* 11. vtkInteractorStyleTrackballCamera类: 轨迹球交互样式
       鼠标操作可以实现视角旋转、平移和缩放。
    */
    vtkNew<vtkInteractorStyleTrackballCamera> style;

    /* 12. vtkRenderWindowInteractor类: 处理用户交互
    */
    vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
    renderWindowInteractor->SetRenderWindow(renderWindow); 
    renderWindowInteractor->SetInteractorStyle(style); 
    renderWindowInteractor->Initialize(); 
    renderWindowInteractor->Start(); 

    return 0;
}
