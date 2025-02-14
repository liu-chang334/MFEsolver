#include <vtkNew.h>
#include <vtkPoints.h>
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

    /* 2. vtkUnstructuredGrid类: 表示非结构化网格数据
    */
    vtkNew<vtkUnstructuredGrid> grid;
    grid->SetPoints(points); 
    vtkIdType pts[8] = {4, 5, 7, 6, 0, 1, 3, 2};
    grid->InsertNextCell(VTK_HEXAHEDRON, 8, pts);

    /* 3. vtkDataSetMapper类: 将数据集映射为几何图形
    */
    vtkNew<vtkDataSetMapper> mapper;
    mapper->SetInputData(grid);

    /* 4. vtkNamedColors类: 用于获取 VTK 中预定义的颜色, 可以通过名称获取颜色值。
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

    /* 5. vtkActor类: 表示可渲染的图形对象
    */
    vtkNew<vtkActor> actor;
    actor->SetMapper(mapper);
    actor->GetProperty()->SetColor(colors->GetColor3d("MyColor").GetData());
    actor->GetProperty()->SetEdgeVisibility(true);
    actor->GetProperty()->SetEdgeColor(1.0, 0.0, 0.0);
    actor->GetProperty()->SetLineWidth(2.0);

    /* 6. vtkRenderer类: 用于渲染图形
    */
    vtkNew<vtkRenderer> renderer;
    renderer->AddActor(actor); 
    renderer->SetBackground(0.5, 0.5, 0.5); 

    /* 7. vtkRenderWindow类: 用于显示渲染结果
    */
    vtkNew<vtkRenderWindow> renderWindow;
    renderWindow->AddRenderer(renderer); 
    renderWindow->SetSize(640, 480); 
    renderWindow->Render(); 

    /* 8. vtkInteractorStyleTrackballCamera类: 轨迹球交互样式
       鼠标操作可以实现视角旋转、平移和缩放。
    */
    vtkNew<vtkInteractorStyleTrackballCamera> style;

    /* 9. vtkRenderWindowInteractor类: 处理用户交互
    */
    vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
    renderWindowInteractor->SetRenderWindow(renderWindow); 
    renderWindowInteractor->SetInteractorStyle(style); 
    renderWindowInteractor->Initialize(); 
    renderWindowInteractor->Start(); 

    return 0;
}
