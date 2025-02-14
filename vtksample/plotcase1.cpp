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

// 初始化 VTK 模块（仅在静态编译时需要）
VTK_MODULE_INIT(vtkRenderingOpenGL2);
VTK_MODULE_INIT(vtkInteractionStyle);

int main()
{
    /* 1. vtkPoints类: 用于存储几何点的位置
       InsertNextPoint(double x, double y, double z): 添加一个点到点集中，点的坐标为(x, y, z)。
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
       SetPoints(vtkPoints *): 设置网格中的点数据。
       InsertNextCell(int cellType, vtkIdType numPts, vtkIdType *pts): 添加一个单元到网格中，指定单元类型和顶点索引。
    */
    vtkNew<vtkUnstructuredGrid> grid;
    grid->SetPoints(points); 
    vtkIdType pts[8] = {4, 5, 7, 6, 0, 1, 3, 2};
    grid->InsertNextCell(VTK_HEXAHEDRON, 8, pts);

    /* 3. vtkDataSetMapper类: 将数据集映射为几何图形
       SetInputData(vtkDataSet *input): 设置输入数据，用于可视化。
    */
    vtkNew<vtkDataSetMapper> mapper;
    mapper->SetInputData(grid);

    /* 4. vtkActor类: 表示可渲染的图形对象
       SetMapper(vtkMapper *): 设置映射器，将几何数据与图形对象关联。
       GetProperty(): 获取图形对象的属性。
       SetColor(double r, double g, double b): 设置图形对象的颜色。
       SetEdgeVisibility(bool visibility): 设置图形对象的边可见性。
       SetEdgeColor(double r, double g, double b): 设置图形对象的边颜色。
       SetLineWidth(double width): 设置图形对象的线宽。
    */
    vtkNew<vtkActor> actor;
    actor->SetMapper(mapper);
    actor->GetProperty()->SetColor(0.0, 0.0, 1.0);
    actor->GetProperty()->SetEdgeVisibility(true);
    actor->GetProperty()->SetEdgeColor(1.0, 0.0, 0.0);
    actor->GetProperty()->SetLineWidth(2.0);

    /* 5. vtkRenderer类: 用于渲染图形
       AddActor(vtkProp *): 添加图形对象到渲染器。
       SetBackground(double r, double g, double b): 设置背景颜色，范围为0.0~1.0。
    */
    vtkNew<vtkRenderer> renderer;
    renderer->AddActor(actor); 
    renderer->SetBackground(0.5, 0.5, 0.5); 

    /* 6. vtkRenderWindow类: 用于显示渲染结果
       AddRenderer(vtkRenderer *): 添加渲染器到窗口。
       SetSize(int width, int height): 设置窗口的宽和高。
       Render(): 执行渲染。
    */
    vtkNew<vtkRenderWindow> renderWindow;
    renderWindow->AddRenderer(renderer); 
    renderWindow->SetSize(640, 480); 
    renderWindow->Render(); 

    /* 7. vtkInteractorStyleTrackballCamera类: 轨迹球交互样式
       鼠标操作可以实现视角旋转、平移和缩放。
    */
    vtkNew<vtkInteractorStyleTrackballCamera> style;

    /* 8. vtkRenderWindowInteractor类: 处理用户交互
       SetRenderWindow(vtkRenderWindow *): 设置关联的渲染窗口。
       SetInteractorStyle(vtkInteractorObserver *): 设置交互样式。
       Initialize(): 初始化交互器。
       Start(): 开始用户交互。
    */
    vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
    renderWindowInteractor->SetRenderWindow(renderWindow); 
    renderWindowInteractor->SetInteractorStyle(style); 
    renderWindowInteractor->Initialize(); 
    renderWindowInteractor->Start(); 

    return 0;
}
