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

// ��ʼ�� VTK ģ�飨���ھ�̬����ʱ��Ҫ��
VTK_MODULE_INIT(vtkRenderingOpenGL2);
VTK_MODULE_INIT(vtkInteractionStyle);

int main()
{
    /* 1. vtkPoints��: ���ڴ洢���ε��λ��
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

    /* 2. vtkUnstructuredGrid��: ��ʾ�ǽṹ����������
    */
    vtkNew<vtkUnstructuredGrid> grid;
    grid->SetPoints(points); 
    vtkIdType pts[8] = {4, 5, 7, 6, 0, 1, 3, 2};
    grid->InsertNextCell(VTK_HEXAHEDRON, 8, pts);

    /* 3. vtkDataSetMapper��: �����ݼ�ӳ��Ϊ����ͼ��
    */
    vtkNew<vtkDataSetMapper> mapper;
    mapper->SetInputData(grid);

    /* 4. vtkNamedColors��: ���ڻ�ȡ VTK ��Ԥ�������ɫ, ����ͨ�����ƻ�ȡ��ɫֵ��
    GetColor3d(const vtkStdString &name): ��ȡָ�����Ƶ���ɫ�� RGB ֵ��
    GetData(): ��ȡ��ɫ�����ָ�롣
    SetColor(const vtkStdString &name, double r, double g, double b): ����ָ�����Ƶ���ɫ�� RGB ֵ��
    ColorExists(const vtkStdString &name): ���ָ�����Ƶ���ɫ�Ƿ���ڡ�
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

    /* 5. vtkActor��: ��ʾ����Ⱦ��ͼ�ζ���
    */
    vtkNew<vtkActor> actor;
    actor->SetMapper(mapper);
    actor->GetProperty()->SetColor(colors->GetColor3d("MyColor").GetData());
    actor->GetProperty()->SetEdgeVisibility(true);
    actor->GetProperty()->SetEdgeColor(1.0, 0.0, 0.0);
    actor->GetProperty()->SetLineWidth(2.0);

    /* 6. vtkRenderer��: ������Ⱦͼ��
    */
    vtkNew<vtkRenderer> renderer;
    renderer->AddActor(actor); 
    renderer->SetBackground(0.5, 0.5, 0.5); 

    /* 7. vtkRenderWindow��: ������ʾ��Ⱦ���
    */
    vtkNew<vtkRenderWindow> renderWindow;
    renderWindow->AddRenderer(renderer); 
    renderWindow->SetSize(640, 480); 
    renderWindow->Render(); 

    /* 8. vtkInteractorStyleTrackballCamera��: �켣�򽻻���ʽ
       ����������ʵ���ӽ���ת��ƽ�ƺ����š�
    */
    vtkNew<vtkInteractorStyleTrackballCamera> style;

    /* 9. vtkRenderWindowInteractor��: �����û�����
    */
    vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
    renderWindowInteractor->SetRenderWindow(renderWindow); 
    renderWindowInteractor->SetInteractorStyle(style); 
    renderWindowInteractor->Initialize(); 
    renderWindowInteractor->Start(); 

    return 0;
}
