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

// ��ʼ�� VTK ģ�飨���ھ�̬����ʱ��Ҫ��
VTK_MODULE_INIT(vtkRenderingOpenGL2);
VTK_MODULE_INIT(vtkInteractionStyle);

int main()
{
    /* 1. vtkPoints��: ���ڴ洢���ε��λ��
       InsertNextPoint(double x, double y, double z): ���һ���㵽�㼯�У��������Ϊ(x, y, z)��
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
       SetPoints(vtkPoints *): ���������еĵ����ݡ�
       InsertNextCell(int cellType, vtkIdType numPts, vtkIdType *pts): ���һ����Ԫ�������У�ָ����Ԫ���ͺͶ���������
    */
    vtkNew<vtkUnstructuredGrid> grid;
    grid->SetPoints(points); 
    vtkIdType pts[8] = {4, 5, 7, 6, 0, 1, 3, 2};
    grid->InsertNextCell(VTK_HEXAHEDRON, 8, pts);

    /* 3. vtkDataSetMapper��: �����ݼ�ӳ��Ϊ����ͼ��
       SetInputData(vtkDataSet *input): �����������ݣ����ڿ��ӻ���
    */
    vtkNew<vtkDataSetMapper> mapper;
    mapper->SetInputData(grid);

    /* 4. vtkActor��: ��ʾ����Ⱦ��ͼ�ζ���
       SetMapper(vtkMapper *): ����ӳ������������������ͼ�ζ��������
       GetProperty(): ��ȡͼ�ζ�������ԡ�
       SetColor(double r, double g, double b): ����ͼ�ζ������ɫ��
       SetEdgeVisibility(bool visibility): ����ͼ�ζ���ı߿ɼ��ԡ�
       SetEdgeColor(double r, double g, double b): ����ͼ�ζ���ı���ɫ��
       SetLineWidth(double width): ����ͼ�ζ�����߿�
    */
    vtkNew<vtkActor> actor;
    actor->SetMapper(mapper);
    actor->GetProperty()->SetColor(0.0, 0.0, 1.0);
    actor->GetProperty()->SetEdgeVisibility(true);
    actor->GetProperty()->SetEdgeColor(1.0, 0.0, 0.0);
    actor->GetProperty()->SetLineWidth(2.0);

    /* 5. vtkRenderer��: ������Ⱦͼ��
       AddActor(vtkProp *): ���ͼ�ζ�����Ⱦ����
       SetBackground(double r, double g, double b): ���ñ�����ɫ����ΧΪ0.0~1.0��
    */
    vtkNew<vtkRenderer> renderer;
    renderer->AddActor(actor); 
    renderer->SetBackground(0.5, 0.5, 0.5); 

    /* 6. vtkRenderWindow��: ������ʾ��Ⱦ���
       AddRenderer(vtkRenderer *): �����Ⱦ�������ڡ�
       SetSize(int width, int height): ���ô��ڵĿ�͸ߡ�
       Render(): ִ����Ⱦ��
    */
    vtkNew<vtkRenderWindow> renderWindow;
    renderWindow->AddRenderer(renderer); 
    renderWindow->SetSize(640, 480); 
    renderWindow->Render(); 

    /* 7. vtkInteractorStyleTrackballCamera��: �켣�򽻻���ʽ
       ����������ʵ���ӽ���ת��ƽ�ƺ����š�
    */
    vtkNew<vtkInteractorStyleTrackballCamera> style;

    /* 8. vtkRenderWindowInteractor��: �����û�����
       SetRenderWindow(vtkRenderWindow *): ���ù�������Ⱦ���ڡ�
       SetInteractorStyle(vtkInteractorObserver *): ���ý�����ʽ��
       Initialize(): ��ʼ����������
       Start(): ��ʼ�û�������
    */
    vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
    renderWindowInteractor->SetRenderWindow(renderWindow); 
    renderWindowInteractor->SetInteractorStyle(style); 
    renderWindowInteractor->Initialize(); 
    renderWindowInteractor->Start(); 

    return 0;
}
