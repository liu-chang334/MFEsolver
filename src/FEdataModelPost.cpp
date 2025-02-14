#include "../include/FEdataModelPost.h"

// constructor for FEDataModelPost
FEDataModelPost::FEDataModelPost(FiniteElementModel feModel) 
{
    Node = feModel.Node;
    Element = feModel.Element;
}

// cool to warm color map
void ApplyParaViewColorMap(vtkMapper* mapper, const double range[2]) {
    vtkNew<vtkColorTransferFunction> colorTransferFunction;
    colorTransferFunction->SetColorSpaceToDiverging();
    colorTransferFunction->AddRGBPoint(range[0], 0.231, 0.298, 0.753); // cool
    colorTransferFunction->AddRGBPoint((range[0] + range[1]) / 2, 0.865, 0.865, 0.865); // white
    colorTransferFunction->AddRGBPoint(range[1], 0.706, 0.016, 0.150); // warm

    mapper->SetLookupTable(colorTransferFunction);
    mapper->UseLookupTableScalarRangeOn();
}
// read result 
void FEDataModelPost::ReadResult(std::string fieldname)
{
    // std::string resultpath = "D://liuchang//FEsolvercxx//FEoutput";
    std::string current_path = std::filesystem::current_path().string();
    std::string resultpath = current_path + "\\.." + "\\.." + "\\FEoutput";

    if (fieldname == "U")
    {
        std::string filename = "U.txt";
        std::string resultfullpath = resultpath + "\\" + filename;
        Eigen::MatrixXd result = loadMatrixFromTXT(resultfullpath);
        // transform the displacement (3*node, 1) to (node, 3) 
        int node = Node.rows();
        int dof = 3;
        Displacement = Eigen::MatrixXd::Zero(node, dof);
        for (int i = 0; i < node; i++)
        {
            Displacement(i, 0) = result(3*i, 0);
            Displacement(i, 1) = result(3*i+1, 0);
            Displacement(i, 2) = result(3*i+2, 0);
        }
    }
    if (fieldname == "S")
    {
        std::cout << "Stress plotting is not implemented yet" << std::endl;
    }
}

// set the points
void FEDataModelPost::FEdataSetPoints()
{
    // set the points
    for (int i = 0; i < Node.rows(); i++)
    {
        points->InsertNextPoint(Node(i, 0), Node(i, 1), Node(i, 2));
    }
}

// set the points and grid
void FEDataModelPost::FEdataSetGrid()
{
    // set the points
    FEdataSetPoints();
    ugrid->SetPoints(points);
    // set the grid
    for (int i = 0; i < Element.rows(); i++)
    {
        vtkIdType pts[8] = {Element(i, 0)-1, Element(i, 1)-1, Element(i, 2)-1, Element(i, 3)-1, 
                Element(i, 4)-1, Element(i, 5)-1, Element(i, 6)-1, Element(i, 7)-1};
        ugrid->InsertNextCell(VTK_HEXAHEDRON, 8, pts);
    }
}

void FEDataModelPost::FEdataSetGridScalar(std::string fieldname)
{
    ReadResult(fieldname);
    // set the points and grid
    FEdataSetGrid();
    // set the scalar
    if (fieldname == "U")
    {
        scalar->SetName("Displacement");
        scalar->SetNumberOfComponents(3);
        for (int i = 0; i < Node.rows(); i++)
        {
            scalar->InsertNextTuple3(Displacement(i, 0), Displacement(i, 1), Displacement(i, 2));
        }
    }
    ugrid->GetPointData()->SetVectors(scalar);
}

void FEDataModelPost::FEdataPlot()
{
    // set the points and grid
    FEdataSetGrid();

    // set the mapper
    vtkNew<vtkDataSetMapper> mapper;
    mapper->SetInputData(ugrid);

    // set the actor
    vtkNew<vtkActor> actor;
    actor->SetMapper(mapper);
    actor->GetProperty()->SetRepresentationToWireframe();
    actor->GetProperty()->SetColor(0, 0, 0);
    actor->GetProperty()->SetLineWidth(2.0); 

    // set the renderer
    vtkNew<vtkRenderer> renderer;
    renderer->AddActor(actor);
    renderer->SetBackground(0.9, 0.9, 0.9);

    // vtkNew<vtkCamera> camera;
    // camera->SetParallelProjection(true);  // Enable orthographic projection
    // renderer->SetActiveCamera(camera);

    // set the render window
    vtkNew<vtkRenderWindow> renderWindow;
    renderWindow->AddRenderer(renderer);
    renderWindow->SetSize(1200, 800);
    renderWindow->Render();

    // set the render window interactor
    vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
    vtkNew<vtkInteractorStyleTrackballCamera> style;
    renderWindowInteractor->SetRenderWindow(renderWindow);
    renderWindowInteractor->SetInteractorStyle(style);
    renderWindowInteractor->Initialize();
    renderWindowInteractor->Start();
}


void FEDataModelPost::FEdataPlotScalar(std::string fieldname, int component)
{
    FEdataSetGridScalar(fieldname);
    vtkNew<vtkArrayCalculator> calculator;
    calculator->SetInputData(ugrid);
    calculator->SetAttributeTypeToPointData();
    if (fieldname == "U")
    {
        calculator->AddVectorArrayName("Displacement");
        if (component == 1)
        {
            calculator->SetFunction("Displacement[0]");
            calculator->SetResultArrayName("DisplacementX");
        }else if (component == 2)
        {
            calculator->SetFunction("Displacement[1]");
            calculator->SetResultArrayName("DisplacementY");
        }else if (component == 3)
        {
            calculator->SetFunction("Displacement[2]");
            calculator->SetResultArrayName("DisplacementZ");
        }
    }
    calculator->Update();

    double range[2];
    vtkDataSet* dataSet = vtkDataSet::SafeDownCast(calculator->GetOutput());
    std::string resultarrayname = calculator->GetResultArrayName();
    dataSet->GetPointData()->GetArray(resultarrayname.c_str())->GetRange(range);

    vtkNew<vtkDataSetMapper> mapper;
    mapper->SetInputConnection(calculator->GetOutputPort());
    mapper->SetScalarModeToUsePointFieldData();
    mapper->SelectColorArray(resultarrayname.c_str());
    mapper->SetScalarRange(range);

    ApplyParaViewColorMap(mapper, range);

    vtkNew<vtkNamedColors> colors;
    vtkNew<vtkActor> surfaceActor;
    surfaceActor->SetMapper(mapper);
    surfaceActor->GetProperty()->SetColor(colors->GetColor3d("Blue").GetData());
    vtkNew<vtkActor> wireframeActor;
    wireframeActor->SetMapper(mapper);
    wireframeActor->GetProperty()->SetRepresentationToWireframe();
    wireframeActor->GetProperty()->SetColor(colors->GetColor3d("Black").GetData());
    wireframeActor->GetProperty()->SetLineWidth(2.0);

    vtkNew<vtkRenderer> renderer;
    renderer->AddActor(surfaceActor);
    renderer->AddActor(wireframeActor);
    renderer->SetBackground(colors->GetColor3d("SlateGray").GetData());

    vtkNew<vtkScalarBarActor> scalarBar;
    scalarBar->SetLookupTable(mapper->GetLookupTable());
    scalarBar->SetTitle(resultarrayname.c_str());
    scalarBar->SetNumberOfLabels(10);
    scalarBar->SetLabelFormat("%.3e");
    renderer->AddActor2D(scalarBar);

    vtkNew<vtkRenderWindow> renderWindow;
    renderWindow->AddRenderer(renderer);
    renderWindow->SetSize(1200, 800);
    renderWindow->Render();
    vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
    vtkNew<vtkInteractorStyleTrackballCamera> style;
    renderWindowInteractor->SetRenderWindow(renderWindow);
    renderWindowInteractor->SetInteractorStyle(style);
    renderWindowInteractor->Initialize();
    renderWindowInteractor->Start();
}