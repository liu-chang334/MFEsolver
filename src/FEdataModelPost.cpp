#include "../include/FEdataModelPost.h"


/**
 * @brief Construct a new FEDataModelPost object
 *
 * @param[in] feModel Finite element model
 */
FEDataModelPost::FEDataModelPost(FiniteElementModel feModel) 
{
    Node = feModel.Node;
    Element = feModel.Element;
}

/**
 * @brief Apply a color map to the mapper
 *
 * @param[in,out] mapper vtkMapper
 * @param[in] range range of the scalar
 * @note The color map is cool to warm color map
 */
void ApplyParaViewColorMap(vtkMapper* mapper, const double range[2]) {
    vtkNew<vtkDiscretizableColorTransferFunction> ctf;
    ctf->SetColorSpaceToDiverging();

    ctf->AddRGBPoint(range[0],               0.231, 0.298, 0.753); 
    // ctf->AddRGBPoint(0.4*(range[0]+range[1]),0.865, 0.865, 0.865);
    ctf->AddRGBPoint(range[1],               0.706, 0.016, 0.150);
    
    ctf->DiscretizeOn();
    ctf->SetNumberOfValues(10); 

    mapper->SetLookupTable(ctf);
    mapper->UseLookupTableScalarRangeOn();
}

/**
 * @brief Read the result from the txt file
 * 
 * @param[in] fieldname field name, "U" for displacement, "S" for stress
 * @note The result is stored in the member variablies: Displacement, Stress, etc.
 *      but only Displacement is implemented now
 */
void FEDataModelPost::ReadResult(std::string fieldname)
{
    std::string current_path = std::filesystem::current_path().string();
    std::string resultpath = current_path + "\\.." + "\\.." + "\\FEoutput";

    if (fieldname == "U")
    {
        std::string filename = "Displacement.txt";
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
    if (fieldname == "E")
    {
        std::string filename = "Strain.txt";
        std::string resultfullpath = resultpath + "\\" + filename;
        Strain = loadMatrixFromTXT(resultfullpath);
        // the strain matrix is already (node, 6)
        // Strain = result;
    }
    if (fieldname == "S")
    {
        std::string filename = "Stress.txt";
        std::string resultfullpath = resultpath + "\\" + filename;
        Stress = loadMatrixFromTXT(resultfullpath);
    }
}

/**
 * @brief Set the points object
 * 
 * @note The points are stored in the member variable: points
 *      The points are stored in the order of nodeID
 */
void FEDataModelPost::FEdataSetPoints()
{
    for (int i = 0; i < Node.rows(); i++)
    {
        points->InsertNextPoint(Node(i, 0), Node(i, 1), Node(i, 2));
    }
}

/**
 * @brief Set the grid object
 *
 * @note The grid is stored in the member variable: ugrid
 *      The grid is stored in the order of elementID
 */
void FEDataModelPost::FEdataSetGrid()
{
    FEdataSetPoints();
    ugrid->SetPoints(points);

    for (int i = 0; i < Element.rows(); i++)
    {
        vtkIdType pts[8] = {Element(i, 0)-1, Element(i, 1)-1, Element(i, 2)-1, Element(i, 3)-1, 
                Element(i, 4)-1, Element(i, 5)-1, Element(i, 6)-1, Element(i, 7)-1};
        ugrid->InsertNextCell(VTK_HEXAHEDRON, 8, pts);
    }
}

/**
 * @brief Set the grid scalar object
 *
 * @param[in] fieldname field name, "U" for displacement, "S" for stress
 * @note The grid scalar is stored in the member variable: scalar
 *      The grid scalar is stored in the order of nodeID
 */
void FEDataModelPost::FEdataSetGridScalar(std::string fieldname)
{
    ReadResult(fieldname);
    FEdataSetGrid();

    if (fieldname == "U")
    {
        scalar->SetName("Displacement");
        scalar->SetNumberOfComponents(3);
        for (int i = 0; i < Node.rows(); i++)
        {
            scalar->InsertNextTuple3(Displacement(i, 0), Displacement(i, 1), Displacement(i, 2));
        }
        ugrid->GetPointData()->SetVectors(scalar);
    }
    if (fieldname == "E")
    {
        scalar->SetName("Strain");
        scalar->SetNumberOfComponents(6);
        for (int i = 0; i < Node.rows(); i++)
        {
            scalar->InsertNextTuple6(Strain(i, 0), Strain(i, 1), Strain(i, 2), Strain(i, 3), Strain(i, 4), Strain(i, 5)); 
        } 
        ugrid->GetPointData()->AddArray(scalar);
    }
    if (fieldname == "S")
    {
        scalar->SetName("Stress");
        scalar->SetNumberOfComponents(6);
        for (int i = 0; i < Node.rows(); i++)
        {
            scalar->InsertNextTuple6(Stress(i, 0), Stress(i, 1), Stress(i, 2), Stress(i, 3), Stress(i, 4), Stress(i, 5)); 
        } 
        ugrid->GetPointData()->AddArray(scalar);
    }
}

/**
 * @brief Plot the grid
 * 
 * @note The grid is plotted without scalar
 */
void FEDataModelPost::FEdataPlot()
{
    FEdataSetGrid();
    vtkNew<vtkDataSetMapper> mapper;
    mapper->SetInputData(ugrid);

    vtkNew<vtkActor> actor;
    actor->SetMapper(mapper);
    actor->GetProperty()->SetRepresentationToWireframe();
    actor->GetProperty()->SetColor(0, 0, 0);
    actor->GetProperty()->SetLineWidth(2.0); 

    vtkNew<vtkRenderer> renderer;
    renderer->AddActor(actor);
    renderer->SetBackground(0.9, 0.9, 0.9);

    // vtkNew<vtkCamera> camera;
    // camera->SetParallelProjection(true);  // Enable orthographic projection
    // renderer->SetActiveCamera(camera);

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

/**
 * @brief Plot the grid with scalar
 *
 * @param[in] fieldname field name, "U" for displacement, "S" for stress
 * @param[in] component component of the scalar, 1 for x, 2 for y, 3 for z
 * @note Only displacement is implemented now
 */
void FEDataModelPost::FEdataPlotScalar(std::string fieldname, int component)
{
    FEdataSetGridScalar(fieldname);
    vtkNew<vtkArrayCalculator> calculator;
    calculator->SetInputData(ugrid);
    calculator->SetAttributeTypeToPointData();
    if (fieldname == "U")
    {
        calculator->AddVectorArrayName("U");
        if (component == 1)
        {
            calculator->SetFunction("U[0]");
            calculator->SetResultArrayName("U1");
        }else if (component == 2)
        {
            calculator->SetFunction("U[1]");
            calculator->SetResultArrayName("U2");
        }else if (component == 3)
        {
            calculator->SetFunction("U[2]");
            calculator->SetResultArrayName("U3");
        }
    }
    if (fieldname == "E")
    {
        calculator->AddScalarVariable("E11", "Strain", 0);
        calculator->AddScalarVariable("E22", "Strain", 1);
        calculator->AddScalarVariable("E33", "Strain", 2);
        calculator->AddScalarVariable("E12", "Strain", 3);
        calculator->AddScalarVariable("E13", "Strain", 4);
        calculator->AddScalarVariable("E23", "Strain", 5);
        if (component == 11)
        {
            calculator->SetFunction("E11");
            calculator->SetResultArrayName("E11");
        }else if (component == 22)
        {
            calculator->SetFunction("E22");
            calculator->SetResultArrayName("E22"); 
        }else if (component == 33)
        {
            calculator->SetFunction("E33");
            calculator->SetResultArrayName("E33");
        }else if (component == 12)
        {
            calculator->SetFunction("E12");
            calculator->SetResultArrayName("E12");
        }else if (component == 13)
        {
            calculator->SetFunction("E13");
            calculator->SetResultArrayName("E13"); 
        }
    }
    if (fieldname == "S")
    {
        calculator->AddScalarVariable("S11", "Stress", 0);
        calculator->AddScalarVariable("S22", "Stress", 1);
        calculator->AddScalarVariable("S33", "Stress", 2);
        calculator->AddScalarVariable("S12", "Stress", 3);
        calculator->AddScalarVariable("S13", "Stress", 4);
        calculator->AddScalarVariable("S23", "Stress", 5);
        if (component == 11)
        {
            calculator->SetFunction("S11");
            calculator->SetResultArrayName("S11");
        }else if (component == 22)
        {
            calculator->SetFunction("S22");
            calculator->SetResultArrayName("S22"); 
        }else if (component == 33)
        {
            calculator->SetFunction("S33");
            calculator->SetResultArrayName("S33");
        }else if (component == 12)
        {
            calculator->SetFunction("S12");
            calculator->SetResultArrayName("S12"); 
        }else if (component == 13)
        {
            calculator->SetFunction("S13");
            calculator->SetResultArrayName("S13");
        }else if (component == 23)
        {
            calculator->SetFunction("S23");
            calculator->SetResultArrayName("S23");
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
    scalarBar->UnconstrainedFontSizeOn(); 
    scalarBar->SetLookupTable(mapper->GetLookupTable());
    scalarBar->SetTitle(resultarrayname.c_str());
    scalarBar->GetTitleTextProperty()->SetFontSize(30);
    scalarBar->SetNumberOfLabels(10);
    scalarBar->SetLabelFormat("%.3e");
    scalarBar->GetLabelTextProperty()->SetFontSize(24);
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