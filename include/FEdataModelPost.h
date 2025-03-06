#ifndef FEDataModelPost_H
#define FEDataModelPost_H

#include <iostream>
#include <filesystem>
#include <vector>
#include <stdexcept>
#include <vtkUnstructuredGrid.h>
#include <vtkDataSetMapper.h>
#include <vtkCellType.h>
#include <vtkProperty.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkAutoInit.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkCylinderSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkCamera.h>
#include <vtkNamedColors.h>
#include <vtkScalarBarActor.h> 
#include <vtkPointData.h> 
#include <vtkArrayCalculator.h>
#include <vtkColorTransferFunction.h>
#include <vtkDiscretizableColorTransferFunction.h>
#include <vtkTextProperty.h>

#include <Eigen/Dense>
#include "FiniteElementModel.h"    
#include "Tools.h"

class FEDataModelPost{
public:
    Eigen::MatrixXd Node;
    Eigen::MatrixXi Element;
    Eigen::MatrixXd Displacement;
    Eigen::MatrixXd Stress;
    Eigen::MatrixXd Strain;
    
    vtkNew<vtkPoints> points;
    vtkNew<vtkUnstructuredGrid> ugrid;
    vtkNew<vtkDoubleArray> scalar;
    // vtkNew<vtkDataSetMapper> mapper;
    // vtkNew<vtkActor> actor;
    // vtkNew<vtkRenderer> renderer;
    // vtkNew<vtkRenderWindow> renderWindow;
    // vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
    // vtkNew<vtkInteractorStyleTrackballCamera> style;


public:
    FEDataModelPost(FiniteElementModel feModel);
    void ReadResult(std::string fieldname);
    void FEdataSetPoints();
    void FEdataSetGrid();
    void FEdataSetGridScalar(std::string fieldname);
    void FEdataPlot();
    void FEdataPlotScalar(std::string fieldname, int component);
};
#endif