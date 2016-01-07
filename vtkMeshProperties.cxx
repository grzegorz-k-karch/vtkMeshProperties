#include "vtkObjectFactory.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkConnectivityFilter.h"
#include "vtkCellData.h"
#include "vtkCellArray.h"
#include "vtkPointData.h"
#include "vtkDataArray.h"
#include "vtkDoubleArray.h"
#include "vtkPoints.h"
#include "vtkDataSetSurfaceFilter.h"
#include "vtkUnstructuredGrid.h"
#include "vtkTriangle.h"
#include "vtkIdTypeArray.h"
#include "vtkMeshProperties.h"

#include <iostream>
#include <vector>
#include <map>

vtkStandardNewMacro(vtkMeshProperties);

namespace
{
  //----------------------------------------------------------------------------
  double dot(double* a, double* b)
  { 
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  }
  //----------------------------------------------------------------------------
  void cross(double* a, double* b, double* c)
  {
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
  }
}

//----------------------------------------------------------------------------
int vtkMeshProperties::FillInputPortInformation(int port, vtkInformation* info)
{
  if (!this->Superclass::FillInputPortInformation(port, info)) {
    return 0;
  }
  if (port == 0) {
      info->Set( vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
      return 1;
  }
  return 0;
}

//----------------------------------------------------------------------------
int vtkMeshProperties::RequestData(vtkInformation *vtkNotUsed(request),
			     vtkInformationVector **inputVector,
			     vtkInformationVector *outputVector)
{
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkPolyData *input = 
    vtkPolyData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkSmartPointer<vtkConnectivityFilter> connectivityFilter = 
    vtkSmartPointer<vtkConnectivityFilter>::New();
  connectivityFilter->SetExtractionModeToAllRegions();
  connectivityFilter->ColorRegionsOn();
  connectivityFilter->SetInputData(input);
  connectivityFilter->Update();

  vtkSmartPointer<vtkDataSetSurfaceFilter> surfaceFilter =
    vtkSmartPointer<vtkDataSetSurfaceFilter>::New();  
  surfaceFilter->SetInputData(connectivityFilter->GetOutput());
  surfaceFilter->Update();
  vtkSmartPointer<vtkPolyData> surface = surfaceFilter->GetOutput();
  vtkSmartPointer<vtkDataArray> regionId = 
    surface->GetCellData()->GetArray("RegionId");

  const int numRegions = connectivityFilter->GetNumberOfExtractedRegions();
  const int numCells = surface->GetNumberOfCells();

  const vtkIdType *cellArrayData = surface->GetPolys()->GetPointer();
  int arrayIdx = 0;

  for (int i = 0; i < numCells; ++i) {
    int npts = cellArrayData[arrayIdx];
    if (npts != 3) {
      vtkErrorMacro("vtkMeshProperties works only with triangle polygons");
      return 0;
    }
    arrayIdx += 4;
  }

  std::vector<bool> closedSurfaces(numRegions, true);

  //--------------------------------------------------------------------------
  // test if meshes are closed
  {
    std::vector<std::vector<vtkIdType>> regions(numRegions);
    for (auto &reg : regions) {
      reg.clear();
    }

    arrayIdx = 0;
    for (int i = 0; i < numCells; ++i) {
      int pid[3] = {cellArrayData[arrayIdx+1],
		    cellArrayData[arrayIdx+2],
		    cellArrayData[arrayIdx+3]};
      int id = regionId->GetComponent(i, 0);

      regions[id].push_back(pid[0]);
      regions[id].push_back(pid[1]);
      regions[id].push_back(pid[2]);

      arrayIdx += 4;
    }

    int regIdx = 0;
    for (const auto &region : regions) {
      std::map<std::pair<vtkIdType,vtkIdType>,int> edges;
      edges.clear();
      
      for (int t = 0; t < region.size()/3; ++t) {
	for (int j = 0; j < 3; ++j) {
	  int e0 = region[t*3+j+0];
	  int e1 = region[t*3+(j+1)%3];
	  edges[std::pair<vtkIdType,vtkIdType>(std::min(e0,e1),std::max(e0,e1))] += 1;
	}
      }

      for (const auto &e : edges) {
	if (e.second != 2) {
	  closedSurfaces[regIdx] = false;
	}
      }
      ++regIdx;
    }
  }
  //--------------------------------------------------------------------------
  
  vtkSmartPointer<vtkDoubleArray> areasArray = 
    vtkSmartPointer<vtkDoubleArray>::New();
  areasArray->SetName("Area");
  areasArray->SetNumberOfComponents(1);
  areasArray->SetNumberOfTuples(numCells);
  std::vector<double> areas(numRegions, 0.0);

  vtkSmartPointer<vtkDoubleArray> volumesArray = 
    vtkSmartPointer<vtkDoubleArray>::New();
  volumesArray->SetName("Volume");
  volumesArray->SetNumberOfComponents(1);
  volumesArray->SetNumberOfTuples(numCells);
  std::vector<double> volumes(numRegions, 0.0);

  vtkSmartPointer<vtkPoints> points = surface->GetPoints();

  arrayIdx = 0;
  for (int i = 0; i < numCells; ++i) {

    double p[3][3];
    points->GetPoint(cellArrayData[arrayIdx+1], p[0]);
    points->GetPoint(cellArrayData[arrayIdx+2], p[1]);
    points->GetPoint(cellArrayData[arrayIdx+3], p[2]);
    int id = regionId->GetComponent(i, 0);
    
    double triaArea = vtkTriangle::TriangleArea(p[0],p[1],p[2]);
    areas[id] += triaArea;

    if (closedSurfaces[id] == true) {
      double cr[3];
      cross(p[1], p[2], cr);
      volumes[id] += dot(p[0], cr)/6.0;
    }

    arrayIdx += 4;
  }
  for (int i = 0; i < numCells; ++i) {
    int id = regionId->GetComponent(i, 0);
    areasArray->SetValue(i, areas[id]);
    volumesArray->SetValue(i, volumes[id]);
  }

  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  vtkPolyData *output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  output->CopyStructure(input);
  output->GetPointData()->PassData(input->GetPointData());
  output->GetCellData()->PassData(input->GetCellData());
  output->GetCellData()->AddArray(areasArray);
  output->GetCellData()->AddArray(volumesArray);

  return 1;
}

////////// External Operators /////////////
void vtkMeshProperties::PrintSelf(ostream &os, vtkIndent indent)
{
}
