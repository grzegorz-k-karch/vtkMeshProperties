#ifndef __vtkMeshProperties_h
#define __vtkMeshProperties_h

#include "vtkPolyDataAlgorithm.h"

class vtkMeshProperties : public vtkPolyDataAlgorithm
{
 public:
  static vtkMeshProperties *New();
  vtkTypeMacro(vtkMeshProperties, vtkPolyDataAlgorithm);
  void PrintSelf(ostream &os, vtkIndent indent);

 protected:
  int FillInputPortInformation(int port, vtkInformation* info);
  int RequestData(vtkInformation*, 
		  vtkInformationVector**, 
		  vtkInformationVector*);
};
#endif
