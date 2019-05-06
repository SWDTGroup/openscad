#pragma once

#include "PolygonCarve.h"
#include "polyset.h"
#include "Polygon2d.h"


//TODO: delete
namespace Look {
	inline vtkSmartPointer<vtkPolyData> loadPolyData(std::string path)
	{
		std::string lowerPath(path);
		transform(path.begin(), path.end(), lowerPath.begin(), tolower);
		std::string::size_type indexNum = lowerPath.rfind('.');
		std::string extStr = lowerPath.substr( indexNum==std::string::npos ? lowerPath.length() : (indexNum+1) );
		if (extStr == "stl")
		{
			vtkSmartPointer<vtkSTLReader>  pSTLReader = vtkSmartPointer<vtkSTLReader>::New();
			//pSTLReader->MergingOff();
			pSTLReader->SetFileName(path.c_str());
			pSTLReader->Update();
			return pSTLReader->GetOutput()  ;
		}
		else if (extStr == "ply")
		{
			vtkSmartPointer<vtkPLYReader>  pPLYReader = vtkSmartPointer<vtkPLYReader>::New();
			//pSTLReader->MergingOff();
			pPLYReader->SetFileName(path.c_str());
			pPLYReader->Update();
			return pPLYReader->GetOutput();
		}
		else if (extStr == "vtk")
		{
			vtkSmartPointer<vtkPolyDataReader>  pVTKReader = vtkSmartPointer<vtkPolyDataReader>::New();
			//pSTLReader->MergingOff();
			pVTKReader->SetFileName(path.c_str());
			pVTKReader->Update();
			return pVTKReader->GetOutput();
		}
		else
		{
			fprintf(stderr, "mesh文件格式无法识别!\n");
			return NULL;
		}
	}



	inline void savePolyData(vtkSmartPointer<vtkPolyData> polyData, std::string path, bool isBinary = true)
	{
		std::string lowerPath(path);
		transform(path.begin(), path.end(), lowerPath.begin(), tolower);
		std::string::size_type indexNum = lowerPath.rfind('.');
		std::string extStr = lowerPath.substr( indexNum==std::string::npos ? lowerPath.length() : (indexNum+1) );
		if (extStr == "stl")
		{
			vtkSTLWriter *stlWriter = vtkSTLWriter::New();
			stlWriter->SetInputData(polyData);
			stlWriter->SetFileName( path.c_str() );
			if (isBinary) stlWriter->SetFileTypeToBinary();
			stlWriter->Update();
			stlWriter->Delete();
		}
		else if (extStr == "ply")
		{
			vtkPLYWriter *plyWriter = vtkPLYWriter::New();
			plyWriter->SetInputData(polyData);
			plyWriter->SetFileName(path.c_str());
			if (isBinary) plyWriter->SetFileTypeToBinary();
			plyWriter->Update();
			plyWriter->Delete();
		}
		else if (extStr == "vtk")
		{
			vtkPolyDataWriter *polydataWriter = vtkPolyDataWriter::New();
			polydataWriter->SetInputData( polyData );
			polydataWriter->SetFileName(path.c_str());
			if (isBinary) polydataWriter->SetFileTypeToBinary();
			polydataWriter->Update();
			polydataWriter->Delete();
		}
		else
		{
			fprintf(stderr, "要保存的文件格式不对！\n");
		}
	}
}


namespace Shapetizer2
{
	vtkSmartPointer<vtkPolyData> polySetToPolyData(const PolySet *polyset); //polyset must be tesslated


	PolySet* polyDataToPolysetPtr(vtkSmartPointer<vtkPolyData> polydata);


	LKPolygon polygon2dToLKPolygonValue(const Polygon2d *polygon2d);

	Polygon2d LKPolygonToPolygon2dValue(LKPolygon *lkPolygon);

}
