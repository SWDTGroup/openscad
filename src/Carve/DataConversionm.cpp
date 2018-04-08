#include "GeometryEvaluator.h"
#include "traverser.h"
#include "Tree.h"
#include "GeometryCache.h"
#include "CGALCache.h"
#include "Polygon2d.h"
#include "module.h"
#include "state.h"
#include "offsetnode.h"
#include "transformnode.h"
#include "polarizationNode.h"//add by Look
#include "decimationNode.h"   //add by zwbrush
#include "alignNode.h"   //add by zwbrush
#include "linearextrudenode.h"
#include "rotateextrudenode.h"
#include "csgnode.h"
#include "cgaladvnode.h"
#include "projectionnode.h"
#include "textnode.h"
#include "CGAL_Nef_polyhedron.h"
#include "cgalutils.h"
#include "rendernode.h"
#include "clipper-utils.h"
#include "polyset-utils.h"
#include "polyset.h"
#include "calc.h"
#include "printutils.h"
#include "svg.h"
#include "calc.h"
#include "dxfdata.h"


#include "DataConversionm.h"




namespace Shapetizer2
{
	vtkSmartPointer<vtkPolyData> polySetToPolyData(const PolySet *polyset)
	{
		vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
		vtkPoints *points = vtkPoints::New();
		polydata->SetPoints(points);
		points->Delete();
		vtkCellArray *cells = vtkCellArray::New();
		polydata->SetPolys(cells);
		cells->Delete();
		//merge points
		points->Allocate(polyset->polygons.size()*3 /2);
		vtkSmartPointer<vtkIncrementalPointLocator> locator = NULL;
		locator.TakeReference(vtkMergePoints::New());
		BoundingBox bbox = polyset->getBoundingBox();
		double bounds[6] = {(bbox.min)()[0], (bbox.max)()[0], (bbox.min)()[1], (bbox.max)()[1], (bbox.min)()[2], (bbox.max)()[2]};
		locator->InitPointInsertion(points, bounds);
		vtkIdType ptIds[3];
		for (int i = 0; i < polyset->polygons.size(); i++)
		{
			const Polygon &plyg = polyset->polygons[i];
			if (plyg.size() != 3)
			{
				printf("polyset's polygon is %i \n", plyg.size());
				throw "polyset must be tesslated";
			}
			for (int j = 0; j < 3; j++)
			{
				ptIds[j] = points->GetNumberOfPoints();
				points->InsertNextPoint((double*)&plyg[j]);
				//locator->InsertUniquePoint((double*)&plyg[j], ptIds[j]);
				//printf("id:%i  pt(%f %f %f)\n", ptIds[j], plyg[j][0], plyg[j][1], plyg[j][2]);
			}
			cells->InsertNextCell(3, ptIds);
		}

		//Look::savePolyData(polydata, "/root/polydataSave.stl");


		return polydata;
	}


	PolySet* polyDataToPolysetPtr(vtkSmartPointer<vtkPolyData> polydata)
	{
		Look::savePolyData(polydata, "/home/zwbrush/polydataSave.stl");

		PolySet *ps = new PolySet(3);
		polydata->BuildCells();
		
		for (vtkIdType cellId = 0; cellId < polydata->GetNumberOfCells(); cellId++)
		{
		  vtkIdType npts, *pts;
			polydata->GetCellPoints(cellId, npts, pts);
   		Polygon plyg;
   		plyg.reserve(3);
			for (int i = 0; i < 3; i++)
			{
				double pt[3];
				polydata->GetPoint(pts[i], pt);
				plyg.push_back(Vector3d(pt[0], pt[1], pt[2]));
			}
			ps->append_poly(plyg);
		}
		
		return ps;
	}


	LKPolygon polygon2dToLKPolygonValue(const Polygon2d *polygon2d)
	{
		LKPolygon lkPolygon;
		const Polygon2d::Outlines2d &outlines = polygon2d->outlines();
		lkPolygon.rings.resize(outlines.size());
		printf("rings' size is %i \n", lkPolygon.rings.size());
		for (int i = 0; i < outlines.size(); i++)
		{
			lkPolygon.rings[i].resize(outlines[i].vertices.size());
			for (int j = 0; j<lkPolygon.rings[i].size(); j++)
			{
				lkPolygon.rings[i][j] = NewLKPoint(outlines[i].vertices[j].x(), 0, outlines[i].vertices[j].y());
				//lkPolygon.rings[i][j].v.y = 0;
				//printf("pt(%f %f %f)\n", lkPolygon.rings[i][j].v.x, lkPolygon.rings[i][j].v.y, lkPolygon.rings[i][j].v.z);
			}
		}

		return lkPolygon;
	}

// 	Polygon2d LKPolygonToPolygon2dValue(LKPolygon *lkPolygon);


}

