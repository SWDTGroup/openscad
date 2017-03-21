#include "polyset-utils.h"
#include "polyset.h"
#include "Polygon2d.h"
#include "printutils.h"
#include "GeometryUtils.h"
#include "Reindexer.h"
#include "grid.h"
#ifdef ENABLE_CGAL
#include "cgalutils.h"
#endif
#include <vtkPolyData.h>
#include <vtkTriangle.h>
#include <vtkCellArray.h>
#include <vtkSmartPointer.h>

#include <boost/foreach.hpp>

namespace PolysetUtils {

	// Project all polygons (also back-facing) into a Polygon2d instance.
  // It's important to select all faces, since filtering by normal vector here
	// will trigger floating point incertainties and cause problems later.
	Polygon2d *project(const PolySet &ps) {
		Polygon2d *poly = new Polygon2d;

		BOOST_FOREACH(const Polygon &p, ps.polygons) {
			Outline2d outline;
			BOOST_FOREACH(const Vector3d &v, p) {
				outline.vertices.push_back(Vector2d(v[0], v[1]));
			}
			poly->addOutline(outline);
		}
		return poly;
	}

/* Tessellation of 3d PolySet faces
	 
	 This code is for tessellating the faces of a 3d PolySet, assuming that
	 the faces are near-planar polygons.
	 
	 The purpose of this code is originally to fix github issue 349. Our CGAL
	 kernel does not accept polygons for Nef_Polyhedron_3 if each of the
	 points is not exactly coplanar. "Near-planar" or "Almost planar" polygons
	 often occur due to rounding issues on, for example, polyhedron() input.
	 By tessellating the 3d polygon into individual smaller tiles that
	 are perfectly coplanar (triangles, for example), we can get CGAL to accept
	 the polyhedron() input.
*/
	
/* Given a 3D PolySet with near planar polygonal faces, tessellate the
	 faces. As of writing, our only tessellation method is triangulation
	 using CGAL's Constrained Delaunay algorithm. This code assumes the input
	 polyset has simple polygon faces with no holes.
	 The tessellation will be robust wrt. degenerate and self-intersecting
*/
	void tessellate_faces(const PolySet &inps, PolySet &outps, float max_edge_length = -1)
	{
		  vtkSmartPointer<vtkPolyData> vtk_polyData =
   			 vtkSmartPointer<vtkPolyData>::New();

 		 vtkSmartPointer<vtkCellArray> vtk_cells =
   			 vtkSmartPointer<vtkCellArray>::New();	
 
		int degeneratePolygons = 0;

		// Build Indexed PolyMesh
		Reindexer<Vector3f> allVertices;
		std::vector<std::vector<IndexedFace> > polygons;

		BOOST_FOREACH(const Polygon &pgon, inps.polygons) {
			if (pgon.size() < 3) {
				degeneratePolygons++;
				continue;
			}
			if (pgon.size() == 3) { // Short-circuit
				
				
				if(max_edge_length > 0)
				{
					int ii=0;
					vtkSmartPointer<vtkTriangle> vtk_triangle =
  						  vtkSmartPointer<vtkTriangle>::New();
					BOOST_FOREACH (const Vector3d &v, pgon) 
					{
						// Create vertex indices and remove consecutive duplicate vertices
						int idx = allVertices.lookup(v.cast<float>());
 					 	vtk_triangle->GetPointIds()->SetId ( ii, idx );
  					}

					vtk_cells->InsertNextCell(vtk_triangle); 

				}
				else
					outps.append_poly(pgon);

				
				continue;
			}
			
			polygons.push_back(std::vector<IndexedFace>());
			std::vector<IndexedFace> &faces = polygons.back();
			faces.push_back(IndexedFace());
			IndexedFace &currface = faces.back();
			BOOST_FOREACH (const Vector3d &v, pgon) {
				// Create vertex indices and remove consecutive duplicate vertices
				int idx = allVertices.lookup(v.cast<float>());
				if (currface.empty() || idx != currface.back()) currface.push_back(idx);
			}
			if (currface.front() == currface.back()) currface.pop_back();
                    if (currface.size() < 3) {
                        faces.pop_back(); // Cull empty triangles
                        if (faces.empty()) polygons.pop_back(); // All faces were culled
                    }
		}

		// Tessellate indexed mesh
		const Vector3f *verts = allVertices.getArray();
		std::vector<IndexedTriangle> allTriangles;
		BOOST_FOREACH(const std::vector<IndexedFace> &faces, polygons) {
			std::vector<IndexedTriangle> triangles;
			if (faces[0].size() == 3) {
				triangles.push_back(IndexedTriangle(faces[0][0], faces[0][1], faces[0][2]));
			}
			else {
				bool err = GeometryUtils::tessellatePolygonWithHoles(verts, faces, triangles, NULL);
				if (!err) {
					BOOST_FOREACH(const IndexedTriangle &t, triangles) {
						
						if(max_edge_length > 0)
						{
							vtkSmartPointer<vtkTriangle> vtk_triangle =
  								  vtkSmartPointer<vtkTriangle>::New();
					 		vtk_triangle->GetPointIds()->SetId ( 0, t[0]);
					 		vtk_triangle->GetPointIds()->SetId ( 1, t[1]);
					 		vtk_triangle->GetPointIds()->SetId ( 2, t[2] );
							vtk_cells->InsertNextCell(vtk_triangle); 
						}
						else
						{
							outps.append_poly();
							outps.append_vertex(verts[t[0]]);
							outps.append_vertex(verts[t[1]]);
							outps.append_vertex(verts[t[2]]);
						}
	
					}
				}
			}
		}
		if (degeneratePolygons > 0) PRINT("WARNING: PolySet has degenerate polygons");
		
		if(max_edge_length > 0)
		{
			vtkSmartPointer<vtkPoints> vtk_points =
   				 vtkSmartPointer<vtkPoints>::New();
			const Vector3f*  e_pts = allVertices.getArray();
			for(unsigned int i=0;i<allVertices.size();i++)
			{
				vtk_points->InsertNextPoint((double)e_pts[i][0], (double)e_pts[i][1], (double)e_pts[i][2]);
			}
			vtk_polyData->SetPoints ( vtk_points);
  			vtk_polyData->SetPolys ( vtk_cells);

			//subdivide ploydata by max_edge_length
			
			 vtkPolyData *vtk_polyData;
			 vtkIdType npts, *pts;
			 for(vtkIdType i=0;i<vtk_polyData->GetNumberOfCells();i++)
			 {
				 vtk_polyData->GetCellPoints(i, npts, pts);
				 outps.append_poly();
				 for (int k = 0;k<3;k++)
				 {
					double* vertex = vtk_polyData->GetPoint(pts[k]);
				 	outps.append_vertex(Vector3f((float)vertex[0],(float)vertex[1],(float)vertex[2]));
				 }
			 }			

		}
	}

	bool is_approximately_convex(const PolySet &ps) {
#ifdef ENABLE_CGAL
		return CGALUtils::is_approximately_convex(ps);
#else
		return false;
#endif
	}

}
