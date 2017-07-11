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


//#include <boost/foreach.hpp>

namespace PolysetUtils {

	

	double getLen2BetweenAB_axis(const double A[3], const double B[3],  subdivide_axis axis)
	{
		switch(axis)
		{

		case subdivide_axis_x:
			return (A[0]-B[0])*(A[0]-B[0]);
		case subdivide_axis_y:
			return (A[1]-B[1])*(A[1]-B[1]);

		case subdivide_axis_z:
			return (A[2]-B[2])*(A[2]-B[2]);

		default:		
		{
			const double *orderPtPtrs[2];
			if (A[0] > B[0] || (A[0]==B[0] && (A[1]>B[1] || (A[1]==B[1]&&A[2]>B[2]))))
			{
				orderPtPtrs[0] = A;
				orderPtPtrs[1] = B;
			}
			else
			{
				orderPtPtrs[0] = B;
				orderPtPtrs[1] = A;
			}
			double v3d[3];
			v3d[0] = orderPtPtrs[0][0] - orderPtPtrs[1][0];
			v3d[1] = orderPtPtrs[0][1] - orderPtPtrs[1][1];
			v3d[2] = orderPtPtrs[0][2] - orderPtPtrs[1][2];
			return v3d[0]*v3d[0] + v3d[1]*v3d[1] + v3d[2]*v3d[2];
		}
		}
	}

	//三角形类

	struct SBTriangle {

		double vertices[3][3];
		double edgeLen2s[3];	// 0--[2-0]	1--[0-1] 2--[1-2]
		int longestEdgeId;

		void install(subdivide_axis axis)
		{
			double longestEdgeLen2 = -1;
			for (int preI = 2, i = 0; i < 3; preI=i,i++)
			{
				edgeLen2s[i] = getLen2BetweenAB_axis(vertices[preI], vertices[i], axis);
				if (edgeLen2s[i] > longestEdgeLen2)
				{
					longestEdgeLen2 = edgeLen2s[i];
					longestEdgeId = i;
				}
			}
		}
	};

	void recurseSubDivideTriangle(std::vector<SBTriangle> &leafTriangles, std::stack<SBTriangle> &waitingSubdivideTris, const SBTriangle &triangle, const double maxLen2,  subdivide_axis axis)
	{

		if (triangle.edgeLen2s[triangle.longestEdgeId] > maxLen2)
		{
			double newPt[3];
			int idB = triangle.longestEdgeId; //AB为最长边
			int idA = (idB + 2) % 3;
			int idC = (idB + 1) % 3;
			const double *ptA = triangle.vertices[idA];
			const double *ptB = triangle.vertices[idB];
			//const double *ptC = triangle.vertices[idC];
			newPt[0] = ptA[0] > ptB[0] ? (ptA[0] + ptB[0]) / 2 : (ptB[0] + ptA[0]) / 2;
			newPt[1] = ptA[1] > ptB[1] ? (ptA[1] + ptB[1]) / 2 : (ptB[1] + ptA[1]) / 2;
			newPt[2] = ptA[2] > ptB[2] ? (ptA[2] + ptB[2]) / 2 : (ptB[2] + ptA[2]) / 2;

			//?阜秩切?

			{

				SBTriangle newTri;
				memcpy(newTri.vertices[0], triangle.vertices[idA], sizeof(double)*3);
				memcpy(newTri.vertices[1], newPt, sizeof(double)*3);
				memcpy(newTri.vertices[2], triangle.vertices[idC], sizeof(double)*3);
				newTri.install(axis);
				waitingSubdivideTris.push(newTri);
			}
			//细分三角形2

			{

				SBTriangle newTri;
				memcpy(newTri.vertices[0], newPt, sizeof(double)*3);
				memcpy(newTri.vertices[1], triangle.vertices[idB], sizeof(double)*3);
				memcpy(newTri.vertices[2], triangle.vertices[idC], sizeof(double)*3);
				newTri.install(axis);
				waitingSubdivideTris.push(newTri);
			}
		}
		else
		{
			leafTriangles.push_back(triangle);
		}
	}

	void polyset_subdivide(const PolySet &inps, PolySet &outps, double max_edge_length,  subdivide_axis axis)
	{
		//printf("polyset_subdivide %d polygons, max edge length %lf\n", inps.polygons.size(), max_edge_length);							

		double maxLen2 = max_edge_length * max_edge_length;

		std::stack<SBTriangle> waitingSubdivideTris;
		///模型读取为三角形
		shared_ptr<PolySet> tess_ps;
		tess_ps.reset( new PolySet(3));
		tessellate_faces(inps,  *tess_ps);

		
		for (const auto  &pgon : tess_ps->polygons) {
			SBTriangle tri;
			int i=0;
for(const auto &v : pgon) {

//			BOOST_FOREACH (const Vector3d &v, pgon) {
				tri.vertices[i][0] = v[0];
				tri.vertices[i][1] = v[1];
				tri.vertices[i][2] = v[2];
				i++;
			}
			tri.install(axis);
			waitingSubdivideTris.push(tri);

		}

		std::vector<SBTriangle> leafTriangles;
		//leafTriangles.reserve(waitingSubdivideTris.size() * 1.5);
		while(!waitingSubdivideTris.empty())
		{
			SBTriangle theTri = waitingSubdivideTris.top();
			waitingSubdivideTris.pop();
			recurseSubDivideTriangle(leafTriangles, waitingSubdivideTris, theTri, maxLen2, axis);
		}

		for (unsigned int i = 0; i < leafTriangles.size(); i++)
		{
			SBTriangle &tri = leafTriangles[i];

			outps.append_poly();
			outps.append_vertex(tri.vertices[0][0], tri.vertices[0][1], tri.vertices[0][2]);
			outps.append_vertex(tri.vertices[1][0], tri.vertices[1][1], tri.vertices[1][2]);
			outps.append_vertex(tri.vertices[2][0], tri.vertices[2][1], tri.vertices[2][2]);
		}

		//printf("polyset_subdivide done. output %d %d polygons\n", outps.polygons.size(), leafTriangles.size());							

	}



	// Project all polygons (also back-facing) into a Polygon2d instance.
  // It's important to select all faces, since filtering by normal vector here
	// will trigger floating point incertainties and cause problems later.
	Polygon2d *project(const PolySet &ps) {
		Polygon2d *poly = new Polygon2d;

		for(const auto &p : ps.polygons) {
			Outline2d outline;
			for(const auto &v : p) {
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
	void tessellate_faces(const PolySet &inps, PolySet &outps)
	{
		int degeneratePolygons = 0;

		// Build Indexed PolyMesh
		Reindexer<Vector3f> allVertices;
		std::vector<std::vector<IndexedFace>> polygons;

		for (const auto &pgon : inps.polygons) {
			if (pgon.size() < 3) {
				degeneratePolygons++;
				continue;
			}
			
			polygons.push_back(std::vector<IndexedFace>());
			std::vector<IndexedFace> &faces = polygons.back();
			faces.push_back(IndexedFace());
			IndexedFace &currface = faces.back();
			for(const auto &v : pgon) {
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
		for(const auto &faces : polygons) {
			std::vector<IndexedTriangle> triangles;
			bool err = false;
			if (faces[0].size() == 3) {
				triangles.push_back(IndexedTriangle(faces[0][0], faces[0][1], faces[0][2]));
			}
			else {
				err = GeometryUtils::tessellatePolygonWithHoles(verts, faces, triangles, NULL);
			}
			if (!err) {
				for(const auto &t : triangles) {
					outps.append_poly();
					outps.append_vertex(verts[t[0]]);
					outps.append_vertex(verts[t[1]]);
					outps.append_vertex(verts[t[2]]);
				}
			}
		}
		if (degeneratePolygons > 0) PRINT("WARNING: PolySet has degenerate polygons");
	}

	bool is_approximately_convex(const PolySet &ps) {
#ifdef ENABLE_CGAL
		return CGALUtils::is_approximately_convex(ps);
#else
		return false;
#endif
	}
}
