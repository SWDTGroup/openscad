#pragma once

class Polygon2d;
class PolySet;

namespace PolysetUtils {
	enum subdivide_axis {
		subdivide_axis_x,
		subdivide_axis_y,
		subdivide_axis_z,
		subdivide_axis_all
	};
	
	Polygon2d *project(const PolySet &ps);
	void tessellate_faces(const PolySet &inps, PolySet &outps);
	bool is_approximately_convex(const PolySet &ps);
	void polyset_subdivide(const PolySet &inps, PolySet &outps, double max_edge_length, subdivide_axis axis = subdivide_axis_all);
	

};
