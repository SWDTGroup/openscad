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

#include "Carve/DataConversionm.h"	//add by Look
#undef max(a, b)	//add by Look
#undef min(a, b)	//add by Look

#include <algorithm>
#include <queue>
#include <boost/foreach.hpp>
#include <boost/heap/binomial_heap.hpp>

#include <CGAL/convex_hull_2.h>
#include <CGAL/Point_2.h>

GeometryEvaluator::GeometryEvaluator(const class Tree &tree):
	tree(tree)
{
}

/*!
	Set allownef to false to force the result to _not_ be a Nef polyhedron
*/
shared_ptr<const Geometry> GeometryEvaluator::evaluateGeometry(const AbstractNode &node, 
																															 bool allownef)
{
	if (!GeometryCache::instance()->contains(this->tree.getIdString(node))) {
		shared_ptr<const CGAL_Nef_polyhedron> N;
		if (CGALCache::instance()->contains(this->tree.getIdString(node))) {
			N = CGALCache::instance()->get(this->tree.getIdString(node));
		}

		// If not found in any caches, we need to evaluate the geometry
		if (N) {
			this->root = N;
		}	
    else {
			Traverser trav(*this, node, Traverser::PRE_AND_POSTFIX);
			trav.execute();
		}

		if (!allownef) {
			if (shared_ptr<const CGAL_Nef_polyhedron> N = dynamic_pointer_cast<const CGAL_Nef_polyhedron>(this->root)) {
				PolySet *ps = new PolySet(3);
				ps->setConvexity(N->getConvexity());
				this->root.reset(ps);
                if (!N->isEmpty()) {
                    bool err = CGALUtils::createPolySetFromNefPolyhedron3(*N->p3, *ps);
                    if (err) {
                        PRINT("ERROR: Nef->PolySet failed");
                    }
                }

				smartCacheInsert(node, this->root);
			}
		}
		return this->root;
	}
	return GeometryCache::instance()->get(this->tree.getIdString(node));
}

GeometryEvaluator::ResultObject GeometryEvaluator::applyToChildren(const AbstractNode &node, OpenSCADOperator op)
{
	unsigned int dim = 0;
	/*
	BOOST_FOREACH(const Geometry::ChildItem &item, this->visitedchildren[node.index()]) {
		if (!item.first->modinst->isBackground() && item.second) {
			if (!dim) dim = item.second->getDimension();
			else if (dim != item.second->getDimension()) {
				PRINT("WARNING: Mixing 2D and 3D objects is not supported.");
				break;
			}
		}
	}
    if (dim == 2) {
        Polygon2d *p2d = applyToChildren2D(node, op);
        assert(p2d);
        return ResultObject(p2d);
    }
    else if (dim == 3) return applyToChildren3D(node, op);
	return ResultObject();
	*/
	if (op != OPENSCAD_CARVE)
	{
		BOOST_FOREACH(const Geometry::ChildItem &item, this->visitedchildren[node.index()]) {
			if (!item.first->modinst->isBackground() && item.second) {
				if (!dim) dim = item.second->getDimension();
				else if (dim != item.second->getDimension()) {
					PRINT("WARNING: Mixing 2D and 3D objects is not supported.");
					break;
				}
			}
		}
		if (dim == 2) {
			Polygon2d *p2d = applyToChildren2D(node, op);
			assert(p2d);
			return ResultObject(p2d);
		}
		else if (dim == 3) return applyToChildren3D(node, op);
		return ResultObject();
	}
	else
	{	//carve()
		bool isCarve = false;
		bool is2and3dimension = false;
		bool isFirst2dim = false;
		BOOST_FOREACH(const Geometry::ChildItem &item, this->visitedchildren[node.index()]) {
			if (!item.first->modinst->isBackground() && item.second) {
				if (!dim) dim = item.second->getDimension();
				else if (dim != item.second->getDimension()) {
					//PRINT("WARNING: Mixing 2D and 3D objects is not supported.");
					if (dim * item.second->getDimension() == 6 && dim + item.second->getDimension()==5)
					{
						is2and3dimension = true;
						if (dim == 2) isFirst2dim = true;
						else isFirst2dim = false;
					}
					else
					{
						PRINT("WARNING: Mixing 2D and 3D objects is not supported.");
					}
					break;
				}
			}
		}
		if (is2and3dimension && this->visitedchildren[node.index()].size()==2)
		{
			isCarve = true;
		}


		if (!isCarve)
		{
			PRINT("ERROR: carve function's arguments were wrong!");
		}
		else
		{
			//carve
			printf("begin carve polygon now!\n");
			shared_ptr<const Geometry> ps;
			shared_ptr<const class Geometry> geom_3d, geom_2d;
			if (!isFirst2dim) {
				geom_2d = this->visitedchildren[node.index()].back().second;
				geom_3d = this->visitedchildren[node.index()].front().second;
			} else {
				geom_2d = this->visitedchildren[node.index()].front().second;
				geom_3d = this->visitedchildren[node.index()].back().second;
			}
			//3d model
			shared_ptr<const PolySet> ps_3d = dynamic_pointer_cast<const PolySet>(geom_3d);
			if (ps_3d==NULL) 
			{

				shared_ptr<const CGAL_Nef_polyhedron> chN = dynamic_pointer_cast<const CGAL_Nef_polyhedron>(geom_3d);
				if (chN) 
				{
					PolySet *ps_new = new PolySet(3);
					bool err = CGALUtils::createPolySetFromNefPolyhedron3(*chN->p3, *ps_new);
					if (err) {
						PRINT("ERROR: Nef->PolySet failed");
					}
					else {
						ps_3d.reset(ps_new);
					}
				}
				else
				{
					PRINT("ERROR: Nef->PolySet failed  node's model data is NULL");
				}
			}
			/* */
			PolySet tessPSValue(3);
			PolysetUtils::tessellate_faces(*ps_3d.get(), tessPSValue);
			PolySet *newPs = new PolySet(3);
			PolysetUtils::polyset_subdivide(tessPSValue, *newPs, 8);
			shared_ptr<PolySet> ps_3d_dynamic;
			ps_3d_dynamic.reset(newPs);
			ps_3d = ps_3d_dynamic;


			//to vtkPolyData
			vtkSmartPointer<vtkPolyData> polyData = Shapetizer2::polySetToPolyData(ps_3d.get());
			//Look::savePolyData(polyData, "/root/polydataSave.stl");
			//to LKPolygon
			Shapetizer2::LKPolygon lkPolygon = Shapetizer2::polygon2dToLKPolygonValue((const Polygon2d*)geom_2d.get());
			// 		Shapetizer2::LKPolygon lkPolygon;
			// 		lkPolygon.rings.clear();
			// 		Shapetizer2::LKRing ring;
			// 		ring.push_back(Shapetizer2::NewLKPoint(0.273621*15, 0.000000, 3.839996));
			// 		ring.push_back(Shapetizer2::NewLKPoint(-0.245361*15, 0.000000, 3.839996));
			// 		ring.push_back(Shapetizer2::NewLKPoint(-0.245361*15, 0.000000, 0.000000));
			// 		ring.push_back(Shapetizer2::NewLKPoint(0.273621*15, 0.000000, 0.000000));
			// 		ring.push_back(Shapetizer2::NewLKPoint(-3, 0, -3));
			// 		ring.push_back(Shapetizer2::NewLKPoint(3, 0, -3));
			// 		ring.push_back(Shapetizer2::NewLKPoint(3, 0, 3));
			// 		ring.push_back(Shapetizer2::NewLKPoint(-3, 0, 3));
			//ring.push_back(Shapetizer2::NewLKPoint(-3, 0, -3));
			//		lkPolygon.rings.push_back(ring);


			//do carve
			CsgNode &csgNode = *(CsgNode *)&node;
			Shapetizer2::PolygonCarve plgCarve;
			plgCarve.setFontPolygon( lkPolygon );
			plgCarve.setDepth(csgNode.carve_depth);
			plgCarve.setSubdivideLen(csgNode.carve_subdivideLen);
			plgCarve.setModelPolyData(polyData);

			plgCarve.update();
			if (plgCarve.getErrCode() != 0)
			{
				fprintf(stderr, "carve error: %i %s\n", plgCarve.getErrCode(), plgCarve.getError());
				ps = geom_3d;
			}
			else
			{
				//Look::savePolyData(plgCarve.getResult(), "/root/newPolygonCarve.stl", true);
				//printf("carving is successed!\n");
				ps = dynamic_pointer_cast<const Geometry>(Shapetizer2::polyDataToPolysetPtr(plgCarve.getResult()));
			}


			return (ps);
		}
	}
	return ResultObject();
}

/*!
	Applies the operator to all child nodes of the given node.
	
	May return NULL or any 3D Geometry object (can be either PolySet or CGAL_Nef_polyhedron)
*/
GeometryEvaluator::ResultObject GeometryEvaluator::applyToChildren3D(const AbstractNode &node, OpenSCADOperator op)
{
	Geometry::ChildList children = collectChildren3D(node);
	if (children.size() == 0) return ResultObject();

	if (op == OPENSCAD_HULL) {
		PolySet *ps = new PolySet(3, true);

		if (CGALUtils::applyHull(children, *ps)) {
			return ps;
		}

		delete ps;
		return ResultObject();
	}
	
	// Only one child -> this is a noop
	if (children.size() == 1) return ResultObject(children.front().second);

	if (op == OPENSCAD_MINKOWSKI) {
		Geometry::ChildList actualchildren;
		BOOST_FOREACH(const Geometry::ChildItem &item, children) {
			if (!item.second->isEmpty()) actualchildren.push_back(item);
		}
		if (actualchildren.empty()) return ResultObject();
		if (actualchildren.size() == 1) return ResultObject(actualchildren.front().second);
		return ResultObject(CGALUtils::applyMinkowski(actualchildren));
	}

	CGAL_Nef_polyhedron *N = CGALUtils::applyOperator(children, op);
	// FIXME: Clarify when we can return NULL and what that means
	if (!N) N = new CGAL_Nef_polyhedron;
	return ResultObject(N);
}



/*!
	Apply 2D hull.

	May return an empty geometry but will not return NULL.
*/
Polygon2d *GeometryEvaluator::applyHull2D(const AbstractNode &node)
{
	std::vector<const Polygon2d *> children = collectChildren2D(node);
	Polygon2d *geometry = new Polygon2d();

	typedef CGAL::Point_2<CGAL::Cartesian<double> > CGALPoint2;
	// Collect point cloud
	std::list<CGALPoint2> points;
	BOOST_FOREACH(const Polygon2d *p, children) {
		BOOST_FOREACH(const Outline2d &o, p->outlines()) {
			BOOST_FOREACH(const Vector2d &v, o.vertices) {
				points.push_back(CGALPoint2(v[0], v[1]));
			}
		}
	}
	if (points.size() > 0) {
		// Apply hull
		std::list<CGALPoint2> result;
		CGAL::convex_hull_2(points.begin(), points.end(), std::back_inserter(result));

		// Construct Polygon2d
		Outline2d outline;
		BOOST_FOREACH(const CGALPoint2 &p, result) {
			outline.vertices.push_back(Vector2d(p[0], p[1]));
		}
		geometry->addOutline(outline);
	}
	return geometry;
}

Geometry *GeometryEvaluator::applyHull3D(const AbstractNode &node)
{
	Geometry::ChildList children = collectChildren3D(node);

	PolySet *P = new PolySet(3);
	if (CGALUtils::applyHull(children, *P)) {
		return P;
	}
	delete P;
	return NULL;
}

//add by zwbrush
Polygon2d *GeometryEvaluator::applyOutline2D(const AbstractNode &node)
{
	std::vector<const Polygon2d *> children = collectChildren2D(node);
	if (!children.empty()) {
		return ClipperUtils::applyOutline2D(children);
	}
	return NULL;
}
//end of add by zwbrush


Polygon2d *GeometryEvaluator::applyMinkowski2D(const AbstractNode &node)
{
	std::vector<const Polygon2d *> children = collectChildren2D(node);
	if (!children.empty()) {
		return ClipperUtils::applyMinkowski(children);
	}
	return NULL;
}

/*!
	Returns a list of Polygon2d children of the given node.
	May return empty Polygon2d object, but not NULL objects
*/
std::vector<const class Polygon2d *> GeometryEvaluator::collectChildren2D(const AbstractNode &node)
{
	std::vector<const Polygon2d *> children;
	BOOST_FOREACH(const Geometry::ChildItem &item, this->visitedchildren[node.index()]) {
		const AbstractNode *chnode = item.first;
		const shared_ptr<const Geometry> &chgeom = item.second;
		// FIXME: Don't use deep access to modinst members
		if (chnode->modinst->isBackground()) continue;

		// NB! We insert into the cache here to ensure that all children of
		// a node is a valid object. If we inserted as we created them, the 
		// cache could have been modified before we reach this point due to a large
		// sibling object. 
		smartCacheInsert(*chnode, chgeom);
		
		if (chgeom) {
			if (chgeom->getDimension() == 2) {
				const Polygon2d *polygons = dynamic_cast<const Polygon2d *>(chgeom.get());
				assert(polygons);
				children.push_back(polygons);
			}
			else {
				PRINT("WARNING: Ignoring 3D child object for 2D operation");
			}
		}
	}
	return children;
}

/*!
	Since we can generate both Nef and non-Nef geometry, we need to insert it into
	the appropriate cache.
	This method inserts the geometry into the appropriate cache if it's not already cached.
*/
void GeometryEvaluator::smartCacheInsert(const AbstractNode &node, 
																				 const shared_ptr<const Geometry> &geom)
{
	const std::string &key = this->tree.getIdString(node);

	shared_ptr<const CGAL_Nef_polyhedron> N = dynamic_pointer_cast<const CGAL_Nef_polyhedron>(geom);
	if (N) {
		if (!CGALCache::instance()->contains(key)) CGALCache::instance()->insert(key, N);
	}
	else {
		if (!GeometryCache::instance()->contains(key)) {
			if (!GeometryCache::instance()->insert(key, geom)) {
				PRINT("WARNING: GeometryEvaluator: Node didn't fit into cache");
			}
		}
	}
}

bool GeometryEvaluator::isSmartCached(const AbstractNode &node)
{
	const std::string &key = this->tree.getIdString(node);
	return (GeometryCache::instance()->contains(key) ||
					CGALCache::instance()->contains(key));
}

shared_ptr<const Geometry> GeometryEvaluator::smartCacheGet(const AbstractNode &node, bool preferNef)
{
	const std::string &key = this->tree.getIdString(node);
	shared_ptr<const Geometry> geom;
	bool hasgeom = GeometryCache::instance()->contains(key);
	bool hascgal = CGALCache::instance()->contains(key);
	if (hascgal && (preferNef || !hasgeom)) geom = CGALCache::instance()->get(key);
	else if (hasgeom) geom = GeometryCache::instance()->get(key);
	return geom;
}

/*!
	Returns a list of 3D Geometry children of the given node.
	May return empty geometries, but not NULL objects
*/
Geometry::ChildList GeometryEvaluator::collectChildren3D(const AbstractNode &node)
{
	Geometry::ChildList children;
	BOOST_FOREACH(const Geometry::ChildItem &item, this->visitedchildren[node.index()]) {
		const AbstractNode *chnode = item.first;
		const shared_ptr<const Geometry> &chgeom = item.second;
		// FIXME: Don't use deep access to modinst members
		if (chnode->modinst->isBackground()) continue;

		// NB! We insert into the cache here to ensure that all children of
		// a node is a valid object. If we inserted as we created them, the 
		// cache could have been modified before we reach this point due to a large
		// sibling object. 
		smartCacheInsert(*chnode, chgeom);
		
		if (chgeom) {
			if (chgeom->getDimension() == 2) {
				PRINT("WARNING: Ignoring 2D child object for 3D operation");
			}
			else if (chgeom->isEmpty() || chgeom->getDimension() == 3) {
				children.push_back(item);
			}
		}
	}
	return children;
}

/*!
	
*/
Polygon2d *GeometryEvaluator::applyToChildren2D(const AbstractNode &node, OpenSCADOperator op)
{
	//add by zwbrush
	if (op == OPENSCAD_OUTLINE) {
		return applyOutline2D(node);
	}
	//end of add by zwbrush

	if (op == OPENSCAD_MINKOWSKI) {
		return applyMinkowski2D(node);
	}
	else if (op == OPENSCAD_HULL) {
		return applyHull2D(node);
	}

	std::vector<const Polygon2d *> children = collectChildren2D(node);

	if (children.empty()) {
		return NULL;
	}

	if (children.size() == 1) {
		return new Polygon2d(*children[0]); // Copy
	}

	ClipperLib::ClipType clipType;
	switch (op) {
	case OPENSCAD_UNION:
		clipType = ClipperLib::ctUnion;
		break;
	case OPENSCAD_INTERSECTION:
		clipType = ClipperLib::ctIntersection;
		break;
	case OPENSCAD_DIFFERENCE:
		clipType = ClipperLib::ctDifference;
		break;
	default:
		PRINTB("Error: Unknown boolean operation %d", int(op));
		return NULL;
		break;
	}

	return ClipperUtils::apply(children, clipType);
}

/*!
	Adds ourself to out parent's list of traversed children.
	Call this for _every_ node which affects output during traversal.
	Usually, this should be called from the postfix stage, but for some nodes, 
	we defer traversal letting other components (e.g. CGAL) render the subgraph, 
	and we'll then call this from prefix and prune further traversal.

	The added geometry can be NULL if it wasn't possible to evaluate it.
*/
void GeometryEvaluator::addToParent(const State &state, 
																		const AbstractNode &node, 
																		const shared_ptr<const Geometry> &geom)
{
	this->visitedchildren.erase(node.index());
	if (state.parent()) {
		this->visitedchildren[state.parent()->index()].push_back(std::make_pair(&node, geom));
	}
	else {
		// Root node, insert into cache
		smartCacheInsert(node, geom);
		this->root = geom;
        assert(this->visitedchildren.empty());
	}
}

/*!
   Custom nodes are handled here => implicit union
*/
Response GeometryEvaluator::visit(State &state, const AbstractNode &node)
{
	if (state.isPrefix()) {
		if (isSmartCached(node)) return PruneTraversal;
		state.setPreferNef(true); // Improve quality of CSG by avoiding conversion loss
	}
	if (state.isPostfix()) {
		shared_ptr<const class Geometry> geom;
		if (!isSmartCached(node)) {
			geom = applyToChildren(node, OPENSCAD_UNION).constptr();
		}
		else {
			geom = smartCacheGet(node, state.preferNef());
		}
		addToParent(state, node, geom);
	}
	return ContinueTraversal;
}

Response GeometryEvaluator::visit(State &state, const OffsetNode &node)
{
	if (state.isPrefix() && isSmartCached(node)) return PruneTraversal;
	if (state.isPostfix()) {
		shared_ptr<const Geometry> geom;
		if (!isSmartCached(node)) {
			const Geometry *geometry = applyToChildren2D(node, OPENSCAD_UNION);
			if (geometry) {
				const Polygon2d *polygon = dynamic_cast<const Polygon2d*>(geometry);
				// ClipperLib documentation: The formula for the number of steps in a full
				// circular arc is ... Pi / acos(1 - arc_tolerance / abs(delta))
                double n = Calc::get_fragments_from_r(fabs(node.delta), node.fn, node.fs, node.fa);
				double arc_tolerance = fabs(node.delta) * (1 - cos(M_PI / n));
				const Polygon2d *result = ClipperUtils::applyOffset(*polygon, node.delta, node.join_type, node.miter_limit, arc_tolerance);
				assert(result);
				geom.reset(result);
				delete geometry;
			}
		}
		else {
			geom = smartCacheGet(node, false);
		}
		addToParent(state, node, geom);
	}
	return ContinueTraversal;
}

//add by zwbrush begin
struct struct_pt 
{
	int x;
	int y;
	double error;
};
	

struct cmp_pt{  
    bool operator()(const struct_pt& a, const struct_pt& b) const{  
        if(a.error > b.error)  return true;  
        return false;  
    }  
};


typedef typename boost::heap::binomial_heap < struct_pt, boost::heap::compare<cmp_pt> > pt_priority_queue;

struct struct_pt_link 
{
	Eigen::Vector2d pt;
	int prev ;
	int next;
	pt_priority_queue::handle_type self;
	bool exist;
};

static double computerError(Eigen::Vector2d pt, Eigen::Vector2d pt_prev, Eigen::Vector2d pt_next)
{
	Eigen::Vector2d e1 = pt  - pt_prev; 
 	Eigen::Vector2d e2 =pt_next  - pt; 
	double e2_long2 = e2[0] * e2[0] + e2[1]*e2[1];
	
	double priority_error;
	if( e2_long2 ==0.0)
		priority_error=0.0;
	else
	{
		//priority_error= e1.dot(e2)/ e2.dot(e2); //by angle.  normal ?
		priority_error= fabs(e1[0]*e2[1] - e1[1]*e2[0]) / sqrt(e2_long2); //by distance
	}
//	printf("priority_error %lf\n", priority_error);
	return priority_error;
}

static Polygon2d  *decimatePolygon2d(const DecimationNode &node, const Polygon2d &poly1)
{
	if(node.target==0)
		return new Polygon2d(poly1);

  if(poly1.outlines().size()==0)
  	return new Polygon2d(poly1);
  	
	const Polygon2d *p_poly = &poly1;

	if(node.keep_main)
		p_poly =  ClipperUtils::findLargest(poly1); 
	const Polygon2d& poly = *p_poly;	
	pt_priority_queue p_queue;
	std::vector<std::vector<struct_pt_link> > pt_vecs;
	pt_vecs.reserve(poly.outlines().size());

	unsigned int pt_num = 0;
	for (unsigned int o_idx = 0;o_idx < poly.outlines().size();o_idx++) {
			const Outline2d &o = poly.outlines()[o_idx];
			pt_num += o.vertices.size();
	}	
	if(pt_num <= node.target)
		return new Polygon2d(poly);

	for (unsigned int o_idx = 0;o_idx < poly.outlines().size();o_idx++) {
			
		
		const Outline2d &o = poly.outlines()[o_idx];
		std::vector<struct_pt_link> pt_vec;
		pt_vec.reserve(o.vertices.size());
		for (unsigned int v_idx = 0;v_idx < o.vertices.size();v_idx++) {

			int prev = v_idx == 0 ?  (o.vertices.size()-1) : (v_idx-1);
			int next = v_idx == o.vertices.size() -1 ?  0 : (v_idx+1);
			struct_pt pt_priority;
			pt_priority.x = o_idx;
			pt_priority.y =	 v_idx;
			pt_priority.error = computerError(o.vertices[v_idx], o.vertices[prev], o.vertices[next]);
	//		printf("##x= %d, y = %d, error = %lf\n", pt_priority.x, pt_priority.y, pt_priority.error);

			pt_priority_queue::handle_type ha =  p_queue.push(pt_priority);
			struct_pt_link link;
			link.exist = true;
			link.self = ha;
			link.pt = o.vertices[v_idx];
			
			if(v_idx>0) 
			{
				link.prev = prev;
				pt_vec[prev].next = v_idx;
			}

			if(v_idx == o.vertices.size() -1)
			{
				pt_vec[0].prev = v_idx;
				link.next = 0;
			}

			pt_vec.push_back(link) ;
			
				
		}
		pt_vecs.push_back(pt_vec);
	}
//	printf("origin size %d  dciamate to %d\n", p_queue.size(), node.target);
	if(node.target>0 && p_queue.size() > node.target)
	{
 		unsigned int decimate_num = p_queue.size() -  node.target ;
		for(unsigned int i =0;i<decimate_num ;i++)
		{
			const struct_pt pt_delete = p_queue.top();
			std::vector<struct_pt_link>& pt_vec=pt_vecs[pt_delete.x];
\
			pt_vec[pt_delete.y].exist = false;
	//		printf("1---x= %d, y = %d, error = %lf\n", pt_delete.x, pt_delete.y, pt_delete.error);


			int prev = pt_vec[pt_delete.y].prev;
			int next = pt_vec[pt_delete.y].next;
			p_queue.pop();
		
			//to recompute error

			if(pt_vec[prev].prev == next ) 
			{
				//if only two left, all erase
				p_queue.erase(pt_vec[prev].self);
				p_queue.erase(pt_vec[next].self);
				i +=2;
			}
			else
			{
		//		printf("2.2\n");
				struct_pt struct_pt_prev, struct_pt_next;
				struct_pt_prev.x = pt_delete.x;
				struct_pt_prev.y = prev;
				struct_pt_prev.error = computerError(pt_vec[prev].pt, pt_vec[pt_vec[prev].prev].pt, pt_vec[pt_vec[prev].next].pt);
	//			printf("2.2.1---x= %d, y = %d, error = %lf\n", struct_pt_prev.x, struct_pt_prev.y, struct_pt_prev.error);


				struct_pt_next.x = pt_delete.x;
				struct_pt_next.y = next;
				struct_pt_next.error = computerError(pt_vec[next].pt, pt_vec[pt_vec[next].prev].pt, pt_vec[pt_vec[next].next].pt);
	//			printf("2.2.2---x= %d, y = %d, error = %lf\n", struct_pt_next.x, struct_pt_next.y, struct_pt_next.error);

				p_queue.update(pt_vec[prev].self, struct_pt_prev);
				p_queue.update(pt_vec[next].self, struct_pt_next);

				pt_vec[next].prev = prev;
				pt_vec[prev].next = next;



			}
			

		}
	}
//	printf("remain size %d\n", p_queue.size());
	

	Polygon2d  new_poly;
	int cc = 0;
	for (unsigned int o_idx = 0;o_idx < pt_vecs.size();o_idx++) {
		
		std::vector<struct_pt_link> &pt_vec = pt_vecs[o_idx];
		
		const Outline2d &o = poly.outlines()[o_idx];
		Outline2d new_outline;
		new_outline.positive = o.positive;
		for (unsigned int v_idx = 0;v_idx < pt_vec.size();v_idx++) {
			if(pt_vec[v_idx].exist)
			{
				new_outline.vertices.push_back(o.vertices[v_idx]);
			}
		}
		if(new_outline.vertices.size()<3)
			continue;
		cc += new_outline.vertices.size();
		new_poly.addOutline(new_outline);
	}
	new_poly.setSanitized(true);
	
	
//	printf("remain size after %d\n", cc);
	

	ClipperLib::Clipper clipper;
	ClipperLib::Paths _paths = ClipperUtils::fromPolygon2d(new_poly);
	ClipperLib::CleanPolygons(_paths);
	clipper.AddPaths(_paths, ClipperLib::ptSubject, true);
	ClipperLib::PolyTree result;
	clipper.Execute(ClipperLib::ctUnion, result, ClipperLib::pftPositive);
	
	return ClipperUtils::toPolygon2d(result);
	
	

	

//	return new_poly;
}


Response GeometryEvaluator::visit(State &state, const DecimationNode &node)
{
	if (state.isPrefix() && isSmartCached(node)) return PruneTraversal;
	if (state.isPostfix()) {
		shared_ptr<const Geometry> geom;
		if (!isSmartCached(node)) {
			const Geometry *geometry = applyToChildren2D(node, OPENSCAD_UNION);
			if (geometry) {
				const Polygon2d *polygon = dynamic_cast<const Polygon2d*>(geometry);
				Polygon2d *result = decimatePolygon2d(node, *polygon );
				
				geom.reset(result);
				delete geometry;
				
			}
		}
		else {
			geom = smartCacheGet(node, false);
		}
		addToParent(state, node, geom);
	}
	return ContinueTraversal;
}


Transform3d GeometryEvaluator::getAlignMartix(const BoundingBox& bbox, const Vector3d& plane)
{
	Transform3d	_matrix = Transform3d::Identity();
	Vector3d translatevec(0,0,0);
	if(plane==Vector3d(0,0,0)) //center
	{
		translatevec = bbox.center() * -1 ;
	}	
	else 
	{
		if(plane.x()!=0)
			translatevec[0]= plane.x() > 0 ? -1 * bbox.min().x() : -1 * bbox.max().x();
		if(plane.y()!=0)
			translatevec[1]= plane.y() > 0 ? -1 * bbox.min().y() : -1 * bbox.max().y();
		if(plane.z()!=0)
			translatevec[2]= plane.z() > 0 ? -1 * bbox.min().z() : -1 * bbox.max().z();
	}
	
	_matrix.translate(translatevec);

	
	return _matrix;
} 



Response GeometryEvaluator::visit(State &state, const AlignNode &node)
{
	if (state.isPrefix() && isSmartCached(node)) return PruneTraversal;
	if (state.isPostfix()) {
		shared_ptr<const Geometry> geom;
		if (!isSmartCached(node)) {
			ResultObject res = applyToChildren(node, OPENSCAD_UNION);
				if ((geom = res.constptr())) {
					if (geom->getDimension() == 2) {
						
						shared_ptr<const Polygon2d> polygons = dynamic_pointer_cast<const Polygon2d>(geom);
						assert(polygons);
						
						// If we got a const object, make a copy
						shared_ptr<Polygon2d> newpoly;
						if (res.isConst()) newpoly.reset(new Polygon2d(*polygons));
						else newpoly = dynamic_pointer_cast<Polygon2d>(res.ptr());
						
						BoundingBox bbox= newpoly->getBoundingBox();
						Transform3d _matrix = getAlignMartix(bbox, node.m_plane);

						Transform2d mat2;
						mat2.matrix() << 
							_matrix(0,0), _matrix(0,1), _matrix(0,3),
							_matrix(1,0), _matrix(1,1), _matrix(1,3),
							_matrix(3,0), _matrix(3,1), _matrix(3,3);
						newpoly->transform(mat2);
					
						geom = newpoly;
						
						
					}
					else if (geom->getDimension() == 3) {
						shared_ptr<const PolySet> ps = dynamic_pointer_cast<const PolySet>(geom);
						if (ps) {
							// If we got a const object, make a copy
							shared_ptr<PolySet> newps;
							if (res.isConst()) newps.reset(new PolySet(*ps));
							else newps = dynamic_pointer_cast<PolySet>(res.ptr());
							
							BoundingBox bbox= newps->getBoundingBox();
							Transform3d _matrix = getAlignMartix(bbox, node.m_plane);

							newps ->transform(_matrix);
							geom = newps;
						}
						else {
							shared_ptr<const CGAL_Nef_polyhedron> N = dynamic_pointer_cast<const CGAL_Nef_polyhedron>(geom);
							assert(N);
							// If we got a const object, make a copy
							shared_ptr<CGAL_Nef_polyhedron> newN;
							if (res.isConst()) newN.reset((CGAL_Nef_polyhedron*)N->copy());
							else newN = dynamic_pointer_cast<CGAL_Nef_polyhedron>(res.ptr());
							
							BoundingBox bbox = newN->getBoundingBox();
							Transform3d _matrix = getAlignMartix(bbox, node.m_plane);
							newN->transform(_matrix);
							geom = newN;
						}
					}
				}
		}
		else {
			geom = smartCacheGet(node, false);
		}
		addToParent(state, node, geom);
	}
	return ContinueTraversal;
}




//add by zwbrush end


/*!
   RenderNodes just pass on convexity
*/
Response GeometryEvaluator::visit(State &state, const RenderNode &node)
{
	if (state.isPrefix()) {
		if (isSmartCached(node)) return PruneTraversal;
		state.setPreferNef(true); // Improve quality of CSG by avoiding conversion loss
	}
	if (state.isPostfix()) {
		shared_ptr<const class Geometry> geom;
		if (!isSmartCached(node)) {
			ResultObject res = applyToChildren(node, OPENSCAD_UNION);

			geom = res.constptr();
			if (shared_ptr<const PolySet> ps = dynamic_pointer_cast<const PolySet>(geom)) {
				// If we got a const object, make a copy
				shared_ptr<PolySet> newps;
				if (res.isConst()) newps.reset(new PolySet(*ps));
				else newps = dynamic_pointer_cast<PolySet>(res.ptr());
				newps->setConvexity(node.convexity);
				geom = newps;
			}
			else if (shared_ptr<const CGAL_Nef_polyhedron> N = dynamic_pointer_cast<const CGAL_Nef_polyhedron>(geom)) {
				// If we got a const object, make a copy
				shared_ptr<CGAL_Nef_polyhedron> newN;
				if (res.isConst()) newN.reset((CGAL_Nef_polyhedron*)N->copy());
				else newN = dynamic_pointer_cast<CGAL_Nef_polyhedron>(res.ptr());
				newN->setConvexity(node.convexity);
				geom = newN;
			}
		}
		else {
			geom = smartCacheGet(node, state.preferNef());
		}
		addToParent(state, node, geom);
	}
	return ContinueTraversal;
}

/*!
	Leaf nodes can create their own geometry, so let them do that

	input: None
	output: PolySet or Polygon2d
*/
Response GeometryEvaluator::visit(State &state, const LeafNode &node)
{
	if (state.isPrefix()) {
		shared_ptr<const Geometry> geom;
		if (!isSmartCached(node)) {
			const Geometry *geometry = node.createGeometry();
            assert(geometry);
			if (const Polygon2d *polygon = dynamic_cast<const Polygon2d*>(geometry)) {
				if (!polygon->isSanitized()) {
					Polygon2d *p = ClipperUtils::sanitize(*polygon);
					delete geometry;
					geometry = p;
				}
			}
            geom.reset(geometry);
		}
		else geom = smartCacheGet(node, state.preferNef());
		addToParent(state, node, geom);
	}
	return PruneTraversal;
}

Response GeometryEvaluator::visit(State &state, const TextNode &node)
{
	if (state.isPrefix()) {
		shared_ptr<const Geometry> geom;
		if (!isSmartCached(node)) {
			std::vector<const Geometry *> geometrylist = node.createGeometryList();
			std::vector<const Polygon2d *> polygonlist;
			BOOST_FOREACH(const Geometry *geometry, geometrylist) {
				const Polygon2d *polygon = dynamic_cast<const Polygon2d*>(geometry);
				assert(polygon);
				polygonlist.push_back(polygon);
			}
			geom.reset(ClipperUtils::apply(polygonlist, ClipperLib::ctUnion));
		}
		else geom = GeometryCache::instance()->get(this->tree.getIdString(node));
		addToParent(state, node, geom);
	}
	return PruneTraversal;
}


/*!
	input: List of 2D or 3D objects (not mixed)
	output: Polygon2d or 3D PolySet
	operation:
	  o Perform csg op on children
 */			
Response GeometryEvaluator::visit(State &state, const CsgNode &node)
{
	if (state.isPrefix()) {
		if (isSmartCached(node)) return PruneTraversal;
		state.setPreferNef(true); // Improve quality of CSG by avoiding conversion loss
	}
	if (state.isPostfix()) {
		shared_ptr<const Geometry> geom;
		if (!isSmartCached(node)) {
			geom = applyToChildren(node, node.type).constptr();
		}
		else {
			geom = smartCacheGet(node, state.preferNef());
		}
		addToParent(state, node, geom);
	}
	return ContinueTraversal;
}

/*!
	input: List of 2D or 3D objects (not mixed)
	output: Polygon2d or 3D PolySet
	operation:
	  o Union all children
	  o Perform transform
 */			
Response GeometryEvaluator::visit(State &state, const TransformNode &node)
{
	if (state.isPrefix() && isSmartCached(node)) return PruneTraversal;
	if (state.isPostfix()) {
		shared_ptr<const class Geometry> geom;
		if (!isSmartCached(node)) {
			if (matrix_contains_infinity(node.matrix) || matrix_contains_nan(node.matrix)) {
				// due to the way parse/eval works we can't currently distinguish between NaN and Inf
				PRINT("WARNING: Transformation matrix contains Not-a-Number and/or Infinity - removing object.");
			}
			else {
				// First union all children
				ResultObject res = applyToChildren(node, OPENSCAD_UNION);
				if ((geom = res.constptr())) {
					if (geom->getDimension() == 2) {
						shared_ptr<const Polygon2d> polygons = dynamic_pointer_cast<const Polygon2d>(geom);
						assert(polygons);
						
						// If we got a const object, make a copy
						shared_ptr<Polygon2d> newpoly;
						if (res.isConst()) newpoly.reset(new Polygon2d(*polygons));
						else newpoly = dynamic_pointer_cast<Polygon2d>(res.ptr());
						
						Transform2d mat2;
						mat2.matrix() << 
							node.matrix(0,0), node.matrix(0,1), node.matrix(0,3),
							node.matrix(1,0), node.matrix(1,1), node.matrix(1,3),
							node.matrix(3,0), node.matrix(3,1), node.matrix(3,3);
						newpoly->transform(mat2);
						// A 2D transformation may flip the winding order of a polygon.
						// If that happens with a sanitized polygon, we need to reverse
						// the winding order for it to be correct.
						if (newpoly->isSanitized() && mat2.matrix().determinant() <= 0) {
							geom.reset(ClipperUtils::sanitize(*newpoly));
						}
					}
					else if (geom->getDimension() == 3) {
						shared_ptr<const PolySet> ps = dynamic_pointer_cast<const PolySet>(geom);
						if (ps) {
							// If we got a const object, make a copy
							shared_ptr<PolySet> newps;
							if (res.isConst()) newps.reset(new PolySet(*ps));
							else newps = dynamic_pointer_cast<PolySet>(res.ptr());
							newps->transform(node.matrix);
							geom = newps;
						}
						else {
							shared_ptr<const CGAL_Nef_polyhedron> N = dynamic_pointer_cast<const CGAL_Nef_polyhedron>(geom);
							assert(N);
							// If we got a const object, make a copy
							shared_ptr<CGAL_Nef_polyhedron> newN;
							if (res.isConst()) newN.reset((CGAL_Nef_polyhedron*)N->copy());
							else newN = dynamic_pointer_cast<CGAL_Nef_polyhedron>(res.ptr());
							newN->transform(node.matrix);
							geom = newN;
						}
					}
				}
			}
		}
		else {
			geom = smartCacheGet(node, state.preferNef());
		}
		addToParent(state, node, geom);
	}
	return ContinueTraversal;
}

//add by Look begin
Response GeometryEvaluator::visit(State &state, const PolarizationNode &node)
{
	if (state.isPrefix() && isSmartCached(node)) return PruneTraversal;
	if (state.isPostfix()) {
		shared_ptr<const class Geometry> geom;
		if (!isSmartCached(node)) {
			if (!"YanZheng") {
				//
			}
			else {
				//printf("polarization - visit \n");
				// First union all children
				ResultObject res = applyToChildren(node, OPENSCAD_UNION);
				if ((geom = res.constptr())) {
					if (geom->getDimension() == 2) {
						PRINT("2d polygon could not use the funcation of polarization!\n");
					}
					else if (geom->getDimension() == 3) {
 						shared_ptr<PolySet> newps;
 						shared_ptr<const PolySet> ps = dynamic_pointer_cast<const PolySet>(geom);
 						if (ps==NULL) 
						{

							shared_ptr<const CGAL_Nef_polyhedron> chN = dynamic_pointer_cast<const CGAL_Nef_polyhedron>(geom);
							if (chN) 
							{
								PolySet *ps_new = new PolySet(3);
								bool err = CGALUtils::createPolySetFromNefPolyhedron3(*chN->p3, *ps_new);
								if (err) {
									PRINT("ERROR: Nef->PolySet failed");
								}
								else {
									newps.reset(ps_new);
								}
							}	
						}
						else 
						{
							// If we got a const object, make a copy
 							if (res.isConst()) newps.reset(new PolySet(*ps));
 							else newps = dynamic_pointer_cast<PolySet>(res.ptr());							
						}

							shared_ptr<PolySet> newps_subdivided;
							
							newps_subdivided.reset( new PolySet(3));
							double max_edge_lendgth = -1.0;

							BoundingBox bbox = newps->getBoundingBox();
							double o_size = bbox.max()[0] - bbox.min()[0];

							
							
							if(node.fn > 0.0)
								max_edge_lendgth =  o_size  / node.fn;
							else if(node.fa > 12)
								max_edge_lendgth =  o_size * node.fa / node.angle;
							else 
								max_edge_lendgth =  node.fs;

							printf("%lf , %lf, %lf, %lf, max_edge_lendgth %lf\n", o_size,node.fn,node.fa,  node.fs, max_edge_lendgth);
							PolysetUtils::polyset_subdivide(*newps, *newps_subdivided, max_edge_lendgth, PolysetUtils::subdivide_axis_x);
 														//注意：强制把极化周长设置为node模型的周长
							//BoundingBox bbox = newps_subdivided->getBoundingBox();
							//double o_size = bbox.max()[0] - bbox.min()[0];
							if(node.exp.empty())
								newps_subdivided->polarization(o_size, node.angle);
							else
								if(!newps_subdivided->lua_exp(node.exp))
									PRINTB("ERROR: lua -- %s" , node.exp.c_str());


 							geom = newps_subdivided; 
					}
				}
			}
		}
		else {
			geom = smartCacheGet(node, state.preferNef());
		}
		addToParent(state, node, geom);
	}
	return ContinueTraversal;
}
//add by Look end

static void translate_PolySet(PolySet &ps, const Vector3d &translation)
{
	BOOST_FOREACH(Polygon &p, ps.polygons) {
		BOOST_FOREACH(Vector3d &v, p) {
			v += translation;
		}
	}
}

static void add_slice(PolySet *ps, const Polygon2d &poly, 
											double rot1, double rot2, 
											double h1, double h2, 
											const Vector2d &scale1,
											const Vector2d &scale2)
{
	Eigen::Affine2d trans1(Eigen::Scaling(scale1) * Eigen::Rotation2D<double>(-rot1*M_PI/180));
	Eigen::Affine2d trans2(Eigen::Scaling(scale2) * Eigen::Rotation2D<double>(-rot2*M_PI/180));
	
	bool splitfirst = sin((rot1 - rot2)*M_PI/180) > 0.0;
	BOOST_FOREACH(const Outline2d &o, poly.outlines()) {
		Vector2d prev1 = trans1 * o.vertices[0];
		Vector2d prev2 = trans2 * o.vertices[0];
		for (size_t i=1;i<=o.vertices.size();i++) {
			Vector2d curr1 = trans1 * o.vertices[i % o.vertices.size()];
			Vector2d curr2 = trans2 * o.vertices[i % o.vertices.size()];
			ps->append_poly();
			
			// Make sure to split negative outlines correctly
			if (splitfirst xor !o.positive) {
				ps->insert_vertex(prev1[0], prev1[1], h1);
				ps->insert_vertex(curr2[0], curr2[1], h2);
				ps->insert_vertex(curr1[0], curr1[1], h1);
				if (scale2[0] > 0 || scale2[1] > 0) {
					ps->append_poly();
					ps->insert_vertex(curr2[0], curr2[1], h2);
					ps->insert_vertex(prev1[0], prev1[1], h1);
					ps->insert_vertex(prev2[0], prev2[1], h2);
				}
			}
			else {
				ps->insert_vertex(prev1[0], prev1[1], h1);
				ps->insert_vertex(prev2[0], prev2[1], h2);
				ps->insert_vertex(curr1[0], curr1[1], h1);
				if (scale2[0] > 0 || scale2[1] > 0) {
					ps->append_poly();
					ps->insert_vertex(prev2[0], prev2[1], h2);
					ps->insert_vertex(curr2[0], curr2[1], h2);
					ps->insert_vertex(curr1[0], curr1[1], h1);
				}
			}
			prev1 = curr1;
			prev2 = curr2;
		}
	}
}

/*!
	Input to extrude should be sanitized. This means non-intersecting, correct winding order
	etc., the input coming from a library like Clipper.
*/
static Geometry *extrudePolygon(const LinearExtrudeNode &node, const Polygon2d &poly)
{
//	PRINTB("brefore simplify %d", poly1.outlines().size());
//	Polygon2d &poly = *ClipperUtils::SimplifyPolygon2d(poly1);
//	PRINTB("after simplify %d", poly.outlines().size());
	bool cvx = poly.is_convex();
	PolySet *ps = new PolySet(3, !cvx ? boost::tribool(false) : node.twist == 0 ? boost::tribool(true) : unknown);
	ps->setConvexity(node.convexity);
	if (node.height <= 0) return ps;

	double h1, h2;

	if (node.center) {
		h1 = -node.height/2.0;
		h2 = +node.height/2.0;
	} else {
		h1 = 0;
		h2 = node.height;
	}

	PolySet *ps_bottom = poly.tessellate(); // bottom
	
	// Flip vertex ordering for bottom polygon
	BOOST_FOREACH(Polygon &p, ps_bottom->polygons) {
		std::reverse(p.begin(), p.end());
	}
	translate_PolySet(*ps_bottom, Vector3d(0,0,h1));
         
	ps->append(*ps_bottom);
	delete ps_bottom;
	
	if (node.scale_x > 0 || node.scale_y > 0) {
		Polygon2d top_poly(poly);
		Eigen::Affine2d trans(Eigen::Scaling(node.scale_x, node.scale_y) *
													 Eigen::Rotation2D<double>(-node.twist*M_PI/180));
		top_poly.transform(trans); // top
		PolySet *ps_top = top_poly.tessellate();
		translate_PolySet(*ps_top, Vector3d(0,0,h2));
		ps->append(*ps_top);
		delete ps_top;
	}

    size_t slices = node.slices;

	for (unsigned int j = 0; j < slices; j++) {
		double rot1 = node.twist*j / slices;
		double rot2 = node.twist*(j+1) / slices;
		double height1 = h1 + (h2-h1)*j / slices;
		double height2 = h1 + (h2-h1)*(j+1) / slices;
		Vector2d scale1(1 - (1-node.scale_x)*j / slices,
										1 - (1-node.scale_y)*j / slices);
		Vector2d scale2(1 - (1-node.scale_x)*(j+1) / slices,
										1 - (1-node.scale_y)*(j+1) / slices);
		add_slice(ps, poly, rot1, rot2, height1, height2, scale1, scale2);
	}

	//delete &poly;
	return ps;
}

/*!
	input: List of 2D objects
	output: 3D PolySet
	operation:
	  o Union all children
	  o Perform extrude
 */			
Response GeometryEvaluator::visit(State &state, const LinearExtrudeNode &node)
{
	if (state.isPrefix() && isSmartCached(node)) return PruneTraversal;
	if (state.isPostfix()) {
		shared_ptr<const Geometry> geom;
		if (!isSmartCached(node)) {
			const Geometry *geometry = NULL;
			if (!node.filename.empty()) {
				DxfData dxf(node.fn, node.fs, node.fa, node.filename, node.layername, node.origin_x, node.origin_y, node.scale_x);

				Polygon2d *p2d = dxf.toPolygon2d();
				if (p2d) geometry = ClipperUtils::sanitize(*p2d);
			}
			else {
				geometry = applyToChildren2D(node, OPENSCAD_UNION);
			}
			if (geometry) {
				const Polygon2d *polygons = dynamic_cast<const Polygon2d*>(geometry);
				Geometry *extruded = extrudePolygon(node, *polygons);
				assert(extruded);
				geom.reset(extruded);
				delete geometry;
			}
		}
		else {
			geom = smartCacheGet(node, false);
		}
		addToParent(state, node, geom);
	}
	return ContinueTraversal;
}

static void fill_ring(std::vector<Vector3d> &ring, const Outline2d &o, double a)
{
	for (unsigned int i=0;i<o.vertices.size();i++) {
		ring[i][0] = o.vertices[i][0] * sin(a);
		ring[i][1] = o.vertices[i][0] * cos(a);
		ring[i][2] = o.vertices[i][1];
	}
}

/*!
	Input to extrude should be clean. This means non-intersecting, correct winding order
	etc., the input coming from a library like Clipper.

	FIXME: We should handle some common corner cases better:
  o 2D polygon having an edge being on the Y axis:
	  In this case, we don't need to generate geometry involving this edge as it
		will be an internal edge.
  o 2D polygon having a vertex touching the Y axis:
    This is more complex as the resulting geometry will (may?) be nonmanifold.
		In any case, the previous case is a specialization of this, so the following
		should be handled for both cases:
		Since the ring associated with this vertex will have a radius of zero, it will
		collapse to one vertex. Any quad using this ring will be collapsed to a triangle.

	Currently, we generate a lot of zero-area triangles

*/
static Geometry *rotatePolygon(const RotateExtrudeNode &node, const Polygon2d &poly)
{
	PolySet *ps = new PolySet(3);
	ps->setConvexity(node.convexity);

	BOOST_FOREACH(const Outline2d &o, poly.outlines()) {
		double min_x = 0;
		double max_x = 0;
		BOOST_FOREACH(const Vector2d &v, o.vertices) {
			min_x = fmin(min_x, v[0]);
			max_x = fmax(max_x, v[0]);

			if ((max_x - min_x) > max_x && (max_x - min_x) > fabs(min_x)) {
				PRINTB("ERROR: all points for rotate_extrude() must have the same X coordinate sign (range is %.2f -> %.2f)", min_x % max_x);
				delete ps;
				return NULL;
			}
		}
		int fragments = Calc::get_fragments_from_r(max_x - min_x, node.fn, node.fs, node.fa);

		std::vector<Vector3d> rings[2];
		rings[0].resize(o.vertices.size());
		rings[1].resize(o.vertices.size());

		fill_ring(rings[0], o, -M_PI/2); // first ring
		for (int j = 0; j < fragments; j++) {
			double a = ((j+1)%fragments*2*M_PI) / fragments - M_PI/2; // start on the X axis
			fill_ring(rings[(j+1)%2], o, a);

			for (size_t i=0;i<o.vertices.size();i++) {
				ps->append_poly();
				ps->insert_vertex(rings[j%2][i]);
				ps->insert_vertex(rings[(j+1)%2][(i+1)%o.vertices.size()]);
				ps->insert_vertex(rings[j%2][(i+1)%o.vertices.size()]);
				ps->append_poly();
				ps->insert_vertex(rings[j%2][i]);
				ps->insert_vertex(rings[(j+1)%2][i]);
				ps->insert_vertex(rings[(j+1)%2][(i+1)%o.vertices.size()]);
			}
		}
	}
	return ps;
}

/*!
	input: List of 2D objects
	output: 3D PolySet
	operation:
	  o Union all children
	  o Perform extrude
 */			
Response GeometryEvaluator::visit(State &state, const RotateExtrudeNode &node)
{
	if (state.isPrefix() && isSmartCached(node)) return PruneTraversal;
	if (state.isPostfix()) {
		shared_ptr<const Geometry> geom;
		if (!isSmartCached(node)) {
			const Geometry *geometry = NULL;
			if (!node.filename.empty()) {
				DxfData dxf(node.fn, node.fs, node.fa, node.filename, node.layername, node.origin_x, node.origin_y, node.scale);
				Polygon2d *p2d = dxf.toPolygon2d();
				if (p2d) geometry = ClipperUtils::sanitize(*p2d);
			}
			else {
				geometry = applyToChildren2D(node, OPENSCAD_UNION);
			}
			if (geometry) {
				const Polygon2d *polygons = dynamic_cast<const Polygon2d*>(geometry);
				Geometry *rotated = rotatePolygon(node, *polygons);
				geom.reset(rotated);
				delete geometry;
			}
		}
		else {
			geom = smartCacheGet(node, false);
		}
		addToParent(state, node, geom);
	}
	return ContinueTraversal;
}

/*!
	Handles non-leaf PolyNodes; projection
*/
Response GeometryEvaluator::visit(State &state, const AbstractPolyNode &node)
{
	assert(false);
	return AbortTraversal;
}

/*!
	input: List of 3D objects
	output: Polygon2d
	operation:
	  o Union all children
		o Perform projection
 */			
Response GeometryEvaluator::visit(State &state, const ProjectionNode &node)
{
	if (state.isPrefix() && isSmartCached(node)) return PruneTraversal;
	if (state.isPostfix()) {
		shared_ptr<const class Geometry> geom;
		if (!isSmartCached(node)) {

			if (!node.cut_mode) {
				ClipperLib::Clipper sumclipper;
				BOOST_FOREACH(const Geometry::ChildItem &item, this->visitedchildren[node.index()]) {
					const AbstractNode *chnode = item.first;
					const shared_ptr<const Geometry> &chgeom = item.second;
					// FIXME: Don't use deep access to modinst members
					if (chnode->modinst->isBackground()) continue;

					const Polygon2d *poly = NULL;

// CGAL version of Geometry projection
// Causes crashes in createNefPolyhedronFromGeometry() for this model:
// projection(cut=false) {
//    cube(10);
//    difference() {
//      sphere(10);
//      cylinder(h=30, r=5, center=true);
//    }
// }
#if 0
					shared_ptr<const PolySet> chPS = dynamic_pointer_cast<const PolySet>(chgeom);
					const PolySet *ps2d = NULL;
					shared_ptr<const CGAL_Nef_polyhedron> chN = dynamic_pointer_cast<const CGAL_Nef_polyhedron>(chgeom);
					if (chN) chPS.reset(chN->convertToPolyset());
					if (chPS) ps2d = PolysetUtils::flatten(*chPS);
					if (ps2d) {
						CGAL_Nef_polyhedron *N2d = CGALUtils::createNefPolyhedronFromGeometry(*ps2d);
						poly = N2d->convertToPolygon2d();
					}
#endif

// Clipper version of Geometry projection
// Clipper doesn't handle meshes very well.
// It's better in V6 but not quite there. FIXME: stand-alone example.
#if 1
					// project chgeom -> polygon2d
					shared_ptr<const PolySet> chPS = dynamic_pointer_cast<const PolySet>(chgeom);
					if (!chPS) {
						shared_ptr<const CGAL_Nef_polyhedron> chN = dynamic_pointer_cast<const CGAL_Nef_polyhedron>(chgeom);
						if (chN) {
							PolySet *ps = new PolySet(3);
							bool err = CGALUtils::createPolySetFromNefPolyhedron3(*chN->p3, *ps);
							if (err) {
								PRINT("ERROR: Nef->PolySet failed");
							}
							else {
								chPS.reset(ps);
							}
						}
					}
					if (chPS) poly = PolysetUtils::project(*chPS);
#endif

					if (poly) {
						ClipperLib::Paths result = ClipperUtils::fromPolygon2d(*poly);
						// Using NonZero ensures that we don't create holes from polygons sharing
						// edges since we're unioning a mesh
						result = ClipperUtils::process(result, 
																					 ClipperLib::ctUnion, 
																					 ClipperLib::pftNonZero);
						// Add correctly winded polygons to the main clipper
						sumclipper.AddPaths(result, ClipperLib::ptSubject, true);
					}
				}
				ClipperLib::PolyTree sumresult;
				// This is key - without StrictlySimple, we tend to get self-intersecting results
				sumclipper.StrictlySimple(true);
				sumclipper.Execute(ClipperLib::ctUnion, sumresult, ClipperLib::pftNonZero, ClipperLib::pftNonZero);
				if (sumresult.Total() > 0) geom.reset(ClipperUtils::toPolygon2d(sumresult));
			}
			else {
				shared_ptr<const Geometry> newgeom = applyToChildren3D(node, OPENSCAD_UNION).constptr();
				if (newgeom) {
					shared_ptr<const CGAL_Nef_polyhedron> Nptr = dynamic_pointer_cast<const CGAL_Nef_polyhedron>(newgeom);
					if (!Nptr) {
						Nptr.reset(CGALUtils::createNefPolyhedronFromGeometry(*newgeom));
					}
					if (!Nptr->isEmpty()) {
						Polygon2d *poly = CGALUtils::project(*Nptr, node.cut_mode);
						if (poly) {
							poly->setConvexity(node.convexity);
							geom.reset(poly);
						}
					}
				}
			}
		}
		else {
			geom = smartCacheGet(node, false);
		}
		addToParent(state, node, geom);
	}
	return ContinueTraversal;
}		

/*!
	input: List of 2D or 3D objects (not mixed)
	output: any Geometry
	operation:
	  o Perform cgal operation
 */			
Response GeometryEvaluator::visit(State &state, const CgaladvNode &node)
{
	if (state.isPrefix() && isSmartCached(node)) return PruneTraversal;
	if (state.isPostfix()) {
		shared_ptr<const Geometry> geom;
		if (!isSmartCached(node)) {
			switch (node.type) {
			//add by zwbrush
			case OUTLINE:
			{
				geom = applyToChildren(node, OPENSCAD_OUTLINE).constptr();
				break;

			}
			//end of add by zwbrush
			case MINKOWSKI: {
				ResultObject res = applyToChildren(node, OPENSCAD_MINKOWSKI);
				geom = res.constptr();
				// If we added convexity, we need to pass it on
				if (geom && geom->getConvexity() != node.convexity) {
					shared_ptr<Geometry> editablegeom;
					// If we got a const object, make a copy
					if (res.isConst()) editablegeom.reset(geom->copy());
					else editablegeom = res.ptr();
					geom = editablegeom;
					editablegeom->setConvexity(node.convexity);
				}
				break;
			}
			case HULL: {
				geom = applyToChildren(node, OPENSCAD_HULL).constptr();
				break;
			}
			case RESIZE: {
				ResultObject res = applyToChildren(node, OPENSCAD_UNION);
				geom = res.constptr();
				if (geom) {
					shared_ptr<Geometry> editablegeom;
					// If we got a const object, make a copy
					if (res.isConst()) editablegeom.reset(geom->copy());
					else editablegeom = res.ptr();
					geom = editablegeom;

					shared_ptr<CGAL_Nef_polyhedron> N = dynamic_pointer_cast<CGAL_Nef_polyhedron>(editablegeom);
					if (N) {
						N->resize(node.newsize, node.autosize);
					}
					else {
						shared_ptr<Polygon2d> poly = dynamic_pointer_cast<Polygon2d>(editablegeom);
						if (poly) {
							poly->resize(Vector2d(node.newsize[0], node.newsize[1]),
													 Eigen::Matrix<bool,2,1>(node.autosize[0], node.autosize[1]));
						}
						else {
							shared_ptr<PolySet> ps = dynamic_pointer_cast<PolySet>(editablegeom);
							if (ps) {
								ps->resize(node.newsize, node.autosize);
							}
							else {
								assert(false);
							}
						}
					}
				}
				break;
			}
			default:
				assert(false && "not implemented");
			}
		}
		else {
			geom = smartCacheGet(node, state.preferNef());
		}
		addToParent(state, node, geom);
	}
	return ContinueTraversal;
}

Response GeometryEvaluator::visit(State &state, const AbstractIntersectionNode &node)
{
	if (state.isPrefix()) {
		if (isSmartCached(node)) return PruneTraversal;
		state.setPreferNef(true); // Improve quality of CSG by avoiding conversion loss
	}
	if (state.isPostfix()) {
		shared_ptr<const class Geometry> geom;
		if (!isSmartCached(node)) {
			geom = applyToChildren(node, OPENSCAD_INTERSECTION).constptr();
		}
		else {
			geom = smartCacheGet(node, state.preferNef());
		}
		addToParent(state, node, geom);
	}
	return ContinueTraversal;
}
