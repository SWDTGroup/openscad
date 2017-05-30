#include "PolygonCarve.h"

// #include <glutess.h>
// #include <tess.h>


#undef _TEST_DEMO


namespace Shapetizer2
{
	// 添加点并返回点的索引, 直接加点，不剔除重复的点
	static inline vtkIdType AddPointReturnIndex(vtkPoints *points, NewLKPoint &point, bool isGrayColor = false, std::vector<unsigned char> *resultPointRGBsPtr = NULL)	// isGrayColor-是否为灰颜色顶点（false为白色）
	{
		points->InsertNextPoint((double *)&point);
		if (isGrayColor)
		{
			int newSize = points->GetNumberOfPoints()*3;
			resultPointRGBsPtr->resize(newSize);
			unsigned char grayColor[3] = RGB_COLOR_UNDER;
			(*resultPointRGBsPtr)[newSize-3] = grayColor[0];
			(*resultPointRGBsPtr)[newSize-2] = grayColor[1];
			(*resultPointRGBsPtr)[newSize-1] = grayColor[2];
		}
		return (points->GetNumberOfPoints()-1);
	}
	// 添加点并返回点的索引， 没有重复点
	static inline vtkIdType AddPointNoRepeatReturnIndex(std::vector<NewLKPoint> &points, NewLKPoint &point)
	{
		for (vtkIdType i = 0; i < points.size(); i++)
		{
			if (points[i].v.x == point.v.x && points[i].v.y == point.v.y && points[i].v.z == point.v.z)
			{
				return i;
			}
		}
		points.push_back(point);
		return (points.size()-1);
	}

	// 判断点是否在三角形内（2D）
	// 返回值代表含义：
	//   0 在多边形外
	//   1 点在顶点上
	//   2 点在边上
	//   3 点在多边形内
	enum PointWithPolygon{
		PointWithPolygonOutside = 0,
		PointWithPolygonInside = 1,
		PointWithPolygonOnVertex = 2,
		PointWithPolygonOnEdge = 3
	};
	static inline PointWithPolygon pointWithWorkTriangle(const NewLKPoint &pt, const STWorkTriangle &workTriangle, size_t &index)
	{
		// 判断点在顶点上
		for (int i=0;i<3; i++)
		{
			const double delta = 2.e-10;
			//if (x == workTriangle.pt[i].x && y == workTriangle.pt[i].z)
			if ( abs(pt.v.x - workTriangle.pt[i].v.x) < delta && abs(pt.v.z - workTriangle.pt[i].v.z) < delta )
			{
				index = i;
				return PointWithPolygonOnVertex;	//点在顶点上， 顶点非常靠近的情况没有处理
			}
		}

		// 判断点在水平边上
		for (int i = 0, j = 2; i < 3; i++)
		{
			const double delta = 1.e-10;
			const STWorkLine &line = workTriangle.line[j];
			//if( workTriangle.pt[i].z == workTriangle.pt[j].z && pt.z == workTriangle.pt[i].z && pt.x > min(workTriangle.pt[j].x, workTriangle.pt[i].x) && pt.x < max(workTriangle.pt[j].x, workTriangle.pt[i].x) )
			if( abs(line.a * pt.v.x + line.b * pt.v.z + line.c) < delta
				&& pt.v.x > line.minX - delta && pt.v.x < line.maxX + delta
				&& pt.v.z > line.minZ - delta && pt.v.z < line.maxZ + delta
				)
			{
				index = j;
				return PointWithPolygonOnEdge;
			}
			j=i;
		}


		bool  oddNodes = false;
		for (int i = 0, j=2; i < 3; i++)
		{
			if((workTriangle.pt[i].v.z < pt.v.z && workTriangle.pt[j].v.z >= pt.v.z
				||   workTriangle.pt[j].v.z < pt.v.z && workTriangle.pt[i].v.z >= pt.v.z)
				&& (workTriangle.pt[i].v.x <= pt.v.x || workTriangle.pt[j].v.x <= pt.v.x))
			{
				/* 有问题 * /
				double expValue = pt.x * workTriangle.line[j].a + pt.z * workTriangle.line[j].b + workTriangle.line[j].c;

				oddNodes ^= (expValue > 0);

				/* */
				//a. 两点顺序交换成同用情况
				//double intersectionXSubX;
				//if (workTriangle.pt[i].x > workTriangle.pt[j].x
				//	|| (workTriangle.pt[i].x == workTriangle.pt[j].x && workTriangle.pt[i].z > workTriangle.pt[j].z)
				//	)
				//{
				//	intersectionXSubX = workTriangle.pt[j].x+(pt.z-workTriangle.pt[j].z)/(workTriangle.pt[i].z-workTriangle.pt[j].z)*(workTriangle.pt[i].x-workTriangle.pt[j].x) - pt.x;
				//}
				//else
				//{
				//	intersectionXSubX = workTriangle.pt[i].x+(pt.z - workTriangle.pt[i].z)/(workTriangle.pt[j].z-workTriangle.pt[i].z)*(workTriangle.pt[j].x-workTriangle.pt[i].x) - pt.x;
				//}

				//b.
				//double intersectionXSubX = workTriangle.pt[i].x+(pt.z-workTriangle.pt[i].z)/(workTriangle.pt[j].z-workTriangle.pt[i].z)*(workTriangle.pt[j].x-workTriangle.pt[i].x) - pt.x;

				// 点与边非常靠近算直接算在多形内
				//if( abs(intersectionXSubX) < 1.e-9 )
				//{
				//	index = j;
				//	return PointWithPolygonOnEdge;
				//}

				//oddNodes ^= (intersectionXSubX < 0);
				oddNodes ^= (workTriangle.pt[i].v.x+(pt.v.z-workTriangle.pt[i].v.z)/(workTriangle.pt[j].v.z-workTriangle.pt[i].v.z)*(workTriangle.pt[j].v.x-workTriangle.pt[i].v.x) < pt.v.x);
				/**/
			}
			j = i;
		}
		return oddNodes ? PointWithPolygonInside : PointWithPolygonOutside; 
	}

	// 求直线方程 a*X + b * Z + c = 0
	static inline void getLinearEquationXZ(const NewLKPoint &pt1, const NewLKPoint &pt2, double &a, double &b, double &c)
	{
		if (pt1.v.z == pt2.v.z)
		{
			a = 0;
			b = 1;
			c = -pt1.v.z;
		}
		else 
		{
			a = 1;
			b = (pt2.v.x - pt1.v.x)/(pt1.v.z - pt2.v.z);
			c = (pt1.v.x * pt2.v.z - pt2.v.x * pt1.v.z) / (pt1.v.z - pt2.v.z);
			double hypotenuse = hypot(a, b);
			a /= hypotenuse;
			b /= hypotenuse;
			c /= hypotenuse;
		}
	}



	// 判断两线段是否相交
	//   交点在端点上判为false
	static inline bool isSegmentIntersectXZ(const STWorkLine &line1, const STWorkLine &line2)
	{
		const double &a = line1.a, &b = line1.b, &c = line1.c;
		const double &A = line2.a, &B = line2.b, &C = line2.c;
		const NewLKPoint &p1 = line1.pt[0],
			&p2 = line1.pt[1],
			&p3 = line2.pt[0],
			&p4 = line2.pt[1];


		// 端点离直线近，且点在包围盒内，则相交
		double exp_p3 = a * p3.v.x + b * p3.v.z + c;
		double exp_p4 = a * p4.v.x + b * p4.v.z + c;

// 		const double delta = 1.e-10;
// 		if (
// 			( abs(exp_p3) < delta
// 			  && p3.x > line1.minX - delta && p3.x < line1.maxX + delta
// 		      && p3.z > line1.minZ - delta && p3.z < line1.maxZ + delta
// 			)
// 			||
// 			( abs(exp_p4) < delta
// 			  && p4.x > line1.minX - delta && p4.x < line1.maxX + delta
// 			  && p4.z > line1.minZ - delta && p4.z < line1.maxZ + delta
// 			)
// 			)
// 		{
// 			return true;
// 		}

		if ( exp_p3 * exp_p4 < 0 )
		{
			double exp_p1 = A * p1.v.x + B * p1.v.z + C;
			double exp_p2 = A * p2.v.x + B * p2.v.z + C;

// 			if (
// 				( abs(exp_p1) < delta
// 				&& p1.x > line2.minX - delta && p1.x < line2.maxX + delta
// 				&& p1.z > line2.minZ - delta && p1.z < line2.maxZ + delta
// 				)
// 				||
// 				( abs(exp_p2) < delta
// 				&& p2.x > line2.minX - delta && p2.x < line2.maxX + delta
// 				&& p2.z > line2.minZ - delta && p2.z < line2.maxZ + delta
// 				)
// 				)
// 			{
// 				return true;
// 			}

			return exp_p1 * exp_p2 < 0 //(注意4点共线的问题)
				&& max(line1.minX, line2.minX) <= min(line1.maxX, line2.maxX)
				&& max(line1.minZ, line2.minZ) <= min(line1.maxZ, line2.maxZ);
		}
		else 
		{
			return false;
		}

	}


	//求两线段的交点 http://blog.csdn.net/dgq8211/article/details/7952825
	static inline double Cross(const NewLKPoint& p1,const NewLKPoint& p2,const NewLKPoint& p3,const NewLKPoint& p4)  
	{  
		return (p2.v.x-p1.v.x)*(p4.v.z-p3.v.z) - (p2.v.z-p1.v.z)*(p4.v.x-p3.v.x);  
	}  

	static inline double Area(const NewLKPoint& p1,const NewLKPoint& p2,const NewLKPoint& p3)  
	{  
		return Cross(p1,p2,p1,p3);  
	}  

	static inline double fArea(const NewLKPoint& p1,const NewLKPoint& p2,const NewLKPoint& p3)  
	{  
		return abs(Area(p1,p2,p3));  
	}  
	static inline NewLKPoint segmentIntersetionXZ(const NewLKPoint& pt1,const NewLKPoint& pt2,const NewLKPoint& p3,const NewLKPoint& p4)  
	{
		// p1 p2设置为通用情况-p1在p2的左下，避免算出的交点不同 （p3, p4）为pt1 pt2 传入顺序不会变
		bool notNeedSwap_pt1_pt2 = pt1.v.x < pt2.v.x || (pt1.v.x==pt2.v.x && pt1.v.z < pt2.v.z);
		const NewLKPoint& p1 = notNeedSwap_pt1_pt2 ? pt1 : pt2;
		const NewLKPoint& p2 = notNeedSwap_pt1_pt2 ? pt2 : pt1;
		//const STPoint& p1 = pt1;
		//const STPoint& p2 = pt2;

		// fArea(p1,p2,p4)可能为0
		//double k = fArea(p1,p2,p3) / fArea(p1,p2,p4);
		//STPoint intersection;
		//intersection.x = (p3.x + k*p4.x)/(1+k);
		//intersection.z = (p3.z + k*p4.z)/(1+k);

		double s1 = fArea(p1,p2,p3);
		double s2 = fArea(p1,p2,p4);
		NewLKPoint intersection;
		intersection.v.x = (p4.v.x*s1+p3.v.x*s2)/(s1+s2);
		intersection.v.z = (p4.v.z*s1+p3.v.z*s2)/(s1+s2);

		// 给交?愕膟赋值
		NewLKPoint vectorP1P2;
		NewLKPoint vectorP1Insert;
		vectorP1P2.v.x = p2.v.x - p1.v.x;
		vectorP1P2.v.z = p2.v.z - p1.v.z;
		vectorP1P2.v.y = p2.v.y - p1.v.y;
		vectorP1Insert.v.x = intersection.v.x - p1.v.x;
		vectorP1Insert.v.z = intersection.v.z - p1.v.z;
		if (vectorP1P2.v.x == 0 && vectorP1P2.v.z == 0)
			intersection.v.y = p1.v.y;
		else
			intersection.v.y = p1.v.y + vectorP1P2.v.y * sqrt((vectorP1Insert.v.x * vectorP1Insert.v.x + vectorP1Insert.v.z * vectorP1Insert.v.z) /(vectorP1P2.v.x * vectorP1P2.v.x + vectorP1P2.v.z * vectorP1P2.v.z));

		return intersection;
	}
	
	//判断线段是否与矩形相交（XZ平面上）
	bool isSegmentCrossRectangeXZ(double pt1[3], double pt2[3], double bounds[6])
	{
		if (   (pt1[0] >= bounds[0] && pt1[0] <= bounds[1] && pt1[2] >= bounds[4] && pt1[2] <= bounds[5])
			|| (pt2[0] >= bounds[0] && pt2[0] <= bounds[1] && pt2[2] >= bounds[4] && pt2[2] <= bounds[5])
			)
			return true;
		else
		{
			STWorkLine line1(*(NewLKPoint *)pt1, *(NewLKPoint *)pt2);
			NewLKPoint pts[4];
			pts[0].v.x = bounds[0];
			pts[0].v.z = bounds[5];
			pts[1].v.x = bounds[1];
			pts[1].v.z = bounds[5];
			pts[2].v.x = bounds[1];
			pts[2].v.z = bounds[4];
			pts[3].v.x = bounds[0];
			pts[3].v.z = bounds[4];
			for (int i = 0, j = 3; i < 4; j = i, i++)
			{
				STWorkLine line2(pts[j], pts[i]);
				if (isSegmentIntersectXZ(line1, line2))
				{
					return true;
				}
			}

			return false;
		}
	}

	// 是否在多边形内
	bool PolygonMesh::pointInsideXZ2D(double pt[3])
	{
		bool  oddNodes = false;
		for (size_t k = 0; k < boundEdgeNodess.size(); k++)
		{
			BoundsEdges &boundEdges = boundEdgeNodess[k];
			for (int i = 0;i<boundEdges.size(); i++)
			{
				NewLKPoint &ptPrev = points[boundEdges[i]->ptIdFirst];
				NewLKPoint &ptCur = points[boundEdges[i]->ptId];
				if((ptCur.v.z< pt[2] && ptPrev.v.z>=pt[2]
				||   ptPrev.v.z<pt[2] && ptCur.v.z>=pt[2])
					&& (ptCur.v.x<=pt[0] || ptPrev.v.x<=pt[0]))
				{
					oddNodes ^= (ptCur.v.x+(pt[2]-ptCur.v.z)/(ptPrev.v.z-ptCur.v.z)*(ptPrev.v.x-ptCur.v.x)<pt[0]);
				}
			}
		}
		return oddNodes; 
	}



	// triangulate 三角化两个点在三角形内部的情况，两点可能在三角形顶点上
	static inline void triangulateFivePoint(STWorkTriangle &workTriangle, STWorkLine line, NewLKPoint &pt1, NewLKPoint &pt2,
		PointWithPolygon &pt1WithTriangle, PointWithPolygon &pt2WithTriangle,
		size_t &pt1OnIndex, size_t &pt2OnIndex,
		std::vector<STWorkTriangle> &tempNewWorkTris,
		std::vector<STWorkTriangle *> &tempEraseTris,
		std::vector<NewLKPoint> &intersections);
	static inline void triangulateCross_1_2(STWorkTriangle &workTriangle, STWorkLine line, NewLKPoint &pt1, NewLKPoint &pt2,
		PointWithPolygon &pt1WithTriangle, PointWithPolygon &pt2WithTriangle,
		size_t &pt1OnIndex, size_t &pt2OnIndex,
		std::vector<STWorkTriangle> &tempNewWorkTris,
		std::vector<STWorkTriangle *> &tempEraseTris,
		std::vector<NewLKPoint> &intersections
		);
	static inline void triangulateCross_2(STWorkTriangle &workTriangle, STWorkLine line, NewLKPoint &pt1, NewLKPoint &pt2,
		PointWithPolygon &pt1WithTriangle, PointWithPolygon &pt2WithTriangle,
		size_t &pt1OnIndex, size_t &pt2OnIndex,
		std::vector<STWorkTriangle> &tempNewWorkTris,
		std::vector<STWorkTriangle *> &tempEraseTris,
		std::vector<NewLKPoint> &intersections
		);




	PolygonCarve::PolygonCarve()
	{
		errCode = 0;
		result = NULL;
		_subdivideLen2 = 0.2 * 0.2;
		fontPolygonY = -1e20;
		bspTree = NULL;
	}


	PolygonCarve::~PolygonCarve(void)
	{


	}


	void PolygonMesh::getFontBounds(double bounds[6])
	{
		if (!points.empty())
		{
			bounds[0] = bounds[1] = points[0].v.x;
			bounds[2] = bounds[3] = points[0].v.y;
			bounds[4] = bounds[5] = points[0].v.z;
		}
		else
		{
			bounds[0] = bounds[1] = bounds[2] = bounds[3] = bounds[4] = bounds[5] = 0.0;
		}

		for (size_t i = 1; i < points.size(); i++)
		{
			NewLKPoint &point = points[i];
			if (bounds[0] > point.v.x) bounds[0] = point.v.x;
			if (bounds[2] > point.v.y) bounds[2] = point.v.y;
			if (bounds[4] > point.v.z) bounds[4] = point.v.z;
			if (bounds[1] < point.v.x) bounds[1] = point.v.x;
			if (bounds[3] < point.v.y) bounds[3] = point.v.y;
			if (bounds[5] < point.v.z) bounds[5] = point.v.z;
		}
	}


	void PolygonCarve::getFontBounds(double bounds[6])
	{
		fontPolygonMesh.getFontBounds(bounds);

		if (wordDepth < 0.0)
		{
			bounds[3] -= wordDepth;
		}
		else 
		{
			bounds[2] -= wordDepth;
		}
	}



	struct SubdivideLine{
		NewLKPoint pts[2];
		double length2;
		SubdivideLine(NewLKPoint pt0, NewLKPoint pt1) {
			pts[0] = pt0; pts[1] = pt1;
		}
		SubdivideLine(const struct SubdivideLine &line) {
			*this = line;
		}
		void calcLen2(){
			double a = pts[1].v.x - pts[0].v.x;
			double b = pts[1].v.z - pts[0].v.z;
			length2 = a * a + b * b;
		}
	};

	bool PolygonCarve::update()
	{
		if (errCode) return false;

		// 变换
		if (modelPolyDataOrigin)
		{
			vtkSmartPointer<vtkTransformPolyDataFilter> tpdFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
			tpdFilter->SetTransform( allTransformToVerticality );
			tpdFilter->SetInputData( modelPolyDataOrigin );
			tpdFilter->Update();
			modelPolyData = tpdFilter->GetOutput();
		}

		resultPoints = vtkSmartPointer<vtkPoints>::New();


		if (!this->isPicCarving)
		{
			//// 删除邻边交叠的情况 —— "w"字母特殊处理
			//for (size_t i = 0; i < fontPolygon.rings.size(); i++)
			//{
			//	STRing &ring = fontPolygon.rings[i];
			//	for (size_t j = 0; j < ring.size(); j++)
			//	{
			//		STPoint pt1 = ring[j];
			//		STPoint pt2 = ring[(j + 1) % ring.size()];
			//		STPoint pt2Next = ring[(j + 2) % ring.size()];
			//		if ( pt1 == pt2Next )
			//		{
			//			//printf("已经处理异常情况：字形有邻边交叠在一起！\n");
			//			ring.erase( find(ring.begin(), ring.end(), pt2) );
			//			ring.erase( find(ring.begin(), ring.end(), pt2Next) );
			//		}
			//	}
			//}






			/// 生成文字底面
			{
#ifdef _TEST_DEMO
				BruceLiu::TimingName timingName("生成文字底面");
				BruceLiu::TimingBlock timingBlock(timingName);
#endif
				buildAndSubdivideUnderSuface();
			}
		}


		// 求出文字平面包围矩形
		if (!isPicCarving)
		{
			double *fontPolygonBounds = fontPolygon.getBounds();
			memcpy(wordBounds, fontPolygonBounds, sizeof(double) * 6);
			wordBounds[2] = wordBounds[3] = this->fontPolygonY;
		}

		/// 文字贴模型表面 + 筛选三角形
		{
#ifdef _TEST_DEMO
			BruceLiu::TimingName timingName("贴模型表面.");
			BruceLiu::TimingBlock timingBlock(timingName);
#endif

			vtkSmartPointer<vtkPoints> bspPoints = vtkSmartPointer<vtkPoints>::New();
			vtkSmartPointer<vtkCellArray> bspCells = vtkSmartPointer<vtkCellArray>::New();
			bspPolyData = vtkSmartPointer<vtkPolyData>::New();
			vtkSmartPointer<vtkIncrementalPointLocator> locator = NULL;
			locator.TakeReference(vtkMergePoints::New());
			locator->InitPointInsertion(bspPoints, modelPolyData->GetPoints()->GetBounds());

			vtkIdType npts, *pts;
			double A[3], B[3], C[3];
			double ABC_Min_X, ABC_Max_X, ABC_Min_Y, ABC_Max_Y, ABC_Min_Z, ABC_Max_Z;
			modelPolyData->BuildCells();
			for (vtkIdType i = 0; i < modelPolyData->GetNumberOfCells(); i++)
			{
				modelPolyData->GetCellPoints(i, npts, pts);
				modelPolyData->GetPoint(pts[0], A);
				modelPolyData->GetPoint(pts[1], B);
				modelPolyData->GetPoint(pts[2], C);
				ABC_Min_X = min( min(A[0], B[0]), C[0]);
				ABC_Max_X = max( max(A[0], B[0]), C[0]);
				//ABC_Min_Y = min( min(A[1], B[1]), C[1]);
				ABC_Max_Y = max( max(A[1], B[1]), C[1]);
				ABC_Min_Z = min( min(A[2], B[2]), C[2]);
				ABC_Max_Z = max( max(A[2], B[2]), C[2]);

				bool needAddToWorks;
				bool inBound = max(ABC_Min_X, wordBounds[0]) <= min(ABC_Max_X, wordBounds[1])
					//&& max(ABC_Min_Y, wordBounds[2]) <= min(ABC_Max_Y, wordBounds[3])
					&& max(ABC_Min_Z, wordBounds[4]) <= min(ABC_Max_Z, wordBounds[5])
					&& ABC_Max_Y > this->fontPolygonY;
				if (inBound)
				{
					//modelPolyData->GetCellPoints(i, npts, pts);
					//bspCells->InsertNextCell(3, pts);
					vtkIdType nodes[3];
					locator->InsertUniquePoint(A, nodes[0]);
					locator->InsertUniquePoint(B, nodes[1]);
					locator->InsertUniquePoint(C, nodes[2]);
					bspCells->InsertNextCell(3, nodes);
				}
				else
				{
					STTriangle triangle;
					for (vtkIdType j = 0; j < 3; j++)
					{
						double *pt = modelPolyData->GetPoint(pts[j]);
						triangle.vId[j] = AddPointReturnIndex(resultPoints, *(NewLKPoint *)pt);
					}
					resultTriangles.push_back(triangle);
				}
			}
			bspPolyData->SetPoints(bspPoints);
			bspPolyData->SetPolys(bspCells);

			// 初始化bsp树
			bspTree = vtkSmartPointer<vtkModifiedBSPTree>::New();
			bspTree->SetDataSet(bspPolyData);

			// 文字贴模型表面
			bool pastSuccessed = polygonPastToModel();
			if (!pastSuccessed)
			{
				errCode = 1;
				error = "有文字超出模型范围！";
				//assert(0);
				return false;
			}
		}


		/// 筛选出需要做拼接处理的三角面片
		{
#ifdef _TEST_DEMO
			BruceLiu::TimingName timingName("筛选三角面片.");
			BruceLiu::TimingBlock timingBlock(timingName);
#endif
			// 求出文字包围盒
			//double *wordBounds = fontPolygon.getBounds();
			//wordBounds[2] = wordBounds[3] = this->fontPolygonY;	//注意，此处很重要
			//if (wordDepth < 0.0)
			//{
			//	wordBounds[2] += wordDepth;
			//	wordBounds[3] -= wordDepth * 2;
			//}
			//else 
			//{
			//	wordBounds[2] -= wordDepth * 2;
			//	wordBounds[3] += wordDepth;
			//}

			/// 筛选要切割的三角形
			// 挑出一个最下面的三角形
			vtkIdType downmostCellId = -1;
			double intersectPtsLowerY = 1e20;
			double minPoint[3];
			minPoint[0] = wordBounds[0];
			minPoint[1] = wordBounds[2];
			minPoint[2] = wordBounds[4];
			{
				NewLKPoint point = *(NewLKPoint *)minPoint;

				double tol = 1e-10;
				double linePt1[3], linePt2[3];
				vtkPoints *intersectPoints = vtkPoints::New();
				vtkIdList *intersectCells = vtkIdList::New();
				linePt1[0] = point.v.x;
				linePt1[1] = point.v.y + 1.e5;
				linePt1[2] = point.v.z;
				linePt2[0] = point.v.x;
				linePt2[1] = point.v.y;
				linePt2[2] = point.v.z;

				//vtkSmartPointer<vtkModifiedBSPTree> underBspTree = vtkSmartPointer<vtkModifiedBSPTree>::New();
				//underBspTree->SetDataSet(modelPolyData);
				//underBspTree->IntersectWithLine(linePt1, linePt2, tol, intersectPoints, intersectCells);
				bspTree->IntersectWithLine(linePt1, linePt2, tol, intersectPoints, intersectCells);

				// 偏移
				if (intersectPoints->GetNumberOfPoints() > 0)
				{
					for (vtkIdType j = 0; j < intersectPoints->GetNumberOfPoints(); j++)
					{
						double intersectPoint[3];
						intersectPoints->GetPoint(j, intersectPoint);
						if ( intersectPtsLowerY > intersectPoint[1] )
						{
							intersectPtsLowerY = intersectPoint[1];
							downmostCellId = intersectCells->GetId(j);
						}
					}
				}
				else
				{
					// 没有交点
				}

				intersectCells->Delete();
				intersectPoints->Delete();
			}

			if (downmostCellId < 0)
			{
				errCode = 1;
				error = "有文字超出模型范围！";
				return false;
			}


			// 递归找邻边，邻边包围盒与矩形相交的，判断法向量朝向，向y轴负方向的添加
			//bool *isCellNeedAddToWorks = (bool *)calloc(modelPolyData->GetNumberOfCells(), sizeof(bool));// 三角面片是否需要添加到切割队列
			//bool *isCellCheckedNeighbers = (bool *)calloc(modelPolyData->GetNumberOfCells(), sizeof(bool));// 三角面片是否需要添加到切割队列
			bool *isCellNeedAddToWorks = (bool *)calloc(bspPolyData->GetNumberOfCells(), sizeof(bool));// 三角面片是否需要添加到切割队列
			bool *isCellCheckedNeighbers = (bool *)calloc(bspPolyData->GetNumberOfCells(), sizeof(bool));// 三角面片是否需要添加到切割队列
			std::stack<vtkIdType> waitingCheckNeighbersCells;

			isCellCheckedNeighbers[downmostCellId] = true;
			isCellNeedAddToWorks[downmostCellId] = true;
			waitingCheckNeighbersCells.push(downmostCellId);


			vtkIdType npts, *pts;
			double A[3], B[3], C[3];
			double ABC_Min_X, ABC_Max_X, ABC_Min_Y, ABC_Max_Y, ABC_Min_Z, ABC_Max_Z;
			//modelPolyData->BuildLinks();
			bspPolyData->BuildLinks();
 			while(!waitingCheckNeighbersCells.empty())
 			{
 				vtkIdType cellId = waitingCheckNeighbersCells.top();
 				waitingCheckNeighbersCells.pop();

				//vtkSmartPointer<vtkCell> cell = modelPolyData->GetCell(cellId);
				vtkSmartPointer<vtkCell> cell = bspPolyData->GetCell(cellId);
 				int edgeNum = cell->GetNumberOfEdges();
 				for (int j = 0; j < edgeNum; j++)
 				{
 					vtkSmartPointer<vtkCell> edgeCell = cell->GetEdge(j);
					vtkSmartPointer<vtkIdList> neighborIdList = vtkSmartPointer<vtkIdList>::New();
					//modelPolyData->GetCellNeighbors(cellId, edgeCell->GetPointIds(), neighborIdList);
					bspPolyData->GetCellNeighbors(cellId, edgeCell->GetPointIds(), neighborIdList);
 					for (int i = 0; i < neighborIdList->GetNumberOfIds(); i++)
 					{
 						vtkIdType nCellId = neighborIdList->GetId(i);
 
 						if (!isCellCheckedNeighbers[nCellId] && !isCellNeedAddToWorks[nCellId])
						{
							//modelPolyData->GetCellPoints(nCellId, npts, pts);
							//modelPolyData->GetPoint(pts[0], A);
							//modelPolyData->GetPoint(pts[1], B);
							//modelPolyData->GetPoint(pts[2], C);
							bspPolyData->GetCellPoints(nCellId, npts, pts);
							bspPolyData->GetPoint(pts[0], A);
							bspPolyData->GetPoint(pts[1], B);
							bspPolyData->GetPoint(pts[2], C);
 							ABC_Min_X = min( min(A[0], B[0]), C[0]);
 							ABC_Max_X = max( max(A[0], B[0]), C[0]);
 							ABC_Min_Y = min( min(A[1], B[1]), C[1]);
 							ABC_Max_Y = max( max(A[1], B[1]), C[1]);
 							ABC_Min_Z = min( min(A[2], B[2]), C[2]);
 							ABC_Max_Z = max( max(A[2], B[2]), C[2]);
 
 							bool needAddToWorks;
 							bool inBound = max(ABC_Min_X, wordBounds[0]) <= min(ABC_Max_X, wordBounds[1])
 								&& max(ABC_Min_Z, wordBounds[4]) <= min(ABC_Max_Z, wordBounds[5])
								&& ABC_Max_Y > this->fontPolygonY;
 							if (!inBound)
 							{
 								needAddToWorks = false;
 							}
 							else
 							{
								// 检测 5 文字框与模型相交
								if (ABC_Min_Y < this->fontPolygonY)
								{
									//出错
									errCode = 5;
									error = "文字框与模型相交";
									return false;
								}

 								// 检测 4 文字粘贴处过于凹凸不平
 								// 验证三角形法向量y值是否大于0
 								NewLKPoint vectorAB, vectorAC;
 								vectorAB.v.x = B[0] - A[0];
 								vectorAB.v.z = B[2] - A[2];
 								vectorAC.v.x = C[0] - A[0];
 								vectorAC.v.z = C[2] - A[2];
 								double normalY = vectorAC.v.x * vectorAB.v.z - vectorAB.v.x * vectorAC.v.z;
 								if (normalY > 0)
 								{
 									// 判断相邻的这条边是否与bounds相交
 									double pt1[3], pt2[3];
 									//modelPolyData->GetPoint(edgeCell->GetPointId(0), pt1);
 									//modelPolyData->GetPoint(edgeCell->GetPointId(1), pt2); 
 									bspPolyData->GetPoint(edgeCell->GetPointId(0), pt1);
 									bspPolyData->GetPoint(edgeCell->GetPointId(1), pt2); 
 									if (isSegmentCrossRectangeXZ(pt1, pt2, wordBounds))
 									{
 										// 出错
 										errCode = 4;
 										error = "刻文字处过于凹凸不平";
 										return false;
 									}
 									else
 									{
 										// 不要的三角形
 										needAddToWorks = false;
 									}
 
 								}
 								else
 								{
 									needAddToWorks = true;
 								}
 							}
 
 							// 将在文字框内的三角面片id添加到检测队列
 							if (needAddToWorks)
 							{
 								isCellNeedAddToWorks[nCellId] = true;
 								if (!isCellCheckedNeighbers[nCellId])
 								{
 									waitingCheckNeighbersCells.push(nCellId);
 								}
 							}
 						}
 					}
 				}
 
 				isCellCheckedNeighbers[cellId] = true;
 			}

 			// 最终筛选三角面片
 			//for (vtkIdType i = 0; i < modelPolyData->GetNumberOfCells(); i++)
 			for (vtkIdType i = 0; i < bspPolyData->GetNumberOfCells(); i++)
 			{
				//modelPolyData->GetCellPoints(i, npts, pts);
				bspPolyData->GetCellPoints(i, npts, pts);
 
 				if (!isCellNeedAddToWorks[i])
 				{
 					STTriangle triangle;
 					for (vtkIdType j = 0; j < 3; j++)
 					{
 						//double *pt = modelPolyData->GetPoint(pts[j]);
 						double *pt = bspPolyData->GetPoint(pts[j]);
 						triangle.vId[j] = AddPointReturnIndex(resultPoints, *(NewLKPoint *)pt);
 					}
 					resultTriangles.push_back(triangle);
 				}
 				else
 				{
 					// 添加到考察三角面片组中
 					STWorkTriangle checkTri;
 					for (vtkIdType j = 0; j < 3; j++)
 					{
						//double *pt = modelPolyData->GetPoint(pts[j]);
						double *pt = bspPolyData->GetPoint(pts[j]);
 						checkTri.pt[j].v.x = pt[0];
 						checkTri.pt[j].v.y = pt[1];
 						checkTri.pt[j].v.z = pt[2];
 						if (fontPolygonMesh.pointInsideXZ2D(pt)) //fontPolygon.pointInsideXZ2D(pt)
 						{
 							checkTri.positionEnum[j] = PositionInside;
 						}
 						else 
 						{
 							checkTri.positionEnum[j] = PositionOutside;
 						}
 					}
 					checkTri.calcBounds();
 					workingTriangles.push_back(checkTri);
 				}
 			}

			free(isCellNeedAddToWorks);
			free(isCellCheckedNeighbers);
 		}


		// 做切割计算
		{
#ifdef _TEST_DEMO
			BruceLiu::TimingName timingName("切割计算");
			BruceLiu::TimingBlock timingBlock(timingName);
#endif

			std::vector<STWorkTriangle> tempNewWorkTris;	// 储存新建三角面片
			std::vector<STWorkTriangle *> tempEraseTris;	// 要删除的三角面片
			std::vector<NewLKPoint> intersections;				// 与三角面片的交点
			std::vector<NewLKPoint> boundEdgeSubdividePts;				// 底面线段细分点
			for (size_t i = 0; i < fontPolygonMesh.boundEdgeNodess.size(); i++)
			{
				BoundsEdges &boundEdges = fontPolygonMesh.boundEdgeNodess[i];
				for (size_t j = 0; j < boundEdges.size(); j++)
				{
					BoundEdge *boundEdge = boundEdges[j];
					NewLKPoint &pt1 = fontPolygonMesh.points[boundEdge->ptIdFirst];
					NewLKPoint &pt2 = fontPolygonMesh.points[boundEdge->ptId];
					STWorkLine line(pt1, pt2);
					double pt1_offsetY = -wordDepth, pt2_offsetY = -wordDepth;
					if (isPicCarving)
					{
						pt1_offsetY = fontPolygonMesh.picture_pointsYOffset[boundEdge->ptIdFirst];
						pt2_offsetY = fontPolygonMesh.picture_pointsYOffset[boundEdge->ptId];
					}

					tempNewWorkTris.resize(0);
					tempEraseTris.resize(0);
					intersections.resize(0);
					boundEdgeSubdividePts.resize(0);

					for (vtkIdType n = 0; n < workingTriangles.size(); n++)
					{
						STWorkTriangle &workTriangle = workingTriangles[n];

						//printf("tri%i bounds %f, %f, %f, %f\n", n, workTriangle.minX, workTriangle.minZ, workTriangle.maxX, workTriangle.maxZ);
						// 判断线段是否与包围盒相交
						if (
							max(line.minX, workTriangle.minX) > min(line.maxX, workTriangle.maxX)
							|| max(line.minZ, workTriangle.minZ) > min(line.maxZ, workTriangle.maxZ)
							)
						{
							continue;
						}


						// 两端点和三角形位置关系
						size_t pt1OnIndex = 100, pt2OnIndex = 100;
						PointWithPolygon pt1WithTriangle = pointWithWorkTriangle(pt1, workTriangle, pt1OnIndex);
						PointWithPolygon pt2WithTriangle = pointWithWorkTriangle(pt2, workTriangle, pt2OnIndex);

						// 一、有点在三角形里
						//   1. 两个点都在三角形里
						if (pt1WithTriangle && pt2WithTriangle)
						{
							//三角化
							triangulateFivePoint(workTriangle, line, pt1, pt2, pt1WithTriangle, pt2WithTriangle, pt1OnIndex, pt2OnIndex, tempNewWorkTris, tempEraseTris, intersections);
							break;
						}

						//   2. 一个点在三角形里
						if ( !pt1WithTriangle && pt2WithTriangle || pt1WithTriangle && !pt2WithTriangle )
						{
							triangulateCross_1_2(workTriangle, line, pt1, pt2, pt1WithTriangle, pt2WithTriangle, pt1OnIndex, pt2OnIndex, tempNewWorkTris, tempEraseTris, intersections);
							continue;
						}

						// 二、没有点在三角形里
						if (!pt1WithTriangle && !pt2WithTriangle)
						{
							triangulateCross_2(workTriangle, line, pt1, pt2, pt1WithTriangle, pt2WithTriangle, pt1OnIndex, pt2OnIndex, tempNewWorkTris, tempEraseTris, intersections);
							continue;
						}

					}

					//把废弃的三角面片清除，添加新生成的三角面片
					//workingTriangles.reserve(workingTriangles.size()+100);
					std::sort(tempEraseTris.begin(), tempEraseTris.end());
					for (int idx = tempEraseTris.size() - 1; idx >= 0; idx--)
					{
						workingTriangles.erase(std::find(workingTriangles.begin(), workingTriangles.end(), *tempEraseTris[idx]));
					}
					workingTriangles.insert(workingTriangles.end(), tempNewWorkTris.begin(), tempNewWorkTris.end());


					/// 生成文字边缘竖直三角面片
					if (intersections.size() < 2)
					{
						fprintf(stderr, "线段两点重叠或过于接近\n");
					}
					else
					{
						NewLKPoint vectorPt1Pt2;
						vectorPt1Pt2.v.x = pt2.v.x - pt1.v.x;
						vectorPt1Pt2.v.z = pt2.v.z - pt1.v.z;
						// 端点加入交点表 这里不再添加，端点会有与三角形顶点合并处理的情况
						//AddPointNoRepeatReturnIndex(intersections, pt1); //把端点加入交点表
						//AddPointNoRepeatReturnIndex(intersections, pt2); //把端点加入交点表
						// 交点排序
						{
							int j, k, h;
							NewLKPoint tempPt; 
							for (h=intersections.size()-1; h>0; h=k) /*循环到没有比较范围*/
							{
								for (j=0, k=0; j<h; j++) /*每次预置k=0，循环扫描后更新k*/
								{
									if ( (intersections[j].v.x - intersections[j+1].v.x) * vectorPt1Pt2.v.x > 1.e-16
										|| ( abs((intersections[j].v.x - intersections[j+1].v.x) * vectorPt1Pt2.v.x) <= 1.e-16
										&& (intersections[j].v.z - intersections[j+1].v.z) * vectorPt1Pt2.v.z > 0 )
										) /*大的放在后面，小的放到前面*/
									{
										tempPt = intersections[j];
										intersections[j] = intersections[j+1];
										intersections[j+1] = tempPt; /*完成交换*/
										k = j; /*保存最后下沉的位置。这样k后面的都是排序排好了的。*/
									}
								}
							}
						}
						// 组成三角形
						NewLKPoint pt1UnderSurface = pt1, pt2UnderSurface = pt2;
						pt1UnderSurface.v.y += pt1_offsetY;
						pt2UnderSurface.v.y += pt2_offsetY;
						// 模型表面轮廓增加顶点的 面片群
						for (int idx = 1; idx < intersections.size(); idx++)
						{
							STTriangle newTriangle;
							newTriangle.vId[0] = AddPointReturnIndex(resultPoints, pt1UnderSurface);
							newTriangle.vId[1] = AddPointReturnIndex(resultPoints, intersections[idx-1]);
							newTriangle.vId[2] = AddPointReturnIndex(resultPoints, intersections[idx]);
							resultTriangles.push_back(newTriangle);
						}
						//// 底面半部分
						//if (!intersections.empty())
						//{
						//	STTriangle newTriangle;
						//	newTriangle.vId[0] = AddPointReturnIndex(resultPoints, pt1UnderSurface);
						//	newTriangle.vId[1] = AddPointReturnIndex(resultPoints, pt2UnderSurface);
						//	newTriangle.vId[2] = AddPointReturnIndex(resultPoints, intersections.back());
						//	resultTriangles.push_back(newTriangle);
						//}
						// 底面轮廓增加顶点的面片群
						if (isPicCarving) {
							STTriangle newTriangle;
							newTriangle.vId[0] = AddPointReturnIndex(resultPoints, intersections.back());
							newTriangle.vId[1] = AddPointReturnIndex(resultPoints, pt2UnderSurface);
							newTriangle.vId[2] = AddPointReturnIndex(resultPoints, pt1UnderSurface);
							resultTriangles.push_back(newTriangle);
						}
						else {
							boundEdgeSubdividePts.push_back(pt1UnderSurface);
							if (fontPolygonMesh.boundEdgeNodess[i][j]->pmTri) // 接收到细分的边
							{
								std::stack<BoundEdge *> waitingRecurseAddEdge;	// 等待加入的edge
								waitingRecurseAddEdge.push(fontPolygonMesh.boundEdgeNodess[i][j]);
								while(!waitingRecurseAddEdge.empty()) {
									BoundEdge *boundEdge = waitingRecurseAddEdge.top();
									waitingRecurseAddEdge.pop();
									if (boundEdge->isLeaf)
									{
										NewLKPoint newPt = fontPolygonMesh.points[boundEdge->ptId];
										newPt.v.y -= wordDepth;
										boundEdgeSubdividePts.push_back(newPt);
									}
									else
									{
										waitingRecurseAddEdge.push(boundEdge->subdivides[1]);
										waitingRecurseAddEdge.push(boundEdge->subdivides[0]);
									}
								}
							}
							else	//没有接收到细分的边，再次细分
							{
								std::stack<SubdivideLine> waitingSubdividelines;
								SubdivideLine sLine(pt1UnderSurface, pt2UnderSurface);
								sLine.calcLen2();
								waitingSubdividelines.push(sLine);
								while(!waitingSubdividelines.empty()) {
									SubdivideLine line = waitingSubdividelines.top();
									waitingSubdividelines.pop();
									if (line.length2 <= _subdivideLen2)
									{
										NewLKPoint newPt = line.pts[1];
										pastedToModelFace(newPt);
										newPt.v.y -= wordDepth;
										boundEdgeSubdividePts.push_back(newPt);
									}
									else
									{
										NewLKPoint midPt;
										midPt.v.x = (line.pts[0].v.x + line.pts[1].v.x) / 2;
										midPt.v.y = (line.pts[0].v.y + line.pts[1].v.y) / 2;
										midPt.v.z = (line.pts[0].v.z + line.pts[1].v.z) / 2;
										SubdivideLine line0(line.pts[0], midPt), line1(midPt, line.pts[1]);
										line0.length2 = line1.length2 = line.length2 / 4;
										waitingSubdividelines.push(line1);
										waitingSubdividelines.push(line0);
									}
								}
							}
							for (int idx = 1; idx < boundEdgeSubdividePts.size(); idx++)
							{
								STTriangle newTriangle;
								newTriangle.vId[0] = AddPointReturnIndex(resultPoints, intersections.back());
								newTriangle.vId[1] = AddPointReturnIndex(resultPoints, boundEdgeSubdividePts[idx]);
								newTriangle.vId[2] = AddPointReturnIndex(resultPoints, boundEdgeSubdividePts[idx-1]);
								resultTriangles.push_back(newTriangle);
							}
						}
					}

				}
			}

		}


		// workingTriangle切割结果合并到三角网中
		for (size_t i = 0; i < workingTriangles.size(); i++)
		{
			STWorkTriangle &triangle = workingTriangles[i];
  			bool needAddTriangle = false;
  			if( triangle.positionEnum[0] == PositionOutside
  				|| triangle.positionEnum[1] == PositionOutside
  				|| triangle.positionEnum[2] == PositionOutside )
  			{
  				needAddTriangle = true;
  			}
  			else if ( triangle.positionEnum[0] == PositionOnEdge
  				&& triangle.positionEnum[1] == PositionOnEdge
  				&& triangle.positionEnum[2] == PositionOnEdge )
  			{
  				double center[3];
  				center[0] = (triangle.pt[0].v.x + triangle.pt[1].v.x + triangle.pt[2].v.x)/3;
  				center[2] = (triangle.pt[0].v.z + triangle.pt[1].v.z + triangle.pt[2].v.z)/3;
  				if (!fontPolygonMesh.pointInsideXZ2D(center)) //fontPolygon
  				{
  					needAddTriangle = true;
  				}
  			}
  			else
  			{
  				//do nothing
  			}
  
  			if (needAddTriangle)
			{
				STTriangle newTriangle;
				for (vtkIdType j = 0; j < 3; j++)
				{
					newTriangle.vId[j] = AddPointReturnIndex(resultPoints, triangle.pt[j]);
				}
				resultTriangles.push_back(newTriangle);
			}
		}


		/// 生成文字底面
		buildUnderSuface();



		/// 生成polydata
		// create hash table and merge points/triangles.
		{
#ifdef _TEST_DEMO
			BruceLiu::TimingName timingName("合并顶点");
			BruceLiu::TimingBlock timingBlock(timingName);
#endif

			//对模型顶点颜色处理
			if (resultPointRGBs.size())
			{
				if (resultPointRGBs.size()/3 < resultPoints->GetNumberOfPoints())
				{
					resultPointRGBs.resize(resultPoints->GetNumberOfPoints()*3);
					//resultPointRGBs->SetNumberOfTuples(resultPoints->GetNumberOfPoints());
					//resultPointRGBs->SetTuple3(resultPoints->GetNumberOfPoints()-1, 1.0, 1.0, 1.0);
				}
				unsigned char colorNormal[3] = RGB_COLOR_NORMAL;
				for (vtkIdType i = 0; i < resultPointRGBs.size(); i=i+3)
				{
					//resultPointRGBs->GetTupleValue(i, color);
					//unsigned char *color = &resultPointRGBs[i];
					if (!(resultPointRGBs[i+0] || resultPointRGBs[i+1] || resultPointRGBs[i+2]))
					{
						//resultPointRGBs->SetTuple3(i, 255, 255, 255);
						//resultPointRGBs->SetTuple3(i, 0, 0, 0);
						resultPointRGBs[i+0] = colorNormal[0];
						resultPointRGBs[i+1] = colorNormal[1];
						resultPointRGBs[i+2] = colorNormal[2];
					}
				}
			}



			result = vtkSmartPointer<vtkPolyData>::New();
			vtkPoints *mergedPts = vtkPoints::New();
			mergedPts->Allocate(resultPoints->GetNumberOfPoints() /2);
			vtkCellArray *mergedPolys = vtkCellArray::New();
			mergedPolys->Allocate(resultTriangles.size());
			vtkSmartPointer<vtkIncrementalPointLocator> locator = NULL;
			locator.TakeReference(vtkMergePoints::New());
			locator->InitPointInsertion(mergedPts, resultPoints->GetBounds());
			//顶点颜色
			vtkUnsignedCharArray *pointRGBs = NULL;
			std::vector<unsigned char> pointRGBVector;
			if (resultPointRGBs.size())
			{
				pointRGBs = vtkUnsignedCharArray::New();
				//pointRGBs->SetNumberOfTuples(resultPoints->GetNumberOfPoints());
				pointRGBs->SetName("RGB");
				pointRGBs->SetNumberOfComponents(3);
				result->GetPointData()->SetScalars(pointRGBs);
				//result->GetPointData()->SetActiveScalars("RGB");
				pointRGBs->Delete();

			}

			double bounds[6];
			resultPoints->GetBounds(bounds);

			for (int i = resultTriangles.size()-1; i >= 0; i--)
			{
				STTriangle &triangle = resultTriangles[i];
				vtkIdType nodes[3];
				for (int i = 0; i < 3; i++)
				{
					double x[3];
					resultPoints->GetPoint(triangle.vId[i], x);
					locator->InsertUniquePoint(x, nodes[i]);
					if (pointRGBs && pointRGBVector.size()/3 == (mergedPts->GetNumberOfPoints()-1)) {
						//resultPointRGBs->GetTupleValue(triangle.vId[i], color);
						//pointRGBs->SetNumberOfTuples(nodes[i]+1);
						//pointRGBs->SetTupleValue(nodes[i], &resultPointRGBs[triangle.vId[i]*3]);
						pointRGBVector.push_back(resultPointRGBs[triangle.vId[i]*3]);
						pointRGBVector.push_back(resultPointRGBs[triangle.vId[i]*3+1]);
						pointRGBVector.push_back(resultPointRGBs[triangle.vId[i]*3+2]);
					}
				}

				if (nodes[0] != nodes[1] &&
					nodes[0] != nodes[2] &&
					nodes[1] != nodes[2])
				{
					mergedPolys->InsertNextCell(3, nodes);
				}
			}

			//颜色
			if (pointRGBs)
			{
				pointRGBs->SetNumberOfTuples(pointRGBVector.size()/3);
				for (int i = 0; i < pointRGBs->GetNumberOfTuples(); i++)
				{
					pointRGBs->SetTupleValue(i, &pointRGBVector[i*3]);
				}
			}


			result->SetPoints(mergedPts);
			mergedPts->Delete();
			result->SetPolys(mergedPolys);
			mergedPolys->Delete();
		}

		// 把模型旋转还原回原状
		if (modelPolyDataOrigin)
		{
			vtkSmartPointer<vtkTransformPolyDataFilter> tpdFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
			allTransformToVerticality->Inverse();
			tpdFilter->SetTransform( allTransformToVerticality );
			tpdFilter->SetInputData( result );
			tpdFilter->Update();
			result = tpdFilter->GetOutput();
		}


#ifdef _TEST_DEMO
		//BruceLiu::timerPtr->print();
#endif


		return true;
	}




	bool PolygonCarve::polygonPastToModel()
	{
		bool reVal = true;

		const double tol = 1e-10;

		// 三角形贴表面 和 偏移到底面（顶面）
		for (size_t i = 0; i < fontPolygonMesh.points.size(); i++)
		{
			NewLKPoint &point = fontPolygonMesh.points[i];

			double linePt1[3], linePt2[3];
			vtkPoints *intersectPoints = vtkPoints::New();
			vtkIdList *intersectCells = vtkIdList::New();
			linePt1[0] = point.v.x;
			linePt1[1] = point.v.y + 1.e5;
			linePt1[2] = point.v.z;
			linePt2[0] = point.v.x;
			linePt2[1] = point.v.y;
			linePt2[2] = point.v.z;

			bspTree->IntersectWithLine(linePt1, linePt2, tol, intersectPoints, intersectCells);

			// 偏移
			if (intersectPoints->GetNumberOfPoints() > 0)
			{
				double intersectPtsLowerY = 1e30;
				for (vtkIdType j = 0; j < intersectPoints->GetNumberOfPoints(); j++)
				{
					double intersectPoint[3];
					intersectPoints->GetPoint(j, intersectPoint);
					if ( intersectPtsLowerY > intersectPoint[1] )
					{
						intersectPtsLowerY = intersectPoint[1];
					}
				}
				point.v.y = intersectPtsLowerY;
				// 偏移
				//point.y  -= wordDepth;
				// 重合点偏移
				if (!isPicCarving && fontPolygonMesh.nDeltaYSubdividePts[i] != 0)
				{
					point.v.y += UNITSubdivideDoublelicationDelta * fontPolygonMesh.nDeltaYSubdividePts[i];
				}
			}
			else
			{
				reVal = false;
			}

			intersectCells->Delete();
			intersectPoints->Delete();
		}

		return reVal;
	}

	bool PolygonCarve::pastedToModelFace(NewLKPoint &point)
	{
		bool reVal = true;

		const double tol = 1e-10;

		double linePt1[3], linePt2[3];
		vtkPoints *intersectPoints = vtkPoints::New();
		vtkIdList *intersectCells = vtkIdList::New();
		linePt1[0] = point.v.x;
		linePt1[1] = point.v.y + 1000000.0;
		linePt1[2] = point.v.z;
		linePt2[0] = point.v.x;
		linePt2[1] = point.v.y - 1000000.0;
		linePt2[2] = point.v.z;

		bspTree->IntersectWithLine(linePt1, linePt2, tol, intersectPoints, intersectCells);

		// 偏移
		if (intersectPoints->GetNumberOfPoints() > 0)
		{
			double intersectPtsLowerY = 1e20;
			for (vtkIdType j = 0; j < intersectPoints->GetNumberOfPoints(); j++)
			{
				double intersectPoint[3];
				intersectPoints->GetPoint(j, intersectPoint);
				if ( intersectPtsLowerY > intersectPoint[1] )
				{
					intersectPtsLowerY = intersectPoint[1];
				}
			}
			point.v.y = intersectPtsLowerY;
		}
		else
		{
			reVal = false;
		}

		intersectCells->Delete();
		intersectPoints->Delete();

		return reVal;
	}


	// 细分三角形添加到结果
	void PolygonCarve::addSubdivideTriToResult(PMTriangle *pmTri)
	{
		if (pmTri->isLeaf)
		{
			underCellsNum++;
			STTriangle newTriangle;
			if (fontPolygonMesh.pointIsGrays.size() == 0)
			{
				newTriangle.vId[0] = AddPointReturnIndex(resultPoints, fontPolygonMesh.points[pmTri->ids[0]]);
				newTriangle.vId[1] = AddPointReturnIndex(resultPoints, fontPolygonMesh.points[pmTri->ids[1]]);
				newTriangle.vId[2] = AddPointReturnIndex(resultPoints, fontPolygonMesh.points[pmTri->ids[2]]);
			}
			else
			{
				newTriangle.vId[0] = AddPointReturnIndex(resultPoints, fontPolygonMesh.points[pmTri->ids[0]], fontPolygonMesh.pointIsGrays[pmTri->ids[0]], &resultPointRGBs);
				newTriangle.vId[1] = AddPointReturnIndex(resultPoints, fontPolygonMesh.points[pmTri->ids[1]], fontPolygonMesh.pointIsGrays[pmTri->ids[1]], &resultPointRGBs);
				newTriangle.vId[2] = AddPointReturnIndex(resultPoints, fontPolygonMesh.points[pmTri->ids[2]], fontPolygonMesh.pointIsGrays[pmTri->ids[2]], &resultPointRGBs);
			}
			resultTriangles.push_back(newTriangle);
		}
		else
		{
			addSubdivideTriToResult(pmTri->subdivides[0]);
			addSubdivideTriToResult(pmTri->subdivides[1]);
		}
	}

	void PolygonCarve::buildUnderSuface()
	{
		if (!isPicCarving)
		{
			for (int i = 0; i < fontPolygonMesh.points.size(); i++)
			{
				fontPolygonMesh.points[i].v.y -= wordDepth;
			}
		}
		else
		{
			for (int i = 0; i < fontPolygonMesh.points.size(); i++)
			{
				fontPolygonMesh.points[i].v.y += fontPolygonMesh.picture_pointsYOffset[i];
			}
		}

		//递归添加三角形
		underCellsNum = 0;
		for (int i = 0; i < fontPolygonMesh.cellNodes.size(); i++)
		{
			addSubdivideTriToResult(fontPolygonMesh.cellNodes[i]);
		}
	}

	void PolygonCarve::subdivideTriangle(PMTriangle *pmTri)
	{
		if ( pmTri->edgeLen2s[pmTri->longestEdgeId] > _subdivideLen2)
		{
			pmTri->isLeaf = false;

			NewLKPoint newPt;
			vtkIdType newPtRealIndex;
			vtkIdType newPtFirstOneIndex;
			int idB = pmTri->longestEdgeId;
			int idA = (idB + 2) % 3; //AB为最长边
			int idC = (idB + 1) % 3;
			NewLKPoint ptA = fontPolygonMesh.points[pmTri->ids[idA]]; //注意数值对应的地址会改变
			NewLKPoint ptB = fontPolygonMesh.points[pmTri->ids[idB]];
			NewLKPoint ptC = fontPolygonMesh.points[pmTri->ids[idC]];

			double len = (ptA.v.x - ptB.v.x) * (ptA.v.x - ptB.v.x) + (ptA.v.z - ptB.v.z) * (ptA.v.z - ptB.v.z);
			double lenBC = (ptC.v.x - ptB.v.x) * (ptC.v.x - ptB.v.x) + (ptC.v.z - ptB.v.z) * (ptC.v.z - ptB.v.z);
			double lenCA = (ptA.v.x - ptC.v.x) * (ptA.v.x - ptC.v.x) + (ptA.v.z - ptC.v.z) * (ptA.v.z - ptC.v.z);
			//if (abs(pmTri->edgeLen2s[pmTri->longestEdgeId] - len) > 1e-4
			//	|| lenBC > len
			//	|| lenCA > len
			//	)
			//{
			//	printf("%i	%i,%i,%i %f \n", pmTri->level, pmTri->ids[idA], pmTri->ids[idB], pmTri->ids[idC], pmTri->edgeLen2s[pmTri->longestEdgeId]);
			//}
			newPt.v.x = ptA.v.x > ptB.v.x ? (ptA.v.x + ptB.v.x) / 2 : (ptB.v.x + ptA.v.x) / 2;
			newPt.v.y = ptA.v.y > ptB.v.y ? (ptA.v.y + ptB.v.y) / 2 : (ptB.v.y + ptA.v.y) / 2;
			newPt.v.z = ptA.v.z > ptB.v.z ? (ptA.v.z + ptB.v.z) / 2 : (ptB.v.z + ptA.v.z) / 2;
			double a = ptC.v.x > newPt.v.x ? (ptC.v.x - newPt.v.x) : (newPt.v.x - ptC.v.x);
			double b = ptC.v.z > newPt.v.z ? (ptC.v.z - newPt.v.z) : (newPt.v.z - ptC.v.z);
			double newLen = a * a + b * b;
			//注意数值对应的地址会改变
			//fontPolygonMesh.points.push_back(newPt);
			//vtkIdType ptsBeforeSize = subdivideMergePts->GetNumberOfPoints();
			fontPolygonMesh.subdivideLocator->InsertUniquePoint((double*)&newPt, newPtFirstOneIndex);
			newPtRealIndex = newPtFirstOneIndex;
			bool isRealNewPt = newPtFirstOneIndex == fontPolygonMesh.points.size();
			//颜色处理
			if (isRealNewPt && fontPolygonMesh.pointIsGrays.size()>0)
			{
				char isGrayA = fontPolygonMesh.pointIsGrays[pmTri->ids[idA]]; //注意数值对应的地址会改变
				char isGrayB = fontPolygonMesh.pointIsGrays[pmTri->ids[idB]];
				char isGrayC = fontPolygonMesh.pointIsGrays[pmTri->ids[idC]];
				char newPtIsGray = isGrayA && isGrayB;
				fontPolygonMesh.pointIsGrays.push_back(newPtIsGray);
			}
			// 该点父线段
			Array<vtkIdType, 2> newPtParentSegment;
			newPtParentSegment.value[0] = min(pmTri->ids[idA], pmTri->ids[idB]);
			newPtParentSegment.value[1] = max(pmTri->ids[idA], pmTri->ids[idB]);
			bool isSamePtButParent = false;
			if (!isRealNewPt)
			{
				if (!(newPtParentSegment == fontPolygonMesh.subdividePtParentSegments[newPtRealIndex]))
				{
					if (fontPolygonMesh.doublelicationMap[newPtFirstOneIndex] == -1 ) // 扩展为多点共位置时更改此处
					{
						fontPolygonMesh.subdivideLocator->InsertNextPoint((double *)&newPt);
						newPtRealIndex = fontPolygonMesh.subdivideMergePts->GetNumberOfPoints() - 1;
						fontPolygonMesh.doublelicationMap[newPtFirstOneIndex] = newPtRealIndex;
						isSamePtButParent = true;
						isRealNewPt = true;
					}
					else
					{
						newPtRealIndex = fontPolygonMesh.doublelicationMap[newPtFirstOneIndex];
						assert(fontPolygonMesh.doublelicationMap[newPtRealIndex] == newPtFirstOneIndex);// 多点共位置会有问题
					}
				}
			}
			if (isRealNewPt)
			{
				fontPolygonMesh.points.push_back(newPt);
				fontPolygonMesh.subdividePtParentSegments.push_back(newPtParentSegment);
				if (isSamePtButParent)
				{
					// 偏移量处理
					int newNDeltaY = 0;
					int &nDeltaY_A = fontPolygonMesh.nDeltaYSubdividePts[pmTri->ids[idA]];
					int &nDeltaY_B = fontPolygonMesh.nDeltaYSubdividePts[pmTri->ids[idB]];
					int minDelta_A_B = min(nDeltaY_A, nDeltaY_B);
					newNDeltaY += minDelta_A_B == 0 ? 1 : minDelta_A_B;
					newNDeltaY += max(nDeltaY_A, nDeltaY_B);
					//this->nDeltaYSubdividePts.push_back(0);
					this->fontPolygonMesh.nDeltaYSubdividePts.push_back(newNDeltaY);
					//this->nDeltaYSubdividePts[newPtFirstOneIndex] = -newNDeltaY;
					fontPolygonMesh.doublelicationMap.push_back(newPtFirstOneIndex);
				}
				else
				{
					this->fontPolygonMesh.nDeltaYSubdividePts.push_back(0);
					fontPolygonMesh.doublelicationMap.push_back(-1);
				}
			}



			{
				PMTriangle *newTri = new PMTriangle();
				newTri->level = pmTri->level + 1;
				pmTri->subdivides[0] = newTri;
				newTri->isLeaf = true;
				newTri->ids[0] = pmTri->ids[idA];
				newTri->edgeLen2s[0] = pmTri->edgeLen2s[idA];
				newTri->ids[1] = newPtRealIndex;
				newTri->edgeLen2s[1] = pmTri->edgeLen2s[idB] / 4;
				newTri->ids[2] = pmTri->ids[idC];
				newTri->edgeLen2s[2] = newLen;
				/// 边界线处理
				newTri->boundEdges[0] = pmTri->boundEdges[idA]; 
				newTri->boundEdges[1] = (BoundEdge *)NULL;
				newTri->boundEdges[2] = (BoundEdge *)NULL;
				if (pmTri->boundEdges[idA])
				{
					pmTri->boundEdges[idA]->pmTri = newTri;
				}
				// 细分边界线
				if (pmTri->boundEdges[idB])
				{
					BoundEdge *boundEdge = pmTri->boundEdges[idB];
					boundEdge->isLeaf = false;
					BoundEdge *newBoundEdge = new BoundEdge();
					newBoundEdge->isLeaf = true;
					newBoundEdge->pmTri = newTri;
					newBoundEdge->sameClockwiseAsTriangle = boundEdge->sameClockwiseAsTriangle;
					newBoundEdge->ptId = boundEdge->sameClockwiseAsTriangle ? newPtRealIndex : pmTri->ids[idA];
					newBoundEdge->edgeIndex = 1;
					//边树枝节的子节点赋值
					boundEdge->subdivides[boundEdge->sameClockwiseAsTriangle?0:1] = newBoundEdge;
					//边界赋予三角形属性
					newTri->boundEdges[1] = newBoundEdge;
				}

				double longestEdgeLen2 = -1;
				for (int i = 0; i < 3; i++)
				{
					if (newTri->edgeLen2s[i] > longestEdgeLen2)
					{
						longestEdgeLen2 = newTri->edgeLen2s[i];
						newTri->longestEdgeId = i;
					}
				}
				//递归细分三角形
				if (longestEdgeLen2 > _subdivideLen2)
				{
					//this->subdivideTriangle(newTri);
					fontPolygonMesh.subdividingTri.push(newTri);
				}
			}
			{
				PMTriangle *newTri = new PMTriangle();
				newTri->level = pmTri->level + 1;
				pmTri->subdivides[1] = newTri;
				newTri->isLeaf = true;
				newTri->ids[0] = newPtRealIndex;
				newTri->edgeLen2s[0] = newLen;
				newTri->ids[1] = pmTri->ids[idB];
				newTri->edgeLen2s[1] = pmTri->edgeLen2s[idB] / 4;
				newTri->ids[2] = pmTri->ids[idC];
				newTri->edgeLen2s[2] = pmTri->edgeLen2s[idC];
				/// 边界线处理
				newTri->boundEdges[0] = (BoundEdge *)NULL; 
				newTri->boundEdges[1] = (BoundEdge *)NULL; 
				newTri->boundEdges[2] = pmTri->boundEdges[idC];
				if (pmTri->boundEdges[idC])
				{
					pmTri->boundEdges[idC]->pmTri = newTri;
				}
				// 细分边界线
				if (pmTri->boundEdges[idB])
				{
					BoundEdge *boundEdge = pmTri->boundEdges[idB];
					boundEdge->isLeaf = false;
					BoundEdge *newBoundEdge = new BoundEdge();
					newBoundEdge->isLeaf = true;
					newBoundEdge->pmTri = newTri;
					newBoundEdge->sameClockwiseAsTriangle = boundEdge->sameClockwiseAsTriangle;
					newBoundEdge->ptId = boundEdge->sameClockwiseAsTriangle ? pmTri->ids[idB] : newPtRealIndex;
					newBoundEdge->edgeIndex = 1;
					//边树枝节?淖咏诘愀持?
					boundEdge->subdivides[boundEdge->sameClockwiseAsTriangle?1:0] = newBoundEdge;
					//边界赋予三角形属性
					newTri->boundEdges[1] = newBoundEdge;
				}

				double longestEdgeLen2 = -1;
				for (int i = 0; i < 3; i++)
				{
					if (newTri->edgeLen2s[i] > longestEdgeLen2)
					{
						longestEdgeLen2 = newTri->edgeLen2s[i];
						newTri->longestEdgeId = i;
					}
				}
				//递归细分三角形
				if (longestEdgeLen2 > _subdivideLen2)
				{
					//this->subdivideTriangle(newTri);
					fontPolygonMesh.subdividingTri.push(newTri);
				}
			}
		}
	}

	void PolygonCarve::addUnderTriangle(NewLKPoint pt0, NewLKPoint pt1, NewLKPoint pt2)
	{
		vtkIdType nodes[3];
		fontPolygonMesh.fb_locator->InsertUniquePoint((double *)&pt0, nodes[0]);
		fontPolygonMesh.fb_locator->InsertUniquePoint((double *)&pt1, nodes[1]);
		fontPolygonMesh.fb_locator->InsertUniquePoint((double *)&pt2, nodes[2]);
		if (nodes[0] != nodes[1] &&
			nodes[0] != nodes[2] &&
			nodes[1] != nodes[2])
		{
			fontPolygonMesh.fb_fontPolydata->GetPolys()->InsertNextCell(3, nodes);
		}
		
	}


	void PolygonMesh::findBoundEdgess()
	{
		vtkIdType ptNum = fb_fontPolydata->GetNumberOfPoints();
		// 初始化 点-选段
		fb_ptid_segmentArrs.reserve(ptNum * 1.01);
		fb_ptid_segmentArrs.resize(ptNum);
		for (int i = 0; i < ptNum; i++)
		{
			fb_ptid_segmentArrs[i].reserve(2);
		}

		fb_fontPolydata->BuildLinks();
		vtkSmartPointer<vtkIdList> neighborIdList = vtkSmartPointer<vtkIdList>::New();
		vtkIdType npts, *pts;
		for (int cellId = 0; cellId < fb_fontPolydata->GetNumberOfCells(); cellId++)
		{
			fb_fontPolydata->GetCellPoints(cellId, npts, pts);
			for (int cur = 0, prev = npts - 1; cur < npts; prev = cur, cur++)
			{
				neighborIdList->Reset();
				fb_fontPolydata->GetCellEdgeNeighbors(cellId, pts[prev], pts[cur], neighborIdList);
				if (neighborIdList->GetNumberOfIds() == 0)
				{
					FBSegment *fbSegment = new FBSegment();
					fbSegment->added = false;
					fbSegment->ptids[0] = pts[prev];
					fbSegment->ptids[1] = pts[cur];
					fb_ptid_segmentArrs[pts[prev]].push_back(fbSegment);
					fb_ptid_segmentArrs[pts[cur]].push_back(fbSegment);
				}
			}
		}

		// 连接边界线段
		int findFromPtid = 0;
		for (; findFromPtid < fb_ptid_segmentArrs.size(); findFromPtid++)
		{
			std::vector<FBSegment *> &fbSegments = fb_ptid_segmentArrs[findFromPtid];
			for (int i = 0; i < fbSegments.size(); i++)
			{
				if (fbSegments[i]->added == false)
				{
					findBoundEdgeRing(fbSegments[i]);
				}
			}
		}
	}

	void PolygonMesh::findBoundEdgeRing(FBSegment *firstSegment)
	{
		unsigned short firstId = boundEdgeNodess.size();
		unsigned short secondeId = 0;
		FBSegment *curSegment = firstSegment;
		vtkIdType curPtid = firstSegment->ptids[0];
		vtkIdType nextPtid = firstSegment->ptids[1];
		bool curSameAsTriangle = curPtid == curSegment->ptids[0];
		PtRingIds ptRingIds;
		ptRingIds.firstId = firstId;
		// 先给第一个点赋值顺序值
		ptRingIds.secondId = -1;	//注意，如果不闭合为-1
		ptRingIds.ptId = curPtid;
		points[curPtid].v.y = *(double *)&ptRingIds;
		//增加ring
		boundEdgeNodess.resize(boundEdgeNodess.size()+1);
		BoundsEdges &boudEdges = boundEdgeNodess.back();
		do 
		{
			nextPtid = curSegment->ptids[curSameAsTriangle ? 1 : 0];


			// 点赋予顺序值
			//ptRingIds.firstId = firstId;
			ptRingIds.secondId = secondeId;
			ptRingIds.ptId = nextPtid;
			points[nextPtid].v.y = *(double *)&ptRingIds;


			BoundEdge *newBoundEdge = new BoundEdge;
			newBoundEdge->isLeaf = true;
			newBoundEdge->ptIdFirst = curPtid;
			newBoundEdge->ptId = nextPtid;
			newBoundEdge->sameClockwiseAsTriangle = curSameAsTriangle;
			newBoundEdge->pmTri = NULL;
			newBoundEdge->edgeIndex = -1;
			boudEdges.push_back(newBoundEdge);
			
			curSegment->added = true;
			curSegment = NULL;


			///// 重复的点添加新顶点  需要改动较大：a）连轮廓线时要    b）三角形顶点索引值
			//int addedNum = 0;
			//std::vector<FBSegment *> &curSegmentArr = fb_ptid_segmentArrs[curPtid];
			//if (curSegmentArr.size() > 2)
			//{
			//	for (int i = 0; i < curSegmentArr.size(); i++)
			//	{
			//		if (curSegmentArr[i]->added)
			//		{
			//			addedNum++;
			//		}
			//	}

			//	// 添加新点，偏移量+1
			//	if (addedNum >= 2)
			//	{
			//		STPoint newPt = points[curPtid];
			//		// 添加新点
			//		vtkIdType newPtid = points.size();
			//		points.push_back(newPt);
			//		newPt.y = 0;
			//		subdivideLocator->InsertNextPoint((double *)&newPt);
			//		//偏移量+1
			//		nDeltaYSubdividePts.push_back(nDeltaYSubdividePts[curPtid]+1);
			//		// 共点映射
			//		doublelicationMap.push_back(curPtid);
			//		doublelicationMap[curPtid] = newPtid;
			//		// 父线段赋值
			//		Array<vtkIdType, 2> pSegment;
			//		pSegment.value[0] = pSegment.value[1] = -1;
			//		subdividePtParentSegments.push_back(pSegment);

			//		//点-线段 表更新
			//		fb_ptid_segmentArrs.push_back(std::vector<FBSegment *>(0));
			//		std::vector<FBSegment *> &newSegmentArr = fb_ptid_segmentArrs.back();
			//		for (int i = 0; i < curSegmentArr.size(); i++)
			//		{
			//			if (!curSegmentArr[i]->added)
			//			{
			//				FBSegment *moveFBSegment = curSegmentArr[i];
			//				for (int j = 0; j < 2; j++)
			//				{
			//					if( curPtid == moveFBSegment->ptids[j] )
			//					{
			//						moveFBSegment->ptids[j] = newPtid;
			//					}
			//				}
			//				newSegmentArr.push_back(moveFBSegment);
			//				curSegmentArr.erase(find(curSegmentArr.begin(), curSegmentArr.end(), moveFBSegment));
			//				i--;
			//			}
			//		}
			//	}
			//}





			curSameAsTriangle = false;
			std::vector<FBSegment *> &nextSegmentArrs = fb_ptid_segmentArrs[nextPtid];
			for (int i = 0; i < nextSegmentArrs.size(); i++)
			{
				if (!nextSegmentArrs[i]->added)
				{
					curSegment = nextSegmentArrs[i];
					if (curSegment->ptids[0] == nextPtid)
					{
						curSameAsTriangle = true;
						break;
					}
				}
			}

			curPtid = nextPtid;
			secondeId++;
		} while (curSegment);
	}


	void PolygonMesh::cleanNarrowTriangles(double minAngleDegree)
	{
		double minAngleRadian = vtkMath::RadiansFromDegrees(minAngleDegree);
		double minCos = cos(minAngleRadian);

		double trianglePts[3][3];
		double vector1[3], vector2[3];
		bool isSmallAngle[3];
		int smallAngleNum, biggerAngleIndex, biggerAngleDot;
		int fromCellId = 0;
		bool *deletedCell = (bool *)malloc(fb_fontPolydata->GetNumberOfCells() * sizeof(bool));
		while (fromCellId < fb_fontPolydata->GetNumberOfCells())
		{
LoopFilter:
			memset(deletedCell, 0, fb_fontPolydata->GetNumberOfCells());
			vtkSmartPointer<vtkCellArray> newCells = vtkSmartPointer<vtkCellArray>::New();
			fb_fontPolydata->BuildCells();
			vtkIdType npts, *pts;
			for (int cellId = fromCellId; cellId < fb_fontPolydata->GetNumberOfCells(); cellId++)
			{
				fb_fontPolydata->GetCellPoints(cellId, npts, pts);
				for (int i = 0; i < 3; i++)
				{
					fb_fontPolydata->GetPoint(pts[i], trianglePts[i]);
				}
				// 求三个角
				smallAngleNum = 0;
				biggerAngleDot = 2;
				for (int i = 0; i < 3; i++)
				{
					int iPrev = (i+2)%3;
					int iNext = (i+1)%3;
					vector1[0] = trianglePts[iPrev][0] - trianglePts[i][0];
					vector1[1] = 0;
					vector1[2] = trianglePts[iPrev][2] - trianglePts[i][2];
					vtkMath::Normalize(vector1);
					vector2[0] = trianglePts[iNext][0] - trianglePts[i][0];
					vector2[1] = 0;
					vector2[2] = trianglePts[iNext][2] - trianglePts[i][2];
					vtkMath::Normalize(vector2);

					double dot = vector1[0] * vector2[0] + vector1[2] * vector2[2];
					isSmallAngle[i] = dot > minCos;
					if (isSmallAngle[i])
					{
						smallAngleNum++;
					}
					else
					{
						biggerAngleIndex = i;
					}
					//if (biggerAngleDot > dot)
					//{
					//	biggerAngleDot = dot;
					//	biggerAngleIndex = i;
					//}
				}

				assert (smallAngleNum != 3);
				if (smallAngleNum > 1)
				{
					deletedCell[cellId] = true;// 删除当前cell
					
					int idA = biggerAngleIndex,
						idB = (biggerAngleIndex+1)%3,
						idC = (biggerAngleIndex+2)%3;// BC边为要切割的边
					int ptIdA = pts[idA];
 					//STWorkLine lineBC(*(STPoint *)trianglePts[idB], *(STPoint *)trianglePts[idC]);
 					//double a, b, c;
 					//a = lineBC.b;
 					//b = -lineBC.a;
 					//c = -( a * trianglePts[idA][0] + b * trianglePts[idA][2] );
 					//double newPtA[3];
 					//newPtA[1] = trianglePts[idA][1];
 					//newPtA[0] = - lineBC.a * lineBC.c - a * c;
 					//newPtA[2] = - lineBC.b * lineBC.c - b * c;
 					//fb_fontPolydata->GetPoints()->InsertPoint(ptIdA, newPtA);
 					//double ptPie[3];
 					//fb_fontPolydata->GetPoint(ptIdA, ptPie);

					// 邻边
					vtkSmartPointer<vtkIdList> neighborIdList = vtkSmartPointer<vtkIdList>::New();
					vtkSmartPointer<vtkIdList> edge = vtkSmartPointer<vtkIdList>::New();
					int ptIdB = pts[idB],
						ptIdC = pts[idC];
					edge->InsertNextId(ptIdB);
					edge->InsertNextId(ptIdC);
					fb_fontPolydata->GetCellNeighbors(cellId, edge, neighborIdList);
					if (neighborIdList->GetNumberOfIds() > 0)
					{
						assert(neighborIdList->GetNumberOfIds() == 1);

						int neighborCellId = neighborIdList->GetId(0);

						deletedCell[neighborCellId] = true;// 删除此相邻cell

						// 更新polydata进入下一轮检测
						int fb_fontPolydataCellNum = fb_fontPolydata->GetNumberOfCells();
						for (int cellId2 = 0; cellId2 < fb_fontPolydataCellNum; cellId2++)
						{
							if (!deletedCell[cellId2])
							{
								fb_fontPolydata->GetCellPoints(cellId2, npts, pts);
								newCells->InsertNextCell(3, pts);
							}
						}

						/// 插入两个新三角形
						vtkIdType neighborNpts, *neighborPts;
						fb_fontPolydata->GetCellPoints(neighborCellId, neighborNpts, neighborPts);
						int nIdA = -1, nIdB, nIdC; // AB为邻边
						for (int i = 0, iPrev = 2; i < 3; iPrev = i, i++)
						{
							if (min(ptIdB, ptIdC) == min(neighborPts[iPrev], neighborPts[i]) &&
								max(ptIdB, ptIdC) == max(neighborPts[iPrev], neighborPts[i])
								)
							{
								nIdA = iPrev;
								nIdB = i;
								nIdC = (i + 1) % 3;
								break;
							}
						}
						// 插入新三角形
						{
							vtkIdType newPts[3];
							newPts[0] = neighborPts[nIdA];
							newPts[1] = ptIdA;
							newPts[2] = neighborPts[nIdC];
							newCells->InsertNextCell(3, newPts);
						}
						{
							vtkIdType newPts[3];
							newPts[0] = ptIdA;
							newPts[1] = neighborPts[nIdB];
							newPts[2] = neighborPts[nIdC];
							newCells->InsertNextCell(3, newPts);
						}

						vtkSmartPointer<vtkPolyData> new_fb_fontPolydata = vtkSmartPointer<vtkPolyData>::New();
						new_fb_fontPolydata->SetPoints(fb_fontPolydata->GetPoints());
						new_fb_fontPolydata->SetPolys(newCells);
						fb_fontPolydata = new_fb_fontPolydata;
						//fb_fontPolydata->SetPolys(newCells);
						fromCellId = cellId;
						if (neighborCellId < cellId)
						{
							fromCellId--;
						}
						goto LoopFilter; // 注意，使用了goto
						//break;
					}
					else
					{
						//printf("no neighbor \n");

						// 更新polydata进入下一轮检测
						int fb_fontPolydataCellNum = fb_fontPolydata->GetNumberOfCells();
						for (int cellId2 = 0; cellId2 < fb_fontPolydataCellNum; cellId2++)
						{
							if (!deletedCell[cellId2])
							{
								fb_fontPolydata->GetCellPoints(cellId2, npts, pts);
								newCells->InsertNextCell(3, pts);
							}
						}
						vtkSmartPointer<vtkPolyData> new_fb_fontPolydata = vtkSmartPointer<vtkPolyData>::New();
						new_fb_fontPolydata->SetPoints(fb_fontPolydata->GetPoints());
						new_fb_fontPolydata->SetPolys(newCells);
						fb_fontPolydata = new_fb_fontPolydata;
						//fb_fontPolydata->SetPolys(newCells);
						fromCellId = cellId;
						goto LoopFilter; // 注意，使用了goto
						//break;
					}
				}
				else
				{
					// do nothing
					//newCells->InsertNextCell(3, pts);
				}
			}

			fromCellId++;
		}

		free(deletedCell);
	}






	////////TODO: need delete begin///////////////
#define FTPoint NewLKPoint
#ifndef FTGL_DOUBLE
#define FTGL_DOUBLE double
#endif
#ifndef CALLBACK
#define CALLBACK
#endif

#if defined __APPLE_CC__ && __APPLE_CC__ < 5465
	typedef GLvoid (*GLUTesselatorFunction) (...);
#elif defined WIN32 && !defined __CYGWIN__
	typedef GLvoid (CALLBACK *GLUTesselatorFunction) ();
#else
	typedef GLvoid (*GLUTesselatorFunction) ();
#endif




	
	/**
	 * FTTesselation captures points that are output by OpenGL's gluTesselator.
	 */
	class FTTesselation
	{
		public:
			/**
			 * Default constructor
			 */
			FTTesselation(GLenum m)
			: meshType(m)
			{
				pointList.reserve(128);
			}

			/**
			 *  Destructor
			 */
			~FTTesselation()
			{
				pointList.clear();
			}

			/**
			 * Add a point to the mesh.
			 */
			void AddPoint(const FTGL_DOUBLE x, const FTGL_DOUBLE y,
						  const FTGL_DOUBLE z)
			{
				pointList.push_back(NewLKPoint(x, y, z));
			}

			/**
			 * The number of points in this mesh
			 */
			size_t PointCount() const { return pointList.size(); }

			/**
			 *
			 */
			const NewLKPoint& Point(unsigned int index) const
			{ return pointList[index]; }

			/**
			 * Return the OpenGL polygon type.
			 */
			GLenum PolygonType() const { return meshType; }

		private:
			/**
			 * Points generated by gluTesselator.
			 */
			typedef std::vector<NewLKPoint> PointVector;
			PointVector pointList;

			/**
			 * OpenGL primitive type from gluTesselator.
			 */
			GLenum meshType;
	};


	/**
	 * FTMesh is a container of FTTesselation's that make up a polygon glyph
	 */
	class FTMesh
	{
			typedef std::vector<FTTesselation*> TesselationVector;
			typedef std::list<NewLKPoint> PointList;

		public:
			/**
			 * Default constructor
			 */
			FTMesh();

			/**
			 *  Destructor
			 */
			~FTMesh();

			/**
			 * Add a point to the mesh
			 */
			void AddPoint(const FTGL_DOUBLE x, const FTGL_DOUBLE y,
						  const FTGL_DOUBLE z);

			/**
			 *  Create a combine point for the gluTesselator
			 */
			const FTGL_DOUBLE* Combine(const FTGL_DOUBLE x, const FTGL_DOUBLE y,
									   const FTGL_DOUBLE z);

			/**
			 * Begin a new polygon
			 */
			void Begin(GLenum meshType);

			/**
			 * End a polygon
			 */
			void End();

			/**
			 * Record a gluTesselation error
			 */
			void Error(GLenum e) { err = e; }

			/**
			 * The number of tesselations in the mesh
			 */
			size_t TesselationCount() const { return tesselationList.size(); }

			/**
			 * Get a tesselation by index
			 */
			const FTTesselation* const Tesselation(size_t index) const;

			/**
			 * Return the temporary point list. For testing only.
			 */
			const PointList& TempPointList() const { return tempPointList; }

			/**
			 * Get the GL ERROR returned by the glu tesselator
			 */
			GLenum Error() const { return err; }

		private:
			/**
			 * The current sub mesh that we are constructing.
			 */
			FTTesselation* currentTesselation;

			/**
			 * Holds each sub mesh that comprises this glyph.
			 */
			TesselationVector tesselationList;

			/**
			 * Holds extra points created by gluTesselator. See ftglCombine.
			 */
			PointList tempPointList;

			/**
			 * GL ERROR returned by the glu tesselator
			 */
			GLenum err;

	};

 	FTMesh::FTMesh()
		: currentTesselation(0),
		err(0)
	{
		tesselationList.reserve(16);
	}


	FTMesh::~FTMesh()
	{
		for(size_t t = 0; t < tesselationList.size(); ++t)
		{
			delete tesselationList[t];
		}

		tesselationList.clear();
	}


	void FTMesh::AddPoint(const FTGL_DOUBLE x, const FTGL_DOUBLE y, const FTGL_DOUBLE z)
	{
		currentTesselation->AddPoint(x, y, z);
	}


	const FTGL_DOUBLE* FTMesh::Combine(const FTGL_DOUBLE x, const FTGL_DOUBLE y, const FTGL_DOUBLE z)
	{
		tempPointList.push_back(FTPoint(x, y,z));
		return static_cast<const FTGL_DOUBLE*>(tempPointList.back());
	}


	void FTMesh::Begin(GLenum meshType)
	{
		currentTesselation = new FTTesselation(meshType);
	}


	void FTMesh::End()
	{
		tesselationList.push_back(currentTesselation);
	}


	const FTTesselation* const FTMesh::Tesselation(size_t index) const
	{
		return (index < tesselationList.size()) ? tesselationList[index] : NULL;
	}




	void CALLBACK ftglError(GLenum errCode, FTMesh* mesh)
	{
		mesh->Error(errCode);
	}


	void CALLBACK ftglVertex(void* data, FTMesh* mesh)
	{
		FTGL_DOUBLE* vertex = static_cast<FTGL_DOUBLE*>(data);
		mesh->AddPoint(vertex[0], vertex[1], vertex[2]);
	}


	void CALLBACK ftglCombine(FTGL_DOUBLE coords[3], void* vertex_data[4], GLfloat weight[4], void** outData, FTMesh* mesh)
	{
		const FTGL_DOUBLE* vertex = static_cast<const FTGL_DOUBLE*>(coords);
		*outData = const_cast<FTGL_DOUBLE*>(mesh->Combine(vertex[0], vertex[1], vertex[2]));
	}

	void CALLBACK ftglBegin(GLenum type, FTMesh* mesh)
	{
		mesh->Begin(type);
	}


	void CALLBACK ftglEnd(FTMesh* mesh)
	{
		mesh->End();
	}
	////////TODO: need delete end///////////////



	void PolygonCarve::buildAndSubdivideUnderSuface()
	{
		FTMesh meshData;
		FTMesh *mesh = &meshData;

		GLUtesselator* tobj = gluNewTess();

		gluTessCallback(tobj, GLU_TESS_BEGIN_DATA,     (GLUTesselatorFunction)ftglBegin);
		gluTessCallback(tobj, GLU_TESS_VERTEX_DATA,    (GLUTesselatorFunction)ftglVertex);
		gluTessCallback(tobj, GLU_TESS_COMBINE_DATA,   (GLUTesselatorFunction)ftglCombine);
		gluTessCallback(tobj, GLU_TESS_END_DATA,       (GLUTesselatorFunction)ftglEnd);
		gluTessCallback(tobj, GLU_TESS_ERROR_DATA,     (GLUTesselatorFunction)ftglError);

		gluTessProperty(tobj, GLU_TESS_TOLERANCE, 0);
		gluTessNormal(tobj, 0.0f, 1.0f, 0.0f);
		gluTessBeginPolygon(tobj, mesh);

		int index = 0;
		for (size_t i = 0; i < fontPolygon.rings.size(); i++)
		{
			LKRing &ring = fontPolygon.rings[i];
			gluTessBeginContour(tobj);
			for (size_t j = 0; j < ring.size(); j++, index++)
			{
				ring[j].v.y = 0;
				gluTessVertex(tobj, (GLdouble *)&ring[j], (GLvoid *)&ring[j]);
			}

			gluTessEndContour(tobj);
		}
		gluTessEndPolygon(tobj);

		gluDeleteTess(tobj);


		// 初始化线段
		fontPolygonMesh.fb_fontPolydata = vtkSmartPointer<vtkPolyData>::New();
		fontPolygonMesh.fb_fontPolydata->SetPolys(vtkSmartPointer<vtkCellArray>::New());
		fontPolygonMesh.fb_fontPolydata->SetPoints(vtkSmartPointer<vtkPoints>::New());
		fontPolygonMesh.fb_locator.TakeReference(vtkMergePoints::New());
		double *polygonBounds = fontPolygon.getBounds();
		fontPolygonMesh.fb_locator->InitPointInsertion(fontPolygonMesh.fb_fontPolydata->GetPoints(), polygonBounds);



		for(size_t j = 0; j < mesh->TesselationCount(); ++j)
		{
			const FTTesselation* subMesh = mesh->Tesselation(j);
			unsigned int polygonType = subMesh->PolygonType();
			// 生成三角形
			switch(polygonType) {
			case GL_TRIANGLES:
				{
					for(unsigned int i = 0; i < subMesh->PointCount(); i = i + 3)
					{
						NewLKPoint pt0 = *(NewLKPoint *)&subMesh->Point(i + 0);
						NewLKPoint pt1 = *(NewLKPoint *)&subMesh->Point(i + 1);
						NewLKPoint pt2 = *(NewLKPoint *)&subMesh->Point(i + 2);

						addUnderTriangle(*(NewLKPoint *)&pt0, *(NewLKPoint *)&pt2, *(NewLKPoint *)&pt1);
					}
				}
				break;
			case GL_TRIANGLE_STRIP:
				for(unsigned int i = 0; i < subMesh->PointCount(); ++i)
				{
					NewLKPoint pt = *(NewLKPoint *)&subMesh->Point(i);
					if(i>=2)
					{
						NewLKPoint pt_prev_2 =  *(NewLKPoint *)&subMesh->Point(i-2);
						NewLKPoint pt_prev_1 =  *(NewLKPoint *)&subMesh->Point(i-1);
						if (i % 2 == 0)
						{
							addUnderTriangle(*(NewLKPoint *)&pt_prev_2, *(NewLKPoint *)&pt, *(NewLKPoint *)&pt_prev_1);
						}
						else {
							addUnderTriangle(*(NewLKPoint *)&pt_prev_1, *(NewLKPoint *)&pt, *(NewLKPoint *)&pt_prev_2);
						}
					}
				}
				break;
			case GL_TRIANGLE_FAN:
				for(unsigned int i = 0; i < subMesh->PointCount(); ++i)
				{
					NewLKPoint pt = *(NewLKPoint *)&subMesh->Point(i);
					if(i>=2)
					{
						NewLKPoint pt_prev_origin =  *(NewLKPoint *)&subMesh->Point(0);
						NewLKPoint pt_prev_1 =  *(NewLKPoint *)&subMesh->Point(i-1);
						addUnderTriangle(*(NewLKPoint *)&pt_prev_origin, *(NewLKPoint *)&pt, *(NewLKPoint *)&pt_prev_1);
					}
				}
				break;
			default:
				{
					fprintf(stderr, "没有考虑的polygon type : %i \n", polygonType);
					assert(0);
				}
			}
		}



		/// 太窄的三角形删除 两个夹角小于
		fontPolygonMesh.cleanNarrowTriangles(0.2);





		/// 浓缩顶点初始化
		fontPolygonMesh.subdivideMergePts = vtkSmartPointer<vtkPoints>::New();
		fontPolygonMesh.subdivideLocator.TakeReference(vtkMergePoints::New());
		double *polygonBounds1 = fontPolygon.getBounds();
		polygonBounds1[2] = 0;
		polygonBounds1[3] = 0;
		fontPolygonMesh.subdivideLocator->InitPointInsertion(fontPolygonMesh.subdivideMergePts, polygonBounds1);


		// 初始化点
		vtkIdType ptNum = fontPolygonMesh.fb_fontPolydata->GetNumberOfPoints();
		//points.reserve(ptNum * 1.01);
		for (int i = 0; i < ptNum; i++)
		{
			double *pt = fontPolygonMesh.fb_fontPolydata->GetPoint(i);
			fontPolygonMesh.points.push_back(*(NewLKPoint *)pt);
		}
		// 添加点到浓缩点阵
		for (int i = 0; i < fontPolygonMesh.points.size(); i++)
		{
			// 添加点
			NewLKPoint insertPt = fontPolygonMesh.points[i];
			insertPt.v.y = 0;
			fontPolygonMesh.subdivideLocator->InsertNextPoint((double *)&insertPt);
			//fontPolygonMesh.points.push_back(insertPt);
			this->fontPolygonMesh.nDeltaYSubdividePts.push_back(0);
			// 共点映射
			fontPolygonMesh.doublelicationMap.push_back(-1);
			// 父线段赋值
			Array<vtkIdType, 2> pSegment;
			pSegment.value[0] = pSegment.value[1] = -1;
			fontPolygonMesh.subdividePtParentSegments.push_back(pSegment);
		}



		// mesh查找并创建边界
		fontPolygonMesh.findBoundEdgess();


		// 添加细分三角形，三角形与边建立联系
		vtkIdType npts, *pts;
		for (vtkIdType cellId = 0; cellId < fontPolygonMesh.fb_fontPolydata->GetNumberOfCells(); cellId++)
		{
			fontPolygonMesh.fb_fontPolydata->GetCellPoints(cellId, npts, pts);

			NewLKPoint vertices[3];
			for (int i = 0; i < 3; i++)
			{
				vertices[i] = fontPolygonMesh.points[pts[i]];
			}

			PMTriangle *pmTri = new PMTriangle();
			pmTri->level = 0;
			pmTri->isLeaf = true;
			double longestEdge2 = -1;
			for (int i = 0; i < 3; i++)
			{
				pmTri->boundEdges[i] = (BoundEdge *)NULL;

				NewLKPoint *pt = vertices+i, *ptPrev = vertices + (i+2)%3;
				PtRingIds *rIds = (PtRingIds *)&(pt->v.y), *rIdsPrev = (PtRingIds *)&(ptPrev->v.y);
				//插入要merge的点
				pmTri->ids[i] = rIds->ptId;
				//vtkIdType ptIndex;
				//STPoint insertPt = *pt;
				//insertPt.y = rIds->ptId;
				//subdivideLocator->InsertUniquePoint((double *)&insertPt, ptIndex);
				//pmTri->ids[i] = ptIndex;

				double a = pt->v.x > ptPrev->v.x ? (pt->v.x - ptPrev->v.x) : (ptPrev->v.x - pt->v.x);
				double b = pt->v.z > ptPrev->v.z ? (pt->v.z - ptPrev->v.z) : (ptPrev->v.z - pt->v.z);
				pmTri->edgeLen2s[i] = a*a + b*b;
				if (longestEdge2 < pmTri->edgeLen2s[i])
				{
					longestEdge2 = pmTri->edgeLen2s[i];
					pmTri->longestEdgeId = i;
				}
				// 多边形边缘检测
				if (rIds->firstId == rIdsPrev->firstId)
				{
					int delta = rIds->secondId - rIdsPrev->secondId;
					BoundEdge *boundEdge;
					if (delta == 1 || delta == -fontPolygonMesh.boundEdgeNodess[rIds->firstId].size()+1)
					{
						boundEdge = fontPolygonMesh.boundEdgeNodess[rIds->firstId][rIds->secondId];
						boundEdge->sameClockwiseAsTriangle = true;
						boundEdge->pmTri = pmTri;
						boundEdge->edgeIndex = i;
						pmTri->boundEdges[i] = boundEdge;
					}
					else if (delta == -1 || delta == fontPolygonMesh.boundEdgeNodess[rIds->firstId].size()-1)
					{
						boundEdge = fontPolygonMesh.boundEdgeNodess[rIdsPrev->firstId][rIdsPrev->secondId];
						boundEdge->sameClockwiseAsTriangle = false;
						boundEdge->pmTri = pmTri;
						boundEdge->edgeIndex = i;
						pmTri->boundEdges[i] = boundEdge;
					}
				}
			}
			fontPolygonMesh.cellNodes.push_back(pmTri);
			// 添加到细化三角形队列
			if (longestEdge2 > _subdivideLen2)
			{
				fontPolygonMesh.subdividingTri.push(pmTri);
				//this->subdivideTriangle(pmTri);
			}
		}




		/// 细分三角形
		while (!fontPolygonMesh.subdividingTri.empty())
		{
			PMTriangle *pmTri = fontPolygonMesh.subdividingTri.top();
			fontPolygonMesh.subdividingTri.pop();
			subdivideTriangle(pmTri);
		}



		///还原顶点y值
		for (size_t i = 0; i < fontPolygonMesh.points.size(); i++)
		{
			fontPolygonMesh.points[i].v.y = this->fontPolygonY;
		}
	}

	// 二、没有点在三角形里
	void triangulateCross_2(STWorkTriangle &workTriangle, STWorkLine line, NewLKPoint &pt1, NewLKPoint &pt2,
		PointWithPolygon &pt1WithTriangle, PointWithPolygon &pt2WithTriangle,
		size_t &pt1OnIndex, size_t &pt2OnIndex,
		std::vector<STWorkTriangle> &tempNewWorkTris,
		std::vector<STWorkTriangle *> &tempEraseTris,
		std::vector<NewLKPoint> &intersections)
	{
		//   判断是否与三条边有交点
		//判断pt1-pt2 与A-B B-C C-A是否有交点
		bool isIntersection[3] = {false, false, false};
		bool isOnEdge[3] = {false, false, false};
		bool isPt1OnEdge[3] = {false, false, false};
		if (pt1WithTriangle == PointWithPolygonOnEdge)
		{
		 	isIntersection[pt1OnIndex] = true;
			isOnEdge[pt1OnIndex] = true;
			isPt1OnEdge[pt1OnIndex] = true;
		}
		if (pt2WithTriangle == PointWithPolygonOnEdge)
		{
		 	isIntersection[pt2OnIndex] = true;
			isOnEdge[pt2OnIndex] = true;
		}
		size_t intersectionCount = 0;
		for (int ii = 0, jj = 2; ii < 3; ii++)
		{
			isIntersection[jj] = isIntersection[jj] || isSegmentIntersectXZ(workTriangle.line[jj], line);
			if ( isIntersection[jj] )
				intersectionCount++;
			jj = ii;
		}

		//   如果有交点生成3个新三角形，剔除一个三角形
		if (intersectionCount == 2)
		{
			// 把变换成通用情况：pt1-pt2 不与 A-B相交
			// 交点pt1Intersect是与B-C的交点， pt2Intersect为C-A
			STWorkTriangle tempTriangle;
			size_t notIntersectEdgeIndex;
			for (size_t ii = 0; ii < 3; ii++)
			{
				if (!isIntersection[ii])
				{
					notIntersectEdgeIndex = ii;
					break;
				}
			}
			tempTriangle.pt[0] = workTriangle.pt[(notIntersectEdgeIndex+0)%3];
			tempTriangle.positionEnum[0] = workTriangle.positionEnum[(notIntersectEdgeIndex+0)%3];
			tempTriangle.pt[1] = workTriangle.pt[(notIntersectEdgeIndex+1)%3];
			tempTriangle.positionEnum[1] = workTriangle.positionEnum[(notIntersectEdgeIndex+1)%3];
			tempTriangle.pt[2] = workTriangle.pt[(notIntersectEdgeIndex+2)%3];
			tempTriangle.positionEnum[2] = workTriangle.positionEnum[(notIntersectEdgeIndex+2)%3];
			bool tempIsOnEdge[3];
			tempIsOnEdge[0] = isOnEdge[(notIntersectEdgeIndex+0)%3];
			tempIsOnEdge[1] = isOnEdge[(notIntersectEdgeIndex+1)%3];
			tempIsOnEdge[2] = isOnEdge[(notIntersectEdgeIndex+2)%3];
			bool tempIsPt1OnEdge[3];
			tempIsPt1OnEdge[0] = isPt1OnEdge[(notIntersectEdgeIndex+0)%3];
			tempIsPt1OnEdge[1] = isPt1OnEdge[(notIntersectEdgeIndex+1)%3];
			tempIsPt1OnEdge[2] = isPt1OnEdge[(notIntersectEdgeIndex+2)%3];

			// 求交点坐标
			NewLKPoint pt1Intersect = tempIsOnEdge[1] ? (isPt1OnEdge[1]?pt1:pt2) : segmentIntersetionXZ(tempTriangle.pt[1], tempTriangle.pt[2], pt1, pt2);
			NewLKPoint pt2Intersect = tempIsOnEdge[2] ? (isPt1OnEdge[2]?pt1:pt2) : segmentIntersetionXZ(tempTriangle.pt[2], tempTriangle.pt[0], pt1, pt2);

			// 添加交点到列表
			AddPointNoRepeatReturnIndex(intersections, pt1Intersect);
			AddPointNoRepeatReturnIndex(intersections, pt2Intersect);

			// 切成3个新三角形
			STWorkTriangle newTriangles[3];

			newTriangles[0].pt[0] = tempTriangle.pt[0];
			newTriangles[0].positionEnum[0] = tempTriangle.positionEnum[0];
			newTriangles[0].pt[1] = tempTriangle.pt[1];
			newTriangles[0].positionEnum[1] = tempTriangle.positionEnum[1];
			newTriangles[0].pt[2] = pt2Intersect;
			newTriangles[0].positionEnum[2] = PositionOnEdge;
			newTriangles[0].calcBounds();

			newTriangles[1].pt[0] = pt2Intersect;
			newTriangles[1].positionEnum[0] = PositionOnEdge;
			newTriangles[1].pt[1] = tempTriangle.pt[1];
			newTriangles[1].positionEnum[1] =  tempTriangle.positionEnum[1];
			newTriangles[1].pt[2] = pt1Intersect;
			newTriangles[1].positionEnum[2] = PositionOnEdge;
			newTriangles[1].calcBounds();

			newTriangles[2].pt[0] = pt2Intersect;
			newTriangles[2].positionEnum[0] = PositionOnEdge;
			newTriangles[2].pt[1] = pt1Intersect;
			newTriangles[2].positionEnum[1] = PositionOnEdge;
			newTriangles[2].pt[2] = tempTriangle.pt[2];
			newTriangles[2].positionEnum[2] = tempTriangle.positionEnum[2];
			newTriangles[2].calcBounds();

			tempEraseTris.push_back(&workTriangle);

			// 添加必须放到后面，否则地址改变
			tempNewWorkTris.push_back(newTriangles[0]);
			tempNewWorkTris.push_back(newTriangles[1]);
			tempNewWorkTris.push_back(newTriangles[2]);
		}
		else if (intersectionCount == 0)
		{
			//do nothing
		}
		else if (intersectionCount == 1)
		{
			fprintf(stderr, "两点在三角形外，交点数是%i，可能直线太靠近三角形某个边或某个顶点 \n", intersectionCount);
			assert(0);
		}
		else if (intersectionCount == 3)
		{
			fprintf(stderr, "两点在三角形外，交点数是%i，可能太靠近三角形某个边 \n", intersectionCount);
			assert(0);
		}

	}


	//   2. 一个点在三角形里
	void triangulateCross_1_2(STWorkTriangle &workTriangle, STWorkLine line, NewLKPoint &pt1, NewLKPoint &pt2,
		PointWithPolygon &pt1WithTriangle, PointWithPolygon &pt2WithTriangle,
		size_t &pt1OnIndex, size_t &pt2OnIndex,
		std::vector<STWorkTriangle> &tempNewWorkTris,
		std::vector<STWorkTriangle *> &tempEraseTris,
		std::vector<NewLKPoint> &intersections)
	{
		if (pt1WithTriangle == PointWithPolygonOnVertex && pt2WithTriangle == PointWithPolygonOnVertex)
		{
			//do nothing 两个点都在三角形顶点上（没有这种情况）
			fprintf(stderr, "两点在三角形内外，两点都在三角形顶点上\n");
			assert(0);
		}
		else if (pt1WithTriangle == PointWithPolygonOnVertex || pt2WithTriangle == PointWithPolygonOnVertex)
		{
			size_t ptOnTriVertexIndex = pt1WithTriangle == PointWithPolygonOnVertex ? pt1OnIndex : pt2OnIndex;
			// 把三角形转化为通用的三角形: pt1Pie在三角形顶点B上
			size_t offset_0 = (ptOnTriVertexIndex+2)%3;
			size_t offset_1 = (ptOnTriVertexIndex+0)%3;
			size_t offset_2 = (ptOnTriVertexIndex+1)%3;

			NewLKPoint pt1Pie = workTriangle.pt[offset_1];

			//判断pt1-pt2 与 C-A相交
			bool isIntersect1 = isSegmentIntersectXZ(workTriangle.line[offset_2], line);

			if (isIntersect1)
			{
				NewLKPoint ptIntersect = segmentIntersetionXZ(workTriangle.pt[offset_2], workTriangle.pt[offset_0], pt1, pt2);

				// 添加交点到列表
				AddPointNoRepeatReturnIndex(intersections, ptIntersect);

				AddPointNoRepeatReturnIndex(intersections, pt1Pie); //把端点加入交点表

				// 添加新生成的三角面片
				STWorkTriangle newTriangles[2];

				newTriangles[0].pt[0] = workTriangle.pt[offset_0];
				newTriangles[0].positionEnum[0] = workTriangle.positionEnum[offset_0];
				newTriangles[0].pt[1] = pt1Pie;
				newTriangles[0].positionEnum[1] = PositionOnEdge;
				newTriangles[0].pt[2] = ptIntersect;
				newTriangles[0].positionEnum[2] = PositionOnEdge;
				newTriangles[0].calcBounds();

				newTriangles[1].pt[0] = ptIntersect;
				newTriangles[1].positionEnum[0] = PositionOnEdge;
				newTriangles[1].pt[1] = pt1Pie;
				newTriangles[1].positionEnum[1] = PositionOnEdge;
				newTriangles[1].pt[2] = workTriangle.pt[offset_2];
				newTriangles[1].positionEnum[2] = workTriangle.positionEnum[offset_2];
				newTriangles[1].calcBounds();

				tempEraseTris.push_back(&workTriangle);

				// 添加必须放到后面，否则地址改变
				tempNewWorkTris.push_back(newTriangles[0]);
				tempNewWorkTris.push_back(newTriangles[1]);

			}
			else
			{
				//do nothing
			}
			

		}
		else 
		{
			//判断pt1-pt2 与A-B B-C C-A是否有交点
			bool isIntersection[3] = {false, false, false};
			bool isOnEdge[3] = {false, false, false};
			bool isPt1OnEdge[3] = {false, false, false};
			if (pt1WithTriangle == PointWithPolygonOnEdge)
			{
				isIntersection[pt1OnIndex] = true;
				isOnEdge[pt1OnIndex] = true;
				isPt1OnEdge[pt1OnIndex] = true;
			}
			if (pt2WithTriangle == PointWithPolygonOnEdge)
			{
				isIntersection[pt2OnIndex] = true;
				isOnEdge[pt2OnIndex] = true;
			}
			size_t intersectionCount = 0;
			for (int ii = 0, jj = 2; ii < 3; ii++)
			{
				isIntersection[jj] = isIntersection[jj] || isSegmentIntersectXZ(workTriangle.line[jj], line);
				if ( isIntersection[jj] )
					intersectionCount++;
				jj = ii;
			}

			if (intersectionCount == 0)
			{
				fprintf(stderr, "两点分别在三角形内、外，交点数为0， 可能原因\n    1.isSegmentIntersectXZ2程序有问题\n    2.判断点在三角形内外有问题  \n");
				assert(0);
			}

			// 如果交点个数为 1个， 切4个三角形
			if (intersectionCount == 1)
			{
				// 把三角形转化为通用的三角形
				//   1.pt1-pt2 与 C-A相交
				//   2. pt1Pie在三角形内
				//   3. ptIntersect为交点
				STWorkTriangle tempTriangle;
				size_t intersectEdgeIndex;
				for (size_t ii = 0; ii < 3; ii++)
				{
					if (isIntersection[ii])
					{
						intersectEdgeIndex = ii;
						break;
					}
				}
				tempTriangle.pt[0] = workTriangle.pt[(intersectEdgeIndex+1)%3];
				tempTriangle.positionEnum[0] = workTriangle.positionEnum[(intersectEdgeIndex+1)%3];
				tempTriangle.pt[1] = workTriangle.pt[(intersectEdgeIndex+2)%3];
				tempTriangle.positionEnum[1] = workTriangle.positionEnum[(intersectEdgeIndex+2)%3];
				tempTriangle.pt[2] = workTriangle.pt[(intersectEdgeIndex+0)%3];
				tempTriangle.positionEnum[2] = workTriangle.positionEnum[(intersectEdgeIndex+0)%3];

				NewLKPoint pt1Pie = pt1WithTriangle != PointWithPolygonOutside ? pt1 : pt2;
				AddPointNoRepeatReturnIndex(intersections, pt1Pie); //把端点加入交点表


				// pt1Pie在C-A边上
				if (pt1WithTriangle == PointWithPolygonOnEdge || pt2WithTriangle == PointWithPolygonOnEdge)
				{// 添加新生成的三角面片
					STWorkTriangle newTriangles[2];

					newTriangles[0].pt[0] = tempTriangle.pt[0];
					newTriangles[0].positionEnum[0] = tempTriangle.positionEnum[0];
					newTriangles[0].pt[1] = tempTriangle.pt[1];
					newTriangles[0].positionEnum[1] = tempTriangle.positionEnum[1];
					newTriangles[0].pt[2] = pt1Pie;
					newTriangles[0].positionEnum[2] = PositionOnEdge;
					newTriangles[0].calcBounds();

					newTriangles[1].pt[0] = pt1Pie;
					newTriangles[1].positionEnum[0] = PositionOnEdge;
					newTriangles[1].pt[1] = tempTriangle.pt[1];
					newTriangles[1].positionEnum[1] = tempTriangle.positionEnum[1];
					newTriangles[1].pt[2] = tempTriangle.pt[2];
					newTriangles[1].positionEnum[2] = tempTriangle.positionEnum[2];
					newTriangles[1].calcBounds();

					tempEraseTris.push_back(&workTriangle);

					// 添加必须放到后面，否则地址改变
					tempNewWorkTris.push_back(newTriangles[0]);
					tempNewWorkTris.push_back(newTriangles[1]);


				}
				// pt1Pie不在C-A边上
				else {
					// 求交点坐标
					NewLKPoint ptIntersect = segmentIntersetionXZ(tempTriangle.pt[2], tempTriangle.pt[0], pt1, pt2);

					// 添加交点到列表
					AddPointNoRepeatReturnIndex(intersections, ptIntersect);

					// 添加新生成的三角面片
					STWorkTriangle newTriangles[4];

					newTriangles[0].pt[0] = tempTriangle.pt[0];
					newTriangles[0].positionEnum[0] = tempTriangle.positionEnum[0];
					newTriangles[0].pt[1] = tempTriangle.pt[1];
					newTriangles[0].positionEnum[1] = tempTriangle.positionEnum[1];
					newTriangles[0].pt[2] = pt1Pie;
					newTriangles[0].positionEnum[2] = PositionOnEdge;
					newTriangles[0].calcBounds();

					newTriangles[1].pt[0] = tempTriangle.pt[0];
					newTriangles[1].positionEnum[0] = tempTriangle.positionEnum[0];
					newTriangles[1].pt[1] = pt1Pie;
					newTriangles[1].positionEnum[1] = PositionOnEdge;
					newTriangles[1].pt[2] = ptIntersect;
					newTriangles[1].positionEnum[2] = PositionOnEdge;
					newTriangles[1].calcBounds();

					newTriangles[2].pt[0] = pt1Pie;
					newTriangles[2].positionEnum[0] = PositionOnEdge;
					newTriangles[2].pt[1] = tempTriangle.pt[1];
					newTriangles[2].positionEnum[1] = tempTriangle.positionEnum[1];
					newTriangles[2].pt[2] = tempTriangle.pt[2];
					newTriangles[2].positionEnum[2] = tempTriangle.positionEnum[2];
					newTriangles[2].calcBounds();

					newTriangles[3].pt[0] = ptIntersect;
					newTriangles[3].positionEnum[0] = PositionOnEdge;
					newTriangles[3].pt[1] = pt1Pie;
					newTriangles[3].positionEnum[1] = PositionOnEdge;
					newTriangles[3].pt[2] = tempTriangle.pt[2];
					newTriangles[3].positionEnum[2] = tempTriangle.positionEnum[2];
					newTriangles[3].calcBounds();

					tempEraseTris.push_back(&workTriangle);

					// 添加必须放到后面，否则地址改变
					tempNewWorkTris.push_back(newTriangles[0]);
					tempNewWorkTris.push_back(newTriangles[1]);
					tempNewWorkTris.push_back(newTriangles[2]);
					tempNewWorkTris.push_back(newTriangles[3]);
				}


			}

			// 如果交点个数为 2个， 肯定三角形内的点在边上
			if (intersectionCount == 2)
			{
				if (pt1WithTriangle != PointWithPolygonOnEdge && pt2WithTriangle != PointWithPolygonOnEdge)
				{
					fprintf(stderr, "两点在三角形内外，没有点在边上，交点为两个，可能直线通过某个顶点\n");
					assert(0);
				}
				else
				{

					// 把变换成通用情况：pt1-pt2 不与 A-B相交
					// 交点pt1Intersect是与B-C的交点， pt2Intersect为C-A
					STWorkTriangle tempTriangle;
					size_t notIntersectEdgeIndex;
					for (size_t ii = 0; ii < 3; ii++)
					{
						if (!isIntersection[ii])
						{
							notIntersectEdgeIndex = ii;
							break;
						}
					}
					tempTriangle.pt[0] = workTriangle.pt[(notIntersectEdgeIndex+0)%3];
					tempTriangle.positionEnum[0] = workTriangle.positionEnum[(notIntersectEdgeIndex+0)%3];
					tempTriangle.pt[1] = workTriangle.pt[(notIntersectEdgeIndex+1)%3];
					tempTriangle.positionEnum[1] = workTriangle.positionEnum[(notIntersectEdgeIndex+1)%3];
					tempTriangle.pt[2] = workTriangle.pt[(notIntersectEdgeIndex+2)%3];
					tempTriangle.positionEnum[2] = workTriangle.positionEnum[(notIntersectEdgeIndex+2)%3];
					bool tempIsOnEdge[3];
					tempIsOnEdge[0] = isOnEdge[(notIntersectEdgeIndex+0)%3];
					tempIsOnEdge[1] = isOnEdge[(notIntersectEdgeIndex+1)%3];
					tempIsOnEdge[2] = isOnEdge[(notIntersectEdgeIndex+2)%3];
					bool tempIsPt1OnEdge[3];
					tempIsPt1OnEdge[0] = isPt1OnEdge[(notIntersectEdgeIndex+0)%3];
					tempIsPt1OnEdge[1] = isPt1OnEdge[(notIntersectEdgeIndex+1)%3];
					tempIsPt1OnEdge[2] = isPt1OnEdge[(notIntersectEdgeIndex+2)%3];

					// 求交点坐标
					NewLKPoint pt1Intersect = tempIsOnEdge[1] ? (tempIsPt1OnEdge[1]?pt1:pt2) : segmentIntersetionXZ(tempTriangle.pt[1], tempTriangle.pt[2], pt1, pt2);
					NewLKPoint pt2Intersect = tempIsOnEdge[2] ? (tempIsPt1OnEdge[2]?pt1:pt2) : segmentIntersetionXZ(tempTriangle.pt[2], tempTriangle.pt[0], pt1, pt2);

					// 添加交点到列表
					AddPointNoRepeatReturnIndex(intersections, pt1Intersect);
					AddPointNoRepeatReturnIndex(intersections, pt2Intersect);

					// 切成3个新三角形
					STWorkTriangle newTriangles[3];

					newTriangles[0].pt[0] = tempTriangle.pt[0];
					newTriangles[0].positionEnum[0] = tempTriangle.positionEnum[0];
					newTriangles[0].pt[1] = tempTriangle.pt[1];
					newTriangles[0].positionEnum[1] = tempTriangle.positionEnum[1];
					newTriangles[0].pt[2] = pt2Intersect;
					newTriangles[0].positionEnum[2] = PositionOnEdge;
					newTriangles[0].calcBounds();

					newTriangles[1].pt[0] = pt2Intersect;
					newTriangles[1].positionEnum[0] = PositionOnEdge;
					newTriangles[1].pt[1] = tempTriangle.pt[1];
					newTriangles[1].positionEnum[1] =  tempTriangle.positionEnum[1];
					newTriangles[1].pt[2] = pt1Intersect;
					newTriangles[1].positionEnum[2] = PositionOnEdge;
					newTriangles[1].calcBounds();

					newTriangles[2].pt[0] = pt2Intersect;
					newTriangles[2].positionEnum[0] = PositionOnEdge;
					newTriangles[2].pt[1] = pt1Intersect;
					newTriangles[2].positionEnum[1] = PositionOnEdge;
					newTriangles[2].pt[2] = tempTriangle.pt[2];
					newTriangles[2].positionEnum[2] = tempTriangle.positionEnum[2];
					newTriangles[2].calcBounds();

					tempEraseTris.push_back(&workTriangle);

					// 添加必须放到后面，否则地址改变
					tempNewWorkTris.push_back(newTriangles[0]);
					tempNewWorkTris.push_back(newTriangles[1]);
					tempNewWorkTris.push_back(newTriangles[2]);
				}
			}

			// 如果交点个数为 3个 错误
			if (intersectionCount == 3)
			{
				fprintf(stderr, "两点分别在三角形内、外，交点数为3个!!\n");
				assert(0);
			}
		}
	}

	//   1. 两个点都在三角形里
	void triangulateFivePoint(STWorkTriangle &workTriangle, STWorkLine line, NewLKPoint &pt1, NewLKPoint &pt2,
		PointWithPolygon &pt1WithTriangle, PointWithPolygon &pt2WithTriangle,
		size_t &pt1OnIndex, size_t &pt2OnIndex,
		std::vector<STWorkTriangle> &tempNewWorkTris,
		std::vector<STWorkTriangle *> &tempEraseTris,
		std::vector<NewLKPoint> &intersections)
	{
		bool point1OnTriangleVertices = false, point2OnTriangleVertices = false;
		size_t onTriangleVertexIndex;
		NewLKPoint pt1Pie, pt2Pie;//重合的点在pt1Pie和tempTriangle.pt[0]上
		if (pt1WithTriangle == PointWithPolygonOnVertex)
		{
			point1OnTriangleVertices = true;
			onTriangleVertexIndex = pt1OnIndex;
			pt1Pie = pt1;
			pt2Pie = pt2;
		}
		if (pt2WithTriangle == PointWithPolygonOnVertex)
		{
			point2OnTriangleVertices = true;
			onTriangleVertexIndex = pt2OnIndex;
			pt1Pie = pt2;
			pt2Pie = pt1;
		}

		//判断pt1-pt2 与A-B B-C C-A是否有交点
		bool isIntersection[3] = {false, false, false};
		bool isOnEdge[3] = {false, false, false};
		bool isPt1OnEdge[3] = {false, false, false};
		if (pt1WithTriangle == PointWithPolygonOnEdge)
		{
			isIntersection[pt1OnIndex] = true;
			isOnEdge[pt1OnIndex] = true;
			isPt1OnEdge[pt1OnIndex] = true;
		}
		if (pt2WithTriangle == PointWithPolygonOnEdge)
		{
			isIntersection[pt2OnIndex] = true;
			isOnEdge[pt2OnIndex] = true;
		}
		size_t intersectionCount = 0;
		for (int ii = 0, jj = 2; ii < 3; ii++)
		{
			if ( isIntersection[jj] )
				intersectionCount++;
			jj = ii;
		}

		if (point1OnTriangleVertices && point2OnTriangleVertices)
		{
			//printf("两点都在三角形顶点上\n");
			if (pt1OnIndex == pt2OnIndex)
			{
				printf("两点在三角形同一顶点（%i）上 \n", pt1OnIndex);
				//AddPointNoRepeatReturnIndex(intersections, pt1); //把端点加入交点表
				//AddPointNoRepeatReturnIndex(intersections, pt2); //把端点加入交点表
			}
			else
			{
				//printf("两点在三角形两顶点(%i)(%i)上 \n", pt1OnIndex, pt2OnIndex);
				AddPointNoRepeatReturnIndex(intersections, pt1); //把端点加入交点表
				AddPointNoRepeatReturnIndex(intersections, pt2); //把端点加入交点表
			}
		}
		else if (point1OnTriangleVertices || point2OnTriangleVertices)
		{
			// 通用情况 pt1Pie在顶点A上
			STWorkTriangle tempTriangle;
			tempTriangle.pt[0] = workTriangle.pt[onTriangleVertexIndex];
			tempTriangle.positionEnum[0] = workTriangle.positionEnum[onTriangleVertexIndex];
			tempTriangle.pt[1] = workTriangle.pt[(onTriangleVertexIndex+1)%3];
			tempTriangle.positionEnum[1] = workTriangle.positionEnum[(onTriangleVertexIndex+1)%3];
			tempTriangle.pt[2] = workTriangle.pt[(onTriangleVertexIndex+2)%3];
			tempTriangle.positionEnum[2] = workTriangle.positionEnum[(onTriangleVertexIndex+2)%3];
			bool tempIsOnEdge[3];
			tempIsOnEdge[0] = isOnEdge[(onTriangleVertexIndex+0)%3];
			tempIsOnEdge[1] = isOnEdge[(onTriangleVertexIndex+1)%3];
			tempIsOnEdge[2] = isOnEdge[(onTriangleVertexIndex+2)%3];

			// 合并相近点时有用
			pt1Pie = tempTriangle.pt[0];

			if (intersectionCount == 0)
			{
				AddPointNoRepeatReturnIndex(intersections, pt1Pie); //把端点加入交点表
				AddPointNoRepeatReturnIndex(intersections, pt2Pie); //把端点加入交点表

				STWorkTriangle newTriangles[3];

				newTriangles[0].pt[0] = tempTriangle.pt[0];
				newTriangles[0].positionEnum[0] = tempTriangle.positionEnum[0];
				newTriangles[0].pt[1] = tempTriangle.pt[1];
				newTriangles[0].positionEnum[1] = tempTriangle.positionEnum[1];
				newTriangles[0].pt[2] = pt2Pie;
				newTriangles[0].positionEnum[2] = PositionOnEdge;
				newTriangles[0].calcBounds();

				newTriangles[1].pt[0] = tempTriangle.pt[0];
				newTriangles[1].positionEnum[0] = tempTriangle.positionEnum[0];
				newTriangles[1].pt[1] = pt2Pie;
				newTriangles[1].positionEnum[1] = PositionOnEdge;
				newTriangles[1].pt[2] = tempTriangle.pt[2];
				newTriangles[1].positionEnum[2] = tempTriangle.positionEnum[2];
				newTriangles[1].calcBounds();

				newTriangles[2].pt[0] = pt2Pie;
				newTriangles[2].positionEnum[0] = PositionOnEdge;
				newTriangles[2].pt[1] = tempTriangle.pt[1];
				newTriangles[2].positionEnum[1] = tempTriangle.positionEnum[1];
				newTriangles[2].pt[2] = tempTriangle.pt[2];
				newTriangles[2].positionEnum[2] = tempTriangle.positionEnum[2];
				newTriangles[2].calcBounds();

				tempEraseTris.push_back(&workTriangle);

				// 添加必须放到后面，否则地址改变
				tempNewWorkTris.push_back(newTriangles[0]);
				tempNewWorkTris.push_back(newTriangles[1]);
				tempNewWorkTris.push_back(newTriangles[2]);
			}
			else if (intersectionCount == 1)
			{
				//printf("情况一.1 \n");
				if (tempIsOnEdge[1])
				{
					printf("    情况一.1.b \n");

					AddPointNoRepeatReturnIndex(intersections, pt1Pie); //把端点加入交点表
					AddPointNoRepeatReturnIndex(intersections, pt2Pie); //把端点加入交点表
 
  					STWorkTriangle newTriangles[2];

					newTriangles[0].pt[0] = pt1Pie;
					newTriangles[0].positionEnum[0] = PositionOnEdge;
					newTriangles[0].pt[1] = tempTriangle.pt[1];
					newTriangles[0].positionEnum[1] = tempTriangle.positionEnum[1];
					newTriangles[0].pt[2] = pt2Pie;
					newTriangles[0].positionEnum[2] = PositionOnEdge;
					newTriangles[0].calcBounds();

					newTriangles[1].pt[0] = pt1Pie;
					newTriangles[1].positionEnum[0] = PositionOnEdge;
					newTriangles[1].pt[1] = pt2Pie;
					newTriangles[1].positionEnum[1] = PositionOnEdge;
					newTriangles[1].pt[2] = tempTriangle.pt[2];
					newTriangles[1].positionEnum[2] = tempTriangle.positionEnum[2];
					newTriangles[1].calcBounds();
  
   					tempEraseTris.push_back(&workTriangle);

  					// 添加必须放到后面，否则地址改变
    				tempNewWorkTris.push_back(newTriangles[0]);
    				tempNewWorkTris.push_back(newTriangles[1]);

				}
				if (tempIsOnEdge[0])
				{
					printf("    情况一.1.c \n");

					AddPointNoRepeatReturnIndex(intersections, pt1Pie); //把端点加入交点表
					AddPointNoRepeatReturnIndex(intersections, pt2Pie); //把端点加入交点表

					STWorkTriangle newTriangles[2];

					newTriangles[0].pt[0] = pt1Pie;
					newTriangles[0].positionEnum[0] = PositionOnEdge;
					newTriangles[0].pt[1] = pt2Pie;
					newTriangles[0].positionEnum[1] = PositionOnEdge;
					newTriangles[0].pt[2] = tempTriangle.pt[2];
					newTriangles[0].positionEnum[2] = tempTriangle.positionEnum[2];
					newTriangles[0].calcBounds();

					newTriangles[1].pt[0] = pt2Pie;
					newTriangles[1].positionEnum[0] = PositionOnEdge;
					newTriangles[1].pt[1] = tempTriangle.pt[1];
					newTriangles[1].positionEnum[1] = tempTriangle.positionEnum[1];
					newTriangles[1].pt[2] = tempTriangle.pt[2];
					newTriangles[1].positionEnum[2] = tempTriangle.positionEnum[2];
					newTriangles[1].calcBounds();

					tempEraseTris.push_back(&workTriangle);

					// 添加必须放到后面，否则地址改变
					tempNewWorkTris.push_back(newTriangles[0]);
					tempNewWorkTris.push_back(newTriangles[1]);
				}

				if (tempIsOnEdge[2])
				{
					printf("    情况一.1.d \n");

					AddPointNoRepeatReturnIndex(intersections, pt1Pie); //把端点加入交点表
					AddPointNoRepeatReturnIndex(intersections, pt2Pie); //把端点加入交点表

					STWorkTriangle newTriangles[2];

					newTriangles[0].pt[0] = pt1Pie;
					newTriangles[0].positionEnum[0] = PositionOnEdge;
					newTriangles[0].pt[1] = tempTriangle.pt[1];
					newTriangles[0].positionEnum[1] = tempTriangle.positionEnum[1];
					newTriangles[0].pt[2] = pt2Pie;
					newTriangles[0].positionEnum[2] = PositionOnEdge;
					newTriangles[0].calcBounds();

					newTriangles[1].pt[0] = pt2Pie;
					newTriangles[1].positionEnum[0] = PositionOnEdge;
					newTriangles[1].pt[1] = tempTriangle.pt[1];
					newTriangles[1].positionEnum[1] = tempTriangle.positionEnum[1];
					newTriangles[1].pt[2] = tempTriangle.pt[2];
					newTriangles[1].positionEnum[2] = tempTriangle.positionEnum[2];
					newTriangles[1].calcBounds();

					tempEraseTris.push_back(&workTriangle);

					// 添加必须放到后面，否则地址改变
					tempNewWorkTris.push_back(newTriangles[0]);
					tempNewWorkTris.push_back(newTriangles[1]);
				}
			}
		}
		else
		{
			//printf("情况一.1/2/3 \n");

			// 两个点都在边上
			if (pt1WithTriangle == PointWithPolygonOnEdge && pt2WithTriangle == PointWithPolygonOnEdge)
			{
				//printf("情况一.2 \n");
				if (pt1OnIndex != pt2OnIndex)
				{
					printf("    情况一.2.a \n");

					// 把变换成通用情况：pt1-pt2 不与 A-B相交
					// 交点pt1Intersect是与B-C的交点， pt2Intersect为C-A
					STWorkTriangle tempTriangle;
					size_t notIntersectEdgeIndex;
					for (size_t ii = 0; ii < 3; ii++)
					{
						if (!isIntersection[ii])
						{
							notIntersectEdgeIndex = ii;
							break;
						}
					}
					tempTriangle.pt[0] = workTriangle.pt[(notIntersectEdgeIndex+0)%3];
					tempTriangle.positionEnum[0] = workTriangle.positionEnum[(notIntersectEdgeIndex+0)%3];
					tempTriangle.pt[1] = workTriangle.pt[(notIntersectEdgeIndex+1)%3];
					tempTriangle.positionEnum[1] = workTriangle.positionEnum[(notIntersectEdgeIndex+1)%3];
					tempTriangle.pt[2] = workTriangle.pt[(notIntersectEdgeIndex+2)%3];
					tempTriangle.positionEnum[2] = workTriangle.positionEnum[(notIntersectEdgeIndex+2)%3];
					bool tempIsOnEdge[3];
					tempIsOnEdge[0] = isOnEdge[(notIntersectEdgeIndex+0)%3];
					tempIsOnEdge[1] = isOnEdge[(notIntersectEdgeIndex+1)%3];
					tempIsOnEdge[2] = isOnEdge[(notIntersectEdgeIndex+2)%3];
					bool tempIsPt1OnEdge[3];
					tempIsPt1OnEdge[0] = isPt1OnEdge[(notIntersectEdgeIndex+0)%3];
					tempIsPt1OnEdge[1] = isPt1OnEdge[(notIntersectEdgeIndex+1)%3];
					tempIsPt1OnEdge[2] = isPt1OnEdge[(notIntersectEdgeIndex+2)%3];

					// 求交点坐标
					NewLKPoint pt1Intersect = tempIsOnEdge[1] ? (isPt1OnEdge[1]?pt1:pt2) : segmentIntersetionXZ(tempTriangle.pt[1], tempTriangle.pt[2], pt1, pt2);
					NewLKPoint pt2Intersect = tempIsOnEdge[2] ? (isPt1OnEdge[2]?pt1:pt2) : segmentIntersetionXZ(tempTriangle.pt[2], tempTriangle.pt[0], pt1, pt2);

					// 添加交点到列表
					AddPointNoRepeatReturnIndex(intersections, pt1Intersect);
					AddPointNoRepeatReturnIndex(intersections, pt2Intersect);

					// 切成3个新三角形
					STWorkTriangle newTriangles[3];

					newTriangles[0].pt[0] = tempTriangle.pt[0];
					newTriangles[0].positionEnum[0] = tempTriangle.positionEnum[0];
					newTriangles[0].pt[1] = tempTriangle.pt[1];
					newTriangles[0].positionEnum[1] = tempTriangle.positionEnum[1];
					newTriangles[0].pt[2] = pt2Intersect;
					newTriangles[0].positionEnum[2] = PositionOnEdge;
					newTriangles[0].calcBounds();

					newTriangles[1].pt[0] = pt2Intersect;
					newTriangles[1].positionEnum[0] = PositionOnEdge;
					newTriangles[1].pt[1] = tempTriangle.pt[1];
					newTriangles[1].positionEnum[1] =  tempTriangle.positionEnum[1];
					newTriangles[1].pt[2] = pt1Intersect;
					newTriangles[1].positionEnum[2] = PositionOnEdge;
					newTriangles[1].calcBounds();

					newTriangles[2].pt[0] = pt2Intersect;
					newTriangles[2].positionEnum[0] = PositionOnEdge;
					newTriangles[2].pt[1] = pt1Intersect;
					newTriangles[2].positionEnum[1] = PositionOnEdge;
					newTriangles[2].pt[2] = tempTriangle.pt[2];
					newTriangles[2].positionEnum[2] = tempTriangle.positionEnum[2];
					newTriangles[2].calcBounds();

					tempEraseTris.push_back(&workTriangle);

					// 添加必须放到后面，否则地址改变
					tempNewWorkTris.push_back(newTriangles[0]);
					tempNewWorkTris.push_back(newTriangles[1]);
					tempNewWorkTris.push_back(newTriangles[2]);

				}
				else if (pt1OnIndex == pt2OnIndex)
				{
					printf("    情况一.2.b \n");

					// 把变换成通用情况：pt1-pt2 在 C-A上，pt1靠近A
					STWorkTriangle tempTriangle;
					size_t intersectEdgeIndex;
					for (size_t ii = 0; ii < 3; ii++)
					{
						if (isIntersection[ii])
						{
							intersectEdgeIndex = ii;
							break;
						}
					}
					tempTriangle.pt[0] = workTriangle.pt[(intersectEdgeIndex+1)%3];
					tempTriangle.positionEnum[0] = workTriangle.positionEnum[(intersectEdgeIndex+1)%3];
					tempTriangle.pt[1] = workTriangle.pt[(intersectEdgeIndex+2)%3];
					tempTriangle.positionEnum[1] = workTriangle.positionEnum[(intersectEdgeIndex+2)%3];
					tempTriangle.pt[2] = workTriangle.pt[(intersectEdgeIndex+0)%3];
					tempTriangle.positionEnum[2] = workTriangle.positionEnum[(intersectEdgeIndex+0)%3];

					// pt1 pt2位置调整
					NewLKPoint ACVector, pt1pt2Vector;
					ACVector.v.x = tempTriangle.pt[2].v.x - tempTriangle.pt[0].v.x;
					ACVector.v.z = tempTriangle.pt[2].v.z - tempTriangle.pt[0].v.z;
					pt1pt2Vector.v.x = pt2.v.x - pt1.v.x;
					pt1pt2Vector.v.z = pt2.v.z - pt1.v.z;
					// ACVector的向量夹角约等于pt1pt2Vector的，不换位置
					// 向量同向检测
					bool notNeedSwap = ACVector.v.x * pt1pt2Vector.v.x > 1.e-16 || (abs(ACVector.v.x * pt1pt2Vector.v.x) <= 1.e-16 && ACVector.v.z * pt1pt2Vector.v.z > 0);
					pt1Pie = notNeedSwap ? pt1 : pt2;
					pt2Pie = notNeedSwap ? pt2 : pt1;

 					AddPointNoRepeatReturnIndex(intersections, pt1Pie); //把端点加入交点表
 					AddPointNoRepeatReturnIndex(intersections, pt2Pie); //把端点加入交点表

					STWorkTriangle newTriangles[3];

					newTriangles[0].pt[0] = tempTriangle.pt[0];
					newTriangles[0].positionEnum[0] = tempTriangle.positionEnum[0];
					newTriangles[0].pt[1] = tempTriangle.pt[1];
					newTriangles[0].positionEnum[1] = tempTriangle.positionEnum[1];
					newTriangles[0].pt[2] = pt1Pie;
					newTriangles[0].positionEnum[2] = PositionOnEdge;
					newTriangles[0].calcBounds();

					newTriangles[1].pt[0] = pt1Pie;
					newTriangles[1].positionEnum[0] = PositionOnEdge;
					newTriangles[1].pt[1] = tempTriangle.pt[1];
					newTriangles[1].positionEnum[1] = tempTriangle.positionEnum[1];
					newTriangles[1].pt[2] = pt2Pie;
					newTriangles[1].positionEnum[2] = PositionOnEdge;
					newTriangles[1].calcBounds();

					newTriangles[2].pt[0] = pt2Pie;
					newTriangles[2].positionEnum[0] = PositionOnEdge;
					newTriangles[2].pt[1] = tempTriangle.pt[1];
					newTriangles[2].positionEnum[1] = tempTriangle.positionEnum[1];
					newTriangles[2].pt[2] = tempTriangle.pt[2];
					newTriangles[2].positionEnum[2] = tempTriangle.positionEnum[2];
					newTriangles[2].calcBounds();

					tempEraseTris.push_back(&workTriangle);

					// 添加必须放到后面，否则地址改变
					tempNewWorkTris.push_back(newTriangles[0]);
					tempNewWorkTris.push_back(newTriangles[1]);
					tempNewWorkTris.push_back(newTriangles[2]);
				}
			}

			// 一个点在边上
			if (pt1WithTriangle == PointWithPolygonOnEdge ^ pt2WithTriangle == PointWithPolygonOnEdge)
			{
				//printf("情况一.1 \n");
				printf("    情况一.1.a \n");

				// 把三角形转化为通用的三角形
				//   1.pt1-pt2 与 C-A相交
				//   2. pt1Pie在三角形内
				//   3. ptIntersect为交点
				STWorkTriangle tempTriangle;
				size_t intersectEdgeIndex;
				for (size_t ii = 0; ii < 3; ii++)
				{
					if (isIntersection[ii])
					{
						intersectEdgeIndex = ii;
						break;
					}
				}
				tempTriangle.pt[0] = workTriangle.pt[(intersectEdgeIndex+1)%3];
				tempTriangle.positionEnum[0] = workTriangle.positionEnum[(intersectEdgeIndex+1)%3];
				tempTriangle.pt[1] = workTriangle.pt[(intersectEdgeIndex+2)%3];
				tempTriangle.positionEnum[1] = workTriangle.positionEnum[(intersectEdgeIndex+2)%3];
				tempTriangle.pt[2] = workTriangle.pt[(intersectEdgeIndex+0)%3];
				tempTriangle.positionEnum[2] = workTriangle.positionEnum[(intersectEdgeIndex+0)%3];

				NewLKPoint pt1Pie = pt1WithTriangle != PointWithPolygonOnEdge ? pt1 : pt2;
				AddPointNoRepeatReturnIndex(intersections, pt1Pie); //把端点加入交点表

				// 求交点坐标
				NewLKPoint ptIntersect = pt1WithTriangle == PointWithPolygonOnEdge ? pt1 : pt2;

				// 添加交点到列表
				AddPointNoRepeatReturnIndex(intersections, ptIntersect);

				// 添加新生成的三角面片
				STWorkTriangle newTriangles[4];

				newTriangles[0].pt[0] = tempTriangle.pt[0];
				newTriangles[0].positionEnum[0] = tempTriangle.positionEnum[0];
				newTriangles[0].pt[1] = tempTriangle.pt[1];
				newTriangles[0].positionEnum[1] = tempTriangle.positionEnum[1];
				newTriangles[0].pt[2] = pt1Pie;
				newTriangles[0].positionEnum[2] = PositionOnEdge;
				newTriangles[0].calcBounds();

				newTriangles[1].pt[0] = tempTriangle.pt[0];
				newTriangles[1].positionEnum[0] = tempTriangle.positionEnum[0];
				newTriangles[1].pt[1] = pt1Pie;
				newTriangles[1].positionEnum[1] = PositionOnEdge;
				newTriangles[1].pt[2] = ptIntersect;
				newTriangles[1].positionEnum[2] = PositionOnEdge;
				newTriangles[1].calcBounds();

				newTriangles[2].pt[0] = pt1Pie;
				newTriangles[2].positionEnum[0] = PositionOnEdge;
				newTriangles[2].pt[1] = tempTriangle.pt[1];
				newTriangles[2].positionEnum[1] = tempTriangle.positionEnum[1];
				newTriangles[2].pt[2] = tempTriangle.pt[2];
				newTriangles[2].positionEnum[2] = tempTriangle.positionEnum[2];
				newTriangles[2].calcBounds();

				newTriangles[3].pt[0] = ptIntersect;
				newTriangles[3].positionEnum[0] = PositionOnEdge;
				newTriangles[3].pt[1] = pt1Pie;
				newTriangles[3].positionEnum[1] = PositionOnEdge;
				newTriangles[3].pt[2] = tempTriangle.pt[2];
				newTriangles[3].positionEnum[2] = tempTriangle.positionEnum[2];
				newTriangles[3].calcBounds();

				tempEraseTris.push_back(&workTriangle);

				// 添加必须放到后面，否则地址改变
				tempNewWorkTris.push_back(newTriangles[0]);
				tempNewWorkTris.push_back(newTriangles[1]);
				tempNewWorkTris.push_back(newTriangles[2]);
				tempNewWorkTris.push_back(newTriangles[3]);

			}

			if (pt1WithTriangle != PointWithPolygonOnEdge && pt2WithTriangle != PointWithPolygonOnEdge)
			{
				// 校正ABC为通用情况， 点A、C在直线pt1 pt2的两端
				// aX + bY + c = 0 直线方程
				double a, b, c;
				getLinearEquationXZ(pt1, pt2, a, b, c);

				size_t ACTwoSidesOfLinePt1Pt2_AIndex = 100;
				for (size_t i = 0; i < 3; i++)
				{
					double lineEquationValueVertex1 = (a * workTriangle.pt[i].v.x + b * workTriangle.pt[i].v.z + c);
					double lineEquationValueVertex2 = (a * workTriangle.pt[(i+2)%3].v.x + b * workTriangle.pt[(i+2)%3].v.z + c);
					double pointACLinearExpValMulti = lineEquationValueVertex1 * lineEquationValueVertex2;
					// 做了处理，待验证 TODO: 有顶点在pt1 pt2的直线上的情况处理
					if ( pointACLinearExpValMulti < 0 && abs(lineEquationValueVertex1) >= 5.e-10 && abs(lineEquationValueVertex2) >= 5.e-10 )
					{
						ACTwoSidesOfLinePt1Pt2_AIndex = i;
						break;
					}
				}
				if (ACTwoSidesOfLinePt1Pt2_AIndex > 3)
				{
					printf("没?姓业搅礁龅阍谥毕?pt1,pt2)两侧\n");
					assert( 0 );
				}
				STWorkTriangle tempTriangle;
				tempTriangle.pt[0] = workTriangle.pt[ACTwoSidesOfLinePt1Pt2_AIndex];
				tempTriangle.positionEnum[0] = workTriangle.positionEnum[ACTwoSidesOfLinePt1Pt2_AIndex];
				tempTriangle.pt[1] = workTriangle.pt[(ACTwoSidesOfLinePt1Pt2_AIndex+1)%3];
				tempTriangle.positionEnum[1] = workTriangle.positionEnum[(ACTwoSidesOfLinePt1Pt2_AIndex+1)%3];
				tempTriangle.pt[2] = workTriangle.pt[(ACTwoSidesOfLinePt1Pt2_AIndex+2)%3];
				tempTriangle.positionEnum[2] = workTriangle.positionEnum[(ACTwoSidesOfLinePt1Pt2_AIndex+2)%3];

				// 校正pt1Pie pt2Pie 为通用情况， pt1Pie距离AC更远
				double distancePt1ToAC, distancePt2ToAC;
				distancePt1ToAC = vtkLine::DistanceToLine((double*)&pt1, (double*)&tempTriangle.pt[0], (double *)&tempTriangle.pt[2]);
				distancePt2ToAC = vtkLine::DistanceToLine((double*)&pt2, (double*)&tempTriangle.pt[0], (double *)&tempTriangle.pt[2]);
				if (distancePt1ToAC > distancePt2ToAC)
				{
					pt1Pie = pt1;
					pt2Pie = pt2;
				}
				else
				{
					pt1Pie = pt2;
					pt2Pie = pt1;
				}

				AddPointNoRepeatReturnIndex(intersections, pt1Pie); //把端点加入交点表
				AddPointNoRepeatReturnIndex(intersections, pt2Pie); //把端点加入交点表


				STWorkTriangle newTriangles[5];

				newTriangles[0].pt[0] = tempTriangle.pt[0];
				newTriangles[0].positionEnum[0] = tempTriangle.positionEnum[0];
				newTriangles[0].pt[1] = tempTriangle.pt[1];
				newTriangles[0].positionEnum[1] = tempTriangle.positionEnum[1];
				newTriangles[0].pt[2] = pt1Pie;
				newTriangles[0].positionEnum[2] = PositionOnEdge;
				newTriangles[0].calcBounds();

				newTriangles[1].pt[0] = tempTriangle.pt[0];
				newTriangles[1].positionEnum[0] = tempTriangle.positionEnum[0];
				newTriangles[1].pt[1] = pt1Pie;
				newTriangles[1].positionEnum[1] = PositionOnEdge;
				newTriangles[1].pt[2] = pt2Pie;
				newTriangles[1].positionEnum[2] = PositionOnEdge;
				newTriangles[1].calcBounds();

				newTriangles[2].pt[0] = tempTriangle.pt[0];
				newTriangles[2].positionEnum[0] = tempTriangle.positionEnum[0];
				newTriangles[2].pt[1] = pt2Pie;
				newTriangles[2].positionEnum[1] = PositionOnEdge;
				newTriangles[2].pt[2] = tempTriangle.pt[2];
				newTriangles[2].positionEnum[2] = tempTriangle.positionEnum[2];
				newTriangles[2].calcBounds();

				newTriangles[3].pt[0] = pt1Pie;
				newTriangles[3].positionEnum[0] = PositionOnEdge;
				newTriangles[3].pt[1] = tempTriangle.pt[1];
				newTriangles[3].positionEnum[1] = tempTriangle.positionEnum[1];
				newTriangles[3].pt[2] = tempTriangle.pt[2];
				newTriangles[3].positionEnum[2] = tempTriangle.positionEnum[2];
				newTriangles[3].calcBounds();

				newTriangles[4].pt[0] = pt2Pie;
				newTriangles[4].positionEnum[0] = PositionOnEdge;
				newTriangles[4].pt[1] = pt1Pie;
				newTriangles[4].positionEnum[1] = PositionOnEdge;
				newTriangles[4].pt[2] = tempTriangle.pt[2];
				newTriangles[4].positionEnum[2] = tempTriangle.positionEnum[2];
				newTriangles[4].calcBounds();


				// 清除原三角面片
				tempEraseTris.push_back(&workTriangle);

				// 添加必须放到后面，否则地址改变.
				tempNewWorkTris.push_back(newTriangles[0]);
				tempNewWorkTris.push_back(newTriangles[1]);
				tempNewWorkTris.push_back(newTriangles[2]);
				tempNewWorkTris.push_back(newTriangles[3]);
				tempNewWorkTris.push_back(newTriangles[4]);
			}
		}
	}

}