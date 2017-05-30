#include "PolygonCarve.h"

// #include <glutess.h>
// #include <tess.h>


#undef _TEST_DEMO


namespace Shapetizer2
{
	// ��ӵ㲢���ص������, ֱ�Ӽӵ㣬���޳��ظ��ĵ�
	static inline vtkIdType AddPointReturnIndex(vtkPoints *points, NewLKPoint &point, bool isGrayColor = false, std::vector<unsigned char> *resultPointRGBsPtr = NULL)	// isGrayColor-�Ƿ�Ϊ����ɫ���㣨falseΪ��ɫ��
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
	// ��ӵ㲢���ص�������� û���ظ���
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

	// �жϵ��Ƿ����������ڣ�2D��
	// ����ֵ�����壺
	//   0 �ڶ������
	//   1 ���ڶ�����
	//   2 ���ڱ���
	//   3 ���ڶ������
	enum PointWithPolygon{
		PointWithPolygonOutside = 0,
		PointWithPolygonInside = 1,
		PointWithPolygonOnVertex = 2,
		PointWithPolygonOnEdge = 3
	};
	static inline PointWithPolygon pointWithWorkTriangle(const NewLKPoint &pt, const STWorkTriangle &workTriangle, size_t &index)
	{
		// �жϵ��ڶ�����
		for (int i=0;i<3; i++)
		{
			const double delta = 2.e-10;
			//if (x == workTriangle.pt[i].x && y == workTriangle.pt[i].z)
			if ( abs(pt.v.x - workTriangle.pt[i].v.x) < delta && abs(pt.v.z - workTriangle.pt[i].v.z) < delta )
			{
				index = i;
				return PointWithPolygonOnVertex;	//���ڶ����ϣ� ����ǳ����������û�д���
			}
		}

		// �жϵ���ˮƽ����
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
				/* ������ * /
				double expValue = pt.x * workTriangle.line[j].a + pt.z * workTriangle.line[j].b + workTriangle.line[j].c;

				oddNodes ^= (expValue > 0);

				/* */
				//a. ����˳�򽻻���ͬ�����
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

				// ����߷ǳ�������ֱ�����ڶ�����
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

	// ��ֱ�߷��� a*X + b * Z + c = 0
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



	// �ж����߶��Ƿ��ཻ
	//   �����ڶ˵�����Ϊfalse
	static inline bool isSegmentIntersectXZ(const STWorkLine &line1, const STWorkLine &line2)
	{
		const double &a = line1.a, &b = line1.b, &c = line1.c;
		const double &A = line2.a, &B = line2.b, &C = line2.c;
		const NewLKPoint &p1 = line1.pt[0],
			&p2 = line1.pt[1],
			&p3 = line2.pt[0],
			&p4 = line2.pt[1];


		// �˵���ֱ�߽����ҵ��ڰ�Χ���ڣ����ཻ
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

			return exp_p1 * exp_p2 < 0 //(ע��4�㹲�ߵ�����)
				&& max(line1.minX, line2.minX) <= min(line1.maxX, line2.maxX)
				&& max(line1.minZ, line2.minZ) <= min(line1.maxZ, line2.maxZ);
		}
		else 
		{
			return false;
		}

	}


	//�����߶εĽ��� http://blog.csdn.net/dgq8211/article/details/7952825
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
		// p1 p2����Ϊͨ�����-p1��p2�����£���������Ľ��㲻ͬ ��p3, p4��Ϊpt1 pt2 ����˳�򲻻��
		bool notNeedSwap_pt1_pt2 = pt1.v.x < pt2.v.x || (pt1.v.x==pt2.v.x && pt1.v.z < pt2.v.z);
		const NewLKPoint& p1 = notNeedSwap_pt1_pt2 ? pt1 : pt2;
		const NewLKPoint& p2 = notNeedSwap_pt1_pt2 ? pt2 : pt1;
		//const STPoint& p1 = pt1;
		//const STPoint& p2 = pt2;

		// fArea(p1,p2,p4)����Ϊ0
		//double k = fArea(p1,p2,p3) / fArea(p1,p2,p4);
		//STPoint intersection;
		//intersection.x = (p3.x + k*p4.x)/(1+k);
		//intersection.z = (p3.z + k*p4.z)/(1+k);

		double s1 = fArea(p1,p2,p3);
		double s2 = fArea(p1,p2,p4);
		NewLKPoint intersection;
		intersection.v.x = (p4.v.x*s1+p3.v.x*s2)/(s1+s2);
		intersection.v.z = (p4.v.z*s1+p3.v.z*s2)/(s1+s2);

		// ����?��y��ֵ
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
	
	//�ж��߶��Ƿ�������ཻ��XZƽ���ϣ�
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

	// �Ƿ��ڶ������
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



	// triangulate ���ǻ����������������ڲ����������������������ζ�����
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

		// �任
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
			//// ɾ���ڱ߽�������� ���� "w"��ĸ���⴦��
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
			//			//printf("�Ѿ������쳣������������ڱ߽�����һ��\n");
			//			ring.erase( find(ring.begin(), ring.end(), pt2) );
			//			ring.erase( find(ring.begin(), ring.end(), pt2Next) );
			//		}
			//	}
			//}






			/// �������ֵ���
			{
#ifdef _TEST_DEMO
				BruceLiu::TimingName timingName("�������ֵ���");
				BruceLiu::TimingBlock timingBlock(timingName);
#endif
				buildAndSubdivideUnderSuface();
			}
		}


		// �������ƽ���Χ����
		if (!isPicCarving)
		{
			double *fontPolygonBounds = fontPolygon.getBounds();
			memcpy(wordBounds, fontPolygonBounds, sizeof(double) * 6);
			wordBounds[2] = wordBounds[3] = this->fontPolygonY;
		}

		/// ������ģ�ͱ��� + ɸѡ������
		{
#ifdef _TEST_DEMO
			BruceLiu::TimingName timingName("��ģ�ͱ���.");
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

			// ��ʼ��bsp��
			bspTree = vtkSmartPointer<vtkModifiedBSPTree>::New();
			bspTree->SetDataSet(bspPolyData);

			// ������ģ�ͱ���
			bool pastSuccessed = polygonPastToModel();
			if (!pastSuccessed)
			{
				errCode = 1;
				error = "�����ֳ���ģ�ͷ�Χ��";
				//assert(0);
				return false;
			}
		}


		/// ɸѡ����Ҫ��ƴ�Ӵ����������Ƭ
		{
#ifdef _TEST_DEMO
			BruceLiu::TimingName timingName("ɸѡ������Ƭ.");
			BruceLiu::TimingBlock timingBlock(timingName);
#endif
			// ������ְ�Χ��
			//double *wordBounds = fontPolygon.getBounds();
			//wordBounds[2] = wordBounds[3] = this->fontPolygonY;	//ע�⣬�˴�����Ҫ
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

			/// ɸѡҪ�и��������
			// ����һ���������������
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

				// ƫ��
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
					// û�н���
				}

				intersectCells->Delete();
				intersectPoints->Delete();
			}

			if (downmostCellId < 0)
			{
				errCode = 1;
				error = "�����ֳ���ģ�ͷ�Χ��";
				return false;
			}


			// �ݹ����ڱߣ��ڱ߰�Χ��������ཻ�ģ��жϷ�����������y�Ḻ��������
			//bool *isCellNeedAddToWorks = (bool *)calloc(modelPolyData->GetNumberOfCells(), sizeof(bool));// ������Ƭ�Ƿ���Ҫ��ӵ��и����
			//bool *isCellCheckedNeighbers = (bool *)calloc(modelPolyData->GetNumberOfCells(), sizeof(bool));// ������Ƭ�Ƿ���Ҫ��ӵ��и����
			bool *isCellNeedAddToWorks = (bool *)calloc(bspPolyData->GetNumberOfCells(), sizeof(bool));// ������Ƭ�Ƿ���Ҫ��ӵ��и����
			bool *isCellCheckedNeighbers = (bool *)calloc(bspPolyData->GetNumberOfCells(), sizeof(bool));// ������Ƭ�Ƿ���Ҫ��ӵ��и����
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
								// ��� 5 ���ֿ���ģ���ཻ
								if (ABC_Min_Y < this->fontPolygonY)
								{
									//����
									errCode = 5;
									error = "���ֿ���ģ���ཻ";
									return false;
								}

 								// ��� 4 ����ճ�������ڰ�͹��ƽ
 								// ��֤�����η�����yֵ�Ƿ����0
 								NewLKPoint vectorAB, vectorAC;
 								vectorAB.v.x = B[0] - A[0];
 								vectorAB.v.z = B[2] - A[2];
 								vectorAC.v.x = C[0] - A[0];
 								vectorAC.v.z = C[2] - A[2];
 								double normalY = vectorAC.v.x * vectorAB.v.z - vectorAB.v.x * vectorAC.v.z;
 								if (normalY > 0)
 								{
 									// �ж����ڵ��������Ƿ���bounds�ཻ
 									double pt1[3], pt2[3];
 									//modelPolyData->GetPoint(edgeCell->GetPointId(0), pt1);
 									//modelPolyData->GetPoint(edgeCell->GetPointId(1), pt2); 
 									bspPolyData->GetPoint(edgeCell->GetPointId(0), pt1);
 									bspPolyData->GetPoint(edgeCell->GetPointId(1), pt2); 
 									if (isSegmentCrossRectangeXZ(pt1, pt2, wordBounds))
 									{
 										// ����
 										errCode = 4;
 										error = "�����ִ����ڰ�͹��ƽ";
 										return false;
 									}
 									else
 									{
 										// ��Ҫ��������
 										needAddToWorks = false;
 									}
 
 								}
 								else
 								{
 									needAddToWorks = true;
 								}
 							}
 
 							// �������ֿ��ڵ�������Ƭid��ӵ�������
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

 			// ����ɸѡ������Ƭ
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
 					// ��ӵ�����������Ƭ����
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


		// ���и����
		{
#ifdef _TEST_DEMO
			BruceLiu::TimingName timingName("�и����");
			BruceLiu::TimingBlock timingBlock(timingName);
#endif

			std::vector<STWorkTriangle> tempNewWorkTris;	// �����½�������Ƭ
			std::vector<STWorkTriangle *> tempEraseTris;	// Ҫɾ����������Ƭ
			std::vector<NewLKPoint> intersections;				// ��������Ƭ�Ľ���
			std::vector<NewLKPoint> boundEdgeSubdividePts;				// �����߶�ϸ�ֵ�
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
						// �ж��߶��Ƿ����Χ���ཻ
						if (
							max(line.minX, workTriangle.minX) > min(line.maxX, workTriangle.maxX)
							|| max(line.minZ, workTriangle.minZ) > min(line.maxZ, workTriangle.maxZ)
							)
						{
							continue;
						}


						// ���˵��������λ�ù�ϵ
						size_t pt1OnIndex = 100, pt2OnIndex = 100;
						PointWithPolygon pt1WithTriangle = pointWithWorkTriangle(pt1, workTriangle, pt1OnIndex);
						PointWithPolygon pt2WithTriangle = pointWithWorkTriangle(pt2, workTriangle, pt2OnIndex);

						// һ���е�����������
						//   1. �����㶼����������
						if (pt1WithTriangle && pt2WithTriangle)
						{
							//���ǻ�
							triangulateFivePoint(workTriangle, line, pt1, pt2, pt1WithTriangle, pt2WithTriangle, pt1OnIndex, pt2OnIndex, tempNewWorkTris, tempEraseTris, intersections);
							break;
						}

						//   2. һ��������������
						if ( !pt1WithTriangle && pt2WithTriangle || pt1WithTriangle && !pt2WithTriangle )
						{
							triangulateCross_1_2(workTriangle, line, pt1, pt2, pt1WithTriangle, pt2WithTriangle, pt1OnIndex, pt2OnIndex, tempNewWorkTris, tempEraseTris, intersections);
							continue;
						}

						// ����û�е�����������
						if (!pt1WithTriangle && !pt2WithTriangle)
						{
							triangulateCross_2(workTriangle, line, pt1, pt2, pt1WithTriangle, pt2WithTriangle, pt1OnIndex, pt2OnIndex, tempNewWorkTris, tempEraseTris, intersections);
							continue;
						}

					}

					//�ѷ�����������Ƭ�������������ɵ�������Ƭ
					//workingTriangles.reserve(workingTriangles.size()+100);
					std::sort(tempEraseTris.begin(), tempEraseTris.end());
					for (int idx = tempEraseTris.size() - 1; idx >= 0; idx--)
					{
						workingTriangles.erase(std::find(workingTriangles.begin(), workingTriangles.end(), *tempEraseTris[idx]));
					}
					workingTriangles.insert(workingTriangles.end(), tempNewWorkTris.begin(), tempNewWorkTris.end());


					/// �������ֱ�Ե��ֱ������Ƭ
					if (intersections.size() < 2)
					{
						fprintf(stderr, "�߶������ص�����ڽӽ�\n");
					}
					else
					{
						NewLKPoint vectorPt1Pt2;
						vectorPt1Pt2.v.x = pt2.v.x - pt1.v.x;
						vectorPt1Pt2.v.z = pt2.v.z - pt1.v.z;
						// �˵���뽻��� ���ﲻ����ӣ��˵�����������ζ���ϲ���������
						//AddPointNoRepeatReturnIndex(intersections, pt1); //�Ѷ˵���뽻���
						//AddPointNoRepeatReturnIndex(intersections, pt2); //�Ѷ˵���뽻���
						// ��������
						{
							int j, k, h;
							NewLKPoint tempPt; 
							for (h=intersections.size()-1; h>0; h=k) /*ѭ����û�бȽϷ�Χ*/
							{
								for (j=0, k=0; j<h; j++) /*ÿ��Ԥ��k=0��ѭ��ɨ������k*/
								{
									if ( (intersections[j].v.x - intersections[j+1].v.x) * vectorPt1Pt2.v.x > 1.e-16
										|| ( abs((intersections[j].v.x - intersections[j+1].v.x) * vectorPt1Pt2.v.x) <= 1.e-16
										&& (intersections[j].v.z - intersections[j+1].v.z) * vectorPt1Pt2.v.z > 0 )
										) /*��ķ��ں��棬С�ķŵ�ǰ��*/
									{
										tempPt = intersections[j];
										intersections[j] = intersections[j+1];
										intersections[j+1] = tempPt; /*��ɽ���*/
										k = j; /*��������³���λ�á�����k����Ķ��������ź��˵ġ�*/
									}
								}
							}
						}
						// ���������
						NewLKPoint pt1UnderSurface = pt1, pt2UnderSurface = pt2;
						pt1UnderSurface.v.y += pt1_offsetY;
						pt2UnderSurface.v.y += pt2_offsetY;
						// ģ�ͱ����������Ӷ���� ��ƬȺ
						for (int idx = 1; idx < intersections.size(); idx++)
						{
							STTriangle newTriangle;
							newTriangle.vId[0] = AddPointReturnIndex(resultPoints, pt1UnderSurface);
							newTriangle.vId[1] = AddPointReturnIndex(resultPoints, intersections[idx-1]);
							newTriangle.vId[2] = AddPointReturnIndex(resultPoints, intersections[idx]);
							resultTriangles.push_back(newTriangle);
						}
						//// ����벿��
						//if (!intersections.empty())
						//{
						//	STTriangle newTriangle;
						//	newTriangle.vId[0] = AddPointReturnIndex(resultPoints, pt1UnderSurface);
						//	newTriangle.vId[1] = AddPointReturnIndex(resultPoints, pt2UnderSurface);
						//	newTriangle.vId[2] = AddPointReturnIndex(resultPoints, intersections.back());
						//	resultTriangles.push_back(newTriangle);
						//}
						// �����������Ӷ������ƬȺ
						if (isPicCarving) {
							STTriangle newTriangle;
							newTriangle.vId[0] = AddPointReturnIndex(resultPoints, intersections.back());
							newTriangle.vId[1] = AddPointReturnIndex(resultPoints, pt2UnderSurface);
							newTriangle.vId[2] = AddPointReturnIndex(resultPoints, pt1UnderSurface);
							resultTriangles.push_back(newTriangle);
						}
						else {
							boundEdgeSubdividePts.push_back(pt1UnderSurface);
							if (fontPolygonMesh.boundEdgeNodess[i][j]->pmTri) // ���յ�ϸ�ֵı�
							{
								std::stack<BoundEdge *> waitingRecurseAddEdge;	// �ȴ������edge
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
							else	//û�н��յ�ϸ�ֵıߣ��ٴ�ϸ��
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


		// workingTriangle�и����ϲ�����������
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


		/// �������ֵ���
		buildUnderSuface();



		/// ����polydata
		// create hash table and merge points/triangles.
		{
#ifdef _TEST_DEMO
			BruceLiu::TimingName timingName("�ϲ�����");
			BruceLiu::TimingBlock timingBlock(timingName);
#endif

			//��ģ�Ͷ�����ɫ����
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
			//������ɫ
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

			//��ɫ
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

		// ��ģ����ת��ԭ��ԭ״
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

		// ������������ �� ƫ�Ƶ����棨���棩
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

			// ƫ��
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
				// ƫ��
				//point.y  -= wordDepth;
				// �غϵ�ƫ��
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

		// ƫ��
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


	// ϸ����������ӵ����
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

		//�ݹ����������
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
			int idA = (idB + 2) % 3; //ABΪ���
			int idC = (idB + 1) % 3;
			NewLKPoint ptA = fontPolygonMesh.points[pmTri->ids[idA]]; //ע����ֵ��Ӧ�ĵ�ַ��ı�
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
			//ע����ֵ��Ӧ�ĵ�ַ��ı�
			//fontPolygonMesh.points.push_back(newPt);
			//vtkIdType ptsBeforeSize = subdivideMergePts->GetNumberOfPoints();
			fontPolygonMesh.subdivideLocator->InsertUniquePoint((double*)&newPt, newPtFirstOneIndex);
			newPtRealIndex = newPtFirstOneIndex;
			bool isRealNewPt = newPtFirstOneIndex == fontPolygonMesh.points.size();
			//��ɫ����
			if (isRealNewPt && fontPolygonMesh.pointIsGrays.size()>0)
			{
				char isGrayA = fontPolygonMesh.pointIsGrays[pmTri->ids[idA]]; //ע����ֵ��Ӧ�ĵ�ַ��ı�
				char isGrayB = fontPolygonMesh.pointIsGrays[pmTri->ids[idB]];
				char isGrayC = fontPolygonMesh.pointIsGrays[pmTri->ids[idC]];
				char newPtIsGray = isGrayA && isGrayB;
				fontPolygonMesh.pointIsGrays.push_back(newPtIsGray);
			}
			// �õ㸸�߶�
			Array<vtkIdType, 2> newPtParentSegment;
			newPtParentSegment.value[0] = min(pmTri->ids[idA], pmTri->ids[idB]);
			newPtParentSegment.value[1] = max(pmTri->ids[idA], pmTri->ids[idB]);
			bool isSamePtButParent = false;
			if (!isRealNewPt)
			{
				if (!(newPtParentSegment == fontPolygonMesh.subdividePtParentSegments[newPtRealIndex]))
				{
					if (fontPolygonMesh.doublelicationMap[newPtFirstOneIndex] == -1 ) // ��չΪ��㹲λ��ʱ���Ĵ˴�
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
						assert(fontPolygonMesh.doublelicationMap[newPtRealIndex] == newPtFirstOneIndex);// ��㹲λ�û�������
					}
				}
			}
			if (isRealNewPt)
			{
				fontPolygonMesh.points.push_back(newPt);
				fontPolygonMesh.subdividePtParentSegments.push_back(newPtParentSegment);
				if (isSamePtButParent)
				{
					// ƫ��������
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
				/// �߽��ߴ���
				newTri->boundEdges[0] = pmTri->boundEdges[idA]; 
				newTri->boundEdges[1] = (BoundEdge *)NULL;
				newTri->boundEdges[2] = (BoundEdge *)NULL;
				if (pmTri->boundEdges[idA])
				{
					pmTri->boundEdges[idA]->pmTri = newTri;
				}
				// ϸ�ֱ߽���
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
					//����֦�ڵ��ӽڵ㸳ֵ
					boundEdge->subdivides[boundEdge->sameClockwiseAsTriangle?0:1] = newBoundEdge;
					//�߽縳������������
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
				//�ݹ�ϸ��������
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
				/// �߽��ߴ���
				newTri->boundEdges[0] = (BoundEdge *)NULL; 
				newTri->boundEdges[1] = (BoundEdge *)NULL; 
				newTri->boundEdges[2] = pmTri->boundEdges[idC];
				if (pmTri->boundEdges[idC])
				{
					pmTri->boundEdges[idC]->pmTri = newTri;
				}
				// ϸ�ֱ߽���
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
					//����֦��?��ӽڵ㸳�?
					boundEdge->subdivides[boundEdge->sameClockwiseAsTriangle?1:0] = newBoundEdge;
					//�߽縳������������
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
				//�ݹ�ϸ��������
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
		// ��ʼ�� ��-ѡ��
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

		// ���ӱ߽��߶�
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
		// �ȸ���һ���㸳ֵ˳��ֵ
		ptRingIds.secondId = -1;	//ע�⣬������պ�Ϊ-1
		ptRingIds.ptId = curPtid;
		points[curPtid].v.y = *(double *)&ptRingIds;
		//����ring
		boundEdgeNodess.resize(boundEdgeNodess.size()+1);
		BoundsEdges &boudEdges = boundEdgeNodess.back();
		do 
		{
			nextPtid = curSegment->ptids[curSameAsTriangle ? 1 : 0];


			// �㸳��˳��ֵ
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


			///// �ظ��ĵ�����¶���  ��Ҫ�Ķ��ϴ�a����������ʱҪ    b�������ζ�������ֵ
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

			//	// ����µ㣬ƫ����+1
			//	if (addedNum >= 2)
			//	{
			//		STPoint newPt = points[curPtid];
			//		// ����µ�
			//		vtkIdType newPtid = points.size();
			//		points.push_back(newPt);
			//		newPt.y = 0;
			//		subdivideLocator->InsertNextPoint((double *)&newPt);
			//		//ƫ����+1
			//		nDeltaYSubdividePts.push_back(nDeltaYSubdividePts[curPtid]+1);
			//		// ����ӳ��
			//		doublelicationMap.push_back(curPtid);
			//		doublelicationMap[curPtid] = newPtid;
			//		// ���߶θ�ֵ
			//		Array<vtkIdType, 2> pSegment;
			//		pSegment.value[0] = pSegment.value[1] = -1;
			//		subdividePtParentSegments.push_back(pSegment);

			//		//��-�߶� �����
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
				// ��������
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
					deletedCell[cellId] = true;// ɾ����ǰcell
					
					int idA = biggerAngleIndex,
						idB = (biggerAngleIndex+1)%3,
						idC = (biggerAngleIndex+2)%3;// BC��ΪҪ�и�ı�
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

					// �ڱ�
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

						deletedCell[neighborCellId] = true;// ɾ��������cell

						// ����polydata������һ�ּ��
						int fb_fontPolydataCellNum = fb_fontPolydata->GetNumberOfCells();
						for (int cellId2 = 0; cellId2 < fb_fontPolydataCellNum; cellId2++)
						{
							if (!deletedCell[cellId2])
							{
								fb_fontPolydata->GetCellPoints(cellId2, npts, pts);
								newCells->InsertNextCell(3, pts);
							}
						}

						/// ����������������
						vtkIdType neighborNpts, *neighborPts;
						fb_fontPolydata->GetCellPoints(neighborCellId, neighborNpts, neighborPts);
						int nIdA = -1, nIdB, nIdC; // ABΪ�ڱ�
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
						// ������������
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
						goto LoopFilter; // ע�⣬ʹ����goto
						//break;
					}
					else
					{
						//printf("no neighbor \n");

						// ����polydata������һ�ּ��
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
						goto LoopFilter; // ע�⣬ʹ����goto
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


		// ��ʼ���߶�
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
			// ����������
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
					fprintf(stderr, "û�п��ǵ�polygon type : %i \n", polygonType);
					assert(0);
				}
			}
		}



		/// ̫խ��������ɾ�� �����н�С��
		fontPolygonMesh.cleanNarrowTriangles(0.2);





		/// Ũ�������ʼ��
		fontPolygonMesh.subdivideMergePts = vtkSmartPointer<vtkPoints>::New();
		fontPolygonMesh.subdivideLocator.TakeReference(vtkMergePoints::New());
		double *polygonBounds1 = fontPolygon.getBounds();
		polygonBounds1[2] = 0;
		polygonBounds1[3] = 0;
		fontPolygonMesh.subdivideLocator->InitPointInsertion(fontPolygonMesh.subdivideMergePts, polygonBounds1);


		// ��ʼ����
		vtkIdType ptNum = fontPolygonMesh.fb_fontPolydata->GetNumberOfPoints();
		//points.reserve(ptNum * 1.01);
		for (int i = 0; i < ptNum; i++)
		{
			double *pt = fontPolygonMesh.fb_fontPolydata->GetPoint(i);
			fontPolygonMesh.points.push_back(*(NewLKPoint *)pt);
		}
		// ��ӵ㵽Ũ������
		for (int i = 0; i < fontPolygonMesh.points.size(); i++)
		{
			// ��ӵ�
			NewLKPoint insertPt = fontPolygonMesh.points[i];
			insertPt.v.y = 0;
			fontPolygonMesh.subdivideLocator->InsertNextPoint((double *)&insertPt);
			//fontPolygonMesh.points.push_back(insertPt);
			this->fontPolygonMesh.nDeltaYSubdividePts.push_back(0);
			// ����ӳ��
			fontPolygonMesh.doublelicationMap.push_back(-1);
			// ���߶θ�ֵ
			Array<vtkIdType, 2> pSegment;
			pSegment.value[0] = pSegment.value[1] = -1;
			fontPolygonMesh.subdividePtParentSegments.push_back(pSegment);
		}



		// mesh���Ҳ������߽�
		fontPolygonMesh.findBoundEdgess();


		// ���ϸ�������Σ���������߽�����ϵ
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
				//����Ҫmerge�ĵ�
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
				// ����α�Ե���
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
			// ��ӵ�ϸ�������ζ���
			if (longestEdge2 > _subdivideLen2)
			{
				fontPolygonMesh.subdividingTri.push(pmTri);
				//this->subdivideTriangle(pmTri);
			}
		}




		/// ϸ��������
		while (!fontPolygonMesh.subdividingTri.empty())
		{
			PMTriangle *pmTri = fontPolygonMesh.subdividingTri.top();
			fontPolygonMesh.subdividingTri.pop();
			subdivideTriangle(pmTri);
		}



		///��ԭ����yֵ
		for (size_t i = 0; i < fontPolygonMesh.points.size(); i++)
		{
			fontPolygonMesh.points[i].v.y = this->fontPolygonY;
		}
	}

	// ����û�е�����������
	void triangulateCross_2(STWorkTriangle &workTriangle, STWorkLine line, NewLKPoint &pt1, NewLKPoint &pt2,
		PointWithPolygon &pt1WithTriangle, PointWithPolygon &pt2WithTriangle,
		size_t &pt1OnIndex, size_t &pt2OnIndex,
		std::vector<STWorkTriangle> &tempNewWorkTris,
		std::vector<STWorkTriangle *> &tempEraseTris,
		std::vector<NewLKPoint> &intersections)
	{
		//   �ж��Ƿ����������н���
		//�ж�pt1-pt2 ��A-B B-C C-A�Ƿ��н���
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

		//   ����н�������3���������Σ��޳�һ��������
		if (intersectionCount == 2)
		{
			// �ѱ任��ͨ�������pt1-pt2 ���� A-B�ཻ
			// ����pt1Intersect����B-C�Ľ��㣬 pt2IntersectΪC-A
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

			// �󽻵�����
			NewLKPoint pt1Intersect = tempIsOnEdge[1] ? (isPt1OnEdge[1]?pt1:pt2) : segmentIntersetionXZ(tempTriangle.pt[1], tempTriangle.pt[2], pt1, pt2);
			NewLKPoint pt2Intersect = tempIsOnEdge[2] ? (isPt1OnEdge[2]?pt1:pt2) : segmentIntersetionXZ(tempTriangle.pt[2], tempTriangle.pt[0], pt1, pt2);

			// ��ӽ��㵽�б�
			AddPointNoRepeatReturnIndex(intersections, pt1Intersect);
			AddPointNoRepeatReturnIndex(intersections, pt2Intersect);

			// �г�3����������
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

			// ��ӱ���ŵ����棬�����ַ�ı�
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
			fprintf(stderr, "�������������⣬��������%i������ֱ��̫����������ĳ���߻�ĳ������ \n", intersectionCount);
			assert(0);
		}
		else if (intersectionCount == 3)
		{
			fprintf(stderr, "�������������⣬��������%i������̫����������ĳ���� \n", intersectionCount);
			assert(0);
		}

	}


	//   2. һ��������������
	void triangulateCross_1_2(STWorkTriangle &workTriangle, STWorkLine line, NewLKPoint &pt1, NewLKPoint &pt2,
		PointWithPolygon &pt1WithTriangle, PointWithPolygon &pt2WithTriangle,
		size_t &pt1OnIndex, size_t &pt2OnIndex,
		std::vector<STWorkTriangle> &tempNewWorkTris,
		std::vector<STWorkTriangle *> &tempEraseTris,
		std::vector<NewLKPoint> &intersections)
	{
		if (pt1WithTriangle == PointWithPolygonOnVertex && pt2WithTriangle == PointWithPolygonOnVertex)
		{
			//do nothing �����㶼�������ζ����ϣ�û�����������
			fprintf(stderr, "���������������⣬���㶼�������ζ�����\n");
			assert(0);
		}
		else if (pt1WithTriangle == PointWithPolygonOnVertex || pt2WithTriangle == PointWithPolygonOnVertex)
		{
			size_t ptOnTriVertexIndex = pt1WithTriangle == PointWithPolygonOnVertex ? pt1OnIndex : pt2OnIndex;
			// ��������ת��Ϊͨ�õ�������: pt1Pie�������ζ���B��
			size_t offset_0 = (ptOnTriVertexIndex+2)%3;
			size_t offset_1 = (ptOnTriVertexIndex+0)%3;
			size_t offset_2 = (ptOnTriVertexIndex+1)%3;

			NewLKPoint pt1Pie = workTriangle.pt[offset_1];

			//�ж�pt1-pt2 �� C-A�ཻ
			bool isIntersect1 = isSegmentIntersectXZ(workTriangle.line[offset_2], line);

			if (isIntersect1)
			{
				NewLKPoint ptIntersect = segmentIntersetionXZ(workTriangle.pt[offset_2], workTriangle.pt[offset_0], pt1, pt2);

				// ��ӽ��㵽�б�
				AddPointNoRepeatReturnIndex(intersections, ptIntersect);

				AddPointNoRepeatReturnIndex(intersections, pt1Pie); //�Ѷ˵���뽻���

				// ��������ɵ�������Ƭ
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

				// ��ӱ���ŵ����棬�����ַ�ı�
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
			//�ж�pt1-pt2 ��A-B B-C C-A�Ƿ��н���
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
				fprintf(stderr, "����ֱ����������ڡ��⣬������Ϊ0�� ����ԭ��\n    1.isSegmentIntersectXZ2����������\n    2.�жϵ�������������������  \n");
				assert(0);
			}

			// ����������Ϊ 1���� ��4��������
			if (intersectionCount == 1)
			{
				// ��������ת��Ϊͨ�õ�������
				//   1.pt1-pt2 �� C-A�ཻ
				//   2. pt1Pie����������
				//   3. ptIntersectΪ����
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
				AddPointNoRepeatReturnIndex(intersections, pt1Pie); //�Ѷ˵���뽻���


				// pt1Pie��C-A����
				if (pt1WithTriangle == PointWithPolygonOnEdge || pt2WithTriangle == PointWithPolygonOnEdge)
				{// ��������ɵ�������Ƭ
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

					// ��ӱ���ŵ����棬�����ַ�ı�
					tempNewWorkTris.push_back(newTriangles[0]);
					tempNewWorkTris.push_back(newTriangles[1]);


				}
				// pt1Pie����C-A����
				else {
					// �󽻵�����
					NewLKPoint ptIntersect = segmentIntersetionXZ(tempTriangle.pt[2], tempTriangle.pt[0], pt1, pt2);

					// ��ӽ��㵽�б�
					AddPointNoRepeatReturnIndex(intersections, ptIntersect);

					// ��������ɵ�������Ƭ
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

					// ��ӱ���ŵ����棬�����ַ�ı�
					tempNewWorkTris.push_back(newTriangles[0]);
					tempNewWorkTris.push_back(newTriangles[1]);
					tempNewWorkTris.push_back(newTriangles[2]);
					tempNewWorkTris.push_back(newTriangles[3]);
				}


			}

			// ����������Ϊ 2���� �϶��������ڵĵ��ڱ���
			if (intersectionCount == 2)
			{
				if (pt1WithTriangle != PointWithPolygonOnEdge && pt2WithTriangle != PointWithPolygonOnEdge)
				{
					fprintf(stderr, "���������������⣬û�е��ڱ��ϣ�����Ϊ����������ֱ��ͨ��ĳ������\n");
					assert(0);
				}
				else
				{

					// �ѱ任��ͨ�������pt1-pt2 ���� A-B�ཻ
					// ����pt1Intersect����B-C�Ľ��㣬 pt2IntersectΪC-A
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

					// �󽻵�����
					NewLKPoint pt1Intersect = tempIsOnEdge[1] ? (tempIsPt1OnEdge[1]?pt1:pt2) : segmentIntersetionXZ(tempTriangle.pt[1], tempTriangle.pt[2], pt1, pt2);
					NewLKPoint pt2Intersect = tempIsOnEdge[2] ? (tempIsPt1OnEdge[2]?pt1:pt2) : segmentIntersetionXZ(tempTriangle.pt[2], tempTriangle.pt[0], pt1, pt2);

					// ��ӽ��㵽�б�
					AddPointNoRepeatReturnIndex(intersections, pt1Intersect);
					AddPointNoRepeatReturnIndex(intersections, pt2Intersect);

					// �г�3����������
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

					// ��ӱ���ŵ����棬�����ַ�ı�
					tempNewWorkTris.push_back(newTriangles[0]);
					tempNewWorkTris.push_back(newTriangles[1]);
					tempNewWorkTris.push_back(newTriangles[2]);
				}
			}

			// ����������Ϊ 3�� ����
			if (intersectionCount == 3)
			{
				fprintf(stderr, "����ֱ����������ڡ��⣬������Ϊ3��!!\n");
				assert(0);
			}
		}
	}

	//   1. �����㶼����������
	void triangulateFivePoint(STWorkTriangle &workTriangle, STWorkLine line, NewLKPoint &pt1, NewLKPoint &pt2,
		PointWithPolygon &pt1WithTriangle, PointWithPolygon &pt2WithTriangle,
		size_t &pt1OnIndex, size_t &pt2OnIndex,
		std::vector<STWorkTriangle> &tempNewWorkTris,
		std::vector<STWorkTriangle *> &tempEraseTris,
		std::vector<NewLKPoint> &intersections)
	{
		bool point1OnTriangleVertices = false, point2OnTriangleVertices = false;
		size_t onTriangleVertexIndex;
		NewLKPoint pt1Pie, pt2Pie;//�غϵĵ���pt1Pie��tempTriangle.pt[0]��
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

		//�ж�pt1-pt2 ��A-B B-C C-A�Ƿ��н���
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
			//printf("���㶼�������ζ�����\n");
			if (pt1OnIndex == pt2OnIndex)
			{
				printf("������������ͬһ���㣨%i���� \n", pt1OnIndex);
				//AddPointNoRepeatReturnIndex(intersections, pt1); //�Ѷ˵���뽻���
				//AddPointNoRepeatReturnIndex(intersections, pt2); //�Ѷ˵���뽻���
			}
			else
			{
				//printf("������������������(%i)(%i)�� \n", pt1OnIndex, pt2OnIndex);
				AddPointNoRepeatReturnIndex(intersections, pt1); //�Ѷ˵���뽻���
				AddPointNoRepeatReturnIndex(intersections, pt2); //�Ѷ˵���뽻���
			}
		}
		else if (point1OnTriangleVertices || point2OnTriangleVertices)
		{
			// ͨ����� pt1Pie�ڶ���A��
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

			// �ϲ������ʱ����
			pt1Pie = tempTriangle.pt[0];

			if (intersectionCount == 0)
			{
				AddPointNoRepeatReturnIndex(intersections, pt1Pie); //�Ѷ˵���뽻���
				AddPointNoRepeatReturnIndex(intersections, pt2Pie); //�Ѷ˵���뽻���

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

				// ��ӱ���ŵ����棬�����ַ�ı�
				tempNewWorkTris.push_back(newTriangles[0]);
				tempNewWorkTris.push_back(newTriangles[1]);
				tempNewWorkTris.push_back(newTriangles[2]);
			}
			else if (intersectionCount == 1)
			{
				//printf("���һ.1 \n");
				if (tempIsOnEdge[1])
				{
					printf("    ���һ.1.b \n");

					AddPointNoRepeatReturnIndex(intersections, pt1Pie); //�Ѷ˵���뽻���
					AddPointNoRepeatReturnIndex(intersections, pt2Pie); //�Ѷ˵���뽻���
 
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

  					// ��ӱ���ŵ����棬�����ַ�ı�
    				tempNewWorkTris.push_back(newTriangles[0]);
    				tempNewWorkTris.push_back(newTriangles[1]);

				}
				if (tempIsOnEdge[0])
				{
					printf("    ���һ.1.c \n");

					AddPointNoRepeatReturnIndex(intersections, pt1Pie); //�Ѷ˵���뽻���
					AddPointNoRepeatReturnIndex(intersections, pt2Pie); //�Ѷ˵���뽻���

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

					// ��ӱ���ŵ����棬�����ַ�ı�
					tempNewWorkTris.push_back(newTriangles[0]);
					tempNewWorkTris.push_back(newTriangles[1]);
				}

				if (tempIsOnEdge[2])
				{
					printf("    ���һ.1.d \n");

					AddPointNoRepeatReturnIndex(intersections, pt1Pie); //�Ѷ˵���뽻���
					AddPointNoRepeatReturnIndex(intersections, pt2Pie); //�Ѷ˵���뽻���

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

					// ��ӱ���ŵ����棬�����ַ�ı�
					tempNewWorkTris.push_back(newTriangles[0]);
					tempNewWorkTris.push_back(newTriangles[1]);
				}
			}
		}
		else
		{
			//printf("���һ.1/2/3 \n");

			// �����㶼�ڱ���
			if (pt1WithTriangle == PointWithPolygonOnEdge && pt2WithTriangle == PointWithPolygonOnEdge)
			{
				//printf("���һ.2 \n");
				if (pt1OnIndex != pt2OnIndex)
				{
					printf("    ���һ.2.a \n");

					// �ѱ任��ͨ�������pt1-pt2 ���� A-B�ཻ
					// ����pt1Intersect����B-C�Ľ��㣬 pt2IntersectΪC-A
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

					// �󽻵�����
					NewLKPoint pt1Intersect = tempIsOnEdge[1] ? (isPt1OnEdge[1]?pt1:pt2) : segmentIntersetionXZ(tempTriangle.pt[1], tempTriangle.pt[2], pt1, pt2);
					NewLKPoint pt2Intersect = tempIsOnEdge[2] ? (isPt1OnEdge[2]?pt1:pt2) : segmentIntersetionXZ(tempTriangle.pt[2], tempTriangle.pt[0], pt1, pt2);

					// ��ӽ��㵽�б�
					AddPointNoRepeatReturnIndex(intersections, pt1Intersect);
					AddPointNoRepeatReturnIndex(intersections, pt2Intersect);

					// �г�3����������
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

					// ��ӱ���ŵ����棬�����ַ�ı�
					tempNewWorkTris.push_back(newTriangles[0]);
					tempNewWorkTris.push_back(newTriangles[1]);
					tempNewWorkTris.push_back(newTriangles[2]);

				}
				else if (pt1OnIndex == pt2OnIndex)
				{
					printf("    ���һ.2.b \n");

					// �ѱ任��ͨ�������pt1-pt2 �� C-A�ϣ�pt1����A
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

					// pt1 pt2λ�õ���
					NewLKPoint ACVector, pt1pt2Vector;
					ACVector.v.x = tempTriangle.pt[2].v.x - tempTriangle.pt[0].v.x;
					ACVector.v.z = tempTriangle.pt[2].v.z - tempTriangle.pt[0].v.z;
					pt1pt2Vector.v.x = pt2.v.x - pt1.v.x;
					pt1pt2Vector.v.z = pt2.v.z - pt1.v.z;
					// ACVector�������н�Լ����pt1pt2Vector�ģ�����λ��
					// ����ͬ����
					bool notNeedSwap = ACVector.v.x * pt1pt2Vector.v.x > 1.e-16 || (abs(ACVector.v.x * pt1pt2Vector.v.x) <= 1.e-16 && ACVector.v.z * pt1pt2Vector.v.z > 0);
					pt1Pie = notNeedSwap ? pt1 : pt2;
					pt2Pie = notNeedSwap ? pt2 : pt1;

 					AddPointNoRepeatReturnIndex(intersections, pt1Pie); //�Ѷ˵���뽻���
 					AddPointNoRepeatReturnIndex(intersections, pt2Pie); //�Ѷ˵���뽻���

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

					// ��ӱ���ŵ����棬�����ַ�ı�
					tempNewWorkTris.push_back(newTriangles[0]);
					tempNewWorkTris.push_back(newTriangles[1]);
					tempNewWorkTris.push_back(newTriangles[2]);
				}
			}

			// һ�����ڱ���
			if (pt1WithTriangle == PointWithPolygonOnEdge ^ pt2WithTriangle == PointWithPolygonOnEdge)
			{
				//printf("���һ.1 \n");
				printf("    ���һ.1.a \n");

				// ��������ת��Ϊͨ�õ�������
				//   1.pt1-pt2 �� C-A�ཻ
				//   2. pt1Pie����������
				//   3. ptIntersectΪ����
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
				AddPointNoRepeatReturnIndex(intersections, pt1Pie); //�Ѷ˵���뽻���

				// �󽻵�����
				NewLKPoint ptIntersect = pt1WithTriangle == PointWithPolygonOnEdge ? pt1 : pt2;

				// ��ӽ��㵽�б�
				AddPointNoRepeatReturnIndex(intersections, ptIntersect);

				// ��������ɵ�������Ƭ
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

				// ��ӱ���ŵ����棬�����ַ�ı�
				tempNewWorkTris.push_back(newTriangles[0]);
				tempNewWorkTris.push_back(newTriangles[1]);
				tempNewWorkTris.push_back(newTriangles[2]);
				tempNewWorkTris.push_back(newTriangles[3]);

			}

			if (pt1WithTriangle != PointWithPolygonOnEdge && pt2WithTriangle != PointWithPolygonOnEdge)
			{
				// У��ABCΪͨ������� ��A��C��ֱ��pt1 pt2������
				// aX + bY + c = 0 ֱ�߷���
				double a, b, c;
				getLinearEquationXZ(pt1, pt2, a, b, c);

				size_t ACTwoSidesOfLinePt1Pt2_AIndex = 100;
				for (size_t i = 0; i < 3; i++)
				{
					double lineEquationValueVertex1 = (a * workTriangle.pt[i].v.x + b * workTriangle.pt[i].v.z + c);
					double lineEquationValueVertex2 = (a * workTriangle.pt[(i+2)%3].v.x + b * workTriangle.pt[(i+2)%3].v.z + c);
					double pointACLinearExpValMulti = lineEquationValueVertex1 * lineEquationValueVertex2;
					// ���˴�������֤ TODO: �ж�����pt1 pt2��ֱ���ϵ��������
					if ( pointACLinearExpValMulti < 0 && abs(lineEquationValueVertex1) >= 5.e-10 && abs(lineEquationValueVertex2) >= 5.e-10 )
					{
						ACTwoSidesOfLinePt1Pt2_AIndex = i;
						break;
					}
				}
				if (ACTwoSidesOfLinePt1Pt2_AIndex > 3)
				{
					printf("û?��ҵ���������ֱ�?pt1,pt2)����\n");
					assert( 0 );
				}
				STWorkTriangle tempTriangle;
				tempTriangle.pt[0] = workTriangle.pt[ACTwoSidesOfLinePt1Pt2_AIndex];
				tempTriangle.positionEnum[0] = workTriangle.positionEnum[ACTwoSidesOfLinePt1Pt2_AIndex];
				tempTriangle.pt[1] = workTriangle.pt[(ACTwoSidesOfLinePt1Pt2_AIndex+1)%3];
				tempTriangle.positionEnum[1] = workTriangle.positionEnum[(ACTwoSidesOfLinePt1Pt2_AIndex+1)%3];
				tempTriangle.pt[2] = workTriangle.pt[(ACTwoSidesOfLinePt1Pt2_AIndex+2)%3];
				tempTriangle.positionEnum[2] = workTriangle.positionEnum[(ACTwoSidesOfLinePt1Pt2_AIndex+2)%3];

				// У��pt1Pie pt2Pie Ϊͨ������� pt1Pie����AC��Զ
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

				AddPointNoRepeatReturnIndex(intersections, pt1Pie); //�Ѷ˵���뽻���
				AddPointNoRepeatReturnIndex(intersections, pt2Pie); //�Ѷ˵���뽻���


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


				// ���ԭ������Ƭ
				tempEraseTris.push_back(&workTriangle);

				// ��ӱ���ŵ����棬�����ַ�ı�.
				tempNewWorkTris.push_back(newTriangles[0]);
				tempNewWorkTris.push_back(newTriangles[1]);
				tempNewWorkTris.push_back(newTriangles[2]);
				tempNewWorkTris.push_back(newTriangles[3]);
				tempNewWorkTris.push_back(newTriangles[4]);
			}
		}
	}

}