#pragma once


#ifdef WIN32

    // Under windows avoid including <windows.h> is overrated.
    // Sure, it can be avoided and "name space pollution" can be
    // avoided, but why? It really doesn't make that much difference
    // these days.
    #define  WIN32_LEAN_AND_MEAN
    #include <windows.h>

    #ifndef __gl_h_
		#ifdef _USE_FTGL_OPENGL
			#include <GL/gl.h>
		#endif
        #include <GL/glu.h>
    #endif

#else

    // Non windows platforms - don't require nonsense as seen above :-)
    #ifndef __gl_h_
        #ifdef SDL_main
			#ifdef _USE_FTGL_OPENGL
				#include "SDL_opengl.h"
			#endif
        #elif __APPLE_CC__
			#ifdef _USE_FTGL_OPENGL
				#include <OpenGL/gl.h>
			#endif
            #include <OpenGL/glu.h>
        #else
			#ifdef _USE_FTGL_OPENGL
				#include <GL/gl.h>
			#endif
            #if defined (__sun__) && !defined (__sparc__)
                #include <mesa/glu.h>
            #else
                #include <GL/glu.h>
            #endif
        #endif

    #endif

    // Required for compatibility with glext.h style function definitions of
    // OpenGL extensions, such as in src/osg/Point.cpp.
    #ifndef APIENTRY
        #define APIENTRY
    #endif
#endif
















#define UNITSubdivideDoublelicationDelta (-2e-5)  // ϸ�������κ��غϵ㵥λƫ��


// vtk include
#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkCommand.h>
#include <vtkSTLReader.h>
#include <vtkSTLWriter.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkPLYReader.h>
#include <vtkPLYWriter.h>
#include <vtkPNGWriter.h>
#include <vtkSphereSource.h>
#include <vtkGraphicsFactory.h>
//#include <vtkImagingFactory.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkPlane.h>
#include <vtkProperty.h>
#include <vtkPolyDataMapper.h>
#include <vtkAxesActor.h>
#include <vtkActor.h>
#include <vtkTriangle.h>
#include <vtkDecimatePro.h>
#include <vtkQuadricDecimation.h>
#include <vtkTransform.h>
#include <vtkMatrix4x4.h>
#include <vtkTriangleFilter.h>
#include <vtkWindowToImageFilter.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkOrientationMarkerWidget.h>
#include <vtkAppendPolyData.h>
#include <vtkClipPolyData.h>
#include <vtkPolyDataNormals.h>
#include <vtkLight.h>
#include <vtkPolyDataAlgorithm.h>
#include <vtkObjectFactory.h>
#include <vtkDataObject.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkCell.h>
#include <vtkMath.h>
#include <vtkLine.h>
#include <vtkTriangle.h>
#include <vtkCleanPolyData.h>
#include <vtkFeatureEdges.h>
#include <vtkModifiedBSPTree.h>
#include <vtkMergePoints.h>
#include <vtkPolyLine.h>
#include <vtkPointData.h>

#include <vtkCamera.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>


#include <vector>
#include <list>
#include <stack>
#include <algorithm>



//#if WIN32
#ifndef max
#define max(a,b)            (((a) > (b)) ? (a) : (b))
#endif

#ifndef min
#define min(a,b)            (((a) < (b)) ? (a) : (b))
#endif
//#endif



//////////////////////////////////���㺯��//////////////////////////////////////////////
namespace Shapetizer2
{






}


//////////////////////////////////������//////////////////////////////////////////////
namespace Shapetizer2
{
	typedef struct {
		vtkIdType vId[3];
	} STTriangle;

	typedef union NewLKPoint {
		double d[3];
		struct {double x, y, z;} v;

		inline NewLKPoint(){}
		inline NewLKPoint(union NewLKPoint &pt){ *this = pt; }
		inline NewLKPoint(double x_, double y_, double z_) {
			this->v.x = x_; this->v.y = y_; this->v.z = z_;
		}
		inline NewLKPoint(double data[3]) {
			d[0] = data[0]; d[1] = data[1]; d[2] = data[2];
		}
		inline NewLKPoint(const NewLKPoint &pt) {
			*this = pt;
		}

		inline operator const double *() const {
			return d;
		}

		inline bool operator == (const union NewLKPoint &pt)
		{
			return v.x == pt.v.x && v.z == pt.v.z && v.y == pt.v.y;
		}
	} NewLKPoint;


	typedef std::vector<NewLKPoint> LKRing;


	class LKPolygon
	{
	public:
		LKPolygon(void){}
		~LKPolygon(void){}

		// return (xmin,xmax, ymin,ymax, zmin,zmax).
		double *getBounds()
		{
			if (!rings.empty() && !rings[0].empty())
			{
				bounds[0] = bounds[1] = rings[0][0].v.x;
				bounds[2] = bounds[3] = rings[0][0].v.y;
				bounds[4] = bounds[5] = rings[0][0].v.z;
			}
			else
			{
				bounds[0] = bounds[1] = bounds[2] = bounds[3] = bounds[4] = bounds[5] = 0.0;
			}

			for (size_t i = 0; i < rings.size(); i++)
			{
				LKRing &ring = rings[i];
				for (size_t j = 0; j < ring.size(); j++)
				{
					NewLKPoint &point = ring[j];
					if (bounds[0] > point.v.x) bounds[0] = point.v.x;
					if (bounds[2] > point.v.y) bounds[2] = point.v.y;
					if (bounds[4] > point.v.z) bounds[4] = point.v.z;
					if (bounds[1] < point.v.x) bounds[1] = point.v.x;
					if (bounds[3] < point.v.y) bounds[3] = point.v.y;
					if (bounds[5] < point.v.z) bounds[5] = point.v.z;
				}
			}
			return bounds;
		}

		void doTransform(double matrix[16]);

	public:

		std::vector<LKRing> rings; 

		double bounds[6];

	public:
		size_t getPointsCount() {
			size_t count = 0;
			for (int i = 0; i < rings.size(); count += rings[i].size(), i++);
			return count;
		}

		// �жϵ�pt�Ƿ��ڶ�����ڣ�XZ����ɵ�2Dƽ��
		bool pointInsideXZ2D(double *pt) {
			bool  oddNodes = false;
			for (size_t k = 0; k < rings.size(); k++)
			{
				LKRing &ring = rings[k];
				for (int i = 0, j = ring.size()-1;i<ring.size(); i++)
				{
					if((ring[i].v.z< pt[2] && ring[j].v.z>=pt[2]
					||   ring[j].v.z<pt[2] && ring[i].v.z>=pt[2])
						&& (ring[i].v.x<=pt[0] || ring[j].v.x<=pt[0]))
					{
						oddNodes ^= (ring[i].v.x+(pt[2]-ring[i].v.z)/(ring[j].v.z-ring[i].v.z)*(ring[j].v.x-ring[i].v.x)<pt[0]);
					}
					j = i;
				}
			}
			return oddNodes; 
		}
	};






	// ��ɸѡ��������
	enum PointPositionWithPolygon {
		PositionInside = 0,
		PositionOnEdge = 1,
		PositionOutside = 2
	};
	//�� ��ֱ�����귽��
	typedef struct STWorkLineStruct {
		NewLKPoint pt[2];
		double a, b, c;//ֱ�����귽�� a^2 + b^2 == 1
		double minX, maxX, minZ, maxZ;
		STWorkLineStruct(){}
		STWorkLineStruct(const NewLKPoint &pt1, const NewLKPoint &pt2)
		{
			pt[0] = pt1;
			pt[1] = pt2;

			// ֱ�����귽�� a^2 + b^2 == 1
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

			// ��Χ��
			minX = min(pt1.v.x, pt2.v.x);
			minZ = min(pt1.v.z, pt2.v.z);
			maxX = max(pt1.v.x, pt2.v.x);
			maxZ = max(pt1.v.z, pt2.v.z);
		}
	} STWorkLine;

	typedef struct STWorkTriangleStruct{
		NewLKPoint pt[3];
		STWorkLine line[3];
		PointPositionWithPolygon positionEnum[3];
		double minX, minZ, maxX, maxZ;

		inline bool operator==(const struct STWorkTriangleStruct &b)
		{
			return this == &b;
		}

		void calcBounds() {
			minX = maxX = pt[0].v.x;
			minZ = maxZ = pt[0].v.z;
			for (size_t i = 1; i < 3; i++)
			{
				if (minX > pt[i].v.x)
					minX = pt[i].v.x;
				if (minZ > pt[i].v.z)
					minZ = pt[i].v.z;
				if (maxX < pt[i].v.x)
					maxX = pt[i].v.x;
				if (maxZ < pt[i].v.z)
					maxZ = pt[i].v.z;
			}

			line[0] = STWorkLine(pt[0], pt[1]);
			line[1] = STWorkLine(pt[1], pt[2]);
			line[2] = STWorkLine(pt[2], pt[0]);
		}
	} STWorkTriangle;


	// �������Ķ���α߽����������
	struct PMTriangle;;
	//struct BoundEdge;
	struct BoundEdge {
		bool isLeaf; //�Ƿ�ΪҶ�ӱ�
		BoundEdge *subdivides[2];
		vtkIdType ptId;	// �߶εڶ����˵�id
		vtkIdType ptIdFirst; // �߶ε�һ���˵�id
		PMTriangle *pmTri;
		char edgeIndex;
		bool sameClockwiseAsTriangle; // �����Ƿ�������εķ�����ͬ
	};
	struct PMTriangle {
		int level;
		bool isLeaf;
		PMTriangle *subdivides[2];
		vtkIdType ids[3];
		BoundEdge *boundEdges[3];
		double edgeLen2s[3];
		char longestEdgeId;
	};


	typedef std::vector<BoundEdge *> BoundsEdges;

	/// Ѱ��mesh�߽�������
	struct FBSegment {
		bool added;
		vtkIdType ptids[2];
	};

	// ��������
	template<typename T, unsigned int C>
	class Array {
	public:
		T value[C];

		inline bool operator == (const Array<T, C> &b)
		{
			bool isEqual = true;
			for (unsigned int i = 0; i < C; i++)
			{
				if ( value[i] != b.value[i] )
				{
					isEqual = false;
					break;
				}
			}
			return isEqual;
		}
	};

}






//////////////////////////////////�����ۺ���//////////////////////////////////////////////
namespace Shapetizer2
{
	struct PolygonMesh {
		std::vector<NewLKPoint> points;
		std::vector<char>	pointIsGrays;	//����Ϊ�յĻ���ʾ����ɫ
		std::vector<double>	picture_pointsYOffset;
		std::vector<PMTriangle *> cellNodes;
		std::vector<BoundsEdges> boundEdgeNodess;

		///Ѱ��mesh�߽�����Find Boundary
		std::vector<std::vector<FBSegment *> > fb_ptid_segmentArrs;
		vtkSmartPointer<vtkPolyData> fb_fontPolydata;
		vtkSmartPointer<vtkIncrementalPointLocator> fb_locator;

		std::stack<PMTriangle *> subdividingTri; //�ȴ�ϸ�ֵ�������
		///���ڼ���ϸ�������ζ���
		vtkSmartPointer<vtkPoints> subdivideMergePts;
		vtkSmartPointer<vtkIncrementalPointLocator> subdivideLocator;
		//��ֹϸ�ֶ����غ��ã�Ѱ���غϵ�
		std::vector<Array<vtkIdType, 2> > subdividePtParentSegments;	// ϸ�ֵ�ĸ��߶ζ˵�����ֵ  a)�����߶θ�ֵ(-1, -1) 
		std::vector<int> doublelicationMap;	// index - (key �� value)������XZƽ����ͶӰ�غ�(value = -1 or index)
		std::vector<int> nDeltaYSubdividePts;	// �غϵ����y��ƫ�Ƶ�λ��������λ1e-5)




		// ��ȡ���������ֵİ�Χ��
		void getFontBounds(double bounds[6]);


		//�Ƿ��ڱ߿���
		bool pointInsideXZ2D(double pt[3]);


		/// ̫խ��������ɾ�� �����н�С��minAngleDegree
		void cleanNarrowTriangles(double minAngleDegree);


		/// ����mesh��Ե 
		void findBoundEdgess();

		// Ѱ�ұ߽绷
		void findBoundEdgeRing(FBSegment *firstSegment);


		PolygonMesh()
		{
		}

		~PolygonMesh()
		{
			// �ͷ��ڴ�
			for (unsigned int i = 0; i < cellNodes.size(); i++)
			{
				deletePMTriangle(cellNodes[i]);
			}
			for (unsigned int i = 0; i < boundEdgeNodess.size(); i++)
			{
				BoundsEdges &boundEdgeNodes = boundEdgeNodess[i];
				for (int j = 0; j < boundEdgeNodes.size(); j++)
				{
					deleteBoundEdge(boundEdgeNodes[j]);
				}
			}
			// �ͷ� ��-�߶�
			for (unsigned int i = 0; i < fb_ptid_segmentArrs.size(); i++)
			{
				std::vector<FBSegment *> &segmentArrs = fb_ptid_segmentArrs[i];
				for (unsigned int j = 0; j < segmentArrs.size(); j++)
				{
					FBSegment *fbSegment = segmentArrs[j];
					if (fbSegment->added)
					{
						fbSegment->added = false;
					}
					else
					{
						delete fbSegment;	
					}
				}
			}
		}
		void deletePMTriangle(PMTriangle *pmTri)
		{
			if (pmTri->isLeaf)
			{
				delete pmTri;
			}
			else
			{
				deletePMTriangle(pmTri->subdivides[0]);
				deletePMTriangle(pmTri->subdivides[1]);
				delete pmTri;
			}
		}
		void deleteBoundEdge(BoundEdge *boundEdge)
		{
			if (boundEdge->isLeaf)
			{
				delete boundEdge;
			}
			else
			{
				deleteBoundEdge(boundEdge->subdivides[0]);
				deleteBoundEdge(boundEdge->subdivides[1]);
				delete boundEdge;
			}
		}
	};
	// ���ڶ���ε�����
	struct PtRingIds {
		short firstId;
		short secondId;
		unsigned int ptId;
	};




	//�̶���ε���
	class PolygonCarve
	{
	public:
		PolygonCarve(void);
		~PolygonCarve(void);


	public:
		// ���ò�����ʽһ
		void setModelPolyData(vtkSmartPointer<vtkPolyData> modelPolyData_){ modelPolyData = modelPolyData_; }
		void setFontPolygon(LKPolygon &fontPolygon_){ isPicCarving = false; fontPolygon = fontPolygon_; double *wordBounds = fontPolygon.getBounds(); fontPolygonY = (wordBounds[2] + wordBounds[3]) / 2; } // ���ͷŴ��ܴ�
		void setDepth(double depth){wordDepth = depth;}
		void setSubdivideLen(double subdivideLen){ _subdivideLen2 = subdivideLen * subdivideLen;} // ��Сϸ�������α߳�����


		
		/**
		 * ��ʼ����
		 *     �������������update()
		 * @return �Ƿ���ֳɹ�
		 */
		bool update();
		
		// ��ȡ����ģ��
		vtkSmartPointer<vtkPolyData> getResult(){return result;}

		// ��ȡ������������Ŀ
		int getUnderCellsNum() {return underCellsNum;}

		// ��ȡ��Χ��
		void getFontBounds(double bounds[6]);

		/**
		  * ��ȡ�������
		  *
		  �����б�
			1 - �����ֳ���ģ�ͷ�Χ
			2 - �����ڵ�����
			3 - �����ļ�����ʧ��
			4 - �����ִ����ڰ�͹��ƽ
			5 - ���ֿ���ģ���ཻ
		  */
		const int getErrCode() { return errCode; }
		// ��ȡ������Ϣ
		const char *getError(){ return error.c_str();}



	protected:
		int errCode;
		std::string error;

		bool isPicCarving;

		vtkSmartPointer<vtkPolyData> modelPolyData;
		vtkSmartPointer<vtkPolyData> modelPolyDataOrigin;
		LKPolygon fontPolygon;
		vtkSmartPointer<vtkPolyData> result;
		double wordDepth;
		double wordBounds[6];

		PolygonMesh fontPolygonMesh;	// ����
		double fontPolygonY;		// ��������yֵ
		double _subdivideLen2;		// ��Сϸ�������α߳����ȵ�ƽ��
		int underCellsNum;			// ���ֵ������
		vtkSmartPointer<vtkModifiedBSPTree> bspTree; //���Խ���������û��ϸ�֣�ϸ�ֺ�Ҫ����������
		vtkSmartPointer<vtkPolyData> bspPolyData;


		vtkSmartPointer<vtkTransform> allTransformToVerticality;

		// ���ؽ���Ķ��㣬��������
		//std::vector<STPoint> resultPoints;
		std::vector<STTriangle> resultTriangles;
		vtkSmartPointer<vtkPoints> resultPoints;
		//vtkCellArray *resultCellArray;
		//vtkSmartPointer<vtkUnsignedCharArray> resultPointRGBs;
		std::vector<unsigned char> resultPointRGBs;
#define RGB_COLOR_UNDER_GRAY 0.61
#define RGB_COLOR_UNDER {RGB_COLOR_UNDER_GRAY*0.73*255, RGB_COLOR_UNDER_GRAY*0.73*255, RGB_COLOR_UNDER_GRAY*0.73*255}	//setting: �ײ���ɫ
#define RGB_COLOR_NORMAL {RGB_COLOR_UNDER_GRAY*0.85*255, RGB_COLOR_UNDER_GRAY*0.85*255, RGB_COLOR_UNDER_GRAY*0.85*255}	//setting: ������ɫ


	protected:
		std::vector<STWorkTriangle> workingTriangles;


		// ������������y������ģ�ͱ��棬��ͶӰ��ģ�ͱ���
		bool polygonPastToModel();
		// ��pt��ģ���ϵ�ͶӰ��
		bool pastedToModelFace(NewLKPoint &point);


		// �������ֵ���
		void buildUnderSuface();
		// ϸ����������ӵ����
		void addSubdivideTriToResult(PMTriangle *pmTri);

		//���ɲ�ϸ������
		void buildAndSubdivideUnderSuface();
		// ��Ӳ�ϸ�������ε�����
		void addUnderTriangle(NewLKPoint pt0, NewLKPoint pt1, NewLKPoint pt2);
		// ϸ��������
		void subdivideTriangle(PMTriangle *pmTri);

		// �����µ��������������
		void addBoundEdgeToNewPolygon(BoundEdge *boundEdge);
	};

}
