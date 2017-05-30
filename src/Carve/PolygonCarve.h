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
















#define UNITSubdivideDoublelicationDelta (-2e-5)  // 细化三角形后重合点单位偏移


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



//////////////////////////////////计算函数//////////////////////////////////////////////
namespace Shapetizer2
{






}


//////////////////////////////////基本类//////////////////////////////////////////////
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

		// 判断点pt是否在多边形内，XZ轴组成的2D平面
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






	// 待筛选的三角面
	enum PointPositionWithPolygon {
		PositionInside = 0,
		PositionOnEdge = 1,
		PositionOutside = 2
	};
	//线 含直线坐标方程
	typedef struct STWorkLineStruct {
		NewLKPoint pt[2];
		double a, b, c;//直线坐标方程 a^2 + b^2 == 1
		double minX, maxX, minZ, maxZ;
		STWorkLineStruct(){}
		STWorkLineStruct(const NewLKPoint &pt1, const NewLKPoint &pt2)
		{
			pt[0] = pt1;
			pt[1] = pt2;

			// 直线坐标方程 a^2 + b^2 == 1
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

			// 包围盒
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


	// 点索引的多边形边界和三角网格
	struct PMTriangle;;
	//struct BoundEdge;
	struct BoundEdge {
		bool isLeaf; //是否为叶子边
		BoundEdge *subdivides[2];
		vtkIdType ptId;	// 线段第二个端点id
		vtkIdType ptIdFirst; // 线段第一个端点id
		PMTriangle *pmTri;
		char edgeIndex;
		bool sameClockwiseAsTriangle; // 方向是否和三角形的方向相同
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

	/// 寻找mesh边界多边形用
	struct FBSegment {
		bool added;
		vtkIdType ptids[2];
	};

	// 定长数组
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






//////////////////////////////////主体综合类//////////////////////////////////////////////
namespace Shapetizer2
{
	struct PolygonMesh {
		std::vector<NewLKPoint> points;
		std::vector<char>	pointIsGrays;	//长度为空的话表示无颜色
		std::vector<double>	picture_pointsYOffset;
		std::vector<PMTriangle *> cellNodes;
		std::vector<BoundsEdges> boundEdgeNodess;

		///寻找mesh边界多边形Find Boundary
		std::vector<std::vector<FBSegment *> > fb_ptid_segmentArrs;
		vtkSmartPointer<vtkPolyData> fb_fontPolydata;
		vtkSmartPointer<vtkIncrementalPointLocator> fb_locator;

		std::stack<PMTriangle *> subdividingTri; //等待细分的三角形
		///用于减少细化三角形顶点
		vtkSmartPointer<vtkPoints> subdivideMergePts;
		vtkSmartPointer<vtkIncrementalPointLocator> subdivideLocator;
		//防止细分顶点重合用，寻找重合点
		std::vector<Array<vtkIdType, 2> > subdividePtParentSegments;	// 细分点的父线段端点索引值  a)轮廓线段赋值(-1, -1) 
		std::vector<int> doublelicationMap;	// index - (key 和 value)的两点XZ平面上投影重合(value = -1 or index)
		std::vector<int> nDeltaYSubdividePts;	// 重合点添加y轴偏移单位个数（单位1e-5)




		// 获取文字立体字的包围盒
		void getFontBounds(double bounds[6]);


		//是否在边框里
		bool pointInsideXZ2D(double pt[3]);


		/// 太窄的三角形删除 两个夹角小于minAngleDegree
		void cleanNarrowTriangles(double minAngleDegree);


		/// 查找mesh边缘 
		void findBoundEdgess();

		// 寻找边界环
		void findBoundEdgeRing(FBSegment *firstSegment);


		PolygonMesh()
		{
		}

		~PolygonMesh()
		{
			// 释放内存
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
			// 释放 点-线段
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
	// 点在多边形的索引
	struct PtRingIds {
		short firstId;
		short secondId;
		unsigned int ptId;
	};




	//刻多边形的类
	class PolygonCarve
	{
	public:
		PolygonCarve(void);
		~PolygonCarve(void);


	public:
		// 设置参数方式一
		void setModelPolyData(vtkSmartPointer<vtkPolyData> modelPolyData_){ modelPolyData = modelPolyData_; }
		void setFontPolygon(LKPolygon &fontPolygon_){ isPicCarving = false; fontPolygon = fontPolygon_; double *wordBounds = fontPolygon.getBounds(); fontPolygonY = (wordBounds[2] + wordBounds[3]) / 2; } // 会释放此能存
		void setDepth(double depth){wordDepth = depth;}
		void setSubdivideLen(double subdivideLen){ _subdivideLen2 = subdivideLen * subdivideLen;} // 最小细分三角形边长长度


		
		/**
		 * 开始刻字
		 *     设置完参数调用update()
		 * @return 是否刻字成功
		 */
		bool update();
		
		// 获取刻字模型
		vtkSmartPointer<vtkPolyData> getResult(){return result;}

		// 获取底面三角形数目
		int getUnderCellsNum() {return underCellsNum;}

		// 获取包围盒
		void getFontBounds(double bounds[6]);

		/**
		  * 获取错误编码
		  *
		  错误列表：
			1 - 有文字超出模型范围
			2 - 不存在的字体
			3 - 字体文件加载失败
			4 - 刻文字处过于凹凸不平
			5 - 文字框与模型相交
		  */
		const int getErrCode() { return errCode; }
		// 获取错误信息
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

		PolygonMesh fontPolygonMesh;	// 底面
		double fontPolygonY;		// 文字轮廓y值
		double _subdivideLen2;		// 最小细分三角形边长长度的平放
		int underCellsNum;			// 文字底面个数
		vtkSmartPointer<vtkModifiedBSPTree> bspTree; //有自交轮廓，边没有细分，细分后要重新贴表面
		vtkSmartPointer<vtkPolyData> bspPolyData;


		vtkSmartPointer<vtkTransform> allTransformToVerticality;

		// 返回结果的顶点，和三角面
		//std::vector<STPoint> resultPoints;
		std::vector<STTriangle> resultTriangles;
		vtkSmartPointer<vtkPoints> resultPoints;
		//vtkCellArray *resultCellArray;
		//vtkSmartPointer<vtkUnsignedCharArray> resultPointRGBs;
		std::vector<unsigned char> resultPointRGBs;
#define RGB_COLOR_UNDER_GRAY 0.61
#define RGB_COLOR_UNDER {RGB_COLOR_UNDER_GRAY*0.73*255, RGB_COLOR_UNDER_GRAY*0.73*255, RGB_COLOR_UNDER_GRAY*0.73*255}	//setting: 底部颜色
#define RGB_COLOR_NORMAL {RGB_COLOR_UNDER_GRAY*0.85*255, RGB_COLOR_UNDER_GRAY*0.85*255, RGB_COLOR_UNDER_GRAY*0.85*255}	//setting: 正常颜色


	protected:
		std::vector<STWorkTriangle> workingTriangles;


		// 文字轮廓沿着y轴贴到模型表面，即投影到模型表面
		bool polygonPastToModel();
		// 点pt在模型上的投影点
		bool pastedToModelFace(NewLKPoint &point);


		// 生成文字底面
		void buildUnderSuface();
		// 细分三角形添加到结果
		void addSubdivideTriToResult(PMTriangle *pmTri);

		//生成并细化底面
		void buildAndSubdivideUnderSuface();
		// 添加并细分三角形到底面
		void addUnderTriangle(NewLKPoint pt0, NewLKPoint pt1, NewLKPoint pt2);
		// 细分三角形
		void subdivideTriangle(PMTriangle *pmTri);

		// 生成新的切碎的文字轮廓
		void addBoundEdgeToNewPolygon(BoundEdge *boundEdge);
	};

}
