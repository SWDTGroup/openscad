VTK_INCLUDEPATH += -I/usr/include/vtk-5.8
LIBS +=  \
-lvtkalglib \
-lvtkCharts \
-lvtkCommon \
-lvtkDICOMParser \
-lvtkexoIIc \
-lvtkFiltering \
-lvtkftgl \
-lvtkGenericFiltering \
-lvtkGeovis \
-lvtkGraphics \
-lvtkHybrid \
-lvtkImaging \
-lvtkInfovis \
-lvtkIO \
-lvtkmetaio \
-lvtkParallel \
-lvtkproj4 \
-lvtkRendering \
-lvtksys \
-lvtkverdict \
-lvtkViews \
-lvtkVolumeRendering \
-lvtkWidgets
QMAKE_CXXFLAGS += $$VTK_INCLUDEPATH
