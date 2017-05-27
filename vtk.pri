QMAKE_CXXFLAGS += -I/usr/local/include/vtk-5.8
QMAKE_LFLAGS +=  -Wl,--copy-dt-needed-entries
LIBS +=  -L/usr/local/lib/vtk-5.8\
-lvtkCharts	\
-lvtkGeovis	\
-lvtkViews	\
-lvtkInfovis	\
-lvtkWidgets	\
-lvtkHybrid	\
-lvtkVolumeRendering	\
-lvtkRendering	\
-lvtkIO	\
-lvtkGenericFiltering	\
-lvtkGraphics	\
-lvtkImaging	\
-lvtkFiltering	\
-lvtkCommon	\
-lvtkftgl	\
-lvtkalglib	\
-lvtkexoIIc	\
-lvtksqlite	\
-lvtkmetaio	\
-lvtkNetCDF_cxx	\
-lvtkNetCDF	\
-lvtkverdict	\
-lMapReduceMPI	\
-lmpistubs	\
-lvtkproj4	\
-lvtkDICOMParser	\
-lvtklibxml2	\
-lvtkfreetype	\
-lvtkexpat	\
-lvtktiff	\
-lvtkpng	\
-lvtkjpeg	\
-lvtkhdf5	\
-lvtkzlib	\
-lvtksys	\






