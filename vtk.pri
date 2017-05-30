
# read environment variables
OPENSCAD_LIBRARIES_DIR = $$(OPENSCAD_LIBRARIES)

!isEmpty(OPENSCAD_LIBRARIES_DIR) {
        VTK6_INCLUDEPATH = $$OPENSCAD_LIBRARIES_DIR/include/vtk-6.3
       VTK6_LIBPATH = $$OPENSCAD_LIBRARIES_DIR/lib
 }
else
{
        VTK6_INCLUDEPATH = /usr/local/include/vtk-6.3
        VTK6_LIBPATH = /usr/local/lib 
}


CONFIG(mingw-cross-env) {
  #message("mingw vtk")
    MXE_TARGET_DIR=$$(MXETARGETDIR)
    #message($$MXE_TARGET_DIR)
       VTK6_INCLUDEPATH = $$MXE_TARGET_DIR/include/vtk-6.3
      VTK6_LIBPATH = $$MXE_TARGET_DIR/lib

}



QMAKE_CXXFLAGS += -I$$VTK6_INCLUDEPATH
QMAKE_LFLAGS +=  -Wl,--copy-dt-needed-entries
LIBS +=  -L$$VTK6_LIBPATH \
-lvtkDomainsChemistry-6.3 \
-lvtkFiltersFlowPaths-6.3 \
-lvtkFiltersGeneric-6.3 \
-lvtkFiltersHyperTree-6.3 \
-lvtkFiltersParallelImaging-6.3 \
-lvtkFiltersProgrammable-6.3 \
-lvtkFiltersSelection-6.3 \
-lvtkFiltersSMP-6.3 \
-lvtkFiltersTexture-6.3 \
-lvtkFiltersVerdict-6.3 \
-lvtkverdict-6.3 \
-lvtkGeovisCore-6.3 \
-lvtkproj4-6.3 \
-lvtkImagingMath-6.3 \
-lvtkImagingMorphological-6.3 \
-lvtkImagingStatistics-6.3 \
-lvtkImagingStencil-6.3 \
-lvtkInteractionImage-6.3 \
-lvtkIOAMR-6.3 \
-lvtkIOEnSight-6.3 \
-lvtkIOExodus-6.3 \
-lvtkIOExport-6.3 \
-lvtkRenderingGL2PS-6.3 \
-lvtkRenderingContextOpenGL-6.3 \
-lvtkIOImport-6.3 \
-lvtkIOInfovis-6.3 \
-lvtklibxml2-6.3 \
-lvtkIOLSDyna-6.3 \
-lvtkIOMINC-6.3 \
-lvtkIOMovie-6.3 \
-lvtkIOParallel-6.3 \
-lvtkIOParallelXML-6.3 \
-lvtkIOPLY-6.3 \
-lvtkIOSQL-6.3 \
-lvtksqlite-6.3 \
-lvtkIOVideo-6.3 \
-lvtkRenderingImage-6.3 \
-lvtkRenderingLIC-6.3 \
-lvtkRenderingLOD-6.3 \
-lvtkRenderingVolumeOpenGL-6.3 \
-lvtkViewsContext2D-6.3 \
-lvtkViewsInfovis-6.3 \
-lvtkFiltersAMR-6.3 \
-lvtkgl2ps-6.3 \
-lvtkexoIIc-6.3 \
-lvtkFiltersParallel-6.3 \
-lvtkIONetCDF-6.3 \
-lvtkNetCDF_cxx-6.3 \
-lvtkNetCDF-6.3 \
-lvtkParallelCore-6.3 \
-lvtkIOXML-6.3 \
-lvtkIOGeometry-6.3 \
-lvtkjsoncpp-6.3 \
-lvtkIOXMLParser-6.3 \
-lvtkexpat-6.3 \
-lvtkIOLegacy-6.3 \
-lvtkRenderingOpenGL-6.3 \
-lvtkChartsCore-6.3 \
-lvtkRenderingContext2D-6.3 \
-lvtkFiltersImaging-6.3 \
-lvtkInfovisLayout-6.3 \
-lvtkInfovisCore-6.3 \
-lvtkViewsCore-6.3 \
-lvtkInteractionWidgets-6.3 \
-lvtkFiltersHybrid-6.3 \
-lvtkImagingGeneral-6.3 \
-lvtkImagingSources-6.3 \
-lvtkFiltersModeling-6.3 \
-lvtkImagingHybrid-6.3 \
-lvtkIOImage-6.3 \
-lvtkDICOMParser-6.3 \
-lvtkIOCore-6.3 \
-lvtkmetaio-6.3 \
-lvtkpng-6.3 \
-lvtktiff-6.3 \
-lvtkjpeg-6.3 \
-lvtkInteractionStyle-6.3 \
-lvtkRenderingAnnotation-6.3 \
-lvtkImagingColor-6.3 \
-lvtkRenderingVolume-6.3 \
-lvtkRenderingLabel-6.3 \
-lvtkRenderingFreeType-6.3 \
-lvtkRenderingCore-6.3 \
-lvtkCommonColor-6.3 \
-lvtkFiltersExtraction-6.3 \
-lvtkFiltersStatistics-6.3 \
-lvtkalglib-6.3 \
-lvtkImagingFourier-6.3 \
-lvtkImagingCore-6.3 \
-lvtkFiltersGeometry-6.3 \
-lvtkFiltersSources-6.3 \
-lvtkFiltersGeneral-6.3 \
-lvtkFiltersCore-6.3 \
-lvtkCommonExecutionModel-6.3 \
-lvtkCommonComputationalGeometry-6.3 \
-lvtkCommonDataModel-6.3 \
-lvtkCommonMisc-6.3 \
-lvtkCommonTransforms-6.3 \
-lvtkCommonMath-6.3 \
-lvtkCommonSystem-6.3 \
-lvtkCommonCore-6.3 \
-lvtksys-6.3 \
-lvtkftgl-6.3 \
-lvtkfreetype-6.3 \
-lvtkzlib-6.3 