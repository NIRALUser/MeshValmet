PROJECT(MeshValmet)

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

FIND_PACKAGE(Qt3 REQUIRED)

FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

INCLUDE_DIRECTORIES(
 ${QT_INCLUDE_DIR}
 ${CMAKE_CURRENT_SOURCE_DIR}
 ${CMAKE_CURRENT_BINARY_DIR}
 ${CMAKE_CURRENT_SOURCE_DIR}/lib3d
 ${CMAKE_CURRENT_SOURCE_DIR}/mesh
 ${CMAKE_CURRENT_SOURCE_DIR}/gui
)

LINK_LIBRARIES(
 ${QT_QT_LIBRARY}
 ${QT_GL_LIBRARY}
 vtkRendering  vtkWidgets
)

SET(MeshValmet_SRCS
 lib3d/3dmodel.h
 lib3d/block_list.h
 lib3d/debug_print.h
 lib3d/geomutils.h
 lib3d/model_in.h
 lib3d/model_in_ply.h
 lib3d/types.h
 mesh/compute_error.h
 mesh/mesh_run.h
 mesh/model_analysis.h
 mesh/reporting.h
 mesh/xalloc.h
 mesh/compute_volume_overlap.h
 gui/MeshValmetControls.h
 gui/vtkQtRenderWindow.h
 gui/vtkQtRenderWindowInteractor.h
 gui/TextWidget.h
 gui/ColormapWidget.h
 lib3d/read_model.h
 
 lib3d/block_list.cxx
 lib3d/geomutils.cxx
 lib3d/model_in.cxx
 lib3d/model_in_ply.cxx
 lib3d/model_in_raw.cxx
 lib3d/model_in_smf.cxx
 lib3d/model_in_vrml_iv.cxx
 
 mesh/xalloc.cxx
 mesh/reporting.cxx
 mesh/model_analysis.cxx
 mesh/compute_error.cxx
 mesh/mesh_run.cxx
 mesh/compute_volume_overlap.cxx
 gui/MeshValmetControls.cxx
 gui/vtkQtRenderWindow.cxx
 gui/vtkQtRenderWindowInteractor.cxx
 gui/TextWidget.cxx
 gui/ColormapWidget.cxx
 MeshValmetMain.cxx
)

SET(MeshValmet_GUI_SRCS 
  gui/MeshValmetGUI.ui
)

SET(MeshValmet_MOC_SRCS
  gui/vtkQtRenderWindowInteractor.h
  gui/ColormapWidget.h
)

IF(QT_WRAP_UI)
  QT_WRAP_UI(MeshValmet MeshValmet_HDRS  MeshValmet_SRCS  ${MeshValmet_GUI_SRCS} )
ENDIF(QT_WRAP_UI)

IF(QT_WRAP_CPP)
  QT_WRAP_CPP(MeshValmet MeshValmet_SRCS ${MeshValmet_MOC_SRCS} )
ENDIF(QT_WRAP_CPP)

ADD_DEFINITIONS(-DQT_DLL)

ADD_DEFINITIONS(-DQT_THREAD_SUPPORT)

ADD_EXECUTABLE(MeshValmet ${MeshValmet_SRCS})
