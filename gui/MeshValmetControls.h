/* Author: Christine Xu, University of North Carolina at Chapel Hill*/

#ifndef MESHVALMETCONTROLS_H
#define MESHVALEMTCONTROLS_H
#include <qdir.h>
#include <qstringlist.h>
#include <qprogressdialog.h>
#include <qmenubar.h>
#include <qpopupmenu.h>

#include "vtkQtRenderWindow.h"
#include "vtkQtRenderWindowInteractor.h"

#include <MeshValmetGUI.h>

#include "vtkActor.h"
#include "vtkBYUReader.h"
#include "vtkCamera.h"
#include "vtkCellArray.h"
#include "vtkColorTransferFunction.h"
#include "vtkCoordinate.h"
#include "vtkDataSet.h"
#include "vtkDoubleArray.h"
#include "vtkFloatArray.h"
#include "vtkLookupTable.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkProperty.h"
#include "vtkProperty2D.h"
#include "vtkRenderer.h"
#include "vtkVRMLExporter.h"
#include "vtkXYPlotActor.h"

#include <3dmodel.h>
#include <read_model.h>
#include <mesh_run.h>
#include <TextWidget.h>



class MeshValmetControls : public MeshValmetGUI
{


public:
  QMenuBar* menubar;
  QPopupMenu *output,*view, *comp;
  int idtext,idhist,idcomphist,idsyn,idoverlay,iderror,idverror,idstat, idreset, idbg;

public:
  QString commonPath;
  QString fileName1;
  QString fileName2;
  
  vtkRenderer* ren1;
  vtkRenderer* ren2;
  vtkRenderer* ren;
  vtkQtRenderWindowInteractor *iren;
  vtkQtRenderWindowInteractor *iren1;
  vtkQtRenderWindowInteractor *iren2;
  
  vtkPolyData* mesh1;
  vtkPolyData* mesh2;
  vtkActor* modelActor1;
  vtkActor* modelActor11;
  vtkActor* modelActor2;
  vtkActor* modelActor22;
  
  vtkColorTransferFunction* lut;

  vtkXYPlotActor * xyplot;
  
  int model1Line;
  int model2Line;
  
  TextWidget* textOut;
  struct outbuf *log;
  
  struct model_error model1,model2;
  int downsampling;
  struct args pargs;
  struct dist_surf_surf_stats stats;
  struct dist_surf_surf_stats stats_rev;
  double abs_sampling_step,abs_sampling_dens;
    

  int *histogram1;
  int *histogram2;
  int *hist;
  int cmap_len;
  double dmax, dmin;
  double middle;
  
  double hausdorff_error;
  double abs_sum_error, abs_sum_sqr_error, mad_error, msd_error, dist_volume;
  //mean_error is vertex mean; sigma_error is vertex sigma(rms).
  double mean_error,sum_error,sum_square_deviations, abs_sum_square_deviations, sigma_error, abs_sigma_error, error_max;
  //face_mean_error is mean of face error; face_rms_eroor is face sigma(rms).
  double face_mean_error, face_rms_error;
  //Median; 95percentile; 68precentile
  double median, percentile95, percentile75, percentile68, percentile25;    
  double dice_coefficient[1], int_union_ratio[1];
  
  bool computed;
  int bgcolor[3];
private:
  struct look
      {
          float R;
          float G;
          float B;
      };
  look lookuptable[256];
public:
  MeshValmetControls( QWidget * parent = 0, const char* name = 0);
  ~MeshValmetControls();
  void CheckTextWindow();
  void OutputRawErrors();
  void OutputVertexErrors();
  void OutputStats();
  void SaveRawErrors(QString outFile);
  void SaveVertexErrors(QString outFile);
  void SaveStats(QString outFile);

  void Model1Open();
  void Model2Open();
  void Model1Switch();
  void Model2Switch();
  void SwitchSync();
  void SwitchOverlay();
  void ResetCameras();
  void SwitchBackground();
  void Compute();
  void SelectStep(const QString& step);
  vtkPolyData* BuildMeshFromModel(model* currModel);
  void SetupWindow(vtkPolyData* mesh, vtkRenderer* ren, vtkActor* actor);
  void MeshSetup(struct args* pargs);
  void SaveModel1();
  void SaveModel2();
  void SaveMeshToIV(vtkQtRenderWindow* win, QString name);
  void drawVertexErrorT(vtkPolyData* mesh, model_error* model, vtkActor* actor, vtkQtRenderWindow* win, int tag = 0);//0: colormap 1: no colormap
  void ChangeSamplingSlider();
  void ColormapCheckbox1Changed(bool checked);
  void ColormapCheckbox2Changed(bool checked);
  void AboutMeshValmet();
  
  void ComputeHistogram();
  void UpdateHistogram();
  void OutputHistogram();
  void SaveHist(QString filename);
  void DoHistogram();
  void UpdateColor();
  void middleSliderMoved(int pos);
  void DoStat();
  void OutputStatInfo();

};

extern "C" {
  void QT_prog(void *out, int p);
}

#endif
