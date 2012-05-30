/* Author: Christine Xu, University of North Carolina at Chapel Hill*/

#define QT_ALTERNATE_QTSMANIP

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

//#include <assert.h>

#include <qdir.h>
#include <qfiledialog.h> 
#include <qstring.h>
#include <qfileinfo.h>
#include <qmessagebox.h>
#include <qcheckbox.h>
#include <qslider.h>
#include <qspinbox.h>
#include <qradiobutton.h>
#include <qtextstream.h>
#include <qlineedit.h>
#include <qfile.h>
#include <qtoolbutton.h>
#include <qlistbox.h>
#include <qdialog.h>

#include <MeshValmetControls.h>
#include <ColormapWidget.h>

#include <vtkTextProperty.h>
#include <vtkAxisActor2D.h> 
#include <vtkInteractorStyleSwitch.h>
#include <vtkImplicitBoolean.h>
//#include <vtkClipPolyData.h>

#include "geomutils.h"
#include "compute_volume_overlap.h"

int compare_doubles (const void * a, const void * b)
{
  double x = *((double *)a);
  double y = *((double *)b);

  if(x>y)
    return 1;
  else if(x<y)
    return -1;
  else
    return 0;
}

MeshValmetControls::MeshValmetControls( QWidget * parent, const char* name )
  :MeshValmetGUI(parent, name)
{
  
  commonPath = "";
  fileName1 = "";
  fileName2 = "";
  model1Line = 0;
  model2Line = 0;
  stepList->setSelected(2,true);
  
  ren1 = vtkRenderer::New();
  renWin1->AddRenderer(ren1);
  ren2 = vtkRenderer::New();
  renWin2->AddRenderer(ren2);
  ren = vtkRenderer::New();  
  renWin->AddRenderer(ren);
  iren = vtkQtRenderWindowInteractor::New();  
  iren1 = vtkQtRenderWindowInteractor::New();
  iren2 = vtkQtRenderWindowInteractor::New();
  vtkInteractorStyleSwitch *trackBall1 = vtkInteractorStyleSwitch::New(); 
  trackBall1->SetCurrentStyleToTrackballCamera(); 
  iren1->SetInteractorStyle(trackBall1);
  vtkInteractorStyleSwitch *trackBall2 = vtkInteractorStyleSwitch::New(); 
  trackBall2->SetCurrentStyleToTrackballCamera(); 
  iren2->SetInteractorStyle(trackBall2); 
  //iren1->SetInteractorStyle(vtkInteractorStyleSwitch::New());
  //iren2->SetInteractorStyle(vtkInteractorStyleSwitch::New());
  iren1->SetRenderWindow(renWin1);
  iren2->SetRenderWindow(renWin2);
  iren->SetRenderWindow(renWin);

  modelActor1 = vtkActor::New();
  modelActor11 = vtkActor::New();
  modelActor2 = vtkActor::New();
  modelActor22 = vtkActor::New();
        
  downsampling = 1;  
  
  /*When displaying model with wireframe representation, some of the edges
  *are displayed as black color, which is invisible when the background color
  *is also black.
  */
  bgcolor[0] = 1;
  bgcolor[1] = 1;
  bgcolor[2] = 1;
  ren1->SetBackground(1,1,1);
  ren2->SetBackground(1,1,1);
  ren->SetBackground(1,1,1);
  ren->SetViewport(0,0,1,1);        
        
  histogram1 = NULL;
  histogram2 = NULL;
  hist = NULL;
  cmap_len = 256;

  mean_error=0;
  sum_error=0;
  sum_square_deviations=0;
  abs_sum_square_deviations=0;
  sigma_error=0;
  error_max=0;
  face_mean_error=0;
  face_rms_error=0;
  //median1=0;
  //median2=0;
  median=0;
  percentile95=0;
  percentile75=0;
  percentile68=0;
  percentile25=0;
  dice_coefficient[0] = 0.0;
  int_union_ratio[0] = 0.0;
  
  computed = false;
  //RecomputeHistButton->setEnabled(false);
  //outputHistButton->setEnabled(false);
  UpdateColorButton->setEnabled(false);
  UpdateHistButton->setEnabled(false);
  middleSlider->setEnabled(false);
  IntervalSpinBox->setEnabled(false);
  pargs.do_symmetric = 0;
    
  ColormapCheckbox1->setEnabled(false);
  ColormapCheckbox2->setEnabled(false);
  
  
  //SynchronizeButton->setToggleButton(true);
    
  middle = 0;
  for (float i=0;i<256;i++) 
  {
      lookuptable[(int)i].R = 255*((-(i*i)/(float)4096)+(i/32));
      if (lookuptable[(int)i].R<0)lookuptable[(int)i].R =0;
      lookuptable[(int)i].G = 255*((-((i-64)*(i-64))/(float)4096)+((i-64)/(float)32));
      if (lookuptable[(int)i].G<0)lookuptable[(int)i].G =0;
      lookuptable[(int)i].B = 255*((-((i-128)*(i-128))/(float)4096)+((i-128)/(float)32));
      if (lookuptable[(int)i].B<0)lookuptable[(int)i].B =0;
  }
  
  menubar = new QMenuBar( this, "menu" );

  QPopupMenu * file = new QPopupMenu( this );
  menubar->insertItem( "File", file );
  output = new QPopupMenu( this );
  output->setCheckable(true);
  menubar->insertItem( "Output", output );
  view = new QPopupMenu( this );
  view->setCheckable(true);
  menubar->insertItem( "View", view );
  comp = new QPopupMenu( this );
  menubar->insertItem( "Compute", comp );
  QPopupMenu * help = new QPopupMenu( this);
  menubar->insertItem("Help", help);
  
  file->insertItem("Load Model A",this,SLOT(Model1Open()));
  file->insertItem("Load Model B",this,SLOT(Model2Open()));
  file->insertSeparator();
  file->insertItem("Exit",this,SLOT(close()));
  
  output->insertItem("Save Model 1",this,SLOT(SaveModel1()));
  output->insertItem("Save Model 2",this,SLOT(SaveModel2()));
  output->insertSeparator();
  idtext = output->insertItem("External TextWindow",this,SLOT(CheckTextWindow()));
  output->setItemChecked(idtext,false);
  output->insertSeparator();
  idhist = output->insertItem("Output Histogram",this,SLOT(OutputHistogram()));
  output->setItemEnabled(idhist,false);
  output->insertSeparator();
  iderror = output->insertItem("Output Sample Errors",this,SLOT(OutputRawErrors()));
  output->setItemEnabled(iderror,false);
  output->insertSeparator();
  idverror = output->insertItem("Output Vertex Errors",this,SLOT(OutputVertexErrors()));
  output->setItemEnabled(idverror,false);
  output->insertSeparator();
  idstat = output->insertItem("Output Statistical Infos",this,SLOT(OutputStats()));
  output->setItemEnabled(idstat,false);
  
  comp->insertItem("Compute Error",this,SLOT(Compute()));
  comp->insertSeparator();
  idcomphist = comp->insertItem("Update Histogram",this,SLOT(ComputeHistogram()));
  comp->setItemEnabled(idcomphist,false);
  
  idbg = view->insertItem("Switch BLK/WHT Background",this,SLOT(SwitchBackground()));
  view->insertSeparator();
  idreset = view->insertItem("Reset Cameras", this, SLOT(ResetCameras()));
  view->insertSeparator();
  idsyn = view->insertItem("Synchronize Viewport",this, SLOT(SwitchSync()));
  view->setItemEnabled(idsyn,false);
  view->insertSeparator();
  idoverlay = view->insertItem("Overlay Two Views",this, SLOT(SwitchOverlay()));
  view->setItemEnabled(idoverlay,false);
  view->insertSeparator();
  view->insertItem("Fill/Line/Point for Model 1",this, SLOT(Model1Switch()));
  view->insertItem("Fill/Line/Point for Model 2",this, SLOT(Model2Switch()));
  
  help->insertItem("About MeshValmet",this, SLOT(AboutMeshValmet()));

}

MeshValmetControls::~MeshValmetControls()
{
  if(renWin1 != NULL)
  {
    renWin1->Delete();
    renWin1 = NULL;
  }
  if(renWin2 != NULL)
  {
    renWin2->Delete();
    renWin2 = NULL;  
  }
  if(renWin != NULL)
  {
    renWin->Delete();
    renWin = NULL;
  }
  
}

void MeshValmetControls::AboutMeshValmet()
{
  QMessageBox::information(this,"About MeshValmet", "MeshValmet 3.0 \n Author: Christine Xu \n Contact: xushun@gmail.com \n Funded by: Guido Gerig and Martin Styner \n Copyright (c) 2010");
}


void MeshValmetControls::CheckTextWindow()
{
  if(output->isItemChecked(idtext))
    output->setItemChecked(idtext,false);
  else
    output->setItemChecked(idtext,true);
}

vtkPolyData* MeshValmetControls::BuildMeshFromModel(model* currModel)
{
                        
  vtkPolyData *mesh = vtkPolyData::New();
  vtkPoints *points = vtkPoints::New();
  vtkCellArray *polys = vtkCellArray::New();

  vtkIdType i;
  float v[3];
  vtkIdType f[3];
  for (i=0; i<currModel->num_vert; i++)
  {
      vertex_t t = currModel->vertices[i];
      v[0] = t.x;
      v[1] = t.y;
      v[2] = t.z;
      points->InsertPoint(i,v);
  }
  for (i=0; i<currModel->num_faces; i++) 
  {
      face_t t = currModel->faces[i];
      f[0] = t.f0;
      f[1] = t.f1;
      f[2] = t.f2;
      polys->InsertNextCell(3,f);
  }
    
  mesh->SetPoints(points);
  points->Delete();
  mesh->SetPolys(polys);
  polys->Delete();
    
  return mesh;
}

void MeshValmetControls::SetupWindow(vtkPolyData* mesh, vtkRenderer* ren, vtkActor* actor)
{
  vtkPolyDataMapper *meshMapper = vtkPolyDataMapper::New();
  meshMapper->SetInput(mesh);
  actor->SetMapper(meshMapper);
        
  ren->AddActor(actor);
        
  vtkCamera *camera = vtkCamera::New();
  camera->SetPosition(1,1,1);
  camera->SetFocalPoint(0,0,0);
        
  ren->SetActiveCamera(camera);
  ren->ResetCamera();
        
}

void MeshValmetControls::Model1Open()
{
  QString model1file;
  if(commonPath == "")
    model1file = QFileDialog::getOpenFileName(
                    ".",
                    "3D Models (*.raw *.byu *.iv)",//m3d model boundaries, ct, mri, analyze image implied surfaces
                    this,
                    "open file dialog",
                    "Choose a file to open" );
  else
    model1file = QFileDialog::getOpenFileName(
                    commonPath,
                    "3D Models (*.raw *.byu *.iv)",
                    this,
                    "open file dialog",
                    "Choose a file to open" );
  if(model1file != "" && model1file != NULL)
  {
    QFileInfo fi(model1file.latin1());
    commonPath = fi.dirPath();
    model1Name->setText(fi.fileName());
    
    /*if(model1file.right(4).lower() == ".byu")
    {
      vtkBYUReader* byureader = vtkBYUReader::New();
      byureader->SetFileName(model1file.latin1());
      mesh1 = byureader->GetOutput();
      
    }
    else */
    if(model1file.right(4).lower() == ".byu" || model1file.right(4).lower() == ".raw" || model1file.right(3).lower() == ".iv")
    {
      struct model* model1 = read_model_file(model1file.latin1());
      mesh1 = BuildMeshFromModel(model1);
    }
    
    SetupWindow(mesh1,ren1,modelActor1);
    modelActor1->GetProperty()->SetColor(0,0,1);
    
    renWin1->Render();
    
    fileName1 = model1file;  
    
    ColormapCheckbox1->setEnabled(false);
    ColormapCheckbox2->setEnabled(false);
    ColormapCheckbox1->setChecked(false);
    //memset(&model1,0,sizeof(*model1));
    
    computed = false;
    //RecomputeHistButton->setEnabled(false);
    comp->setItemEnabled(idcomphist,false);
    //outputHistButton->setEnabled(false);
    output->setItemEnabled(idhist,false);
    output->setItemEnabled(iderror,false);
    output->setItemEnabled(idverror,false);
    output->setItemEnabled(idstat,false);
    
    if(model1Name->text() == ""  || model2Name->text() == "" )
    {
    	view->setItemEnabled(idsyn,false);
    	view->setItemEnabled(idoverlay,false);
    }
    else
    {
    	view->setItemEnabled(idsyn,true);
    	view->setItemEnabled(idoverlay,true);
    }
    UpdateColorButton->setEnabled(false);
    UpdateHistButton->setEnabled(false);
    middleSlider->setEnabled(false);
    IntervalSpinBox->setEnabled(false);
  }
}

void MeshValmetControls::Model2Open()
{
  QString model2file;
  if(commonPath == "")
    model2file = QFileDialog::getOpenFileName(
                    ".",
                    "3D Models (*.raw *.byu *.iv)",
                    this,
                    "open file dialog",
                    "Choose a file to open" );
  else
    model2file = QFileDialog::getOpenFileName(
		    commonPath,
                    "3D Models (*.raw *.byu *.iv)",
                    this,
                    "open file dialog",
                    "Choose a file to open" );
                    
  if(model2file != "" && model2file != NULL)
  {
    QFileInfo fi(model2file.latin1());
    commonPath = fi.dirPath();
    model2Name->setText(fi.fileName());
    
    if(model2file.right(4).lower() == ".byu")
    {
      vtkBYUReader* byureader = vtkBYUReader::New();
      byureader->SetFileName(model2file.latin1());
      mesh2 = byureader->GetOutput();
    }
    else if(model2file.right(4).lower() == ".raw" || model2file.right(3).lower() == ".iv")
    {
      struct model* model2 = read_model_file(model2file.latin1());
      mesh2 = BuildMeshFromModel(model2);
    }
    SetupWindow(mesh2,ren2,modelActor2);
    modelActor2->GetProperty()->SetColor(1,0,0);
    renWin2->Render();
    
    fileName2 = model2file;
    
    ColormapCheckbox1->setEnabled(false);
    ColormapCheckbox2->setEnabled(false);
    ColormapCheckbox2->setChecked(false);
    //memset(&model2,0,sizeof(*model2));
    
    computed = false;
    //RecomputeHistButton->setEnabled(false);
    comp->setItemEnabled(idcomphist,false);
    //outputHistButton->setEnabled(false);
    output->setItemEnabled(idhist,false);
    output->setItemEnabled(iderror,false);
    output->setItemEnabled(idverror,false);
    output->setItemEnabled(idstat,false);

    if(model1Name->text() == ""  || model2Name->text() == "" )
    {
    	view->setItemEnabled(idsyn,false);
    	view->setItemEnabled(idoverlay,false);
    }
    else
    {
    	view->setItemEnabled(idsyn,true);
    	view->setItemEnabled(idoverlay,true);
    }
    UpdateColorButton->setEnabled(false);
    UpdateHistButton->setEnabled(false);
    middleSlider->setEnabled(false);
    IntervalSpinBox->setEnabled(false);
  }
}



void MeshValmetControls::Model1Switch()
{
  model1Line++;
  model1Line = model1Line % 3;
  if(model1Line == 0)
    modelActor1->GetProperty()->SetRepresentationToSurface();
  else if(model1Line == 1)
    modelActor1->GetProperty()->SetRepresentationToWireframe();
  else if(model1Line == 2)
  {
    modelActor1->GetProperty()->SetRepresentationToPoints();
    modelActor1->GetProperty()->SetPointSize(3.0);
  }
    
  renWin1->Render();
}

void MeshValmetControls::Model2Switch()
{
  model2Line++;
  model2Line = model2Line % 3;
  if(model2Line == 0)
    modelActor2->GetProperty()->SetRepresentationToSurface();
  else if(model2Line == 1)
    modelActor2->GetProperty()->SetRepresentationToWireframe();
  else if(model2Line == 2)
  {
    modelActor2->GetProperty()->SetRepresentationToPoints();
    modelActor2->GetProperty()->SetPointSize(3.0);
  }
    
  renWin2->Render();
}

void MeshValmetControls::SwitchSync()
{
  if(view->isItemEnabled(idsyn))
  {
    if(view->isItemChecked(idsyn))
    {
      view->setItemChecked(idsyn,false);
      renWin1->SetInteractorSync(NULL);
      renWin2->SetInteractorSync(NULL);
    }
    else
    {
      view->setItemChecked(idsyn,true);
      
      vtkCamera *camera1 = vtkCamera::New();
      camera1->SetPosition(1,1,1);
      camera1->SetFocalPoint(0,0,0);
      vtkCamera *camera2 = vtkCamera::New();
      camera2->SetPosition(1,1,1);
      camera2->SetFocalPoint(0,0,0);
      
      ren1->SetActiveCamera(camera1);
      ren2->SetActiveCamera(camera2);
      ren1->ResetCamera();
      ren2->ResetCamera();
      
      renWin1->SetInteractorSync(iren2);
      renWin2->SetInteractorSync(iren1);
      
      renWin1->Render();
      renWin2->Render();
      
    }
  }
}

void MeshValmetControls::ResetCameras()
{
  vtkCamera *camera1 = vtkCamera::New();
  camera1->SetPosition(1,1,1);
  camera1->SetFocalPoint(0,0,0);
  vtkCamera *camera2 = vtkCamera::New();
  camera2->SetPosition(1,1,1);
  camera2->SetFocalPoint(0,0,0);
  
  ren1->SetActiveCamera(camera1);
  ren2->SetActiveCamera(camera2);
  ren1->ResetCamera();
  ren2->ResetCamera();
	
	renWin1->Render();
  renWin2->Render();
}

void MeshValmetControls::SwitchBackground()
{
  // Currently I only allow two colors: BLK or WHT
  if(bgcolor[0] == 0)
  {
    bgcolor[0] = 1;
    bgcolor[1] = 1;
    bgcolor[2] = 1;
  }
  else
  {
    bgcolor[0] = 0;
    bgcolor[1] = 0;
    bgcolor[2] = 0;
  }
  ren1->SetBackground(bgcolor[0],bgcolor[1],bgcolor[2]);
  ren2->SetBackground(bgcolor[0],bgcolor[1],bgcolor[2]);
  renWin1->Render();
  renWin2->Render();
}

void MeshValmetControls::SwitchOverlay()
{
  if(view->isItemEnabled(idoverlay))
  {
    if(view->isItemChecked(idoverlay))
    {
      view->setItemChecked(idoverlay,false);
      ren1->RemoveActor(modelActor22);
      ren2->RemoveActor(modelActor11);
      
      renWin1->Render();
    	renWin2->Render();
    }
    else
    {
      view->setItemChecked(idoverlay,true);
      SetupWindow(mesh1,ren2,modelActor11);
    	modelActor11->GetProperty()->SetColor(1,1,0);
    	modelActor11->GetProperty()->SetRepresentationToPoints();
    	
    	SetupWindow(mesh2,ren1,modelActor22);
    	modelActor22->GetProperty()->SetColor(1,1,0);
    	modelActor22->GetProperty()->SetRepresentationToPoints();
    	
    	renWin1->Render();
    	renWin2->Render();
    }
  }
	
}

void MeshValmetControls::DoStat()
{
  double *serror, *allserror;
  int i,n;
  
  abs_sum_error = 0;
  abs_sum_sqr_error = 0;
  sum_error = 0;
  error_max = -100000;
  sum_square_deviations = 0; 
  abs_sum_square_deviations = 0; 

  double percent95= 0.95;
  double percent68 = 0.68;
  double percent75 = 0.75;
  double percent50 = 0.5; //median
  double percent25 = 0.25;
  percentile95 = 0;
  percentile75 = 0;
  percentile68 = 0;
  percentile25 = 0;
  median = 0; // percentile50
  
  if(pargs.do_symmetric == 1)
  {
    n = model1.n_samples;
    serror = model1.fe[0].serror;
    // We calculate the mean and the STD (sigma) of sample error
    for (i=0; i<n; i++) 
    { 
       sum_error += serror[i];
       abs_sum_error += fabs(serror[i]);
       abs_sum_sqr_error += serror[i]*serror[i];
       if(error_max < serror[i])
    	   error_max = serror[i];
    }
    if(pargs.signeddist)
      mean_error = sum_error / n;
    else
      mean_error = abs_sum_error / n;
    mad_error = abs_sum_error / n;
    msd_error = abs_sum_sqr_error /n;
      
    for (i=0; i<n; i++)
    {
    	sum_square_deviations += (serror[i]-mean_error)*(serror[i]-mean_error);
    	abs_sum_square_deviations += (fabs(serror[i])-mad_error)*(fabs(serror[i])-mad_error);
    }

    if(pargs.signeddist)
      sigma_error = sqrt(sum_square_deviations/(n-1));
    else
      sigma_error = sqrt(abs_sum_square_deviations/(n-1));
    abs_sigma_error = sqrt(abs_sum_square_deviations/(n-1));

    // Calculate median, and percentiles

    // Previous computation was wrong
    // Official def: The 95th percentile is the smallest number that is 
    //               greater that 95% of the numbers in a given set.   

    allserror = new double[n];
    if(pargs.signeddist)
      memcpy(allserror, serror, n*sizeof(double));
    else {
      for (i=0; i<n; i++) 
	allserror[i] = fabs(serror[i]);
    }

    qsort((void *)allserror,n,sizeof(double),compare_doubles);

    percentile25 = allserror[(int)(n*percent25+1)];
    median = allserror[(int)(n*percent50+1)];
    percentile68 = allserror[(int)(n*percent68+1)];
    percentile75 = allserror[(int)(n*percent75+1)];
    percentile95 = allserror[(int)(n*percent95+1)];

    delete [] allserror;
      
  }
  else if(pargs.do_symmetric == 2)
  {
    n = model2.n_samples;
    serror = model2.fe[0].serror;
    // We calculate the mean and the STD (sigma) error
   for (i=0; i<n; i++) 
   { 
       sum_error += serror[i];
       abs_sum_error += fabs(serror[i]);
       abs_sum_sqr_error += serror[i]*serror[i];
       if(error_max < serror[i])
    	   error_max = serror[i];
    }
    if(pargs.signeddist)
      mean_error = sum_error / n;
    else
      mean_error = abs_sum_error / n;
    mad_error = abs_sum_error / n;
    msd_error = abs_sum_sqr_error /n;
      
    for (i=0; i<n; i++)
    {
    	sum_square_deviations += (serror[i]-mean_error)*(serror[i]-mean_error);
    	abs_sum_square_deviations += (fabs(serror[i])-mad_error)*(fabs(serror[i])-mad_error);
    }

    if(pargs.signeddist)
      sigma_error = sqrt(sum_square_deviations/(n-1));
    else
      sigma_error = sqrt(abs_sum_square_deviations/(n-1));
    abs_sigma_error = sqrt(abs_sum_square_deviations/(n-1));

    allserror = new double[n];
    if(pargs.signeddist)
      memcpy(allserror, serror, n*sizeof(double));
    else {
      for (i=0; i<n; i++) 
	allserror[i] = fabs(serror[i]);
    }

    qsort((void *)allserror,n,sizeof(double),compare_doubles);

    percentile25 = allserror[(int)(n*percent25+1)];
    median = allserror[(int)(n*percent50+1)];
    percentile68 = allserror[(int)(n*percent68+1)];
    percentile75 = allserror[(int)(n*percent75+1)];
    percentile95 = allserror[(int)(n*percent95+1)];

    delete [] allserror;
  }
  else if(pargs.do_symmetric == 3)
  {
    int n1, n2;
    n1 = model1.n_samples;
    serror = model1.fe[0].serror;
    // We calculate the mean and the SD (sigma) error
    for (i=0; i<n1; i++) 
    { 
      sum_error += serror[i];
      abs_sum_error += fabs(serror[i]);
      abs_sum_sqr_error += serror[i]*serror[i];
      if(error_max < serror[i])
     	 error_max = serror[i];
    }
    n2 = model2.n_samples;
    serror = model2.fe[0].serror;
    
    for (i=0; i<n2; i++) 
    {
      sum_error += serror[i];
      abs_sum_error += fabs(serror[i]);
      abs_sum_sqr_error += serror[i]*serror[i];
      if(error_max < serror[i])
    	  error_max = serror[i];
    }

    n = n1 + n2;
    if(pargs.signeddist)
      mean_error = sum_error / n;
    else
      mean_error = abs_sum_error / n;
    mad_error = abs_sum_error / n;
    msd_error = abs_sum_sqr_error /n;
        
    serror = model1.fe[0].serror;
    for (i=0; i<n1; i++) 
    {
        sum_square_deviations += (serror[i]-mean_error)*(serror[i]-mean_error);
        abs_sum_square_deviations += (fabs(serror[i])-mad_error)*(fabs(serror[i])-mad_error);
    }
    serror = model2.fe[0].serror;
    for (i=0; i<n2; i++) 
    {
      sum_square_deviations += (serror[i]-mean_error)*(serror[i]-mean_error);
      abs_sum_square_deviations += (fabs(serror[i])-mad_error)*(fabs(serror[i])-mad_error);
    }
    if(pargs.signeddist)
      sigma_error = sqrt(sum_square_deviations/(n-1));
    else
      sigma_error = sqrt(abs_sum_square_deviations/(n-1));
    abs_sigma_error = sqrt(abs_sum_square_deviations/(n-1));

    allserror = new double[n];
    if(pargs.signeddist) {
      memcpy(allserror, model1.fe[0].serror, n1*sizeof(double));
      memcpy(&allserror[n1], model2.fe[0].serror, n2*sizeof(double));
    }
    else {
      serror = model1.fe[0].serror;
      for (i=0; i<n1; i++) 
	allserror[i] = fabs(serror[i]);

      serror = model2.fe[0].serror;
      for (i=n1; i<n1+n2; i++) 
	allserror[i] = fabs(serror[i-n1]);
    }

    qsort((void *)allserror,n,sizeof(double),compare_doubles);

    percentile25 = allserror[(int)(n*percent25+1)];
    median = allserror[(int)(n*percent50+1)];
    percentile68 = allserror[(int)(n*percent68+1)];
    percentile75 = allserror[(int)(n*percent75+1)];
    percentile95 = allserror[(int)(n*percent95+1)];

    delete [] allserror;
  }
  
  //Calculate face mean and face rms
  if(pargs.signeddist)
  {
    if(pargs.do_symmetric == 1)
    {
      face_mean_error = stats.mean_dist;
      dist_volume = stats.abs_mean_tot;
      face_rms_error = stats.rms_dist;
    }
    else if(pargs.do_symmetric == 2)
    {
      face_mean_error = stats_rev.mean_dist;
      dist_volume = stats_rev.abs_mean_tot;
      face_rms_error = stats_rev.rms_dist;
    }
    else if(pargs.do_symmetric == 3)
    {
      face_mean_error = (stats.mean_tot+stats_rev.mean_tot)/(stats.st_m1_area+stats_rev.st_m1_area);
      dist_volume = (stats.abs_mean_tot+stats_rev.abs_mean_tot)/2.0;
      face_rms_error = sqrt((stats.rms_tot+stats_rev.rms_tot)/(stats.st_m1_area+stats_rev.st_m1_area));
    }
  }
  else
  {
    if(pargs.do_symmetric == 1)
    {
      face_mean_error = stats.abs_mean_dist;
      dist_volume = stats.abs_mean_tot;
      face_rms_error = stats.abs_rms_dist;
    }
    else if(pargs.do_symmetric == 2)
    {
      face_mean_error = stats_rev.abs_mean_dist;
      dist_volume = stats_rev.abs_mean_tot;
      face_rms_error = stats_rev.abs_rms_dist;
    }
    else if(pargs.do_symmetric == 3)
    {
      face_mean_error = (stats.abs_mean_tot+stats_rev.abs_mean_tot)/(stats.st_m1_area+stats_rev.st_m1_area);
      dist_volume = (stats.abs_mean_tot+stats_rev.abs_mean_tot)/2.0;
      face_rms_error = sqrt((stats.abs_rms_tot+stats_rev.abs_rms_tot)/(stats.st_m1_area+stats_rev.st_m1_area));
    }	
  }
  
}


void MeshValmetControls::OutputStatInfo()
{
  if(pargs.signeddist)
  {
    StdMeanEdit->setText(QString::number(mean_error));
    SigmaEdit->setText(QString::number(sigma_error));
  }
  else
  {
    StdMeanEdit->setText(QString::number(mad_error));
    SigmaEdit->setText(QString::number(abs_sigma_error));
  }
  
  
  MeanEdit->setText(QString::number(face_mean_error));
  RMSEdit->setText(QString::number(face_rms_error));  
  
  MedianEdit->setText(QString::number(median));  
  Percentile95Edit->setText(QString::number(percentile95));
  Percentile75Edit->setText(QString::number(percentile75));
  Percentile68Edit->setText(QString::number(percentile68));
  Percentile25Edit->setText(QString::number(percentile25));

  DiceEdit->setText(QString::number(dice_coefficient[0]));
  IntUnionEdit->setText(QString::number(int_union_ratio[0]));
    
  /* Christine: Hausdorff distance is defined as maximum of absolute distance here.*/
  hausdorff_error = fabs(dmax)>fabs(dmin)?fabs(dmax):fabs(dmin);  
  HausdorffEdit->setText(QString::number(hausdorff_error));
  MSDEdit->setText(QString::number(msd_error));
  MADEdit->setText(QString::number(mad_error));
  //DistVolEdit->setText(QString::number(dist_volume));
  if(pargs.do_symmetric == 1)
  {
  	if(pargs.signeddist)
  	{
    	MinEdit->setText(QString::number(model1.min_error));
    	MaxEdit->setText(QString::number(model1.max_error));    
    }
    else
    {
    	MinEdit->setText(QString::number(model1.abs_min_error));
    	MaxEdit->setText(QString::number(model1.abs_max_error)); 
    }
  }
  else if(pargs.do_symmetric == 2)
  {
    if(pargs.signeddist)
  	{
    	MinEdit->setText(QString::number(model2.min_error));
    	MaxEdit->setText(QString::number(model2.max_error));    
    }
    else
    {
    	MinEdit->setText(QString::number(model2.abs_min_error));
    	MaxEdit->setText(QString::number(model2.abs_max_error)); 
    } 
  }
  else if(pargs.do_symmetric == 3)
  {
    MinEdit->setText(QString::number(dmin));
    MaxEdit->setText(QString::number(dmax));    
  }
}

void MeshValmetControls::Compute()
{
  QString file1 = model1Name->text();
  QString file2 = model2Name->text();
  if(file1 == "" || file2 == "")
  {
    QMessageBox::warning(this,"Warning","Please input two model files!",QMessageBox::Ok,QMessageBox::NoButton,QMessageBox::NoButton);
    return;
  }
  
    QProgressDialog *qProg;    
    struct prog_reporter pr;
    
    qProg = NULL;
    memset(&model1,0,sizeof(model1));
    memset(&model2,0,sizeof(model2));
    memset(&pr,0,sizeof(pr));
    log = NULL;
    
    qProg = new QProgressDialog("Calculating distance",0,100);
    qProg->setMinimumDuration(1500);
    pr.prog = QT_prog;
    pr.cb_out = qProg;
    
    MeshSetup(&pargs);
    mesh_run(&pargs, &model1, &model2, log, &pr, &stats, &stats_rev, &abs_sampling_step, &abs_sampling_dens);

    int num_vert1 = model1.mesh->num_vert;
    int num_vert2 = model2.mesh->num_vert;
    int num_faces1 = model1.mesh->num_faces;
    int num_faces2 = model2.mesh->num_faces;

    double * L1 = new double[3*num_vert1];
    double * L2 = new double[3*num_vert2];
    int * T1 = new int[3*num_faces1];
    int * T2 = new int[3*num_faces2];

    vertex_t * vert_list = model1.mesh->vertices;
    for(int i=0; i<num_vert1; i++) {
      L1[3*i+0] = vert_list[i].x;
      L1[3*i+1] = vert_list[i].y;
      L1[3*i+2] = vert_list[i].z;
    }
    vert_list = model2.mesh->vertices;
    for(int i=0; i<num_vert2; i++) {
      L2[3*i+0] = vert_list[i].x;
      L2[3*i+1] = vert_list[i].y;
      L2[3*i+2] = vert_list[i].z;
    }

    face_t * face_list = model1.mesh->faces;
    for(int i=0; i<num_faces1; i++) {
      T1[3*i+0] = face_list[i].f0+1;
      T1[3*i+1] = face_list[i].f1+1;
      T1[3*i+2] = face_list[i].f2+1;
    }
    face_list = model2.mesh->faces;
    for(int i=0; i<num_faces2; i++) {
      T2[3*i+0] = face_list[i].f0+1;
      T2[3*i+1] = face_list[i].f1+1;
      T2[3*i+2] = face_list[i].f2+1;
    } 

    ComputeRobustVolumeOverlap(L1,L2,num_vert1,num_vert2,T1,T2,num_faces1,num_faces2,dice_coefficient,int_union_ratio);

    delete [] L1;
    delete [] L2;
    delete [] T1;
    delete [] T2;

    int t = pargs.do_symmetric;
    bool signedd = pargs.signeddist;
    
    // dmin dmax are updated according to signedd
    if(t == 1)
    {
      if(signedd)
      {
	dmax = model1.max_error; dmin = model1.min_error;
      }
      else{
	dmax = model1.abs_max_error; dmin = model1.abs_min_error;
      }
    }
    else if(t == 2)
    {
      if(signedd)
      {
	dmax = model2.max_error; dmin = model2.min_error;
      }
      else
      {
	dmax = model2.abs_max_error; dmin = model2.abs_min_error;
      }
    }
    else if(t == 3)
    {
      if(signedd)
      {
	dmax = (model1.max_error > model2.max_error)?model1.max_error:model2.max_error;
	dmin = (model1.min_error > model2.min_error)?model2.min_error:model1.min_error;
      }
      else
      {
	dmax = (model1.abs_max_error > model2.abs_max_error)?model1.abs_max_error:model2.abs_max_error;
	dmin = (model1.abs_min_error > model2.abs_min_error)?model2.abs_min_error:model1.abs_min_error;
      }
    }     		
    
    downsampling = DownsamplingSlider->value();
    if(pargs.do_symmetric == 1 || pargs.do_symmetric == 3)
    {
      drawVertexErrorT(mesh1, &model1, modelActor1, renWin1);
      ColormapCheckbox1->setEnabled(true);
      ColormapCheckbox1->setChecked(TRUE);
    }
    if(pargs.do_symmetric == 2 || pargs.do_symmetric == 3)
    {
      drawVertexErrorT(mesh2, &model2, modelActor2, renWin2);
      ColormapCheckbox2->setEnabled(true);
      ColormapCheckbox2->setChecked(TRUE);
    }    
    
  ColormapMin->setText(QString::number(dmin));
  ColormapMax->setText(QString::number(dmax));
  
  //compute sample error mean and std, face error mean and root mean square error
  DoStat();

  OutputStatInfo();
 
  computed = true;
  //RecomputeHistButton->setEnabled(true);
  comp->setItemEnabled(idcomphist,true);
  //outputHistButton->setEnabled(true);
  output->setItemEnabled(idhist,true);
  output->setItemEnabled(iderror,true);
  output->setItemEnabled(idverror,true);
  output->setItemEnabled(idstat,true);
  view->setItemEnabled(idsyn,true);
  UpdateColorButton->setEnabled(true);
  UpdateHistButton->setEnabled(true);
  middleSlider->setEnabled(true);
  IntervalSpinBox->setEnabled(true);
  
  ComputeHistogram();

  // output results to the standard output
  fprintf(stdout,"\n\nStatistical results:\n");
  fprintf(stdout,"\t Min:\t\t\t%f\n", dmin);
  fprintf(stdout,"\t Max:\t\t\t%f\n", dmax);
  fprintf(stdout,"\t VertexMean:\t\t%f\n", mean_error);
  fprintf(stdout,"\t VertexStd:\t\t%f\n", sigma_error);
  fprintf(stdout,"\t FaceMean:\t\t%f\n", face_mean_error);
  fprintf(stdout,"\t FaceRMS:\t\t%f\n", face_rms_error);
  fprintf(stdout,"\t Hausdorff:\t\t%f\n", hausdorff_error);
  fprintf(stdout,"\t MSD:\t\t\t%f\n", msd_error);
  fprintf(stdout,"\t MAD:\t\t\t%f\n", mad_error);
  fprintf(stdout,"\t Median:\t\t%f\n", median);
  //fprintf(stdout,"\t DistVolume:\t\t%f\n", dist_volume);
  fprintf(stdout,"\t 25percentile:\t\t%f\n", percentile25);
  fprintf(stdout,"\t 68percentile:\t\t%f\n", percentile68);
  fprintf(stdout,"\t 75percentile:\t\t%f\n", percentile75);
  fprintf(stdout,"\t 95percentile:\t\t%f\n", percentile95);

  // fprintf(stdout,"\t Min\tMax\tVertexMean\tVertexSigma\tFaceMean\tFaceRMS\tHausdorff\tMSD\tMAD\tMedian\tDistVolume\t25percentile\t68percentile\t75percentile\t95percentile\n");
  //fprintf(stdout,"\t %f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",dmin, dmax, mean_error, sigma_error, face_mean_error, face_rms_error, hausdorff_error, msd_error, mad_error, median, dist_volume, percentile25, percentile68, percentile75,percentile95);

}

void MeshValmetControls::drawVertexErrorT(vtkPolyData* mesh, model_error* model, vtkActor* actor, vtkQtRenderWindow* win, int tag)
{
  //mesh->Delete();
  mesh = vtkPolyData::New();
  vtkPoints *points = vtkPoints::New();
  vtkCellArray *polys = vtkCellArray::New();
  vtkDoubleArray *scalars = vtkDoubleArray::New();
  
  int k,i,j,jmax,n;
  vertex_t u,v;
  vertex_t a,b,c;
  face_t *cur_face;
  int i0,i1,i2,i3;
  int j0,j1,j2,j3;
  int l0,l1,l2,l3;
  vertex_t v0,v1,v2,v3;
  
  vtkIdType index = 0;
  double vertex[3];
  vtkIdType f[3];
  
  for (k=0; k<model->mesh->num_faces; k++) 
  {
    n = model->fe[k].sample_freq;
    cur_face = &(model->mesh->faces[k]);
    if (n == 1 && downsampling == 1) 
    {
      // displaying only at triangle vertices + center
      a = model->mesh->vertices[cur_face->f0];
      b = model->mesh->vertices[cur_face->f1];
      c = model->mesh->vertices[cur_face->f2];
      v3.x = 1/3.0*(a.x+b.x+c.x);
      v3.y = 1/3.0*(a.y+b.y+c.y);
      v3.z = 1/3.0*(a.z+b.z+c.z);
  
      vertex[0] = a.x;
      vertex[1] = a.y;
      vertex[2] = a.z;
      points->InsertPoint(index,vertex);
      scalars->InsertTuple1(index++, model->verror[cur_face->f0]);        
      vertex[0] = b.x;
      vertex[1] = b.y;
      vertex[2] = b.z;
      points->InsertPoint(index,vertex);
      scalars->InsertTuple1(index++, model->verror[cur_face->f1]);
      vertex[0] = v3.x;
      vertex[1] = v3.y;
      vertex[2] = v3.z;
      points->InsertPoint(index,vertex);
      scalars->InsertTuple1(index++, model->fe[k].serror[0]);//?
      f[0] = index - 3;
      f[1] = index - 2;
      f[2] = index - 1;
      polys->InsertNextCell(3,f);
      
      vertex[0] = a.x;
      vertex[1] = a.y;
      vertex[2] = a.z;
      points->InsertPoint(index,vertex);
      scalars->InsertTuple1(index++, model->verror[cur_face->f0]);        
      vertex[0] = v3.x;
      vertex[1] = v3.y;
      vertex[2] = v3.z;
      points->InsertPoint(index,vertex);
      scalars->InsertTuple1(index++, model->fe[k].serror[0]);//?
      vertex[0] = c.x;
      vertex[1] = c.y;
      vertex[2] = c.z;
      points->InsertPoint(index,vertex);
      scalars->InsertTuple1(index++, model->verror[cur_face->f2]);        
      f[0] = index - 3;
      f[1] = index - 2;
      f[2] = index - 1;
      polys->InsertNextCell(3,f);
        
      vertex[0] = b.x;
      vertex[1] = b.y;
      vertex[2] = b.z;
      points->InsertPoint(index,vertex);
      scalars->InsertTuple1(index++, model->verror[cur_face->f1]);
      vertex[0] = c.x;
      vertex[1] = c.y;
      vertex[2] = c.z;
      points->InsertPoint(index,vertex);
      scalars->InsertTuple1(index++, model->verror[cur_face->f2]);        
      vertex[0] = v3.x;
      vertex[1] = v3.y;
      vertex[2] = v3.z;
      points->InsertPoint(index,vertex);
      scalars->InsertTuple1(index++, model->fe[k].serror[0]);//?                
      f[0] = index - 3;
      f[1] = index - 2;
      f[2] = index - 1;
      polys->InsertNextCell(3,f);

    } 
    else if (downsampling >= n) 
    {
      //displaying only at triangle vertices 
      vertex[0] = model->mesh->vertices[cur_face->f0].x;
      vertex[1] = model->mesh->vertices[cur_face->f0].y;
      vertex[2] = model->mesh->vertices[cur_face->f0].z;
      points->InsertPoint(index,vertex);
      scalars->InsertTuple1(index++, model->verror[cur_face->f0]);
    
      vertex[0] = model->mesh->vertices[cur_face->f1].x;
      vertex[1] = model->mesh->vertices[cur_face->f1].y;
      vertex[2] = model->mesh->vertices[cur_face->f1].z;
      points->InsertPoint(index,vertex);
      scalars->InsertTuple1(index++, model->verror[cur_face->f1]);
      
      vertex[0] = model->mesh->vertices[cur_face->f2].x;
      vertex[1] = model->mesh->vertices[cur_face->f2].y;
      vertex[2] = model->mesh->vertices[cur_face->f2].z;
      points->InsertPoint(index,vertex);
      scalars->InsertTuple1(index++, model->verror[cur_face->f2]);
      
      f[0] = index - 3;
      f[1] = index - 2;
      f[2] = index - 1;
      polys->InsertNextCell(3,f);
    }
    else
    {
      /* displaying at error samples and triangle vertices */
      //assert(n > 1);
      a = model->mesh->vertices[cur_face->f0];
      b = model->mesh->vertices[cur_face->f1];
      c = model->mesh->vertices[cur_face->f2];
      substract_v(&b,&a,&u);
      substract_v(&c,&a,&v);
      prod_v(1/(float)(n-1),&u,&u);
      prod_v(1/(float)(n-1),&v,&v);
      for (i=0; i<n-1; i+=downsampling) 
      {
	i2 = (i+downsampling < n) ? i+downsampling : n-1;
	for (j=0, jmax=n-i-1; j<jmax; j+=downsampling) 
	{
	  if (i+j+downsampling < n) 
	  {
	    i0 = i;
	    j0 = j;
	    i1 = i+downsampling;
	    j1 = j;
	    i2 = i;
	    j2 = j+downsampling;
	    i3 = i1;
	    j3 = j2;
	  } 
	  else 
	  {
	    i2 = i;
	    j2 = j;
	    i0 = (i+downsampling < n) ? i+downsampling : n-1;
	    j0 = (j>0) ? j-downsampling : j;
	    //assert(j0 >= 0);
	    i1 = i0;
	    j1 = n-1-i1;
	    i3 = i;
	    j3 = n-1-i3;
	    //assert(j3 >= 0);
	  }
	  l0 = j0+i0*(2*n-i0+1)/2;
	  l1 = j1+i1*(2*n-i1+1)/2;
	  l2 = j2+i2*(2*n-i2+1)/2;
	  v0.x = a.x+i0*u.x+j0*v.x;
	  v0.y = a.y+i0*u.y+j0*v.y;
	  v0.z = a.z+i0*u.z+j0*v.z;
	  v1.x = a.x+i1*u.x+j1*v.x;
	  v1.y = a.y+i1*u.y+j1*v.y;
	  v1.z = a.z+i1*u.z+j1*v.z;
	  v2.x = a.x+i2*u.x+j2*v.x;
	  v2.y = a.y+i2*u.y+j2*v.y;
	  v2.z = a.z+i2*u.z+j2*v.z;
	  if (i0 != i1 || j0 != j1) /* avoid possible degenerate */
	  {                  
	    vertex[0] = v0.x;
	    vertex[1] = v0.y;
	    vertex[2] = v0.z;
	    points->InsertPoint(index,vertex);
	    scalars->InsertTuple1(index++, model->fe[k].serror[l0]);//?
                    
	    vertex[0] = v1.x;
	    vertex[1] = v1.y;
	    vertex[2] = v1.z;
	    points->InsertPoint(index,vertex);
	    scalars->InsertTuple1(index++, model->fe[k].serror[l1]);//?
          
	    vertex[0] = v2.x;
	    vertex[1] = v2.y;
	    vertex[2] = v2.z;
	    points->InsertPoint(index,vertex);
	    scalars->InsertTuple1(index++, model->fe[k].serror[l2]);//?

	    f[0] = index - 3;
	    f[1] = index - 2;
	    f[2] = index - 1;
	    polys->InsertNextCell(3,f);
	  }
	  if (i3+j3 < n) 
	  {
	    l3 = j3+i3*(2*n-i3+1)/2;
	    v3.x = a.x+i3*u.x+j3*v.x;
	    v3.y = a.y+i3*u.y+j3*v.y;
	    v3.z = a.z+i3*u.z+j3*v.z;
          
	    vertex[0] = v3.x;
	    vertex[1] = v3.y;
	    vertex[2] = v3.z;
	    points->InsertPoint(index,vertex);
	    scalars->InsertTuple1(index++, model->fe[k].serror[l3]);//?
                    
	    vertex[0] = v2.x;
	    vertex[1] = v2.y;
	    vertex[2] = v2.z;
	    points->InsertPoint(index,vertex);
	    scalars->InsertTuple1(index++, model->fe[k].serror[l2]);
          
	    vertex[0] = v1.x;
	    vertex[1] = v1.y;
	    vertex[2] = v1.z;
	    points->InsertPoint(index,vertex);
	    scalars->InsertTuple1(index++, model->fe[k].serror[l1]);
                    
	    f[0] = index - 3;
	    f[1] = index - 2;
	    f[2] = index - 1;
	    polys->InsertNextCell(3,f);
	  }
	}
      }
    }
  }

  mesh->SetPoints(points);
  points->Delete();
  mesh->SetPolys(polys);
  polys->Delete();
  if(tag == 0)  mesh->GetPointData()->SetScalars(scalars);  
  scalars->Delete();
  
  vtkPolyDataMapper *meshMapper = vtkPolyDataMapper::New(); 
  meshMapper->SetInput(mesh); 
  
  double mmax,mmin;
  int inter = IntervalSpinBox->value();

  mmax = dmax;
  mmin = dmin;
  
  if(tag == 0)
  {
    lut = vtkColorTransferFunction::New();
    
    //Begin seting up my own lookup table
    //The Middle point always points to zero distance
    lut->AddRGBPoint(middle,0,1,0);
    
    //lut->AddRGBPoint(0.02,lookuptable[58].R/(float)255,lookuptable[58].G/(float)255,lookuptable[58].B/(float)255);
    //lut->AddRGBPoint(-0.02,lookuptable[198].R/(float)255,lookuptable[198].G/(float)255,lookuptable[198].B/(float)255);
    
    if(mmax-middle>0 && mmin-middle>=0)
    {
      int i;
      for(i = 1; i <=inter; i++)
        lut->AddRGBPoint((float)(middle+(mmax-middle)*i/(double)inter),lookuptable[128-i*22].R/(float)255,lookuptable[128-i*22].G/(float)255,lookuptable[128-i*22].B/(float)255);

    }
    else if(mmax-middle<=0 && mmin-middle<0)
    {
      int i;
      for(i = 1; i <=inter; i++)
        lut->AddRGBPoint((float)(middle+(mmin-middle)*i/(double)inter),lookuptable[128+i*22].R/(float)255,lookuptable[128+i*22].G/(float)255,lookuptable[128+i*22].B/(float)255);
    }
    else if(mmax-middle>0 && mmin-middle<0)
    {
      //The Positive outside distance
      int i;
      for(i = 1; i <=inter; i++)
      {
        lut->AddRGBPoint((float)(middle+(mmax-middle)*i/(double)inter),lookuptable[128-i*22].R/(float)255,lookuptable[128-i*22].G/(float)255,lookuptable[128-i*22].B/(float)255);
        //printf("%f\t%f\t%f\t%f\n",(float)(mmax*i/(double)5),lookuptable[128-i].R/(float)255,lookuptable[128-i].G/(float)255,lookuptable[128-i].B/(float)255);
      }
        
      //The Negative inside distance
      for(i = 1; i <=inter; i++)
        lut->AddRGBPoint((float)(middle+(mmin-middle)*i/(double)inter),lookuptable[128+i*22].R/(float)255,lookuptable[128+i*22].G/(float)255,lookuptable[128+i*22].B/(float)255);
  
    }
        
    meshMapper->SetLookupTable(lut);
  }
  actor->SetMapper(meshMapper);
  
  win->updateGL();

}

void MeshValmetControls::MeshSetup(struct args* pargs) {
  pargs->m1_fname =  (char *)fileName1.latin1();
  pargs->m2_fname =  (char *)fileName2.latin1();
  /* The sampling step, as fraction of the bounding box diagonal of model 2. */
  pargs->sampling_step = atof((char*)samplingStep->text().latin1())/100;
  /* Minimum sampling frequency to enforce on each triangle */
  pargs->min_sample_freq = atoi((char*)minSamplingFq->text().latin1());
  pargs->do_wlog = output->isItemChecked(idtext);
  pargs->verb_analysis = false;
  pargs->do_texture = false;
  pargs->signeddist = true;
  
  if(a2bRadio->isChecked())
    pargs->do_symmetric = 1;
  else if(b2aRadio->isChecked())
    pargs->do_symmetric = 2;
  else if(symmetricRadio->isChecked())
    pargs->do_symmetric = 3;
    
  if(SignedRadio->isChecked())
  	pargs->signeddist = true;
  else if(AbsRadio->isChecked())
  	pargs->signeddist = false;

  if (!pargs->do_wlog) {
    log = outbuf_new(stdio_puts,stdout);
  } else {
    textOut = new TextWidget();
    log = outbuf_new(TextWidget_puts,textOut);
    textOut->show();
  }
}

void MeshValmetControls::SelectStep(const QString& step)
{
  samplingStep->setText(step);  
}

void MeshValmetControls::SaveModel1()
{
  if(ren1->GetNumberOfPropsRendered() != 0 )
  {
    QString outFile;
    if(commonPath == "")
      outFile = QFileDialog::getSaveFileName(
                      ".",
                      "Mesh Info (*.iv)",
                      this,
                      "save file dialog",
                      "Choose a filename to save under" );
          else
            outFile = QFileDialog::getSaveFileName(
                      commonPath,
                      "Mesh Info (*.iv)",
                      this,
                      "save file dialog",
                      "Choose a filename to save under" );
  
    if(outFile != 0 )
    {
      if(outFile.right(3) != ".iv")
        outFile = outFile + ".iv";      
      
      if ( QFile::exists( outFile ) )
      {
        QMessageBox mb( "MeshValmet",
                        "Saving the file will overwrite the old file on disk.\n"
                        "Do you really want to save?",
                        QMessageBox::Information,
                        QMessageBox::Yes | QMessageBox::Default,
                        QMessageBox::No,
                        QMessageBox::Cancel | QMessageBox::Escape );
            mb.setButtonText( QMessageBox::Yes, "Save" );
            mb.setButtonText( QMessageBox::No, "Don't Save" );
            switch( mb.exec() ) 
            {
                case QMessageBox::Yes:
		  ren1->SetBackground(0,0,0);
		  SaveMeshToIV(renWin1, outFile);
		  ren1->SetBackground(bgcolor[0],bgcolor[1],bgcolor[2]);
		  break;
                case QMessageBox::No:
		  break;
                case QMessageBox::Cancel:
		  break;
            }
      }
      else
      {
        ren1->SetBackground(0,0,0);
        SaveMeshToIV(renWin1, outFile);
        ren1->SetBackground(bgcolor[0],bgcolor[1],bgcolor[2]);
      }
              
      QFileInfo fi(outFile);
      commonPath = fi.dirPath();

    }
  }
}

void MeshValmetControls::SaveModel2()
{
  if(ren2->GetNumberOfPropsRendered() != 0 )
  {
    QString outFile;
    if(commonPath == "")
      outFile = QFileDialog::getSaveFileName(
                      ".",
                      "Mesh Info (*.iv)",
                      this,
                      "save file dialog",
                      "Choose a filename to save under" );
          else
            outFile = QFileDialog::getSaveFileName(
                      commonPath,
                      "Mesh Info (*.iv)",
                      this,
                      "save file dialog",
                      "Choose a filename to save under" );
  
    if(outFile != 0 )
    {
      if(outFile.right(3) != ".iv")
        outFile = outFile + ".iv";
      
      if ( QFile::exists( outFile ) )
      {
        QMessageBox mb( "MeshValmet",
                        "Saving the file will overwrite the old file on disk.\n"
                        "Do you really want to save?",
                        QMessageBox::Information,
                        QMessageBox::Yes | QMessageBox::Default,
                        QMessageBox::No,
                        QMessageBox::Cancel | QMessageBox::Escape );
            mb.setButtonText( QMessageBox::Yes, "Save" );
            mb.setButtonText( QMessageBox::No, "Don't Save" );
            switch( mb.exec() ) 
            {
                case QMessageBox::Yes:
                  ren2->SetBackground(0,0,0);
		  SaveMeshToIV(renWin2, outFile);
		  ren2->SetBackground(1,1,1);
		  break;
                case QMessageBox::No:
		  break;
	        case QMessageBox::Cancel:
		  break;
            }
      }
      else
      {
        ren2->SetBackground(0,0,0);
        SaveMeshToIV(renWin2, outFile);
        ren2->SetBackground(1,1,1);
      }
      
      QFileInfo fi(outFile);
      commonPath = fi.dirPath();    
    }
  }
}

void MeshValmetControls::SaveMeshToIV(vtkQtRenderWindow* win, QString name)
{
  vtkVRMLExporter* exporter = vtkVRMLExporter::New();
  exporter->SetRenderWindow(win);
  exporter->SetFileName(name.latin1());
  exporter->Write();
  exporter->Delete();
}

void MeshValmetControls::ChangeSamplingSlider()
{
  if(!computed)
    return;
  downsampling = DownsamplingSlider->value();
  if(pargs.do_symmetric == 1 || pargs.do_symmetric == 3)
    if(ColormapCheckbox1->isChecked())
          drawVertexErrorT(mesh1, &model1, modelActor1, renWin1);
        else
          drawVertexErrorT(mesh1, &model1, modelActor1, renWin1, 1);
    if(pargs.do_symmetric == 2 || pargs.do_symmetric == 3)
      if(ColormapCheckbox2->isChecked())
          drawVertexErrorT(mesh2, &model2, modelActor2, renWin2);
        else
          drawVertexErrorT(mesh2, &model2, modelActor2, renWin2, 1);
}

void MeshValmetControls::ColormapCheckbox1Changed(bool checked)
{
  if(!computed)
    return;
  if(checked)
    drawVertexErrorT(mesh1, &model1, modelActor1, renWin1, 0);
  else
    drawVertexErrorT(mesh1, &model1, modelActor1, renWin1, 1);
}

void MeshValmetControls::ColormapCheckbox2Changed(bool checked)
{
  if(!computed)
    return;
  if(checked)
    drawVertexErrorT(mesh2, &model2, modelActor2, renWin2, 0);
  else
    drawVertexErrorT(mesh2, &model2, modelActor2, renWin2, 1);
}


void MeshValmetControls::ComputeHistogram()
{  
  if(computed == false)
    return;
  cmap_len = BinsSpinBox->value();
  DoHistogram();
  
  //draw histogram
  xyplot = vtkXYPlotActor::New();
  xyplot->PlotPointsOn();
  xyplot->GetProperty()->SetPointSize(3);
  xyplot->GetProperty()->SetColor(0,0,1);
  xyplot->GetTitleTextProperty()->SetColor(0,0,1);
  xyplot->SetTitle("Histograms of Mesh Errors");
  xyplot->SetNumberOfXLabels(6);
  xyplot->GetAxisTitleTextProperty()->SetColor(0,0,1);
  xyplot->SetXTitle("Distance");
  xyplot->SetYTitle("# of sample points");
  xyplot->GetPositionCoordinate()->SetValue(0, 0, 0);
  xyplot->GetPosition2Coordinate()->SetValue(1, 1, 0);//relative to Position
  xyplot->SetXValuesToValue();
  xyplot->LegendOn();
  xyplot->GetAxisLabelTextProperty()->SetColor(0,0,1);
  
    
  vtkPolyData *mesh1 = vtkPolyData::New();
  vtkPolyData *mesh2 = vtkPolyData::New();
  vtkPolyData *mesh3 = vtkPolyData::New();
  
  if(pargs.do_symmetric == 1 || pargs.do_symmetric == 3)
  {
    vtkPoints *points = vtkPoints::New();
    vtkDoubleArray *scalars = vtkDoubleArray::New();
    
    int i;
    double v[3];
    for (i=0; i<cmap_len; i++)
    {
      v[0] = dmin + (double)(dmax-dmin)*i / (double)cmap_len;
      v[1] = dmin + (double)(dmax-dmin)*i / (double)cmap_len;
      v[2] = dmin + (double)(dmax-dmin)*i / (double)cmap_len;
      points->InsertPoint(i,v);
      scalars->InsertTuple1(i,(double)histogram1[i]);
    }
    mesh1->SetPoints(points);
    points->Delete();
    mesh1->GetPointData()->SetScalars(scalars);
    scalars->Delete();        
    xyplot->AddInput((vtkDataSet *)mesh1);
    
    xyplot->SetPlotColor(0,0,0,1);
    xyplot->SetPlotLabel(0, "A 2 B");

  }
  if(pargs.do_symmetric == 2 || pargs.do_symmetric == 3)
  {
    vtkPoints *points = vtkPoints::New();
    vtkDoubleArray *scalars = vtkDoubleArray::New();
    
    int i;
    double v[3];
    for (i=0; i<cmap_len; i++)
    {
      v[0] = dmin + (double)(dmax-dmin)*i / (double)cmap_len;
      v[1] = dmin + (double)(dmax-dmin)*i / (double)cmap_len;
      v[2] = dmin + (double)(dmax-dmin)*i / (double)cmap_len;
      points->InsertPoint(i,v);
      scalars->InsertTuple1(i,(double)histogram2[i]);
    }
    mesh2->SetPoints(points);
    points->Delete();
    mesh2->GetPointData()->SetScalars(scalars);
    scalars->Delete();        
    xyplot->AddInput((vtkDataSet *)mesh2);
        
    if(pargs.do_symmetric == 2)
    {
      xyplot->SetPlotColor(0,1,0,0);
      xyplot->SetPlotLabel(0, "B 2 A");
    }
    else
    {
      xyplot->SetPlotColor(1,1,0,0);
      xyplot->SetPlotLabel(1, "B 2 A");
    }
  }
  if(pargs.do_symmetric == 3)
  {
    vtkPoints *points = vtkPoints::New();
    vtkDoubleArray *scalars = vtkDoubleArray::New();
    
    int i;
    double v[3];
    for (i=0; i<cmap_len; i++)
    {
        v[0] = dmin + (double)(dmax-dmin)*i / (double)cmap_len;
        v[1] = dmin + (double)(dmax-dmin)*i / (double)cmap_len;
        v[2] = dmin + (double)(dmax-dmin)*i / (double)cmap_len;
        points->InsertPoint(i,v);
        scalars->InsertTuple1(i,(double)histogram1[i]+histogram2[i]);
    }
    mesh3->SetPoints(points);
    points->Delete();
    mesh3->GetPointData()->SetScalars(scalars);
    scalars->Delete();        
    xyplot->AddInput((vtkDataSet *)mesh3);
      
    xyplot->SetPlotColor(2,0,1,0);
    xyplot->SetPlotLabel(2, "A and B");      
  }
  
  ren->RemoveAllViewProps();
  ren->AddActor2D(xyplot);
  renWin->update();
  
  /*  double* range = xyplot->GetXRange();
  HistXMin->setText(QString::number(range[0]));
  HistXMax->setText(QString::number(range[1]));
  double* range2 = xyplot->GetYRange();
  HistYMin->setText(QString::number(range2[0]));
  HistYMax->setText(QString::number(range2[1]));*/

}

void MeshValmetControls::DoHistogram()
{
    double drange,off;
    double *serror;
    int i,bin_idx,n;
    int len = cmap_len;

    // This is a potentially slow operation
    QApplication::setOverrideCursor(Qt::waitCursor);
    
    delete [] histogram1;
    delete [] histogram2;
    delete [] hist;
    histogram1 = new int[len];
    histogram2 = new int[len];
    hist = new int[len];
    memset(histogram1, 0, sizeof(*histogram1)*len);
    memset(histogram2, 0, sizeof(*histogram2)*len);
    memset(hist, 0, sizeof(*hist)*len);
  
  if(pargs.do_symmetric == 1 || pargs.do_symmetric == 3)
  {
    n = model1.n_samples;
    if(pargs.signeddist == true)
    {
      drange = model1.max_error-model1.min_error;
      off = model1.min_error;
      serror = model1.fe[0].serror;
      for (i=0; i<n; i++) 
      {
         bin_idx = (int) ((serror[i]-off)/drange*len);
          if (bin_idx >= len) bin_idx = len-1;
          histogram1[bin_idx]++;
      }
    }
    else
    {
    	drange = model1.abs_max_error-model1.abs_min_error;
      off = model1.abs_min_error;
      serror = model1.fe[0].serror;
      for (i=0; i<n; i++) 
      {
         bin_idx = (int) ((fabs(serror[i])-off)/drange*len);
          if (bin_idx >= len) bin_idx = len-1;
          histogram1[bin_idx]++;
      }
    }
  }
  if(pargs.do_symmetric == 2 || pargs.do_symmetric == 3)
  {
    n = model2.n_samples;
    if(pargs.signeddist == true)
    {
      drange = model2.max_error-model2.min_error;
      off = model2.min_error;
      serror = model2.fe[0].serror;
      for (i=0; i<n; i++) 
      {
         bin_idx = (int) ((serror[i]-off)/drange*len);
          if (bin_idx >= len) bin_idx = len-1;
          histogram2[bin_idx]++;
      }
    }
    else
    {
    	drange = model2.abs_max_error-model2.abs_min_error;
      off = model2.abs_min_error;
      serror = model2.fe[0].serror;
      for (i=0; i<n; i++) 
      {
         bin_idx = (int) ((fabs(serror[i])-off)/drange*len);
          if (bin_idx >= len) bin_idx = len-1;
          histogram2[bin_idx]++;
      }
    }
  }
  if(pargs.do_symmetric == 3)
    for(i=0; i<len; i++)
    {
      hist[i] = histogram1[i] + histogram2[i];  
      //cout<<hist[i]<<"\t";
    }
  
  QApplication::restoreOverrideCursor();
}

void MeshValmetControls::UpdateHistogram()
{
  xyplot->SetXRange(HistXMin->text().toDouble(),HistXMax->text().toDouble());
  xyplot->SetYRange(HistYMin->text().toDouble(),HistYMax->text().toDouble());
  renWin->update();
    
}

void MeshValmetControls::OutputVertexErrors()
{
  QString outFile;
  if(commonPath == "")
    outFile = QFileDialog::getSaveFileName(
                      ".",
                      "Raw Errors Info (*.verrs)",
                      this,
                      "save file dialog",
                      "Choose a filename to save under" );
  else
          outFile = QFileDialog::getSaveFileName(
                      commonPath,
                      "Raw Errors  Info (*.verrs)",
                      this,
                      "save file dialog",
                      "Choose a filename to save under" );
  if(outFile != 0 )
  {
    if(outFile.right(6) != ".verrs")
      outFile = outFile + ".verrs";      
      
    if ( QFile::exists( outFile ) )
    {
      QMessageBox mb( "MeshValmet",
                        "Saving the file will overwrite the old file on disk.\n"
                        "Do you really want to save?",
                        QMessageBox::Information,
                        QMessageBox::Yes | QMessageBox::Default,
                        QMessageBox::No,
                        QMessageBox::Cancel | QMessageBox::Escape );
          mb.setButtonText( QMessageBox::Yes, "Save" );
          mb.setButtonText( QMessageBox::No, "Don't Save" );
          switch( mb.exec() ) 
          {
              case QMessageBox::Yes:
                    SaveVertexErrors(outFile);
                    break;
              case QMessageBox::No:
                    break;
              case QMessageBox::Cancel:
                    break;
          }
    }
    else
      SaveVertexErrors(outFile);
              
    QFileInfo fi(outFile);
    commonPath = fi.dirPath();

  }
}

void MeshValmetControls::SaveVertexErrors(QString outFile)
{
  int i,n;
  float *verror;
  FILE* fp;
  fp = fopen(outFile.latin1(),"w");
  if(pargs.do_symmetric == 1||pargs.do_symmetric == 3)
  {
    fprintf(fp,"%s->%s\n",fileName1.latin1(),fileName2.latin1());
    
    n = model1.mesh->num_vert;
    fprintf(fp,"%d\n",n);
    verror = model1.verror;
    
    if(pargs.signeddist == true)
    {
	    for (i=0; i<n; i++) 
	      fprintf(fp,"%f\n",verror[i]);
	    fprintf(fp,"\n\n");
	  }
	  else
	  {
	  	for (i=0; i<n; i++) 
	      fprintf(fp,"%f\n",fabs(verror[i]));
	    fprintf(fp,"\n\n");
	  }
  }
  if(pargs.do_symmetric == 2||pargs.do_symmetric == 3)
  {
    fprintf(fp,"%s->%s\n",fileName2.latin1(),fileName1.latin1());
    
    n = model2.mesh->num_vert;
    fprintf(fp,"%d\n",n);
    verror = model2.verror;
    
    if(pargs.signeddist == true)
    {
	    for (i=0; i<n; i++) 
	      fprintf(fp,"%f\n",verror[i]);
	    fprintf(fp,"\n\n");
	  }
	  else
	  {
	  	for (i=0; i<n; i++) 
	      fprintf(fp,"%f\n",fabs(verror[i]));
	    fprintf(fp,"\n\n");
	  }
  }
  fclose(fp);  
	
}

void MeshValmetControls::OutputRawErrors()
{
  QString outFile;
  if(commonPath == "")
    outFile = QFileDialog::getSaveFileName(
                      ".",
                      "Raw Errors Info (*.errs)",
                      this,
                      "save file dialog",
                      "Choose a filename to save under" );
  else
          outFile = QFileDialog::getSaveFileName(
                      commonPath,
                      "Raw Errors  Info (*.errs)",
                      this,
                      "save file dialog",
                      "Choose a filename to save under" );
  if(outFile != 0 )
  {
    if(outFile.right(5) != ".errs")
      outFile = outFile + ".errs";      
      
    if ( QFile::exists( outFile ) )
    {
      QMessageBox mb( "MeshValmet",
                        "Saving the file will overwrite the old file on disk.\n"
                        "Do you really want to save?",
                        QMessageBox::Information,
                        QMessageBox::Yes | QMessageBox::Default,
                        QMessageBox::No,
                        QMessageBox::Cancel | QMessageBox::Escape );
          mb.setButtonText( QMessageBox::Yes, "Save" );
          mb.setButtonText( QMessageBox::No, "Don't Save" );
          switch( mb.exec() ) 
          {
              case QMessageBox::Yes:
                    SaveRawErrors(outFile);
                    break;
              case QMessageBox::No:
                    break;
              case QMessageBox::Cancel:
                    break;
          }
    }
    else
      SaveRawErrors(outFile);
              
    QFileInfo fi(outFile);
    commonPath = fi.dirPath();

  }
}

void MeshValmetControls::SaveRawErrors(QString outFile)
{
  int i,n;
  double * serror;
  FILE* fp;
  fp = fopen(outFile.latin1(),"w");
  if(pargs.do_symmetric == 1||pargs.do_symmetric == 3)
  {
    fprintf(fp,"%s->%s\n",fileName1.latin1(),fileName2.latin1());
    n = model1.n_samples;
    fprintf(fp,"%d\n",n);
    serror = model1.fe[0].serror;
    
    if(pargs.signeddist == true)
    {
	    for (i=0; i<n; i++) 
	      fprintf(fp,"%f\n",serror[i]);
	    fprintf(fp,"\n\n");
	  }
	  else
	  {
	  	for (i=0; i<n; i++) 
	      fprintf(fp,"%f\n",fabs(serror[i]));
	    fprintf(fp,"\n\n");
	  }
  }
  if(pargs.do_symmetric == 2||pargs.do_symmetric == 3)
  {
    fprintf(fp,"%s->%s\n",fileName2.latin1(),fileName1.latin1());
    n = model2.n_samples;
    fprintf(fp,"%d\n",n);
    serror = model2.fe[0].serror;
    
    if(pargs.signeddist == true)
    {
	    for (i=0; i<n; i++) 
	      fprintf(fp,"%f\n",serror[i]);
	    fprintf(fp,"\n\n");
	  }
	  else
	  {
	  	for (i=0; i<n; i++) 
	      fprintf(fp,"%f\n",fabs(serror[i]));
	    fprintf(fp,"\n\n");
	  }
  }
  fclose(fp);  
}

void MeshValmetControls::SaveStats(QString outFile)
{
  
  FILE* fp = fopen(outFile.latin1(),"w");

  fprintf(fp,"Min\tMax\tVertexMean\tVertexSigma\tFaceMean\tFaceRMS\tHausdorff\tMSD\tMAD\tMedian\t25percentile\t68percentile\t75percentile\t95percentile\n");
  fprintf(fp,"%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",dmin, dmax, mean_error, sigma_error, face_mean_error, face_rms_error, hausdorff_error, mad_error, mad_error, median, percentile25, percentile68, percentile75,percentile95);

  //cout<<error_min<<":"<<error_max<<":"<<mean_error<<":"<<sigma_error<<":"<<face_mean_error<<":"<<face_rms_error<<":"<<hausdorff_error<<":"<<mad_error<<":"<<mad_error<<":"<<median<<":"<<dist_volume<<":"<<percentile25<<":"<<percentile68<<":"<<percentile75<<":"<<percentile95<<endl;
  fclose(fp);
}

void MeshValmetControls::OutputStats()
{
  QString outFile;
  if(commonPath == "")
    outFile = QFileDialog::getSaveFileName(
                      ".",
                      "Statistics Info (*.stat)",
                      this,
                      "save file dialog",
                      "Choose a filename to save under" );
  else
          outFile = QFileDialog::getSaveFileName(
                      commonPath,
                      "Statistics Info (*.stat)",
                      this,
                      "save file dialog",
                      "Choose a filename to save under" );
  if(outFile != 0 )
  {
    if(outFile.right(5) != ".stat")
      outFile = outFile + ".stat";      
      
    if ( QFile::exists( outFile ) )
    {
      QMessageBox mb( "MeshValmet",
                        "Saving the file will overwrite the old file on disk.\n"
                        "Do you really want to save?",
                        QMessageBox::Information,
                        QMessageBox::Yes | QMessageBox::Default,
                        QMessageBox::No,
                        QMessageBox::Cancel | QMessageBox::Escape );
          mb.setButtonText( QMessageBox::Yes, "Save" );
          mb.setButtonText( QMessageBox::No, "Don't Save" );
          switch( mb.exec() ) 
          {
              case QMessageBox::Yes:
                    SaveStats(outFile);
                    break;
              case QMessageBox::No:
                    break;
              case QMessageBox::Cancel:
                    break;
          }
    }
    else
      SaveStats(outFile);
              
    QFileInfo fi(outFile);
    commonPath = fi.dirPath();

  }
}

void MeshValmetControls::OutputHistogram()
{
  if( (pargs.do_symmetric == 1 || pargs.do_symmetric == 3) && histogram1 == NULL || (pargs.do_symmetric == 2 || pargs.do_symmetric == 3) && histogram2 == NULL || pargs.do_symmetric == 0)
    return;
  QString outFile;
  if(commonPath == "")
    outFile = QFileDialog::getSaveFileName(
                      ".",
                      "Histogram Info (*.hist)",
                      this,
                      "save file dialog",
                      "Choose a filename to save under" );
  else
          outFile = QFileDialog::getSaveFileName(
                      commonPath,
                      "Histogram Info (*.hist)",
                      this,
                      "save file dialog",
                      "Choose a filename to save under" );
  if(outFile != 0 )
  {
    if(outFile.right(5) != ".hist")
      outFile = outFile + ".hist";      
      
    if ( QFile::exists( outFile ) )
    {
      QMessageBox mb( "MeshValmet",
                        "Saving the file will overwrite the old file on disk.\n"
                        "Do you really want to save?",
                        QMessageBox::Information,
                        QMessageBox::Yes | QMessageBox::Default,
                        QMessageBox::No,
                        QMessageBox::Cancel | QMessageBox::Escape );
          mb.setButtonText( QMessageBox::Yes, "Save" );
          mb.setButtonText( QMessageBox::No, "Don't Save" );
          switch( mb.exec() ) 
          {
              case QMessageBox::Yes:
                    SaveHist(outFile);
                    break;
              case QMessageBox::No:
                    break;
              case QMessageBox::Cancel:
                    break;
          }
    }
    else
      SaveHist(outFile);
              
    QFileInfo fi(outFile);
    commonPath = fi.dirPath();

  }
}

void MeshValmetControls::SaveHist(QString filename)
{  
  QFile file(filename);
  if ( file.open( IO_WriteOnly ) ) 
  {
          QTextStream stream( &file );   
          if(pargs.do_symmetric == 1)
            stream << "#Histogram from Model1 \"" << pargs.m1_fname
                   << "\" to Model2 \"" << pargs.m2_fname << "\"\n";      
          else if(pargs.do_symmetric == 2)
            stream << "#Histogram from Model2 \"" << pargs.m2_fname
                   << "\" to Model1 \"" << pargs.m1_fname << "\"\n";
        else if(pargs.do_symmetric == 3)
          stream << "#Histgram of symmetric errors between Model1 \""
                 << pargs.m1_fname << "\" and Model2 \""
                 << pargs.m2_fname << "\"\n";
    stream << "#MAX: " << dmax << "\n";
    stream << "#MIN: " << dmin << "\n";
    stream << "#BINS: " << cmap_len << "\n";
    stream << "#STEP: " << (double)(dmax - dmin)/(double)cmap_len << "\n";
    stream << "#VERTEXMEAN: " << mean_error << "\n";
    stream << "#SIGMA: " << sigma_error << "\n";
    stream << "#MEADIAN: " << median << "\n";
    stream << "#25PERCENTILE: " << percentile25 <<  "\n";
    stream << "#68PERCENTILE: " << percentile68 <<  "\n";
    stream << "#75PERCENTILE: " << percentile75 <<  "\n";
    stream << "#95PERCENTILE: " << percentile95 <<  "\n";
    stream << "#FACEMEAN: " << face_mean_error << "\n";
    stream << "#RMS: " << face_rms_error << "\n";
    for(int i = 0; i < cmap_len; i++)
      if(pargs.do_symmetric == 1)
        stream << histogram1[i] << "\n";
      else if(pargs.do_symmetric == 2)
        stream << histogram2[i] << "\n";
      else if(pargs.do_symmetric == 3)
        stream << histogram1[i]+histogram2[i] << "\n";
          file.close();
      }
      else
        cout << "Cant open file for writing!" << endl;

}

void MeshValmetControls::UpdateColor()
{
  double mmax,mmin;
  if(ColormapMax->text() != "")
    mmax = ColormapMax->text().toDouble();
  else
      mmax = dmax;
  if(ColormapMin->text() != "")
      mmin = ColormapMin->text().toDouble();
  else
      mmin = dmin;
  int inter = IntervalSpinBox->value();
  
  lut = vtkColorTransferFunction::New();
    
    //Begin seting up my own lookup table
    //The Middle point always points to zero distance
    lut->AddRGBPoint(middle,0,1,0);
    //lut->AddRGBPoint(0.02,lookuptable[58].R/(float)255,lookuptable[58].G/(float)255,lookuptable[58].B/(float)255);
    //lut->AddRGBPoint(-0.02,lookuptable[198].R/(float)255,lookuptable[198].G/(float)255,lookuptable[198].B/(float)255);

    if(mmax-middle>0 && mmin-middle>=0)
    {
      int i;
      for(i = 1; i <=inter; i++)
        lut->AddRGBPoint((float)(middle+(mmax-middle)*i/(double)inter),lookuptable[128-i*22].R/(float)255,lookuptable[128-i*22].G/(float)255,lookuptable[128-i*22].B/(float)255);

    }
    else if(mmax-middle<=0 && mmin-middle<0)
    {
      int i;
      for(i = 1; i <=inter; i++)
        lut->AddRGBPoint((float)(middle+(mmin-middle)*i/(double)inter),lookuptable[128+i*22].R/(float)255,lookuptable[128+i*22].G/(float)255,lookuptable[128+i*22].B/(float)255);
    }
    else if(mmax-middle>0 && mmin-middle<0)
    {
      //The Positive outside distance
      int i;
      for(i = 1; i <=inter; i++)
      {
        lut->AddRGBPoint((float)(middle+(mmax-middle)*i/(double)inter),lookuptable[128-i*22].R/(float)255,lookuptable[128-i*22].G/(float)255,lookuptable[128-i*22].B/(float)255);
        //printf("%f\t%f\t%f\t%f\n",(float)(mmax*i/(double)5),lookuptable[128-i].R/(float)255,lookuptable[128-i].G/(float)255,lookuptable[128-i].B/(float)255);
      }
        
      //The Negative inside distance
      for(i = 1; i <=inter; i++)
        lut->AddRGBPoint((float)(middle+(mmin-middle)*i/(double)inter),lookuptable[128+i*22].R/(float)255,lookuptable[128+i*22].G/(float)255,lookuptable[128+i*22].B/(float)255);
  
    }
    
    vtkPolyDataMapper* mapper;
    
    if(pargs.do_symmetric == 1 || pargs.do_symmetric == 3)
    {
      mapper = (vtkPolyDataMapper *)modelActor1->GetMapper();
      mapper->SetLookupTable(lut);
      renWin1->update();
    }
    if(pargs.do_symmetric == 2 || pargs.do_symmetric == 3)
    {
      mapper = (vtkPolyDataMapper *)modelActor2->GetMapper();
      mapper->SetLookupTable(lut);
      renWin2->update();
    }
    
}

void MeshValmetControls::middleSliderMoved(int pos)
{
  double range = fabs(dmax)>fabs(dmin)?fabs(dmax):fabs(dmin);
  middle = range*(double)middleSlider->value()/110;
  TextLabel3->setText(QString::number(middle));
}

void QT_prog(void *out, int p) {
  QProgressDialog *qpd;
  qpd = (QProgressDialog*)out;
  qpd->setProgress(p<0 ? qpd->totalSteps() : p );
  qApp->processEvents();
}

