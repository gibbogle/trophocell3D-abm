#include "ui_SimpleView3DUI.h"
#include "SimpleView3DUI.h"
 
#include <vtkPolyDataMapper.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkSphereSource.h>
 
#include "vtkSmartPointer.h"

// From Hedgehog
#include <vtkVersion.h>
#include "vtkSmartPointer.h"
#include "vtkCamera.h"
#include "vtkFloatArray.h"
//#include "vtkHedgeHog.h"
#include "vtkMath.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyDataMapper.h"
#include "vtkProperty.h"
#include "vtkRenderWindow.h"
#include "vtkRenderer.h"
#include "vtkStructuredGrid.h"

#include <vtkArrowSource.h>
#include <vtkGlyph3D.h>
#include"vtkLight.h"
#include "log.h"
#include "ImageSave.h"

LOG_USE();

// From Hedgehog...

//------------------------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------------------------
SimpleView3D::SimpleView3D()
{
    this->ui = new Ui_SimpleView3D;
    this->ui->setupUi(this);
    LOG_QMSG("SimpleView3D");
}

//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
void SimpleView3D::create()
{
//    sprintf(msg,"Slice axis: %d fraction: %f use_strength: %d",axis,fraction,use_strength);
//    LOG_MSG(msg);
    int isgrid;
    vtkSmartPointer<vtkStructuredGrid> sgrid;
    vtkSmartPointer<vtkArrowSource> arrowSource;
//    vtkSmartPointer<vtkGlyph2D> glyphFilter;
    vtkSmartPointer<vtkGlyph3D> glyphFilter;
    vtkSmartPointer<vtkPolyDataMapper> sgridMapper;
    vtkSmartPointer<vtkActor> sgridActor;

    max_chemo = 2;
    // Create the structured grids.
    for (isgrid=0; isgrid<max_chemo; isgrid++) {
        sgrid_array[isgrid] = vtkSmartPointer<vtkStructuredGrid>::New();
    }
    float gmaxx, gmax[4], scaling;
    CreateGradientData(sgrid_array, chemo_select, gmax);
    gmaxx = 0;
    for (int i=0; i<max_chemo; i++) {
        if (gmax[i] > gmaxx)
            gmaxx = gmax[i];
    }
    if (scale == 0)
        scaling = 2.0/gmaxx;
    else
        scaling = 0.2*scale;
    // Create the usual rendering stuff
    renderer = vtkSmartPointer<vtkRenderer>::New();

    for (isgrid=0; isgrid<max_chemo; isgrid++) {
        sgrid = sgrid_array[isgrid];
        // We create a simple pipeline to display the data.
          // Setup the arrows
        arrowSource = vtkSmartPointer<vtkArrowSource>::New();
        arrowSource->SetShaftResolution(12);
        arrowSource->SetTipResolution(12);
        arrowSource->Update();

//        glyphFilter = vtkSmartPointer<vtkGlyph2D>::New();
        glyphFilter = vtkSmartPointer<vtkGlyph3D>::New();
        glyphFilter->SetSourceConnection(arrowSource->GetOutputPort());
        glyphFilter->OrientOn();
        glyphFilter->SetVectorModeToUseVector();
        glyphFilter->SetScaleFactor(scaling);
        glyphFilter->SetScaleModeToScaleByVector();
    //    glyphFilter->SetColorModeToColorByVector();
        glyphFilter->SetColorModeToColorByScalar();
        glyphFilter->SetInputConnection(sgrid->GetProducerPort());
        glyphFilter->Update();

        sgridMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        sgridMapper->SetInputConnection(glyphFilter->GetOutputPort());

        sgridActor_array[isgrid] = vtkSmartPointer<vtkActor>::New();
        sgridActor = sgridActor_array[isgrid];
        sgridActor->SetMapper(sgridMapper);
        if (isgrid == 0) {
            sgridActor->GetProperty()->SetColor(1.0,0,0);
        } else if (isgrid == 1) {
            sgridActor->GetProperty()->SetColor(0,1.0,0);
        } else if (isgrid == 2){
            sgridActor->GetProperty()->SetColor(0,0,1.0);
        } else {
            sgridActor->GetProperty()->SetColor(1.0,1.0,0);
        }
    }
    displayFields();
    renderer->SetBackground(1,1,1);
  //  renderer->SetBackground(0,0,0);
    renderer->ResetCamera();
    renderer->GetActiveCamera()->Elevation(60.0);
    renderer->GetActiveCamera()->Azimuth(30.0);
    renderer->GetActiveCamera()->Zoom(1.25);
  // VTK/Qt wedded
    renWin = this->ui->qvtkWidget_gradient->GetRenderWindow();
    renWin->AddRenderer(renderer);

  // Set up action signals and slots
  connect(this->ui->actionExit, SIGNAL(triggered()), this, SLOT(slotExit()));
};

/*
QString name = "checkBox_S1P";
QCheckBox *cb = this->ui->centralwidget->findChild<QCheckBox *>(name);
if (cb) {
    LOG_QMSG("Found a checkbox called: " + cb->objectName());
}
LOG_QMSG("All checkboxes");
QList<QCheckBox *> allCheckBoxes = this->ui->centralwidget->findChildren<QCheckBox *>();
for (int i=0; i<allCheckBoxes.size(); i++) {
    name = allCheckBoxes[i]->objectName();
    LOG_QMSG(name);
}
*/

//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
void SimpleView3D::slotExit()
{
  qApp->exit();
}

//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
void SimpleView3D::resizeEvent(QResizeEvent *event)
{
    int rwsize[2];

    QSize qsize = this->size();
    rwsize[0] = qsize.width();
    rwsize[1] = qsize.height();
    GetRenderWindow()->SetSize(rwsize);
}

//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
void SimpleView3D::saveImage(void)
{
    QString fname = "";
    ImageSave *myImageSave = new ImageSave(GetRenderWindow());
    myImageSave->save(fname);
}

//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
void SimpleView3D::setParameters()
{
    chemo_select[0] = 1;
    chemo_select[1] = 0;
    chemo_displayed[0] = false;
    chemo_displayed[1] = false;
    scale = 0;
    use_strength = 0;
}
//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
void SimpleView3D::chooseParameters(void)
{
    chemo_select[0] = 1;
    chemo_select[1] = 0;
    chemo_displayed[0] = false;
    chemo_displayed[1] = false;

    bool ok;
    QString text;
    scale = QInputDialog::getDouble(this, tr("QInputDialog::getDouble()"),
                                       tr("Scaling factor: (0 normalizes scale)"), 0.0, 0, 100, 2, &ok);
    if (scale == 0) {
        text = "Normalized vectors";
    } else {
        sprintf(msg,"Vector scaling: %5.2f",scale);
        text = QString(msg);
    }
    this->ui->label_scaling->setText(text);

    int ret = QMessageBox::question(this, tr("Chemokine relative strength"),
                                        tr("Do you want to multiply gradients by the chemokine strength?"),
                                        QMessageBox::Yes | QMessageBox::No, QMessageBox::No);
    if (ret == QMessageBox::Yes) {
        use_strength = 1;
        text = "Using relative strength";
    } else {
        use_strength = 0;
        text = "Not using relative strength";
    }
    this->ui->label_strength->setText(text);
}

//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
void SimpleView3D::aimCamera(void)
{
    double cp[3], fp[3], up[3];
    vtkSmartPointer<vtkCamera> camera;

    axis = 3;
//    return;     // perhaps not needed

    camera = renderer->GetActiveCamera();
    camera->GetFocalPoint(fp);
    camera->GetPosition(cp);
    camera->GetViewUp(up);
    sprintf(msg,"AimCamera:\n  position: %f %f %f\n  focus: %f %f %f\n  up: %f %f %f",cp[0],cp[1],cp[2],fp[0],fp[1],fp[2],up[0],up[1],up[2]);
    LOG_MSG(msg);
    if (axis == 1) {
        cp[0] = 2*fp[0];
        cp[1] = fp[1];
        cp[2] = fp[2];
        camera->SetPosition(cp);
        up[0] = 0;
        up[1] = 1;
        up[2] = 0;
        camera->SetViewUp(up);
    } else if (axis == 2) {
        cp[0] = fp[0];
        cp[1] = 3*fp[1];
        cp[2] = fp[2];
        camera->SetPosition(cp);
        up[0] = 0;
        up[1] = 0;
        up[2] = -1;
        camera->SetViewUp(up);
    } else if (axis == 3) {
        cp[0] = fp[0];
        cp[1] = fp[1];
        cp[2] = 2*fp[2];
        camera->SetPosition(cp);
        up[0] = 0;
        up[1] = 1;
        up[2] = 0;
        camera->SetViewUp(up);
    }
// This suppresses rotation
//    vtkSmartPointer<vtkInteractorStyleImage> imageStyle = vtkSmartPointer<vtkInteractorStyleImage>::New();
//    renWin->GetInteractor()->SetInteractorStyle(imageStyle);
    renderer->ResetCamera();
    renderer->Render();
}

//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
void SimpleView3D::displayFields(void)
{
    int ichemo;
    for (ichemo=0; ichemo<max_chemo; ichemo++) {
        if (chemo_select[ichemo] == 0) {
            if (chemo_displayed[ichemo]) {
                renderer->RemoveActor(sgridActor_array[ichemo]);
                chemo_displayed[ichemo] = false;
                LOG_QMSG("Removed actor");
            }
        } else {
            if (!chemo_displayed[ichemo]) {
                renderer->AddActor(sgridActor_array[ichemo]);
                chemo_displayed[ichemo] = true;
                LOG_QMSG("Added actor");
            }
        }
    }
    iren = this->ui->qvtkWidget_gradient->GetInteractor();
    iren->Render();
    // Set up the lighting.
    vtkSmartPointer<vtkLight> light0 = vtkSmartPointer<vtkLight>::New();
    light0->SetPosition(100,0,0);
    renderer->AddLight(light0);
    vtkSmartPointer<vtkLight> light1 = vtkSmartPointer<vtkLight>::New();
    light1->SetPosition(-100,0,0);
    renderer->AddLight(light1);
    vtkSmartPointer<vtkLight> light2 = vtkSmartPointer<vtkLight>::New();
    light2->SetPosition(0,100,0);
    renderer->AddLight(light2);
    vtkSmartPointer<vtkLight> light3 = vtkSmartPointer<vtkLight>::New();
    light3->SetPosition(0,-100,0);
    renderer->AddLight(light3);
    vtkSmartPointer<vtkLight> light4 = vtkSmartPointer<vtkLight>::New();
    light4->SetPosition(0,0,100);
    renderer->AddLight(light4);
    vtkSmartPointer<vtkLight> light5 = vtkSmartPointer<vtkLight>::New();
    light5->SetPosition(0,0,-100);
    renderer->AddLight(light5);
}

//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
void SimpleView3D::stateChanged_CheckBox_CCL3(void)
{
    if (!chemo_used[0]) return;
    if (chemo_select[0] == 1)
        chemo_select[0] = 0;
    else
        chemo_select[0] = 1;
    sprintf(msg,"CCL3 select is now: %d\n",chemo_select[0]);
    LOG_MSG(msg);
    displayFields();
}

//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
void SimpleView3D::stateChanged_CheckBox_other(void)
{
    if (!chemo_used[1]) return;
    if (chemo_select[1] == 1)
        chemo_select[1] = 0;
    else
        chemo_select[1] = 1;
    sprintf(msg,"Other select is now: %d\n",chemo_select[1]);
    LOG_MSG(msg);
    displayFields();
}

//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
vtkSmartPointer<vtkRenderWindow> SimpleView3D::GetRenderWindow(void)
{
    vtkSmartPointer<vtkRenderWindow> renWin = this->ui->qvtkWidget_gradient->GetRenderWindow();
    return renWin;
}

//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
void SimpleView3D::ShowSize(int *size)
{
    sprintf(msg,"Window size: %d %d\n",size[0],size[1]);
    LOG_MSG(msg);
}

//------------------------------------------------------------------------------------------------
// Note that the amount of data retrieved from the DLL depends on nchem_used and nsites, therefore
// these values must be retrieved before gradient_array is allocated.
// We need a way to determine (for display and selection) the chemokine names corresponding to
// the identifying numbers.  For now these can be hard-coded:
// 0 = CCL3/4
// 1 = other
// The array chem_used[] conveys which of these are in use.  (This info is also available from
// the UI).
//------------------------------------------------------------------------------------------------
void SimpleView3D::CreateGradientData(vtkSmartPointer<vtkStructuredGrid> sgrid_array[], int chemo_select[], float gmax[])
{
    int k, iga;
    float x[3], v[3], g;
    static int dims[3]={50,50,50};
    int nchemo_used, ndata;
    int chemo_simulated[4], nsites;
    float *gradient_array;
    int ichemo;

    for (ichemo=0; ichemo<max_chemo; ichemo++) {
        sgrid_array[ichemo]->SetDimensions(dims);
    }
//  get_gradient_info(chemo_simulated, &nsites);
  nchemo_used = 0;
  for (ichemo=0; ichemo<max_chemo; ichemo++) {
      if (chemo_simulated[ichemo] == 1) {
          nchemo_used++;
          chemo_used[ichemo] = true;
      } else {
          chemo_used[ichemo] = false;
      }
  }
  if (nchemo_used == 0) return;
  ndata = nsites*(3 + 3*max_chemo);
  sprintf(msg,"nchem_used: %d nsites: %d ndata: %d",nchemo_used,nsites,ndata);
  LOG_MSG(msg);
  gradient_array = (float *)malloc(ndata*sizeof(float));
//  get_gradients(chemo_simulated, &nsites, gradient_array, &use_strength);

  vtkSmartPointer<vtkFloatArray> vectors;
  vtkSmartPointer<vtkPoints> points;
  for (ichemo=0; ichemo<max_chemo; ichemo++) {
      if (!chemo_used[ichemo]) continue;

      // We create the points and vectors.
      vectors = vtkSmartPointer<vtkFloatArray>::New();
      vectors->SetNumberOfComponents(3);
      vectors->SetNumberOfTuples(nsites);

      points = vtkSmartPointer<vtkPoints>::New();
      points->Allocate(nsites);
      gmax[ichemo] = 0;
      for (k=0; k<nsites; k++) {
          iga = k*(3 + 3*max_chemo);
          x[0] = gradient_array[iga];
          x[1] = gradient_array[iga+1];
          x[2] = gradient_array[iga+2];
          iga += 3 + 3*ichemo;
          v[0] = gradient_array[iga];
          v[1] = gradient_array[iga+1];
          v[2] = gradient_array[iga+2];
          g = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
          if (g > gmax[ichemo]) gmax[ichemo] = g;
          points->InsertPoint(k,x);
          vectors->InsertTuple(k,v);
    //      sprintf(msg,"site: %d  %d %d %d v: %f %f %f",k,int(x[0]),int(x[1]),int(x[2]),v[0],v[1],v[2]);
    //      LOG_MSG(msg);
      }
      sgrid_array[ichemo]->SetPoints(points);
      sgrid_array[ichemo]->GetPointData()->SetVectors(vectors);
      sprintf(msg,"ichemo: %d gmax: %f\n",ichemo,gmax[ichemo]);
      LOG_MSG(msg);
  }
  free(gradient_array);
}

//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
void SimpleView3D::CreateTestData(vtkStructuredGrid* sgrid)
{
  int i, j, k, kOffset, jOffset, offset, myoffset;
  float x[3], v[3], rMin=0.5, rMax=1.0, deltaRad, deltaZ;
  float radius, theta;
  static int dims[3]={13,11,11};

  sgrid->SetDimensions(dims);

  // We also create the points and vectors. The points
  // form a hemi-cylinder of data.
  vtkSmartPointer<vtkFloatArray> vectors = vtkSmartPointer<vtkFloatArray>::New();
  vectors->SetNumberOfComponents(3);
  vectors->SetNumberOfTuples(dims[0]*dims[1]*dims[2]);

  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  points->Allocate(dims[0]*dims[1]*dims[2]);

  myoffset = 0;
  deltaZ = 2.0 / (dims[2]-1);
  deltaRad = (rMax-rMin) / (dims[1]-1);
  v[2]=0.0;
    for (j=0; j<dims[1]; j++)
      {
      radius = rMin + j*deltaRad;
      jOffset = j * dims[0];
      for (i=0; i<dims[0]; i++)
        {
  for ( k=0; k<dims[2]; k++)
    {
    x[2] = -1.0 + k*deltaZ;
    kOffset = k * dims[0] * dims[1];
        theta = i * vtkMath::RadiansFromDegrees(15.0);
        x[0] = radius * cos(theta);
        x[1] = radius * sin(theta);
        v[0] = -x[1];
        v[1] = x[0];
        offset = i + jOffset + kOffset;
        points->InsertPoint(myoffset,x);
        vectors->InsertTuple(myoffset,v);
        myoffset++;
        }
      }
    }
  sgrid->SetPoints(points);

  sgrid->GetPointData()->SetVectors(vectors);
}

//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
void SimpleView3D::makeFrame(int i)
{
    setParameters();
    create();
    show();
    aimCamera();
    char fname_str[64];
    sprintf(fname_str,"D:/DCvisits/execution/image/frame%04d.png",i);
    LOG_MSG(fname_str);
    QString fname = QString(fname_str);
    ImageSave *is = new ImageSave(GetRenderWindow());
    is->save(fname);
    delete is;
}
