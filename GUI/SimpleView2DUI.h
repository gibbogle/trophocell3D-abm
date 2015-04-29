#ifndef SimpleView2DUI_H
#define SimpleView2DUI_H
 
#include "vtkSmartPointer.h"
#include <QMainWindow>
#include "vtkArrowSource.h"
#include "vtkStructuredGrid.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkInteractorStyleImage.h"
#include "vtkActor.h"
#include "log.h"
#include <QinputDialog>
#include <QMessageBox>

#ifdef __cplusplus
extern "C" {
#endif
//    void get_gradient2d_info(int *, int *, int *, float *);
//    void get_gradients2d(int *, int *, float *, int *, float *, int *);
#ifdef __cplusplus
}
#endif

// Forward Qt class declarations
class Ui_SimpleView2D;
 
class SimpleView2D : public QMainWindow
{
  Q_OBJECT
public:
 
  // Constructor/Destructor
  SimpleView2D();
  ~SimpleView2D() {};

  vtkSmartPointer<vtkRenderWindow> GetRenderWindow();
  void ShowSize(int *);
  void displayFields(void);
  void aimCamera(void);
  void create();
  void chooseParameters();
  void setParameters();
  void makeFrame(int i);

  int max_chemo;
  char msg[1024];
  int axis;
  float fraction;
  double scale;
  int use_strength;
  int chemo_select[2];  // Currently hard-coded for 2 chemokines
  bool chemo_displayed[2];
  bool chemo_used[2];
  vtkSmartPointer<vtkRenderer> renderer;
  vtkSmartPointer<vtkRenderWindow> renWin;
  vtkRenderWindowInteractor * iren;
  vtkSmartPointer<vtkArrowSource> arrowSource_array[2];
  vtkSmartPointer<vtkStructuredGrid> sgrid_array[2];
  vtkSmartPointer<vtkActor> sgridActor_array[2];

public slots:
 
  virtual void slotExit();
  void stateChanged_CheckBox_CCL3();
  void stateChanged_CheckBox_other();
  void saveImage();

protected:
 
protected slots:
 
private:

  void CreateGradientData(vtkSmartPointer<vtkStructuredGrid> sgrid_array[], int chemo_select[], float gmax[]);
  void CreateTestData(vtkStructuredGrid* sgrid);
  void resizeEvent(QResizeEvent *);

  // Designer form
  Ui_SimpleView2D *ui;
};
 
#endif // SimpleView2DUI_H
