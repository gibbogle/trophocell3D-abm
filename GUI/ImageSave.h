#ifndef IMAGESAVE_H
#define IMAGESAVE_H


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
#include <vtkPNGWriter.h>
#include <vtkJPEGWriter.h>
#include <vtkTIFFWriter.h>
#include <vtkBMPWriter.h>
#include <vtkWindowToImageFilter.h>
#include <vtkSmartPointer.h>

class ImageSave : public QMainWindow
{
  Q_OBJECT
public:

  // Constructor/Destructor
  ImageSave(vtkSmartPointer<vtkRenderWindow> renWin);
  ~ImageSave() {};


  char msg[1024];
  vtkSmartPointer<vtkRenderWindow> renWin;
  void save(QString filename);

public slots:

protected:

protected slots:

private:


};


#endif // IMAGESAVE_H
