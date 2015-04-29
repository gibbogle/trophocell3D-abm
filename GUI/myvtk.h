// myvtk.h
#ifndef MYVTK_H
#define MYVTK_H

#include <QtGui>
#include <QtCore>
#include <QIODevice>
#include <QVTKWidget.h>
#include <vtkRenderer.h> 
#include <vtkRenderWindow.h>
#include "vtkSphereSource.h"
#include "vtkCylinderSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkProperty.h"
#include "vtkCamera.h"
//#include <vtkMPEG2Writer.h>
#include <vtkPNGWriter.h>
#include <vtkJPEGWriter.h>
#include <vtkTIFFWriter.h>
#include <vtkBMPWriter.h>
#include <vtkWindowToImageFilter.h>
#include <vtkSmartPointer.h>
#include <vtkImageCast.h>

//#include <vtkConfigure.h>

using namespace std;

struct cell_pos {
	int tag;
	int x, y, z;
	double diameter;
	double state;
};
typedef cell_pos CELL_POS;

struct bond_pos {
	int TCtag;
	int DCtag;
};
typedef bond_pos BOND_POS;

class MyVTK
{
public:
    MyVTK(QWidget *, QWidget *);
	~MyVTK();

    void test_canvas(QWidget *test_page);
	void read_cell_positions(QString, QString, bool);
	void get_cell_positions(bool fast);
	void init();
	void cleanup();
    void renderCells(bool);
    void process_Tcells();
    void process_Dcells(bool);
    void process_bonds();
    void makeLines();
	void getTCColor(int state, double *r, double *g, double *b);
	void setColor(double *r, double *g, double *b, int col[]);
	bool startPlayer(QString, QTimer *, bool);
	bool nextFrame();
	void pause();
	void playon();
	void saveSnapshot(QString, QString);
    void record(QString, int number);
	void stop();

	QList<CELL_POS > TCpos_list;
	QList<CELL_POS > DCpos_list;
	QList<BOND_POS > bondpos_list;
	QList<vtkActor *> T_Actor_list;
	QList<vtkActor *> D_Actor_list;
	QList<vtkActor *> B_Actor_list;

	QVTKWidget* qvtkWidget;
	vtkRenderWindow *renWin;	
	vtkRenderer* ren;
	vtkRenderWindowInteractor * iren;
	vtkPolyDataMapper *TcellMapper;
	vtkPolyDataMapper *DcellMapper;
	vtkPolyDataMapper *bondMapper;
    vtkPolyDataMapper *cylMapper;
//	vtkMPEG2Writer *mpg;
//	vtkSmartPointer<vtkPNGWriter> writer;
//	vtkSmartPointer<vtkBMPWriter> writer;
	vtkSmartPointer<vtkJPEGWriter> writer;
//	vtkSmartPointer<vtkTIFFWriter> writer;
	vtkSmartPointer<vtkImageCast> castFilter;
	vtkWindowToImageFilter *w2i;
	vtkWindowToImageFilter *w2img;
//	vtkSmartPointer<vtkPNGWriter> pngwriter;
//	vtkSmartPointer<vtkJPEGWriter> jpgwriter;

	char msg[2048];
	double zoomlevel;
	double Pi;
	bool DCmotion;
	bool DCfade;
	bool first_VTK;
	bool playing;
	bool paused;
	bool save_image;
	QString casename;
	int framenum;
	QTimer *timer;
	QString infile;
	QFile *playerData;
	QTextStream *playerStream;

//	double *TCColor, *DCColor;
};

#endif
