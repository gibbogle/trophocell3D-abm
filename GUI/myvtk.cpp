// myvtk.cpp

#include <vtkCamera.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkObjectFactory.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkTransform.h>

#ifdef _WIN32
#include "windows.h"
#endif
#include "myvtk.h"
#include "log.h"
#include "global.h"

LOG_USE();

// Define interaction style
class MouseInteractorStyle4 : public vtkInteractorStyleTrackballCamera
{
  public:
	static MouseInteractorStyle4* New();
	vtkTypeMacro(MouseInteractorStyle4, vtkInteractorStyleTrackballCamera);

	virtual void OnLeftButtonDown()
	{
      Global::leftb = true;
	  // Forward events
	  vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
	}

	virtual void OnMiddleButtonDown()
	{
//	  std::cout << "Pressed middle mouse button." << std::endl;
	  // Forward events
	  vtkInteractorStyleTrackballCamera::OnMiddleButtonDown();
	}

	virtual void OnRightButtonDown()
	{
//	  std::cout << "Pressed right mouse button." << std::endl;
	  // Forward events
	  vtkInteractorStyleTrackballCamera::OnRightButtonDown();
	}

	virtual void OnLeftButtonUp()
	{
//	  std::cout << "Released left mouse button." << std::endl;
//	  LOG_QMSG("Released left mouse button.");
      Global::leftb = false;
	  // Forward events
	  vtkInteractorStyleTrackballCamera::OnLeftButtonUp();
	}

	virtual void OnMiddleButtonUp()
	{
//	  std::cout << "Released middle mouse button." << std::endl;
	  // Forward events
	  vtkInteractorStyleTrackballCamera::OnMiddleButtonUp();
	}

	virtual void OnRightButtonUp()
	{
//	  std::cout << "Released right mouse button." << std::endl;
	  // Forward events
	  vtkInteractorStyleTrackballCamera::OnRightButtonUp();
	}

};

vtkStandardNewMacro(MouseInteractorStyle4);

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
MyVTK::MyVTK(QWidget *page, QWidget *test_page, double radius, double length)
{
    zoomlevel = 0.4;
	double backgroundColor[] = {0.0,0.0,0.0};

    test_canvas(test_page);
	Pi = 4*atan(1.0);
    Global::leftb = false;
	qvtkWidget = new QVTKWidget(page,QFlag(0));
	LOG_MSG("Created a new QVTKWidget");
	QVBoxLayout *layout = new QVBoxLayout;
    layout->addWidget(qvtkWidget);

	// Associate the layout with page_VTK
    page->setLayout(layout);

	// Create a renderer, and add it to qvtkWidget's render window.
	// The renderer renders into the render window. 
	ren = vtkRenderer::New();     
    renWin = qvtkWidget->GetRenderWindow();
    renWin->AddRenderer(ren);
	ren->SetBackground(backgroundColor);
//	ren->SetBackground(0.1, 0.2, 0.4);		// backgroundColor
	ren->ResetCamera();
	iren = qvtkWidget->GetInteractor();

	vtkSmartPointer<MouseInteractorStyle4> style = vtkSmartPointer<MouseInteractorStyle4>::New();
	iren->SetInteractorStyle( style );

	iren->Initialize();

	// Create mappers
	vtkSphereSource *Tcell = vtkSphereSource::New();
    Tcell->SetThetaResolution(12);
    Tcell->SetPhiResolution(12);
    Tcell->SetRadius(0.5);
	TcellMapper = vtkPolyDataMapper::New();

	TcellMapper->SetInputConnection(Tcell->GetOutputPort());

    vtkCylinderSource *cyl = vtkCylinderSource::New();
    cyl->SetCenter(0.0, 0.0, 0.0);
    cyl->SetRadius(radius);
    cyl->SetHeight(length);
    cyl->SetResolution(100);
    cyl->CappingOff();

    vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
    transform->RotateWXYZ(90., 1., 0., 0.);
    transform->Translate(0., length/2, 0.);
    vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
    transformFilter->SetTransform(transform);
    transformFilter->SetInputConnection(cyl->GetOutputPort());
    transformFilter->Update();

    cylMapper = vtkPolyDataMapper::New();
    cylMapper->SetInputConnection(transformFilter->GetOutputPort());
    cylactor = vtkSmartPointer<vtkActor>::New();
    cylactor->SetMapper(cylMapper);
//    cylactor->GetProperty()->SetOpacity(0.5); // translucent !!!
    cylactor->GetProperty()->SetColor(0, 1, 0);
    ren->AddActor(cylactor);
//    ren->SetBackground(.1, .3,.2); // Background color

	// Create image filter for save Snapshot()
	w2img = vtkWindowToImageFilter::New();

	first_VTK = true;
	playing = false;
	paused = false;
    framenum = 0;

    ren->GetActiveCamera()->SetFocalPoint(0.,0.,300.);
    ren->GetActiveCamera()->SetPosition(0.,-3000, 0.);
//    ren->GetActiveCamera()->SetPosition(0.,-15*radius, 0.);
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
MyVTK::~MyVTK()
{
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void MyVTK::test_canvas(QWidget *test_page)
{
    QGraphicsScene* scene = new QGraphicsScene(QRect(0, 0, 200, 300));
    QBrush brush;
    QGraphicsTextItem *text;
    brush.setColor(QColor(255,0,0));
    brush.setStyle(Qt::SolidPattern);
    scene->addRect(10,10,20,20,Qt::NoPen, brush);
    text = scene->addText("Red square");
    text->setPos(35, 10);
    brush.setColor(QColor(0,255,0));
    scene->addEllipse(10,40,20,20,Qt::NoPen, brush);
    text = scene->addText("Green circle");
    text->setPos(35, 40);
    QGraphicsView* view = new QGraphicsView(test_page);
    view->setScene(scene);
    view->setGeometry(QRect(0, 0, 220, 320));
    view->show();
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void MyVTK::setWallOpacity(double opacity)
{
    cylactor->GetProperty()->SetOpacity(opacity);

}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void MyVTK::get_cell_positions(bool dummy)
{
    double TC_diam = 10.0;
	TCpos_list.clear();
	DCpos_list.clear();
	bondpos_list.clear();
    sprintf(msg,"nTC_list: %d TC_diam: %f",Global::nTC_list,TC_diam);
    LOG_MSG(msg);
    for (int i=0; i<Global::nTC_list; i++) {
		int j = 5*i;
		CELL_POS cp;
        cp.tag = Global::TC_list[j];
        cp.x = Global::TC_list[j+1];
        cp.y = Global::TC_list[j+2];
        cp.z = Global::TC_list[j+3];
        cp.state = Global::TC_list[j+4];
        cp.diameter = TC_diam;
		TCpos_list.append(cp);
        sprintf(msg,"cell pos: i: %d x,y,z: %d %d %d state: %d diam: %f",i,cp.x,cp.y,cp.z,cp.state,cp.diameter);
        LOG_MSG(msg);
	}
}

//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------
void MyVTK::init()
{
    T_Actor_list.clear();
}

//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------
void MyVTK::cleanup()
{
	int i;
	vtkActor *actor;
	LOG_MSG("VTK cleanup");
	for (i = 0; i<T_Actor_list.length(); i++) {
		actor = T_Actor_list[i];
        ren->RemoveActor(actor);
	}
    T_Actor_list.clear();
	first_VTK = true;	
}

//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------
void MyVTK::renderCells(bool dummy)
{
    process_Tcells();
	if (first_VTK) {
		LOG_MSG("Initializing the renderer");
	}
    iren->Render();
	first_VTK = false;	
}

//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------
void MyVTK::makeLines()
{
    double x, y, z, y0, z0, theta, beta, pos[3], s[3], PI=3.14159;
    double r=0, g=1, b=0;
    vtkActor *actor;

    x = Global::tubeLength/2.0;
    y0 = Global::tubeRadius + 1.5;
    z0 = Global::tubeRadius + 1.5;
    beta = 90.0;
    s[1] = Global::tubeLength;
    s[0] = 1;
    s[2] = 2;
    for (int i=0; i<Global::nLines; i++) {
        theta = i*2*PI/Global::nLines;
        y = y0 + Global::tubeRadius*cos(theta);
        z = z0 + Global::tubeRadius*sin(theta);
        actor = vtkActor::New();
        actor->SetMapper(cylMapper);
        pos[0] = x;
        pos[1] = y;
        pos[2] = z;
        actor->SetPosition(pos);
        actor->GetProperty()->SetColor(r, g, b);
        actor->SetScale(s);
        actor->RotateWXYZ(beta,0,0,1);
        ren->AddActor(actor);
//        printf("line: %d %f %f %f\n",i,x,y,z);
    }
}

//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------
void MyVTK::process_Tcells()
{
	int i, tag;
	double r, g, b;
	CELL_POS cp;
	vtkActor *actor;

    int na = T_Actor_list.length();
    int np = TCpos_list.length();
    int n = na;
	for (i=0; i<np; i++) {
        cp = TCpos_list[i];
        tag = cp.tag;
        n = max(tag+1,n);
	}
    bool *active;
	active = new bool[n];
	for (i=0; i<n; i++)
		active[i] = false;
	for (i=0; i<np; i++) {
        cp = TCpos_list[i];
        tag = cp.tag;
        active[tag] = true;
		if (tag >= na) {   // need to add actor, and possibly fill gaps
			if (tag > na) {
                for (int j=na; j<tag; j++)	//j in range(na,tag):
                    T_Actor_list.append(0);
			}
			actor = vtkActor::New();
            actor->SetMapper(TcellMapper);
 //           actor->GetProperty()->SetColor(TCColor);
            ren->AddActor(actor);
            T_Actor_list.append(actor);
            na = tag + 1;
		}
		getTCColor(cp.state,&r,&g,&b);
        actor = T_Actor_list[tag];
        actor->GetProperty()->SetColor(r, g, b);
        actor->SetPosition(cp.x, cp.y, cp.z);
        actor->SetScale(cp.diameter);
		if (actor != 0) 
			actor->SetPosition(cp.x, cp.y, cp.z);
		else {
			sprintf(msg,"T_actor = 0: %d",tag);
			LOG_MSG(msg);
			exit(1);
		}
	}

	for (int k=0; k<na; k++) {	// k in range(0,na):
		if (T_Actor_list[k] != 0 && !active[k]) {     // need to remove actor from list
            actor = T_Actor_list[k];
            ren->RemoveActor(actor);
            T_Actor_list[k] = 0;
		}
	}
}

//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------
void MyVTK::getTCColor(int state, double *r, double *g, double *b)
{
//	bool CD4;
//	int deep_blue[]   = {30,20,255};
//	int deep_green[]  = {0,150,0};
//	int light_blue[]  = {0,200,255};
//	int light_green[] = {50,255,150};
//	int purple[]      = {200,30,255};
//	int yellow[]      = {255,255,30};
  int red[]         = {255,0,0};

    setColor(r,g,b,red);
}

//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------
void MyVTK::setColor(double *r, double *g, double *b, int col[])
{
	*r = col[0]/255.;
	*g = col[1]/255.;
	*b = col[2]/255.;
}

//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------
void MyVTK::process_Dcells(bool redo)
{
	int i, tag;
	CELL_POS cp;
	double antigen, color[3];
	vtkActor *actor;
	double minlevel = 0.3;

	double DCColor[] = {1.0,0.0,0.0};
    int na = D_Actor_list.length();
    int np = DCpos_list.length();
    int n = na;
	for (i=0; i<np; i++) {
        cp = DCpos_list[i];
        tag = cp.tag;
        n = max(tag+1,n);
	}
    bool *active;
	active = new bool[n];
	for (i=0; i<n; i++)
		active[i] = false;

	for (i=0; i<np; i++) {
        cp = DCpos_list[i];
        tag = cp.tag;
        active[tag] = true;
        bool newDC = false;
		if (tag >= na) {   // need to add actor, and possibly fill gaps
			if (tag > na) {
                for (int j=na; j<tag; j++)	//j in range(na,tag):
                    D_Actor_list.append(0);
			}
			actor = vtkActor::New();
            actor->SetMapper(DcellMapper);
            actor->GetProperty()->SetColor(DCColor);

            ren->AddActor(actor);
            D_Actor_list.append(actor);
            na = tag + 1;
            newDC = true;
		} else {
			actor = D_Actor_list[tag];
		}
		if (redo || newDC || DCmotion || DCfade) {
			if (actor != 0) {
				if (DCfade) {
					antigen = cp.state;
					color[0] = (minlevel + (1-minlevel)*antigen)*DCColor[0];
					color[1] = (minlevel + (1-minlevel)*antigen)*DCColor[1];
					color[2] = (minlevel + (1-minlevel)*antigen)*DCColor[2];
					actor->GetProperty()->SetColor(color);
				}
				actor->SetPosition(cp.x, cp.y, cp.z);
			} else {
				sprintf(msg,"D_actor = 0: %d",tag);
				LOG_MSG(msg);
				exit(1);
			}
		}
	}

	for (int k=0; k<na; k++) {	// k in range(0,na):
		if (D_Actor_list[k] != 0 && !active[k]) {     // need to remove actor from list
            actor = D_Actor_list[k];
            ren->RemoveActor(actor);
            D_Actor_list[k] = 0;
		}
	}
}


//---------------------------------------------------------------------------------------------
// A cylinder is created orientated along the y-axis, i.e. along b = (0,1,0)
// To change the orientation to the vector v, we first create a vector r
// normal to both b and v: r = bxv, this will be the axis of rotation.
// We now need to rotate the cylinder by theta about r, where theta is the angle between
// b and v, i.e. sin(theta) = |r|/(|b||v|) = |r|/|v|
// We can now use actor.RotateWXYZ(theta,r[0],r[1],r[2]) where theta is in degrees
// What is bxv when b = (0,1,0) and v = (v0,v1,v2)?
// r = [v[2],0,-v[0]]
//---------------------------------------------------------------------------------------------
void MyVTK::process_bonds()
{
	int i, j;
	BOND_POS bp;
	vtkActor *actor, *T_actor, *D_actor;
	double bpos[3], v[3];
	double Pi = 3.15159;
	double *tcpos, *dcpos;
	double bondColor[] = {0.5,0.0,0.0};

    int na = B_Actor_list.length();
    int np = bondpos_list.length();

    // First remove all old bonds (strictly speaking we should remove only those not in the new list)

	for (int k=0; k<na; k++) {
        ren->RemoveActor(B_Actor_list[k]);
	}

    B_Actor_list.clear();    

	for (i=0; i<np; i++) {
        bp = bondpos_list[i];
		actor = vtkActor::New();
        actor->SetMapper(bondMapper);
		actor->GetProperty()->SetColor(bondColor);
		T_actor = T_Actor_list[bp.TCtag];
		if (T_actor != 0)
	        tcpos = T_actor->GetPosition();
		else {
			sprintf(msg,"T_actor = 0 in bond: %d %d",i,bp.TCtag);
			LOG_MSG(msg);
			exit(1);
		}
		D_actor = D_Actor_list[bp.DCtag];
		if (D_actor != 0)
	        dcpos = D_actor->GetPosition();
		else {
			sprintf(msg,"D_actor = 0 in bond: %d %d",i,bp.DCtag);
			LOG_MSG(msg);
			exit(1);
		}
	
		for (j=0; j<3; j++) {
            bpos[j] = (tcpos[j] + dcpos[j])/2;
            v[j] = tcpos[j] - dcpos[j];
		}
        double v_mod = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
		double s[] = {1, v_mod, 1};
        actor->SetScale(s);
        for (j=0; j<3; j++)
            v[j] = v[j]/v_mod;
            
        double sina = sqrt(v[0]*v[0] + v[2]*v[2]);
        double cosa = v[1];
        double theta = asin(sina)*(180.0/Pi);
		if (cosa < 0) 
            theta = 180 - theta;
		
        actor->SetPosition(bpos);
        actor->RotateWXYZ(theta,v[2],0,-v[0]);
        ren->AddActor(actor);
        B_Actor_list.append(actor);
	}
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
bool MyVTK::startPlayer(QString posfile, QTimer *theTimer, bool save)
{
	save_image = save;
	LOG_QMSG(posfile);
	timer = theTimer;
	playerData = new QFile(posfile);
	if (!playerData->open(QFile::ReadOnly)) {
		LOG_MSG("Open failure on VTK file");
		return false;
	}
	playerStream = new QTextStream(playerData);
	if (!first_VTK) {
        cleanup();
	}
	playing = true;
	paused = false;

	if (save_image) {
		w2i = vtkWindowToImageFilter::New();
		w2i->SetInput(renWin);	//the render window
//		writer = vtkSmartPointer<vtkPNGWriter>::New();
		writer = vtkSmartPointer<vtkJPEGWriter>::New();
		writer->SetInputConnection(w2i->GetOutputPort()); 
		framenum = 0;
		LOG_MSG("set up writer");
	}
	LOG_MSG("playing");
	return true;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
bool MyVTK::nextFrame()
{
	LOG_MSG("VTK: nextFrame");
	if (!playing)
		return false;
	if (paused)
		return true;
	if (playerStream->atEnd()) {
		LOG_MSG("nextFrame: no more data");
		stop();
		return false;
	}
	TCpos_list.clear();
	DCpos_list.clear();
	bondpos_list.clear();
	int k = 0;
	QString line;
	do {
		line = playerStream->readLine();
		if (line.length() > 0) {
			k++;
			QStringList s = line.split(" ",QString::SkipEmptyParts);
			if (s[0].compare("T") == 0) {
				CELL_POS cp;
				cp.tag = s[1].toInt();
				cp.x = s[2].toInt();
				cp.y = s[3].toInt();
				cp.z = s[4].toInt();
				cp.diameter = s[5].toDouble();
				cp.state = s[6].toDouble();
				TCpos_list.append(cp);
			} else if (s[0].compare("E") == 0) {
				break;
			}
		}
	} while (true);

	bool redo = false;
	if (first_VTK) {
		redo = true;
	}
    renderCells(redo);
    char numstr[6];
    sprintf(numstr,"%05d",framenum);
	if (save_image) {
		w2i->Modified();	//importante 
		writer->SetFileName((casename + numstr + ".jpg").toStdString().c_str()); 
		writer->Write(); 
	}
	framenum++;
	return true;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void MyVTK::saveSnapshot(QString fileName, QString imgType)
{
//    LOG_QMSG("saveSnapshot");
    w2img->SetInput(renWin);
    w2img->Modified();	//important
    if (imgType.compare("png") == 0) {
		vtkSmartPointer<vtkPNGWriter> pngwriter = vtkPNGWriter::New();
		pngwriter->SetInputConnection(w2img->GetOutputPort()); 
		pngwriter->SetFileName((fileName.toStdString()).c_str()); 
		pngwriter->Write();
	} else if (imgType.compare("jpg") == 0) {
		vtkJPEGWriter *jpgwriter = vtkJPEGWriter::New();
		jpgwriter->SetInputConnection(w2img->GetOutputPort()); 
		jpgwriter->SetFileName((fileName.toStdString()).c_str()); 
		jpgwriter->Write();
	} else if (imgType.compare("tif") == 0) {
		vtkTIFFWriter *tifwriter = vtkTIFFWriter::New();
		tifwriter->SetInputConnection(w2img->GetOutputPort()); 
		tifwriter->SetFileName((fileName.toStdString()).c_str()); 
		tifwriter->Write();
	} else if (imgType.compare("bmp") == 0) {
		vtkBMPWriter *bmpwriter = vtkBMPWriter::New();
		bmpwriter->SetInputConnection(w2img->GetOutputPort()); 
		bmpwriter->SetFileName((fileName.toStdString()).c_str()); 
		bmpwriter->Write();
	}
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void MyVTK::record(QString basename, int number)
{
    char numstr[6];
    sprintf(numstr,"%05d",number);
    QString imgType = "png";
    QString fileName = basename + numstr + ".png";
    saveSnapshot(fileName,imgType);
    framenum++;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void MyVTK::pause()
{
	paused = true;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void MyVTK::playon()
{
	paused = false;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void MyVTK::stop()
{
	if (save_image) {
		writer->Delete();
		w2i->Delete();
	}
	delete playerStream;
	playerData->close();
	delete playerData;
	timer->stop();
	playing = false;
	paused = false;
}

