#include "plot.h"
#include "log.h"

LOG_USE();

#define USE_LEGEND false

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
QString getImageFile()
{
	QString fileName = QFileDialog::getSaveFileName(0,"Select image file", ".", 
		"Image files (*.png *.jpg *.tif *.bmp)");
	return fileName;
}


//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
Plot::Plot(QString aname, QString acasename, QWidget *parent)
	: QwtPlot(parent)
{
	name = aname;
	for (int i=0; i<ncmax; i++) {
		curve[i] = 0;
	}
	ncurves = 0;
	yscale = 0;
//	setAxisTitle(QwtPlot::xBottom, "Time (hours)");
	if (name.compare("") != 0) {
		curve[0] = new QwtPlotCurve(acasename);
		curve[0]->attach(this);
		ncurves = 1;
		replot();
	}
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
Plot::~Plot()
{
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void Plot::mousePressEvent (QMouseEvent *event) {
	event->accept();
	if (event->button() == Qt::RightButton) {
		int w = this->width();
		int h = this->height();
		QPixmap pixmap(w, h);
        pixmap.fill(Qt::white);

		QwtPlotPrintFilter filter;
		int options = QwtPlotPrintFilter::PrintAll;
		options &= ~QwtPlotPrintFilter::PrintBackground;
		options |= QwtPlotPrintFilter::PrintFrameWithScales;
		filter.setOptions(options);

		this->print(pixmap, filter);

		QString fileName = getImageFile();
		if (fileName.isEmpty()) {
			return;
		}
		pixmap.save(fileName,0,-1);

        // Testing conversion to QImage, for qvideooutput
        QImage image = pixmap.toImage();
	}
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void Plot::addCurve(QString name)
{
	for (int k=0; k<ncmax; k++) {
		if (curve[k] == 0) {
			curve[k] = new QwtPlotCurve(name);
			curve[k]->attach(this);
			ncurves++;
			replot();
			break;
		}
	}
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void Plot::removeCurve(QString name)
{
	for (int i=0; i<ncmax; i++) {
		if (curve[i] != 0) {
			if (name.compare(curve[i]->title().text()) == 0) {
				curve[i]->detach();
				delete curve[i];
				curve[i] = 0;
				ncurves--;
			}
		}
	}
	replot();
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void Plot::removeAllCurves()
{
	for (int i=0; i<ncmax; i++) {
		if (curve[i] != 0) {
			curve[i]->detach();
			delete curve[i];
			curve[i] = 0;
			ncurves--;
		}
	}
	replot();
}
	
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void Plot::setYScale(double maxval)
{
    yscale = calc_yscale_ts(maxval);
	setAxisScale(QwtPlot::yLeft, 0, yscale, 0);
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void Plot::redraw(double *x, double *y, int n, QString name, QString tag, double fixed_yscale, bool profile)
{
    QwtLegend *legend;
//    if (n == 1) { // That is, this is the first plotting instance.
//        if (fixed_yscale == 0) {
//            yscale = max(yscale,calc_yscale_ts(y[0]));
//        } else {
//            yscale = fixed_yscale;
//        }
//        setAxisScale(QwtPlot::yLeft, 0, yscale, 0);
        if (USE_LEGEND){
            legend = this->legend();
            if (legend == NULL) {
                legend = new QwtLegend();
                this->insertLegend(legend, QwtPlot::RightLegend);
            }
        }
//    }
    // Note: Number of pen colors should match ncmax
    QColor pencolor[] = {Qt::black, Qt::red, Qt::blue, Qt::darkGreen, Qt::magenta, Qt::darkCyan };
    QPen *pen = new QPen();
    for (int k=0; k<ncmax; k++) {
        if (curve[k] == 0) continue;
//        if (profile) LOG_QMSG("profile redraw "+name);
        if (name.compare(curve[k]->title().text()) == 0) {
//            LOG_QMSG("redraw " + tag);
            // Just in case someone set ncmax > # of pen colors (currently = 6)
            if (k < 6) {
                pen->setColor(pencolor[k]);
            } else {
                pen->setColor(pencolor[0]);
            }
            curve[k]->setPen(*pen);
            curve[k]->setData(x, y, n);
//            legend = new QwtLegend();
//            this->insertLegend(legend, QwtPlot::RightLegend);
//            LOG_MSG("did insertLegend");
//            if (!profile) {
                if (fixed_yscale == 0) {
                    double ylast = y[n-1];
                    if (ylast > yscale) {
                        yscale = max(yscale,calc_yscale_ts(ylast));
//                        sprintf(msg,"ylast: %f yscale: %f",ylast,yscale);
//                        LOG_MSG(msg);
                    }
                } else {
                    yscale = fixed_yscale;
                }
                setAxisScale(QwtPlot::yLeft, 0, yscale, 0);
//            }
            replot();
        }
    }
    delete pen;
}
//-----------------------------------------------------------------------------------------
// This is to plot the total number of cognate cells (y2 = ytotal), and the number of 
// "seed" cells (y1 = yseed)
//-----------------------------------------------------------------------------------------
void Plot::redraw2(double *x1, double *y1, double *x2, double *y2, int n1, int n2)
{
	LOG_MSG("redraw2");
	if (n1 == 1) { // That is, this is the first plotting instance.
        yscale = max(yscale,calc_yscale_ts(y1[0]));
        yscale = max(yscale,calc_yscale_ts(y2[0]));
		setAxisScale(QwtPlot::yLeft, 0, yscale, 0);
	}
    curve[0]->setData(x1, y1, n1);
    curve[1]->setData(x2, y2, n2);
	QPen *pen = new QPen();
	pen->setColor(Qt::red);
	curve[1]->setPen(*pen);
	delete pen;
	this->insertLegend(&QwtLegend(), QwtPlot::RightLegend);
	double ylast = y1[n1-1];
	double ylast2 = y2[n2-1];
	if (ylast2 > ylast) {
		ylast = ylast2;
	}
	if (ylast > yscale) {
        yscale = calc_yscale_ts(ylast);
		setAxisScale(QwtPlot::yLeft, 0, yscale, 0);
	}
    replot();
}

//-----------------------------------------------------------------------------------------
// This is for a one-off plot of y1 vs x1 and y2 vs x2, where there are n1 data points
// for y1 and n2 for y2, and n1 and n2 may differ, as may x1 and x2.
// curve is the current simulation output, curve2 is the previous run
// Note: assumes that data are hourly.
//-----------------------------------------------------------------------------------------
void Plot::draw2(double *x1, double *y1, double *x2, double *y2, int n1, int n2)
{
    curve[0]->setData(x1, y1, n1);
    curve[1]->setData(x2, y2, n2);
	QPen *pen = new QPen();
	pen->setColor(Qt::red);
	curve[1]->setPen(*pen);
	delete pen;
	this->insertLegend(&QwtLegend(), QwtPlot::RightLegend);
	double ylast = 0;
	for (int i=0; i<n1; i++)
		ylast = max(ylast,y1[i]);
	for (int i=0; i<n2; i++)
		ylast = max(ylast,y2[i]);
    yscale = calc_yscale_ts(ylast);
	setAxisScale(QwtPlot::yLeft, 0, yscale, 0);
    replot();
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
double Plot::calc_yscale_ts(double yval)
{
    int v;
    double yscale;
    if (yval > 10) {
        v = int(1.3*yval);
        yscale = double(v);
    } else if (yval > 1) {
        v = int(13*yval);
        yscale = v/10;
    } else if (yval > 0.1) {
        v = int(130*yval);
        yscale = v/100.;
    } else {
        yscale = 0.1;
    }
//    sprintf(msg,"calc_yscale_ts: %f %f",yval,yscale);
//    LOG_MSG(msg);
    return yscale;
}
