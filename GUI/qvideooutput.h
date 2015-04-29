////////////////////////////////////////////////////////////////////////////////
//! @file   : qvideooutput.h
//! @date   : feb 2013
//!
//! @brief  : Declares a QT style video output class
//!
//! The Basement Lab for Computer Vision and Personal Robotics
//! Copyright (C) 2013 - All Rights Reserved
////////////////////////////////////////////////////////////////////////////////
#ifndef QVIDEOOUTPUT_H
#define QVIDEOOUTPUT_H
////////////////////////////////////////////////////////////////////////////////
// Includes
////////////////////////////////////////////////////////////////////////////////
#include "QObject"
#include "QImage"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
extern "C" {
#include "libavutil/mathematics.h"
#include "libavformat/avformat.h"
#include "libswscale/swscale.h"
}

#include <QPainter>
#include <QWidget>
#include <QFileDialog>
#include <QTemporaryFile>
#include <vtkSmartPointer.h>
#include <vtkRenderWindow.h>
#include <vtkWindowToImageFilter.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <qwt_plot.h>
#include <qwt_painter.h>
#include <qwt_plot_canvas.h>

#define VTK_SOURCE 0
#define QWT_SOURCE 1

////////////////////////////////////////////////////////////////////////////////
//! @class : QVideoOutput
//!
//! @brief : Implements a Qt style wrapper class for some functions
//!          from the FFmpeg library
//!
////////////////////////////////////////////////////////////////////////////////
class QVideoOutput : public QObject
{
public:
   QVideoOutput(QObject * parent=0, int source=0, vtkRenderWindow *VTKrenWin=0, QwtPlot *qwtplot=0);
   virtual ~QVideoOutput();
   bool openMediaFile(int width,
                      int height,
                      const QString & filename);
   bool closeMediaFile();
   void setResolution(int width, int height);
   bool newVtkFrame(vtkImageData * imageData);
   bool isOpen();
   void recorder();
   bool record;
public slots:
   bool newFrame(const QImage & image);
   void startRecorder(QString aviFileName, QString fileFormat, QString codec, int nframes);
   void stopRecorder();
protected:
   // Protected methods ////////////////////////////////////////////////////////
   AVStream * addStream(AVFormatContext *inFormatContext,
                        AVCodec **codec,
                        AVCodecID codecId)const;
   bool openVideo(AVCodec *codec, AVStream *st);
   bool writeVideoFrame(const AVPicture &src,
                        int srcWidth,
                        int srcHeight,
                        AVFormatContext *inFormatContext,
                        AVStream *st);
   void closeVideo(AVStream *st);
   bool flushVideo(AVStream *stream);
   // Protected members ////////////////////////////////////////////////////////
   AVFormatContext * formatContext;
   AVOutputFormat  * outputFormat;
   AVStream        * videoStream;
   AVCodec         * videoCodec;
   SwsContext      * swsContext;
   AVFrame         * frame;
   AVPicture srcPicture;
   AVPicture dstPicture;
   AVPixelFormat streamPixFmt;
   int swsFlags;
   int streamFrameRate;
   int width;
   int height;
   int frameCount;
   int source;  // 0 = VTK, 1 = FACS

//   QWidget *dlgWidget;
   bool openedMediaFile;
   vtkSmartPointer<vtkPNGWriter> pngwriter;
   vtkRenderWindow *renWin;
   vtkSmartPointer<vtkWindowToImageFilter> w2i;
   QwtPlot *qp;

   QTemporaryFile *tempFile;
   int record_nframes, record_it, framenum;
   QString record_basename;
   QString record_fileName;
   QString record_fileFormat;
   QString record_codec;
   char msg[1024];
};
#endif // QVIDEOOUTPUT_H
