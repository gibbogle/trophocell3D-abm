////////////////////////////////////////////////////////////////////////////////
//! @file   : qvideooutput.cpp
//! @date   : feb 2013
//!
//! @brief  : Defines a QT style video output class
//!
//! The Basement Lab for Computer Vision and Personal Robotics
//! Copyright (C) 2013 - All Rights Reserved
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Includes
////////////////////////////////////////////////////////////////////////////////
#include "qvideooutput.h"
#include "stdio.h"
#include "log.h"

LOG_USE();

////////////////////////////////////////////////////////////////////////////////
// Macro redefinition since original does not compile in c++
///////////////////////////////////////////////////////////////////////////////
#undef av_err2str
#define av_err2str(errnum) \
        av_make_error_string(reinterpret_cast<char*>(_alloca(AV_ERROR_MAX_STRING_SIZE)),\
                             AV_ERROR_MAX_STRING_SIZE, errnum)
////////////////////////////////////////////////////////////////////////////////
//  QVideoOutput::QVideoOutput
//!
//! @brief Constructor
//!
//! @param[in] parent : A parent object.
//!
////////////////////////////////////////////////////////////////////////////////
QVideoOutput::QVideoOutput(QObject *parent, int imageSource, vtkRenderWindow *VTKrenWin, QwtPlot *qwtplot)
: QObject(parent)
, swsContext(0x0)
, formatContext(0x0)
, outputFormat(0x0)
, videoStream(0x0)
, videoCodec(0x0)
, frame(0x0)
, swsFlags(SWS_BICUBIC)
, streamPixFmt(AV_PIX_FMT_YUV420P) // default pix_fmt
, streamFrameRate(25)              // 25 images/s
, width(640)
, height(480)
, w2i(0x0)
, openedMediaFile(false)
{
   source = imageSource;
   if (source == VTK_SOURCE) {
       // Set renWin
       renWin = VTKrenWin;
       LOG_MSG("created videoVTK");
   } else if (source == QWT_SOURCE) {
       qp = qwtplot;
       LOG_MSG("created videoFACS");
   }
   record_it = 0;
   record = false;
   // Init FFmpeg
   av_register_all();
}
////////////////////////////////////////////////////////////////////////////////
//  QVideoOutput::~QVideoOutput
//!
//! @brief Destructor
//!
////////////////////////////////////////////////////////////////////////////////
QVideoOutput::~QVideoOutput()
{
}
////////////////////////////////////////////////////////////////////////////////
//  QVideoOutput::setResolution
//!
//! @brief Sets resolution.
//!
//! @param[in] inWidth  : Specifies width. Must be a multiple of two
//! @param[in] inHeight : Specifies height. Must be a multiple of two
//!
////////////////////////////////////////////////////////////////////////////////
void QVideoOutput::setResolution(int inWidth, int inHeight)
{
   Q_ASSERT(inWidth%2  == 0);
   Q_ASSERT(inHeight%2 == 0);
   width  = inWidth;
   height = inHeight;
}
////////////////////////////////////////////////////////////////////////////////
//  QVideoOutput::addStream
//!
//! @brief Adds stream
//!
//! @param[in]  inFormatContext
//! @param[out] codec
//! @param[in]  codecId
//!
////////////////////////////////////////////////////////////////////////////////
AVStream *QVideoOutput::addStream(AVFormatContext * inFormatContext,
                                  AVCodec **codec,
                                  AVCodecID codecId)const
{
    AVCodecContext *c;
    AVStream *st;
    // find the encoder
    *codec = avcodec_find_encoder(codecId);
    if (!(*codec))
    {
        fprintf(stderr, "Could not find encoder for '%s'\n",
                avcodec_get_name(codecId));
        return 0x0;
    }
    st = avformat_new_stream(inFormatContext, *codec);
    if (!st)
    {
        fprintf(stderr, "Could not allocate stream\n");
        return 0x0;
    }
    st->id = inFormatContext->nb_streams-1;
    c = st->codec;
    switch ((*codec)->type)
    {
    case AVMEDIA_TYPE_AUDIO:
        st->id = 1;
        c->sample_fmt  = AV_SAMPLE_FMT_S16;
        c->bit_rate    = 64000;
        c->sample_rate = 44100;
        c->channels    = 2;
        break;
    case AVMEDIA_TYPE_VIDEO:
        avcodec_get_context_defaults3(c, *codec);
        c->codec_id = codecId;
        c->bit_rate = 400000;
        c->width    = width;
        c->height   = height;
        // timebase: This is the fundamental unit of time (in seconds) in terms
        // of which frame timestamps are represented. For fixed-fps content,
        // timebase should be 1/framerate and timestamp increments should be
        // identical to 1.
        c->time_base.den = streamFrameRate;
        c->time_base.num = 1;
        c->gop_size      = 12; // emit one intra frame every twelve frames at most
        c->pix_fmt       = streamPixFmt;
        if (c->codec_id == AV_CODEC_ID_MPEG2VIDEO)
        {
            // just for testing, we also add B frames
            c->max_b_frames = 2;
        }
        if (c->codec_id == AV_CODEC_ID_MPEG1VIDEO)
        {
            // Needed to avoid using macroblocks in which some coeffs overflow.
            // This does not happen with normal video, it just happens here as
            // the motion of the chroma plane does not match the luma plane.
            c->mb_decision = 2;
        }
    break;
    default:
        break;
    }
    // Some formats want stream headers to be separate.
    if (inFormatContext->oformat->flags & AVFMT_GLOBALHEADER)
        c->flags |= CODEC_FLAG_GLOBAL_HEADER;
    return st;
}
////////////////////////////////////////////////////////////////////////////////
//  QVideoOutput::openVideo
//!
//! @brief Opens video
//!
//! @param[in]  codec
//! @param[in]  stream
//!
////////////////////////////////////////////////////////////////////////////////
bool QVideoOutput::openVideo(AVCodec *codec, AVStream *stream)
{
    int ret;
    AVCodecContext *c = stream->codec;
    // open the codec
    ret = avcodec_open2(c, codec, NULL);
    if (ret < 0)
    {
        sprintf(msg, "Could not open video codec: %s\n", av_err2str(ret));
        LOG_MSG(msg);
       fprintf(stderr, "Could not open video codec: %s\n", av_err2str(ret));
       return false;
    }
    // allocate and init a re-usable frame
    frame = avcodec_alloc_frame();
    if (!frame)
    {
        sprintf(msg, "Could not allocate video frame\n");
        LOG_MSG(msg);
       fprintf(stderr, "Could not allocate video frame\n");
       return false;
    }
    // Allocate the encoded raw picture.
    ret = avpicture_alloc(&dstPicture, c->pix_fmt, c->width, c->height);
    if (ret < 0)
    {
        sprintf(msg, "Could not allocate picture: %s\n", av_err2str(ret));
        LOG_MSG(msg);
        fprintf(stderr, "Could not allocate picture: %s\n", av_err2str(ret));
        return false;
    }
    // copy data and linesize picture pointers to frame
    *((AVPicture *)frame) = dstPicture;
    return true;
}
////////////////////////////////////////////////////////////////////////////////
//  QVideoOutput::writeVideoFrame
//!
//! @brief Writes video frame
//!
//! @param[in]  src
//! @param[in]  srcWidth
//! @param[in]  srcHeight
//! @param[in]  inFormatContext
//! @param[in]  stream
//!
////////////////////////////////////////////////////////////////////////////////
bool QVideoOutput::writeVideoFrame(const AVPicture &src,
                                   int srcWidth,
                                   int srcHeight,
                                   AVFormatContext * inFormatContext,
                                   AVStream *stream)
{
    int ret;

    AVCodecContext *c = stream->codec;
    if (c->pix_fmt != AV_PIX_FMT_RGBA)
    {
       // as we only use RGBA picture, we must convert it
       // to the codec pixel format if needed
      if (!swsContext)
      {
           swsContext = sws_getContext(srcWidth,
                                       srcHeight,
                                       AV_PIX_FMT_BGRA,
                                       c->width,
                                       c->height,
                                       c->pix_fmt,
                                       swsFlags,
                                       NULL,
                                       NULL,
                                       NULL);
           if (!swsContext)
           {
               fprintf(stderr, "Could not initialize the conversion context\n");
               return false;
           }
       }
       sws_scale(swsContext,
                 (const uint8_t * const *)src.data,
                 src.linesize,
                 0,
                 c->height,
                 dstPicture.data,
                 dstPicture.linesize);
    }
    if (inFormatContext->oformat->flags & AVFMT_RAWPICTURE)
    {
        // Raw video case - directly store the picture in the packet
        AVPacket pkt;
        av_init_packet(&pkt);
        pkt.flags        |= AV_PKT_FLAG_KEY;
        pkt.stream_index  = stream->index;
        pkt.data          = dstPicture.data[0];
        pkt.size          = sizeof(AVPicture);
        ret = av_interleaved_write_frame(inFormatContext, &pkt);
    }
    else
    {
        // encode the image
        AVPacket pkt;
        int gotOutput;
        av_init_packet(&pkt);
        pkt.data = NULL;    // packet data will be allocated by the encoder
        pkt.size = 0;
        ret = avcodec_encode_video2(c, &pkt, frame, &gotOutput);
        if (ret < 0)
        {
            LOG_QMSG("Error encoding video frame");
            fprintf(stderr, "Error encoding video frame: %s\n",
                    av_err2str(ret));
            return false;
        }
        // If size is zero, it means the image was buffered.
        if (gotOutput)
        {
            if (c->coded_frame->key_frame)
                pkt.flags |= AV_PKT_FLAG_KEY;
            pkt.stream_index = stream->index;
            // Write the compressed frame to the media file.
//            LOG_QMSG("av_interleaved_write_frame");
            ret = av_interleaved_write_frame(inFormatContext, &pkt);
        }
        else
        {
            ret = 0;
        }
    }
    if (ret != 0)
    {
        fprintf(stderr, "Error while writing video frame: %s\n",
                av_err2str(ret));
        return false;
    }
    frameCount++;
    return true;
}

////////////////////////////////////////////////////////////////////////////////
//  QVideoOutput::openMediaFile
//!
//! @brief Opens media file
//!
//! @param[in]  width    :
//! @param[in]  height   :
//! @param[in]  filename :
//!
////////////////////////////////////////////////////////////////////////////////
bool QVideoOutput::openMediaFile(int imwidth, int imheight, const QString & filename)
{
    openedMediaFile = false;
    width = imwidth;
    height = imheight;
    sprintf(msg,"width: %d  height: %d",width,height);
    LOG_MSG(msg);
    // allocate the output media context
    if (record_codec.contains("h264")) {
        avformat_alloc_output_context2(&formatContext, NULL, "h264", filename.toAscii().data());
    } else if (record_codec.contains("mpeg4")) {
        avformat_alloc_output_context2(&formatContext, NULL, "mpeg4", filename.toAscii().data());
    } else {
        avformat_alloc_output_context2(&formatContext, NULL, NULL, filename.toAscii().data());     // => default = mpeg4
    }
   if (!formatContext)
   {
       printf("Could not deduce output format from file extension: using MPEG.\n");
       LOG_QMSG("Could not deduce output format from file extension: using MPEG.");
       avformat_alloc_output_context2(&formatContext, NULL, "mpeg", filename.toAscii().data());
   }
   if (!formatContext)
   {
       LOG_QMSG("No formatContext found!")
       return false;
   }
   outputFormat = formatContext->oformat;
   // Add the video streams using the default format codecs
   // and initialize the codecs.
   sprintf(msg,"outputFormat: %s %d",outputFormat->name,outputFormat->video_codec);
   LOG_MSG(msg);
   videoStream = NULL;
   if (outputFormat->video_codec != AV_CODEC_ID_NONE)
   {
       videoStream = addStream(formatContext, &videoCodec, outputFormat->video_codec);
   }
   // Now that all the parameters are set, we can open the audio and
   // video codecs and allocate the necessary encode buffers.
   if (videoStream) {
       LOG_QMSG("addStream successful");
       bool success = openVideo(videoCodec, videoStream);
       if (!success) {
           LOG_QMSG("openVideo failed");
           return false;
       }
   } else {
       LOG_QMSG("addStream failed");
       return false;
   }
   av_dump_format(formatContext, 0, filename.toAscii().data(), 1);
   int ret = 0;
   // open the output file, if needed
   if (!(outputFormat->flags & AVFMT_NOFILE))
   {
       ret = avio_open(&formatContext->pb, filename.toAscii().data(), AVIO_FLAG_WRITE);
       if (ret < 0)
       {
           fprintf(stderr, "Could not open '%s': %s\n", filename, av_err2str(ret));
           return false;
       }
       LOG_QMSG("Opened output file");
   }
   // Write the stream header, if any.
   ret = avformat_write_header(formatContext, NULL);
   if (ret < 0)
   {
       fprintf(stderr, "Error occurred when opening output file: %s\n", av_err2str(ret));
       return false;
   }
   if (frame)
       frame->pts = 0;
   bool success = (avpicture_alloc(&srcPicture, AV_PIX_FMT_RGBA, imwidth, imheight) >= 0);
   sprintf(msg,"openMediaFile: width,height: %d %d",imwidth,imheight);
   LOG_MSG(msg);
   openedMediaFile = success;
   return success;
}

////////////////////////////////////////////////////////////////////////////////
//  QVideoOutput::newFrame
//!
//! @brief Adds new frame to output stream
//!
//! @param[in]  image :
//!

/*
// To create a QImage of the QwtPlot:

virtual void YourPlot::drawCanvas( QPainter *painter )
{
    QImage image( canvas()->size(), QImage::QImage::Format_RGB32 );
    image.fill( QColor( Qt::white ).rgb() ); // guess you don't need this line

    QPainter p( &image );
    QwtPlot::drawCanvas( &p );
    painter->drawImage( 0, 0, image );
}

// but ... look at plot::mousePressEvent(), which uses QPixmap:
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

// In QPixmap class there is a method QImage QPixmap::toImage()
// to convert a pixmap to image

// For individual points (flow cytometry data, CFSE):
    QwtPlotMarker* m = new QwtPlotMarker();
    m->setSymbol( QwtSymbol( QwtSymbol::Diamond, Qt::red, Qt::NoPen, QSize( 10, 10 ) ) );
    m->setValue( QPointF( 1.5, 2.2 ) );
    m->attach( plot );
*/

////////////////////////////////////////////////////////////////////////////////
bool QVideoOutput::newFrame(const QImage & image)
{
//    LOG_QMSG("newFrame");
   const int width  = image.width();
   const int height = image.height();
   // write video frame
   for (int y = 0; y < height; y++)
   {
      const uint8_t * scanline = image.scanLine(y);
      for (int x = 0; x < width*4; x++)
      {
         srcPicture.data[0][y * srcPicture.linesize[0] + x] = scanline[x];
      }
   }
   bool success = writeVideoFrame(srcPicture, width, height, formatContext, videoStream);
   if (!success) {
       LOG_QMSG("ERROR: writeVideoFrame failed");
       return false;
   }
   frame->pts += av_rescale_q(1,
                              videoStream->codec->time_base,
                              videoStream->time_base);
   return true;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
bool QVideoOutput::newVtkFrame(vtkImageData * imageData)
{
//    LOG_QMSG("newVtkFrame");
    int width = imageData->GetDimensions()[0];
    int height = imageData->GetDimensions()[1];
    QImage image(width, height, QImage::Format_RGB32);
//    sprintf(msg,"newVtkFrame: width,height: %d %d image.width,height: %d %d",width,height,image.width(),image.height());
//    LOG_MSG(msg);
    QRgb* rgbPtr = reinterpret_cast<QRgb*>(image.bits()) + width * (height-1);
    unsigned char* colorsPtr = reinterpret_cast<unsigned char *>(imageData->GetScalarPointer());
     // mirror vertically
    for (int row = 0; row < height; ++row) {
       for (int col = 0; col < width; ++col) {
         // Swap rgb
         *(rgbPtr++) = QColor(colorsPtr[0], colorsPtr[1], colorsPtr[2]).rgb();
         colorsPtr += 3;
       }
       rgbPtr -= width * 2;
     }
    bool success = newFrame(image);
    if (!success) {
        LOG_QMSG("ERROR: newFrame failed");
        return false;
    }
    return true;
}

////////////////////////////////////////////////////////////////////////////////
//  QVideoOutput::closeMediaFile
//!
//! @brief Closes media file
//!
////////////////////////////////////////////////////////////////////////////////
bool QVideoOutput::closeMediaFile()
{
    if (videoStream)
        flushVideo(videoStream);
    av_free(srcPicture.data[0]);
   // Write the trailer, if any. The trailer must be written before you
   // close the CodecContexts open when you wrote the header; otherwise
   // av_write_trailer() may try to use memory that was freed on
   // av_codec_close().
   av_write_trailer(formatContext);
   // Close each codec.
   if (videoStream)
       closeVideo(videoStream);
   if (swsContext)
   {
      sws_freeContext(swsContext);
      swsContext = 0x0;
   }
   // Free the streams.
   for (unsigned int i = 0; i < formatContext->nb_streams; i++)
   {
       av_freep(&formatContext->streams[i]->codec);
       av_freep(&formatContext->streams[i]);
   }
   if (!(outputFormat->flags & AVFMT_NOFILE))
   {
      // Close the output file.
      avio_close(formatContext->pb);
   }
   // free the stream
   av_free(formatContext);
   openedMediaFile = false;
   return true;
}

////////////////////////////////////////////////////////////////////////////////
//  QVideoOutput::closeVideo
//!
//! @brief Closes video
//!
//! @param[in]  stream :
//!
////////////////////////////////////////////////////////////////////////////////
void QVideoOutput::closeVideo(AVStream *stream)
{
    avcodec_close(stream->codec);
    av_free(dstPicture.data[0]);
    av_free(frame);
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
bool QVideoOutput::isOpen()
{
    return openedMediaFile;
}

bool QVideoOutput::flushVideo(AVStream *stream)
{
    int ret, gotOutput, i=0;
    AVCodecContext *c = stream->codec;
    for (gotOutput = 1; gotOutput; i++) {
        AVPacket pkt;
        av_init_packet(&pkt);
        pkt.data = NULL;    // packet data will be allocated by the encoder
        pkt.size = 0;
        ret = avcodec_encode_video2(c, &pkt, NULL, &gotOutput);
        if (ret < 0) {
            LOG_QMSG("Error encoding video frame");
            return false;
        }
        // If size is zero, it means the image was buffered.
        if (gotOutput) {
            if (c->coded_frame->key_frame)
                pkt.flags |= AV_PKT_FLAG_KEY;
            pkt.stream_index = stream->index;
            // Write the compressed frame to the media file.
            ret = av_interleaved_write_frame(formatContext, &pkt);
        }
        else
        {
            ret = 0;
        }
        if (ret != 0)
        {
            LOG_QMSG("Error flushing video frame");
            return false;
        }
    }
    return true;
}

//-----------------------------------------------------------------------------------------
// Uses: w2i, renWin, pngwriter,
//       record, record_basename, record_nframes, record_it, framenum
//-----------------------------------------------------------------------------------------
void QVideoOutput::startRecorder(QString videoFileName, QString fileFormat, QString codec, int nframes)
{
    if (source == VTK_SOURCE) {
        if (w2i == 0) {
            w2i = vtkWindowToImageFilter::New();
            w2i->SetInput(renWin);	//the render window
        }
    }
    record = true;
    record_fileName = videoFileName;
    record_fileFormat = fileFormat;
    record_codec = codec;
    record_nframes = nframes;
    record_it = 0;

    LOG_MSG("Started recording");
}

//-----------------------------------------------------------------------------------------
// Uses: w2i, videoOutput, pngwriter,
//       record, record_basename, record_nframes, record_it, framenum
//-----------------------------------------------------------------------------------------
void QVideoOutput::recorder()
{
    int imwidth, imheight;
    vtkImageData *id;
    QImage im;

    sprintf(msg,"recorder: record_it: %d",record_it);
    LOG_MSG(msg);
    if (!record) return;
    if (source == VTK_SOURCE) {
        id = vtkImageData::New();
        id = w2i->GetOutput();
        w2i->Modified();	//important
        id->Update();
        imwidth = id->GetDimensions()[0];
        imheight = id->GetDimensions()[1];
        if (imwidth == 0) {
            LOG_QMSG("ERROR: recorder: vtkImageData dimension = 0");
            exit(1);
        }
    } else if (source == QWT_SOURCE) {
        // Create an image
//        QImage image( qp->canvas()->size(), QImage::Format_RGB32 );
        QImage image( qp->size(), QImage::Format_RGB32 );
        image.fill( QColor( Qt::white ).rgb() ); // guess you don't need this line
//        QPainter p( &image );
//        qp->drawCanvas( &p );
//        p.drawImage( 0, 0, image );
        qp->print(image);
        im = image;
        imwidth = im.width();
        imheight = im.height();
    }
    record_it++;
    if (!isOpen()) {
        // Generate temporary filename
        tempFile = new QTemporaryFile("qt_temp.XXXXXX.avi");
        if (tempFile->open())
        {
           // Open media file and prepare for recording
           QString fileName = tempFile->fileName();
           bool recording = openMediaFile(imwidth, imheight, fileName.toAscii().data());
            if (!recording) {
                LOG_QMSG("ERROR: openMediaFile failed");
                record = false;
                return;
            }
        }
    }
    bool success;
    if (source == VTK_SOURCE) {
        success = newVtkFrame(id);
        if (!success) {
            LOG_QMSG("ERROR: newVtkFrame failed");
            record = false;
            exit(1);
        }
    } else if (source == QWT_SOURCE) {
        success = newFrame(im);
        if (!success) {
            LOG_QMSG("ERROR: newFrame failed for QWT_SOURCE");
            record = false;
            exit(1);
        }
    }
    if (record_it == record_nframes) {
        stopRecorder();
        return;
    }
}

//-----------------------------------------------------------------------------------------
// Uses: record
//-----------------------------------------------------------------------------------------
void QVideoOutput::stopRecorder()
{
    LOG_QMSG("stopRecorder");
    if (!record) return;
    record = false;
    closeMediaFile();
    LOG_QMSG("Closed media file");
    QString fileName = record_fileName;
    LOG_QMSG(fileName);
//  if (record_fileName.contains(".mov") && record_codec.contains("h264")) {
    if (record_fileName.contains(".mov") || record_codec.contains("h264")) {
        char cmd[1024];
        sprintf(cmd,"ffmpeg -i %s -vcodec h264 -y %s",tempFile->fileName().toStdString().c_str(), fileName.toStdString().c_str());
        LOG_MSG(cmd);
        int res = system(cmd);
        sprintf(msg,"Result code: %d",res);
        LOG_MSG(msg);
    } else if (fileName.isNull() == false) {
       QFile::copy(tempFile->fileName(), fileName);
    }
    delete tempFile;
    tempFile = 0x0;
    LOG_MSG("Stopped recording");
}

