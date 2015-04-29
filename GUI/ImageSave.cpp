#include "ImageSave.h"

#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkSphereSource.h>

#include "vtkSmartPointer.h"

#include "log.h"
#include <QFileDialog>

LOG_USE();


// Constructor
ImageSave::ImageSave(vtkSmartPointer<vtkRenderWindow> renderWin)
{
    renWin = renderWin;
};

void ImageSave::save(QString fname)
{
    QString fileName;
    bool checkDone;
    if (fname.isEmpty()) {
        fileName = QFileDialog::getSaveFileName(this, tr("Select image file"), ".",
            tr("Image files (*.png *.jpg *.tif *.bmp)"));
        if (fileName.compare("") == 0) {
            return;
        }
        checkDone = false;
    } else {
        fileName = fname;
        checkDone = false;
    }
    QFileInfo fi(fileName);
    QString imgType = fi.suffix();
    LOG_QMSG("ImageSave: imgType: " + imgType);
    vtkSmartPointer<vtkWindowToImageFilter> w2img = vtkWindowToImageFilter::New();
    w2img->SetInput(renWin);
    if (imgType.compare("png") == 0) {
        vtkSmartPointer<vtkPNGWriter> pngwriter = vtkPNGWriter::New();
        LOG_QMSG("new w2img");
        pngwriter->SetInputConnection(w2img->GetOutputPort());
        LOG_QMSG("SetInputConnection")
        w2img->Modified();
        LOG_QMSG("Modified");
        pngwriter->SetFileName((fileName.toStdString()).c_str());
        LOG_QMSG("SetFilename");
        pngwriter->Write();
        LOG_QMSG("Write");
    } else if (imgType.compare("jpg") == 0) {
        vtkJPEGWriter *jpgwriter = vtkJPEGWriter::New();
        jpgwriter->SetInputConnection(w2img->GetOutputPort());
        w2img->Modified();
        jpgwriter->SetFileName((fileName.toStdString()).c_str());
        jpgwriter->Write();
    } else if (imgType.compare("tif") == 0) {
        vtkTIFFWriter *tifwriter = vtkTIFFWriter::New();
        tifwriter->SetInputConnection(w2img->GetOutputPort());
        w2img->Modified();
        tifwriter->SetFileName((fileName.toStdString()).c_str());
        tifwriter->Write();
    } else if (imgType.compare("bmp") == 0) {
        vtkBMPWriter *bmpwriter = vtkBMPWriter::New();
        bmpwriter->SetInputConnection(w2img->GetOutputPort());
        w2img->Modified();
        bmpwriter->SetFileName((fileName.toStdString()).c_str());
        bmpwriter->Write();
    }
    if (checkDone) {
        QMessageBox::question(this, tr("Completion"),
                                        tr("Done?"),
                                        QMessageBox::Yes | QMessageBox::No, QMessageBox::No);
    }
}


