#ifndef QMYCHECKBOX_H
#define QMYCHECKBOX_H

#include <qcheckbox.h>
#include <QMouseEvent>
#include "log.h"
LOG_USE();

class QMyCheckBox: public QCheckBox
{
    Q_OBJECT

public:
    QMyCheckBox(QWidget *parent = 0);
	
signals:
    void checkBoxClicked(QString text);

private:

	void mousePressEvent (QMouseEvent *event);

public:
    QString description;

};

#endif

