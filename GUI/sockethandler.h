#ifndef MYSOCKET_H
#define MYSOCKET_H

#include <QThread>
#include <QTcpServer>

class SocketHandler : public QThread
 {
    Q_OBJECT

public:
    SocketHandler(int newport, QObject *parent = 0);
    ~SocketHandler();
    void run();
	void end();
	bool quitMessage(QString);

	QTcpServer *tcpServer;
	int port;
	bool exiting;
	quint16 blockSize;
	QTcpSocket *socket;
	char msg[2048];
	static const int CPORT0 = 5000;
	static const int CPORT1 = 5001;

private slots:
	 void processor();
signals:
	 void sh_connected();
	 void sh_disconnected();
	 void sh_output(QString);
	 void sh_error(QString);
};

class ExecThread: public QThread
{
	Q_OBJECT
	QString inputFile;
	QString dll_path;
public:
//	ExecThread::ExecThread(QString, QString);
        ExecThread(QString, QString);
        void run();
};


#endif
