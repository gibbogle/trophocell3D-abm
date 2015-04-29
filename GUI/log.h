//log.h

#ifndef __LOG_H
#define __LOG_H
//It is bad to include files within other include files, but this is an exception to simplify the code
#include <stdio.h>

//must be put in global scope and used once in all files
#define LOG_DECLARE Log __log;
//used in any file need to access the log files, except the file that has LOG_DECLARE statment
#define LOG_USE(); extern Log __log;
//init log files, called once in program initialization
#define LOG_INIT(fileName) do{if(fileName!= NULL)__log.init(fileName); else printf("error init log file %s", fileName);}while(0);

//here are the different levels of logging
#define LOG_VERBOSE(msg) __log.write("verbose", msg, __FILE__, __LINE__); //detailed info
#define LOG_MSG(msg) __log.write("msg", msg, __FILE__, __LINE__); //brief info
#define LOG_WARN(msg) __log.write("warn", msg, __FILE__, __LINE__); //warning
#define LOG_ERR(msg) __log.write("error", msg, __FILE__, __LINE__); //error
#define LOG_FATAL(msg) __log.write("fatal", msg, __FILE__, __LINE__); //fatal error

#define LOG_QMSG(qmsg) __log.qwrite("qmsg", qmsg, __FILE__, __LINE__); //brief info

class Log
{
FILE *fp;
bool logOk;

public:

Log()
{fp = NULL;}
~Log() //the destructor closes the file
{close();}

void init(char *pfileName)
{ 
	if (pfileName != NULL) {
//Gib		fp = fopen(pfileName, "a+"); 
		fp = fopen(pfileName, "w"); 
		if (fp != NULL) 
			fseek(fp, 0, SEEK_END);
	} 
}

void close()
{if(fp != NULL) fclose(fp); fp = NULL;}

//FIXME: A critical section is required to protect the file writing function in multithreading programms
void write(char *pType, char *pMsg, char *pFileName, int lineNo)
{ if(fp != NULL) {fprintf(fp, "%s\t%s\t%s\t%s\t%d\t%s\n", pType, __DATE__, __TIME__, pFileName, lineNo, pMsg); fflush(fp);} }

void qwrite(char *pType, QString qMsg, char *pFileName, int lineNo)
{ 
	if (fp != NULL) {
		char *pMsg = (qMsg.toAscii()).data();
		fprintf(fp, "%s\t%s\t%s\t%s\t%d\t%s\n", pType, __DATE__, __TIME__, pFileName, lineNo, pMsg); 
		fflush(fp);
	} 
}

};

#endif
