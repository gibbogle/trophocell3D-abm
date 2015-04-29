/* $Id: clientc.c,v 1.1.1.1 2003/03/06 19:38:16 maine Exp $ */
/* $Name:  $ */

/* clientc.c
 * open a tcp connection for a getData client.
 *
 * The read/write/close routines are used in both clients and servers
 * and are in a separate file.
 *
 * This file has wrappers for the system routines that are
 * not directly or easily callable from Fortran.
 * Specific to system Fortran and C calling conventions.
 *
 * Version for MS Windows MS C and CVF.
 * Differences from nag version are:
 *   Different system include files.
 *   Different procedure name mangling convention in the defines
 *   Specify __stdcall in procedure headers.
 *   Add extern decl for write_error_sub (to get the __stdcall).
 *   Different placement of implicit string length argument.
 *   Have to call Win-specific WSAStartup or nothing works.
 *   Cast *port (perhaps appropriate for nag version also).
 *   Used SOCKET_ERROR and INVALID_SOCKET to test for errors per MS docs,
 *   though the bsd-style test did appear to work.
 *
 * 11 Oct 91, Richard Maine: Version 1.0
 * 17 Jul 01, Richard Maine: MS Windows port.
 */

//#include <sys/types.h>
//#include <winsock.h>
#include <stdio.h>
#include <stdlib.h>
#ifdef _WIN32
   #include <Winsock2.h>
#else
    #include <sys/types.h>
    #include <sys/socket.h>
    #include <unistd.h>
    #include <netinet/in.h>
    #include <netdb.h>

    #define closesocket close
#endif

/* Defines to make routines fortran-callable */
#if defined(__INTEL_COMPILER) || defined(_MSC_VER)
#define tcp_connect TCP_CONNECT
#define tcp_init TCP_INIT
#else
#define tcp_connect tcp_connect_
#define tcp_init tcp_init_
//extern void  write_error_msg (char * msg, int n);
/* Wrapper to call fortran error message routine. */
//#define write_error_msg write_error_sub_
#endif
/*
 * Would probably be better to return an error message to the
 * calling routine instead of calling write_error_msg from in here.
 */

static int started = 0;

/* We want to avoid standard c i/o. */
void printerrormsg(msg)
  char *msg;
{
//    printf("%s\n",msg);
//  write_error_msg(msg,strlen(msg));
}


/* Connect to a specified host and port.
 * Return a socket number in sock. */
/* This version works with gfortran, which places hidden char len arguments at the end */
void tcp_connect(sock, host_name, host_name_len, port, error)
  int *sock, *port, *error;
  char *host_name;
  int *host_name_len;
{
  char cstring[129];
  int i;
  struct hostent *host_ent;
  struct sockaddr_in host;
  int keepalive;
  char *from, *to;

#ifdef _WIN32
  WORD wver;
  WSADATA wsaData;
  /* Start up Windows sockets */
  if (started == 0) {
    started = 1;
    wver = MAKEWORD(2,0);
    i = WSAStartup(wver, &wsaData);
  }
#endif

  *error = 1;
//printf("host_name_len: %d\n",*host_name_len);
//printf("port: %d\n",*port);
  /* Make host_name into a c string */
  for (i=0; (i<*host_name_len) && (i<127) && (host_name[i] != ' '); i++)
    cstring[i] = host_name[i];
  cstring[i] = '\0';
  printerrormsg("In tcp_connect");
  printerrormsg("host_name:");
  printerrormsg(cstring);

  /* Find the host IP address */
  host.sin_family = AF_INET;
  host_ent = gethostbyname(cstring);
  if (host_ent==0) {printerrormsg("Can't find host address."); return;}
  /* Avoid bcopy/memcpy dependence. */
  from = (char *) host_ent->h_addr;
  to = (char *) &host.sin_addr;
  for (i=0; i<host_ent->h_length; i++) *to++ = *from++;

  /* Set the port */
  host.sin_port = htons((unsigned short int) *port);

  /* Create a socket */
  *sock = socket(AF_INET, SOCK_STREAM, IPPROTO_TCP);
//  if (*sock==INVALID_SOCKET) {printerrormsg("Can't create socket.");
  if (*sock < 0) {printerrormsg("Can't create socket.");
    return;
  }

  /* Attempt connection */
//  if (connect(*sock, (struct sockaddr *) &host, sizeof(host)) == SOCKET_ERROR)
  if (connect(*sock, (struct sockaddr *) &host, sizeof(host)) < 0)
    {printerrormsg("Connect failed."); closesocket(*sock); return;}

  /* Enable keepalive packets to check socket connection. */
  keepalive = 1;
  setsockopt(*sock, SOL_SOCKET, SO_KEEPALIVE, (char *) &keepalive, 4);

  *error = 0;
}

void tcp_init(int *res)
{
#ifdef _WIN32
  /* Start up Windows sockets */
    WORD wver;
    WSADATA wsaData;

    *res = 0;
    if (started == 0) {
        wver = MAKEWORD(2,0);
        *res = WSAStartup(wver, &wsaData);
    }
#else
	*res = 0;
#endif
    started = 1;
}

