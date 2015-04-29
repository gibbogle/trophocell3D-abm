/* $Id: tcpc.c,v 1.1.1.1 2003/03/06 19:38:16 maine Exp $ */
/* $Name:  $ */

/* tcpc.c
 * tcp communication routines for getData.
 *
 * The read/write/close here are used in both clients and servers.
 * These routines must not directly or indirectly reference the
 * fortran write_error_msg routine because some versions of it
 * may call these routines.
 *
 * The routines for opening a tcp connection are different for
 * clients and servers.  Furthermore, the client one may
 * call write_error_msg.  Thus those routines are in separate files.
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
 *   Use closesocket instead of close.
 *   Used SOCKET_ERROR to test for errors per MS docs,
 *     though the bsd-style test did appear to work.
 *
 * 11 Oct 91, Richard Maine: Version 1.0
 * 17 Jul 01, Richard Maine: MS Windows port.
 */


#include <sys/types.h>
#include <stdio.h>
#ifdef _WIN32		// was __WIN32
#include <winsock.h>
#define sleep(n) Sleep(1000 * n)
#else
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <arpa/inet.h>
#include <netinet/in.h>
#include <sys/socket.h>
#include <netdb.h>

#define closesocket close
#define INVALID_SOCKET 1
#define SOCKET_ERROR 2
#endif

/* Defines to make routines fortran-callable */
#if defined(__INTEL_COMPILER) || defined(_MSC_VER)
#define tcp_close TCP_CLOSE
#define tcp_read TCP_READ
#define tcp_write TCP_WRITE
#define sleeper SLEEPER	
#else
#define tcp_close tcp_close_
#define tcp_read tcp_read_
#define tcp_write tcp_write_
#define sleeper sleeper_
#endif



void tcp_close(sock)
  int *sock;
{
  shutdown(*sock,2); /* clears any pending i/o. */
  closesocket(*sock);
}

/* Write to a tcp socket. */
void tcp_write(sock, buf, buflen, error)
  int *sock, *buflen, *error;
  char *buf;
{
  int i;

  i = send(*sock,buf,*buflen,0);
//  *error = (i == SOCKET_ERROR);
  if (i == *buflen) {
    *error = 0;
  } else {
    *error = i;
  }
}

/* Read from a tcp socket. */
void tcp_read(sock, buf, buflen, error)
  int *sock, *buflen, *error;
  char *buf;
{
  int nleft, nread;
  char *bufptr;

  bufptr = buf;
  nleft = *buflen;
  *error = 1;

//while (nleft > 0)
//    {
      nread = recv(*sock, bufptr, nleft, 0);
      if (nread == SOCKET_ERROR) return;
//     nleft -= nread;
//     bufptr += nread;
//    }

  *error = 0;
}

void sleeper(int *nsecs)
{
    unsigned int nt;
    nt = *nsecs;
    sleep(nt);
}
