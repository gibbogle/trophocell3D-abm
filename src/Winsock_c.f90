module winsock
implicit none

TYPE winsockport
    LOGICAL           :: is_open
    INTEGER           :: handle
    character*(64)    :: host_name
    integer           :: ip_addr    ! server IP address u_long: wp%ip_addr.S_addr = inet_addr(ip_address)
    INTEGER           :: ip_port    ! IP port on server
    INTEGER           :: protocol   ! TCP or UDP for ethernet
END TYPE winsockport

integer :: sock = 0
integer, parameter :: nflog_ws=21
integer, parameter :: IPPROTO_TCP = 1

contains

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
logical FUNCTION Winsock_init (startstop)
IMPLICIT NONE
INTEGER, INTENT(IN)     :: startstop
INTEGER                 :: res

res = 0
!   Initialize Winsock
IF (startstop == 1) THEN
    write(nflog_ws,*) 'do tcp_init'
    call tcp_init(res)
    write(nflog_ws,*) 'did tcp_init: res: ',res
ELSE
    call tcp_close(sock)
END IF
Winsock_init = (res == 0)
END FUNCTION Winsock_init

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
SUBROUTINE Set_Winsock_Port (wp,error)
IMPLICIT NONE
TYPE(winsockport), INTENT(INOUT) :: wp
integer :: error, host_name_len

error = 0
IF (wp%is_open) then
!    call tcp_close (wp%handle)
!    wp%is_open = .false.
    write(nflog_ws,*) 'port is already open'
endif
host_name_len = len_trim(wp%host_name)
write(nflog_ws,*) 'call tcp_connect: ',wp%ip_port,' ',wp%host_name
call tcp_connect(sock, wp%host_name, host_name_len, wp%ip_port, error)
write(nflog_ws,*) 'port is opened: ',wp%ip_port
wp%handle = sock
wp%is_open = .true.
END SUBROUTINE

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
SUBROUTINE Winsock_send (wp, outbuf, ncout,error)
IMPLICIT NONE
TYPE(winsockport), INTENT(IN)   :: wp
CHARACTER(LEN=*), INTENT(IN)    :: outbuf
INTEGER, INTENT(IN)             :: ncout
INTEGER	                        :: error
!CHARACTER(LEN=128)              :: msg

error = 0
!msg = outbuf(1:ncout)
IF (.NOT.wp%is_open) then
    write(nflog_ws,*) 'Winsock_send: ',outbuf(1:ncout)	!msg
    write(nflog_ws,*) 'port not open: ',wp%ip_port
    error = 1
    RETURN
endif

!call tcp_write(wp%handle,msg,ncout,error)
call tcp_write(wp%handle,outbuf(1:ncout),ncout,error)

end subroutine

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
SUBROUTINE winsock_close (wp)
IMPLICIT NONE
TYPE(winsockport), INTENT(INOUT)    :: wp
call tcp_close (wp%handle)
wp%is_open = .false.
end subroutine

end module
