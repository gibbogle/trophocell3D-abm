!-----------------------------------------------------------------------------------------
! Runs tropho-abm
!-----------------------------------------------------------------------------------------
PROGRAM tropho_main
use tropho_mod
integer :: ncpu, res, summarydata(100)
character*(128) :: infile,outfile
character*(64) :: travelfile = 'travel_time_dist.out'
integer :: status, nlen, cnt, i, inbuflen, outbuflen
integer :: jstep, hour, Nsteps_tmp
character*(128) :: b, c, progname

call process_command_line(ncpu,infile,outfile)

outfile = 'tropho_main.out'
inbuflen = len(infile)
outbuflen = len(outfile)

call get_command (b, nlen, status)
if (status .ne. 0) then
    write (*,*) 'get_command failed with status = ', status
    stop
end if
!write (*,*) 'command line = ', b(1:len)
call get_command_argument (0, c, nlen, status)
if (status .ne. 0) then
    write (*,*) 'Getting command name failed with status = ', status
    stop
end if
!write (*,*) 'command name = ', c(1:len)
progname = c(1:nlen)
cnt = command_argument_count ()
!write (*,*) 'number of command arguments = ', cnt
if (cnt < 2) then
    write(*,*) 'Use: ',trim(progname),' num_cpu input_file'
    stop
endif

do i = 1, cnt
    call get_command_argument (i, c, nlen, status)
    if (status .ne. 0) then
        write (*,*) 'get_command_argument failed: status = ', status, ' arg = ', i
        stop
    end if
    if (i == 1) then
!        read(c(1:len),'(i)') ncpu
        read(c(1:nlen),*) ncpu															! --> ncpu
        write(*,*) 'Requested threads: ',ncpu
    elseif (i == 2) then
        infile = c(1:nlen)																! --> infile
        write(*,*) 'Input file: ',infile
    elseif (i == 3) then
        outfile = c(1:nlen)																! --> outfile
        write(*,*) 'Output file: ',outfile
!    elseif (i == 4) then
!        resfile = c(1:len)																! --> resfile
!        write(*,*) 'Result file: ',resfile
    endif
end do

!ncpu = 1
!infile = 'expta3.inp'

!runfile = 'running.out'

write(*,*) 'call execute!:'
call execute(ncpu,infile,inbuflen,outfile,outbuflen)
write(*,*) 'Nsteps: ',Nsteps
call get_Nsteps(Nsteps_tmp)
write(*,*) 'Nsteps_tmp: ',Nsteps_tmp
!call get_dimensions(NX,NY,NZ,Nsteps)
do jstep = 1,Nsteps_tmp
	call simulate_step(res)
	if (res < 0) then
		write(*,*) 'Error exit'
		stop
	endif
	if (res > 0) then
		write(*,*) 'Successful execution'
		exit
	endif
	if (mod(jstep,240) == 0) then
		call get_summary(summarydata)
!		hour = summaryData(1)
!		ntot = summaryData(4)
!		inflow = summaryData(10)
!		exits = summaryData(11)
!		write(*,'(5(a,i8))') 'Hour: ',hour,' ncells: ',ntot,' ncog: ',ncog(1),' inflow: ',inflow,' nexits: ', exits
	endif
enddo
call terminate_run(0)
end

