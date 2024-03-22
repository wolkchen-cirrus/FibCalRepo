program final_glmt

implicit none

integer jobID
character(100) jobIDString
real beamWidthX
real beamWidthY
real minRadius
real maxRadius
integer steps
real, allocatable, dimension(:) :: glmt_dataset
integer i, j, k, n, io
real, allocatable, dimension(:) :: radii
character(100) radiusString
real, allocatable, dimension(:) :: theta
real junk, thetaIn, p11In
real, allocatable, dimension(:) :: glmt_scattcross1, glmt_scattcross2
real dtheta
real pi

pi = 3.14159265
print*,'start final_glmt'

! need to read in several parameters
! job ID, so we know which directory to read and write files in
! min beam dimension, so we know where to fit glmt to Mie

! general code outline:
! read parameters
! read in glmt datasets
! trim to 0 to 180 degree phase function
! compute scattering cross section
! read in mie cross sections
! do some fitting to adjust the normalisation of glmt
! output scatt cross for mie and glmt


! read job ID
open(unit = 10, file = "jobID.txt")
	read(10,*) jobID
close(10)
print*,'job ID:',jobID

! make job ID string
write(jobIDString,*)jobID
call StripSpaces(jobIDString)
print*,"Directory: job"//trim(jobIDString)//"/"

! read job parameters
open(unit = 10, file = "Parameters of incident beam.txt")
	read(10,*) beamWidthX,beamWidthY
close(10)
open(unit = 10, file = "Parameters for particle.txt")
	read(10,*) minRadius,maxRadius,steps
close(10)
print*,'beamWidthX',beamWidthX
print*,'beamWidthY',beamWidthY
print*,'Smallest beam dimension: ',min(beamWidthX,beamWidthY)
print*,'min radius: ',minRadius
print*,'max radius: ',maxRadius
print*,'number of glmt datasets: ',steps

allocate(radii(steps))
open(unit = 10, file = "input.txt")

	do i = 1, steps
		read(10,*)radii(i)
	end do
	print*,'particle radii:', radii
close(10)

! all glmt datasets have 3601 vals
! forward direction starts at 1802

allocate(glmt_dataset(1:1801))
allocate(theta(1:1801))
allocate(glmt_scattcross1(1:steps)) ! glmt scatt cross before adjustment
allocate(glmt_scattcross2(1:steps)) ! glmt scatt cross after adjustment

j = 0
do i = 1, steps
	write(radiusString,'(f10.6)')radii(i)
	call StripSpaces(radiusString)
	print*,"radius string:"//trim(radiusString)
	! print*,'trying to open file:',"job"//trim(jobIDString)//"/LAM_Far field scattering_radius_"//trim(radiusString)//".xls"
	open(unit = 10, file = "job"//trim(jobIDString)//"/LAM_Far field scattering_radius_"//trim(radiusString)//".xls")
	! print*,'opened file at location: ',"job"//"jobID/LAM_Far field scattering_radius_"//trim(radiusString)//".xls"
	j = 0 ! counter for number of lines
    do  ! scan through the lines in the file...
        read(10,*,iostat=io) thetaIn, junk, p11In
        if (io/=0) exit
        ! print*,'thetaIn',thetaIn,"p11In",p11In
        if(thetaIn .ge. 0) then
        	j = j + 1
        	theta(j) = thetaIn
        	glmt_dataset(j) = p11In
        end if
    end do
    ! do k = 1,1801
    ! 	print*,'theta',theta(k),' p11:',glmt_dataset(k)
    ! end do
    ! print*,'finished reading this file. number of theta vals >= 0 :',j
    close(10)

    ! we've just read in the i^th glmt dataset, now we compute scattering cross section and add to an array

    glmt_scattcross1(i) = 0 ! set starting cross section to 0
    do k = 2, 1801
    	dtheta = theta(k) - theta(k-1)
    	glmt_scattcross1(i) = glmt_scattcross1(i) + 0.5*glmt_dataset(k)*sin(theta(k-1)*pi/180)*dtheta*pi/180
    end do
    print*,'scatt cross: ',glmt_scattcross1(i)
end do

open(unit = 10, file = "job"//trim(jobIDString)//"/CSCAT_GLMT")

	write(10,*)"AA		CSCAT"

	do i = 1, steps
		write(10,*) radii(i), glmt_scattcross1(i)
	end do

close(10)
















contains

subroutine StripSpaces(string) ! taken from: https://stackoverflow.com/questions/27179549/removing-whitespace-in-string
character(len=*) :: string
integer :: stringLen 
integer :: last, actual

stringLen = len (string)
last = 1
actual = 1

do while (actual < stringLen)
    if (string(last:last) == ' ') then
        actual = actual + 1
        string(last:last) = string(actual:actual)
        string(actual:actual) = ' '
    else
        last = last + 1
        if (actual < last) &
            actual = last
    endif
end do

end subroutine

end program final_glmt ! program rt containing all subroutines
