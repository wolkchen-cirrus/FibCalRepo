program init_glmt

use ifport !  includes functions and subroutines that ease porting of code to or from a PC

implicit none

integer(8) i, ste
integer(8) diameter_steps
real(8) diameter_start, diameter_end
real(8) rbi, ibi
real(8) beam_y,beam_x, lambda
logical exists ! for checking if directories or files exist
logical result ! logical for navigating directories
character(100) cwd ! current working directory
integer jobID
character(100) jobIDString

print*,'initialising'

open(unit = 10, file = "Parameters for particle.txt")

read(10,*) diameter_start, diameter_end, diameter_steps
read(10,*) rbi, ibi

close(10)

print*,'diameter_start',diameter_start
print*,'diameter_end',diameter_end
print*,'diameter_steps',diameter_steps

open(unit = 10, file = "input.txt")

do i = 0, diameter_steps-1

	write(10,*) diameter_start + i*(diameter_end - diameter_start)/(diameter_steps-1)

end do

close(10)

open(unit = 10, file = "Parameters of incident beam.txt")

	read(10,*) beam_x, beam_y, lambda

close(10)

print*,'beam width x:',beam_x
print*,'beam width y:',beam_y
print*,'lambda:',lambda

open(unit = 10, file = "mie_params.txt")

	write(10,*)lambda
	write(10,*)rbi, ibi

close(10)

result = getcwd(cwd)

jobID = 1
write(jobIDString,*)jobID
call StripSpaces(jobIDString)
print*,"job"//trim(jobIDString)
exists = .false.
inquire(directory = trim(cwd)//"/"//"job"//trim(jobIDString), exist = exists)
print*,'outside sr, dir exists?',exists

do while (exists .eq. .true.)
	jobID = jobID + 1
	write(jobIDString,*)jobID
	call StripSpaces(jobIDString)
	inquire(directory = trim(cwd)//"/"//"job"//trim(jobIDString), exist = exists)
end do

call makeDir("job"//trim(jobIDString))

open(unit = 10, file = "jobID.txt")

	write(10,*) jobID

close(10)

contains

subroutine makeDir(dirName)
  character(len=*) dirName
  logical exists ! for checking if directories or files exist
  logical result ! true if subdirectory was made, false if subdirectory was not made

inquire(directory = trim(cwd)//"/"//dirName, exist = exists)
if(exists .eq. .false.) then 
  print*,'Creating new subdirectory at location:',trim(cwd)//"/"//dirName
  result = makedirqq(trim(cwd)//"/"//dirName)
  print*,'success?',result
end if
end subroutine makeDir

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
    end if
end do

end subroutine

end program init_glmt ! program rt containing all subroutines