module err_exit
  ! Set all Global values and routines to private by default
  ! and then explicitly set their exposure
  !----------------------------------------------------------
  implicit none
  private

  public :: iulog
  public :: endrun

  integer, parameter::iulog=6

  contains
    !=====================================================================
    subroutine endrun(Message)
      character(len=*) Message
      write(iulog,*) Message
      stop
    end subroutine
    !=====================================================================

end module err_exit
