#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


module derivative_mod

  ! Useful modules
  !---------------
  use SE_Constants  ,only: real_kind, longdouble_kind
  use SE_Constants  ,only: np, nc, npdg, nep
  use SE_Constants  ,only: rrearth 
  use err_exit      ,only: endrun
  use quadrature_mod,only: quadrature_t, gauss, gausslobatto,legendre, jacobi
  use element_mod   ,only: element_t

  ! Set all Global values and routines to private by default
  ! and then explicitly set their exposure
  !-----------------------------------------------------------
  implicit none
  private

  public :: derivative_t
  public :: derivative_stag_t

  public :: derivinit
  public :: deriv_print
  private:: dmatinit
  private:: dpvinit
  public :: v2pinit
  private:: dvvinit
  public :: divergence
  private:: divergence_stag
  private:: divergence_nonstag
  public :: gradient_wk
  private:: gradient_wk_stag
  private:: gradient_wk_nonstag
  public :: gradient
  private:: gradient_str_stag
  private:: gradient_str_nonstag
  public :: vorticity
  public :: interpolate_gll2fvm_points
  public :: interpolate_gll2spelt_points
  public :: interpolate_gll2fvm_corners
  public :: remap_phys2gll

  ! these routines compute spherical differential operators as opposed to
  ! the gnomonic coordinate operators above.  Vectors (input or output)
  ! are always expressed in lat-lon coordinates
  !
  ! note that weak derivatives (integrated by parts form) can be defined using
  ! contra or co-variant test functions, so 
  !-------------------------------------------------------------------------
  public :: gradient_sphere
  public :: curl_sphere_wk_testcov
  public :: gradient_sphere_wk_testcov
  public :: gradient_sphere_wk_testcontra   ! only used for debugging
  public :: ugradv_sphere
  public :: curl_sphere
!  public :: curl_sphere_wk_testcontra  ! not coded
  public :: divergence_sphere_wk
  public :: element_boundary_integral
  public :: edge_flux_u_cg
  public :: vorticity_sphere
  public :: vorticity_sphere_diag
  public :: divergence_sphere
  public :: laplace_sphere_wk
  public :: vlaplace_sphere_wk
  private:: vlaplace_sphere_wk_cartesian
  private:: vlaplace_sphere_wk_contra

  public :: gll_to_dgmodal
  public :: dgmodal_to_gll
  public :: subcell_integration
  private:: allocate_subcell_integration_matrix

  ! PFC: Hack to set global option values
  !--------------------------------------
  public :: init_derivative_mod
  private:: hypervis_scaling
  private:: hypervis_power
  real(kind=real_kind):: hypervis_scaling
  real(kind=real_kind):: hypervis_power

  ! Public Interfaces
  !-------------------
  interface divergence
    module procedure divergence_nonstag
    module procedure divergence_stag
  end interface

  interface gradient
    module procedure gradient_str_nonstag
    module procedure gradient_str_stag
  end interface

  interface gradient_wk
    module procedure gradient_wk_nonstag
    module procedure gradient_wk_stag
  end interface

  ! Type Definitions
  !--------------------
  type derivative_t
    real(kind=real_kind):: Dvv     (np,np)
    real(kind=real_kind):: Dvv_diag(np,np)
    real(kind=real_kind):: Dvv_twt (np,np)
    real(kind=real_kind):: Mvv_twt (np,np)  ! diagonal matrix of GLL weights
    real(kind=real_kind):: vvtemp  (np,np)
    real(kind=real_kind):: vvtempt (np,np,2)
    real(kind=real_kind):: Mfvm    (np,nc+1)
    real(kind=real_kind):: Cfvm    (np,nc)
    real(kind=real_kind):: Sfvm    (np,nep)
    real(kind=real_kind):: legdg   (np,np)
  end type derivative_t

  type derivative_stag_t
    real(kind=real_kind):: D     (np,np)
    real(kind=real_kind):: M     (np,np)
    real(kind=real_kind):: Dpv   (np,np)
    real(kind=real_kind):: D_twt (np,np)
    real(kind=real_kind):: M_twt (np,np)
    real(kind=real_kind):: M_t   (np,np)
    real(kind=real_kind):: vtemp (np,np,2)
    real(kind=real_kind):: vtempt(np,np,2)
  end type derivative_stag_t

  ! Global Data
  !---------------------
  real(kind=real_kind),allocatable:: integration_matrix(:,:)

contains
  !==================================================================
  subroutine init_derivative_mod(I_SEopt)
    ! PFC: Hack to set global option values
    !--------------------------------------
    use SE_Options,only: SEoptions_t
    ! Passed Variables
    !------------------
    type(SEoptions_t),intent(in):: I_SEopt

    hypervis_scaling = I_SEopt%hypervis_scaling
    hypervis_power   = I_SEopt%hypervis_power

    ! End Routine
    !-------------
    return
  end subroutine init_derivative_mod
  !==================================================================


  !==================================================================
  subroutine derivinit(deriv,fvm_corners, fvm_points, spelt_refnep)
    ! derivinit: Initialize the matrices for taking derivatives 
    !            and interpolating
    !================================================================
    !
    ! Passed Variables
    !-----------------
    type(derivative_t)                 :: deriv
    real(kind=longdouble_kind),optional:: fvm_corners (nc+1)
    real(kind=longdouble_kind),optional:: fvm_points  (nc)
    real(kind=longdouble_kind),optional:: spelt_refnep(nep)
!    real (kind=longdouble_kind),optional:: phys_points(:)
    !
    ! Local Values
    !-------------
    type(quadrature_t)        :: gp   ! Quadrature points and weights on pressure grid
    real(kind=longdouble_kind):: dmat(np,np)
    real(kind=longdouble_kind):: dpv(np,np)
    real(kind=longdouble_kind):: v2p(np,np)
    real(kind=longdouble_kind):: p2v(np,np)
    real(kind=longdouble_kind):: dvv(np,np)
    real(kind=longdouble_kind):: dvv_diag(np,np)
    real(kind=longdouble_kind):: v2v(np,np)
    real(kind=longdouble_kind):: xnorm

    integer ii,jj

    ! initialize matrices in longdouble_kind precision
    ! and transfer results into real_kind floating point precision
    !--------------------------------------------------------------
    gp=gausslobatto(np)

    ! Legendre polynomials of degree npdg-1, on the np GLL grid:
    !-----------------------------------------------------------
    if(npdg>np) call endrun( 'FATAL ERROR: npdg>np')

    if((npdg>0).and.(npdg<np)) then
      ! in this case, we use a DG basis of Legendre polynomials
      ! stored at the GLL points and normalize.
      !--------------------------------------------------------
      do ii=1,np
        deriv%legdg(1:npdg,ii) = legendre(gp%points(ii),npdg-1)
      end do
      do jj=1,npdg
        xnorm=sqrt(sum(deriv%legdg(jj,:)*deriv%legdg(jj,:)*gp%weights(:)))
        deriv%legdg(jj,:)=deriv%legdg(jj,:)/xnorm
      end do
    endif

    call dvvinit(dvv,gp)
!DIAG******************************************
!    print *,' dvv(1)=',dvv(:,1)
!    print *,' dvv(2)=',dvv(:,2)
!    print *,' dvv(3)=',dvv(:,3)
!    print *,' dvv(4)=',dvv(:,4)
!DIAG******************************************
    deriv%Dvv(:,:)=dvv(:,:)

    do ii=1,np
    do jj=1,np
      if(ii.eq.jj) then
        deriv%dvv_diag(ii,jj) = dvv(ii,jj)
      else
        deriv%dvv_diag(ii,jj) = 0.0D0
      endif 
    end do
    end do

    v2v = 0.0D0
    do ii=1,np
      v2v(ii,ii) = gp%weights(ii)
    end do

    do ii=1,np
    do jj=1,np
      dvv(jj,ii) = dvv(jj,ii)*gp%weights(ii)
    end do
    end do

    deriv%Dvv_twt = transpose(dvv)
    deriv%Mvv_twt = v2v
    if(present(fvm_corners )) call v2pinit(deriv%Mfvm,gp%points,fvm_corners ,np,nc+1)
    if(present(fvm_points  )) call v2pinit(deriv%Cfvm,gp%points,fvm_points  ,np,nc  )
    if(present(spelt_refnep)) call v2pinit(deriv%Sfvm,gp%points,spelt_refnep,np,nep )
         
    ! notice we deallocate this memory here even though it was allocated 
    ! by the call to gausslobatto.
    !-------------------------------------------------------------------
    deallocate(gp%points)
    deallocate(gp%weights)

    ! End Routine
    !-------------
    return
  end subroutine derivinit
  !==================================================================


  !==================================================================
  subroutine deriv_print(deriv)
    ! deriv_print: 
    !
    !================================================================
    ! 
    ! Passed Variables
    !-----------------
    type(derivative_t):: deriv
    !
    ! Local Values
    !---------------
    integer jj

    print *,"Derivative Matrix Dvv"
    do jj=1,np
      write(6,*) deriv%Dvv(:,jj)
    end do

    print *,"Weak Derivative Matrix Dvv_twt"
    do jj=1,np
      write(6,*) deriv%Dvv_twt(:,jj)
    end do

    ! End Routine
    !-------------
    return
  end subroutine deriv_print
  !==================================================================


  !==================================================================
  subroutine dmatinit(dmat)
    ! dmatinit: Compute rectangular v->p derivative matrix (dmat)
    !=============================================================
    !
    ! Passed Variables
    !------------------
    real(kind=longdouble_kind):: dmat(np,np)
    !
    ! Local Values
    !-------------
    type(quadrature_t)        :: gll
    type(quadrature_t)        :: gs
    real(kind=longdouble_kind):: fact,f1,f2
    real(kind=longdouble_kind):: func0,func1
    real(kind=longdouble_kind):: dis,c0,c1
    real(kind=longdouble_kind):: leg(np,np)
    real(kind=longdouble_kind)::  jac(0:np-1)
    real(kind=longdouble_kind):: djac(0:np-1)

    integer ii,jj

    c0 = 0.0_longdouble_kind
    c1 = 1.0_longdouble_kind

    gll= gausslobatto(np)
    gs = gauss       (np)

    ! Compute Legendre polynomials on Gauss-Lobatto grid (velocity)
    !--------------------------------------------------------------
    do ii=1,np
      leg(:,ii) = legendre(gll%points(ii),np-1)
    end do

    !  Derivatives of velocity cardinal functions on pressure grid
    !  d(i,j) = D(j,i) = D' (D-transpose) since D(i,j) = dh_j(x_i)/dx
    !-----------------------------------------------------------------
    fact = np*(np-1)
    do jj=1,np
      call jacobi(np-1,gs%points(jj),c0,c0,jac(0:np-1),djac(0:np-1))
      func0 =  jac(np-1)
      func1 = djac(np-1)
      f1 = fact*func0
      f2 = (c1 - gs%points(jj))*(c1 + gs%points(jj))*func1
      do ii=1,np
        if( gs%points(jj) /= gll%points(ii) ) then
          dis = gs%points(jj) - gll%points(ii)
          dmat(ii,jj) = func0/(leg(np,ii)*dis ) + f2/(fact*leg(np,ii)*dis*dis)
!!! OTHER           dmat(ii,jj) = (1.0D0/(fact*leg(np,ii)*dis*dis))* (func0*fact*dis + f2)
        else
          dmat(ii,jj) = c0
        endif
      end do
    end do

    ! Clean up
    !----------
    deallocate(gll%points )
    deallocate(gll%weights)
    deallocate( gs%points )
    deallocate( gs%weights)

    ! End Routine
    !-------------
    return
  end subroutine dmatinit
  !==================================================================


  !==================================================================
  subroutine dpvinit(dmat)
    ! dpvinit: Compute rectangular p->v derivative matrix (dmat) 
    !          for strong gradients
    !==============================================================
    !
    ! Passed Variables
    !------------------
    real(kind=longdouble_kind):: dmat(np,np)
    !
    ! Local Values
    !-------------
    type(quadrature_t)        :: gll
    type(quadrature_t)        :: gs
    real(kind=longdouble_kind):: dis,c0,c1
    real(kind=longdouble_kind)::  legv(0:np,np)
    real(kind=longdouble_kind):: dlegv(0:np,np)
    real(kind=longdouble_kind)::  leg(0:np)
    real(kind=longdouble_kind):: dleg(0:np)

    integer ii,jj

    c0 = 0.0_longdouble_kind
    c1 = 1.0_longdouble_kind

    gll= gausslobatto(np)
    gs = gauss(np)

    ! Compute Legendre polynomials on Gauss-Lobatto grid (velocity)
    !--------------------------------------------------------------
    do ii=1,np
      call jacobi(np,gll%points(ii),c0,c0,legv(0:np,ii),dlegv(0:np,ii))
    end do

    !  Derivatives of velocity cardinal functions on pressure grid
    !  d(i,j) = D(j,i) = D' (D-transpose) since D(i,j) = dh_j(x_i)/dx
    !----------------------------------------------------------------
    do jj=1,np
      call jacobi(np,gs%points(jj),c0,c0,leg(0:np),dleg(0:np))
      do ii=1,np
        if(gs%points(jj) /= gll%points(ii)) then
          dis = gll%points(ii) - gs%points(jj)
          dmat(jj,ii) = dlegv(np,ii)/(dleg(np)*dis) -  legv(np,ii)/(dleg(np)*dis*dis)
        else
          dmat(jj,ii) = c0
        endif
      end do
    end do

    ! Clean up
    !----------
    deallocate(gll%points )
    deallocate(gll%weights)
    deallocate( gs%points )
    deallocate( gs%weights)

    ! End Routine
    !-------------
    return
  end subroutine dpvinit
  !==================================================================


  !==================================================================
  subroutine v2pinit(v2p,gll,gs,n1,n2)
    ! v2pinit: Compute interpolation matrix from gll(1:n1) -> gs(1:n2)
    !================================================================
    ! 
    ! Passed Variables
    !-------------------
    real(kind=real_kind)      :: v2p(n1,n2)
    real(kind=longdouble_kind):: gll(n1),gs(n2)
    integer                   :: n1,n2
    !
    ! Local Values
    !-----------------
    real(kind=real_kind)      :: v2p_new(n1,n2)
    real(kind=longdouble_kind):: fact,f1, sum
    real(kind=longdouble_kind):: func0,func1
    real(kind=longdouble_kind):: leg(n1,n1)
    real(kind=longdouble_kind)::  jac(0:n1-1)
    real(kind=longdouble_kind):: djac(0:n1-1)
    real(kind=longdouble_kind):: c0,c1
    type(quadrature_t)        :: gll_pts
    real(kind=longdouble_kind):: leg_out(n1,n2)
    real(kind=longdouble_kind):: gamma(n1)

    integer ii,jj,kk,mm,ll

    c0 = 0.0_longdouble_kind
    c1 = 1.0_longdouble_kind

    ! Compute Legendre polynomials on Gauss-Lobatto grid (velocity)
    !--------------------------------------------------------------
    fact = -n1*(n1-1)
    do ii=1,n1
      leg( :,ii) = legendre(gll(ii),n1-1)
      leg(n1,ii) = fact*leg(n1,ii)
    end do

    ! Velocity cardinal functions on pressure grid
    !-----------------------------------------------
#if 0
    do jj=1,n2
      call jacobi(n1-1,gs(jj),c0,c0,jac(0:n1-1),djac(0:n1-1))
      func0 =  jac(n1-1)
      func1 = djac(n1-1)
      f1 = (c1 - gs(jj)**2)*func1
      do ii=1,n1
        if(gs(jj) /= gll(ii)) then
          v2p(ii,jj) = f1/(leg(n1,ii)*(gs(jj)-gll(ii)))
        else
          v2p(ii,jj) = c1
        endif
      end do
    end do
#endif
    ! NEW VERSION, with no division by (gs(j)-gll(i)):
    ! compute legendre polynomials at output points:
    !------------------------------------------------
    gll_pts = gausslobatto(n1)
    fact    = -n1*(n1-1)
    do ii=1,n2
      leg_out( :,ii) = legendre(gs(ii),n1-1)
      leg_out(n1,ii) = fact*leg_out(n1,ii)
    end do

    ! compute gamma: (normalization factor for inv(leg)
    !---------------------------------------------------
    do mm=1,n1
      gamma(mm)=0.0_longdouble_kind
      do ii=1,n1
        gamma(mm)=gamma(mm)+leg(mm,ii)*leg(mm,ii)*gll_pts%weights(ii) 
      end do
      gamma(mm)=1.0_longdouble_kind/gamma(mm)
    end do

    ! compute product of leg_out * inv(leg):
    !----------------------------------------
    do jj=1,n2    ! this should be fvm points
    do ll=1,n1    ! this should be GLL points
      sum=0.0_longdouble_kind
      do kk=1,n1  ! number of polynomials = number of GLL points
        sum=sum + leg_out(kk,jj)*gamma(kk)*leg(kk,ll)
      end do
      v2p_new(ll,jj) = gll_pts%weights(ll)*sum
    end do
    end do
#if 0
    do jj=1,n2   ! this should be fvm points
    do ll=1,n1   ! this should be GLL points
      print *,ll,jj,v2p_new(ll,jj),v2p(ll,jj)
    end do
    end do
    print *,'max error: ',maxval(abs(v2p_new-v2p))
#endif
    v2p(:,:)=v2p_new(:,:)

    ! Clean up
    !---------------
    deallocate(gll_pts%points )
    deallocate(gll_pts%weights)

    ! End Routine
    !-------------
    return
  end subroutine v2pinit
  !==================================================================


  !==================================================================
  subroutine dvvinit(dvv,gll)
    ! dvvinit: Compute rectangular v->v derivative matrix (dvv)
    !==========================================================
    !
    ! Passed Variables
    !-----------------
    real(kind=longdouble_kind):: dvv(np,np)
    type(quadrature_t)        :: gll
    !
    ! Local Values
    !-----------------
    real(kind=longdouble_kind):: leg(np,np)
    real(kind=longdouble_kind):: c0,c1,c4

    integer ii,jj

    c0 = 0.0_longdouble_kind
    c1 = 1.0_longdouble_kind
    c4 = 4.0_longdouble_kind

    do ii=1,np
      leg(:,ii) = legendre(gll%points(ii),np-1)
    end do

    dvv(:,:) = c0
    do jj=1,np
      do ii=1,jj-1
        dvv(jj,ii) = (c1/(gll%points(ii)-gll%points(jj)))*leg(np,ii)/leg(np,jj)
      end do
      dvv(jj,jj) = c0
      do ii=jj+1,np
        dvv(jj,ii) = (c1/(gll%points(ii)-gll%points(jj)))*leg(np,ii)/leg(np,jj)
      end do
    end do

    dvv(np,np) = + np*(np-1)/c4
    dvv(1,1)   = - np*(np-1)/c4

    ! End Routine
    !-------------
    return
  end subroutine dvvinit
  !==================================================================


  !==================================================================
  function divergence_stag(v,deriv) result(div)
    ! divergence_stag: Compute divergence (maps v grid -> p grid)
    !=============================================================
    !
    ! Passed Variables
    !-------------------
    real(kind=real_kind),intent(in):: v(np,np,2)
    type(derivative_stag_t)        :: deriv
    real(kind=real_kind)           :: div(np,np)
    !
    ! Local Values
    !---------------
    logical,parameter   :: UseUnroll = .TRUE.
    real(kind=real_kind)::  sumx00,sumx01
    real(kind=real_kind)::  sumy00,sumy01
    real(kind=real_kind)::  sumx10,sumx11
    real(kind=real_kind)::  sumy10,sumy11

    integer ii,jj,ll

#ifdef DEBUG
    print *, "divergence_stag"
#endif

    if((modulo(np,2) == 0).and.(UseUnroll)) then 
      !JMD====================================
      !JMD  2*np*np*np Flops
      !JMD====================================
      do jj=1,np,2
      do ll=1,np,2
        sumx00=0.0d0
        sumx01=0.0d0
        sumx10=0.0d0
        sumx11=0.0d0
        sumy00=0.0d0
        sumy01=0.0d0
        sumy10=0.0d0
        sumy11=0.0d0
        do ii=1,np
          sumx00 = sumx00 + deriv%D(ii,ll  )*v(ii,jj  ,1)
          sumx01 = sumx01 + deriv%D(ii,ll+1)*v(ii,jj  ,1)
          sumx10 = sumx10 + deriv%D(ii,ll  )*v(ii,jj+1,1)
          sumx11 = sumx11 + deriv%D(ii,ll+1)*v(ii,jj+1,1)
          sumy00 = sumy00 + deriv%M(ii,ll  )*v(ii,jj  ,2)
          sumy01 = sumy01 + deriv%M(ii,ll+1)*v(ii,jj  ,2)
          sumy10 = sumy10 + deriv%M(ii,ll  )*v(ii,jj+1,2)
          sumy11 = sumy11 + deriv%M(ii,ll+1)*v(ii,jj+1,2)
        end do
        deriv%vtemp(jj  ,ll  ,1) = sumx00
        deriv%vtemp(jj  ,ll+1,1) = sumx01
        deriv%vtemp(jj+1,ll  ,1) = sumx10
        deriv%vtemp(jj+1,ll+1,1) = sumx11
        deriv%vtemp(jj  ,ll  ,2) = sumy00
        deriv%vtemp(jj  ,ll+1,2) = sumy01
        deriv%vtemp(jj+1,ll  ,2) = sumy10
        deriv%vtemp(jj+1,ll+1,2) = sumy11
      end do
      end do

      !JMD====================================
      !JMD  2*np*np*np Flops
      !JMD====================================
      do jj=1,np,2
      do ii=1,np,2
        sumx00=0.0d0
        sumx01=0.0d0
        sumx10=0.0d0
        sumx11=0.0d0
        sumy00=0.0d0
        sumy01=0.0d0
        sumy10=0.0d0
        sumy11=0.0d0
        do ll=1,np
          sumx00 = sumx00 +  deriv%M(ll,jj  )*deriv%vtemp(ll,ii  ,1)
          sumx01 = sumx01 +  deriv%M(ll,jj+1)*deriv%vtemp(ll,ii  ,1)
          sumx10 = sumx10 +  deriv%M(ll,jj  )*deriv%vtemp(ll,ii+1,1)
          sumx11 = sumx11 +  deriv%M(ll,jj+1)*deriv%vtemp(ll,ii+1,1)
          sumy00 = sumy00 +  deriv%D(ll,jj  )*deriv%vtemp(ll,ii  ,2)
          sumy01 = sumy01 +  deriv%D(ll,jj+1)*deriv%vtemp(ll,ii  ,2)
          sumy10 = sumy10 +  deriv%D(ll,jj  )*deriv%vtemp(ll,ii+1,2)
          sumy11 = sumy11 +  deriv%D(ll,jj+1)*deriv%vtemp(ll,ii+1,2)
        end do
        div(ii  ,jj  ) = sumx00 + sumy00
        div(ii  ,jj+1) = sumx01 + sumy01
        div(ii+1,jj  ) = sumx10 + sumy10
        div(ii+1,jj+1) = sumx11 + sumy11
      end do
      end do
    else
      do jj=1,np
      do ll=1,np
        sumx00=0.0d0
        sumy00=0.0d0
        do ii=1,np
          sumx00 = sumx00 + deriv%D(ii,ll)*v(ii,jj,1)
          sumy00 = sumy00 + deriv%M(ii,ll)*v(ii,jj,2)
        end do
        deriv%vtemp(jj,ll,1) = sumx00
        deriv%vtemp(jj,ll,2) = sumy00
      end do
      end do

      do jj=1,np
      do ii=1,np
        sumx00=0.0d0
        sumy00=0.0d0
        do ll=1,np
          sumx00 = sumx00 +  deriv%M(ll,jj)*deriv%vtemp(ll,ii,1)
          sumy00 = sumy00 +  deriv%D(ll,jj)*deriv%vtemp(ll,ii,2)
        end do
        div(ii,jj) = sumx00 + sumy00
      end do
      end do
    endif

    ! End Function
    !-------------
    return
  end function divergence_stag
  !==================================================================


  !==================================================================
  function divergence_nonstag(v,deriv) result(div)
    !  divergence_nonstag: Compute divergence (maps v->v)
    !=====================================================
    !
    ! Passed Variables
    !------------------
    real(kind=real_kind),intent(in):: v(np,np,2)
    type(derivative_t)             :: deriv
    real(kind=real_kind)           :: div(np,np)
    !
    ! Local Values
    !-------------
    logical,parameter   :: UseUnroll = .TRUE.
    real(kind=real_kind):: dudx00,dudx01
    real(kind=real_kind):: dudx10,dudx11
    real(kind=real_kind):: dvdy00,dvdy01
    real(kind=real_kind):: dvdy10,dvdy11

    integer ii,jj,ll

    if((modulo(np,2).eq.0).and.(UseUnroll)) then
      ! this is just loop unrolling - a good compiler should do it for you jpe
      !------------------------------------------------------------------------
      do jj=1,np,2
      do ll=1,np,2
        dudx00=0.0d0
        dudx01=0.0d0
        dudx10=0.0d0
        dudx11=0.0d0
        dvdy00=0.0d0
        dvdy01=0.0d0
        dvdy10=0.0d0
        dvdy11=0.0d0
        do ii=1,np
          dudx00 = dudx00 + deriv%Dvv(ii,ll  )*v(ii,jj  ,1)
          dudx01 = dudx01 + deriv%Dvv(ii,ll+1)*v(ii,jj  ,1)
          dudx10 = dudx10 + deriv%Dvv(ii,ll  )*v(ii,jj+1,1)
          dudx11 = dudx11 + deriv%Dvv(ii,ll+1)*v(ii,jj+1,1)
          dvdy00 = dvdy00 + deriv%Dvv(ii,ll  )*v(jj  ,ii,2)
          dvdy01 = dvdy01 + deriv%Dvv(ii,ll+1)*v(jj  ,ii,2)
          dvdy10 = dvdy10 + deriv%Dvv(ii,ll  )*v(jj+1,ii,2)
          dvdy11 = dvdy11 + deriv%Dvv(ii,ll+1)*v(jj+1,ii,2)
        end do
        div(ll  ,jj  ) = dudx00
        div(ll+1,jj  ) = dudx01
        div(ll  ,jj+1) = dudx10
        div(ll+1,jj+1) = dudx11
        deriv%vvtemp(jj  ,ll  ) = dvdy00
        deriv%vvtemp(jj  ,ll+1) = dvdy01
        deriv%vvtemp(jj+1,ll  ) = dvdy10
        deriv%vvtemp(jj+1,ll+1) = dvdy11
      end do
      end do
    else
      do jj=1,np
      do ll=1,np
        dudx00=0.0d0
        dvdy00=0.0d0
        do ii=1,np
          dudx00 = dudx00 + deriv%Dvv(ii,ll)*v(ii,jj,1)
          dvdy00 = dvdy00 + deriv%Dvv(ii,ll)*v(jj,ii,2)
        end do
        div(ll,jj)          = dudx00
        deriv%vvtemp(jj,ll) = dvdy00
      end do
      end do
    endif

    do jj=1,np
    do ii=1,np
      div(ii,jj)=div(ii,jj)+deriv%vvtemp(ii,jj)
    end do
    end do

    ! End Function
    !-------------
    return
  end function divergence_nonstag
  !==================================================================


  !==================================================================
  function gradient_wk_stag(p,deriv) result(dp)
    !  gradient_wk_stag: Compute the weak form gradient: maps scalar field on 
    !                    the pressure grid to the velocity grid
    !========================================================================
    !
    ! Passed Variables
    !---------------------
    real(kind=real_kind),intent(in):: p(np,np)
    type(derivative_stag_t)        :: deriv
    real(kind=real_kind)           :: dp(np,np,2)
    !
    ! Local Values
    !----------------
    logical,parameter   :: UseUnroll = .TRUE.
    real(kind=real_kind):: sumx00,sumx01
    real(kind=real_kind):: sumy00,sumy01
    real(kind=real_kind):: sumx10,sumx11
    real(kind=real_kind):: sumy10,sumy11

    integer ii,jj,ll

#ifdef DEBUG
    print *, "gradient_wk_stag"
#endif

    if((modulo(np,2) == 0).and.(UseUnroll)) then 
      !JMD ================================
      !JMD 2*np*np*np Flops 
      !JMD ================================
      do jj=1,np,2
      do ll=1,np,2
        sumx00=0.0d0
        sumx01=0.0d0
        sumx10=0.0d0
        sumx11=0.0d0
        sumy00=0.0d0
        sumy01=0.0d0
        sumy10=0.0d0
        sumy11=0.0d0
        do ii=1,np
          sumx00 = sumx00 + deriv%D_twt(ii,ll  )*p(ii,jj  )
          sumx01 = sumx01 + deriv%D_twt(ii,ll+1)*p(ii,jj  )
          sumx10 = sumx10 + deriv%D_twt(ii,ll  )*p(ii,jj+1)
          sumx11 = sumx11 + deriv%D_twt(ii,ll+1)*p(ii,jj+1)
          sumy00 = sumy00 + deriv%M_twt(ii,ll  )*p(ii,jj  )
          sumy01 = sumy01 + deriv%M_twt(ii,ll+1)*p(ii,jj  )
          sumy10 = sumy10 + deriv%M_twt(ii,ll  )*p(ii,jj+1)
          sumy11 = sumy11 + deriv%M_twt(ii,ll+1)*p(ii,jj+1)
        end do
        deriv%vtempt(jj  ,ll  ,1) = sumx00
        deriv%vtempt(jj  ,ll+1,1) = sumx01
        deriv%vtempt(jj+1,ll  ,1) = sumx10
        deriv%vtempt(jj+1,ll+1,1) = sumx11
        deriv%vtempt(jj  ,ll  ,2) = sumy00
        deriv%vtempt(jj  ,ll+1,2) = sumy01
        deriv%vtempt(jj+1,ll  ,2) = sumy10
        deriv%vtempt(jj+1,ll+1,2) = sumy11
      end do
      end do

      !JMD ================================
      !JMD 2*np*np*np Flops 
      !JMD ================================
      do jj=1,np,2
      do ii=1,np,2
        sumx00=0.0d0
        sumx01=0.0d0
        sumx10=0.0d0
        sumx11=0.0d0
        sumy00=0.0d0
        sumy01=0.0d0
        sumy10=0.0d0
        sumy11=0.0d0
        do ll=1,np
          sumx00 = sumx00 +  deriv%M_twt(ll,jj  )*deriv%vtempt(ll,ii  ,1)
          sumx01 = sumx01 +  deriv%M_twt(ll,jj+1)*deriv%vtempt(ll,ii  ,1)
          sumx10 = sumx10 +  deriv%M_twt(ll,jj  )*deriv%vtempt(ll,ii+1,1)
          sumx11 = sumx11 +  deriv%M_twt(ll,jj+1)*deriv%vtempt(ll,ii+1,1)
          sumy00 = sumy00 +  deriv%D_twt(ll,jj  )*deriv%vtempt(ll,ii  ,2)
          sumy01 = sumy01 +  deriv%D_twt(ll,jj+1)*deriv%vtempt(ll,ii  ,2)
          sumy10 = sumy10 +  deriv%D_twt(ll,jj  )*deriv%vtempt(ll,ii+1,2)
          sumy11 = sumy11 +  deriv%D_twt(ll,jj+1)*deriv%vtempt(ll,ii+1,2)
        end do
        dp(ii  ,jj  ,1) = sumx00
        dp(ii  ,jj+1,1) = sumx01
        dp(ii+1,jj  ,1) = sumx10
        dp(ii+1,jj+1,1) = sumx11
        dp(ii  ,jj  ,2) = sumy00
        dp(ii  ,jj+1,2) = sumy01
        dp(ii+1,jj  ,2) = sumy10
        dp(ii+1,jj+1,2) = sumy11
      end do
      end do
    else
      do jj=1,np
      do ll=1,np
        sumx00=0.0d0
        sumy00=0.0d0
        do ii=1,np
          sumx00 = sumx00 + deriv%D_twt(ii,ll)*p(ii,jj)
          sumy00 = sumy00 + deriv%M_twt(ii,ll)*p(ii,jj)
        end do
        deriv%vtempt(jj,ll,1) = sumx00
        deriv%vtempt(jj,ll,2) = sumy00
      end do
      end do

      do jj=1,np
      do ii=1,np
        sumx00=0.0d0
        sumy00=0.0d0
        do ll=1,np
          sumx00 = sumx00 +  deriv%M_twt(ll,jj)*deriv%vtempt(ll,ii,1)
          sumy00 = sumy00 +  deriv%D_twt(ll,jj)*deriv%vtempt(ll,ii,2)
        end do
        dp(ii,jj,1) = sumx00
        dp(ii,jj,2) = sumy00
      end do
      end do
    endif

    ! End Function
    !-------------
    return
  end function gradient_wk_stag
  !==================================================================


  !==================================================================
  function gradient_wk_nonstag(p,deriv) result(dp)
    !  gradient_wk_nonstag: Compute the weak form gradient: maps scalar field 
    !                       on the Gauss-Lobatto grid to the weak gradient on 
    !                       the Gauss-Lobbatto grid
    !=======================================================================
    !
    ! Passed Variables
    !------------------
    real(kind=real_kind),intent(in):: p(np,np)
    type(derivative_t)             :: deriv
    real(kind=real_kind)           :: dp(np,np,2)
    !
    ! Local Values
    !--------------
    logical,parameter   :: UseUnroll = .TRUE.
    real(kind=real_kind):: sumx00,sumx01
    real(kind=real_kind):: sumy00,sumy01
    real(kind=real_kind):: sumx10,sumx11
    real(kind=real_kind):: sumy10,sumy11

    integer ii,jj,ll

!   print *, "gradient_wk_nonstag"

    !JMD ================================
    !JMD 2*np*np*np Flops 
    !JMD ================================

    if(modulo(np,2) .eq. 0 .and. UseUnroll) then
      ! this is just loop unrolling - a good compiler should do it for you jpe
      !------------------------------------------------------------------------
      do jj=1,np,2
      do ll=1,np,2
        sumx00=0.0d0
        sumx01=0.0d0
        sumx10=0.0d0
        sumx11=0.0d0
        sumy00=0.0d0
        sumy01=0.0d0
        sumy10=0.0d0
        sumy11=0.0d0
        do ii=1,np
          sumx00 = sumx00 + deriv%Dvv_twt(ii,ll  )*p(ii,jj  )
          sumx01 = sumx01 + deriv%Dvv_twt(ii,ll+1)*p(ii,jj  )
          sumx10 = sumx10 + deriv%Dvv_twt(ii,ll  )*p(ii,jj+1)
          sumx11 = sumx11 + deriv%Dvv_twt(ii,ll+1)*p(ii,jj+1)
          sumy00 = sumy00 + deriv%Mvv_twt(ii,ll  )*p(ii,jj  )
          sumy01 = sumy01 + deriv%Mvv_twt(ii,ll+1)*p(ii,jj  )
          sumy10 = sumy10 + deriv%Mvv_twt(ii,ll  )*p(ii,jj+1)
          sumy11 = sumy11 + deriv%Mvv_twt(ii,ll+1)*p(ii,jj+1)
        end do
        deriv%vvtempt(jj  ,ll  ,1) = sumx00
        deriv%vvtempt(jj  ,ll+1,1) = sumx01
        deriv%vvtempt(jj+1,ll  ,1) = sumx10
        deriv%vvtempt(jj+1,ll+1,1) = sumx11
        deriv%vvtempt(jj  ,ll  ,2) = sumy00
        deriv%vvtempt(jj  ,ll+1,2) = sumy01
        deriv%vvtempt(jj+1,ll  ,2) = sumy10
        deriv%vvtempt(jj+1,ll+1,2) = sumy11
      end do
      end do

      ! vvtempt1 = p'*Dvv_twt
      ! vvtempt2 = p'*Mvv_twt
      ! dp1 = dy*Mvv_twt*vvtempt1' = dy*Mvv_twt*(p'*Dvv_twt)' = dy*Mvv_twt*Dvv_twt'*p
      ! dp2 = dx*Dvv_twt*vvtempt2' = dx*Dvv_twt*(p'*Mvv_twt)' = dx*Dvv_twt*Mvv_twt'*p
      !     New formulation 
      ! dp1 = dy*MvvDvvt*p
      ! dp2 = dx*DvvMvvt*p
      ! MvvDvvt = Mvv_twt*Dvv_twt'
      ! DvvMvvt = Dvv_twt*Mvv_twt'

      !JMD ================================
      !JMD 2*np*np*np Flops 
      !JMD ================================
      do jj=1,np,2
      do ii=1,np,2
        sumx00=0.0d0
        sumx01=0.0d0
        sumx10=0.0d0
        sumx11=0.0d0
        sumy00=0.0d0
        sumy01=0.0d0
        sumy10=0.0d0
        sumy11=0.0d0
        do ll=1,np
          sumx00 = sumx00 +  deriv%Mvv_twt(ll,jj  )*deriv%vvtempt(ll,ii  ,1)
          sumx01 = sumx01 +  deriv%Mvv_twt(ll,jj+1)*deriv%vvtempt(ll,ii  ,1)
          sumx10 = sumx10 +  deriv%Mvv_twt(ll,jj  )*deriv%vvtempt(ll,ii+1,1)
          sumx11 = sumx11 +  deriv%Mvv_twt(ll,jj+1)*deriv%vvtempt(ll,ii+1,1)
          sumy00 = sumy00 +  deriv%Dvv_twt(ll,jj  )*deriv%vvtempt(ll,ii  ,2)
          sumy01 = sumy01 +  deriv%Dvv_twt(ll,jj+1)*deriv%vvtempt(ll,ii  ,2)
          sumy10 = sumy10 +  deriv%Dvv_twt(ll,jj  )*deriv%vvtempt(ll,ii+1,2)
          sumy11 = sumy11 +  deriv%Dvv_twt(ll,jj+1)*deriv%vvtempt(ll,ii+1,2)
        end do
        dp(ii  ,jj  ,1) = sumx00
        dp(ii  ,jj+1,1) = sumx01
        dp(ii+1,jj  ,1) = sumx10
        dp(ii+1,jj+1,1) = sumx11
        dp(ii  ,jj  ,2) = sumy00
        dp(ii  ,jj+1,2) = sumy01
        dp(ii+1,jj  ,2) = sumy10
        dp(ii+1,jj+1,2) = sumy11
      end do
      end do
    else
      do jj=1,np
      do ll=1,np
        sumx00=0.0d0
        sumy00=0.0d0
        do ii=1,np
          sumx00 = sumx00 + deriv%Dvv_twt(ii,ll)*p(ii,jj)
          sumy00 = sumy00 + deriv%Mvv_twt(ii,ll)*p(ii,jj)
        end do
        deriv%vvtempt(jj,ll,1) = sumx00
        deriv%vvtempt(jj,ll,2) = sumy00
      end do
      end do

      !JMD ================================
      !JMD 2*np*np*np Flops 
      !JMD ================================

      do jj=1,np
      do ii=1,np
        sumx00=0.0d0
        sumy00=0.0d0
        do ll=1,np
          sumx00 = sumx00 +  deriv%Mvv_twt(ll,jj)*deriv%vvtempt(ll,ii,1)
          sumy00 = sumy00 +  deriv%Dvv_twt(ll,jj)*deriv%vvtempt(ll,ii,2)
        end do
        dp(ii,jj,1) = sumx00
        dp(ii,jj,2) = sumy00
      end do
      end do
    endif

    ! End Function
    !-------------
    return
  end function gradient_wk_nonstag
  !==================================================================


  !==================================================================
  function gradient_str_stag(p,deriv) result(dp)
    !  gradient_str_stag: Compute the *strong* form gradient: maps scalar field 
    !                     on the pressure grid to the velocity grid
    !======================================================================
    !
    ! Passed Variables
    !-------------------
    real(kind=real_kind),intent(in):: p(np,np)
    type(derivative_stag_t)        :: deriv
    real(kind=real_kind)           :: dp(np,np,2)
    !
    ! Local Values
    !---------------
    logical,parameter   :: UseUnroll=.TRUE.
    real(kind=real_kind):: sumx00,sumx01
    real(kind=real_kind):: sumy00,sumy01
    real(kind=real_kind):: sumx10,sumx11
    real(kind=real_kind):: sumy10,sumy11

    integer ii,jj,ll

#ifdef DEBUG
    print *, "gradient_str_stag"
#endif

    if((modulo(np,2) == 0).and.(UseUnroll)) then 
      do jj=1,np,2
      do ll=1,np,2
        sumx00=0.0d0
        sumx01=0.0d0
        sumx10=0.0d0
        sumx11=0.0d0
        sumy00=0.0d0
        sumy01=0.0d0
        sumy10=0.0d0
        sumy11=0.0d0
        do ii=1,np
          sumx00 = sumx00 + deriv%Dpv(ii,ll  )*p(ii,jj  )
          sumx01 = sumx01 + deriv%Dpv(ii,ll+1)*p(ii,jj  )
          sumx10 = sumx10 + deriv%Dpv(ii,ll  )*p(ii,jj+1)
          sumx11 = sumx11 + deriv%Dpv(ii,ll+1)*p(ii,jj+1)
          sumy00 = sumy00 + deriv%M_t(ii,ll  )*p(ii,jj  )
          sumy01 = sumy01 + deriv%M_t(ii,ll+1)*p(ii,jj  )
          sumy10 = sumy10 + deriv%M_t(ii,ll  )*p(ii,jj+1)
          sumy11 = sumy11 + deriv%M_t(ii,ll+1)*p(ii,jj+1)
        end do
        deriv%vtempt(jj  ,ll  ,1) = sumx00
        deriv%vtempt(jj  ,ll+1,1) = sumx01
        deriv%vtempt(jj+1,ll  ,1) = sumx10
        deriv%vtempt(jj+1,ll+1,1) = sumx11
        deriv%vtempt(jj  ,ll  ,2) = sumy00
        deriv%vtempt(jj  ,ll+1,2) = sumy01
        deriv%vtempt(jj+1,ll  ,2) = sumy10
        deriv%vtempt(jj+1,ll+1,2) = sumy11
      end do
      end do

      do jj=1,np,2
      do ii=1,np,2
        sumx00=0.0d0
        sumx01=0.0d0
        sumx10=0.0d0
        sumx11=0.0d0
        sumy00=0.0d0
        sumy01=0.0d0
        sumy10=0.0d0
        sumy11=0.0d0
        do ll=1,np
          sumx00 = sumx00 +  deriv%M_t(ll,jj  )*deriv%vtempt(ll,ii  ,1)
          sumx01 = sumx01 +  deriv%M_t(ll,jj+1)*deriv%vtempt(ll,ii  ,1)
          sumx10 = sumx10 +  deriv%M_t(ll,jj  )*deriv%vtempt(ll,ii+1,1)
          sumx11 = sumx11 +  deriv%M_t(ll,jj+1)*deriv%vtempt(ll,ii+1,1)
          sumy00 = sumy00 +  deriv%Dpv(ll,jj  )*deriv%vtempt(ll,ii  ,2)
          sumy01 = sumy01 +  deriv%Dpv(ll,jj+1)*deriv%vtempt(ll,ii  ,2)
          sumy10 = sumy10 +  deriv%Dpv(ll,jj  )*deriv%vtempt(ll,ii+1,2)
          sumy11 = sumy11 +  deriv%Dpv(ll,jj+1)*deriv%vtempt(ll,ii+1,2)
        end do
        dp(ii  ,jj  ,1) = sumx00
        dp(ii  ,jj+1,1) = sumx01
        dp(ii+1,jj  ,1) = sumx10
        dp(ii+1,jj+1,1) = sumx11
        dp(ii  ,jj  ,2) = sumy00
        dp(ii  ,jj+1,2) = sumy01
        dp(ii+1,jj  ,2) = sumy10
        dp(ii+1,jj+1,2) = sumy11
      end do
      end do
    else
      do jj=1,np
      do ll=1,np
        sumx00=0.0d0
        sumy00=0.0d0
        do ii=1,np
          sumx00 = sumx00 + deriv%Dpv(ii,ll)*p(ii,jj)
          sumy00 = sumy00 + deriv%M_t(ii,ll)*p(ii,jj)
        end do
        deriv%vtempt(jj,ll,1) = sumx00
        deriv%vtempt(jj,ll,2) = sumy00
      end do
      end do

      do jj=1,np
      do ii=1,np
        sumx00=0.0d0
        sumy00=0.0d0
        do ll=1,np
          sumx00 = sumx00 +  deriv%M_t(ll,jj)*deriv%vtempt(ll,ii,1)
          sumy00 = sumy00 +  deriv%Dpv(ll,jj)*deriv%vtempt(ll,ii,2)
        end do
        dp(ii,jj,1) = sumx00
        dp(ii,jj,2) = sumy00
      end do
      end do
    endif

    ! End Function
    !-------------
    return
  end function gradient_str_stag
  !==================================================================


  !==================================================================
  function gradient_str_nonstag(s,deriv) result(ds)
    !  gradient_str_nonstag: Compute the *strong* gradient on the velocity grid
    !                        of a scalar field on the velocity grid
    !====================================================================
    !
    ! Passed Variables
    !---------------------
    real(kind=real_kind),intent(in):: s(np,np)
    type(derivative_t)             :: deriv
    real(kind=real_kind)           :: ds(np,np,2)
    !
    ! Local Values
    !----------------
    logical, parameter  :: UseUnroll = .TRUE.
    real(kind=real_kind):: dsdx00,dsdx01
    real(kind=real_kind):: dsdx10,dsdx11
    real(kind=real_kind):: dsdy00,dsdy01
    real(kind=real_kind):: dsdy10,dsdy11

    integer ii,jj,ll

#ifdef DEBUG
    print *, "gradient_str_nonstag"
!   write(17) np,s,deriv
#endif

    if((modulo(np,2).eq.0).and.(UseUnroll)) then
       do jj=1,np,2
       do ll=1,np,2
         dsdx00=0.0d0
         dsdx01=0.0d0
         dsdx10=0.0d0
         dsdx11=0.0d0
         dsdy00=0.0d0
         dsdy01=0.0d0
         dsdy10=0.0d0
         dsdy11=0.0d0
         do ii=1,np
           dsdx00 = dsdx00 + deriv%Dvv(ii,ll  )*s(ii,jj  )
           dsdx01 = dsdx01 + deriv%Dvv(ii,ll+1)*s(ii,jj  )
           dsdx10 = dsdx10 + deriv%Dvv(ii,ll  )*s(ii,jj+1)
           dsdx11 = dsdx11 + deriv%Dvv(ii,ll+1)*s(ii,jj+1)
           dsdy00 = dsdy00 + deriv%Dvv(ii,ll  )*s(jj  ,ii)
           dsdy01 = dsdy01 + deriv%Dvv(ii,ll+1)*s(jj  ,ii)
           dsdy10 = dsdy10 + deriv%Dvv(ii,ll  )*s(jj+1,ii)
           dsdy11 = dsdy11 + deriv%Dvv(ii,ll+1)*s(jj+1,ii)
         end do
#ifdef DEBUG
         if(jj.eq.3.and.ll.eq.1) then
           print *, dsdx00
         endif
#endif
         ds(ll  ,jj  ,1) = dsdx00
         ds(ll+1,jj  ,1) = dsdx01
         ds(ll  ,jj+1,1) = dsdx10
         ds(ll+1,jj+1,1) = dsdx11
         ds(jj  ,ll  ,2) = dsdy00
         ds(jj  ,ll+1,2) = dsdy01
         ds(jj+1,ll  ,2) = dsdy10
         ds(jj+1,ll+1,2) = dsdy11
       end do
       end do
     else
       do jj=1,np
       do ll=1,np
         dsdx00=0.0d0
         dsdy00=0.0d0
         do ii=1,np
           dsdx00 = dsdx00 + deriv%Dvv(ii,ll)*s(ii,jj)
           dsdy00 = dsdy00 + deriv%Dvv(ii,ll)*s(jj,ii)
         end do
         ds(ll,jj,1) = dsdx00
         ds(jj,ll,2) = dsdy00
       end do
       end do
     endif

    ! End Function
    !-------------
    return
  end function gradient_str_nonstag
  !==================================================================


  !==================================================================
  function vorticity(v,deriv) result(vort)
    !  vorticity: Compute the vorticity of the velocity field on the
    !             velocity grid
    !=================================================================
    !
    ! Passed Variables
    !------------------
    real(kind=real_kind),intent(in):: v(np,np,2)
    type(derivative_t)             :: deriv
    real(kind=real_kind)            :: vort(np,np)
    !
    ! Local Values
    !--------------
    logical, parameter  :: UseUnroll = .TRUE.
    real(kind=real_kind):: dvdx00,dvdx01
    real(kind=real_kind):: dvdx10,dvdx11
    real(kind=real_kind):: dudy00,dudy01
    real(kind=real_kind):: dudy10,dudy11

    integer ii,jj,ll
    
    if((modulo(np,2) == 0).and.(UseUnroll)) then 
      do jj=1,np,2
      do ll=1,np,2
        dudy00=0.0d0
        dudy01=0.0d0
        dudy10=0.0d0
        dudy11=0.0d0
        dvdx00=0.0d0
        dvdx01=0.0d0
        dvdx10=0.0d0
        dvdx11=0.0d0
        do ii=1,np
          dvdx00 = dvdx00 + deriv%Dvv(ii,ll  )*v(ii,jj  ,2)
          dvdx01 = dvdx01 + deriv%Dvv(ii,ll+1)*v(ii,jj  ,2)
          dvdx10 = dvdx10 + deriv%Dvv(ii,ll  )*v(ii,jj+1,2)
          dvdx11 = dvdx11 + deriv%Dvv(ii,ll+1)*v(ii,jj+1,2)
          dudy00 = dudy00 + deriv%Dvv(ii,ll  )*v(jj  ,ii,1)
          dudy01 = dudy01 + deriv%Dvv(ii,ll+1)*v(jj  ,ii,1)
          dudy10 = dudy10 + deriv%Dvv(ii,ll  )*v(jj+1,ii,1)
          dudy11 = dudy11 + deriv%Dvv(ii,ll+1)*v(jj+1,ii,1)
        end do
        vort(ll  ,jj  ) = dvdx00
        vort(ll+1,jj  ) = dvdx01
        vort(ll  ,jj+1) = dvdx10
        vort(ll+1,jj+1) = dvdx11
        deriv%vvtemp(jj  ,ll  ) = dudy00
        deriv%vvtemp(jj  ,ll+1) = dudy01
        deriv%vvtemp(jj+1,ll  ) = dudy10
        deriv%vvtemp(jj+1,ll+1) = dudy11
      end do
      end do
    else
      do jj=1,np
      do ll=1,np
        dudy00=0.0d0
        dvdx00=0.0d0
        do ii=1,np
          dvdx00 = dvdx00 + deriv%Dvv(ii,ll)*v(ii,jj,2)
          dudy00 = dudy00 + deriv%Dvv(ii,ll)*v(jj,ii,1)
        end do
        vort(ll,jj) = dvdx00
        deriv%vvtemp(jj,ll) = dudy00
      end do
      end do
    endif

    do jj=1,np
    do ii=1,np
      vort(ii,jj)=vort(ii,jj)-deriv%vvtemp(ii,jj)
    end do
    end do

    ! End Function
    !-------------
    return
  end function vorticity
  !==================================================================


  !==================================================================
  function interpolate_gll2fvm_points(v,deriv) result(p)
    !  interpolate_gll2fvm_points: shape funtion interpolation from data on 
    !                              GLL grid to cellcenters on physics grid
    !                              Author: Christoph Erath
    !==========================================================================
    ! 
    ! Passed Variables
    !------------------
    real(kind=real_kind),intent(in):: v(np,np)
    type(derivative_t)             :: deriv
    real(kind=real_kind)           :: p(nc,nc)
    !
    ! Local Values
    !--------------
    real(kind=real_kind):: sumx00,sumx01
    real(kind=real_kind):: sumx10,sumx11
    real(kind=real_kind):: vtemp(np,nc)

    integer ii,jj,ll

    do jj=1,np
    do ll=1,nc
      sumx00=0.0d0
      do ii=1,np
        sumx00 = sumx00 + deriv%Cfvm(ii,ll)*v(ii,jj)
      end do
      vtemp(jj,ll) = sumx00
    end do
    end do

    do jj=1,nc
    do ii=1,nc
      sumx00=0.0d0
      do ll=1,np
        sumx00 = sumx00 + deriv%Cfvm(ll,jj)*vtemp(ll,ii)
      end do
      p(ii,jj) = sumx00
    end do
    end do

    ! End Function
    !-------------
    return
  end function interpolate_gll2fvm_points
  !==================================================================


  !==================================================================
  function interpolate_gll2spelt_points(v,deriv) result(p)
    !  interpolate_gll2spelt_points: shape function interpolation from data 
    !                                on GLL grid the spelt grid
    !                                Author: Christoph Erath
    !=====================================================================
    !
    ! Passed Variables
    !-------------------
    real(kind=real_kind),intent(in):: v(np,np)
    type(derivative_t)             :: deriv
    real(kind=real_kind)           :: p(nep,nep)
    !
    ! Local Values
    !-------------
    real(kind=real_kind):: sumx00,sumx01
    real(kind=real_kind):: sumx10,sumx11
    real(kind=real_kind):: vtemp(np,nep)

    integer ii,jj,ll

    do jj=1,np
    do ll=1,nep
      sumx00=0.0d0
      do ii=1,np
        sumx00 = sumx00 + deriv%Sfvm(ii,ll)*v(ii,jj)
      end do
      vtemp(jj,ll) = sumx00
    end do
    end do

    do jj=1,nep
    do ii=1,nep
      sumx00=0.0d0
      do ll=1,np
        sumx00 = sumx00 + deriv%Sfvm(ll,jj)*vtemp(ll,ii)
      end do
      p(ii,jj) = sumx00
    end do
    end do

    ! End Function
    !-------------
    return
  end function interpolate_gll2spelt_points
  !==================================================================


  !==================================================================
  function interpolate_gll2fvm_corners(v,deriv) result(p)
    !  interpolate_gll2fvm_corners: shape funtion interpolation from data 
    !                               on GLL grid to physics grid
    !
    !====================================================================
    !
    ! Passed Variables
    !------------------
    real(kind=real_kind),intent(in):: v(np,np)
    type(derivative_t)             :: deriv
    real(kind=real_kind)           :: p(nc+1,nc+1)
    !
    ! Local Values
    !--------------
    real(kind=real_kind):: sumx00,sumx01
    real(kind=real_kind):: sumx10,sumx11
    real(kind=real_kind):: vtemp(np,nc+1)

    integer ii,jj,ll

    do jj=1,np
    do ll=1,nc+1
      sumx00=0.0d0
      do ii=1,np
        sumx00 = sumx00 + deriv%Mfvm(ii,ll)*v(ii,jj)
      end do
      vtemp(jj,ll) = sumx00
    end do
    end do

    do jj=1,nc+1
    do ii=1,nc+1
      sumx00=0.0d0
      do ll=1,np
        sumx00 = sumx00 + deriv%Mfvm(ll,jj)*vtemp(ll,ii)
      end do
      p(ii,jj) = sumx00
    end do
    end do

    ! End Function
    !-------------
    return
  end function interpolate_gll2fvm_corners
  !==================================================================


  !==================================================================
  function remap_phys2gll(pin,nphys) result(pout)
    !  remap_phys2gll: interpolate to an equally spaced (in reference element 
    !                  coordinate system) "physics" grid to the GLL grid
    !
    !                  1st order, monotone, conservative
    !=======================================================================
    !
    ! Passed Variables
    !------------------
    real(kind=real_kind),intent(in):: pin(nphys*nphys)
    integer                        :: nphys
    real(kind=real_kind)           :: pout(np,np)
    !
    ! Local Values
    !---------------
    integer,save                     :: nphys_init=0
    integer,save                     :: nintersect
    real(kind=real_kind),save,pointer:: acell(:)  ! arrivial cell index of i'th intersection
    real(kind=real_kind),save,pointer:: dcell(:)  ! departure cell index of i'th intersection
    real(kind=real_kind),save,pointer:: delta(:)  ! length of i'th intersection
    real(kind=real_kind),save,pointer:: delta_a(:) ! length of arrival cells

    logical                   :: found
    real(kind=real_kind)      :: tol=1e-13
    real(kind=real_kind)      :: weight,x1,x2,dx
    real(kind=longdouble_kind):: gll_edges(np+1),phys_edges(nphys+1)
    type(quadrature_t)        :: gll_pts

    integer in_i,in_j,ia,ja,id,jd,count,ii,jj

    ! setup (most be done on masterthread only) since all data is static
    !---------------------------------------------------------------------
#if (! defined ELEMENT_OPENMP)
!OMP MASTER
#endif
    if(nphys_init/=nphys) then

      ! find number of intersections
      !------------------------------
      nintersect = np+nphys-1  ! max number of possible intersections
      allocate(acell(nintersect))
      allocate(dcell(nintersect))
      allocate(delta(nintersect))
      allocate(delta_a(np))

      ! compute phys grid cell edges on [-1,1]
      !---------------------------------------
      do ii=1,nphys+1
        dx = 2d0/nphys
        phys_edges(ii)=-1 + (ii-1)*dx
      end do

      ! compute GLL cell edges on [-1,1]
      !---------------------------------
      gll_pts = gausslobatto(np)
      gll_edges(1)=-1
      do ii=2,np
        gll_edges(ii) = gll_edges(ii-1) + gll_pts%weights(ii-1)
      end do
      gll_edges(np+1)=1
      delta_a=gll_pts%weights
      deallocate(gll_pts%points )
      deallocate(gll_pts%weights)

      count=0
      x1=-1
      do while ( abs(x1-1) > tol )
        ! find point x2 closet to x1 and x2>x1:
        !---------------------------------------
        x2=1.1
        do ia=2,np+1
          if(gll_edges(ia)>x1) then
            if((gll_edges(ia)-x1) < (x2-x1)) then
              x2=gll_edges(ia)
            endif
          endif
        end do
        do id=2,nphys+1
          if(phys_edges(id)>x1) then
            if((phys_edges(id)-x1) < (x2-x1)) then
              x2=phys_edges(id)
            endif
          endif
        end do
        if(x2>1+tol) call endrun('ERROR: did not find next intersection point')
             
        count=count+1
        delta(count)=x2-x1
          
        found=.false.
        do ia=1,np
          if((gll_edges(ia) <= (x1+tol)).and.((x2-tol) <= gll_edges(ia+1))) then
            found=.true.
            acell(count)=ia
          endif
        end do
        if(.not. found) call endrun('ERROR: interval search problem')
       
        found=.false.
        do id=1,nphys
          if((phys_edges(id) <= (x1+tol)).and.((x2-tol) <= phys_edges(id+1))) then
            found=.true.
            dcell(count)=id
          endif
        end do
        if(.not. found) call endrun('ERROR: interval search problem')
        x1=x2
      end do ! while ( abs(x1-1) > tol )

      if(count>nintersect) call endrun('ERROR: nintersect was too small')
      nintersect=count
#if 0
      print *,'gll->phys conservative monotone remap algorithm:'
      print *,'np,nphys,nintersect',np,nphys,nintersect
      print *,'i   [x1,x2]   [acell]   [dcell]'
      x1=-1
      do in_i=1,nintersect
        ia=acell(in_i)
        id=dcell(in_i)
        write(*,'(i3,a,2f10.6,a,a,2f10.6,a,a,2f10.6,a)') in_i,&
                '[',x1,x1+delta(in_i),']',                    &
                '[',gll_edges(ia),gll_edges(ia+1),']',        &
                '[',phys_edges(id),phys_edges(id+1),']'
        x1=x1+delta(in_i)
      end do

      pout=0
      do in_i = 1,nintersect
      do in_j = 1,nintersect
        ia = acell(in_i)
        ja = acell(in_j)
        id = dcell(in_i)
        jd = dcell(in_j)
        weight = (  delta(in_i)*delta(in_j) ) / ( delta_a(ia)*delta_a(ja))
        pout(ia,ja) = pout(ia,ja) + weight
      end do
      end do
      print *,'sum of weights: ',pout(:,:)
      call endrun(__FILE__)
#endif
      ! "nphys" initializaton done.
      !---------------------------
      nphys_init=nphys
    endif

#if (! defined ELEMENT_OPENMP)
    !OMP END MASTER
    !OMP BARRIER
#endif

    pout=0
    do in_i = 1,nintersect
    do in_j = 1,nintersect
      ia = acell(in_i)
      ja = acell(in_j)
      id = dcell(in_i)
      jd = dcell(in_j)

      ! mass in intersection region:  value*area_intersect
      ! value_arrival = value*area_intersect/area_arrival
      !----------------------------------------------------
      weight = (delta(in_i)*delta(in_j)) / (delta_a(ia)*delta_a(ja))

      ! accumulate contribution from each intersection region:
      !-------------------------------------------------------
      pout(ia,ja) = pout(ia,ja) + weight*pin(id+(jd-1)*nphys)
    end do
    end do
    
    ! End Function
    !-------------
    return
  end function remap_phys2gll
  !==================================================================

!----------------------------------------------------------------
    
  !==================================================================
  function gradient_sphere(s,deriv,Dinv) result(ds)
    !  gradient_sphere:
    !                    input s:  scalar
    !                    output  ds: spherical gradient of s, lat-lon coordinates
    !----------------------------------------------------------------
    !
    ! Passed variables
    !-----------------
    real(kind=real_kind), intent(in):: s(np,np)
    type(derivative_t)              :: deriv
    real(kind=real_kind),intent(in) :: Dinv(2,2,np,np)
    real(kind=real_kind)            :: ds(np,np,2)
    !
    ! Local Values
    !--------------
    real(kind=real_kind):: dsdx00
    real(kind=real_kind):: dsdy00
    real(kind=real_kind):: v1(np,np),v2(np,np)

    integer ii,jj,ll

    do jj=1,np
    do ll=1,np
      dsdx00=0.0d0
      dsdy00=0.0d0
      do ii=1,np
        dsdx00 = dsdx00 + deriv%Dvv(ii,ll)*s(ii,jj)
        dsdy00 = dsdy00 + deriv%Dvv(ii,ll)*s(jj,ii)
      end do
      v1(ll,jj) = dsdx00*rrearth
      v2(jj,ll) = dsdy00*rrearth
    end do
    end do

    ! convert covarient to latlon
    !-----------------------------
    do jj=1,np
    do ii=1,np
      ds(ii,jj,1)=Dinv(1,1,ii,jj)*v1(ii,jj) + Dinv(2,1,ii,jj)*v2(ii,jj)
      ds(ii,jj,2)=Dinv(1,2,ii,jj)*v1(ii,jj) + Dinv(2,2,ii,jj)*v2(ii,jj)
    end do
    end do

    ! End Function
    !-------------
    return
  end function gradient_sphere
  !==================================================================


  !==================================================================
  function curl_sphere_wk_testcov(s,deriv,elem) result(ds)
    !  curl_sphere_wk_testcov:
    !
    !   integrated-by-parts gradient, w.r.t. COVARIANT test functions
    !   input s:  scalar  (assumed to be s*khat)
    !   output  ds: weak curl, lat/lon coordinates
    !   
    ! starting with: 
    !   PHIcov1 = (PHI,0)  covariant vector 
    !   PHIcov2 = (0,PHI)  covariant vector 
    !
    !   ds1 = integral[ PHIcov1 dot curl(s*khat) ] 
    !   ds2 = integral[ PHIcov2 dot curl(s*khat) ] 
    ! integrate by parts: 
    !   ds1 = integral[ vor(PHIcov1) * s ]       
    !   ds2 = integral[ vor(PHIcov1) * s ]
    !
    !     PHIcov1 = (PHI^mn,0)   
    !     PHIcov2 = (0,PHI^mn)
    !  vorticity() acts on covariant vectors:
    !   ds1 = sum wij g  s_ij 1/g (  (PHIcov1_2)_x  - (PHIcov1_1)_y ) 
    !       = -sum wij s_ij  d/dy (PHI^mn )
    ! for d/dy component, only sum over i=m
    !       = -sum  w_mj s_mj   d( PHI^n)(j)
    !           j
    !
    !   ds2 = sum wij g  s_ij 1/g (  (PHIcov2_2)_x  - (PHIcov2_1)_y ) 
    !       = +sum wij s_ij  d/dx (PHI^mn )
    ! for d/dx component, only sum over j=n
    !       = +sum  w_in s_in  d( PHI^m)(i)
    !           i
    !===========================================================================
    !
    ! Passed Varibales
    !-------------------
    real(kind=real_kind),intent(in):: s(np,np)
    type(derivative_t)             :: deriv
    type(element_t)                :: elem
    real(kind=real_kind)           :: ds(np,np,2)
    !
    ! Local Values
    !---------------
    real(kind=real_kind)::  dscontra(np,np,2)

    integer ii,jj,ll,mm,nn

    dscontra=0
    do nn=1,np
    do mm=1,np
    do jj=1,np
      ! phi(n)_y  sum over second index,    1st index fixed at m
      ! phi(m)_x  sum over  first index, second index fixed at n
      !-------------------------------------------------------
      dscontra(mm,nn,1)=dscontra(mm,nn,1)-(elem%mp(mm,jj)*s(mm,jj)*deriv%Dvv(nn,jj))*rrearth
      dscontra(mm,nn,2)=dscontra(mm,nn,2)+(elem%mp(jj,nn)*s(jj,nn)*deriv%Dvv(mm,jj))*rrearth
    end do
    end do
    end do

    ! convert contra -> latlon 
    !--------------------------
    do jj=1,np
    do ii=1,np
      ds(ii,jj,1)=( elem%D(1,1,ii,jj)*dscontra(ii,jj,1) &
                   +elem%D(1,2,ii,jj)*dscontra(ii,jj,2) )
      ds(ii,jj,2)=( elem%D(2,1,ii,jj)*dscontra(ii,jj,1) &
                  + elem%D(2,2,ii,jj)*dscontra(ii,jj,2) )
    end do
    end do

    ! End Function
    !-------------
    return
  end function curl_sphere_wk_testcov
  !==================================================================


  !==================================================================
  function gradient_sphere_wk_testcov(s,deriv,elem) result(ds)
    !  gradient_sphere_wk_testcov:
    !
    !   integrated-by-parts gradient, w.r.t. COVARIANT test functions
    !   input s:  scalar
    !   output  ds: weak gradient, lat/lon coordinates
    !   ds = - integral[ div(PHIcov) s ]
    !
    !     PHIcov1 = (PHI^mn,0)   
    !     PHIcov2 = (0,PHI^mn)
    !   div() acts on contra components, so convert test function to contra: 
    !     PHIcontra1 =  metinv PHIcov1  = (a^mn,b^mn)*PHI^mn   
    !                                     a = metinv(1,1)  b=metinv(2,1)
    !
    !   ds1 = sum wij g  s_ij 1/g ( g a PHI^mn)_x  + ( g b PHI^mn)_y ) 
    !       = sum  wij s_ij  ag(m,n)  d/dx( PHI^mn ) + bg(m,n) d/dy( PHI^mn)
    !          i,j 
    ! for d/dx component, only sum over j=n
    !       = sum  w_in s_in  ag(m,n)  d( PHI^m)(i)
    !          i
    ! for d/dy component, only sum over i=m
    !       = sum  w_mj s_mj  bg(m,n)  d( PHI^n)(j)
    !          j
    !  
    !
    ! This formula is identical to gradient_sphere_wk_testcontra, except that
    !    g(m,n) is replaced by a(m,n)*g(m,n)   
    !  and we have two terms for each componet of ds 
    !
    !========================================================================
    !
    ! Passed Variables
    !------------------
    real(kind=real_kind),intent(in):: s(np,np)
    type(derivative_t)             :: deriv
    type(element_t)                :: elem
    real(kind=real_kind)           :: ds(np,np,2)
    !
    ! Local Values
    !--------------
    real(kind=real_kind)::  dscontra(np,np,2)

    integer ii,jj,ll,mm,nn

    dscontra=0
    do nn=1,np
    do mm=1,np
    do jj=1,np
      dscontra(mm,nn,1)=dscontra(mm,nn,1)                                              &
                       -rrearth*( ( elem%mp        (jj,nn)                             &
                                   *elem%metinv(1,1,mm,nn)                             &
                                   *elem%metdet    (mm,nn)*s(jj,nn)*deriv%Dvv(mm,jj) ) &
                                 +( elem%mp        (mm,jj)                             &
                                   *elem%metinv(2,1,mm,nn)                             &
                                   *elem%metdet    (mm,nn)*s(mm,jj)*deriv%Dvv(nn,jj) ) ) 
      dscontra(mm,nn,2)=dscontra(mm,nn,2)                                              &
                       -rrearth*( ( elem%mp        (jj,nn)                             &
                                   *elem%metinv(1,2,mm,nn)                             &
                                   *elem%metdet    (mm,nn)*s(jj,nn)*deriv%Dvv(mm,jj) ) &
                                 +(elem%mp         (mm,jj)                             &
                                   *elem%metinv(2,2,mm,nn)                             &
                                   *elem%metdet    (mm,nn)*s(mm,jj)*deriv%Dvv(nn,jj) ) ) 
    end do
    end do
    end do

    ! convert contra -> latlon 
    !--------------------------
    do jj=1,np
    do ii=1,np
      ds(ii,jj,1)=( elem%D(1,1,ii,jj)*dscontra(ii,jj,1) &
                   +elem%D(1,2,ii,jj)*dscontra(ii,jj,2) )
      ds(ii,jj,2)=( elem%D(2,1,ii,jj)*dscontra(ii,jj,1) &
                   +elem%D(2,2,ii,jj)*dscontra(ii,jj,2) )
    end do
    end do

    ! End Function
    !-------------
    return
  end function gradient_sphere_wk_testcov
  !==================================================================


  !==================================================================
  function gradient_sphere_wk_testcontra(s,deriv,elem) result(ds)
    !  gradient_sphere_wk_testcontra:
    !
    !   integrated-by-parts gradient, w.r.t. CONTRA test functions
    !   input s:  scalar
    !   output  ds: weak gradient, lat/lon coordinates
    !
    !   integral[ div(phivec) s ] = sum  spheremp()* divergence_sphere(phivec) *s
    !   ds1 = above formual with phivec=(PHI,0) in CONTRA coordinates
    !   ds2 = above formual with phivec=(0,PHI) in CONTRA coordinates
    !   
    ! PHI = (phi,0)
    !   s1 =  sum w_ij s_ij g_ij 1/g_ij ( g_ij PHI^mn )x  
    !      =  sum w_ij s_ij g_mn dx(PHI^mn)_ij 
    !         ij
    ! because x derivative is zero for j<>n, only have to sum over j=n
    !   s1(m,n)  =  sum w_i,n g_mn dx(PHI^m)_i,n s_i,n
    !                i
    !================================================================
    !
    ! Passed Variables
    !--------------------
    real(kind=real_kind),intent(in):: s(np,np)
    type(derivative_t)             :: deriv
    type(element_t)                :: elem
    real(kind=real_kind)           :: ds(np,np,2)
    !
    ! Local Values
    !-------------
    real(kind=real_kind):: dscov(np,np,2)

    integer ii,jj,ll,mm,nn

    ! debug: 
    real(kind=real_kind):: vcontra(np,np,2)
    real(kind=real_kind):: v(np,np,2)
    real(kind=real_kind):: div(np,np)

    dscov=0
    do nn=1,np
    do mm=1,np
    do jj=1,np
      ! phi(m)_x  sum over  first index, second index fixed at n
      ! phi(n)_y  sum over second index,    1st index fixed at m
      !---------------------------------------------------------------
      dscov(mm,nn,1)=dscov(mm,nn,1)                                         &
                    -rrearth*( elem%mp    (jj,nn)                           &
                              *elem%metdet(mm,nn)*s(jj,nn)*deriv%Dvv(mm,jj) )
      dscov(mm,nn,2)=dscov(mm,nn,2)                                         &
                    -rrearth*( elem%mp    (mm,jj)                           &
                              *elem%metdet(mm,nn)*s(mm,jj)*deriv%Dvv(nn,jj) )
    end do
    end do
    end do

#if 0
    ! slow form, for debugging
    !------------------------
    do mm=1,np
    do nn=1,np
      vcontra=0
      vcontra(mm,nn,1)=1

      ! contra->latlon:
      !----------------
      v(:,:,1)=(elem%D(1,1,:,:)*vcontra(:,:,1) + elem%D(1,2,:,:)*vcontra(:,:,2))
      v(:,:,2)=(elem%D(2,1,:,:)*vcontra(:,:,1) + elem%D(2,2,:,:)*vcontra(:,:,2))

      ! compute div(metdet phivec) * s
      !--------------------------------
      div = divergence_sphere(v,deriv,elem)

      ! compute integral[ div(phi) * s ]
      !---------------------------------
      ds(mm,nn,1)=0
      do ii=1,np
      do jj=1,np
        ds(mm,nn,1)=ds(mm,nn,1) + div(ii,jj)*s(ii,jj)*elem%spheremp(ii,jj)
      end do
      end do

      vcontra=0
      vcontra(mm,nn,2)=1

      ! contra->latlon:
      !----------------
      v(:,:,1)=(elem%D(1,1,:,:)*vcontra(:,:,1) + elem%D(1,2,:,:)*vcontra(:,:,2))
      v(:,:,2)=(elem%D(2,1,:,:)*vcontra(:,:,1) + elem%D(2,2,:,:)*vcontra(:,:,2))

      ! compute div(metdet phivec) * s
      !----------------------------------
      div = divergence_sphere(v,deriv,elem)

      ! compute integral[ div(phi) * s ]
      !-----------------------------------
      ds(mm,nn,2)=0
      do ii=1,np
      do jj=1,np
        ds(mm,nn,2)=ds(mm,nn,2) + div(ii,jj)*s(ii,jj)*elem%spheremp(ii,jj)
      end do
      end do
    end do
    end do

    ! change sign 
    ds=-ds
    print *,'ds,dscov:1 ',ds(1,1,1),dscov(1,1,1),ds(1,1,1)/dscov(1,1,1)
    print *,'ds,dscov:2 ',ds(1,1,2),dscov(1,1,2),ds(1,1,2)/dscov(1,1,2)

    dscov=ds
#endif

    ! convert covariant -> latlon 
    !-----------------------------
    ds(:,:,1)=elem%Dinv(1,1,:,:)*dscov(:,:,1) + elem%Dinv(2,1,:,:)*dscov(:,:,2)
    ds(:,:,2)=elem%Dinv(1,2,:,:)*dscov(:,:,1) + elem%Dinv(2,2,:,:)*dscov(:,:,2)

    ! End Function
    !-------------
    return
  end function gradient_sphere_wk_testcontra
  !==================================================================


  !==================================================================
  function ugradv_sphere(u,v,deriv,elem) result(ugradv)
    ! ugradv_sphere:
    !
    !   input:  vectors u and v  (latlon coordinates)
    !   output: vector  [ u dot grad ] v  (latlon coordinates)
    !================================================================
    !
    ! Passed Variables
    !------------------
    real(kind=real_kind),intent(in):: u(np,np,2)
    real(kind=real_kind),intent(in):: v(np,np,2)
    type(derivative_t)             :: deriv
    type(element_t)                :: elem
    real(kind=real_kind)           :: ugradv(np,np,2)
    !
    ! Local Values
    !---------------
    real(kind=real_kind):: dum_cart(np,np,3)
    integer             :: component

    ! latlon -> cartesian
    !---------------------
    do component=1,3
       ! Summing along the third dimension is a sum over components for each point.
       ! (This is just a faster way of doing a dot product for each grid point,
       ! since reindexing the inputs to use the intrinsic effectively would be
       ! just asking for trouble.)
       !-----------------------------------------------------------------------------
       dum_cart(:,:,component)=sum( elem%vec_sphere2cart(:,:,component,:)*v(:,:,:) ,3)
    end do

    ! Do ugradv on the cartesian components.
    !---------------------------------------
    do component=1,3
       ! Dot u with the gradient of each component
       !-------------------------------------------
       dum_cart(:,:,component) = sum( u(:,:,:) * &
               gradient_sphere(dum_cart(:,:,component),deriv,elem%Dinv) ,3)
    end do

    ! cartesian -> latlon
    !-------------------------
    do component=1,2
       ! vec_sphere2cart is its own pseudoinverse.
       !--------------------------------------------
       ugradv(:,:,component)=sum( dum_cart(:,:,:)*elem%vec_sphere2cart(:,:,:,component) ,3)
    end do

    ! End Function
    !-------------
    return
  end function ugradv_sphere
  !==================================================================


  !==================================================================
  function curl_sphere(s,deriv,elem) result(ds)
    ! curl_sphere:
    !
    !   input s:  scalar  (assumed to be  s khat)
    !   output  curl(s khat) vector in lat-lon coordinates
    ! 
    !   This subroutine can be used to compute divergence free velocity fields,
    !   since div(ds)=0
    !
    !    first compute:  
    !    curl(s khat) = (1/jacobian) ( ds/dy, -ds/dx ) in contra-variant coordinates
    !    then map to lat-lon
    !================================================================================
    !
    ! Passed Variables
    !-------------------
    real(kind=real_kind),intent(in):: s(np,np)
    type(derivative_t)             :: deriv
    type(element_t)                :: elem
    real(kind=real_kind) :: ds(np,np,2)
    ! 
    ! Local Values
    !---------------
    real(kind=real_kind):: dsdx00
    real(kind=real_kind):: dsdy00
    real(kind=real_kind):: v1(np,np),v2(np,np)
    
    integer ii,jj,ll

    do jj=1,np
    do ll=1,np
      dsdx00=0.0d0
      dsdy00=0.0d0
      do ii=1,np
        dsdx00 = dsdx00 + deriv%Dvv(ii,ll)*s(ii,jj)
        dsdy00 = dsdy00 + deriv%Dvv(ii,ll)*s(jj,ii)
      end do
      v2(ll,jj) = -dsdx00*rrearth
      v1(jj,ll) =  dsdy00*rrearth
    end do
    end do

    ! convert contra -> latlon *and* divide by jacobian
    !--------------------------------------------------
    do jj=1,np
    do ii=1,np
      ds(ii,jj,1)=( elem%D(1,1,ii,jj)*v1(ii,jj)                   &
                   +elem%D(1,2,ii,jj)*v2(ii,jj))/elem%metdet(ii,jj)
      ds(ii,jj,2)=( elem%D(2,1,ii,jj)*v1(ii,jj)                   &
                   +elem%D(2,2,ii,jj)*v2(ii,jj))/elem%metdet(ii,jj)
    end do
    end do
 
    ! End Function
    !-------------
    return
  end function curl_sphere
  !==================================================================

!--------------------------------------------------------------------------

  !==================================================================
  function divergence_sphere_wk(v,deriv,elem) result(div)
    ! divergence_sphere_wk:
    !
    !   input:  v = velocity in lat-lon coordinates
    !   ouput:  div(v)  spherical divergence of v, integrated by parts
    !
    !   Computes  -< grad(psi) dot v > 
    !   (the integrated by parts version of < psi div(v) > )
    !
    !   note: after DSS, divergence_sphere() and divergence_sphere_wk() 
    !   are identical to roundoff, as theory predicts.
    !=========================================================================
    !
    ! Passed Variables
    !------------------
    real(kind=real_kind),intent(in):: v(np,np,2)  ! in lat-lon coordinates
    type(derivative_t)             :: deriv
    type(element_t)                :: elem
    real(kind=real_kind)           :: div(np,np)
    !
    ! Local Values
    !-----------------
    real(kind=real_kind)::  vtemp(np,np,2)
    real(kind=real_kind):: ggtemp(np,np,2)
    real(kind=real_kind)::  gtemp(np,np,2)
    real(kind=real_kind)::    psi(np,np)
    real(kind=real_kind):: xtmp

    integer ii,jj,mm,nn

    ! latlon- > contra
    !--------------------
    do jj=1,np
    do ii=1,np
      vtemp(ii,jj,1)=(elem%Dinv(1,1,ii,jj)*v(ii,jj,1) + elem%Dinv(1,2,ii,jj)*v(ii,jj,2))
      vtemp(ii,jj,2)=(elem%Dinv(2,1,ii,jj)*v(ii,jj,1) + elem%Dinv(2,2,ii,jj)*v(ii,jj,2))
    end do
    end do

    do nn=1,np
    do mm=1,np
      div(mm,nn)=0
      do jj=1,np
        div(mm,nn)=div(mm,nn)                                                     &
                  -rrearth*( elem%spheremp(jj,nn)*vtemp(jj,nn,1)*deriv%Dvv(mm,jj) &
                            +elem%spheremp(mm,jj)*vtemp(mm,jj,2)*deriv%Dvv(nn,jj) ) 
      end do
#if 0
      ! debug the above formula using the N^4 slow formulation:
      !---------------------------------------------------------
      psi=0
      psi(mm,nn)=1
      ggtemp=gradient_sphere(psi,deriv,elem%Dinv)

      ! latlon -> covarient
      !--------------------
      do jj=1,np
      do ii=1,np
        gtemp(ii,jj,1)=( elem%D(1,1,ii,jj)*ggtemp(ii,jj,1) &
                        +elem%D(2,1,ii,jj)*ggtemp(ii,jj,2) )
        gtemp(ii,jj,2)=( elem%D(1,2,ii,jj)*ggtemp(ii,jj,1) &
                        +elem%D(2,2,ii,jj)*ggtemp(ii,jj,2) )
      end do
      end do

      ! grad(psi) dot v:
      !-----------------
      xtmp=0
      do jj=1,np
      do ii=1,np
        xtmp=xtmp-elem%spheremv(ii,jj)*( vtemp(ii,jj,1)*gtemp(ii,jj,1) &
                                        +vtemp(ii,jj,2)*gtemp(ii,jj,2) )
      end do
      end do
      if(abs(xtmp-div(mm,nn)) > 3e-17) then
        print *,mm,nn,xtmp,div(mm,nn),xtmp-div(mm,nn)
      endif
#endif          
    end do
    end do
    
    ! End Function
    !-------------
    return
  end function divergence_sphere_wk
  !==================================================================


  !==================================================================
  function element_boundary_integral(v,deriv,elem) result(result)
    ! element_boundary_integral:
    !
    !   input:  v = velocity in lat-lon coordinates
    !   ouput:  result(i,j) = contour integral of PHI_ij * v dot normal
    !           where PHI_ij = cardinal function at i,j GLL point 
    !
    !   this routine is used just to check spectral element integration 
    !   by parts identities
    !=========================================================================
    ! 
    ! Passed Variables
    !-------------------
    real(kind=real_kind),intent(in):: v(np,np,2)  ! in lat-lon coordinates
    type(derivative_t)             :: deriv
    type(element_t)                :: elem
    real(kind=real_kind)           :: result(np,np)
    !
    ! Local Values
    !-----------------
    real(kind=real_kind):: ucontra(np,np,2)  ! in lat-lon coordinates

    integer ii,jj

    ! latlon->contra
    !----------------
    do jj=1,np
    do ii=1,np
      ucontra(ii,jj,1)=( elem%Dinv(1,1,ii,jj)*v(ii,jj,1) &
                        +elem%Dinv(1,2,ii,jj)*v(ii,jj,2) )
      ucontra(ii,jj,2)=( elem%Dinv(2,1,ii,jj)*v(ii,jj,1) &
                        +elem%Dinv(2,2,ii,jj)*v(ii,jj,2) )
    end do
    end do

    ! note: GLL weights  weight(i) = Mvv_twt(i,i)
    !---------------------------------------------
    result=0
    jj=1
    do ii=1,np
      result(ii,jj)=result(ii,jj)-deriv%Mvv_twt(ii,ii)*elem%metdet(ii,jj) &
                                                          *ucontra(ii,jj,2)*rrearth
    end do
    
    jj=np
    do ii=1,np
      result(ii,jj)=result(ii,jj)+deriv%Mvv_twt(ii,ii)*elem%metdet(ii,jj) &
                                                          *ucontra(ii,jj,2)*rrearth
    end do
    
    ii=1
    do jj=1,np
      result(ii,jj)=result(ii,jj)-deriv%Mvv_twt(jj,jj)*elem%metdet(ii,jj) &
                                                          *ucontra(ii,jj,1)*rrearth
    end do
    
    ii=np
    do jj=1,np
       result(ii,jj)=result(ii,jj)+deriv%Mvv_twt(jj,jj)*elem%metdet(ii,jj) &
                                                           *ucontra(ii,jj,1)*rrearth
    end do

    ! End Function
    !-------------
    return
  end function element_boundary_integral
  !==================================================================


  !==================================================================
  function edge_flux_u_cg( v,p,pedges, deriv, elem, u_is_contra) result(result)
    ! edge_flux_u_cg:
    !
    !
    !   input:  v = velocity in contra or lat-lon coordinates (CONTINUIOUS)
    !           p      = scalar on this element
    !           pedges = scalar edge data from neighbor elements
    !
    !   ouput:  result(i,j) = contour integral of PHI_ij * pstar * v dot normal
    !           where PHI_ij = cardinal function at i,j GLL point 
    !           pstar = centered or other flux
    !========================================================================
    !
    ! Passed Variables
    !------------------
    real(kind=real_kind),intent(in):: v(np,np,2) 
    real(kind=real_kind),intent(in):: p(np,np) 
    real(kind=real_kind),intent(in):: pedges(0:np+1,0:np+1) 
    type(derivative_t)             :: deriv
    type(element_t)                :: elem
    logical                        :: u_is_contra
    real(kind=real_kind)           :: result(np,np)
    !
    ! Local Values
    !--------------
    real(kind=real_kind):: ucontra(np,np,2)  ! in lat-lon coordinates
    real(kind=real_kind):: flux,pstar

    integer ii,jj

    result=0
    if(u_is_contra) then
      ucontra=v
    else
      ! latlon->contra
      !---------------
      do jj=1,np
      do ii=1,np
        ucontra(ii,jj,1)=( elem%Dinv(1,1,ii,jj)*v(ii,jj,1) &
                          +elem%Dinv(1,2,ii,jj)*v(ii,jj,2) )
        ucontra(ii,jj,2)=( elem%Dinv(2,1,ii,jj)*v(ii,jj,1) &
                          +elem%Dinv(2,2,ii,jj)*v(ii,jj,2) )
      end do
      end do
    endif

#if 0
    ! centered
    !---------
    do ii=1,np
      jj=1
      pstar=(pedges(ii,0) + p(ii,jj))/2
      flux = -pstar*ucontra(ii,jj,2)*(deriv%Mvv_twt(ii,ii)*elem%metdet(ii,jj)*rrearth)
      result(ii,jj)=result(ii,jj)+flux
      jj=np
      pstar=(pedges(ii,np+1) + p(ii,jj))/2
      flux =  pstar*ucontra(ii,jj,2)*(deriv%Mvv_twt(ii,ii)*elem%metdet(ii,jj)*rrearth)
      result(ii,jj)=result(ii,jj)+flux
    end do
    
    do jj=1,np
      ii=1
      pstar=(pedges(0,jj) + p(ii,jj))/2
      flux = -pstar*ucontra(ii,jj,1)*(deriv%Mvv_twt(jj,jj)*elem%metdet(ii,jj)*rrearth)
      result(ii,jj)=result(ii,jj)+flux
      ii=np  
      pstar=(pedges(np+1,jj) + p(ii,jj))/2
      flux =  pstar*ucontra(ii,jj,1)*(deriv%Mvv_twt(jj,jj)*elem%metdet(ii,jj)*rrearth)
      result(ii,jj)=result(ii,jj)+flux
    end do
#else
    ! upwind
    !---------
    do ii=1,np
      jj=1
      pstar=p(ii,jj)
      if(ucontra(ii,jj,2)>0) pstar=pedges(ii,0)
      flux = -pstar*ucontra(ii,jj,2)*(deriv%Mvv_twt(ii,ii)*elem%metdet(ii,jj)*rrearth)
      result(ii,jj)=result(ii,jj)+flux
      jj=np
      pstar=p(ii,jj)
      if(ucontra(ii,jj,2)<0) pstar=pedges(ii,np+1)
      flux =  pstar*ucontra(ii,jj,2)*(deriv%Mvv_twt(ii,ii)*elem%metdet(ii,jj)*rrearth)
      result(ii,jj)=result(ii,jj)+flux
    end do
    
    do jj=1,np
      ii=1
      pstar=p(ii,jj)
      if(ucontra(ii,jj,1)>0) pstar=pedges(0,jj)
      flux = -pstar*ucontra(ii,jj,1)*(deriv%Mvv_twt(jj,jj)*elem%metdet(ii,jj)*rrearth)
      result(ii,jj)=result(ii,jj)+flux
      ii=np  
      pstar=p(ii,jj)
      if(ucontra(ii,jj,1)<0) pstar=pedges(np+1,jj)
      flux =  pstar*ucontra(ii,jj,1)*(deriv%Mvv_twt(jj,jj)*elem%metdet(ii,jj)*rrearth)
      result(ii,jj)=result(ii,jj)+flux
    end do
#endif    

    ! End Function
    !-------------
    return
  end function edge_flux_u_cg
  !==================================================================

    
  !==================================================================
  function vorticity_sphere(v,deriv,elem) result(vort)
    ! vorticity_sphere:
    !
    !   input:  v = velocity in lat-lon coordinates
    !   ouput:  spherical vorticity of v
    !================================================================
    !
    ! Passed Variables
    !---------------------
    real(kind=real_kind),intent(in):: v(np,np,2)
    type(derivative_t)             :: deriv
    type(element_t)                :: elem
    real(kind=real_kind)           :: vort(np,np)
    !
    ! Local Values
    !--------------
    real(kind=real_kind):: dvdx00
    real(kind=real_kind):: dudy00
    real(kind=real_kind):: vco(np,np,2)
    real(kind=real_kind):: vtemp(np,np)

    integer ii,jj,ll
    
    ! convert to covariant form
    !--------------------------
    do jj=1,np
    do ii=1,np
      vco(ii,jj,1)=(elem%D(1,1,ii,jj)*v(ii,jj,1) + elem%D(2,1,ii,jj)*v(ii,jj,2))
      vco(ii,jj,2)=(elem%D(1,2,ii,jj)*v(ii,jj,1) + elem%D(2,2,ii,jj)*v(ii,jj,2))
    end do
    end do

    do jj=1,np
    do ll=1,np
      dudy00=0.0d0
      dvdx00=0.0d0
      do ii=1,np
        dvdx00 = dvdx00 + deriv%Dvv(ii,ll)*vco(ii,jj,2)
        dudy00 = dudy00 + deriv%Dvv(ii,ll)*vco(jj,ii,1)
      end do
      vort (ll,jj) = dvdx00
      vtemp(jj,ll) = dudy00
    end do
    end do

    do jj=1,np
    do ii=1,np
      vort(ii,jj)=(vort(ii,jj)-vtemp(ii,jj))*(elem%rmetdet(ii,jj)*rrearth)
    end do
    end do

    ! End Function
    !-------------
    return
  end function vorticity_sphere
  !==================================================================


  !==================================================================
  function vorticity_sphere_diag(v,deriv,elem) result(vort)
    ! vorticity_sphere_diag:
    !
    !   input:  v = velocity in lat-lon coordinates
    !   ouput:  diagonal component of spherical vorticity of v
    !=============================================================
    !
    ! Passed Variables
    !----------------------
    real(kind=real_kind),intent(in):: v(np,np,2)
    type(derivative_t)             :: deriv
    type(element_t)                :: elem
    real(kind=real_kind)           :: vort(np,np)
    !
    ! Local Values
    !----------------
    real(kind=real_kind):: dvdx00
    real(kind=real_kind):: dudy00
    real(kind=real_kind)::   vco(np,np,2)
    real(kind=real_kind):: vtemp(np,np)
    real(kind=real_kind):: rdx
    real(kind=real_kind):: rdy

    integer ii,jj,ll

    ! convert to covariant form
    !---------------------------
    do jj=1,np
    do ii=1,np
      vco(ii,jj,1)=(elem%D(1,1,ii,jj)*v(ii,jj,1) + elem%D(2,1,ii,jj)*v(ii,jj,2))
      vco(ii,jj,2)=(elem%D(1,2,ii,jj)*v(ii,jj,1) + elem%D(2,2,ii,jj)*v(ii,jj,2))
    end do
    end do

    do jj=1,np
    do ll=1,np
      dudy00=0.0d0
      dvdx00=0.0d0
      do ii=1,np
        dvdx00 = dvdx00 + deriv%Dvv_diag(ii,ll)*vco(ii,jj,2)
        dudy00 = dudy00 + deriv%Dvv_diag(ii,ll)*vco(jj,ii,1)
      enddo 
      vort (ll,jj) = dvdx00 
      vtemp(jj,ll) = dudy00
    end do
    end do

    do jj=1,np
    do ii=1,np 
      vort(ii,jj)=(vort(ii,jj)-vtemp(ii,jj))*(elem%rmetdet(ii,jj)*rrearth)
    end do 
    end do 
     
    ! End Function
    !-------------
    return
  end function vorticity_sphere_diag
  !==================================================================


  !==================================================================
  function divergence_sphere(v,deriv,elem) result(div)
    ! divergence_sphere:
    !
    !   input:  v = velocity in lat-lon coordinates
    !   ouput:  div(v)  spherical divergence of v
    !================================================================
    !
    ! Passed Variables
    !------------------
    real(kind=real_kind),intent(in):: v(np,np,2)  ! in lat-lon coordinates
    type(derivative_t)             :: deriv
    type(element_t)                :: elem
    real(kind=real_kind)           :: div(np,np)
    !
    ! Local Values
    !----------------
    real(kind=real_kind):: dudx00
    real(kind=real_kind):: dvdy00
    real(kind=real_kind):: gv(np,np,2),vvtemp(np,np)

    integer ii,jj,ll

    ! convert to contra variant form and multiply by g
    !--------------------------------------------------
    do jj=1,np
    do ii=1,np
      gv(ii,jj,1)=elem%metdet(ii,jj)*( elem%Dinv(1,1,ii,jj)*v(ii,jj,1) &
                                      +elem%Dinv(1,2,ii,jj)*v(ii,jj,2) )
      gv(ii,jj,2)=elem%metdet(ii,jj)*( elem%Dinv(2,1,ii,jj)*v(ii,jj,1) &
                                      +elem%Dinv(2,2,ii,jj)*v(ii,jj,2) )
    end do
    end do

    ! compute d/dx and d/dy         
    !-----------------------
    do jj=1,np
    do ll=1,np
      dudx00=0.0d0
      dvdy00=0.0d0
      do ii=1,np
        dudx00 = dudx00 + deriv%Dvv(ii,ll)*gv(ii,jj,1)
        dvdy00 = dvdy00 + deriv%Dvv(ii,ll)*gv(jj,ii,2)
      end do
      div   (ll,jj) = dudx00
      vvtemp(jj,ll) = dvdy00
    end do
    end do

    do jj=1,np
    do ii=1,np
      div(ii,jj)=(div(ii,jj)+vvtemp(ii,jj))*(elem%rmetdet(ii,jj)*rrearth)
    end do
    end do
    
    ! End Function
    !-------------
    return
  end function divergence_sphere
  !==================================================================


  !==================================================================
  function laplace_sphere_wk(s,deriv,elem,var_coef) result(laplace)
    ! laplace_sphere_wk:
    !
    !   input:  s = scalar
    !   ouput:  -< grad(PHI), grad(s) >   = weak divergence of grad(s)
    !     note: for this form of the operator, grad(s) does not need to be made C0
    !==============================================================================
    !
    ! Passed Variables
    !-------------------
    real(kind=real_kind),intent(in):: s(np,np) 
    type(derivative_t)             :: deriv
    type(element_t)                :: elem
    logical                        :: var_coef
    real(kind=real_kind)           :: laplace(np,np)
    !
    ! Local Values
    !---------------
    real(kind=real_kind):: laplace2(np,np)
    real(kind=real_kind):: grads(np,np,2), oldgrads(np,np,2)

    integer ii,jj

    grads=gradient_sphere(s,deriv,elem%Dinv)
 
    if(var_coef) then
      if(hypervis_power/=0 ) then
        ! scalar viscosity with variable coefficient
        !--------------------------------------------
        grads(:,:,1) = grads(:,:,1)*elem%variable_hyperviscosity(:,:)
        grads(:,:,2) = grads(:,:,2)*elem%variable_hyperviscosity(:,:)
      elseif(hypervis_scaling /=0 ) then
        ! tensor hv, (3)
        !-----------------
        oldgrads=grads
        do jj=1,np
        do ii=1,np
          grads(ii,jj,1) = sum(oldgrads(ii,jj,:)*elem%tensorVisc(1,:,ii,jj))
          grads(ii,jj,2) = sum(oldgrads(ii,jj,:)*elem%tensorVisc(2,:,ii,jj))
        end do
        end do
      else
        ! do nothing: constant coefficient viscsoity
        !--------------------------------------------
      endif
    endif

    ! note: divergnece_sphere and divergence_sphere_wk are identical *after* 
    !       bndry_exchange if input is C_0.  Here input is not C_0, so we 
    !       should use divergence_sphere_wk().  
    !-----------------------------------------------------------------------
    laplace=divergence_sphere_wk(grads,deriv,elem)

    ! End Function
    !-------------
    return
  end function laplace_sphere_wk
  !==================================================================


  !==================================================================
  function vlaplace_sphere_wk(v,deriv,elem,var_coef,nu_ratio) result(laplace)
    ! vlaplace_sphere_wk:
    !
    !   input:  v = vector in lat-lon coordinates
    !   ouput:  weak laplacian of v, in lat-lon coordinates
    !
    !   logic:
    !      tensorHV:     requires cartesian
    !      nu_div/=nu:   requires contra formulatino
    !
    !   One combination NOT supported:  tensorHV and nu_div/=nu then endrun
    !======================================================================
    !
    ! Passed Variables
    !-------------------
    real(kind=real_kind),intent(in):: v(np,np,2) 
    type(derivative_t)             :: deriv
    type(element_t)                :: elem
    logical                        :: var_coef
    real(kind=real_kind),optional  :: nu_ratio
    real(kind=real_kind)           :: laplace(np,np,2)

    if((hypervis_scaling/=0).and.(var_coef)) then
      ! tensorHV is turned on - requires cartesian formulation
      !---------------------------------------------------------
      if(present(nu_ratio)) then
        if(nu_ratio /= 1) then
          call endrun('ERROR: tensorHV can not be used with nu_div/=nu')
        endif
      endif
      laplace=vlaplace_sphere_wk_cartesian(v,deriv,elem,var_coef)
    else  
      ! all other cases, use contra formulation:
      !-----------------------------------------
      laplace=vlaplace_sphere_wk_contra(v,deriv,elem,var_coef,nu_ratio)
    endif

    ! End Function
    !-------------
    return
  end function vlaplace_sphere_wk
  !==================================================================


  !==================================================================
  function vlaplace_sphere_wk_cartesian(v,deriv,elem,var_coef) result(laplace)
    ! vlaplace_sphere_wk_cartesian:
    !
    !   input:  v = vector in lat-lon coordinates
    !   ouput:  weak laplacian of v, in lat-lon coordinates
    !====================================================================
    !
    ! Passed Variables
    !--------------------
    real(kind=real_kind),intent(in):: v(np,np,2) 
    type(derivative_t)             :: deriv
    type(element_t)                :: elem
    logical                        :: var_coef
    real(kind=real_kind)           :: laplace(np,np,2)
    !
    ! Local Values
    !------------------
    real(kind=real_kind) :: dum_cart(np,np,3)
    integer component

    ! latlon -> cartesian
    !----------------------
    do component=1,3
      dum_cart(:,:,component)=sum( elem%vec_sphere2cart(:,:,component,:)*v(:,:,:) ,3)
    end do

    ! Do laplace on cartesian comps
    !-------------------------------
    do component=1,3
      dum_cart(:,:,component) = laplace_sphere_wk(dum_cart(:,:,component),deriv,elem,var_coef)
    end do

    ! cartesian -> latlon
    !---------------------
    do component=1,2
      ! vec_sphere2cart is its own pseudoinverse.
      !-------------------------------------------
      laplace(:,:,component)=sum( dum_cart(:,:,:)*elem%vec_sphere2cart(:,:,:,component) ,3)
    end do 

    ! End Function
    !-------------
    return
  end function vlaplace_sphere_wk_cartesian
  !==================================================================


  !==================================================================
  function vlaplace_sphere_wk_contra(v,deriv,elem,var_coef,nu_ratio) result(laplace)
    ! vlaplace_sphere_wk_contra:
    !
    !   input:  v = vector in lat-lon coordinates
    !   ouput:  weak laplacian of v, in lat-lon coordinates
    !
    !   du/dt = laplace(u) = grad(div) - curl(vor)
    !   < PHI du/dt > = < PHI laplace(u) >        PHI = covariant, u = contravariant
    !                 = < PHI grad(div) >  - < PHI curl(vor) >
    !                 = grad_wk(div) - curl_wk(vor)               
    !================================================================================
    !
    ! Passed Variables
    !--------------------
    real(kind=real_kind),intent(in):: v(np,np,2) 
    type(derivative_t)             :: deriv
    type(element_t)                :: elem
    logical                        :: var_coef
    real(kind=real_kind),optional  :: nu_ratio
    real(kind=real_kind)           :: laplace(np,np,2)
    !
    ! Local Values
    !------------------
    real(kind=real_kind):: vor(np,np),div(np,np)
    real(kind=real_kind):: v1,v2,div1,div2,vor1,vor2,phi_x,phi_y

    integer ii,jj,ll,mm,nn

    div=divergence_sphere(v,deriv,elem)
    vor= vorticity_sphere(v,deriv,elem)

    if((var_coef).and.(hypervis_power/=0)) then
      ! scalar viscosity with variable coefficient
      !---------------------------------------------
      div = div*elem%variable_hyperviscosity(:,:)
      vor = vor*elem%variable_hyperviscosity(:,:)
    endif

    if(present(nu_ratio)) div = nu_ratio*div

    laplace = gradient_sphere_wk_testcov(div,deriv,elem) &
                 -curl_sphere_wk_testcov(vor,deriv,elem)

#define UNDAMPRR
#ifdef UNDAMPRR
    do nn=1,np
    do mm=1,np
      ! add in correction so we dont damp rigid rotation
      !-------------------------------------------------
      laplace(mm,nn,1)=laplace(mm,nn,1) + 2*elem%spheremp(mm,nn)*v(mm,nn,1)*(rrearth**2)
      laplace(mm,nn,2)=laplace(mm,nn,2) + 2*elem%spheremp(mm,nn)*v(mm,nn,2)*(rrearth**2)
    end do
    end do
#endif

    ! End Function
    !-------------
    return
  end function vlaplace_sphere_wk_contra
  !==================================================================

!-----------------------------------------------------------------------------------

  !==================================================================
  function gll_to_dgmodal(p,deriv) result(phat)
    ! gll_to_dgmodal:
    !
    !   input:  v = velocity in lat-lon coordinates
    !   ouput:  phat = Legendre coefficients
    !
    !   Computes  < g dot p  > = SUM  g(i,j) p(i,j) w(i) w(j)
    !   (the quadrature approximation on the *reference element* of the integral of p against
    !    all Legendre polynomials up to degree npdg
    !
    !   for npdg < np, this routine gives the (exact) modal expansion of p/spheremp()
    !================================================================================
    !
    ! Passed Variables
    !---------------------
    real(kind=real_kind),intent(in):: p(np,np) 
    type(derivative_t)             :: deriv
    real(kind=real_kind)           :: phat(npdg,npdg)
    !
    ! Local Values
    !------------------
    real(kind=real_kind):: A(np,npdg)

    integer ii,jj,mm,nn

    A   =0
    phat=0

    ! N^3 tensor product formulation:
    !---------------------------------
    do mm=1,npdg
    do jj=1,np
    do ii=1,np
      A(jj,mm)=A(jj,mm)+(p(ii,jj)*deriv%Mvv_twt(ii,ii)                   &
                                 *deriv%Mvv_twt(jj,jj))*deriv%legdg(mm,ii)
    end do
    end do
    end do

    do nn=1,npdg
    do mm=1,npdg
    do jj=1,np
      phat(mm,nn)=phat(mm,nn)+A(jj,mm)*deriv%legdg(nn,jj)
    end do
    end do
    end do
    
#if 0
    do mm=1,npdg
    do nn=1,npdg
    do jj=1,np
    do ii=1,np
      gmn = deriv%legdg(mm,ii)*deriv%legdg(nn,jj)                  ! basis function
      phat(mm,nn)=phat(mm,nn)+gmn*p(ii,jj)*deriv%Mvv_twt(ii,ii)*deriv%Mvv_twt(jj,jj)
    end do
    end do
    end do
    end do
#endif

    ! End Function
    !-------------
    return
  end function gll_to_dgmodal
  !==================================================================


  !==================================================================
  function dgmodal_to_gll(phat,deriv) result(p)
    ! dgmodal_to_gll:
    !
    !   input:  phat = coefficients of Legendre expansion
    !   ouput:  p    = sum expansion to evaluate phat at GLL points
    !==================================================================
    !
    ! Passed variables
    !-------------------
    real(kind=real_kind):: phat(npdg,npdg)
    type(derivative_t)  :: deriv
    real(kind=real_kind):: p(np,np) 
    !
    ! Local Values
    !---------------
    real(kind=real_kind):: A(npdg,np)

    integer ii,jj,mm,nn

    p(:,:)=0
    A     =0

    ! tensor product version
    !-----------------------
    do ii=1,np
    do nn=1,npdg
    do mm=1,npdg
      A(nn,ii)=A(nn,ii)+phat(mm,nn)*deriv%legdg(mm,ii)
    end do
    end do
    end do

    do jj=1,np
    do ii=1,np
    do nn=1,npdg
      p(ii,jj)=p(ii,jj)+A(nn,ii)*deriv%legdg(nn,jj)
    end do
    end do
    end do

#if 0
    do jj=1,np
    do ii=1,np
    do mm=1,npdg
    do nn=1,npdg
      p(ii,jj)=p(ii,jj)+phat(mm,nn)*deriv%legdg(mm,ii)*deriv%legdg(nn,jj) 
    end do
    end do
    end do
    end do
#endif

    ! End Function
    !-------------
    return
  end function dgmodal_to_gll
  !==================================================================


  !==================================================================
  function subcell_integration(sampled_val, metdet, np, intervals) result(values)
    ! subcell_integration:
    !
    ! Given a field defined on the unit element, [-1,1]x[-1,1]
    ! sample values, sampled_val, and integration weights, metdet,
    ! at a number, np, of Gauss-Lobatto-Legendre points. Divide
    ! the square up into intervals by intervals sub-squares so that
    ! there are now intervals**2 sub-cells.  Integrate the 
    ! function defined by sampled_val and metdet over each of these
    ! sub-cells and return the integrated values as an 
    ! intervals by intervals matrix.
    !
    ! Efficiency is obtained by computing and caching the appropriate
    ! integration matrix the first time the function is called.
    !============================================================================
    !
    ! Passed Variables
    !------------------
    real(kind=real_kind),intent(in):: sampled_val(np,np)
    real(kind=real_kind),intent(in):: metdet     (np,np)
    integer             ,intent(in):: np
    integer             ,intent(in):: intervals
    real(kind=real_kind)           :: values(intervals,intervals)
    !
    ! Local Values
    !----------------
    real(kind=real_kind):: V(np,np)

    integer ii,jj

    V  = sampled_val * metdet

    if((.not.allocated(integration_matrix)               ).or. &
       (          size(integration_matrix,1).ne.intervals).or. &
       (          size(integration_matrix,2).ne.np       )     ) then
      call allocate_subcell_integration_matrix(np,intervals)
    endif

    ! Multiply the sampled values by the weighted jacobians.  
    ! Symmetry allows us to write this as J^t V J
    ! where J is a vector.  
    !-------------------------------------------------------
    values = matmul(integration_matrix,matmul(V,transpose(integration_matrix)))

    ! End Function
    !-------------
    return
  end function subcell_integration
  !==================================================================


  !==================================================================
  subroutine allocate_subcell_integration_matrix(np, intervals)
    ! allocate_subcell_integration_matrix:
    !
    ! Helper subroutine that will fill in a matrix needed to 
    ! integrate a function defined on the GLL points of a unit
    ! square on sub-cells.  So np is the number of integration
    ! GLL points defined on the unit square (actually [-1,1]x[-1,1])
    ! and intervals is the number to cut it up into, say a 3 by 3
    ! set of uniform sub-cells.  This function will fill the 
    ! subcell_integration matrix with the correct coefficients
    ! to integrate over each subcell.  
    !==================================================================
    !
    ! Passed Variables
    !-------------------
    integer,intent(in):: np
    integer,intent(in):: intervals
    !
    ! Local Values
    !---------------
    real(kind=real_kind),parameter:: zero=0.0D0, one=1.0D0, two=2.0D0
    real(kind=real_kind)          :: values(intervals,intervals)
    real(kind=real_kind)          :: sub_gll(intervals,np)
    real(kind=real_kind)          :: Lagrange_interp(intervals,np,np)
    type(quadrature_t)            :: gll 
    real(kind=real_kind)          :: legrange_div(np)
    real(kind=real_kind)          :: aa,bb,xx,yy, x_j, x_i 
    real(kind=real_kind)          :: rr(1) 

    integer ii,jj,nn,mm

    if(allocated(integration_matrix)) deallocate(integration_matrix)
    allocate(integration_matrix(intervals,np))

    gll = gausslobatto(np)
 
    ! The GLL (Gauss-Lobatto-Legendre) points are from [-1,1], 
    ! we have a bunch of sub-intervals defined by intervals that 
    ! go from [a,b] so we need to linearly map [-1,1] -> [a,b] 
    ! all the  GLL points by  y = (a/2)(1-x) + (b/2)(1+x)
    !-----------------------------------------------------------
    do ii=1,intervals
      aa = -one + (ii-one)*two/intervals   
      bb = -one +  ii     *two/intervals  
      sub_gll(ii,:) = (aa+bb)/two + gll%points(:)/intervals
    end do

    ! Now to interpolate from the values at the input GLL
    ! points to the sub-GLL points.  Do this by Lagrange
    ! interpolation.  The jth Lagrange interpolating polynomial
    ! for points x_i is 
    !              \prod_{i\ne j} (x-x_i)/(x_j-x_i)
    ! These are then multiplied by the sampled values y_i 
    ! and summed. 
    
    ! Save some time by pre-computing the denominitor. I think 
    ! this is OK since all the points are of order 1 so should
    ! be well behaved.
    !------------------------------------------------------------
    do nn = 1,np
      x_j = gll%points(nn)
      xx  = one 
      do mm = 1,np 
        if(mm.ne.nn) then
          x_i = gll%points(mm)
          xx  = xx*(x_j - x_i)
        endif
      end do
      legrange_div(nn)= xx
    end do 

    do ii=1,intervals
    do nn=1,np
      xx = sub_gll(ii,nn)
      do jj = 1,np
        yy = one
        do mm = 1,np
          if(mm.ne.jj) then
            x_i = gll%points(mm)
            yy = yy*(xx - x_i)
          endif
        end do
        Lagrange_interp(ii,nn,jj) = yy/legrange_div(jj)
      end do
    end do
    end do

    ! Integration is the GLL weights times Jacobians times
    ! the interpolated values:
    !                   w^t I Y I^t w 
    ! where  
    ! w is GLL weights and Jacobians, 
    ! I is the Lagrange_interp matrix, and
    ! Y is the coefficient matrix, sampled_val.
    ! This can be written  J Y J^t where
    !                       J = w^t I
    ! J is integration_matrix
    !-------------------------------------------------------
    do ii=1,intervals
      integration_matrix(ii,:) = matmul(gll%weights(:),Lagrange_interp(ii,:,:))
    end do

    ! There is still the Jacobian to consider.  We are 
    ! integrating over [a,b] x [c,d] where 
    !        |b-a| = |d-c| = 2/Intervals
    ! Multiply the weights appropriately given that 
    ! they are defined for a 2x2 square
    !------------------------------------------------------
    integration_matrix = integration_matrix/intervals

    ! End Routine
    !-------------
    return
  end subroutine allocate_subcell_integration_matrix
  !==================================================================

end module derivative_mod
