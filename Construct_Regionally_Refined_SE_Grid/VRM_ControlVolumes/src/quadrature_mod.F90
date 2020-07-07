#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#undef _GAUSS_TABLE
#undef _QUAD_DBG
module quadrature_mod

  ! Useful modules
  !----------------
  use SE_Constants, only : longdouble_kind

  ! Set all Global values and routines to private by default
  ! and then explicitly set their exposure
  !-----------------------------------------------------------
  implicit none
  private

  public :: quadrature_t

  public :: gauss
  private:: gauss_pts
  private:: gauss_wts
  public :: test_gauss
  public :: gausslobatto
  private:: gausslobatto_pts
  private:: gausslobatto_wts
  public :: test_gausslobatto
  public :: jacobi
  private:: jacobi_polynomials
  private:: jacobi_derivatives
  public :: legendre
  public :: quad_norm
  private:: trapN
  public :: trapezoid
  public :: simpsons
  public :: gaussian_int

  ! Type Definitions
  !--------------------
  type quadrature_t
     real (kind=longdouble_kind), dimension(:), pointer :: points
     real (kind=longdouble_kind), dimension(:), pointer :: weights
  end type quadrature_t


contains

  !==================================================================
  function gauss(npts) result(gs)
    ! gauss: Find the Gauss collocation points and the corresponding weights.
    !
    !==============================================================
    !
    ! Passed Variables
    !----------------------
    integer,intent(in):: npts
    type(quadrature_t):: gs

    allocate(gs%points (npts))
    allocate(gs%weights(npts))

    gs%points =gauss_pts(npts)
    gs%weights=gauss_wts(npts,gs%points)

    ! End Function
    !-------------
    return
  end function gauss
  !==================================================================


#if defined(_GAUSS_TABLE)
  !==================================================================
  function gauss_pts(npts) result(pts)
    ! gauss_pts:
    !
    !==============================================================
    !
    ! Passed Variables
    !--------------------
    integer,intent(in)        :: npts
    real(kind=longdouble_kind):: pts(npts)

    pts(1) = -0.93246951420315202781d0
    pts(2) = -0.66120938646626451366d0
    pts(3) = -0.23861918608319690863d0
    pts(4) = -pts(3)
    pts(5) = -pts(2)
    pts(6) = -pts(1)

    ! End Function
    !-------------
    return
  end function gauss_pts
  !==================================================================


  !==================================================================
  function gauss_wts(npts,pts) result(wts)
    ! gauss_wts:
    !
    !===============================================================
    !
    ! Passed Variables
    !----------------------
    integer,intent(in)        :: npts
    real(kind=longdouble_kind):: pts(npts)
    real(kind=longdouble_kind):: wts(npts)

    wts(1)  =  0.17132449237917034504d0
    wts(2)  =  0.36076157304813860756d0
    wts(3)  =  0.46791393457269104738d0
    wts(4)  =  wts(3)
    wts(5)  =  wts(2)
    wts(6)  =  wts(1)

    ! End Function
    !-------------
    return
  end function gauss_wts
  !==================================================================
#else
  !==================================================================
  function gauss_pts(np1) result(pts)
    ! gauss_pts: Compute the Gauss Collocation points for Jacobi Polynomials
    !
    !==============================================================
    use SE_Constants,only: QQ_PI
    !
    ! Passed Variables
    !-----------------
    integer,intent(in)        :: np1        ! Number of velocity grid points
    real(kind=longdouble_kind):: pts(np1)
    !
    ! Local Values
    !---------------
    real(kind=longdouble_kind),parameter:: convthresh = 10 ! convergence threshold relative to machine epsilon
    real(kind=longdouble_kind)          :: eps             ! machine epsilon
    integer                             :: prec            ! number of mantissa bits
    integer                  ,parameter :: kstop = 30      ! max iterations for polynomial deflation

    real(kind=longdouble_kind):: alpha,beta
    real(kind=longdouble_kind):: xjac(0:np1-1)
    real(kind=longdouble_kind)::  jac(0:np1)
    real(kind=longdouble_kind):: djac(0:np1)
    real(kind=longdouble_kind):: poly
    real(kind=longdouble_kind):: pder
    real(kind=longdouble_kind):: recsum,thresh
    real(kind=longdouble_kind):: dth
    real(kind=longdouble_kind):: xx
    real(kind=longdouble_kind):: delx
    real(kind=longdouble_kind):: c0,c1,c2,c10

    integer ii,jj,kk
    integer nn,nh

    nn   = np1 - 1
    c0   = 0.0_longdouble_kind
    c1   = 1.0_longdouble_kind
    c2   = 2.0_longdouble_kind
    c10  = 10.0_longdouble_kind
    alpha= c0
    beta = c0

    ! compute machine precision and set the convergence
    ! threshold thresh to 10 times that level
    !----------------------------------------------------
    prec   = precision(c10)
    eps    = c10**(-prec)
    thresh = convthresh*eps

    ! Compute first half of the roots by "polynomial deflation".
    !------------------------------------------------------------
    dth = QQ_PI/(2*nn+2)
    nh  = (nn+1)/2
    do jj=0,nh-1
      xx=COS((c2*jj+1)*dth)          ! first guess at root
      kk=0
      delx=c1
      do while((kk < kstop).and.(ABS(delx) > thresh))
        call jacobi(nn+1,xx,alpha,beta,jac(0:nn+1),djac(0:nn+1))
        poly =  jac(nn+1)
        pder = djac(nn+1)
        recsum=c0
        do ii=0,jj-1
          recsum = recsum + c1/(xx-xjac(ii))
        end do
        delx = -poly/(pder-recsum*poly)
        xx = xx + delx
        kk = kk + 1
      end do ! while((kk < kstop).and.(ABS(delx) > thresh))
      xjac(jj)=xx
    end do

    ! compute the second half of the roots by symmetry
    !--------------------------------------------------
    do jj=0,nh
      xjac(nn-jj) = -xjac(jj)
    end do
    if(MODULO(nn,2)==0) xjac(nh)=c0

    ! Reverse the sign of everything so that indexing
    ! increases with position
    !-------------------------------------------------
    do jj=0,nn
      pts(jj+1) = -xjac(jj)
    end do

    ! End Function
    !-------------
    return
  end function gauss_pts
  !==================================================================


  !==================================================================
  function gauss_wts(np1, gpts) result(wts)
    ! gauss_wts: Gauss Legendre Weights
    !
    !=====================================================
    !
    ! Passed Variables
    !----------------------
    integer,intent(in)                   :: np1
    real(kind=longdouble_kind),intent(in):: gpts(np1)  ! Gauss-Legendre points
    real(kind=longdouble_kind)           :: wts(np1)   ! Gauss-Legendre weights
    !
    ! Local Values
    !--------------
    real(kind=longdouble_kind):: c0,c1,c2
    real(kind=longdouble_kind):: alpha
    real(kind=longdouble_kind):: beta
    real(kind=longdouble_kind):: djac(np1)

    integer ii,nn

    c0    = 0.0_longdouble_kind
    c1    = 1.0_longdouble_kind
    c2    = 2.0_longdouble_kind
    alpha = c0
    beta  = c0
    nn    = np1-1

    djac=jacobi_derivatives(np1,alpha,beta,np1,gpts)     

    do ii=1,np1
      wts(ii)=c2/((c1-gpts(ii)**2)*djac(ii)*djac(ii))
    end do

    ! End Function
    !-------------
    return
  end function gauss_wts
  !==================================================================
#endif


  !==================================================================
  subroutine test_gauss(npts)
    ! test_gauss: Unit Tester for Gaussian Points, Weights
    !==============================================================
    use SE_Constants,only: real_kind
    !
    ! Passed Variables
    !----------------------
    integer,intent(in):: npts
    ! 
    ! Local Values
    !--------------
    type(quadrature_t)  :: gs
    real(kind=real_kind):: gssum

    integer ii

    print *,' npts=',npts

    gs=gauss(npts)

    print *
    print *,"============================================"
    print *,"        Testing Gaussian Quadrature..."
    print *
    print *,"          points              weights"
    print *,"============================================"
    do ii=1,npts
      print *,ii,gs%points(ii),gs%weights(ii)
    end do

    print *,"============================================"
    gssum=SUM(gs%weights(:))
    print *,"sum of Gaussian weights=",gssum
    print *,"============================================"

    deallocate(gs%points )
    deallocate(gs%weights)

    ! End Routine
    !-------------
    return
  end subroutine test_gauss
  !==================================================================


  !==================================================================
  function gausslobatto(npts) result(gll)
    ! gausslobatto: Find the Gauss-Lobatto Legendre collocation points xgl(i) 
    !               and the corresponding weights.
    !==============================================================
    !
    ! Passed Variables
    !-------------------
    integer,intent(in):: npts
    type(quadrature_t):: gll

    allocate(gll%points (npts))
    allocate(gll%weights(npts))

    gll%points =gausslobatto_pts(npts)
    gll%weights=gausslobatto_wts(npts,gll%points)

    ! End Function
    !-------------
    return
  end function gausslobatto
  !==================================================================


  !==================================================================
  function gausslobatto_pts(np1) result(pts)
    ! gausslobatto_pts: Compute the Gauss-Lobatto Collocation points 
    !                   for Jacobi Polynomials
    !==============================================================
    use SE_Constants, only : QQ_PI
    !
    ! Passed Variables
    !--------------------
    integer,intent(in)        :: np1        ! Number of velocity grid points
    real(kind=longdouble_kind):: pts(np1)
    !
    ! Local Values
    !---------------
    real(kind=longdouble_kind),parameter:: convthresh = 10 ! convergence threshold relative to machine epsilon 
    real(kind=longdouble_kind)          :: eps             ! machine epsilon
    integer                             :: prec            ! number of mantissa bits 
    integer                   ,parameter:: kstop = 30      ! max iterations for polynomial deflation

    real(kind=longdouble_kind):: alpha,beta
    real(kind=longdouble_kind)::  xjac(0:np1-1)
    real(kind=longdouble_kind)::   jac(0:np1)
    real(kind=longdouble_kind):: jacm1(0:np1)
    real(kind=longdouble_kind)::  djac(0:np1)
    real(kind=longdouble_kind):: aa,bb,det
    real(kind=longdouble_kind):: poly
    real(kind=longdouble_kind):: pder
    real(kind=longdouble_kind):: recsum,thresh
    real(kind=longdouble_kind):: dth,cd,sd,cs,ss,cstmp
    real(kind=longdouble_kind):: xx
    real(kind=longdouble_kind):: delx
    real(kind=longdouble_kind):: c0,c1,c2,c10

    integer ii,jj,kk
    integer nn,nh

    nn   = np1 - 1
    c0   = 0.0_longdouble_kind
    c1   = 1.0_longdouble_kind
    c2   = 2.0_longdouble_kind
    c10  = 10.0_longdouble_kind
    alpha= c0
    beta = c0

    ! compute machine precision and set the convergence
    ! threshold thresh to 10 times that level 
    !---------------------------------------------------
    prec   = PRECISION(c10)
    eps    = c10**(-prec)
    thresh = convthresh*eps

    ! initialize the end points
    !-----------------------------
    xjac(0) =  c1
    xjac(nn)= -c1

    !===========================================================
    ! Compute first half of the roots by "polynomial deflation".
    !===========================================================

    ! compute the parameters in the polynomial whose 
    ! roots are desired...
    !--------------------------------------------------
    call jacobi(nn+1, c1,alpha,beta,  jac(0:nn+1),djac(0:nn+1))
    call jacobi(nn+1,-c1,alpha,beta,jacm1(0:nn+1),djac(0:nn+1))

    det =   jac(nn  )*jacm1(nn-1)-jacm1(nn  )*jac(nn-1)
    aa  = -(jac(nn+1)*jacm1(nn-1)-jacm1(nn+1)*jac(nn-1))/det
    bb  = -(jac(nn  )*jacm1(nn+1)-jacm1(nn  )*jac(nn+1))/det

    dth = QQ_PI/(2*nn+1)
    cd  = COS(c2*dth)
    sd  = SIN(c2*dth)
    cs  = COS(dth)
    ss  = SIN(dth)
    nh  = (nn+1)/2

    do jj=1,nh-1
      xx=cs          ! first guess at root 
      kk=0
      delx=c1
      do while((kk < kstop).and.(ABS(delx) > thresh))
        call jacobi(nn+1,xx,alpha,beta,jac(0:nn+1),djac(0:nn+1))
        poly =  jac(nn+1)+aa* jac(nn)+bb* jac(nn-1)
        pder = djac(nn+1)+aa*djac(nn)+bb*djac(nn-1)
        recsum=c0
        do ii=0,jj-1
          recsum = recsum + c1/(xx-xjac(ii))
        end do
        delx = -poly/(pder-recsum*poly)
        xx = xx + delx
        kk = kk + 1
      end do ! while((kk < kstop).and.(ABS(delx) > thresh))
      xjac(jj)=xx

      !  compute the guesses for the roots
      !  for the next points, i.e :
      !
      !  ss = sn(theta) => sin(theta+2*dth)
      !  cs = cs(theta) => cs(theta+2*dth)
      !-------------------------------------
      cstmp=cs*cd-ss*sd      
      ss=cs*sd+ss*cd    
      cs=cstmp          
    end do

    ! compute the second half of the roots by symmetry
    !-------------------------------------------------
    do jj=1,nh 
      xjac(nn-jj) = -xjac(jj)
    end do
    if(MODULO(nn,2)==0) xjac(nh)=c0

    ! Reverse the sign of everything so that indexing
    ! increases with position          
    !----------------------------------------------------
    do jj=0,nn
      pts(jj+1) = -xjac(jj)
    end do

    ! End Function
    !-------------
    return
  end function gausslobatto_pts
  !==================================================================


  !==================================================================
  function gausslobatto_wts(np1, glpts) result(wts)
    ! gausslobatto_wts: Gauss Lobatto Legendre Weights   
    !==========================================================
    !
    ! Passed Variables
    !------------------
    integer,intent(in)                   :: np1
    real(kind=longdouble_kind),intent(in):: glpts(np1)
    real(kind=longdouble_kind)           ::   wts(np1)
    !
    ! Local Values
    !----------------
    real(kind=longdouble_kind):: c0,c2
    real(kind=longdouble_kind):: alpha
    real(kind=longdouble_kind):: beta
    real(kind=longdouble_kind):: jac(np1)

    integer ii,nn

    c0    = 0.0_longdouble_kind
    c2    = 2.0_longdouble_kind
    alpha = c0
    beta  = c0
    nn    = np1-1

    jac=jacobi_polynomials(nn,alpha,beta,np1,glpts)

    do ii=1,np1
      wts(ii)=c2/(nn*(nn+1)*jac(ii)*jac(ii))
    end do

    ! End Function
    !-------------
    return
  end function gausslobatto_wts
  !==================================================================


  !==================================================================
  subroutine test_gausslobatto(npts)
    ! test_gausslobatto: Unit Tester for Gaussian Lobatto Quadrature...
    !
    !==============================================================
    use SE_Constants,only: real_kind
    !
    ! Passed Variables
    !-------------------
    integer,intent(in):: npts
    !
    ! Local Values
    !----------------
    type(quadrature_t)  :: gll
    real(kind=real_kind):: gllsum

    integer ii

    gll=gausslobatto(npts)

    print *
    print *,"============================================"
    print *,"      Testing Gauss-Lobatto Quadrature..."
    print *
    print *,"          points              weights"
    print *,"============================================"
    do ii=1,npts
      print *,ii,gll%points(ii),gll%weights(ii)
    end do

    print *,"============================================"
    gllsum=SUM(gll%weights(:))
    print *,"sum of Gauss-Lobatto weights=",gllsum
    print *,"============================================"

    deallocate(gll%points )
    deallocate(gll%weights)

    ! End Routine
    !-------------
    return
  end subroutine test_gausslobatto
  !==================================================================


  !==================================================================
  subroutine jacobi(nn,xx,alpha,beta,jac,djac)
    ! subroutine jacobi:
    !
    !  Computes the Jacobi Polynomials (jac) and their
    !  first derivatives up to and including degree n 
    !  at point x on the interval (-1,1).
    !
    !    See for example the recurrence relations 
    !    in equation 2.5.4 (page 70) in 
    !
    !    "Spectral Methods in Fluid Dynamics",
    !    by C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A.Zang
    !    Springer-Verlag, 1988.
    !=================================================================
    !
    ! Passes Variables
    !------------------
    integer                   ,intent(in):: nn
    real(kind=longdouble_kind),intent(in):: xx
    real(kind=longdouble_kind),intent(in):: alpha
    real(kind=longdouble_kind),intent(in):: beta
    real(kind=longdouble_kind)           :: jac(0:nn)
    real(kind=longdouble_kind)           :: djac(0:nn)
    !
    ! Local Values
    !----------------
    real(kind=longdouble_kind):: a1k
    real(kind=longdouble_kind):: a2k
    real(kind=longdouble_kind):: a3k
    real(kind=longdouble_kind):: da2kdx
    real(kind=longdouble_kind):: c2,c1,c0

    integer ::  kk

    c0 = 0.0_longdouble_kind
    c1 = 1.0_longdouble_kind
    c2 = 2.0_longdouble_kind

    jac (0)= c1
    jac (1)=(c1 + alpha)*xx
    djac(0)= c0
    djac(1)=(c1 + alpha)

    do kk=1,nn-1
      a1k       =  c2*(kk + c1)*(kk + alpha + beta + c1)*(c2*kk + alpha + beta)
      da2kdx    = (c2*(kk + c1) + alpha + beta)*(c2*kk + alpha + beta + c1) &
                                               *(c2*kk + alpha + beta)
      a2k       = (c2* kk + alpha + beta + c1)*(alpha*alpha - beta*beta) + xx*da2kdx
      a3k       =  c2*(kk + alpha)*(kk + beta)*(c2*kk + alpha + beta + c2)
      jac (kk+1)= (a2k* jac(kk)                  - a3k* jac(kk-1))/a1k
      djac(kk+1)= (a2k*djac(kk) + da2kdx*jac(kk) - a3k*djac(kk-1))/a1k          
    end do

    ! End Routine
    !-------------
    return
  end subroutine jacobi
  !==================================================================


  !==================================================================
  function jacobi_polynomials(nn,alpha,beta,npoints,xx) result(jac)
    ! jacobi_polynomials:
    !
    ! This routine computes the Nth order Jacobi Polynomials 
    ! (jac) for a vector of positions x on the interval (-1,1),
    ! of length npoints.
    !
    !    See for example the recurrence relations 
    !    in equation 2.5.4 (page 70) in 
    !
    !     "Spectral Methods in Fluid Dynamics",
    !     by C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A.Zang
    !     Springer-Verlag, 1988.
    !
    !================================================================
    !
    ! Passed Variables
    !---------------------
    integer,intent(in)        :: nn         ! order of the Jacobi Polynomial
    real(kind=longdouble_kind):: alpha 
    real(kind=longdouble_kind):: beta
    integer, intent(in)       :: npoints
    real(kind=longdouble_kind)::  xx(npoints)
    real(kind=longdouble_kind):: jac(npoints)
    !
    ! Local Values
    !--------------
    real(kind=longdouble_kind):: a1k
    real(kind=longdouble_kind):: a2k
    real(kind=longdouble_kind):: a3k
    real(kind=longdouble_kind):: da2kdx
    real(kind=longdouble_kind):: jacp1
    real(kind=longdouble_kind):: jacm1
    real(kind=longdouble_kind):: jac0
    real(kind=longdouble_kind):: xtmp
    real(kind=longdouble_kind):: c2,c1,c0

    integer jj,kk

    c0 = 0.0_longdouble_kind
    c1 = 1.0_longdouble_kind
    c2 = 2.0_longdouble_kind

    do jj=1,npoints
      xtmp =xx(jj)
      jacm1=c1
      jac0 =(c1+alpha)*xtmp
      do kk=1,nn-1
        a1k   = c2*(kk+c1)*(kk+alpha+beta+c1)*(c2*kk+alpha+beta)
        da2kdx=(c2*kk+alpha+beta+c2)*(c2*kk+alpha+beta+c1)*(c2*kk+alpha+beta)
        a2k   =(c2*kk+alpha+beta+c1)*(alpha*alpha-beta*beta) + xtmp*da2kdx
        a3k   = c2*(kk+alpha)*(kk+beta)*(c2*kk+alpha+beta+c2)
        jacp1 =(a2k*jac0-a3k*jacm1)/a1k
        jacm1 =jac0
        jac0  =jacp1
      end do
      if (nn==0)jac0=jacm1
      jac(jj)=jac0
    end do

    ! End Function
    !-------------
    return
  end function jacobi_polynomials
  !==================================================================


  !==================================================================
  function jacobi_derivatives(nn,alpha,beta,npoints,xx) result(djac)
    ! jacobi_derivatives:
    !
    ! This routine computes the first derivatives of Nth
    ! order Jacobi Polynomials (djac) for a vector of 
    ! positions x on the interval (-1,1), of length npoints.
    !
    ! See for example the recurrence relations 
    ! in equation 2.5.4 (page 70) in 
    !
    ! "Spectral Methods in Fluid Dynamics",
    ! by C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A.Zang
    ! Springer-Verlag, 1988.
    !
    !===============================================================
    !
    ! Passed Variables
    !-------------------
    integer                   ,intent(in):: nn         ! order of the Jacobi Polynomial
    real(kind=longdouble_kind),intent(in):: alpha 
    real(kind=longdouble_kind),intent(in):: beta
    integer                   ,intent(in):: npoints
    real(kind=longdouble_kind),intent(in)::   xx(npoints)
    real(kind=longdouble_kind)           :: djac(npoints)
    !
    ! Local Values
    !--------------
    real(kind=longdouble_kind):: a1k
    real(kind=longdouble_kind):: a2k
    real(kind=longdouble_kind):: a3k
    real(kind=longdouble_kind):: da2kdx
    real(kind=longdouble_kind):: jacp1
    real(kind=longdouble_kind):: jacm1
    real(kind=longdouble_kind):: jac0
    real(kind=longdouble_kind):: djacp1
    real(kind=longdouble_kind):: djacm1
    real(kind=longdouble_kind):: djac0
    real(kind=longdouble_kind):: xtmp
    real(kind=longdouble_kind):: c2,c1,c0

    integer jj,kk

    c0 = 0.0_longdouble_kind
    c1 = 1.0_longdouble_kind
    c2 = 2.0_longdouble_kind

    do jj=1,npoints
      xtmp  =xx(jj)
      jacm1 = c1
      jac0  =(c1+alpha)*xtmp
      djacm1= c0
      djac0 =(c1+alpha)
      do kk=1,nn-1
        a1k   = c2*(kk+c1)*(kk+alpha+beta+c1)*(c2*kk+alpha+beta)
        da2kdx=(c2*kk+alpha+beta+c2)*(c2*kk+alpha+beta+c1)*(c2*kk+alpha+beta)
        a2k   =(c2*kk+alpha+beta+c1)*(alpha*alpha-beta*beta) + xtmp*da2kdx
        a3k   = c2*(kk+alpha)*(kk+beta)*(c2*kk+alpha+beta+c2)
        jacp1 =(a2k*jac0-a3k*jacm1)/a1k
        djacp1=(a2k*djac0+da2kdx*jac0-a3k*djacm1)/a1k
        jacm1 =jac0
        jac0  =jacp1
        djacm1=djac0
        djac0 =djacp1
      end do
      if(nn==0) djac0=djacm1
      djac(jj)=djac0
    end do

    ! End Function
    !-------------
    return
  end function jacobi_derivatives
  !==================================================================


  !==================================================================
  function legendre(xx,nn) result(leg)
    ! legendre:
    !
    ! Compute the legendre polynomials using
    ! the recurrence relationship.
    ! return leg(m+1) = P_N(x) for m=0..N
    ! p_3 = Legendre polynomial of degree N
    ! p_2 = Legendre polynomial of degree N-1 at x
    ! p_1 = Legendre polynomial of degree N-2 at x
    !
    !===============================================================
    !
    ! Passed Variables
    !------------------
    real(kind=longdouble_kind):: xx
    integer                   :: nn
    real(kind=longdouble_kind):: leg(nn+1)
    !
    ! Local Values
    !-------------
    real(kind=longdouble_kind):: p_1, p_2, p_3

    integer kk

    p_3   = 1.0_longdouble_kind
    leg(1)=p_3
    if(nn.ne.0) then
      p_2   = p_3
      p_3   = xx
      leg(2)= p_3
      do kk=2,nn
        p_1      = p_2
        p_2      = p_3
        p_3      = ((2*kk-1)*xx*p_2 - (kk-1)*p_1)/kk
        leg(kk+1)= p_3
      end do
    endif

    ! End Function
    !-------------
    return
  end function legendre
  !==================================================================


  !==================================================================
  function quad_norm(gquad,nn) result(gamma)
    ! quad_norm: compute normalization constants for k=1,N order 
    !            Legendre polynomials
    !            e.g. gamma(k) in Canuto, page 58.
    !===============================================================
    !
    ! Passed Variables
    !------------------
    type(quadrature_t),intent(in):: gquad
    integer           ,intent(in):: nn
    real(kind=longdouble_kind)   :: gamma(nn)
    !
    ! Local Values
    !---------------
    real(kind=longdouble_kind):: leg(nn)

    integer ii,kk

    gamma(:)=0.0_longdouble_kind
    do ii=1,nn
      leg=legendre(gquad%points(ii),nn-1)
      do kk=1,nn
        gamma(kk)= gamma(kk)+leg(kk)*leg(kk)*gquad%weights(ii)
      end do
    end do

    ! End Function
    !-------------
    return
  end function quad_norm
  !==================================================================


  !==================================================================
  subroutine trapN(ff,aa,bb,nn,it,ss)
    ! TrapN: Numerical recipes
    !
    !============================================
    use SE_Constants,only: real_kind
    !
    ! Passed Variables
    !------------------
    INTERFACE
       FUNCTION ff(xx) RESULT(f_x)   ! Function to be integrated
         use SE_Constants,only: real_kind
         real(kind=real_kind),intent(in):: xx
         real(kind=real_kind)           :: f_x
       END FUNCTION ff
    END INTERFACE
    real(kind=real_kind),intent(in   ):: aa,bb   ! The integral bounds
    integer             ,intent(in   ):: nn
    integer             ,intent(inout):: it
    real(kind=real_kind),intent(inout):: ss
    !
    ! Local Values
    !-----------------
    real(kind=real_kind):: ssum
    real(kind=real_kind):: del
    real(kind=real_kind):: rtnm
    real(kind=real_kind):: xx

    integer jj

    if(nn==1) then
      ss = 0.5d0*(bb-aa)*(ff(aa) + ff(bb))
      it = 1
    else
      ssum = 0.0d0
      rtnm = 1.0d0/it
      del  = (bb-aa)*rtnm
      xx   = aa+0.5*del
      do jj=1,it
        ssum= ssum + ff(xx)
        xx  = xx + del
      end do
      ss = 0.5d0*(ss + del*ssum)
      it = 2*it  
    endif

    ! End Routine
    !-------------
    return
  end subroutine trapN
  !==================================================================


  !==================================================================
  function trapezoid(ff,aa,bb,eps) result(Integral)
    ! trapezoid:
    !
    ! Trapezoid Rule for integrating functions 
    ! from a to b with residual error eps
    !=========================================================
    use SE_Constants,only: real_kind
    integer,parameter:: NMAX = 25  ! At most 2^NMAX + 1 points in integral
    !
    ! Passed Variables
    !---------------------
    INTERFACE
       FUNCTION ff(xx) RESULT(f_x)   ! Function to be integrated
         use SE_Constants,only: real_kind
         real(kind=real_kind),intent(in):: xx
         real(kind=real_kind)           :: f_x
       END FUNCTION ff
    END INTERFACE
    real(kind=real_kind),intent(in):: aa,bb     ! The integral bounds
    real(kind=real_kind),intent(in):: eps       ! relative error bound for integral
    real(kind=real_kind)           :: Integral  ! the integral result (within eps)
    !
    ! Local Values
    !--------------
    real(kind=real_kind):: ss        ! Integral approximation
    real(kind=real_kind):: sold      ! previous integral approx

    integer nn,it

    ! Calculate I here using trapezoid rule using f and a DO loop...
    !-----------------------------------------------------------------
    ss  = 1.0D30
    sold= 0.0D0
    nn  = 1
    it  = 0
    do while((nn<=NMAX).and.(ABS(ss-sold)>eps*ABS(sold)))
      sold=ss
      call trapN(ff,aa,bb,nn,it,ss)
#ifdef _QUAD_DBG
      print *,"N=",nn," ABS(s-sold)",ABS(ss-sold)," threshold=",ABS(sold)*eps
#endif
      nn=nn+1
    end do
    Integral = ss

    ! End Function
    !-------------
    return
  end function trapezoid
  !==================================================================


  !==================================================================
  function simpsons(ff,aa,bb,eps) result(Integral)
    ! simpsons:
    !
    ! Simpsons Rule for integrating functions 
    ! from a to b with residual error eps
    !=============================================================
    use SE_Constants, only : real_kind
    integer,parameter:: NMAX = 25  ! At most 2^NMAX + 1 points in integral
    !
    ! Passed Variables
    !------------------
    INTERFACE
       FUNCTION ff(xx) RESULT(f_x)   ! Function to be integrated
         use SE_Constants,only: real_kind
         real(kind=real_kind),intent(in):: xx
         real(kind=real_kind)           :: f_x
       END FUNCTION ff
    END INTERFACE
    real(kind=real_kind),intent(in):: aa,bb     ! The integral bounds
    real(kind=real_kind),intent(in):: eps       ! relative error bound for integral
    real(kind=real_kind)           :: Integral  ! the integral result (within eps)
    !
    ! Local Values
    !---------------
    real(kind=real_kind):: ss        ! Integral approximation
    real(kind=real_kind):: os        ! previous integral approx
    real(kind=real_kind):: st        ! Integral approximation
    real(kind=real_kind):: ost       ! previous integral approx

    integer nn,it

    ! Calculate I here using trapezoid rule using f and a DO loop...
    !-----------------------------------------------------------------
    ost= 0.0D0
    ss = 1.0D30
    os = 0.0D0
    nn = 1
    it = 0
    do while(((nn<=NMAX).and.(ABS(ss-os)>eps*ABS(os))).or.nn<=2)
      os = ss
      call trapN(ff,aa,bb,nn,it,st)
      ss = (4.0D0*st-ost)/3.0D0
#ifdef _QUAD_DBG
      print *,"N=",nn," ABS(s-os)=",ABS(ss-os)," threshold=",ABS(os)*eps
#endif
      ost=st
      nn =nn+1
    end do
    Integral = ss

    ! End Function
    !-------------
    return
  end function simpsons
  !==================================================================


  !==================================================================
  function gaussian_int(ff,aa,bb,gs) result(Integral)
    ! gaussian_int: Gaussian Quadrature Rule for integrating function f 
    !               from a to b  with gs weights and points with 
    !               precomputed gaussian quadrature and weights.
    !===================================================================
    use SE_Constants, only : real_kind
    integer,parameter:: NMAX = 10  ! At most 2^NMAX + 1 points in integral
    !
    ! Passed Variables
    !------------------
    INTERFACE
       FUNCTION ff(xx) RESULT(f_x)   ! Function to be integrated
         use SE_Constants,only: real_kind
         real(kind=real_kind),intent(in):: xx
         real(kind=real_kind)           :: f_x
       END FUNCTION ff
    END INTERFACE
    real(kind=real_kind),intent(in):: aa,bb     ! The integral bounds
    type(quadrature_t)  ,intent(in):: gs        ! gaussian points/wts
    real(kind=real_kind)           :: Integral  ! the integral result (within eps)
    !
    ! Local Values
    !---------------
    real(kind=real_kind):: ss,xx

    integer ii

    ! Calculate I = S f(x)dx here using gaussian quadrature
    !------------------------------------------------------
    ss = 0.0D0
    do ii=1,SIZE(gs%points)
      xx = 0.50D0*((bb-aa)*gs%points(ii) + (bb+aa))
      ss = ss + gs%weights(ii)*ff(xx)
    end do
    Integral = ss*(0.5D0*(bb-aa))

    ! End Function
    !-------------
    return
  end function gaussian_int
  !==================================================================

end module quadrature_mod
