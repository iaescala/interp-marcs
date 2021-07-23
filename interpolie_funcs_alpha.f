****************************************************************************
* interpolation of model atmosphere
* parameter space of interpolation {Teff,logg,z}
* 8 MARCS binary models as input required:
* ! the order of input models matters !
* ! only standard composition has been tested, no guaranty for peculiar composition !
* turbospectrum/babsma format compatible output

* Interpolation scheme of model atmosphere(@)
* in {Teff,logg,z} stellar parameter space
*
*         |_ _ _ _ _
*         |        /|:
*        /|       /
*       _ |_ _ _ /  |
*       | |      |
*         |         |
*       | |    @ |
*         |--------------> logg
*       |/       | /
*       /_ _ _ _ ./
*      /
*     /
*    Teff
*
* Each structure component of each model (T,Pe,Pg,xit,optical depth, kappaross)
* is first resampled on a common optical depth basis tau(5000 or ross) (see resample routine).
*
* Then each component of the atmosphere(@) (T,Pe,Pg,xit,kappa,geom depth +tau5000 and tauross))
* is interpolated individually along each dimension between each input model value(*)
* (physically it is impossible to ensure a simple relationship between models in more than one dimension at the same time,
* that's why the input models MUST form a "cube" in the stellar parameter space {Teff,logg,z})
* The interpolation is successively done at each optical depth (tau)
* + interpolation point weighted by an empirical value (see TM thesis + manual)
*
*
*                 ^ T,Pe,Pg,xit or radius
*                 |
*                 |            *
*                 |           /|
*                 |          @
*                 |         /| |
*                 |   /    /
*                 |  /    /  | |
*         ^       | /    *
*         |       |/     |   | |
*         |       -----------------------> Teff, logg and z
*         |      /      low ref up
*         |     /*      /   / /
*         |    / |\
*         |   /    \  /   / /
*         |  /   |  \
*         | /       /@  / /
*         |/     |   |\
*         / _  _ /_ _ /*_ _ _ _ _
*        /      low ref up
*       /
*   tau{5000 or Ross}
*
***************************************************************************
* TM 07/2003

c  07/2004 resampling of each model on a common optical depth basis
c  06/2006 works for spherical geometry models
c  09/2007 new calibration of free parameter alpha
c           + modified to read Uppsala ascii models
c           + kappa interpolation
c           + rhox calculation
c           + 2 outputs (babsma and ATLAS/MOOG)
c  10/2007 error estimates
c  10/2011 two non crucial bugs fixed
c             -unformatted->formatted read for ascii models because there is a couple of trouble makers in the grid
c             -dimension of taubas was not matching tauR (emo@astro.indiana.edu)
c  04/2012 unformatted reading of ascii models reinstated (troublemakers hopefully fixed) /BE
****************************************************************************
!  compile with Fortran 90 or 95

c---------------------------------------------------------------------------------
      subroutine blend_103 (r,s,t,u,
     & x0000,x0010,x0100,x0110,x1000,x1010,x1100,x1110,
     & x0001,x0011,x0101,x0111,x1001,x1011,x1101,x1111, x )
!
!*******************************************************************************
!from http://www.math.iastate.edu/burkardt/f_src/
!
!! BLEND_103 extends scalar point data into a cube.
!
!
!  Diagram:
!
!    011--------------111
!      |               |
!      |               |
!      |               |
!      |               |
!      |               |
!    001--------------101
!
!
!      *---------------*
!      |               |
!      |               |
!      |      rst      |
!      |               |
!      |               |
!      *---------------*
!
!
!    010--------------110
!      |               |
!      |               |
!      |               |
!      |               |
!      |               |
!    000--------------100
!
!
!  Formula:
!
!    Written as a polynomial in R, S and T, the interpolation map has the
!    form:
!
!      X(R,S,T) =
!        1         * ( + x000 )
!      + r         * ( - x000 + x100 )
!      +     s     * ( - x000        + x010 )
!      +         t * ( - x000               + x001 )
!      + r * s     * ( + x000 - x100 - x010                       + x110 )
!      + r     * t * ( + x000 - x100        - x001        + x101 )
!      +     s * t * ( + x000        - x010 - x001 + x011 )
!      + r * s * t * ( - x000 + x100 + x010 + x001 - x011 - x101 - x110 + x111 )
!
!  Reference:
!
!    William Gordon,
!    Blending-Function Methods of Bivariate and Multivariate Interpolation
!      and Approximation,
!    SIAM Journal on Numerical Analysis,
!    Volume 8, Number 1, March 1971, pages 158-177.
!
!    William Gordon and Charles Hall,
!    Transfinite Element Methods: Blending-Function Interpolation over
!      Arbitrary Curved Element Domains,
!    Numerische Mathematik,
!    Volume 21, Number 1, 1973, pages 109-129.
!
!    William Gordon and Charles Hall,
!    Construction of Curvilinear Coordinate Systems and Application to
!      Mesh Generation,
!    International Journal of Numerical Methods in Engineering,
!    Volume 7, 1973, pages 461-477.
!
!    Joe Thompson, Bharat Soni, Nigel Weatherill,
!    Handbook of Grid Generation,
!    CRC Press, 1999.
!
!  Modified:
!
!    11 October 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R, S, T, the coordinates where an interpolated value
!    is desired.
!
!    Input, real X000, X001, X010, X011, X100, X101, X110, X111, the
!    data values at the corners.
!
!    Output, real X, the interpolated data value at (R,S,T).
!
      implicit none
!
      real r,s,t,u,x
      real x0000,x0010,x0100,x0110
      real x1000,x1010,x1100,x1110
      real x0001,x0011,x0101,x0111
      real x1001,x1011,x1101,x1111
!
!  Interpolate the interior point.
!
      x =
     & 1.0E+00         * ( + x0000 )
     & + r             * ( - x0000 + x1000 )
     & + s             * ( - x0000 + x0100 )
     & + t             * ( - x0000 + x0010 )
     & + u             * ( - x0000 + x0001 )
     & + r * s         * ( + x0000 - x1000 - x0100 + x1100 )
     & + r * t         * ( + x0000 - x1000 - x0010 + x1010 )
     & + r * u         * ( + x0000 - x1000 - x0001 + x1001 )
     & + s * t         * ( + x0000 - x0100 - x0010 + x0110 )
     & + s * u         * ( + x0000 - x0100 - x0001 + x0101 )
     & + t * u         * ( + x0000 - x0010 - x0001 + x0011 )
     & + r * s * t     * ( - x0000 + x1000 + x0100 + x0010
     &                             - x0110 - x1010 + x1100 + x1110 )
     & + r * t * u     * ( - x0000 + x1000 + x0010 + x0001
     &                             - x0011 - x1001 - x1010 + x1011 )
     & + r * s * u     * ( - x0000 + x1000 + x0100 + x0001
     &                             - x1100 - x1001 - x0101 + x1101 )
     & + s * t * u     * ( - x0000 + x0100 + x0010 + x0001
     &                             - x0110 - x0101 - x0011 + x0111 )
     & + r * s * t * u * ( - x0000 - x1000 - x0100 - x0010 - x0001
     &                             + x1100 + x1010 + x1001
     &                             + x0110 + x0101 + x0011
     &                             - x1110 - x1101 - x1011 - x0111
     &                             + x1111 )

      return
      end

c---------------------------------------------------------------------------------
      subroutine extract_bin(FILE,TEFF,grav,metal,ndepth,xlr_ref,tau5,
     &                tauR,temp,prese,presg,xit,rad,sph,kappa)

!extracted from osplot.f 07/2003, to get tau,T,Pe,Pg,microturb from a model
      implicit none
      integer :: ndp,ndepth,k,nlp,nlb
      parameter(ndp=200)
      CHARACTER*117 ADUM
      CHARACTER*256 FILE,file2
      CHARACTER*30 COMMENT
      real metal,grav,xlr_ref,TEFF,GG,radius
      logical :: sph
      real :: tau(ndp),t(ndp),z(ndp),ptot(ndp),prad(ndp),
     &  pg(ndp),pturb(ndp),pe(ndp),ro(ndp),xmass(ndp),xkapr(ndp)
      real :: gradp(ndp),gravity(ndp),pcheck(ndp),rr(ndp),
     &   xit(ndp),geff(ndp),gradptur(ndp),dp(ndp),taus(ndp),xlr(30),
     &   coldens(ndp)
      real :: tau5(ndp),tauR(ndp),temp(ndp),prese(ndp),
     &   presg(ndp),rad(ndp),kappa(ndp)
      real :: xlb(155000),w(155000),fluxme(155000)
      real :: presmo(30,ndp),ptio(ndp)
      real :: bPPR(NDP),bPPT(NDP),bPP(NDP),bGG(NDP),
     & bZZ(NDP),bDD(NDP),
     & bVV(NDP),bFFC(NDP),bPPE(NDP),bTT(NDP),
     & bTAULN(NDP),erad(ndp)
      integer :: NbTAU,IbTER
      common /struct/ tau,t,z,ptot,prad,pg,pturb,pe,ro,rr,taus,xlr,
     &                nlp,xkapr
      common /spectr/ nlb,xlb,w,fluxme
      common /pressure/ presmo,ptio
      common /binstruc/bPPR,bPPT,bPP,bGG,
     & bZZ,bDD,
     & bVV,bFFC,bPPE,bTT,
     & bTAULN,NbTAU,IbTER,erad
      common /radius/ radius
      OPEN(UNIT=10,FILE=FILE,STATUS='OLD',FORM='UNFORMATTED',
     &     convert='big_endian')
ccc     &     RECL=152600)
*
      CALL READMO(10,NDEPTH,TEFF,GG,metal,sph)
c         open(21,file=file2,status='unknown')
      do k=1,ndepth
        tau5(k)=log10(taus(k))
        tauR(k)=log10(tau(k))
         temp(k)=t(k)
         prese(k)=log10(pe(k))
         presg(k)= log10(pg(k))
         kappa(k)=log10(xkapr(k))
      end do
         xit=2.0
         xlr_ref=xlr(nlp)
         grav=log10(GG)
*        rad=rr
         if(.not.sph) radius=0.0
         rad=radius-z
c         write(21,1966) ndepth,xlr(nlp),log10(GG)
c1966     format('''INTERPOL''',1x,i3,f8.0,2x,f4.2,1x,'0 0.00')
c         do k=1,ndepth
c           write(21,1965) log10(taus(k)),t(k),log10(pe(k)),
c     &                    log10(pg(k)), xit
c1965       format(f8.4,1x,f8.2,3(x,f8.4))
c         enddo
c      close(21)
      close(10)
      END
C
      SUBROUTINE READMO(IARCH,JTAU,TEFF,G,metal,spherical)
C        THIS ROUTINE READS ONE MODEL, TO GET INFO ON PRAD
C             ( All features taken from listmo )
      PARAMETER (NDP=200)
C
      DIMENSION ABUND(16),TKORRM(NDP),FCORR(NDP),TAU(NDP),TAUS(NDP),
     *T(NDP),PE(NDP),PG(NDP),PRAD(NDP),PTURB(NDP),XKAPR(NDP),RO(NDP),
     *CP(NDP),CV(NDP),AGRAD(NDP),Q(NDP),U(NDP),V(NDP),ANCONV(NDP),
     *PRESMO(30,NDP),FCONV(NDP),RR(NDP),Z(NDP),EMU(NDP),HNIC(NDP)
     *,NJ(16),XLR(30),IEL(16),ANJON(16,5),PART(16,5),PROV(50,20+1),
     *ABSKA(50),SPRIDA(50),XLB(155000),FLUXME(155000),FLUMAG(155000),
     & PEP(16),
     * ABNAME(50),SOURCE(50),PTOT(NDP)
      DIMENSION W(155000),UW(12),BW(21),VW(25)
      CHARACTER*10 DAG,NAME,NAMEP,KLOCK
      CHARACTER*8 ABNAME,SOURCE
      DIMENSION WAVFLX(10)
      dimension PTIO(NDP)
      real*8 dluminosity
      real abSc,abTi,abV,abMn,abCo,metal
      logical spherical
      common /binstruc/ dummy(11*ndp+2),erad(ndp)
      common /struct/ tau,t,z,ptot,prad,pg,pturb,pe,ro,rr,taus,xlr,
     &                nlp,xkapr
      common /spectr/ nlb,xlb,w,fluxme
      common /pressure/ presmo,ptio
      common /radius/ radius
      DATA UW/0.145,0.436,0.910,1.385,1.843,2.126,2.305,2.241,1.270,
     *0.360,0.128,0.028/,BW/0.003,0.026,0.179,0.612,1.903,2.615,2.912,
     *3.005,2.990,2.876,2.681,2.388,2.058,1.725,1.416,1.135,0.840,0.568,
     *0.318,0.126,0.019/,VW/0.006,0.077,0.434,1.455,2.207,2.703,2.872,
     *2.738,2.505,2.219,1.890,1.567,1.233,0.918,0.680,0.474,0.312,0.200,
     *0.132,0.096,0.069,0.053,0.037,0.022,0.012/
      DATA NAME/'LOCAL'/,NAMEP/'PARSONS'/
      DATA A,B/.34785485,.65214515/
      IREAD=5
C
      REWIND IARCH
      READ(IARCH)
      erad=-1.e30
      READ(IARCH) INORD,DAG,KLOCK
      READ(IARCH) TEFF,FLUX,G,PALFA,PNY,PY,PBETA,ILINE,ISTRAL,
     &                MIHAL,IDRAB1,IDRAB2,IDRAB3,IDRAB4,IDRAB5,IDRAB6,
     &            ITMAX,NEL,(ABUND(I),I=1,NEL),abSc,abTi,abV,abMn,abCo
      GLOG=ALOG10(G)
      FNORD=0.1*INORD
C        CONVERT TO 'PHYSICAL FLUX'
      FLUX=3.14159*FLUX
      DO 2 I=1,NEL
    2 ABUND(I)=ALOG10(ABUND(I))+12.
      metal=abund(15)-7.50
      READ(IARCH)JTAU,NCORE,DIFLOG,TAUM,RADIUS,(RR(K),K=1,JTAU)
      if (jtau.gt.ndp) then
        print*, 'ERROR !!! Jtau (number of depths of model) = ',jtau
        print*, ' is larger than ndp!! Increase NDP.'
        stop
      endif
      if (radius.le.2.) then
         spherical = .false.
      else
         spherical = .true.
         rr=radius-rr
      endif

      READ(IARCH)JTAU,(TKORRM(I),I=1,JTAU),(FCORR(K),K=1,JTAU)
      NTPO=0
      DO 3 K=1,JTAU
        READ(IARCH) KR,TAU(K),TAUS(K),Z(K),T(K),PE(K),PG(K),PRAD(K),
     &              PTURB(K),XKAPR(K),RO(K),EMU(K),CP(K),CV(K),
     &              AGRAD(K),Q(K),U(K),V(K),ANCONV(K),HNIC(K),NMOL,
     &              (PRESMO(J,K),J=1,NMOL),ptio(k)
        TAUK=ALOG10(TAU(K))+10.01
        KTAU=TAUK
        IF(ABS(TAUK-KTAU).GT.0.02) GO TO 31
        IF(KTAU.EQ.10) K0=K
        NTPO=NTPO+1
   31   CONTINUE
    3 CONTINUE
c      Z0=Z(K0)
c      DO 5 I=1,JTAU
c        Z(I)=Z(I)-Z0
c        i1=min(i+1,jtau)
c        PTOT(I)=PG(I)+PRAD(I)+0.5*(pturb(i)+pturb(i1))
    5 CONTINUE

      READ(IARCH)(NJ(I),I=1,NEL),NLP,(XLR(I),I=1,NLP)
     & ,NPROV,NPROVA,NPROVS,(ABNAME(KP),SOURCE(KP),KP=1,NPROV)
c      DO 22 KTAU=1,NTPO
c      DO 20 IE=1,NEL
c      NJP=NJ(IE)
c      READ(IARCH) KR,TAUI,TI,PEI,IEL(IE),ABUND(IE),
c     &            (ANJON(IE,JJ),JJ=1,NJP),(PART(IE,JJ),JJ=1,NJP)
c   20 CONTINUE
c      DO 21 KLAM=1,NLP
c      READ(IARCH) KR,TAUIL,(PROV(J,KLAM),J=1,NPROV),
c     &            ABSKA(KLAM),SPRIDA(KLAM)
   21 CONTINUE
   22 continue
c      READ(IARCH) NLB,(XLB(J),FLUXME(J),J=1,NLB),(W(J),J=1,NLB)
C CONVERT TO 'PHYSICAL' FLUXES
c      DO 24 J=1,NLB
c   24 FLUXME(J)=3.14159*FLUXME(J)

c      dluminosity=0.
c      do 25 j=1,nlb
c       dluminosity=dluminosity+fluxme(j)*w(j)
25    continue
c      dluminosity=dluminosity*4.*3.14159*radius**2/3.82d33
c      ddddd=real(dluminosity)
c      print*,'luminosity: ',dluminosity*3.82d33,' erg/s  = ',ddddd,
c     &  ' solar luminosities'

      RETURN
         END

c---------------------------------------------------------------------------------
      subroutine extract_ascii(FILE,TEFF,grav,metal,alpha,ndepth,
     &xlr_ref,tau5,tauR,temp,prese,presg,xit,rad,sph,xkapr)
c     adapted from P. DeLaverny
      implicit none
      integer :: ndp,k
      parameter(ndp=200)
      integer :: imod,idum,ndepth,stat
      CHARACTER*117 ADUM
      CHARACTER*256 FILE,file2,line
      CHARACTER*30 COMMENT
      CHARACTER*50 MOCODE,blabla
      real metal,radius,mass,grav,xlr_ref,TEFF,GG,xic,alpha
      logical :: sph
      real :: tau(ndp),t(ndp),z(ndp),ptot(ndp),prad(ndp),
     &  pg(ndp),pturb(ndp),pe(ndp),ro(ndp),xmass(ndp),xkapr(ndp)
      real :: gradp(ndp),gravity(ndp),pcheck(ndp),rr(ndp),
     &  xit(ndp),geff(ndp),gradptur(ndp),dp(ndp),taus(ndp),xlr(30),
     &   coldens(ndp)
      real :: tau5(ndp),tauR(ndp),temp(ndp),prese(ndp),
     &   presg(ndp),rad(ndp),emu(ndp),vconv(ndp),fconv(ndp)
      real :: dimension xlb(155000),w(155000),fluxme(155000)

          blabla=''
          imod =10
          OPEN(UNIT=imod,FILE=FILE,STATUS='OLD')
          read(imod,'(a)') mocode
c          print*,mocode,' = mocode'
          if (mocode(1:1).eq.'p' .or. mocode(1:3).eq.'sun') then
            print*,' this model is PLANE PARALLEL'
          else if (mocode(1:1).eq.'s') then
            print*,' this model is SPHERICALLY SYMMETRIC'
          else
            print*,' This model may not be a NewMARCS model!'
          endif
          sph=(mocode(1:1).eq.'s')
          xlr_ref=5000
          read(imod,*)TEFF
          read(imod,*)
          read(imod,*)grav
          grav=log10(grav)
          read(imod,*)xic
          read(imod,*)mass
          read(imod,*)metal,alpha
          read(imod,*)radius
          do while (blabla.ne.'Model structure')
            read(imod,'(a)') blabla
          enddo
          backspace(imod)
          backspace(imod)
          read(imod,*)ndepth
          read(imod,*)
          read(imod,*)
          do k=1,ndepth
            read(imod,*) idum,tauR(k),tau5(k),rad(k),temp(k),
     &                   Pe(k),Pg(k)
            prese(k)  = log10(Pe(k))
            presg(k)  = log10(Pg(k))
            xit(k) = xic
            rad(k)=radius-rad(k)
          enddo
          read(imod,*)

          do k=1,ndepth

            read(imod,*,iostat=stat) idum,tauR(k),xkapr(k),ro(k),emu(k),
     &                   Vconv(k),Fconv(k)

c           Added by I.Escala to handle MARCS file format errors
            if (stat.eq.5010) then
              line=''
              backspace(imod)
              read(imod,'(a)') line
              read(line(39:48), '(f10.3)') Vconv(k)
              read(line(49:56), '(f8.5)') Fconv(k)
c              print*,tauR(k),xkapr(k),ro(k),emu(k),Vconv(k),Fconv(k)
            endif

            xkapr(k) = log10(xkapr(k))

          enddo

         close(imod)

      END

c---------------------------------------------------------------------------------

      subroutine common_depth_basis(tauresample,tau,nlinemod,nfile)
      implicit none
      integer :: file,nlinemod,nfile
      real,dimension(nlinemod,nfile) :: tau
      real,dimension(nlinemod) :: tauresample

      tauresample=0
c initialize the common tau(5000) with min depth = max of the min depth of the models
c                                  and max depth = min of the max depth of the models
c essential for  the resampling with cubic spline
      tauresample(1)=tau(1,1)
      tauresample(nlinemod)=tau(nlinemod,1)
      do file=2,nfile
         if (tauresample(1).lt.tau(1,file)) then
            tauresample(1)=tau(1,file)
          end if
         if (tauresample(nlinemod).gt.tau(nlinemod,file)) then
            tauresample(nlinemod)=tau(nlinemod,file)
         end if
      end do
      call blend_i_0d1 ( tauresample, nlinemod )
      end

!*******************************************************************************

      subroutine blend_i_0d1 ( x, m )
!
!
!! BLEND_I_0D1 extends indexed scalar data at endpoints along a line.
!http://orion.math.iastate.edu/burkardt/f_src/f_src.html
!
!  Diagram:
!
!    ( X1, ..., ..., ..., ..., ..., XM )
!
!  Reference:
!
!    William Gordon,
!    Blending-Function Methods of Bivariate and Multivariate Interpolation
!      and Approximation,
!    SIAM Journal on Numerical Analysis,
!    Volume 8, Number 1, March 1971, pages 158-177.
!
!    William Gordon and Charles Hall,
!    Transfinite Element Methods: Blending-Function Interpolation over
!      Arbitrary Curved Element Domains,
!    Numerische Mathematik,
!    Volume 21, Number 1, 1973, pages 109-129.
!
!    William Gordon and Charles Hall,
!    Construction of Curvilinear Coordinate Systems and Application to
!      Mesh Generation,
!    International Journal of Numerical Methods in Engineering,
!    Volume 7, 1973, pages 461-477.
!
!    Joe Thompson, Bharat Soni, Nigel Weatherill,
!    Handbook of Grid Generation,
!    CRC Press, 1999.
!
!  Modified:
!
!    15 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real X(M).
!
!    On input, X(1) and X(M) contain scalar values which are to be
!    interpolated through the entries X(2) through X(M).  It is assumed
!    that the dependence of the data is linear in the vector index I.
!
!    On output, X(2) through X(M-1) have been assigned interpolated
!    values.
!
!    Input, integer M, the number of entries in X.

      implicit none
!
      integer m
!
      integer i
      real r
      real x(m)
!
      do i = 2, m - 1

        r = real ( i - 1 ) / real ( m - 1 )

        call blend_101 ( r, x(1), x(m), x(i) )

      end do

      return
      end

!*******************************************************************************

      subroutine blend_101 ( r, x0, x1, x )
!
!
!! BLEND_101 extends scalar endpoint data to a line.
!http://orion.math.iastate.edu/burkardt/f_src/f_src.html
!
!  Diagram:
!
!    0-----r-----1
!
!  Reference:
!
!    William Gordon,
!    Blending-Function Methods of Bivariate and Multivariate Interpolation
!      and Approximation,
!    SIAM Journal on Numerical Analysis,
!    Volume 8, Number 1, March 1971, pages 158-177.
!
!    William Gordon and Charles Hall,
!    Transfinite Element Methods: Blending-Function Interpolation over
!      Arbitrary Curved Element Domains,
!    Numerische Mathematik,
!    Volume 21, Number 1, 1973, pages 109-129.
!
!    William Gordon and Charles Hall,
!    Construction of Curvilinear Coordinate Systems and Application to
!      Mesh Generation,
!    International Journal of Numerical Methods in Engineering,
!    Volume 7, 1973, pages 461-477.
!
!    Joe Thompson, Bharat Soni, Nigel Weatherill,
!    Handbook of Grid Generation,
!    CRC Press, 1999.
!
!  Modified:
!
!    14 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R, the coordinate where an interpolated value is desired.
!
!    Input, real X0, X1, the data values at the ends of the line.
!
!    Output, real X, the interpolated data value at (R).
!
      implicit none
!
      real r
      real x
      real x0
      real x1
!
      x = ( 1.0E+00 - r ) * x0 + r * x1

      return
      end

c*******************************************************************************

      subroutine calcrhox(tau,kappa,ndepth,rhox)
c     2 ways to calculate rhox : int(ro*dx) or int(1/kappa*dtau)
c      A&A 387, 595-604 (2002)
      implicit none
      integer :: i,ndepth
      real :: first
      real, dimension(ndepth) :: tau,kappa,rhox,f,x
      real :: tot
      real, external :: rinteg


      f=(1/10**kappa)
      x=10**(tau)
      first = x(1)*f(1)
      tot=rinteg(x,f,rhox,ndepth,first)
      do i=2,ndepth
           rhox(i) = rhox(i-1) + rhox(i)
      enddo
      end

      real function rinteg(x,f,fint,n,start)
c******************************************************************************
c     This routine is from ATLAS6
c******************************************************************************
      implicit none
      integer :: n,i,n1
      real x(5000), f(5000), fint(5000)
      real a(5000), b(5000), c(5000)
      real :: start

      call parcoe (f,x,a,b,c,n)
      fint(1) = start
      rinteg = start
      n1 = n - 1
      do 10 i=1,n1
         fint(i+1)= (a(i)+b(i)/2.*(x(i+1)+x(i))+
     .     c(i)/3.*((x(i+1)+x(i))*x(i+1)+x(i)*x(i)))*(x(i+1)-x(i))
10    rinteg = rinteg + fint(i+1)

      return
      end

c******************************************************************************

      subroutine parcoe(f,x,a,b,c,n)

      implicit none
      integer :: n,n1,j,j1
      real f(5000), x(5000), a(5000), b(5000), c(5000)
      real :: d,wt

      c(1)=0.
      b(1)=(f(2)-f(1))/(x(2)-x(1))
      a(1)=f(1)-x(1)*b(1)
      n1=n-1
      c(n)=0.
      b(n)=(f(n)-f(n1))/(x(n)-x(n1))
      a(n)=f(n)-x(n)*b(n)
      if(n.eq.2)return
      do 1 j=2,n1
      j1=j-1
      d=(f(j)-f(j1))/(x(j)-x(j1))
      c(j)=f(j+1)/((x(j+1)-x(j))*(x(j+1)-x(j1)))-f(j)/((x(j)-x(j1))*
     1(x(j+1)-x(j)))+f(j1)/((x(j)-x(j1))*(x(j+1)-x(j1)))
      b(j)=d-(x(j)+x(j1))*c(j)
    1 a(j)=f(j1)-x(j1)*d+x(j)*x(j1)*c(j)
      c(2)=0.
      b(2)=(f(3)-f(2))/(x(3)-x(2))
      a(2)=f(2)-x(2)*b(2)
      c(3)=0.
      b(3)=(f(4)-f(3))/(x(4)-x(3))
      a(3)=f(3)-x(3)*b(3)
      if(c(j).eq.0.)go to 2
      j1=j+1
      wt=abs(c(j1))/(abs(c(j1))+abs(c(j)))
      a(j)=a(j1)+wt*(a(j)-a(j1))
      b(j)=b(j1)+wt*(b(j)-b(j1))
      c(j)=c(j1)+wt*(c(j)-c(j1))
    2 continue
      a(n1)=a(n)
      b(n1)=b(n)
      c(n1)=c(n)
      return
      end

***********************************************************************************
      subroutine calc_error(xinf,xsup,yinf,ysup,zinf,zsup,teff_ref,
     &   logg_ref,z_ref)
      implicit none
      real :: xinf,xsup,yinf,ysup,zinf,zsup,teff_ref,logg_ref,z_ref
      real ::  error_T,error_Pe,error_Pg,error_kappa
      real :: errorTeffT,errorloggT,errorzT,
     &        errorTeffPe,errorloggPe,errorzPe,
     &         errorTeffPg,errorloggPg,errorzPg,
     &         errorTeffkappa,errorloggkappa,errorzkappa

! values read out of the figures of the manual and scaled down o the according step
              errorTeffT=0.055/32
              errorloggT=0.008/5
              errorzT=0.015/4
              errorTeffPe=0.65/32
              errorloggPe=0.4/5
              errorzPe=0.38/4
              errorTeffPg=0.25/32
              errorloggPg=0.23/5
              errorzPg=0.35/4
              errorTeffkappa=0.8/32
              errorloggkappa=0.36/5
              errorzkappa=0.38/4



      error_T=min((xsup-teff_ref),(teff_ref-xinf))/100*errorTeffT +
     &         min((ysup-logg_ref),(logg_ref-yinf))*errorloggT +
     &         min((zsup-z_ref),(z_ref-zinf))*errorzT
      write(*,1409) 'estimated max error on T =',error_T*100,'%'


      error_Pe=min((xsup-teff_ref),(teff_ref-xinf))/100*errorTeffPe +
     &         min((ysup-logg_ref),(logg_ref-yinf))*errorloggPe +
     &         min((zsup-z_ref),(z_ref-zinf))*errorzPe
       write(*,1409) 'estimated max error on Pe =',error_Pe*100,'%'


      error_Pg=min((xsup-teff_ref),(teff_ref-xinf))/100*errorTeffPg +
     &         min((ysup-logg_ref),(logg_ref-yinf))*errorloggPg +
     &         min((zsup-z_ref),(z_ref-zinf))*errorzPg
       write(*,1409) 'estimated max error on Pg =',error_Pg*100,'%'


      error_kappa=min((xsup-teff_ref),(teff_ref-xinf))/100
     &                                                *errorTeffkappa +
     &         min((ysup-logg_ref),(logg_ref-yinf))*errorloggkappa +
     &         min((zsup-z_ref),(z_ref-zinf))*errorzkappa
      write(*,1409) 'estimated max error on kappa =',error_kappa*100,'%'

 1409 format(a30,f5.1,a2)
      end

c*******************************************************************************

      REAL FUNCTION SevalSingle(u,x,y,n)

      REAL,INTENT(IN) :: u
      REAL,INTENT(IN),DIMENSION(:) :: x
      REAL,INTENT(IN),DIMENSION(:) :: y
      integer, intent(in) :: n

! for f2py wrapping output arguments cannot be of assumed SIZE or allocatable

c      REAL,DIMENSION(:),allocatable :: b,c,d
      real, dimension(n) :: b,c,d

      INTEGER, SAVE :: i=1
      INTEGER :: j, k
      REAL:: dx

c      INTERFACE
c       subroutine FMMsplineSingle(x, y, b, c, d, n)
c          REAL,INTENT(IN),DIMENSION(:) :: x
c          REAL,INTENT(IN),DIMENSION(:) :: y
c          integer, intent(in) :: n
c          REAL,INTENT(OUT),DIMENSION(:) :: b,c,d
c          real, intent(out), dimension(n) :: b,c,d
c        end subroutine FMMsplineSingle
c      END INTERFACE

c      n=SIZE(x)
c      allocate(b(n),c(n),d(n))
      call FMMsplineSingle(x, y, b, c, d, n)

      IF (  (i<1) .OR. (i >= n) ) i=1
      IF ( (u < x(i))  .OR.  (u >= x(i+1)) ) THEN
          i=1   ! binary search
          j=n+1

          DO
            k=(i+j)/2
            IF (u < x(k)) THEN
              j=k
            ELSE
              i=k
            END IF
            IF (j <= i+1) EXIT
          END DO
      END IF

      dx=u-x(i)   ! evaluate the spline
      SevalSingle=y(i)+dx*(b(i)+dx*(c(i)+dx*d(i)))

      RETURN
c      deallocate(b,c,d)
      END Function SevalSingle

c*******************************************************************************

      subroutine FMMsplineSingle(x, y, b, c, d, n)

c      REAL,DIMENSION(:), INTENT(IN)  :: x
c      REAL,DIMENSION(:), INTENT(IN)  :: y
c      REAL,DIMENSION(:), INTENT(OUT) :: b,c,d

      integer :: n
      real, dimension(n) :: x,y
      real, dimension(n) :: b,c,d

      REAL:: t,aux
      INTEGER:: k
      REAL,PARAMETER:: ZERO=0.0, TWO=2.0, THREE=3.0

c      n=SIZE(x)

       IF (n < 3) THEN
         b(1)=ZERO
         IF (n == 2) b(1)=(y(2)-y(1))/(x(2)-x(1))
         c(1)=ZERO
         d(1)=ZERO
         IF (n < 2) RETURN
         b(2)=b(1)
         c(2)=ZERO
         d(2)=ZERO
         RETURN
       END IF

        d(1)=x(2)-x(1)
        c(2)=(y(2)-y(1))/d(1)
       DO k=2,n-1
         d(k)=x(k+1)-x(k)
         b(k)=TWO*(d(k-1)+d(k))
         c(k+1)=(y(k+1)-y(k))/d(k)
         c(k)=c(k+1)-c(k)
       END DO

       b(1)=-d(1)
       b(n)=-d(n-1)
       c(1)=ZERO
       c(n)=ZERO
       IF (n > 3) THEN
         c(1)=c(3)/(x(4)-x(2))-c(2)/(x(3)-x(1))
         c(n)=c(n-1)/(x(n)-x(n-2))-c(n-2)/(x(n-1)-x(n-3))
         c(1)=c(1)*d(1)*d(1)/(x(4)-x(1))
         c(n)=-c(n)*d(n-1)*d(n-1)/(x(n)-x(n-3))
       END IF

       DO k=2,n
         t=d(k-1)/b(k-1)
         b(k)=b(k)-t*d(k-1)
         c(k)=c(k)-t*c(k-1)
       END DO

       c(n)=c(n)/b(n)
       DO k=n-1,1,-1
         c(k)=(c(k)-d(k)*c(k+1))/b(k)
       END DO

       b(n)=(y(n)-y(n-1))/d(n-1)+d(n-1)*(c(n-1)+c(n)+c(n))
       DO k=1,n-1
         b(k)=(y(k+1)-y(k))/d(k)-d(k)*(c(k+1)+c(k)+c(k))
         d(k)=(c(k+1)-c(k))/d(k)
         c(k)=THREE*c(k)
       END DO
       c(n)=THREE*c(n)
       d(n)=d(n-1)

       RETURN
       end subroutine FMMsplineSingle

c---------------------------------------------------------------------------------

      subroutine resample(taus,tauR,T,Pe,Pg,xit,rr,xkapref)
      implicit  none
      integer :: nlinemod,file,nfile,k,i
      real,dimension(:,:) :: taus,tauR,T,Pe,Pg,xit,rr,xkapref
      real,dimension(:,:),allocatable :: taubas
      real,dimension(:),allocatable :: tauresample,taustemp,tauRtemp,
     & Ttemp,Petemp,Pgtemp,xittemp,rrtemp,xkapreftemp

      INTERFACE
      function SevalSingle(u,x,y,n)
      REAL,INTENT(IN) :: u
      REAL,INTENT(IN),DIMENSION(:) :: x
      REAL,INTENT(IN),DIMENSION(:):: y
      integer, intent(in) :: n
      real SevalSingle
      end
      END INTERFACE

      nlinemod=size(taus,1)
      nfile=size(taus,2)-3

      allocate(tauresample(nlinemod),taustemp(nlinemod),
     &  tauRtemp(nlinemod),Ttemp(nlinemod)
     & ,Petemp(nlinemod),Pgtemp(nlinemod),xittemp(nlinemod),
     & rrtemp(nlinemod),xkapreftemp(nlinemod),
     &   taubas(nlinemod,size(taus,2)))

      taubas=tauR
      write(*,*) 'resample models on common depth basis: tauRoss'

      call common_depth_basis(tauresample,taubas,nlinemod,nfile)
c now do the resampling with the common tau
       do file=1,nfile
         do k=1,nlinemod

       taustemp(k)=SevalSingle(tauresample(k),taubas(:,file),
     &             taus(:,file), size(taubas(:,file)))

       tauRtemp(k)=SevalSingle(tauresample(k),taubas(:,file),
     &             tauR(:,file), size(taubas(:,file)))
       Ttemp(k)=SevalSingle(tauresample(k),taubas(:,file),T(:,file),
     &                      size(taubas(:,file)))
       Petemp(k)=SevalSingle(tauresample(k),taubas(:,file),Pe(:,file),
     &                      size(taubas(:,file)))
       Pgtemp(k)=SevalSingle(tauresample(k),taubas(:,file),Pg(:,file),
     &                      size(taubas(:,file)))
       xittemp(k)=SevalSingle(tauresample(k),taubas(:,file),xit(:,file),
     &                      size(taubas(:,file)))
       rrtemp(k)=SevalSingle(tauresample(k),taubas(:,file),rr(:,file),
     &                      size(taubas(:,file)))
       xkapreftemp(k)=SevalSingle(tauresample(k),taubas(:,file),
     &                 xkapref(:,file), size(taubas(:,file)))
         end do

         taus(:,file)=taustemp
         tauR(:,file)=tauRtemp
         T(:,file)=Ttemp
         Pe(:,file)=Petemp
         Pg(:,file)=Pgtemp
         xit(:,file)=xittemp
         rr(:,file)=rrtemp
         xkapref(:,file)=xkapreftemp
       end do

       deallocate(tauresample,taustemp,tauRtemp,Ttemp,Petemp
     &,Pgtemp,xittemp,rrtemp,xkapreftemp,taubas)

       end
