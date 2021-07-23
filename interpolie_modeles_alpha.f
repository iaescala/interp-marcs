      program interpol_modeles
c     Modified by I. Escala (Carnegie 2020-2021) to be compatible with
c     Python wrapping. Also modified to do 4D interpolation in alpha/Fe

      implicit none
      character*256 :: model1,model2,model3,model4,model5,
     &model6,model7,model8,model9,model10,model11,model12,
     &model13,model14,model15,model16,model1out, model2out
      real :: temp_ref, logg_ref, z_ref, a_ref
      logical :: optimize

c      model1 = 'Testwebformat/p5500:g+4.0:m0.0:t01:ap:z-0.50:'//
c     &'a+0.00:c+0.00:n+0.00:o+0.00:r+0.00:s+0.00.mod'
c      model2 = 'Testwebformat/p5500:g+4.0:m0.0:t01:st:z+0.00:'//
c     &'a+0.00:c+0.00:n+0.00:o+0.00:r+0.00:s+0.00.mod'
c      model3 = 'Testwebformat/p5500:g+5.0:m0.0:t01:ap:z-0.50:'//
c     &'a+0.00:c+0.00:n+0.00:o+0.00:r+0.00:s+0.00.mod'
c      model4 = 'Testwebformat/p5500:g+5.0:m0.0:t01:st:z+0.00:'//
c     &'a+0.00:c+0.00:n+0.00:o+0.00:r+0.00:s+0.00.mod'
c      model5 = 'Testwebformat/p6000:g+4.0:m0.0:t01:ap:z-0.50:'//
c     &'a+0.00:c+0.00:n+0.00:o+0.00:r+0.00:s+0.00.mod'
c      model6='Testwebformat/p6000:g+4.0:m0.0:t01:st:z+0.00:'//
c     &'a+0.00:c+0.00:n+0.00:o+0.00:r+0.00:s+0.00.mod'
c      model7='Testwebformat/p6000:g+5.0:m0.0:t01:ap:z-0.50:'//
c     &'a+0.00:c+0.00:n+0.00:o+0.00:r+0.00:s+0.00.mod'
c      model8='Testwebformat/p6000:g+5.0:m0.0:t01:st:z+0.00:'//
c     &'a+0.00:c+0.00:n+0.00:o+0.00:r+0.00:s+0.00.mod'
c      model1out='Testoutpy/p5750:g+4.5:m0.0:t01:ap:z-0.25:'//
c     &'a+0.00:c+0.00:n+0.00:o+0.00:r+0.00:s+0.00.mod'
c      model2out='Testoutpy/p5750:g+4.5:m0.0:t01:ap:z-0.25:'//
c     &'a+0.00:c+0.00:n+0.00:o+0.00:r+0.00:s+0.00.alt'

c      temp_ref = 5750
c      logg_ref = 4.5
c      z_ref = -0.25

      call interpolie(model1,model2,model3,model4,model5,model6,
     &model7,model8,model9,model10,model11,model12,model13,model14,
     &model15,model16,model1out,model2out,temp_ref,logg_ref,z_ref,
     &a_ref,optimize)
      end

c*****************************************************************************

      subroutine interpolie(model1,model2,model3,model4,model5,
     &model6,model7,model8,model9,model10,model11,model12,model13,
     &model14,model15,model16,model1out,model2out,temp_ref,logg_ref,
     &z_ref,a_ref,optimize)

c      include 'interpolie_funcs.f'

      implicit none
      integer :: file,nfile,k,ndp,ndepth_ref,out,nlinemod
      parameter (ndp=200)
      parameter (nfile=18)
      logical :: verif,check,test,extrapol,binary,optimize
      real :: lambda_ref,temp_ref,logg_ref,z_ref,a_ref,x,y,z,w,xinf,
     &        xsup,yinf,ysup,zinf,zsup,winf,wsup,teffpoint,loggpoint,
     &        metpoint,alphapoint
      character*256 :: model1,model2,model3,model4,model5,model6,
     &model7,model8,model9,model10,model11,model12,model13,model14,
     &model15,model16,model1out,model2out
      character*256, dimension (nfile) :: FILE_IN
      real, dimension (:,:), allocatable:: taus,tauR,T,Pe,Pg,xit,rr,
     &xkapref,rhox,taus_aux,tauR_aux,T_aux,Pe_aux,Pg_aux,xit_aux,rr_aux,
     &xkapref_aux
      integer, dimension (nfile) :: ndepth
      real, dimension (nfile) :: xlr,teff,logg,metal,alpha
      logical, dimension (nfile) :: sph
      external :: blend_103
      real, external :: inf,sup
      real, dimension(8,3) :: lin_dif,power

Cf2py intent(in) model1,model2,model3,model4,model5,model6,model7,model8,model9,model10,model11,model12,model13,model14,model15,model16,model1out,model2out,temp_ref,logg_ref,z_ref,a_ref,optimize

      INTERFACE reec
        subroutine resample(taus,tauR,T,Pe,Pg,xit,rr,xkapref)
        real, dimension (:,:) :: taus,tauR,T,Pe,Pg,xit,rr,xkapref
        end
      END INTERFACE reec

      FILE_IN = (/model1,model2,model3,model4,model5,model6,
     &model7,model8,model9,model10,model11,model12,model13,
     &model14,model15,model16,model1out,model2out/)

      write(*,*) '*****************************'
      write(*,*) '* begining of interpolation *'
      write(*,*) '*****************************'

******* you can choose here to switch of the "optimization" and prefer simple linear interpolation

c      optimize = .true.

******  read 8 models, put in tables,
****** check number of layer, reference optical depth, and geometry compatibility ******

      out=17
      write(*,*) 'Interpolation between :'
      do file=1,18
         write(*,*) FILE_IN(file)
      end do

      test=.false.
      verif=.true.
      check=.true.
      nlinemod=ndp

      allocate(taus_aux(nlinemod,nfile),tauR_aux(nlinemod,nfile),
     & T_aux(nlinemod,nfile),Pe_aux(nlinemod,nfile),
     & Pg_aux(nlinemod,nfile),xit_aux(nlinemod,nfile),
     & rr_aux(nlinemod,nfile),xkapref_aux(nlinemod,nfile))

      binary=.false.
      if (binary) then
      do file=1,16
      call extract_bin(FILE_IN(file),teff(file),logg(file),metal(file),
     & ndepth(file),xlr(file),taus_aux(:,file),tauR_aux(:,file),
     & T_aux(:,file),Pe_aux(:,file),Pg_aux(:,file),xit_aux(:,file),
     & rr_aux(:,file),sph(file),xkapref_aux(:,file))
      verif=verif.and.(ndepth(file).eq.ndepth(1))
      check=check.and.(xlr(file).eq.xlr(1))
      if (.not.(((sph(1).and.sph(file))).or.
     &    ((.not.(sph(1))).and.(.not.(sph(file)))))) then
      write(*,*) 'geometry compatibility problem with'
       write(*,78) file,teff(file),logg(file),metal(file)
      stop
      endif
      write(*,78) file,teff(file),logg(file),metal(file)
      end do

      else
          do file=1,16
      call extract_ascii(FILE_IN(file),teff(file),logg(file),
     & metal(file),alpha(file),ndepth(file),xlr(file),taus_aux(:,file),
     & tauR_aux(:,file),T_aux(:,file),Pe_aux(:,file),Pg_aux(:,file),
     & xit_aux(:,file),rr_aux(:,file),sph(file),xkapref_aux(:,file))
      verif=verif.and.(ndepth(file).eq.ndepth(1))
      check=check.and.(xlr(file).eq.xlr(1))
      if (.not.(((sph(1).and.sph(file))).or.
     &    ((.not.(sph(1))).and.(.not.(sph(file)))))) then
      write(*,*) 'geometry compatibility problem with'
       write(*,78) file,teff(file),logg(file),metal(file),alpha(file)
      stop
      end if
      write(*,78) file,teff(file),logg(file),metal(file),alpha(file)
      end do
      endif
 78   format('model',i2,'  Teff=',f8.0,'  logg=',f5.2,'  z=',f6.2,
     &       ' a=',f6.2)


      write(*,79) temp_ref,logg_ref,z_ref,a_ref
 79   format('Interpolation point : Teff=',f8.0,'  logg=',f5.2,
     &                            '  z=',f6.2,' a=',f6.2)
************** check if files are length and depth ref compatible *******
      if (.not.(check)) then
         write(*,*) 'All the models do not have the same'
          write(*,*) 'lambda ref'
          write(*,*) 'no interpolation done'
          stop

       else
      if (.not.(verif)) then
         write(*,*) 'WARNING : All the models do not have the same'
         write(*,*) 'number of layers, resampling to',
     &                 ndepth(1),'layers'
      end if
      ndepth_ref=ndepth(1)
      lambda_ref=xlr(1)
********* calculation of the interpolation point(x,y,z) in {Teff,logg,z} space******************

       allocate(taus(ndepth_ref,nfile),tauR(ndepth_ref,nfile),
     & T(ndepth_ref,nfile),Pe(ndepth_ref,nfile),Pg(ndepth_ref,nfile),
     & xit(ndepth_ref,nfile),rr(ndepth_ref,nfile),
     &  xkapref(ndepth_ref,nfile))
       taus=taus_aux(1:ndepth_ref,:)
       tauR=tauR_aux(1:ndepth_ref,:)
       T=T_aux(1:ndepth_ref,:)
       Pe=Pe_aux(1:ndepth_ref,:)
       Pg=Pg_aux(1:ndepth_ref,:)
       xit=xit_aux(1:ndepth_ref,:)
       rr=rr_aux(1:ndepth_ref,:)
       xkapref=xkapref_aux(1:ndepth_ref,:)

         xinf=inf(teff)
         yinf=inf(logg)
         zinf=inf(metal)
         winf=inf(alpha)
         xsup=sup(teff)
         ysup=sup(logg)
         zsup=sup(metal)
         wsup=sup(alpha)
         if (xsup.eq.xinf) then
            teffpoint=0
         else
          teffpoint=(temp_ref-xinf)/(xsup-xinf)
         end if
         if (ysup.eq.yinf) then
            loggpoint=0
         else
          loggpoint=(logg_ref-yinf)/(ysup-yinf)
         end if
         if (zsup.eq.zinf) then
            metpoint=0
         else
          metpoint=(z_ref-zinf)/(zsup-zinf)
         end if
         if (wsup.eq.winf) then
            alphapoint=0
         else
            alphapoint=(a_ref-winf)/(wsup-winf)
         endif
         extrapol=((teffpoint.lt.0).or.(teffpoint.gt.1)
     &         .or.(loggpoint.lt.0).or.(loggpoint.gt.1)
     &         .or.(metpoint.lt.0).or.(metpoint.gt.1)
     &         .or.(alphapoint.lt.0).or.(alphapoint.gt.1))
         if (extrapol) then
         write(*,*) '!!!  WARNING : extrapolation  !!!'
         end if

*
*******resample each layer of each input model on a common depth basis(tau5000 or tauRoss, see resample routine)*****************
!if you don't want to resample all the model to the same depth scale, just comment the following line
         call resample(taus,tauR,T,Pe,Pg,xit,rr,xkapref)


****** initialisation of empirical constants for optimized interpolation (see TM thesis)*************

         lin_dif(1,1)=0                          !tau5000 vs Teff
         lin_dif(1,2)=0                          ! ...    vs logg
         lin_dif(1,3)=0                          ! ...    vs z
         lin_dif(2,1)=0                          !tauross vs Teff
         lin_dif(2,2)=0                          ! ...    vs logg
         lin_dif(2,3)=0                          ! ...    vs z
         lin_dif(3,1)=0.15                       !T       vs Teff
         lin_dif(3,2)=0.3                        ! ...    vs logg
         lin_dif(3,3)=1-(temp_ref/4000)**2.0     ! ...    vs z
         lin_dif(4,1)=0.15                       !logPe   vs Teff
         lin_dif(4,2)=0.06                       ! ...    vs logg
         lin_dif(4,3)=1-(temp_ref/3500)**2.5     ! ...    vs z
         lin_dif(5,1)=-0.4                       !logPg   vs Teff
         lin_dif(5,2)=0.06                       ! ...    vs logg
         lin_dif(5,3)=1-(temp_ref/4100)**4       ! ...    vs z
         lin_dif(6,1)=0                          !xit     vs Teff
         lin_dif(6,2)=0                          ! ...    vs logg
         lin_dif(6,3)=0                          ! ...    vs z
         lin_dif(7,1)=0                          !rr      vs Teff
         lin_dif(7,2)=0                          ! ...    vs logg
         lin_dif(7,3)=0                          ! ...    vs z
         lin_dif(8,1)=-0.15                      !logxkapref vs Teff
         lin_dif(8,2)=-0.12                      ! ...    vs logg
         lin_dif(8,3)=1-(temp_ref/3700)**3.5     ! ...    vs z

         if (optimize) then
          write(*,*) 'optimized interpolation applied for standard compo
     &sition models'
         else
            lin_dif=0.
             write(*,*) 'linear interpolation applied'
         end if
!these constants are calibrated on a broad range of stellar parameters; scale them now to the present one.
            power(:,1)= 1-(lin_dif(:,1)*(abs(xsup-xinf)/(7000-3800)))
            power(:,2)= 1-(lin_dif(:,2)*(abs(ysup-yinf)/(5-0.0)))
            power(:,3)= 1-(lin_dif(:,3)*(abs(zsup-zinf)/(0-(-4))))

****** interpolation of each component of the atmosphere (taus,teff,Pe,Pg,microt,rr) and at each layer *****************

        do k=1,ndepth_ref
          x=(teffpoint)**power(1,1)
          y=(loggpoint)**power(1,2)
          z=(metpoint)**power(1,3)
          w=alphapoint
          call blend_103(x,y,z,w,taus(k,1),taus(k,2),
     &     taus(k,3),taus(k,4),taus(k,5),taus(k,6),taus(k,7),taus(k,8)
     &     ,taus(k,9),taus(k,10),taus(k,11),taus(k,12),taus(k,13)
     &     ,taus(k,14),taus(k,15),taus(k,16),taus(k,out))

          x=(teffpoint)**power(2,1)
          y=(loggpoint)**power(2,2)
          z=(metpoint)**power(2,3)
          w=alphapoint
          call blend_103(x,y,z,w,tauR(k,1),tauR(k,2),
     &     tauR(k,3),tauR(k,4),tauR(k,5),tauR(k,6),tauR(k,7),tauR(k,8)
     &     ,tauR(k,9),tauR(k,10),tauR(k,11),tauR(k,12),tauR(k,13)
     &     ,tauR(k,14),tauR(k,15),tauR(k,16),tauR(k,out))

          x=(teffpoint)**power(3,1)
          y=(loggpoint)**power(3,2)
          z=(metpoint)**power(3,3)
          w=alphapoint
          call blend_103(x,y,z,w,T(k,1),T(k,2),T(k,3),T(k,4)
     &     ,T(k,5),T(k,6),T(k,7),T(k,8),T(k,9),T(k,10)
     &     ,T(k,11),T(k,12),T(k,13),T(k,14),T(k,15)
     &     ,T(k,16),T(k,out))

          x=(teffpoint)**power(4,1)
          y=(loggpoint)**power(4,2)
          z=(metpoint)**power(4,3)
          w=alphapoint
          call blend_103(x,y,z,w,Pe(k,1),Pe(k,2),Pe(k,3),Pe(k,4)
     &     ,Pe(k,5),Pe(k,6),Pe(k,7),Pe(k,8),Pe(k,9),Pe(k,10)
     &     ,Pe(k,11),Pe(k,12),Pe(k,13),Pe(k,14),Pe(k,15),Pe(k,16)
     &     ,Pe(k,out))

          x=(teffpoint)**power(5,1)
          y=(loggpoint)**power(5,2)
          z=(metpoint)**power(5,3)
          w=alphapoint
          call blend_103(x,y,z,w,Pg(k,1),Pg(k,2),Pg(k,3),Pg(k,4)
     &     ,Pg(k,5),Pg(k,6),Pg(k,7),Pg(k,8),Pg(k,9),Pg(k,10)
     &     ,Pg(k,11),Pg(k,12),Pg(k,13),Pg(k,14),Pg(k,15)
     &     ,Pg(k,16),Pg(k,out))

          x=(teffpoint)**power(6,1)
          y=(loggpoint)**power(6,2)
          z=(metpoint)**power(6,3)
          w=alphapoint
          call blend_103(x,y,z,w,xit(k,1),xit(k,2),
     &     xit(k,3),xit(k,4),xit(k,5),xit(k,6),xit(k,7),xit(k,8)
     &     ,xit(k,9),xit(k,10),xit(k,11),xit(k,12),xit(k,13)
     &     ,xit(k,14),xit(k,15),xit(k,16),xit(k,out))

          x=(teffpoint)**power(7,1)
          y=(loggpoint)**power(7,2)
          z=(metpoint)**power(7,3)
          w=alphapoint
          call blend_103(x,y,z,w,rr(k,1),rr(k,2),
     &     rr(k,3),rr(k,4),rr(k,5),rr(k,6),rr(k,7),
     &     rr(k,8),rr(k,9),rr(k,10),rr(k,11),rr(k,12),
     &     rr(k,13),rr(k,14),rr(k,15),rr(k,16),
     &     rr(k,out))

          x=(teffpoint)**power(8,1)
          y=(loggpoint)**power(8,2)
          z=(metpoint)**power(8,3)
          w=alphapoint
          call blend_103(x,y,z,w,xkapref(k,1),xkapref(k,2),
     &     xkapref(k,3),xkapref(k,4),xkapref(k,5),xkapref(k,6),
     &     xkapref(k,7),xkapref(k,8),xkapref(k,9),xkapref(k,10),
     &     xkapref(k,11),xkapref(k,12),xkapref(k,13),
     &     xkapref(k,14),xkapref(k,15),xkapref(k,16),
     &     xkapref(k,out))
       end do
       ndepth(out)=ndepth_ref
       xlr(out)=lambda_ref


**********now calculate rhox*****************
      write(*,*) 'now calculate rhox'
      allocate(rhox(ndepth_ref,nfile))
      do file=1,17
      call calcrhox(tauR(:,file),xkapref(:,file),ndepth_ref,
     &                                                 rhox(:,file))
      enddo

**********calculate estimated error********
      if (optimize) then
         write(*,*) 'now calculate error'
      call calc_error(xinf,xsup,yinf,ysup,zinf,zsup,temp_ref,
     & logg_ref,z_ref)
      endif

******** write interpolated model in file nber out (basma compatible format)***********
      write(*,*) 'now write result'
       open(unit=23,file=FILE_IN(out))
       open(unit=25,file=FILE_IN(out+1))

       if (sph(1)) then
        write(*,*) 'spherical models'
        write(23,1967) ndepth_ref,xlr(out),logg_ref
 1967   format('''sphINTERPOL''',1x,i3,f8.0,2x,f4.2,1x,'0 0.00')
         do k=1,ndepth_ref
           write(23,1968) taus(k,out),T(k,out),Pe(k,out),
     &                    Pg(k,out),xit(k,out),rr(k,out),taur(k,out)
         enddo
 1968   format(f8.4,1x,f8.2,3(1x,f8.4),1x,e15.6,1x,f8.4)
        write(25,19671) ndepth_ref,xlr(out)
19671   format('sphINTERPOL',1x,i3,f8.0)
        write(25,19672)
19672   format('    k    log(tau)  T    log(Pe)   log(Pg)    rhox')
         do k=1,ndepth_ref
           write(25,19681) k,taus(k,out),T(k,out),Pe(k,out),Pg(k,out),
     &     rhox(k,out)
19681      format(i5,1x,f8.4,1x,f8.2,2(1x,f8.4),1x,e15.6)
          enddo



        else
        write(*,*) 'plane parallel models'
        write(23,1966) ndepth_ref,xlr(out),logg_ref
1966     format('''ppINTERPOL''',1x,i3,f8.0,2x,f4.2,1x,'0 0.00')
         do k=1,ndepth_ref
           write(23,1965) taus(k,out),T(k,out),Pe(k,out),
     &                    Pg(k,out),xit(k,out),rr(k,out),taur(k,out)
1965       format(f8.4,1x,f8.2,3(1x,f8.4),1x,e15.6,1x,f8.4)
          enddo
        write(25,19661) xlr(out)
19661   format('ppINTERPOL',f8.0)
        write(25,19672)
         do k=1,ndepth_ref
           write(25,19681) k,taus(k,out),T(k,out),Pe(k,out),Pg(k,out),
     &                   rhox(k,out)
         enddo
       end if


       write(23,*) (FILE_IN(file),file=1,16)
       write(25,*) (FILE_IN(file),file=1,16)

       close(23)
       close(25)

      if (extrapol) then
          write (*,*) 'extrapolation done'
          else
      write (*,*) 'interpolation done'
      end if
      end if

      deallocate(taus,tauR,T,Pe,Pg,xit,taus_aux,tauR_aux,T_aux,Pe_aux,
     & Pg_aux,xit_aux,rr_aux,rr,rhox)

      end

c*****************************************************************************

      real function inf(tab)
      implicit none
      integer :: n
      real,dimension(17) :: tab
      inf=tab(1)
      do n=2,16
         if (tab(n).lt.inf) then
            inf=tab(n)
         end if
      end do
      end function

      real function sup(tab)
      implicit none
      integer :: n
      real,dimension(17) :: tab
      sup=tab(1)
      do n=2,16
         if (tab(n).gt.sup) then
            sup=tab(n)
         end if
      end do
      end function
