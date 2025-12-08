!     Last change:  MG   05 Dec 2025    
PROGRAM mc

  ! Purpose:
  ! --------
  !   To simulate the light field in an aquatic environment.
  !
  ! Description
  ! -----------
  !   This is the main program. It calls various subroutines
  !   and functions that simulate the fate of the photon during
  !   its residence time in the aquatic system.                        
  !
  ! Record and revisions
  ! --------------------
  !     Date     Programmer           Description of change
  !     ====     ==========           =====================
  !  3/26/2003   Manuel Gimond        Original code
  !  11/20/2025  Manuel Gimond		    Added routine to generate a non-surfer formatted grid output file
  !  12/02/2025  Manuel Gimond        Added missing PHI column in the rad.out output
  !  12/02/2025  Manuel Gimond        Allow for more than 4 PHI columns in rad.out output
  !  12/05/2025  Manuel Gimond        Update license
  !  12/05/2025  Manuel Gimond        Replaced custom RNG with intrinsic RNG in modules.f90
  !  12/07/2025  Manuel Gimond		  Change recording of output aop.out and rad.out
  !
  ! License/Disclaimer
  ! ------------------
  ! The AOMC source code is open source and available under the MIT License.
  ! Copyright (c) 2003-2025 Manuel Gimond
  !
  ! THE AUTHOR MAKES NO WARRANTY (EXPRESSED OR IMPLIED) AS TO THE
  ! USEFULNESS OF THE SOFTWARE AND ITS DOCUMENTION AND ASSUMES NO
  ! RESPONSIBILITY FOR THE (MIS)USE OF THE SOFTWARE.
  !
  !                      =====
  !                       TOC
  !                      =====
  !
  ! 1.0  Initialization
  ! 2.0  Reading general parameters from 'amc.inp'
  !      2.1 Post initialization
  ! 3.0  Allocation of memory
  ! 4.0  Reading of input files
  !      4.1  reading of 'conc.inp'
  !      4.2  reading of 'abs.inp'
  !      4.3  reading of 'scat.inp'
  !      4.4  Deallocation of some memory
  !      4.5  Reading of 'spf.inp'
  !      4.6  Reading of 'lambbot.inp'
  !      4.7  Reading of 'difcol.inp'
  ! 5.0  Calculation of IOP
  !      5.1  Beam attenuation coefficient
  !      5.2  Scattering albedo
  !      5.3  Write IOP's to file 'iop.out'
  !    * 5.4  Write IOP's to database file(s)
  ! 6.0  Run the model
  !      6.1  Free up some memory no longer needed
  ! 7.0  Compute AOP
  !      7.1  Allocate arrays
  !      7.2  Calculate radiance (L)
  !      7.3  Calculate planar irradiances (Eu, Ed)
  !      7.4  Caclulate scalar irradiances (Eou, Eod)
  !      7.5  Output AOP's to file
  ! 8.0  Free up remaining allocated memory
  !   
  !


  USE rand_global 
  USE randmod
  USE const_global
  USE physical_global
  USE log_global
  USE light_global
  USE math_global
  USE propagation_global
  USE aop_global

  IMPLICIT NONE

  LOGICAL   :: erro           !Error of 0 indicates no error
  CHARACTER(LEN=24) :: fmt    !On the fly format string
  CHARACTER(LEN=20) :: fname  !On the fly filename creation
  CHARACTER(LEN=5)  :: cdepth !Conversion of NUM to CHAR
  CHARACTER(LEN=7)  :: clamb  !Conversion of NUM to CHAR
  CHARACTER(LEN=15)  :: cformat !Conversion of NUM to CHAR
  CHARACTER(LEN=64) :: rowfmt
  CHARACTER(LEN=64) :: headfmt
  INTEGER, PARAMETER :: p8i = SELECTED_INT_KIND(14)    ! Allows for values greater than 32 bit
  INTEGER   :: ios,records
  INTEGER   :: ii,jj,kk,ll !Loop index
  INTEGER   :: astat       !Status of 0 indicates successful alloc/dealloc
  INTEGER   :: iter        !Number of photons to simulate (from user input)
  INTEGER(p8i):: titer     !True number of photons to simulate (different
  INTEGER   :: lam         !Wavelength counter used in loops
  !from iterations if tabular data is used)
  INTEGER   :: botout      !Log to file photons absorbed by the bottom
  INTEGER   :: absout      !Log to file photons absorbed within water column
  INTEGER   :: sideout     !Log to file photons absorbed by side boundaries
  INTEGER   :: normh2o     !Normalize to below water (1) or above water (0)
  INTEGER   :: spfang      !Number of scaterring phase function angles to read in
  INTEGER   :: surfer      !output to surfer format
  INTEGER   :: numskyelem  !number of sky elements
  INTEGER   :: numskywave  !number of wavelengths used
  INTEGER   :: action      !Type of interaction which ended iteration
  INTEGER   :: sidebound   !number of grids cells about the origin
  INTEGER   :: total       !Incremental summation placeholder
  INTEGER   :: wavematch   !Wavelength of sky-element data that most closely matches
  INTEGER   :: wavecontrol !Swith to indicate if a new set of sky dist. is used
  INTEGER   :: switch      !Indicator for transmittance or reflectance across layer
                           !that of the simulation
  INTEGER, PARAMETER :: p8r = SELECTED_INT_KIND(18)!Allows for values greater than 32 bit
  REAL(p8r),ALLOCATABLE,DIMENSION(:) :: cositer    !Cosine of angle of impinging
                                                   !light
  REAL,ALLOCATABLE,DIMENSION(:,:)::conc     !Concentration of each constituent in 
                                            !each layer
  REAL,ALLOCATABLE,DIMENSION(:,:)::specabs  !Specific absorption coefficient 
                                            !of each constituent
  REAL,ALLOCATABLE,DIMENSION(:,:)::specscat !Specific scattering coefficient of 
                                            !each constituent
  REAL,ALLOCATABLE,DIMENSION(:)  :: norma   !Value used to normalize the irradiance
                                            ! 'quartet'
  REAL,ALLOCATABLE,DIMENSION(:)::wavelength !Wavelength associated with photon
  REAL,ALLOCATABLE,DIMENSION(:,:,:,:):: polar_ld, polar_lu

  REAL      :: scratch     !used to read unused data
  REAL      :: randradius  !Randomly generated radius
  REAL      :: randangle   !Randomly generated angle
  REAL      :: sumfrac     !cumulitive summary of cumputed fracscat()
  REAL      :: ai          !Angle in rad
  REAL      :: sr          !Solid angle
  REAL      :: R_aop, mu_up_aop, mu_down_aop, mu_aop, rd_aop, rup_aop, pld_aop, plu_aop
  REAL      :: ed_sum_aop, eu_sum_aop, eou_sum_aop, eod_sum_aop, efwd_plus_eback_aop
  REAL      :: R_val, Lu_val, Lua_val, ed_sum_0, norma_val
  REAL, ALLOCATABLE, DIMENSION(:) :: rad_vals_row
  REAL      :: dmu


  WRITE(*,'(7(A60,/))')                                             &
       '===================================================',&
       '       Aquatic Optics Monte Carlo Model            ',&
       '       A       O      M     C                      ',&
       '                                                   ',&
       '       Version 1.2                                 ',&
       '                                                   ',&
       '       Original Author: Manuel Gimond              ',&
       '==================================================='


  ! ----------------------------------------------------------
  !                          - 1.0 -
  !                       Initialization
  ! ----------------------------------------------------------

  ! ----------------------------------------------------------
  !     Initialize the seeds to be used with the random number
  !     generator. This only needs to be done once for this
  !     simulation.
  ! ----------------------------------------------------------

  CALL SYSTEM_CLOCK (xseed,yseed,zseed)
  wseed = ABS(xseed - 19027983)

  ! ----------------------------------------------------------
  !                 Initialize variables.
  ! ----------------------------------------------------------

  erro         = .FALSE.
  const        = 1
  numskyelem   = 0
  x(:)         = 0.
  y(:)         = 0.
  z(:)         = 0.
  sumfrac      = 0.
  sidebound    = 0

  ! ----------------------------------------------------------
  !                          - 2.0 -
  !            Reading general parameters from 'amc.inp'
  ! ----------------------------------------------------------

  ! ----------------------------------------------------------
  !         Open 'amc.inp' and read main input parameters.
  ! ----------------------------------------------------------

  OPEN(1,iostat=ios, file='amc.inp', status='old')

  SELECT CASE(ios)

  CASE(0)

     WRITE(*,*) 'File ''amc.inp''... successfully opened'

     READ(1,'(48X,I12)')iter ;  WRITE(*,'(A30,I12)')'Iteration =',iter
     CALL check_range_int(iter,1,999999999,erro)

     READ(1,'(48X,I1)')lsource ;  WRITE(*,'(A30,I12)')'Light source =',lsource
     CALL check_range_int(lsource,0,1,erro)

     READ(1,'(48X,I8)')azimuth ;  WRITE(*,'(A30,I12)')'Azimuth =',azimuth
     CALL check_range_int(azimuth,0,359,erro)

     READ(1,'(48X,I8)')zenith ;  WRITE(*,'(A30,I12)')' Zenith =',zenith
     CALL check_range_int(zenith,0,89,erro)

     READ(1,'(48X,I8)')const ;  WRITE(*,'(A30,I12)')'Constituents =',const
     CALL check_range_int(const,1,20,erro)

     READ(1,'(48X,F12.2)')depthb ;  WRITE(*,'(A30,F12.2)')'Depth =',depthb
     CALL check_range_real(depthb,0.,10000.,erro)

     READ(1,'(48X,F12.2)')sideb ;  WRITE(*,'(A30,F12.2)')'Side =',sideb
     CALL check_range_real(sideb,0.,100000.,erro)

     READ(1,'(48X,I8)')numlogly ;  WRITE(*,'(A30,I12)')'Logging layers =',numlogly
     CALL check_range_int(numlogly,1,500,erro)

     READ(1,'(48X,F12.4)')intlogly ;  WRITE(*,'(A30,F12.2)')'Thickness of logging layers =',intlogly
     CALL check_range_real(intlogly,-1.,500.,erro)

     READ(1,'(48X,I8)')alphaint ;  WRITE(*,'(A30,I12)')'Alpha intervals =',alphaint
     CALL check_range_int(alphaint,1,180,erro)

     READ(1,'(48X,I8)')phiint ;  WRITE(*,'(A30,I12)')'Phi intervals =',phiint
     CALL check_range_int(phiint,1,360,erro)

     READ(1,'(48X,I1)')angint ;  WRITE(*,'(A30,I12)')'Angular or cosine intervals =',angint
     CALL check_range_int(angint,0,1,erro)

     READ(1,'(48X,I1)')absout ;  WRITE(*,'(A30,I12)')'Absorption file =',absout
     CALL check_range_int(absout,0,1,erro)

     READ(1,'(48X,I1)')botout ;  WRITE(*,'(A30,I12)')'Bottom file =',botout
     CALL check_range_int(botout,0,1,erro)

     READ(1,'(48X,I1)')sideout ;  WRITE(*,'(A30,I12)')'Side file =',sideout
     CALL check_range_int(sideout,0,1,erro)

     READ(1,'(48X,I1)')spftyp ;  WRITE(*,'(A30,I12)')'Source of SPF =',spftyp
     CALL check_range_int(spftyp,0,1,erro)

     READ(1,'(48X,I8)')spfang ;  WRITE(*,'(A30,I12)')'Number of SPF angles =',spfang
     CALL check_range_int(spfang,0,360,erro)

     READ(1,'(48X,I8)')numlycst ;  WRITE(*,'(A30,I12)')'Depth dependent layer =',numlycst
     CALL check_range_int(numlycst,1,100,erro)

     READ(1,'(48X,I8)')numgrid ;  WRITE(*,'(A30,I12)')'Number of cells =',numgrid
     CALL check_range_int(numgrid,0,101,erro)
     CALL check_odd(numgrid,erro)

     READ(1,'(48X,F12.6)')cellsize ;  WRITE(*,'(A30,F12.2)')'Cell size =',cellsize
     CALL check_range_real(cellsize,0.,100000.,erro)

     READ(1,'(48X,I1)')surfer ;  WRITE(*,'(A30,I12)')'Surfer output =',surfer
     CALL check_range_int(surfer,0,1,erro)

     READ(1,'(48X,I1)')lgttyp ;  WRITE(*,'(A30,I12)')'Light footprint =',lgttyp
     CALL check_range_int(lgttyp,1,3,erro)

     READ(1,'(48X,F12.6)')lgtdiam ;  WRITE(*,'(A30,F12.6)')'Circle diameter =',lgtdiam
     CALL check_range_real(lgtdiam,0.,10000000.,erro)

     READ(1,'(48X,F12.6)')lgtx ;  WRITE(*,'(A30,F12.6)')'X size of light footprint =',lgtx
     CALL check_range_real(lgtx,0.,10000000.,erro)

     READ(1,'(48X,F12.6)')lgty ;  WRITE(*,'(A30,F12.6)')'Y size of light footprint =',lgty
     CALL check_range_real(lgty,0.,10000000.,erro)

     READ(1,'(48X,I1)')targbot ;  WRITE(*,'(A30,I12)')'Bottom target =',targbot
     CALL check_range_int(targbot,0,1,erro)

     READ(1,'(48X,F12.6)')targx ;  WRITE(*,'(A30,F12.6)')'X size of target =',targx
     CALL check_range_real(targx,0.,10000000.,erro)

     READ(1,'(48X,F12.6)')targy ;  WRITE(*,'(A30,F12.6)')'Y size of target =',targy
     CALL check_range_real(targy,0.,10000000.,erro)

     READ(1,'(48X,I1)')normh2o ;  WRITE(*,'(A30,I12)')'Normalize to =',normh2o
     CALL check_range_int(normh2o,0,1,erro)

     READ(1,'(48X,I8)')lambda ;  WRITE(*,'(A30,I12)')'Number of wavebands =',lambda
     CALL check_range_int(lambda,1,500,erro)

     WRITE (*,*)'----------------------------------------------------------------'

  CASE DEFAULT

     WRITE(*,*) 'Failed to open ''amc.inp'' '
     WRITE(*,*) 'IOSTAT error code = ', ios
     erro = .TRUE.

  END SELECT

  CLOSE(1)

  ! ----------------------------------------------------------
  !                          - 2.1 -
  !                    Post initialization
  ! ----------------------------------------------------------

  nalpha       = PI / REAL( alphaint )
  nphi         = 2 * PI /phiint
  muu          = phiint / ( (REAL(alphaint) /2. -1 )* phiint +1 )
  mum          = muu / phiint

  ! Use polar cap if the cosine of the polar angle is used as equal intervals


  ! ----------------------------------------------------------
  !                          - 3.0 -
  !                    Allocation of memory
  ! ----------------------------------------------------------

  ! Allocated memory for N

  IF (numgrid == 0) THEN

     ALLOCATE (n(lambda,0:0,0:0,-1:numlogly,alphaint,phiint) , STAT=astat)

  ELSE

     sidebound = NINT( (numgrid -1 ) / 2.)
     ALLOCATE (n( lambda,-sidebound : sidebound,-sidebound : sidebound,  &
          -1:numlogly,alphaint,phiint) , STAT=astat)
  END IF

  IF(astat/=0) THEN

     WRITE(*,*) 'Error: array ''N'' not allocated! Do you have enough memory?'
     erro = .TRUE.

  ELSE

     n(:,:,:,:,:,:) = 0

  END IF

  ! Allocated memory for layval

  ALLOCATE (layval(-1:numlogly), STAT=astat)

  IF(astat/=0) THEN

     WRITE(*,*) 'Error: array ''layval'' not allocated! Do you have enough memory?'
     erro = .TRUE.

  ELSE

     layval(-1:0) = 0.0

     IF (intlogly > -1.0) THEN

        DO ii = 1 , (numlogly - 1)

           layval(ii) = layval(ii-1)+intlogly

        END DO

        layval(numlogly) = depthb

     ELSE

        DO ii = 1 , numlogly

           layval(ii) = (depthb / numlogly) * ii

        END DO

     END IF

  END IF

  ! Allocated memory for conc

  ALLOCATE (conc(numlycst,const), STAT=astat)

  IF(astat/=0) THEN

     WRITE(*,*) 'Error: array ''conc'' not allocated! Do you have enough memory?'
     erro = .TRUE.

  ELSE

     conc(:,:) = 0.0

  END IF

  ! Allocated memory for specabs

  ALLOCATE (specabs(lambda,const), STAT=astat)

  IF(astat/=0) THEN

     WRITE(*,*) 'Error: array ''specabs'' not allocated! Do you have enough memory?'
     erro = .TRUE.

  ELSE

     specabs (:,:) = 0.0

  END IF

  ! Allocated memory for absorb

  ALLOCATE (absorb(lambda,numlycst), STAT=astat)

  IF(astat/=0) THEN

     WRITE(*,*) 'Error: array ''absorb'' not allocated! Do you have enough memory?'
     erro = .TRUE.

  ELSE

     absorb (:,:) = 0.0

  END IF

  ! Allocated memory for specscat

  ALLOCATE (specscat(lambda,const), STAT=astat)

  IF(astat/=0) THEN

     WRITE(*,*) 'Error: array ''specscat'' not allocated! Do you have enough memory?'
     erro = .TRUE.

  ELSE

     specscat(:,:) = 0.0

  END IF

  ! Allocated memory for scatter

  ALLOCATE (scatter(lambda,numlycst), STAT=astat)

  IF(astat/=0) THEN

     WRITE(*,*) 'Error: array ''scatter'' not allocated! Do you have enough memory?'
     erro = .TRUE.

  ELSE

     scatter (:,:) = 0.0

  END IF

  ! Allocated memory for fracscat

  ALLOCATE (fracscat(lambda,numlycst,const), STAT=astat)

  IF(astat/=0) THEN

     WRITE(*,*) 'Error: array ''fracscat'' not allocated! Do you have enough memory?'
     erro = .TRUE.

  ELSE

     fracscat (:,:,:) = 0.0

  END IF

  ! Allocated memory for atten

  ALLOCATE (atten(lambda,numlycst), STAT=astat)

  IF(astat/=0) THEN

     WRITE(*,*) 'Error: array ''atten'' not allocated! Do you have enough memory?'
     erro = .TRUE.

  ELSE

     atten (:,:) = 0.0

  END IF

  ! Allocated memory for salbedo

  ALLOCATE (salbedo(lambda,numlycst), STAT=astat)

  IF(astat/=0) THEN

     WRITE(*,*) 'Error: array ''salbedo'' not allocated! Do you have enough memory?'
     erro = .TRUE.

  ELSE

     salbedo (:,:) = 0.0

  END IF

  ! Allocated memory for spf

  ALLOCATE (spf(spfang,0:const), STAT=astat)

  IF(astat/=0) THEN

     WRITE(*,*) 'Error: array ''spf'' not allocated! Do you have enough memory?'
     erro = .TRUE.

  ELSE

     spf (:,:)= 0.0

  END IF

  ! Allocated memory for slopespf

  ALLOCATE (slopespf(2:spfang,1:const), STAT=astat)

  IF(astat/=0) THEN

     WRITE(*,*) 'Error: array ''slopespf'' not allocated! Do you have enough memory?'
     erro = .TRUE.

  ELSE

     slopespf (:,:)= 0.0

  END IF

  ! Allocated memory for refr

  ALLOCATE (refr(0:numlycst), STAT=astat)

  IF(astat/=0) THEN

     WRITE(*,*) 'Error: array ''refr'' not allocated! Do you have enough memory?'
     erro = .TRUE.

  ELSE

     refr (:)= 1.0

  END IF

  ! Allocated memory for dptlycst

  ALLOCATE (dptlycst(numlycst), STAT=astat)

  IF(astat/=0) THEN

     WRITE(*,*) 'Error: array ''dptlycst'' not allocated! Do you have enough memory?'
     erro = .TRUE.

  ELSE

     dptlycst (:)= 1.0

  END IF

  ! Allocated memory for wavelength

  ALLOCATE (wavelength(lambda), STAT=astat)

  IF(astat/=0) THEN

     WRITE(*,*) 'Error: array ''wavelength'' not allocated! Do you have enough memory?'
     erro = .TRUE.

  ELSE

     wavelength (:)= 0.0

  END IF

  ! Allocated memory for bottomr

  ALLOCATE (bottomr(lambda), STAT=astat)

  IF(astat/=0) THEN

     WRITE(*,*) 'Error: array ''bottomr'' not allocated! Do you have enough memory?'
     erro = .TRUE.

  ELSE

     bottomr (:)= 0.0

  END IF

  ! Allocated memory for bottomspc

  ALLOCATE (bottomspc(lambda), STAT=astat)

  IF(astat/=0) THEN

     WRITE(*,*) 'Error: array ''bottomspc'' not allocated! Do you have enough memory?'
     erro = .TRUE.

  ELSE

     bottomspc (:)= 0.0

  END IF

  SELECT CASE (targbot)

  CASE (1)

     ! Allocated memory for targref

     ALLOCATE (targref(lambda), STAT=astat)

     IF(astat/=0) THEN

        WRITE(*,*) 'Error: array ''targref'' not allocated! Do you have enough memory?'
        erro = .TRUE.

     ELSE

        targref (:)= 0.0

     END IF

     ! Allocated memory for targspec

     ALLOCATE (targspc(lambda), STAT=astat)

     IF(astat/=0) THEN

        WRITE(*,*) 'Error: array ''targspc'' not allocated! Do you have enough memory?'
        erro = .TRUE.

     ELSE

        targspc (:)= 0.0

     END IF

  CASE DEFAULT

  END SELECT

  IF (lsource == 0 )THEN

     ! Allocated memory for direct

     ALLOCATE (direct(lambda), STAT=astat)

     IF(astat/=0) THEN

        WRITE(*,*) 'Error: array ''direct'' not allocated! Do you have enough memory?'
        erro = .TRUE.

     ELSE

        direct (:)= 0.0

     END IF

  END IF

  ! Allocated memory for intensity

  ALLOCATE (intensity(lambda), STAT=astat)

  IF(astat/=0) THEN

     WRITE(*,*) 'Error: array ''intensity'' not allocated! Do you have enough memory?'
     erro = .TRUE.

  ELSE

     intensity (:)= 0.0

  END IF

  ! Allocated memory for fdown

  ALLOCATE (fdown(lambda,-1:numlogly), STAT=astat)

  IF(astat/=0) THEN

     WRITE(*,*) 'Error: array ''fdown'' not allocated! Do you have enough memory?'
     erro = .TRUE.

  ELSE

     fdown (:,:)= 0.0

  END IF

  ! Allocated memory for fup

  ALLOCATE (fup(lambda,-1:numlogly), STAT=astat)

  IF(astat/=0) THEN

     WRITE(*,*) 'Error: array ''fup'' not allocated! Do you have enough memory?'
     erro = .TRUE.

  ELSE

     fup (:,:)= 0.0

  END IF

  ! Allocated memory for bdown

  ALLOCATE (bdown(lambda,-1:numlogly), STAT=astat)

  IF(astat/=0) THEN

     WRITE(*,*) 'Error: array ''bdown'' not allocated! Do you have enough memory?'
     erro = .TRUE.

  ELSE

     bdown (:,:)= 0.0

  END IF

  ! Allocated memory for bup

  ALLOCATE (bup(lambda,-1:numlogly), STAT=astat)

  IF(astat/=0) THEN

     WRITE(*,*) 'Error: array ''bup'' not allocated! Do you have enough memory?'
     erro = .TRUE.

  ELSE

     bup (:,:)= 0.0

  END IF

  ! Allocated memory for totalback

  ALLOCATE (totalback(lambda,-1:numlogly), STAT=astat)

  IF(astat/=0) THEN

     WRITE(*,*) 'Error: array ''totalback'' not allocated! Do you have enough memory?'
     erro = .TRUE.

  ELSE

     totalback (:,:)= 0.0

  END IF

  ! Allocated memory for totalfwd

  ALLOCATE (totalfwd(lambda,-1:numlogly), STAT=astat)

  IF(astat/=0) THEN

     WRITE(*,*) 'Error: array ''totalfwd'' not allocated! Do you have enough memory?'
     erro = .TRUE.

  ELSE

     totalfwd (:,:)= 0.0

  END IF

  ! Allocated memory for cositer

  ALLOCATE (cositer(lambda), STAT=astat)

  IF(astat/=0) THEN

     WRITE(*,*) 'Error: array ''cositer'' not allocated! Do you have enough memory?'
     erro = .TRUE.

  ELSE

     cositer (:)= 0.0

  END IF


  ! ----------------------------------------------------------
  !                          - 4.0 -
  !                  Reading of input files
  ! ----------------------------------------------------------


  ! ----------------------------------------------------------
  !                          - 4.1 -                        
  ! Open 'conc.inp' and read depth dependent concentration data
  ! ----------------------------------------------------------

  OPEN(1,iostat=ios, file='conc.inp', status='old')

  SELECT CASE(ios)

  CASE(0)

     WRITE(*,*) 'File ''conc.inp''... successfully opened'

     ! First value in file reprents the number of lines (records) to read

     READ(1,*)records

     IF (records == numlycst ) THEN

        DO ii=1,numlycst

           READ(1,*)(conc(ii,jj),jj=1,const),refr(ii),dptlycst(ii)

           ! Check that dptlycst does not exceed maximum depth

           IF (dptlycst(ii) > depthb) THEN

              WRITE(*,*) 'ERROR!! boundary of layer',ii,'is greater than depth.'
              erro = .TRUE.

           END IF

        END DO

        ! Make sure that the top boundary of the top layer is 0

        dptlycst(1)=0.0

     ELSE

        WRITE(*,*)'number of records in ''conc.inp'' do not match number of layers'
        erro = .TRUE.

     END IF

  CASE DEFAULT

     WRITE(*,*) 'Failed to open ''conc.inp'' '
     WRITE(*,*) 'IOSTAT error code = ', ios
     erro = .TRUE.

  END SELECT

  CLOSE(1)

  ! ----------------------------------------------------------
  !                          - 4.2 -
  !  Open 'abs.inp' and read specific absorption coefficients
  ! ----------------------------------------------------------

  OPEN(2,iostat=ios, file='abs.inp', status='old')

  SELECT CASE(ios)

  CASE(0)

     WRITE(*,*) 'File ''abs.inp''... successfully opened'

     ! First value in file reprents the number of lines (records) to read

     READ(2,*)records

     IF (records == lambda ) THEN

        ! Read specific absorption coefficient data

        DO ii=1,lambda

           READ(2,*)(specabs(ii,jj),jj=1,const)

        END DO

        ! Calculate the bulk absorption coefficient

        absorb(:,:)  = 0.0

        DO jj=1,numlycst

           DO ii=1,lambda

              DO kk=1,const

                 absorb(ii,jj)=absorb(ii,jj)+specabs(ii,kk)*conc(jj,kk) 

                 ! The preceding line can be replaced with the following
                 ! for simulating Mobley est al.'s (1993) base problem 3 
                 !
                 ! absorb(ii,jj)=absorb(ii,jj)+specabs(ii,kk)*conc(jj,kk)**0.602 
                 !

              END DO

           END DO

        END DO

     ELSE

        WRITE(*,*)'number of records in ''abs.inp'' do not match lambda'
        erro = .TRUE.

     END IF

  CASE DEFAULT

     WRITE(*,*) 'Failed to open ''abs.inp'' '
     WRITE(*,*) 'IOSTAT error code = ', ios
     erro = .TRUE.

  END SELECT

  CLOSE(2)

  ! ----------------------------------------------------------
  !                          - 4.3 -
  ! Open 'scat.inp' and read specific scattering coefficients
  ! ----------------------------------------------------------

  OPEN(3,iostat=ios, file='scat.inp', status='old')

  SELECT CASE(ios)

  CASE(0)

     WRITE(*,*) 'File ''scat.inp''... successfully opened'

     ! First value in file represents the number of lines (records) to read

     READ(3,*)records

     IF (records == lambda) THEN

        ! Read specific scattering coefficient data

        DO ii=1,records

           READ(3,*)(specscat(ii,jj),jj=1,const)

        END DO

        ! Calculate the bulk scattering coefficient

        DO jj=1,numlycst

           DO ii=1,lambda

              DO kk=1,const

                 scatter(ii,jj)=scatter(ii,jj)+specscat(ii,kk)*conc(jj,kk)

                 ! The preceding line can be replaced with the following
                 ! for simulating Mobley est al.'s (1993) base problem 3 
                 !
                 ! scatter(ii,jj)=scatter(ii,jj)+specscat(ii,kk) * conc(jj,kk)**0.620 
                 ! 

              END DO

           END DO

        END DO

        ! Calculate the probability of scattering for
        ! each constituent at each wavelength and optical
        ! depth

        DO jj=1,numlycst

           DO ii=1,lambda

              sumfrac = 0.0

              DO kk=1,const

                 ! Replace the following first uncommented line with the following
                 ! commented line to simulated Mobley et al.'s (1993) problem 3
                 ! fracscat(ii,jj,kk) = specscat(ii,kk)*conc(jj,kk)**0.620 / &
                 fracscat(ii,jj,kk) = specscat(ii,kk)*conc(jj,kk)/ &  
                      scatter(ii,jj) + sumfrac
                 sumfrac = fracscat(ii,jj,kk)  !!ORIGINAL

              END DO

           END DO

        END DO

     ELSE

        WRITE(*,*)'number of records in ''scat.inp'' do not match lambda'
        erro = .TRUE.

     END IF

  CASE DEFAULT

     WRITE(*,*) 'Failed to open ''scat.inp'' '
     WRITE(*,*) 'IOSTAT error code = ', ios
     erro = .TRUE.

  END SELECT

  CLOSE(3)

  ! ----------------------------------------------------------
  !                          - 4.4 -
  !                 Deallocation of some memory
  ! ----------------------------------------------------------

  DEALLOCATE (specabs,STAT=astat)

  IF (astat /= 0) THEN

     WRITE(*,*)'ERROR: problem in array ''specabs'' deallocation'
     erro = .TRUE.

  END IF

  DEALLOCATE (specscat,STAT=astat)

  IF (astat /= 0) THEN

     WRITE(*,*)'Warning: problem in array ''specscat'' deallocation'
     erro = .TRUE.

  END IF

  DEALLOCATE (conc,STAT=astat)

  IF (astat /= 0) THEN

     WRITE(*,*)'ERROR: problem in array ''conc'' deallocation'
     erro = .TRUE.

  END IF

  ! ----------------------------------------------------------
  !                          - 4.5 -
  !  Open 'spf.inp' and read scattering phase function data
  ! ----------------------------------------------------------

  OPEN(4,iostat=ios, file='spf.inp', status='old')

  SELECT CASE(ios)

  CASE(0)

     WRITE(*,*) 'File ''spf.inp''... successfully opened'
     READ(4,*)records

     IF (records == spfang) THEN

        DO ii=1,spfang

           READ(4,*)(spf(ii,jj),jj=0,const)

        END DO

        ! Calculate the slope between each value

        DO ii=2,spfang

           slopespf(ii,1:const) = ( spf(ii,1:const) - spf(ii-1,1:const) ) /  &
                ( spf(ii,0) - spf(ii-1,0) )
        END DO

     ELSE

        WRITE(*,*)'number of records in ''spf.inp'' do not match number of angles.'
        erro = .TRUE.

     END IF

  CASE DEFAULT

     WRITE(*,*) 'Failed to open ''spf.inp'' '
     WRITE(*,*) 'IOSTAT error code = ', ios
     erro = .TRUE.

  END SELECT

  CLOSE(4)

  ! ----------------------------------------------------------
  !                          - 4.6 -
  ! Open 'lambbot.inp' and read wavelength and bottom data
  ! ----------------------------------------------------------

  OPEN(8,iostat=ios, file='lambbot.inp', status='old')

  SELECT CASE(ios)

  CASE(0)

     WRITE(*,*) 'File ''lambbot.inp''... successfully opened'
     READ(8,*)records

     IF (records == lambda) THEN

        SELECT CASE (targbot)

        CASE (0)

           DO ii=1,lambda

              READ(8,*)wavelength(ii),bottomr(ii),bottomspc(ii)

           END DO

        CASE DEFAULT

           DO ii=1,lambda

              READ(8,*)wavelength(ii),bottomr(ii),bottomspc(ii),&
                   targref(ii),targspc(ii)

           END DO

        END SELECT

     ELSE

        WRITE(*,*)'number of records in ''lambbot.inp'' do not match number of wavelengths.'
        erro = .TRUE.

     END IF

  CASE DEFAULT

     WRITE(*,*) 'Failed to open ''lambbot.inp'' '
     WRITE(*,*) 'IOSTAT error code = ', ios
     erro = .TRUE.

  END SELECT

  CLOSE(8)

  ! ----------------------------------------------------------
  !                          - 4.7 -
  ! Open 'difcol.inp' and read the diffuse/colimated
  ! probability. Only read intensity() data  if tabulated data
  ! is not used.
  ! If the user chooses to input external skylight distribution
  ! data, this step is skipped
  ! ----------------------------------------------------------



  SELECT CASE (lsource)

  CASE (0) ! Photon source is internally computed

     OPEN(9,iostat=ios, file='difcol.inp', status='old')

     SELECT CASE(ios)

     CASE(0)

        WRITE(*,*) 'File ''difcol.inp''... successfully opened'
        READ(9,*)records

        IF (records == lambda) THEN

           DO ii=1,lambda

              READ(9,*)direct(ii),intensity(ii)

           END DO

        ELSE

           WRITE(*,*)'Number of records in ''difcol.inp'' do not match number of wavelengths.'
           erro = .TRUE.

        END IF

     CASE DEFAULT

        WRITE(*,*) 'Failed to open ''difcol.inp'' '
        WRITE(*,*) 'IOSTAT error code = ', ios
        erro = .TRUE.

     END SELECT

     CLOSE(9)

  CASE DEFAULT ! Read data from tabulated skylight distribution file

     OPEN(10,iostat=ios, file='skydist.inp', status='old')

     SELECT CASE(ios)

     CASE(0)

        WRITE(*,*) 'File ''skydist.inp''... successfully opened'
        READ(10,*) numskyelem,sky_zenith,sky_azimuth,numskywave
        ALLOCATE (sky_wave(numskywave), STAT=astat)

        IF(astat/=0) THEN

           WRITE(*,*) 'Error: array ''sky_wave'' not allocated!'
           erro = .TRUE.

        ELSE

           READ(10,*)(sky_wave(ii),ii=1,numskywave)

        END IF

        ALLOCATE (sky_int(numskywave,sky_zenith,0:sky_azimuth), STAT=astat)

        IF(astat/=0) THEN

           WRITE(*,*) 'Error: array ''directions'' not allocated!'
           erro = .TRUE.

        ELSE

           sky_int(:,:,:) = 0

           DO ii=1,numskyelem

              READ(10,*)jj,kk,ll,sky_int(jj,kk,ll)

           END DO

        END IF

        ! Now convert radiance values to actual number of photons
        ! for each sky element

        ALLOCATE (sky_nint(numskywave,sky_zenith,0:sky_azimuth), STAT=astat)

        IF(astat/=0) THEN

           WRITE(*,*) 'Error: array ''sky_nint'' not allocated!'
           erro = .TRUE.

        ELSE

           sky_nint(:,:,:) = 0

           DO ii = 1 , numskywave

              DO jj = 1 , sky_zenith

                 ai = jj * PI / 180  

                 IF (jj == 1 .OR. jj == 89) THEN

                    sky_int(ii,jj,:)  =  sky_int( ii,jj,: ) * SIN ( ai ) * &  
                         ( 0.0261799) *                    &  
                         ( 0.01745329)*                    &
                         COS( ai )

                 ELSE 

                    sky_int(ii,jj,:)  =  sky_int( ii,jj,: ) * SIN ( ai ) * &   
                         ( 0.01745329) *                   &   
                         ( 0.01745329) *                   &
                         COS( ai )


                 END IF

              END DO

              ! Normalize directional component of incident photons

              sky_int(ii,:,:) = NINT(sky_int(ii,:,:) / SUM (sky_int(ii,:,:))* &
                   iter)
              total = 0

              ! Compute cumulative photons for each wavelength

              DO jj = 1 , sky_zenith

                 DO kk = 0 , sky_azimuth

                    sky_nint(ii,jj,kk) = NINT(sky_int(ii,jj,kk)+ total)
                    total =  sky_nint(ii,jj,kk)

                 END DO

              END DO

           END DO

        END IF

        DEALLOCATE( sky_int,STAT = astat )

        IF (astat /= 0) THEN

           WRITE(*,*)'ERROR: problem in array ''direction'' deallocation'
           erro = .TRUE.

        END IF

     CASE DEFAULT

        WRITE(*,*) 'Failed to open ''skydist.inp'' '
        WRITE(*,*) 'IOSTAT error code = ', ios
        erro = .TRUE.

     END SELECT

     CLOSE(10)

  END SELECT


  ! ----------------------------------------------------------
  !                          - 5.0 -
  !                  Calculation of IOP
  ! ----------------------------------------------------------

  ! ----------------------------------------------------------
  !                          - 5.1 -  
  !        Calculate the bulk beam attenuation coefficient
  ! ----------------------------------------------------------

  atten(:,:) = absorb(:,:) + scatter(:,:)

  ! ----------------------------------------------------------
  !                          - 5.2 -
  !        Calculate the scattering albeldo
  ! ----------------------------------------------------------

  salbedo(:,:) = scatter(:,:) /  atten(:,:)

  ! ----------------------------------------------------------
  !                          - 5.3 -
  !      Write calculated IOP's to file 'iop.out'
  ! ----------------------------------------------------------

  OPEN(6,iostat=ios, file='iop.out', status='unknown')

  SELECT CASE(ios)

  CASE(0)

     DO jj=1,numlycst

        WRITE(6,'(a10,2x,i11,2x,a11,2x,f11.8,1x,2x,A10,a1)') &
             'Layer = ', jj ,'Starting at:', dptlycst(jj), 'm'

        WRITE(6,'(a10,4(2x,a11))')'Wavelength','absorption',&
             'scattering','attenuation', 'scat albedo'

        DO ii=1,lambda

           WRITE(6,'(f10.3,4(2x,f11.7))') wavelength(ii), &
                absorb(ii,jj),scatter(ii,jj),atten(ii,jj), &
                salbedo(ii,jj)

        END DO

     END DO

  CASE default

     WRITE(*,*) 'Failed to open ''attenuation.out'' '
     WRITE(*,*) '  Could not calculate the bulk attenuation coefficient'
     WRITE(*,*) '  IOSTAT error code = ', ios
     erro = .TRUE.

  END SELECT

  CLOSE(6)

  WRITE (*,*)'----------------------------------------------------------------'

  ! ----------------------------------------------------------
  !                          - 6.0 -
  !                         Run model
  ! ----------------------------------------------------------

  WRITE(*,*)'Running model .....'

  IF ( absout == 1 ) THEN

     OPEN(10,iostat=ios, file='absorption.out', status='unknown')

     SELECT CASE(ios)

     CASE(0)

        WRITE(*,*) 'Succefully opened absorption.out file for output'

     CASE default

        WRITE(*,*) 'Failed to open ''absorption.out'' '
        erro = .TRUE.

     END SELECT

  END IF

  SELECT CASE (erro)

  CASE (.FALSE.)    ! No errors occurred, simulation proceeded

     LAMBDA_LOOP: DO lam= 1 , lambda

        WRITE(*,*)'Running simulation for wavelength = ',wavelength(lam)

        SELECT CASE (lsource)

        CASE (0)

           titer = NINT(intensity(lam) * iter)

        CASE DEFAULT  ! Determine which wavelength to use from sky dist. data and
                      ! compute the number of iterations to run

           wavecontrol = 1

           DO jj = 1 , numskywave

              IF ( jj == numskywave) THEN

                 titer = MAXVAL(sky_nint(jj,:,:))
                 wavematch = jj
                 EXIT

              ELSEIF ( wavelength(lam)< ( sky_wave(jj) + sky_wave(jj+1))/2. ) THEN

                 titer = MAXVAL(sky_nint(jj,:,:))
                 wavematch = jj
                 EXIT

              END IF

           END DO

           WRITE(*,*) 'titer = ', titer,'at wavelength', sky_wave(wavematch)

        END SELECT

        totalop = 0.0

        INTENSITY_LOOP: DO jj = 1 , titer

           ! Determine direction of travel of incident light

           CALL light(wavematch,lam,jj,wavecontrol)

           ! Determine the (x,y,z) location of the photon's interaction
           ! with the air water interface

           IF (lgttyp == 1) THEN  ! Point location

              x(1) = 0.0
              y(1) = 0.0

           ELSEIF (lgttyp == 2) THEN   ! Circle foot print

              randradius = SQRT(rand()) * lgtdiam / 2.
              randangle  = rand() * 2. * PI
              x(1) = COS(randangle) * randradius
              y(1) = SIN(randangle) * randradius

           ELSEIF (lgttyp == 3) THEN   ! Square foot print

              x(1) = ( rand() * lgtx ) - lgtx / 2.
              y(1) = ( rand() * lgty ) - lgty / 2.

           END IF

           z(1) = 0.0

           ! Log all photons striking the surface  

           CALL logbin(x(1),y(1),-1,lam)

           ! Log information about incident photon (used for normalizing to above water irradiance)

           cositer(lam) = cositer(lam) + 1 / COS( theta )

           ! Determine outcome of air water interface

           CALL interface_sub(refr(0),refr(1),switch)

           ! Log photon activity at the interface and proceed with
           ! simulation

           IF (theta < ( PI / 2.) ) THEN  ! Photon enters water

              CALL logbin(x(1),y(1),0,lam)
              CALL water(lam,action)

              IF (action == 1 .AND. absout == 1) THEN

                 WRITE(10,'(3(F12.4,2x))')x(1),y(1),z(1)

              END IF

           ELSE

              CALL logbin(x(1),y(1),-1,lam) ! Photon reflects off of water

           END IF

        END DO INTENSITY_LOOP

        WRITE(*,'(A28,2x,F6.2)')' % scattered forward  = ', 100 * REAL( SUM( totalfwd(lam,:) ) /    &
             REAL( SUM( totalfwd(lam,:) ) + SUM( totalback(lam,:) )))
        WRITE(*,'(A28,2x,F6.2)')' % scattered backward = ', 100 * REAL( SUM( totalback(lam,:) ) /   &
             REAL( SUM( totalfwd(lam,:) ) + SUM( totalback(lam,:) )))
        WRITE(*,'(A28,2x,F6.2,A22)')' mean optical pathlength = ', REAL( totalop / REAL( titer )),' (still experimental)'

     END DO LAMBDA_LOOP

     IF (absout == 1) THEN

        CLOSE(10)

     END IF

  CASE DEFAULT  ! Errors occurred, simulation did not proceed

     WRITE(*,*) 'Errors occurred... Could not continue with simulation'

  END SELECT

  WRITE(*,*) 'Done with simulation. Now computing/writing output...'

  ! ----------------------------------------------------------
  !                          - 6.1 -
  !                 Free up some memory
  ! ----------------------------------------------------------

  DEALLOCATE(absorb , STAT=astat)

  IF (astat /= 0) THEN

     WRITE(*,*)'ERROR: problem in array ''absorb'' deallocation'
     erro = .TRUE.

  END IF

  DEALLOCATE(scatter , STAT=astat)

  IF (astat /= 0) THEN

     WRITE(*,*)'ERROR: problem in array ''scatter'' deallocation'
     erro = .TRUE.

  END IF

  DEALLOCATE(fracscat , STAT=astat)

  IF (astat /= 0) THEN

     WRITE(*,*)'ERROR: problem in array ''fracscat'' deallocation'
     erro = .TRUE.

  END IF

  DEALLOCATE(atten , STAT=astat)

  IF (astat /= 0) THEN

     WRITE(*,*)'ERROR: problem in array ''atten'' deallocation'
     erro = .TRUE.

  END IF

  DEALLOCATE(salbedo , STAT=astat)

  IF (astat /= 0) THEN

     WRITE(*,*)'ERROR: problem in array ''salbedo'' deallocation'
     erro = .TRUE.

  END IF

  DEALLOCATE(dptlycst , STAT=astat)

  IF (astat /= 0) THEN

     WRITE(*,*)'ERROR: problem in array ''dptlycst'' deallocation'
     erro = .TRUE.

  END IF

  DEALLOCATE(spf , STAT=astat)

  IF (astat /= 0) THEN

     WRITE(*,*)'ERROR: problem in array ''spf'' deallocation'
     erro = .TRUE.

  END IF

  DEALLOCATE(slopespf , STAT=astat)

  IF (astat /= 0) THEN

     WRITE(*,*)'ERROR: problem in array ''slopespf'' deallocation'
     erro = .TRUE.

  END IF

  DEALLOCATE (bottomr, STAT=astat)

  IF(astat/=0) THEN

     WRITE(*,*) 'Error: array ''bottomr'' not deallocated!'
     erro = .TRUE.

  END IF

  DEALLOCATE (bottomspc, STAT=astat)

  IF(astat/=0) THEN
     WRITE(*,*) 'Error: array ''bottomspc'' not deallocated!'
     erro = .TRUE.

  END IF

  SELECT CASE (targbot)

  CASE (1)

     DEALLOCATE (targref, STAT=astat)

     IF(astat/=0) THEN

        WRITE(*,*) 'Error: array ''targref'' not deallocated!'
        erro = .TRUE.

     END IF

     DEALLOCATE (targspc, STAT=astat)

     IF(astat/=0) THEN

        WRITE(*,*) 'Error: array ''targspc'' not deallocated!'
        erro = .TRUE.

     END IF

  CASE DEFAULT

  END SELECT

  IF (lsource == 0) THEN

     DEALLOCATE (direct, STAT=astat)

     IF(astat/=0) THEN

        WRITE(*,*) 'Error: array ''direct'' not deallocated!'
        erro = .TRUE.

     END IF

     DEALLOCATE (intensity, STAT=astat)

     IF(astat/=0) THEN

        WRITE(*,*) 'Error: array ''intensity'' not deallocated!'
        erro = .TRUE.

     END IF

  ELSE

     DEALLOCATE(sky_nint,STAT=astat)

     IF (astat /= 0) THEN

        WRITE(*,*)'ERROR: problem in array ''ndirection'' deallocation'
        erro = .TRUE.

     END IF

     DEALLOCATE(sky_wave,STAT=astat)

     IF (astat /= 0) THEN

        WRITE(*,*)'ERROR: problem in array ''sky_wave'' deallocation'
        erro = .TRUE.

     END IF

  END IF

  ! ----------------------------------------------------------
  !                          - 7.0 -
  !                        Compute AOP
  ! ----------------------------------------------------------

  SELECT CASE (erro)

  CASE (.FALSE.)    ! No errors occurred, compute AOP

     ! ----------------------------------------------------------
     !                          - 7.1 -
     !                       Allocate arrays
     ! ----------------------------------------------------------

     ! - L -

     IF (numgrid == 0) THEN

        ALLOCATE (rad(lambda,0:0,0:0,-1:numlogly,alphaint,phiint) , STAT=astat)

     ELSE

        sidebound = NINT( (numgrid -1 ) / 2.)
        ALLOCATE (rad( lambda,-sidebound : sidebound,-sidebound : sidebound,  &
             -1:numlogly,alphaint,phiint) , STAT=astat)

     END IF

     IF(astat/=0) THEN

        WRITE(*,*) 'Error: array ''rad'' not allocated!'
        erro = .TRUE.

     ELSE

        rad(:,:,:,:,:,:) = 0

     END IF

     ! - Eu -

     IF (numgrid == 0) THEN

        ALLOCATE (eu(lambda,0:0,0:0,-1:numlogly) , STAT=astat)
     ELSE

        sidebound = NINT( (numgrid -1 ) / 2.)
        ALLOCATE (eu( lambda,-sidebound : sidebound,-sidebound : sidebound,  &
             -1:numlogly) , STAT=astat)

     END IF

     IF(astat/=0) THEN

        WRITE(*,*) 'Error: array ''eu'' not allocated!'
        erro = .TRUE.

     ELSE

        eu(:,:,:,:) = 0

     END IF

     ! - Ed -

     IF (numgrid == 0) THEN

        ALLOCATE (ed(lambda,0:0,0:0,-1:numlogly) , STAT=astat)

     ELSE

        sidebound = NINT( (numgrid -1 ) / 2.)
        ALLOCATE (ed( lambda,-sidebound : sidebound,-sidebound : sidebound,  &
             -1:numlogly) , STAT=astat)

     END IF

     IF(astat/=0) THEN

        WRITE(*,*) 'Error: array ''ed'' not allocated!'
        erro = .TRUE.

     ELSE

        ed(:,:,:,:) = 0

     END IF

     ! - Eou-

     IF (numgrid == 0) THEN

        ALLOCATE (eou(lambda,0:0,0:0,-1:numlogly) , STAT=astat)

     ELSE

        sidebound = NINT( (numgrid -1 ) / 2.)
        ALLOCATE (eou( lambda,-sidebound : sidebound,-sidebound : sidebound,  &
             -1:numlogly) , STAT=astat)

     END IF

     IF(astat/=0) THEN

        WRITE(*,*) 'Error: array ''eou'' not allocated!'
        erro = .TRUE.

     ELSE

        eou(:,:,:,:) = 0

     END IF

     ! - Eod -

     IF (numgrid == 0) THEN

        ALLOCATE (eod(lambda,0:0,0:0,-1:numlogly) , STAT=astat)

     ELSE

        sidebound = NINT( (numgrid -1 ) / 2.)
        ALLOCATE (eod( lambda,-sidebound : sidebound,-sidebound : sidebound,  &
             -1:numlogly) , STAT=astat)

     END IF

     IF(astat/=0) THEN

        WRITE(*,*) 'Error: array ''eod'' not allocated!'
        erro = .TRUE.

     ELSE

        eod(:,:,:,:) = 0

     END IF

     ! - polar_ld -

     IF (numgrid == 0) THEN

        ALLOCATE (polar_ld(lambda,0:0,0:0,-1:numlogly) , STAT=astat)
     ELSE

        sidebound = NINT( (numgrid -1 ) / 2.)
        ALLOCATE (polar_ld( lambda,-sidebound : sidebound,-sidebound : sidebound,  &
             -1:numlogly) , STAT=astat)

     END IF

     IF(astat/=0) THEN

        WRITE(*,*) 'Error: array ''polar_ld'' not allocated!'
        erro = .TRUE.

     ELSE

        polar_ld(:,:,:,:) = 0

     END IF

     ! - polar_lu -

     IF (numgrid == 0) THEN

        ALLOCATE (polar_lu(lambda,0:0,0:0,-1:numlogly) , STAT=astat)
     ELSE

        sidebound = NINT( (numgrid -1 ) / 2.)
        ALLOCATE (polar_lu( lambda,-sidebound : sidebound,-sidebound : sidebound,  &
             -1:numlogly) , STAT=astat)

     END IF

     IF(astat/=0) THEN

        WRITE(*,*) 'Error: array ''polar_lu'' not allocated!'
        erro = .TRUE.

     ELSE

        polar_lu(:,:,:,:) = 0

     END IF

     ALLOCATE (norma(lambda), STAT=astat)

     IF(astat/=0) THEN

        WRITE(*,*) 'Error: array ''norma'' not allocated!'
        erro = .TRUE.

     ELSE

        norma(:)= 0.0

     END IF


     ! Open files for output

     OPEN(11,iostat=ios, file='aop.out', status='unknown')
     OPEN(12,iostat=ios, file='rad.out', status='unknown')

     ! If model is run in grid mode, data will also be written out to grid

     IF ( numgrid > 0 .AND. surfer == 0 ) THEN  !Grid files are not in Surfer ASCII grid format

        OPEN(14,iostat=ios, file='gridedo.out', status='unknown')
        OPEN(15,iostat=ios, file='grided.out', status='unknown')
		OPEN(16,iostat=ios, file='grideu.out', status='replace')
		OPEN(17,iostat=ios, file='grideuo.out', status='replace')


     END IF

     ! Start wavelength initialization

     LAMBDA_LOOP2: DO lam = 1, lambda

        ! ------------------------------------------------------------
        !                          - 7.2 -
        !                          Radiance
        ! ------------------------------------------------------------

        !INITIALIZE

        ai=0.
        sr=0.

        ! ALL PHOTONS THAT WERE BINNED AT BIN# 1 AND ALPHAINT WILL BE GROUPED
        ! INTO A SINGLE BIN. THIS WILL BE THE POLAR BINS (NORTH AND SOUTH)

        ALPHA_LOOP: DO  jj = 1, ( alphaint )

           IF  (angint == 1 ) THEN

              ai = ( PI * REAL(jj - 1) / alphaint ) + nalpha / 2.
              sr = SIN( ai ) * nalpha * nphi * ABS( COS( ai ) )
              rad( :,:,:,:,jj,: ) = REAL( n(:,:,:,:,jj,:) ) / sr

           ELSE IF ( jj /= 1 .AND. jj /= alphaint) THEN  

              ai =  ( ACOS(-(jj - REAL(alphaint)/2. -1) * muu) +   &
                   ACOS(-(jj - REAL(alphaint)/2.) * muu)) /2.
              sr = 2 * PI / ( ( alphaint / 2 - 1 )*phiint + 1)
              rad( :,:,:,:,jj,: ) = REAL( n( :,:,:,:,jj,:) ) / (sr * ABS( COS( ai))) 

           END IF

        END DO ALPHA_LOOP

        DO ii = -sidebound , sidebound
           DO jj =  -sidebound , sidebound
              DO kk = -1 , numlogly 

                 IF(angint == 1) THEN
                    polar_ld(lam,ii,jj,kk) = REAL( SUM( n(lam,ii,jj,kk,1,:)) ) / ( 2 * pi * (1 - COS( nalpha )) )
                    polar_lu(lam,ii,jj,kk) = REAL( SUM( n(lam,ii,jj,kk,alphaint,:)) ) / ( 2 * pi * (1 - COS( nalpha )) )

                    rad( lam,ii,jj,kk,1,: ) = REAL( n(lam,ii,jj,kk,1,:) ) / ( nphi * (1 - COS( nalpha )) )
                    rad( lam,ii,jj,kk,alphaint,: ) = REAL( n(lam,ii,jj,kk,alphaint,:) ) / ( nphi * (1 - COS( nalpha )) )
                 ELSE
                    ! For angint=0 (equal cosine), dmu = 2.0 / alphaint.
                    ! Solid angle of the entire polar cap = 2 * PI * dmu.
                    polar_ld(lam,ii,jj,kk) = REAL( SUM( n(lam,ii,jj,kk,1,:)) ) / (2.0 * PI * (2.0 / REAL(alphaint)))
                    polar_lu(lam,ii,jj,kk) = REAL( SUM( n(lam,ii,jj,kk,alphaint,:)) ) / (2.0 * PI * (2.0 / REAL(alphaint)))

                    ! Solid angle of a single bin in the polar cap = dmu * dphi = (2/alphaint) * (2*PI/phiint)
                    rad( lam,ii,jj,kk,1,: ) = REAL( n(lam,ii,jj,kk,1,:) ) / &
                         ((2.0 / REAL(alphaint)) * (2.0 * PI / REAL(phiint)))
                    rad( lam,ii,jj,kk,alphaint,: ) = REAL( n(lam,ii,jj,kk,alphaint,:) ) / &
                         ((2.0 / REAL(alphaint)) * (2.0 * PI / REAL(phiint)))
                 END IF

              END DO
           END DO
        END DO

        ! ------------------------------------------------------------
        !                          - 7.3 -
        !                         Eu and Ed 
        ! ------------------------------------------------------------

        DO ii = -sidebound , sidebound

           DO jj = -sidebound , sidebound

              DO kk = -1 , numlogly

                 eu(lam,ii,jj,kk) = SUM( n(lam,ii,jj,kk,(alphaint/2+1):alphaint,:) )
                 ed(lam,ii,jj,kk) = SUM( n(lam,ii,jj,kk,1:(alphaint/2),:) )

              END DO

           END DO

        END DO

        ! Compute value used for normalization

        IF (normh2o == 1)  THEN ! Normalize to below water irradiance

           norma(lam) = SUM( ed(lam,:,:,0) )

        ELSE

           norma(lam) = cositer(lam)

        END IF

        ! Normalize Irradiance

        eu(lam,:,:,:) = eu(lam,:,:,:)  / norma(lam)
        ed(lam,:,:,:) = ed(lam,:,:,:)  / norma(lam)

        ! ------------------------------------------------------------
        !                          - 7.4 -
        !                        Eou and Eod 
        ! ------------------------------------------------------------

        DO ii = -sidebound , sidebound

           DO jj = -sidebound , sidebound

              DO kk = -1, numlogly

                 DO ll = 1, alphaint

                    IF (angint == 1 ) THEN

                       ai=(PI * REAL(ll - 1) / alphaint ) + nalpha / 2.

                    ELSE IF (ll >= ( REAL(alphaint)/2. + 1.) .AND.      &
                         ( ll /= 1 .AND. ll /= alphaint )) THEN

                       ai = ACOS( ( (-(ll - REAL(alphaint)/2. -1) * muu) +    &
                            (-(ll - REAL(alphaint)/2.) * muu)) /2. )

                    ELSE IF ( ll /= 1 .AND. ll /= alphaint ) THEN

                       ai = ACOS( (((1-mum) - (ll-2)*muu) + ((1-mum) - (ll-1)*muu)) /2.)

                    END IF

                    IF (ll >= ( REAL(alphaint)/2. + 1.)) THEN

                       IF (ll == alphaint) THEN

                          eou(lam,ii,jj,kk) = SUM( n(lam,ii,jj,kk,ll,:) )       &
                               + eou(lam,ii,jj,kk)

                       ELSE

                          eou(lam,ii,jj,kk) = SUM( n(lam,ii,jj,kk,ll,:) ) /      &
                               ABS(COS( ai ))  + eou(lam,ii,jj,kk)

                       END IF

                    ELSE

                       IF (ll == 1) THEN

                          eod(lam,ii,jj,kk) = SUM( n(lam,ii,jj,kk,ll,:) )       &
                               + eod(lam,ii,jj,kk)

                       ELSE

                          eod(lam,ii,jj,kk) = SUM( n(lam,ii,jj,kk,ll,:) ) /     &
                               ABS(COS(ai))  + eod(lam,ii,jj,kk)

                       END IF

                    END IF

                 END DO

              END DO

           END DO

        END DO

        ! Now normalize Eou and Eod to Ed(0)

        eou(lam,:,:,:) = eou(lam,:,:,:)  / norma(lam)
        eod(lam,:,:,:) = eod(lam,:,:,:)  / norma(lam)

        ! ------------------------------------------------------------
        !                          - 7.5 -
        !                   Write AOP out to file
        ! ------------------------------------------------------------

                SELECT CASE(ios)
        
                CASE(0)
        
                   WRITE(11,'(14(A12,2x))')'depth (m)',' Eu ',' Ed ',' E ', ' R* ',' Eou ',     &
                        ' Eod ' , ' meancosUP ',' meancosDown',              &
                        ' meancos ', 'rd ',' rup', 'Ld', 'Lu'
        
           DO kk = -1 , numlogly
              ! Safely calculate AOPs to avoid division by zero
              eu_sum_aop = SUM( eu(lam,:,:,kk) )
              ed_sum_aop = SUM( ed(lam,:,:,kk) )
              eou_sum_aop = SUM( eou(lam,:,:,kk) )
              eod_sum_aop = SUM( eod(lam,:,:,kk) )
              efwd_plus_eback_aop = REAL(totalfwd(lam,kk) + totalback(lam,kk))

              ! R
              IF (ed_sum_aop > 0.0) THEN
                 R_aop = eu_sum_aop / ed_sum_aop
              ELSE
                 R_aop = 0.0
              END IF
              
              ! meancosUP
              IF (eou_sum_aop > 0.0) THEN
                 mu_up_aop = eu_sum_aop / eou_sum_aop
              ELSE
                 mu_up_aop = 0.0
              END IF
              
              ! meancosDown
              IF (eod_sum_aop > 0.0) THEN
                 mu_down_aop = ed_sum_aop / eod_sum_aop
              ELSE
                 mu_down_aop = 0.0
              END IF

              ! meancos
              IF ((eou_sum_aop + eod_sum_aop) > 0.0) THEN
                 mu_aop = (ed_sum_aop - eu_sum_aop) / (eou_sum_aop + eod_sum_aop)
              ELSE
                 mu_aop = 0.0
              END IF
              
              ! rd
              IF (fdown(lam,kk) > 0.0 .AND. efwd_plus_eback_aop > 0.0) THEN
                 rd_aop = (bdown(lam,kk)/fdown(lam,kk)) / &
                          (REAL(totalback(lam,kk)) / efwd_plus_eback_aop)
              ELSE
                 rd_aop = 0.0
              END IF

              ! rup
              IF (fup(lam,kk) > 0.0 .AND. efwd_plus_eback_aop > 0.0) THEN
                 rup_aop = (bup(lam,kk)/fup(lam,kk)) / &
                           (REAL(totalback(lam,kk)) / efwd_plus_eback_aop)
              ELSE
                 rup_aop = 0.0
              END IF
              
              ! polar_ld/u
              IF (norma(lam) > 0.0) THEN
                 pld_aop = SUM(polar_ld(lam,:,:,kk)) / norma(lam)
                 plu_aop = SUM(polar_lu(lam,:,:,kk)) / norma(lam)
              ELSE
                 pld_aop = 0.0
                 plu_aop = 0.0
              END IF

              WRITE(11,'(14(F12.6,2x))') layval(kk), eu_sum_aop, ed_sum_aop, (ed_sum_aop - eu_sum_aop), &
                R_aop, eou_sum_aop, eod_sum_aop, mu_up_aop, mu_down_aop, mu_aop, rd_aop, rup_aop, pld_aop, plu_aop

           END DO
           WRITE(11,*)'* The above water reflectance is actually an albedo calculation'

        CASE DEFAULT

           WRITE(*,*) 'Failed to open ''aop.out'' '
           WRITE(*,*) '  IOSTAT error code = ', ios
           erro = .TRUE.

        END SELECT

        ! ------------------------------------------------------------
        !                 Write radiance out to file
        ! ------------------------------------------------------------

                SELECT CASE(ios)

        

                CASE(0)

                

                   ALLOCATE(rad_vals_row(phiint), STAT=astat)

                   IF (astat /= 0) THEN

                      erro = .TRUE.

                      WRITE(*,*) "Error allocating rad_vals_row"

                      CYCLE LAMBDA_LOOP2

                   END IF

        

                   ! Create format string for header (phi angles)

                   WRITE(fmt,'(a,I0,a)') '(A12,2x,',phiint,'(F6.2,6x))'

                   ! Write header

                   WRITE(12,fmt) 'Angle(deg)', ( (nphi * REAL(ii-1) * 180.0 / PI), ii=1, phiint)

                   

                   ! Create format string for data rows

                   WRITE(fmt,'(a,I0,a)') '(F12.6,2x,',phiint,'(F11.7,1x))'

        

                   DO kk = -1 , numlogly

        

                      WRITE(12,'(A,F10.4,A)') 'Depth: ', layval(kk),' m'

        

                      DO ii = 1, alphaint

                         

                                          ! Calculate zenith angle for the current row

                         

                                          IF (angint == 1) THEN

                         

                                             ai = (nalpha * (REAL(ii - 1)) + nalpha/2.0) * 180.0/PI

                         

                                          ELSE

                         

                                             IF (ii == 1) THEN

                         

                                                ai = (ACOS(1.0 - 0.5 * (2.0/REAL(alphaint)))) * 180.0/PI

                         

                                             ELSE IF (ii == alphaint) THEN

                         

                                                ai = (ACOS(-1.0 + 0.5 * (2.0/REAL(alphaint)))) * 180.0/PI

                         

                                             ELSE IF (ii >= (REAL(alphaint)/2.0 + 1.0)) THEN

                         

                                                ai = ( ACOS(-(ii - REAL(alphaint)/2.0 - 1.0) * muu) + &

                         

                                                       ACOS(-(ii - REAL(alphaint)/2.0) * muu) ) / 2.0 * 180.0/PI

                         

                                             ELSE

                         

                                                ai = ( ACOS(((1.0-mum) - (ii-2.0)*muu) + ((1.0-mum) - (ii-1.0)*muu)) /2.0 )*180/PI

                         

                                             END IF

                         

                                          END IF

                         

                                          

                         

                                          ! Safely calculate radiance values for the row

                         

                                          IF (norma(lam) > 0.0) THEN

                         

                                             DO jj = 1, phiint

                         

                                                rad_vals_row(jj) = SUM(rad(lam,:,:,kk,ii,jj)) / norma(lam)

                         

                                             END DO

                         

                                          ELSE

                         

                                             rad_vals_row(:) = 0.0

                         

                                          END IF

                         

                         

                         

                                          ! Write the full row

                         

                                          WRITE(12,fmt) ai, rad_vals_row

        

                      END DO

        

                   END DO

        

                   DEALLOCATE(rad_vals_row, STAT=astat)

                   IF (astat /= 0) THEN

                      WRITE(*,*) "Error deallocating rad_vals_row"

                   END IF

        

                CASE DEFAULT

           WRITE(*,*) 'Failed to open ''aop.out'' '
           WRITE(*,*) '  IOSTAT error code = ', ios
           erro = .TRUE.

        END SELECT

        ! ------------------------------------------------------------
        !                 Write gridded data out to file
        ! ------------------------------------------------------------

        IF ( numgrid > 0) THEN

           IF (surfer == 1 ) THEN !Output in Surfer ASCII grid format

              WRITE(clamb,'(I7.7)') INT(wavelength(lam) * 1000)

              DO kk = -1 , numlogly

                 ! Open new file for each depth and each wavelength

                 IF ( kk == -1 ) THEN

                    cdepth = '0000a_'

                 ELSE IF ( kk == 0 ) THEN

                    cdepth = '0000w_'

                 ELSE

                    WRITE( cdepth , '(I5.5)' ) kk

                 END IF

                 fname = 'Eu_'//cdepth//'_'//clamb//'.grd'  
                 OPEN(14,iostat=ios, file=fname, status='unknown',RECL=1414)

                 fname = 'Ed_'//cdepth//'_'//clamb//'.grd' 
                 OPEN(15,iostat=ios, file=fname, status='unknown',RECL=1414)


                 WRITE(14,'(A4)')'DSAA'
                 WRITE(15,'(A4)')'DSAA'

                 ! number of grid lines along the X axis

                 WRITE(14,'(2(I4,2x))')numgrid , numgrid
                 WRITE(15,'(2(I4,2x))')numgrid , numgrid

                 ! min and max x value of grid

                 WRITE(14,'(2(F12.6,2x))') (- cellsize * sidebound) , (cellsize * sidebound)
                 WRITE(15,'(2(F12.6,2x))') (- cellsize * sidebound) , (cellsize * sidebound)

                 ! min and max y value of grid

                 WRITE(14,'(2(F12.6,2x))') (- cellsize * sidebound) , (cellsize * sidebound)
                 WRITE(15,'(2(F12.6,2x))') (- cellsize * sidebound) , (cellsize * sidebound)

                 ! min and max value for  entire grid

                 WRITE(14,'(2(F12.6,2x))') MINVAL( eu(lam,:,:,kk)) ,MAXVAL(eu(lam,:,:,kk))
                 WRITE(15,'(2(F12.6,2x))') MINVAL( ed(lam,:,:,kk)) ,MAXVAL(ed(lam,:,:,kk))

                 ! Write format info to string

                 WRITE(cformat,'(A1,I3,A11)') '(' , numgrid , '(E12.4,2x))'

                 DO jj = -sidebound , sidebound

                    WRITE(14,cformat) (eu(lam,ii,jj,kk), ii=-sidebound , sidebound)
                    WRITE(15,cformat) (ed(lam,ii,jj,kk) , ii=-sidebound , sidebound)

                 END DO

                 CLOSE(14)
                 CLOSE(15)


              END DO

           ELSE
                
			! Plain ASCII matrix output to already-open units 14 (Eod) and 15 (Ed)
			  ! One grid per depth per wavelength; include small headers for human readability.

			  ! Build a row format with numgrid values per line (engineering format for safety)
			  WRITE(rowfmt,'(A,I0,A)') '(', numgrid, '(ES16.8,1X))'
			  headfmt = '(A,1X,I0,1X,A,1X,F12.5,1X,A)'

			  ! Loop over depths and write both Eod and Ed grids
			  DO kk = -1 , numlogly
				
				
				! ---- Eod grid to unit 14 ----
				WRITE(14,'(A,1X,I0,1X,A,1X,F10.3,1X,A,1X,A,1X,F12.6,1X,A)') &
					 'GRID',  numgrid,                        &
					 'LAMBDA', wavelength(lam), 'nm',         &
					 'DEPTH', layval(kk),        'm'
				WRITE(14,'(A,1X,F12.6,1X,F12.6)') 'XMIN XMAX', -cellsize*sidebound, cellsize*sidebound
				WRITE(14,'(A,1X,F12.6,1X,F12.6)') 'YMIN YMAX', -cellsize*sidebound, cellsize*sidebound
				WRITE(14,'(A,1X,ES16.8,1X,ES16.8)') 'ZMIN ZMAX', MINVAL(eod(lam,:,:,kk)), MAXVAL(eod(lam,:,:,kk))
				DO jj = -sidebound , sidebound
				  WRITE(14,rowfmt) ( eod(lam,ii,jj,kk), ii = -sidebound , sidebound )
				END DO
			
				! ---- Ed grid to unit 15 ----
				WRITE(15,'(A,1X,I0,1X,A,1X,F10.3,1X,A,1X,A,1X,F12.6,1X,A)') &
					 'GRID',  numgrid,                        &
					 'LAMBDA', wavelength(lam), 'nm',         &
					 'DEPTH', layval(kk),        'm'
				WRITE(15,'(A,1X,F12.6,1X,F12.6)') 'XMIN XMAX', -cellsize*sidebound, cellsize*sidebound
				WRITE(15,'(A,1X,F12.6,1X,F12.6)') 'YMIN YMAX', -cellsize*sidebound, cellsize*sidebound
				WRITE(15,'(A,1X,ES16.8,1X,ES16.8)') 'ZMIN ZMAX', MINVAL(ed(lam,:,:,kk)), MAXVAL(ed(lam,:,:,kk))
				DO jj = -sidebound , sidebound
				  WRITE(15,rowfmt) ( ed(lam,ii,jj,kk), ii = -sidebound , sidebound )
				END DO
  
				! ---- Eu grid to unit 16 ----
				  WRITE(16,'(A,1X,I0,1X,A,1X,F10.3,1X,A,1X,A,1X,F12.6,1X,A)') &
					   'GRID', numgrid, 'LAMBDA', wavelength(lam), 'nm', 'DEPTH', layval(kk), 'm'
				  WRITE(16,'(A,1X,F12.6,1X,F12.6)') 'XMIN XMAX', -cellsize*sidebound, cellsize*sidebound
				  WRITE(16,'(A,1X,F12.6,1X,F12.6)') 'YMIN YMAX', -cellsize*sidebound, cellsize*sidebound
				  WRITE(16,'(A,1X,ES16.8,1X,ES16.8)') 'ZMIN ZMAX', MINVAL(eu(lam,:,:,kk)), MAXVAL(eu(lam,:,:,kk))
				  DO jj = -sidebound , sidebound
					WRITE(16,rowfmt) ( eu(lam,ii,jj,kk), ii = -sidebound , sidebound )
				  END DO

				  ! ---- Euo grid to unit 17 ----
				  WRITE(17,'(A,1X,I0,1X,A,1X,F10.3,1X,A,1X,A,1X,F12.6,1X,A)') &
					   'GRID', numgrid, 'LAMBDA', wavelength(lam), 'nm', 'DEPTH', layval(kk), 'm'
				  WRITE(17,'(A,1X,F12.6,1X,F12.6)') 'XMIN XMAX', -cellsize*sidebound, cellsize*sidebound
				  WRITE(17,'(A,1X,F12.6,1X,F12.6)') 'YMIN YMAX', -cellsize*sidebound, cellsize*sidebound
				  WRITE(17,'(A,1X,ES16.8,1X,ES16.8)') 'ZMIN ZMAX', MINVAL(eou(lam,:,:,kk)), MAXVAL(eou(lam,:,:,kk))
				  DO jj = -sidebound , sidebound
					WRITE(17,rowfmt) ( eou(lam,ii,jj,kk), ii = -sidebound , sidebound )
				  END DO


			  END DO

           END IF

        END IF

        ! End of bulk data output at each wavelength


     END DO LAMBDA_LOOP2

     CLOSE(11)
     CLOSE(12)

     IF ( numgrid > 0 .AND. surfer == 0) THEN

        CLOSE(14)
        CLOSE(15)
		CLOSE(16)
		CLOSE(17)

     END IF

     ! ------------------------------------------------------------
     !    Write data to file as a function of wavelength 
     ! ------------------------------------------------------------

     IF (lambda > 1) THEN

        OPEN(14,iostat=ios,FILE='wave.out',status = 'unknown')

        SELECT CASE(ios)

        CASE(0)

           WRITE(14,'(7(A12,2x))')'wavelength',' Eu ', ' Ed ', ' R ', &
                ' Lu ',' Ed(a) ', ' Lu(a) '

           DO lam = 1, lambda
              
              ed_sum_0 = SUM(ed(lam,:,:,0))
              norma_val = norma(lam)

              IF (ed_sum_0 > 0.0) THEN
                 R_val = SUM(eu(lam,:,:,0)) / ed_sum_0
              ELSE
                 R_val = 0.0
              END IF

              IF (norma_val > 0.0) THEN
                 Lu_val = SUM(polar_lu(lam,:,:,0)) / norma_val
                 Lua_val = SUM(polar_lu(lam,:,:,-1)) / norma_val
              ELSE
                 Lu_val = 0.0
                 Lua_val = 0.0
              END IF
              
              WRITE(14,'(F12.6, 6(2X,F13.7))') wavelength(lam), &
                   SUM( eu(lam,:,:,0)), ed_sum_0,      &
                   R_val,       &
                   Lu_val,  &
                   SUM( ed(lam,:,:,-1)),                         &
                   Lua_val
           END DO

        CASE DEFAULT

           WRITE(*,*)'Failed to open ''wave.out'' '
           erro = .TRUE.

        END SELECT

     END IF

     ! ------------------------------------------------------------
     !                   END of AOP computation
     ! ------------------------------------------------------------

  CASE DEFAULT      ! Errors occurred, simulation did not proceed

     WRITE(*,*) 'Errors occurred... Could not compute AOP''s'

  END SELECT

  ! ----------------------------------------------------------
  !                          - 8.0 -
  !                    Free up some memory
  ! ----------------------------------------------------------

  DEALLOCATE(n,STAT=astat)

  IF (astat /= 0) THEN

     WRITE(*,*)'ERROR: problem in array ''N'' deallocation'
     erro = .TRUE.

  END IF

  DEALLOCATE(rad,STAT=astat)

  IF (astat /= 0) THEN

     WRITE(*,*)'ERROR: problem in array ''rad'' deallocation'
     erro = .TRUE.

  END IF

  DEALLOCATE(eu,STAT=astat)

  IF (astat /= 0) THEN

     WRITE(*,*)'ERROR: problem in array ''eu'' deallocation'
     erro = .TRUE.

  END IF

  DEALLOCATE(ed,STAT=astat)

  IF (astat /= 0) THEN

     WRITE(*,*)'ERROR: problem in array ''ed'' deallocation'
     erro = .TRUE.

  END IF

  DEALLOCATE(polar_ld,STAT=astat)

  IF (astat /= 0) THEN

     WRITE(*,*)'ERROR: problem in array ''polar_ld'' deallocation'
     erro = .TRUE.

  END IF

  DEALLOCATE(polar_lu,STAT=astat)

  IF (astat /= 0) THEN

     WRITE(*,*)'ERROR: problem in array ''polar_lu'' deallocation'
     erro = .TRUE.

  END IF

  DEALLOCATE(wavelength,STAT=astat)

  IF (astat /= 0) THEN

     WRITE(*,*)'ERROR: problem in array ''wavelength'' deallocation'
     erro = .TRUE.

  END IF

  DEALLOCATE(layval,STAT=astat)

  IF (astat /= 0) THEN

     WRITE(*,*)'ERROR: problem in array ''layval'' deallocation'
     erro = .TRUE.

  END IF

  DEALLOCATE(fdown,STAT=astat)

  IF (astat /= 0) THEN

     WRITE(*,*)'ERROR: problem in array ''fdown'' deallocation'
     erro = .TRUE.

  END IF

  DEALLOCATE(fup,STAT=astat)

  IF (astat /= 0) THEN

     WRITE(*,*)'ERROR: problem in array ''fup'' deallocation'
     erro = .TRUE.

  END IF

  DEALLOCATE(bdown,STAT=astat)

  IF (astat /= 0) THEN

     WRITE(*,*)'ERROR: problem in array ''bdown'' deallocation'
     erro = .TRUE.

  END IF

  DEALLOCATE(bup,STAT=astat)

  IF (astat /= 0) THEN

     WRITE(*,*)'ERROR: problem in array ''bup'' deallocation'
     erro = .TRUE.

  END IF

  DEALLOCATE(totalback,STAT=astat)

  IF (astat /= 0) THEN

     WRITE(*,*)'ERROR: problem in array ''totalback'' deallocation'
     erro = .TRUE.

  END IF

  DEALLOCATE(totalfwd,STAT=astat)

  IF (astat /= 0) THEN

     WRITE(*,*)'ERROR: problem in array ''totalfwd'' deallocation'
     erro = .TRUE.

  END IF

  DEALLOCATE(cositer,STAT=astat)

  IF (astat /= 0) THEN

     WRITE(*,*)'ERROR: problem in array ''cositer'' deallocation'
     erro = .TRUE.

  END IF

  SELECT CASE (erro)

  CASE (.FALSE.)

     WRITE(*,*) 'AOMC model run complete!'

  CASE (.TRUE.)

     WRITE(*,*) '* Errors occured, simulation may not have properly terminated *'

  END SELECT


END PROGRAM mc

! ********************************************************************
! ********************************************************************
!
!      Subroutine to check the proper range of the INTEGER value
!
! ********************************************************************
! ********************************************************************

SUBROUTINE check_range_int(value,low,high,erro)

  IMPLICIT NONE

  INTEGER,INTENT(IN)  ::value,low,high
  LOGICAL,INTENT(INOUT) ::erro

  IF ( value<low .OR. value>high ) THEN

     erro = .TRUE.
     WRITE(*,*)"!!! Value is out of range !!!"

  END IF

END SUBROUTINE  check_range_int

! ********************************************************************
! ********************************************************************
!
!      Subroutine to check the proper range of the REAL value
!
! ********************************************************************
! ********************************************************************

SUBROUTINE check_range_real(value,low,high,erro)

  IMPLICIT NONE

  REAL,INTENT(IN)     ::value,low,high
  LOGICAL,INTENT(INOUT) ::erro

  IF ( value<low .OR. value>high ) THEN

     erro = .TRUE.
     WRITE(*,*)"!!! Value is out of range !!!"

  END IF

END SUBROUTINE  check_range_real

! ********************************************************************
! ********************************************************************
!
!      Subroutine to check if the value is and ODD number
!
! ********************************************************************
! ********************************************************************

SUBROUTINE check_odd(value,erro)

  IMPLICIT NONE

  INTEGER,INTENT(IN)     ::value
  LOGICAL,INTENT(INOUT)    ::erro

  IF (value > 0) THEN

     IF ( ABS(MOD(value,2)) < 0.01 ) THEN

        erro = .TRUE.
        WRITE(*,*)'!!! Value must be an ODD number !!!'

     END IF

  END IF

END SUBROUTINE  check_odd

! ********************************************************************
! ********************************************************************
!
!      Subroutine to check if the value is and EVEN number
!
! ********************************************************************
! ********************************************************************

SUBROUTINE check_even(value,erro)

  IMPLICIT NONE

  INTEGER,INTENT(IN)       ::value
  LOGICAL,INTENT(INOUT)    ::erro

  IF (value > 0) THEN

     IF ( ABS(MOD(value,2)) > 0.01 ) THEN

        erro = .TRUE.
        WRITE(*,*)'!!! Value must be an EVEN number !!!'

     END IF

  END IF

END SUBROUTINE  check_even

