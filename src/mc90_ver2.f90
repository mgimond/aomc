!     Last change:  MG   18 Dec 2025
PROGRAM mc
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
   REAL      :: randradius  !Randomly generated radius
   REAL      :: randangle   !Randomly generated angle
   REAL      :: sumfrac     !cumulitive summary of cumputed fracscat()
   REAL      :: ai          !Angle in rad
   REAL      :: alpha1, alpha2 ! Upper and lower bound of bins
   REAL      :: dmu
   REAL      :: mu_low, mu_high ! Bin parameters for angint=0
   REAL      :: sr          !Solid angle
   REAL      :: R_aop, mu_up_aop, mu_down_aop, mu_aop, rd_aop, rup_aop, pld_aop, plu_aop
   REAL      :: ed_sum_aop, eu_sum_aop, eou_sum_aop, eod_sum_aop, efwd_plus_eback_aop
   REAL      :: R_val, Lu_val, Lua_val, ed_sum_0, norma_val
   REAL, ALLOCATABLE, DIMENSION(:) :: rad_vals_row

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
   CALL SYSTEM_CLOCK (xseed,yseed,zseed)
   wseed = ABS(xseed - 19027983)

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
   nphi         = 2 * PI / REAL(phiint)
   muu          = phiint / ( (REAL(alphaint) /2. -1 )* phiint +1 )
   mum          = muu / phiint

   ! ----------------------------------------------------------
   !                          - 3.0 -
   !                    Allocation of memory
   ! ----------------------------------------------------------
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

   ALLOCATE (conc(numlycst,const), STAT=astat)
   IF(astat/=0) THEN
      WRITE(*,*) 'Error: array ''conc'' not allocated! Do you have enough memory?'
      erro = .TRUE.
   ELSE
      conc(:,:) = 0.0
   END IF

   ALLOCATE (specabs(lambda,const), STAT=astat)
   IF(astat/=0) THEN
      WRITE(*,*) 'Error: array ''specabs'' not allocated! Do you have enough memory?'
      erro = .TRUE.
   ELSE
      specabs (:,:) = 0.0
   END IF

   ALLOCATE (absorb(lambda,numlycst), STAT=astat)
   IF(astat/=0) THEN
      WRITE(*,*) 'Error: array ''absorb'' not allocated! Do you have enough memory?'
      erro = .TRUE.
   ELSE
      absorb (:,:) = 0.0
   END IF

   ALLOCATE (specscat(lambda,const), STAT=astat)
   IF(astat/=0) THEN
      WRITE(*,*) 'Error: array ''specscat'' not allocated! Do you have enough memory?'
      erro = .TRUE.
   ELSE
      specscat(:,:) = 0.0
   END IF

   ALLOCATE (scatter(lambda,numlycst), STAT=astat)
   IF(astat/=0) THEN
      WRITE(*,*) 'Error: array ''scatter'' not allocated! Do you have enough memory?'
      erro = .TRUE.
   ELSE
      scatter (:,:) = 0.0
   END IF

   ALLOCATE (fracscat(lambda,numlycst,const), STAT=astat)
   IF(astat/=0) THEN
      WRITE(*,*) 'Error: array ''fracscat'' not allocated! Do you have enough memory?'
      erro = .TRUE.
   ELSE
      fracscat (:,:,:) = 0.0
   END IF

   ALLOCATE (atten(lambda,numlycst), STAT=astat)
   IF(astat/=0) THEN
      WRITE(*,*) 'Error: array ''atten'' not allocated! Do you have enough memory?'
      erro = .TRUE.
   ELSE
      atten (:,:) = 0.0
   END IF

   ALLOCATE (salbedo(lambda,numlycst), STAT=astat)
   IF(astat/=0) THEN
      WRITE(*,*) 'Error: array ''salbedo'' not allocated! Do you have enough memory?'
      erro = .TRUE.
   ELSE
      salbedo (:,:) = 0.0
   END IF

   ALLOCATE (spf(spfang,0:const), STAT=astat)
   IF(astat/=0) THEN
      WRITE(*,*) 'Error: array ''spf'' not allocated! Do you have enough memory?'
      erro = .TRUE.
   ELSE
      spf (:,:)= 0.0
   END IF

   ALLOCATE (slopespf(2:spfang,1:const), STAT=astat)
   IF(astat/=0) THEN
      WRITE(*,*) 'Error: array ''slopespf'' not allocated! Do you have enough memory?'
      erro = .TRUE.
   ELSE
      slopespf (:,:)= 0.0
   END IF

   ALLOCATE (refr(0:numlycst), STAT=astat)
   IF(astat/=0) THEN
      WRITE(*,*) 'Error: array ''refr'' not allocated! Do you have enough memory?'
      erro = .TRUE.
   ELSE
      refr (:)= 1.0
   END IF

   ALLOCATE (dptlycst(numlycst), STAT=astat)
   IF(astat/=0) THEN
      WRITE(*,*) 'Error: array ''dptlycst'' not allocated! Do you have enough memory?'
      erro = .TRUE.
   ELSE
      dptlycst (:)= 1.0
   END IF

   ALLOCATE (wavelength(lambda), STAT=astat)
   IF(astat/=0) THEN
      WRITE(*,*) 'Error: array ''wavelength'' not allocated! Do you have enough memory?'
      erro = .TRUE.
   ELSE
      wavelength (:)= 0.0
   END IF

   ALLOCATE (bottomr(lambda), STAT=astat)
   IF(astat/=0) THEN
      WRITE(*,*) 'Error: array ''bottomr'' not allocated! Do you have enough memory?'
      erro = .TRUE.
   ELSE
      bottomr (:)= 0.0
   END IF

   ALLOCATE (bottomspc(lambda), STAT=astat)
   IF(astat/=0) THEN
      WRITE(*,*) 'Error: array ''bottomspc'' not allocated! Do you have enough memory?'
      erro = .TRUE.
   ELSE
      bottomspc (:)= 0.0
   END IF
   SELECT CASE (targbot)
   CASE (1)
      ALLOCATE (targref(lambda), STAT=astat)
      IF(astat/=0) THEN
         WRITE(*,*) 'Error: array ''targref'' not allocated! Do you have enough memory?'
         erro = .TRUE.
      ELSE
         targref (:)= 0.0
      END IF
      ALLOCATE (targspc(lambda), STAT=astat)
      IF(astat/=0) THEN
         WRITE(*,*) 'Error: array ''targspc'' not allocated! Do you have enough memory?'
         erro = .TRUE.
      ELSE
         targspc (:)= 0.0
      END IF
   END SELECT
   IF (lsource == 0 )THEN
      ALLOCATE (direct(lambda), STAT=astat)
      IF(astat/=0) THEN
         WRITE(*,*) 'Error: array ''direct'' not allocated! Do you have enough memory?'
         erro = .TRUE.
      ELSE
         direct (:)= 0.0
      END IF
   END IF

   ALLOCATE (intensity(lambda), STAT=astat)
   IF(astat/=0) THEN
      WRITE(*,*) 'Error: array ''intensity'' not allocated! Do you have enough memory?'
      erro = .TRUE.
   ELSE
      intensity (:)= 0.0
   END IF

   ALLOCATE (fdown(lambda,-1:numlogly), STAT=astat)
   IF(astat/=0) THEN
      WRITE(*,*) 'Error: array ''fdown'' not allocated! Do you have enough memory?'
      erro = .TRUE.
   ELSE
      fdown (:,:)= 0.0
   END IF

   ALLOCATE (fup(lambda,-1:numlogly), STAT=astat)
   IF(astat/=0) THEN
      WRITE(*,*) 'Error: array ''fup'' not allocated! Do you have enough memory?'
      erro = .TRUE.
   ELSE
      fup (:,:)= 0.0
   END IF

   ALLOCATE (bdown(lambda,-1:numlogly), STAT=astat)
   IF(astat/=0) THEN
      WRITE(*,*) 'Error: array ''bdown'' not allocated! Do you have enough memory?'
      erro = .TRUE.
   ELSE
      bdown (:,:)= 0.0
   END IF

   ALLOCATE (bup(lambda,-1:numlogly), STAT=astat)
   IF(astat/=0) THEN
      WRITE(*,*) 'Error: array ''bup'' not allocated! Do you have enough memory?'
      erro = .TRUE.
   ELSE
      bup (:,:)= 0.0
   END IF

   ALLOCATE (totalback(lambda,-1:numlogly), STAT=astat)
   IF(astat/=0) THEN
      WRITE(*,*) 'Error: array ''totalback'' not allocated! Do you have enough memory?'
      erro = .TRUE.
   ELSE
      totalback (:,:)= 0.0
   END IF

   ALLOCATE (totalfwd(lambda,-1:numlogly), STAT=astat)
   IF(astat/=0) THEN
      WRITE(*,*) 'Error: array ''totalfwd'' not allocated! Do you have enough memory?'
      erro = .TRUE.
   ELSE
      totalfwd (:,:)= 0.0
   END IF

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
   OPEN(1,iostat=ios, file='conc.inp', status='old')
   SELECT CASE(ios)
   CASE(0)
      WRITE(*,*) 'File ''conc.inp''... successfully opened'
      READ(1,*)records
      IF (records == numlycst ) THEN
         DO ii=1,numlycst
            READ(1,*)(conc(ii,jj),jj=1,const),refr(ii),dptlycst(ii)
            IF (dptlycst(ii) > depthb) THEN
               WRITE(*,*) 'ERROR!! boundary of layer',ii,'is greater than depth.'
               erro = .TRUE.
            END IF
         END DO
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

   OPEN(2,iostat=ios, file='abs.inp', status='old')
   SELECT CASE(ios)
   CASE(0)
      WRITE(*,*) 'File ''abs.inp''... successfully opened'
      READ(2,*)records
      IF (records == lambda ) THEN
         DO ii=1,lambda
            READ(2,*)(specabs(ii,jj),jj=1,const)
         END DO
         absorb(:,:)  = 0.0
         DO jj=1,numlycst
            DO ii=1,lambda
               DO kk=1,const
                  absorb(ii,jj)=absorb(ii,jj)+specabs(ii,kk)*conc(jj,kk) 
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

   OPEN(3,iostat=ios, file='scat.inp', status='old')
   SELECT CASE(ios)
   CASE(0)
      WRITE(*,*) 'File ''scat.inp''... successfully opened'
      READ(3,*)records
      IF (records == lambda) THEN
         DO ii=1,records
            READ(3,*)(specscat(ii,jj),jj=1,const)
         END DO
         DO jj=1,numlycst
            DO ii=1,lambda
               DO kk=1,const
                  scatter(ii,jj)=scatter(ii,jj)+specscat(ii,kk)*conc(jj,kk)
               END DO
            END DO
         END DO
         DO jj=1,numlycst
            DO ii=1,lambda
               sumfrac = 0.0
               DO kk=1,const
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

   OPEN(4,iostat=ios, file='spf.inp', status='old')
   SELECT CASE(ios)
   CASE(0)
      WRITE(*,*) 'File ''spf.inp''... successfully opened'
      READ(4,*)records
      IF (records == spfang) THEN
         DO ii=1,spfang
            READ(4,*)(spf(ii,jj),jj=0,const)
         END DO
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
               sky_int(ii,:,:) = NINT(sky_int(ii,:,:) / SUM (sky_int(ii,:,:))* &
                    iter)
               total = 0
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
   atten(:,:) = absorb(:,:) + scatter(:,:)
   salbedo(:,:) = scatter(:,:) /  atten(:,:)
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
         CASE DEFAULT
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
            CALL light(wavematch,lam,jj,wavecontrol)
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
            CALL logbin(x(1),y(1),-1,lam)
            cositer(lam) = cositer(lam) + 1 / COS( theta )
            CALL interface_sub(refr(0),refr(1),switch)
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
   CASE DEFAULT
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

      OPEN(11,iostat=ios, file='aop.out', status='unknown')
      OPEN(12,iostat=ios, file='rad.out', status='unknown')
      IF ( numgrid > 0 .AND. surfer == 0 ) THEN  !Grid files are not in Surfer ASCII grid format
         OPEN(14,iostat=ios, file='gridedo.out', status='unknown')
         OPEN(15,iostat=ios, file='grided.out', status='unknown')
         OPEN(16,iostat=ios, file='grideu.out', status='replace')
         OPEN(17,iostat=ios, file='grideuo.out', status='replace')
      END IF
      LAMBDA_LOOP2: DO lam = 1, lambda
         ! ------------------------------------------------------------
         !                          - 7.2 -
         !                          Radiance
         ! ------------------------------------------------------------
         ai=0.
         sr=0.
         ALPHA_LOOP: DO  jj = 1, ( alphaint )
            IF  (angint == 1 ) THEN
               alpha1 = (jj-1) * nalpha 
               alpha2 = jj * nalpha   
            ELSE  ! angint == 0
               mu_high = 1.0 - (REAL(jj-1) * 2.0 / REAL(alphaint))
               mu_low = 1.0 - (REAL(jj) * 2.0 / REAL(alphaint))
               alpha1 = ACOS(mu_high)
               alpha2 = ACOS(mu_low)
            END IF  

            IF (alpha2 <= PI/2.0) THEN
               sr = 0.5 * ( SIN(alpha2)**2 - SIN(alpha1)**2 ) * nphi
            ELSEIF (alpha1 >= PI/2.0) THEN
               sr = 0.5 * ( SIN(alpha1)**2 - SIN(alpha2)**2 ) * nphi
            ELSE
               sr = 0.5 * ( 2.0 - SIN(alpha1)**2 - SIN(alpha2)**2 ) * nphi
            END IF
            rad( :,:,:,:,jj,: ) = REAL( n(:,:,:,:,jj,:) ) / sr   
         END DO ALPHA_LOOP
         DO ii = -sidebound , sidebound
            DO jj =  -sidebound , sidebound
               DO kk = -1 , numlogly 
                  polar_ld(lam,ii,jj,kk) = SUM( rad(lam,ii,jj,kk,1,:) ) / REAL(phiint)
                  polar_lu(lam,ii,jj,kk) = SUM( rad(lam,ii,jj,kk,alphaint,:) ) / REAL(phiint)
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

         IF (normh2o == 1)  THEN ! Normalize to below water irradiance
            norma(lam) = SUM( ed(lam,:,:,0) )
         ELSE
            norma(lam) = cositer(lam)
         END IF
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