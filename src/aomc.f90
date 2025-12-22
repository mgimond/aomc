!     Last change:  MG   11 Apr 2003    3:07 pm
SUBROUTINE geom2 (thetascat, phiscat, pathlength)

  ! Program to calculate new position of photon based
  ! on the scatter angle
  ! 
  ! The following lists the input parameters:
  ! 
  !      X0,Y0,Z0 = previous position of photon
  !      X1,Y1,Z1 = current position of photon
  !      *THETA   = current angle between trajectory and Z-axis
  !                 (this value is calculated in this subroutine)
  !      *PHI     = current angle between trajectory along XY-plane
  !                 (this value is calculated in this subroutine)
  !      THETASCAT    = angle of scatter about current trajectory
  !      PHISCAT  = angle of rotation of scattered light
  !                 about trajectory
  !      PATHLENGTH = length of new trajectory
  ! 
  ! The following lists the output parameters:
  ! 
  !      X,Y,Z    = new loaction of photon
  !      THETA    = angle between new trajectory and z-axis
  ! 
  ! The routine works as follows:
  !
  !      1) translate last 2 coordinate positions of the photon
  !        ((X0,Y0,Z0) and(X1,Y1,Z1)) such that
  !        (X0,Y0,Z0) is at the origin of the master coordinate
  !        system. The 2 resulting points are P0 and P1.
  !
  !                | 0 |
  !          P0 =  | 0 |
  !                | 0 |
  ! 
  !               | (X1-X0) |
  !          P1 = | (Y1-Y0) |
  !               | (Z1-Z0) |
  ! 
  !     2) use existing theta angle (angle between line
  !        segment P0 and P1 and z-axis of master coordinate
  !        system). Add the scatter angle thetascat to theta.
  !        At this point, we are assuming that phiscat (the
  !        angle of rotation about the current path trajectory)
  !        is 0.
  !
  !     3) Calculate the interim photon position (P2) based on the
  !        new angle theta (at this point we still have not
  !        included the angle about the trajectory axis).
  !
  !     4) rotate the interim postion P2 about the trajectory
  !        axis by the angle phiscat. The resulting position is P3.
  !
  !     5) translate the position P3 back by adding P0. This
  !        gives us our new point location of the photon.
  ! 
  ! 
  !                                      | (P3(1)+ X0) |
  !         New position of photon  =    | (P3(2)+ Y0) |
  !                                      | (P3(3)+ Z0) |
  ! 
  ! WARNING: (X0,Y0,Z0) must not equal (X1,Y1,Z1).
  ! 
  ! 
  ! NOTE: This routine was developed with the assumption that the
  ! PHISCAT angle distribution is uniform across all 360 degrees. This
  ! routine was not thouroughly tested for a non-uniform phiscat
  ! distribution function.


  USE math_global
  USE propagation_global

  IMPLICIT NONE

  REAL,INTENT(IN)     :: phiscat,thetascat
  REAL,INTENT(IN)     :: pathlength
  REAL,DIMENSION(1:3) :: p1,p2,p3
  REAL                :: h,s,dx,dy,dz
  REAL                :: xt,yt,zt

     ! First, lets translate the photon's last point ( P(0)=0 )

     p1(1) = x(1) - x(0)
     p1(2) = y(1) - y(0)
     p1(3) = z(1) - z(0)

     theta = theta + thetascat

     !  Now calculate interim postion P2 given the new theta angle
     !  The angle phiscat is set to 0.0 at this point

     p2(1) = pathlength * SIN( theta ) * COS( phi ) + p1(1)
     p2(2) = pathlength * SIN( theta ) * SIN( phi ) + p1(2)
     p2(3) = pathlength * COS( theta ) + p1(3)

     ! Now we will rotate this new point about the arbitrary axis

     IF ( ( ABS(p1(1) - p2(1)) < 0.00001).AND.( ABS(p1(2) - p2(2)) < 0.00001) &
          .AND. (ABS(p1(3) - p2(3)) < 0.00001) ) THEN

        p3(1) = p2(1)
        p3(2) = p2(2)
        p3(3) = p2(3)

     ELSE

        CALL rot3d (p1, phiscat, p2, p3)

     END IF

     xt = p3(1) + x(0)
     yt = p3(2) + y(0)
     zt = p3(3) + z(0)
     dx = ( xt - x(1))
     dy = ( yt - y(1))
     dz = ( zt - z(1))
     x(0) = x(1)
     y(0) = y(1)
     z(0) = z(1)
     x(1) = xt
     y(1) = yt
     z(1) = zt

     h = SQRT( dx * dx + dy * dy + dz * dz )
     s = SQRT( dx * dx + dy * dy )
     theta = ACOS( dz / h )

     IF (s /= 0.0) THEN
        IF (dx >= 0.0) THEN
           phi = ASIN( dy / s )
        ELSE
           phi= PI - ASIN( dy / s )
        END IF
     ELSE
        phi = 0.0
     END IF

END SUBROUTINE geom2


SUBROUTINE rot3d (axis, aangle, vv, ww)

  ! ***********************************************************************
  ! 
  ! ROTATION_AXIS_VECTOR_3D rotates a vector around an axis vector in 3D.
  !
  ! This subroutine was originally written by John Burkardt (7/31/1999).
  ! http://www.psc.edu/~burkardt/src/geometry/geometry.html
  ! 
  ! The following lists the input parameters:
  !            AXIS(3) = the axis vector for the rotation.
  !            ANGLE   = the angle, in radians, of the rotation.
  !            V(3)    = the vector to be rotated.
  ! 
  ! The following lists the output parameters:
  !            W(3)    = the rotated vector.
  ! ***********************************************************************
  !
  IMPLICIT NONE

  REAL,DIMENSION(3),INTENT(IN)    :: axis,vv
  REAL,INTENT(IN)                 :: aangle
  REAL,DIMENSION(3),INTENT(INOUT) :: ww
  REAL :: dot,norm,norm2,xa,xn
  REAL :: xn2,xp,xr,ya,yn,yn2,yp,yr,za,zn,zn2,zp,zr
  ! 
  !      Compute the length of the rotation axis.
  ! 
  xa = axis(1)
  ya = axis(2)
  za = axis(3)

  norm = SQRT( xa * xa + ya * ya + za * za)


  IF1: IF (norm /= 0.0 ) THEN

     xa = xa / norm
     ya = ya / norm
     za = za / norm

     ! 
     !      Compute the dot product of the vector and the rotation axis.
     ! 

     dot = vv(1) * xa + vv(2) * ya + vv(3) * za

     ! 
     !      Compute the parallel component of the vector.
     ! 

     xp = dot * xa
     yp = dot * ya
     zp = dot * za

     ! 
     !      Compute the normal component of the vector.
     ! 

     xn = vv(1) - xp
     yn = vv(2) - yp
     zn = vv(3) - zp

     norm2 = SQRT(xn * xn + yn * yn + zn * zn)


     IF (norm2 == 0.0) THEN

        ww(1) = xp
        ww(2) = yp
        ww(3) = zp


     ELSE


        xn= xn / norm2
        yn= yn / norm2
        zn= zn / norm2

        ! 
        !      Compute a second vector, lying in the plane, perpendicular
        !      to V, and forming a right-handed system.
        ! 
        !    The cross product in 3D can be regarded as the determinant of the
        !    symbolic matrix:
        ! 
        !          |  i  j  k |
        !      det | x1 y1 z1 |
        !          | x2 y2 z2 |
        ! 
        !      = ( y1 * z2 - z1 * y2 ) * i
        !      + ( z1 * x2 - x1 * z2 ) * j
        !      + ( x1 * y2 - y1 * x2 ) * k

        xn2 = ya * zn - za * yn
        yn2 = za * xn - xa * zn
        zn2 = xa * yn - ya * xn

        norm = SQRT( xn2 * xn2 + yn2 * yn2 + zn2 * zn2 )

        xn2 = xn2 / norm
        yn2 = yn2 / norm
        zn2 = zn2 / norm

        ! 
        !      Rotate the normal component by the angle.
        ! 

        xr = norm2 * ( COS( aangle ) * xn + SIN( aangle ) * xn2 )
        yr = norm2 * ( COS( aangle ) * yn + SIN( aangle ) * yn2 )
        zr = norm2 * ( COS( aangle ) * zn + SIN( aangle ) * zn2 )

        ! 
        !      The rotated vector is the parallel component plus the 
        !      rotated com
        ! 

        ww(1) = xp + xr
        ww(2) = yp + yr
        ww(3) = zp + zr

     END IF

  END IF IF1

END SUBROUTINE rot3d

!     Last change:  MG   20 Mar 2003    3:09 pm

! *******************************************************
! **            Mathematical constants                 **
! *******************************************************

MODULE math_global
  IMPLICIT NONE
  REAL, PARAMETER                :: PI = 3.14159265
END MODULE math_global

! *******************************************************
! **     Random number generator global variables      **
! *******************************************************

MODULE rand_global
  IMPLICIT NONE
  INTEGER                        :: xseed,yseed,zseed,wseed
END MODULE rand_global

! *******************************************************
! **     System constituents characteristics           **
! *******************************************************

MODULE const_global

  IMPLICIT NONE
  INTEGER                         :: lambda  !Total number of wavebands

  ! The following parameters are wavelength dependent

  REAL,ALLOCATABLE,DIMENSION(:,:) :: absorb  !Bulk absorption coefficient in each layer
  REAL,ALLOCATABLE,DIMENSION(:,:) :: scatter !Bulk scattering coefficient in each layer
  REAL,ALLOCATABLE,DIMENSION(:,:) :: atten   !Bulk beam attenuation coefficient
                                             !atten = absorb + scatter
  REAL,ALLOCATABLE,DIMENSION(:,:)  ::salbedo !Scattering albedo
  REAL,ALLOCATABLE,DIMENSION(:,:,:)::fracscat!Fraction of scattering probability associated
                                             !with each constituent at each wavelength and layer
                                             !Sum at each wavelength and layer must equal 1
  REAL,ALLOCATABLE,DIMENSION(:)   :: refr    !Refractive index in each layer
  INTEGER                         :: spftyp  !SPF internally generated (0) or from file(1)
  REAL,ALLOCATABLE,DIMENSION(:,:) :: spf     !Cumulative scattering phase function for
                                             !each constituent. The first element (0) of
                                             !the 2nd dimension is reserved for the angle (rad)
  REAL,ALLOCATABLE,DIMENSION(:,:) :: slopespf!Slope between each data point of SPF							   
  INTEGER                         :: const   !Number of constituents in the system

END MODULE const_global

! *******************************************************
! **          System's physical characteristics        **
! *******************************************************

MODULE physical_global

  IMPLICIT NONE
  REAL                             ::depthb    !Depth of bottom boundary
  REAL,ALLOCATABLE,DIMENSION(:)    ::bottomr   !Bottom reflectance
  REAL,ALLOCATABLE,DIMENSION(:)    ::bottomspc !Specular fraction of reflecting bottom
  REAL                             ::sideb     !Width of side boundary
  INTEGER                          ::numlycst  !Number of depth dependent layers
  REAL,ALLOCATABLE,DIMENSION(:)    ::dptlycst  !Depth of each depth dependent layers
  INTEGER                          ::targbot   !Bottom target switch
  REAL,ALLOCATABLE,DIMENSION(:)    ::targref   !Reflection fraction of bottom target
  REAL,ALLOCATABLE,DIMENSION(:)    ::targspc   !Specular fraction of reflecting target
  REAL                             ::targx     !Width (x) of target
  REAL                             ::targy     !Width (y) of target

END MODULE physical_global

! *******************************************************
! **          Light boundary characteristics           **
! *******************************************************

MODULE light_global

  IMPLICIT NONE
  INTEGER                      :: lsource   !Skylight distribution (0=model, 1=tabular)
  REAL,ALLOCATABLE,DIMENSION(:):: intensity !Intensity of incident light at a wavelength
  REAL,ALLOCATABLE,DIMENSION(:):: direct    !Fraction of light collimated or direct
                                            !(non-tabular data)

  ! Following arrays are used if tabulated sky light distribution is used
  ! for simulation

  REAL,ALLOCATABLE,DIMENSION(:)        :: sky_wave
  REAL,ALLOCATABLE,DIMENSION(:,:,:)    :: sky_int
  INTEGER,ALLOCATABLE,DIMENSION(:,:,:) :: sky_nint    !Cumulative incident photon distribution
  INTEGER                              :: sky_zenith  !Number of zenith intervals
  INTEGER                              :: sky_azimuth !Number of azimuth intervals

  ! Following variables are used if internal sky light distribution is used
  ! for simulation

  INTEGER                      :: zenith  !Zenith angle of point source (for non-tabular data)
  INTEGER                      :: azimuth !Azimuth angle of point source (for non-tabular data)
  INTEGER                      :: lgttyp  !Type of light footprint (point,circle or rectangle)
  REAL                         :: lgtdiam !Diameter of light footprint if circular source
  REAL                         :: lgtx    !Width (x) of light footprint if rectangular source
  REAL                         :: lgty    !Width (y) of light footprint if rectangular source

END MODULE light_global

! *******************************************************
! **           Logging parameters                     **
! *******************************************************

MODULE log_global

  IMPLICIT NONE
  INTEGER         :: numlogly  !Number of logging (recording) layers
  REAL            :: intlogly  !Thickness of each logging (recording)
                               !layers. A value of -1.0 indicates that
                               !the interval is internally calculated
  INTEGER         :: alphaint  !Number of logging bins in polar direction
  REAL            :: nalpha    !Size of alpha intervals
  INTEGER         :: phiint    !Number of logging bins in azimuthal direction
  REAL            :: nphi      !Size of phi intervals
  REAL            :: muu       !mu spacing in polar direction
  REAL            :: mum       !mu spacing for polar caps
  INTEGER         :: polar     !Value of 1 if polar cap is to be used for binning
  INTEGER         :: angint    !1= equal angular intervals, 0 = equal meancosine intervals 
  INTEGER         :: numgrid   !Number of grid cells in x and y direction
                               !A value of 0 indicates that output is not
                               !gridded
  REAL            :: cellsize  !Cell size (x=y) of grid

END MODULE log_global

! *******************************************************
! **           Logging environment                     **
! *******************************************************

MODULE propagation_global

  IMPLICIT NONE
  REAL                :: theta     !Polar angle of propagation of photon
  REAL                :: phi       !azimuth angle of propagation of photon
  REAL,DIMENSION(0:1) :: x,y,z     !x, y and z position of photon
  REAL,ALLOCATABLE,DIMENSION(:)              :: layval ! Depth value for each layer interface
  INTEGER, PARAMETER :: k8i = SELECTED_INT_KIND(14)    ! Allows for values greater than 32 bit
  INTEGER, PARAMETER :: k8r = SELECTED_REAL_KIND(18)   ! Allows for values greater than 32 bit
  INTEGER(K8i),ALLOCATABLE,DIMENSION(:,:,:,:,:,:) :: n ! Counter for photon's crossing
                                                       ! of an interface
  INTEGER(k8i),ALLOCATABLE,DIMENSION(:,:) :: totalfwd  ! Total number of photons scattered
                                                       ! in the forward direction
  INTEGER(K8i),ALLOCATABLE,DIMENSION(:,:) :: totalback ! Total number of photons scattered
                                                       ! in the backward direction
  REAL(K8r)                               :: totalop   !Total optical pathlength traveled
  REAL(K8r),ALLOCATABLE,DIMENSION(:,:)    :: fdown,fup,bdown,bup   ! Shape factors

END MODULE  propagation_global

! *******************************************************
! **                  AOP quantities                   **
! *******************************************************

MODULE aop_global

  IMPLICIT NONE
  REAL,ALLOCATABLE,DIMENSION(:,:,:,:,:,:) :: rad    ! Computed radiance
  REAL,ALLOCATABLE,DIMENSION(:,:,:,:)     :: eu     ! Upwelling irradiance
  REAL,ALLOCATABLE,DIMENSION(:,:,:,:)     :: ed     ! Downwelling irradiance
  REAL,ALLOCATABLE,DIMENSION(:,:,:,:)     :: eou    ! Upwelling scalar irradiance
  REAL,ALLOCATABLE,DIMENSION(:,:,:,:)     :: eod    ! Downwelling scalar irradiance

END MODULE aop_global
!     Last change:  MG   20 Mar 2003    3:06 pm
SUBROUTINE INTERFACE_SUB(na,nb,switch)
  USE math_global
  USE rand_global
  USE randmod
  USE anglesum_mod
  USE const_global
  USE propagation_global
  IMPLICIT NONE

  ! Purpose:
  ! --------
  !   To determine the outcome of the photon's interaction with
  !   the air-water interface.
  !
  ! Description
  ! -----------
  ! THE FOLLOWING ROUTINE WILL DETERMINE IF THE PHOTON  WILL
  ! REFLECT OFF THE AIR-WATER INTERACE. IF IT DOES NOT REFLECT OFF OF
  ! THE INTERACE, THIS SUBROUTINE WILL CALCULATE THE CHANGE
  ! IN POLAR ANGLE OF THE PHOTON ONCE IT PASSES THE INTERFACE
  ! If direction =0, the photon is traveling from the
  ! air medium to the water medium.
  !
  ! Record and revisions
  ! --------------------
  !     Date     Programmer           Description of change
  !     ====     ==========           =====================
  !  3/20/2003  Manuel Gimond         Original code


  INTEGER,INTENT(OUT) :: switch !0 = not reflected, 1 = reflected
  REAL,INTENT(IN)     :: na,nb
  REAL                :: term1,term2,term3,term4
  REAL                :: r
  REAL                :: theta_refl
  REAL                :: anglemax

  switch = 0

  IF ( ( nb - na ) > (0.00001) ) THEN    ! Photon is traveling from low density to high density

     ! THE FOLLOWING DETERMINES THE REFRACTED ANGLE OF THE PHOTON IN THE WATER
     ! THE ONLY ANGLE THAT WILL CHANGE IS THE POLAR ANGLE

     theta_refl = ASIN((na/nb)*SIN(theta))

     ! THE FOLLOWING EQUATION WILL DETERMINE IF THE INCIDENT PHOTON
     ! PENETRATES THE  AIR/WATER INTERFACE.
     ! IF THE ANGLE OF INCIDENCE IS EQUAL TO 0, THE REFLECTANCE
     ! VALUE IS SET TO 2% (SINCE THE ANGLE VALUE 0 IS BEYOND
     ! THE DOMAIN OF FRESNEL'S EQUATION)

     term1=(SIN( theta - theta_refl ))**2
     term2=(SIN( theta + theta_refl ))**2
     term3=(TAN( theta - theta_refl ))**2
     term4=(TAN( theta + theta_refl ))**2

     IF (nb /= 1.0 .AND. theta == 0.) THEN

        r = 0.02

     ELSE IF ( ABS(nb - na ) < 0.000001) THEN

        r = 0.00

     ELSE

        r = 0.5 * ( term1 / term2 ) + 0.5 * ( term3  /term4 )

     END IF

     ! THE VALUE R JUST CALCULATED IS THE PERCENT PROBABILITY OF THE INCIDENT
     ! REFLECTANCE OFF OF THE WATER SURFACE.  A CALL TO RAND() WILL BE MADE TO
     ! A RANDOM NUMBER WHICH WILL BE USED IN CONJUNCTION WITH R IN DETERMININ
     ! PHOTON IS TO REFLECT OFF OF THE SURFACE OR TRANSMIT THROUGH

     IF ( rand() > r ) THEN  ! photon enters water

        theta = theta_refl
        phi = anglesum(1,phi,PI)

     ELSE

        theta = PI - theta
        phi = anglesum(1,phi,PI)
        switch = 1

     END IF

  ELSE IF ( ( na - nb ) > (0.00001) ) THEN        ! The photon is going from water to air

     ! Check for the critical angle of reflection (i.e. the angle at which
     ! the photon will scatter of the interface 100% of the time)

     theta = PI - theta  ! the reflected photon's angle needs to be adjusted

     anglemax = ASIN( refr(0) /refr(1) )

     IF ( ( theta ) < anglemax ) THEN

        ! THE FOLLOWING DETERMINES THE REFRACTED ANGLE OF THE PHOTON IN THE AIR
        ! THE ONLY ANGLE THAT WILL CHANGE IS THE POLAR ANGLE

        theta_refl = ASIN((refr(1) / refr(0)) * SIN( theta ))

        ! THE FOLLOWING EQUATION WILL DETERMINE IF THE INCIDENT PHOTON
        ! PENETRATES THE  AIR/WATER INTERFACE.
        ! IF THE ANGLE OF INCIDENCE IS EQUAL TO 0, THE REFLECTANCE
        ! VALUE IS SET TO 2% (SINCE THE ANGLE VALUE 0 IS BEYOND
        ! THE DOMAIN OF FRESNEL'S EQUATION)

        TERM1 = ( SIN( theta - theta_refl ))**2
        TERM2 = ( SIN( theta + theta_refl ))**2
        TERM3 = ( TAN( theta - theta_refl ))**2
        TERM4 = ( TAN( theta + theta_refl ))**2

        IF (theta /= 0) THEN

           R = 0.5 * (TERM1 / TERM2) + 0.5 * (TERM3 / TERM4)

        ELSE

           R = 0.02

        END IF

        ! THE VALUE R JUST CALCULATED IS THE PERCENT PROBABILITY OF THE INCIDENT
        ! REFLECTANCE OFF OF THE WATER SURFACE.  A CALL TO RAND() WILL BE MADE TO
        ! A RANDOM NUMBER WHICH WILL BE USED IN CONJUNCTION WITH R IN DETERMININ
        ! PHOTON IS TO REFLECT OFF OF THE SURFACE OR TRANSMIT THROUGH

        IF (rand() > r) THEN        ! photon is transmitted through interface

           theta = PI - theta_refl   ! Relative to below water coordinate system

        ELSE

           switch = 1

        END IF


     END IF

  ELSE ! Index of refractions are equal, light crosses interface

     theta = theta
     phi = anglesum(1,phi,PI)

  END IF

END  SUBROUTINE interface_sub
!     Last change:  MG   26 Mar 2003   10:55 am
SUBROUTINE light(wavenum,wavel,iteration,control)

  ! Purpose:
  ! --------
  !   To determine the angle of incidence of the photon.
  !
  ! Description
  ! -----------
  !   This subroutine will determine if the photon's angle of
  !   incidence is randomly generated (diffuse skylight) or
  !   given the zenith and azimuth angle of the collimated
  !   (or direct) source.
  !
  ! Record and revisions
  ! --------------------
  !     Date     Programmer           Description of change
  !     ====     ==========           =====================
  !   3/30/2003  Manuel Gimond        Original code

  USE rand_global
  USE randmod
  USE math_global
  USE propagation_global
  USE light_global

  IMPLICIT NONE
  INTEGER,INTENT(IN)       :: wavenum      ! Wavelength ID
  INTEGER,INTENT(IN)       :: iteration    ! Iteration number
  INTEGER,INTENT(IN)       :: wavel        ! Wavelength if photon direction is
                                           ! internally generated
  INTEGER,INTENT(INOUT)    :: control      ! Check for change in wavelength
  INTEGER                  :: ii,jj        ! Loop index
  INTEGER,SAVE             :: starti = 1, startj = 0 ! Start of Loop index

  IF (lsource == 1 ) THEN

     IF (control == 1) THEN

        control = 0
        starti = 1
        startj = 0

     END IF

  END IF

  SELECT CASE (lsource)

  CASE (0)  ! photon direction is internally computed

     IF (rand() < direct(wavel)) THEN  ! photon is from direct/collimated source

        theta = zenith  * PI / 180
        phi   = azimuth * PI / 180

     ELSE

        phi   = rand() * 2 * PI
        theta = PI / 2 - ASIN( rand() )  

     END IF

  CASE (1)  ! photon direction is read from file

     ZEN_LOOP:  DO ii = starti, sky_zenith
        AZ_LOOP: DO jj = startj, sky_azimuth

           IF( sky_nint(wavenum,ii,jj )>= iteration) THEN

              theta = REAL(ii) * PI / 180
              phi = REAL(jj) * PI /180
              starti = ii
              startj = jj

              EXIT ZEN_LOOP

           ELSE 

              startj = 0

           END IF

        END DO AZ_LOOP

     END DO ZEN_LOOP

  END SELECT

END SUBROUTINE light
!     Last change:  MG   18 Dec 2025   
SUBROUTINE logbin(xint,yint,k,lam)

  USE math_global
  USE propagation_global
  USE log_global
  USE physical_global

  IMPLICIT NONE

  ! Purpose:
  ! --------
  !   To log the crossing of a photon across a layer interface
  !
  ! Description
  ! -----------
  !  If a photon crosses a layer interface where a recording event
  !  is to occur, the photons angles of propagation (theta and phi)
  !  are recorded along with the grid cell (i,j,k) being crossed.
  !
  ! Record and revisions
  ! --------------------
  !     Date     Programmer           Description of change
  !     ====     ==========           =====================
  !  12/31/2002  Manuel Gimond        Original code
  !  12/18/2025  Manuel Gimond        Simplified binning method for angint=0
  
  INTEGER,INTENT(IN)   :: k    ! ID of layer interface
  REAL,INTENT(IN)      :: xint ! X position on interval k
  REAL,INTENT(IN)      :: yint ! Y position on interval k
  INTEGER, INTENT(IN)  :: lam  ! Wavelength
  INTEGER              :: i,j,t,p,limit

  limit = 0

  ! Determine polar and azimuthal bins

  IF (theta < PI) THEN

     limit = 1

  END IF

  ! ---------------------------- SELECT OPTION -------------------------------
  IF ( angint == 1) THEN
     ! For angint=1, bins are of equal theta angle size.
     t = INT( theta / nalpha ) + limit
  ELSE  ! angint == 0
     ! For angint=0, bins are of equal cosine(theta) size.
     ! Bin is determined by mapping cos(theta) from [1, -1] to [1, alphaint].
     IF (theta >= PI) THEN
        t = alphaint
     ELSE
        t = INT( (1.0 - COS(theta)) * (REAL(alphaint) / 2.0) ) + 1
     END IF
  END IF

  p = INT( phi / nphi ) + 1        ! since phi is never 2PI, no need to check 'limit'


  ! ------------------------ END OF SELECT OPTION ----------------------------

  IF (numgrid == 0) THEN

     IF ( (ABS( xint ) < ABS ( sideb / 2)) .AND. &
          (ABS( yint ) < ABS ( sideb / 2)) ) THEN

        IF ( phi < 0. .OR. phi >= ( 2* PI) ) THEN

           WRITE(*,*) theta, phi , k

        END IF

        n(lam,0,0,k,t,p) = n(lam,0,0,k,t,p) + 1

     END IF

  ELSE

     IF ( ((ABS(xint)) < (numgrid * cellsize / 2. )) .AND.   &   ! Only log event occuring within boundary
          ((ABS(yint)) < (numgrid * cellsize / 2. )) ) THEN

        i = NINT( xint / cellsize )
        j = NINT( yint / cellsize )
        n(lam,i,j,k,t,p) = n(lam,i,j,k,t,p) + 1

     END IF

  END IF

END SUBROUTINE logbin
!     Last change:  MG   19 Dec 2025    
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
  !  12/19/2025  Manuel Gimond		  Changed radiance calculations, reformat aop.out precision
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
       '       Version 1.2.1                                 ',&
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
  nphi         = 2 * PI / REAL(phiint)
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
     
        ! Compute polar radiances Lu and Ld
        DO ii = -sidebound , sidebound
           DO jj =  -sidebound , sidebound
              DO kk = -1 , numlogly 
!               polar_ld(lam,ii,jj,kk) = SUM( rad(lam,ii,jj,kk,1,:) ) / REAL(phiint)
!				 polar_lu(lam,ii,jj,kk) = SUM( rad(lam,ii,jj,kk,alphaint,:) ) / REAL(phiint)
          IF(angint == 1) THEN                                                     
             polar_ld(lam,ii,jj,kk) = REAL( SUM( n(lam,ii,jj,kk,1,:)) ) / ( 2 * pi * (1.0 - COS( nalpha )) )                                                           
             polar_lu(lam,ii,jj,kk) = REAL( SUM( n(lam,ii,jj,kk,alphaint,:)) ) / ( 2 * pi * (1.0 - COS( nalpha )) )                                                    
                                                                                   
!             rad( lam,ii,jj,kk,1,: ) = REAL( n(lam,ii,jj,kk,1,:) ) / ( nphi * (1.0 - COS( nalpha )) )                                                                  
!             rad( lam,ii,jj,kk,alphaint,: ) = REAL( n(lam,ii,jj,kk,alphaint,:) ) / ( nphi * (1.0 - COS( nalpha )) )                                                    
          ELSE                                                                     
             ! For angint=0 (equal cosine), dmu = 2.0 / alphaint. Solid angle of the entire polar cap = 2 * PI * dmu.                 
             polar_ld(lam,ii,jj,kk) = REAL( SUM( n(lam,ii,jj,kk,1,:)) ) / (2.0 * PI * (2.0 / REAL(alphaint)))                                                        
             polar_lu(lam,ii,jj,kk) = REAL( SUM( n(lam,ii,jj,kk,alphaint,:)) ) / ( 2.0 * PI * (2.0 / REAL(alphaint)))                                                  
                                                                                   
             ! Solid angle of a single bin in the polar cap = dmu * dphi = (2/alphaint) * (2*PI/phiint)                                                        
 !            rad( lam,ii,jj,kk,1,: ) = REAL( n(lam,ii,jj,kk,1,:) ) / &             
 !                 ((2.0 / REAL(alphaint)) * (2.0 * PI / REAL(phiint)))             
 !            rad( lam,ii,jj,kk,alphaint,: ) = REAL( n(lam,ii,jj,kk,alphaint,:) ) / &                                                                                   
 !                 ((2.0 / REAL(alphaint)) * (2.0 * PI / REAL(phiint)))             
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

              WRITE(11,'(14(F12.8,2x))') layval(kk), eu_sum_aop, ed_sum_aop, (ed_sum_aop - eu_sum_aop), &
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

                   WRITE(fmt,'(a,I0,a)') '(F12.6,2x,',phiint,'(ES11.4,1x))'

        

                   DO kk = -1 , numlogly

        

                      WRITE(12,'(A,F10.4,A)') 'Depth: ', layval(kk),' m'

        

                      DO ii = 1, alphaint

                         

                                          ! Calculate zenith angle for the current row

                         

                                          IF (angint == 1) THEN
                                          
										     ai = (nalpha * (REAL(ii - 1)) + nalpha/2.0) * 180.0/PI
                      
                                          ELSE
                                             dmu = 2.0/REAL(alphaint)
                                             ai = ACOS(1.0 - (REAL(ii - 1) + 0.5) * dmu) * 180.0/PI

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
                  DO ll = 1, alphaint!     Last change:  MG   20 Mar 2003    2:44 pm

! *******************************************************
! **                    MODULE randmod                **
! *******************************************************

MODULE randmod 

! Purpose:
! --------
!   To generate a random number between 0 and 1.
!
! Description
! -----------
!   This  random number generator combines:  
!   (1) The congruential generator x(n)=69069*x(n-1)+1327217885, period 2^32.
!   (2) A 3-shift shift-register generator, period 2^32-1,                
!   (3) Two 16-bit multiply-with-carry generators, period 597273182964842497>2^59
!       Overall period>2^123;  Default seeds x,y,z,w.                        
!
! Record and revisions
! --------------------
!     Date     Programmer           Description of change
!     ====     ==========           =====================
!     1993     George Marsaglia     Original code in F77
!     2000     Mike Metcalf         Converted code to F90
!   1/10/2003  Manuel Gimond        Modified code to ouptut
!                                   values between 0 and 1
!     

IMPLICIT NONE
PRIVATE
PUBLIC :: rand

CONTAINS  
                       
   FUNCTION rand ()

      REAL :: rand

! This  random number generator combines:  
! (1) The congruential generator x(n)=69069*x(n-1)+1327217885, period 2^32.
! (2) A 3-shift shift-register generator, period 2^32-1,                
! (3) Two 16-bit multiply-with-carry generators, period 597273182964842497>2^59
!  Overall period>2^123;  Default seeds x,y,z,w.                        
!  Set your own seeds with statement i=kisset(ix,iy,iz,iw).             
!     
!      xseed = 69069 * xseed + 1327217885 
!      yseed = m (m (m (yseed, 13), - 17), 5) 
!      zseed = 18000 * IAND (zseed, 65535) + ISHFT (zseed, - 16) 
!      wseed = 30903 * IAND (wseed, 65535) + ISHFT (wseed, - 16) 
!      rand = ABS(xseed + yseed + ISHFT (zseed, 16) + wseed)/2147483648. !modified by MG
      CALL random_number(rand)

!   CONTAINS
!
!      FUNCTION m(k, n)
!         IMPLICIT NONE
!         INTEGER :: m, k, n
!
!         m = IEOR (k, ISHFT (k, n) )
!
!      END FUNCTION m

   END FUNCTION rand 

END MODULE randmod


! *******************************************************
! **                    MODULE anglesum                **
! *******************************************************

MODULE anglesum_mod

! Purpose:
! --------
!   To calculate the resulting angle between [0,PI] or [0,2PI[
!   from the summation of two angles.
!
! Description
! -----------
!   If the angles to be summed are polar, the output will result
!   in an angle ranging between 0 (downward) and PI (upward)
!   Else, the resulting output will range between 0 and 2PI.
!   Input angles : radian
!   Output angle : radian
!
!
! Record and revisions
! --------------------
!     Date     Programmer           Description of change
!     ====     ==========           =====================
! 12/27/2002   Manuel Gimond        Original code

IMPLICIT NONE
PRIVATE
PUBLIC :: anglesum

CONTAINS
  FUNCTION anglesum(switch,orig,add)
    REAL    :: anglesum
    INTEGER,INTENT(IN) :: switch ! 0=> poleward , 1=>horizontal
    REAL,INTENT(IN)    :: orig   ! original angle(in radian)
    REAL,INTENT(IN)    :: add    ! angle to add  (in radian)
    REAL, PARAMETER    :: PI = 3.14159265
    INTEGER :: tip


    SELECT CASE (switch)

    CASE (0)  !Polar angle (range [0,PI])

       anglesum = ACOS(COS(orig + add))

    CASE DEFAULT ! (range [0,2PI[

       IF ( ( ABS((2. * PI) - (orig+add))  < 0.0000001 ) .AND. ((orig+add) >= 0.)) THEN

          anglesum = orig + add
          tip=1

       ELSE IF ( ABS( COS(orig + add) - 1) < 0.0000001 ) THEN

          anglesum = 0.0
          tip=2

       ELSE IF ( SIN(orig + add) < 0) THEN

          anglesum = 2 * PI - ACOS( COS(orig + add) )
          tip=3

       ELSE IF ( SIN(orig + add) > 0 ) THEN

          anglesum = ACOS( COS(orig + add) )
          tip=4

       ELSE IF ( ABS( COS(orig + add) + 1) < 0.001 )  THEN

          anglesum = PI
          tip=5

       ELSE

          anglesum = 0.
          tip=6

       END IF

       IF( (anglesum - 2 * PI) > 0.0) THEN
          
          ! If sum of angle in the azimuth direction is greater than 2PI, 
          ! then let the user know. Angle should alwayd be less than 2PI.
          WRITE(*,*) 'BAD SUMMATION IN anglesum ROUTINE ', 'orig = ',orig, 'new = ', anglesum, tip
          anglesum = anglesum - ABS(anglesum - 2 * PI)  !! angle must no be greater than 2PI
      
       END IF

       IF( ABS(anglesum - 2 * PI) < 0.0001) THEN
          anglesum = 0.0
       END IF

    END SELECT

  END FUNCTION anglesum
END MODULE anglesum_mod


!     Last change:  MG   20 Mar 2003    2:46 pm
SUBROUTINE vsf(lam,thetascat,coptlayer)

  ! Purpose:
  ! --------
  !   To randomly generate an angle based on the scattering phase function.
  !
  ! Description
  ! -----------
  !
  !   The subroutine will generate a random polar angle for a scattered 
  !   photon based on a tabulated scattering phase function tabulated
  !   in the input file "spf.inp".
  !   The input file must be setup such that it is read as:
  !   "Probability of scatter up to, and including, angle thetascat"
  !
  ! Record and revisions
  ! --------------------
  !     Date     Programmer           Description of change
  !     ====     ==========           =====================
  !  3/10/2003   Manuel Gimond        Original code 


  USE rand_global
  USE randmod
  USE const_global

  IMPLICIT NONE

  INTEGER             :: ii,jj
  INTEGER, INTENT(IN) :: lam
  INTEGER, INTENT(IN) :: coptlayer
  REAL,INTENT(OUT)    :: thetascat
  REAL :: random
  REAL :: dthetascat

  thetascat = 0.0
  jj = 2

  TYPE_IF:  IF (spftyp == 1) THEN     ! Use SPF data read from file

     ! Determine which constituent is scattering for this interaction

     random = rand()

     FRACSCAT_DO: DO ii = 1, const

        IF (random <= fracscat(lam, coptlayer, ii)) THEN

           EXIT FRACSCAT_DO

        END IF

     END DO FRACSCAT_DO

     random = rand()

     IF (random <= spf(1,ii) ) THEN

        dthetascat = spf(1,0) 
        thetascat = rand() * dthetascat

     ELSE

        ANGLE_LOOP: DO

           IF (random > spf(jj - 1, ii) .AND. random <= spf(jj,ii) ) THEN

              thetascat = spf(jj,0) 

              ! A random value between 2 discrete SPF values can be generated.
              ! Uncomment the follwoing 2 lines if this option is desired
              !
              !           dthetascat = spf(jj, 0) - spf(jj - 1, 0)
              !           thetascat  = rand() * dthetascat + spf(jj - 1, 0)
              !
              ! End of option

              EXIT  ANGLE_LOOP

           ELSE

              jj = jj + 1

           END IF

        END DO ANGLE_LOOP

     END IF

  ELSE

     thetascat = ACOS( rand() * 2 - 1 )

  END IF TYPE_IF

END SUBROUTINE vsf
!     Last change:  MG   25 Nov 2025
!     Last change:  MG   11 Apr 2003    2:18 pm
SUBROUTINE water(lam,action)

  ! Purpose:
  ! --------
  !   To trace the path taken by a photon within a fluid medium
  !
  ! Description
  ! -----------
  ! this subroutine will trace the path of a single photon
  ! through the medium. At each interaction (with a
  ! boundary or a constituent of the medium), the x,y and z
  ! position of the photon is tracked.
  !
  ! Record and revisions
  ! --------------------
  !     Date     Programmer           Description of change
  !     ====     ==========           =====================
  !  3/12/2003  Manuel Gimond        Original code

  USE rand_global
  USE randmod
  USE math_global
  USE const_global
  USE anglesum_mod
  USE propagation_global
  USE physical_global
  USE log_global

  IMPLICIT NONE

  INTEGER,INTENT(IN)   :: lam     ! Wavelength
  INTEGER,INTENT(OUT)  :: action  ! Type of interaction
  INTEGER  :: optlayer            ! Optical layer within which the photon is traveling
  INTEGER  :: layer               ! Layer interface
  INTEGER  :: ii
  INTEGER  :: switch              ! Status at air-water interface (not needed in this sub)
  INTEGER  :: coptlayer           !actual optical layer location of photon 
  
  CHARACTER(LEN=11)   :: state

  REAL     :: op  ! Optical pathlength between each event.
                  ! (In this subroutine, an event can be (1) an
                  ! interaction with a physical boundary, (2)
                  ! absorption, (3) scattering, or (4) a change in 
                  ! IOP of the layer)
  REAL     :: intx, inty ! x,y,z values of intersection with a boundary
  REAL     :: random    ! Temporary palceholder for random number
  REAL     :: xb,yb     ! point location of photon's interaction with bottom
  REAL     :: thetascat ! scattering angle determined from VSF
  REAL     :: phiscat   ! scattering angle randomly generated
  REAL     :: xs,ys     ! X,Y,Z location on the surface
  REAL     :: b,b2
  REAL     :: theta_orig

  ! initialize variables

  action = 0
  x(0)=x(1)
  y(0)=y(1)
  z(0)=z(1)
  op = 0.0
  layer = 1
  optlayer = 1
  coptlayer = 1
  thetascat = 0.0

  ! The following will determine the max optical pathlength
  ! traveled by the photon before interacting with a particle.

  op = LOG( rand() ) / atten ( lam, coptlayer ) * ( -1 )

  ! The following will randomly pick a distance along the max OP.
  ! the resulting location of the photon is then calculated. this
  ! is not part of the loop. it is for use only with the first 
  ! appearance of the photon in the water body.

  x(1) =  x(0) + (op * SIN( theta ) * COS( phi )) 
  y(1) =  y(0) + (op * SIN( theta ) * SIN( phi )) 
  z(1) =  (op * COS( theta ) )

  ! The following loop tracks the photon until it is gone
  ! from the system. the position of the photon in space is
  ! also calculated. 

  ! ********************* depth intervals **************************
  ! Has the photon crossed any of the depth intervals (markers) 
  ! specified by the user? if so, the angles theta and phi and their
  ! associated optical pathlength will have to be recorded. Also, if
  ! the model is running in heterogeneous mode, the optical
  ! pathlength will have to be adjusted for the new inherent optical
  ! properties associated with that depth.
  ! ****************************************************************

  PHOTONLIFE_LOOP: DO

     state = 'attenuation'     !Default

     ! Does the photon cross a layer of differing  optical property

     IF (numlycst > 1) THEN

        IF ( theta >= PI/2.) THEN  !Photon traveling up

           OPTIC_LAYER_UP2: DO optlayer = numlycst  ,2 ,-1

              IF ( z(1) < dptlycst(optlayer) .AND. dptlycst(optlayer) < z(0) ) THEN

                 CALL olayer(lam,optlayer,op)

              END IF

           END DO OPTIC_LAYER_UP2

        ELSE  ! Photon traveling down

           OPTIC_LAYER_DN2: DO optlayer = 1 , (numlycst - 1)

              IF ( z(1) > dptlycst(optlayer+1) .AND. dptlycst(optlayer+1) > z(0) ) THEN

                 CALL olayer(lam,optlayer,op)

              END IF

           END DO OPTIC_LAYER_DN2

        END IF

        ! Determine the location of the photon vis-a-vis the optical layers

        DO ii = numlycst , 2, -1

           IF ( z(1) > dptlycst(ii) ) THEN

              coptlayer = ii
              EXIT

           END IF

        END DO

     END IF

     ! Check that all logging is conducted as needed

     LAYER_LOOP: DO

        ! Downward traveling photon

        LAYER_IF: IF ( z(1) > layval(layer) .AND.  layval(layer)> z(0) ) THEN

           b = ABS( layval( layer ) - z(0)) * TAN( theta )
           intx = x(0) + b * COS(phi)
           inty = y(0) + b * SIN(phi)

           ! Does the photon reach one of the side boundaries

           IF ( (ABS( intx) > ( sideb / 2.)) .OR. (ABS( inty) > ( sideb / 2.)) ) THEN
              state = 'side'
              EXIT PHOTONLIFE_LOOP
           END IF

           ! Now log photon crossing

           CALL logbin(intx,inty,layer,lam)

           IF (layer >= numlogly) THEN

              IF (z(1) >= depthb) THEN

                 state = 'bottom'

              END IF

              EXIT LAYER_LOOP

           END IF

           layer = layer + 1

           ! -------------- Upward traveling photon -------------------------

        ELSE IF  ( (z(1) < layval(layer)) .AND. (layval(layer) <= z(0)) ) THEN

           b = ABS( layval( layer ) - z(0)) * TAN(PI - theta)
           intx = x(0) + b * COS(phi)
           inty = y(0) + b * SIN(phi)

           ! Does the photon reach one of the side boundaries

           IF ( (ABS( intx) > ( sideb / 2.)) .OR. (ABS( inty) > ( sideb / 2.)) ) THEN

              state = 'side'
              EXIT PHOTONLIFE_LOOP

           END IF

           ! Now log photon crossing

           IF (layer > 0 ) THEN

              CALL logbin(intx,inty,layer,lam)
              layer = layer - 1

           ELSE  ! layer == 0

              state = 'surface'
              EXIT LAYER_LOOP

           END IF

           ! -------------- No more logging interfaces crossed ------------

        ELSE 

           IF ( (ABS( x(1)) > ( sideb / 2.)) .OR. (ABS( y(1)) > ( sideb / 2.)) ) THEN

              state = 'side'

           END IF

           IF (( theta >= PI/2.) .AND. (z(1) <= 0.) ) THEN

              state = 'surface'

           END IF

           IF (( theta < PI/2.) .AND. (z(1) >= depthb) ) THEN

              state = 'bottom'

           END IF

           EXIT LAYER_LOOP

        END IF LAYER_IF

     END DO LAYER_LOOP

     ! Does the photon reach the bottom

     IF (z(1) > depthb) THEN

        state = 'bottom'

     END IF

     ! Sum the OP traveled by the photon

     totalop = totalop + op

     ! ----------------------------------------------------------
     !              Select type of interaction between photon
     !                      and medium/boundary
     ! ----------------------------------------------------------

     SELECT CASE (state)

        ! ----------------------------------------------------------
        !                  Bottom boundary is reached
        ! ----------------------------------------------------------

     CASE ('bottom')

        ! Calculate x,y and z location where photon collides with bottom

        b = ABS( depthb - z(0) ) * TAN( theta )
        xb = x(0) + b * COS( phi )
        yb = y(0) + b * SIN( phi )

        ! Calculate the remaining optical pathlength

        op = op - SQRT( B**(2) + (depthb - z(0))**(2) )

        ! Determine outcome of photon's collision with bottom
        
        IF ( (targbot == 1) .AND. ( ABS(xb) < targx/2.) .AND. &
           ( ABS(yb) < targy/2.) ) THEN 
           
           IF (rand() <= targref(lam) ) THEN    !photon is reflected off target
              
              IF (rand() > targspc(lam) ) THEN  !photon is reflected isotropically
                 
                 theta = ACOS( rand() * 2. - 1 ) / 2. + PI / 2.  !new reflected angle
                 phi = anglesum(1, (rand() * 2 * PI) , 0.)
                 
                 x(0) = xb
                 y(0) = yb
                 z(0) = depthb
                 
                 x(1) =  (op * SIN( theta ) * COS( phi )) + xb
                 y(1) =  (op * SIN( theta ) * SIN( phi )) + yb
                 z(1) =  depthb + (op * COS( theta ) )
                 
              ELSE                                 ! Photon is reflected specularly
                 
                 theta = PI - theta
                 
                 x(0) = xb
                 y(0) = yb
                 z(0) = depthb
                 
                 x(1) =  (op * SIN( theta ) * COS( phi )) + xb
                 y(1) =  (op * SIN( theta ) * SIN( phi )) + yb
                 z(1) =  depthb + (op * COS( theta ) )
                 
              END IF
              
           ELSE                                    ! Photon is absorbed by target
              
              x(1) = xb
              y(1) = yb
              z(1) = depthb
              totalop = totalop - op  ! Remove untraveled path from total OP

              EXIT PHOTONLIFE_LOOP

           END IF
           
        ELSE

           IF (rand() <= bottomr(lam) ) THEN      !photon is reflected off bottom

              IF (rand() > bottomspc(lam) ) THEN  !photon is reflected isotropically

                 theta = ACOS( rand() * 2. - 1 ) / 2. + PI / 2.  !new reflected angle
                 phi = anglesum(1, (rand() * 2 * PI) , 0.)

                 x(0) = xb
                 y(0) = yb
                 z(0) = depthb

                 x(1) =  (op * SIN( theta ) * COS( phi )) + xb
                 y(1) =  (op * SIN( theta ) * SIN( phi )) + yb
                 z(1) =  depthb + (op * COS( theta ) )

              ELSE                                 ! Photon is reflected specularly

                 theta = PI - theta

                 x(0) = xb
                 y(0) = yb
                 z(0) = depthb

                 x(1) =  (op * SIN( theta ) * COS( phi )) + xb
                 y(1) =  (op * SIN( theta ) * SIN( phi )) + yb
                 z(1) =  depthb + (op * COS( theta ) )

              END IF

           ELSE                                    ! Photon is absorbed by bottom

              x(1) = xb
              y(1) = yb
              z(1) = depthb
              totalop = totalop - op  ! Remove untraveled path from total OP
              EXIT PHOTONLIFE_LOOP
           END IF
           
        END IF

        ! ----------------------------------------------------------
        !                  Surface boundary is reached
        ! ----------------------------------------------------------

     CASE ('surface')

        b  = ABS( z(0) ) * TAN( PI - theta )
        xs = x(0) + b * COS( phi )
        ys = y(0) + b * SIN( phi )
        op = op - SQRT( z(0)**2 + b**2) ! remaining OP if photon remains in water

        CALL logbin(xs,ys,0,lam)    ! Log the photon as it 'hits' the interface
        CALL interface_sub(refr(1),refr(0),switch)


        IF (theta >= PI/2.) THEN     !photon is not reflected back into the medium

           z(1) = 0.0
           x(1) = xs
           y(1) = ys
           CALL logbin(x(1),y(1),-1,lam)
           totalop = totalop - op    !Adjust total OP
           EXIT PHOTONLIFE_LOOP      !photon has left the medium

        ELSE         ! photon is reflected back into medium

           b2 = (ABS( z(1) ) + z(0) )  * TAN( theta )
           x(0)= xs
           y(0)= ys
           z(0)= 0.

           x(1) = (b2 - b) * COS( phi ) + xs
           y(1) = (b2 - b) * SIN( phi ) + ys
           z(1)= (ABS(COS( theta )) * op )

           CALL logbin(x(0),y(0),0,lam)
           layer = 1

        END IF

        ! ----------------------------------------------------------
        !                  Side boundary is reached
        ! ----------------------------------------------------------

     CASE ('side')

        EXIT PHOTONLIFE_LOOP

        ! ----------------------------------------------------------
        !         Photon is absorbed or scattered by the medium
        ! ----------------------------------------------------------

     CASE ('attenuation')

        ! Is the photon scattered or absorbed

        IF (rand() <= salbedo(lam,coptlayer)) THEN

           ! ******************* scatter *******************
           !

           op = LOG(rand()) / atten(lam,coptlayer) * (-1)  !Determine new optical pathlength
           theta_orig = theta
           CALL vsf(lam,thetascat,coptlayer)          !Determine angle of scatter
                                                      !relative to incident direction.
           IF (COS( theta_orig ) > 0.) THEN

              IF(thetascat < (PI / 2.)) totalfwd(lam, layer - 1) = &
              totalfwd(lam, layer - 1) + 1   !log % sf
              IF(thetascat > (PI / 2.)) totalback(lam, layer - 1) = &
              totalback(lam, layer - 1) + 1 !log % sb

           ELSE

              IF(thetascat < (PI / 2.)) totalfwd(lam,layer) = &
              totalfwd(lam,layer) + 1   !log % bf
              IF(thetascat > (PI / 2.)) totalback(lam,layer) = &
              totalback(lam,layer) + 1 !log % bb

           END IF


           PHI_LOOP: DO

              random = rand()

              IF (random /= 1.0) THEN

                 phiscat = random * 2 * PI
                 EXIT PHI_LOOP

              END IF

           END DO PHI_LOOP

           !    Determine angle of scatter relative to the parent
           !    coordinate system.

           IF( op > 0.0001) THEN  ! If pathlength is very small,
                                  ! assume no movement from the photon
            CALL geom2(thetascat,phiscat,op)
            phi = anglesum(1, phi, 0.)

           END IF

           !    Add the contribution of the scattering event to the shape factors


           IF (COS( theta_orig ) > 0.) THEN

              fdown( lam, layer - 1 ) = fdown( lam, layer -1 ) + 1

           ELSE

              fup( lam,layer )  = fup( lam,layer ) + 1

           END IF

           IF ((COS( theta ) * COS( theta_orig )) < 0) THEN 
              
              ! Photon is traveling in the opposite stream

              IF( COS( theta ) > 0 ) THEN    
                 
                 ! Photon is scattered in downward direction

                 bup( lam,layer )  = bup( lam,layer )   + 1

              ELSE                                    
                 
                 ! Photon is scattered in upward direction

                 bdown( lam,layer - 1 )= bdown( lam,layer - 1 ) + 1

              ENDIF

           END IF

           ! Adjust layer location if  photon trajectory switched from downward
           ! to upward

           IF ( (theta_orig < ( PI / 2.)) .AND. (theta >= ( PI / 2.)) ) THEN

              layer = layer - 1

           END IF


           IF ( (theta_orig >= ( PI / 2.)) .AND. (theta < ( PI / 2.)) ) THEN

              layer = layer + 1

           END IF

           ! ****************** absorb *********************

        ELSE

           action = 1
           EXIT PHOTONLIFE_LOOP

           ! *************** END of all possible outcome ***

        END IF

     END SELECT

  END DO PHOTONLIFE_LOOP

END SUBROUTINE water


! ********************************************************************
! ********************************************************************
!
!      Subroutine to adjust the photon's trajectory in an optically
!      different medium.
!
! ********************************************************************
! ********************************************************************

SUBROUTINE olayer(lam,optlayer, op)

  USE rand_global
  USE randmod
  USE const_global
  USE propagation_global
  USE physical_global

  IMPLICIT NONE

  INTEGER,INTENT(IN)     :: lam
  INTEGER,INTENT(INOUT)  :: optlayer
  REAL,INTENT(INOUT)     :: op
  REAL                   :: opl1,opl2

  ! ******** for downwelling light ************************************

  IF(z(1) > z(0))THEN

     ! Calculate the part of the OP traveled in the new optical layer

     opl1 = ABS( dptlycst( optlayer +1) - z(1) ) / COS( theta )

     ! Now calculate the new OP for the new optical layer using that
     ! layer's IOP and the fraction of the original untraveled photon

     opl2 = atten( lam , (optlayer) ) / atten(lam, (optlayer+1)) * opl1 
     op   = (op - opl1) + opl2

     ! Recalculate the end-point of the photon's OP in the main
     ! Coordinate system

     x(1) = x(0) + ( op * SIN( theta ) * COS( phi ))
     y(1) = y(0) + ( op * SIN( theta ) * SIN( phi ))
     z(1) = z(0) + ( op * COS( theta ))

     ! ******** for upwelling light ************************************

  ELSE

     !      Calculate the coordinate where the photon interacts with the 
     !      layer interface crossed

     opl1 = ABS( dptlycst( optlayer ) - z(1)) / ABS(COS( theta ))

     opl2 = atten( lam , (optlayer) ) / atten(lam, (optlayer-1)) * opl1

     op   = (op - opl1) + opl2

     x(1) = x(0) + ( op * SIN( theta ) * COS( phi ))
     y(1) = y(0) + ( op * SIN( theta ) * SIN( phi ))
     z(1) = z(0) + ( op * COS( theta ))

  END IF

END SUBROUTINE olayer

