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
