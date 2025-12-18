!     Last change:  MG   20 Mar 2003    3:03 pm
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
  ! The following logging technique applies to a spherical grid
  ! whose delta alpha are equal over the entire sphere

  IF ( angint == 1) THEN
     t = INT( theta / nalpha ) + limit                 

     ! The following logging technique applies to a spherical grid
     ! whose delta cosine are equal over the entire sphere (hence, differing
     ! delta alpha)

  ELSE

     IF ( COS(theta) > (1 - mum) ) THEN  ! Downward polar cap

        t = 1

     ELSE IF (COS(theta) < (mum -1 )) THEN ! Upward polar cap

        t = alphaint

     ELSE IF (COS(theta) >= 0 ) THEN

        DO t = 2 , alphaint/2

           IF ( COS(theta) < ( (1-mum) - (t-2)*muu) .AND. COS(theta) > ( (1-mum) - (t-1)*muu) ) THEN

              EXIT

           END IF

        END DO

     ELSE

        DO t =  (alphaint/2 +1) , (alphaint -1) 

           IF( COS(theta) < (-(t - REAL(alphaint)/2. -1) * muu) .AND. COS(theta) > (-(t - REAL(alphaint)/2.) * muu)) THEN

              EXIT

           END IF

        END DO

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
