!     Last change:  MG   23 December 2025
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
  !  12/23/2025  Manuel Gimond        Improved interpolation scheme, added power-law
  !                                   cumulative distribution function for forward-most
  !                                   angle


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
  REAL :: fraction

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

     IF (random <= spf(1,ii) ) THEN ! The forward-most scattering direction

!        dthetascat = spf(1,0) 
!        thetascat = rand() * dthetascat  ! uniform interpolation

        ! The following adopts Mobley et al.'s  power-law extrapolation (1993)
        ! thetascat = spf(1,0) * (random / spf(1,ii))**(1.0 / (2.0 - 1.346))
        thetascat = spf(1,0) * (random / spf(1,ii))**(1.529052)
        
     ELSE

        ANGLE_LOOP: DO

           IF (random > spf(jj - 1, ii) .AND. random <= spf(jj,ii) ) THEN

              ! Discrete sampling technique (quickest, but least accurate)
              ! thetascat = spf(jj,0)  

              ! A random value between 2 discrete SPF values can be generated.
              ! Uncomment the follwoing 2 lines if this option is desired
              !
               fraction  = (random - spf(jj - 1, ii)) / (spf(jj, ii) - spf(jj - 1, ii))
               thetascat = spf(jj - 1, 0) + fraction * (spf(jj, 0) - spf(jj - 1, 0))


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
