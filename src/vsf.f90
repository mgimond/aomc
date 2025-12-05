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
