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
