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
