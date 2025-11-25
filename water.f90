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

