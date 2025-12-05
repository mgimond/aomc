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

