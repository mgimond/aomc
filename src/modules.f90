!     Last change:  MG   20 Mar 2003    2:44 pm

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


