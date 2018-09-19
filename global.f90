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
