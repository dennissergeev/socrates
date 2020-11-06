! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Set the droplet and ice crystal dimension in the Socrates cloud type
!
!------------------------------------------------------------------------------
module socrates_set_cld_dim
implicit none
character(len=*), parameter, private :: ModuleName = 'SOCRATES_SET_CLD_DIM'
contains

subroutine set_cld_dim(cld, control, dimen, spectrum, atm, &
  liq_nc, liq_conv_nc, &
  liq_dim, ice_dim, liq_conv_dim, ice_conv_dim, &
  liq_nc_1d, liq_conv_nc_1d, &
  liq_dim_1d, ice_dim_1d, liq_conv_dim_1d, ice_conv_dim_1d, &
  l_invert, l_debug, i_profile_debug)

use def_cld,      only: StrCld
use def_control,  only: StrCtrl
use def_dimen,    only: StrDim
use def_spectrum, only: StrSpecData
use def_atm,      only: StrAtm
use realtype_rd,  only: RealK
use rad_pcf,      only: &
  ip_cloud_split_homogen, ip_cloud_split_ice_water, &
  ip_clcmp_st_water, ip_clcmp_st_ice, ip_clcmp_cnv_water, ip_clcmp_cnv_ice, &
  ip_re_external, &
  i_normal, i_err_fatal
use ereport_mod,  only: ereport
use errormessagelength_mod, only: errormessagelength

implicit none


! Cloud properties:
type(StrCld),      intent(inout) :: cld

! Control options:
type(StrCtrl),     intent(in)  :: control

! Dimensions:
type(StrDim),      intent(in)  :: dimen

! Spectral data:
type(StrSpecData), intent(in)  :: spectrum

! Atmospheric properties:
type(StrAtm),      intent(in)  :: atm

real(RealK), intent(in), dimension(:, :), optional :: &
  liq_nc, liq_conv_nc, &
  liq_dim, ice_dim, liq_conv_dim, ice_conv_dim
real(RealK), intent(in), dimension(:), optional :: &
  liq_nc_1d, liq_conv_nc_1d, &
  liq_dim_1d, ice_dim_1d, liq_conv_dim_1d, ice_conv_dim_1d
!   Liquid number concentration, liquid and ice effective dimensions

logical, intent(in), optional :: l_invert
!   Flag to invert fields in the vertical

logical, intent(in), optional :: l_debug
integer, intent(in), optional :: i_profile_debug
!   Options for outputting debugging information


! Local variables
integer :: i, k, l
!   Loop variables
integer :: i_param_type
!   Working variable
logical :: l_inv
!   Local flag to invert fields in the vertical
logical :: l_split
!   Split cloud into optically thick and thin regions

integer                      :: ierr = i_normal
character (len=*), parameter :: RoutineName = 'SET_CLD_DIM'
character (len=errormessagelength) :: cmessage


if (.not.control%l_cloud) then
  return
end if

if (present(l_invert)) then
  l_inv = l_invert
else
  l_inv = .false.
end if


select case (control%i_cloud_representation)
case (ip_cloud_split_homogen, ip_cloud_split_ice_water)
  ! Stratiform and convective cloud is combined and the
  ! total cloud is split into optically thick and thin regions
  l_split = .true.
case default
  l_split = .false.
end select


do i=1, cld%n_condensed
  if (cld%type_condensed(i) == ip_clcmp_st_water .or. &
     (cld%type_condensed(i) == ip_clcmp_cnv_water .and. l_split)) then
    if (control%i_drop_re == ip_re_external .and. &
        (present(liq_dim).or.present(liq_dim_1d))) then
      ! Set the droplet effective radius
      call set_cld_field(cld%condensed_dim_char(:, :, i), &
                         liq_dim, liq_dim_1d)
    else
      call set_liq_dim(liq_nc, liq_nc_1d)
    end if
    i_param_type = control%i_st_water

  else if (cld%type_condensed(i) == ip_clcmp_cnv_water) then
    if (control%i_drop_re == ip_re_external .and. &
        (present(liq_conv_dim).or.present(liq_conv_dim_1d))) then
      call set_cld_field(cld%condensed_dim_char(:, :, i), &
                         liq_conv_dim, liq_conv_dim_1d)
    else
      call set_liq_dim(liq_conv_nc, liq_conv_nc_1d)
    end if
    i_param_type = control%i_cnv_water

  else if (cld%type_condensed(i) == ip_clcmp_st_ice .or. &
          (cld%type_condensed(i) == ip_clcmp_cnv_ice .and. l_split)) then
    if (present(ice_dim).or.present(ice_dim_1d)) then
      ! Set the ice-crystal effective dimension
      call set_cld_field(cld%condensed_dim_char(:, :, i), &
                         ice_dim, ice_dim_1d)
    else
      call set_ice_dim()
    end if
    i_param_type = control%i_st_ice

  else if (cld%type_condensed(i) == ip_clcmp_cnv_ice) then
    if (present(ice_conv_dim).or.present(ice_conv_dim_1d)) then
      call set_cld_field(cld%condensed_dim_char(:, :, i), &
                         ice_conv_dim, ice_conv_dim_1d)
    else
      call set_ice_dim()
    end if
    i_param_type = control%i_cnv_ice

  end if

  select case (cld%type_condensed(i))
  case (ip_clcmp_st_water, ip_clcmp_cnv_water)
    ! Constrain effective radius to be within parametrisation bounds
    do k = dimen%id_cloud_top, atm%n_layer
      do l = 1, atm%n_profile
        cld%condensed_dim_char(l, k, i) = min( max( &
          cld%condensed_dim_char(l, k, i), &
          spectrum%drop%parm_min_dim(i_param_type) ), &
          spectrum%drop%parm_max_dim(i_param_type) )
      end do
    end do
  case (ip_clcmp_st_ice, ip_clcmp_cnv_ice)
    ! Constrain effective dimension to be within parametrisation bounds
    do k = dimen%id_cloud_top, atm%n_layer
      do l = 1, atm%n_profile
        cld%condensed_dim_char(l, k, i) = min( max( &
          cld%condensed_dim_char(l, k, i), &
          spectrum%ice%parm_min_dim(i_param_type) ), &
          spectrum%ice%parm_max_dim(i_param_type) )
      end do
    end do
  end select

end do ! over condensed components



contains

  subroutine set_cld_field(out_field, full_field, oned_field)
    implicit none

    real(RealK), intent(inout) :: out_field(:, dimen%id_cloud_top:)
!     Output field
    real(RealK), intent(in), optional :: full_field(:, :)
!     Full field variable
    real(RealK), intent(in), optional :: oned_field(:)
!     One-dimensional variable

    if (present(full_field)) then
      if (l_inv) then
        do k = dimen%id_cloud_top, atm%n_layer
          do l=1, atm%n_profile
            out_field(l, k) = full_field(l, atm%n_layer+1-k)
          end do
        end do
      else
        do k = dimen%id_cloud_top, atm%n_layer
          do l=1, atm%n_profile
            out_field(l, k) = full_field(l, k)
          end do
        end do
      end if
    else if (present(oned_field)) then
      if (l_inv) then
        do k = dimen%id_cloud_top, atm%n_layer
          do l=1, atm%n_profile
            out_field(l, k) = oned_field(atm%n_layer+1-k)
          end do
        end do
      else
        do k = dimen%id_cloud_top, atm%n_layer
          do l=1, atm%n_profile
            out_field(l, k) = oned_field(k)
          end do
        end do
      end if
    else
      cmessage = 'The required cloud fields have not been provided.'
      ierr=i_err_fatal
      call ereport(ModuleName//':'//RoutineName, ierr, cmessage)      
    end if
  end subroutine set_cld_field


  subroutine set_liq_dim(full_nc, oned_nc)

    use rad_pcf, only: ip_re_liu
    use rad_ccf, only: pi, rho_water
    implicit none

    real(RealK), intent(in), optional :: full_nc(:, :)
!     Full field cloud droplet number concentration
    real(RealK), intent(in), optional :: oned_nc(:)
!     One-dimensional cloud droplet number concentration

    real(RealK) :: cdnc(dimen%nd_profile, dimen%id_cloud_top:dimen%nd_layer)
!     Working cloud droplet number concentration
    real(RealK), parameter :: eps = epsilon(1.0_RealK)

    ! Parameters for Liu spectral dispersion
    real(RealK), parameter :: aparam = 0.07_RealK
    real(RealK), parameter :: bparam = -0.14_RealK
    real(RealK) :: beta

    select case (control%i_drop_re)
    case (ip_re_liu)
      call set_cld_field(cdnc, full_nc, oned_nc)
      ! Apply Liu spectral dispersion
      do k = dimen%id_cloud_top, atm%n_layer
        do l=1, atm%n_profile
          beta = aparam * ((MAX(eps, &
            cld%condensed_mix_ratio(l, k, i)) * atm%density(l, k) &
            * 1.0e-3_RealK / (cdnc(l, k) * 1e-06_RealK))**(bparam))
          cld%condensed_dim_char(l, k, i) = MAX(0.0_RealK, &
            3.0_RealK * cld%condensed_mix_ratio(l, k, i) * atm%density(l, k) &
            / (4.0_RealK * pi * rho_water * (beta**(-3)) * cdnc(l, k))) &
            **(1.0_RealK/3.0_RealK)
        end do
      end do
    case default
      ! A default value of 7-microns is assumed.
      cld%condensed_dim_char(:, :, i) = 7.0e-6_RealK
    end select

  end subroutine set_liq_dim


  subroutine set_ice_dim()

    use rad_pcf, only: &
      ip_ice_adt, ip_ice_agg_de, ip_ice_agg_de_sun
    implicit none

    ! Parameters for the aggregate parametrization.
    real (RealK), parameter :: a0_agg_cold = 7.5094588e-04_RealK
    real (RealK), parameter :: b0_agg_cold = 5.0830326e-07_RealK
    real (RealK), parameter :: a0_agg_warm = 1.3505403e-04_RealK
    real (RealK), parameter :: b0_agg_warm = 2.6517429e-05_RealK
    real (RealK), parameter :: t_switch    = 216.208_RealK
    real (RealK), parameter :: t0_agg      = 279.5_RealK
    real (RealK), parameter :: s0_agg      = 0.05_RealK
    real (RealK), parameter :: dge2de      = &
      (3.0_RealK/2.0_RealK)*(3.0_RealK/(2.0_RealK*sqrt(3.0_RealK)))

    select case (cld%i_condensed_param(i))

    case (ip_ice_adt)
      ! This parametrization is based on the mean maximum
      ! dimension of the crystal, determined as a function of
      ! the local temperature. The size is limited to its value
      ! at the freezing level.
      do k = dimen%id_cloud_top, atm%n_layer
        do l = 1, atm%n_profile
          cld%condensed_dim_char(l, k, i) = min( 7.198755e-04_RealK, &
            exp(5.522e-02_RealK * (atm%t(l, k) - 2.7965e+02_RealK)) &
            / 9.702e+02_RealK )
        end do
      end do

    case (ip_ice_agg_de, ip_ice_agg_de_sun)
      ! Aggregate parametrization based on effective dimension.
      ! The fit provided here is based on Stephan Havemann's fit of
      ! Dge with temperature, consistent with David Mitchell's treatment
      ! of the variation of the size distribution with temperature. The
      ! parametrization of the optical properties is based on De
      ! (=(3/2)volume/projected area), whereas Stephan's fit gives Dge
      ! (=(2*SQRT(3)/3)*volume/projected area), which explains the
      ! conversion factor. The fit to Dge is in two sections, because
      ! Mitchell's relationship predicts a cusp at 216.208 K. Limits
      ! of 8 and 124 microns are imposed on Dge: these are based on this
      ! relationship and should be reviewed if it is changed. Note also
      ! that the relationship given here is for polycrystals only.
      do k = dimen%id_cloud_top, atm%n_layer
        do l = 1, atm%n_profile
          ! Preliminary calculation of Dge.
          if (atm%t(l, k) < t_switch) then
            cld%condensed_dim_char(l, k, i) = &
              a0_agg_cold*exp(s0_agg*(atm%t(l, k)-t0_agg))+b0_agg_cold
          else
            cld%condensed_dim_char(l, k, i) = &
              a0_agg_warm*exp(s0_agg*(atm%t(l, k)-t0_agg))+b0_agg_warm
          end if
          ! Limit and convert to De.
          cld%condensed_dim_char(l, k, i) = dge2de * min( 1.24e-04_RealK, &
            max(8.0e-06_RealK, cld%condensed_dim_char(l, k, i)) )
        end do
      end do

    case default
      ! A default value of 30-microns is assumed.
      cld%condensed_dim_char(:, :, i) = 30.0e-6_RealK

    end select

  end subroutine set_ice_dim

end subroutine set_cld_dim
end module socrates_set_cld_dim
