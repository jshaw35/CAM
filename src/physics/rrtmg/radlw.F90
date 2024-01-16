
module radlw
!----------------------------------------------------------------------- 
! 
! Purpose: Longwave radiation calculations.
!
!-----------------------------------------------------------------------
use shr_kind_mod,      only: r8 => shr_kind_r8
use ppgrid,            only: pcols, pver, pverp
use scamMod,           only: single_column, scm_crm_mode
use parrrtm,           only: nbndlw, ngptlw
use rrtmg_lw_init,     only: rrtmg_lw_ini
use rrtmg_lw_rad,      only: rrtmg_lw
use spmd_utils,        only: masterproc
use perf_mod,          only: t_startf, t_stopf
use cam_logfile,       only: iulog
use cam_abortutils,    only: endrun
use radconstants,      only: nlwbands

implicit none

private
save

! Public methods

public ::&
   radlw_init,   &! initialize constants
   rad_rrtmg_lw   ! driver for longwave radiation code
   
! Private data
integer :: ntoplw    ! top level to solve for longwave cooling

! Flag for cloud overlap method
! 0=clear, 1=random, 2=maximum/random, 3=maximum
integer, parameter :: icld = 2
                              

!===============================================================================
CONTAINS
!===============================================================================

subroutine rad_rrtmg_lw(lchnk   ,ncol      ,rrtmg_levs,r_state,       &
                        pmid    ,aer_lw_abs,cld       ,tauc_lw,       &
                        qrl     ,qrlc      ,                          &
                        flns    ,flnt      ,flnsc     ,flntc  ,flwds, &
                        flut    ,flutc     ,fnl       ,fcnl   ,fldsc, &
                        lu      ,ld        )

!-----------------------------------------------------------------------
   use cam_history,         only: outfld
   use mcica_subcol_gen_lw, only: mcica_subcol_lw
   use physconst,           only: cpair
   use rrtmg_state,         only: rrtmg_state_t
   
   use wv_saturation,    only: qsat, qsat_water ! JKS add subroutine to determine WV saturation

!------------------------------Arguments--------------------------------
!
! Input arguments
!
   integer, intent(in) :: lchnk                 ! chunk identifier
   integer, intent(in) :: ncol                  ! number of atmospheric columns
   integer, intent(in) :: rrtmg_levs            ! number of levels rad is applied

!
! Input arguments which are only passed to other routines
!
    type(rrtmg_state_t), intent(in) :: r_state

   real(r8), intent(in) :: pmid(pcols,pver)     ! Level pressure (Pascals)

   real(r8), intent(in) :: aer_lw_abs (pcols,pver,nbndlw) ! aerosol absorption optics depth (LW)

   real(r8), intent(in) :: cld(pcols,pver)      ! Cloud cover
   real(r8), intent(in) :: tauc_lw(nbndlw,pcols,pver)   ! Cloud longwave optical depth by band

!
! Output arguments
!
   real(r8), intent(out) :: qrl (pcols,pver)     ! Longwave heating rate
   real(r8), intent(out) :: qrlc(pcols,pver)     ! Clearsky longwave heating rate
   real(r8), intent(out) :: flns(pcols)          ! Surface cooling flux
   real(r8), intent(out) :: flnt(pcols)          ! Net outgoing flux
   real(r8), intent(out) :: flut(pcols)          ! Upward flux at top of model
   real(r8), intent(out) :: flnsc(pcols)         ! Clear sky surface cooing
   real(r8), intent(out) :: flntc(pcols)         ! Net clear sky outgoing flux
   real(r8), intent(out) :: flutc(pcols)         ! Upward clear-sky flux at top of model
   real(r8), intent(out) :: flwds(pcols)         ! Down longwave flux at surface
   real(r8), intent(out) :: fldsc(pcols)         ! Down longwave clear flux at surface
   real(r8), intent(out) :: fcnl(pcols,pverp)    ! clear sky net flux at interfaces
   real(r8), intent(out) :: fnl(pcols,pverp)     ! net flux at interfaces

   real(r8), pointer, dimension(:,:,:) :: lu ! longwave spectral flux up
   real(r8), pointer, dimension(:,:,:) :: ld ! longwave spectral flux down
   
!
!---------------------------Local variables-----------------------------
!
   integer :: i, k, kk, nbnd         ! indices

   real(r8) :: ful(pcols,pverp)     ! Total upwards longwave flux
   real(r8) :: fsul(pcols,pverp)    ! Clear sky upwards longwave flux
   real(r8) :: fdl(pcols,pverp)     ! Total downwards longwave flux
   real(r8) :: fsdl(pcols,pverp)    ! Clear sky downwards longwv flux

   real(r8) :: tsfc(pcols)          ! surface temperature
   real(r8) :: emis(pcols,nbndlw)   ! surface emissivity

   real(r8) :: taua_lw(pcols,rrtmg_levs-1,nbndlw)     ! aerosol optical depth by band

   real(r8), parameter :: dps = 1._r8/86400._r8 ! Inverse of seconds per day
   
   real(r8), parameter :: amdw = 1.607793_r8    ! JKS - Molecular weight of dry air / water vapor

   ! Cloud arrays for McICA 
   integer, parameter :: nsubclw = ngptlw       ! rrtmg_lw g-point (quadrature point) dimension
   integer :: permuteseed                       ! permute seed for sub-column generator

   real(r8) :: cicewp(pcols,rrtmg_levs-1)   ! in-cloud cloud ice water path
   real(r8) :: cliqwp(pcols,rrtmg_levs-1)   ! in-cloud cloud liquid water path
   real(r8) :: rei(pcols,rrtmg_levs-1)      ! ice particle effective radius (microns)
   real(r8) :: rel(pcols,rrtmg_levs-1)      ! liquid particle radius (micron)

   real(r8) :: cld_stolw(nsubclw, pcols, rrtmg_levs-1)     ! cloud fraction (mcica)
   real(r8) :: cicewp_stolw(nsubclw, pcols, rrtmg_levs-1)  ! cloud ice water path (mcica)
   real(r8) :: cliqwp_stolw(nsubclw, pcols, rrtmg_levs-1)  ! cloud liquid water path (mcica)
   real(r8) :: rei_stolw(pcols,rrtmg_levs-1)               ! ice particle size (mcica)
   real(r8) :: rel_stolw(pcols,rrtmg_levs-1)               ! liquid particle size (mcica)
   real(r8) :: tauc_stolw(nsubclw, pcols, rrtmg_levs-1)    ! cloud optical depth (mcica - optional)

   ! Includes extra layer above model top
   real(r8) :: uflx(pcols,rrtmg_levs+1)  ! Total upwards longwave flux
   real(r8) :: uflxc(pcols,rrtmg_levs+1) ! Clear sky upwards longwave flux
   real(r8) :: dflx(pcols,rrtmg_levs+1)  ! Total downwards longwave flux
   real(r8) :: dflxc(pcols,rrtmg_levs+1) ! Clear sky downwards longwv flux
   real(r8) :: hr(pcols,rrtmg_levs)      ! Longwave heating rate (K/d)
   real(r8) :: hrc(pcols,rrtmg_levs)     ! Clear sky longwave heating rate (K/d)
   real(r8) lwuflxs(nbndlw,pcols,pverp+1)  ! Longwave spectral flux up
   real(r8) lwdflxs(nbndlw,pcols,pverp+1)  ! Longwave spectral flux down
   
   ! JKS - duplicate RRTMG_LW outputs:
   real(r8) :: uflx_temp(pcols,rrtmg_levs+1)  ! Total upwards longwave flux
   real(r8) :: uflxc_temp(pcols,rrtmg_levs+1) ! Clear sky upwards longwave flux
   real(r8) :: dflx_temp(pcols,rrtmg_levs+1)  ! Total downwards longwave flux
   real(r8) :: dflxc_temp(pcols,rrtmg_levs+1) ! Clear sky downwards longwv flux
   real(r8) :: hr_temp(pcols,rrtmg_levs)      ! Longwave heating rate (K/d)
   real(r8) :: hrc_temp(pcols,rrtmg_levs)     ! Clear sky longwave heating rate (K/d)
   real(r8) lwuflxs_temp(nbndlw,pcols,pverp+1)  ! Longwave spectral flux up
   real(r8) lwdflxs_temp(nbndlw,pcols,pverp+1)  ! Longwave spectral flux down   
   
   ! JKS - For clear-sky WV repartitioning
   real(r8) :: es_temp(pcols,rrtmg_levs)         ! Saturation vapor pressure
   real(r8) :: qs_temp(pcols,rrtmg_levs)         ! Saturation specific humidity (kg/kg)
   real(r8) :: h2ovmr_sat(pcols,rrtmg_levs)      ! WV saturation volume mixing ratio
   real(r8) :: h2ovmr_cld(pcols,rrtmg_levs)      ! WV cloudy-sky volume mixing ratio on subcolumns
   real(r8) :: h2ovmr_clr(pcols,rrtmg_levs)      ! WV clear-sky volume mixing ratio on subcolumns
   
   !-----------------------------------------------------------------------

   ! mji/rrtmg

   ! Calculate cloud optical properties here if using CAM method, or if using one of the
   ! methods in RRTMG_LW, then pass in cloud physical properties and zero out cloud optical 
   ! properties here
   
   ! Zero optional cloud optical depth input array tauc_lw, 
   ! if inputting cloud physical properties into RRTMG_LW
   !          tauc_lw(:,:,:) = 0.
   ! Or, pass in CAM cloud longwave optical depth to RRTMG_LW
   ! do nbnd = 1, nbndlw
   !    tauc_lw(nbnd,:ncol,:pver) = cldtau(:ncol,:pver)
   ! end do

   ! Call mcica sub-column generator for RRTMG_LW

   ! Call sub-column generator for McICA in radiation
   call t_startf('mcica_subcol_lw')

   ! Set permute seed (must be offset between LW and SW by at least 140 to insure 
   ! effective randomization)
   permuteseed = 150

   ! These fields are no longer supplied by CAM.
   cicewp = 0.0_r8
   cliqwp = 0.0_r8
   rei = 0.0_r8
   rel = 0.0_r8

   call mcica_subcol_lw(lchnk, ncol, rrtmg_levs-1, icld, permuteseed, pmid(:, pverp-rrtmg_levs+1:pverp-1), &
      cld(:, pverp-rrtmg_levs+1:pverp-1), cicewp, cliqwp, rei, rel, tauc_lw(:, :ncol, pverp-rrtmg_levs+1:pverp-1), &
      cld_stolw, cicewp_stolw, cliqwp_stolw, rei_stolw, rel_stolw, tauc_stolw)

   call t_stopf('mcica_subcol_lw')

   
   call t_startf('rrtmg_lw')

   ! Convert incoming water amounts from specific humidity to vmr as needed;
   ! Convert other incoming molecular amounts from mmr to vmr as needed;
   ! Convert pressures from Pa to hPa;
   ! Set surface emissivity to 1.0 here, this is treated in land surface model;
   ! Set surface temperature
   ! Set aerosol optical depth to zero for now

   emis(:ncol,:nbndlw) = 1._r8
   tsfc(:ncol) = r_state%tlev(:ncol,rrtmg_levs+1)
   taua_lw(:ncol, 1:rrtmg_levs-1, :nbndlw) = aer_lw_abs(:ncol,pverp-rrtmg_levs+1:pverp-1,:nbndlw)

   if (associated(lu)) lu(1:ncol,:,:) = 0.0_r8
   if (associated(ld)) ld(1:ncol,:,:) = 0.0_r8   
   
! Original call to rrtmg_lw:
!   call rrtmg_lw(lchnk  ,ncol ,rrtmg_levs    ,icld    ,                 &
!        r_state%pmidmb  ,r_state%pintmb  ,r_state%tlay    ,r_state%tlev    ,tsfc    ,r_state%h2ovmr, &
!        r_state%o3vmr   ,r_state%co2vmr  ,r_state%ch4vmr  ,r_state%o2vmr   ,r_state%n2ovmr  ,r_state%cfc11vmr,r_state%cfc12vmr, &
!        r_state%cfc22vmr,r_state%ccl4vmr ,emis    ,&
!        cld_stolw,tauc_stolw,cicewp_stolw,cliqwp_stolw ,rei, rel, &
!        taua_lw, &
!        uflx    ,dflx    ,hr      ,uflxc   ,dflxc   ,hrc, &
!        lwuflxs, lwdflxs)

!  Identical call to rrtmg_lw but with new arguments repeated
   call rrtmg_lw(lchnk  ,ncol ,rrtmg_levs    ,icld    ,                 &
        r_state%pmidmb  ,r_state%pintmb  ,r_state%tlay    ,r_state%tlev    ,tsfc ,r_state%h2ovmr, r_state%h2ovmr, r_state%h2ovmr, &
        r_state%o3vmr   ,r_state%co2vmr  ,r_state%ch4vmr  ,r_state%o2vmr   ,r_state%n2ovmr  ,r_state%cfc11vmr,r_state%cfc12vmr, &
        r_state%cfc22vmr,r_state%ccl4vmr ,emis    ,&
        cld_stolw,tauc_stolw,cicewp_stolw,cliqwp_stolw ,rei, rel, &
        taua_lw, &
        .false.,  &
        uflx    ,dflx    ,hr      ,uflxc   ,dflxc   ,hrc, & ! JKS - outputs here and below
!        uflx_temp    ,dflx_temp    ,hr_temp      ,uflxc_temp   ,dflxc_temp   ,hrc_temp, &
        lwuflxs, lwdflxs)   

   write(iulog,*) 'Ran new rrtmg_lw subroutine with repeated old values'
   write(iulog,*) 'dflx_temp (normal):  ',dflx_temp(1,:)
   write(iulog,*) 'dflxc_temp (normal):  ',dflxc_temp(1,:)
   write(iulog,*) 'uflx_temp (normal):  ',uflx_temp(1,:)
   write(iulog,*) 'uflxc_temp (normal):  ',uflxc_temp(1,:)

   ! JKS - order of operations.
   !   1. Compute the WV saturation specific humidity over the column levels.
   !   2. Convert to a saturation vmr, equivalent to the in-cloud H2O vmr
   !   3. Compute the residual clear-sky H2O vmr
   !   4. Call rrtmg_lw with the clear-sky H2O values
       
   call t_startf('rrtmg_lw_clr')
   
   ! Repartition WV into both cloudy and clear components for each column.
   do i = 1,ncol ! Iterate over model gridcells. Note: ncol is the number of columns a chunk is actually using, pcols is maximum number
      do k = 1,rrtmg_levs ! Iterate over model levels
         kk = max(k + (pverp-rrtmg_levs)-1,1) ! Convert between the vertical coordinates of rstate (k) and pstate (kk)
         if (cld(i,kk) .lt. 0.01 .or. cld(i,kk) .gt. 0.99) then ! There is no repartitioning in a totally clear or cloudy layer.
            h2ovmr_clr(i,k) = r_state%h2ovmr(i,k)         
            h2ovmr_cld(i,k) = r_state%h2ovmr(i,k)
            cycle ! skip rest of the do loop, no need to repartition
         end if
         !   1. Compute the WV saturation specific humidity over the column levels. Using r_state levels (k)
         ! Use level midpoint rather than interfaces
         call qsat_water(r_state%tlay(i,k),r_state%pmidmb(i,k)* 1.e2_r8,es_temp(i,k),qs_temp(i,k))
         !   2. Convert to a saturation vmr, equivalent to the in-cloud H2O vmr. Using r_state levels (k).
         h2ovmr_sat(i,k) = qs_temp(i,k) / (1._r8 - qs_temp(i,k)) * amdw
         if (r_state%h2ovmr(i,k) .gt. h2ovmr_sat(i,k)) then ! If the column is super-saturated, do not repartition.
            h2ovmr_clr(i,k) = r_state%h2ovmr(i,k)         
            h2ovmr_cld(i,k) = r_state%h2ovmr(i,k)
            cycle ! skip rest of the do loop, no need to repartition
         end if         
         !   3. Compute the cloudy and clear sky WV vmr.
         ! If there is not enough WV in the column to saturate the cloudy fraction, instead use all of the remaining WV for the cloudy column.
         ! Cloudy-first version (newer)
         h2ovmr_cld(i,k) = min(h2ovmr_sat(i,k),r_state%h2ovmr(i,k) / cld(i,kk))
         ! Repartition any remaining WV into the clear column. Arbitrary lower limit here.
         h2ovmr_clr(i,k) = max((r_state%h2ovmr(i,k) - cld(i,kk) * h2ovmr_cld(i,k)) / (1 - cld(i,kk)),1.0e-8_r8) ! arbitrary lower limit of 1e-8
         
         ! Clear-first version (but this will make the clear-sky artificial moist if the column is super-saturated w.r.t. water.
         !h2ovmr_clr(i,k) = max((r_state%h2ovmr(i,k) - (h2ovmr_sat(i,k) * cld(i,kk))) / (1 - cld(i,kk)),1.0e-8_r8) ! arbitrary lower limit of 1e-8
      end do
   end do   

   if (masterproc) then
!      write(iulog,*) 'r_state%tlay(1,:) =',r_state%tlay(1,:)
!      write(iulog,*) 'r_state%pmidmb(1,:) =',r_state%pmidmb(1,:)
!      write(iulog,*) 'es_temp(1,:) =',es_temp(1,:)
!      write(iulog,*) 'qs_temp(1,:) (sp. hum.) =',qs_temp(1,:)
      write(iulog,*) 'h2ovmr_sat(1,:) =',h2ovmr_sat(1,:)
      write(iulog,*) 'h2ovmr_clr(1,:) =',h2ovmr_clr(1,:)
      write(iulog,*) 'h2ovmr_cld(1,:) =',h2ovmr_cld(1,:)
      write(iulog,*) 'r_state%h2ovmr(1,:) =',r_state%h2ovmr(1,:)
      write(iulog,*) 'r_state%h2ovmr(1,:) / h2ovmr_sat(1,:) =',r_state%h2ovmr(1,:) / h2ovmr_sat(1,:)
      write(iulog,*) 'cld(1,:) =',cld(1,:)
   endif 

   !   4. Call rrtmg_lw now passing the new clear- and cloudy-sky h2ovmrs
   ! Remove "_temp" values for clear-sky variables to make interactive
  if (masterproc) then
     call rrtmg_lw(lchnk  ,ncol ,rrtmg_levs    ,icld    ,                 &
          r_state%pmidmb  ,r_state%pintmb  ,r_state%tlay    ,r_state%tlev    ,tsfc    ,r_state%h2ovmr, h2ovmr_clr, h2ovmr_cld, &
          r_state%o3vmr   ,r_state%co2vmr  ,r_state%ch4vmr  ,r_state%o2vmr   ,r_state%n2ovmr  ,r_state%cfc11vmr,r_state%cfc12vmr, &
          r_state%cfc22vmr,r_state%ccl4vmr ,emis    ,&
          cld_stolw,tauc_stolw,cicewp_stolw,cliqwp_stolw ,rei, rel, &
          taua_lw, &
          .true.,  &
          uflx    ,dflx    ,hr      ,uflxc   ,dflxc   ,hrc, & ! JKS - outputs here and below
         !  uflx_temp    ,dflx_temp    ,hr_temp      ,uflxc_temp   ,dflxc_temp   ,hrc_temp, &
          lwuflxs, lwdflxs)
  else
     call rrtmg_lw(lchnk  ,ncol ,rrtmg_levs    ,icld    ,                 &
          r_state%pmidmb  ,r_state%pintmb  ,r_state%tlay    ,r_state%tlev    ,tsfc    ,r_state%h2ovmr, h2ovmr_clr, h2ovmr_cld, &
          r_state%o3vmr   ,r_state%co2vmr  ,r_state%ch4vmr  ,r_state%o2vmr   ,r_state%n2ovmr  ,r_state%cfc11vmr,r_state%cfc12vmr, &
          r_state%cfc22vmr,r_state%ccl4vmr ,emis    ,&
          cld_stolw,tauc_stolw,cicewp_stolw,cliqwp_stolw ,rei, rel, &
          taua_lw, &
          .false.,  &
          uflx    ,dflx    ,hr      ,uflxc   ,dflxc   ,hrc, & ! JKS - outputs here and below
         !  uflx_temp    ,dflx_temp    ,hr_temp      ,uflxc_temp   ,dflxc_temp   ,hrc_temp, &
          lwuflxs, lwdflxs)   
  end if

!   if (masterproc) then
!      write(iulog,*) 'Ran new rrtmg_lw subroutine with repartitioned values non-interactive'
!      write(iulog,*) 'dflx_temp:  ',dflx_temp(1,:)
!      write(iulog,*) 'dflxc_temp:  ',dflxc_temp(1,:)
!      write(iulog,*) 'uflx_temp:  ',uflx_temp(1,:)
!      write(iulog,*) 'uflxc_temp:  ',uflxc_temp(1,:)
!   endif 
        
   call t_stopf('rrtmg_lw_clr')

   !
   !----------------------------------------------------------------------
   ! All longitudes: store history tape quantities
   ! Flux units are in W/m2 on output from rrtmg_lw and contain output for
   ! extra layer above model top with vertical indexing from bottom to top.
   ! Heating units are in K/d on output from RRTMG and contain output for
   ! extra layer above model top with vertical indexing from bottom to top.
   ! Heating units are converted to J/kg/s below for use in CAM. 

   flwds(:ncol) = dflx (:ncol,1)
   fldsc(:ncol) = dflxc(:ncol,1)
   flns(:ncol)  = uflx (:ncol,1) - dflx (:ncol,1)
   flnsc(:ncol) = uflxc(:ncol,1) - dflxc(:ncol,1)
   flnt(:ncol)  = uflx (:ncol,rrtmg_levs) - dflx (:ncol,rrtmg_levs)
   flntc(:ncol) = uflxc(:ncol,rrtmg_levs) - dflxc(:ncol,rrtmg_levs)
   flut(:ncol)  = uflx (:ncol,rrtmg_levs)
   flutc(:ncol) = uflxc(:ncol,rrtmg_levs)

   !
   ! Reverse vertical indexing here for CAM arrays to go from top to bottom.
   !
   ful = 0._r8
   fdl = 0._r8
   fsul = 0._r8
   fsdl = 0._r8
   ful (:ncol,pverp-rrtmg_levs+1:pverp)= uflx(:ncol,rrtmg_levs:1:-1)
   fdl (:ncol,pverp-rrtmg_levs+1:pverp)= dflx(:ncol,rrtmg_levs:1:-1)
   fsul(:ncol,pverp-rrtmg_levs+1:pverp)=uflxc(:ncol,rrtmg_levs:1:-1)
   fsdl(:ncol,pverp-rrtmg_levs+1:pverp)=dflxc(:ncol,rrtmg_levs:1:-1)

   if (single_column.and.scm_crm_mode) then
      call outfld('FUL     ',ful,pcols,lchnk)
      call outfld('FDL     ',fdl,pcols,lchnk)
      call outfld('FULC    ',fsul,pcols,lchnk)
      call outfld('FDLC    ',fsdl,pcols,lchnk)
   endif
   
   fnl(:ncol,:) = ful(:ncol,:) - fdl(:ncol,:)
   ! mji/ cam excluded this?
   fcnl(:ncol,:) = fsul(:ncol,:) - fsdl(:ncol,:)

   ! Pass longwave heating to CAM arrays and convert from K/d to J/kg/s
   qrl = 0._r8
   qrlc = 0._r8
   qrl (:ncol,pverp-rrtmg_levs+1:pver)=hr (:ncol,rrtmg_levs-1:1:-1)*cpair*dps
   qrlc(:ncol,pverp-rrtmg_levs+1:pver)=hrc(:ncol,rrtmg_levs-1:1:-1)*cpair*dps

   ! Return 0 above solution domain
   if ( ntoplw > 1 )then
      qrl(:ncol,:ntoplw-1) = 0._r8
      qrlc(:ncol,:ntoplw-1) = 0._r8
   end if

   ! Pass spectral fluxes, reverse layering
   ! order=(/3,1,2/) maps the first index of lwuflxs to the third index of lu.
   if (associated(lu)) then
      lu(:ncol,pverp-rrtmg_levs+1:pverp,:) = reshape(lwuflxs(:,:ncol,rrtmg_levs:1:-1), &
           (/ncol,rrtmg_levs,nbndlw/), order=(/3,1,2/))
   end if
   
   if (associated(ld)) then
      ld(:ncol,pverp-rrtmg_levs+1:pverp,:) = reshape(lwdflxs(:,:ncol,rrtmg_levs:1:-1), &
           (/ncol,rrtmg_levs,nbndlw/), order=(/3,1,2/))
   end if
   
   call t_stopf('rrtmg_lw')

end subroutine rad_rrtmg_lw

!-------------------------------------------------------------------------------

subroutine radlw_init()
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Initialize various constants for radiation scheme.
!
!-----------------------------------------------------------------------

   use ref_pres, only : pref_mid

   integer :: k

   ! If the top model level is above ~90 km (0.1 Pa), set the top level to compute
   ! longwave cooling to about 80 km (1 Pa)
   if (pref_mid(1) .lt. 0.1_r8) then
      do k = 1, pver
         if (pref_mid(k) .lt. 1._r8) ntoplw  = k
      end do
   else
      ntoplw  = 1
   end if
   if (masterproc) then
      write(iulog,*) 'radlw_init: ntoplw =',ntoplw
   endif

   call rrtmg_lw_ini

end subroutine radlw_init

!-------------------------------------------------------------------------------

end module radlw
