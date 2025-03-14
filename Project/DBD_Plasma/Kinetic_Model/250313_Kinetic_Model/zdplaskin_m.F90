!
! ZDPLASKIN version 2.0a
! (c) 2008, Sergey Pancheshnyi (pancheshnyi@gmail.com)
!
! BOLSIG+
! (c) 2005, Gerjan Hagelaar (gerjan.hagelaar@laplace.univ-tlse.fr)
!
! http://www.zdplaskin.laplace.univ-tlse.fr/
! This software is provided "as is" without warranty and non-commercial use is freely
! granted provided proper reference is made in publications resulting from its use.
! Use of ZDPlasKin in commerical software requires a license.
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
! Fri Mar 14 10:20:01 2025
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
! ELEMENTS: E C H 
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
! MODULE ZDPlasKin
!
!-----------------------------------------------------------------------------------------------------------------------------------
module ZDPlasKin
  use dvode_f90_m, only : vode_opts
  implicit none
  public
!
! config
!
  integer, parameter :: species_max = 39, species_electrons = 1, species_length = 9, reactions_max = 114, reactions_length = 30
  double precision                          :: density(species_max)
  integer                                   :: species_charge(species_max)
  character(species_length)                 :: species_name(species_max)
  character(reactions_length)               :: reaction_sign(reactions_max)
  logical                                   :: lreaction_block(reactions_max)
!
! internal config
!
  double precision, parameter, private      :: vode_atol = 1.00D-10, vode_rtol = 1.00D-05, cfg_rtol = 1.00D-01
  integer, parameter, private               :: vode_neq = species_max + 1
  type (vode_opts), private                 :: vode_options
  integer, private                          :: vode_itask, vode_istate, ifile_unit = 5
  double precision, private                 :: stat_dens(species_max), stat_src(species_max), stat_rrt(reactions_max), stat_time, &
                                               dens_loc(vode_neq,0:3), rrt_loc(reactions_max), tsav = -huge(tsav), &
                                               mach_accur, mach_tiny
  double precision                          :: rrt(reactions_max), mrtm(species_max, reactions_max), ZDPlasKin_cfg(14)
  logical, private                          :: lZDPlasKin_init = .false., lprint, lstat_accum, ldensity_constant, &
                                               density_constant(species_max), lgas_heating
!
! qtplaskin config
!
  logical, private                          :: lqtplaskin, lqtplaskin_first = .true.
  double precision, parameter, private      :: qtplaskin_atol = 1.00D+00, qtplaskin_rtol = 1.00D-02
  character(32), allocatable                :: qtplaskin_user_names(:)
  double precision, allocatable             :: qtplaskin_user_data(:)
!
! physical constants
!
  double precision, parameter, private      :: eV_to_K = 1.16045052d4, q_elem = 1.60217662d-19, k_B = 1.38064852d-23
!
! bolsig+ config
!
  double precision, parameter, private      :: bolsig_rtol = 1.00D-03, bolsig_rtol_half = 3.16D-02, &
                                               bolsig_field_min = 1.00D-01, bolsig_field_max = 1.00D+03, &
                                               bolsig_eecol_frac_def = 1.00D-05
  double precision, private                 :: bolsig_eecol_frac
  integer, parameter, private               :: bolsig_species_max = 28, bolsig_species_length = 9, bolsig_rates_max = 74 
  character(*), parameter, private          :: bolsigfile = "bolsigdb.dat"
  integer                                   :: bolsig_pointer(bolsig_rates_max) = -1
  integer, private                          :: bolsig_species_index(bolsig_species_max) = -1, bolsig_collisions_max = 0 
  logical, private                          :: lbolsig_ignore_gas_temp, lbolsig_Maxwell_EEDF
  double precision, allocatable             :: bolsig_rates(:)
  character(bolsig_species_length), private :: bolsig_species(bolsig_species_max)
  interface
    subroutine ZDPlasKin_bolsig_Init(a)
      character(*), intent(in) :: a
    end subroutine ZDPlasKin_bolsig_Init
    subroutine ZDPlasKin_bolsig_ReadCollisions(a)
      character(*), intent(in) :: a
    end subroutine ZDPlasKin_bolsig_ReadCollisions
    subroutine ZDPlasKin_bolsig_GetCollisions(i,j)
      integer, intent(out) :: i, j
    end subroutine ZDPlasKin_bolsig_GetCollisions
    subroutine ZDPlasKin_bolsig_GetSpeciesName(a,i)
      integer, intent(in) :: i
      character(*), intent(out) :: a
    end subroutine ZDPlasKin_bolsig_GetSpeciesName
    subroutine ZDPlasKin_bolsig_GetReactionName(a,i)
      integer, intent(in) :: i
      character(*), intent(out) :: a
    end subroutine ZDPlasKin_bolsig_GetReactionName
    subroutine ZDPlasKin_bolsig_SolveBoltzmann(i,a,j,b)
      integer, intent(in) :: i, j
      double precision, intent(in)  :: a(i)
      double precision, intent(out) :: b(j)
    end subroutine ZDPlasKin_bolsig_SolveBoltzmann
    subroutine ZDPlasKin_bolsig_GetEEDF(i,a,b)
      integer, intent(in) :: i
      double precision, intent(out) :: a,b
    end subroutine ZDPlasKin_bolsig_GetEEDF
  end interface
!
! data section
!
  data species_charge(1:species_max) &
  /-1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0/
  data species_name(1:species_max) &
  /"E        ","C        ","CH       ","CH2      ","CH3      ","CH3^+    ","CH4      ","CH4(V13) ","CH4(V24) ","CH4^+    ",&
   "CH5^+    ","C2H2     ","C2H2(V13)","C2H2(V2) ","C2H2(V5) ","C2H2^+   ","C2H3     ","C2H4     ","C2H4(V1) ","C2H4(V2) ",&
   "C2H4^+   ","C2H5     ","C2H5^+   ","C2H6     ","C2H6(V13)","C2H6(V24)","C2H6^+   ","C3H6     ","C3H6(V)  ","C3H6^+   ",&
   "C3H7     ","C3H8     ","C3H8(V1) ","C3H8(V2) ","C3H8^+   ","C4H9H    ","C5H12    ","H        ","H2       "/
  data reaction_sign(1:54) &
  /"bolsig:CH4->CH4(V24)          ","bolsig:CH4->CH4(V13)          ","bolsig:CH4(V24)->CH4          ",&
   "bolsig:CH4(V13)->CH4          ","bolsig:C2H6->C2H6(V24)        ","bolsig:C2H6->C2H6(V13)        ",&
   "bolsig:C2H6(V24)->C2H6        ","bolsig:C2H6(V13)->C2H6        ","bolsig:C2H4->C2H4(V2)         ",&
   "bolsig:C2H4->C2H4(V1)         ","bolsig:C2H4(V2)->C2H4         ","bolsig:C2H4(V1)->C2H4         ",&
   "bolsig:C2H2->C2H2(V2)         ","bolsig:C2H2->C2H2(V5)         ","bolsig:C2H2->C2H2(V13)        ",&
   "bolsig:C2H2(V2)->C2H2         ","bolsig:C2H2(V5)->C2H2         ","bolsig:C2H2(V13)->C2H2        ",&
   "bolsig:C3H8->C3H8(V2)         ","bolsig:C3H8->C3H8(V1)         ","bolsig:C3H8(V2)->C3H8         ",&
   "bolsig:C3H8(V1)->C3H8         ","bolsig:C3H6->C3H6(V)          ","bolsig:C3H6(V)->C3H6          ",&
   "bolsig:CH3->CH3^+             ","bolsig:CH4->CH4^+             ","bolsig:CH4(V24)->CH4^+        ",&
   "bolsig:CH4(V13)->CH4^+        ","bolsig:C2H6->C2H6^+           ","bolsig:C2H6(V24)->C2H6^+      ",&
   "bolsig:C2H6(V13)->C2H6^+      ","bolsig:C2H4->C2H4^+           ","bolsig:C2H4(V2)->C2H4^+       ",&
   "bolsig:C2H4(V1)->C2H4^+       ","bolsig:C2H2->C2H2^+           ","bolsig:C2H2(V2)->C2H2^+       ",&
   "bolsig:C2H2(V5)->C2H2^+       ","bolsig:C2H2(V13)->C2H2^+      ","bolsig:C3H8->C3H8^+           ",&
   "bolsig:C3H8(V2)->C3H8^+       ","bolsig:C3H8(V1)->C3H8^+       ","bolsig:C3H6->C3H6^+           ",&
   "bolsig:C3H6(V)->C3H6^+        ","bolsig:CH3^+->CH3             ","bolsig:CH4^+->CH4             ",&
   "bolsig:C2H6^+->C2H6           ","bolsig:C2H4^+->C2H4           ","bolsig:C2H2^+->C2H2           ",&
   "bolsig:C3H8^+->C3H8           ","bolsig:C3H6^+->C3H6           ","bolsig:CH4->CH3H              ",&
   "bolsig:CH4(V24)->CH3H         ","bolsig:CH4(V13)->CH3H         ","bolsig:CH4->CH2H2             "/
  data reaction_sign(55:108) &
  /"bolsig:CH4(V24)->CH2H2        ","bolsig:CH4(V13)->CH2H2        ","bolsig:C2H6->C2H4H2           ",&
   "bolsig:C2H6(V24)->C2H4H2      ","bolsig:C2H6(V13)->C2H4H2      ","E+C2H5^+=>C2H3+H+H            ",&
   "bolsig:C2H4->C2H2H2           ","bolsig:C2H4(V1)->C2H2H2       ","bolsig:C2H4(V2)->C2H2H2       ",&
   "E+C2H5^+=>C2H2+H2+H           ","bolsig:C3H6->C2H2CH4          ","bolsig:C3H6(V)->C2H2CH4       ",&
   "E+C2H5^+=>C2H2+H+H+H          ","bolsig:C2H4->C2H3H            ","bolsig:C2H4(V1)->C2H3H        ",&
   "bolsig:C2H4(V2)->C2H3H        ","bolsig:C3H8->C3H6H2           ","bolsig:C3H8(V1)->C3H6H2       ",&
   "bolsig:C3H8(V2)->C3H6H2       ","bolsig:CH4->CHH2H             ","bolsig:CH4(V24)->CHH2H        ",&
   "bolsig:CH4(V13)->CHH2H        ","bolsig:CH->CH                 ","CH3+H=>CH4                    ",&
   "CH3+CH3=>C2H6                 ","C2H4+H=>C2H5                  ","C2H4(V1)+H=>C2H5              ",&
   "C2H4(V2)+H=>C2H5              ","CH+CH4=>C2H4+H                ","CH+CH4(V24)=>C2H4+H           ",&
   "CH+CH4(V13)=>C2H4+H           ","C2H5+H=>C2H6                  ","C2H5+C2H5=>C4H9H              ",&
   "CH2+H=>CH+H2                  ","CH2+CH3=>C2H4+H               ","C2H5+H=>CH3+CH3               ",&
   "CH4+CH4^+=>CH3+CH5^+          ","CH4(V24)+CH4^+=>CH3+CH5^+     ","CH4(V13)+CH4^+=>CH3+CH5^+     ",&
   "CH3+C2H5=>C3H8                ","CH3+C2H3=>C3H6                ","C3H6+H=>C3H7                  ",&
   "C3H6(V)+H=>C3H7               ","C4H9H+CH2=>C5H12              ","C3H7+H2=>C3H8+H               ",&
   "C2H6+CH5^+=>CH4+H2+C2H5^+     ","C2H6(V24)+CH5^+=>CH4+H2+C2H5^+","C2H6(V13)+CH5^+=>CH4+H2+C2H5^+",&
   "C2H4+CH5^+=>CH4+C2H5^+        ","C2H4(V1)+CH5^+=>CH4+C2H5^+    ","C2H4(V2)+CH5^+=>CH4+C2H5^+    ",&
   "C3H7+H=>C3H8                  ","C2H3+H=>C2H2+H2               ","C3H8+CH2=>C4H9H               "/
  data reaction_sign(109:114) &
  /"C3H8(V1)+CH2=>C4H9H           ","C3H8(V2)+CH2=>C4H9H           ","CH4+CH3^+=>H2+C2H5^+          ",&
   "CH4(V24)+CH3^+=>H2+C2H5^+     ","CH4(V13)+CH3^+=>H2+C2H5^+     ","H+H=>H2                       "/
  data bolsig_species(1:bolsig_species_max) &
  /"CH       ","CH3      ","CH3^+    ","CH4      ","CH4(V13) ","CH4(V24) ","CH4^+    ","C2H2     ","C2H2(V13)","C2H2(V2) ",&
   "C2H2(V5) ","C2H2^+   ","C2H4     ","C2H4(V1) ","C2H4(V2) ","C2H4^+   ","C2H5^+   ","C2H6     ","C2H6(V13)","C2H6(V24)",&
   "C2H6^+   ","C3H6     ","C3H6(V)  ","C3H6^+   ","C3H8     ","C3H8(V1) ","C3H8(V2) ","C3H8^+   "/
contains
!-----------------------------------------------------------------------------------------------------------------------------------
!
! initialization
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_init()
  implicit none
  character(256) :: string
  integer :: i, j, k
  write(*,"(/,A)") "ZDPlasKin (version " // "2.0a" // ") INIT:"
  if( lZDPlasKin_init ) call ZDPlasKin_stop("   ERROR: the ZDPlasKin library has been initialized")
  write(string,*) species_max
  write(*,"(2x,A)")  "species        ... " // trim(adjustl(string))
  write(string,*) reactions_max
  write(*,"(2x,A)")  "reactions      ... " // trim(adjustl(string))
  if(species_max<=0 .or. reactions_max<=0) call ZDPlasKin_stop("   ERROR: wrong preconfig data")
  write(*,"(2x,A,$)") "BOLSIG+ loader ... " // trim(adjustl(bolsigfile)) // " : "
  call ZDPlasKin_bolsig_Init(bolsigfile)
  do i = 1, bolsig_species_max
    call ZDPlasKin_bolsig_ReadCollisions(trim(bolsig_species(i)))
    j = bolsig_collisions_max
    call ZDPlasKin_bolsig_GetCollisions(k,bolsig_collisions_max)
    if(bolsig_collisions_max <= j) then
      write(*,*)
      call ZDPlasKin_stop("ERROR: wrong file or missing data " // &
                        "(" // trim(adjustl(bolsigfile)) // ": <" // trim(bolsig_species(i)) // ">).")
    endif
  enddo
  if(bolsig_species_max /= k) then
    write(*,*)
    call ZDPlasKin_stop("ERROR: internal error in BOLSIG+ loader")
  endif
  write(string,*) bolsig_species_max
  write(*,"(A,$)") trim(adjustl(string)) // " species & "
  write(string,*) bolsig_collisions_max
  write(*,"(A)")   trim(adjustl(string)) // " collisions"
  write(*,"(2x,A,$)") "species  link  ... "
  j = 0
  do i = 1, bolsig_species_max
    j = j + 1
    k = 1
    do while(k<=species_max .and. bolsig_species_index(i)<=0)
      call ZDPlasKin_bolsig_GetSpeciesName(string,i)
      if(trim(species_name(k)) == trim(string)) then
        bolsig_species_index(i) = k
      else
        k = k + 1
      endif
    enddo
    if(bolsig_species_index(i) <= 0) call ZDPlasKin_stop("cannot find species link for <" // trim(string) // ">")
  enddo
  write(string,*) j
  write(*,"(A)") trim(adjustl(string))
  write(*,"(2x,A,$)") "process  link  ... "
  i = 1
  j = 1
  do while(i<=reactions_max .and. j<=bolsig_rates_max)
    if(reaction_sign(i)(1:7) == "bolsig:") then
      k = 1
      do while(k<=bolsig_collisions_max .and. bolsig_pointer(j)<=0)
        call ZDPlasKin_bolsig_GetReactionName(string,k)
        if(trim(string) == trim(reaction_sign(i)(8:))) then
          bolsig_pointer(j) = k
        else
          k = k + 1
        endif
      enddo
      if(bolsig_pointer(j) <= 0) call ZDPlasKin_stop("cannot find processes link for <" // trim(reaction_sign(i)) // ">")
      j = j + 1
    endif
    i = i + 1
  enddo
  if(j <= bolsig_rates_max) then
    call ZDPlasKin_stop("internal error")
  else
    write(string,*) bolsig_rates_max
    write(*,"(A)") trim(adjustl(string))
  endif
  i = 0
  do while((1.0d0+10.0d0**(i-1)) /= 1.0d0)
    i = i - 1
  enddo
  mach_accur = 10.0d0**i
  mach_tiny  = sqrt( tiny(mach_tiny) )
  lZDPlasKin_init = .true.
  call ZDPlasKin_reset()
  write(*,"(A,/)") "ZDPlasKin INIT DONE"
  return
end subroutine ZDPlasKin_init
!-----------------------------------------------------------------------------------------------------------------------------------
!
! timestep integration using implicit solver dvode_f90
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_timestep(time,dtime)
  use dvode_f90_m, only : dvode_f90
  implicit none
  double precision, intent(in)    ::  time
  double precision, intent(inout) :: dtime
  double precision, save :: densav(vode_neq) = 0.0d0, cfgsav(3) = 0.0d0
  double precision :: tout
  if(.not. lZDPlasKin_init) call ZDPlasKin_init()
  if( lqtplaskin .and. lqtplaskin_first ) call ZDPlasKin_write_qtplaskin(time)
  if(time < tsav) vode_istate = 1
  tsav = time
  if(dtime > 0.0d0) then
    vode_itask = 1
    tout = time + dtime
    if(dtime < mach_accur*abs(tout)) &
      call ZDPlasKin_stop("ZDPlasKin ERROR: dtime parameter is too small (subroutine ZDPlasKin_timestep)")
  else
    vode_itask = 2
    tout = ( 1.0d0 + mach_accur ) * time + mach_tiny
  endif
  dens_loc(1:species_max,0) = density(:)
  dens_loc(1:species_max,1) = 0.5d0 * ( density(:) + abs( density(:) ) )
  if(any(dens_loc(1:species_max,1) /= densav(1:species_max))) vode_istate = 1
  densav(1:species_max) = dens_loc(1:species_max,1)
  if(vode_istate /= 1 .and. any( abs(cfgsav(:)-ZDPlasKin_cfg(1:3)) > cfg_rtol*abs(cfgsav(:)+ZDPlasKin_cfg(1:3)) )) vode_istate = 1
  cfgsav(:) = ZDPlasKin_cfg(1:3)
  if( lgas_heating ) then
    if(ZDPlasKin_cfg(1) /= densav(species_max+1)) vode_istate = 1
    densav(species_max+1) = ZDPlasKin_cfg(1)
  endif
  call dvode_f90(ZDPlasKin_fex,vode_neq,densav,tsav,tout,vode_itask,vode_istate,vode_options,j_fcn=ZDPlasKin_jex)
  if(vode_istate < 0) then
    write(*,"(A,1pd11.4)") "Tgas   =", ZDPlasKin_cfg(1)
    write(*,"(A,1pd11.4)") "    Te =", ZDPlasKin_cfg(4)
    call ZDPlasKin_stop("ZDPlasKin ERROR: DVODE solver issued an error (subroutine ZDPlasKin_timestep)")
  endif
  if( lgas_heating ) ZDPlasKin_cfg(1) = densav(species_max+1)
  density(:) = dens_loc(1:species_max,0) - dens_loc(1:species_max,1) + densav(1:species_max)
  if(dtime <= 0.0d0) dtime = tsav - time
  if( lstat_accum ) then
    call ZDPlasKin_get_rates(SOURCE_TERMS=dens_loc(1:species_max,0),REACTION_RATES=rrt_loc)
    stat_dens(:) = stat_dens(:) + dtime * density(:)
    stat_src(:)  = stat_src(:)  + dtime * dens_loc(1:species_max,0)
    stat_rrt(:)  = stat_rrt(:)  + dtime * rrt_loc(:)
    stat_time    = stat_time    + dtime
  endif
  if( lqtplaskin ) call ZDPlasKin_write_qtplaskin(time+dtime)
  return
end subroutine ZDPlasKin_timestep
!-----------------------------------------------------------------------------------------------------------------------------------
!
! timestep integration using explicit Euler method
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_timestep_explicit(time,dtime,rtol_loc,atol_loc,switch_implicit)
  implicit none
  double precision, intent(in) ::  time, rtol_loc, atol_loc
  double precision, intent(inout) :: dtime
  double precision, optional, intent(in) :: switch_implicit
  double precision :: time_loc, time_end, dtime_loc, dtime_max
  logical, save :: lwarn = .true.
  if(.not. lZDPlasKin_init) call ZDPlasKin_init()
  if( lqtplaskin .and. lqtplaskin_first ) call ZDPlasKin_write_qtplaskin(time)
  if(rtol_loc <= 0.0d0) call ZDPlasKin_stop("ZDPlasKin ERROR: rtol_loc must be positive (subroutine ZDPlasKin_timestep_explicit)")
  if(atol_loc <= 0.0d0) call ZDPlasKin_stop("ZDPlasKin ERROR: atol_loc must be positive (subroutine ZDPlasKin_timestep_explicit)")
  tsav     = time
  time_loc = 0.0d0
  time_end = 0.5d0 * ( dtime + abs(dtime) ) + mach_tiny
  do while(time_loc < time_end)
    dens_loc(1:species_max,0) = density(:)
    if( lgas_heating ) dens_loc(species_max+1,0) = ZDPlasKin_cfg(1)
    dens_loc(:,1) = 0.5d0 * ( dens_loc(:,0) + abs( dens_loc(:,0) ) )
    call ZDPlasKin_fex(vode_neq,tsav,dens_loc(:,1),dens_loc(:,2))
    where(dens_loc(:,2) >= 0.0d0)
      dens_loc(:,3) = + dens_loc(:,2) /    ( rtol_loc * dens_loc(:,1) + atol_loc ) 
    elsewhere
      dens_loc(:,3) = - dens_loc(:,2) / min( rtol_loc * dens_loc(:,1) + atol_loc , dens_loc(:,1) + mach_tiny )
    endwhere
    dtime_loc = 1.0d0 / ( maxval( dens_loc(:,3) ) + mach_tiny )
    if(dtime > 0.0d0) then
      dtime_max = dtime - time_loc
      dtime_loc = min( dtime_loc , dtime_max )
      if( present(switch_implicit) ) then
        if(dtime_loc*switch_implicit < dtime_max) then
          if(lprint .and. lwarn) then
            write(*,"(A,/,A,1pd9.2,A)") "ZDPlasKin INFO: low efficiency of Euler method (subroutine ZDPlasKin_timestep_explicit)", &
                        "                ZDPlasKin_timestep subroutine will be used in similar conditions (", switch_implicit, ")"
            lwarn = .false.
          endif
          time_loc = tsav
          density(:) = dens_loc(1:species_max,0)
          if( lgas_heating ) ZDPlasKin_cfg(1) = dens_loc(species_max+1,0)
          call ZDPlasKin_timestep(time_loc,dtime_max)
          return
        endif
      endif
    else
      dtime = dtime_loc
    endif
    time_loc = time_loc + dtime_loc
    tsav     = time     +  time_loc
    density(:) = dens_loc(1:species_max,0) + dtime_loc * dens_loc(1:species_max,2)
    if( lgas_heating ) ZDPlasKin_cfg(1) = dens_loc(species_max+1,0) + dtime_loc * dens_loc(species_max+1,2)
  enddo
  if( lstat_accum ) then
    call ZDPlasKin_get_rates(SOURCE_TERMS=dens_loc(1:species_max,0),REACTION_RATES=rrt_loc)
    stat_dens(:) = stat_dens(:) + dtime * density(:)
    stat_src(:)  = stat_src(:)  + dtime * dens_loc(1:species_max,0)
    stat_rrt(:)  = stat_rrt(:)  + dtime * rrt_loc(:)
    stat_time    = stat_time    + dtime
  endif
  if( lqtplaskin ) call ZDPlasKin_write_qtplaskin(time+dtime)
  return
end subroutine ZDPlasKin_timestep_explicit
!-----------------------------------------------------------------------------------------------------------------------------------
!
! update BOLSIG+ solution and get electron parameters
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_bolsig_rates(lbolsig_force)
  implicit none
  logical, optional, intent(in) :: lbolsig_force
  logical :: lforce
  integer :: i, j, k
  integer, save :: bolsig_points_max
  logical, save :: lfirst = .true., leecol = .true.
  double precision :: error, density_loc, cfg_loc(6+bolsig_species_max)
  double precision, save :: low_density_limit = bolsig_rtol, bolsig_mesh_a, bolsig_mesh_b
  double precision, save, allocatable :: bolsig_cfg(:,:), bolsig_reslt(:,:)
  if( lfirst ) then
    if(.not. lZDPlasKin_init) call ZDPlasKin_init()
    bolsig_mesh_a = 1.0d0 / log( 1.0d0 + bolsig_rtol )
    bolsig_mesh_b = bolsig_mesh_a * log( bolsig_field_min ) - 0.5d0
    bolsig_points_max = int( bolsig_mesh_a * log( bolsig_field_max ) - bolsig_mesh_b )
    allocate(bolsig_rates(bolsig_collisions_max), &
             bolsig_cfg(6+bolsig_species_max,0:bolsig_points_max), &
             bolsig_reslt(10+bolsig_collisions_max,0:bolsig_points_max),stat=i)
    if(i /= 0) call ZDPlasKin_stop("ZDPlasKin ERROR: memory allocation error (subroutine ZDPlasKin_bolsig_rates)")
    bolsig_cfg(:,:) = 0.0d0
    lfirst = .false.
  endif
  if( present(lbolsig_force) ) then
    lforce = lbolsig_force
  else
    lforce = .false.
  endif
  if(ZDPlasKin_cfg(1) <= 0.0d0) &
    call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined GAS_TEMPERATURE (subroutine ZDPlasKin_bolsig_rates)")
  if(.not. lbolsig_Maxwell_EEDF ) then
    if(ZDPlasKin_cfg(2) < 0.0d0) &
      call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined REDUCED_FREQUENCY (subroutine ZDPlasKin_bolsig_rates)")
    if(ZDPlasKin_cfg(3) < 0.0d0) &
      call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined REDUCED_FIELD (subroutine ZDPlasKin_bolsig_rates)")
    ZDPlasKin_cfg(4) = 0.0d0
  else
    if(ZDPlasKin_cfg(4) <= 0.0d0) then
      call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined ELECTRON_TEMPERATURE (subroutine ZDPlasKin_bolsig_rates)")
    elseif(lprint .and. ZDPlasKin_cfg(4) < ZDPlasKin_cfg(1)) then
      write(*,"(A)") "ZDPlasKin INFO: ELECTRON_TEMPERATURE is below GAS_TEMPERATURE (subroutine ZDPlasKin_bolsig_rates)"
    endif
    ZDPlasKin_cfg(2:3) = 0.0d0
  endif
  density_loc = 0.5d0 * ( sum(density(bolsig_species_index(:))) + sum(abs(density(bolsig_species_index(:)))) )
  if(density_loc <= mach_tiny) then
    call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined densities configured for BOLSIG+ solver " // &
                                                                     "(subroutine ZDPlasKin_bolsig_rates)")
  elseif( lprint ) then
    call ZDPlasKin_get_density_total(ALL_NEUTRAL=error)
    error = abs( 1.0d0 - density_loc / error )
    if(error > low_density_limit) then
      write(*,"(A,1pd9.2)") "ZDPlasKin INFO: the density of species not configured for BOLSIG+ solver exceeds", error
      low_density_limit = sqrt( low_density_limit )
    endif
  endif
  cfg_loc(1:4) = ZDPlasKin_cfg(1:4)
  cfg_loc(2)   = cfg_loc(2) * 1.0d-6
  cfg_loc(6)   = 0.5d0 * ( density(species_electrons) + abs(density(species_electrons)) )
  cfg_loc(5)   = cfg_loc(6) * 1.0d6
  cfg_loc(6)   = cfg_loc(6) / density_loc
  if(cfg_loc(6) < bolsig_eecol_frac) then
    cfg_loc(6) = 0.0d0
  elseif(lprint .and. leecol) then
    write(*,"(A)") "ZDPlasKin INFO: set electron-electron collisions ON ..."
    leecol = .false.
  endif
  cfg_loc(7:) = 0.5d0 * ( density(bolsig_species_index(:)) + abs(density(bolsig_species_index(:))) ) / density_loc
  if(lbolsig_Maxwell_EEDF .or. bolsig_points_max==0) then
    lforce = .true.
    i = 0
  else
    error = min( max( cfg_loc(3) , bolsig_field_min ) , bolsig_field_max )
    i = int( bolsig_mesh_a * log( error ) - bolsig_mesh_b )
    i = max(0,min(i,bolsig_points_max))
  endif
  if( lforce ) then
    error = 2.0d0
  else
    if(lbolsig_ignore_gas_temp .and. bolsig_cfg(1,i)>0.0d0) then
      cfg_loc(1) = bolsig_cfg(1,i)
      error = 0.0d0
    else
      error = abs( ( cfg_loc(1) - bolsig_cfg(1,i) ) / ( 0.5d0 * ( cfg_loc(1) + bolsig_cfg(1,i) ) + mach_tiny) ) / bolsig_rtol_half
    endif
    if(error <= 1.0d0) then
      error = abs( ( cfg_loc(2) - bolsig_cfg(2,i) ) / ( 0.5d0 * ( cfg_loc(2) + bolsig_cfg(2,i) ) + mach_tiny) ) / bolsig_rtol_half
      if(error <= 1.0d0) then
        error = abs( ( cfg_loc(3) - bolsig_cfg(3,i) ) / ( 0.5d0 * ( cfg_loc(3) + bolsig_cfg(3,i) ) + mach_tiny) ) / bolsig_rtol
        if(error <= 1.0d0) then
          error = abs( ( max(cfg_loc(6),bolsig_eecol_frac) - max(bolsig_cfg(6,i),bolsig_eecol_frac) ) &
           / ( 0.5d0 * ( max(cfg_loc(6),bolsig_eecol_frac) + max(bolsig_cfg(6,i),bolsig_eecol_frac) ) + mach_tiny) ) &
           / bolsig_rtol_half
          if(error <= 1.0d0) error = maxval( abs( cfg_loc(7:) - bolsig_cfg(7:,i) ) ) &
                                     / ( 0.5d0 * maxval( cfg_loc(7:) + bolsig_cfg(7:,i) ) ) / bolsig_rtol
        endif
      endif
    endif
  endif
  if(error > 1.0d0) then
    j = 6 + bolsig_species_max
    k = 10 + bolsig_collisions_max
    bolsig_cfg(:,i) = cfg_loc(:)
    call ZDPlasKin_bolsig_SolveBoltzmann(j,bolsig_cfg(1:j,i),k,bolsig_reslt(1:k,i))
    if(.not. lbolsig_Maxwell_EEDF) then
      bolsig_reslt(2,i) = bolsig_reslt(2, i) * eV_to_K / 1.5d0
    else
      bolsig_reslt(2,i) = cfg_loc(4)
    endif
    bolsig_reslt(3, i) = bolsig_reslt(3, i) * 1.0d-2
    bolsig_reslt(4, i) = bolsig_reslt(4, i) * 1.0d-2 / density_loc
    bolsig_reslt(5, i) = bolsig_reslt(5, i) * 1.0d-2
    bolsig_reslt(6, i) = bolsig_reslt(6, i) * 1.0d-2
    bolsig_reslt(7:,i) = bolsig_reslt(7:,i) * 1.0d6
  endif
  ZDPlasKin_cfg(3:12) = bolsig_reslt(:10,i)
  bolsig_rates(:)     = bolsig_reslt(11:,i)
  return
end subroutine ZDPlasKin_bolsig_rates
!-----------------------------------------------------------------------------------------------------------------------------------
!
! get index of species
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_get_species_index(str,i)
  implicit none
  character(*), intent(in ) :: str
  integer,      intent(out) :: i
  character(species_length) :: string
  integer :: j, istr
  if(.not. lZDPlasKin_init) call ZDPlasKin_init()
  string = trim(adjustl(str))
  istr   = len_trim(string)
  do i = 1, istr
    j = iachar(string(i:i))
    if(j>=97 .and. j<=122) string(i:i) = achar(j-32)
  enddo
  i = 0
  j = 0
  do while(i==0 .and. j<species_max)
    j = j + 1
    if(string(1:istr) == trim(species_name(j))) i = j
  enddo
  if(i <= 0) &
    call ZDPlasKin_stop("ZDPlasKin ERROR: cannot identify species <"//trim(str)//"> (subroutine ZDPlasKin_get_species_index)")
  return
end subroutine ZDPlasKin_get_species_index
!-----------------------------------------------------------------------------------------------------------------------------------
!
! set/get density for species
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_set_density(string,DENS,LDENS_CONST)
  implicit none
  character(*), intent(in) :: string
  logical, optional, intent(in) :: LDENS_CONST
  double precision, optional, intent(in) :: DENS
  integer :: i
  call ZDPlasKin_get_species_index(string,i)
  if( present(DENS) ) density(i) = DENS
  if( present(LDENS_CONST) ) then
    density_constant(i) = LDENS_CONST
    ldensity_constant   = any( density_constant(:) )
  endif
  return
end subroutine ZDPlasKin_set_density
subroutine ZDPlasKin_get_density(string,DENS,LDENS_CONST)
  implicit none
  character(*), intent(in ) :: string
  logical, optional, intent(out) :: LDENS_CONST
  double precision, optional, intent(out) :: DENS
  integer :: i
  call ZDPlasKin_get_species_index(string,i)
  if( present( DENS)       )  DENS       = density(i)
  if( present(LDENS_CONST) ) LDENS_CONST = density_constant(i)
  return
end subroutine ZDPlasKin_get_density
!-----------------------------------------------------------------------------------------------------------------------------------
!
! get total densities
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_get_density_total(ALL_SPECIES,ALL_NEUTRAL,ALL_ION_POSITIVE,ALL_ION_NEGATIVE,ALL_CHARGE)
  double precision, optional, intent(out) :: ALL_SPECIES, ALL_NEUTRAL, ALL_ION_POSITIVE, ALL_ION_NEGATIVE, ALL_CHARGE
  if(.not. lZDPlasKin_init) call ZDPlasKin_init()
  if( present(ALL_SPECIES)      ) ALL_SPECIES      = sum(density(:))
  if( present(ALL_NEUTRAL)      ) ALL_NEUTRAL      = sum(density(:), mask = species_charge(:)==0)
  if( present(ALL_ION_POSITIVE) ) ALL_ION_POSITIVE = sum(density(:), mask = species_charge(:)>0)
  if( present(ALL_ION_NEGATIVE) ) ALL_ION_NEGATIVE = sum(density(:), mask = species_charge(:)<0) - density(species_electrons)
  if( present(ALL_CHARGE)       ) ALL_CHARGE       = sum(density(:) * dble(species_charge(:)))
  return
end subroutine ZDPlasKin_get_density_total
!-----------------------------------------------------------------------------------------------------------------------------------
!
! get species source terms & reaction rates
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_get_rates(SOURCE_TERMS,REACTION_RATES,SOURCE_TERMS_MATRIX,MEAN_DENSITY, &
                               MEAN_SOURCE_TERMS,MEAN_REACTION_RATES,MEAN_SOURCE_TERMS_MATRIX)
  double precision, optional, intent(out) :: SOURCE_TERMS(species_max), REACTION_RATES(reactions_max), &
                                             SOURCE_TERMS_MATRIX(species_max,reactions_max), MEAN_DENSITY(species_max), &
                                             MEAN_SOURCE_TERMS(species_max), MEAN_REACTION_RATES(reactions_max), &
                                             MEAN_SOURCE_TERMS_MATRIX(species_max,reactions_max)
  if(.not. lZDPlasKin_init) call ZDPlasKin_init()
  if(present(SOURCE_TERMS) .or. present(REACTION_RATES) .or. present(SOURCE_TERMS_MATRIX)) then
    dens_loc(1:species_max,0) = 0.5d0 * ( density(:) + abs( density(:) ) )
    if( lgas_heating ) dens_loc(species_max+1,0) = ZDPlasKin_cfg(1)
    call ZDPlasKin_fex(vode_neq,tsav,dens_loc(:,0),dens_loc(:,1))
    if( present(SOURCE_TERMS)                 ) SOURCE_TERMS(:)   = dens_loc(1:species_max,1)
    if( present(REACTION_RATES)               ) REACTION_RATES(:) = rrt(:)
    if( present(SOURCE_TERMS_MATRIX)          ) call ZDPlasKin_reac_source_matrix(rrt(:),SOURCE_TERMS_MATRIX(:,:))
  endif
  if(present(MEAN_DENSITY)        .or. present(MEAN_SOURCE_TERMS) .or. &
     present(MEAN_REACTION_RATES) .or. present(MEAN_SOURCE_TERMS_MATRIX)) then
    if( lstat_accum ) then
      if(stat_time > 0.0d0) then
        if( present(MEAN_DENSITY)             ) MEAN_DENSITY        = stat_dens(:) / stat_time
        if( present(MEAN_SOURCE_TERMS)        ) MEAN_SOURCE_TERMS   = stat_src(:)  / stat_time
        if( present(MEAN_REACTION_RATES)      ) MEAN_REACTION_RATES = stat_rrt(:)  / stat_time
        if( present(MEAN_SOURCE_TERMS_MATRIX) ) then
          call ZDPlasKin_reac_source_matrix(stat_rrt(:),MEAN_SOURCE_TERMS_MATRIX(:,:))
          MEAN_SOURCE_TERMS_MATRIX(:,:)  =  MEAN_SOURCE_TERMS_MATRIX(:,:) / stat_time
        endif
      else
        dens_loc(1:species_max,0) = 0.5d0 * ( density(:) + abs( density(:) ) )
        if( lgas_heating ) dens_loc(species_max+1,0) = ZDPlasKin_cfg(1)
        call ZDPlasKin_fex(vode_neq,tsav,dens_loc(:,0),dens_loc(:,1))
        if( present(MEAN_DENSITY)             ) MEAN_DENSITY        = density(:)
        if( present(MEAN_SOURCE_TERMS)        ) MEAN_SOURCE_TERMS   = dens_loc(1:species_max,1)
        if( present(MEAN_REACTION_RATES)      ) MEAN_REACTION_RATES = rrt(:)
        if( present(MEAN_SOURCE_TERMS_MATRIX) ) call ZDPlasKin_reac_source_matrix(rrt(:),MEAN_SOURCE_TERMS_MATRIX(:,:))
      endif
    else
      call ZDPlasKin_stop("ZDPlasKin ERROR: set statistics acquisition ON before (subroutine ZDPlasKin_get_rates)")
    endif
  endif
  return
end subroutine ZDPlasKin_get_rates
!-----------------------------------------------------------------------------------------------------------------------------------
!
! set config
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_set_config(ATOL,RTOL,SILENCE_MODE,STAT_ACCUM,QTPLASKIN_SAVE,BOLSIG_EE_FRAC,BOLSIG_IGNORE_GAS_TEMPERATURE)
  use dvode_f90_m, only : set_intermediate_opts
  implicit none
  logical, optional, intent(in) :: SILENCE_MODE, STAT_ACCUM, QTPLASKIN_SAVE, BOLSIG_IGNORE_GAS_TEMPERATURE
  double precision, optional, intent(in) :: ATOL, RTOL, BOLSIG_EE_FRAC
  integer :: i
  logical, save :: lfirst = .true.
  integer, save :: bounded_components(vode_neq)
  double precision :: atol_loc, rtol_loc
  double precision, save :: atol_save = -1.0d0, rtol_save = -1.0d0
  if( lfirst ) then
    if(.not. lZDPlasKin_init) call ZDPlasKin_init()
    do i = 1, vode_neq
      bounded_components(i) = i
    enddo
    lfirst = .false.
  endif
  if( present(SILENCE_MODE) ) lprint = ( .not. SILENCE_MODE )
  if( present(BOLSIG_EE_FRAC) ) bolsig_eecol_frac = 0.5d0 * ( BOLSIG_EE_FRAC + abs(BOLSIG_EE_FRAC) )
  if( present(BOLSIG_IGNORE_GAS_TEMPERATURE) ) lbolsig_ignore_gas_temp = BOLSIG_IGNORE_GAS_TEMPERATURE
  if( present(STAT_ACCUM) ) then
    if( lprint ) then
      if(lstat_accum .neqv. STAT_ACCUM) then
        if( STAT_ACCUM ) then
          write(*,"(A)") "ZDPlasKin INFO: set statistic acquisition ON ..."
        else
          write(*,"(A)") "ZDPlasKin INFO: set statistic acquisition OFF ..."
        endif
      elseif( STAT_ACCUM ) then
          write(*,"(A)") "ZDPlasKin INFO: reset statistic acquisition data ..."
      endif
    endif
    stat_dens(:) = 0.0d0
    stat_src(:)  = 0.0d0
    stat_rrt(:)  = 0.0d0
    stat_time    = 0.0d0
    lstat_accum  = STAT_ACCUM
  endif
  if( present(QTPLASKIN_SAVE) ) then
    if( lprint ) then
      if(lqtplaskin .neqv. QTPLASKIN_SAVE) then
        if( QTPLASKIN_SAVE ) then
          write(*,"(A)") "ZDPlasKin INFO: set autosave in QTplaskin format ON ..."
        else
          write(*,"(A)") "ZDPlasKin INFO: set autosave in QTplaskin format OFF ..."
        endif
      endif
    endif
    lqtplaskin = QTPLASKIN_SAVE
  endif
  if( present(ATOL) ) then
    atol_loc = ATOL
  else
    atol_loc = atol_save
  endif
  if( present(RTOL) ) then
    rtol_loc = RTOL
  else
    rtol_loc = rtol_save
  endif
  if(min(atol_loc,rtol_loc)<0.0d0 .or. max(atol_loc,rtol_loc)<=0.0d0) &
    call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined ATOL/RTOL (ZDPlasKin_set_config)")
  if(atol_loc/=atol_save .or. rtol_loc/=rtol_save) then
    atol_save = atol_loc
    rtol_save = rtol_loc
    if( lprint ) write(*,"(2(A,1pd9.2),A)") "ZDPlasKin INFO: set accuracy", atol_save, " (absolute) &", rtol_save, " (relative)"
    dens_loc(:,0) = 0.0d0
    dens_loc(:,1) = huge(dens_loc)
    vode_options  = set_intermediate_opts(abserr=atol_save,relerr=rtol_save, &
                                          dense_j=.true.,user_supplied_jacobian=.true., &
                                          constrained=bounded_components(:),clower=dens_loc(:,0),cupper=dens_loc(:,1))
    if(vode_istate /= 1) vode_istate = 3
  endif
  return
end subroutine ZDPlasKin_set_config
!-----------------------------------------------------------------------------------------------------------------------------------
!
! set/get conditions
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_set_conditions(GAS_TEMPERATURE,REDUCED_FREQUENCY,REDUCED_FIELD, &
                                    ELEC_TEMPERATURE,GAS_HEATING,SPEC_HEAT_RATIO,HEAT_SOURCE,SOFT_RESET)
  implicit none
  double precision, optional, intent(in) :: GAS_TEMPERATURE, REDUCED_FREQUENCY, REDUCED_FIELD, ELEC_TEMPERATURE, &
                                            SPEC_HEAT_RATIO, HEAT_SOURCE
  logical,          optional, intent(in) :: GAS_HEATING, SOFT_RESET
  if(.not. lZDPlasKin_init) call ZDPlasKin_init()
  if( present(GAS_TEMPERATURE) ) then
    if(GAS_TEMPERATURE <= 0.0d0) &
      call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined GAS_TEMPERATURE (subroutine ZDPlasKin_set_conditions)")
    ZDPlasKin_cfg(1) = GAS_TEMPERATURE
  endif
  if( present(REDUCED_FREQUENCY) ) then
    if(REDUCED_FREQUENCY < 0.0d0) &
      call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined REDUCED_FREQUENCY (subroutine ZDPlasKin_set_conditions)")
    ZDPlasKin_cfg(2) = REDUCED_FREQUENCY
  endif
  if( present(REDUCED_FIELD) ) then
    if(REDUCED_FIELD < 0.0d0) &
      call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined REDUCED_FIELD (subroutine ZDPlasKin_set_conditions)")
    ZDPlasKin_cfg(3) = REDUCED_FIELD
  endif
  if( present(SOFT_RESET) ) then
    if( SOFT_RESET ) vode_istate = 1
  endif
  if( present(ELEC_TEMPERATURE) ) then
    if(ELEC_TEMPERATURE < 0.0d0) &
      call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined ELEC_TEMPERATURE (subroutine ZDPlasKin_set_conditions)")
    if(ELEC_TEMPERATURE > 0.0d0) then
      lbolsig_Maxwell_EEDF = .true.
    else
      lbolsig_Maxwell_EEDF = .false.
    endif
    ZDPlasKin_cfg(4) = ELEC_TEMPERATURE
  endif
  if( present(GAS_HEATING) ) then
    if(lgas_heating .neqv. GAS_HEATING) then
      if( GAS_HEATING ) then
        if(present(SPEC_HEAT_RATIO) .or. ZDPlasKin_cfg(13)>0.0d0) then
          if( lprint ) write(*,"(A)") "ZDPlasKin INFO: set gas heating ON ..."
        else
          ZDPlasKin_cfg(13) = 2.0d0/3.0d0
          if( lprint ) write(*,"(A,1pd9.2)") "ZDPlasKin INFO: set gas heating ON; specific heat ratio =", ZDPlasKin_cfg(13) + 1.0d0
        endif
      else
        if( lprint ) write(*,"(A)") "ZDPlasKin INFO: set gas heating OFF ..."
      endif
      lgas_heating = GAS_HEATING
    endif
  endif
  if( present(SPEC_HEAT_RATIO) ) then
    if(SPEC_HEAT_RATIO > 1.0d0) then
      ZDPlasKin_cfg(13) = SPEC_HEAT_RATIO - 1.0d0
      if( lprint ) write(*,"(A,1pd9.2)") "ZDPlasKin INFO: set specific heat ratio =", ZDPlasKin_cfg(13) + 1.0d0
    else
      call ZDPlasKin_stop("ZDPlasKin ERROR: wrong value of SPEC_HEAT_RATIO (subroutine ZDPlasKin_set_conditions)")
    endif
  endif
  if( present(HEAT_SOURCE) ) then
    ZDPlasKin_cfg(14) = HEAT_SOURCE
    if( lprint ) write(*,"(A,1pd9.2,A)") "ZDPlasKin INFO: set heat source =", ZDPlasKin_cfg(14), " W/cm3"
  endif
end subroutine ZDPlasKin_set_conditions
subroutine ZDPlasKin_get_conditions(GAS_TEMPERATURE,REDUCED_FREQUENCY,REDUCED_FIELD, &
                                    ELEC_TEMPERATURE,ELEC_DRIFT_VELOCITY,ELEC_DIFF_COEFF,ELEC_MOBILITY_N, &
                                    ELEC_MU_EPS_N,ELEC_DIFF_EPS_N,ELEC_FREQUENCY_N, &
                                    ELEC_POWER_N,ELEC_POWER_ELASTIC_N,ELEC_POWER_INELASTIC_N,ELEC_EEDF)
  implicit none
  double precision, optional, intent(out) :: GAS_TEMPERATURE, REDUCED_FREQUENCY, REDUCED_FIELD, &
                                             ELEC_TEMPERATURE, ELEC_DRIFT_VELOCITY, ELEC_DIFF_COEFF, ELEC_MOBILITY_N, &
                                             ELEC_MU_EPS_N, ELEC_DIFF_EPS_N, ELEC_FREQUENCY_N, &
                                             ELEC_POWER_N, ELEC_POWER_ELASTIC_N, ELEC_POWER_INELASTIC_N
  double precision, optional, dimension(:,:), intent(out) :: ELEC_EEDF
  integer :: i
  double precision :: x,y
  if(.not. lZDPlasKin_init) call ZDPlasKin_init()
  if( present(ELEC_EEDF) ) then
    call ZDPlasKin_bolsig_rates(lbolsig_force=.true.)
  else
    call ZDPlasKin_bolsig_rates()
  endif
  if( present(GAS_TEMPERATURE)        ) GAS_TEMPERATURE        = ZDPlasKin_cfg(1)
  if( present(REDUCED_FREQUENCY)      ) REDUCED_FREQUENCY      = ZDPlasKin_cfg(2)
  if( present(REDUCED_FIELD)          ) REDUCED_FIELD          = ZDPlasKin_cfg(3)
  if( present(ELEC_TEMPERATURE)       ) ELEC_TEMPERATURE       = ZDPlasKin_cfg(4)
  if( present(ELEC_DRIFT_VELOCITY)    ) ELEC_DRIFT_VELOCITY    = ZDPlasKin_cfg(5) * ZDPlasKin_cfg(3) * 1.0d-17
  if( present(ELEC_DIFF_COEFF)        ) ELEC_DIFF_COEFF        = ZDPlasKin_cfg(6)
  if( present(ELEC_MOBILITY_N)        ) ELEC_MOBILITY_N        = ZDPlasKin_cfg(5)
  if( present(ELEC_MU_EPS_N)          ) ELEC_MU_EPS_N          = ZDPlasKin_cfg(7)
  if( present(ELEC_DIFF_EPS_N)        ) ELEC_DIFF_EPS_N        = ZDPlasKin_cfg(8)
  if( present(ELEC_FREQUENCY_N)       ) ELEC_FREQUENCY_N       = ZDPlasKin_cfg(9)
  if( present(ELEC_POWER_N)           ) ELEC_POWER_N           = ZDPlasKin_cfg(10)
  if( present(ELEC_POWER_ELASTIC_N)   ) ELEC_POWER_ELASTIC_N   = ZDPlasKin_cfg(11)
  if( present(ELEC_POWER_INELASTIC_N) ) ELEC_POWER_INELASTIC_N = ZDPlasKin_cfg(12)
  if( present(ELEC_EEDF) ) then
    ELEC_EEDF = 0d0
  	 if( size(ELEC_EEDF,dim=1) < 2 ) then
      if(lprint) write(*,"(A)") &
  	     "ZDPlasKin WARNING: insufficient first dimention of array ELEC_EEDF (subroutine ZDPlasKin_get_conditions)"
  	  else
  		y = 1.0d0
  		do i = 1, size(ELEC_EEDF,dim=2)
  		  call ZDPlasKin_bolsig_GetEEDF(i,x,y)
  		  if( x >= 0d0 .and. y > 0d0) then
  			ELEC_EEDF(1,i) = x
  			ELEC_EEDF(2,i) = y
  		  else
  			exit
  		  endif
  		enddo
  		if(lprint .and. y>0d0) write(*,"(A)") &
  		  "ZDPlasKin WARNING: insufficient second dimention of array ELEC_EEDF (subroutine ZDPlasKin_get_conditions)"
     endif
  endif
end subroutine ZDPlasKin_get_conditions
!-----------------------------------------------------------------------------------------------------------------------------------
!
! reset
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_reset()
  implicit none
  vode_istate         =  1
  density(:)          =  0.0d0
  ZDPlasKin_cfg(:)    =  0.0d0
  ldensity_constant   = .false.
  density_constant(:) = .false.
  lreaction_block(:)  = .false.
  lprint              = .true.
  lstat_accum         = .false.
  lqtplaskin          = .false.
  lgas_heating        = .false.
  bolsig_eecol_frac       = bolsig_eecol_frac_def
  lbolsig_ignore_gas_temp = .false.
  lbolsig_Maxwell_EEDF    = .false.
  write(*,"(A)") "ZDPlasKin INFO: reset data and configuration"
  call ZDPlasKin_set_config(ATOL=vode_atol,RTOL=vode_rtol)
  return
end subroutine ZDPlasKin_reset
!-----------------------------------------------------------------------------------------------------------------------------------
!
! stop
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_stop(string)
  implicit none
  character(*), intent(in) :: string
  if(string /= "") write(*,"(A)") trim(string)
  write(*,"(A,$)") "PRESS ENTER TO EXIT ... "
  read(*,*)
  stop
end subroutine ZDPlasKin_stop
!-----------------------------------------------------------------------------------------------------------------------------------
!
! save data to file
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_write_file(FILE_SPECIES,FILE_REACTIONS,FILE_SOURCE_MATRIX,FILE_UNIT)
  implicit none
  character(*), optional, intent(in) :: FILE_SPECIES, FILE_REACTIONS, FILE_SOURCE_MATRIX
  integer, optional, intent(in) :: FILE_UNIT
  logical :: lerror
  integer :: i
  if( present(FILE_UNIT) ) ifile_unit = FILE_UNIT
  if( present(FILE_SPECIES) ) then
    lerror = .true.
    open(ifile_unit,file=trim(adjustl(FILE_SPECIES)),action="write",err=100)
    do i = 1, species_max
      write(ifile_unit,111,err=100) i, species_name(i)
    enddo
    lerror = .false.
100 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file <" &
                                    // trim(adjustl(FILE_SPECIES)) // "> (subroutine ZDPlasKin_write_file)")
    close(ifile_unit)
111 format(i2,1x,A9)
  endif
  if( present(FILE_REACTIONS) ) then
    lerror = .true.
    open(ifile_unit,file=trim(adjustl(FILE_REACTIONS)),action="write",err=200)
    do i = 1, reactions_max
      write(ifile_unit,211,err=200) i, reaction_sign(i)
    enddo
    lerror = .false.
200 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file <" &
                                    // trim(adjustl(FILE_REACTIONS)) // "> (subroutine ZDPlasKin_write_file)")
    close(ifile_unit)
211 format(i3,1x,A30)
  endif
  if( present(FILE_SOURCE_MATRIX) ) then
    if( lstat_accum ) then
      call ZDPlasKin_reac_source_matrix(stat_rrt(:),mrtm(:,:))
      if(stat_time > 0.0d0) mrtm(:,:) = mrtm(:,:) / stat_time
    else
      call ZDPlasKin_stop("ZDPlasKin ERROR: set statistics acquisition ON before (subroutine ZDPlasKin_write_file)")
    endif
    lerror = .true.
    open(ifile_unit,file=trim(adjustl(FILE_SOURCE_MATRIX)),action="write",err=300)
    write(ifile_unit,311,err=300) ( i, i = 1, species_max )
    write(ifile_unit,312,err=300) "N", "reaction", ( trim(species_name(i)), i = 1, species_max )
    do i = 1, reactions_max
      write(ifile_unit,313,err=300) i, reaction_sign(i), mrtm(:,i)
    enddo
    lerror = .false.
300 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file <" &
                                    // trim(adjustl(FILE_SOURCE_MATRIX)) // "> (subroutine ZDPlasKin_write_file)")
    close(ifile_unit)
311 format(351x,39(1x,i9))
312 format(A3,1x,A30,1x,39(1x,A9))
313 format(i3,1x,A30,1x,39(1x,1pd9.2))
  endif
  return
end subroutine ZDPlasKin_write_file
!-----------------------------------------------------------------------------------------------------------------------------------
!
! save data in qtplaskin format
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_write_qtplaskin(time,LFORCE_WRITE)
  implicit none
  double precision, intent(in) :: time
  logical, optional, intent(in) :: LFORCE_WRITE
  integer, parameter :: idef_data = 5
  character(24), parameter :: qtplaskin_names(idef_data) = (/ "Reduced field [Td]      ", "Gas temperature [K]     ", &
                                  "Electron temperature [K]", "Current density [A/cm2] ", "Power density [W/cm3]   " /)
  double precision, save :: densav(0:species_max,2) = -huge(densav)
  double precision :: rtol, cond(idef_data)
  logical, save :: lfirst = .true.
  logical :: lerror
  integer, save :: iuser_data = 0
  integer :: i
  if( time < densav(0,1) ) lfirst = .true.
  if( lfirst ) then
    call ZDPlasKin_write_file(FILE_SPECIES="qt_species_list.txt",FILE_REACTIONS="qt_reactions_list.txt")
    if( allocated(qtplaskin_user_data) ) then
      iuser_data = size(qtplaskin_user_data)
      iuser_data = min(iuser_data,90)
      if( iuser_data > 0 ) then
        if( allocated(qtplaskin_user_names) ) then
          if( size(qtplaskin_user_names) /= iuser_data ) deallocate(qtplaskin_user_names)
        endif
        if( .not. allocated(qtplaskin_user_names) ) then
          allocate(qtplaskin_user_names(iuser_data))
          do i = 1, iuser_data
            write(qtplaskin_user_names(i),"(A,i2.2)") "user defined #", i
          enddo
        endif
      endif
    endif
    lerror = .true.
    open(ifile_unit,file="qt_conditions_list.txt",action="write",err=100)
    do i = 1, idef_data
      write(ifile_unit,"(i3,1x,A)",err=100) i, trim(adjustl(qtplaskin_names(i)))
    enddo
    if( iuser_data > 0 ) then
      do i = 1, iuser_data
        write(ifile_unit,"(i3,1x,A)",err=100) (i+idef_data), trim(adjustl(qtplaskin_user_names(i)))
      enddo
    endif
    lerror = .false.
100 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file " // &
                                     "<qt_conditions_list.txt> (subroutine writer_save_qtplaskin)")
    close(ifile_unit)
    rrt(:) = 1.0d0
    call ZDPlasKin_reac_source_matrix(rrt(:),mrtm(:,:))
    open(ifile_unit,file="qt_matrix.txt",action="write",err=200)
    do i = 1, species_max
      write(ifile_unit,"(114(i3))",err=200) int(mrtm(i,:))
    enddo
    lerror = .false.
200 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file " // &
                                     "<qt_matrix.txt> (subroutine writer_save_qtplaskin)")
    close(ifile_unit)
    open(ifile_unit,file="qt_densities.txt",action="write",err=300)
    write(ifile_unit,"(1x,A14,39(111x,i2.2))",err=300) "Time_s", ( i, i = 1, species_max )
    lerror = .false.
300 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file " // &
                                     "<qt_densities.txt> (subroutine writer_save_qtplaskin)")
    close(ifile_unit)
    open(ifile_unit,file="qt_conditions.txt",action="write",err=400)
    write(ifile_unit,"(1x,A12,$)",err=400) "Time_s"
    do i = 1, idef_data + iuser_data
      write(ifile_unit,"(11x,i2.2,$)",err=400) i
    enddo
    write(ifile_unit,*,err=400)
    lerror = .false.
400 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file " // &
                                     "<qt_conditions.txt> (subroutine writer_save_qtplaskin)")
    close(ifile_unit)
    open(ifile_unit,file="qt_rates.txt",action="write",err=500)
    write(ifile_unit,"(1x,A12,114(101x,i3.3))",err=500) "Time_s", ( i, i = 1, reactions_max )
    lerror = .false.
500 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file " // &
                                     "<qt_rates.txt> (subroutine writer_save_qtplaskin)")
    close(ifile_unit)
  endif
  if( present(LFORCE_WRITE) ) then
    if( LFORCE_WRITE ) lfirst = .true.
  endif
  rtol = 10.0d0 ** ( floor( log10( abs(densav(0,1)) + tiny(rtol) ) ) - 6 )
  if( ( time - densav(0,1) ) >= rtol .or. lfirst ) then
    densav(0,2) = time
    densav(1:species_max,2) = density(:)
    where( densav(:,2) < 1.0d-99 ) densav(:,2) = 0.0d0
    if( time > 2.0d0 * densav(0,1) ) then
      rtol = huge(rtol)
    else
      rtol = maxval( abs(densav(1:,1)-densav(1:,2)) / ( abs(densav(1:,1)+densav(1:,2))/2.0d0 + qtplaskin_atol ) )
    endif
    if( rtol > qtplaskin_rtol .or. lfirst ) then
      open(ifile_unit,file="qt_densities.txt",access="append")
      write(ifile_unit,"(1pe15.6,39(1pe13.4))") densav(0,2), densav(1:,2)
      close(ifile_unit)
      open(ifile_unit,file="qt_conditions.txt",access="append")
      cond(1) = ZDPlasKin_cfg(3)
      cond(2) = ZDPlasKin_cfg(1)
      cond(3) = ZDPlasKin_cfg(4)
      cond(4) = q_elem * density(species_electrons) * ZDPlasKin_cfg(5) * ZDPlasKin_cfg(3) * 1.0d-17
      call ZDPlasKin_get_density_total(ALL_NEUTRAL=cond(5))
      cond(5) = cond(4) * cond(5) * ZDPlasKin_cfg(3) * 1.0d-17
      where( abs(cond(:)) < 1.0d-99 ) cond(:) = 0.0d0
      write(ifile_unit,"(6(1pe13.4),$)") densav(0,2), cond(:)
      if( iuser_data > 0 ) then
        where( abs(qtplaskin_user_data(1:iuser_data)) < 1.0d-99 ) qtplaskin_user_data(1:iuser_data) = 0.0d0
        write(ifile_unit,"(90(1pe13.4))") qtplaskin_user_data(1:iuser_data)
      else
        write(ifile_unit,*)
      endif
      close(ifile_unit)
      call ZDPlasKin_get_rates(REACTION_RATES=rrt_loc)
      where( abs(rrt_loc(:)) < 1.0d-99 ) rrt_loc(:) = 0.0d0
      open(ifile_unit,file="qt_rates.txt",access="append")
      write(ifile_unit,"(115(1pe13.4))") densav(0,2), rrt_loc(:)
      close(ifile_unit)
      densav(:,1) = densav(:,2)
    endif
  endif
  lfirst = .false.
  lqtplaskin_first = .false.
  return
end subroutine ZDPlasKin_write_qtplaskin
!-----------------------------------------------------------------------------------------------------------------------------------
!
! reaction sensitivity acquisition
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_reac_source_matrix(reac_rate_local,reac_source_local)
  implicit none
  double precision, intent(in)  :: reac_rate_local(reactions_max)
  double precision, intent(out) :: reac_source_local(species_max,reactions_max)
  reac_source_local(:,:) = 0.0d0
  reac_source_local(07,001) = - reac_rate_local(001) 
  reac_source_local(09,001) = + reac_rate_local(001) 
  reac_source_local(07,002) = - reac_rate_local(002) 
  reac_source_local(08,002) = + reac_rate_local(002) 
  reac_source_local(07,003) = + reac_rate_local(003) 
  reac_source_local(09,003) = - reac_rate_local(003) 
  reac_source_local(07,004) = + reac_rate_local(004) 
  reac_source_local(08,004) = - reac_rate_local(004) 
  reac_source_local(24,005) = - reac_rate_local(005) 
  reac_source_local(26,005) = + reac_rate_local(005) 
  reac_source_local(24,006) = - reac_rate_local(006) 
  reac_source_local(25,006) = + reac_rate_local(006) 
  reac_source_local(24,007) = + reac_rate_local(007) 
  reac_source_local(26,007) = - reac_rate_local(007) 
  reac_source_local(24,008) = + reac_rate_local(008) 
  reac_source_local(25,008) = - reac_rate_local(008) 
  reac_source_local(18,009) = - reac_rate_local(009) 
  reac_source_local(20,009) = + reac_rate_local(009) 
  reac_source_local(18,010) = - reac_rate_local(010) 
  reac_source_local(19,010) = + reac_rate_local(010) 
  reac_source_local(18,011) = + reac_rate_local(011) 
  reac_source_local(20,011) = - reac_rate_local(011) 
  reac_source_local(18,012) = + reac_rate_local(012) 
  reac_source_local(19,012) = - reac_rate_local(012) 
  reac_source_local(12,013) = - reac_rate_local(013) 
  reac_source_local(14,013) = + reac_rate_local(013) 
  reac_source_local(12,014) = - reac_rate_local(014) 
  reac_source_local(15,014) = + reac_rate_local(014) 
  reac_source_local(12,015) = - reac_rate_local(015) 
  reac_source_local(13,015) = + reac_rate_local(015) 
  reac_source_local(12,016) = + reac_rate_local(016) 
  reac_source_local(14,016) = - reac_rate_local(016) 
  reac_source_local(12,017) = + reac_rate_local(017) 
  reac_source_local(15,017) = - reac_rate_local(017) 
  reac_source_local(12,018) = + reac_rate_local(018) 
  reac_source_local(13,018) = - reac_rate_local(018) 
  reac_source_local(32,019) = - reac_rate_local(019) 
  reac_source_local(34,019) = + reac_rate_local(019) 
  reac_source_local(32,020) = - reac_rate_local(020) 
  reac_source_local(33,020) = + reac_rate_local(020) 
  reac_source_local(32,021) = + reac_rate_local(021) 
  reac_source_local(34,021) = - reac_rate_local(021) 
  reac_source_local(32,022) = + reac_rate_local(022) 
  reac_source_local(33,022) = - reac_rate_local(022) 
  reac_source_local(28,023) = - reac_rate_local(023) 
  reac_source_local(29,023) = + reac_rate_local(023) 
  reac_source_local(28,024) = + reac_rate_local(024) 
  reac_source_local(29,024) = - reac_rate_local(024) 
  reac_source_local(01,025) = + reac_rate_local(025) 
  reac_source_local(05,025) = - reac_rate_local(025) 
  reac_source_local(06,025) = + reac_rate_local(025) 
  reac_source_local(01,026) = + reac_rate_local(026) 
  reac_source_local(07,026) = - reac_rate_local(026) 
  reac_source_local(10,026) = + reac_rate_local(026) 
  reac_source_local(01,027) = + reac_rate_local(027) 
  reac_source_local(09,027) = - reac_rate_local(027) 
  reac_source_local(10,027) = + reac_rate_local(027) 
  reac_source_local(01,028) = + reac_rate_local(028) 
  reac_source_local(08,028) = - reac_rate_local(028) 
  reac_source_local(10,028) = + reac_rate_local(028) 
  reac_source_local(01,029) = + reac_rate_local(029) 
  reac_source_local(24,029) = - reac_rate_local(029) 
  reac_source_local(27,029) = + reac_rate_local(029) 
  reac_source_local(01,030) = + reac_rate_local(030) 
  reac_source_local(26,030) = - reac_rate_local(030) 
  reac_source_local(27,030) = + reac_rate_local(030) 
  reac_source_local(01,031) = + reac_rate_local(031) 
  reac_source_local(25,031) = - reac_rate_local(031) 
  reac_source_local(27,031) = + reac_rate_local(031) 
  reac_source_local(01,032) = + reac_rate_local(032) 
  reac_source_local(18,032) = - reac_rate_local(032) 
  reac_source_local(21,032) = + reac_rate_local(032) 
  reac_source_local(01,033) = + reac_rate_local(033) 
  reac_source_local(20,033) = - reac_rate_local(033) 
  reac_source_local(21,033) = + reac_rate_local(033) 
  reac_source_local(01,034) = + reac_rate_local(034) 
  reac_source_local(19,034) = - reac_rate_local(034) 
  reac_source_local(21,034) = + reac_rate_local(034) 
  reac_source_local(01,035) = + reac_rate_local(035) 
  reac_source_local(12,035) = - reac_rate_local(035) 
  reac_source_local(16,035) = + reac_rate_local(035) 
  reac_source_local(01,036) = + reac_rate_local(036) 
  reac_source_local(14,036) = - reac_rate_local(036) 
  reac_source_local(16,036) = + reac_rate_local(036) 
  reac_source_local(01,037) = + reac_rate_local(037) 
  reac_source_local(15,037) = - reac_rate_local(037) 
  reac_source_local(16,037) = + reac_rate_local(037) 
  reac_source_local(01,038) = + reac_rate_local(038) 
  reac_source_local(13,038) = - reac_rate_local(038) 
  reac_source_local(16,038) = + reac_rate_local(038) 
  reac_source_local(01,039) = + reac_rate_local(039) 
  reac_source_local(32,039) = - reac_rate_local(039) 
  reac_source_local(35,039) = + reac_rate_local(039) 
  reac_source_local(01,040) = + reac_rate_local(040) 
  reac_source_local(34,040) = - reac_rate_local(040) 
  reac_source_local(35,040) = + reac_rate_local(040) 
  reac_source_local(01,041) = + reac_rate_local(041) 
  reac_source_local(33,041) = - reac_rate_local(041) 
  reac_source_local(35,041) = + reac_rate_local(041) 
  reac_source_local(01,042) = + reac_rate_local(042) 
  reac_source_local(28,042) = - reac_rate_local(042) 
  reac_source_local(30,042) = + reac_rate_local(042) 
  reac_source_local(01,043) = + reac_rate_local(043) 
  reac_source_local(29,043) = - reac_rate_local(043) 
  reac_source_local(30,043) = + reac_rate_local(043) 
  reac_source_local(01,044) = - reac_rate_local(044) 
  reac_source_local(05,044) = + reac_rate_local(044) 
  reac_source_local(06,044) = - reac_rate_local(044) 
  reac_source_local(01,045) = - reac_rate_local(045) 
  reac_source_local(07,045) = + reac_rate_local(045) 
  reac_source_local(10,045) = - reac_rate_local(045) 
  reac_source_local(01,046) = - reac_rate_local(046) 
  reac_source_local(24,046) = + reac_rate_local(046) 
  reac_source_local(27,046) = - reac_rate_local(046) 
  reac_source_local(01,047) = - reac_rate_local(047) 
  reac_source_local(18,047) = + reac_rate_local(047) 
  reac_source_local(21,047) = - reac_rate_local(047) 
  reac_source_local(01,048) = - reac_rate_local(048) 
  reac_source_local(12,048) = + reac_rate_local(048) 
  reac_source_local(16,048) = - reac_rate_local(048) 
  reac_source_local(01,049) = - reac_rate_local(049) 
  reac_source_local(32,049) = + reac_rate_local(049) 
  reac_source_local(35,049) = - reac_rate_local(049) 
  reac_source_local(01,050) = - reac_rate_local(050) 
  reac_source_local(28,050) = + reac_rate_local(050) 
  reac_source_local(30,050) = - reac_rate_local(050) 
  reac_source_local(05,051) = + reac_rate_local(051) 
  reac_source_local(07,051) = - reac_rate_local(051) 
  reac_source_local(38,051) = + reac_rate_local(051) 
  reac_source_local(05,052) = + reac_rate_local(052) 
  reac_source_local(09,052) = - reac_rate_local(052) 
  reac_source_local(38,052) = + reac_rate_local(052) 
  reac_source_local(05,053) = + reac_rate_local(053) 
  reac_source_local(08,053) = - reac_rate_local(053) 
  reac_source_local(38,053) = + reac_rate_local(053) 
  reac_source_local(04,054) = + reac_rate_local(054) 
  reac_source_local(07,054) = - reac_rate_local(054) 
  reac_source_local(39,054) = + reac_rate_local(054) 
  reac_source_local(04,055) = + reac_rate_local(055) 
  reac_source_local(09,055) = - reac_rate_local(055) 
  reac_source_local(39,055) = + reac_rate_local(055) 
  reac_source_local(04,056) = + reac_rate_local(056) 
  reac_source_local(08,056) = - reac_rate_local(056) 
  reac_source_local(39,056) = + reac_rate_local(056) 
  reac_source_local(18,057) = + reac_rate_local(057) 
  reac_source_local(24,057) = - reac_rate_local(057) 
  reac_source_local(39,057) = + reac_rate_local(057) 
  reac_source_local(18,058) = + reac_rate_local(058) 
  reac_source_local(26,058) = - reac_rate_local(058) 
  reac_source_local(39,058) = + reac_rate_local(058) 
  reac_source_local(18,059) = + reac_rate_local(059) 
  reac_source_local(25,059) = - reac_rate_local(059) 
  reac_source_local(39,059) = + reac_rate_local(059) 
  reac_source_local(01,060) = - reac_rate_local(060) 
  reac_source_local(17,060) = + reac_rate_local(060) 
  reac_source_local(23,060) = - reac_rate_local(060) 
  reac_source_local(38,060) = + reac_rate_local(060) * 2.d0
  reac_source_local(12,061) = + reac_rate_local(061) 
  reac_source_local(18,061) = - reac_rate_local(061) 
  reac_source_local(39,061) = + reac_rate_local(061) 
  reac_source_local(12,062) = + reac_rate_local(062) 
  reac_source_local(19,062) = - reac_rate_local(062) 
  reac_source_local(39,062) = + reac_rate_local(062) 
  reac_source_local(12,063) = + reac_rate_local(063) 
  reac_source_local(20,063) = - reac_rate_local(063) 
  reac_source_local(39,063) = + reac_rate_local(063) 
  reac_source_local(01,064) = - reac_rate_local(064) 
  reac_source_local(12,064) = + reac_rate_local(064) 
  reac_source_local(23,064) = - reac_rate_local(064) 
  reac_source_local(38,064) = + reac_rate_local(064) 
  reac_source_local(39,064) = + reac_rate_local(064) 
  reac_source_local(07,065) = + reac_rate_local(065) 
  reac_source_local(12,065) = + reac_rate_local(065) 
  reac_source_local(28,065) = - reac_rate_local(065) 
  reac_source_local(07,066) = + reac_rate_local(066) 
  reac_source_local(12,066) = + reac_rate_local(066) 
  reac_source_local(29,066) = - reac_rate_local(066) 
  reac_source_local(01,067) = - reac_rate_local(067) 
  reac_source_local(12,067) = + reac_rate_local(067) 
  reac_source_local(23,067) = - reac_rate_local(067) 
  reac_source_local(38,067) = + reac_rate_local(067) * 3.d0
  reac_source_local(17,068) = + reac_rate_local(068) 
  reac_source_local(18,068) = - reac_rate_local(068) 
  reac_source_local(38,068) = + reac_rate_local(068) 
  reac_source_local(17,069) = + reac_rate_local(069) 
  reac_source_local(19,069) = - reac_rate_local(069) 
  reac_source_local(38,069) = + reac_rate_local(069) 
  reac_source_local(17,070) = + reac_rate_local(070) 
  reac_source_local(20,070) = - reac_rate_local(070) 
  reac_source_local(38,070) = + reac_rate_local(070) 
  reac_source_local(28,071) = + reac_rate_local(071) 
  reac_source_local(32,071) = - reac_rate_local(071) 
  reac_source_local(39,071) = + reac_rate_local(071) 
  reac_source_local(28,072) = + reac_rate_local(072) 
  reac_source_local(33,072) = - reac_rate_local(072) 
  reac_source_local(39,072) = + reac_rate_local(072) 
  reac_source_local(28,073) = + reac_rate_local(073) 
  reac_source_local(34,073) = - reac_rate_local(073) 
  reac_source_local(39,073) = + reac_rate_local(073) 
  reac_source_local(03,074) = + reac_rate_local(074) 
  reac_source_local(07,074) = - reac_rate_local(074) 
  reac_source_local(38,074) = + reac_rate_local(074) 
  reac_source_local(39,074) = + reac_rate_local(074) 
  reac_source_local(03,075) = + reac_rate_local(075) 
  reac_source_local(09,075) = - reac_rate_local(075) 
  reac_source_local(38,075) = + reac_rate_local(075) 
  reac_source_local(39,075) = + reac_rate_local(075) 
  reac_source_local(03,076) = + reac_rate_local(076) 
  reac_source_local(08,076) = - reac_rate_local(076) 
  reac_source_local(38,076) = + reac_rate_local(076) 
  reac_source_local(39,076) = + reac_rate_local(076) 
  reac_source_local(02,077) = + reac_rate_local(077) 
  reac_source_local(03,077) = - reac_rate_local(077) 
  reac_source_local(38,077) = + reac_rate_local(077) 
  reac_source_local(05,078) = - reac_rate_local(078) 
  reac_source_local(07,078) = + reac_rate_local(078) 
  reac_source_local(38,078) = - reac_rate_local(078) 
  reac_source_local(05,079) = - reac_rate_local(079) * 2.d0
  reac_source_local(24,079) = + reac_rate_local(079) 
  reac_source_local(18,080) = - reac_rate_local(080) 
  reac_source_local(22,080) = + reac_rate_local(080) 
  reac_source_local(38,080) = - reac_rate_local(080) 
  reac_source_local(19,081) = - reac_rate_local(081) 
  reac_source_local(22,081) = + reac_rate_local(081) 
  reac_source_local(38,081) = - reac_rate_local(081) 
  reac_source_local(20,082) = - reac_rate_local(082) 
  reac_source_local(22,082) = + reac_rate_local(082) 
  reac_source_local(38,082) = - reac_rate_local(082) 
  reac_source_local(03,083) = - reac_rate_local(083) 
  reac_source_local(07,083) = - reac_rate_local(083) 
  reac_source_local(18,083) = + reac_rate_local(083) 
  reac_source_local(38,083) = + reac_rate_local(083) 
  reac_source_local(03,084) = - reac_rate_local(084) 
  reac_source_local(09,084) = - reac_rate_local(084) 
  reac_source_local(18,084) = + reac_rate_local(084) 
  reac_source_local(38,084) = + reac_rate_local(084) 
  reac_source_local(03,085) = - reac_rate_local(085) 
  reac_source_local(08,085) = - reac_rate_local(085) 
  reac_source_local(18,085) = + reac_rate_local(085) 
  reac_source_local(38,085) = + reac_rate_local(085) 
  reac_source_local(22,086) = - reac_rate_local(086) 
  reac_source_local(24,086) = + reac_rate_local(086) 
  reac_source_local(38,086) = - reac_rate_local(086) 
  reac_source_local(22,087) = - reac_rate_local(087) * 2.d0
  reac_source_local(36,087) = + reac_rate_local(087) 
  reac_source_local(03,088) = + reac_rate_local(088) 
  reac_source_local(04,088) = - reac_rate_local(088) 
  reac_source_local(38,088) = - reac_rate_local(088) 
  reac_source_local(39,088) = + reac_rate_local(088) 
  reac_source_local(04,089) = - reac_rate_local(089) 
  reac_source_local(05,089) = - reac_rate_local(089) 
  reac_source_local(18,089) = + reac_rate_local(089) 
  reac_source_local(38,089) = + reac_rate_local(089) 
  reac_source_local(05,090) = + reac_rate_local(090) * 2.d0
  reac_source_local(22,090) = - reac_rate_local(090) 
  reac_source_local(38,090) = - reac_rate_local(090) 
  reac_source_local(05,091) = + reac_rate_local(091) 
  reac_source_local(07,091) = - reac_rate_local(091) 
  reac_source_local(10,091) = - reac_rate_local(091) 
  reac_source_local(11,091) = + reac_rate_local(091) 
  reac_source_local(05,092) = + reac_rate_local(092) 
  reac_source_local(09,092) = - reac_rate_local(092) 
  reac_source_local(10,092) = - reac_rate_local(092) 
  reac_source_local(11,092) = + reac_rate_local(092) 
  reac_source_local(05,093) = + reac_rate_local(093) 
  reac_source_local(08,093) = - reac_rate_local(093) 
  reac_source_local(10,093) = - reac_rate_local(093) 
  reac_source_local(11,093) = + reac_rate_local(093) 
  reac_source_local(05,094) = - reac_rate_local(094) 
  reac_source_local(22,094) = - reac_rate_local(094) 
  reac_source_local(32,094) = + reac_rate_local(094) 
  reac_source_local(05,095) = - reac_rate_local(095) 
  reac_source_local(17,095) = - reac_rate_local(095) 
  reac_source_local(28,095) = + reac_rate_local(095) 
  reac_source_local(28,096) = - reac_rate_local(096) 
  reac_source_local(31,096) = + reac_rate_local(096) 
  reac_source_local(38,096) = - reac_rate_local(096) 
  reac_source_local(29,097) = - reac_rate_local(097) 
  reac_source_local(31,097) = + reac_rate_local(097) 
  reac_source_local(38,097) = - reac_rate_local(097) 
  reac_source_local(04,098) = - reac_rate_local(098) 
  reac_source_local(36,098) = - reac_rate_local(098) 
  reac_source_local(37,098) = + reac_rate_local(098) 
  reac_source_local(31,099) = - reac_rate_local(099) 
  reac_source_local(32,099) = + reac_rate_local(099) 
  reac_source_local(38,099) = + reac_rate_local(099) 
  reac_source_local(39,099) = - reac_rate_local(099) 
  reac_source_local(07,100) = + reac_rate_local(100) 
  reac_source_local(11,100) = - reac_rate_local(100) 
  reac_source_local(23,100) = + reac_rate_local(100) 
  reac_source_local(24,100) = - reac_rate_local(100) 
  reac_source_local(39,100) = + reac_rate_local(100) 
  reac_source_local(07,101) = + reac_rate_local(101) 
  reac_source_local(11,101) = - reac_rate_local(101) 
  reac_source_local(23,101) = + reac_rate_local(101) 
  reac_source_local(26,101) = - reac_rate_local(101) 
  reac_source_local(39,101) = + reac_rate_local(101) 
  reac_source_local(07,102) = + reac_rate_local(102) 
  reac_source_local(11,102) = - reac_rate_local(102) 
  reac_source_local(23,102) = + reac_rate_local(102) 
  reac_source_local(25,102) = - reac_rate_local(102) 
  reac_source_local(39,102) = + reac_rate_local(102) 
  reac_source_local(07,103) = + reac_rate_local(103) 
  reac_source_local(11,103) = - reac_rate_local(103) 
  reac_source_local(18,103) = - reac_rate_local(103) 
  reac_source_local(23,103) = + reac_rate_local(103) 
  reac_source_local(07,104) = + reac_rate_local(104) 
  reac_source_local(11,104) = - reac_rate_local(104) 
  reac_source_local(19,104) = - reac_rate_local(104) 
  reac_source_local(23,104) = + reac_rate_local(104) 
  reac_source_local(07,105) = + reac_rate_local(105) 
  reac_source_local(11,105) = - reac_rate_local(105) 
  reac_source_local(20,105) = - reac_rate_local(105) 
  reac_source_local(23,105) = + reac_rate_local(105) 
  reac_source_local(31,106) = - reac_rate_local(106) 
  reac_source_local(32,106) = + reac_rate_local(106) 
  reac_source_local(38,106) = - reac_rate_local(106) 
  reac_source_local(12,107) = + reac_rate_local(107) 
  reac_source_local(17,107) = - reac_rate_local(107) 
  reac_source_local(38,107) = - reac_rate_local(107) 
  reac_source_local(39,107) = + reac_rate_local(107) 
  reac_source_local(04,108) = - reac_rate_local(108) 
  reac_source_local(32,108) = - reac_rate_local(108) 
  reac_source_local(36,108) = + reac_rate_local(108) 
  reac_source_local(04,109) = - reac_rate_local(109) 
  reac_source_local(33,109) = - reac_rate_local(109) 
  reac_source_local(36,109) = + reac_rate_local(109) 
  reac_source_local(04,110) = - reac_rate_local(110) 
  reac_source_local(34,110) = - reac_rate_local(110) 
  reac_source_local(36,110) = + reac_rate_local(110) 
  reac_source_local(06,111) = - reac_rate_local(111) 
  reac_source_local(07,111) = - reac_rate_local(111) 
  reac_source_local(23,111) = + reac_rate_local(111) 
  reac_source_local(39,111) = + reac_rate_local(111) 
  reac_source_local(06,112) = - reac_rate_local(112) 
  reac_source_local(09,112) = - reac_rate_local(112) 
  reac_source_local(23,112) = + reac_rate_local(112) 
  reac_source_local(39,112) = + reac_rate_local(112) 
  reac_source_local(06,113) = - reac_rate_local(113) 
  reac_source_local(08,113) = - reac_rate_local(113) 
  reac_source_local(23,113) = + reac_rate_local(113) 
  reac_source_local(39,113) = + reac_rate_local(113) 
  reac_source_local(38,114) = - reac_rate_local(114) * 2.d0
  reac_source_local(39,114) = + reac_rate_local(114) 
  return
end subroutine ZDPlasKin_reac_source_matrix
!-----------------------------------------------------------------------------------------------------------------------------------
!
! reaction source terms
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_fex(neq,t,y,ydot)
  implicit none
  integer,          intent(in)  :: neq
  double precision, intent(in)  :: t, y(neq)
  double precision, intent(out) :: ydot(neq)
  if( lgas_heating ) ZDPlasKin_cfg(1) = y(40)
  density(:) = y(1:species_max)
  call ZDPlasKin_reac_rates(t)
  rrt(001) = rrt(001) * density(01) * density(07) 
  rrt(002) = rrt(002) * density(01) * density(07) 
  rrt(003) = rrt(003) * density(01) * density(09) 
  rrt(004) = rrt(004) * density(01) * density(08) 
  rrt(005) = rrt(005) * density(01) * density(24) 
  rrt(006) = rrt(006) * density(01) * density(24) 
  rrt(007) = rrt(007) * density(01) * density(26) 
  rrt(008) = rrt(008) * density(01) * density(25) 
  rrt(009) = rrt(009) * density(01) * density(18) 
  rrt(010) = rrt(010) * density(01) * density(18) 
  rrt(011) = rrt(011) * density(01) * density(20) 
  rrt(012) = rrt(012) * density(01) * density(19) 
  rrt(013) = rrt(013) * density(01) * density(12) 
  rrt(014) = rrt(014) * density(01) * density(12) 
  rrt(015) = rrt(015) * density(01) * density(12) 
  rrt(016) = rrt(016) * density(01) * density(14) 
  rrt(017) = rrt(017) * density(01) * density(15) 
  rrt(018) = rrt(018) * density(01) * density(13) 
  rrt(019) = rrt(019) * density(01) * density(32) 
  rrt(020) = rrt(020) * density(01) * density(32) 
  rrt(021) = rrt(021) * density(01) * density(34) 
  rrt(022) = rrt(022) * density(01) * density(33) 
  rrt(023) = rrt(023) * density(01) * density(28) 
  rrt(024) = rrt(024) * density(01) * density(29) 
  rrt(025) = rrt(025) * density(01) * density(05) 
  rrt(026) = rrt(026) * density(01) * density(07) 
  rrt(027) = rrt(027) * density(01) * density(09) 
  rrt(028) = rrt(028) * density(01) * density(08) 
  rrt(029) = rrt(029) * density(01) * density(24) 
  rrt(030) = rrt(030) * density(01) * density(26) 
  rrt(031) = rrt(031) * density(01) * density(25) 
  rrt(032) = rrt(032) * density(01) * density(18) 
  rrt(033) = rrt(033) * density(01) * density(20) 
  rrt(034) = rrt(034) * density(01) * density(19) 
  rrt(035) = rrt(035) * density(01) * density(12) 
  rrt(036) = rrt(036) * density(01) * density(14) 
  rrt(037) = rrt(037) * density(01) * density(15) 
  rrt(038) = rrt(038) * density(01) * density(13) 
  rrt(039) = rrt(039) * density(01) * density(32) 
  rrt(040) = rrt(040) * density(01) * density(34) 
  rrt(041) = rrt(041) * density(01) * density(33) 
  rrt(042) = rrt(042) * density(01) * density(28) 
  rrt(043) = rrt(043) * density(01) * density(29) 
  rrt(044) = rrt(044) * density(01)**2 * density(06) 
  rrt(045) = rrt(045) * density(01)**2 * density(10) 
  rrt(046) = rrt(046) * density(01)**2 * density(27) 
  rrt(047) = rrt(047) * density(01)**2 * density(21) 
  rrt(048) = rrt(048) * density(01)**2 * density(16) 
  rrt(049) = rrt(049) * density(01)**2 * density(35) 
  rrt(050) = rrt(050) * density(01)**2 * density(30) 
  rrt(051) = rrt(051) * density(01) * density(07) 
  rrt(052) = rrt(052) * density(01) * density(09) 
  rrt(053) = rrt(053) * density(01) * density(08) 
  rrt(054) = rrt(054) * density(01) * density(07) 
  rrt(055) = rrt(055) * density(01) * density(09) 
  rrt(056) = rrt(056) * density(01) * density(08) 
  rrt(057) = rrt(057) * density(01) * density(24) 
  rrt(058) = rrt(058) * density(01) * density(26) 
  rrt(059) = rrt(059) * density(01) * density(25) 
  rrt(060) = rrt(060) * density(01) * density(23) 
  rrt(061) = rrt(061) * density(01) * density(18) 
  rrt(062) = rrt(062) * density(01) * density(19) 
  rrt(063) = rrt(063) * density(01) * density(20) 
  rrt(064) = rrt(064) * density(01) * density(23) 
  rrt(065) = rrt(065) * density(01) * density(28) 
  rrt(066) = rrt(066) * density(01) * density(29) 
  rrt(067) = rrt(067) * density(01) * density(23) 
  rrt(068) = rrt(068) * density(01) * density(18) 
  rrt(069) = rrt(069) * density(01) * density(19) 
  rrt(070) = rrt(070) * density(01) * density(20) 
  rrt(071) = rrt(071) * density(01) * density(32) 
  rrt(072) = rrt(072) * density(01) * density(33) 
  rrt(073) = rrt(073) * density(01) * density(34) 
  rrt(074) = rrt(074) * density(01) * density(07) 
  rrt(075) = rrt(075) * density(01) * density(09) 
  rrt(076) = rrt(076) * density(01) * density(08) 
  rrt(077) = rrt(077) * density(01) * density(03) 
  rrt(078) = rrt(078) * density(05) * density(38) 
  rrt(079) = rrt(079) * density(05)**2 
  rrt(080) = rrt(080) * density(18) * density(38) 
  rrt(081) = rrt(081) * density(19) * density(38) 
  rrt(082) = rrt(082) * density(20) * density(38) 
  rrt(083) = rrt(083) * density(03) * density(07) 
  rrt(084) = rrt(084) * density(03) * density(09) 
  rrt(085) = rrt(085) * density(03) * density(08) 
  rrt(086) = rrt(086) * density(22) * density(38) 
  rrt(087) = rrt(087) * density(22)**2 
  rrt(088) = rrt(088) * density(04) * density(38) 
  rrt(089) = rrt(089) * density(04) * density(05) 
  rrt(090) = rrt(090) * density(22) * density(38) 
  rrt(091) = rrt(091) * density(07) * density(10) 
  rrt(092) = rrt(092) * density(09) * density(10) 
  rrt(093) = rrt(093) * density(08) * density(10) 
  rrt(094) = rrt(094) * density(05) * density(22) 
  rrt(095) = rrt(095) * density(05) * density(17) 
  rrt(096) = rrt(096) * density(28) * density(38) 
  rrt(097) = rrt(097) * density(29) * density(38) 
  rrt(098) = rrt(098) * density(04) * density(36) 
  rrt(099) = rrt(099) * density(31) * density(39) 
  rrt(100) = rrt(100) * density(11) * density(24) 
  rrt(101) = rrt(101) * density(11) * density(26) 
  rrt(102) = rrt(102) * density(11) * density(25) 
  rrt(103) = rrt(103) * density(11) * density(18) 
  rrt(104) = rrt(104) * density(11) * density(19) 
  rrt(105) = rrt(105) * density(11) * density(20) 
  rrt(106) = rrt(106) * density(31) * density(38) 
  rrt(107) = rrt(107) * density(17) * density(38) 
  rrt(108) = rrt(108) * density(04) * density(32) 
  rrt(109) = rrt(109) * density(04) * density(33) 
  rrt(110) = rrt(110) * density(04) * density(34) 
  rrt(111) = rrt(111) * density(06) * density(07) 
  rrt(112) = rrt(112) * density(06) * density(09) 
  rrt(113) = rrt(113) * density(06) * density(08) 
  rrt(114) = rrt(114) * density(38)**2 
  ydot(01) = +rrt(025)+rrt(026)+rrt(027)+rrt(028)+rrt(029)+rrt(030)+rrt(031)+rrt(032)+rrt(033)+rrt(034)+rrt(035)+rrt(036)+rrt(037)&
             +rrt(038)+rrt(039)+rrt(040)+rrt(041)+rrt(042)+rrt(043)-rrt(044)-rrt(045)-rrt(046)-rrt(047)-rrt(048)-rrt(049)-rrt(050)&
             -rrt(060)-rrt(064)-rrt(067) 
  ydot(02) = +rrt(077) 
  ydot(03) = +rrt(074)+rrt(075)+rrt(076)-rrt(077)-rrt(083)-rrt(084)-rrt(085)+rrt(088) 
  ydot(04) = +rrt(054)+rrt(055)+rrt(056)-rrt(088)-rrt(089)-rrt(098)-rrt(108)-rrt(109)-rrt(110) 
  ydot(05) = -rrt(025)+rrt(044)+rrt(051)+rrt(052)+rrt(053)-rrt(078)-  2.d0 * rrt(079)-rrt(089)+  2.d0 * rrt(090)+rrt(091)+rrt(092)&
             +rrt(093)-rrt(094)-rrt(095) 
  ydot(06) = +rrt(025)-rrt(044)-rrt(111)-rrt(112)-rrt(113) 
  ydot(07) = -rrt(001)-rrt(002)+rrt(003)+rrt(004)-rrt(026)+rrt(045)-rrt(051)-rrt(054)+rrt(065)+rrt(066)-rrt(074)+rrt(078)-rrt(083)&
             -rrt(091)+rrt(100)+rrt(101)+rrt(102)+rrt(103)+rrt(104)+rrt(105)-rrt(111) 
  ydot(08) = +rrt(002)-rrt(004)-rrt(028)-rrt(053)-rrt(056)-rrt(076)-rrt(085)-rrt(093)-rrt(113) 
  ydot(09) = +rrt(001)-rrt(003)-rrt(027)-rrt(052)-rrt(055)-rrt(075)-rrt(084)-rrt(092)-rrt(112) 
  ydot(10) = +rrt(026)+rrt(027)+rrt(028)-rrt(045)-rrt(091)-rrt(092)-rrt(093) 
  ydot(11) = +rrt(091)+rrt(092)+rrt(093)-rrt(100)-rrt(101)-rrt(102)-rrt(103)-rrt(104)-rrt(105) 
  ydot(12) = -rrt(013)-rrt(014)-rrt(015)+rrt(016)+rrt(017)+rrt(018)-rrt(035)+rrt(048)+rrt(061)+rrt(062)+rrt(063)+rrt(064)+rrt(065)&
             +rrt(066)+rrt(067)+rrt(107) 
  ydot(13) = +rrt(015)-rrt(018)-rrt(038) 
  ydot(14) = +rrt(013)-rrt(016)-rrt(036) 
  ydot(15) = +rrt(014)-rrt(017)-rrt(037) 
  ydot(16) = +rrt(035)+rrt(036)+rrt(037)+rrt(038)-rrt(048) 
  ydot(17) = +rrt(060)+rrt(068)+rrt(069)+rrt(070)-rrt(095)-rrt(107) 
  ydot(18) = -rrt(009)-rrt(010)+rrt(011)+rrt(012)-rrt(032)+rrt(047)+rrt(057)+rrt(058)+rrt(059)-rrt(061)-rrt(068)-rrt(080)+rrt(083)&
             +rrt(084)+rrt(085)+rrt(089)-rrt(103) 
  ydot(19) = +rrt(010)-rrt(012)-rrt(034)-rrt(062)-rrt(069)-rrt(081)-rrt(104) 
  ydot(20) = +rrt(009)-rrt(011)-rrt(033)-rrt(063)-rrt(070)-rrt(082)-rrt(105) 
  ydot(21) = +rrt(032)+rrt(033)+rrt(034)-rrt(047) 
  ydot(22) = +rrt(080)+rrt(081)+rrt(082)-rrt(086)-  2.d0 * rrt(087)-rrt(090)-rrt(094) 
  ydot(23) = -rrt(060)-rrt(064)-rrt(067)+rrt(100)+rrt(101)+rrt(102)+rrt(103)+rrt(104)+rrt(105)+rrt(111)+rrt(112)+rrt(113) 
  ydot(24) = -rrt(005)-rrt(006)+rrt(007)+rrt(008)-rrt(029)+rrt(046)-rrt(057)+rrt(079)+rrt(086)-rrt(100) 
  ydot(25) = +rrt(006)-rrt(008)-rrt(031)-rrt(059)-rrt(102) 
  ydot(26) = +rrt(005)-rrt(007)-rrt(030)-rrt(058)-rrt(101) 
  ydot(27) = +rrt(029)+rrt(030)+rrt(031)-rrt(046) 
  ydot(28) = -rrt(023)+rrt(024)-rrt(042)+rrt(050)-rrt(065)+rrt(071)+rrt(072)+rrt(073)+rrt(095)-rrt(096) 
  ydot(29) = +rrt(023)-rrt(024)-rrt(043)-rrt(066)-rrt(097) 
  ydot(30) = +rrt(042)+rrt(043)-rrt(050) 
  ydot(31) = +rrt(096)+rrt(097)-rrt(099)-rrt(106) 
  ydot(32) = -rrt(019)-rrt(020)+rrt(021)+rrt(022)-rrt(039)+rrt(049)-rrt(071)+rrt(094)+rrt(099)+rrt(106)-rrt(108) 
  ydot(33) = +rrt(020)-rrt(022)-rrt(041)-rrt(072)-rrt(109) 
  ydot(34) = +rrt(019)-rrt(021)-rrt(040)-rrt(073)-rrt(110) 
  ydot(35) = +rrt(039)+rrt(040)+rrt(041)-rrt(049) 
  ydot(36) = +rrt(087)-rrt(098)+rrt(108)+rrt(109)+rrt(110) 
  ydot(37) = +rrt(098) 
  ydot(38) = +rrt(051)+rrt(052)+rrt(053)+  2.d0 * rrt(060)+rrt(064)+  3.d0 * rrt(067)+rrt(068)+rrt(069)+rrt(070)+rrt(074)+rrt(075)&
             +rrt(076)+rrt(077)-rrt(078)-rrt(080)-rrt(081)-rrt(082)+rrt(083)+rrt(084)+rrt(085)-rrt(086)-rrt(088)+rrt(089)-rrt(090)&
             -rrt(096)-rrt(097)+rrt(099)-rrt(106)-rrt(107)-  2.d0 * rrt(114) 
  ydot(39) = +rrt(054)+rrt(055)+rrt(056)+rrt(057)+rrt(058)+rrt(059)+rrt(061)+rrt(062)+rrt(063)+rrt(064)+rrt(071)+rrt(072)+rrt(073)&
             +rrt(074)+rrt(075)+rrt(076)+rrt(088)-rrt(099)+rrt(100)+rrt(101)+rrt(102)+rrt(107)+rrt(111)+rrt(112)+rrt(113)+rrt(114) 
  if( ldensity_constant ) where( density_constant(:) ) ydot(1:species_max) = 0.0d0
  ydot(40) = 0.0d0
  if( lgas_heating ) then
    ydot(40) = ( ZDPlasKin_cfg(14)/k_B + ydot(40) ) / ( sum(density(1:species_max)) - density(species_electrons) ) &
            + eV_to_K * ZDPlasKin_cfg(11) * density(species_electrons)
    ydot(40) = ydot(40) * ZDPlasKin_cfg(13)
  endif
  return
end subroutine ZDPlasKin_fex
!-----------------------------------------------------------------------------------------------------------------------------------
!
! reaction jacobian
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_jex(neq,t,y,ml,mu,pd,nrpd)
  implicit none
  integer,          intent(in)  :: neq, ml, mu, nrpd
  double precision, intent(in)  :: t, y(neq)
  double precision, intent(out) :: pd(nrpd,neq)
  integer                       :: i
  if( lgas_heating ) ZDPlasKin_cfg(1) = y(40)
  density(:) = y(1:species_max)
  call ZDPlasKin_reac_rates(t)
  pd(07,01) = pd(07,01) - rrt(001) * density(07) 
  pd(07,07) = pd(07,07) - rrt(001) * density(01) 
  pd(09,01) = pd(09,01) + rrt(001) * density(07) 
  pd(09,07) = pd(09,07) + rrt(001) * density(01) 
  pd(07,01) = pd(07,01) - rrt(002) * density(07) 
  pd(07,07) = pd(07,07) - rrt(002) * density(01) 
  pd(08,01) = pd(08,01) + rrt(002) * density(07) 
  pd(08,07) = pd(08,07) + rrt(002) * density(01) 
  pd(07,01) = pd(07,01) + rrt(003) * density(09) 
  pd(07,09) = pd(07,09) + rrt(003) * density(01) 
  pd(09,01) = pd(09,01) - rrt(003) * density(09) 
  pd(09,09) = pd(09,09) - rrt(003) * density(01) 
  pd(07,01) = pd(07,01) + rrt(004) * density(08) 
  pd(07,08) = pd(07,08) + rrt(004) * density(01) 
  pd(08,01) = pd(08,01) - rrt(004) * density(08) 
  pd(08,08) = pd(08,08) - rrt(004) * density(01) 
  pd(24,01) = pd(24,01) - rrt(005) * density(24) 
  pd(24,24) = pd(24,24) - rrt(005) * density(01) 
  pd(26,01) = pd(26,01) + rrt(005) * density(24) 
  pd(26,24) = pd(26,24) + rrt(005) * density(01) 
  pd(24,01) = pd(24,01) - rrt(006) * density(24) 
  pd(24,24) = pd(24,24) - rrt(006) * density(01) 
  pd(25,01) = pd(25,01) + rrt(006) * density(24) 
  pd(25,24) = pd(25,24) + rrt(006) * density(01) 
  pd(24,01) = pd(24,01) + rrt(007) * density(26) 
  pd(24,26) = pd(24,26) + rrt(007) * density(01) 
  pd(26,01) = pd(26,01) - rrt(007) * density(26) 
  pd(26,26) = pd(26,26) - rrt(007) * density(01) 
  pd(24,01) = pd(24,01) + rrt(008) * density(25) 
  pd(24,25) = pd(24,25) + rrt(008) * density(01) 
  pd(25,01) = pd(25,01) - rrt(008) * density(25) 
  pd(25,25) = pd(25,25) - rrt(008) * density(01) 
  pd(18,01) = pd(18,01) - rrt(009) * density(18) 
  pd(18,18) = pd(18,18) - rrt(009) * density(01) 
  pd(20,01) = pd(20,01) + rrt(009) * density(18) 
  pd(20,18) = pd(20,18) + rrt(009) * density(01) 
  pd(18,01) = pd(18,01) - rrt(010) * density(18) 
  pd(18,18) = pd(18,18) - rrt(010) * density(01) 
  pd(19,01) = pd(19,01) + rrt(010) * density(18) 
  pd(19,18) = pd(19,18) + rrt(010) * density(01) 
  pd(18,01) = pd(18,01) + rrt(011) * density(20) 
  pd(18,20) = pd(18,20) + rrt(011) * density(01) 
  pd(20,01) = pd(20,01) - rrt(011) * density(20) 
  pd(20,20) = pd(20,20) - rrt(011) * density(01) 
  pd(18,01) = pd(18,01) + rrt(012) * density(19) 
  pd(18,19) = pd(18,19) + rrt(012) * density(01) 
  pd(19,01) = pd(19,01) - rrt(012) * density(19) 
  pd(19,19) = pd(19,19) - rrt(012) * density(01) 
  pd(12,01) = pd(12,01) - rrt(013) * density(12) 
  pd(12,12) = pd(12,12) - rrt(013) * density(01) 
  pd(14,01) = pd(14,01) + rrt(013) * density(12) 
  pd(14,12) = pd(14,12) + rrt(013) * density(01) 
  pd(12,01) = pd(12,01) - rrt(014) * density(12) 
  pd(12,12) = pd(12,12) - rrt(014) * density(01) 
  pd(15,01) = pd(15,01) + rrt(014) * density(12) 
  pd(15,12) = pd(15,12) + rrt(014) * density(01) 
  pd(12,01) = pd(12,01) - rrt(015) * density(12) 
  pd(12,12) = pd(12,12) - rrt(015) * density(01) 
  pd(13,01) = pd(13,01) + rrt(015) * density(12) 
  pd(13,12) = pd(13,12) + rrt(015) * density(01) 
  pd(12,01) = pd(12,01) + rrt(016) * density(14) 
  pd(12,14) = pd(12,14) + rrt(016) * density(01) 
  pd(14,01) = pd(14,01) - rrt(016) * density(14) 
  pd(14,14) = pd(14,14) - rrt(016) * density(01) 
  pd(12,01) = pd(12,01) + rrt(017) * density(15) 
  pd(12,15) = pd(12,15) + rrt(017) * density(01) 
  pd(15,01) = pd(15,01) - rrt(017) * density(15) 
  pd(15,15) = pd(15,15) - rrt(017) * density(01) 
  pd(12,01) = pd(12,01) + rrt(018) * density(13) 
  pd(12,13) = pd(12,13) + rrt(018) * density(01) 
  pd(13,01) = pd(13,01) - rrt(018) * density(13) 
  pd(13,13) = pd(13,13) - rrt(018) * density(01) 
  pd(32,01) = pd(32,01) - rrt(019) * density(32) 
  pd(32,32) = pd(32,32) - rrt(019) * density(01) 
  pd(34,01) = pd(34,01) + rrt(019) * density(32) 
  pd(34,32) = pd(34,32) + rrt(019) * density(01) 
  pd(32,01) = pd(32,01) - rrt(020) * density(32) 
  pd(32,32) = pd(32,32) - rrt(020) * density(01) 
  pd(33,01) = pd(33,01) + rrt(020) * density(32) 
  pd(33,32) = pd(33,32) + rrt(020) * density(01) 
  pd(32,01) = pd(32,01) + rrt(021) * density(34) 
  pd(32,34) = pd(32,34) + rrt(021) * density(01) 
  pd(34,01) = pd(34,01) - rrt(021) * density(34) 
  pd(34,34) = pd(34,34) - rrt(021) * density(01) 
  pd(32,01) = pd(32,01) + rrt(022) * density(33) 
  pd(32,33) = pd(32,33) + rrt(022) * density(01) 
  pd(33,01) = pd(33,01) - rrt(022) * density(33) 
  pd(33,33) = pd(33,33) - rrt(022) * density(01) 
  pd(28,01) = pd(28,01) - rrt(023) * density(28) 
  pd(28,28) = pd(28,28) - rrt(023) * density(01) 
  pd(29,01) = pd(29,01) + rrt(023) * density(28) 
  pd(29,28) = pd(29,28) + rrt(023) * density(01) 
  pd(28,01) = pd(28,01) + rrt(024) * density(29) 
  pd(28,29) = pd(28,29) + rrt(024) * density(01) 
  pd(29,01) = pd(29,01) - rrt(024) * density(29) 
  pd(29,29) = pd(29,29) - rrt(024) * density(01) 
  pd(01,01) = pd(01,01) + rrt(025) * density(05) 
  pd(01,05) = pd(01,05) + rrt(025) * density(01) 
  pd(05,01) = pd(05,01) - rrt(025) * density(05) 
  pd(05,05) = pd(05,05) - rrt(025) * density(01) 
  pd(06,01) = pd(06,01) + rrt(025) * density(05) 
  pd(06,05) = pd(06,05) + rrt(025) * density(01) 
  pd(01,01) = pd(01,01) + rrt(026) * density(07) 
  pd(01,07) = pd(01,07) + rrt(026) * density(01) 
  pd(07,01) = pd(07,01) - rrt(026) * density(07) 
  pd(07,07) = pd(07,07) - rrt(026) * density(01) 
  pd(10,01) = pd(10,01) + rrt(026) * density(07) 
  pd(10,07) = pd(10,07) + rrt(026) * density(01) 
  pd(01,01) = pd(01,01) + rrt(027) * density(09) 
  pd(01,09) = pd(01,09) + rrt(027) * density(01) 
  pd(09,01) = pd(09,01) - rrt(027) * density(09) 
  pd(09,09) = pd(09,09) - rrt(027) * density(01) 
  pd(10,01) = pd(10,01) + rrt(027) * density(09) 
  pd(10,09) = pd(10,09) + rrt(027) * density(01) 
  pd(01,01) = pd(01,01) + rrt(028) * density(08) 
  pd(01,08) = pd(01,08) + rrt(028) * density(01) 
  pd(08,01) = pd(08,01) - rrt(028) * density(08) 
  pd(08,08) = pd(08,08) - rrt(028) * density(01) 
  pd(10,01) = pd(10,01) + rrt(028) * density(08) 
  pd(10,08) = pd(10,08) + rrt(028) * density(01) 
  pd(01,01) = pd(01,01) + rrt(029) * density(24) 
  pd(01,24) = pd(01,24) + rrt(029) * density(01) 
  pd(24,01) = pd(24,01) - rrt(029) * density(24) 
  pd(24,24) = pd(24,24) - rrt(029) * density(01) 
  pd(27,01) = pd(27,01) + rrt(029) * density(24) 
  pd(27,24) = pd(27,24) + rrt(029) * density(01) 
  pd(01,01) = pd(01,01) + rrt(030) * density(26) 
  pd(01,26) = pd(01,26) + rrt(030) * density(01) 
  pd(26,01) = pd(26,01) - rrt(030) * density(26) 
  pd(26,26) = pd(26,26) - rrt(030) * density(01) 
  pd(27,01) = pd(27,01) + rrt(030) * density(26) 
  pd(27,26) = pd(27,26) + rrt(030) * density(01) 
  pd(01,01) = pd(01,01) + rrt(031) * density(25) 
  pd(01,25) = pd(01,25) + rrt(031) * density(01) 
  pd(25,01) = pd(25,01) - rrt(031) * density(25) 
  pd(25,25) = pd(25,25) - rrt(031) * density(01) 
  pd(27,01) = pd(27,01) + rrt(031) * density(25) 
  pd(27,25) = pd(27,25) + rrt(031) * density(01) 
  pd(01,01) = pd(01,01) + rrt(032) * density(18) 
  pd(01,18) = pd(01,18) + rrt(032) * density(01) 
  pd(18,01) = pd(18,01) - rrt(032) * density(18) 
  pd(18,18) = pd(18,18) - rrt(032) * density(01) 
  pd(21,01) = pd(21,01) + rrt(032) * density(18) 
  pd(21,18) = pd(21,18) + rrt(032) * density(01) 
  pd(01,01) = pd(01,01) + rrt(033) * density(20) 
  pd(01,20) = pd(01,20) + rrt(033) * density(01) 
  pd(20,01) = pd(20,01) - rrt(033) * density(20) 
  pd(20,20) = pd(20,20) - rrt(033) * density(01) 
  pd(21,01) = pd(21,01) + rrt(033) * density(20) 
  pd(21,20) = pd(21,20) + rrt(033) * density(01) 
  pd(01,01) = pd(01,01) + rrt(034) * density(19) 
  pd(01,19) = pd(01,19) + rrt(034) * density(01) 
  pd(19,01) = pd(19,01) - rrt(034) * density(19) 
  pd(19,19) = pd(19,19) - rrt(034) * density(01) 
  pd(21,01) = pd(21,01) + rrt(034) * density(19) 
  pd(21,19) = pd(21,19) + rrt(034) * density(01) 
  pd(01,01) = pd(01,01) + rrt(035) * density(12) 
  pd(01,12) = pd(01,12) + rrt(035) * density(01) 
  pd(12,01) = pd(12,01) - rrt(035) * density(12) 
  pd(12,12) = pd(12,12) - rrt(035) * density(01) 
  pd(16,01) = pd(16,01) + rrt(035) * density(12) 
  pd(16,12) = pd(16,12) + rrt(035) * density(01) 
  pd(01,01) = pd(01,01) + rrt(036) * density(14) 
  pd(01,14) = pd(01,14) + rrt(036) * density(01) 
  pd(14,01) = pd(14,01) - rrt(036) * density(14) 
  pd(14,14) = pd(14,14) - rrt(036) * density(01) 
  pd(16,01) = pd(16,01) + rrt(036) * density(14) 
  pd(16,14) = pd(16,14) + rrt(036) * density(01) 
  pd(01,01) = pd(01,01) + rrt(037) * density(15) 
  pd(01,15) = pd(01,15) + rrt(037) * density(01) 
  pd(15,01) = pd(15,01) - rrt(037) * density(15) 
  pd(15,15) = pd(15,15) - rrt(037) * density(01) 
  pd(16,01) = pd(16,01) + rrt(037) * density(15) 
  pd(16,15) = pd(16,15) + rrt(037) * density(01) 
  pd(01,01) = pd(01,01) + rrt(038) * density(13) 
  pd(01,13) = pd(01,13) + rrt(038) * density(01) 
  pd(13,01) = pd(13,01) - rrt(038) * density(13) 
  pd(13,13) = pd(13,13) - rrt(038) * density(01) 
  pd(16,01) = pd(16,01) + rrt(038) * density(13) 
  pd(16,13) = pd(16,13) + rrt(038) * density(01) 
  pd(01,01) = pd(01,01) + rrt(039) * density(32) 
  pd(01,32) = pd(01,32) + rrt(039) * density(01) 
  pd(32,01) = pd(32,01) - rrt(039) * density(32) 
  pd(32,32) = pd(32,32) - rrt(039) * density(01) 
  pd(35,01) = pd(35,01) + rrt(039) * density(32) 
  pd(35,32) = pd(35,32) + rrt(039) * density(01) 
  pd(01,01) = pd(01,01) + rrt(040) * density(34) 
  pd(01,34) = pd(01,34) + rrt(040) * density(01) 
  pd(34,01) = pd(34,01) - rrt(040) * density(34) 
  pd(34,34) = pd(34,34) - rrt(040) * density(01) 
  pd(35,01) = pd(35,01) + rrt(040) * density(34) 
  pd(35,34) = pd(35,34) + rrt(040) * density(01) 
  pd(01,01) = pd(01,01) + rrt(041) * density(33) 
  pd(01,33) = pd(01,33) + rrt(041) * density(01) 
  pd(33,01) = pd(33,01) - rrt(041) * density(33) 
  pd(33,33) = pd(33,33) - rrt(041) * density(01) 
  pd(35,01) = pd(35,01) + rrt(041) * density(33) 
  pd(35,33) = pd(35,33) + rrt(041) * density(01) 
  pd(01,01) = pd(01,01) + rrt(042) * density(28) 
  pd(01,28) = pd(01,28) + rrt(042) * density(01) 
  pd(28,01) = pd(28,01) - rrt(042) * density(28) 
  pd(28,28) = pd(28,28) - rrt(042) * density(01) 
  pd(30,01) = pd(30,01) + rrt(042) * density(28) 
  pd(30,28) = pd(30,28) + rrt(042) * density(01) 
  pd(01,01) = pd(01,01) + rrt(043) * density(29) 
  pd(01,29) = pd(01,29) + rrt(043) * density(01) 
  pd(29,01) = pd(29,01) - rrt(043) * density(29) 
  pd(29,29) = pd(29,29) - rrt(043) * density(01) 
  pd(30,01) = pd(30,01) + rrt(043) * density(29) 
  pd(30,29) = pd(30,29) + rrt(043) * density(01) 
  pd(01,01) = pd(01,01) - rrt(044) * density(01) * density(06) * 2.0d0
  pd(01,06) = pd(01,06) - rrt(044) * density(01)**2 
  pd(05,01) = pd(05,01) + rrt(044) * density(01) * density(06) * 2.0d0
  pd(05,06) = pd(05,06) + rrt(044) * density(01)**2 
  pd(06,01) = pd(06,01) - rrt(044) * density(01) * density(06) * 2.0d0
  pd(06,06) = pd(06,06) - rrt(044) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(045) * density(01) * density(10) * 2.0d0
  pd(01,10) = pd(01,10) - rrt(045) * density(01)**2 
  pd(07,01) = pd(07,01) + rrt(045) * density(01) * density(10) * 2.0d0
  pd(07,10) = pd(07,10) + rrt(045) * density(01)**2 
  pd(10,01) = pd(10,01) - rrt(045) * density(01) * density(10) * 2.0d0
  pd(10,10) = pd(10,10) - rrt(045) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(046) * density(01) * density(27) * 2.0d0
  pd(01,27) = pd(01,27) - rrt(046) * density(01)**2 
  pd(24,01) = pd(24,01) + rrt(046) * density(01) * density(27) * 2.0d0
  pd(24,27) = pd(24,27) + rrt(046) * density(01)**2 
  pd(27,01) = pd(27,01) - rrt(046) * density(01) * density(27) * 2.0d0
  pd(27,27) = pd(27,27) - rrt(046) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(047) * density(01) * density(21) * 2.0d0
  pd(01,21) = pd(01,21) - rrt(047) * density(01)**2 
  pd(18,01) = pd(18,01) + rrt(047) * density(01) * density(21) * 2.0d0
  pd(18,21) = pd(18,21) + rrt(047) * density(01)**2 
  pd(21,01) = pd(21,01) - rrt(047) * density(01) * density(21) * 2.0d0
  pd(21,21) = pd(21,21) - rrt(047) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(048) * density(01) * density(16) * 2.0d0
  pd(01,16) = pd(01,16) - rrt(048) * density(01)**2 
  pd(12,01) = pd(12,01) + rrt(048) * density(01) * density(16) * 2.0d0
  pd(12,16) = pd(12,16) + rrt(048) * density(01)**2 
  pd(16,01) = pd(16,01) - rrt(048) * density(01) * density(16) * 2.0d0
  pd(16,16) = pd(16,16) - rrt(048) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(049) * density(01) * density(35) * 2.0d0
  pd(01,35) = pd(01,35) - rrt(049) * density(01)**2 
  pd(32,01) = pd(32,01) + rrt(049) * density(01) * density(35) * 2.0d0
  pd(32,35) = pd(32,35) + rrt(049) * density(01)**2 
  pd(35,01) = pd(35,01) - rrt(049) * density(01) * density(35) * 2.0d0
  pd(35,35) = pd(35,35) - rrt(049) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(050) * density(01) * density(30) * 2.0d0
  pd(01,30) = pd(01,30) - rrt(050) * density(01)**2 
  pd(28,01) = pd(28,01) + rrt(050) * density(01) * density(30) * 2.0d0
  pd(28,30) = pd(28,30) + rrt(050) * density(01)**2 
  pd(30,01) = pd(30,01) - rrt(050) * density(01) * density(30) * 2.0d0
  pd(30,30) = pd(30,30) - rrt(050) * density(01)**2 
  pd(05,01) = pd(05,01) + rrt(051) * density(07) 
  pd(05,07) = pd(05,07) + rrt(051) * density(01) 
  pd(07,01) = pd(07,01) - rrt(051) * density(07) 
  pd(07,07) = pd(07,07) - rrt(051) * density(01) 
  pd(38,01) = pd(38,01) + rrt(051) * density(07) 
  pd(38,07) = pd(38,07) + rrt(051) * density(01) 
  pd(05,01) = pd(05,01) + rrt(052) * density(09) 
  pd(05,09) = pd(05,09) + rrt(052) * density(01) 
  pd(09,01) = pd(09,01) - rrt(052) * density(09) 
  pd(09,09) = pd(09,09) - rrt(052) * density(01) 
  pd(38,01) = pd(38,01) + rrt(052) * density(09) 
  pd(38,09) = pd(38,09) + rrt(052) * density(01) 
  pd(05,01) = pd(05,01) + rrt(053) * density(08) 
  pd(05,08) = pd(05,08) + rrt(053) * density(01) 
  pd(08,01) = pd(08,01) - rrt(053) * density(08) 
  pd(08,08) = pd(08,08) - rrt(053) * density(01) 
  pd(38,01) = pd(38,01) + rrt(053) * density(08) 
  pd(38,08) = pd(38,08) + rrt(053) * density(01) 
  pd(04,01) = pd(04,01) + rrt(054) * density(07) 
  pd(04,07) = pd(04,07) + rrt(054) * density(01) 
  pd(07,01) = pd(07,01) - rrt(054) * density(07) 
  pd(07,07) = pd(07,07) - rrt(054) * density(01) 
  pd(39,01) = pd(39,01) + rrt(054) * density(07) 
  pd(39,07) = pd(39,07) + rrt(054) * density(01) 
  pd(04,01) = pd(04,01) + rrt(055) * density(09) 
  pd(04,09) = pd(04,09) + rrt(055) * density(01) 
  pd(09,01) = pd(09,01) - rrt(055) * density(09) 
  pd(09,09) = pd(09,09) - rrt(055) * density(01) 
  pd(39,01) = pd(39,01) + rrt(055) * density(09) 
  pd(39,09) = pd(39,09) + rrt(055) * density(01) 
  pd(04,01) = pd(04,01) + rrt(056) * density(08) 
  pd(04,08) = pd(04,08) + rrt(056) * density(01) 
  pd(08,01) = pd(08,01) - rrt(056) * density(08) 
  pd(08,08) = pd(08,08) - rrt(056) * density(01) 
  pd(39,01) = pd(39,01) + rrt(056) * density(08) 
  pd(39,08) = pd(39,08) + rrt(056) * density(01) 
  pd(18,01) = pd(18,01) + rrt(057) * density(24) 
  pd(18,24) = pd(18,24) + rrt(057) * density(01) 
  pd(24,01) = pd(24,01) - rrt(057) * density(24) 
  pd(24,24) = pd(24,24) - rrt(057) * density(01) 
  pd(39,01) = pd(39,01) + rrt(057) * density(24) 
  pd(39,24) = pd(39,24) + rrt(057) * density(01) 
  pd(18,01) = pd(18,01) + rrt(058) * density(26) 
  pd(18,26) = pd(18,26) + rrt(058) * density(01) 
  pd(26,01) = pd(26,01) - rrt(058) * density(26) 
  pd(26,26) = pd(26,26) - rrt(058) * density(01) 
  pd(39,01) = pd(39,01) + rrt(058) * density(26) 
  pd(39,26) = pd(39,26) + rrt(058) * density(01) 
  pd(18,01) = pd(18,01) + rrt(059) * density(25) 
  pd(18,25) = pd(18,25) + rrt(059) * density(01) 
  pd(25,01) = pd(25,01) - rrt(059) * density(25) 
  pd(25,25) = pd(25,25) - rrt(059) * density(01) 
  pd(39,01) = pd(39,01) + rrt(059) * density(25) 
  pd(39,25) = pd(39,25) + rrt(059) * density(01) 
  pd(01,01) = pd(01,01) - rrt(060) * density(23) 
  pd(01,23) = pd(01,23) - rrt(060) * density(01) 
  pd(17,01) = pd(17,01) + rrt(060) * density(23) 
  pd(17,23) = pd(17,23) + rrt(060) * density(01) 
  pd(23,01) = pd(23,01) - rrt(060) * density(23) 
  pd(23,23) = pd(23,23) - rrt(060) * density(01) 
  pd(38,01) = pd(38,01) + rrt(060) * density(23) * 2.0d0
  pd(38,23) = pd(38,23) + rrt(060) * density(01) * 2.0d0
  pd(12,01) = pd(12,01) + rrt(061) * density(18) 
  pd(12,18) = pd(12,18) + rrt(061) * density(01) 
  pd(18,01) = pd(18,01) - rrt(061) * density(18) 
  pd(18,18) = pd(18,18) - rrt(061) * density(01) 
  pd(39,01) = pd(39,01) + rrt(061) * density(18) 
  pd(39,18) = pd(39,18) + rrt(061) * density(01) 
  pd(12,01) = pd(12,01) + rrt(062) * density(19) 
  pd(12,19) = pd(12,19) + rrt(062) * density(01) 
  pd(19,01) = pd(19,01) - rrt(062) * density(19) 
  pd(19,19) = pd(19,19) - rrt(062) * density(01) 
  pd(39,01) = pd(39,01) + rrt(062) * density(19) 
  pd(39,19) = pd(39,19) + rrt(062) * density(01) 
  pd(12,01) = pd(12,01) + rrt(063) * density(20) 
  pd(12,20) = pd(12,20) + rrt(063) * density(01) 
  pd(20,01) = pd(20,01) - rrt(063) * density(20) 
  pd(20,20) = pd(20,20) - rrt(063) * density(01) 
  pd(39,01) = pd(39,01) + rrt(063) * density(20) 
  pd(39,20) = pd(39,20) + rrt(063) * density(01) 
  pd(01,01) = pd(01,01) - rrt(064) * density(23) 
  pd(01,23) = pd(01,23) - rrt(064) * density(01) 
  pd(12,01) = pd(12,01) + rrt(064) * density(23) 
  pd(12,23) = pd(12,23) + rrt(064) * density(01) 
  pd(23,01) = pd(23,01) - rrt(064) * density(23) 
  pd(23,23) = pd(23,23) - rrt(064) * density(01) 
  pd(38,01) = pd(38,01) + rrt(064) * density(23) 
  pd(38,23) = pd(38,23) + rrt(064) * density(01) 
  pd(39,01) = pd(39,01) + rrt(064) * density(23) 
  pd(39,23) = pd(39,23) + rrt(064) * density(01) 
  pd(07,01) = pd(07,01) + rrt(065) * density(28) 
  pd(07,28) = pd(07,28) + rrt(065) * density(01) 
  pd(12,01) = pd(12,01) + rrt(065) * density(28) 
  pd(12,28) = pd(12,28) + rrt(065) * density(01) 
  pd(28,01) = pd(28,01) - rrt(065) * density(28) 
  pd(28,28) = pd(28,28) - rrt(065) * density(01) 
  pd(07,01) = pd(07,01) + rrt(066) * density(29) 
  pd(07,29) = pd(07,29) + rrt(066) * density(01) 
  pd(12,01) = pd(12,01) + rrt(066) * density(29) 
  pd(12,29) = pd(12,29) + rrt(066) * density(01) 
  pd(29,01) = pd(29,01) - rrt(066) * density(29) 
  pd(29,29) = pd(29,29) - rrt(066) * density(01) 
  pd(01,01) = pd(01,01) - rrt(067) * density(23) 
  pd(01,23) = pd(01,23) - rrt(067) * density(01) 
  pd(12,01) = pd(12,01) + rrt(067) * density(23) 
  pd(12,23) = pd(12,23) + rrt(067) * density(01) 
  pd(23,01) = pd(23,01) - rrt(067) * density(23) 
  pd(23,23) = pd(23,23) - rrt(067) * density(01) 
  pd(38,01) = pd(38,01) + rrt(067) * density(23) * 3.0d0
  pd(38,23) = pd(38,23) + rrt(067) * density(01) * 3.0d0
  pd(17,01) = pd(17,01) + rrt(068) * density(18) 
  pd(17,18) = pd(17,18) + rrt(068) * density(01) 
  pd(18,01) = pd(18,01) - rrt(068) * density(18) 
  pd(18,18) = pd(18,18) - rrt(068) * density(01) 
  pd(38,01) = pd(38,01) + rrt(068) * density(18) 
  pd(38,18) = pd(38,18) + rrt(068) * density(01) 
  pd(17,01) = pd(17,01) + rrt(069) * density(19) 
  pd(17,19) = pd(17,19) + rrt(069) * density(01) 
  pd(19,01) = pd(19,01) - rrt(069) * density(19) 
  pd(19,19) = pd(19,19) - rrt(069) * density(01) 
  pd(38,01) = pd(38,01) + rrt(069) * density(19) 
  pd(38,19) = pd(38,19) + rrt(069) * density(01) 
  pd(17,01) = pd(17,01) + rrt(070) * density(20) 
  pd(17,20) = pd(17,20) + rrt(070) * density(01) 
  pd(20,01) = pd(20,01) - rrt(070) * density(20) 
  pd(20,20) = pd(20,20) - rrt(070) * density(01) 
  pd(38,01) = pd(38,01) + rrt(070) * density(20) 
  pd(38,20) = pd(38,20) + rrt(070) * density(01) 
  pd(28,01) = pd(28,01) + rrt(071) * density(32) 
  pd(28,32) = pd(28,32) + rrt(071) * density(01) 
  pd(32,01) = pd(32,01) - rrt(071) * density(32) 
  pd(32,32) = pd(32,32) - rrt(071) * density(01) 
  pd(39,01) = pd(39,01) + rrt(071) * density(32) 
  pd(39,32) = pd(39,32) + rrt(071) * density(01) 
  pd(28,01) = pd(28,01) + rrt(072) * density(33) 
  pd(28,33) = pd(28,33) + rrt(072) * density(01) 
  pd(33,01) = pd(33,01) - rrt(072) * density(33) 
  pd(33,33) = pd(33,33) - rrt(072) * density(01) 
  pd(39,01) = pd(39,01) + rrt(072) * density(33) 
  pd(39,33) = pd(39,33) + rrt(072) * density(01) 
  pd(28,01) = pd(28,01) + rrt(073) * density(34) 
  pd(28,34) = pd(28,34) + rrt(073) * density(01) 
  pd(34,01) = pd(34,01) - rrt(073) * density(34) 
  pd(34,34) = pd(34,34) - rrt(073) * density(01) 
  pd(39,01) = pd(39,01) + rrt(073) * density(34) 
  pd(39,34) = pd(39,34) + rrt(073) * density(01) 
  pd(03,01) = pd(03,01) + rrt(074) * density(07) 
  pd(03,07) = pd(03,07) + rrt(074) * density(01) 
  pd(07,01) = pd(07,01) - rrt(074) * density(07) 
  pd(07,07) = pd(07,07) - rrt(074) * density(01) 
  pd(38,01) = pd(38,01) + rrt(074) * density(07) 
  pd(38,07) = pd(38,07) + rrt(074) * density(01) 
  pd(39,01) = pd(39,01) + rrt(074) * density(07) 
  pd(39,07) = pd(39,07) + rrt(074) * density(01) 
  pd(03,01) = pd(03,01) + rrt(075) * density(09) 
  pd(03,09) = pd(03,09) + rrt(075) * density(01) 
  pd(09,01) = pd(09,01) - rrt(075) * density(09) 
  pd(09,09) = pd(09,09) - rrt(075) * density(01) 
  pd(38,01) = pd(38,01) + rrt(075) * density(09) 
  pd(38,09) = pd(38,09) + rrt(075) * density(01) 
  pd(39,01) = pd(39,01) + rrt(075) * density(09) 
  pd(39,09) = pd(39,09) + rrt(075) * density(01) 
  pd(03,01) = pd(03,01) + rrt(076) * density(08) 
  pd(03,08) = pd(03,08) + rrt(076) * density(01) 
  pd(08,01) = pd(08,01) - rrt(076) * density(08) 
  pd(08,08) = pd(08,08) - rrt(076) * density(01) 
  pd(38,01) = pd(38,01) + rrt(076) * density(08) 
  pd(38,08) = pd(38,08) + rrt(076) * density(01) 
  pd(39,01) = pd(39,01) + rrt(076) * density(08) 
  pd(39,08) = pd(39,08) + rrt(076) * density(01) 
  pd(02,01) = pd(02,01) + rrt(077) * density(03) 
  pd(02,03) = pd(02,03) + rrt(077) * density(01) 
  pd(03,01) = pd(03,01) - rrt(077) * density(03) 
  pd(03,03) = pd(03,03) - rrt(077) * density(01) 
  pd(38,01) = pd(38,01) + rrt(077) * density(03) 
  pd(38,03) = pd(38,03) + rrt(077) * density(01) 
  pd(05,05) = pd(05,05) - rrt(078) * density(38) 
  pd(05,38) = pd(05,38) - rrt(078) * density(05) 
  pd(07,05) = pd(07,05) + rrt(078) * density(38) 
  pd(07,38) = pd(07,38) + rrt(078) * density(05) 
  pd(38,05) = pd(38,05) - rrt(078) * density(38) 
  pd(38,38) = pd(38,38) - rrt(078) * density(05) 
  pd(05,05) = pd(05,05) - rrt(079) * density(05) * 4.0d0
  pd(24,05) = pd(24,05) + rrt(079) * density(05) * 2.0d0
  pd(18,18) = pd(18,18) - rrt(080) * density(38) 
  pd(18,38) = pd(18,38) - rrt(080) * density(18) 
  pd(22,18) = pd(22,18) + rrt(080) * density(38) 
  pd(22,38) = pd(22,38) + rrt(080) * density(18) 
  pd(38,18) = pd(38,18) - rrt(080) * density(38) 
  pd(38,38) = pd(38,38) - rrt(080) * density(18) 
  pd(19,19) = pd(19,19) - rrt(081) * density(38) 
  pd(19,38) = pd(19,38) - rrt(081) * density(19) 
  pd(22,19) = pd(22,19) + rrt(081) * density(38) 
  pd(22,38) = pd(22,38) + rrt(081) * density(19) 
  pd(38,19) = pd(38,19) - rrt(081) * density(38) 
  pd(38,38) = pd(38,38) - rrt(081) * density(19) 
  pd(20,20) = pd(20,20) - rrt(082) * density(38) 
  pd(20,38) = pd(20,38) - rrt(082) * density(20) 
  pd(22,20) = pd(22,20) + rrt(082) * density(38) 
  pd(22,38) = pd(22,38) + rrt(082) * density(20) 
  pd(38,20) = pd(38,20) - rrt(082) * density(38) 
  pd(38,38) = pd(38,38) - rrt(082) * density(20) 
  pd(03,03) = pd(03,03) - rrt(083) * density(07) 
  pd(03,07) = pd(03,07) - rrt(083) * density(03) 
  pd(07,03) = pd(07,03) - rrt(083) * density(07) 
  pd(07,07) = pd(07,07) - rrt(083) * density(03) 
  pd(18,03) = pd(18,03) + rrt(083) * density(07) 
  pd(18,07) = pd(18,07) + rrt(083) * density(03) 
  pd(38,03) = pd(38,03) + rrt(083) * density(07) 
  pd(38,07) = pd(38,07) + rrt(083) * density(03) 
  pd(03,03) = pd(03,03) - rrt(084) * density(09) 
  pd(03,09) = pd(03,09) - rrt(084) * density(03) 
  pd(09,03) = pd(09,03) - rrt(084) * density(09) 
  pd(09,09) = pd(09,09) - rrt(084) * density(03) 
  pd(18,03) = pd(18,03) + rrt(084) * density(09) 
  pd(18,09) = pd(18,09) + rrt(084) * density(03) 
  pd(38,03) = pd(38,03) + rrt(084) * density(09) 
  pd(38,09) = pd(38,09) + rrt(084) * density(03) 
  pd(03,03) = pd(03,03) - rrt(085) * density(08) 
  pd(03,08) = pd(03,08) - rrt(085) * density(03) 
  pd(08,03) = pd(08,03) - rrt(085) * density(08) 
  pd(08,08) = pd(08,08) - rrt(085) * density(03) 
  pd(18,03) = pd(18,03) + rrt(085) * density(08) 
  pd(18,08) = pd(18,08) + rrt(085) * density(03) 
  pd(38,03) = pd(38,03) + rrt(085) * density(08) 
  pd(38,08) = pd(38,08) + rrt(085) * density(03) 
  pd(22,22) = pd(22,22) - rrt(086) * density(38) 
  pd(22,38) = pd(22,38) - rrt(086) * density(22) 
  pd(24,22) = pd(24,22) + rrt(086) * density(38) 
  pd(24,38) = pd(24,38) + rrt(086) * density(22) 
  pd(38,22) = pd(38,22) - rrt(086) * density(38) 
  pd(38,38) = pd(38,38) - rrt(086) * density(22) 
  pd(22,22) = pd(22,22) - rrt(087) * density(22) * 4.0d0
  pd(36,22) = pd(36,22) + rrt(087) * density(22) * 2.0d0
  pd(03,04) = pd(03,04) + rrt(088) * density(38) 
  pd(03,38) = pd(03,38) + rrt(088) * density(04) 
  pd(04,04) = pd(04,04) - rrt(088) * density(38) 
  pd(04,38) = pd(04,38) - rrt(088) * density(04) 
  pd(38,04) = pd(38,04) - rrt(088) * density(38) 
  pd(38,38) = pd(38,38) - rrt(088) * density(04) 
  pd(39,04) = pd(39,04) + rrt(088) * density(38) 
  pd(39,38) = pd(39,38) + rrt(088) * density(04) 
  pd(04,04) = pd(04,04) - rrt(089) * density(05) 
  pd(04,05) = pd(04,05) - rrt(089) * density(04) 
  pd(05,04) = pd(05,04) - rrt(089) * density(05) 
  pd(05,05) = pd(05,05) - rrt(089) * density(04) 
  pd(18,04) = pd(18,04) + rrt(089) * density(05) 
  pd(18,05) = pd(18,05) + rrt(089) * density(04) 
  pd(38,04) = pd(38,04) + rrt(089) * density(05) 
  pd(38,05) = pd(38,05) + rrt(089) * density(04) 
  pd(05,22) = pd(05,22) + rrt(090) * density(38) * 2.0d0
  pd(05,38) = pd(05,38) + rrt(090) * density(22) * 2.0d0
  pd(22,22) = pd(22,22) - rrt(090) * density(38) 
  pd(22,38) = pd(22,38) - rrt(090) * density(22) 
  pd(38,22) = pd(38,22) - rrt(090) * density(38) 
  pd(38,38) = pd(38,38) - rrt(090) * density(22) 
  pd(05,07) = pd(05,07) + rrt(091) * density(10) 
  pd(05,10) = pd(05,10) + rrt(091) * density(07) 
  pd(07,07) = pd(07,07) - rrt(091) * density(10) 
  pd(07,10) = pd(07,10) - rrt(091) * density(07) 
  pd(10,07) = pd(10,07) - rrt(091) * density(10) 
  pd(10,10) = pd(10,10) - rrt(091) * density(07) 
  pd(11,07) = pd(11,07) + rrt(091) * density(10) 
  pd(11,10) = pd(11,10) + rrt(091) * density(07) 
  pd(05,09) = pd(05,09) + rrt(092) * density(10) 
  pd(05,10) = pd(05,10) + rrt(092) * density(09) 
  pd(09,09) = pd(09,09) - rrt(092) * density(10) 
  pd(09,10) = pd(09,10) - rrt(092) * density(09) 
  pd(10,09) = pd(10,09) - rrt(092) * density(10) 
  pd(10,10) = pd(10,10) - rrt(092) * density(09) 
  pd(11,09) = pd(11,09) + rrt(092) * density(10) 
  pd(11,10) = pd(11,10) + rrt(092) * density(09) 
  pd(05,08) = pd(05,08) + rrt(093) * density(10) 
  pd(05,10) = pd(05,10) + rrt(093) * density(08) 
  pd(08,08) = pd(08,08) - rrt(093) * density(10) 
  pd(08,10) = pd(08,10) - rrt(093) * density(08) 
  pd(10,08) = pd(10,08) - rrt(093) * density(10) 
  pd(10,10) = pd(10,10) - rrt(093) * density(08) 
  pd(11,08) = pd(11,08) + rrt(093) * density(10) 
  pd(11,10) = pd(11,10) + rrt(093) * density(08) 
  pd(05,05) = pd(05,05) - rrt(094) * density(22) 
  pd(05,22) = pd(05,22) - rrt(094) * density(05) 
  pd(22,05) = pd(22,05) - rrt(094) * density(22) 
  pd(22,22) = pd(22,22) - rrt(094) * density(05) 
  pd(32,05) = pd(32,05) + rrt(094) * density(22) 
  pd(32,22) = pd(32,22) + rrt(094) * density(05) 
  pd(05,05) = pd(05,05) - rrt(095) * density(17) 
  pd(05,17) = pd(05,17) - rrt(095) * density(05) 
  pd(17,05) = pd(17,05) - rrt(095) * density(17) 
  pd(17,17) = pd(17,17) - rrt(095) * density(05) 
  pd(28,05) = pd(28,05) + rrt(095) * density(17) 
  pd(28,17) = pd(28,17) + rrt(095) * density(05) 
  pd(28,28) = pd(28,28) - rrt(096) * density(38) 
  pd(28,38) = pd(28,38) - rrt(096) * density(28) 
  pd(31,28) = pd(31,28) + rrt(096) * density(38) 
  pd(31,38) = pd(31,38) + rrt(096) * density(28) 
  pd(38,28) = pd(38,28) - rrt(096) * density(38) 
  pd(38,38) = pd(38,38) - rrt(096) * density(28) 
  pd(29,29) = pd(29,29) - rrt(097) * density(38) 
  pd(29,38) = pd(29,38) - rrt(097) * density(29) 
  pd(31,29) = pd(31,29) + rrt(097) * density(38) 
  pd(31,38) = pd(31,38) + rrt(097) * density(29) 
  pd(38,29) = pd(38,29) - rrt(097) * density(38) 
  pd(38,38) = pd(38,38) - rrt(097) * density(29) 
  pd(04,04) = pd(04,04) - rrt(098) * density(36) 
  pd(04,36) = pd(04,36) - rrt(098) * density(04) 
  pd(36,04) = pd(36,04) - rrt(098) * density(36) 
  pd(36,36) = pd(36,36) - rrt(098) * density(04) 
  pd(37,04) = pd(37,04) + rrt(098) * density(36) 
  pd(37,36) = pd(37,36) + rrt(098) * density(04) 
  pd(31,31) = pd(31,31) - rrt(099) * density(39) 
  pd(31,39) = pd(31,39) - rrt(099) * density(31) 
  pd(32,31) = pd(32,31) + rrt(099) * density(39) 
  pd(32,39) = pd(32,39) + rrt(099) * density(31) 
  pd(38,31) = pd(38,31) + rrt(099) * density(39) 
  pd(38,39) = pd(38,39) + rrt(099) * density(31) 
  pd(39,31) = pd(39,31) - rrt(099) * density(39) 
  pd(39,39) = pd(39,39) - rrt(099) * density(31) 
  pd(07,11) = pd(07,11) + rrt(100) * density(24) 
  pd(07,24) = pd(07,24) + rrt(100) * density(11) 
  pd(11,11) = pd(11,11) - rrt(100) * density(24) 
  pd(11,24) = pd(11,24) - rrt(100) * density(11) 
  pd(23,11) = pd(23,11) + rrt(100) * density(24) 
  pd(23,24) = pd(23,24) + rrt(100) * density(11) 
  pd(24,11) = pd(24,11) - rrt(100) * density(24) 
  pd(24,24) = pd(24,24) - rrt(100) * density(11) 
  pd(39,11) = pd(39,11) + rrt(100) * density(24) 
  pd(39,24) = pd(39,24) + rrt(100) * density(11) 
  pd(07,11) = pd(07,11) + rrt(101) * density(26) 
  pd(07,26) = pd(07,26) + rrt(101) * density(11) 
  pd(11,11) = pd(11,11) - rrt(101) * density(26) 
  pd(11,26) = pd(11,26) - rrt(101) * density(11) 
  pd(23,11) = pd(23,11) + rrt(101) * density(26) 
  pd(23,26) = pd(23,26) + rrt(101) * density(11) 
  pd(26,11) = pd(26,11) - rrt(101) * density(26) 
  pd(26,26) = pd(26,26) - rrt(101) * density(11) 
  pd(39,11) = pd(39,11) + rrt(101) * density(26) 
  pd(39,26) = pd(39,26) + rrt(101) * density(11) 
  pd(07,11) = pd(07,11) + rrt(102) * density(25) 
  pd(07,25) = pd(07,25) + rrt(102) * density(11) 
  pd(11,11) = pd(11,11) - rrt(102) * density(25) 
  pd(11,25) = pd(11,25) - rrt(102) * density(11) 
  pd(23,11) = pd(23,11) + rrt(102) * density(25) 
  pd(23,25) = pd(23,25) + rrt(102) * density(11) 
  pd(25,11) = pd(25,11) - rrt(102) * density(25) 
  pd(25,25) = pd(25,25) - rrt(102) * density(11) 
  pd(39,11) = pd(39,11) + rrt(102) * density(25) 
  pd(39,25) = pd(39,25) + rrt(102) * density(11) 
  pd(07,11) = pd(07,11) + rrt(103) * density(18) 
  pd(07,18) = pd(07,18) + rrt(103) * density(11) 
  pd(11,11) = pd(11,11) - rrt(103) * density(18) 
  pd(11,18) = pd(11,18) - rrt(103) * density(11) 
  pd(18,11) = pd(18,11) - rrt(103) * density(18) 
  pd(18,18) = pd(18,18) - rrt(103) * density(11) 
  pd(23,11) = pd(23,11) + rrt(103) * density(18) 
  pd(23,18) = pd(23,18) + rrt(103) * density(11) 
  pd(07,11) = pd(07,11) + rrt(104) * density(19) 
  pd(07,19) = pd(07,19) + rrt(104) * density(11) 
  pd(11,11) = pd(11,11) - rrt(104) * density(19) 
  pd(11,19) = pd(11,19) - rrt(104) * density(11) 
  pd(19,11) = pd(19,11) - rrt(104) * density(19) 
  pd(19,19) = pd(19,19) - rrt(104) * density(11) 
  pd(23,11) = pd(23,11) + rrt(104) * density(19) 
  pd(23,19) = pd(23,19) + rrt(104) * density(11) 
  pd(07,11) = pd(07,11) + rrt(105) * density(20) 
  pd(07,20) = pd(07,20) + rrt(105) * density(11) 
  pd(11,11) = pd(11,11) - rrt(105) * density(20) 
  pd(11,20) = pd(11,20) - rrt(105) * density(11) 
  pd(20,11) = pd(20,11) - rrt(105) * density(20) 
  pd(20,20) = pd(20,20) - rrt(105) * density(11) 
  pd(23,11) = pd(23,11) + rrt(105) * density(20) 
  pd(23,20) = pd(23,20) + rrt(105) * density(11) 
  pd(31,31) = pd(31,31) - rrt(106) * density(38) 
  pd(31,38) = pd(31,38) - rrt(106) * density(31) 
  pd(32,31) = pd(32,31) + rrt(106) * density(38) 
  pd(32,38) = pd(32,38) + rrt(106) * density(31) 
  pd(38,31) = pd(38,31) - rrt(106) * density(38) 
  pd(38,38) = pd(38,38) - rrt(106) * density(31) 
  pd(12,17) = pd(12,17) + rrt(107) * density(38) 
  pd(12,38) = pd(12,38) + rrt(107) * density(17) 
  pd(17,17) = pd(17,17) - rrt(107) * density(38) 
  pd(17,38) = pd(17,38) - rrt(107) * density(17) 
  pd(38,17) = pd(38,17) - rrt(107) * density(38) 
  pd(38,38) = pd(38,38) - rrt(107) * density(17) 
  pd(39,17) = pd(39,17) + rrt(107) * density(38) 
  pd(39,38) = pd(39,38) + rrt(107) * density(17) 
  pd(04,04) = pd(04,04) - rrt(108) * density(32) 
  pd(04,32) = pd(04,32) - rrt(108) * density(04) 
  pd(32,04) = pd(32,04) - rrt(108) * density(32) 
  pd(32,32) = pd(32,32) - rrt(108) * density(04) 
  pd(36,04) = pd(36,04) + rrt(108) * density(32) 
  pd(36,32) = pd(36,32) + rrt(108) * density(04) 
  pd(04,04) = pd(04,04) - rrt(109) * density(33) 
  pd(04,33) = pd(04,33) - rrt(109) * density(04) 
  pd(33,04) = pd(33,04) - rrt(109) * density(33) 
  pd(33,33) = pd(33,33) - rrt(109) * density(04) 
  pd(36,04) = pd(36,04) + rrt(109) * density(33) 
  pd(36,33) = pd(36,33) + rrt(109) * density(04) 
  pd(04,04) = pd(04,04) - rrt(110) * density(34) 
  pd(04,34) = pd(04,34) - rrt(110) * density(04) 
  pd(34,04) = pd(34,04) - rrt(110) * density(34) 
  pd(34,34) = pd(34,34) - rrt(110) * density(04) 
  pd(36,04) = pd(36,04) + rrt(110) * density(34) 
  pd(36,34) = pd(36,34) + rrt(110) * density(04) 
  pd(06,06) = pd(06,06) - rrt(111) * density(07) 
  pd(06,07) = pd(06,07) - rrt(111) * density(06) 
  pd(07,06) = pd(07,06) - rrt(111) * density(07) 
  pd(07,07) = pd(07,07) - rrt(111) * density(06) 
  pd(23,06) = pd(23,06) + rrt(111) * density(07) 
  pd(23,07) = pd(23,07) + rrt(111) * density(06) 
  pd(39,06) = pd(39,06) + rrt(111) * density(07) 
  pd(39,07) = pd(39,07) + rrt(111) * density(06) 
  pd(06,06) = pd(06,06) - rrt(112) * density(09) 
  pd(06,09) = pd(06,09) - rrt(112) * density(06) 
  pd(09,06) = pd(09,06) - rrt(112) * density(09) 
  pd(09,09) = pd(09,09) - rrt(112) * density(06) 
  pd(23,06) = pd(23,06) + rrt(112) * density(09) 
  pd(23,09) = pd(23,09) + rrt(112) * density(06) 
  pd(39,06) = pd(39,06) + rrt(112) * density(09) 
  pd(39,09) = pd(39,09) + rrt(112) * density(06) 
  pd(06,06) = pd(06,06) - rrt(113) * density(08) 
  pd(06,08) = pd(06,08) - rrt(113) * density(06) 
  pd(08,06) = pd(08,06) - rrt(113) * density(08) 
  pd(08,08) = pd(08,08) - rrt(113) * density(06) 
  pd(23,06) = pd(23,06) + rrt(113) * density(08) 
  pd(23,08) = pd(23,08) + rrt(113) * density(06) 
  pd(39,06) = pd(39,06) + rrt(113) * density(08) 
  pd(39,08) = pd(39,08) + rrt(113) * density(06) 
  pd(38,38) = pd(38,38) - rrt(114) * density(38) * 4.0d0
  pd(39,38) = pd(39,38) + rrt(114) * density(38) * 2.0d0
  if( ldensity_constant ) then
    do i = 1, species_max
      if( density_constant(i) ) pd(i,:) = 0.0d0
    enddo
  endif
  if( lgas_heating ) then
    pd(40,1) = eV_to_K * ZDPlasKin_cfg(11)
    pd(40,:) = pd(40,:) * ZDPlasKin_cfg(13)
  endif
  return
end subroutine ZDPlasKin_jex
end module ZDPlasKin
!-----------------------------------------------------------------------------------------------------------------------------------
!
! reaction constant rates
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_reac_rates(Time)
  use ZDPlasKin, only : ZDPlasKin_bolsig_rates, bolsig_rates, bolsig_pointer, ZDPlasKin_cfg, ZDPlasKin_get_density_total, &
                        lreaction_block, rrt
  implicit none
  double precision, intent(in) :: Time
  double precision :: Tgas
  double precision :: Te
  DOUBLE PRECISION, PARAMETER :: R = 8.314D-3
  DOUBLE PRECISION, PARAMETER :: F0 = 1.5702D-06
  DOUBLE PRECISION, PARAMETER :: F1 = 9.0885D-16
  DOUBLE PRECISION, PARAMETER :: F2 = 3.8183D+05
  DOUBLE PRECISION, PARAMETER :: F3 = 4.0458D-04
  DOUBLE PRECISION, PARAMETER :: F4 = 1.1691D-02
  DOUBLE PRECISION, PARAMETER :: F5 = 9.2209D-04
  DOUBLE PRECISION, PARAMETER :: F6 = 3.6707D-02
  DOUBLE PRECISION, PARAMETER :: F7 = 8.7219D-01
  DOUBLE PRECISION, PARAMETER :: F8 = 7.7141D-01
  DOUBLE PRECISION, PARAMETER :: F9 = 7.6561D+00
  DOUBLE PRECISION, PARAMETER :: F10 = 4.7236D-03
  DOUBLE PRECISION, PARAMETER :: F11 = 7.3528D-04
  DOUBLE PRECISION, PARAMETER :: F12 = 1.4673D-02
  DOUBLE PRECISION, PARAMETER :: F13 = 3.8013D-03
  DOUBLE PRECISION, PARAMETER :: F14 = 9.8675D-03
  DOUBLE PRECISION, PARAMETER :: F15 = 4.0357D-03
  DOUBLE PRECISION, PARAMETER :: F16 = 2.2334D+09
  DOUBLE PRECISION, PARAMETER :: F17 = 1.0874D-02
  DOUBLE PRECISION, PARAMETER :: F18 = 1.2339D+14
  DOUBLE PRECISION, PARAMETER :: F19 = 1.0769D-04
  DOUBLE PRECISION, PARAMETER :: F20 = 1.7496D-01
  DOUBLE PRECISION, PARAMETER :: F21 = 4.4750D-02
  DOUBLE PRECISION, PARAMETER :: F22 = 4.7353D-03
  DOUBLE PRECISION, PARAMETER :: F23 = 9.9563D-03
  DOUBLE PRECISION, PARAMETER :: F24 = 3.0270D-03
  DOUBLE PRECISION, PARAMETER :: F25 = 3.7200D-02
  DOUBLE PRECISION, PARAMETER :: F26 = 4.7195D-03
  DOUBLE PRECISION, PARAMETER :: F27 = 1.4808D-02
  DOUBLE PRECISION, PARAMETER :: F28 = 1.8243D-02
  DOUBLE PRECISION, PARAMETER :: F29 = 7.0354D-03
  DOUBLE PRECISION, PARAMETER :: F30 = 4.6848D-03
  DOUBLE PRECISION, PARAMETER :: F31 = 1.0441D+01
  DOUBLE PRECISION, PARAMETER :: F32 = 1.0160D-02
  DOUBLE PRECISION, PARAMETER :: F33 = 3.9466D-14
  DOUBLE PRECISION, PARAMETER :: F34 = 5.1643D+00
  DOUBLE PRECISION, PARAMETER :: F35 = 1.1300D+09
  DOUBLE PRECISION, PARAMETER :: F36 = 5.2790D-01
  call ZDPlasKin_bolsig_rates()
  Tgas = ZDPlasKin_cfg(1)
  Te  = ZDPlasKin_cfg(4)
  rrt(001) = F0*bolsig_rates(bolsig_pointer(1))
  rrt(002) = F0*bolsig_rates(bolsig_pointer(2))
  rrt(003) = F0*bolsig_rates(bolsig_pointer(3))
  rrt(004) = F0*bolsig_rates(bolsig_pointer(4))
  rrt(005) = F0*bolsig_rates(bolsig_pointer(5))
  rrt(006) = F0*bolsig_rates(bolsig_pointer(6))
  rrt(007) = F0*bolsig_rates(bolsig_pointer(7))
  rrt(008) = F0*bolsig_rates(bolsig_pointer(8))
  rrt(009) = F0*bolsig_rates(bolsig_pointer(9))
  rrt(010) = F0*bolsig_rates(bolsig_pointer(10))
  rrt(011) = F0*bolsig_rates(bolsig_pointer(11))
  rrt(012) = F0*bolsig_rates(bolsig_pointer(12))
  rrt(013) = F0*bolsig_rates(bolsig_pointer(13))
  rrt(014) = F0*bolsig_rates(bolsig_pointer(14))
  rrt(015) = F0*bolsig_rates(bolsig_pointer(15))
  rrt(016) = F0*bolsig_rates(bolsig_pointer(16))
  rrt(017) = F0*bolsig_rates(bolsig_pointer(17))
  rrt(018) = F0*bolsig_rates(bolsig_pointer(18))
  rrt(019) = F0*bolsig_rates(bolsig_pointer(19))
  rrt(020) = F0*bolsig_rates(bolsig_pointer(20))
  rrt(021) = F0*bolsig_rates(bolsig_pointer(21))
  rrt(022) = F0*bolsig_rates(bolsig_pointer(22))
  rrt(023) = F0*bolsig_rates(bolsig_pointer(23))
  rrt(024) = F0*bolsig_rates(bolsig_pointer(24))
  rrt(025) = F1*bolsig_rates(bolsig_pointer(25))
  rrt(026) = F1*bolsig_rates(bolsig_pointer(26))
  rrt(027) = F1*bolsig_rates(bolsig_pointer(27))
  rrt(028) = F1*bolsig_rates(bolsig_pointer(28))
  rrt(029) = F1*bolsig_rates(bolsig_pointer(29))
  rrt(030) = F1*bolsig_rates(bolsig_pointer(30))
  rrt(031) = F1*bolsig_rates(bolsig_pointer(31))
  rrt(032) = F1*bolsig_rates(bolsig_pointer(32))
  rrt(033) = F1*bolsig_rates(bolsig_pointer(33))
  rrt(034) = F1*bolsig_rates(bolsig_pointer(34))
  rrt(035) = F1*bolsig_rates(bolsig_pointer(35))
  rrt(036) = F1*bolsig_rates(bolsig_pointer(36))
  rrt(037) = F1*bolsig_rates(bolsig_pointer(37))
  rrt(038) = F1*bolsig_rates(bolsig_pointer(38))
  rrt(039) = F1*bolsig_rates(bolsig_pointer(39))
  rrt(040) = F1*bolsig_rates(bolsig_pointer(40))
  rrt(041) = F1*bolsig_rates(bolsig_pointer(41))
  rrt(042) = F1*bolsig_rates(bolsig_pointer(42))
  rrt(043) = F1*bolsig_rates(bolsig_pointer(43))
  rrt(044) = F1*F2*bolsig_rates(bolsig_pointer(44))
  rrt(045) = F1*F2*bolsig_rates(bolsig_pointer(45))
  rrt(046) = F1*F2*bolsig_rates(bolsig_pointer(46))
  rrt(047) = F1*F2*bolsig_rates(bolsig_pointer(47))
  rrt(048) = F1*F2*bolsig_rates(bolsig_pointer(48))
  rrt(049) = F1*F2*bolsig_rates(bolsig_pointer(49))
  rrt(050) = F1*F2*bolsig_rates(bolsig_pointer(50))
  rrt(051) = F3*bolsig_rates(bolsig_pointer(51))
  rrt(052) = F3*bolsig_rates(bolsig_pointer(52))
  rrt(053) = F3*bolsig_rates(bolsig_pointer(53))
  rrt(054) = F7*bolsig_rates(bolsig_pointer(54))
  rrt(055) = F7*bolsig_rates(bolsig_pointer(55))
  rrt(056) = F7*bolsig_rates(bolsig_pointer(56))
  rrt(057) = F8*bolsig_rates(bolsig_pointer(57))
  rrt(058) = F8*bolsig_rates(bolsig_pointer(58))
  rrt(059) = F8*bolsig_rates(bolsig_pointer(59))
  rrt(060) = F17*1.92D-08*(300./TGAS)**0.71
  rrt(061) = F22*bolsig_rates(bolsig_pointer(60))
  rrt(062) = F22*bolsig_rates(bolsig_pointer(61))
  rrt(063) = F22*bolsig_rates(bolsig_pointer(62))
  rrt(064) = F25*1.60D-08*(300./TGAS)**0.71
  rrt(065) = F26*bolsig_rates(bolsig_pointer(63))
  rrt(066) = F26*bolsig_rates(bolsig_pointer(64))
  rrt(067) = F29*8.98D-09*(300./TGAS)**0.71
  rrt(068) = F31*bolsig_rates(bolsig_pointer(65))
  rrt(069) = F31*bolsig_rates(bolsig_pointer(66))
  rrt(070) = F31*bolsig_rates(bolsig_pointer(67))
  rrt(071) = F33*bolsig_rates(bolsig_pointer(68))
  rrt(072) = F33*bolsig_rates(bolsig_pointer(69))
  rrt(073) = F33*bolsig_rates(bolsig_pointer(70))
  rrt(074) = F34*bolsig_rates(bolsig_pointer(71))
  rrt(075) = F34*bolsig_rates(bolsig_pointer(72))
  rrt(076) = F34*bolsig_rates(bolsig_pointer(73))
  rrt(077) = F35*bolsig_rates(bolsig_pointer(74))
  rrt(078) = F4*1.44D-14*TGAS**(-4.76)*EXP(-1227.98/TGAS)
  rrt(079) = F5*1.87D-6*TGAS**(-7.03)*EXP(-1390.54/TGAS)
  rrt(080) = F6*9.68D-12*(TGAS/298.)**1.28*EXP(-5.40/(R*TGAS))
  rrt(081) = rrt(80)
  rrt(082) = rrt(80)
  rrt(083) = F9*9.97D-11
  rrt(084) = rrt(83)
  rrt(085) = rrt(83)
  rrt(086) = F10*8.65D-7*TGAS**(-0.99)*EXP(-795.17/TGAS)
  rrt(087) = F11*9.55D-12
  rrt(088) = F12*1.00D-11*EXP(7.48/(R*TGAS))
  rrt(089) = F13*7.10D-11
  rrt(090) = F14*5.99D-11
  rrt(091) = F15*1.50D-9
  rrt(092) = rrt(91)
  rrt(093) = rrt(91)
  rrt(094) = F16*1.49D27*TGAS**(-16.82)*EXP(-6575.24/TGAS)
  rrt(095) = F18*3.8D-29
  rrt(096) = F19*1.29D-11*(TGAS/298.)**0.51*EXP(-5.15/(R*TGAS))
  rrt(097) = rrt(96)
  rrt(098) = F20*9.61D-13
  rrt(099) = F21*6.00D-11
  rrt(100) = F23*2.25D-10
  rrt(101) = rrt(100)
  rrt(102) = rrt(100)
  rrt(103) = F24*1.50D-9
  rrt(104) = rrt(103)
  rrt(105) = rrt(103)
  rrt(106) = F27*6.00D-11
  rrt(107) = F28*1.60D-10
  rrt(108) = F30*4.42D-12
  rrt(109) = rrt(108)
  rrt(110) = rrt(108)
  rrt(111) = F32*1.20D-9
  rrt(112) = rrt(111)
  rrt(113) = rrt(111)
  rrt(114) = F36*5.52D-30*TGAS**(-1.00)
  where( lreaction_block(:) ) rrt(:) = 0.0d0
  return
end subroutine ZDPlasKin_reac_rates
!-----------------------------------------------------------------------------------------------------------------------------------
!
! END OF FILE
!
!-----------------------------------------------------------------------------------------------------------------------------------
