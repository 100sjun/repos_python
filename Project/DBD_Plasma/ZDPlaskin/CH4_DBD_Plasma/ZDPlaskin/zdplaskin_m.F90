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
! Thu Aug 22 21:12:23 2024
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
  integer, parameter :: species_max = 45, species_electrons = 1, species_length = 9, reactions_max = 239, reactions_length = 28
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
  integer, parameter, private               :: bolsig_species_max = 26, bolsig_species_length = 9, bolsig_rates_max = 239 
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
  /-1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1,&
    1, 0, 0/
  data species_name(1:species_max) &
  /"E        ","C2H4(V2) ","CH2^+    ","C3H8(V1) ","C2H3     ","C3H4     ","H^+      ","CH2      ","C2H6(V13)","H        ",&
   "C3H4^+   ","C2H6^+   ","C3H6(V)  ","CH4^+    ","C2H4(V1) ","CH3      ","CH3^+    ","C2H6     ","C3H6^+   ","C3H8     ",&
   "C3H5     ","C3H5^+   ","C3H7^+   ","C3H8(V2) ","C3H8^+   ","C2H2(V5) ","C3H7     ","C2H2(V13)","CH       ","C2H5     ",&
   "C3H6     ","C2H2     ","C2H4^+   ","H2       ","C2H6(V24)","C2H2(V2) ","CH4(V13) ","C2H3^+   ","CH4(V24) ","C2H2^+   ",&
   "C2H5^+   ","CH^+     ","H2^+     ","CH4      ","C2H4     "/
  data reaction_sign(1:72) &
  /"bolsig:CH4->CH4(V24)        ","bolsig:CH4->CH4(V13)        ","bolsig:C2H6->C2H6(V24)      ","bolsig:C2H6->C2H6(V13)      ",&
   "bolsig:C2H4->C2H4(V1)       ","bolsig:C2H4->C2H4(V2)       ","bolsig:C2H2->C2H2(V5)       ","bolsig:C2H2->C2H2(V2)       ",&
   "bolsig:C2H2->C2H2(V13)      ","bolsig:C3H8->C3H8(V1)       ","bolsig:C3H8->C3H8(V2)       ","bolsig:C3H6->C3H6(V)        ",&
   "bolsig:CH4->CH3H            ","bolsig:CH4->CH2H2           ","bolsig:CH4->CHH2H           ","bolsig:CH3->CH2H            ",&
   "bolsig:CH3->CHH2            ","bolsig:CH2->CHH             ","bolsig:CH4->CH4^+           ","bolsig:CH4->HCH3^+          ",&
   "bolsig:CH4->H2CH2^+         ","bolsig:CH4->H2HCH^+         ","bolsig:CH3->CH3^+           ","bolsig:CH3->HCH2^+          ",&
   "bolsig:CH3->H2CH^+          ","bolsig:CH2->CH2^+           ","bolsig:CH2->HCH^+           ","bolsig:CH->CH^+             ",&
   "bolsig:C2H6->C2H5H          ","bolsig:C2H6->C2H4H2         ","bolsig:C2H6->C2H3H2H        ","bolsig:C2H6->C2H2H2H2       ",&
   "bolsig:C2H6->CH4CH2         ","bolsig:C2H6->CH3CH3         ","bolsig:C2H5->C2H4H          ","bolsig:C2H5->C2H3H2         ",&
   "bolsig:C2H5->C2H3HH         ","bolsig:C2H5->C2H2H2H        ","bolsig:C2H5->CH4CH          ","bolsig:C2H5->CH3CH2         ",&
   "bolsig:C2H4->C2H3H          ","bolsig:C2H4->C2H2H2         ","bolsig:C2H4->C2H2HH         ","bolsig:C2H4->CH3CH          ",&
   "bolsig:C2H4->CH2CH2         ","bolsig:C2H3->C2H2H          ","bolsig:C2H3->CH2CH          ","bolsig:C2H2->CHCH           ",&
   "bolsig:C3H8->C3H7H          ","bolsig:C3H8->C3H6H2         ","bolsig:C3H8->C3H4H2H2       ","bolsig:C3H8->C2H6CH2        ",&
   "bolsig:C3H8->C2H5CH3        ","bolsig:C3H8->C2H4CH4        ","bolsig:C3H7->C3H6H          ","bolsig:C3H7->C3H5H2         ",&
   "bolsig:C3H7->C3H4H2H        ","bolsig:C3H7->C2H4CH3        ","bolsig:C3H7->C2H3CH4        ","bolsig:C3H6->C3H5H          ",&
   "bolsig:C3H6->C3H4H2         ","bolsig:C3H6->C2H4CH2        ","bolsig:C3H6->C2H3CH3        ","bolsig:C3H6->C2H2CH4        ",&
   "bolsig:C3H5->C3H4H          ","bolsig:C3H5->C2H2CH3        ","bolsig:C3H4->C2H3CH         ","bolsig:C3H4->C2H2CH2        ",&
   "bolsig:C2H6->C2H6^+         ","bolsig:C2H6->HC2H5^+        ","bolsig:C2H6->H2C2H4^+       ","bolsig:C2H6->H2HC2H3^+      "/
  data reaction_sign(73:144) &
  /"bolsig:C2H6->H2H2C2H2^+     ","bolsig:C2H6->CH3CH3^+       ","bolsig:C2H6->CH4CH2^+       ","bolsig:C2H5->C2H5^+         ",&
   "bolsig:C2H5->HC2H4^+        ","bolsig:C2H5->H2C2H3^+       ","bolsig:C2H5->H2HC2H2^+      ","bolsig:C2H5->CH2CH3^+       ",&
   "bolsig:C2H5->CH3CH2^+       ","bolsig:C2H5->CH4CH^+        ","bolsig:C2H4->C2H4^+         ","bolsig:C2H4->HC2H3^+        ",&
   "bolsig:C2H4->CHCH3^+        ","bolsig:C2H4->CH2CH2^+       ","bolsig:C2H4->CH3CH^+        ","bolsig:C2H3->C2H3^+         ",&
   "bolsig:C2H3->HC2H2^+        ","bolsig:C2H3->CHCH2^+        ","bolsig:C2H3->CH2CH^+        ","bolsig:C2H3->C2H2H^+        ",&
   "bolsig:C2H2->C2H2^+         ","bolsig:C2H2->CHCH^+         ","bolsig:C3H8->C3H8^+         ","bolsig:C3H8->HC3H7^+        ",&
   "bolsig:C3H8->H2C3H6^+       ","bolsig:C3H8->H2HC3H5^+      ","bolsig:C3H8->H2H2C3H4^+     ","bolsig:C3H8->CH3C2H5^+      ",&
   "bolsig:C3H8->CH4C2H4^+      ","bolsig:C3H8->C2H5CH3^+      ","bolsig:C3H8->C2H6CH2^+      ","bolsig:C3H7->C3H7^+         ",&
   "bolsig:C3H7->HC3H6^+        ","bolsig:C3H7->H2C3H5^+       ","bolsig:C3H7->H2HC3H4^+      ","bolsig:C3H7->CH2C2H5^+      ",&
   "bolsig:C3H7->CH3C2H4^+      ","bolsig:C3H7->CH4C2H3^+      ","bolsig:C3H7->C2H3CH4^+      ","bolsig:C3H7->C2H4CH3^+      ",&
   "bolsig:C3H7->C2H5CH2^+      ","bolsig:C3H7->C2H6CH^+       ","bolsig:C3H6->C3H6^+         ","bolsig:C3H6->HC3H5^+        ",&
   "bolsig:C3H6->H2C3H4^+       ","bolsig:C3H6->CHC2H5^+       ","bolsig:C3H6->CH2C2H4^+      ","bolsig:C3H6->CH3C2H3^+      ",&
   "bolsig:C3H6->CH4C2H2^+      ","bolsig:C3H6->C2H2CH4^+      ","bolsig:C3H6->C2H3CH3^+      ","bolsig:C3H6->C2H4CH2^+      ",&
   "bolsig:C3H6->C2H5CH^+       ","bolsig:C3H5->C3H5^+         ","bolsig:C3H5->HC3H4^+        ","bolsig:C3H5->CHC2H4^+       ",&
   "bolsig:C3H5->CH2C2H3^+      ","bolsig:C3H5->CH3C2H2^+      ","bolsig:C3H5->C2H2CH3^+      ","bolsig:C3H5->C2H3CH2^+      ",&
   "bolsig:C3H5->C2H4CH^+       ","bolsig:C3H4->C3H4^+         ","bolsig:C3H4->CHC2H3^+       ","bolsig:C3H4->CH2C2H2^+      ",&
   "bolsig:C3H4->C2H2CH2^+      ","bolsig:C3H4->C2H3CH^+       ","bolsig:CH4(V24)->CH3H       ","bolsig:CH4(V13)->CH3H       ",&
   "bolsig:CH4(V24)->CH2H2      ","bolsig:CH4(V13)->CH2H2      ","bolsig:CH4(V24)->CHH2H      ","bolsig:CH4(V13)->CHH2H      "/
  data reaction_sign(145:216) &
  /"bolsig:CH4(V24)->CH4^+      ","bolsig:CH4(V13)->CH4^+      ","bolsig:CH4(V24)->HCH3^+     ","bolsig:CH4(V13)->HCH3^+     ",&
   "bolsig:CH4(V24)->H2CH2^+    ","bolsig:CH4(V13)->H2CH2^+    ","bolsig:CH4(V24)->H2HCH^+    ","bolsig:CH4(V13)->H2HCH^+    ",&
   "bolsig:C2H6(V24)->C2H5H     ","bolsig:C2H6(V13)->C2H5H     ","bolsig:C2H6(V24)->C2H4H2    ","bolsig:C2H6(V13)->C2H4H2    ",&
   "bolsig:C2H6(V24)->C2H3H2H   ","bolsig:C2H6(V13)->C2H3H2H   ","bolsig:C2H6(V24)->C2H2H2H2  ","bolsig:C2H6(V13)->C2H2H2H2  ",&
   "bolsig:C2H6(V24)->CH4CH2    ","bolsig:C2H6(V13)->CH4CH2    ","bolsig:C2H6(V24)->CH3CH3    ","bolsig:C2H6(V13)->CH3CH3    ",&
   "bolsig:C2H4(V1)->C2H3H      ","bolsig:C2H4(V2)->C2H3H      ","bolsig:C2H4(V1)->C2H2H2     ","bolsig:C2H4(V2)->C2H2H2     ",&
   "bolsig:C2H4(V1)->C2H2HH     ","bolsig:C2H4(V2)->C2H2HH     ","bolsig:C2H4(V1)->CH3CH      ","bolsig:C2H4(V2)->CH3CH      ",&
   "bolsig:C2H4(V1)->CH2CH2     ","bolsig:C2H4(V2)->CH2CH2     ","bolsig:C2H2(V5)->CHCH       ","bolsig:C2H2(V2)->CHCH       ",&
   "bolsig:C2H2(V13)->CHCH      ","bolsig:C3H8(V1)->C3H7H      ","bolsig:C3H8(V2)->C3H7H      ","bolsig:C3H8(V1)->C3H6H2     ",&
   "bolsig:C3H8(V2)->C3H6H2     ","bolsig:C3H8(V1)->C3H4H2H2   ","bolsig:C3H8(V2)->C3H4H2H2   ","bolsig:C3H8(V1)->C2H6CH2    ",&
   "bolsig:C3H8(V2)->C2H6CH2    ","bolsig:C3H8(V1)->C2H5CH3    ","bolsig:C3H8(V2)->C2H5CH3    ","bolsig:C3H8(V1)->C2H4CH4    ",&
   "bolsig:C3H8(V2)->C2H4CH4    ","bolsig:C2H6(V24)->C2H6^+    ","bolsig:C2H6(V13)->C2H6^+    ","bolsig:C2H6(V24)->HC2H5^+   ",&
   "bolsig:C2H6(V13)->HC2H5^+   ","bolsig:C2H6(V24)->H2C2H4^+  ","bolsig:C2H6(V13)->H2C2H4^+  ","bolsig:C2H6(V24)->H2HC2H3^+ ",&
   "bolsig:C2H6(V13)->H2HC2H3^+ ","bolsig:C2H6(V24)->H2H2C2H2^+","bolsig:C2H6(V13)->H2H2C2H2^+","bolsig:C2H6(V24)->CH3CH3^+  ",&
   "bolsig:C2H6(V13)->CH3CH3^+  ","bolsig:C2H6(V24)->CH4CH2^+  ","bolsig:C2H6(V13)->CH4CH2^+  ","bolsig:C2H4(V1)->C2H4^+     ",&
   "bolsig:C2H4(V2)->C2H4^+     ","bolsig:C2H4(V1)->HC2H3^+    ","bolsig:C2H4(V2)->HC2H3^+    ","bolsig:C2H4(V1)->CHCH3^+    ",&
   "bolsig:C2H4(V2)->CHCH3^+    ","bolsig:C2H4(V1)->CH2CH2^+   ","bolsig:C2H4(V2)->CH2CH2^+   ","bolsig:C2H4(V1)->CH3CH^+    ",&
   "bolsig:C2H4(V2)->CH3CH^+    ","bolsig:C2H2(V5)->C2H2^+     ","bolsig:C2H2(V2)->C2H2^+     ","bolsig:C2H2(V13)->C2H2^+    "/
  data reaction_sign(217:239) &
  /"bolsig:C2H2(V5)->CHCH^+     ","bolsig:C2H2(V2)->CHCH^+     ","bolsig:C2H2(V13)->CHCH^+    ","bolsig:C3H8(V1)->C3H8^+     ",&
   "bolsig:C3H8(V2)->C3H8^+     ","bolsig:C3H8(V1)->HC3H7^+    ","bolsig:C3H8(V2)->HC3H7^+    ","bolsig:C3H8(V1)->H2C3H6^+   ",&
   "bolsig:C3H8(V2)->H2C3H6^+   ","bolsig:C3H8(V1)->H2HC3H5^+  ","bolsig:C3H8(V2)->H2HC3H5^+  ","bolsig:C3H8(V1)->H2H2C3H4^+ ",&
   "bolsig:C3H8(V2)->H2H2C3H4^+ ","bolsig:C3H8(V1)->CH3C2H5^+  ","bolsig:C3H8(V2)->CH3C2H5^+  ","bolsig:C3H8(V1)->CH4C2H4^+  ",&
   "bolsig:C3H8(V2)->CH4C2H4^+  ","bolsig:C3H8(V1)->C2H5CH3^+  ","bolsig:C3H8(V2)->C2H5CH3^+  ","bolsig:C3H8(V1)->C2H6CH2^+  ",&
   "bolsig:C3H8(V2)->C2H6CH2^+  ","bolsig:H2->HH               ","bolsig:H2->H2^+             "/
  data bolsig_species(1:bolsig_species_max) &
  /"C2H4(V2) ","C2H3     ","C3H8(V1) ","C3H4     ","CH2      ","C2H6(V13)","C2H4(V1) ","CH3      ","C2H6     ","C3H8     ",&
   "C3H5     ","C3H8(V2) ","C2H2(V5) ","C3H7     ","C2H5     ","CH       ","C2H2(V13)","C3H6     ","C2H2     ","H2       ",&
   "C2H6(V24)","C2H2(V2) ","CH4(V13) ","CH4(V24) ","CH4      ","C2H4     "/
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
211 format(i3,1x,A28)
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
311 format(331x,45(1x,i9))
312 format(A3,1x,A28,1x,45(1x,A9))
313 format(i3,1x,A28,1x,45(1x,1pd9.2))
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
      write(ifile_unit,"(239(i3))",err=200) int(mrtm(i,:))
    enddo
    lerror = .false.
200 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file " // &
                                     "<qt_matrix.txt> (subroutine writer_save_qtplaskin)")
    close(ifile_unit)
    open(ifile_unit,file="qt_densities.txt",action="write",err=300)
    write(ifile_unit,"(1x,A14,45(111x,i2.2))",err=300) "Time_s", ( i, i = 1, species_max )
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
    write(ifile_unit,"(1x,A12,239(101x,i3.3))",err=500) "Time_s", ( i, i = 1, reactions_max )
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
      write(ifile_unit,"(1pe15.6,45(1pe13.4))") densav(0,2), densav(1:,2)
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
      write(ifile_unit,"(240(1pe13.4))") densav(0,2), rrt_loc(:)
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
  reac_source_local(39,001) = + reac_rate_local(001) 
  reac_source_local(44,001) = - reac_rate_local(001) 
  reac_source_local(37,002) = + reac_rate_local(002) 
  reac_source_local(44,002) = - reac_rate_local(002) 
  reac_source_local(18,003) = - reac_rate_local(003) 
  reac_source_local(35,003) = + reac_rate_local(003) 
  reac_source_local(09,004) = + reac_rate_local(004) 
  reac_source_local(18,004) = - reac_rate_local(004) 
  reac_source_local(15,005) = + reac_rate_local(005) 
  reac_source_local(45,005) = - reac_rate_local(005) 
  reac_source_local(02,006) = + reac_rate_local(006) 
  reac_source_local(45,006) = - reac_rate_local(006) 
  reac_source_local(26,007) = + reac_rate_local(007) 
  reac_source_local(32,007) = - reac_rate_local(007) 
  reac_source_local(32,008) = - reac_rate_local(008) 
  reac_source_local(36,008) = + reac_rate_local(008) 
  reac_source_local(28,009) = + reac_rate_local(009) 
  reac_source_local(32,009) = - reac_rate_local(009) 
  reac_source_local(04,010) = + reac_rate_local(010) 
  reac_source_local(20,010) = - reac_rate_local(010) 
  reac_source_local(20,011) = - reac_rate_local(011) 
  reac_source_local(24,011) = + reac_rate_local(011) 
  reac_source_local(13,012) = + reac_rate_local(012) 
  reac_source_local(31,012) = - reac_rate_local(012) 
  reac_source_local(10,013) = + reac_rate_local(013) 
  reac_source_local(16,013) = + reac_rate_local(013) 
  reac_source_local(44,013) = - reac_rate_local(013) 
  reac_source_local(08,014) = + reac_rate_local(014) 
  reac_source_local(34,014) = + reac_rate_local(014) 
  reac_source_local(44,014) = - reac_rate_local(014) 
  reac_source_local(10,015) = + reac_rate_local(015) 
  reac_source_local(29,015) = + reac_rate_local(015) 
  reac_source_local(34,015) = + reac_rate_local(015) 
  reac_source_local(44,015) = - reac_rate_local(015) 
  reac_source_local(08,016) = + reac_rate_local(016) 
  reac_source_local(10,016) = + reac_rate_local(016) 
  reac_source_local(16,016) = - reac_rate_local(016) 
  reac_source_local(16,017) = - reac_rate_local(017) 
  reac_source_local(29,017) = + reac_rate_local(017) 
  reac_source_local(34,017) = + reac_rate_local(017) 
  reac_source_local(08,018) = - reac_rate_local(018) 
  reac_source_local(10,018) = + reac_rate_local(018) 
  reac_source_local(29,018) = + reac_rate_local(018) 
  reac_source_local(01,019) = + reac_rate_local(019) 
  reac_source_local(14,019) = + reac_rate_local(019) 
  reac_source_local(44,019) = - reac_rate_local(019) 
  reac_source_local(01,020) = + reac_rate_local(020) 
  reac_source_local(10,020) = + reac_rate_local(020) 
  reac_source_local(17,020) = + reac_rate_local(020) 
  reac_source_local(44,020) = - reac_rate_local(020) 
  reac_source_local(01,021) = + reac_rate_local(021) 
  reac_source_local(03,021) = + reac_rate_local(021) 
  reac_source_local(34,021) = + reac_rate_local(021) 
  reac_source_local(44,021) = - reac_rate_local(021) 
  reac_source_local(01,022) = + reac_rate_local(022) 
  reac_source_local(10,022) = + reac_rate_local(022) 
  reac_source_local(34,022) = + reac_rate_local(022) 
  reac_source_local(42,022) = + reac_rate_local(022) 
  reac_source_local(44,022) = - reac_rate_local(022) 
  reac_source_local(01,023) = + reac_rate_local(023) 
  reac_source_local(16,023) = - reac_rate_local(023) 
  reac_source_local(17,023) = + reac_rate_local(023) 
  reac_source_local(01,024) = + reac_rate_local(024) 
  reac_source_local(03,024) = + reac_rate_local(024) 
  reac_source_local(10,024) = + reac_rate_local(024) 
  reac_source_local(16,024) = - reac_rate_local(024) 
  reac_source_local(01,025) = + reac_rate_local(025) 
  reac_source_local(16,025) = - reac_rate_local(025) 
  reac_source_local(34,025) = + reac_rate_local(025) 
  reac_source_local(42,025) = + reac_rate_local(025) 
  reac_source_local(01,026) = + reac_rate_local(026) 
  reac_source_local(03,026) = + reac_rate_local(026) 
  reac_source_local(08,026) = - reac_rate_local(026) 
  reac_source_local(01,027) = + reac_rate_local(027) 
  reac_source_local(08,027) = - reac_rate_local(027) 
  reac_source_local(10,027) = + reac_rate_local(027) 
  reac_source_local(42,027) = + reac_rate_local(027) 
  reac_source_local(01,028) = + reac_rate_local(028) 
  reac_source_local(29,028) = - reac_rate_local(028) 
  reac_source_local(42,028) = + reac_rate_local(028) 
  reac_source_local(10,029) = + reac_rate_local(029) 
  reac_source_local(18,029) = - reac_rate_local(029) 
  reac_source_local(30,029) = + reac_rate_local(029) 
  reac_source_local(18,030) = - reac_rate_local(030) 
  reac_source_local(34,030) = + reac_rate_local(030) 
  reac_source_local(45,030) = + reac_rate_local(030) 
  reac_source_local(05,031) = + reac_rate_local(031) 
  reac_source_local(10,031) = + reac_rate_local(031) 
  reac_source_local(18,031) = - reac_rate_local(031) 
  reac_source_local(34,031) = + reac_rate_local(031) 
  reac_source_local(18,032) = - reac_rate_local(032) 
  reac_source_local(32,032) = + reac_rate_local(032) 
  reac_source_local(34,032) = + reac_rate_local(032) * 2.d0
  reac_source_local(08,033) = + reac_rate_local(033) 
  reac_source_local(18,033) = - reac_rate_local(033) 
  reac_source_local(44,033) = + reac_rate_local(033) 
  reac_source_local(16,034) = + reac_rate_local(034) * 2.d0
  reac_source_local(18,034) = - reac_rate_local(034) 
  reac_source_local(10,035) = + reac_rate_local(035) 
  reac_source_local(30,035) = - reac_rate_local(035) 
  reac_source_local(45,035) = + reac_rate_local(035) 
  reac_source_local(05,036) = + reac_rate_local(036) 
  reac_source_local(30,036) = - reac_rate_local(036) 
  reac_source_local(34,036) = + reac_rate_local(036) 
  reac_source_local(05,037) = + reac_rate_local(037) 
  reac_source_local(10,037) = + reac_rate_local(037) * 2.d0
  reac_source_local(30,037) = - reac_rate_local(037) 
  reac_source_local(10,038) = + reac_rate_local(038) 
  reac_source_local(30,038) = - reac_rate_local(038) 
  reac_source_local(32,038) = + reac_rate_local(038) 
  reac_source_local(34,038) = + reac_rate_local(038) 
  reac_source_local(29,039) = + reac_rate_local(039) 
  reac_source_local(30,039) = - reac_rate_local(039) 
  reac_source_local(44,039) = + reac_rate_local(039) 
  reac_source_local(08,040) = + reac_rate_local(040) 
  reac_source_local(16,040) = + reac_rate_local(040) 
  reac_source_local(30,040) = - reac_rate_local(040) 
  reac_source_local(05,041) = + reac_rate_local(041) 
  reac_source_local(10,041) = + reac_rate_local(041) 
  reac_source_local(45,041) = - reac_rate_local(041) 
  reac_source_local(32,042) = + reac_rate_local(042) 
  reac_source_local(34,042) = + reac_rate_local(042) 
  reac_source_local(45,042) = - reac_rate_local(042) 
  reac_source_local(10,043) = + reac_rate_local(043) * 2.d0
  reac_source_local(32,043) = + reac_rate_local(043) 
  reac_source_local(45,043) = - reac_rate_local(043) 
  reac_source_local(16,044) = + reac_rate_local(044) 
  reac_source_local(29,044) = + reac_rate_local(044) 
  reac_source_local(45,044) = - reac_rate_local(044) 
  reac_source_local(08,045) = + reac_rate_local(045) * 2.d0
  reac_source_local(45,045) = - reac_rate_local(045) 
  reac_source_local(05,046) = - reac_rate_local(046) 
  reac_source_local(10,046) = + reac_rate_local(046) 
  reac_source_local(32,046) = + reac_rate_local(046) 
  reac_source_local(05,047) = - reac_rate_local(047) 
  reac_source_local(08,047) = + reac_rate_local(047) 
  reac_source_local(29,047) = + reac_rate_local(047) 
  reac_source_local(29,048) = + reac_rate_local(048) * 2.d0
  reac_source_local(32,048) = - reac_rate_local(048) 
  reac_source_local(10,049) = + reac_rate_local(049) 
  reac_source_local(20,049) = - reac_rate_local(049) 
  reac_source_local(27,049) = + reac_rate_local(049) 
  reac_source_local(20,050) = - reac_rate_local(050) 
  reac_source_local(31,050) = + reac_rate_local(050) 
  reac_source_local(34,050) = + reac_rate_local(050) 
  reac_source_local(06,051) = + reac_rate_local(051) 
  reac_source_local(20,051) = - reac_rate_local(051) 
  reac_source_local(34,051) = + reac_rate_local(051) * 2.d0
  reac_source_local(08,052) = + reac_rate_local(052) 
  reac_source_local(18,052) = + reac_rate_local(052) 
  reac_source_local(20,052) = - reac_rate_local(052) 
  reac_source_local(16,053) = + reac_rate_local(053) 
  reac_source_local(20,053) = - reac_rate_local(053) 
  reac_source_local(30,053) = + reac_rate_local(053) 
  reac_source_local(20,054) = - reac_rate_local(054) 
  reac_source_local(44,054) = + reac_rate_local(054) 
  reac_source_local(45,054) = + reac_rate_local(054) 
  reac_source_local(10,055) = + reac_rate_local(055) 
  reac_source_local(27,055) = - reac_rate_local(055) 
  reac_source_local(31,055) = + reac_rate_local(055) 
  reac_source_local(21,056) = + reac_rate_local(056) 
  reac_source_local(27,056) = - reac_rate_local(056) 
  reac_source_local(34,056) = + reac_rate_local(056) 
  reac_source_local(06,057) = + reac_rate_local(057) 
  reac_source_local(10,057) = + reac_rate_local(057) 
  reac_source_local(27,057) = - reac_rate_local(057) 
  reac_source_local(34,057) = + reac_rate_local(057) 
  reac_source_local(16,058) = + reac_rate_local(058) 
  reac_source_local(27,058) = - reac_rate_local(058) 
  reac_source_local(45,058) = + reac_rate_local(058) 
  reac_source_local(05,059) = + reac_rate_local(059) 
  reac_source_local(27,059) = - reac_rate_local(059) 
  reac_source_local(44,059) = + reac_rate_local(059) 
  reac_source_local(10,060) = + reac_rate_local(060) 
  reac_source_local(21,060) = + reac_rate_local(060) 
  reac_source_local(31,060) = - reac_rate_local(060) 
  reac_source_local(06,061) = + reac_rate_local(061) 
  reac_source_local(31,061) = - reac_rate_local(061) 
  reac_source_local(34,061) = + reac_rate_local(061) 
  reac_source_local(08,062) = + reac_rate_local(062) 
  reac_source_local(31,062) = - reac_rate_local(062) 
  reac_source_local(45,062) = + reac_rate_local(062) 
  reac_source_local(05,063) = + reac_rate_local(063) 
  reac_source_local(16,063) = + reac_rate_local(063) 
  reac_source_local(31,063) = - reac_rate_local(063) 
  reac_source_local(31,064) = - reac_rate_local(064) 
  reac_source_local(32,064) = + reac_rate_local(064) 
  reac_source_local(44,064) = + reac_rate_local(064) 
  reac_source_local(06,065) = + reac_rate_local(065) 
  reac_source_local(10,065) = + reac_rate_local(065) 
  reac_source_local(21,065) = - reac_rate_local(065) 
  reac_source_local(16,066) = + reac_rate_local(066) 
  reac_source_local(21,066) = - reac_rate_local(066) 
  reac_source_local(32,066) = + reac_rate_local(066) 
  reac_source_local(05,067) = + reac_rate_local(067) 
  reac_source_local(06,067) = - reac_rate_local(067) 
  reac_source_local(29,067) = + reac_rate_local(067) 
  reac_source_local(06,068) = - reac_rate_local(068) 
  reac_source_local(08,068) = + reac_rate_local(068) 
  reac_source_local(32,068) = + reac_rate_local(068) 
  reac_source_local(01,069) = + reac_rate_local(069) 
  reac_source_local(12,069) = + reac_rate_local(069) 
  reac_source_local(18,069) = - reac_rate_local(069) 
  reac_source_local(01,070) = + reac_rate_local(070) 
  reac_source_local(10,070) = + reac_rate_local(070) 
  reac_source_local(18,070) = - reac_rate_local(070) 
  reac_source_local(41,070) = + reac_rate_local(070) 
  reac_source_local(01,071) = + reac_rate_local(071) 
  reac_source_local(18,071) = - reac_rate_local(071) 
  reac_source_local(33,071) = + reac_rate_local(071) 
  reac_source_local(34,071) = + reac_rate_local(071) 
  reac_source_local(01,072) = + reac_rate_local(072) 
  reac_source_local(10,072) = + reac_rate_local(072) 
  reac_source_local(18,072) = - reac_rate_local(072) 
  reac_source_local(34,072) = + reac_rate_local(072) 
  reac_source_local(38,072) = + reac_rate_local(072) 
  reac_source_local(01,073) = + reac_rate_local(073) 
  reac_source_local(18,073) = - reac_rate_local(073) 
  reac_source_local(34,073) = + reac_rate_local(073) * 2.d0
  reac_source_local(40,073) = + reac_rate_local(073) 
  reac_source_local(01,074) = + reac_rate_local(074) 
  reac_source_local(16,074) = + reac_rate_local(074) 
  reac_source_local(17,074) = + reac_rate_local(074) 
  reac_source_local(18,074) = - reac_rate_local(074) 
  reac_source_local(01,075) = + reac_rate_local(075) 
  reac_source_local(03,075) = + reac_rate_local(075) 
  reac_source_local(18,075) = - reac_rate_local(075) 
  reac_source_local(44,075) = + reac_rate_local(075) 
  reac_source_local(01,076) = + reac_rate_local(076) 
  reac_source_local(30,076) = - reac_rate_local(076) 
  reac_source_local(41,076) = + reac_rate_local(076) 
  reac_source_local(01,077) = + reac_rate_local(077) 
  reac_source_local(10,077) = + reac_rate_local(077) 
  reac_source_local(30,077) = - reac_rate_local(077) 
  reac_source_local(33,077) = + reac_rate_local(077) 
  reac_source_local(01,078) = + reac_rate_local(078) 
  reac_source_local(30,078) = - reac_rate_local(078) 
  reac_source_local(34,078) = + reac_rate_local(078) 
  reac_source_local(38,078) = + reac_rate_local(078) 
  reac_source_local(01,079) = + reac_rate_local(079) 
  reac_source_local(10,079) = + reac_rate_local(079) 
  reac_source_local(30,079) = - reac_rate_local(079) 
  reac_source_local(34,079) = + reac_rate_local(079) 
  reac_source_local(40,079) = + reac_rate_local(079) 
  reac_source_local(01,080) = + reac_rate_local(080) 
  reac_source_local(08,080) = + reac_rate_local(080) 
  reac_source_local(17,080) = + reac_rate_local(080) 
  reac_source_local(30,080) = - reac_rate_local(080) 
  reac_source_local(01,081) = + reac_rate_local(081) 
  reac_source_local(03,081) = + reac_rate_local(081) 
  reac_source_local(16,081) = + reac_rate_local(081) 
  reac_source_local(30,081) = - reac_rate_local(081) 
  reac_source_local(01,082) = + reac_rate_local(082) 
  reac_source_local(30,082) = - reac_rate_local(082) 
  reac_source_local(42,082) = + reac_rate_local(082) 
  reac_source_local(44,082) = + reac_rate_local(082) 
  reac_source_local(01,083) = + reac_rate_local(083) 
  reac_source_local(33,083) = + reac_rate_local(083) 
  reac_source_local(45,083) = - reac_rate_local(083) 
  reac_source_local(01,084) = + reac_rate_local(084) 
  reac_source_local(10,084) = + reac_rate_local(084) 
  reac_source_local(38,084) = + reac_rate_local(084) 
  reac_source_local(45,084) = - reac_rate_local(084) 
  reac_source_local(01,085) = + reac_rate_local(085) 
  reac_source_local(17,085) = + reac_rate_local(085) 
  reac_source_local(29,085) = + reac_rate_local(085) 
  reac_source_local(45,085) = - reac_rate_local(085) 
  reac_source_local(01,086) = + reac_rate_local(086) 
  reac_source_local(03,086) = + reac_rate_local(086) 
  reac_source_local(08,086) = + reac_rate_local(086) 
  reac_source_local(45,086) = - reac_rate_local(086) 
  reac_source_local(01,087) = + reac_rate_local(087) 
  reac_source_local(16,087) = + reac_rate_local(087) 
  reac_source_local(42,087) = + reac_rate_local(087) 
  reac_source_local(45,087) = - reac_rate_local(087) 
  reac_source_local(01,088) = + reac_rate_local(088) 
  reac_source_local(05,088) = - reac_rate_local(088) 
  reac_source_local(38,088) = + reac_rate_local(088) 
  reac_source_local(01,089) = + reac_rate_local(089) 
  reac_source_local(05,089) = - reac_rate_local(089) 
  reac_source_local(10,089) = + reac_rate_local(089) 
  reac_source_local(40,089) = + reac_rate_local(089) 
  reac_source_local(01,090) = + reac_rate_local(090) 
  reac_source_local(03,090) = + reac_rate_local(090) 
  reac_source_local(05,090) = - reac_rate_local(090) 
  reac_source_local(29,090) = + reac_rate_local(090) 
  reac_source_local(01,091) = + reac_rate_local(091) 
  reac_source_local(05,091) = - reac_rate_local(091) 
  reac_source_local(08,091) = + reac_rate_local(091) 
  reac_source_local(42,091) = + reac_rate_local(091) 
  reac_source_local(01,092) = + reac_rate_local(092) 
  reac_source_local(05,092) = - reac_rate_local(092) 
  reac_source_local(07,092) = + reac_rate_local(092) 
  reac_source_local(32,092) = + reac_rate_local(092) 
  reac_source_local(01,093) = + reac_rate_local(093) 
  reac_source_local(32,093) = - reac_rate_local(093) 
  reac_source_local(40,093) = + reac_rate_local(093) 
  reac_source_local(01,094) = + reac_rate_local(094) 
  reac_source_local(29,094) = + reac_rate_local(094) 
  reac_source_local(32,094) = - reac_rate_local(094) 
  reac_source_local(42,094) = + reac_rate_local(094) 
  reac_source_local(01,095) = + reac_rate_local(095) 
  reac_source_local(20,095) = - reac_rate_local(095) 
  reac_source_local(25,095) = + reac_rate_local(095) 
  reac_source_local(01,096) = + reac_rate_local(096) 
  reac_source_local(10,096) = + reac_rate_local(096) 
  reac_source_local(20,096) = - reac_rate_local(096) 
  reac_source_local(23,096) = + reac_rate_local(096) 
  reac_source_local(01,097) = + reac_rate_local(097) 
  reac_source_local(19,097) = + reac_rate_local(097) 
  reac_source_local(20,097) = - reac_rate_local(097) 
  reac_source_local(34,097) = + reac_rate_local(097) 
  reac_source_local(01,098) = + reac_rate_local(098) 
  reac_source_local(10,098) = + reac_rate_local(098) 
  reac_source_local(20,098) = - reac_rate_local(098) 
  reac_source_local(22,098) = + reac_rate_local(098) 
  reac_source_local(34,098) = + reac_rate_local(098) 
  reac_source_local(01,099) = + reac_rate_local(099) 
  reac_source_local(11,099) = + reac_rate_local(099) 
  reac_source_local(20,099) = - reac_rate_local(099) 
  reac_source_local(34,099) = + reac_rate_local(099) * 2.d0
  reac_source_local(01,100) = + reac_rate_local(100) 
  reac_source_local(16,100) = + reac_rate_local(100) 
  reac_source_local(20,100) = - reac_rate_local(100) 
  reac_source_local(41,100) = + reac_rate_local(100) 
  reac_source_local(01,101) = + reac_rate_local(101) 
  reac_source_local(20,101) = - reac_rate_local(101) 
  reac_source_local(33,101) = + reac_rate_local(101) 
  reac_source_local(44,101) = + reac_rate_local(101) 
  reac_source_local(01,102) = + reac_rate_local(102) 
  reac_source_local(17,102) = + reac_rate_local(102) 
  reac_source_local(20,102) = - reac_rate_local(102) 
  reac_source_local(30,102) = + reac_rate_local(102) 
  reac_source_local(01,103) = + reac_rate_local(103) 
  reac_source_local(03,103) = + reac_rate_local(103) 
  reac_source_local(18,103) = + reac_rate_local(103) 
  reac_source_local(20,103) = - reac_rate_local(103) 
  reac_source_local(01,104) = + reac_rate_local(104) 
  reac_source_local(23,104) = + reac_rate_local(104) 
  reac_source_local(27,104) = - reac_rate_local(104) 
  reac_source_local(01,105) = + reac_rate_local(105) 
  reac_source_local(10,105) = + reac_rate_local(105) 
  reac_source_local(19,105) = + reac_rate_local(105) 
  reac_source_local(27,105) = - reac_rate_local(105) 
  reac_source_local(01,106) = + reac_rate_local(106) 
  reac_source_local(22,106) = + reac_rate_local(106) 
  reac_source_local(27,106) = - reac_rate_local(106) 
  reac_source_local(34,106) = + reac_rate_local(106) 
  reac_source_local(01,107) = + reac_rate_local(107) 
  reac_source_local(10,107) = + reac_rate_local(107) 
  reac_source_local(11,107) = + reac_rate_local(107) 
  reac_source_local(27,107) = - reac_rate_local(107) 
  reac_source_local(34,107) = + reac_rate_local(107) 
  reac_source_local(01,108) = + reac_rate_local(108) 
  reac_source_local(08,108) = + reac_rate_local(108) 
  reac_source_local(27,108) = - reac_rate_local(108) 
  reac_source_local(41,108) = + reac_rate_local(108) 
  reac_source_local(01,109) = + reac_rate_local(109) 
  reac_source_local(16,109) = + reac_rate_local(109) 
  reac_source_local(27,109) = - reac_rate_local(109) 
  reac_source_local(33,109) = + reac_rate_local(109) 
  reac_source_local(01,110) = + reac_rate_local(110) 
  reac_source_local(27,110) = - reac_rate_local(110) 
  reac_source_local(38,110) = + reac_rate_local(110) 
  reac_source_local(44,110) = + reac_rate_local(110) 
  reac_source_local(01,111) = + reac_rate_local(111) 
  reac_source_local(05,111) = + reac_rate_local(111) 
  reac_source_local(14,111) = + reac_rate_local(111) 
  reac_source_local(27,111) = - reac_rate_local(111) 
  reac_source_local(01,112) = + reac_rate_local(112) 
  reac_source_local(17,112) = + reac_rate_local(112) 
  reac_source_local(27,112) = - reac_rate_local(112) 
  reac_source_local(45,112) = + reac_rate_local(112) 
  reac_source_local(01,113) = + reac_rate_local(113) 
  reac_source_local(03,113) = + reac_rate_local(113) 
  reac_source_local(27,113) = - reac_rate_local(113) 
  reac_source_local(30,113) = + reac_rate_local(113) 
  reac_source_local(01,114) = + reac_rate_local(114) 
  reac_source_local(18,114) = + reac_rate_local(114) 
  reac_source_local(27,114) = - reac_rate_local(114) 
  reac_source_local(42,114) = + reac_rate_local(114) 
  reac_source_local(01,115) = + reac_rate_local(115) 
  reac_source_local(19,115) = + reac_rate_local(115) 
  reac_source_local(31,115) = - reac_rate_local(115) 
  reac_source_local(01,116) = + reac_rate_local(116) 
  reac_source_local(10,116) = + reac_rate_local(116) 
  reac_source_local(22,116) = + reac_rate_local(116) 
  reac_source_local(31,116) = - reac_rate_local(116) 
  reac_source_local(01,117) = + reac_rate_local(117) 
  reac_source_local(11,117) = + reac_rate_local(117) 
  reac_source_local(31,117) = - reac_rate_local(117) 
  reac_source_local(34,117) = + reac_rate_local(117) 
  reac_source_local(01,118) = + reac_rate_local(118) 
  reac_source_local(29,118) = + reac_rate_local(118) 
  reac_source_local(31,118) = - reac_rate_local(118) 
  reac_source_local(41,118) = + reac_rate_local(118) 
  reac_source_local(01,119) = + reac_rate_local(119) 
  reac_source_local(08,119) = + reac_rate_local(119) 
  reac_source_local(31,119) = - reac_rate_local(119) 
  reac_source_local(33,119) = + reac_rate_local(119) 
  reac_source_local(01,120) = + reac_rate_local(120) 
  reac_source_local(16,120) = + reac_rate_local(120) 
  reac_source_local(31,120) = - reac_rate_local(120) 
  reac_source_local(38,120) = + reac_rate_local(120) 
  reac_source_local(01,121) = + reac_rate_local(121) 
  reac_source_local(31,121) = - reac_rate_local(121) 
  reac_source_local(40,121) = + reac_rate_local(121) 
  reac_source_local(44,121) = + reac_rate_local(121) 
  reac_source_local(01,122) = + reac_rate_local(122) 
  reac_source_local(14,122) = + reac_rate_local(122) 
  reac_source_local(31,122) = - reac_rate_local(122) 
  reac_source_local(32,122) = + reac_rate_local(122) 
  reac_source_local(01,123) = + reac_rate_local(123) 
  reac_source_local(05,123) = + reac_rate_local(123) 
  reac_source_local(17,123) = + reac_rate_local(123) 
  reac_source_local(31,123) = - reac_rate_local(123) 
  reac_source_local(01,124) = + reac_rate_local(124) 
  reac_source_local(03,124) = + reac_rate_local(124) 
  reac_source_local(31,124) = - reac_rate_local(124) 
  reac_source_local(45,124) = + reac_rate_local(124) 
  reac_source_local(01,125) = + reac_rate_local(125) 
  reac_source_local(30,125) = + reac_rate_local(125) 
  reac_source_local(31,125) = - reac_rate_local(125) 
  reac_source_local(42,125) = + reac_rate_local(125) 
  reac_source_local(01,126) = + reac_rate_local(126) 
  reac_source_local(21,126) = - reac_rate_local(126) 
  reac_source_local(22,126) = + reac_rate_local(126) 
  reac_source_local(01,127) = + reac_rate_local(127) 
  reac_source_local(10,127) = + reac_rate_local(127) 
  reac_source_local(11,127) = + reac_rate_local(127) 
  reac_source_local(21,127) = - reac_rate_local(127) 
  reac_source_local(01,128) = + reac_rate_local(128) 
  reac_source_local(21,128) = - reac_rate_local(128) 
  reac_source_local(29,128) = + reac_rate_local(128) 
  reac_source_local(33,128) = + reac_rate_local(128) 
  reac_source_local(01,129) = + reac_rate_local(129) 
  reac_source_local(08,129) = + reac_rate_local(129) 
  reac_source_local(21,129) = - reac_rate_local(129) 
  reac_source_local(38,129) = + reac_rate_local(129) 
  reac_source_local(01,130) = + reac_rate_local(130) 
  reac_source_local(16,130) = + reac_rate_local(130) 
  reac_source_local(21,130) = - reac_rate_local(130) 
  reac_source_local(40,130) = + reac_rate_local(130) 
  reac_source_local(01,131) = + reac_rate_local(131) 
  reac_source_local(17,131) = + reac_rate_local(131) 
  reac_source_local(21,131) = - reac_rate_local(131) 
  reac_source_local(32,131) = + reac_rate_local(131) 
  reac_source_local(01,132) = + reac_rate_local(132) 
  reac_source_local(03,132) = + reac_rate_local(132) 
  reac_source_local(05,132) = + reac_rate_local(132) 
  reac_source_local(21,132) = - reac_rate_local(132) 
  reac_source_local(01,133) = + reac_rate_local(133) 
  reac_source_local(21,133) = - reac_rate_local(133) 
  reac_source_local(42,133) = + reac_rate_local(133) 
  reac_source_local(45,133) = + reac_rate_local(133) 
  reac_source_local(01,134) = + reac_rate_local(134) 
  reac_source_local(06,134) = - reac_rate_local(134) 
  reac_source_local(11,134) = + reac_rate_local(134) 
  reac_source_local(01,135) = + reac_rate_local(135) 
  reac_source_local(06,135) = - reac_rate_local(135) 
  reac_source_local(29,135) = + reac_rate_local(135) 
  reac_source_local(38,135) = + reac_rate_local(135) 
  reac_source_local(01,136) = + reac_rate_local(136) 
  reac_source_local(06,136) = - reac_rate_local(136) 
  reac_source_local(08,136) = + reac_rate_local(136) 
  reac_source_local(40,136) = + reac_rate_local(136) 
  reac_source_local(01,137) = + reac_rate_local(137) 
  reac_source_local(03,137) = + reac_rate_local(137) 
  reac_source_local(06,137) = - reac_rate_local(137) 
  reac_source_local(32,137) = + reac_rate_local(137) 
  reac_source_local(01,138) = + reac_rate_local(138) 
  reac_source_local(05,138) = + reac_rate_local(138) 
  reac_source_local(06,138) = - reac_rate_local(138) 
  reac_source_local(42,138) = + reac_rate_local(138) 
  reac_source_local(10,139) = + reac_rate_local(139) 
  reac_source_local(16,139) = + reac_rate_local(139) 
  reac_source_local(39,139) = - reac_rate_local(139) 
  reac_source_local(10,140) = + reac_rate_local(140) 
  reac_source_local(16,140) = + reac_rate_local(140) 
  reac_source_local(37,140) = - reac_rate_local(140) 
  reac_source_local(08,141) = + reac_rate_local(141) 
  reac_source_local(34,141) = + reac_rate_local(141) 
  reac_source_local(39,141) = - reac_rate_local(141) 
  reac_source_local(08,142) = + reac_rate_local(142) 
  reac_source_local(34,142) = + reac_rate_local(142) 
  reac_source_local(37,142) = - reac_rate_local(142) 
  reac_source_local(10,143) = + reac_rate_local(143) 
  reac_source_local(29,143) = + reac_rate_local(143) 
  reac_source_local(34,143) = + reac_rate_local(143) 
  reac_source_local(39,143) = - reac_rate_local(143) 
  reac_source_local(10,144) = + reac_rate_local(144) 
  reac_source_local(29,144) = + reac_rate_local(144) 
  reac_source_local(34,144) = + reac_rate_local(144) 
  reac_source_local(37,144) = - reac_rate_local(144) 
  reac_source_local(01,145) = + reac_rate_local(145) 
  reac_source_local(14,145) = + reac_rate_local(145) 
  reac_source_local(39,145) = - reac_rate_local(145) 
  reac_source_local(01,146) = + reac_rate_local(146) 
  reac_source_local(14,146) = + reac_rate_local(146) 
  reac_source_local(37,146) = - reac_rate_local(146) 
  reac_source_local(01,147) = + reac_rate_local(147) 
  reac_source_local(10,147) = + reac_rate_local(147) 
  reac_source_local(17,147) = + reac_rate_local(147) 
  reac_source_local(39,147) = - reac_rate_local(147) 
  reac_source_local(01,148) = + reac_rate_local(148) 
  reac_source_local(10,148) = + reac_rate_local(148) 
  reac_source_local(17,148) = + reac_rate_local(148) 
  reac_source_local(37,148) = - reac_rate_local(148) 
  reac_source_local(01,149) = + reac_rate_local(149) 
  reac_source_local(03,149) = + reac_rate_local(149) 
  reac_source_local(34,149) = + reac_rate_local(149) 
  reac_source_local(39,149) = - reac_rate_local(149) 
  reac_source_local(01,150) = + reac_rate_local(150) 
  reac_source_local(03,150) = + reac_rate_local(150) 
  reac_source_local(34,150) = + reac_rate_local(150) 
  reac_source_local(37,150) = - reac_rate_local(150) 
  reac_source_local(01,151) = + reac_rate_local(151) 
  reac_source_local(10,151) = + reac_rate_local(151) 
  reac_source_local(34,151) = + reac_rate_local(151) 
  reac_source_local(39,151) = - reac_rate_local(151) 
  reac_source_local(42,151) = + reac_rate_local(151) 
  reac_source_local(01,152) = + reac_rate_local(152) 
  reac_source_local(10,152) = + reac_rate_local(152) 
  reac_source_local(34,152) = + reac_rate_local(152) 
  reac_source_local(37,152) = - reac_rate_local(152) 
  reac_source_local(42,152) = + reac_rate_local(152) 
  reac_source_local(10,153) = + reac_rate_local(153) 
  reac_source_local(30,153) = + reac_rate_local(153) 
  reac_source_local(35,153) = - reac_rate_local(153) 
  reac_source_local(09,154) = - reac_rate_local(154) 
  reac_source_local(10,154) = + reac_rate_local(154) 
  reac_source_local(30,154) = + reac_rate_local(154) 
  reac_source_local(34,155) = + reac_rate_local(155) 
  reac_source_local(35,155) = - reac_rate_local(155) 
  reac_source_local(45,155) = + reac_rate_local(155) 
  reac_source_local(09,156) = - reac_rate_local(156) 
  reac_source_local(34,156) = + reac_rate_local(156) 
  reac_source_local(45,156) = + reac_rate_local(156) 
  reac_source_local(05,157) = + reac_rate_local(157) 
  reac_source_local(10,157) = + reac_rate_local(157) 
  reac_source_local(34,157) = + reac_rate_local(157) 
  reac_source_local(35,157) = - reac_rate_local(157) 
  reac_source_local(05,158) = + reac_rate_local(158) 
  reac_source_local(09,158) = - reac_rate_local(158) 
  reac_source_local(10,158) = + reac_rate_local(158) 
  reac_source_local(34,158) = + reac_rate_local(158) 
  reac_source_local(32,159) = + reac_rate_local(159) 
  reac_source_local(34,159) = + reac_rate_local(159) * 2.d0
  reac_source_local(35,159) = - reac_rate_local(159) 
  reac_source_local(09,160) = - reac_rate_local(160) 
  reac_source_local(32,160) = + reac_rate_local(160) 
  reac_source_local(34,160) = + reac_rate_local(160) * 2.d0
  reac_source_local(08,161) = + reac_rate_local(161) 
  reac_source_local(35,161) = - reac_rate_local(161) 
  reac_source_local(44,161) = + reac_rate_local(161) 
  reac_source_local(08,162) = + reac_rate_local(162) 
  reac_source_local(09,162) = - reac_rate_local(162) 
  reac_source_local(44,162) = + reac_rate_local(162) 
  reac_source_local(16,163) = + reac_rate_local(163) * 2.d0
  reac_source_local(35,163) = - reac_rate_local(163) 
  reac_source_local(09,164) = - reac_rate_local(164) 
  reac_source_local(16,164) = + reac_rate_local(164) * 2.d0
  reac_source_local(05,165) = + reac_rate_local(165) 
  reac_source_local(10,165) = + reac_rate_local(165) 
  reac_source_local(15,165) = - reac_rate_local(165) 
  reac_source_local(02,166) = - reac_rate_local(166) 
  reac_source_local(05,166) = + reac_rate_local(166) 
  reac_source_local(10,166) = + reac_rate_local(166) 
  reac_source_local(15,167) = - reac_rate_local(167) 
  reac_source_local(32,167) = + reac_rate_local(167) 
  reac_source_local(34,167) = + reac_rate_local(167) 
  reac_source_local(02,168) = - reac_rate_local(168) 
  reac_source_local(32,168) = + reac_rate_local(168) 
  reac_source_local(34,168) = + reac_rate_local(168) 
  reac_source_local(10,169) = + reac_rate_local(169) * 2.d0
  reac_source_local(15,169) = - reac_rate_local(169) 
  reac_source_local(32,169) = + reac_rate_local(169) 
  reac_source_local(02,170) = - reac_rate_local(170) 
  reac_source_local(10,170) = + reac_rate_local(170) * 2.d0
  reac_source_local(32,170) = + reac_rate_local(170) 
  reac_source_local(15,171) = - reac_rate_local(171) 
  reac_source_local(16,171) = + reac_rate_local(171) 
  reac_source_local(29,171) = + reac_rate_local(171) 
  reac_source_local(02,172) = - reac_rate_local(172) 
  reac_source_local(16,172) = + reac_rate_local(172) 
  reac_source_local(29,172) = + reac_rate_local(172) 
  reac_source_local(08,173) = + reac_rate_local(173) * 2.d0
  reac_source_local(15,173) = - reac_rate_local(173) 
  reac_source_local(02,174) = - reac_rate_local(174) 
  reac_source_local(08,174) = + reac_rate_local(174) * 2.d0
  reac_source_local(26,175) = - reac_rate_local(175) 
  reac_source_local(29,175) = + reac_rate_local(175) * 2.d0
  reac_source_local(29,176) = + reac_rate_local(176) * 2.d0
  reac_source_local(36,176) = - reac_rate_local(176) 
  reac_source_local(28,177) = - reac_rate_local(177) 
  reac_source_local(29,177) = + reac_rate_local(177) * 2.d0
  reac_source_local(04,178) = - reac_rate_local(178) 
  reac_source_local(10,178) = + reac_rate_local(178) 
  reac_source_local(27,178) = + reac_rate_local(178) 
  reac_source_local(10,179) = + reac_rate_local(179) 
  reac_source_local(24,179) = - reac_rate_local(179) 
  reac_source_local(27,179) = + reac_rate_local(179) 
  reac_source_local(04,180) = - reac_rate_local(180) 
  reac_source_local(31,180) = + reac_rate_local(180) 
  reac_source_local(34,180) = + reac_rate_local(180) 
  reac_source_local(24,181) = - reac_rate_local(181) 
  reac_source_local(31,181) = + reac_rate_local(181) 
  reac_source_local(34,181) = + reac_rate_local(181) 
  reac_source_local(04,182) = - reac_rate_local(182) 
  reac_source_local(06,182) = + reac_rate_local(182) 
  reac_source_local(34,182) = + reac_rate_local(182) * 2.d0
  reac_source_local(06,183) = + reac_rate_local(183) 
  reac_source_local(24,183) = - reac_rate_local(183) 
  reac_source_local(34,183) = + reac_rate_local(183) * 2.d0
  reac_source_local(04,184) = - reac_rate_local(184) 
  reac_source_local(08,184) = + reac_rate_local(184) 
  reac_source_local(18,184) = + reac_rate_local(184) 
  reac_source_local(08,185) = + reac_rate_local(185) 
  reac_source_local(18,185) = + reac_rate_local(185) 
  reac_source_local(24,185) = - reac_rate_local(185) 
  reac_source_local(04,186) = - reac_rate_local(186) 
  reac_source_local(16,186) = + reac_rate_local(186) 
  reac_source_local(30,186) = + reac_rate_local(186) 
  reac_source_local(16,187) = + reac_rate_local(187) 
  reac_source_local(24,187) = - reac_rate_local(187) 
  reac_source_local(30,187) = + reac_rate_local(187) 
  reac_source_local(04,188) = - reac_rate_local(188) 
  reac_source_local(44,188) = + reac_rate_local(188) 
  reac_source_local(45,188) = + reac_rate_local(188) 
  reac_source_local(24,189) = - reac_rate_local(189) 
  reac_source_local(44,189) = + reac_rate_local(189) 
  reac_source_local(45,189) = + reac_rate_local(189) 
  reac_source_local(01,190) = + reac_rate_local(190) 
  reac_source_local(12,190) = + reac_rate_local(190) 
  reac_source_local(35,190) = - reac_rate_local(190) 
  reac_source_local(01,191) = + reac_rate_local(191) 
  reac_source_local(09,191) = - reac_rate_local(191) 
  reac_source_local(12,191) = + reac_rate_local(191) 
  reac_source_local(01,192) = + reac_rate_local(192) 
  reac_source_local(10,192) = + reac_rate_local(192) 
  reac_source_local(35,192) = - reac_rate_local(192) 
  reac_source_local(41,192) = + reac_rate_local(192) 
  reac_source_local(01,193) = + reac_rate_local(193) 
  reac_source_local(09,193) = - reac_rate_local(193) 
  reac_source_local(10,193) = + reac_rate_local(193) 
  reac_source_local(41,193) = + reac_rate_local(193) 
  reac_source_local(01,194) = + reac_rate_local(194) 
  reac_source_local(33,194) = + reac_rate_local(194) 
  reac_source_local(34,194) = + reac_rate_local(194) 
  reac_source_local(35,194) = - reac_rate_local(194) 
  reac_source_local(01,195) = + reac_rate_local(195) 
  reac_source_local(09,195) = - reac_rate_local(195) 
  reac_source_local(33,195) = + reac_rate_local(195) 
  reac_source_local(34,195) = + reac_rate_local(195) 
  reac_source_local(01,196) = + reac_rate_local(196) 
  reac_source_local(10,196) = + reac_rate_local(196) 
  reac_source_local(34,196) = + reac_rate_local(196) 
  reac_source_local(35,196) = - reac_rate_local(196) 
  reac_source_local(38,196) = + reac_rate_local(196) 
  reac_source_local(01,197) = + reac_rate_local(197) 
  reac_source_local(09,197) = - reac_rate_local(197) 
  reac_source_local(10,197) = + reac_rate_local(197) 
  reac_source_local(34,197) = + reac_rate_local(197) 
  reac_source_local(38,197) = + reac_rate_local(197) 
  reac_source_local(01,198) = + reac_rate_local(198) 
  reac_source_local(34,198) = + reac_rate_local(198) * 2.d0
  reac_source_local(35,198) = - reac_rate_local(198) 
  reac_source_local(40,198) = + reac_rate_local(198) 
  reac_source_local(01,199) = + reac_rate_local(199) 
  reac_source_local(09,199) = - reac_rate_local(199) 
  reac_source_local(34,199) = + reac_rate_local(199) * 2.d0
  reac_source_local(40,199) = + reac_rate_local(199) 
  reac_source_local(01,200) = + reac_rate_local(200) 
  reac_source_local(16,200) = + reac_rate_local(200) 
  reac_source_local(17,200) = + reac_rate_local(200) 
  reac_source_local(35,200) = - reac_rate_local(200) 
  reac_source_local(01,201) = + reac_rate_local(201) 
  reac_source_local(09,201) = - reac_rate_local(201) 
  reac_source_local(16,201) = + reac_rate_local(201) 
  reac_source_local(17,201) = + reac_rate_local(201) 
  reac_source_local(01,202) = + reac_rate_local(202) 
  reac_source_local(03,202) = + reac_rate_local(202) 
  reac_source_local(35,202) = - reac_rate_local(202) 
  reac_source_local(44,202) = + reac_rate_local(202) 
  reac_source_local(01,203) = + reac_rate_local(203) 
  reac_source_local(03,203) = + reac_rate_local(203) 
  reac_source_local(09,203) = - reac_rate_local(203) 
  reac_source_local(44,203) = + reac_rate_local(203) 
  reac_source_local(01,204) = + reac_rate_local(204) 
  reac_source_local(15,204) = - reac_rate_local(204) 
  reac_source_local(33,204) = + reac_rate_local(204) 
  reac_source_local(01,205) = + reac_rate_local(205) 
  reac_source_local(02,205) = - reac_rate_local(205) 
  reac_source_local(33,205) = + reac_rate_local(205) 
  reac_source_local(01,206) = + reac_rate_local(206) 
  reac_source_local(10,206) = + reac_rate_local(206) 
  reac_source_local(15,206) = - reac_rate_local(206) 
  reac_source_local(38,206) = + reac_rate_local(206) 
  reac_source_local(01,207) = + reac_rate_local(207) 
  reac_source_local(02,207) = - reac_rate_local(207) 
  reac_source_local(10,207) = + reac_rate_local(207) 
  reac_source_local(38,207) = + reac_rate_local(207) 
  reac_source_local(01,208) = + reac_rate_local(208) 
  reac_source_local(15,208) = - reac_rate_local(208) 
  reac_source_local(17,208) = + reac_rate_local(208) 
  reac_source_local(29,208) = + reac_rate_local(208) 
  reac_source_local(01,209) = + reac_rate_local(209) 
  reac_source_local(02,209) = - reac_rate_local(209) 
  reac_source_local(17,209) = + reac_rate_local(209) 
  reac_source_local(29,209) = + reac_rate_local(209) 
  reac_source_local(01,210) = + reac_rate_local(210) 
  reac_source_local(03,210) = + reac_rate_local(210) 
  reac_source_local(08,210) = + reac_rate_local(210) 
  reac_source_local(15,210) = - reac_rate_local(210) 
  reac_source_local(01,211) = + reac_rate_local(211) 
  reac_source_local(02,211) = - reac_rate_local(211) 
  reac_source_local(03,211) = + reac_rate_local(211) 
  reac_source_local(08,211) = + reac_rate_local(211) 
  reac_source_local(01,212) = + reac_rate_local(212) 
  reac_source_local(15,212) = - reac_rate_local(212) 
  reac_source_local(16,212) = + reac_rate_local(212) 
  reac_source_local(42,212) = + reac_rate_local(212) 
  reac_source_local(01,213) = + reac_rate_local(213) 
  reac_source_local(02,213) = - reac_rate_local(213) 
  reac_source_local(16,213) = + reac_rate_local(213) 
  reac_source_local(42,213) = + reac_rate_local(213) 
  reac_source_local(01,214) = + reac_rate_local(214) 
  reac_source_local(26,214) = - reac_rate_local(214) 
  reac_source_local(40,214) = + reac_rate_local(214) 
  reac_source_local(01,215) = + reac_rate_local(215) 
  reac_source_local(36,215) = - reac_rate_local(215) 
  reac_source_local(40,215) = + reac_rate_local(215) 
  reac_source_local(01,216) = + reac_rate_local(216) 
  reac_source_local(28,216) = - reac_rate_local(216) 
  reac_source_local(40,216) = + reac_rate_local(216) 
  reac_source_local(01,217) = + reac_rate_local(217) 
  reac_source_local(26,217) = - reac_rate_local(217) 
  reac_source_local(29,217) = + reac_rate_local(217) 
  reac_source_local(42,217) = + reac_rate_local(217) 
  reac_source_local(01,218) = + reac_rate_local(218) 
  reac_source_local(29,218) = + reac_rate_local(218) 
  reac_source_local(36,218) = - reac_rate_local(218) 
  reac_source_local(42,218) = + reac_rate_local(218) 
  reac_source_local(01,219) = + reac_rate_local(219) 
  reac_source_local(28,219) = - reac_rate_local(219) 
  reac_source_local(29,219) = + reac_rate_local(219) 
  reac_source_local(42,219) = + reac_rate_local(219) 
  reac_source_local(01,220) = + reac_rate_local(220) 
  reac_source_local(04,220) = - reac_rate_local(220) 
  reac_source_local(25,220) = + reac_rate_local(220) 
  reac_source_local(01,221) = + reac_rate_local(221) 
  reac_source_local(24,221) = - reac_rate_local(221) 
  reac_source_local(25,221) = + reac_rate_local(221) 
  reac_source_local(01,222) = + reac_rate_local(222) 
  reac_source_local(04,222) = - reac_rate_local(222) 
  reac_source_local(10,222) = + reac_rate_local(222) 
  reac_source_local(23,222) = + reac_rate_local(222) 
  reac_source_local(01,223) = + reac_rate_local(223) 
  reac_source_local(10,223) = + reac_rate_local(223) 
  reac_source_local(23,223) = + reac_rate_local(223) 
  reac_source_local(24,223) = - reac_rate_local(223) 
  reac_source_local(01,224) = + reac_rate_local(224) 
  reac_source_local(04,224) = - reac_rate_local(224) 
  reac_source_local(19,224) = + reac_rate_local(224) 
  reac_source_local(34,224) = + reac_rate_local(224) 
  reac_source_local(01,225) = + reac_rate_local(225) 
  reac_source_local(19,225) = + reac_rate_local(225) 
  reac_source_local(24,225) = - reac_rate_local(225) 
  reac_source_local(34,225) = + reac_rate_local(225) 
  reac_source_local(01,226) = + reac_rate_local(226) 
  reac_source_local(04,226) = - reac_rate_local(226) 
  reac_source_local(10,226) = + reac_rate_local(226) 
  reac_source_local(22,226) = + reac_rate_local(226) 
  reac_source_local(34,226) = + reac_rate_local(226) 
  reac_source_local(01,227) = + reac_rate_local(227) 
  reac_source_local(10,227) = + reac_rate_local(227) 
  reac_source_local(22,227) = + reac_rate_local(227) 
  reac_source_local(24,227) = - reac_rate_local(227) 
  reac_source_local(34,227) = + reac_rate_local(227) 
  reac_source_local(01,228) = + reac_rate_local(228) 
  reac_source_local(04,228) = - reac_rate_local(228) 
  reac_source_local(11,228) = + reac_rate_local(228) 
  reac_source_local(34,228) = + reac_rate_local(228) * 2.d0
  reac_source_local(01,229) = + reac_rate_local(229) 
  reac_source_local(11,229) = + reac_rate_local(229) 
  reac_source_local(24,229) = - reac_rate_local(229) 
  reac_source_local(34,229) = + reac_rate_local(229) * 2.d0
  reac_source_local(01,230) = + reac_rate_local(230) 
  reac_source_local(04,230) = - reac_rate_local(230) 
  reac_source_local(16,230) = + reac_rate_local(230) 
  reac_source_local(41,230) = + reac_rate_local(230) 
  reac_source_local(01,231) = + reac_rate_local(231) 
  reac_source_local(16,231) = + reac_rate_local(231) 
  reac_source_local(24,231) = - reac_rate_local(231) 
  reac_source_local(41,231) = + reac_rate_local(231) 
  reac_source_local(01,232) = + reac_rate_local(232) 
  reac_source_local(04,232) = - reac_rate_local(232) 
  reac_source_local(33,232) = + reac_rate_local(232) 
  reac_source_local(44,232) = + reac_rate_local(232) 
  reac_source_local(01,233) = + reac_rate_local(233) 
  reac_source_local(24,233) = - reac_rate_local(233) 
  reac_source_local(33,233) = + reac_rate_local(233) 
  reac_source_local(44,233) = + reac_rate_local(233) 
  reac_source_local(01,234) = + reac_rate_local(234) 
  reac_source_local(04,234) = - reac_rate_local(234) 
  reac_source_local(17,234) = + reac_rate_local(234) 
  reac_source_local(30,234) = + reac_rate_local(234) 
  reac_source_local(01,235) = + reac_rate_local(235) 
  reac_source_local(17,235) = + reac_rate_local(235) 
  reac_source_local(24,235) = - reac_rate_local(235) 
  reac_source_local(30,235) = + reac_rate_local(235) 
  reac_source_local(01,236) = + reac_rate_local(236) 
  reac_source_local(03,236) = + reac_rate_local(236) 
  reac_source_local(04,236) = - reac_rate_local(236) 
  reac_source_local(18,236) = + reac_rate_local(236) 
  reac_source_local(01,237) = + reac_rate_local(237) 
  reac_source_local(03,237) = + reac_rate_local(237) 
  reac_source_local(18,237) = + reac_rate_local(237) 
  reac_source_local(24,237) = - reac_rate_local(237) 
  reac_source_local(10,238) = + reac_rate_local(238) * 2.d0
  reac_source_local(34,238) = - reac_rate_local(238) 
  reac_source_local(01,239) = + reac_rate_local(239) 
  reac_source_local(34,239) = - reac_rate_local(239) 
  reac_source_local(43,239) = + reac_rate_local(239) 
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
  if( lgas_heating ) ZDPlasKin_cfg(1) = y(46)
  density(:) = y(1:species_max)
  call ZDPlasKin_reac_rates(t)
  rrt(001) = rrt(001) * density(01) * density(44) 
  rrt(002) = rrt(002) * density(01) * density(44) 
  rrt(003) = rrt(003) * density(01) * density(18) 
  rrt(004) = rrt(004) * density(01) * density(18) 
  rrt(005) = rrt(005) * density(01) * density(45) 
  rrt(006) = rrt(006) * density(01) * density(45) 
  rrt(007) = rrt(007) * density(01) * density(32) 
  rrt(008) = rrt(008) * density(01) * density(32) 
  rrt(009) = rrt(009) * density(01) * density(32) 
  rrt(010) = rrt(010) * density(01) * density(20) 
  rrt(011) = rrt(011) * density(01) * density(20) 
  rrt(012) = rrt(012) * density(01) * density(31) 
  rrt(013) = rrt(013) * density(01) * density(44) 
  rrt(014) = rrt(014) * density(01) * density(44) 
  rrt(015) = rrt(015) * density(01) * density(44) 
  rrt(016) = rrt(016) * density(01) * density(16) 
  rrt(017) = rrt(017) * density(01) * density(16) 
  rrt(018) = rrt(018) * density(01) * density(08) 
  rrt(019) = rrt(019) * density(01) * density(44) 
  rrt(020) = rrt(020) * density(01) * density(44) 
  rrt(021) = rrt(021) * density(01) * density(44) 
  rrt(022) = rrt(022) * density(01) * density(44) 
  rrt(023) = rrt(023) * density(01) * density(16) 
  rrt(024) = rrt(024) * density(01) * density(16) 
  rrt(025) = rrt(025) * density(01) * density(16) 
  rrt(026) = rrt(026) * density(01) * density(08) 
  rrt(027) = rrt(027) * density(01) * density(08) 
  rrt(028) = rrt(028) * density(01) * density(29) 
  rrt(029) = rrt(029) * density(01) * density(18) 
  rrt(030) = rrt(030) * density(01) * density(18) 
  rrt(031) = rrt(031) * density(01) * density(18) 
  rrt(032) = rrt(032) * density(01) * density(18) 
  rrt(033) = rrt(033) * density(01) * density(18) 
  rrt(034) = rrt(034) * density(01) * density(18) 
  rrt(035) = rrt(035) * density(01) * density(30) 
  rrt(036) = rrt(036) * density(01) * density(30) 
  rrt(037) = rrt(037) * density(01) * density(30) 
  rrt(038) = rrt(038) * density(01) * density(30) 
  rrt(039) = rrt(039) * density(01) * density(30) 
  rrt(040) = rrt(040) * density(01) * density(30) 
  rrt(041) = rrt(041) * density(01) * density(45) 
  rrt(042) = rrt(042) * density(01) * density(45) 
  rrt(043) = rrt(043) * density(01) * density(45) 
  rrt(044) = rrt(044) * density(01) * density(45) 
  rrt(045) = rrt(045) * density(01) * density(45) 
  rrt(046) = rrt(046) * density(01) * density(05) 
  rrt(047) = rrt(047) * density(01) * density(05) 
  rrt(048) = rrt(048) * density(01) * density(32) 
  rrt(049) = rrt(049) * density(01) * density(20) 
  rrt(050) = rrt(050) * density(01) * density(20) 
  rrt(051) = rrt(051) * density(01) * density(20) 
  rrt(052) = rrt(052) * density(01) * density(20) 
  rrt(053) = rrt(053) * density(01) * density(20) 
  rrt(054) = rrt(054) * density(01) * density(20) 
  rrt(055) = rrt(055) * density(01) * density(27) 
  rrt(056) = rrt(056) * density(01) * density(27) 
  rrt(057) = rrt(057) * density(01) * density(27) 
  rrt(058) = rrt(058) * density(01) * density(27) 
  rrt(059) = rrt(059) * density(01) * density(27) 
  rrt(060) = rrt(060) * density(01) * density(31) 
  rrt(061) = rrt(061) * density(01) * density(31) 
  rrt(062) = rrt(062) * density(01) * density(31) 
  rrt(063) = rrt(063) * density(01) * density(31) 
  rrt(064) = rrt(064) * density(01) * density(31) 
  rrt(065) = rrt(065) * density(01) * density(21) 
  rrt(066) = rrt(066) * density(01) * density(21) 
  rrt(067) = rrt(067) * density(01) * density(06) 
  rrt(068) = rrt(068) * density(01) * density(06) 
  rrt(069) = rrt(069) * density(01) * density(18) 
  rrt(070) = rrt(070) * density(01) * density(18) 
  rrt(071) = rrt(071) * density(01) * density(18) 
  rrt(072) = rrt(072) * density(01) * density(18) 
  rrt(073) = rrt(073) * density(01) * density(18) 
  rrt(074) = rrt(074) * density(01) * density(18) 
  rrt(075) = rrt(075) * density(01) * density(18) 
  rrt(076) = rrt(076) * density(01) * density(30) 
  rrt(077) = rrt(077) * density(01) * density(30) 
  rrt(078) = rrt(078) * density(01) * density(30) 
  rrt(079) = rrt(079) * density(01) * density(30) 
  rrt(080) = rrt(080) * density(01) * density(30) 
  rrt(081) = rrt(081) * density(01) * density(30) 
  rrt(082) = rrt(082) * density(01) * density(30) 
  rrt(083) = rrt(083) * density(01) * density(45) 
  rrt(084) = rrt(084) * density(01) * density(45) 
  rrt(085) = rrt(085) * density(01) * density(45) 
  rrt(086) = rrt(086) * density(01) * density(45) 
  rrt(087) = rrt(087) * density(01) * density(45) 
  rrt(088) = rrt(088) * density(01) * density(05) 
  rrt(089) = rrt(089) * density(01) * density(05) 
  rrt(090) = rrt(090) * density(01) * density(05) 
  rrt(091) = rrt(091) * density(01) * density(05) 
  rrt(092) = rrt(092) * density(01) * density(05) 
  rrt(093) = rrt(093) * density(01) * density(32) 
  rrt(094) = rrt(094) * density(01) * density(32) 
  rrt(095) = rrt(095) * density(01) * density(20) 
  rrt(096) = rrt(096) * density(01) * density(20) 
  rrt(097) = rrt(097) * density(01) * density(20) 
  rrt(098) = rrt(098) * density(01) * density(20) 
  rrt(099) = rrt(099) * density(01) * density(20) 
  rrt(100) = rrt(100) * density(01) * density(20) 
  rrt(101) = rrt(101) * density(01) * density(20) 
  rrt(102) = rrt(102) * density(01) * density(20) 
  rrt(103) = rrt(103) * density(01) * density(20) 
  rrt(104) = rrt(104) * density(01) * density(27) 
  rrt(105) = rrt(105) * density(01) * density(27) 
  rrt(106) = rrt(106) * density(01) * density(27) 
  rrt(107) = rrt(107) * density(01) * density(27) 
  rrt(108) = rrt(108) * density(01) * density(27) 
  rrt(109) = rrt(109) * density(01) * density(27) 
  rrt(110) = rrt(110) * density(01) * density(27) 
  rrt(111) = rrt(111) * density(01) * density(27) 
  rrt(112) = rrt(112) * density(01) * density(27) 
  rrt(113) = rrt(113) * density(01) * density(27) 
  rrt(114) = rrt(114) * density(01) * density(27) 
  rrt(115) = rrt(115) * density(01) * density(31) 
  rrt(116) = rrt(116) * density(01) * density(31) 
  rrt(117) = rrt(117) * density(01) * density(31) 
  rrt(118) = rrt(118) * density(01) * density(31) 
  rrt(119) = rrt(119) * density(01) * density(31) 
  rrt(120) = rrt(120) * density(01) * density(31) 
  rrt(121) = rrt(121) * density(01) * density(31) 
  rrt(122) = rrt(122) * density(01) * density(31) 
  rrt(123) = rrt(123) * density(01) * density(31) 
  rrt(124) = rrt(124) * density(01) * density(31) 
  rrt(125) = rrt(125) * density(01) * density(31) 
  rrt(126) = rrt(126) * density(01) * density(21) 
  rrt(127) = rrt(127) * density(01) * density(21) 
  rrt(128) = rrt(128) * density(01) * density(21) 
  rrt(129) = rrt(129) * density(01) * density(21) 
  rrt(130) = rrt(130) * density(01) * density(21) 
  rrt(131) = rrt(131) * density(01) * density(21) 
  rrt(132) = rrt(132) * density(01) * density(21) 
  rrt(133) = rrt(133) * density(01) * density(21) 
  rrt(134) = rrt(134) * density(01) * density(06) 
  rrt(135) = rrt(135) * density(01) * density(06) 
  rrt(136) = rrt(136) * density(01) * density(06) 
  rrt(137) = rrt(137) * density(01) * density(06) 
  rrt(138) = rrt(138) * density(01) * density(06) 
  rrt(139) = rrt(139) * density(01) * density(39) 
  rrt(140) = rrt(140) * density(01) * density(37) 
  rrt(141) = rrt(141) * density(01) * density(39) 
  rrt(142) = rrt(142) * density(01) * density(37) 
  rrt(143) = rrt(143) * density(01) * density(39) 
  rrt(144) = rrt(144) * density(01) * density(37) 
  rrt(145) = rrt(145) * density(01) * density(39) 
  rrt(146) = rrt(146) * density(01) * density(37) 
  rrt(147) = rrt(147) * density(01) * density(39) 
  rrt(148) = rrt(148) * density(01) * density(37) 
  rrt(149) = rrt(149) * density(01) * density(39) 
  rrt(150) = rrt(150) * density(01) * density(37) 
  rrt(151) = rrt(151) * density(01) * density(39) 
  rrt(152) = rrt(152) * density(01) * density(37) 
  rrt(153) = rrt(153) * density(01) * density(35) 
  rrt(154) = rrt(154) * density(01) * density(09) 
  rrt(155) = rrt(155) * density(01) * density(35) 
  rrt(156) = rrt(156) * density(01) * density(09) 
  rrt(157) = rrt(157) * density(01) * density(35) 
  rrt(158) = rrt(158) * density(01) * density(09) 
  rrt(159) = rrt(159) * density(01) * density(35) 
  rrt(160) = rrt(160) * density(01) * density(09) 
  rrt(161) = rrt(161) * density(01) * density(35) 
  rrt(162) = rrt(162) * density(01) * density(09) 
  rrt(163) = rrt(163) * density(01) * density(35) 
  rrt(164) = rrt(164) * density(01) * density(09) 
  rrt(165) = rrt(165) * density(01) * density(15) 
  rrt(166) = rrt(166) * density(01) * density(02) 
  rrt(167) = rrt(167) * density(01) * density(15) 
  rrt(168) = rrt(168) * density(01) * density(02) 
  rrt(169) = rrt(169) * density(01) * density(15) 
  rrt(170) = rrt(170) * density(01) * density(02) 
  rrt(171) = rrt(171) * density(01) * density(15) 
  rrt(172) = rrt(172) * density(01) * density(02) 
  rrt(173) = rrt(173) * density(01) * density(15) 
  rrt(174) = rrt(174) * density(01) * density(02) 
  rrt(175) = rrt(175) * density(01) * density(26) 
  rrt(176) = rrt(176) * density(01) * density(36) 
  rrt(177) = rrt(177) * density(01) * density(28) 
  rrt(178) = rrt(178) * density(01) * density(04) 
  rrt(179) = rrt(179) * density(01) * density(24) 
  rrt(180) = rrt(180) * density(01) * density(04) 
  rrt(181) = rrt(181) * density(01) * density(24) 
  rrt(182) = rrt(182) * density(01) * density(04) 
  rrt(183) = rrt(183) * density(01) * density(24) 
  rrt(184) = rrt(184) * density(01) * density(04) 
  rrt(185) = rrt(185) * density(01) * density(24) 
  rrt(186) = rrt(186) * density(01) * density(04) 
  rrt(187) = rrt(187) * density(01) * density(24) 
  rrt(188) = rrt(188) * density(01) * density(04) 
  rrt(189) = rrt(189) * density(01) * density(24) 
  rrt(190) = rrt(190) * density(01) * density(35) 
  rrt(191) = rrt(191) * density(01) * density(09) 
  rrt(192) = rrt(192) * density(01) * density(35) 
  rrt(193) = rrt(193) * density(01) * density(09) 
  rrt(194) = rrt(194) * density(01) * density(35) 
  rrt(195) = rrt(195) * density(01) * density(09) 
  rrt(196) = rrt(196) * density(01) * density(35) 
  rrt(197) = rrt(197) * density(01) * density(09) 
  rrt(198) = rrt(198) * density(01) * density(35) 
  rrt(199) = rrt(199) * density(01) * density(09) 
  rrt(200) = rrt(200) * density(01) * density(35) 
  rrt(201) = rrt(201) * density(01) * density(09) 
  rrt(202) = rrt(202) * density(01) * density(35) 
  rrt(203) = rrt(203) * density(01) * density(09) 
  rrt(204) = rrt(204) * density(01) * density(15) 
  rrt(205) = rrt(205) * density(01) * density(02) 
  rrt(206) = rrt(206) * density(01) * density(15) 
  rrt(207) = rrt(207) * density(01) * density(02) 
  rrt(208) = rrt(208) * density(01) * density(15) 
  rrt(209) = rrt(209) * density(01) * density(02) 
  rrt(210) = rrt(210) * density(01) * density(15) 
  rrt(211) = rrt(211) * density(01) * density(02) 
  rrt(212) = rrt(212) * density(01) * density(15) 
  rrt(213) = rrt(213) * density(01) * density(02) 
  rrt(214) = rrt(214) * density(01) * density(26) 
  rrt(215) = rrt(215) * density(01) * density(36) 
  rrt(216) = rrt(216) * density(01) * density(28) 
  rrt(217) = rrt(217) * density(01) * density(26) 
  rrt(218) = rrt(218) * density(01) * density(36) 
  rrt(219) = rrt(219) * density(01) * density(28) 
  rrt(220) = rrt(220) * density(01) * density(04) 
  rrt(221) = rrt(221) * density(01) * density(24) 
  rrt(222) = rrt(222) * density(01) * density(04) 
  rrt(223) = rrt(223) * density(01) * density(24) 
  rrt(224) = rrt(224) * density(01) * density(04) 
  rrt(225) = rrt(225) * density(01) * density(24) 
  rrt(226) = rrt(226) * density(01) * density(04) 
  rrt(227) = rrt(227) * density(01) * density(24) 
  rrt(228) = rrt(228) * density(01) * density(04) 
  rrt(229) = rrt(229) * density(01) * density(24) 
  rrt(230) = rrt(230) * density(01) * density(04) 
  rrt(231) = rrt(231) * density(01) * density(24) 
  rrt(232) = rrt(232) * density(01) * density(04) 
  rrt(233) = rrt(233) * density(01) * density(24) 
  rrt(234) = rrt(234) * density(01) * density(04) 
  rrt(235) = rrt(235) * density(01) * density(24) 
  rrt(236) = rrt(236) * density(01) * density(04) 
  rrt(237) = rrt(237) * density(01) * density(24) 
  rrt(238) = rrt(238) * density(01) * density(34) 
  rrt(239) = rrt(239) * density(01) * density(34) 
  ydot(01) = +rrt(019)+rrt(020)+rrt(021)+rrt(022)+rrt(023)+rrt(024)+rrt(025)+rrt(026)+rrt(027)+rrt(028)+rrt(069)+rrt(070)+rrt(071)&
             +rrt(072)+rrt(073)+rrt(074)+rrt(075)+rrt(076)+rrt(077)+rrt(078)+rrt(079)+rrt(080)+rrt(081)+rrt(082)+rrt(083)+rrt(084)&
             +rrt(085)+rrt(086)+rrt(087)+rrt(088)+rrt(089)+rrt(090)+rrt(091)+rrt(092)+rrt(093)+rrt(094)+rrt(095)+rrt(096)+rrt(097)&
             +rrt(098)+rrt(099)+rrt(100)+rrt(101)+rrt(102)+rrt(103)+rrt(104)+rrt(105)+rrt(106)+rrt(107)+rrt(108)+rrt(109)+rrt(110)&
             +rrt(111)+rrt(112)+rrt(113)+rrt(114)+rrt(115)+rrt(116)+rrt(117)+rrt(118)+rrt(119)+rrt(120)+rrt(121)+rrt(122)+rrt(123)&
             +rrt(124)+rrt(125)+rrt(126)+rrt(127)+rrt(128)+rrt(129)+rrt(130)+rrt(131)+rrt(132)+rrt(133)+rrt(134)+rrt(135)+rrt(136)&
             +rrt(137)+rrt(138)+rrt(145)+rrt(146)+rrt(147)+rrt(148)+rrt(149)+rrt(150)+rrt(151)+rrt(152)+rrt(190)+rrt(191)+rrt(192)&
             +rrt(193)+rrt(194)+rrt(195)+rrt(196)+rrt(197)+rrt(198)+rrt(199)+rrt(200)+rrt(201)+rrt(202)+rrt(203)+rrt(204)+rrt(205)&
             +rrt(206)+rrt(207)+rrt(208)+rrt(209)+rrt(210)+rrt(211)+rrt(212)+rrt(213)+rrt(214)+rrt(215)+rrt(216)+rrt(217)+rrt(218)&
             +rrt(219)+rrt(220)+rrt(221)+rrt(222)+rrt(223)+rrt(224)+rrt(225)+rrt(226)+rrt(227)+rrt(228)+rrt(229)+rrt(230)+rrt(231)&
             +rrt(232)+rrt(233)+rrt(234)+rrt(235)+rrt(236)+rrt(237)+rrt(239) 
  ydot(02) = +rrt(006)-rrt(166)-rrt(168)-rrt(170)-rrt(172)-rrt(174)-rrt(205)-rrt(207)-rrt(209)-rrt(211)-rrt(213) 
  ydot(03) = +rrt(021)+rrt(024)+rrt(026)+rrt(075)+rrt(081)+rrt(086)+rrt(090)+rrt(103)+rrt(113)+rrt(124)+rrt(132)+rrt(137)+rrt(149)&
             +rrt(150)+rrt(202)+rrt(203)+rrt(210)+rrt(211)+rrt(236)+rrt(237) 
  ydot(04) = +rrt(010)-rrt(178)-rrt(180)-rrt(182)-rrt(184)-rrt(186)-rrt(188)-rrt(220)-rrt(222)-rrt(224)-rrt(226)-rrt(228)-rrt(230)&
             -rrt(232)-rrt(234)-rrt(236) 
  ydot(05) = +rrt(031)+rrt(036)+rrt(037)+rrt(041)-rrt(046)-rrt(047)+rrt(059)+rrt(063)+rrt(067)-rrt(088)-rrt(089)-rrt(090)-rrt(091)&
             -rrt(092)+rrt(111)+rrt(123)+rrt(132)+rrt(138)+rrt(157)+rrt(158)+rrt(165)+rrt(166) 
  ydot(06) = +rrt(051)+rrt(057)+rrt(061)+rrt(065)-rrt(067)-rrt(068)-rrt(134)-rrt(135)-rrt(136)-rrt(137)-rrt(138)+rrt(182)+rrt(183) 
  ydot(07) = +rrt(092) 
  ydot(08) = +rrt(014)+rrt(016)-rrt(018)-rrt(026)-rrt(027)+rrt(033)+rrt(040)+  2.d0 * rrt(045)+rrt(047)+rrt(052)+rrt(062)+rrt(068)&
             +rrt(080)+rrt(086)+rrt(091)+rrt(108)+rrt(119)+rrt(129)+rrt(136)+rrt(141)+rrt(142)+rrt(161)+rrt(162)+  2.d0 * rrt(173)&
             +  2.d0 * rrt(174)+rrt(184)+rrt(185)+rrt(210)+rrt(211) 
  ydot(09) = +rrt(004)-rrt(154)-rrt(156)-rrt(158)-rrt(160)-rrt(162)-rrt(164)-rrt(191)-rrt(193)-rrt(195)-rrt(197)-rrt(199)-rrt(201)&
             -rrt(203) 
  ydot(10) = +rrt(013)+rrt(015)+rrt(016)+rrt(018)+rrt(020)+rrt(022)+rrt(024)+rrt(027)+rrt(029)+rrt(031)+rrt(035)+  2.d0 * rrt(037)&
             +rrt(038)+rrt(041)+  2.d0 * rrt(043)+rrt(046)+rrt(049)+rrt(055)+rrt(057)+rrt(060)+rrt(065)+rrt(070)+rrt(072)+rrt(077)&
             +rrt(079)+rrt(084)+rrt(089)+rrt(096)+rrt(098)+rrt(105)+rrt(107)+rrt(116)+rrt(127)+rrt(139)+rrt(140)+rrt(143)+rrt(144)&
             +rrt(147)+rrt(148)+rrt(151)+rrt(152)+rrt(153)+rrt(154)+rrt(157)+rrt(158)+rrt(165)+rrt(166)+  2.d0 * rrt(169)&
             +  2.d0 * rrt(170)+rrt(178)+rrt(179)+rrt(192)+rrt(193)+rrt(196)+rrt(197)+rrt(206)+rrt(207)+rrt(222)+rrt(223)+rrt(226)&
             +rrt(227)+  2.d0 * rrt(238) 
  ydot(11) = +rrt(099)+rrt(107)+rrt(117)+rrt(127)+rrt(134)+rrt(228)+rrt(229) 
  ydot(12) = +rrt(069)+rrt(190)+rrt(191) 
  ydot(13) = +rrt(012) 
  ydot(14) = +rrt(019)+rrt(111)+rrt(122)+rrt(145)+rrt(146) 
  ydot(15) = +rrt(005)-rrt(165)-rrt(167)-rrt(169)-rrt(171)-rrt(173)-rrt(204)-rrt(206)-rrt(208)-rrt(210)-rrt(212) 
  ydot(16) = +rrt(013)-rrt(016)-rrt(017)-rrt(023)-rrt(024)-rrt(025)+  2.d0 * rrt(034)+rrt(040)+rrt(044)+rrt(053)+rrt(058)+rrt(063)&
             +rrt(066)+rrt(074)+rrt(081)+rrt(087)+rrt(100)+rrt(109)+rrt(120)+rrt(130)+rrt(139)+rrt(140)+  2.d0 * rrt(163)&
             +  2.d0 * rrt(164)+rrt(171)+rrt(172)+rrt(186)+rrt(187)+rrt(200)+rrt(201)+rrt(212)+rrt(213)+rrt(230)+rrt(231) 
  ydot(17) = +rrt(020)+rrt(023)+rrt(074)+rrt(080)+rrt(085)+rrt(102)+rrt(112)+rrt(123)+rrt(131)+rrt(147)+rrt(148)+rrt(200)+rrt(201)&
             +rrt(208)+rrt(209)+rrt(234)+rrt(235) 
  ydot(18) = -rrt(003)-rrt(004)-rrt(029)-rrt(030)-rrt(031)-rrt(032)-rrt(033)-rrt(034)+rrt(052)-rrt(069)-rrt(070)-rrt(071)-rrt(072)&
             -rrt(073)-rrt(074)-rrt(075)+rrt(103)+rrt(114)+rrt(184)+rrt(185)+rrt(236)+rrt(237) 
  ydot(19) = +rrt(097)+rrt(105)+rrt(115)+rrt(224)+rrt(225) 
  ydot(20) = -rrt(010)-rrt(011)-rrt(049)-rrt(050)-rrt(051)-rrt(052)-rrt(053)-rrt(054)-rrt(095)-rrt(096)-rrt(097)-rrt(098)-rrt(099)&
             -rrt(100)-rrt(101)-rrt(102)-rrt(103) 
  ydot(21) = +rrt(056)+rrt(060)-rrt(065)-rrt(066)-rrt(126)-rrt(127)-rrt(128)-rrt(129)-rrt(130)-rrt(131)-rrt(132)-rrt(133) 
  ydot(22) = +rrt(098)+rrt(106)+rrt(116)+rrt(126)+rrt(226)+rrt(227) 
  ydot(23) = +rrt(096)+rrt(104)+rrt(222)+rrt(223) 
  ydot(24) = +rrt(011)-rrt(179)-rrt(181)-rrt(183)-rrt(185)-rrt(187)-rrt(189)-rrt(221)-rrt(223)-rrt(225)-rrt(227)-rrt(229)-rrt(231)&
             -rrt(233)-rrt(235)-rrt(237) 
  ydot(25) = +rrt(095)+rrt(220)+rrt(221) 
  ydot(26) = +rrt(007)-rrt(175)-rrt(214)-rrt(217) 
  ydot(27) = +rrt(049)-rrt(055)-rrt(056)-rrt(057)-rrt(058)-rrt(059)-rrt(104)-rrt(105)-rrt(106)-rrt(107)-rrt(108)-rrt(109)-rrt(110)&
             -rrt(111)-rrt(112)-rrt(113)-rrt(114)+rrt(178)+rrt(179) 
  ydot(28) = +rrt(009)-rrt(177)-rrt(216)-rrt(219) 
  ydot(29) = +rrt(015)+rrt(017)+rrt(018)-rrt(028)+rrt(039)+rrt(044)+rrt(047)+  2.d0 * rrt(048)+rrt(067)+rrt(085)+rrt(090)+rrt(094)&
             +rrt(118)+rrt(128)+rrt(135)+rrt(143)+rrt(144)+rrt(171)+rrt(172)+  2.d0 * rrt(175)+  2.d0 * rrt(176)+  2.d0 * rrt(177)&
             +rrt(208)+rrt(209)+rrt(217)+rrt(218)+rrt(219) 
  ydot(30) = +rrt(029)-rrt(035)-rrt(036)-rrt(037)-rrt(038)-rrt(039)-rrt(040)+rrt(053)-rrt(076)-rrt(077)-rrt(078)-rrt(079)-rrt(080)&
             -rrt(081)-rrt(082)+rrt(102)+rrt(113)+rrt(125)+rrt(153)+rrt(154)+rrt(186)+rrt(187)+rrt(234)+rrt(235) 
  ydot(31) = -rrt(012)+rrt(050)+rrt(055)-rrt(060)-rrt(061)-rrt(062)-rrt(063)-rrt(064)-rrt(115)-rrt(116)-rrt(117)-rrt(118)-rrt(119)&
             -rrt(120)-rrt(121)-rrt(122)-rrt(123)-rrt(124)-rrt(125)+rrt(180)+rrt(181) 
  ydot(32) = -rrt(007)-rrt(008)-rrt(009)+rrt(032)+rrt(038)+rrt(042)+rrt(043)+rrt(046)-rrt(048)+rrt(064)+rrt(066)+rrt(068)+rrt(092)&
             -rrt(093)-rrt(094)+rrt(122)+rrt(131)+rrt(137)+rrt(159)+rrt(160)+rrt(167)+rrt(168)+rrt(169)+rrt(170) 
  ydot(33) = +rrt(071)+rrt(077)+rrt(083)+rrt(101)+rrt(109)+rrt(119)+rrt(128)+rrt(194)+rrt(195)+rrt(204)+rrt(205)+rrt(232)+rrt(233) 
  ydot(34) = +rrt(014)+rrt(015)+rrt(017)+rrt(021)+rrt(022)+rrt(025)+rrt(030)+rrt(031)+  2.d0 * rrt(032)+rrt(036)+rrt(038)+rrt(042)&
             +rrt(050)+  2.d0 * rrt(051)+rrt(056)+rrt(057)+rrt(061)+rrt(071)+rrt(072)+  2.d0 * rrt(073)+rrt(078)+rrt(079)+rrt(097)&
             +rrt(098)+  2.d0 * rrt(099)+rrt(106)+rrt(107)+rrt(117)+rrt(141)+rrt(142)+rrt(143)+rrt(144)+rrt(149)+rrt(150)+rrt(151)&
             +rrt(152)+rrt(155)+rrt(156)+rrt(157)+rrt(158)+  2.d0 * rrt(159)+  2.d0 * rrt(160)+rrt(167)+rrt(168)+rrt(180)+rrt(181)&
             +  2.d0 * rrt(182)+  2.d0 * rrt(183)+rrt(194)+rrt(195)+rrt(196)+rrt(197)+  2.d0 * rrt(198)+  2.d0 * rrt(199)+rrt(224)&
             +rrt(225)+rrt(226)+rrt(227)+  2.d0 * rrt(228)+  2.d0 * rrt(229)-rrt(238)-rrt(239) 
  ydot(35) = +rrt(003)-rrt(153)-rrt(155)-rrt(157)-rrt(159)-rrt(161)-rrt(163)-rrt(190)-rrt(192)-rrt(194)-rrt(196)-rrt(198)-rrt(200)&
             -rrt(202) 
  ydot(36) = +rrt(008)-rrt(176)-rrt(215)-rrt(218) 
  ydot(37) = +rrt(002)-rrt(140)-rrt(142)-rrt(144)-rrt(146)-rrt(148)-rrt(150)-rrt(152) 
  ydot(38) = +rrt(072)+rrt(078)+rrt(084)+rrt(088)+rrt(110)+rrt(120)+rrt(129)+rrt(135)+rrt(196)+rrt(197)+rrt(206)+rrt(207) 
  ydot(39) = +rrt(001)-rrt(139)-rrt(141)-rrt(143)-rrt(145)-rrt(147)-rrt(149)-rrt(151) 
  ydot(40) = +rrt(073)+rrt(079)+rrt(089)+rrt(093)+rrt(121)+rrt(130)+rrt(136)+rrt(198)+rrt(199)+rrt(214)+rrt(215)+rrt(216) 
  ydot(41) = +rrt(070)+rrt(076)+rrt(100)+rrt(108)+rrt(118)+rrt(192)+rrt(193)+rrt(230)+rrt(231) 
  ydot(42) = +rrt(022)+rrt(025)+rrt(027)+rrt(028)+rrt(082)+rrt(087)+rrt(091)+rrt(094)+rrt(114)+rrt(125)+rrt(133)+rrt(138)+rrt(151)&
             +rrt(152)+rrt(212)+rrt(213)+rrt(217)+rrt(218)+rrt(219) 
  ydot(43) = +rrt(239) 
  ydot(44) = -rrt(001)-rrt(002)-rrt(013)-rrt(014)-rrt(015)-rrt(019)-rrt(020)-rrt(021)-rrt(022)+rrt(033)+rrt(039)+rrt(054)+rrt(059)&
             +rrt(064)+rrt(075)+rrt(082)+rrt(101)+rrt(110)+rrt(121)+rrt(161)+rrt(162)+rrt(188)+rrt(189)+rrt(202)+rrt(203)+rrt(232)&
             +rrt(233) 
  ydot(45) = -rrt(005)-rrt(006)+rrt(030)+rrt(035)-rrt(041)-rrt(042)-rrt(043)-rrt(044)-rrt(045)+rrt(054)+rrt(058)+rrt(062)-rrt(083)&
             -rrt(084)-rrt(085)-rrt(086)-rrt(087)+rrt(112)+rrt(124)+rrt(133)+rrt(155)+rrt(156)+rrt(188)+rrt(189) 
  if( ldensity_constant ) where( density_constant(:) ) ydot(1:species_max) = 0.0d0
  ydot(46) = 0.0d0
  if( lgas_heating ) then
    ydot(46) = ( ZDPlasKin_cfg(14)/k_B + ydot(46) ) / ( sum(density(1:species_max)) - density(species_electrons) ) &
            + eV_to_K * ZDPlasKin_cfg(11) * density(species_electrons)
    ydot(46) = ydot(46) * ZDPlasKin_cfg(13)
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
  if( lgas_heating ) ZDPlasKin_cfg(1) = y(46)
  density(:) = y(1:species_max)
  call ZDPlasKin_reac_rates(t)
  pd(39,01) = pd(39,01) + rrt(001) * density(44) 
  pd(39,44) = pd(39,44) + rrt(001) * density(01) 
  pd(44,01) = pd(44,01) - rrt(001) * density(44) 
  pd(44,44) = pd(44,44) - rrt(001) * density(01) 
  pd(37,01) = pd(37,01) + rrt(002) * density(44) 
  pd(37,44) = pd(37,44) + rrt(002) * density(01) 
  pd(44,01) = pd(44,01) - rrt(002) * density(44) 
  pd(44,44) = pd(44,44) - rrt(002) * density(01) 
  pd(18,01) = pd(18,01) - rrt(003) * density(18) 
  pd(18,18) = pd(18,18) - rrt(003) * density(01) 
  pd(35,01) = pd(35,01) + rrt(003) * density(18) 
  pd(35,18) = pd(35,18) + rrt(003) * density(01) 
  pd(09,01) = pd(09,01) + rrt(004) * density(18) 
  pd(09,18) = pd(09,18) + rrt(004) * density(01) 
  pd(18,01) = pd(18,01) - rrt(004) * density(18) 
  pd(18,18) = pd(18,18) - rrt(004) * density(01) 
  pd(15,01) = pd(15,01) + rrt(005) * density(45) 
  pd(15,45) = pd(15,45) + rrt(005) * density(01) 
  pd(45,01) = pd(45,01) - rrt(005) * density(45) 
  pd(45,45) = pd(45,45) - rrt(005) * density(01) 
  pd(02,01) = pd(02,01) + rrt(006) * density(45) 
  pd(02,45) = pd(02,45) + rrt(006) * density(01) 
  pd(45,01) = pd(45,01) - rrt(006) * density(45) 
  pd(45,45) = pd(45,45) - rrt(006) * density(01) 
  pd(26,01) = pd(26,01) + rrt(007) * density(32) 
  pd(26,32) = pd(26,32) + rrt(007) * density(01) 
  pd(32,01) = pd(32,01) - rrt(007) * density(32) 
  pd(32,32) = pd(32,32) - rrt(007) * density(01) 
  pd(32,01) = pd(32,01) - rrt(008) * density(32) 
  pd(32,32) = pd(32,32) - rrt(008) * density(01) 
  pd(36,01) = pd(36,01) + rrt(008) * density(32) 
  pd(36,32) = pd(36,32) + rrt(008) * density(01) 
  pd(28,01) = pd(28,01) + rrt(009) * density(32) 
  pd(28,32) = pd(28,32) + rrt(009) * density(01) 
  pd(32,01) = pd(32,01) - rrt(009) * density(32) 
  pd(32,32) = pd(32,32) - rrt(009) * density(01) 
  pd(04,01) = pd(04,01) + rrt(010) * density(20) 
  pd(04,20) = pd(04,20) + rrt(010) * density(01) 
  pd(20,01) = pd(20,01) - rrt(010) * density(20) 
  pd(20,20) = pd(20,20) - rrt(010) * density(01) 
  pd(20,01) = pd(20,01) - rrt(011) * density(20) 
  pd(20,20) = pd(20,20) - rrt(011) * density(01) 
  pd(24,01) = pd(24,01) + rrt(011) * density(20) 
  pd(24,20) = pd(24,20) + rrt(011) * density(01) 
  pd(13,01) = pd(13,01) + rrt(012) * density(31) 
  pd(13,31) = pd(13,31) + rrt(012) * density(01) 
  pd(31,01) = pd(31,01) - rrt(012) * density(31) 
  pd(31,31) = pd(31,31) - rrt(012) * density(01) 
  pd(10,01) = pd(10,01) + rrt(013) * density(44) 
  pd(10,44) = pd(10,44) + rrt(013) * density(01) 
  pd(16,01) = pd(16,01) + rrt(013) * density(44) 
  pd(16,44) = pd(16,44) + rrt(013) * density(01) 
  pd(44,01) = pd(44,01) - rrt(013) * density(44) 
  pd(44,44) = pd(44,44) - rrt(013) * density(01) 
  pd(08,01) = pd(08,01) + rrt(014) * density(44) 
  pd(08,44) = pd(08,44) + rrt(014) * density(01) 
  pd(34,01) = pd(34,01) + rrt(014) * density(44) 
  pd(34,44) = pd(34,44) + rrt(014) * density(01) 
  pd(44,01) = pd(44,01) - rrt(014) * density(44) 
  pd(44,44) = pd(44,44) - rrt(014) * density(01) 
  pd(10,01) = pd(10,01) + rrt(015) * density(44) 
  pd(10,44) = pd(10,44) + rrt(015) * density(01) 
  pd(29,01) = pd(29,01) + rrt(015) * density(44) 
  pd(29,44) = pd(29,44) + rrt(015) * density(01) 
  pd(34,01) = pd(34,01) + rrt(015) * density(44) 
  pd(34,44) = pd(34,44) + rrt(015) * density(01) 
  pd(44,01) = pd(44,01) - rrt(015) * density(44) 
  pd(44,44) = pd(44,44) - rrt(015) * density(01) 
  pd(08,01) = pd(08,01) + rrt(016) * density(16) 
  pd(08,16) = pd(08,16) + rrt(016) * density(01) 
  pd(10,01) = pd(10,01) + rrt(016) * density(16) 
  pd(10,16) = pd(10,16) + rrt(016) * density(01) 
  pd(16,01) = pd(16,01) - rrt(016) * density(16) 
  pd(16,16) = pd(16,16) - rrt(016) * density(01) 
  pd(16,01) = pd(16,01) - rrt(017) * density(16) 
  pd(16,16) = pd(16,16) - rrt(017) * density(01) 
  pd(29,01) = pd(29,01) + rrt(017) * density(16) 
  pd(29,16) = pd(29,16) + rrt(017) * density(01) 
  pd(34,01) = pd(34,01) + rrt(017) * density(16) 
  pd(34,16) = pd(34,16) + rrt(017) * density(01) 
  pd(08,01) = pd(08,01) - rrt(018) * density(08) 
  pd(08,08) = pd(08,08) - rrt(018) * density(01) 
  pd(10,01) = pd(10,01) + rrt(018) * density(08) 
  pd(10,08) = pd(10,08) + rrt(018) * density(01) 
  pd(29,01) = pd(29,01) + rrt(018) * density(08) 
  pd(29,08) = pd(29,08) + rrt(018) * density(01) 
  pd(01,01) = pd(01,01) + rrt(019) * density(44) 
  pd(01,44) = pd(01,44) + rrt(019) * density(01) 
  pd(14,01) = pd(14,01) + rrt(019) * density(44) 
  pd(14,44) = pd(14,44) + rrt(019) * density(01) 
  pd(44,01) = pd(44,01) - rrt(019) * density(44) 
  pd(44,44) = pd(44,44) - rrt(019) * density(01) 
  pd(01,01) = pd(01,01) + rrt(020) * density(44) 
  pd(01,44) = pd(01,44) + rrt(020) * density(01) 
  pd(10,01) = pd(10,01) + rrt(020) * density(44) 
  pd(10,44) = pd(10,44) + rrt(020) * density(01) 
  pd(17,01) = pd(17,01) + rrt(020) * density(44) 
  pd(17,44) = pd(17,44) + rrt(020) * density(01) 
  pd(44,01) = pd(44,01) - rrt(020) * density(44) 
  pd(44,44) = pd(44,44) - rrt(020) * density(01) 
  pd(01,01) = pd(01,01) + rrt(021) * density(44) 
  pd(01,44) = pd(01,44) + rrt(021) * density(01) 
  pd(03,01) = pd(03,01) + rrt(021) * density(44) 
  pd(03,44) = pd(03,44) + rrt(021) * density(01) 
  pd(34,01) = pd(34,01) + rrt(021) * density(44) 
  pd(34,44) = pd(34,44) + rrt(021) * density(01) 
  pd(44,01) = pd(44,01) - rrt(021) * density(44) 
  pd(44,44) = pd(44,44) - rrt(021) * density(01) 
  pd(01,01) = pd(01,01) + rrt(022) * density(44) 
  pd(01,44) = pd(01,44) + rrt(022) * density(01) 
  pd(10,01) = pd(10,01) + rrt(022) * density(44) 
  pd(10,44) = pd(10,44) + rrt(022) * density(01) 
  pd(34,01) = pd(34,01) + rrt(022) * density(44) 
  pd(34,44) = pd(34,44) + rrt(022) * density(01) 
  pd(42,01) = pd(42,01) + rrt(022) * density(44) 
  pd(42,44) = pd(42,44) + rrt(022) * density(01) 
  pd(44,01) = pd(44,01) - rrt(022) * density(44) 
  pd(44,44) = pd(44,44) - rrt(022) * density(01) 
  pd(01,01) = pd(01,01) + rrt(023) * density(16) 
  pd(01,16) = pd(01,16) + rrt(023) * density(01) 
  pd(16,01) = pd(16,01) - rrt(023) * density(16) 
  pd(16,16) = pd(16,16) - rrt(023) * density(01) 
  pd(17,01) = pd(17,01) + rrt(023) * density(16) 
  pd(17,16) = pd(17,16) + rrt(023) * density(01) 
  pd(01,01) = pd(01,01) + rrt(024) * density(16) 
  pd(01,16) = pd(01,16) + rrt(024) * density(01) 
  pd(03,01) = pd(03,01) + rrt(024) * density(16) 
  pd(03,16) = pd(03,16) + rrt(024) * density(01) 
  pd(10,01) = pd(10,01) + rrt(024) * density(16) 
  pd(10,16) = pd(10,16) + rrt(024) * density(01) 
  pd(16,01) = pd(16,01) - rrt(024) * density(16) 
  pd(16,16) = pd(16,16) - rrt(024) * density(01) 
  pd(01,01) = pd(01,01) + rrt(025) * density(16) 
  pd(01,16) = pd(01,16) + rrt(025) * density(01) 
  pd(16,01) = pd(16,01) - rrt(025) * density(16) 
  pd(16,16) = pd(16,16) - rrt(025) * density(01) 
  pd(34,01) = pd(34,01) + rrt(025) * density(16) 
  pd(34,16) = pd(34,16) + rrt(025) * density(01) 
  pd(42,01) = pd(42,01) + rrt(025) * density(16) 
  pd(42,16) = pd(42,16) + rrt(025) * density(01) 
  pd(01,01) = pd(01,01) + rrt(026) * density(08) 
  pd(01,08) = pd(01,08) + rrt(026) * density(01) 
  pd(03,01) = pd(03,01) + rrt(026) * density(08) 
  pd(03,08) = pd(03,08) + rrt(026) * density(01) 
  pd(08,01) = pd(08,01) - rrt(026) * density(08) 
  pd(08,08) = pd(08,08) - rrt(026) * density(01) 
  pd(01,01) = pd(01,01) + rrt(027) * density(08) 
  pd(01,08) = pd(01,08) + rrt(027) * density(01) 
  pd(08,01) = pd(08,01) - rrt(027) * density(08) 
  pd(08,08) = pd(08,08) - rrt(027) * density(01) 
  pd(10,01) = pd(10,01) + rrt(027) * density(08) 
  pd(10,08) = pd(10,08) + rrt(027) * density(01) 
  pd(42,01) = pd(42,01) + rrt(027) * density(08) 
  pd(42,08) = pd(42,08) + rrt(027) * density(01) 
  pd(01,01) = pd(01,01) + rrt(028) * density(29) 
  pd(01,29) = pd(01,29) + rrt(028) * density(01) 
  pd(29,01) = pd(29,01) - rrt(028) * density(29) 
  pd(29,29) = pd(29,29) - rrt(028) * density(01) 
  pd(42,01) = pd(42,01) + rrt(028) * density(29) 
  pd(42,29) = pd(42,29) + rrt(028) * density(01) 
  pd(10,01) = pd(10,01) + rrt(029) * density(18) 
  pd(10,18) = pd(10,18) + rrt(029) * density(01) 
  pd(18,01) = pd(18,01) - rrt(029) * density(18) 
  pd(18,18) = pd(18,18) - rrt(029) * density(01) 
  pd(30,01) = pd(30,01) + rrt(029) * density(18) 
  pd(30,18) = pd(30,18) + rrt(029) * density(01) 
  pd(18,01) = pd(18,01) - rrt(030) * density(18) 
  pd(18,18) = pd(18,18) - rrt(030) * density(01) 
  pd(34,01) = pd(34,01) + rrt(030) * density(18) 
  pd(34,18) = pd(34,18) + rrt(030) * density(01) 
  pd(45,01) = pd(45,01) + rrt(030) * density(18) 
  pd(45,18) = pd(45,18) + rrt(030) * density(01) 
  pd(05,01) = pd(05,01) + rrt(031) * density(18) 
  pd(05,18) = pd(05,18) + rrt(031) * density(01) 
  pd(10,01) = pd(10,01) + rrt(031) * density(18) 
  pd(10,18) = pd(10,18) + rrt(031) * density(01) 
  pd(18,01) = pd(18,01) - rrt(031) * density(18) 
  pd(18,18) = pd(18,18) - rrt(031) * density(01) 
  pd(34,01) = pd(34,01) + rrt(031) * density(18) 
  pd(34,18) = pd(34,18) + rrt(031) * density(01) 
  pd(18,01) = pd(18,01) - rrt(032) * density(18) 
  pd(18,18) = pd(18,18) - rrt(032) * density(01) 
  pd(32,01) = pd(32,01) + rrt(032) * density(18) 
  pd(32,18) = pd(32,18) + rrt(032) * density(01) 
  pd(34,01) = pd(34,01) + rrt(032) * density(18) * 2.0d0
  pd(34,18) = pd(34,18) + rrt(032) * density(01) * 2.0d0
  pd(08,01) = pd(08,01) + rrt(033) * density(18) 
  pd(08,18) = pd(08,18) + rrt(033) * density(01) 
  pd(18,01) = pd(18,01) - rrt(033) * density(18) 
  pd(18,18) = pd(18,18) - rrt(033) * density(01) 
  pd(44,01) = pd(44,01) + rrt(033) * density(18) 
  pd(44,18) = pd(44,18) + rrt(033) * density(01) 
  pd(16,01) = pd(16,01) + rrt(034) * density(18) * 2.0d0
  pd(16,18) = pd(16,18) + rrt(034) * density(01) * 2.0d0
  pd(18,01) = pd(18,01) - rrt(034) * density(18) 
  pd(18,18) = pd(18,18) - rrt(034) * density(01) 
  pd(10,01) = pd(10,01) + rrt(035) * density(30) 
  pd(10,30) = pd(10,30) + rrt(035) * density(01) 
  pd(30,01) = pd(30,01) - rrt(035) * density(30) 
  pd(30,30) = pd(30,30) - rrt(035) * density(01) 
  pd(45,01) = pd(45,01) + rrt(035) * density(30) 
  pd(45,30) = pd(45,30) + rrt(035) * density(01) 
  pd(05,01) = pd(05,01) + rrt(036) * density(30) 
  pd(05,30) = pd(05,30) + rrt(036) * density(01) 
  pd(30,01) = pd(30,01) - rrt(036) * density(30) 
  pd(30,30) = pd(30,30) - rrt(036) * density(01) 
  pd(34,01) = pd(34,01) + rrt(036) * density(30) 
  pd(34,30) = pd(34,30) + rrt(036) * density(01) 
  pd(05,01) = pd(05,01) + rrt(037) * density(30) 
  pd(05,30) = pd(05,30) + rrt(037) * density(01) 
  pd(10,01) = pd(10,01) + rrt(037) * density(30) * 2.0d0
  pd(10,30) = pd(10,30) + rrt(037) * density(01) * 2.0d0
  pd(30,01) = pd(30,01) - rrt(037) * density(30) 
  pd(30,30) = pd(30,30) - rrt(037) * density(01) 
  pd(10,01) = pd(10,01) + rrt(038) * density(30) 
  pd(10,30) = pd(10,30) + rrt(038) * density(01) 
  pd(30,01) = pd(30,01) - rrt(038) * density(30) 
  pd(30,30) = pd(30,30) - rrt(038) * density(01) 
  pd(32,01) = pd(32,01) + rrt(038) * density(30) 
  pd(32,30) = pd(32,30) + rrt(038) * density(01) 
  pd(34,01) = pd(34,01) + rrt(038) * density(30) 
  pd(34,30) = pd(34,30) + rrt(038) * density(01) 
  pd(29,01) = pd(29,01) + rrt(039) * density(30) 
  pd(29,30) = pd(29,30) + rrt(039) * density(01) 
  pd(30,01) = pd(30,01) - rrt(039) * density(30) 
  pd(30,30) = pd(30,30) - rrt(039) * density(01) 
  pd(44,01) = pd(44,01) + rrt(039) * density(30) 
  pd(44,30) = pd(44,30) + rrt(039) * density(01) 
  pd(08,01) = pd(08,01) + rrt(040) * density(30) 
  pd(08,30) = pd(08,30) + rrt(040) * density(01) 
  pd(16,01) = pd(16,01) + rrt(040) * density(30) 
  pd(16,30) = pd(16,30) + rrt(040) * density(01) 
  pd(30,01) = pd(30,01) - rrt(040) * density(30) 
  pd(30,30) = pd(30,30) - rrt(040) * density(01) 
  pd(05,01) = pd(05,01) + rrt(041) * density(45) 
  pd(05,45) = pd(05,45) + rrt(041) * density(01) 
  pd(10,01) = pd(10,01) + rrt(041) * density(45) 
  pd(10,45) = pd(10,45) + rrt(041) * density(01) 
  pd(45,01) = pd(45,01) - rrt(041) * density(45) 
  pd(45,45) = pd(45,45) - rrt(041) * density(01) 
  pd(32,01) = pd(32,01) + rrt(042) * density(45) 
  pd(32,45) = pd(32,45) + rrt(042) * density(01) 
  pd(34,01) = pd(34,01) + rrt(042) * density(45) 
  pd(34,45) = pd(34,45) + rrt(042) * density(01) 
  pd(45,01) = pd(45,01) - rrt(042) * density(45) 
  pd(45,45) = pd(45,45) - rrt(042) * density(01) 
  pd(10,01) = pd(10,01) + rrt(043) * density(45) * 2.0d0
  pd(10,45) = pd(10,45) + rrt(043) * density(01) * 2.0d0
  pd(32,01) = pd(32,01) + rrt(043) * density(45) 
  pd(32,45) = pd(32,45) + rrt(043) * density(01) 
  pd(45,01) = pd(45,01) - rrt(043) * density(45) 
  pd(45,45) = pd(45,45) - rrt(043) * density(01) 
  pd(16,01) = pd(16,01) + rrt(044) * density(45) 
  pd(16,45) = pd(16,45) + rrt(044) * density(01) 
  pd(29,01) = pd(29,01) + rrt(044) * density(45) 
  pd(29,45) = pd(29,45) + rrt(044) * density(01) 
  pd(45,01) = pd(45,01) - rrt(044) * density(45) 
  pd(45,45) = pd(45,45) - rrt(044) * density(01) 
  pd(08,01) = pd(08,01) + rrt(045) * density(45) * 2.0d0
  pd(08,45) = pd(08,45) + rrt(045) * density(01) * 2.0d0
  pd(45,01) = pd(45,01) - rrt(045) * density(45) 
  pd(45,45) = pd(45,45) - rrt(045) * density(01) 
  pd(05,01) = pd(05,01) - rrt(046) * density(05) 
  pd(05,05) = pd(05,05) - rrt(046) * density(01) 
  pd(10,01) = pd(10,01) + rrt(046) * density(05) 
  pd(10,05) = pd(10,05) + rrt(046) * density(01) 
  pd(32,01) = pd(32,01) + rrt(046) * density(05) 
  pd(32,05) = pd(32,05) + rrt(046) * density(01) 
  pd(05,01) = pd(05,01) - rrt(047) * density(05) 
  pd(05,05) = pd(05,05) - rrt(047) * density(01) 
  pd(08,01) = pd(08,01) + rrt(047) * density(05) 
  pd(08,05) = pd(08,05) + rrt(047) * density(01) 
  pd(29,01) = pd(29,01) + rrt(047) * density(05) 
  pd(29,05) = pd(29,05) + rrt(047) * density(01) 
  pd(29,01) = pd(29,01) + rrt(048) * density(32) * 2.0d0
  pd(29,32) = pd(29,32) + rrt(048) * density(01) * 2.0d0
  pd(32,01) = pd(32,01) - rrt(048) * density(32) 
  pd(32,32) = pd(32,32) - rrt(048) * density(01) 
  pd(10,01) = pd(10,01) + rrt(049) * density(20) 
  pd(10,20) = pd(10,20) + rrt(049) * density(01) 
  pd(20,01) = pd(20,01) - rrt(049) * density(20) 
  pd(20,20) = pd(20,20) - rrt(049) * density(01) 
  pd(27,01) = pd(27,01) + rrt(049) * density(20) 
  pd(27,20) = pd(27,20) + rrt(049) * density(01) 
  pd(20,01) = pd(20,01) - rrt(050) * density(20) 
  pd(20,20) = pd(20,20) - rrt(050) * density(01) 
  pd(31,01) = pd(31,01) + rrt(050) * density(20) 
  pd(31,20) = pd(31,20) + rrt(050) * density(01) 
  pd(34,01) = pd(34,01) + rrt(050) * density(20) 
  pd(34,20) = pd(34,20) + rrt(050) * density(01) 
  pd(06,01) = pd(06,01) + rrt(051) * density(20) 
  pd(06,20) = pd(06,20) + rrt(051) * density(01) 
  pd(20,01) = pd(20,01) - rrt(051) * density(20) 
  pd(20,20) = pd(20,20) - rrt(051) * density(01) 
  pd(34,01) = pd(34,01) + rrt(051) * density(20) * 2.0d0
  pd(34,20) = pd(34,20) + rrt(051) * density(01) * 2.0d0
  pd(08,01) = pd(08,01) + rrt(052) * density(20) 
  pd(08,20) = pd(08,20) + rrt(052) * density(01) 
  pd(18,01) = pd(18,01) + rrt(052) * density(20) 
  pd(18,20) = pd(18,20) + rrt(052) * density(01) 
  pd(20,01) = pd(20,01) - rrt(052) * density(20) 
  pd(20,20) = pd(20,20) - rrt(052) * density(01) 
  pd(16,01) = pd(16,01) + rrt(053) * density(20) 
  pd(16,20) = pd(16,20) + rrt(053) * density(01) 
  pd(20,01) = pd(20,01) - rrt(053) * density(20) 
  pd(20,20) = pd(20,20) - rrt(053) * density(01) 
  pd(30,01) = pd(30,01) + rrt(053) * density(20) 
  pd(30,20) = pd(30,20) + rrt(053) * density(01) 
  pd(20,01) = pd(20,01) - rrt(054) * density(20) 
  pd(20,20) = pd(20,20) - rrt(054) * density(01) 
  pd(44,01) = pd(44,01) + rrt(054) * density(20) 
  pd(44,20) = pd(44,20) + rrt(054) * density(01) 
  pd(45,01) = pd(45,01) + rrt(054) * density(20) 
  pd(45,20) = pd(45,20) + rrt(054) * density(01) 
  pd(10,01) = pd(10,01) + rrt(055) * density(27) 
  pd(10,27) = pd(10,27) + rrt(055) * density(01) 
  pd(27,01) = pd(27,01) - rrt(055) * density(27) 
  pd(27,27) = pd(27,27) - rrt(055) * density(01) 
  pd(31,01) = pd(31,01) + rrt(055) * density(27) 
  pd(31,27) = pd(31,27) + rrt(055) * density(01) 
  pd(21,01) = pd(21,01) + rrt(056) * density(27) 
  pd(21,27) = pd(21,27) + rrt(056) * density(01) 
  pd(27,01) = pd(27,01) - rrt(056) * density(27) 
  pd(27,27) = pd(27,27) - rrt(056) * density(01) 
  pd(34,01) = pd(34,01) + rrt(056) * density(27) 
  pd(34,27) = pd(34,27) + rrt(056) * density(01) 
  pd(06,01) = pd(06,01) + rrt(057) * density(27) 
  pd(06,27) = pd(06,27) + rrt(057) * density(01) 
  pd(10,01) = pd(10,01) + rrt(057) * density(27) 
  pd(10,27) = pd(10,27) + rrt(057) * density(01) 
  pd(27,01) = pd(27,01) - rrt(057) * density(27) 
  pd(27,27) = pd(27,27) - rrt(057) * density(01) 
  pd(34,01) = pd(34,01) + rrt(057) * density(27) 
  pd(34,27) = pd(34,27) + rrt(057) * density(01) 
  pd(16,01) = pd(16,01) + rrt(058) * density(27) 
  pd(16,27) = pd(16,27) + rrt(058) * density(01) 
  pd(27,01) = pd(27,01) - rrt(058) * density(27) 
  pd(27,27) = pd(27,27) - rrt(058) * density(01) 
  pd(45,01) = pd(45,01) + rrt(058) * density(27) 
  pd(45,27) = pd(45,27) + rrt(058) * density(01) 
  pd(05,01) = pd(05,01) + rrt(059) * density(27) 
  pd(05,27) = pd(05,27) + rrt(059) * density(01) 
  pd(27,01) = pd(27,01) - rrt(059) * density(27) 
  pd(27,27) = pd(27,27) - rrt(059) * density(01) 
  pd(44,01) = pd(44,01) + rrt(059) * density(27) 
  pd(44,27) = pd(44,27) + rrt(059) * density(01) 
  pd(10,01) = pd(10,01) + rrt(060) * density(31) 
  pd(10,31) = pd(10,31) + rrt(060) * density(01) 
  pd(21,01) = pd(21,01) + rrt(060) * density(31) 
  pd(21,31) = pd(21,31) + rrt(060) * density(01) 
  pd(31,01) = pd(31,01) - rrt(060) * density(31) 
  pd(31,31) = pd(31,31) - rrt(060) * density(01) 
  pd(06,01) = pd(06,01) + rrt(061) * density(31) 
  pd(06,31) = pd(06,31) + rrt(061) * density(01) 
  pd(31,01) = pd(31,01) - rrt(061) * density(31) 
  pd(31,31) = pd(31,31) - rrt(061) * density(01) 
  pd(34,01) = pd(34,01) + rrt(061) * density(31) 
  pd(34,31) = pd(34,31) + rrt(061) * density(01) 
  pd(08,01) = pd(08,01) + rrt(062) * density(31) 
  pd(08,31) = pd(08,31) + rrt(062) * density(01) 
  pd(31,01) = pd(31,01) - rrt(062) * density(31) 
  pd(31,31) = pd(31,31) - rrt(062) * density(01) 
  pd(45,01) = pd(45,01) + rrt(062) * density(31) 
  pd(45,31) = pd(45,31) + rrt(062) * density(01) 
  pd(05,01) = pd(05,01) + rrt(063) * density(31) 
  pd(05,31) = pd(05,31) + rrt(063) * density(01) 
  pd(16,01) = pd(16,01) + rrt(063) * density(31) 
  pd(16,31) = pd(16,31) + rrt(063) * density(01) 
  pd(31,01) = pd(31,01) - rrt(063) * density(31) 
  pd(31,31) = pd(31,31) - rrt(063) * density(01) 
  pd(31,01) = pd(31,01) - rrt(064) * density(31) 
  pd(31,31) = pd(31,31) - rrt(064) * density(01) 
  pd(32,01) = pd(32,01) + rrt(064) * density(31) 
  pd(32,31) = pd(32,31) + rrt(064) * density(01) 
  pd(44,01) = pd(44,01) + rrt(064) * density(31) 
  pd(44,31) = pd(44,31) + rrt(064) * density(01) 
  pd(06,01) = pd(06,01) + rrt(065) * density(21) 
  pd(06,21) = pd(06,21) + rrt(065) * density(01) 
  pd(10,01) = pd(10,01) + rrt(065) * density(21) 
  pd(10,21) = pd(10,21) + rrt(065) * density(01) 
  pd(21,01) = pd(21,01) - rrt(065) * density(21) 
  pd(21,21) = pd(21,21) - rrt(065) * density(01) 
  pd(16,01) = pd(16,01) + rrt(066) * density(21) 
  pd(16,21) = pd(16,21) + rrt(066) * density(01) 
  pd(21,01) = pd(21,01) - rrt(066) * density(21) 
  pd(21,21) = pd(21,21) - rrt(066) * density(01) 
  pd(32,01) = pd(32,01) + rrt(066) * density(21) 
  pd(32,21) = pd(32,21) + rrt(066) * density(01) 
  pd(05,01) = pd(05,01) + rrt(067) * density(06) 
  pd(05,06) = pd(05,06) + rrt(067) * density(01) 
  pd(06,01) = pd(06,01) - rrt(067) * density(06) 
  pd(06,06) = pd(06,06) - rrt(067) * density(01) 
  pd(29,01) = pd(29,01) + rrt(067) * density(06) 
  pd(29,06) = pd(29,06) + rrt(067) * density(01) 
  pd(06,01) = pd(06,01) - rrt(068) * density(06) 
  pd(06,06) = pd(06,06) - rrt(068) * density(01) 
  pd(08,01) = pd(08,01) + rrt(068) * density(06) 
  pd(08,06) = pd(08,06) + rrt(068) * density(01) 
  pd(32,01) = pd(32,01) + rrt(068) * density(06) 
  pd(32,06) = pd(32,06) + rrt(068) * density(01) 
  pd(01,01) = pd(01,01) + rrt(069) * density(18) 
  pd(01,18) = pd(01,18) + rrt(069) * density(01) 
  pd(12,01) = pd(12,01) + rrt(069) * density(18) 
  pd(12,18) = pd(12,18) + rrt(069) * density(01) 
  pd(18,01) = pd(18,01) - rrt(069) * density(18) 
  pd(18,18) = pd(18,18) - rrt(069) * density(01) 
  pd(01,01) = pd(01,01) + rrt(070) * density(18) 
  pd(01,18) = pd(01,18) + rrt(070) * density(01) 
  pd(10,01) = pd(10,01) + rrt(070) * density(18) 
  pd(10,18) = pd(10,18) + rrt(070) * density(01) 
  pd(18,01) = pd(18,01) - rrt(070) * density(18) 
  pd(18,18) = pd(18,18) - rrt(070) * density(01) 
  pd(41,01) = pd(41,01) + rrt(070) * density(18) 
  pd(41,18) = pd(41,18) + rrt(070) * density(01) 
  pd(01,01) = pd(01,01) + rrt(071) * density(18) 
  pd(01,18) = pd(01,18) + rrt(071) * density(01) 
  pd(18,01) = pd(18,01) - rrt(071) * density(18) 
  pd(18,18) = pd(18,18) - rrt(071) * density(01) 
  pd(33,01) = pd(33,01) + rrt(071) * density(18) 
  pd(33,18) = pd(33,18) + rrt(071) * density(01) 
  pd(34,01) = pd(34,01) + rrt(071) * density(18) 
  pd(34,18) = pd(34,18) + rrt(071) * density(01) 
  pd(01,01) = pd(01,01) + rrt(072) * density(18) 
  pd(01,18) = pd(01,18) + rrt(072) * density(01) 
  pd(10,01) = pd(10,01) + rrt(072) * density(18) 
  pd(10,18) = pd(10,18) + rrt(072) * density(01) 
  pd(18,01) = pd(18,01) - rrt(072) * density(18) 
  pd(18,18) = pd(18,18) - rrt(072) * density(01) 
  pd(34,01) = pd(34,01) + rrt(072) * density(18) 
  pd(34,18) = pd(34,18) + rrt(072) * density(01) 
  pd(38,01) = pd(38,01) + rrt(072) * density(18) 
  pd(38,18) = pd(38,18) + rrt(072) * density(01) 
  pd(01,01) = pd(01,01) + rrt(073) * density(18) 
  pd(01,18) = pd(01,18) + rrt(073) * density(01) 
  pd(18,01) = pd(18,01) - rrt(073) * density(18) 
  pd(18,18) = pd(18,18) - rrt(073) * density(01) 
  pd(34,01) = pd(34,01) + rrt(073) * density(18) * 2.0d0
  pd(34,18) = pd(34,18) + rrt(073) * density(01) * 2.0d0
  pd(40,01) = pd(40,01) + rrt(073) * density(18) 
  pd(40,18) = pd(40,18) + rrt(073) * density(01) 
  pd(01,01) = pd(01,01) + rrt(074) * density(18) 
  pd(01,18) = pd(01,18) + rrt(074) * density(01) 
  pd(16,01) = pd(16,01) + rrt(074) * density(18) 
  pd(16,18) = pd(16,18) + rrt(074) * density(01) 
  pd(17,01) = pd(17,01) + rrt(074) * density(18) 
  pd(17,18) = pd(17,18) + rrt(074) * density(01) 
  pd(18,01) = pd(18,01) - rrt(074) * density(18) 
  pd(18,18) = pd(18,18) - rrt(074) * density(01) 
  pd(01,01) = pd(01,01) + rrt(075) * density(18) 
  pd(01,18) = pd(01,18) + rrt(075) * density(01) 
  pd(03,01) = pd(03,01) + rrt(075) * density(18) 
  pd(03,18) = pd(03,18) + rrt(075) * density(01) 
  pd(18,01) = pd(18,01) - rrt(075) * density(18) 
  pd(18,18) = pd(18,18) - rrt(075) * density(01) 
  pd(44,01) = pd(44,01) + rrt(075) * density(18) 
  pd(44,18) = pd(44,18) + rrt(075) * density(01) 
  pd(01,01) = pd(01,01) + rrt(076) * density(30) 
  pd(01,30) = pd(01,30) + rrt(076) * density(01) 
  pd(30,01) = pd(30,01) - rrt(076) * density(30) 
  pd(30,30) = pd(30,30) - rrt(076) * density(01) 
  pd(41,01) = pd(41,01) + rrt(076) * density(30) 
  pd(41,30) = pd(41,30) + rrt(076) * density(01) 
  pd(01,01) = pd(01,01) + rrt(077) * density(30) 
  pd(01,30) = pd(01,30) + rrt(077) * density(01) 
  pd(10,01) = pd(10,01) + rrt(077) * density(30) 
  pd(10,30) = pd(10,30) + rrt(077) * density(01) 
  pd(30,01) = pd(30,01) - rrt(077) * density(30) 
  pd(30,30) = pd(30,30) - rrt(077) * density(01) 
  pd(33,01) = pd(33,01) + rrt(077) * density(30) 
  pd(33,30) = pd(33,30) + rrt(077) * density(01) 
  pd(01,01) = pd(01,01) + rrt(078) * density(30) 
  pd(01,30) = pd(01,30) + rrt(078) * density(01) 
  pd(30,01) = pd(30,01) - rrt(078) * density(30) 
  pd(30,30) = pd(30,30) - rrt(078) * density(01) 
  pd(34,01) = pd(34,01) + rrt(078) * density(30) 
  pd(34,30) = pd(34,30) + rrt(078) * density(01) 
  pd(38,01) = pd(38,01) + rrt(078) * density(30) 
  pd(38,30) = pd(38,30) + rrt(078) * density(01) 
  pd(01,01) = pd(01,01) + rrt(079) * density(30) 
  pd(01,30) = pd(01,30) + rrt(079) * density(01) 
  pd(10,01) = pd(10,01) + rrt(079) * density(30) 
  pd(10,30) = pd(10,30) + rrt(079) * density(01) 
  pd(30,01) = pd(30,01) - rrt(079) * density(30) 
  pd(30,30) = pd(30,30) - rrt(079) * density(01) 
  pd(34,01) = pd(34,01) + rrt(079) * density(30) 
  pd(34,30) = pd(34,30) + rrt(079) * density(01) 
  pd(40,01) = pd(40,01) + rrt(079) * density(30) 
  pd(40,30) = pd(40,30) + rrt(079) * density(01) 
  pd(01,01) = pd(01,01) + rrt(080) * density(30) 
  pd(01,30) = pd(01,30) + rrt(080) * density(01) 
  pd(08,01) = pd(08,01) + rrt(080) * density(30) 
  pd(08,30) = pd(08,30) + rrt(080) * density(01) 
  pd(17,01) = pd(17,01) + rrt(080) * density(30) 
  pd(17,30) = pd(17,30) + rrt(080) * density(01) 
  pd(30,01) = pd(30,01) - rrt(080) * density(30) 
  pd(30,30) = pd(30,30) - rrt(080) * density(01) 
  pd(01,01) = pd(01,01) + rrt(081) * density(30) 
  pd(01,30) = pd(01,30) + rrt(081) * density(01) 
  pd(03,01) = pd(03,01) + rrt(081) * density(30) 
  pd(03,30) = pd(03,30) + rrt(081) * density(01) 
  pd(16,01) = pd(16,01) + rrt(081) * density(30) 
  pd(16,30) = pd(16,30) + rrt(081) * density(01) 
  pd(30,01) = pd(30,01) - rrt(081) * density(30) 
  pd(30,30) = pd(30,30) - rrt(081) * density(01) 
  pd(01,01) = pd(01,01) + rrt(082) * density(30) 
  pd(01,30) = pd(01,30) + rrt(082) * density(01) 
  pd(30,01) = pd(30,01) - rrt(082) * density(30) 
  pd(30,30) = pd(30,30) - rrt(082) * density(01) 
  pd(42,01) = pd(42,01) + rrt(082) * density(30) 
  pd(42,30) = pd(42,30) + rrt(082) * density(01) 
  pd(44,01) = pd(44,01) + rrt(082) * density(30) 
  pd(44,30) = pd(44,30) + rrt(082) * density(01) 
  pd(01,01) = pd(01,01) + rrt(083) * density(45) 
  pd(01,45) = pd(01,45) + rrt(083) * density(01) 
  pd(33,01) = pd(33,01) + rrt(083) * density(45) 
  pd(33,45) = pd(33,45) + rrt(083) * density(01) 
  pd(45,01) = pd(45,01) - rrt(083) * density(45) 
  pd(45,45) = pd(45,45) - rrt(083) * density(01) 
  pd(01,01) = pd(01,01) + rrt(084) * density(45) 
  pd(01,45) = pd(01,45) + rrt(084) * density(01) 
  pd(10,01) = pd(10,01) + rrt(084) * density(45) 
  pd(10,45) = pd(10,45) + rrt(084) * density(01) 
  pd(38,01) = pd(38,01) + rrt(084) * density(45) 
  pd(38,45) = pd(38,45) + rrt(084) * density(01) 
  pd(45,01) = pd(45,01) - rrt(084) * density(45) 
  pd(45,45) = pd(45,45) - rrt(084) * density(01) 
  pd(01,01) = pd(01,01) + rrt(085) * density(45) 
  pd(01,45) = pd(01,45) + rrt(085) * density(01) 
  pd(17,01) = pd(17,01) + rrt(085) * density(45) 
  pd(17,45) = pd(17,45) + rrt(085) * density(01) 
  pd(29,01) = pd(29,01) + rrt(085) * density(45) 
  pd(29,45) = pd(29,45) + rrt(085) * density(01) 
  pd(45,01) = pd(45,01) - rrt(085) * density(45) 
  pd(45,45) = pd(45,45) - rrt(085) * density(01) 
  pd(01,01) = pd(01,01) + rrt(086) * density(45) 
  pd(01,45) = pd(01,45) + rrt(086) * density(01) 
  pd(03,01) = pd(03,01) + rrt(086) * density(45) 
  pd(03,45) = pd(03,45) + rrt(086) * density(01) 
  pd(08,01) = pd(08,01) + rrt(086) * density(45) 
  pd(08,45) = pd(08,45) + rrt(086) * density(01) 
  pd(45,01) = pd(45,01) - rrt(086) * density(45) 
  pd(45,45) = pd(45,45) - rrt(086) * density(01) 
  pd(01,01) = pd(01,01) + rrt(087) * density(45) 
  pd(01,45) = pd(01,45) + rrt(087) * density(01) 
  pd(16,01) = pd(16,01) + rrt(087) * density(45) 
  pd(16,45) = pd(16,45) + rrt(087) * density(01) 
  pd(42,01) = pd(42,01) + rrt(087) * density(45) 
  pd(42,45) = pd(42,45) + rrt(087) * density(01) 
  pd(45,01) = pd(45,01) - rrt(087) * density(45) 
  pd(45,45) = pd(45,45) - rrt(087) * density(01) 
  pd(01,01) = pd(01,01) + rrt(088) * density(05) 
  pd(01,05) = pd(01,05) + rrt(088) * density(01) 
  pd(05,01) = pd(05,01) - rrt(088) * density(05) 
  pd(05,05) = pd(05,05) - rrt(088) * density(01) 
  pd(38,01) = pd(38,01) + rrt(088) * density(05) 
  pd(38,05) = pd(38,05) + rrt(088) * density(01) 
  pd(01,01) = pd(01,01) + rrt(089) * density(05) 
  pd(01,05) = pd(01,05) + rrt(089) * density(01) 
  pd(05,01) = pd(05,01) - rrt(089) * density(05) 
  pd(05,05) = pd(05,05) - rrt(089) * density(01) 
  pd(10,01) = pd(10,01) + rrt(089) * density(05) 
  pd(10,05) = pd(10,05) + rrt(089) * density(01) 
  pd(40,01) = pd(40,01) + rrt(089) * density(05) 
  pd(40,05) = pd(40,05) + rrt(089) * density(01) 
  pd(01,01) = pd(01,01) + rrt(090) * density(05) 
  pd(01,05) = pd(01,05) + rrt(090) * density(01) 
  pd(03,01) = pd(03,01) + rrt(090) * density(05) 
  pd(03,05) = pd(03,05) + rrt(090) * density(01) 
  pd(05,01) = pd(05,01) - rrt(090) * density(05) 
  pd(05,05) = pd(05,05) - rrt(090) * density(01) 
  pd(29,01) = pd(29,01) + rrt(090) * density(05) 
  pd(29,05) = pd(29,05) + rrt(090) * density(01) 
  pd(01,01) = pd(01,01) + rrt(091) * density(05) 
  pd(01,05) = pd(01,05) + rrt(091) * density(01) 
  pd(05,01) = pd(05,01) - rrt(091) * density(05) 
  pd(05,05) = pd(05,05) - rrt(091) * density(01) 
  pd(08,01) = pd(08,01) + rrt(091) * density(05) 
  pd(08,05) = pd(08,05) + rrt(091) * density(01) 
  pd(42,01) = pd(42,01) + rrt(091) * density(05) 
  pd(42,05) = pd(42,05) + rrt(091) * density(01) 
  pd(01,01) = pd(01,01) + rrt(092) * density(05) 
  pd(01,05) = pd(01,05) + rrt(092) * density(01) 
  pd(05,01) = pd(05,01) - rrt(092) * density(05) 
  pd(05,05) = pd(05,05) - rrt(092) * density(01) 
  pd(07,01) = pd(07,01) + rrt(092) * density(05) 
  pd(07,05) = pd(07,05) + rrt(092) * density(01) 
  pd(32,01) = pd(32,01) + rrt(092) * density(05) 
  pd(32,05) = pd(32,05) + rrt(092) * density(01) 
  pd(01,01) = pd(01,01) + rrt(093) * density(32) 
  pd(01,32) = pd(01,32) + rrt(093) * density(01) 
  pd(32,01) = pd(32,01) - rrt(093) * density(32) 
  pd(32,32) = pd(32,32) - rrt(093) * density(01) 
  pd(40,01) = pd(40,01) + rrt(093) * density(32) 
  pd(40,32) = pd(40,32) + rrt(093) * density(01) 
  pd(01,01) = pd(01,01) + rrt(094) * density(32) 
  pd(01,32) = pd(01,32) + rrt(094) * density(01) 
  pd(29,01) = pd(29,01) + rrt(094) * density(32) 
  pd(29,32) = pd(29,32) + rrt(094) * density(01) 
  pd(32,01) = pd(32,01) - rrt(094) * density(32) 
  pd(32,32) = pd(32,32) - rrt(094) * density(01) 
  pd(42,01) = pd(42,01) + rrt(094) * density(32) 
  pd(42,32) = pd(42,32) + rrt(094) * density(01) 
  pd(01,01) = pd(01,01) + rrt(095) * density(20) 
  pd(01,20) = pd(01,20) + rrt(095) * density(01) 
  pd(20,01) = pd(20,01) - rrt(095) * density(20) 
  pd(20,20) = pd(20,20) - rrt(095) * density(01) 
  pd(25,01) = pd(25,01) + rrt(095) * density(20) 
  pd(25,20) = pd(25,20) + rrt(095) * density(01) 
  pd(01,01) = pd(01,01) + rrt(096) * density(20) 
  pd(01,20) = pd(01,20) + rrt(096) * density(01) 
  pd(10,01) = pd(10,01) + rrt(096) * density(20) 
  pd(10,20) = pd(10,20) + rrt(096) * density(01) 
  pd(20,01) = pd(20,01) - rrt(096) * density(20) 
  pd(20,20) = pd(20,20) - rrt(096) * density(01) 
  pd(23,01) = pd(23,01) + rrt(096) * density(20) 
  pd(23,20) = pd(23,20) + rrt(096) * density(01) 
  pd(01,01) = pd(01,01) + rrt(097) * density(20) 
  pd(01,20) = pd(01,20) + rrt(097) * density(01) 
  pd(19,01) = pd(19,01) + rrt(097) * density(20) 
  pd(19,20) = pd(19,20) + rrt(097) * density(01) 
  pd(20,01) = pd(20,01) - rrt(097) * density(20) 
  pd(20,20) = pd(20,20) - rrt(097) * density(01) 
  pd(34,01) = pd(34,01) + rrt(097) * density(20) 
  pd(34,20) = pd(34,20) + rrt(097) * density(01) 
  pd(01,01) = pd(01,01) + rrt(098) * density(20) 
  pd(01,20) = pd(01,20) + rrt(098) * density(01) 
  pd(10,01) = pd(10,01) + rrt(098) * density(20) 
  pd(10,20) = pd(10,20) + rrt(098) * density(01) 
  pd(20,01) = pd(20,01) - rrt(098) * density(20) 
  pd(20,20) = pd(20,20) - rrt(098) * density(01) 
  pd(22,01) = pd(22,01) + rrt(098) * density(20) 
  pd(22,20) = pd(22,20) + rrt(098) * density(01) 
  pd(34,01) = pd(34,01) + rrt(098) * density(20) 
  pd(34,20) = pd(34,20) + rrt(098) * density(01) 
  pd(01,01) = pd(01,01) + rrt(099) * density(20) 
  pd(01,20) = pd(01,20) + rrt(099) * density(01) 
  pd(11,01) = pd(11,01) + rrt(099) * density(20) 
  pd(11,20) = pd(11,20) + rrt(099) * density(01) 
  pd(20,01) = pd(20,01) - rrt(099) * density(20) 
  pd(20,20) = pd(20,20) - rrt(099) * density(01) 
  pd(34,01) = pd(34,01) + rrt(099) * density(20) * 2.0d0
  pd(34,20) = pd(34,20) + rrt(099) * density(01) * 2.0d0
  pd(01,01) = pd(01,01) + rrt(100) * density(20) 
  pd(01,20) = pd(01,20) + rrt(100) * density(01) 
  pd(16,01) = pd(16,01) + rrt(100) * density(20) 
  pd(16,20) = pd(16,20) + rrt(100) * density(01) 
  pd(20,01) = pd(20,01) - rrt(100) * density(20) 
  pd(20,20) = pd(20,20) - rrt(100) * density(01) 
  pd(41,01) = pd(41,01) + rrt(100) * density(20) 
  pd(41,20) = pd(41,20) + rrt(100) * density(01) 
  pd(01,01) = pd(01,01) + rrt(101) * density(20) 
  pd(01,20) = pd(01,20) + rrt(101) * density(01) 
  pd(20,01) = pd(20,01) - rrt(101) * density(20) 
  pd(20,20) = pd(20,20) - rrt(101) * density(01) 
  pd(33,01) = pd(33,01) + rrt(101) * density(20) 
  pd(33,20) = pd(33,20) + rrt(101) * density(01) 
  pd(44,01) = pd(44,01) + rrt(101) * density(20) 
  pd(44,20) = pd(44,20) + rrt(101) * density(01) 
  pd(01,01) = pd(01,01) + rrt(102) * density(20) 
  pd(01,20) = pd(01,20) + rrt(102) * density(01) 
  pd(17,01) = pd(17,01) + rrt(102) * density(20) 
  pd(17,20) = pd(17,20) + rrt(102) * density(01) 
  pd(20,01) = pd(20,01) - rrt(102) * density(20) 
  pd(20,20) = pd(20,20) - rrt(102) * density(01) 
  pd(30,01) = pd(30,01) + rrt(102) * density(20) 
  pd(30,20) = pd(30,20) + rrt(102) * density(01) 
  pd(01,01) = pd(01,01) + rrt(103) * density(20) 
  pd(01,20) = pd(01,20) + rrt(103) * density(01) 
  pd(03,01) = pd(03,01) + rrt(103) * density(20) 
  pd(03,20) = pd(03,20) + rrt(103) * density(01) 
  pd(18,01) = pd(18,01) + rrt(103) * density(20) 
  pd(18,20) = pd(18,20) + rrt(103) * density(01) 
  pd(20,01) = pd(20,01) - rrt(103) * density(20) 
  pd(20,20) = pd(20,20) - rrt(103) * density(01) 
  pd(01,01) = pd(01,01) + rrt(104) * density(27) 
  pd(01,27) = pd(01,27) + rrt(104) * density(01) 
  pd(23,01) = pd(23,01) + rrt(104) * density(27) 
  pd(23,27) = pd(23,27) + rrt(104) * density(01) 
  pd(27,01) = pd(27,01) - rrt(104) * density(27) 
  pd(27,27) = pd(27,27) - rrt(104) * density(01) 
  pd(01,01) = pd(01,01) + rrt(105) * density(27) 
  pd(01,27) = pd(01,27) + rrt(105) * density(01) 
  pd(10,01) = pd(10,01) + rrt(105) * density(27) 
  pd(10,27) = pd(10,27) + rrt(105) * density(01) 
  pd(19,01) = pd(19,01) + rrt(105) * density(27) 
  pd(19,27) = pd(19,27) + rrt(105) * density(01) 
  pd(27,01) = pd(27,01) - rrt(105) * density(27) 
  pd(27,27) = pd(27,27) - rrt(105) * density(01) 
  pd(01,01) = pd(01,01) + rrt(106) * density(27) 
  pd(01,27) = pd(01,27) + rrt(106) * density(01) 
  pd(22,01) = pd(22,01) + rrt(106) * density(27) 
  pd(22,27) = pd(22,27) + rrt(106) * density(01) 
  pd(27,01) = pd(27,01) - rrt(106) * density(27) 
  pd(27,27) = pd(27,27) - rrt(106) * density(01) 
  pd(34,01) = pd(34,01) + rrt(106) * density(27) 
  pd(34,27) = pd(34,27) + rrt(106) * density(01) 
  pd(01,01) = pd(01,01) + rrt(107) * density(27) 
  pd(01,27) = pd(01,27) + rrt(107) * density(01) 
  pd(10,01) = pd(10,01) + rrt(107) * density(27) 
  pd(10,27) = pd(10,27) + rrt(107) * density(01) 
  pd(11,01) = pd(11,01) + rrt(107) * density(27) 
  pd(11,27) = pd(11,27) + rrt(107) * density(01) 
  pd(27,01) = pd(27,01) - rrt(107) * density(27) 
  pd(27,27) = pd(27,27) - rrt(107) * density(01) 
  pd(34,01) = pd(34,01) + rrt(107) * density(27) 
  pd(34,27) = pd(34,27) + rrt(107) * density(01) 
  pd(01,01) = pd(01,01) + rrt(108) * density(27) 
  pd(01,27) = pd(01,27) + rrt(108) * density(01) 
  pd(08,01) = pd(08,01) + rrt(108) * density(27) 
  pd(08,27) = pd(08,27) + rrt(108) * density(01) 
  pd(27,01) = pd(27,01) - rrt(108) * density(27) 
  pd(27,27) = pd(27,27) - rrt(108) * density(01) 
  pd(41,01) = pd(41,01) + rrt(108) * density(27) 
  pd(41,27) = pd(41,27) + rrt(108) * density(01) 
  pd(01,01) = pd(01,01) + rrt(109) * density(27) 
  pd(01,27) = pd(01,27) + rrt(109) * density(01) 
  pd(16,01) = pd(16,01) + rrt(109) * density(27) 
  pd(16,27) = pd(16,27) + rrt(109) * density(01) 
  pd(27,01) = pd(27,01) - rrt(109) * density(27) 
  pd(27,27) = pd(27,27) - rrt(109) * density(01) 
  pd(33,01) = pd(33,01) + rrt(109) * density(27) 
  pd(33,27) = pd(33,27) + rrt(109) * density(01) 
  pd(01,01) = pd(01,01) + rrt(110) * density(27) 
  pd(01,27) = pd(01,27) + rrt(110) * density(01) 
  pd(27,01) = pd(27,01) - rrt(110) * density(27) 
  pd(27,27) = pd(27,27) - rrt(110) * density(01) 
  pd(38,01) = pd(38,01) + rrt(110) * density(27) 
  pd(38,27) = pd(38,27) + rrt(110) * density(01) 
  pd(44,01) = pd(44,01) + rrt(110) * density(27) 
  pd(44,27) = pd(44,27) + rrt(110) * density(01) 
  pd(01,01) = pd(01,01) + rrt(111) * density(27) 
  pd(01,27) = pd(01,27) + rrt(111) * density(01) 
  pd(05,01) = pd(05,01) + rrt(111) * density(27) 
  pd(05,27) = pd(05,27) + rrt(111) * density(01) 
  pd(14,01) = pd(14,01) + rrt(111) * density(27) 
  pd(14,27) = pd(14,27) + rrt(111) * density(01) 
  pd(27,01) = pd(27,01) - rrt(111) * density(27) 
  pd(27,27) = pd(27,27) - rrt(111) * density(01) 
  pd(01,01) = pd(01,01) + rrt(112) * density(27) 
  pd(01,27) = pd(01,27) + rrt(112) * density(01) 
  pd(17,01) = pd(17,01) + rrt(112) * density(27) 
  pd(17,27) = pd(17,27) + rrt(112) * density(01) 
  pd(27,01) = pd(27,01) - rrt(112) * density(27) 
  pd(27,27) = pd(27,27) - rrt(112) * density(01) 
  pd(45,01) = pd(45,01) + rrt(112) * density(27) 
  pd(45,27) = pd(45,27) + rrt(112) * density(01) 
  pd(01,01) = pd(01,01) + rrt(113) * density(27) 
  pd(01,27) = pd(01,27) + rrt(113) * density(01) 
  pd(03,01) = pd(03,01) + rrt(113) * density(27) 
  pd(03,27) = pd(03,27) + rrt(113) * density(01) 
  pd(27,01) = pd(27,01) - rrt(113) * density(27) 
  pd(27,27) = pd(27,27) - rrt(113) * density(01) 
  pd(30,01) = pd(30,01) + rrt(113) * density(27) 
  pd(30,27) = pd(30,27) + rrt(113) * density(01) 
  pd(01,01) = pd(01,01) + rrt(114) * density(27) 
  pd(01,27) = pd(01,27) + rrt(114) * density(01) 
  pd(18,01) = pd(18,01) + rrt(114) * density(27) 
  pd(18,27) = pd(18,27) + rrt(114) * density(01) 
  pd(27,01) = pd(27,01) - rrt(114) * density(27) 
  pd(27,27) = pd(27,27) - rrt(114) * density(01) 
  pd(42,01) = pd(42,01) + rrt(114) * density(27) 
  pd(42,27) = pd(42,27) + rrt(114) * density(01) 
  pd(01,01) = pd(01,01) + rrt(115) * density(31) 
  pd(01,31) = pd(01,31) + rrt(115) * density(01) 
  pd(19,01) = pd(19,01) + rrt(115) * density(31) 
  pd(19,31) = pd(19,31) + rrt(115) * density(01) 
  pd(31,01) = pd(31,01) - rrt(115) * density(31) 
  pd(31,31) = pd(31,31) - rrt(115) * density(01) 
  pd(01,01) = pd(01,01) + rrt(116) * density(31) 
  pd(01,31) = pd(01,31) + rrt(116) * density(01) 
  pd(10,01) = pd(10,01) + rrt(116) * density(31) 
  pd(10,31) = pd(10,31) + rrt(116) * density(01) 
  pd(22,01) = pd(22,01) + rrt(116) * density(31) 
  pd(22,31) = pd(22,31) + rrt(116) * density(01) 
  pd(31,01) = pd(31,01) - rrt(116) * density(31) 
  pd(31,31) = pd(31,31) - rrt(116) * density(01) 
  pd(01,01) = pd(01,01) + rrt(117) * density(31) 
  pd(01,31) = pd(01,31) + rrt(117) * density(01) 
  pd(11,01) = pd(11,01) + rrt(117) * density(31) 
  pd(11,31) = pd(11,31) + rrt(117) * density(01) 
  pd(31,01) = pd(31,01) - rrt(117) * density(31) 
  pd(31,31) = pd(31,31) - rrt(117) * density(01) 
  pd(34,01) = pd(34,01) + rrt(117) * density(31) 
  pd(34,31) = pd(34,31) + rrt(117) * density(01) 
  pd(01,01) = pd(01,01) + rrt(118) * density(31) 
  pd(01,31) = pd(01,31) + rrt(118) * density(01) 
  pd(29,01) = pd(29,01) + rrt(118) * density(31) 
  pd(29,31) = pd(29,31) + rrt(118) * density(01) 
  pd(31,01) = pd(31,01) - rrt(118) * density(31) 
  pd(31,31) = pd(31,31) - rrt(118) * density(01) 
  pd(41,01) = pd(41,01) + rrt(118) * density(31) 
  pd(41,31) = pd(41,31) + rrt(118) * density(01) 
  pd(01,01) = pd(01,01) + rrt(119) * density(31) 
  pd(01,31) = pd(01,31) + rrt(119) * density(01) 
  pd(08,01) = pd(08,01) + rrt(119) * density(31) 
  pd(08,31) = pd(08,31) + rrt(119) * density(01) 
  pd(31,01) = pd(31,01) - rrt(119) * density(31) 
  pd(31,31) = pd(31,31) - rrt(119) * density(01) 
  pd(33,01) = pd(33,01) + rrt(119) * density(31) 
  pd(33,31) = pd(33,31) + rrt(119) * density(01) 
  pd(01,01) = pd(01,01) + rrt(120) * density(31) 
  pd(01,31) = pd(01,31) + rrt(120) * density(01) 
  pd(16,01) = pd(16,01) + rrt(120) * density(31) 
  pd(16,31) = pd(16,31) + rrt(120) * density(01) 
  pd(31,01) = pd(31,01) - rrt(120) * density(31) 
  pd(31,31) = pd(31,31) - rrt(120) * density(01) 
  pd(38,01) = pd(38,01) + rrt(120) * density(31) 
  pd(38,31) = pd(38,31) + rrt(120) * density(01) 
  pd(01,01) = pd(01,01) + rrt(121) * density(31) 
  pd(01,31) = pd(01,31) + rrt(121) * density(01) 
  pd(31,01) = pd(31,01) - rrt(121) * density(31) 
  pd(31,31) = pd(31,31) - rrt(121) * density(01) 
  pd(40,01) = pd(40,01) + rrt(121) * density(31) 
  pd(40,31) = pd(40,31) + rrt(121) * density(01) 
  pd(44,01) = pd(44,01) + rrt(121) * density(31) 
  pd(44,31) = pd(44,31) + rrt(121) * density(01) 
  pd(01,01) = pd(01,01) + rrt(122) * density(31) 
  pd(01,31) = pd(01,31) + rrt(122) * density(01) 
  pd(14,01) = pd(14,01) + rrt(122) * density(31) 
  pd(14,31) = pd(14,31) + rrt(122) * density(01) 
  pd(31,01) = pd(31,01) - rrt(122) * density(31) 
  pd(31,31) = pd(31,31) - rrt(122) * density(01) 
  pd(32,01) = pd(32,01) + rrt(122) * density(31) 
  pd(32,31) = pd(32,31) + rrt(122) * density(01) 
  pd(01,01) = pd(01,01) + rrt(123) * density(31) 
  pd(01,31) = pd(01,31) + rrt(123) * density(01) 
  pd(05,01) = pd(05,01) + rrt(123) * density(31) 
  pd(05,31) = pd(05,31) + rrt(123) * density(01) 
  pd(17,01) = pd(17,01) + rrt(123) * density(31) 
  pd(17,31) = pd(17,31) + rrt(123) * density(01) 
  pd(31,01) = pd(31,01) - rrt(123) * density(31) 
  pd(31,31) = pd(31,31) - rrt(123) * density(01) 
  pd(01,01) = pd(01,01) + rrt(124) * density(31) 
  pd(01,31) = pd(01,31) + rrt(124) * density(01) 
  pd(03,01) = pd(03,01) + rrt(124) * density(31) 
  pd(03,31) = pd(03,31) + rrt(124) * density(01) 
  pd(31,01) = pd(31,01) - rrt(124) * density(31) 
  pd(31,31) = pd(31,31) - rrt(124) * density(01) 
  pd(45,01) = pd(45,01) + rrt(124) * density(31) 
  pd(45,31) = pd(45,31) + rrt(124) * density(01) 
  pd(01,01) = pd(01,01) + rrt(125) * density(31) 
  pd(01,31) = pd(01,31) + rrt(125) * density(01) 
  pd(30,01) = pd(30,01) + rrt(125) * density(31) 
  pd(30,31) = pd(30,31) + rrt(125) * density(01) 
  pd(31,01) = pd(31,01) - rrt(125) * density(31) 
  pd(31,31) = pd(31,31) - rrt(125) * density(01) 
  pd(42,01) = pd(42,01) + rrt(125) * density(31) 
  pd(42,31) = pd(42,31) + rrt(125) * density(01) 
  pd(01,01) = pd(01,01) + rrt(126) * density(21) 
  pd(01,21) = pd(01,21) + rrt(126) * density(01) 
  pd(21,01) = pd(21,01) - rrt(126) * density(21) 
  pd(21,21) = pd(21,21) - rrt(126) * density(01) 
  pd(22,01) = pd(22,01) + rrt(126) * density(21) 
  pd(22,21) = pd(22,21) + rrt(126) * density(01) 
  pd(01,01) = pd(01,01) + rrt(127) * density(21) 
  pd(01,21) = pd(01,21) + rrt(127) * density(01) 
  pd(10,01) = pd(10,01) + rrt(127) * density(21) 
  pd(10,21) = pd(10,21) + rrt(127) * density(01) 
  pd(11,01) = pd(11,01) + rrt(127) * density(21) 
  pd(11,21) = pd(11,21) + rrt(127) * density(01) 
  pd(21,01) = pd(21,01) - rrt(127) * density(21) 
  pd(21,21) = pd(21,21) - rrt(127) * density(01) 
  pd(01,01) = pd(01,01) + rrt(128) * density(21) 
  pd(01,21) = pd(01,21) + rrt(128) * density(01) 
  pd(21,01) = pd(21,01) - rrt(128) * density(21) 
  pd(21,21) = pd(21,21) - rrt(128) * density(01) 
  pd(29,01) = pd(29,01) + rrt(128) * density(21) 
  pd(29,21) = pd(29,21) + rrt(128) * density(01) 
  pd(33,01) = pd(33,01) + rrt(128) * density(21) 
  pd(33,21) = pd(33,21) + rrt(128) * density(01) 
  pd(01,01) = pd(01,01) + rrt(129) * density(21) 
  pd(01,21) = pd(01,21) + rrt(129) * density(01) 
  pd(08,01) = pd(08,01) + rrt(129) * density(21) 
  pd(08,21) = pd(08,21) + rrt(129) * density(01) 
  pd(21,01) = pd(21,01) - rrt(129) * density(21) 
  pd(21,21) = pd(21,21) - rrt(129) * density(01) 
  pd(38,01) = pd(38,01) + rrt(129) * density(21) 
  pd(38,21) = pd(38,21) + rrt(129) * density(01) 
  pd(01,01) = pd(01,01) + rrt(130) * density(21) 
  pd(01,21) = pd(01,21) + rrt(130) * density(01) 
  pd(16,01) = pd(16,01) + rrt(130) * density(21) 
  pd(16,21) = pd(16,21) + rrt(130) * density(01) 
  pd(21,01) = pd(21,01) - rrt(130) * density(21) 
  pd(21,21) = pd(21,21) - rrt(130) * density(01) 
  pd(40,01) = pd(40,01) + rrt(130) * density(21) 
  pd(40,21) = pd(40,21) + rrt(130) * density(01) 
  pd(01,01) = pd(01,01) + rrt(131) * density(21) 
  pd(01,21) = pd(01,21) + rrt(131) * density(01) 
  pd(17,01) = pd(17,01) + rrt(131) * density(21) 
  pd(17,21) = pd(17,21) + rrt(131) * density(01) 
  pd(21,01) = pd(21,01) - rrt(131) * density(21) 
  pd(21,21) = pd(21,21) - rrt(131) * density(01) 
  pd(32,01) = pd(32,01) + rrt(131) * density(21) 
  pd(32,21) = pd(32,21) + rrt(131) * density(01) 
  pd(01,01) = pd(01,01) + rrt(132) * density(21) 
  pd(01,21) = pd(01,21) + rrt(132) * density(01) 
  pd(03,01) = pd(03,01) + rrt(132) * density(21) 
  pd(03,21) = pd(03,21) + rrt(132) * density(01) 
  pd(05,01) = pd(05,01) + rrt(132) * density(21) 
  pd(05,21) = pd(05,21) + rrt(132) * density(01) 
  pd(21,01) = pd(21,01) - rrt(132) * density(21) 
  pd(21,21) = pd(21,21) - rrt(132) * density(01) 
  pd(01,01) = pd(01,01) + rrt(133) * density(21) 
  pd(01,21) = pd(01,21) + rrt(133) * density(01) 
  pd(21,01) = pd(21,01) - rrt(133) * density(21) 
  pd(21,21) = pd(21,21) - rrt(133) * density(01) 
  pd(42,01) = pd(42,01) + rrt(133) * density(21) 
  pd(42,21) = pd(42,21) + rrt(133) * density(01) 
  pd(45,01) = pd(45,01) + rrt(133) * density(21) 
  pd(45,21) = pd(45,21) + rrt(133) * density(01) 
  pd(01,01) = pd(01,01) + rrt(134) * density(06) 
  pd(01,06) = pd(01,06) + rrt(134) * density(01) 
  pd(06,01) = pd(06,01) - rrt(134) * density(06) 
  pd(06,06) = pd(06,06) - rrt(134) * density(01) 
  pd(11,01) = pd(11,01) + rrt(134) * density(06) 
  pd(11,06) = pd(11,06) + rrt(134) * density(01) 
  pd(01,01) = pd(01,01) + rrt(135) * density(06) 
  pd(01,06) = pd(01,06) + rrt(135) * density(01) 
  pd(06,01) = pd(06,01) - rrt(135) * density(06) 
  pd(06,06) = pd(06,06) - rrt(135) * density(01) 
  pd(29,01) = pd(29,01) + rrt(135) * density(06) 
  pd(29,06) = pd(29,06) + rrt(135) * density(01) 
  pd(38,01) = pd(38,01) + rrt(135) * density(06) 
  pd(38,06) = pd(38,06) + rrt(135) * density(01) 
  pd(01,01) = pd(01,01) + rrt(136) * density(06) 
  pd(01,06) = pd(01,06) + rrt(136) * density(01) 
  pd(06,01) = pd(06,01) - rrt(136) * density(06) 
  pd(06,06) = pd(06,06) - rrt(136) * density(01) 
  pd(08,01) = pd(08,01) + rrt(136) * density(06) 
  pd(08,06) = pd(08,06) + rrt(136) * density(01) 
  pd(40,01) = pd(40,01) + rrt(136) * density(06) 
  pd(40,06) = pd(40,06) + rrt(136) * density(01) 
  pd(01,01) = pd(01,01) + rrt(137) * density(06) 
  pd(01,06) = pd(01,06) + rrt(137) * density(01) 
  pd(03,01) = pd(03,01) + rrt(137) * density(06) 
  pd(03,06) = pd(03,06) + rrt(137) * density(01) 
  pd(06,01) = pd(06,01) - rrt(137) * density(06) 
  pd(06,06) = pd(06,06) - rrt(137) * density(01) 
  pd(32,01) = pd(32,01) + rrt(137) * density(06) 
  pd(32,06) = pd(32,06) + rrt(137) * density(01) 
  pd(01,01) = pd(01,01) + rrt(138) * density(06) 
  pd(01,06) = pd(01,06) + rrt(138) * density(01) 
  pd(05,01) = pd(05,01) + rrt(138) * density(06) 
  pd(05,06) = pd(05,06) + rrt(138) * density(01) 
  pd(06,01) = pd(06,01) - rrt(138) * density(06) 
  pd(06,06) = pd(06,06) - rrt(138) * density(01) 
  pd(42,01) = pd(42,01) + rrt(138) * density(06) 
  pd(42,06) = pd(42,06) + rrt(138) * density(01) 
  pd(10,01) = pd(10,01) + rrt(139) * density(39) 
  pd(10,39) = pd(10,39) + rrt(139) * density(01) 
  pd(16,01) = pd(16,01) + rrt(139) * density(39) 
  pd(16,39) = pd(16,39) + rrt(139) * density(01) 
  pd(39,01) = pd(39,01) - rrt(139) * density(39) 
  pd(39,39) = pd(39,39) - rrt(139) * density(01) 
  pd(10,01) = pd(10,01) + rrt(140) * density(37) 
  pd(10,37) = pd(10,37) + rrt(140) * density(01) 
  pd(16,01) = pd(16,01) + rrt(140) * density(37) 
  pd(16,37) = pd(16,37) + rrt(140) * density(01) 
  pd(37,01) = pd(37,01) - rrt(140) * density(37) 
  pd(37,37) = pd(37,37) - rrt(140) * density(01) 
  pd(08,01) = pd(08,01) + rrt(141) * density(39) 
  pd(08,39) = pd(08,39) + rrt(141) * density(01) 
  pd(34,01) = pd(34,01) + rrt(141) * density(39) 
  pd(34,39) = pd(34,39) + rrt(141) * density(01) 
  pd(39,01) = pd(39,01) - rrt(141) * density(39) 
  pd(39,39) = pd(39,39) - rrt(141) * density(01) 
  pd(08,01) = pd(08,01) + rrt(142) * density(37) 
  pd(08,37) = pd(08,37) + rrt(142) * density(01) 
  pd(34,01) = pd(34,01) + rrt(142) * density(37) 
  pd(34,37) = pd(34,37) + rrt(142) * density(01) 
  pd(37,01) = pd(37,01) - rrt(142) * density(37) 
  pd(37,37) = pd(37,37) - rrt(142) * density(01) 
  pd(10,01) = pd(10,01) + rrt(143) * density(39) 
  pd(10,39) = pd(10,39) + rrt(143) * density(01) 
  pd(29,01) = pd(29,01) + rrt(143) * density(39) 
  pd(29,39) = pd(29,39) + rrt(143) * density(01) 
  pd(34,01) = pd(34,01) + rrt(143) * density(39) 
  pd(34,39) = pd(34,39) + rrt(143) * density(01) 
  pd(39,01) = pd(39,01) - rrt(143) * density(39) 
  pd(39,39) = pd(39,39) - rrt(143) * density(01) 
  pd(10,01) = pd(10,01) + rrt(144) * density(37) 
  pd(10,37) = pd(10,37) + rrt(144) * density(01) 
  pd(29,01) = pd(29,01) + rrt(144) * density(37) 
  pd(29,37) = pd(29,37) + rrt(144) * density(01) 
  pd(34,01) = pd(34,01) + rrt(144) * density(37) 
  pd(34,37) = pd(34,37) + rrt(144) * density(01) 
  pd(37,01) = pd(37,01) - rrt(144) * density(37) 
  pd(37,37) = pd(37,37) - rrt(144) * density(01) 
  pd(01,01) = pd(01,01) + rrt(145) * density(39) 
  pd(01,39) = pd(01,39) + rrt(145) * density(01) 
  pd(14,01) = pd(14,01) + rrt(145) * density(39) 
  pd(14,39) = pd(14,39) + rrt(145) * density(01) 
  pd(39,01) = pd(39,01) - rrt(145) * density(39) 
  pd(39,39) = pd(39,39) - rrt(145) * density(01) 
  pd(01,01) = pd(01,01) + rrt(146) * density(37) 
  pd(01,37) = pd(01,37) + rrt(146) * density(01) 
  pd(14,01) = pd(14,01) + rrt(146) * density(37) 
  pd(14,37) = pd(14,37) + rrt(146) * density(01) 
  pd(37,01) = pd(37,01) - rrt(146) * density(37) 
  pd(37,37) = pd(37,37) - rrt(146) * density(01) 
  pd(01,01) = pd(01,01) + rrt(147) * density(39) 
  pd(01,39) = pd(01,39) + rrt(147) * density(01) 
  pd(10,01) = pd(10,01) + rrt(147) * density(39) 
  pd(10,39) = pd(10,39) + rrt(147) * density(01) 
  pd(17,01) = pd(17,01) + rrt(147) * density(39) 
  pd(17,39) = pd(17,39) + rrt(147) * density(01) 
  pd(39,01) = pd(39,01) - rrt(147) * density(39) 
  pd(39,39) = pd(39,39) - rrt(147) * density(01) 
  pd(01,01) = pd(01,01) + rrt(148) * density(37) 
  pd(01,37) = pd(01,37) + rrt(148) * density(01) 
  pd(10,01) = pd(10,01) + rrt(148) * density(37) 
  pd(10,37) = pd(10,37) + rrt(148) * density(01) 
  pd(17,01) = pd(17,01) + rrt(148) * density(37) 
  pd(17,37) = pd(17,37) + rrt(148) * density(01) 
  pd(37,01) = pd(37,01) - rrt(148) * density(37) 
  pd(37,37) = pd(37,37) - rrt(148) * density(01) 
  pd(01,01) = pd(01,01) + rrt(149) * density(39) 
  pd(01,39) = pd(01,39) + rrt(149) * density(01) 
  pd(03,01) = pd(03,01) + rrt(149) * density(39) 
  pd(03,39) = pd(03,39) + rrt(149) * density(01) 
  pd(34,01) = pd(34,01) + rrt(149) * density(39) 
  pd(34,39) = pd(34,39) + rrt(149) * density(01) 
  pd(39,01) = pd(39,01) - rrt(149) * density(39) 
  pd(39,39) = pd(39,39) - rrt(149) * density(01) 
  pd(01,01) = pd(01,01) + rrt(150) * density(37) 
  pd(01,37) = pd(01,37) + rrt(150) * density(01) 
  pd(03,01) = pd(03,01) + rrt(150) * density(37) 
  pd(03,37) = pd(03,37) + rrt(150) * density(01) 
  pd(34,01) = pd(34,01) + rrt(150) * density(37) 
  pd(34,37) = pd(34,37) + rrt(150) * density(01) 
  pd(37,01) = pd(37,01) - rrt(150) * density(37) 
  pd(37,37) = pd(37,37) - rrt(150) * density(01) 
  pd(01,01) = pd(01,01) + rrt(151) * density(39) 
  pd(01,39) = pd(01,39) + rrt(151) * density(01) 
  pd(10,01) = pd(10,01) + rrt(151) * density(39) 
  pd(10,39) = pd(10,39) + rrt(151) * density(01) 
  pd(34,01) = pd(34,01) + rrt(151) * density(39) 
  pd(34,39) = pd(34,39) + rrt(151) * density(01) 
  pd(39,01) = pd(39,01) - rrt(151) * density(39) 
  pd(39,39) = pd(39,39) - rrt(151) * density(01) 
  pd(42,01) = pd(42,01) + rrt(151) * density(39) 
  pd(42,39) = pd(42,39) + rrt(151) * density(01) 
  pd(01,01) = pd(01,01) + rrt(152) * density(37) 
  pd(01,37) = pd(01,37) + rrt(152) * density(01) 
  pd(10,01) = pd(10,01) + rrt(152) * density(37) 
  pd(10,37) = pd(10,37) + rrt(152) * density(01) 
  pd(34,01) = pd(34,01) + rrt(152) * density(37) 
  pd(34,37) = pd(34,37) + rrt(152) * density(01) 
  pd(37,01) = pd(37,01) - rrt(152) * density(37) 
  pd(37,37) = pd(37,37) - rrt(152) * density(01) 
  pd(42,01) = pd(42,01) + rrt(152) * density(37) 
  pd(42,37) = pd(42,37) + rrt(152) * density(01) 
  pd(10,01) = pd(10,01) + rrt(153) * density(35) 
  pd(10,35) = pd(10,35) + rrt(153) * density(01) 
  pd(30,01) = pd(30,01) + rrt(153) * density(35) 
  pd(30,35) = pd(30,35) + rrt(153) * density(01) 
  pd(35,01) = pd(35,01) - rrt(153) * density(35) 
  pd(35,35) = pd(35,35) - rrt(153) * density(01) 
  pd(09,01) = pd(09,01) - rrt(154) * density(09) 
  pd(09,09) = pd(09,09) - rrt(154) * density(01) 
  pd(10,01) = pd(10,01) + rrt(154) * density(09) 
  pd(10,09) = pd(10,09) + rrt(154) * density(01) 
  pd(30,01) = pd(30,01) + rrt(154) * density(09) 
  pd(30,09) = pd(30,09) + rrt(154) * density(01) 
  pd(34,01) = pd(34,01) + rrt(155) * density(35) 
  pd(34,35) = pd(34,35) + rrt(155) * density(01) 
  pd(35,01) = pd(35,01) - rrt(155) * density(35) 
  pd(35,35) = pd(35,35) - rrt(155) * density(01) 
  pd(45,01) = pd(45,01) + rrt(155) * density(35) 
  pd(45,35) = pd(45,35) + rrt(155) * density(01) 
  pd(09,01) = pd(09,01) - rrt(156) * density(09) 
  pd(09,09) = pd(09,09) - rrt(156) * density(01) 
  pd(34,01) = pd(34,01) + rrt(156) * density(09) 
  pd(34,09) = pd(34,09) + rrt(156) * density(01) 
  pd(45,01) = pd(45,01) + rrt(156) * density(09) 
  pd(45,09) = pd(45,09) + rrt(156) * density(01) 
  pd(05,01) = pd(05,01) + rrt(157) * density(35) 
  pd(05,35) = pd(05,35) + rrt(157) * density(01) 
  pd(10,01) = pd(10,01) + rrt(157) * density(35) 
  pd(10,35) = pd(10,35) + rrt(157) * density(01) 
  pd(34,01) = pd(34,01) + rrt(157) * density(35) 
  pd(34,35) = pd(34,35) + rrt(157) * density(01) 
  pd(35,01) = pd(35,01) - rrt(157) * density(35) 
  pd(35,35) = pd(35,35) - rrt(157) * density(01) 
  pd(05,01) = pd(05,01) + rrt(158) * density(09) 
  pd(05,09) = pd(05,09) + rrt(158) * density(01) 
  pd(09,01) = pd(09,01) - rrt(158) * density(09) 
  pd(09,09) = pd(09,09) - rrt(158) * density(01) 
  pd(10,01) = pd(10,01) + rrt(158) * density(09) 
  pd(10,09) = pd(10,09) + rrt(158) * density(01) 
  pd(34,01) = pd(34,01) + rrt(158) * density(09) 
  pd(34,09) = pd(34,09) + rrt(158) * density(01) 
  pd(32,01) = pd(32,01) + rrt(159) * density(35) 
  pd(32,35) = pd(32,35) + rrt(159) * density(01) 
  pd(34,01) = pd(34,01) + rrt(159) * density(35) * 2.0d0
  pd(34,35) = pd(34,35) + rrt(159) * density(01) * 2.0d0
  pd(35,01) = pd(35,01) - rrt(159) * density(35) 
  pd(35,35) = pd(35,35) - rrt(159) * density(01) 
  pd(09,01) = pd(09,01) - rrt(160) * density(09) 
  pd(09,09) = pd(09,09) - rrt(160) * density(01) 
  pd(32,01) = pd(32,01) + rrt(160) * density(09) 
  pd(32,09) = pd(32,09) + rrt(160) * density(01) 
  pd(34,01) = pd(34,01) + rrt(160) * density(09) * 2.0d0
  pd(34,09) = pd(34,09) + rrt(160) * density(01) * 2.0d0
  pd(08,01) = pd(08,01) + rrt(161) * density(35) 
  pd(08,35) = pd(08,35) + rrt(161) * density(01) 
  pd(35,01) = pd(35,01) - rrt(161) * density(35) 
  pd(35,35) = pd(35,35) - rrt(161) * density(01) 
  pd(44,01) = pd(44,01) + rrt(161) * density(35) 
  pd(44,35) = pd(44,35) + rrt(161) * density(01) 
  pd(08,01) = pd(08,01) + rrt(162) * density(09) 
  pd(08,09) = pd(08,09) + rrt(162) * density(01) 
  pd(09,01) = pd(09,01) - rrt(162) * density(09) 
  pd(09,09) = pd(09,09) - rrt(162) * density(01) 
  pd(44,01) = pd(44,01) + rrt(162) * density(09) 
  pd(44,09) = pd(44,09) + rrt(162) * density(01) 
  pd(16,01) = pd(16,01) + rrt(163) * density(35) * 2.0d0
  pd(16,35) = pd(16,35) + rrt(163) * density(01) * 2.0d0
  pd(35,01) = pd(35,01) - rrt(163) * density(35) 
  pd(35,35) = pd(35,35) - rrt(163) * density(01) 
  pd(09,01) = pd(09,01) - rrt(164) * density(09) 
  pd(09,09) = pd(09,09) - rrt(164) * density(01) 
  pd(16,01) = pd(16,01) + rrt(164) * density(09) * 2.0d0
  pd(16,09) = pd(16,09) + rrt(164) * density(01) * 2.0d0
  pd(05,01) = pd(05,01) + rrt(165) * density(15) 
  pd(05,15) = pd(05,15) + rrt(165) * density(01) 
  pd(10,01) = pd(10,01) + rrt(165) * density(15) 
  pd(10,15) = pd(10,15) + rrt(165) * density(01) 
  pd(15,01) = pd(15,01) - rrt(165) * density(15) 
  pd(15,15) = pd(15,15) - rrt(165) * density(01) 
  pd(02,01) = pd(02,01) - rrt(166) * density(02) 
  pd(02,02) = pd(02,02) - rrt(166) * density(01) 
  pd(05,01) = pd(05,01) + rrt(166) * density(02) 
  pd(05,02) = pd(05,02) + rrt(166) * density(01) 
  pd(10,01) = pd(10,01) + rrt(166) * density(02) 
  pd(10,02) = pd(10,02) + rrt(166) * density(01) 
  pd(15,01) = pd(15,01) - rrt(167) * density(15) 
  pd(15,15) = pd(15,15) - rrt(167) * density(01) 
  pd(32,01) = pd(32,01) + rrt(167) * density(15) 
  pd(32,15) = pd(32,15) + rrt(167) * density(01) 
  pd(34,01) = pd(34,01) + rrt(167) * density(15) 
  pd(34,15) = pd(34,15) + rrt(167) * density(01) 
  pd(02,01) = pd(02,01) - rrt(168) * density(02) 
  pd(02,02) = pd(02,02) - rrt(168) * density(01) 
  pd(32,01) = pd(32,01) + rrt(168) * density(02) 
  pd(32,02) = pd(32,02) + rrt(168) * density(01) 
  pd(34,01) = pd(34,01) + rrt(168) * density(02) 
  pd(34,02) = pd(34,02) + rrt(168) * density(01) 
  pd(10,01) = pd(10,01) + rrt(169) * density(15) * 2.0d0
  pd(10,15) = pd(10,15) + rrt(169) * density(01) * 2.0d0
  pd(15,01) = pd(15,01) - rrt(169) * density(15) 
  pd(15,15) = pd(15,15) - rrt(169) * density(01) 
  pd(32,01) = pd(32,01) + rrt(169) * density(15) 
  pd(32,15) = pd(32,15) + rrt(169) * density(01) 
  pd(02,01) = pd(02,01) - rrt(170) * density(02) 
  pd(02,02) = pd(02,02) - rrt(170) * density(01) 
  pd(10,01) = pd(10,01) + rrt(170) * density(02) * 2.0d0
  pd(10,02) = pd(10,02) + rrt(170) * density(01) * 2.0d0
  pd(32,01) = pd(32,01) + rrt(170) * density(02) 
  pd(32,02) = pd(32,02) + rrt(170) * density(01) 
  pd(15,01) = pd(15,01) - rrt(171) * density(15) 
  pd(15,15) = pd(15,15) - rrt(171) * density(01) 
  pd(16,01) = pd(16,01) + rrt(171) * density(15) 
  pd(16,15) = pd(16,15) + rrt(171) * density(01) 
  pd(29,01) = pd(29,01) + rrt(171) * density(15) 
  pd(29,15) = pd(29,15) + rrt(171) * density(01) 
  pd(02,01) = pd(02,01) - rrt(172) * density(02) 
  pd(02,02) = pd(02,02) - rrt(172) * density(01) 
  pd(16,01) = pd(16,01) + rrt(172) * density(02) 
  pd(16,02) = pd(16,02) + rrt(172) * density(01) 
  pd(29,01) = pd(29,01) + rrt(172) * density(02) 
  pd(29,02) = pd(29,02) + rrt(172) * density(01) 
  pd(08,01) = pd(08,01) + rrt(173) * density(15) * 2.0d0
  pd(08,15) = pd(08,15) + rrt(173) * density(01) * 2.0d0
  pd(15,01) = pd(15,01) - rrt(173) * density(15) 
  pd(15,15) = pd(15,15) - rrt(173) * density(01) 
  pd(02,01) = pd(02,01) - rrt(174) * density(02) 
  pd(02,02) = pd(02,02) - rrt(174) * density(01) 
  pd(08,01) = pd(08,01) + rrt(174) * density(02) * 2.0d0
  pd(08,02) = pd(08,02) + rrt(174) * density(01) * 2.0d0
  pd(26,01) = pd(26,01) - rrt(175) * density(26) 
  pd(26,26) = pd(26,26) - rrt(175) * density(01) 
  pd(29,01) = pd(29,01) + rrt(175) * density(26) * 2.0d0
  pd(29,26) = pd(29,26) + rrt(175) * density(01) * 2.0d0
  pd(29,01) = pd(29,01) + rrt(176) * density(36) * 2.0d0
  pd(29,36) = pd(29,36) + rrt(176) * density(01) * 2.0d0
  pd(36,01) = pd(36,01) - rrt(176) * density(36) 
  pd(36,36) = pd(36,36) - rrt(176) * density(01) 
  pd(28,01) = pd(28,01) - rrt(177) * density(28) 
  pd(28,28) = pd(28,28) - rrt(177) * density(01) 
  pd(29,01) = pd(29,01) + rrt(177) * density(28) * 2.0d0
  pd(29,28) = pd(29,28) + rrt(177) * density(01) * 2.0d0
  pd(04,01) = pd(04,01) - rrt(178) * density(04) 
  pd(04,04) = pd(04,04) - rrt(178) * density(01) 
  pd(10,01) = pd(10,01) + rrt(178) * density(04) 
  pd(10,04) = pd(10,04) + rrt(178) * density(01) 
  pd(27,01) = pd(27,01) + rrt(178) * density(04) 
  pd(27,04) = pd(27,04) + rrt(178) * density(01) 
  pd(10,01) = pd(10,01) + rrt(179) * density(24) 
  pd(10,24) = pd(10,24) + rrt(179) * density(01) 
  pd(24,01) = pd(24,01) - rrt(179) * density(24) 
  pd(24,24) = pd(24,24) - rrt(179) * density(01) 
  pd(27,01) = pd(27,01) + rrt(179) * density(24) 
  pd(27,24) = pd(27,24) + rrt(179) * density(01) 
  pd(04,01) = pd(04,01) - rrt(180) * density(04) 
  pd(04,04) = pd(04,04) - rrt(180) * density(01) 
  pd(31,01) = pd(31,01) + rrt(180) * density(04) 
  pd(31,04) = pd(31,04) + rrt(180) * density(01) 
  pd(34,01) = pd(34,01) + rrt(180) * density(04) 
  pd(34,04) = pd(34,04) + rrt(180) * density(01) 
  pd(24,01) = pd(24,01) - rrt(181) * density(24) 
  pd(24,24) = pd(24,24) - rrt(181) * density(01) 
  pd(31,01) = pd(31,01) + rrt(181) * density(24) 
  pd(31,24) = pd(31,24) + rrt(181) * density(01) 
  pd(34,01) = pd(34,01) + rrt(181) * density(24) 
  pd(34,24) = pd(34,24) + rrt(181) * density(01) 
  pd(04,01) = pd(04,01) - rrt(182) * density(04) 
  pd(04,04) = pd(04,04) - rrt(182) * density(01) 
  pd(06,01) = pd(06,01) + rrt(182) * density(04) 
  pd(06,04) = pd(06,04) + rrt(182) * density(01) 
  pd(34,01) = pd(34,01) + rrt(182) * density(04) * 2.0d0
  pd(34,04) = pd(34,04) + rrt(182) * density(01) * 2.0d0
  pd(06,01) = pd(06,01) + rrt(183) * density(24) 
  pd(06,24) = pd(06,24) + rrt(183) * density(01) 
  pd(24,01) = pd(24,01) - rrt(183) * density(24) 
  pd(24,24) = pd(24,24) - rrt(183) * density(01) 
  pd(34,01) = pd(34,01) + rrt(183) * density(24) * 2.0d0
  pd(34,24) = pd(34,24) + rrt(183) * density(01) * 2.0d0
  pd(04,01) = pd(04,01) - rrt(184) * density(04) 
  pd(04,04) = pd(04,04) - rrt(184) * density(01) 
  pd(08,01) = pd(08,01) + rrt(184) * density(04) 
  pd(08,04) = pd(08,04) + rrt(184) * density(01) 
  pd(18,01) = pd(18,01) + rrt(184) * density(04) 
  pd(18,04) = pd(18,04) + rrt(184) * density(01) 
  pd(08,01) = pd(08,01) + rrt(185) * density(24) 
  pd(08,24) = pd(08,24) + rrt(185) * density(01) 
  pd(18,01) = pd(18,01) + rrt(185) * density(24) 
  pd(18,24) = pd(18,24) + rrt(185) * density(01) 
  pd(24,01) = pd(24,01) - rrt(185) * density(24) 
  pd(24,24) = pd(24,24) - rrt(185) * density(01) 
  pd(04,01) = pd(04,01) - rrt(186) * density(04) 
  pd(04,04) = pd(04,04) - rrt(186) * density(01) 
  pd(16,01) = pd(16,01) + rrt(186) * density(04) 
  pd(16,04) = pd(16,04) + rrt(186) * density(01) 
  pd(30,01) = pd(30,01) + rrt(186) * density(04) 
  pd(30,04) = pd(30,04) + rrt(186) * density(01) 
  pd(16,01) = pd(16,01) + rrt(187) * density(24) 
  pd(16,24) = pd(16,24) + rrt(187) * density(01) 
  pd(24,01) = pd(24,01) - rrt(187) * density(24) 
  pd(24,24) = pd(24,24) - rrt(187) * density(01) 
  pd(30,01) = pd(30,01) + rrt(187) * density(24) 
  pd(30,24) = pd(30,24) + rrt(187) * density(01) 
  pd(04,01) = pd(04,01) - rrt(188) * density(04) 
  pd(04,04) = pd(04,04) - rrt(188) * density(01) 
  pd(44,01) = pd(44,01) + rrt(188) * density(04) 
  pd(44,04) = pd(44,04) + rrt(188) * density(01) 
  pd(45,01) = pd(45,01) + rrt(188) * density(04) 
  pd(45,04) = pd(45,04) + rrt(188) * density(01) 
  pd(24,01) = pd(24,01) - rrt(189) * density(24) 
  pd(24,24) = pd(24,24) - rrt(189) * density(01) 
  pd(44,01) = pd(44,01) + rrt(189) * density(24) 
  pd(44,24) = pd(44,24) + rrt(189) * density(01) 
  pd(45,01) = pd(45,01) + rrt(189) * density(24) 
  pd(45,24) = pd(45,24) + rrt(189) * density(01) 
  pd(01,01) = pd(01,01) + rrt(190) * density(35) 
  pd(01,35) = pd(01,35) + rrt(190) * density(01) 
  pd(12,01) = pd(12,01) + rrt(190) * density(35) 
  pd(12,35) = pd(12,35) + rrt(190) * density(01) 
  pd(35,01) = pd(35,01) - rrt(190) * density(35) 
  pd(35,35) = pd(35,35) - rrt(190) * density(01) 
  pd(01,01) = pd(01,01) + rrt(191) * density(09) 
  pd(01,09) = pd(01,09) + rrt(191) * density(01) 
  pd(09,01) = pd(09,01) - rrt(191) * density(09) 
  pd(09,09) = pd(09,09) - rrt(191) * density(01) 
  pd(12,01) = pd(12,01) + rrt(191) * density(09) 
  pd(12,09) = pd(12,09) + rrt(191) * density(01) 
  pd(01,01) = pd(01,01) + rrt(192) * density(35) 
  pd(01,35) = pd(01,35) + rrt(192) * density(01) 
  pd(10,01) = pd(10,01) + rrt(192) * density(35) 
  pd(10,35) = pd(10,35) + rrt(192) * density(01) 
  pd(35,01) = pd(35,01) - rrt(192) * density(35) 
  pd(35,35) = pd(35,35) - rrt(192) * density(01) 
  pd(41,01) = pd(41,01) + rrt(192) * density(35) 
  pd(41,35) = pd(41,35) + rrt(192) * density(01) 
  pd(01,01) = pd(01,01) + rrt(193) * density(09) 
  pd(01,09) = pd(01,09) + rrt(193) * density(01) 
  pd(09,01) = pd(09,01) - rrt(193) * density(09) 
  pd(09,09) = pd(09,09) - rrt(193) * density(01) 
  pd(10,01) = pd(10,01) + rrt(193) * density(09) 
  pd(10,09) = pd(10,09) + rrt(193) * density(01) 
  pd(41,01) = pd(41,01) + rrt(193) * density(09) 
  pd(41,09) = pd(41,09) + rrt(193) * density(01) 
  pd(01,01) = pd(01,01) + rrt(194) * density(35) 
  pd(01,35) = pd(01,35) + rrt(194) * density(01) 
  pd(33,01) = pd(33,01) + rrt(194) * density(35) 
  pd(33,35) = pd(33,35) + rrt(194) * density(01) 
  pd(34,01) = pd(34,01) + rrt(194) * density(35) 
  pd(34,35) = pd(34,35) + rrt(194) * density(01) 
  pd(35,01) = pd(35,01) - rrt(194) * density(35) 
  pd(35,35) = pd(35,35) - rrt(194) * density(01) 
  pd(01,01) = pd(01,01) + rrt(195) * density(09) 
  pd(01,09) = pd(01,09) + rrt(195) * density(01) 
  pd(09,01) = pd(09,01) - rrt(195) * density(09) 
  pd(09,09) = pd(09,09) - rrt(195) * density(01) 
  pd(33,01) = pd(33,01) + rrt(195) * density(09) 
  pd(33,09) = pd(33,09) + rrt(195) * density(01) 
  pd(34,01) = pd(34,01) + rrt(195) * density(09) 
  pd(34,09) = pd(34,09) + rrt(195) * density(01) 
  pd(01,01) = pd(01,01) + rrt(196) * density(35) 
  pd(01,35) = pd(01,35) + rrt(196) * density(01) 
  pd(10,01) = pd(10,01) + rrt(196) * density(35) 
  pd(10,35) = pd(10,35) + rrt(196) * density(01) 
  pd(34,01) = pd(34,01) + rrt(196) * density(35) 
  pd(34,35) = pd(34,35) + rrt(196) * density(01) 
  pd(35,01) = pd(35,01) - rrt(196) * density(35) 
  pd(35,35) = pd(35,35) - rrt(196) * density(01) 
  pd(38,01) = pd(38,01) + rrt(196) * density(35) 
  pd(38,35) = pd(38,35) + rrt(196) * density(01) 
  pd(01,01) = pd(01,01) + rrt(197) * density(09) 
  pd(01,09) = pd(01,09) + rrt(197) * density(01) 
  pd(09,01) = pd(09,01) - rrt(197) * density(09) 
  pd(09,09) = pd(09,09) - rrt(197) * density(01) 
  pd(10,01) = pd(10,01) + rrt(197) * density(09) 
  pd(10,09) = pd(10,09) + rrt(197) * density(01) 
  pd(34,01) = pd(34,01) + rrt(197) * density(09) 
  pd(34,09) = pd(34,09) + rrt(197) * density(01) 
  pd(38,01) = pd(38,01) + rrt(197) * density(09) 
  pd(38,09) = pd(38,09) + rrt(197) * density(01) 
  pd(01,01) = pd(01,01) + rrt(198) * density(35) 
  pd(01,35) = pd(01,35) + rrt(198) * density(01) 
  pd(34,01) = pd(34,01) + rrt(198) * density(35) * 2.0d0
  pd(34,35) = pd(34,35) + rrt(198) * density(01) * 2.0d0
  pd(35,01) = pd(35,01) - rrt(198) * density(35) 
  pd(35,35) = pd(35,35) - rrt(198) * density(01) 
  pd(40,01) = pd(40,01) + rrt(198) * density(35) 
  pd(40,35) = pd(40,35) + rrt(198) * density(01) 
  pd(01,01) = pd(01,01) + rrt(199) * density(09) 
  pd(01,09) = pd(01,09) + rrt(199) * density(01) 
  pd(09,01) = pd(09,01) - rrt(199) * density(09) 
  pd(09,09) = pd(09,09) - rrt(199) * density(01) 
  pd(34,01) = pd(34,01) + rrt(199) * density(09) * 2.0d0
  pd(34,09) = pd(34,09) + rrt(199) * density(01) * 2.0d0
  pd(40,01) = pd(40,01) + rrt(199) * density(09) 
  pd(40,09) = pd(40,09) + rrt(199) * density(01) 
  pd(01,01) = pd(01,01) + rrt(200) * density(35) 
  pd(01,35) = pd(01,35) + rrt(200) * density(01) 
  pd(16,01) = pd(16,01) + rrt(200) * density(35) 
  pd(16,35) = pd(16,35) + rrt(200) * density(01) 
  pd(17,01) = pd(17,01) + rrt(200) * density(35) 
  pd(17,35) = pd(17,35) + rrt(200) * density(01) 
  pd(35,01) = pd(35,01) - rrt(200) * density(35) 
  pd(35,35) = pd(35,35) - rrt(200) * density(01) 
  pd(01,01) = pd(01,01) + rrt(201) * density(09) 
  pd(01,09) = pd(01,09) + rrt(201) * density(01) 
  pd(09,01) = pd(09,01) - rrt(201) * density(09) 
  pd(09,09) = pd(09,09) - rrt(201) * density(01) 
  pd(16,01) = pd(16,01) + rrt(201) * density(09) 
  pd(16,09) = pd(16,09) + rrt(201) * density(01) 
  pd(17,01) = pd(17,01) + rrt(201) * density(09) 
  pd(17,09) = pd(17,09) + rrt(201) * density(01) 
  pd(01,01) = pd(01,01) + rrt(202) * density(35) 
  pd(01,35) = pd(01,35) + rrt(202) * density(01) 
  pd(03,01) = pd(03,01) + rrt(202) * density(35) 
  pd(03,35) = pd(03,35) + rrt(202) * density(01) 
  pd(35,01) = pd(35,01) - rrt(202) * density(35) 
  pd(35,35) = pd(35,35) - rrt(202) * density(01) 
  pd(44,01) = pd(44,01) + rrt(202) * density(35) 
  pd(44,35) = pd(44,35) + rrt(202) * density(01) 
  pd(01,01) = pd(01,01) + rrt(203) * density(09) 
  pd(01,09) = pd(01,09) + rrt(203) * density(01) 
  pd(03,01) = pd(03,01) + rrt(203) * density(09) 
  pd(03,09) = pd(03,09) + rrt(203) * density(01) 
  pd(09,01) = pd(09,01) - rrt(203) * density(09) 
  pd(09,09) = pd(09,09) - rrt(203) * density(01) 
  pd(44,01) = pd(44,01) + rrt(203) * density(09) 
  pd(44,09) = pd(44,09) + rrt(203) * density(01) 
  pd(01,01) = pd(01,01) + rrt(204) * density(15) 
  pd(01,15) = pd(01,15) + rrt(204) * density(01) 
  pd(15,01) = pd(15,01) - rrt(204) * density(15) 
  pd(15,15) = pd(15,15) - rrt(204) * density(01) 
  pd(33,01) = pd(33,01) + rrt(204) * density(15) 
  pd(33,15) = pd(33,15) + rrt(204) * density(01) 
  pd(01,01) = pd(01,01) + rrt(205) * density(02) 
  pd(01,02) = pd(01,02) + rrt(205) * density(01) 
  pd(02,01) = pd(02,01) - rrt(205) * density(02) 
  pd(02,02) = pd(02,02) - rrt(205) * density(01) 
  pd(33,01) = pd(33,01) + rrt(205) * density(02) 
  pd(33,02) = pd(33,02) + rrt(205) * density(01) 
  pd(01,01) = pd(01,01) + rrt(206) * density(15) 
  pd(01,15) = pd(01,15) + rrt(206) * density(01) 
  pd(10,01) = pd(10,01) + rrt(206) * density(15) 
  pd(10,15) = pd(10,15) + rrt(206) * density(01) 
  pd(15,01) = pd(15,01) - rrt(206) * density(15) 
  pd(15,15) = pd(15,15) - rrt(206) * density(01) 
  pd(38,01) = pd(38,01) + rrt(206) * density(15) 
  pd(38,15) = pd(38,15) + rrt(206) * density(01) 
  pd(01,01) = pd(01,01) + rrt(207) * density(02) 
  pd(01,02) = pd(01,02) + rrt(207) * density(01) 
  pd(02,01) = pd(02,01) - rrt(207) * density(02) 
  pd(02,02) = pd(02,02) - rrt(207) * density(01) 
  pd(10,01) = pd(10,01) + rrt(207) * density(02) 
  pd(10,02) = pd(10,02) + rrt(207) * density(01) 
  pd(38,01) = pd(38,01) + rrt(207) * density(02) 
  pd(38,02) = pd(38,02) + rrt(207) * density(01) 
  pd(01,01) = pd(01,01) + rrt(208) * density(15) 
  pd(01,15) = pd(01,15) + rrt(208) * density(01) 
  pd(15,01) = pd(15,01) - rrt(208) * density(15) 
  pd(15,15) = pd(15,15) - rrt(208) * density(01) 
  pd(17,01) = pd(17,01) + rrt(208) * density(15) 
  pd(17,15) = pd(17,15) + rrt(208) * density(01) 
  pd(29,01) = pd(29,01) + rrt(208) * density(15) 
  pd(29,15) = pd(29,15) + rrt(208) * density(01) 
  pd(01,01) = pd(01,01) + rrt(209) * density(02) 
  pd(01,02) = pd(01,02) + rrt(209) * density(01) 
  pd(02,01) = pd(02,01) - rrt(209) * density(02) 
  pd(02,02) = pd(02,02) - rrt(209) * density(01) 
  pd(17,01) = pd(17,01) + rrt(209) * density(02) 
  pd(17,02) = pd(17,02) + rrt(209) * density(01) 
  pd(29,01) = pd(29,01) + rrt(209) * density(02) 
  pd(29,02) = pd(29,02) + rrt(209) * density(01) 
  pd(01,01) = pd(01,01) + rrt(210) * density(15) 
  pd(01,15) = pd(01,15) + rrt(210) * density(01) 
  pd(03,01) = pd(03,01) + rrt(210) * density(15) 
  pd(03,15) = pd(03,15) + rrt(210) * density(01) 
  pd(08,01) = pd(08,01) + rrt(210) * density(15) 
  pd(08,15) = pd(08,15) + rrt(210) * density(01) 
  pd(15,01) = pd(15,01) - rrt(210) * density(15) 
  pd(15,15) = pd(15,15) - rrt(210) * density(01) 
  pd(01,01) = pd(01,01) + rrt(211) * density(02) 
  pd(01,02) = pd(01,02) + rrt(211) * density(01) 
  pd(02,01) = pd(02,01) - rrt(211) * density(02) 
  pd(02,02) = pd(02,02) - rrt(211) * density(01) 
  pd(03,01) = pd(03,01) + rrt(211) * density(02) 
  pd(03,02) = pd(03,02) + rrt(211) * density(01) 
  pd(08,01) = pd(08,01) + rrt(211) * density(02) 
  pd(08,02) = pd(08,02) + rrt(211) * density(01) 
  pd(01,01) = pd(01,01) + rrt(212) * density(15) 
  pd(01,15) = pd(01,15) + rrt(212) * density(01) 
  pd(15,01) = pd(15,01) - rrt(212) * density(15) 
  pd(15,15) = pd(15,15) - rrt(212) * density(01) 
  pd(16,01) = pd(16,01) + rrt(212) * density(15) 
  pd(16,15) = pd(16,15) + rrt(212) * density(01) 
  pd(42,01) = pd(42,01) + rrt(212) * density(15) 
  pd(42,15) = pd(42,15) + rrt(212) * density(01) 
  pd(01,01) = pd(01,01) + rrt(213) * density(02) 
  pd(01,02) = pd(01,02) + rrt(213) * density(01) 
  pd(02,01) = pd(02,01) - rrt(213) * density(02) 
  pd(02,02) = pd(02,02) - rrt(213) * density(01) 
  pd(16,01) = pd(16,01) + rrt(213) * density(02) 
  pd(16,02) = pd(16,02) + rrt(213) * density(01) 
  pd(42,01) = pd(42,01) + rrt(213) * density(02) 
  pd(42,02) = pd(42,02) + rrt(213) * density(01) 
  pd(01,01) = pd(01,01) + rrt(214) * density(26) 
  pd(01,26) = pd(01,26) + rrt(214) * density(01) 
  pd(26,01) = pd(26,01) - rrt(214) * density(26) 
  pd(26,26) = pd(26,26) - rrt(214) * density(01) 
  pd(40,01) = pd(40,01) + rrt(214) * density(26) 
  pd(40,26) = pd(40,26) + rrt(214) * density(01) 
  pd(01,01) = pd(01,01) + rrt(215) * density(36) 
  pd(01,36) = pd(01,36) + rrt(215) * density(01) 
  pd(36,01) = pd(36,01) - rrt(215) * density(36) 
  pd(36,36) = pd(36,36) - rrt(215) * density(01) 
  pd(40,01) = pd(40,01) + rrt(215) * density(36) 
  pd(40,36) = pd(40,36) + rrt(215) * density(01) 
  pd(01,01) = pd(01,01) + rrt(216) * density(28) 
  pd(01,28) = pd(01,28) + rrt(216) * density(01) 
  pd(28,01) = pd(28,01) - rrt(216) * density(28) 
  pd(28,28) = pd(28,28) - rrt(216) * density(01) 
  pd(40,01) = pd(40,01) + rrt(216) * density(28) 
  pd(40,28) = pd(40,28) + rrt(216) * density(01) 
  pd(01,01) = pd(01,01) + rrt(217) * density(26) 
  pd(01,26) = pd(01,26) + rrt(217) * density(01) 
  pd(26,01) = pd(26,01) - rrt(217) * density(26) 
  pd(26,26) = pd(26,26) - rrt(217) * density(01) 
  pd(29,01) = pd(29,01) + rrt(217) * density(26) 
  pd(29,26) = pd(29,26) + rrt(217) * density(01) 
  pd(42,01) = pd(42,01) + rrt(217) * density(26) 
  pd(42,26) = pd(42,26) + rrt(217) * density(01) 
  pd(01,01) = pd(01,01) + rrt(218) * density(36) 
  pd(01,36) = pd(01,36) + rrt(218) * density(01) 
  pd(29,01) = pd(29,01) + rrt(218) * density(36) 
  pd(29,36) = pd(29,36) + rrt(218) * density(01) 
  pd(36,01) = pd(36,01) - rrt(218) * density(36) 
  pd(36,36) = pd(36,36) - rrt(218) * density(01) 
  pd(42,01) = pd(42,01) + rrt(218) * density(36) 
  pd(42,36) = pd(42,36) + rrt(218) * density(01) 
  pd(01,01) = pd(01,01) + rrt(219) * density(28) 
  pd(01,28) = pd(01,28) + rrt(219) * density(01) 
  pd(28,01) = pd(28,01) - rrt(219) * density(28) 
  pd(28,28) = pd(28,28) - rrt(219) * density(01) 
  pd(29,01) = pd(29,01) + rrt(219) * density(28) 
  pd(29,28) = pd(29,28) + rrt(219) * density(01) 
  pd(42,01) = pd(42,01) + rrt(219) * density(28) 
  pd(42,28) = pd(42,28) + rrt(219) * density(01) 
  pd(01,01) = pd(01,01) + rrt(220) * density(04) 
  pd(01,04) = pd(01,04) + rrt(220) * density(01) 
  pd(04,01) = pd(04,01) - rrt(220) * density(04) 
  pd(04,04) = pd(04,04) - rrt(220) * density(01) 
  pd(25,01) = pd(25,01) + rrt(220) * density(04) 
  pd(25,04) = pd(25,04) + rrt(220) * density(01) 
  pd(01,01) = pd(01,01) + rrt(221) * density(24) 
  pd(01,24) = pd(01,24) + rrt(221) * density(01) 
  pd(24,01) = pd(24,01) - rrt(221) * density(24) 
  pd(24,24) = pd(24,24) - rrt(221) * density(01) 
  pd(25,01) = pd(25,01) + rrt(221) * density(24) 
  pd(25,24) = pd(25,24) + rrt(221) * density(01) 
  pd(01,01) = pd(01,01) + rrt(222) * density(04) 
  pd(01,04) = pd(01,04) + rrt(222) * density(01) 
  pd(04,01) = pd(04,01) - rrt(222) * density(04) 
  pd(04,04) = pd(04,04) - rrt(222) * density(01) 
  pd(10,01) = pd(10,01) + rrt(222) * density(04) 
  pd(10,04) = pd(10,04) + rrt(222) * density(01) 
  pd(23,01) = pd(23,01) + rrt(222) * density(04) 
  pd(23,04) = pd(23,04) + rrt(222) * density(01) 
  pd(01,01) = pd(01,01) + rrt(223) * density(24) 
  pd(01,24) = pd(01,24) + rrt(223) * density(01) 
  pd(10,01) = pd(10,01) + rrt(223) * density(24) 
  pd(10,24) = pd(10,24) + rrt(223) * density(01) 
  pd(23,01) = pd(23,01) + rrt(223) * density(24) 
  pd(23,24) = pd(23,24) + rrt(223) * density(01) 
  pd(24,01) = pd(24,01) - rrt(223) * density(24) 
  pd(24,24) = pd(24,24) - rrt(223) * density(01) 
  pd(01,01) = pd(01,01) + rrt(224) * density(04) 
  pd(01,04) = pd(01,04) + rrt(224) * density(01) 
  pd(04,01) = pd(04,01) - rrt(224) * density(04) 
  pd(04,04) = pd(04,04) - rrt(224) * density(01) 
  pd(19,01) = pd(19,01) + rrt(224) * density(04) 
  pd(19,04) = pd(19,04) + rrt(224) * density(01) 
  pd(34,01) = pd(34,01) + rrt(224) * density(04) 
  pd(34,04) = pd(34,04) + rrt(224) * density(01) 
  pd(01,01) = pd(01,01) + rrt(225) * density(24) 
  pd(01,24) = pd(01,24) + rrt(225) * density(01) 
  pd(19,01) = pd(19,01) + rrt(225) * density(24) 
  pd(19,24) = pd(19,24) + rrt(225) * density(01) 
  pd(24,01) = pd(24,01) - rrt(225) * density(24) 
  pd(24,24) = pd(24,24) - rrt(225) * density(01) 
  pd(34,01) = pd(34,01) + rrt(225) * density(24) 
  pd(34,24) = pd(34,24) + rrt(225) * density(01) 
  pd(01,01) = pd(01,01) + rrt(226) * density(04) 
  pd(01,04) = pd(01,04) + rrt(226) * density(01) 
  pd(04,01) = pd(04,01) - rrt(226) * density(04) 
  pd(04,04) = pd(04,04) - rrt(226) * density(01) 
  pd(10,01) = pd(10,01) + rrt(226) * density(04) 
  pd(10,04) = pd(10,04) + rrt(226) * density(01) 
  pd(22,01) = pd(22,01) + rrt(226) * density(04) 
  pd(22,04) = pd(22,04) + rrt(226) * density(01) 
  pd(34,01) = pd(34,01) + rrt(226) * density(04) 
  pd(34,04) = pd(34,04) + rrt(226) * density(01) 
  pd(01,01) = pd(01,01) + rrt(227) * density(24) 
  pd(01,24) = pd(01,24) + rrt(227) * density(01) 
  pd(10,01) = pd(10,01) + rrt(227) * density(24) 
  pd(10,24) = pd(10,24) + rrt(227) * density(01) 
  pd(22,01) = pd(22,01) + rrt(227) * density(24) 
  pd(22,24) = pd(22,24) + rrt(227) * density(01) 
  pd(24,01) = pd(24,01) - rrt(227) * density(24) 
  pd(24,24) = pd(24,24) - rrt(227) * density(01) 
  pd(34,01) = pd(34,01) + rrt(227) * density(24) 
  pd(34,24) = pd(34,24) + rrt(227) * density(01) 
  pd(01,01) = pd(01,01) + rrt(228) * density(04) 
  pd(01,04) = pd(01,04) + rrt(228) * density(01) 
  pd(04,01) = pd(04,01) - rrt(228) * density(04) 
  pd(04,04) = pd(04,04) - rrt(228) * density(01) 
  pd(11,01) = pd(11,01) + rrt(228) * density(04) 
  pd(11,04) = pd(11,04) + rrt(228) * density(01) 
  pd(34,01) = pd(34,01) + rrt(228) * density(04) * 2.0d0
  pd(34,04) = pd(34,04) + rrt(228) * density(01) * 2.0d0
  pd(01,01) = pd(01,01) + rrt(229) * density(24) 
  pd(01,24) = pd(01,24) + rrt(229) * density(01) 
  pd(11,01) = pd(11,01) + rrt(229) * density(24) 
  pd(11,24) = pd(11,24) + rrt(229) * density(01) 
  pd(24,01) = pd(24,01) - rrt(229) * density(24) 
  pd(24,24) = pd(24,24) - rrt(229) * density(01) 
  pd(34,01) = pd(34,01) + rrt(229) * density(24) * 2.0d0
  pd(34,24) = pd(34,24) + rrt(229) * density(01) * 2.0d0
  pd(01,01) = pd(01,01) + rrt(230) * density(04) 
  pd(01,04) = pd(01,04) + rrt(230) * density(01) 
  pd(04,01) = pd(04,01) - rrt(230) * density(04) 
  pd(04,04) = pd(04,04) - rrt(230) * density(01) 
  pd(16,01) = pd(16,01) + rrt(230) * density(04) 
  pd(16,04) = pd(16,04) + rrt(230) * density(01) 
  pd(41,01) = pd(41,01) + rrt(230) * density(04) 
  pd(41,04) = pd(41,04) + rrt(230) * density(01) 
  pd(01,01) = pd(01,01) + rrt(231) * density(24) 
  pd(01,24) = pd(01,24) + rrt(231) * density(01) 
  pd(16,01) = pd(16,01) + rrt(231) * density(24) 
  pd(16,24) = pd(16,24) + rrt(231) * density(01) 
  pd(24,01) = pd(24,01) - rrt(231) * density(24) 
  pd(24,24) = pd(24,24) - rrt(231) * density(01) 
  pd(41,01) = pd(41,01) + rrt(231) * density(24) 
  pd(41,24) = pd(41,24) + rrt(231) * density(01) 
  pd(01,01) = pd(01,01) + rrt(232) * density(04) 
  pd(01,04) = pd(01,04) + rrt(232) * density(01) 
  pd(04,01) = pd(04,01) - rrt(232) * density(04) 
  pd(04,04) = pd(04,04) - rrt(232) * density(01) 
  pd(33,01) = pd(33,01) + rrt(232) * density(04) 
  pd(33,04) = pd(33,04) + rrt(232) * density(01) 
  pd(44,01) = pd(44,01) + rrt(232) * density(04) 
  pd(44,04) = pd(44,04) + rrt(232) * density(01) 
  pd(01,01) = pd(01,01) + rrt(233) * density(24) 
  pd(01,24) = pd(01,24) + rrt(233) * density(01) 
  pd(24,01) = pd(24,01) - rrt(233) * density(24) 
  pd(24,24) = pd(24,24) - rrt(233) * density(01) 
  pd(33,01) = pd(33,01) + rrt(233) * density(24) 
  pd(33,24) = pd(33,24) + rrt(233) * density(01) 
  pd(44,01) = pd(44,01) + rrt(233) * density(24) 
  pd(44,24) = pd(44,24) + rrt(233) * density(01) 
  pd(01,01) = pd(01,01) + rrt(234) * density(04) 
  pd(01,04) = pd(01,04) + rrt(234) * density(01) 
  pd(04,01) = pd(04,01) - rrt(234) * density(04) 
  pd(04,04) = pd(04,04) - rrt(234) * density(01) 
  pd(17,01) = pd(17,01) + rrt(234) * density(04) 
  pd(17,04) = pd(17,04) + rrt(234) * density(01) 
  pd(30,01) = pd(30,01) + rrt(234) * density(04) 
  pd(30,04) = pd(30,04) + rrt(234) * density(01) 
  pd(01,01) = pd(01,01) + rrt(235) * density(24) 
  pd(01,24) = pd(01,24) + rrt(235) * density(01) 
  pd(17,01) = pd(17,01) + rrt(235) * density(24) 
  pd(17,24) = pd(17,24) + rrt(235) * density(01) 
  pd(24,01) = pd(24,01) - rrt(235) * density(24) 
  pd(24,24) = pd(24,24) - rrt(235) * density(01) 
  pd(30,01) = pd(30,01) + rrt(235) * density(24) 
  pd(30,24) = pd(30,24) + rrt(235) * density(01) 
  pd(01,01) = pd(01,01) + rrt(236) * density(04) 
  pd(01,04) = pd(01,04) + rrt(236) * density(01) 
  pd(03,01) = pd(03,01) + rrt(236) * density(04) 
  pd(03,04) = pd(03,04) + rrt(236) * density(01) 
  pd(04,01) = pd(04,01) - rrt(236) * density(04) 
  pd(04,04) = pd(04,04) - rrt(236) * density(01) 
  pd(18,01) = pd(18,01) + rrt(236) * density(04) 
  pd(18,04) = pd(18,04) + rrt(236) * density(01) 
  pd(01,01) = pd(01,01) + rrt(237) * density(24) 
  pd(01,24) = pd(01,24) + rrt(237) * density(01) 
  pd(03,01) = pd(03,01) + rrt(237) * density(24) 
  pd(03,24) = pd(03,24) + rrt(237) * density(01) 
  pd(18,01) = pd(18,01) + rrt(237) * density(24) 
  pd(18,24) = pd(18,24) + rrt(237) * density(01) 
  pd(24,01) = pd(24,01) - rrt(237) * density(24) 
  pd(24,24) = pd(24,24) - rrt(237) * density(01) 
  pd(10,01) = pd(10,01) + rrt(238) * density(34) * 2.0d0
  pd(10,34) = pd(10,34) + rrt(238) * density(01) * 2.0d0
  pd(34,01) = pd(34,01) - rrt(238) * density(34) 
  pd(34,34) = pd(34,34) - rrt(238) * density(01) 
  pd(01,01) = pd(01,01) + rrt(239) * density(34) 
  pd(01,34) = pd(01,34) + rrt(239) * density(01) 
  pd(34,01) = pd(34,01) - rrt(239) * density(34) 
  pd(34,34) = pd(34,34) - rrt(239) * density(01) 
  pd(43,01) = pd(43,01) + rrt(239) * density(34) 
  pd(43,34) = pd(43,34) + rrt(239) * density(01) 
  if( ldensity_constant ) then
    do i = 1, species_max
      if( density_constant(i) ) pd(i,:) = 0.0d0
    enddo
  endif
  if( lgas_heating ) then
    pd(46,1) = eV_to_K * ZDPlasKin_cfg(11)
    pd(46,:) = pd(46,:) * ZDPlasKin_cfg(13)
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
  call ZDPlasKin_bolsig_rates()
  rrt(001) = bolsig_rates(bolsig_pointer(1))
  rrt(002) = bolsig_rates(bolsig_pointer(2))
  rrt(003) = bolsig_rates(bolsig_pointer(3))
  rrt(004) = bolsig_rates(bolsig_pointer(4))
  rrt(005) = bolsig_rates(bolsig_pointer(5))
  rrt(006) = bolsig_rates(bolsig_pointer(6))
  rrt(007) = bolsig_rates(bolsig_pointer(7))
  rrt(008) = bolsig_rates(bolsig_pointer(8))
  rrt(009) = bolsig_rates(bolsig_pointer(9))
  rrt(010) = bolsig_rates(bolsig_pointer(10))
  rrt(011) = bolsig_rates(bolsig_pointer(11))
  rrt(012) = bolsig_rates(bolsig_pointer(12))
  rrt(013) = bolsig_rates(bolsig_pointer(13))
  rrt(014) = bolsig_rates(bolsig_pointer(14))
  rrt(015) = bolsig_rates(bolsig_pointer(15))
  rrt(016) = bolsig_rates(bolsig_pointer(16))
  rrt(017) = bolsig_rates(bolsig_pointer(17))
  rrt(018) = bolsig_rates(bolsig_pointer(18))
  rrt(019) = bolsig_rates(bolsig_pointer(19))
  rrt(020) = bolsig_rates(bolsig_pointer(20))
  rrt(021) = bolsig_rates(bolsig_pointer(21))
  rrt(022) = bolsig_rates(bolsig_pointer(22))
  rrt(023) = bolsig_rates(bolsig_pointer(23))
  rrt(024) = bolsig_rates(bolsig_pointer(24))
  rrt(025) = bolsig_rates(bolsig_pointer(25))
  rrt(026) = bolsig_rates(bolsig_pointer(26))
  rrt(027) = bolsig_rates(bolsig_pointer(27))
  rrt(028) = bolsig_rates(bolsig_pointer(28))
  rrt(029) = bolsig_rates(bolsig_pointer(29))
  rrt(030) = bolsig_rates(bolsig_pointer(30))
  rrt(031) = bolsig_rates(bolsig_pointer(31))
  rrt(032) = bolsig_rates(bolsig_pointer(32))
  rrt(033) = bolsig_rates(bolsig_pointer(33))
  rrt(034) = bolsig_rates(bolsig_pointer(34))
  rrt(035) = bolsig_rates(bolsig_pointer(35))
  rrt(036) = bolsig_rates(bolsig_pointer(36))
  rrt(037) = bolsig_rates(bolsig_pointer(37))
  rrt(038) = bolsig_rates(bolsig_pointer(38))
  rrt(039) = bolsig_rates(bolsig_pointer(39))
  rrt(040) = bolsig_rates(bolsig_pointer(40))
  rrt(041) = bolsig_rates(bolsig_pointer(41))
  rrt(042) = bolsig_rates(bolsig_pointer(42))
  rrt(043) = bolsig_rates(bolsig_pointer(43))
  rrt(044) = bolsig_rates(bolsig_pointer(44))
  rrt(045) = bolsig_rates(bolsig_pointer(45))
  rrt(046) = bolsig_rates(bolsig_pointer(46))
  rrt(047) = bolsig_rates(bolsig_pointer(47))
  rrt(048) = bolsig_rates(bolsig_pointer(48))
  rrt(049) = bolsig_rates(bolsig_pointer(49))
  rrt(050) = bolsig_rates(bolsig_pointer(50))
  rrt(051) = bolsig_rates(bolsig_pointer(51))
  rrt(052) = bolsig_rates(bolsig_pointer(52))
  rrt(053) = bolsig_rates(bolsig_pointer(53))
  rrt(054) = bolsig_rates(bolsig_pointer(54))
  rrt(055) = bolsig_rates(bolsig_pointer(55))
  rrt(056) = bolsig_rates(bolsig_pointer(56))
  rrt(057) = bolsig_rates(bolsig_pointer(57))
  rrt(058) = bolsig_rates(bolsig_pointer(58))
  rrt(059) = bolsig_rates(bolsig_pointer(59))
  rrt(060) = bolsig_rates(bolsig_pointer(60))
  rrt(061) = bolsig_rates(bolsig_pointer(61))
  rrt(062) = bolsig_rates(bolsig_pointer(62))
  rrt(063) = bolsig_rates(bolsig_pointer(63))
  rrt(064) = bolsig_rates(bolsig_pointer(64))
  rrt(065) = bolsig_rates(bolsig_pointer(65))
  rrt(066) = bolsig_rates(bolsig_pointer(66))
  rrt(067) = bolsig_rates(bolsig_pointer(67))
  rrt(068) = bolsig_rates(bolsig_pointer(68))
  rrt(069) = bolsig_rates(bolsig_pointer(69))
  rrt(070) = bolsig_rates(bolsig_pointer(70))
  rrt(071) = bolsig_rates(bolsig_pointer(71))
  rrt(072) = bolsig_rates(bolsig_pointer(72))
  rrt(073) = bolsig_rates(bolsig_pointer(73))
  rrt(074) = bolsig_rates(bolsig_pointer(74))
  rrt(075) = bolsig_rates(bolsig_pointer(75))
  rrt(076) = bolsig_rates(bolsig_pointer(76))
  rrt(077) = bolsig_rates(bolsig_pointer(77))
  rrt(078) = bolsig_rates(bolsig_pointer(78))
  rrt(079) = bolsig_rates(bolsig_pointer(79))
  rrt(080) = bolsig_rates(bolsig_pointer(80))
  rrt(081) = bolsig_rates(bolsig_pointer(81))
  rrt(082) = bolsig_rates(bolsig_pointer(82))
  rrt(083) = bolsig_rates(bolsig_pointer(83))
  rrt(084) = bolsig_rates(bolsig_pointer(84))
  rrt(085) = bolsig_rates(bolsig_pointer(85))
  rrt(086) = bolsig_rates(bolsig_pointer(86))
  rrt(087) = bolsig_rates(bolsig_pointer(87))
  rrt(088) = bolsig_rates(bolsig_pointer(88))
  rrt(089) = bolsig_rates(bolsig_pointer(89))
  rrt(090) = bolsig_rates(bolsig_pointer(90))
  rrt(091) = bolsig_rates(bolsig_pointer(91))
  rrt(092) = bolsig_rates(bolsig_pointer(92))
  rrt(093) = bolsig_rates(bolsig_pointer(93))
  rrt(094) = bolsig_rates(bolsig_pointer(94))
  rrt(095) = bolsig_rates(bolsig_pointer(95))
  rrt(096) = bolsig_rates(bolsig_pointer(96))
  rrt(097) = bolsig_rates(bolsig_pointer(97))
  rrt(098) = bolsig_rates(bolsig_pointer(98))
  rrt(099) = bolsig_rates(bolsig_pointer(99))
  rrt(100) = bolsig_rates(bolsig_pointer(100))
  rrt(101) = bolsig_rates(bolsig_pointer(101))
  rrt(102) = bolsig_rates(bolsig_pointer(102))
  rrt(103) = bolsig_rates(bolsig_pointer(103))
  rrt(104) = bolsig_rates(bolsig_pointer(104))
  rrt(105) = bolsig_rates(bolsig_pointer(105))
  rrt(106) = bolsig_rates(bolsig_pointer(106))
  rrt(107) = bolsig_rates(bolsig_pointer(107))
  rrt(108) = bolsig_rates(bolsig_pointer(108))
  rrt(109) = bolsig_rates(bolsig_pointer(109))
  rrt(110) = bolsig_rates(bolsig_pointer(110))
  rrt(111) = bolsig_rates(bolsig_pointer(111))
  rrt(112) = bolsig_rates(bolsig_pointer(112))
  rrt(113) = bolsig_rates(bolsig_pointer(113))
  rrt(114) = bolsig_rates(bolsig_pointer(114))
  rrt(115) = bolsig_rates(bolsig_pointer(115))
  rrt(116) = bolsig_rates(bolsig_pointer(116))
  rrt(117) = bolsig_rates(bolsig_pointer(117))
  rrt(118) = bolsig_rates(bolsig_pointer(118))
  rrt(119) = bolsig_rates(bolsig_pointer(119))
  rrt(120) = bolsig_rates(bolsig_pointer(120))
  rrt(121) = bolsig_rates(bolsig_pointer(121))
  rrt(122) = bolsig_rates(bolsig_pointer(122))
  rrt(123) = bolsig_rates(bolsig_pointer(123))
  rrt(124) = bolsig_rates(bolsig_pointer(124))
  rrt(125) = bolsig_rates(bolsig_pointer(125))
  rrt(126) = bolsig_rates(bolsig_pointer(126))
  rrt(127) = bolsig_rates(bolsig_pointer(127))
  rrt(128) = bolsig_rates(bolsig_pointer(128))
  rrt(129) = bolsig_rates(bolsig_pointer(129))
  rrt(130) = bolsig_rates(bolsig_pointer(130))
  rrt(131) = bolsig_rates(bolsig_pointer(131))
  rrt(132) = bolsig_rates(bolsig_pointer(132))
  rrt(133) = bolsig_rates(bolsig_pointer(133))
  rrt(134) = bolsig_rates(bolsig_pointer(134))
  rrt(135) = bolsig_rates(bolsig_pointer(135))
  rrt(136) = bolsig_rates(bolsig_pointer(136))
  rrt(137) = bolsig_rates(bolsig_pointer(137))
  rrt(138) = bolsig_rates(bolsig_pointer(138))
  rrt(139) = bolsig_rates(bolsig_pointer(139))
  rrt(140) = bolsig_rates(bolsig_pointer(140))
  rrt(141) = bolsig_rates(bolsig_pointer(141))
  rrt(142) = bolsig_rates(bolsig_pointer(142))
  rrt(143) = bolsig_rates(bolsig_pointer(143))
  rrt(144) = bolsig_rates(bolsig_pointer(144))
  rrt(145) = bolsig_rates(bolsig_pointer(145))
  rrt(146) = bolsig_rates(bolsig_pointer(146))
  rrt(147) = bolsig_rates(bolsig_pointer(147))
  rrt(148) = bolsig_rates(bolsig_pointer(148))
  rrt(149) = bolsig_rates(bolsig_pointer(149))
  rrt(150) = bolsig_rates(bolsig_pointer(150))
  rrt(151) = bolsig_rates(bolsig_pointer(151))
  rrt(152) = bolsig_rates(bolsig_pointer(152))
  rrt(153) = bolsig_rates(bolsig_pointer(153))
  rrt(154) = bolsig_rates(bolsig_pointer(154))
  rrt(155) = bolsig_rates(bolsig_pointer(155))
  rrt(156) = bolsig_rates(bolsig_pointer(156))
  rrt(157) = bolsig_rates(bolsig_pointer(157))
  rrt(158) = bolsig_rates(bolsig_pointer(158))
  rrt(159) = bolsig_rates(bolsig_pointer(159))
  rrt(160) = bolsig_rates(bolsig_pointer(160))
  rrt(161) = bolsig_rates(bolsig_pointer(161))
  rrt(162) = bolsig_rates(bolsig_pointer(162))
  rrt(163) = bolsig_rates(bolsig_pointer(163))
  rrt(164) = bolsig_rates(bolsig_pointer(164))
  rrt(165) = bolsig_rates(bolsig_pointer(165))
  rrt(166) = bolsig_rates(bolsig_pointer(166))
  rrt(167) = bolsig_rates(bolsig_pointer(167))
  rrt(168) = bolsig_rates(bolsig_pointer(168))
  rrt(169) = bolsig_rates(bolsig_pointer(169))
  rrt(170) = bolsig_rates(bolsig_pointer(170))
  rrt(171) = bolsig_rates(bolsig_pointer(171))
  rrt(172) = bolsig_rates(bolsig_pointer(172))
  rrt(173) = bolsig_rates(bolsig_pointer(173))
  rrt(174) = bolsig_rates(bolsig_pointer(174))
  rrt(175) = bolsig_rates(bolsig_pointer(175))
  rrt(176) = bolsig_rates(bolsig_pointer(176))
  rrt(177) = bolsig_rates(bolsig_pointer(177))
  rrt(178) = bolsig_rates(bolsig_pointer(178))
  rrt(179) = bolsig_rates(bolsig_pointer(179))
  rrt(180) = bolsig_rates(bolsig_pointer(180))
  rrt(181) = bolsig_rates(bolsig_pointer(181))
  rrt(182) = bolsig_rates(bolsig_pointer(182))
  rrt(183) = bolsig_rates(bolsig_pointer(183))
  rrt(184) = bolsig_rates(bolsig_pointer(184))
  rrt(185) = bolsig_rates(bolsig_pointer(185))
  rrt(186) = bolsig_rates(bolsig_pointer(186))
  rrt(187) = bolsig_rates(bolsig_pointer(187))
  rrt(188) = bolsig_rates(bolsig_pointer(188))
  rrt(189) = bolsig_rates(bolsig_pointer(189))
  rrt(190) = bolsig_rates(bolsig_pointer(190))
  rrt(191) = bolsig_rates(bolsig_pointer(191))
  rrt(192) = bolsig_rates(bolsig_pointer(192))
  rrt(193) = bolsig_rates(bolsig_pointer(193))
  rrt(194) = bolsig_rates(bolsig_pointer(194))
  rrt(195) = bolsig_rates(bolsig_pointer(195))
  rrt(196) = bolsig_rates(bolsig_pointer(196))
  rrt(197) = bolsig_rates(bolsig_pointer(197))
  rrt(198) = bolsig_rates(bolsig_pointer(198))
  rrt(199) = bolsig_rates(bolsig_pointer(199))
  rrt(200) = bolsig_rates(bolsig_pointer(200))
  rrt(201) = bolsig_rates(bolsig_pointer(201))
  rrt(202) = bolsig_rates(bolsig_pointer(202))
  rrt(203) = bolsig_rates(bolsig_pointer(203))
  rrt(204) = bolsig_rates(bolsig_pointer(204))
  rrt(205) = bolsig_rates(bolsig_pointer(205))
  rrt(206) = bolsig_rates(bolsig_pointer(206))
  rrt(207) = bolsig_rates(bolsig_pointer(207))
  rrt(208) = bolsig_rates(bolsig_pointer(208))
  rrt(209) = bolsig_rates(bolsig_pointer(209))
  rrt(210) = bolsig_rates(bolsig_pointer(210))
  rrt(211) = bolsig_rates(bolsig_pointer(211))
  rrt(212) = bolsig_rates(bolsig_pointer(212))
  rrt(213) = bolsig_rates(bolsig_pointer(213))
  rrt(214) = bolsig_rates(bolsig_pointer(214))
  rrt(215) = bolsig_rates(bolsig_pointer(215))
  rrt(216) = bolsig_rates(bolsig_pointer(216))
  rrt(217) = bolsig_rates(bolsig_pointer(217))
  rrt(218) = bolsig_rates(bolsig_pointer(218))
  rrt(219) = bolsig_rates(bolsig_pointer(219))
  rrt(220) = bolsig_rates(bolsig_pointer(220))
  rrt(221) = bolsig_rates(bolsig_pointer(221))
  rrt(222) = bolsig_rates(bolsig_pointer(222))
  rrt(223) = bolsig_rates(bolsig_pointer(223))
  rrt(224) = bolsig_rates(bolsig_pointer(224))
  rrt(225) = bolsig_rates(bolsig_pointer(225))
  rrt(226) = bolsig_rates(bolsig_pointer(226))
  rrt(227) = bolsig_rates(bolsig_pointer(227))
  rrt(228) = bolsig_rates(bolsig_pointer(228))
  rrt(229) = bolsig_rates(bolsig_pointer(229))
  rrt(230) = bolsig_rates(bolsig_pointer(230))
  rrt(231) = bolsig_rates(bolsig_pointer(231))
  rrt(232) = bolsig_rates(bolsig_pointer(232))
  rrt(233) = bolsig_rates(bolsig_pointer(233))
  rrt(234) = bolsig_rates(bolsig_pointer(234))
  rrt(235) = bolsig_rates(bolsig_pointer(235))
  rrt(236) = bolsig_rates(bolsig_pointer(236))
  rrt(237) = bolsig_rates(bolsig_pointer(237))
  rrt(238) = bolsig_rates(bolsig_pointer(238))
  rrt(239) = bolsig_rates(bolsig_pointer(239))
  where( lreaction_block(:) ) rrt(:) = 0.0d0
  return
end subroutine ZDPlasKin_reac_rates
!-----------------------------------------------------------------------------------------------------------------------------------
!
! END OF FILE
!
!-----------------------------------------------------------------------------------------------------------------------------------
