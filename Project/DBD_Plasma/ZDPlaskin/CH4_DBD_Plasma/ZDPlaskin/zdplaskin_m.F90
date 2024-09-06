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
! Fri Sep  6 18:19:24 2024
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
  integer, parameter :: species_max = 52, species_electrons = 1, species_length = 9, reactions_max = 507, reactions_length = 28
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
  /-1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0,&
    0, 1, 0, 1, 1, 0, 0, 1, 0, 0/
  data species_name(1:species_max) &
  /"E        ","C2H4(V2) ","C3H7^+   ","C2H5^+   ","C3H6^+   ","C2H      ","C2H3^+   ","CH       ","C2H2(V2) ","CH4^+    ",&
   "H2       ","H3^+     ","H        ","C3H8^+   ","C3H6(V)  ","CH3^+    ","C2H2^+   ","C3H7     ","C3H4^+   ","C3H5     ",&
   "C2H2     ","CH4(V24) ","CH3      ","C3H8(V2) ","CH4(V13) ","C3H5^+   ","H2^+     ","CH2      ","C2H3     ","H^+      ",&
   "C4H9     ","C3H8(V1) ","C2H4(V1) ","C2H5     ","C3H4     ","C2H^+    ","C3H6     ","C2H2(V13)","C2H6(V24)","C2H6^+   ",&
   "C2H6(V13)","C2H4     ","C2H6     ","CH^+     ","C3H8     ","C2H4^+   ","CH5^+    ","C5H12    ","C4H9H    ","CH2^+    ",&
   "CH4      ","C2H2(V5) "/
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
  data reaction_sign(217:288) &
  /"bolsig:C2H2(V5)->CHCH^+     ","bolsig:C2H2(V2)->CHCH^+     ","bolsig:C2H2(V13)->CHCH^+    ","bolsig:C3H8(V1)->C3H8^+     ",&
   "bolsig:C3H8(V2)->C3H8^+     ","bolsig:C3H8(V1)->HC3H7^+    ","bolsig:C3H8(V2)->HC3H7^+    ","bolsig:C3H8(V1)->H2C3H6^+   ",&
   "bolsig:C3H8(V2)->H2C3H6^+   ","bolsig:C3H8(V1)->H2HC3H5^+  ","bolsig:C3H8(V2)->H2HC3H5^+  ","bolsig:C3H8(V1)->H2H2C3H4^+ ",&
   "bolsig:C3H8(V2)->H2H2C3H4^+ ","bolsig:C3H8(V1)->CH3C2H5^+  ","bolsig:C3H8(V2)->CH3C2H5^+  ","bolsig:C3H8(V1)->CH4C2H4^+  ",&
   "bolsig:C3H8(V2)->CH4C2H4^+  ","bolsig:C3H8(V1)->C2H5CH3^+  ","bolsig:C3H8(V2)->C2H5CH3^+  ","bolsig:C3H8(V1)->C2H6CH2^+  ",&
   "bolsig:C3H8(V2)->C2H6CH2^+  ","bolsig:H2->HH               ","bolsig:H2->H2^+             ","E+E+CH4^+=>E+CH4            ",&
   "E+E+CH3^+=>E+CH3            ","E+E+CH2^+=>E+CH2            ","E+E+CH^+=>E+CH              ","E+E+C2H6^+=>E+C2H6          ",&
   "E+E+C2H5^+=>E+C2H5          ","E+E+C2H4^+=>E+C2H4          ","E+E+C2H3^+=>E+C2H3          ","E+E+C2H2^+=>E+C2H2          ",&
   "E+E+C3H8^+=>E+C3H8          ","E+E+C3H7^+=>E+C3H7          ","E+E+C3H6^+=>E+C3H6          ","E+E+C3H5^+=>E+C3H5          ",&
   "E+E+C3H4^+=>E+C3H4          ","E+E+H^+=>E+H                ","E+E+H2^+=>E+H2              ","E+CH5^+=>CH3+H+H            ",&
   "E+CH5^+=>CH2+H2+H           ","E+CH4^+=>CH3+H              ","E+CH4^+=>CH2+H+H            ","E+CH4^+=>CH+H2+H            ",&
   "E+CH3^+=>CH2+H              ","E+CH3^+=>CH+H2              ","E+CH2^+=>CH+H               ","E+C2H6^+=>C2H5+H            ",&
   "E+C2H6^+=>C2H4+H+H          ","E+C2H5^+=>C2H4+H            ","E+C2H5^+=>C2H3+H+H          ","E+C2H5^+=>C2H2+H2+H         ",&
   "E+C2H5^+=>C2H2+H+H+H        ","E+C2H5^+=>CH3+CH2           ","E+C2H4^+=>C2H3+H            ","E+C2H4^+=>C2H2+H+H          ",&
   "E+C2H3^+=>C2H2+H            ","E+C2H2^+=>CH+CH             ","E+H2^+=>H+H                 ","E+H3^+=>H+H2                ",&
   "H2+H2^+=>H+H3^+             ","CH2+CH5^+=>CH4+CH3^+        ","CH+CH5^+=>CH4+CH2^+         ","C2H6+CH5^+=>CH4+H2+C2H5^+   ",&
   "C2H4+CH5^+=>CH4+C2H5^+      ","C2H2+CH5^+=>CH4+C2H3^+      ","H+CH5^+=>H2+CH4^+           ","CH4+CH4^+=>CH3+CH5^+        ",&
   "C2H6+CH4^+=>CH4+H2+C2H4^+   ","C2H4+CH4^+=>CH3+C2H5^+      ","C2H4+CH4^+=>CH4+C2H4^+      ","C2H2+CH4^+=>CH3+C2H3^+      "/
  data reaction_sign(289:360) &
  /"C2H2+CH4^+=>CH4+C2H2^+      ","H2+CH4^+=>H+CH5^+           ","H+CH4^+=>H2+CH3^+           ","CH4+CH3^+=>CH3+CH4^+        ",&
   "CH4+CH3^+=>H2+C2H5^+        ","CH2+CH3^+=>H2+C2H3^+        ","CH+CH3^+=>H2+C2H2^+         ","C2H6+CH3^+=>CH4+C2H5^+      ",&
   "C2H4+CH3^+=>CH4+C2H3^+      ","C2H3+CH3^+=>CH3+C2H3^+      ","CH4+CH2^+=>CH3+CH3^+        ","CH4+CH2^+=>H+C2H5^+         ",&
   "CH4+CH2^+=>H2+C2H4^+        ","CH4+CH2^+=>H+H2+C2H3^+      ","CH4+CH2^+=>H2+H2+C2H2^+     ","H2+CH2^+=>H+CH3^+           ",&
   "CH4+CH^+=>H+C2H4^+          ","CH4+CH^+=>H2+C2H3^+         ","CH4+CH^+=>H2+H+C2H2^+       ","H2+CH^+=>H+CH2^+            ",&
   "C2H4+C2H6^+=>C2H6+C2H4^+    ","C2H2+C2H6^+=>C2H3+C2H5^+    ","H+C2H6^+=>H2+C2H5^+         ","H+C2H5^+=>H2+C2H4^+         ",&
   "C2H3+C2H4^+=>C2H2+C2H5^+    ","C2H3+C2H4^+=>C2H4+C2H3^+    ","H+C2H4^+=>H2+C2H3^+         ","C2H6+C2H3^+=>C2H4+C2H5^+    ",&
   "C2H4+C2H3^+=>C2H2+C2H5^+    ","C2H+C2H3^+=>C2H2+C2H2^+     ","H+C2H3^+=>H2+C2H2^+         ","CH4+C2H2^+=>CH3+C2H3^+      ",&
   "C2H6+C2H2^+=>C2H3+C2H5^+    ","C2H6+C2H2^+=>C2H4+C2H4^+    ","C2H4+C2H2^+=>C2H2+C2H4^+    ","C2H3+C2H2^+=>C2H2+C2H3^+    ",&
   "H2+C2H2^+=>H+C2H3^+         ","CH4+C2H^+=>CH3+C2H2^+       ","CH4+H3^+=>H2+CH5^+          ","CH3+H3^+=>H2+CH4^+          ",&
   "CH2+H3^+=>H2+CH3^+          ","CH+H3^+=>H2+CH2^+           ","C2H6+H3^+=>H2+H2+C2H5^+     ","C2H5+H3^+=>H2+C2H6^+        ",&
   "C2H4+H3^+=>H2+C2H5^+        ","C2H4+H3^+=>H2+H2+C2H3^+     ","C2H3+H3^+=>H2+C2H4^+        ","C2H+H3^+=>H2+C2H2^+         ",&
   "C2H2+H3^+=>H2+C2H3^+        ","CH4+H2^+=>H+CH5^+           ","CH4+H2^+=>H2+CH4^+          ","CH4+H2^+=>H2+H+CH3^+        ",&
   "CH2+H2^+=>H+CH3^+           ","CH2+H2^+=>H2+CH2^+          ","CH+H2^+=>H+CH2^+            ","CH+H2^+=>H2+CH^+            ",&
   "C2H6+H2^+=>H2+C2H6^+        ","C2H6+H2^+=>H2+H+C2H5^+      ","C2H6+H2^+=>H2+H2+C2H4^+     ","C2H6+H2^+=>H2+H2+H+C2H3^+   ",&
   "C2H6+H2^+=>H2+H2+H2+C2H2^+  ","C2H4+H2^+=>H2+C2H4^+        ","C2H4+H2^+=>H2+H+C2H3^+      ","C2H4+H2^+=>H2+H2+C2H2^+     ",&
   "C2H2+H2^+=>H+C2H3^+         ","C2H2+H2^+=>H2+C2H2^+        ","H+H2^+=>H3^+                ","H+H2^+=>H2+H^+              ",&
   "CH4+H^+=>H+CH4^+            ","CH4+H^+=>H2+CH3^+           ","CH3+H^+=>H+CH3^+            ","CH2+H^+=>H+CH2^+            "/
  data reaction_sign(361:432) &
  /"CH2+H^+=>H2+CH^+            ","CH+H^+=>H+CH^+              ","C2H6+H^+=>H2+C2H5^+         ","C2H6+H^+=>H2+H+C2H4^+       ",&
   "C2H6+H^+=>H2+H2+C2H3^+      ","C2H5+H^+=>H2+C2H4^+         ","C2H5+H^+=>H2+H+C2H3^+       ","C2H4+H^+=>H+C2H4^+          ",&
   "C2H4+H^+=>H2+C2H3^+         ","C2H4+H^+=>H2+H+C2H2^+       ","C2H3+H^+=>H+C2H3^+          ","C2H3+H^+=>H2+C2H2^+         ",&
   "C2H2+H^+=>H+C2H2^+          ","CH4+CH2=>CH3+CH3            ","CH4+CH=>C2H4+H              ","CH4+C2H5=>C2H6+CH3          ",&
   "CH4+C2H3=>C2H4+CH3          ","CH4+C2H=>C2H2+CH3           ","CH4+C3H7=>C3H8+CH3          ","CH4+C3H5=>C3H6+CH3          ",&
   "CH4+H=>CH3+H2               ","CH4+CH3=>H+C2H6             ","CH4+C4H9=>C4H9H+CH3         ","CH4+CH2=>C2H6               ",&
   "CH4=>CH3+H                  ","CH3=>CH2+H                  ","CH3=>CH+H2                  ","CH3+C2H5=>C2H6+CH2          ",&
   "CH2+CH2=>C2H2+H+H           ","CH2+C2H5=>C2H4+CH3          ","CH2+C2H3=>C2H2+CH3          ","CH2+C2H=>C2H2+CH            ",&
   "CH2+C3H8=>C3H7+CH3          ","CH2+C3H7=>C2H4+C2H5         ","CH2+C3H7=>C3H6+CH3          ","CH2+C3H6=>C3H5+CH3          ",&
   "CH2+H2=>CH3+H               ","CH2+H=>CH+H2                ","CH2=>CH+H                   ","CH2+CH2=>C2H2+H2            ",&
   "CH2+H=>CH3                  ","CH+C2H6=>C3H6+H             ","CH+C2H6=>C3H7               ","CH+H2=>CH2+H                ",&
   "CH+CH3=>C2H3+H              ","CH+CH2=>C2H2+H              ","CH+H2=>CH3                  ","CH+C2H3=>CH2+C2H2           ",&
   "C2H6+C2H3=>C2H5+C2H4        ","C2H6+C3H7=>C3H8+C2H5        ","C2H6+C3H5=>C3H6+C2H5        ","C2H6+H=>C2H5+H2             ",&
   "C2H6+H=>CH4+CH3             ","C2H6=>CH3+CH3               ","C2H6+C4H9=>C4H9H+C2H5       ","C2H6+CH=>C2H4+CH3           ",&
   "C2H6+CH2=>C2H5+CH3          ","C2H5+C2H5=>C2H6+C2H4        ","C2H5+C2H4=>C2H6+C2H3        ","C2H5+C2H2=>C2H6+C2H         ",&
   "C2H5+C2H=>C2H4+C2H2         ","C2H5+C3H8=>C2H6+C3H7        ","C2H5+C3H7=>C3H8+C2H4        ","C2H5+C3H7=>C3H6+C2H6        ",&
   "C2H5+C3H6=>C3H5+C2H6        ","C2H5+H2=>C2H6+H             ","C2H5+H=>CH3+CH3             ","C2H5+H=>C2H4+H2             ",&
   "C2H5+H=>C2H6                ","C2H5=>C2H4+H                ","C2H5+C2H5=>C4H9H            ","C2H5+C4H9=>C2H4+C4H9H       "/
  data reaction_sign(433:504) &
  /"C2H5+C2H3=>C2H4+C2H4        ","C2H4+H=>C2H3+H2             ","C2H4+H=>C2H5                ","C2H4+H2=>C2H5+H             ",&
   "C2H4=>C2H3+H                ","C2H4+C3H6=>C3H5+C2H5        ","C2H4+C2H2=>C2H3+C2H3        ","C2H4+C3H6=>C2H3+C3H7        ",&
   "C2H4+C2H4=>C2H5+C2H3        ","C2H4+CH3=>C3H7              ","C2H4=>C2H2+H2               ","C2H4+C2H5=>C4H9             ",&
   "C2H4+H2=>C2H6               ","C2H4+CH2=>C3H6              ","C2H4+C4H9=>C3H6+C3H7        ","C2H3+C2H3=>C2H4+C2H2        ",&
   "C2H3+C3H8=>C2H4+C3H7        ","C2H3+C3H7=>C3H8+C2H2        ","C2H3+C3H7=>C3H6+C2H4        ","C2H3+C3H6=>C3H5+C2H4        ",&
   "C2H3+C3H5=>C3H6+C2H2        ","C2H3+H2=>C2H4+H             ","C2H3+H=>C2H2+H2             ","C2H3+H=>C2H4                ",&
   "C2H3=>C2H2+H                ","C2H3+C4H9=>C2H2+C4H9H       ","C2H2+H=>C2H3                ","C2H2+H2=>C2H4               ",&
   "C2H2+H2=>C2H3+H             ","C2H2+CH3=>C3H5              ","C2H2+C4H9=>C3H6+C3H5        ","C3H8+C3H5=>C3H6+C3H7        ",&
   "C3H8+H=>C3H7+H2             ","C3H8=>C2H5+CH3              ","C3H8+C4H9=>C4H9H+C3H7       ","C3H8+CH2=>C4H9H             ",&
   "C3H7+C3H7=>C3H6+C3H8        ","C3H7+C3H6=>C3H5+C3H8        ","C3H7+C3H5=>C3H6+C3H6        ","C3H7+H2=>C3H8+H             ",&
   "C3H7+H=>C3H6+H2             ","C3H7+H=>C3H8                ","C3H7+H=>CH3+C2H5            ","C3H7=>C3H6+H                ",&
   "C3H7=>C2H4+CH3              ","C3H7+C4H9=>C4H9H+C3H6       ","C3H6+C2H2=>C2H3+C3H5        ","C3H6+C3H6=>C3H7+C3H5        ",&
   "C3H6=>C3H5+H                ","C3H6+H=>C3H5+H2             ","C3H6+H=>C3H7                ","C3H6=>CH3+C2H3              ",&
   "C3H6+CH3=>C4H9              ","C3H6+C4H9=>C4H9H+C3H5       ","C3H5+H2=>C3H6+H             ","C3H5+H=>C3H6                ",&
   "C3H5=>C2H2+CH3              ","C4H9=>C2H4+C2H5             ","C4H9+CH2=>C2H4+C3H7         ","C4H9=>C3H6+CH3              ",&
   "C4H9+H2=>C4H9H+H            ","C4H9H+CH3=>CH4+C4H9         ","C4H9H=>C3H7+CH3             ","C4H9H=>C2H5+C2H5            ",&
   "C4H9H+H=>C4H9+H2            ","C4H9H+CH2=>C4H9+CH3         ","C4H9H+C2H3=>C2H4+C4H9       ","C4H9H+C3H7=>C3H8+C4H9       ",&
   "C4H9H+C2H=>C2H2+C4H9        ","C4H9H+C2H5=>C2H6+C4H9       ","C4H9H+C3H5=>C3H6+C4H9       ","C4H9H+CH2=>C5H12            "/
  data reaction_sign(505:507) &
  /"C5H12=>CH3+C4H9             ","H2=>H+H                     ","H+H=>H2                     "/
  data bolsig_species(1:bolsig_species_max) &
  /"C2H4(V2) ","CH       ","C2H2(V2) ","H2       ","C3H7     ","C3H5     ","C2H2     ","CH4(V24) ","CH3      ","C3H8(V2) ",&
   "CH4(V13) ","CH2      ","C2H3     ","C3H8(V1) ","C2H4(V1) ","C2H5     ","C3H4     ","C3H6     ","C2H2(V13)","C2H6(V24)",&
   "C2H6(V13)","C2H4     ","C2H6     ","C3H8     ","CH4      ","C2H2(V5) "/
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
311 format(331x,52(1x,i9))
312 format(A3,1x,A28,1x,52(1x,A9))
313 format(i3,1x,A28,1x,52(1x,1pd9.2))
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
      write(ifile_unit,"(507(i3))",err=200) int(mrtm(i,:))
    enddo
    lerror = .false.
200 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file " // &
                                     "<qt_matrix.txt> (subroutine writer_save_qtplaskin)")
    close(ifile_unit)
    open(ifile_unit,file="qt_densities.txt",action="write",err=300)
    write(ifile_unit,"(1x,A14,52(111x,i2.2))",err=300) "Time_s", ( i, i = 1, species_max )
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
    write(ifile_unit,"(1x,A12,507(101x,i3.3))",err=500) "Time_s", ( i, i = 1, reactions_max )
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
      write(ifile_unit,"(1pe15.6,52(1pe13.4))") densav(0,2), densav(1:,2)
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
      write(ifile_unit,"(508(1pe13.4))") densav(0,2), rrt_loc(:)
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
  reac_source_local(22,001) = + reac_rate_local(001) 
  reac_source_local(51,001) = - reac_rate_local(001) 
  reac_source_local(25,002) = + reac_rate_local(002) 
  reac_source_local(51,002) = - reac_rate_local(002) 
  reac_source_local(39,003) = + reac_rate_local(003) 
  reac_source_local(43,003) = - reac_rate_local(003) 
  reac_source_local(41,004) = + reac_rate_local(004) 
  reac_source_local(43,004) = - reac_rate_local(004) 
  reac_source_local(33,005) = + reac_rate_local(005) 
  reac_source_local(42,005) = - reac_rate_local(005) 
  reac_source_local(02,006) = + reac_rate_local(006) 
  reac_source_local(42,006) = - reac_rate_local(006) 
  reac_source_local(21,007) = - reac_rate_local(007) 
  reac_source_local(52,007) = + reac_rate_local(007) 
  reac_source_local(09,008) = + reac_rate_local(008) 
  reac_source_local(21,008) = - reac_rate_local(008) 
  reac_source_local(21,009) = - reac_rate_local(009) 
  reac_source_local(38,009) = + reac_rate_local(009) 
  reac_source_local(32,010) = + reac_rate_local(010) 
  reac_source_local(45,010) = - reac_rate_local(010) 
  reac_source_local(24,011) = + reac_rate_local(011) 
  reac_source_local(45,011) = - reac_rate_local(011) 
  reac_source_local(15,012) = + reac_rate_local(012) 
  reac_source_local(37,012) = - reac_rate_local(012) 
  reac_source_local(13,013) = + reac_rate_local(013) 
  reac_source_local(23,013) = + reac_rate_local(013) 
  reac_source_local(51,013) = - reac_rate_local(013) 
  reac_source_local(11,014) = + reac_rate_local(014) 
  reac_source_local(28,014) = + reac_rate_local(014) 
  reac_source_local(51,014) = - reac_rate_local(014) 
  reac_source_local(08,015) = + reac_rate_local(015) 
  reac_source_local(11,015) = + reac_rate_local(015) 
  reac_source_local(13,015) = + reac_rate_local(015) 
  reac_source_local(51,015) = - reac_rate_local(015) 
  reac_source_local(13,016) = + reac_rate_local(016) 
  reac_source_local(23,016) = - reac_rate_local(016) 
  reac_source_local(28,016) = + reac_rate_local(016) 
  reac_source_local(08,017) = + reac_rate_local(017) 
  reac_source_local(11,017) = + reac_rate_local(017) 
  reac_source_local(23,017) = - reac_rate_local(017) 
  reac_source_local(08,018) = + reac_rate_local(018) 
  reac_source_local(13,018) = + reac_rate_local(018) 
  reac_source_local(28,018) = - reac_rate_local(018) 
  reac_source_local(01,019) = + reac_rate_local(019) 
  reac_source_local(10,019) = + reac_rate_local(019) 
  reac_source_local(51,019) = - reac_rate_local(019) 
  reac_source_local(01,020) = + reac_rate_local(020) 
  reac_source_local(13,020) = + reac_rate_local(020) 
  reac_source_local(16,020) = + reac_rate_local(020) 
  reac_source_local(51,020) = - reac_rate_local(020) 
  reac_source_local(01,021) = + reac_rate_local(021) 
  reac_source_local(11,021) = + reac_rate_local(021) 
  reac_source_local(50,021) = + reac_rate_local(021) 
  reac_source_local(51,021) = - reac_rate_local(021) 
  reac_source_local(01,022) = + reac_rate_local(022) 
  reac_source_local(11,022) = + reac_rate_local(022) 
  reac_source_local(13,022) = + reac_rate_local(022) 
  reac_source_local(44,022) = + reac_rate_local(022) 
  reac_source_local(51,022) = - reac_rate_local(022) 
  reac_source_local(01,023) = + reac_rate_local(023) 
  reac_source_local(16,023) = + reac_rate_local(023) 
  reac_source_local(23,023) = - reac_rate_local(023) 
  reac_source_local(01,024) = + reac_rate_local(024) 
  reac_source_local(13,024) = + reac_rate_local(024) 
  reac_source_local(23,024) = - reac_rate_local(024) 
  reac_source_local(50,024) = + reac_rate_local(024) 
  reac_source_local(01,025) = + reac_rate_local(025) 
  reac_source_local(11,025) = + reac_rate_local(025) 
  reac_source_local(23,025) = - reac_rate_local(025) 
  reac_source_local(44,025) = + reac_rate_local(025) 
  reac_source_local(01,026) = + reac_rate_local(026) 
  reac_source_local(28,026) = - reac_rate_local(026) 
  reac_source_local(50,026) = + reac_rate_local(026) 
  reac_source_local(01,027) = + reac_rate_local(027) 
  reac_source_local(13,027) = + reac_rate_local(027) 
  reac_source_local(28,027) = - reac_rate_local(027) 
  reac_source_local(44,027) = + reac_rate_local(027) 
  reac_source_local(01,028) = + reac_rate_local(028) 
  reac_source_local(08,028) = - reac_rate_local(028) 
  reac_source_local(44,028) = + reac_rate_local(028) 
  reac_source_local(13,029) = + reac_rate_local(029) 
  reac_source_local(34,029) = + reac_rate_local(029) 
  reac_source_local(43,029) = - reac_rate_local(029) 
  reac_source_local(11,030) = + reac_rate_local(030) 
  reac_source_local(42,030) = + reac_rate_local(030) 
  reac_source_local(43,030) = - reac_rate_local(030) 
  reac_source_local(11,031) = + reac_rate_local(031) 
  reac_source_local(13,031) = + reac_rate_local(031) 
  reac_source_local(29,031) = + reac_rate_local(031) 
  reac_source_local(43,031) = - reac_rate_local(031) 
  reac_source_local(11,032) = + reac_rate_local(032) * 2.d0
  reac_source_local(21,032) = + reac_rate_local(032) 
  reac_source_local(43,032) = - reac_rate_local(032) 
  reac_source_local(28,033) = + reac_rate_local(033) 
  reac_source_local(43,033) = - reac_rate_local(033) 
  reac_source_local(51,033) = + reac_rate_local(033) 
  reac_source_local(23,034) = + reac_rate_local(034) * 2.d0
  reac_source_local(43,034) = - reac_rate_local(034) 
  reac_source_local(13,035) = + reac_rate_local(035) 
  reac_source_local(34,035) = - reac_rate_local(035) 
  reac_source_local(42,035) = + reac_rate_local(035) 
  reac_source_local(11,036) = + reac_rate_local(036) 
  reac_source_local(29,036) = + reac_rate_local(036) 
  reac_source_local(34,036) = - reac_rate_local(036) 
  reac_source_local(13,037) = + reac_rate_local(037) * 2.d0
  reac_source_local(29,037) = + reac_rate_local(037) 
  reac_source_local(34,037) = - reac_rate_local(037) 
  reac_source_local(11,038) = + reac_rate_local(038) 
  reac_source_local(13,038) = + reac_rate_local(038) 
  reac_source_local(21,038) = + reac_rate_local(038) 
  reac_source_local(34,038) = - reac_rate_local(038) 
  reac_source_local(08,039) = + reac_rate_local(039) 
  reac_source_local(34,039) = - reac_rate_local(039) 
  reac_source_local(51,039) = + reac_rate_local(039) 
  reac_source_local(23,040) = + reac_rate_local(040) 
  reac_source_local(28,040) = + reac_rate_local(040) 
  reac_source_local(34,040) = - reac_rate_local(040) 
  reac_source_local(13,041) = + reac_rate_local(041) 
  reac_source_local(29,041) = + reac_rate_local(041) 
  reac_source_local(42,041) = - reac_rate_local(041) 
  reac_source_local(11,042) = + reac_rate_local(042) 
  reac_source_local(21,042) = + reac_rate_local(042) 
  reac_source_local(42,042) = - reac_rate_local(042) 
  reac_source_local(13,043) = + reac_rate_local(043) * 2.d0
  reac_source_local(21,043) = + reac_rate_local(043) 
  reac_source_local(42,043) = - reac_rate_local(043) 
  reac_source_local(08,044) = + reac_rate_local(044) 
  reac_source_local(23,044) = + reac_rate_local(044) 
  reac_source_local(42,044) = - reac_rate_local(044) 
  reac_source_local(28,045) = + reac_rate_local(045) * 2.d0
  reac_source_local(42,045) = - reac_rate_local(045) 
  reac_source_local(13,046) = + reac_rate_local(046) 
  reac_source_local(21,046) = + reac_rate_local(046) 
  reac_source_local(29,046) = - reac_rate_local(046) 
  reac_source_local(08,047) = + reac_rate_local(047) 
  reac_source_local(28,047) = + reac_rate_local(047) 
  reac_source_local(29,047) = - reac_rate_local(047) 
  reac_source_local(08,048) = + reac_rate_local(048) * 2.d0
  reac_source_local(21,048) = - reac_rate_local(048) 
  reac_source_local(13,049) = + reac_rate_local(049) 
  reac_source_local(18,049) = + reac_rate_local(049) 
  reac_source_local(45,049) = - reac_rate_local(049) 
  reac_source_local(11,050) = + reac_rate_local(050) 
  reac_source_local(37,050) = + reac_rate_local(050) 
  reac_source_local(45,050) = - reac_rate_local(050) 
  reac_source_local(11,051) = + reac_rate_local(051) * 2.d0
  reac_source_local(35,051) = + reac_rate_local(051) 
  reac_source_local(45,051) = - reac_rate_local(051) 
  reac_source_local(28,052) = + reac_rate_local(052) 
  reac_source_local(43,052) = + reac_rate_local(052) 
  reac_source_local(45,052) = - reac_rate_local(052) 
  reac_source_local(23,053) = + reac_rate_local(053) 
  reac_source_local(34,053) = + reac_rate_local(053) 
  reac_source_local(45,053) = - reac_rate_local(053) 
  reac_source_local(42,054) = + reac_rate_local(054) 
  reac_source_local(45,054) = - reac_rate_local(054) 
  reac_source_local(51,054) = + reac_rate_local(054) 
  reac_source_local(13,055) = + reac_rate_local(055) 
  reac_source_local(18,055) = - reac_rate_local(055) 
  reac_source_local(37,055) = + reac_rate_local(055) 
  reac_source_local(11,056) = + reac_rate_local(056) 
  reac_source_local(18,056) = - reac_rate_local(056) 
  reac_source_local(20,056) = + reac_rate_local(056) 
  reac_source_local(11,057) = + reac_rate_local(057) 
  reac_source_local(13,057) = + reac_rate_local(057) 
  reac_source_local(18,057) = - reac_rate_local(057) 
  reac_source_local(35,057) = + reac_rate_local(057) 
  reac_source_local(18,058) = - reac_rate_local(058) 
  reac_source_local(23,058) = + reac_rate_local(058) 
  reac_source_local(42,058) = + reac_rate_local(058) 
  reac_source_local(18,059) = - reac_rate_local(059) 
  reac_source_local(29,059) = + reac_rate_local(059) 
  reac_source_local(51,059) = + reac_rate_local(059) 
  reac_source_local(13,060) = + reac_rate_local(060) 
  reac_source_local(20,060) = + reac_rate_local(060) 
  reac_source_local(37,060) = - reac_rate_local(060) 
  reac_source_local(11,061) = + reac_rate_local(061) 
  reac_source_local(35,061) = + reac_rate_local(061) 
  reac_source_local(37,061) = - reac_rate_local(061) 
  reac_source_local(28,062) = + reac_rate_local(062) 
  reac_source_local(37,062) = - reac_rate_local(062) 
  reac_source_local(42,062) = + reac_rate_local(062) 
  reac_source_local(23,063) = + reac_rate_local(063) 
  reac_source_local(29,063) = + reac_rate_local(063) 
  reac_source_local(37,063) = - reac_rate_local(063) 
  reac_source_local(21,064) = + reac_rate_local(064) 
  reac_source_local(37,064) = - reac_rate_local(064) 
  reac_source_local(51,064) = + reac_rate_local(064) 
  reac_source_local(13,065) = + reac_rate_local(065) 
  reac_source_local(20,065) = - reac_rate_local(065) 
  reac_source_local(35,065) = + reac_rate_local(065) 
  reac_source_local(20,066) = - reac_rate_local(066) 
  reac_source_local(21,066) = + reac_rate_local(066) 
  reac_source_local(23,066) = + reac_rate_local(066) 
  reac_source_local(08,067) = + reac_rate_local(067) 
  reac_source_local(29,067) = + reac_rate_local(067) 
  reac_source_local(35,067) = - reac_rate_local(067) 
  reac_source_local(21,068) = + reac_rate_local(068) 
  reac_source_local(28,068) = + reac_rate_local(068) 
  reac_source_local(35,068) = - reac_rate_local(068) 
  reac_source_local(01,069) = + reac_rate_local(069) 
  reac_source_local(40,069) = + reac_rate_local(069) 
  reac_source_local(43,069) = - reac_rate_local(069) 
  reac_source_local(01,070) = + reac_rate_local(070) 
  reac_source_local(04,070) = + reac_rate_local(070) 
  reac_source_local(13,070) = + reac_rate_local(070) 
  reac_source_local(43,070) = - reac_rate_local(070) 
  reac_source_local(01,071) = + reac_rate_local(071) 
  reac_source_local(11,071) = + reac_rate_local(071) 
  reac_source_local(43,071) = - reac_rate_local(071) 
  reac_source_local(46,071) = + reac_rate_local(071) 
  reac_source_local(01,072) = + reac_rate_local(072) 
  reac_source_local(07,072) = + reac_rate_local(072) 
  reac_source_local(11,072) = + reac_rate_local(072) 
  reac_source_local(13,072) = + reac_rate_local(072) 
  reac_source_local(43,072) = - reac_rate_local(072) 
  reac_source_local(01,073) = + reac_rate_local(073) 
  reac_source_local(11,073) = + reac_rate_local(073) * 2.d0
  reac_source_local(17,073) = + reac_rate_local(073) 
  reac_source_local(43,073) = - reac_rate_local(073) 
  reac_source_local(01,074) = + reac_rate_local(074) 
  reac_source_local(16,074) = + reac_rate_local(074) 
  reac_source_local(23,074) = + reac_rate_local(074) 
  reac_source_local(43,074) = - reac_rate_local(074) 
  reac_source_local(01,075) = + reac_rate_local(075) 
  reac_source_local(43,075) = - reac_rate_local(075) 
  reac_source_local(50,075) = + reac_rate_local(075) 
  reac_source_local(51,075) = + reac_rate_local(075) 
  reac_source_local(01,076) = + reac_rate_local(076) 
  reac_source_local(04,076) = + reac_rate_local(076) 
  reac_source_local(34,076) = - reac_rate_local(076) 
  reac_source_local(01,077) = + reac_rate_local(077) 
  reac_source_local(13,077) = + reac_rate_local(077) 
  reac_source_local(34,077) = - reac_rate_local(077) 
  reac_source_local(46,077) = + reac_rate_local(077) 
  reac_source_local(01,078) = + reac_rate_local(078) 
  reac_source_local(07,078) = + reac_rate_local(078) 
  reac_source_local(11,078) = + reac_rate_local(078) 
  reac_source_local(34,078) = - reac_rate_local(078) 
  reac_source_local(01,079) = + reac_rate_local(079) 
  reac_source_local(11,079) = + reac_rate_local(079) 
  reac_source_local(13,079) = + reac_rate_local(079) 
  reac_source_local(17,079) = + reac_rate_local(079) 
  reac_source_local(34,079) = - reac_rate_local(079) 
  reac_source_local(01,080) = + reac_rate_local(080) 
  reac_source_local(16,080) = + reac_rate_local(080) 
  reac_source_local(28,080) = + reac_rate_local(080) 
  reac_source_local(34,080) = - reac_rate_local(080) 
  reac_source_local(01,081) = + reac_rate_local(081) 
  reac_source_local(23,081) = + reac_rate_local(081) 
  reac_source_local(34,081) = - reac_rate_local(081) 
  reac_source_local(50,081) = + reac_rate_local(081) 
  reac_source_local(01,082) = + reac_rate_local(082) 
  reac_source_local(34,082) = - reac_rate_local(082) 
  reac_source_local(44,082) = + reac_rate_local(082) 
  reac_source_local(51,082) = + reac_rate_local(082) 
  reac_source_local(01,083) = + reac_rate_local(083) 
  reac_source_local(42,083) = - reac_rate_local(083) 
  reac_source_local(46,083) = + reac_rate_local(083) 
  reac_source_local(01,084) = + reac_rate_local(084) 
  reac_source_local(07,084) = + reac_rate_local(084) 
  reac_source_local(13,084) = + reac_rate_local(084) 
  reac_source_local(42,084) = - reac_rate_local(084) 
  reac_source_local(01,085) = + reac_rate_local(085) 
  reac_source_local(08,085) = + reac_rate_local(085) 
  reac_source_local(16,085) = + reac_rate_local(085) 
  reac_source_local(42,085) = - reac_rate_local(085) 
  reac_source_local(01,086) = + reac_rate_local(086) 
  reac_source_local(28,086) = + reac_rate_local(086) 
  reac_source_local(42,086) = - reac_rate_local(086) 
  reac_source_local(50,086) = + reac_rate_local(086) 
  reac_source_local(01,087) = + reac_rate_local(087) 
  reac_source_local(23,087) = + reac_rate_local(087) 
  reac_source_local(42,087) = - reac_rate_local(087) 
  reac_source_local(44,087) = + reac_rate_local(087) 
  reac_source_local(01,088) = + reac_rate_local(088) 
  reac_source_local(07,088) = + reac_rate_local(088) 
  reac_source_local(29,088) = - reac_rate_local(088) 
  reac_source_local(01,089) = + reac_rate_local(089) 
  reac_source_local(13,089) = + reac_rate_local(089) 
  reac_source_local(17,089) = + reac_rate_local(089) 
  reac_source_local(29,089) = - reac_rate_local(089) 
  reac_source_local(01,090) = + reac_rate_local(090) 
  reac_source_local(08,090) = + reac_rate_local(090) 
  reac_source_local(29,090) = - reac_rate_local(090) 
  reac_source_local(50,090) = + reac_rate_local(090) 
  reac_source_local(01,091) = + reac_rate_local(091) 
  reac_source_local(28,091) = + reac_rate_local(091) 
  reac_source_local(29,091) = - reac_rate_local(091) 
  reac_source_local(44,091) = + reac_rate_local(091) 
  reac_source_local(01,092) = + reac_rate_local(092) 
  reac_source_local(21,092) = + reac_rate_local(092) 
  reac_source_local(29,092) = - reac_rate_local(092) 
  reac_source_local(30,092) = + reac_rate_local(092) 
  reac_source_local(01,093) = + reac_rate_local(093) 
  reac_source_local(17,093) = + reac_rate_local(093) 
  reac_source_local(21,093) = - reac_rate_local(093) 
  reac_source_local(01,094) = + reac_rate_local(094) 
  reac_source_local(08,094) = + reac_rate_local(094) 
  reac_source_local(21,094) = - reac_rate_local(094) 
  reac_source_local(44,094) = + reac_rate_local(094) 
  reac_source_local(01,095) = + reac_rate_local(095) 
  reac_source_local(14,095) = + reac_rate_local(095) 
  reac_source_local(45,095) = - reac_rate_local(095) 
  reac_source_local(01,096) = + reac_rate_local(096) 
  reac_source_local(03,096) = + reac_rate_local(096) 
  reac_source_local(13,096) = + reac_rate_local(096) 
  reac_source_local(45,096) = - reac_rate_local(096) 
  reac_source_local(01,097) = + reac_rate_local(097) 
  reac_source_local(05,097) = + reac_rate_local(097) 
  reac_source_local(11,097) = + reac_rate_local(097) 
  reac_source_local(45,097) = - reac_rate_local(097) 
  reac_source_local(01,098) = + reac_rate_local(098) 
  reac_source_local(11,098) = + reac_rate_local(098) 
  reac_source_local(13,098) = + reac_rate_local(098) 
  reac_source_local(26,098) = + reac_rate_local(098) 
  reac_source_local(45,098) = - reac_rate_local(098) 
  reac_source_local(01,099) = + reac_rate_local(099) 
  reac_source_local(11,099) = + reac_rate_local(099) * 2.d0
  reac_source_local(19,099) = + reac_rate_local(099) 
  reac_source_local(45,099) = - reac_rate_local(099) 
  reac_source_local(01,100) = + reac_rate_local(100) 
  reac_source_local(04,100) = + reac_rate_local(100) 
  reac_source_local(23,100) = + reac_rate_local(100) 
  reac_source_local(45,100) = - reac_rate_local(100) 
  reac_source_local(01,101) = + reac_rate_local(101) 
  reac_source_local(45,101) = - reac_rate_local(101) 
  reac_source_local(46,101) = + reac_rate_local(101) 
  reac_source_local(51,101) = + reac_rate_local(101) 
  reac_source_local(01,102) = + reac_rate_local(102) 
  reac_source_local(16,102) = + reac_rate_local(102) 
  reac_source_local(34,102) = + reac_rate_local(102) 
  reac_source_local(45,102) = - reac_rate_local(102) 
  reac_source_local(01,103) = + reac_rate_local(103) 
  reac_source_local(43,103) = + reac_rate_local(103) 
  reac_source_local(45,103) = - reac_rate_local(103) 
  reac_source_local(50,103) = + reac_rate_local(103) 
  reac_source_local(01,104) = + reac_rate_local(104) 
  reac_source_local(03,104) = + reac_rate_local(104) 
  reac_source_local(18,104) = - reac_rate_local(104) 
  reac_source_local(01,105) = + reac_rate_local(105) 
  reac_source_local(05,105) = + reac_rate_local(105) 
  reac_source_local(13,105) = + reac_rate_local(105) 
  reac_source_local(18,105) = - reac_rate_local(105) 
  reac_source_local(01,106) = + reac_rate_local(106) 
  reac_source_local(11,106) = + reac_rate_local(106) 
  reac_source_local(18,106) = - reac_rate_local(106) 
  reac_source_local(26,106) = + reac_rate_local(106) 
  reac_source_local(01,107) = + reac_rate_local(107) 
  reac_source_local(11,107) = + reac_rate_local(107) 
  reac_source_local(13,107) = + reac_rate_local(107) 
  reac_source_local(18,107) = - reac_rate_local(107) 
  reac_source_local(19,107) = + reac_rate_local(107) 
  reac_source_local(01,108) = + reac_rate_local(108) 
  reac_source_local(04,108) = + reac_rate_local(108) 
  reac_source_local(18,108) = - reac_rate_local(108) 
  reac_source_local(28,108) = + reac_rate_local(108) 
  reac_source_local(01,109) = + reac_rate_local(109) 
  reac_source_local(18,109) = - reac_rate_local(109) 
  reac_source_local(23,109) = + reac_rate_local(109) 
  reac_source_local(46,109) = + reac_rate_local(109) 
  reac_source_local(01,110) = + reac_rate_local(110) 
  reac_source_local(07,110) = + reac_rate_local(110) 
  reac_source_local(18,110) = - reac_rate_local(110) 
  reac_source_local(51,110) = + reac_rate_local(110) 
  reac_source_local(01,111) = + reac_rate_local(111) 
  reac_source_local(10,111) = + reac_rate_local(111) 
  reac_source_local(18,111) = - reac_rate_local(111) 
  reac_source_local(29,111) = + reac_rate_local(111) 
  reac_source_local(01,112) = + reac_rate_local(112) 
  reac_source_local(16,112) = + reac_rate_local(112) 
  reac_source_local(18,112) = - reac_rate_local(112) 
  reac_source_local(42,112) = + reac_rate_local(112) 
  reac_source_local(01,113) = + reac_rate_local(113) 
  reac_source_local(18,113) = - reac_rate_local(113) 
  reac_source_local(34,113) = + reac_rate_local(113) 
  reac_source_local(50,113) = + reac_rate_local(113) 
  reac_source_local(01,114) = + reac_rate_local(114) 
  reac_source_local(18,114) = - reac_rate_local(114) 
  reac_source_local(43,114) = + reac_rate_local(114) 
  reac_source_local(44,114) = + reac_rate_local(114) 
  reac_source_local(01,115) = + reac_rate_local(115) 
  reac_source_local(05,115) = + reac_rate_local(115) 
  reac_source_local(37,115) = - reac_rate_local(115) 
  reac_source_local(01,116) = + reac_rate_local(116) 
  reac_source_local(13,116) = + reac_rate_local(116) 
  reac_source_local(26,116) = + reac_rate_local(116) 
  reac_source_local(37,116) = - reac_rate_local(116) 
  reac_source_local(01,117) = + reac_rate_local(117) 
  reac_source_local(11,117) = + reac_rate_local(117) 
  reac_source_local(19,117) = + reac_rate_local(117) 
  reac_source_local(37,117) = - reac_rate_local(117) 
  reac_source_local(01,118) = + reac_rate_local(118) 
  reac_source_local(04,118) = + reac_rate_local(118) 
  reac_source_local(08,118) = + reac_rate_local(118) 
  reac_source_local(37,118) = - reac_rate_local(118) 
  reac_source_local(01,119) = + reac_rate_local(119) 
  reac_source_local(28,119) = + reac_rate_local(119) 
  reac_source_local(37,119) = - reac_rate_local(119) 
  reac_source_local(46,119) = + reac_rate_local(119) 
  reac_source_local(01,120) = + reac_rate_local(120) 
  reac_source_local(07,120) = + reac_rate_local(120) 
  reac_source_local(23,120) = + reac_rate_local(120) 
  reac_source_local(37,120) = - reac_rate_local(120) 
  reac_source_local(01,121) = + reac_rate_local(121) 
  reac_source_local(17,121) = + reac_rate_local(121) 
  reac_source_local(37,121) = - reac_rate_local(121) 
  reac_source_local(51,121) = + reac_rate_local(121) 
  reac_source_local(01,122) = + reac_rate_local(122) 
  reac_source_local(10,122) = + reac_rate_local(122) 
  reac_source_local(21,122) = + reac_rate_local(122) 
  reac_source_local(37,122) = - reac_rate_local(122) 
  reac_source_local(01,123) = + reac_rate_local(123) 
  reac_source_local(16,123) = + reac_rate_local(123) 
  reac_source_local(29,123) = + reac_rate_local(123) 
  reac_source_local(37,123) = - reac_rate_local(123) 
  reac_source_local(01,124) = + reac_rate_local(124) 
  reac_source_local(37,124) = - reac_rate_local(124) 
  reac_source_local(42,124) = + reac_rate_local(124) 
  reac_source_local(50,124) = + reac_rate_local(124) 
  reac_source_local(01,125) = + reac_rate_local(125) 
  reac_source_local(34,125) = + reac_rate_local(125) 
  reac_source_local(37,125) = - reac_rate_local(125) 
  reac_source_local(44,125) = + reac_rate_local(125) 
  reac_source_local(01,126) = + reac_rate_local(126) 
  reac_source_local(20,126) = - reac_rate_local(126) 
  reac_source_local(26,126) = + reac_rate_local(126) 
  reac_source_local(01,127) = + reac_rate_local(127) 
  reac_source_local(13,127) = + reac_rate_local(127) 
  reac_source_local(19,127) = + reac_rate_local(127) 
  reac_source_local(20,127) = - reac_rate_local(127) 
  reac_source_local(01,128) = + reac_rate_local(128) 
  reac_source_local(08,128) = + reac_rate_local(128) 
  reac_source_local(20,128) = - reac_rate_local(128) 
  reac_source_local(46,128) = + reac_rate_local(128) 
  reac_source_local(01,129) = + reac_rate_local(129) 
  reac_source_local(07,129) = + reac_rate_local(129) 
  reac_source_local(20,129) = - reac_rate_local(129) 
  reac_source_local(28,129) = + reac_rate_local(129) 
  reac_source_local(01,130) = + reac_rate_local(130) 
  reac_source_local(17,130) = + reac_rate_local(130) 
  reac_source_local(20,130) = - reac_rate_local(130) 
  reac_source_local(23,130) = + reac_rate_local(130) 
  reac_source_local(01,131) = + reac_rate_local(131) 
  reac_source_local(16,131) = + reac_rate_local(131) 
  reac_source_local(20,131) = - reac_rate_local(131) 
  reac_source_local(21,131) = + reac_rate_local(131) 
  reac_source_local(01,132) = + reac_rate_local(132) 
  reac_source_local(20,132) = - reac_rate_local(132) 
  reac_source_local(29,132) = + reac_rate_local(132) 
  reac_source_local(50,132) = + reac_rate_local(132) 
  reac_source_local(01,133) = + reac_rate_local(133) 
  reac_source_local(20,133) = - reac_rate_local(133) 
  reac_source_local(42,133) = + reac_rate_local(133) 
  reac_source_local(44,133) = + reac_rate_local(133) 
  reac_source_local(01,134) = + reac_rate_local(134) 
  reac_source_local(19,134) = + reac_rate_local(134) 
  reac_source_local(35,134) = - reac_rate_local(134) 
  reac_source_local(01,135) = + reac_rate_local(135) 
  reac_source_local(07,135) = + reac_rate_local(135) 
  reac_source_local(08,135) = + reac_rate_local(135) 
  reac_source_local(35,135) = - reac_rate_local(135) 
  reac_source_local(01,136) = + reac_rate_local(136) 
  reac_source_local(17,136) = + reac_rate_local(136) 
  reac_source_local(28,136) = + reac_rate_local(136) 
  reac_source_local(35,136) = - reac_rate_local(136) 
  reac_source_local(01,137) = + reac_rate_local(137) 
  reac_source_local(21,137) = + reac_rate_local(137) 
  reac_source_local(35,137) = - reac_rate_local(137) 
  reac_source_local(50,137) = + reac_rate_local(137) 
  reac_source_local(01,138) = + reac_rate_local(138) 
  reac_source_local(29,138) = + reac_rate_local(138) 
  reac_source_local(35,138) = - reac_rate_local(138) 
  reac_source_local(44,138) = + reac_rate_local(138) 
  reac_source_local(13,139) = + reac_rate_local(139) 
  reac_source_local(22,139) = - reac_rate_local(139) 
  reac_source_local(23,139) = + reac_rate_local(139) 
  reac_source_local(13,140) = + reac_rate_local(140) 
  reac_source_local(23,140) = + reac_rate_local(140) 
  reac_source_local(25,140) = - reac_rate_local(140) 
  reac_source_local(11,141) = + reac_rate_local(141) 
  reac_source_local(22,141) = - reac_rate_local(141) 
  reac_source_local(28,141) = + reac_rate_local(141) 
  reac_source_local(11,142) = + reac_rate_local(142) 
  reac_source_local(25,142) = - reac_rate_local(142) 
  reac_source_local(28,142) = + reac_rate_local(142) 
  reac_source_local(08,143) = + reac_rate_local(143) 
  reac_source_local(11,143) = + reac_rate_local(143) 
  reac_source_local(13,143) = + reac_rate_local(143) 
  reac_source_local(22,143) = - reac_rate_local(143) 
  reac_source_local(08,144) = + reac_rate_local(144) 
  reac_source_local(11,144) = + reac_rate_local(144) 
  reac_source_local(13,144) = + reac_rate_local(144) 
  reac_source_local(25,144) = - reac_rate_local(144) 
  reac_source_local(01,145) = + reac_rate_local(145) 
  reac_source_local(10,145) = + reac_rate_local(145) 
  reac_source_local(22,145) = - reac_rate_local(145) 
  reac_source_local(01,146) = + reac_rate_local(146) 
  reac_source_local(10,146) = + reac_rate_local(146) 
  reac_source_local(25,146) = - reac_rate_local(146) 
  reac_source_local(01,147) = + reac_rate_local(147) 
  reac_source_local(13,147) = + reac_rate_local(147) 
  reac_source_local(16,147) = + reac_rate_local(147) 
  reac_source_local(22,147) = - reac_rate_local(147) 
  reac_source_local(01,148) = + reac_rate_local(148) 
  reac_source_local(13,148) = + reac_rate_local(148) 
  reac_source_local(16,148) = + reac_rate_local(148) 
  reac_source_local(25,148) = - reac_rate_local(148) 
  reac_source_local(01,149) = + reac_rate_local(149) 
  reac_source_local(11,149) = + reac_rate_local(149) 
  reac_source_local(22,149) = - reac_rate_local(149) 
  reac_source_local(50,149) = + reac_rate_local(149) 
  reac_source_local(01,150) = + reac_rate_local(150) 
  reac_source_local(11,150) = + reac_rate_local(150) 
  reac_source_local(25,150) = - reac_rate_local(150) 
  reac_source_local(50,150) = + reac_rate_local(150) 
  reac_source_local(01,151) = + reac_rate_local(151) 
  reac_source_local(11,151) = + reac_rate_local(151) 
  reac_source_local(13,151) = + reac_rate_local(151) 
  reac_source_local(22,151) = - reac_rate_local(151) 
  reac_source_local(44,151) = + reac_rate_local(151) 
  reac_source_local(01,152) = + reac_rate_local(152) 
  reac_source_local(11,152) = + reac_rate_local(152) 
  reac_source_local(13,152) = + reac_rate_local(152) 
  reac_source_local(25,152) = - reac_rate_local(152) 
  reac_source_local(44,152) = + reac_rate_local(152) 
  reac_source_local(13,153) = + reac_rate_local(153) 
  reac_source_local(34,153) = + reac_rate_local(153) 
  reac_source_local(39,153) = - reac_rate_local(153) 
  reac_source_local(13,154) = + reac_rate_local(154) 
  reac_source_local(34,154) = + reac_rate_local(154) 
  reac_source_local(41,154) = - reac_rate_local(154) 
  reac_source_local(11,155) = + reac_rate_local(155) 
  reac_source_local(39,155) = - reac_rate_local(155) 
  reac_source_local(42,155) = + reac_rate_local(155) 
  reac_source_local(11,156) = + reac_rate_local(156) 
  reac_source_local(41,156) = - reac_rate_local(156) 
  reac_source_local(42,156) = + reac_rate_local(156) 
  reac_source_local(11,157) = + reac_rate_local(157) 
  reac_source_local(13,157) = + reac_rate_local(157) 
  reac_source_local(29,157) = + reac_rate_local(157) 
  reac_source_local(39,157) = - reac_rate_local(157) 
  reac_source_local(11,158) = + reac_rate_local(158) 
  reac_source_local(13,158) = + reac_rate_local(158) 
  reac_source_local(29,158) = + reac_rate_local(158) 
  reac_source_local(41,158) = - reac_rate_local(158) 
  reac_source_local(11,159) = + reac_rate_local(159) * 2.d0
  reac_source_local(21,159) = + reac_rate_local(159) 
  reac_source_local(39,159) = - reac_rate_local(159) 
  reac_source_local(11,160) = + reac_rate_local(160) * 2.d0
  reac_source_local(21,160) = + reac_rate_local(160) 
  reac_source_local(41,160) = - reac_rate_local(160) 
  reac_source_local(28,161) = + reac_rate_local(161) 
  reac_source_local(39,161) = - reac_rate_local(161) 
  reac_source_local(51,161) = + reac_rate_local(161) 
  reac_source_local(28,162) = + reac_rate_local(162) 
  reac_source_local(41,162) = - reac_rate_local(162) 
  reac_source_local(51,162) = + reac_rate_local(162) 
  reac_source_local(23,163) = + reac_rate_local(163) * 2.d0
  reac_source_local(39,163) = - reac_rate_local(163) 
  reac_source_local(23,164) = + reac_rate_local(164) * 2.d0
  reac_source_local(41,164) = - reac_rate_local(164) 
  reac_source_local(13,165) = + reac_rate_local(165) 
  reac_source_local(29,165) = + reac_rate_local(165) 
  reac_source_local(33,165) = - reac_rate_local(165) 
  reac_source_local(02,166) = - reac_rate_local(166) 
  reac_source_local(13,166) = + reac_rate_local(166) 
  reac_source_local(29,166) = + reac_rate_local(166) 
  reac_source_local(11,167) = + reac_rate_local(167) 
  reac_source_local(21,167) = + reac_rate_local(167) 
  reac_source_local(33,167) = - reac_rate_local(167) 
  reac_source_local(02,168) = - reac_rate_local(168) 
  reac_source_local(11,168) = + reac_rate_local(168) 
  reac_source_local(21,168) = + reac_rate_local(168) 
  reac_source_local(13,169) = + reac_rate_local(169) * 2.d0
  reac_source_local(21,169) = + reac_rate_local(169) 
  reac_source_local(33,169) = - reac_rate_local(169) 
  reac_source_local(02,170) = - reac_rate_local(170) 
  reac_source_local(13,170) = + reac_rate_local(170) * 2.d0
  reac_source_local(21,170) = + reac_rate_local(170) 
  reac_source_local(08,171) = + reac_rate_local(171) 
  reac_source_local(23,171) = + reac_rate_local(171) 
  reac_source_local(33,171) = - reac_rate_local(171) 
  reac_source_local(02,172) = - reac_rate_local(172) 
  reac_source_local(08,172) = + reac_rate_local(172) 
  reac_source_local(23,172) = + reac_rate_local(172) 
  reac_source_local(28,173) = + reac_rate_local(173) * 2.d0
  reac_source_local(33,173) = - reac_rate_local(173) 
  reac_source_local(02,174) = - reac_rate_local(174) 
  reac_source_local(28,174) = + reac_rate_local(174) * 2.d0
  reac_source_local(08,175) = + reac_rate_local(175) * 2.d0
  reac_source_local(52,175) = - reac_rate_local(175) 
  reac_source_local(08,176) = + reac_rate_local(176) * 2.d0
  reac_source_local(09,176) = - reac_rate_local(176) 
  reac_source_local(08,177) = + reac_rate_local(177) * 2.d0
  reac_source_local(38,177) = - reac_rate_local(177) 
  reac_source_local(13,178) = + reac_rate_local(178) 
  reac_source_local(18,178) = + reac_rate_local(178) 
  reac_source_local(32,178) = - reac_rate_local(178) 
  reac_source_local(13,179) = + reac_rate_local(179) 
  reac_source_local(18,179) = + reac_rate_local(179) 
  reac_source_local(24,179) = - reac_rate_local(179) 
  reac_source_local(11,180) = + reac_rate_local(180) 
  reac_source_local(32,180) = - reac_rate_local(180) 
  reac_source_local(37,180) = + reac_rate_local(180) 
  reac_source_local(11,181) = + reac_rate_local(181) 
  reac_source_local(24,181) = - reac_rate_local(181) 
  reac_source_local(37,181) = + reac_rate_local(181) 
  reac_source_local(11,182) = + reac_rate_local(182) * 2.d0
  reac_source_local(32,182) = - reac_rate_local(182) 
  reac_source_local(35,182) = + reac_rate_local(182) 
  reac_source_local(11,183) = + reac_rate_local(183) * 2.d0
  reac_source_local(24,183) = - reac_rate_local(183) 
  reac_source_local(35,183) = + reac_rate_local(183) 
  reac_source_local(28,184) = + reac_rate_local(184) 
  reac_source_local(32,184) = - reac_rate_local(184) 
  reac_source_local(43,184) = + reac_rate_local(184) 
  reac_source_local(24,185) = - reac_rate_local(185) 
  reac_source_local(28,185) = + reac_rate_local(185) 
  reac_source_local(43,185) = + reac_rate_local(185) 
  reac_source_local(23,186) = + reac_rate_local(186) 
  reac_source_local(32,186) = - reac_rate_local(186) 
  reac_source_local(34,186) = + reac_rate_local(186) 
  reac_source_local(23,187) = + reac_rate_local(187) 
  reac_source_local(24,187) = - reac_rate_local(187) 
  reac_source_local(34,187) = + reac_rate_local(187) 
  reac_source_local(32,188) = - reac_rate_local(188) 
  reac_source_local(42,188) = + reac_rate_local(188) 
  reac_source_local(51,188) = + reac_rate_local(188) 
  reac_source_local(24,189) = - reac_rate_local(189) 
  reac_source_local(42,189) = + reac_rate_local(189) 
  reac_source_local(51,189) = + reac_rate_local(189) 
  reac_source_local(01,190) = + reac_rate_local(190) 
  reac_source_local(39,190) = - reac_rate_local(190) 
  reac_source_local(40,190) = + reac_rate_local(190) 
  reac_source_local(01,191) = + reac_rate_local(191) 
  reac_source_local(40,191) = + reac_rate_local(191) 
  reac_source_local(41,191) = - reac_rate_local(191) 
  reac_source_local(01,192) = + reac_rate_local(192) 
  reac_source_local(04,192) = + reac_rate_local(192) 
  reac_source_local(13,192) = + reac_rate_local(192) 
  reac_source_local(39,192) = - reac_rate_local(192) 
  reac_source_local(01,193) = + reac_rate_local(193) 
  reac_source_local(04,193) = + reac_rate_local(193) 
  reac_source_local(13,193) = + reac_rate_local(193) 
  reac_source_local(41,193) = - reac_rate_local(193) 
  reac_source_local(01,194) = + reac_rate_local(194) 
  reac_source_local(11,194) = + reac_rate_local(194) 
  reac_source_local(39,194) = - reac_rate_local(194) 
  reac_source_local(46,194) = + reac_rate_local(194) 
  reac_source_local(01,195) = + reac_rate_local(195) 
  reac_source_local(11,195) = + reac_rate_local(195) 
  reac_source_local(41,195) = - reac_rate_local(195) 
  reac_source_local(46,195) = + reac_rate_local(195) 
  reac_source_local(01,196) = + reac_rate_local(196) 
  reac_source_local(07,196) = + reac_rate_local(196) 
  reac_source_local(11,196) = + reac_rate_local(196) 
  reac_source_local(13,196) = + reac_rate_local(196) 
  reac_source_local(39,196) = - reac_rate_local(196) 
  reac_source_local(01,197) = + reac_rate_local(197) 
  reac_source_local(07,197) = + reac_rate_local(197) 
  reac_source_local(11,197) = + reac_rate_local(197) 
  reac_source_local(13,197) = + reac_rate_local(197) 
  reac_source_local(41,197) = - reac_rate_local(197) 
  reac_source_local(01,198) = + reac_rate_local(198) 
  reac_source_local(11,198) = + reac_rate_local(198) * 2.d0
  reac_source_local(17,198) = + reac_rate_local(198) 
  reac_source_local(39,198) = - reac_rate_local(198) 
  reac_source_local(01,199) = + reac_rate_local(199) 
  reac_source_local(11,199) = + reac_rate_local(199) * 2.d0
  reac_source_local(17,199) = + reac_rate_local(199) 
  reac_source_local(41,199) = - reac_rate_local(199) 
  reac_source_local(01,200) = + reac_rate_local(200) 
  reac_source_local(16,200) = + reac_rate_local(200) 
  reac_source_local(23,200) = + reac_rate_local(200) 
  reac_source_local(39,200) = - reac_rate_local(200) 
  reac_source_local(01,201) = + reac_rate_local(201) 
  reac_source_local(16,201) = + reac_rate_local(201) 
  reac_source_local(23,201) = + reac_rate_local(201) 
  reac_source_local(41,201) = - reac_rate_local(201) 
  reac_source_local(01,202) = + reac_rate_local(202) 
  reac_source_local(39,202) = - reac_rate_local(202) 
  reac_source_local(50,202) = + reac_rate_local(202) 
  reac_source_local(51,202) = + reac_rate_local(202) 
  reac_source_local(01,203) = + reac_rate_local(203) 
  reac_source_local(41,203) = - reac_rate_local(203) 
  reac_source_local(50,203) = + reac_rate_local(203) 
  reac_source_local(51,203) = + reac_rate_local(203) 
  reac_source_local(01,204) = + reac_rate_local(204) 
  reac_source_local(33,204) = - reac_rate_local(204) 
  reac_source_local(46,204) = + reac_rate_local(204) 
  reac_source_local(01,205) = + reac_rate_local(205) 
  reac_source_local(02,205) = - reac_rate_local(205) 
  reac_source_local(46,205) = + reac_rate_local(205) 
  reac_source_local(01,206) = + reac_rate_local(206) 
  reac_source_local(07,206) = + reac_rate_local(206) 
  reac_source_local(13,206) = + reac_rate_local(206) 
  reac_source_local(33,206) = - reac_rate_local(206) 
  reac_source_local(01,207) = + reac_rate_local(207) 
  reac_source_local(02,207) = - reac_rate_local(207) 
  reac_source_local(07,207) = + reac_rate_local(207) 
  reac_source_local(13,207) = + reac_rate_local(207) 
  reac_source_local(01,208) = + reac_rate_local(208) 
  reac_source_local(08,208) = + reac_rate_local(208) 
  reac_source_local(16,208) = + reac_rate_local(208) 
  reac_source_local(33,208) = - reac_rate_local(208) 
  reac_source_local(01,209) = + reac_rate_local(209) 
  reac_source_local(02,209) = - reac_rate_local(209) 
  reac_source_local(08,209) = + reac_rate_local(209) 
  reac_source_local(16,209) = + reac_rate_local(209) 
  reac_source_local(01,210) = + reac_rate_local(210) 
  reac_source_local(28,210) = + reac_rate_local(210) 
  reac_source_local(33,210) = - reac_rate_local(210) 
  reac_source_local(50,210) = + reac_rate_local(210) 
  reac_source_local(01,211) = + reac_rate_local(211) 
  reac_source_local(02,211) = - reac_rate_local(211) 
  reac_source_local(28,211) = + reac_rate_local(211) 
  reac_source_local(50,211) = + reac_rate_local(211) 
  reac_source_local(01,212) = + reac_rate_local(212) 
  reac_source_local(23,212) = + reac_rate_local(212) 
  reac_source_local(33,212) = - reac_rate_local(212) 
  reac_source_local(44,212) = + reac_rate_local(212) 
  reac_source_local(01,213) = + reac_rate_local(213) 
  reac_source_local(02,213) = - reac_rate_local(213) 
  reac_source_local(23,213) = + reac_rate_local(213) 
  reac_source_local(44,213) = + reac_rate_local(213) 
  reac_source_local(01,214) = + reac_rate_local(214) 
  reac_source_local(17,214) = + reac_rate_local(214) 
  reac_source_local(52,214) = - reac_rate_local(214) 
  reac_source_local(01,215) = + reac_rate_local(215) 
  reac_source_local(09,215) = - reac_rate_local(215) 
  reac_source_local(17,215) = + reac_rate_local(215) 
  reac_source_local(01,216) = + reac_rate_local(216) 
  reac_source_local(17,216) = + reac_rate_local(216) 
  reac_source_local(38,216) = - reac_rate_local(216) 
  reac_source_local(01,217) = + reac_rate_local(217) 
  reac_source_local(08,217) = + reac_rate_local(217) 
  reac_source_local(44,217) = + reac_rate_local(217) 
  reac_source_local(52,217) = - reac_rate_local(217) 
  reac_source_local(01,218) = + reac_rate_local(218) 
  reac_source_local(08,218) = + reac_rate_local(218) 
  reac_source_local(09,218) = - reac_rate_local(218) 
  reac_source_local(44,218) = + reac_rate_local(218) 
  reac_source_local(01,219) = + reac_rate_local(219) 
  reac_source_local(08,219) = + reac_rate_local(219) 
  reac_source_local(38,219) = - reac_rate_local(219) 
  reac_source_local(44,219) = + reac_rate_local(219) 
  reac_source_local(01,220) = + reac_rate_local(220) 
  reac_source_local(14,220) = + reac_rate_local(220) 
  reac_source_local(32,220) = - reac_rate_local(220) 
  reac_source_local(01,221) = + reac_rate_local(221) 
  reac_source_local(14,221) = + reac_rate_local(221) 
  reac_source_local(24,221) = - reac_rate_local(221) 
  reac_source_local(01,222) = + reac_rate_local(222) 
  reac_source_local(03,222) = + reac_rate_local(222) 
  reac_source_local(13,222) = + reac_rate_local(222) 
  reac_source_local(32,222) = - reac_rate_local(222) 
  reac_source_local(01,223) = + reac_rate_local(223) 
  reac_source_local(03,223) = + reac_rate_local(223) 
  reac_source_local(13,223) = + reac_rate_local(223) 
  reac_source_local(24,223) = - reac_rate_local(223) 
  reac_source_local(01,224) = + reac_rate_local(224) 
  reac_source_local(05,224) = + reac_rate_local(224) 
  reac_source_local(11,224) = + reac_rate_local(224) 
  reac_source_local(32,224) = - reac_rate_local(224) 
  reac_source_local(01,225) = + reac_rate_local(225) 
  reac_source_local(05,225) = + reac_rate_local(225) 
  reac_source_local(11,225) = + reac_rate_local(225) 
  reac_source_local(24,225) = - reac_rate_local(225) 
  reac_source_local(01,226) = + reac_rate_local(226) 
  reac_source_local(11,226) = + reac_rate_local(226) 
  reac_source_local(13,226) = + reac_rate_local(226) 
  reac_source_local(26,226) = + reac_rate_local(226) 
  reac_source_local(32,226) = - reac_rate_local(226) 
  reac_source_local(01,227) = + reac_rate_local(227) 
  reac_source_local(11,227) = + reac_rate_local(227) 
  reac_source_local(13,227) = + reac_rate_local(227) 
  reac_source_local(24,227) = - reac_rate_local(227) 
  reac_source_local(26,227) = + reac_rate_local(227) 
  reac_source_local(01,228) = + reac_rate_local(228) 
  reac_source_local(11,228) = + reac_rate_local(228) * 2.d0
  reac_source_local(19,228) = + reac_rate_local(228) 
  reac_source_local(32,228) = - reac_rate_local(228) 
  reac_source_local(01,229) = + reac_rate_local(229) 
  reac_source_local(11,229) = + reac_rate_local(229) * 2.d0
  reac_source_local(19,229) = + reac_rate_local(229) 
  reac_source_local(24,229) = - reac_rate_local(229) 
  reac_source_local(01,230) = + reac_rate_local(230) 
  reac_source_local(04,230) = + reac_rate_local(230) 
  reac_source_local(23,230) = + reac_rate_local(230) 
  reac_source_local(32,230) = - reac_rate_local(230) 
  reac_source_local(01,231) = + reac_rate_local(231) 
  reac_source_local(04,231) = + reac_rate_local(231) 
  reac_source_local(23,231) = + reac_rate_local(231) 
  reac_source_local(24,231) = - reac_rate_local(231) 
  reac_source_local(01,232) = + reac_rate_local(232) 
  reac_source_local(32,232) = - reac_rate_local(232) 
  reac_source_local(46,232) = + reac_rate_local(232) 
  reac_source_local(51,232) = + reac_rate_local(232) 
  reac_source_local(01,233) = + reac_rate_local(233) 
  reac_source_local(24,233) = - reac_rate_local(233) 
  reac_source_local(46,233) = + reac_rate_local(233) 
  reac_source_local(51,233) = + reac_rate_local(233) 
  reac_source_local(01,234) = + reac_rate_local(234) 
  reac_source_local(16,234) = + reac_rate_local(234) 
  reac_source_local(32,234) = - reac_rate_local(234) 
  reac_source_local(34,234) = + reac_rate_local(234) 
  reac_source_local(01,235) = + reac_rate_local(235) 
  reac_source_local(16,235) = + reac_rate_local(235) 
  reac_source_local(24,235) = - reac_rate_local(235) 
  reac_source_local(34,235) = + reac_rate_local(235) 
  reac_source_local(01,236) = + reac_rate_local(236) 
  reac_source_local(32,236) = - reac_rate_local(236) 
  reac_source_local(43,236) = + reac_rate_local(236) 
  reac_source_local(50,236) = + reac_rate_local(236) 
  reac_source_local(01,237) = + reac_rate_local(237) 
  reac_source_local(24,237) = - reac_rate_local(237) 
  reac_source_local(43,237) = + reac_rate_local(237) 
  reac_source_local(50,237) = + reac_rate_local(237) 
  reac_source_local(11,238) = - reac_rate_local(238) 
  reac_source_local(13,238) = + reac_rate_local(238) * 2.d0
  reac_source_local(01,239) = + reac_rate_local(239) 
  reac_source_local(11,239) = - reac_rate_local(239) 
  reac_source_local(27,239) = + reac_rate_local(239) 
  reac_source_local(01,240) = - reac_rate_local(240) 
  reac_source_local(10,240) = - reac_rate_local(240) 
  reac_source_local(51,240) = + reac_rate_local(240) 
  reac_source_local(01,241) = - reac_rate_local(241) 
  reac_source_local(16,241) = - reac_rate_local(241) 
  reac_source_local(23,241) = + reac_rate_local(241) 
  reac_source_local(01,242) = - reac_rate_local(242) 
  reac_source_local(28,242) = + reac_rate_local(242) 
  reac_source_local(50,242) = - reac_rate_local(242) 
  reac_source_local(01,243) = - reac_rate_local(243) 
  reac_source_local(08,243) = + reac_rate_local(243) 
  reac_source_local(44,243) = - reac_rate_local(243) 
  reac_source_local(01,244) = - reac_rate_local(244) 
  reac_source_local(40,244) = - reac_rate_local(244) 
  reac_source_local(43,244) = + reac_rate_local(244) 
  reac_source_local(01,245) = - reac_rate_local(245) 
  reac_source_local(04,245) = - reac_rate_local(245) 
  reac_source_local(34,245) = + reac_rate_local(245) 
  reac_source_local(01,246) = - reac_rate_local(246) 
  reac_source_local(42,246) = + reac_rate_local(246) 
  reac_source_local(46,246) = - reac_rate_local(246) 
  reac_source_local(01,247) = - reac_rate_local(247) 
  reac_source_local(07,247) = - reac_rate_local(247) 
  reac_source_local(29,247) = + reac_rate_local(247) 
  reac_source_local(01,248) = - reac_rate_local(248) 
  reac_source_local(17,248) = - reac_rate_local(248) 
  reac_source_local(21,248) = + reac_rate_local(248) 
  reac_source_local(01,249) = - reac_rate_local(249) 
  reac_source_local(14,249) = - reac_rate_local(249) 
  reac_source_local(45,249) = + reac_rate_local(249) 
  reac_source_local(01,250) = - reac_rate_local(250) 
  reac_source_local(03,250) = - reac_rate_local(250) 
  reac_source_local(18,250) = + reac_rate_local(250) 
  reac_source_local(01,251) = - reac_rate_local(251) 
  reac_source_local(05,251) = - reac_rate_local(251) 
  reac_source_local(37,251) = + reac_rate_local(251) 
  reac_source_local(01,252) = - reac_rate_local(252) 
  reac_source_local(20,252) = + reac_rate_local(252) 
  reac_source_local(26,252) = - reac_rate_local(252) 
  reac_source_local(01,253) = - reac_rate_local(253) 
  reac_source_local(19,253) = - reac_rate_local(253) 
  reac_source_local(35,253) = + reac_rate_local(253) 
  reac_source_local(01,254) = - reac_rate_local(254) 
  reac_source_local(13,254) = + reac_rate_local(254) 
  reac_source_local(30,254) = - reac_rate_local(254) 
  reac_source_local(01,255) = - reac_rate_local(255) 
  reac_source_local(11,255) = + reac_rate_local(255) 
  reac_source_local(27,255) = - reac_rate_local(255) 
  reac_source_local(01,256) = - reac_rate_local(256) 
  reac_source_local(13,256) = + reac_rate_local(256) * 2.d0
  reac_source_local(23,256) = + reac_rate_local(256) 
  reac_source_local(47,256) = - reac_rate_local(256) 
  reac_source_local(01,257) = - reac_rate_local(257) 
  reac_source_local(11,257) = + reac_rate_local(257) 
  reac_source_local(13,257) = + reac_rate_local(257) 
  reac_source_local(28,257) = + reac_rate_local(257) 
  reac_source_local(47,257) = - reac_rate_local(257) 
  reac_source_local(01,258) = - reac_rate_local(258) 
  reac_source_local(10,258) = - reac_rate_local(258) 
  reac_source_local(13,258) = + reac_rate_local(258) 
  reac_source_local(23,258) = + reac_rate_local(258) 
  reac_source_local(01,259) = - reac_rate_local(259) 
  reac_source_local(10,259) = - reac_rate_local(259) 
  reac_source_local(13,259) = + reac_rate_local(259) * 2.d0
  reac_source_local(28,259) = + reac_rate_local(259) 
  reac_source_local(01,260) = - reac_rate_local(260) 
  reac_source_local(08,260) = + reac_rate_local(260) 
  reac_source_local(10,260) = - reac_rate_local(260) 
  reac_source_local(11,260) = + reac_rate_local(260) 
  reac_source_local(13,260) = + reac_rate_local(260) 
  reac_source_local(01,261) = - reac_rate_local(261) 
  reac_source_local(13,261) = + reac_rate_local(261) 
  reac_source_local(16,261) = - reac_rate_local(261) 
  reac_source_local(28,261) = + reac_rate_local(261) 
  reac_source_local(01,262) = - reac_rate_local(262) 
  reac_source_local(08,262) = + reac_rate_local(262) 
  reac_source_local(11,262) = + reac_rate_local(262) 
  reac_source_local(16,262) = - reac_rate_local(262) 
  reac_source_local(01,263) = - reac_rate_local(263) 
  reac_source_local(08,263) = + reac_rate_local(263) 
  reac_source_local(13,263) = + reac_rate_local(263) 
  reac_source_local(50,263) = - reac_rate_local(263) 
  reac_source_local(01,264) = - reac_rate_local(264) 
  reac_source_local(13,264) = + reac_rate_local(264) 
  reac_source_local(34,264) = + reac_rate_local(264) 
  reac_source_local(40,264) = - reac_rate_local(264) 
  reac_source_local(01,265) = - reac_rate_local(265) 
  reac_source_local(13,265) = + reac_rate_local(265) * 2.d0
  reac_source_local(40,265) = - reac_rate_local(265) 
  reac_source_local(42,265) = + reac_rate_local(265) 
  reac_source_local(01,266) = - reac_rate_local(266) 
  reac_source_local(04,266) = - reac_rate_local(266) 
  reac_source_local(13,266) = + reac_rate_local(266) 
  reac_source_local(42,266) = + reac_rate_local(266) 
  reac_source_local(01,267) = - reac_rate_local(267) 
  reac_source_local(04,267) = - reac_rate_local(267) 
  reac_source_local(13,267) = + reac_rate_local(267) * 2.d0
  reac_source_local(29,267) = + reac_rate_local(267) 
  reac_source_local(01,268) = - reac_rate_local(268) 
  reac_source_local(04,268) = - reac_rate_local(268) 
  reac_source_local(11,268) = + reac_rate_local(268) 
  reac_source_local(13,268) = + reac_rate_local(268) 
  reac_source_local(21,268) = + reac_rate_local(268) 
  reac_source_local(01,269) = - reac_rate_local(269) 
  reac_source_local(04,269) = - reac_rate_local(269) 
  reac_source_local(13,269) = + reac_rate_local(269) * 3.d0
  reac_source_local(21,269) = + reac_rate_local(269) 
  reac_source_local(01,270) = - reac_rate_local(270) 
  reac_source_local(04,270) = - reac_rate_local(270) 
  reac_source_local(23,270) = + reac_rate_local(270) 
  reac_source_local(28,270) = + reac_rate_local(270) 
  reac_source_local(01,271) = - reac_rate_local(271) 
  reac_source_local(13,271) = + reac_rate_local(271) 
  reac_source_local(29,271) = + reac_rate_local(271) 
  reac_source_local(46,271) = - reac_rate_local(271) 
  reac_source_local(01,272) = - reac_rate_local(272) 
  reac_source_local(13,272) = + reac_rate_local(272) * 2.d0
  reac_source_local(21,272) = + reac_rate_local(272) 
  reac_source_local(46,272) = - reac_rate_local(272) 
  reac_source_local(01,273) = - reac_rate_local(273) 
  reac_source_local(07,273) = - reac_rate_local(273) 
  reac_source_local(13,273) = + reac_rate_local(273) 
  reac_source_local(21,273) = + reac_rate_local(273) 
  reac_source_local(01,274) = - reac_rate_local(274) 
  reac_source_local(08,274) = + reac_rate_local(274) * 2.d0
  reac_source_local(17,274) = - reac_rate_local(274) 
  reac_source_local(01,275) = - reac_rate_local(275) 
  reac_source_local(13,275) = + reac_rate_local(275) * 2.d0
  reac_source_local(27,275) = - reac_rate_local(275) 
  reac_source_local(01,276) = - reac_rate_local(276) 
  reac_source_local(11,276) = + reac_rate_local(276) 
  reac_source_local(12,276) = - reac_rate_local(276) 
  reac_source_local(13,276) = + reac_rate_local(276) 
  reac_source_local(11,277) = - reac_rate_local(277) 
  reac_source_local(12,277) = + reac_rate_local(277) 
  reac_source_local(13,277) = + reac_rate_local(277) 
  reac_source_local(27,277) = - reac_rate_local(277) 
  reac_source_local(16,278) = + reac_rate_local(278) 
  reac_source_local(28,278) = - reac_rate_local(278) 
  reac_source_local(47,278) = - reac_rate_local(278) 
  reac_source_local(51,278) = + reac_rate_local(278) 
  reac_source_local(08,279) = - reac_rate_local(279) 
  reac_source_local(47,279) = - reac_rate_local(279) 
  reac_source_local(50,279) = + reac_rate_local(279) 
  reac_source_local(51,279) = + reac_rate_local(279) 
  reac_source_local(04,280) = + reac_rate_local(280) 
  reac_source_local(11,280) = + reac_rate_local(280) 
  reac_source_local(43,280) = - reac_rate_local(280) 
  reac_source_local(47,280) = - reac_rate_local(280) 
  reac_source_local(51,280) = + reac_rate_local(280) 
  reac_source_local(04,281) = + reac_rate_local(281) 
  reac_source_local(42,281) = - reac_rate_local(281) 
  reac_source_local(47,281) = - reac_rate_local(281) 
  reac_source_local(51,281) = + reac_rate_local(281) 
  reac_source_local(07,282) = + reac_rate_local(282) 
  reac_source_local(21,282) = - reac_rate_local(282) 
  reac_source_local(47,282) = - reac_rate_local(282) 
  reac_source_local(51,282) = + reac_rate_local(282) 
  reac_source_local(10,283) = + reac_rate_local(283) 
  reac_source_local(11,283) = + reac_rate_local(283) 
  reac_source_local(13,283) = - reac_rate_local(283) 
  reac_source_local(47,283) = - reac_rate_local(283) 
  reac_source_local(10,284) = - reac_rate_local(284) 
  reac_source_local(23,284) = + reac_rate_local(284) 
  reac_source_local(47,284) = + reac_rate_local(284) 
  reac_source_local(51,284) = - reac_rate_local(284) 
  reac_source_local(10,285) = - reac_rate_local(285) 
  reac_source_local(11,285) = + reac_rate_local(285) 
  reac_source_local(43,285) = - reac_rate_local(285) 
  reac_source_local(46,285) = + reac_rate_local(285) 
  reac_source_local(51,285) = + reac_rate_local(285) 
  reac_source_local(04,286) = + reac_rate_local(286) 
  reac_source_local(10,286) = - reac_rate_local(286) 
  reac_source_local(23,286) = + reac_rate_local(286) 
  reac_source_local(42,286) = - reac_rate_local(286) 
  reac_source_local(10,287) = - reac_rate_local(287) 
  reac_source_local(42,287) = - reac_rate_local(287) 
  reac_source_local(46,287) = + reac_rate_local(287) 
  reac_source_local(51,287) = + reac_rate_local(287) 
  reac_source_local(07,288) = + reac_rate_local(288) 
  reac_source_local(10,288) = - reac_rate_local(288) 
  reac_source_local(21,288) = - reac_rate_local(288) 
  reac_source_local(23,288) = + reac_rate_local(288) 
  reac_source_local(10,289) = - reac_rate_local(289) 
  reac_source_local(17,289) = + reac_rate_local(289) 
  reac_source_local(21,289) = - reac_rate_local(289) 
  reac_source_local(51,289) = + reac_rate_local(289) 
  reac_source_local(10,290) = - reac_rate_local(290) 
  reac_source_local(11,290) = - reac_rate_local(290) 
  reac_source_local(13,290) = + reac_rate_local(290) 
  reac_source_local(47,290) = + reac_rate_local(290) 
  reac_source_local(10,291) = - reac_rate_local(291) 
  reac_source_local(11,291) = + reac_rate_local(291) 
  reac_source_local(13,291) = - reac_rate_local(291) 
  reac_source_local(16,291) = + reac_rate_local(291) 
  reac_source_local(10,292) = + reac_rate_local(292) 
  reac_source_local(16,292) = - reac_rate_local(292) 
  reac_source_local(23,292) = + reac_rate_local(292) 
  reac_source_local(51,292) = - reac_rate_local(292) 
  reac_source_local(04,293) = + reac_rate_local(293) 
  reac_source_local(11,293) = + reac_rate_local(293) 
  reac_source_local(16,293) = - reac_rate_local(293) 
  reac_source_local(51,293) = - reac_rate_local(293) 
  reac_source_local(07,294) = + reac_rate_local(294) 
  reac_source_local(11,294) = + reac_rate_local(294) 
  reac_source_local(16,294) = - reac_rate_local(294) 
  reac_source_local(28,294) = - reac_rate_local(294) 
  reac_source_local(08,295) = - reac_rate_local(295) 
  reac_source_local(11,295) = + reac_rate_local(295) 
  reac_source_local(16,295) = - reac_rate_local(295) 
  reac_source_local(17,295) = + reac_rate_local(295) 
  reac_source_local(04,296) = + reac_rate_local(296) 
  reac_source_local(16,296) = - reac_rate_local(296) 
  reac_source_local(43,296) = - reac_rate_local(296) 
  reac_source_local(51,296) = + reac_rate_local(296) 
  reac_source_local(07,297) = + reac_rate_local(297) 
  reac_source_local(16,297) = - reac_rate_local(297) 
  reac_source_local(42,297) = - reac_rate_local(297) 
  reac_source_local(51,297) = + reac_rate_local(297) 
  reac_source_local(07,298) = + reac_rate_local(298) 
  reac_source_local(16,298) = - reac_rate_local(298) 
  reac_source_local(23,298) = + reac_rate_local(298) 
  reac_source_local(29,298) = - reac_rate_local(298) 
  reac_source_local(16,299) = + reac_rate_local(299) 
  reac_source_local(23,299) = + reac_rate_local(299) 
  reac_source_local(50,299) = - reac_rate_local(299) 
  reac_source_local(51,299) = - reac_rate_local(299) 
  reac_source_local(04,300) = + reac_rate_local(300) 
  reac_source_local(13,300) = + reac_rate_local(300) 
  reac_source_local(50,300) = - reac_rate_local(300) 
  reac_source_local(51,300) = - reac_rate_local(300) 
  reac_source_local(11,301) = + reac_rate_local(301) 
  reac_source_local(46,301) = + reac_rate_local(301) 
  reac_source_local(50,301) = - reac_rate_local(301) 
  reac_source_local(51,301) = - reac_rate_local(301) 
  reac_source_local(07,302) = + reac_rate_local(302) 
  reac_source_local(11,302) = + reac_rate_local(302) 
  reac_source_local(13,302) = + reac_rate_local(302) 
  reac_source_local(50,302) = - reac_rate_local(302) 
  reac_source_local(51,302) = - reac_rate_local(302) 
  reac_source_local(11,303) = + reac_rate_local(303) * 2.d0
  reac_source_local(17,303) = + reac_rate_local(303) 
  reac_source_local(50,303) = - reac_rate_local(303) 
  reac_source_local(51,303) = - reac_rate_local(303) 
  reac_source_local(11,304) = - reac_rate_local(304) 
  reac_source_local(13,304) = + reac_rate_local(304) 
  reac_source_local(16,304) = + reac_rate_local(304) 
  reac_source_local(50,304) = - reac_rate_local(304) 
  reac_source_local(13,305) = + reac_rate_local(305) 
  reac_source_local(44,305) = - reac_rate_local(305) 
  reac_source_local(46,305) = + reac_rate_local(305) 
  reac_source_local(51,305) = - reac_rate_local(305) 
  reac_source_local(07,306) = + reac_rate_local(306) 
  reac_source_local(11,306) = + reac_rate_local(306) 
  reac_source_local(44,306) = - reac_rate_local(306) 
  reac_source_local(51,306) = - reac_rate_local(306) 
  reac_source_local(11,307) = + reac_rate_local(307) 
  reac_source_local(13,307) = + reac_rate_local(307) 
  reac_source_local(17,307) = + reac_rate_local(307) 
  reac_source_local(44,307) = - reac_rate_local(307) 
  reac_source_local(51,307) = - reac_rate_local(307) 
  reac_source_local(11,308) = - reac_rate_local(308) 
  reac_source_local(13,308) = + reac_rate_local(308) 
  reac_source_local(44,308) = - reac_rate_local(308) 
  reac_source_local(50,308) = + reac_rate_local(308) 
  reac_source_local(40,309) = - reac_rate_local(309) 
  reac_source_local(42,309) = - reac_rate_local(309) 
  reac_source_local(43,309) = + reac_rate_local(309) 
  reac_source_local(46,309) = + reac_rate_local(309) 
  reac_source_local(04,310) = + reac_rate_local(310) 
  reac_source_local(21,310) = - reac_rate_local(310) 
  reac_source_local(29,310) = + reac_rate_local(310) 
  reac_source_local(40,310) = - reac_rate_local(310) 
  reac_source_local(04,311) = + reac_rate_local(311) 
  reac_source_local(11,311) = + reac_rate_local(311) 
  reac_source_local(13,311) = - reac_rate_local(311) 
  reac_source_local(40,311) = - reac_rate_local(311) 
  reac_source_local(04,312) = - reac_rate_local(312) 
  reac_source_local(11,312) = + reac_rate_local(312) 
  reac_source_local(13,312) = - reac_rate_local(312) 
  reac_source_local(46,312) = + reac_rate_local(312) 
  reac_source_local(04,313) = + reac_rate_local(313) 
  reac_source_local(21,313) = + reac_rate_local(313) 
  reac_source_local(29,313) = - reac_rate_local(313) 
  reac_source_local(46,313) = - reac_rate_local(313) 
  reac_source_local(07,314) = + reac_rate_local(314) 
  reac_source_local(29,314) = - reac_rate_local(314) 
  reac_source_local(42,314) = + reac_rate_local(314) 
  reac_source_local(46,314) = - reac_rate_local(314) 
  reac_source_local(07,315) = + reac_rate_local(315) 
  reac_source_local(11,315) = + reac_rate_local(315) 
  reac_source_local(13,315) = - reac_rate_local(315) 
  reac_source_local(46,315) = - reac_rate_local(315) 
  reac_source_local(04,316) = + reac_rate_local(316) 
  reac_source_local(07,316) = - reac_rate_local(316) 
  reac_source_local(42,316) = + reac_rate_local(316) 
  reac_source_local(43,316) = - reac_rate_local(316) 
  reac_source_local(04,317) = + reac_rate_local(317) 
  reac_source_local(07,317) = - reac_rate_local(317) 
  reac_source_local(21,317) = + reac_rate_local(317) 
  reac_source_local(42,317) = - reac_rate_local(317) 
  reac_source_local(06,318) = - reac_rate_local(318) 
  reac_source_local(07,318) = - reac_rate_local(318) 
  reac_source_local(17,318) = + reac_rate_local(318) 
  reac_source_local(21,318) = + reac_rate_local(318) 
  reac_source_local(07,319) = - reac_rate_local(319) 
  reac_source_local(11,319) = + reac_rate_local(319) 
  reac_source_local(13,319) = - reac_rate_local(319) 
  reac_source_local(17,319) = + reac_rate_local(319) 
  reac_source_local(07,320) = + reac_rate_local(320) 
  reac_source_local(17,320) = - reac_rate_local(320) 
  reac_source_local(23,320) = + reac_rate_local(320) 
  reac_source_local(51,320) = - reac_rate_local(320) 
  reac_source_local(04,321) = + reac_rate_local(321) 
  reac_source_local(17,321) = - reac_rate_local(321) 
  reac_source_local(29,321) = + reac_rate_local(321) 
  reac_source_local(43,321) = - reac_rate_local(321) 
  reac_source_local(17,322) = - reac_rate_local(322) 
  reac_source_local(42,322) = + reac_rate_local(322) 
  reac_source_local(43,322) = - reac_rate_local(322) 
  reac_source_local(46,322) = + reac_rate_local(322) 
  reac_source_local(17,323) = - reac_rate_local(323) 
  reac_source_local(21,323) = + reac_rate_local(323) 
  reac_source_local(42,323) = - reac_rate_local(323) 
  reac_source_local(46,323) = + reac_rate_local(323) 
  reac_source_local(07,324) = + reac_rate_local(324) 
  reac_source_local(17,324) = - reac_rate_local(324) 
  reac_source_local(21,324) = + reac_rate_local(324) 
  reac_source_local(29,324) = - reac_rate_local(324) 
  reac_source_local(07,325) = + reac_rate_local(325) 
  reac_source_local(11,325) = - reac_rate_local(325) 
  reac_source_local(13,325) = + reac_rate_local(325) 
  reac_source_local(17,325) = - reac_rate_local(325) 
  reac_source_local(17,326) = + reac_rate_local(326) 
  reac_source_local(23,326) = + reac_rate_local(326) 
  reac_source_local(36,326) = - reac_rate_local(326) 
  reac_source_local(51,326) = - reac_rate_local(326) 
  reac_source_local(11,327) = + reac_rate_local(327) 
  reac_source_local(12,327) = - reac_rate_local(327) 
  reac_source_local(47,327) = + reac_rate_local(327) 
  reac_source_local(51,327) = - reac_rate_local(327) 
  reac_source_local(10,328) = + reac_rate_local(328) 
  reac_source_local(11,328) = + reac_rate_local(328) 
  reac_source_local(12,328) = - reac_rate_local(328) 
  reac_source_local(23,328) = - reac_rate_local(328) 
  reac_source_local(11,329) = + reac_rate_local(329) 
  reac_source_local(12,329) = - reac_rate_local(329) 
  reac_source_local(16,329) = + reac_rate_local(329) 
  reac_source_local(28,329) = - reac_rate_local(329) 
  reac_source_local(08,330) = - reac_rate_local(330) 
  reac_source_local(11,330) = + reac_rate_local(330) 
  reac_source_local(12,330) = - reac_rate_local(330) 
  reac_source_local(50,330) = + reac_rate_local(330) 
  reac_source_local(04,331) = + reac_rate_local(331) 
  reac_source_local(11,331) = + reac_rate_local(331) * 2.d0
  reac_source_local(12,331) = - reac_rate_local(331) 
  reac_source_local(43,331) = - reac_rate_local(331) 
  reac_source_local(11,332) = + reac_rate_local(332) 
  reac_source_local(12,332) = - reac_rate_local(332) 
  reac_source_local(34,332) = - reac_rate_local(332) 
  reac_source_local(40,332) = + reac_rate_local(332) 
  reac_source_local(04,333) = + reac_rate_local(333) 
  reac_source_local(11,333) = + reac_rate_local(333) 
  reac_source_local(12,333) = - reac_rate_local(333) 
  reac_source_local(42,333) = - reac_rate_local(333) 
  reac_source_local(07,334) = + reac_rate_local(334) 
  reac_source_local(11,334) = + reac_rate_local(334) * 2.d0
  reac_source_local(12,334) = - reac_rate_local(334) 
  reac_source_local(42,334) = - reac_rate_local(334) 
  reac_source_local(11,335) = + reac_rate_local(335) 
  reac_source_local(12,335) = - reac_rate_local(335) 
  reac_source_local(29,335) = - reac_rate_local(335) 
  reac_source_local(46,335) = + reac_rate_local(335) 
  reac_source_local(06,336) = - reac_rate_local(336) 
  reac_source_local(11,336) = + reac_rate_local(336) 
  reac_source_local(12,336) = - reac_rate_local(336) 
  reac_source_local(17,336) = + reac_rate_local(336) 
  reac_source_local(07,337) = + reac_rate_local(337) 
  reac_source_local(11,337) = + reac_rate_local(337) 
  reac_source_local(12,337) = - reac_rate_local(337) 
  reac_source_local(21,337) = - reac_rate_local(337) 
  reac_source_local(13,338) = + reac_rate_local(338) 
  reac_source_local(27,338) = - reac_rate_local(338) 
  reac_source_local(47,338) = + reac_rate_local(338) 
  reac_source_local(51,338) = - reac_rate_local(338) 
  reac_source_local(10,339) = + reac_rate_local(339) 
  reac_source_local(11,339) = + reac_rate_local(339) 
  reac_source_local(27,339) = - reac_rate_local(339) 
  reac_source_local(51,339) = - reac_rate_local(339) 
  reac_source_local(11,340) = + reac_rate_local(340) 
  reac_source_local(13,340) = + reac_rate_local(340) 
  reac_source_local(16,340) = + reac_rate_local(340) 
  reac_source_local(27,340) = - reac_rate_local(340) 
  reac_source_local(51,340) = - reac_rate_local(340) 
  reac_source_local(13,341) = + reac_rate_local(341) 
  reac_source_local(16,341) = + reac_rate_local(341) 
  reac_source_local(27,341) = - reac_rate_local(341) 
  reac_source_local(28,341) = - reac_rate_local(341) 
  reac_source_local(11,342) = + reac_rate_local(342) 
  reac_source_local(27,342) = - reac_rate_local(342) 
  reac_source_local(28,342) = - reac_rate_local(342) 
  reac_source_local(50,342) = + reac_rate_local(342) 
  reac_source_local(08,343) = - reac_rate_local(343) 
  reac_source_local(13,343) = + reac_rate_local(343) 
  reac_source_local(27,343) = - reac_rate_local(343) 
  reac_source_local(50,343) = + reac_rate_local(343) 
  reac_source_local(08,344) = - reac_rate_local(344) 
  reac_source_local(11,344) = + reac_rate_local(344) 
  reac_source_local(27,344) = - reac_rate_local(344) 
  reac_source_local(44,344) = + reac_rate_local(344) 
  reac_source_local(11,345) = + reac_rate_local(345) 
  reac_source_local(27,345) = - reac_rate_local(345) 
  reac_source_local(40,345) = + reac_rate_local(345) 
  reac_source_local(43,345) = - reac_rate_local(345) 
  reac_source_local(04,346) = + reac_rate_local(346) 
  reac_source_local(11,346) = + reac_rate_local(346) 
  reac_source_local(13,346) = + reac_rate_local(346) 
  reac_source_local(27,346) = - reac_rate_local(346) 
  reac_source_local(43,346) = - reac_rate_local(346) 
  reac_source_local(11,347) = + reac_rate_local(347) * 2.d0
  reac_source_local(27,347) = - reac_rate_local(347) 
  reac_source_local(43,347) = - reac_rate_local(347) 
  reac_source_local(46,347) = + reac_rate_local(347) 
  reac_source_local(07,348) = + reac_rate_local(348) 
  reac_source_local(11,348) = + reac_rate_local(348) * 2.d0
  reac_source_local(13,348) = + reac_rate_local(348) 
  reac_source_local(27,348) = - reac_rate_local(348) 
  reac_source_local(43,348) = - reac_rate_local(348) 
  reac_source_local(11,349) = + reac_rate_local(349) * 3.d0
  reac_source_local(17,349) = + reac_rate_local(349) 
  reac_source_local(27,349) = - reac_rate_local(349) 
  reac_source_local(43,349) = - reac_rate_local(349) 
  reac_source_local(11,350) = + reac_rate_local(350) 
  reac_source_local(27,350) = - reac_rate_local(350) 
  reac_source_local(42,350) = - reac_rate_local(350) 
  reac_source_local(46,350) = + reac_rate_local(350) 
  reac_source_local(07,351) = + reac_rate_local(351) 
  reac_source_local(11,351) = + reac_rate_local(351) 
  reac_source_local(13,351) = + reac_rate_local(351) 
  reac_source_local(27,351) = - reac_rate_local(351) 
  reac_source_local(42,351) = - reac_rate_local(351) 
  reac_source_local(11,352) = + reac_rate_local(352) * 2.d0
  reac_source_local(17,352) = + reac_rate_local(352) 
  reac_source_local(27,352) = - reac_rate_local(352) 
  reac_source_local(42,352) = - reac_rate_local(352) 
  reac_source_local(07,353) = + reac_rate_local(353) 
  reac_source_local(13,353) = + reac_rate_local(353) 
  reac_source_local(21,353) = - reac_rate_local(353) 
  reac_source_local(27,353) = - reac_rate_local(353) 
  reac_source_local(11,354) = + reac_rate_local(354) 
  reac_source_local(17,354) = + reac_rate_local(354) 
  reac_source_local(21,354) = - reac_rate_local(354) 
  reac_source_local(27,354) = - reac_rate_local(354) 
  reac_source_local(12,355) = + reac_rate_local(355) 
  reac_source_local(13,355) = - reac_rate_local(355) 
  reac_source_local(27,355) = - reac_rate_local(355) 
  reac_source_local(11,356) = + reac_rate_local(356) 
  reac_source_local(13,356) = - reac_rate_local(356) 
  reac_source_local(27,356) = - reac_rate_local(356) 
  reac_source_local(30,356) = + reac_rate_local(356) 
  reac_source_local(10,357) = + reac_rate_local(357) 
  reac_source_local(13,357) = + reac_rate_local(357) 
  reac_source_local(30,357) = - reac_rate_local(357) 
  reac_source_local(51,357) = - reac_rate_local(357) 
  reac_source_local(11,358) = + reac_rate_local(358) 
  reac_source_local(16,358) = + reac_rate_local(358) 
  reac_source_local(30,358) = - reac_rate_local(358) 
  reac_source_local(51,358) = - reac_rate_local(358) 
  reac_source_local(13,359) = + reac_rate_local(359) 
  reac_source_local(16,359) = + reac_rate_local(359) 
  reac_source_local(23,359) = - reac_rate_local(359) 
  reac_source_local(30,359) = - reac_rate_local(359) 
  reac_source_local(13,360) = + reac_rate_local(360) 
  reac_source_local(28,360) = - reac_rate_local(360) 
  reac_source_local(30,360) = - reac_rate_local(360) 
  reac_source_local(50,360) = + reac_rate_local(360) 
  reac_source_local(11,361) = + reac_rate_local(361) 
  reac_source_local(28,361) = - reac_rate_local(361) 
  reac_source_local(30,361) = - reac_rate_local(361) 
  reac_source_local(44,361) = + reac_rate_local(361) 
  reac_source_local(08,362) = - reac_rate_local(362) 
  reac_source_local(13,362) = + reac_rate_local(362) 
  reac_source_local(30,362) = - reac_rate_local(362) 
  reac_source_local(44,362) = + reac_rate_local(362) 
  reac_source_local(04,363) = + reac_rate_local(363) 
  reac_source_local(11,363) = + reac_rate_local(363) 
  reac_source_local(30,363) = - reac_rate_local(363) 
  reac_source_local(43,363) = - reac_rate_local(363) 
  reac_source_local(11,364) = + reac_rate_local(364) 
  reac_source_local(13,364) = + reac_rate_local(364) 
  reac_source_local(30,364) = - reac_rate_local(364) 
  reac_source_local(43,364) = - reac_rate_local(364) 
  reac_source_local(46,364) = + reac_rate_local(364) 
  reac_source_local(07,365) = + reac_rate_local(365) 
  reac_source_local(11,365) = + reac_rate_local(365) * 2.d0
  reac_source_local(30,365) = - reac_rate_local(365) 
  reac_source_local(43,365) = - reac_rate_local(365) 
  reac_source_local(11,366) = + reac_rate_local(366) 
  reac_source_local(30,366) = - reac_rate_local(366) 
  reac_source_local(34,366) = - reac_rate_local(366) 
  reac_source_local(46,366) = + reac_rate_local(366) 
  reac_source_local(07,367) = + reac_rate_local(367) 
  reac_source_local(11,367) = + reac_rate_local(367) 
  reac_source_local(13,367) = + reac_rate_local(367) 
  reac_source_local(30,367) = - reac_rate_local(367) 
  reac_source_local(34,367) = - reac_rate_local(367) 
  reac_source_local(13,368) = + reac_rate_local(368) 
  reac_source_local(30,368) = - reac_rate_local(368) 
  reac_source_local(42,368) = - reac_rate_local(368) 
  reac_source_local(46,368) = + reac_rate_local(368) 
  reac_source_local(07,369) = + reac_rate_local(369) 
  reac_source_local(11,369) = + reac_rate_local(369) 
  reac_source_local(30,369) = - reac_rate_local(369) 
  reac_source_local(42,369) = - reac_rate_local(369) 
  reac_source_local(11,370) = + reac_rate_local(370) 
  reac_source_local(13,370) = + reac_rate_local(370) 
  reac_source_local(17,370) = + reac_rate_local(370) 
  reac_source_local(30,370) = - reac_rate_local(370) 
  reac_source_local(42,370) = - reac_rate_local(370) 
  reac_source_local(07,371) = + reac_rate_local(371) 
  reac_source_local(13,371) = + reac_rate_local(371) 
  reac_source_local(29,371) = - reac_rate_local(371) 
  reac_source_local(30,371) = - reac_rate_local(371) 
  reac_source_local(11,372) = + reac_rate_local(372) 
  reac_source_local(17,372) = + reac_rate_local(372) 
  reac_source_local(29,372) = - reac_rate_local(372) 
  reac_source_local(30,372) = - reac_rate_local(372) 
  reac_source_local(13,373) = + reac_rate_local(373) 
  reac_source_local(17,373) = + reac_rate_local(373) 
  reac_source_local(21,373) = - reac_rate_local(373) 
  reac_source_local(30,373) = - reac_rate_local(373) 
  reac_source_local(23,374) = + reac_rate_local(374) * 2.d0
  reac_source_local(28,374) = - reac_rate_local(374) 
  reac_source_local(51,374) = - reac_rate_local(374) 
  reac_source_local(08,375) = - reac_rate_local(375) 
  reac_source_local(13,375) = + reac_rate_local(375) 
  reac_source_local(42,375) = + reac_rate_local(375) 
  reac_source_local(51,375) = - reac_rate_local(375) 
  reac_source_local(23,376) = + reac_rate_local(376) 
  reac_source_local(34,376) = - reac_rate_local(376) 
  reac_source_local(43,376) = + reac_rate_local(376) 
  reac_source_local(51,376) = - reac_rate_local(376) 
  reac_source_local(23,377) = + reac_rate_local(377) 
  reac_source_local(29,377) = - reac_rate_local(377) 
  reac_source_local(42,377) = + reac_rate_local(377) 
  reac_source_local(51,377) = - reac_rate_local(377) 
  reac_source_local(06,378) = - reac_rate_local(378) 
  reac_source_local(21,378) = + reac_rate_local(378) 
  reac_source_local(23,378) = + reac_rate_local(378) 
  reac_source_local(51,378) = - reac_rate_local(378) 
  reac_source_local(18,379) = - reac_rate_local(379) 
  reac_source_local(23,379) = + reac_rate_local(379) 
  reac_source_local(45,379) = + reac_rate_local(379) 
  reac_source_local(51,379) = - reac_rate_local(379) 
  reac_source_local(20,380) = - reac_rate_local(380) 
  reac_source_local(23,380) = + reac_rate_local(380) 
  reac_source_local(37,380) = + reac_rate_local(380) 
  reac_source_local(51,380) = - reac_rate_local(380) 
  reac_source_local(11,381) = + reac_rate_local(381) 
  reac_source_local(13,381) = - reac_rate_local(381) 
  reac_source_local(23,381) = + reac_rate_local(381) 
  reac_source_local(51,381) = - reac_rate_local(381) 
  reac_source_local(13,382) = + reac_rate_local(382) 
  reac_source_local(23,382) = - reac_rate_local(382) 
  reac_source_local(43,382) = + reac_rate_local(382) 
  reac_source_local(51,382) = - reac_rate_local(382) 
  reac_source_local(23,383) = + reac_rate_local(383) 
  reac_source_local(31,383) = - reac_rate_local(383) 
  reac_source_local(49,383) = + reac_rate_local(383) 
  reac_source_local(51,383) = - reac_rate_local(383) 
  reac_source_local(28,384) = - reac_rate_local(384) 
  reac_source_local(43,384) = + reac_rate_local(384) 
  reac_source_local(51,384) = - reac_rate_local(384) 
  reac_source_local(13,385) = + reac_rate_local(385) 
  reac_source_local(23,385) = + reac_rate_local(385) 
  reac_source_local(51,385) = - reac_rate_local(385) 
  reac_source_local(13,386) = + reac_rate_local(386) 
  reac_source_local(23,386) = - reac_rate_local(386) 
  reac_source_local(28,386) = + reac_rate_local(386) 
  reac_source_local(08,387) = + reac_rate_local(387) 
  reac_source_local(11,387) = + reac_rate_local(387) 
  reac_source_local(23,387) = - reac_rate_local(387) 
  reac_source_local(23,388) = - reac_rate_local(388) 
  reac_source_local(28,388) = + reac_rate_local(388) 
  reac_source_local(34,388) = - reac_rate_local(388) 
  reac_source_local(43,388) = + reac_rate_local(388) 
  reac_source_local(13,389) = + reac_rate_local(389) * 2.d0
  reac_source_local(21,389) = + reac_rate_local(389) 
  reac_source_local(28,389) = - reac_rate_local(389) * 2.d0
  reac_source_local(23,390) = + reac_rate_local(390) 
  reac_source_local(28,390) = - reac_rate_local(390) 
  reac_source_local(34,390) = - reac_rate_local(390) 
  reac_source_local(42,390) = + reac_rate_local(390) 
  reac_source_local(21,391) = + reac_rate_local(391) 
  reac_source_local(23,391) = + reac_rate_local(391) 
  reac_source_local(28,391) = - reac_rate_local(391) 
  reac_source_local(29,391) = - reac_rate_local(391) 
  reac_source_local(06,392) = - reac_rate_local(392) 
  reac_source_local(08,392) = + reac_rate_local(392) 
  reac_source_local(21,392) = + reac_rate_local(392) 
  reac_source_local(28,392) = - reac_rate_local(392) 
  reac_source_local(18,393) = + reac_rate_local(393) 
  reac_source_local(23,393) = + reac_rate_local(393) 
  reac_source_local(28,393) = - reac_rate_local(393) 
  reac_source_local(45,393) = - reac_rate_local(393) 
  reac_source_local(18,394) = - reac_rate_local(394) 
  reac_source_local(28,394) = - reac_rate_local(394) 
  reac_source_local(34,394) = + reac_rate_local(394) 
  reac_source_local(42,394) = + reac_rate_local(394) 
  reac_source_local(18,395) = - reac_rate_local(395) 
  reac_source_local(23,395) = + reac_rate_local(395) 
  reac_source_local(28,395) = - reac_rate_local(395) 
  reac_source_local(37,395) = + reac_rate_local(395) 
  reac_source_local(20,396) = + reac_rate_local(396) 
  reac_source_local(23,396) = + reac_rate_local(396) 
  reac_source_local(28,396) = - reac_rate_local(396) 
  reac_source_local(37,396) = - reac_rate_local(396) 
  reac_source_local(11,397) = - reac_rate_local(397) 
  reac_source_local(13,397) = + reac_rate_local(397) 
  reac_source_local(23,397) = + reac_rate_local(397) 
  reac_source_local(28,397) = - reac_rate_local(397) 
  reac_source_local(08,398) = + reac_rate_local(398) 
  reac_source_local(11,398) = + reac_rate_local(398) 
  reac_source_local(13,398) = - reac_rate_local(398) 
  reac_source_local(28,398) = - reac_rate_local(398) 
  reac_source_local(08,399) = + reac_rate_local(399) 
  reac_source_local(13,399) = + reac_rate_local(399) 
  reac_source_local(28,399) = - reac_rate_local(399) 
  reac_source_local(11,400) = + reac_rate_local(400) 
  reac_source_local(21,400) = + reac_rate_local(400) 
  reac_source_local(28,400) = - reac_rate_local(400) * 2.d0
  reac_source_local(13,401) = - reac_rate_local(401) 
  reac_source_local(23,401) = + reac_rate_local(401) 
  reac_source_local(28,401) = - reac_rate_local(401) 
  reac_source_local(08,402) = - reac_rate_local(402) 
  reac_source_local(13,402) = + reac_rate_local(402) 
  reac_source_local(37,402) = + reac_rate_local(402) 
  reac_source_local(43,402) = - reac_rate_local(402) 
  reac_source_local(08,403) = - reac_rate_local(403) 
  reac_source_local(18,403) = + reac_rate_local(403) 
  reac_source_local(43,403) = - reac_rate_local(403) 
  reac_source_local(08,404) = - reac_rate_local(404) 
  reac_source_local(11,404) = - reac_rate_local(404) 
  reac_source_local(13,404) = + reac_rate_local(404) 
  reac_source_local(28,404) = + reac_rate_local(404) 
  reac_source_local(08,405) = - reac_rate_local(405) 
  reac_source_local(13,405) = + reac_rate_local(405) 
  reac_source_local(23,405) = - reac_rate_local(405) 
  reac_source_local(29,405) = + reac_rate_local(405) 
  reac_source_local(08,406) = - reac_rate_local(406) 
  reac_source_local(13,406) = + reac_rate_local(406) 
  reac_source_local(21,406) = + reac_rate_local(406) 
  reac_source_local(28,406) = - reac_rate_local(406) 
  reac_source_local(08,407) = - reac_rate_local(407) 
  reac_source_local(11,407) = - reac_rate_local(407) 
  reac_source_local(23,407) = + reac_rate_local(407) 
  reac_source_local(08,408) = - reac_rate_local(408) 
  reac_source_local(21,408) = + reac_rate_local(408) 
  reac_source_local(28,408) = + reac_rate_local(408) 
  reac_source_local(29,408) = - reac_rate_local(408) 
  reac_source_local(29,409) = - reac_rate_local(409) 
  reac_source_local(34,409) = + reac_rate_local(409) 
  reac_source_local(42,409) = + reac_rate_local(409) 
  reac_source_local(43,409) = - reac_rate_local(409) 
  reac_source_local(18,410) = - reac_rate_local(410) 
  reac_source_local(34,410) = + reac_rate_local(410) 
  reac_source_local(43,410) = - reac_rate_local(410) 
  reac_source_local(45,410) = + reac_rate_local(410) 
  reac_source_local(20,411) = - reac_rate_local(411) 
  reac_source_local(34,411) = + reac_rate_local(411) 
  reac_source_local(37,411) = + reac_rate_local(411) 
  reac_source_local(43,411) = - reac_rate_local(411) 
  reac_source_local(11,412) = + reac_rate_local(412) 
  reac_source_local(13,412) = - reac_rate_local(412) 
  reac_source_local(34,412) = + reac_rate_local(412) 
  reac_source_local(43,412) = - reac_rate_local(412) 
  reac_source_local(13,413) = - reac_rate_local(413) 
  reac_source_local(23,413) = + reac_rate_local(413) 
  reac_source_local(43,413) = - reac_rate_local(413) 
  reac_source_local(51,413) = + reac_rate_local(413) 
  reac_source_local(23,414) = + reac_rate_local(414) * 2.d0
  reac_source_local(43,414) = - reac_rate_local(414) 
  reac_source_local(31,415) = - reac_rate_local(415) 
  reac_source_local(34,415) = + reac_rate_local(415) 
  reac_source_local(43,415) = - reac_rate_local(415) 
  reac_source_local(49,415) = + reac_rate_local(415) 
  reac_source_local(08,416) = - reac_rate_local(416) 
  reac_source_local(23,416) = + reac_rate_local(416) 
  reac_source_local(42,416) = + reac_rate_local(416) 
  reac_source_local(43,416) = - reac_rate_local(416) 
  reac_source_local(23,417) = + reac_rate_local(417) 
  reac_source_local(28,417) = - reac_rate_local(417) 
  reac_source_local(34,417) = + reac_rate_local(417) 
  reac_source_local(43,417) = - reac_rate_local(417) 
  reac_source_local(34,418) = - reac_rate_local(418) * 2.d0
  reac_source_local(42,418) = + reac_rate_local(418) 
  reac_source_local(43,418) = + reac_rate_local(418) 
  reac_source_local(29,419) = + reac_rate_local(419) 
  reac_source_local(34,419) = - reac_rate_local(419) 
  reac_source_local(42,419) = - reac_rate_local(419) 
  reac_source_local(43,419) = + reac_rate_local(419) 
  reac_source_local(06,420) = + reac_rate_local(420) 
  reac_source_local(21,420) = - reac_rate_local(420) 
  reac_source_local(34,420) = - reac_rate_local(420) 
  reac_source_local(43,420) = + reac_rate_local(420) 
  reac_source_local(06,421) = - reac_rate_local(421) 
  reac_source_local(21,421) = + reac_rate_local(421) 
  reac_source_local(34,421) = - reac_rate_local(421) 
  reac_source_local(42,421) = + reac_rate_local(421) 
  reac_source_local(18,422) = + reac_rate_local(422) 
  reac_source_local(34,422) = - reac_rate_local(422) 
  reac_source_local(43,422) = + reac_rate_local(422) 
  reac_source_local(45,422) = - reac_rate_local(422) 
  reac_source_local(18,423) = - reac_rate_local(423) 
  reac_source_local(34,423) = - reac_rate_local(423) 
  reac_source_local(42,423) = + reac_rate_local(423) 
  reac_source_local(45,423) = + reac_rate_local(423) 
  reac_source_local(18,424) = - reac_rate_local(424) 
  reac_source_local(34,424) = - reac_rate_local(424) 
  reac_source_local(37,424) = + reac_rate_local(424) 
  reac_source_local(43,424) = + reac_rate_local(424) 
  reac_source_local(20,425) = + reac_rate_local(425) 
  reac_source_local(34,425) = - reac_rate_local(425) 
  reac_source_local(37,425) = - reac_rate_local(425) 
  reac_source_local(43,425) = + reac_rate_local(425) 
  reac_source_local(11,426) = - reac_rate_local(426) 
  reac_source_local(13,426) = + reac_rate_local(426) 
  reac_source_local(34,426) = - reac_rate_local(426) 
  reac_source_local(43,426) = + reac_rate_local(426) 
  reac_source_local(13,427) = - reac_rate_local(427) 
  reac_source_local(23,427) = + reac_rate_local(427) * 2.d0
  reac_source_local(34,427) = - reac_rate_local(427) 
  reac_source_local(11,428) = + reac_rate_local(428) 
  reac_source_local(13,428) = - reac_rate_local(428) 
  reac_source_local(34,428) = - reac_rate_local(428) 
  reac_source_local(42,428) = + reac_rate_local(428) 
  reac_source_local(13,429) = - reac_rate_local(429) 
  reac_source_local(34,429) = - reac_rate_local(429) 
  reac_source_local(43,429) = + reac_rate_local(429) 
  reac_source_local(13,430) = + reac_rate_local(430) 
  reac_source_local(34,430) = - reac_rate_local(430) 
  reac_source_local(42,430) = + reac_rate_local(430) 
  reac_source_local(34,431) = - reac_rate_local(431) * 2.d0
  reac_source_local(49,431) = + reac_rate_local(431) 
  reac_source_local(31,432) = - reac_rate_local(432) 
  reac_source_local(34,432) = - reac_rate_local(432) 
  reac_source_local(42,432) = + reac_rate_local(432) 
  reac_source_local(49,432) = + reac_rate_local(432) 
  reac_source_local(29,433) = - reac_rate_local(433) 
  reac_source_local(34,433) = - reac_rate_local(433) 
  reac_source_local(42,433) = + reac_rate_local(433) * 2.d0
  reac_source_local(11,434) = + reac_rate_local(434) 
  reac_source_local(13,434) = - reac_rate_local(434) 
  reac_source_local(29,434) = + reac_rate_local(434) 
  reac_source_local(42,434) = - reac_rate_local(434) 
  reac_source_local(13,435) = - reac_rate_local(435) 
  reac_source_local(34,435) = + reac_rate_local(435) 
  reac_source_local(42,435) = - reac_rate_local(435) 
  reac_source_local(11,436) = - reac_rate_local(436) 
  reac_source_local(13,436) = + reac_rate_local(436) 
  reac_source_local(34,436) = + reac_rate_local(436) 
  reac_source_local(42,436) = - reac_rate_local(436) 
  reac_source_local(13,437) = + reac_rate_local(437) 
  reac_source_local(29,437) = + reac_rate_local(437) 
  reac_source_local(42,437) = - reac_rate_local(437) 
  reac_source_local(20,438) = + reac_rate_local(438) 
  reac_source_local(34,438) = + reac_rate_local(438) 
  reac_source_local(37,438) = - reac_rate_local(438) 
  reac_source_local(42,438) = - reac_rate_local(438) 
  reac_source_local(21,439) = - reac_rate_local(439) 
  reac_source_local(29,439) = + reac_rate_local(439) * 2.d0
  reac_source_local(42,439) = - reac_rate_local(439) 
  reac_source_local(18,440) = + reac_rate_local(440) 
  reac_source_local(29,440) = + reac_rate_local(440) 
  reac_source_local(37,440) = - reac_rate_local(440) 
  reac_source_local(42,440) = - reac_rate_local(440) 
  reac_source_local(29,441) = + reac_rate_local(441) 
  reac_source_local(34,441) = + reac_rate_local(441) 
  reac_source_local(42,441) = - reac_rate_local(441) * 2.d0
  reac_source_local(18,442) = + reac_rate_local(442) 
  reac_source_local(23,442) = - reac_rate_local(442) 
  reac_source_local(42,442) = - reac_rate_local(442) 
  reac_source_local(11,443) = + reac_rate_local(443) 
  reac_source_local(21,443) = + reac_rate_local(443) 
  reac_source_local(42,443) = - reac_rate_local(443) 
  reac_source_local(31,444) = + reac_rate_local(444) 
  reac_source_local(34,444) = - reac_rate_local(444) 
  reac_source_local(42,444) = - reac_rate_local(444) 
  reac_source_local(11,445) = - reac_rate_local(445) 
  reac_source_local(42,445) = - reac_rate_local(445) 
  reac_source_local(43,445) = + reac_rate_local(445) 
  reac_source_local(28,446) = - reac_rate_local(446) 
  reac_source_local(37,446) = + reac_rate_local(446) 
  reac_source_local(42,446) = - reac_rate_local(446) 
  reac_source_local(18,447) = + reac_rate_local(447) 
  reac_source_local(31,447) = - reac_rate_local(447) 
  reac_source_local(37,447) = + reac_rate_local(447) 
  reac_source_local(42,447) = - reac_rate_local(447) 
  reac_source_local(21,448) = + reac_rate_local(448) 
  reac_source_local(29,448) = - reac_rate_local(448) * 2.d0
  reac_source_local(42,448) = + reac_rate_local(448) 
  reac_source_local(18,449) = + reac_rate_local(449) 
  reac_source_local(29,449) = - reac_rate_local(449) 
  reac_source_local(42,449) = + reac_rate_local(449) 
  reac_source_local(45,449) = - reac_rate_local(449) 
  reac_source_local(18,450) = - reac_rate_local(450) 
  reac_source_local(21,450) = + reac_rate_local(450) 
  reac_source_local(29,450) = - reac_rate_local(450) 
  reac_source_local(45,450) = + reac_rate_local(450) 
  reac_source_local(18,451) = - reac_rate_local(451) 
  reac_source_local(29,451) = - reac_rate_local(451) 
  reac_source_local(37,451) = + reac_rate_local(451) 
  reac_source_local(42,451) = + reac_rate_local(451) 
  reac_source_local(20,452) = + reac_rate_local(452) 
  reac_source_local(29,452) = - reac_rate_local(452) 
  reac_source_local(37,452) = - reac_rate_local(452) 
  reac_source_local(42,452) = + reac_rate_local(452) 
  reac_source_local(20,453) = - reac_rate_local(453) 
  reac_source_local(21,453) = + reac_rate_local(453) 
  reac_source_local(29,453) = - reac_rate_local(453) 
  reac_source_local(37,453) = + reac_rate_local(453) 
  reac_source_local(11,454) = - reac_rate_local(454) 
  reac_source_local(13,454) = + reac_rate_local(454) 
  reac_source_local(29,454) = - reac_rate_local(454) 
  reac_source_local(42,454) = + reac_rate_local(454) 
  reac_source_local(11,455) = + reac_rate_local(455) 
  reac_source_local(13,455) = - reac_rate_local(455) 
  reac_source_local(21,455) = + reac_rate_local(455) 
  reac_source_local(29,455) = - reac_rate_local(455) 
  reac_source_local(13,456) = - reac_rate_local(456) 
  reac_source_local(29,456) = - reac_rate_local(456) 
  reac_source_local(42,456) = + reac_rate_local(456) 
  reac_source_local(13,457) = + reac_rate_local(457) 
  reac_source_local(21,457) = + reac_rate_local(457) 
  reac_source_local(29,457) = - reac_rate_local(457) 
  reac_source_local(21,458) = + reac_rate_local(458) 
  reac_source_local(29,458) = - reac_rate_local(458) 
  reac_source_local(31,458) = - reac_rate_local(458) 
  reac_source_local(49,458) = + reac_rate_local(458) 
  reac_source_local(13,459) = - reac_rate_local(459) 
  reac_source_local(21,459) = - reac_rate_local(459) 
  reac_source_local(29,459) = + reac_rate_local(459) 
  reac_source_local(11,460) = - reac_rate_local(460) 
  reac_source_local(21,460) = - reac_rate_local(460) 
  reac_source_local(42,460) = + reac_rate_local(460) 
  reac_source_local(11,461) = - reac_rate_local(461) 
  reac_source_local(13,461) = + reac_rate_local(461) 
  reac_source_local(21,461) = - reac_rate_local(461) 
  reac_source_local(29,461) = + reac_rate_local(461) 
  reac_source_local(20,462) = + reac_rate_local(462) 
  reac_source_local(21,462) = - reac_rate_local(462) 
  reac_source_local(23,462) = - reac_rate_local(462) 
  reac_source_local(20,463) = + reac_rate_local(463) 
  reac_source_local(21,463) = - reac_rate_local(463) 
  reac_source_local(31,463) = - reac_rate_local(463) 
  reac_source_local(37,463) = + reac_rate_local(463) 
  reac_source_local(18,464) = + reac_rate_local(464) 
  reac_source_local(20,464) = - reac_rate_local(464) 
  reac_source_local(37,464) = + reac_rate_local(464) 
  reac_source_local(45,464) = - reac_rate_local(464) 
  reac_source_local(11,465) = + reac_rate_local(465) 
  reac_source_local(13,465) = - reac_rate_local(465) 
  reac_source_local(18,465) = + reac_rate_local(465) 
  reac_source_local(45,465) = - reac_rate_local(465) 
  reac_source_local(23,466) = + reac_rate_local(466) 
  reac_source_local(34,466) = + reac_rate_local(466) 
  reac_source_local(45,466) = - reac_rate_local(466) 
  reac_source_local(18,467) = + reac_rate_local(467) 
  reac_source_local(31,467) = - reac_rate_local(467) 
  reac_source_local(45,467) = - reac_rate_local(467) 
  reac_source_local(49,467) = + reac_rate_local(467) 
  reac_source_local(28,468) = - reac_rate_local(468) 
  reac_source_local(45,468) = - reac_rate_local(468) 
  reac_source_local(49,468) = + reac_rate_local(468) 
  reac_source_local(18,469) = - reac_rate_local(469) * 2.d0
  reac_source_local(37,469) = + reac_rate_local(469) 
  reac_source_local(45,469) = + reac_rate_local(469) 
  reac_source_local(18,470) = - reac_rate_local(470) 
  reac_source_local(20,470) = + reac_rate_local(470) 
  reac_source_local(37,470) = - reac_rate_local(470) 
  reac_source_local(45,470) = + reac_rate_local(470) 
  reac_source_local(18,471) = - reac_rate_local(471) 
  reac_source_local(20,471) = - reac_rate_local(471) 
  reac_source_local(37,471) = + reac_rate_local(471) * 2.d0
  reac_source_local(11,472) = - reac_rate_local(472) 
  reac_source_local(13,472) = + reac_rate_local(472) 
  reac_source_local(18,472) = - reac_rate_local(472) 
  reac_source_local(45,472) = + reac_rate_local(472) 
  reac_source_local(11,473) = + reac_rate_local(473) 
  reac_source_local(13,473) = - reac_rate_local(473) 
  reac_source_local(18,473) = - reac_rate_local(473) 
  reac_source_local(37,473) = + reac_rate_local(473) 
  reac_source_local(13,474) = - reac_rate_local(474) 
  reac_source_local(18,474) = - reac_rate_local(474) 
  reac_source_local(45,474) = + reac_rate_local(474) 
  reac_source_local(13,475) = - reac_rate_local(475) 
  reac_source_local(18,475) = - reac_rate_local(475) 
  reac_source_local(23,475) = + reac_rate_local(475) 
  reac_source_local(34,475) = + reac_rate_local(475) 
  reac_source_local(13,476) = + reac_rate_local(476) 
  reac_source_local(18,476) = - reac_rate_local(476) 
  reac_source_local(37,476) = + reac_rate_local(476) 
  reac_source_local(18,477) = - reac_rate_local(477) 
  reac_source_local(23,477) = + reac_rate_local(477) 
  reac_source_local(42,477) = + reac_rate_local(477) 
  reac_source_local(18,478) = - reac_rate_local(478) 
  reac_source_local(31,478) = - reac_rate_local(478) 
  reac_source_local(37,478) = + reac_rate_local(478) 
  reac_source_local(49,478) = + reac_rate_local(478) 
  reac_source_local(20,479) = + reac_rate_local(479) 
  reac_source_local(21,479) = - reac_rate_local(479) 
  reac_source_local(29,479) = + reac_rate_local(479) 
  reac_source_local(37,479) = - reac_rate_local(479) 
  reac_source_local(18,480) = + reac_rate_local(480) 
  reac_source_local(20,480) = + reac_rate_local(480) 
  reac_source_local(37,480) = - reac_rate_local(480) * 2.d0
  reac_source_local(13,481) = + reac_rate_local(481) 
  reac_source_local(20,481) = + reac_rate_local(481) 
  reac_source_local(37,481) = - reac_rate_local(481) 
  reac_source_local(11,482) = + reac_rate_local(482) 
  reac_source_local(13,482) = - reac_rate_local(482) 
  reac_source_local(20,482) = + reac_rate_local(482) 
  reac_source_local(37,482) = - reac_rate_local(482) 
  reac_source_local(13,483) = - reac_rate_local(483) 
  reac_source_local(18,483) = + reac_rate_local(483) 
  reac_source_local(37,483) = - reac_rate_local(483) 
  reac_source_local(23,484) = + reac_rate_local(484) 
  reac_source_local(29,484) = + reac_rate_local(484) 
  reac_source_local(37,484) = - reac_rate_local(484) 
  reac_source_local(23,485) = - reac_rate_local(485) 
  reac_source_local(31,485) = + reac_rate_local(485) 
  reac_source_local(37,485) = - reac_rate_local(485) 
  reac_source_local(20,486) = + reac_rate_local(486) 
  reac_source_local(31,486) = - reac_rate_local(486) 
  reac_source_local(37,486) = - reac_rate_local(486) 
  reac_source_local(49,486) = + reac_rate_local(486) 
  reac_source_local(11,487) = - reac_rate_local(487) 
  reac_source_local(13,487) = + reac_rate_local(487) 
  reac_source_local(20,487) = - reac_rate_local(487) 
  reac_source_local(37,487) = + reac_rate_local(487) 
  reac_source_local(13,488) = - reac_rate_local(488) 
  reac_source_local(20,488) = - reac_rate_local(488) 
  reac_source_local(37,488) = + reac_rate_local(488) 
  reac_source_local(20,489) = - reac_rate_local(489) 
  reac_source_local(21,489) = + reac_rate_local(489) 
  reac_source_local(23,489) = + reac_rate_local(489) 
  reac_source_local(31,490) = - reac_rate_local(490) 
  reac_source_local(34,490) = + reac_rate_local(490) 
  reac_source_local(42,490) = + reac_rate_local(490) 
  reac_source_local(18,491) = + reac_rate_local(491) 
  reac_source_local(28,491) = - reac_rate_local(491) 
  reac_source_local(31,491) = - reac_rate_local(491) 
  reac_source_local(42,491) = + reac_rate_local(491) 
  reac_source_local(23,492) = + reac_rate_local(492) 
  reac_source_local(31,492) = - reac_rate_local(492) 
  reac_source_local(37,492) = + reac_rate_local(492) 
  reac_source_local(11,493) = - reac_rate_local(493) 
  reac_source_local(13,493) = + reac_rate_local(493) 
  reac_source_local(31,493) = - reac_rate_local(493) 
  reac_source_local(49,493) = + reac_rate_local(493) 
  reac_source_local(23,494) = - reac_rate_local(494) 
  reac_source_local(31,494) = + reac_rate_local(494) 
  reac_source_local(49,494) = - reac_rate_local(494) 
  reac_source_local(51,494) = + reac_rate_local(494) 
  reac_source_local(18,495) = + reac_rate_local(495) 
  reac_source_local(23,495) = + reac_rate_local(495) 
  reac_source_local(49,495) = - reac_rate_local(495) 
  reac_source_local(34,496) = + reac_rate_local(496) * 2.d0
  reac_source_local(49,496) = - reac_rate_local(496) 
  reac_source_local(11,497) = + reac_rate_local(497) 
  reac_source_local(13,497) = - reac_rate_local(497) 
  reac_source_local(31,497) = + reac_rate_local(497) 
  reac_source_local(49,497) = - reac_rate_local(497) 
  reac_source_local(23,498) = + reac_rate_local(498) 
  reac_source_local(28,498) = - reac_rate_local(498) 
  reac_source_local(31,498) = + reac_rate_local(498) 
  reac_source_local(49,498) = - reac_rate_local(498) 
  reac_source_local(29,499) = - reac_rate_local(499) 
  reac_source_local(31,499) = + reac_rate_local(499) 
  reac_source_local(42,499) = + reac_rate_local(499) 
  reac_source_local(49,499) = - reac_rate_local(499) 
  reac_source_local(18,500) = - reac_rate_local(500) 
  reac_source_local(31,500) = + reac_rate_local(500) 
  reac_source_local(45,500) = + reac_rate_local(500) 
  reac_source_local(49,500) = - reac_rate_local(500) 
  reac_source_local(06,501) = - reac_rate_local(501) 
  reac_source_local(21,501) = + reac_rate_local(501) 
  reac_source_local(31,501) = + reac_rate_local(501) 
  reac_source_local(49,501) = - reac_rate_local(501) 
  reac_source_local(31,502) = + reac_rate_local(502) 
  reac_source_local(34,502) = - reac_rate_local(502) 
  reac_source_local(43,502) = + reac_rate_local(502) 
  reac_source_local(49,502) = - reac_rate_local(502) 
  reac_source_local(20,503) = - reac_rate_local(503) 
  reac_source_local(31,503) = + reac_rate_local(503) 
  reac_source_local(37,503) = + reac_rate_local(503) 
  reac_source_local(49,503) = - reac_rate_local(503) 
  reac_source_local(28,504) = - reac_rate_local(504) 
  reac_source_local(48,504) = + reac_rate_local(504) 
  reac_source_local(49,504) = - reac_rate_local(504) 
  reac_source_local(23,505) = + reac_rate_local(505) 
  reac_source_local(31,505) = + reac_rate_local(505) 
  reac_source_local(48,505) = - reac_rate_local(505) 
  reac_source_local(11,506) = - reac_rate_local(506) 
  reac_source_local(13,506) = + reac_rate_local(506) * 2.d0
  reac_source_local(11,507) = + reac_rate_local(507) 
  reac_source_local(13,507) = - reac_rate_local(507) * 2.d0
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
  if( lgas_heating ) ZDPlasKin_cfg(1) = y(53)
  density(:) = y(1:species_max)
  call ZDPlasKin_reac_rates(t)
  rrt(001) = rrt(001) * density(01) * density(51) 
  rrt(002) = rrt(002) * density(01) * density(51) 
  rrt(003) = rrt(003) * density(01) * density(43) 
  rrt(004) = rrt(004) * density(01) * density(43) 
  rrt(005) = rrt(005) * density(01) * density(42) 
  rrt(006) = rrt(006) * density(01) * density(42) 
  rrt(007) = rrt(007) * density(01) * density(21) 
  rrt(008) = rrt(008) * density(01) * density(21) 
  rrt(009) = rrt(009) * density(01) * density(21) 
  rrt(010) = rrt(010) * density(01) * density(45) 
  rrt(011) = rrt(011) * density(01) * density(45) 
  rrt(012) = rrt(012) * density(01) * density(37) 
  rrt(013) = rrt(013) * density(01) * density(51) 
  rrt(014) = rrt(014) * density(01) * density(51) 
  rrt(015) = rrt(015) * density(01) * density(51) 
  rrt(016) = rrt(016) * density(01) * density(23) 
  rrt(017) = rrt(017) * density(01) * density(23) 
  rrt(018) = rrt(018) * density(01) * density(28) 
  rrt(019) = rrt(019) * density(01) * density(51) 
  rrt(020) = rrt(020) * density(01) * density(51) 
  rrt(021) = rrt(021) * density(01) * density(51) 
  rrt(022) = rrt(022) * density(01) * density(51) 
  rrt(023) = rrt(023) * density(01) * density(23) 
  rrt(024) = rrt(024) * density(01) * density(23) 
  rrt(025) = rrt(025) * density(01) * density(23) 
  rrt(026) = rrt(026) * density(01) * density(28) 
  rrt(027) = rrt(027) * density(01) * density(28) 
  rrt(028) = rrt(028) * density(01) * density(08) 
  rrt(029) = rrt(029) * density(01) * density(43) 
  rrt(030) = rrt(030) * density(01) * density(43) 
  rrt(031) = rrt(031) * density(01) * density(43) 
  rrt(032) = rrt(032) * density(01) * density(43) 
  rrt(033) = rrt(033) * density(01) * density(43) 
  rrt(034) = rrt(034) * density(01) * density(43) 
  rrt(035) = rrt(035) * density(01) * density(34) 
  rrt(036) = rrt(036) * density(01) * density(34) 
  rrt(037) = rrt(037) * density(01) * density(34) 
  rrt(038) = rrt(038) * density(01) * density(34) 
  rrt(039) = rrt(039) * density(01) * density(34) 
  rrt(040) = rrt(040) * density(01) * density(34) 
  rrt(041) = rrt(041) * density(01) * density(42) 
  rrt(042) = rrt(042) * density(01) * density(42) 
  rrt(043) = rrt(043) * density(01) * density(42) 
  rrt(044) = rrt(044) * density(01) * density(42) 
  rrt(045) = rrt(045) * density(01) * density(42) 
  rrt(046) = rrt(046) * density(01) * density(29) 
  rrt(047) = rrt(047) * density(01) * density(29) 
  rrt(048) = rrt(048) * density(01) * density(21) 
  rrt(049) = rrt(049) * density(01) * density(45) 
  rrt(050) = rrt(050) * density(01) * density(45) 
  rrt(051) = rrt(051) * density(01) * density(45) 
  rrt(052) = rrt(052) * density(01) * density(45) 
  rrt(053) = rrt(053) * density(01) * density(45) 
  rrt(054) = rrt(054) * density(01) * density(45) 
  rrt(055) = rrt(055) * density(01) * density(18) 
  rrt(056) = rrt(056) * density(01) * density(18) 
  rrt(057) = rrt(057) * density(01) * density(18) 
  rrt(058) = rrt(058) * density(01) * density(18) 
  rrt(059) = rrt(059) * density(01) * density(18) 
  rrt(060) = rrt(060) * density(01) * density(37) 
  rrt(061) = rrt(061) * density(01) * density(37) 
  rrt(062) = rrt(062) * density(01) * density(37) 
  rrt(063) = rrt(063) * density(01) * density(37) 
  rrt(064) = rrt(064) * density(01) * density(37) 
  rrt(065) = rrt(065) * density(01) * density(20) 
  rrt(066) = rrt(066) * density(01) * density(20) 
  rrt(067) = rrt(067) * density(01) * density(35) 
  rrt(068) = rrt(068) * density(01) * density(35) 
  rrt(069) = rrt(069) * density(01) * density(43) 
  rrt(070) = rrt(070) * density(01) * density(43) 
  rrt(071) = rrt(071) * density(01) * density(43) 
  rrt(072) = rrt(072) * density(01) * density(43) 
  rrt(073) = rrt(073) * density(01) * density(43) 
  rrt(074) = rrt(074) * density(01) * density(43) 
  rrt(075) = rrt(075) * density(01) * density(43) 
  rrt(076) = rrt(076) * density(01) * density(34) 
  rrt(077) = rrt(077) * density(01) * density(34) 
  rrt(078) = rrt(078) * density(01) * density(34) 
  rrt(079) = rrt(079) * density(01) * density(34) 
  rrt(080) = rrt(080) * density(01) * density(34) 
  rrt(081) = rrt(081) * density(01) * density(34) 
  rrt(082) = rrt(082) * density(01) * density(34) 
  rrt(083) = rrt(083) * density(01) * density(42) 
  rrt(084) = rrt(084) * density(01) * density(42) 
  rrt(085) = rrt(085) * density(01) * density(42) 
  rrt(086) = rrt(086) * density(01) * density(42) 
  rrt(087) = rrt(087) * density(01) * density(42) 
  rrt(088) = rrt(088) * density(01) * density(29) 
  rrt(089) = rrt(089) * density(01) * density(29) 
  rrt(090) = rrt(090) * density(01) * density(29) 
  rrt(091) = rrt(091) * density(01) * density(29) 
  rrt(092) = rrt(092) * density(01) * density(29) 
  rrt(093) = rrt(093) * density(01) * density(21) 
  rrt(094) = rrt(094) * density(01) * density(21) 
  rrt(095) = rrt(095) * density(01) * density(45) 
  rrt(096) = rrt(096) * density(01) * density(45) 
  rrt(097) = rrt(097) * density(01) * density(45) 
  rrt(098) = rrt(098) * density(01) * density(45) 
  rrt(099) = rrt(099) * density(01) * density(45) 
  rrt(100) = rrt(100) * density(01) * density(45) 
  rrt(101) = rrt(101) * density(01) * density(45) 
  rrt(102) = rrt(102) * density(01) * density(45) 
  rrt(103) = rrt(103) * density(01) * density(45) 
  rrt(104) = rrt(104) * density(01) * density(18) 
  rrt(105) = rrt(105) * density(01) * density(18) 
  rrt(106) = rrt(106) * density(01) * density(18) 
  rrt(107) = rrt(107) * density(01) * density(18) 
  rrt(108) = rrt(108) * density(01) * density(18) 
  rrt(109) = rrt(109) * density(01) * density(18) 
  rrt(110) = rrt(110) * density(01) * density(18) 
  rrt(111) = rrt(111) * density(01) * density(18) 
  rrt(112) = rrt(112) * density(01) * density(18) 
  rrt(113) = rrt(113) * density(01) * density(18) 
  rrt(114) = rrt(114) * density(01) * density(18) 
  rrt(115) = rrt(115) * density(01) * density(37) 
  rrt(116) = rrt(116) * density(01) * density(37) 
  rrt(117) = rrt(117) * density(01) * density(37) 
  rrt(118) = rrt(118) * density(01) * density(37) 
  rrt(119) = rrt(119) * density(01) * density(37) 
  rrt(120) = rrt(120) * density(01) * density(37) 
  rrt(121) = rrt(121) * density(01) * density(37) 
  rrt(122) = rrt(122) * density(01) * density(37) 
  rrt(123) = rrt(123) * density(01) * density(37) 
  rrt(124) = rrt(124) * density(01) * density(37) 
  rrt(125) = rrt(125) * density(01) * density(37) 
  rrt(126) = rrt(126) * density(01) * density(20) 
  rrt(127) = rrt(127) * density(01) * density(20) 
  rrt(128) = rrt(128) * density(01) * density(20) 
  rrt(129) = rrt(129) * density(01) * density(20) 
  rrt(130) = rrt(130) * density(01) * density(20) 
  rrt(131) = rrt(131) * density(01) * density(20) 
  rrt(132) = rrt(132) * density(01) * density(20) 
  rrt(133) = rrt(133) * density(01) * density(20) 
  rrt(134) = rrt(134) * density(01) * density(35) 
  rrt(135) = rrt(135) * density(01) * density(35) 
  rrt(136) = rrt(136) * density(01) * density(35) 
  rrt(137) = rrt(137) * density(01) * density(35) 
  rrt(138) = rrt(138) * density(01) * density(35) 
  rrt(139) = rrt(139) * density(01) * density(22) 
  rrt(140) = rrt(140) * density(01) * density(25) 
  rrt(141) = rrt(141) * density(01) * density(22) 
  rrt(142) = rrt(142) * density(01) * density(25) 
  rrt(143) = rrt(143) * density(01) * density(22) 
  rrt(144) = rrt(144) * density(01) * density(25) 
  rrt(145) = rrt(145) * density(01) * density(22) 
  rrt(146) = rrt(146) * density(01) * density(25) 
  rrt(147) = rrt(147) * density(01) * density(22) 
  rrt(148) = rrt(148) * density(01) * density(25) 
  rrt(149) = rrt(149) * density(01) * density(22) 
  rrt(150) = rrt(150) * density(01) * density(25) 
  rrt(151) = rrt(151) * density(01) * density(22) 
  rrt(152) = rrt(152) * density(01) * density(25) 
  rrt(153) = rrt(153) * density(01) * density(39) 
  rrt(154) = rrt(154) * density(01) * density(41) 
  rrt(155) = rrt(155) * density(01) * density(39) 
  rrt(156) = rrt(156) * density(01) * density(41) 
  rrt(157) = rrt(157) * density(01) * density(39) 
  rrt(158) = rrt(158) * density(01) * density(41) 
  rrt(159) = rrt(159) * density(01) * density(39) 
  rrt(160) = rrt(160) * density(01) * density(41) 
  rrt(161) = rrt(161) * density(01) * density(39) 
  rrt(162) = rrt(162) * density(01) * density(41) 
  rrt(163) = rrt(163) * density(01) * density(39) 
  rrt(164) = rrt(164) * density(01) * density(41) 
  rrt(165) = rrt(165) * density(01) * density(33) 
  rrt(166) = rrt(166) * density(01) * density(02) 
  rrt(167) = rrt(167) * density(01) * density(33) 
  rrt(168) = rrt(168) * density(01) * density(02) 
  rrt(169) = rrt(169) * density(01) * density(33) 
  rrt(170) = rrt(170) * density(01) * density(02) 
  rrt(171) = rrt(171) * density(01) * density(33) 
  rrt(172) = rrt(172) * density(01) * density(02) 
  rrt(173) = rrt(173) * density(01) * density(33) 
  rrt(174) = rrt(174) * density(01) * density(02) 
  rrt(175) = rrt(175) * density(01) * density(52) 
  rrt(176) = rrt(176) * density(01) * density(09) 
  rrt(177) = rrt(177) * density(01) * density(38) 
  rrt(178) = rrt(178) * density(01) * density(32) 
  rrt(179) = rrt(179) * density(01) * density(24) 
  rrt(180) = rrt(180) * density(01) * density(32) 
  rrt(181) = rrt(181) * density(01) * density(24) 
  rrt(182) = rrt(182) * density(01) * density(32) 
  rrt(183) = rrt(183) * density(01) * density(24) 
  rrt(184) = rrt(184) * density(01) * density(32) 
  rrt(185) = rrt(185) * density(01) * density(24) 
  rrt(186) = rrt(186) * density(01) * density(32) 
  rrt(187) = rrt(187) * density(01) * density(24) 
  rrt(188) = rrt(188) * density(01) * density(32) 
  rrt(189) = rrt(189) * density(01) * density(24) 
  rrt(190) = rrt(190) * density(01) * density(39) 
  rrt(191) = rrt(191) * density(01) * density(41) 
  rrt(192) = rrt(192) * density(01) * density(39) 
  rrt(193) = rrt(193) * density(01) * density(41) 
  rrt(194) = rrt(194) * density(01) * density(39) 
  rrt(195) = rrt(195) * density(01) * density(41) 
  rrt(196) = rrt(196) * density(01) * density(39) 
  rrt(197) = rrt(197) * density(01) * density(41) 
  rrt(198) = rrt(198) * density(01) * density(39) 
  rrt(199) = rrt(199) * density(01) * density(41) 
  rrt(200) = rrt(200) * density(01) * density(39) 
  rrt(201) = rrt(201) * density(01) * density(41) 
  rrt(202) = rrt(202) * density(01) * density(39) 
  rrt(203) = rrt(203) * density(01) * density(41) 
  rrt(204) = rrt(204) * density(01) * density(33) 
  rrt(205) = rrt(205) * density(01) * density(02) 
  rrt(206) = rrt(206) * density(01) * density(33) 
  rrt(207) = rrt(207) * density(01) * density(02) 
  rrt(208) = rrt(208) * density(01) * density(33) 
  rrt(209) = rrt(209) * density(01) * density(02) 
  rrt(210) = rrt(210) * density(01) * density(33) 
  rrt(211) = rrt(211) * density(01) * density(02) 
  rrt(212) = rrt(212) * density(01) * density(33) 
  rrt(213) = rrt(213) * density(01) * density(02) 
  rrt(214) = rrt(214) * density(01) * density(52) 
  rrt(215) = rrt(215) * density(01) * density(09) 
  rrt(216) = rrt(216) * density(01) * density(38) 
  rrt(217) = rrt(217) * density(01) * density(52) 
  rrt(218) = rrt(218) * density(01) * density(09) 
  rrt(219) = rrt(219) * density(01) * density(38) 
  rrt(220) = rrt(220) * density(01) * density(32) 
  rrt(221) = rrt(221) * density(01) * density(24) 
  rrt(222) = rrt(222) * density(01) * density(32) 
  rrt(223) = rrt(223) * density(01) * density(24) 
  rrt(224) = rrt(224) * density(01) * density(32) 
  rrt(225) = rrt(225) * density(01) * density(24) 
  rrt(226) = rrt(226) * density(01) * density(32) 
  rrt(227) = rrt(227) * density(01) * density(24) 
  rrt(228) = rrt(228) * density(01) * density(32) 
  rrt(229) = rrt(229) * density(01) * density(24) 
  rrt(230) = rrt(230) * density(01) * density(32) 
  rrt(231) = rrt(231) * density(01) * density(24) 
  rrt(232) = rrt(232) * density(01) * density(32) 
  rrt(233) = rrt(233) * density(01) * density(24) 
  rrt(234) = rrt(234) * density(01) * density(32) 
  rrt(235) = rrt(235) * density(01) * density(24) 
  rrt(236) = rrt(236) * density(01) * density(32) 
  rrt(237) = rrt(237) * density(01) * density(24) 
  rrt(238) = rrt(238) * density(01) * density(11) 
  rrt(239) = rrt(239) * density(01) * density(11) 
  rrt(240) = rrt(240) * density(01)**2 * density(10) 
  rrt(241) = rrt(241) * density(01)**2 * density(16) 
  rrt(242) = rrt(242) * density(01)**2 * density(50) 
  rrt(243) = rrt(243) * density(01)**2 * density(44) 
  rrt(244) = rrt(244) * density(01)**2 * density(40) 
  rrt(245) = rrt(245) * density(01)**2 * density(04) 
  rrt(246) = rrt(246) * density(01)**2 * density(46) 
  rrt(247) = rrt(247) * density(01)**2 * density(07) 
  rrt(248) = rrt(248) * density(01)**2 * density(17) 
  rrt(249) = rrt(249) * density(01)**2 * density(14) 
  rrt(250) = rrt(250) * density(01)**2 * density(03) 
  rrt(251) = rrt(251) * density(01)**2 * density(05) 
  rrt(252) = rrt(252) * density(01)**2 * density(26) 
  rrt(253) = rrt(253) * density(01)**2 * density(19) 
  rrt(254) = rrt(254) * density(01)**2 * density(30) 
  rrt(255) = rrt(255) * density(01)**2 * density(27) 
  rrt(256) = rrt(256) * density(01) * density(47) 
  rrt(257) = rrt(257) * density(01) * density(47) 
  rrt(258) = rrt(258) * density(01) * density(10) 
  rrt(259) = rrt(259) * density(01) * density(10) 
  rrt(260) = rrt(260) * density(01) * density(10) 
  rrt(261) = rrt(261) * density(01) * density(16) 
  rrt(262) = rrt(262) * density(01) * density(16) 
  rrt(263) = rrt(263) * density(01) * density(50) 
  rrt(264) = rrt(264) * density(01) * density(40) 
  rrt(265) = rrt(265) * density(01) * density(40) 
  rrt(266) = rrt(266) * density(01) * density(04) 
  rrt(267) = rrt(267) * density(01) * density(04) 
  rrt(268) = rrt(268) * density(01) * density(04) 
  rrt(269) = rrt(269) * density(01) * density(04) 
  rrt(270) = rrt(270) * density(01) * density(04) 
  rrt(271) = rrt(271) * density(01) * density(46) 
  rrt(272) = rrt(272) * density(01) * density(46) 
  rrt(273) = rrt(273) * density(01) * density(07) 
  rrt(274) = rrt(274) * density(01) * density(17) 
  rrt(275) = rrt(275) * density(01) * density(27) 
  rrt(276) = rrt(276) * density(01) * density(12) 
  rrt(277) = rrt(277) * density(11) * density(27) 
  rrt(278) = rrt(278) * density(28) * density(47) 
  rrt(279) = rrt(279) * density(08) * density(47) 
  rrt(280) = rrt(280) * density(43) * density(47) 
  rrt(281) = rrt(281) * density(42) * density(47) 
  rrt(282) = rrt(282) * density(21) * density(47) 
  rrt(283) = rrt(283) * density(13) * density(47) 
  rrt(284) = rrt(284) * density(10) * density(51) 
  rrt(285) = rrt(285) * density(10) * density(43) 
  rrt(286) = rrt(286) * density(10) * density(42) 
  rrt(287) = rrt(287) * density(10) * density(42) 
  rrt(288) = rrt(288) * density(10) * density(21) 
  rrt(289) = rrt(289) * density(10) * density(21) 
  rrt(290) = rrt(290) * density(10) * density(11) 
  rrt(291) = rrt(291) * density(10) * density(13) 
  rrt(292) = rrt(292) * density(16) * density(51) 
  rrt(293) = rrt(293) * density(16) * density(51) 
  rrt(294) = rrt(294) * density(16) * density(28) 
  rrt(295) = rrt(295) * density(08) * density(16) 
  rrt(296) = rrt(296) * density(16) * density(43) 
  rrt(297) = rrt(297) * density(16) * density(42) 
  rrt(298) = rrt(298) * density(16) * density(29) 
  rrt(299) = rrt(299) * density(50) * density(51) 
  rrt(300) = rrt(300) * density(50) * density(51) 
  rrt(301) = rrt(301) * density(50) * density(51) 
  rrt(302) = rrt(302) * density(50) * density(51) 
  rrt(303) = rrt(303) * density(50) * density(51) 
  rrt(304) = rrt(304) * density(11) * density(50) 
  rrt(305) = rrt(305) * density(44) * density(51) 
  rrt(306) = rrt(306) * density(44) * density(51) 
  rrt(307) = rrt(307) * density(44) * density(51) 
  rrt(308) = rrt(308) * density(11) * density(44) 
  rrt(309) = rrt(309) * density(40) * density(42) 
  rrt(310) = rrt(310) * density(21) * density(40) 
  rrt(311) = rrt(311) * density(13) * density(40) 
  rrt(312) = rrt(312) * density(04) * density(13) 
  rrt(313) = rrt(313) * density(29) * density(46) 
  rrt(314) = rrt(314) * density(29) * density(46) 
  rrt(315) = rrt(315) * density(13) * density(46) 
  rrt(316) = rrt(316) * density(07) * density(43) 
  rrt(317) = rrt(317) * density(07) * density(42) 
  rrt(318) = rrt(318) * density(06) * density(07) 
  rrt(319) = rrt(319) * density(07) * density(13) 
  rrt(320) = rrt(320) * density(17) * density(51) 
  rrt(321) = rrt(321) * density(17) * density(43) 
  rrt(322) = rrt(322) * density(17) * density(43) 
  rrt(323) = rrt(323) * density(17) * density(42) 
  rrt(324) = rrt(324) * density(17) * density(29) 
  rrt(325) = rrt(325) * density(11) * density(17) 
  rrt(326) = rrt(326) * density(36) * density(51) 
  rrt(327) = rrt(327) * density(12) * density(51) 
  rrt(328) = rrt(328) * density(12) * density(23) 
  rrt(329) = rrt(329) * density(12) * density(28) 
  rrt(330) = rrt(330) * density(08) * density(12) 
  rrt(331) = rrt(331) * density(12) * density(43) 
  rrt(332) = rrt(332) * density(12) * density(34) 
  rrt(333) = rrt(333) * density(12) * density(42) 
  rrt(334) = rrt(334) * density(12) * density(42) 
  rrt(335) = rrt(335) * density(12) * density(29) 
  rrt(336) = rrt(336) * density(06) * density(12) 
  rrt(337) = rrt(337) * density(12) * density(21) 
  rrt(338) = rrt(338) * density(27) * density(51) 
  rrt(339) = rrt(339) * density(27) * density(51) 
  rrt(340) = rrt(340) * density(27) * density(51) 
  rrt(341) = rrt(341) * density(27) * density(28) 
  rrt(342) = rrt(342) * density(27) * density(28) 
  rrt(343) = rrt(343) * density(08) * density(27) 
  rrt(344) = rrt(344) * density(08) * density(27) 
  rrt(345) = rrt(345) * density(27) * density(43) 
  rrt(346) = rrt(346) * density(27) * density(43) 
  rrt(347) = rrt(347) * density(27) * density(43) 
  rrt(348) = rrt(348) * density(27) * density(43) 
  rrt(349) = rrt(349) * density(27) * density(43) 
  rrt(350) = rrt(350) * density(27) * density(42) 
  rrt(351) = rrt(351) * density(27) * density(42) 
  rrt(352) = rrt(352) * density(27) * density(42) 
  rrt(353) = rrt(353) * density(21) * density(27) 
  rrt(354) = rrt(354) * density(21) * density(27) 
  rrt(355) = rrt(355) * density(13) * density(27) 
  rrt(356) = rrt(356) * density(13) * density(27) 
  rrt(357) = rrt(357) * density(30) * density(51) 
  rrt(358) = rrt(358) * density(30) * density(51) 
  rrt(359) = rrt(359) * density(23) * density(30) 
  rrt(360) = rrt(360) * density(28) * density(30) 
  rrt(361) = rrt(361) * density(28) * density(30) 
  rrt(362) = rrt(362) * density(08) * density(30) 
  rrt(363) = rrt(363) * density(30) * density(43) 
  rrt(364) = rrt(364) * density(30) * density(43) 
  rrt(365) = rrt(365) * density(30) * density(43) 
  rrt(366) = rrt(366) * density(30) * density(34) 
  rrt(367) = rrt(367) * density(30) * density(34) 
  rrt(368) = rrt(368) * density(30) * density(42) 
  rrt(369) = rrt(369) * density(30) * density(42) 
  rrt(370) = rrt(370) * density(30) * density(42) 
  rrt(371) = rrt(371) * density(29) * density(30) 
  rrt(372) = rrt(372) * density(29) * density(30) 
  rrt(373) = rrt(373) * density(21) * density(30) 
  rrt(374) = rrt(374) * density(28) * density(51) 
  rrt(375) = rrt(375) * density(08) * density(51) 
  rrt(376) = rrt(376) * density(34) * density(51) 
  rrt(377) = rrt(377) * density(29) * density(51) 
  rrt(378) = rrt(378) * density(06) * density(51) 
  rrt(379) = rrt(379) * density(18) * density(51) 
  rrt(380) = rrt(380) * density(20) * density(51) 
  rrt(381) = rrt(381) * density(13) * density(51) 
  rrt(382) = rrt(382) * density(23) * density(51) 
  rrt(383) = rrt(383) * density(31) * density(51) 
  rrt(384) = rrt(384) * density(28) * density(51) 
  rrt(385) = rrt(385) * density(51) 
  rrt(386) = rrt(386) * density(23) 
  rrt(387) = rrt(387) * density(23) 
  rrt(388) = rrt(388) * density(23) * density(34) 
  rrt(389) = rrt(389) * density(28)**2 
  rrt(390) = rrt(390) * density(28) * density(34) 
  rrt(391) = rrt(391) * density(28) * density(29) 
  rrt(392) = rrt(392) * density(06) * density(28) 
  rrt(393) = rrt(393) * density(28) * density(45) 
  rrt(394) = rrt(394) * density(18) * density(28) 
  rrt(395) = rrt(395) * density(18) * density(28) 
  rrt(396) = rrt(396) * density(28) * density(37) 
  rrt(397) = rrt(397) * density(11) * density(28) 
  rrt(398) = rrt(398) * density(13) * density(28) 
  rrt(399) = rrt(399) * density(28) 
  rrt(400) = rrt(400) * density(28)**2 
  rrt(401) = rrt(401) * density(13) * density(28) 
  rrt(402) = rrt(402) * density(08) * density(43) 
  rrt(403) = rrt(403) * density(08) * density(43) 
  rrt(404) = rrt(404) * density(08) * density(11) 
  rrt(405) = rrt(405) * density(08) * density(23) 
  rrt(406) = rrt(406) * density(08) * density(28) 
  rrt(407) = rrt(407) * density(08) * density(11) 
  rrt(408) = rrt(408) * density(08) * density(29) 
  rrt(409) = rrt(409) * density(29) * density(43) 
  rrt(410) = rrt(410) * density(18) * density(43) 
  rrt(411) = rrt(411) * density(20) * density(43) 
  rrt(412) = rrt(412) * density(13) * density(43) 
  rrt(413) = rrt(413) * density(13) * density(43) 
  rrt(414) = rrt(414) * density(43) 
  rrt(415) = rrt(415) * density(31) * density(43) 
  rrt(416) = rrt(416) * density(08) * density(43) 
  rrt(417) = rrt(417) * density(28) * density(43) 
  rrt(418) = rrt(418) * density(34)**2 
  rrt(419) = rrt(419) * density(34) * density(42) 
  rrt(420) = rrt(420) * density(21) * density(34) 
  rrt(421) = rrt(421) * density(06) * density(34) 
  rrt(422) = rrt(422) * density(34) * density(45) 
  rrt(423) = rrt(423) * density(18) * density(34) 
  rrt(424) = rrt(424) * density(18) * density(34) 
  rrt(425) = rrt(425) * density(34) * density(37) 
  rrt(426) = rrt(426) * density(11) * density(34) 
  rrt(427) = rrt(427) * density(13) * density(34) 
  rrt(428) = rrt(428) * density(13) * density(34) 
  rrt(429) = rrt(429) * density(13) * density(34) 
  rrt(430) = rrt(430) * density(34) 
  rrt(431) = rrt(431) * density(34)**2 
  rrt(432) = rrt(432) * density(31) * density(34) 
  rrt(433) = rrt(433) * density(29) * density(34) 
  rrt(434) = rrt(434) * density(13) * density(42) 
  rrt(435) = rrt(435) * density(13) * density(42) 
  rrt(436) = rrt(436) * density(11) * density(42) 
  rrt(437) = rrt(437) * density(42) 
  rrt(438) = rrt(438) * density(37) * density(42) 
  rrt(439) = rrt(439) * density(21) * density(42) 
  rrt(440) = rrt(440) * density(37) * density(42) 
  rrt(441) = rrt(441) * density(42)**2 
  rrt(442) = rrt(442) * density(23) * density(42) 
  rrt(443) = rrt(443) * density(42) 
  rrt(444) = rrt(444) * density(34) * density(42) 
  rrt(445) = rrt(445) * density(11) * density(42) 
  rrt(446) = rrt(446) * density(28) * density(42) 
  rrt(447) = rrt(447) * density(31) * density(42) 
  rrt(448) = rrt(448) * density(29)**2 
  rrt(449) = rrt(449) * density(29) * density(45) 
  rrt(450) = rrt(450) * density(18) * density(29) 
  rrt(451) = rrt(451) * density(18) * density(29) 
  rrt(452) = rrt(452) * density(29) * density(37) 
  rrt(453) = rrt(453) * density(20) * density(29) 
  rrt(454) = rrt(454) * density(11) * density(29) 
  rrt(455) = rrt(455) * density(13) * density(29) 
  rrt(456) = rrt(456) * density(13) * density(29) 
  rrt(457) = rrt(457) * density(29) 
  rrt(458) = rrt(458) * density(29) * density(31) 
  rrt(459) = rrt(459) * density(13) * density(21) 
  rrt(460) = rrt(460) * density(11) * density(21) 
  rrt(461) = rrt(461) * density(11) * density(21) 
  rrt(462) = rrt(462) * density(21) * density(23) 
  rrt(463) = rrt(463) * density(21) * density(31) 
  rrt(464) = rrt(464) * density(20) * density(45) 
  rrt(465) = rrt(465) * density(13) * density(45) 
  rrt(466) = rrt(466) * density(45) 
  rrt(467) = rrt(467) * density(31) * density(45) 
  rrt(468) = rrt(468) * density(28) * density(45) 
  rrt(469) = rrt(469) * density(18)**2 
  rrt(470) = rrt(470) * density(18) * density(37) 
  rrt(471) = rrt(471) * density(18) * density(20) 
  rrt(472) = rrt(472) * density(11) * density(18) 
  rrt(473) = rrt(473) * density(13) * density(18) 
  rrt(474) = rrt(474) * density(13) * density(18) 
  rrt(475) = rrt(475) * density(13) * density(18) 
  rrt(476) = rrt(476) * density(18) 
  rrt(477) = rrt(477) * density(18) 
  rrt(478) = rrt(478) * density(18) * density(31) 
  rrt(479) = rrt(479) * density(21) * density(37) 
  rrt(480) = rrt(480) * density(37)**2 
  rrt(481) = rrt(481) * density(37) 
  rrt(482) = rrt(482) * density(13) * density(37) 
  rrt(483) = rrt(483) * density(13) * density(37) 
  rrt(484) = rrt(484) * density(37) 
  rrt(485) = rrt(485) * density(23) * density(37) 
  rrt(486) = rrt(486) * density(31) * density(37) 
  rrt(487) = rrt(487) * density(11) * density(20) 
  rrt(488) = rrt(488) * density(13) * density(20) 
  rrt(489) = rrt(489) * density(20) 
  rrt(490) = rrt(490) * density(31) 
  rrt(491) = rrt(491) * density(28) * density(31) 
  rrt(492) = rrt(492) * density(31) 
  rrt(493) = rrt(493) * density(11) * density(31) 
  rrt(494) = rrt(494) * density(23) * density(49) 
  rrt(495) = rrt(495) * density(49) 
  rrt(496) = rrt(496) * density(49) 
  rrt(497) = rrt(497) * density(13) * density(49) 
  rrt(498) = rrt(498) * density(28) * density(49) 
  rrt(499) = rrt(499) * density(29) * density(49) 
  rrt(500) = rrt(500) * density(18) * density(49) 
  rrt(501) = rrt(501) * density(06) * density(49) 
  rrt(502) = rrt(502) * density(34) * density(49) 
  rrt(503) = rrt(503) * density(20) * density(49) 
  rrt(504) = rrt(504) * density(28) * density(49) 
  rrt(505) = rrt(505) * density(48) 
  rrt(506) = rrt(506) * density(11) 
  rrt(507) = rrt(507) * density(13)**2 
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
             +rrt(232)+rrt(233)+rrt(234)+rrt(235)+rrt(236)+rrt(237)+rrt(239)-rrt(240)-rrt(241)-rrt(242)-rrt(243)-rrt(244)-rrt(245)&
             -rrt(246)-rrt(247)-rrt(248)-rrt(249)-rrt(250)-rrt(251)-rrt(252)-rrt(253)-rrt(254)-rrt(255)-rrt(256)-rrt(257)-rrt(258)&
             -rrt(259)-rrt(260)-rrt(261)-rrt(262)-rrt(263)-rrt(264)-rrt(265)-rrt(266)-rrt(267)-rrt(268)-rrt(269)-rrt(270)-rrt(271)&
             -rrt(272)-rrt(273)-rrt(274)-rrt(275)-rrt(276) 
  ydot(02) = +rrt(006)-rrt(166)-rrt(168)-rrt(170)-rrt(172)-rrt(174)-rrt(205)-rrt(207)-rrt(209)-rrt(211)-rrt(213) 
  ydot(03) = +rrt(096)+rrt(104)+rrt(222)+rrt(223)-rrt(250) 
  ydot(04) = +rrt(070)+rrt(076)+rrt(100)+rrt(108)+rrt(118)+rrt(192)+rrt(193)+rrt(230)+rrt(231)-rrt(245)-rrt(266)-rrt(267)-rrt(268)&
             -rrt(269)-rrt(270)+rrt(280)+rrt(281)+rrt(286)+rrt(293)+rrt(296)+rrt(300)+rrt(310)+rrt(311)-rrt(312)+rrt(313)+rrt(316)&
             +rrt(317)+rrt(321)+rrt(331)+rrt(333)+rrt(346)+rrt(363) 
  ydot(05) = +rrt(097)+rrt(105)+rrt(115)+rrt(224)+rrt(225)-rrt(251) 
  ydot(06) = -rrt(318)-rrt(336)-rrt(378)-rrt(392)+rrt(420)-rrt(421)-rrt(501) 
  ydot(07) = +rrt(072)+rrt(078)+rrt(084)+rrt(088)+rrt(110)+rrt(120)+rrt(129)+rrt(135)+rrt(196)+rrt(197)+rrt(206)+rrt(207)-rrt(247)&
             -rrt(273)+rrt(282)+rrt(288)+rrt(294)+rrt(297)+rrt(298)+rrt(302)+rrt(306)+rrt(314)+rrt(315)-rrt(316)-rrt(317)-rrt(318)&
             -rrt(319)+rrt(320)+rrt(324)+rrt(325)+rrt(334)+rrt(337)+rrt(348)+rrt(351)+rrt(353)+rrt(365)+rrt(367)+rrt(369)+rrt(371) 
  ydot(08) = +rrt(015)+rrt(017)+rrt(018)-rrt(028)+rrt(039)+rrt(044)+rrt(047)+  2.d0 * rrt(048)+rrt(067)+rrt(085)+rrt(090)+rrt(094)&
             +rrt(118)+rrt(128)+rrt(135)+rrt(143)+rrt(144)+rrt(171)+rrt(172)+  2.d0 * rrt(175)+  2.d0 * rrt(176)+  2.d0 * rrt(177)&
             +rrt(208)+rrt(209)+rrt(217)+rrt(218)+rrt(219)+rrt(243)+rrt(260)+rrt(262)+rrt(263)+  2.d0 * rrt(274)-rrt(279)-rrt(295)&
             -rrt(330)-rrt(343)-rrt(344)-rrt(362)-rrt(375)+rrt(387)+rrt(392)+rrt(398)+rrt(399)-rrt(402)-rrt(403)-rrt(404)-rrt(405)&
             -rrt(406)-rrt(407)-rrt(408)-rrt(416) 
  ydot(09) = +rrt(008)-rrt(176)-rrt(215)-rrt(218) 
  ydot(10) = +rrt(019)+rrt(111)+rrt(122)+rrt(145)+rrt(146)-rrt(240)-rrt(258)-rrt(259)-rrt(260)+rrt(283)-rrt(284)-rrt(285)-rrt(286)&
             -rrt(287)-rrt(288)-rrt(289)-rrt(290)-rrt(291)+rrt(292)+rrt(328)+rrt(339)+rrt(357) 
  ydot(11) = +rrt(014)+rrt(015)+rrt(017)+rrt(021)+rrt(022)+rrt(025)+rrt(030)+rrt(031)+  2.d0 * rrt(032)+rrt(036)+rrt(038)+rrt(042)&
             +rrt(050)+  2.d0 * rrt(051)+rrt(056)+rrt(057)+rrt(061)+rrt(071)+rrt(072)+  2.d0 * rrt(073)+rrt(078)+rrt(079)+rrt(097)&
             +rrt(098)+  2.d0 * rrt(099)+rrt(106)+rrt(107)+rrt(117)+rrt(141)+rrt(142)+rrt(143)+rrt(144)+rrt(149)+rrt(150)+rrt(151)&
             +rrt(152)+rrt(155)+rrt(156)+rrt(157)+rrt(158)+  2.d0 * rrt(159)+  2.d0 * rrt(160)+rrt(167)+rrt(168)+rrt(180)+rrt(181)&
             +  2.d0 * rrt(182)+  2.d0 * rrt(183)+rrt(194)+rrt(195)+rrt(196)+rrt(197)+  2.d0 * rrt(198)+  2.d0 * rrt(199)+rrt(224)&
             +rrt(225)+rrt(226)+rrt(227)+  2.d0 * rrt(228)+  2.d0 * rrt(229)-rrt(238)-rrt(239)+rrt(255)+rrt(257)+rrt(260)+rrt(262)&
             +rrt(268)+rrt(276)-rrt(277)+rrt(280)+rrt(283)+rrt(285)-rrt(290)+rrt(291)+rrt(293)+rrt(294)+rrt(295)+rrt(301)+rrt(302)&
             +  2.d0 * rrt(303)-rrt(304)+rrt(306)+rrt(307)-rrt(308)+rrt(311)+rrt(312)+rrt(315)+rrt(319)-rrt(325)+rrt(327)+rrt(328)&
             +rrt(329)+rrt(330)+  2.d0 * rrt(331)+rrt(332)+rrt(333)+  2.d0 * rrt(334)+rrt(335)+rrt(336)+rrt(337)+rrt(339)+rrt(340)&
             +rrt(342)+rrt(344)+rrt(345)+rrt(346)+  2.d0 * rrt(347)+  2.d0 * rrt(348)+  3.d0 * rrt(349)+rrt(350)+rrt(351)&
             +  2.d0 * rrt(352)+rrt(354)+rrt(356)+rrt(358)+rrt(361)+rrt(363)+rrt(364)+  2.d0 * rrt(365)+rrt(366)+rrt(367)+rrt(369)&
             +rrt(370)+rrt(372)+rrt(381)+rrt(387)-rrt(397)+rrt(398)+rrt(400)-rrt(404)-rrt(407)+rrt(412)-rrt(426)+rrt(428)+rrt(434)&
             -rrt(436)+rrt(443)-rrt(445)-rrt(454)+rrt(455)-rrt(460)-rrt(461)+rrt(465)-rrt(472)+rrt(473)+rrt(482)-rrt(487)-rrt(493)&
             +rrt(497)-rrt(506)+rrt(507) 
  ydot(12) = -rrt(276)+rrt(277)-rrt(327)-rrt(328)-rrt(329)-rrt(330)-rrt(331)-rrt(332)-rrt(333)-rrt(334)-rrt(335)-rrt(336)-rrt(337)&
             +rrt(355) 
  ydot(13) = +rrt(013)+rrt(015)+rrt(016)+rrt(018)+rrt(020)+rrt(022)+rrt(024)+rrt(027)+rrt(029)+rrt(031)+rrt(035)+  2.d0 * rrt(037)&
             +rrt(038)+rrt(041)+  2.d0 * rrt(043)+rrt(046)+rrt(049)+rrt(055)+rrt(057)+rrt(060)+rrt(065)+rrt(070)+rrt(072)+rrt(077)&
             +rrt(079)+rrt(084)+rrt(089)+rrt(096)+rrt(098)+rrt(105)+rrt(107)+rrt(116)+rrt(127)+rrt(139)+rrt(140)+rrt(143)+rrt(144)&
             +rrt(147)+rrt(148)+rrt(151)+rrt(152)+rrt(153)+rrt(154)+rrt(157)+rrt(158)+rrt(165)+rrt(166)+  2.d0 * rrt(169)&
             +  2.d0 * rrt(170)+rrt(178)+rrt(179)+rrt(192)+rrt(193)+rrt(196)+rrt(197)+rrt(206)+rrt(207)+rrt(222)+rrt(223)+rrt(226)&
             +rrt(227)+  2.d0 * rrt(238)+rrt(254)+  2.d0 * rrt(256)+rrt(257)+rrt(258)+  2.d0 * rrt(259)+rrt(260)+rrt(261)+rrt(263)&
             +rrt(264)+  2.d0 * rrt(265)+rrt(266)+  2.d0 * rrt(267)+rrt(268)+  3.d0 * rrt(269)+rrt(271)+  2.d0 * rrt(272)+rrt(273)&
             +  2.d0 * rrt(275)+rrt(276)+rrt(277)-rrt(283)+rrt(290)-rrt(291)+rrt(300)+rrt(302)+rrt(304)+rrt(305)+rrt(307)+rrt(308)&
             -rrt(311)-rrt(312)-rrt(315)-rrt(319)+rrt(325)+rrt(338)+rrt(340)+rrt(341)+rrt(343)+rrt(346)+rrt(348)+rrt(351)+rrt(353)&
             -rrt(355)-rrt(356)+rrt(357)+rrt(359)+rrt(360)+rrt(362)+rrt(364)+rrt(367)+rrt(368)+rrt(370)+rrt(371)+rrt(373)+rrt(375)&
             -rrt(381)+rrt(382)+rrt(385)+rrt(386)+  2.d0 * rrt(389)+rrt(397)-rrt(398)+rrt(399)-rrt(401)+rrt(402)+rrt(404)+rrt(405)&
             +rrt(406)-rrt(412)-rrt(413)+rrt(426)-rrt(427)-rrt(428)-rrt(429)+rrt(430)-rrt(434)-rrt(435)+rrt(436)+rrt(437)+rrt(454)&
             -rrt(455)-rrt(456)+rrt(457)-rrt(459)+rrt(461)-rrt(465)+rrt(472)-rrt(473)-rrt(474)-rrt(475)+rrt(476)+rrt(481)-rrt(482)&
             -rrt(483)+rrt(487)-rrt(488)+rrt(493)-rrt(497)+  2.d0 * rrt(506)-  2.d0 * rrt(507) 
  ydot(14) = +rrt(095)+rrt(220)+rrt(221)-rrt(249) 
  ydot(15) = +rrt(012) 
  ydot(16) = +rrt(020)+rrt(023)+rrt(074)+rrt(080)+rrt(085)+rrt(102)+rrt(112)+rrt(123)+rrt(131)+rrt(147)+rrt(148)+rrt(200)+rrt(201)&
             +rrt(208)+rrt(209)+rrt(234)+rrt(235)-rrt(241)-rrt(261)-rrt(262)+rrt(278)+rrt(291)-rrt(292)-rrt(293)-rrt(294)-rrt(295)&
             -rrt(296)-rrt(297)-rrt(298)+rrt(299)+rrt(304)+rrt(329)+rrt(340)+rrt(341)+rrt(358)+rrt(359) 
  ydot(17) = +rrt(073)+rrt(079)+rrt(089)+rrt(093)+rrt(121)+rrt(130)+rrt(136)+rrt(198)+rrt(199)+rrt(214)+rrt(215)+rrt(216)-rrt(248)&
             -rrt(274)+rrt(289)+rrt(295)+rrt(303)+rrt(307)+rrt(318)+rrt(319)-rrt(320)-rrt(321)-rrt(322)-rrt(323)-rrt(324)-rrt(325)&
             +rrt(326)+rrt(336)+rrt(349)+rrt(352)+rrt(354)+rrt(370)+rrt(372)+rrt(373) 
  ydot(18) = +rrt(049)-rrt(055)-rrt(056)-rrt(057)-rrt(058)-rrt(059)-rrt(104)-rrt(105)-rrt(106)-rrt(107)-rrt(108)-rrt(109)-rrt(110)&
             -rrt(111)-rrt(112)-rrt(113)-rrt(114)+rrt(178)+rrt(179)+rrt(250)-rrt(379)+rrt(393)-rrt(394)-rrt(395)+rrt(403)-rrt(410)&
             +rrt(422)-rrt(423)-rrt(424)+rrt(440)+rrt(442)+rrt(447)+rrt(449)-rrt(450)-rrt(451)+rrt(464)+rrt(465)+rrt(467)&
             -  2.d0 * rrt(469)-rrt(470)-rrt(471)-rrt(472)-rrt(473)-rrt(474)-rrt(475)-rrt(476)-rrt(477)-rrt(478)+rrt(480)+rrt(483)&
             +rrt(491)+rrt(495)-rrt(500) 
  ydot(19) = +rrt(099)+rrt(107)+rrt(117)+rrt(127)+rrt(134)+rrt(228)+rrt(229)-rrt(253) 
  ydot(20) = +rrt(056)+rrt(060)-rrt(065)-rrt(066)-rrt(126)-rrt(127)-rrt(128)-rrt(129)-rrt(130)-rrt(131)-rrt(132)-rrt(133)+rrt(252)&
             -rrt(380)+rrt(396)-rrt(411)+rrt(425)+rrt(438)+rrt(452)-rrt(453)+rrt(462)+rrt(463)-rrt(464)+rrt(470)-rrt(471)+rrt(479)&
             +rrt(480)+rrt(481)+rrt(482)+rrt(486)-rrt(487)-rrt(488)-rrt(489)-rrt(503) 
  ydot(21) = -rrt(007)-rrt(008)-rrt(009)+rrt(032)+rrt(038)+rrt(042)+rrt(043)+rrt(046)-rrt(048)+rrt(064)+rrt(066)+rrt(068)+rrt(092)&
             -rrt(093)-rrt(094)+rrt(122)+rrt(131)+rrt(137)+rrt(159)+rrt(160)+rrt(167)+rrt(168)+rrt(169)+rrt(170)+rrt(248)+rrt(268)&
             +rrt(269)+rrt(272)+rrt(273)-rrt(282)-rrt(288)-rrt(289)-rrt(310)+rrt(313)+rrt(317)+rrt(318)+rrt(323)+rrt(324)-rrt(337)&
             -rrt(353)-rrt(354)-rrt(373)+rrt(378)+rrt(389)+rrt(391)+rrt(392)+rrt(400)+rrt(406)+rrt(408)-rrt(420)+rrt(421)-rrt(439)&
             +rrt(443)+rrt(448)+rrt(450)+rrt(453)+rrt(455)+rrt(457)+rrt(458)-rrt(459)-rrt(460)-rrt(461)-rrt(462)-rrt(463)-rrt(479)&
             +rrt(489)+rrt(501) 
  ydot(22) = +rrt(001)-rrt(139)-rrt(141)-rrt(143)-rrt(145)-rrt(147)-rrt(149)-rrt(151) 
  ydot(23) = +rrt(013)-rrt(016)-rrt(017)-rrt(023)-rrt(024)-rrt(025)+  2.d0 * rrt(034)+rrt(040)+rrt(044)+rrt(053)+rrt(058)+rrt(063)&
             +rrt(066)+rrt(074)+rrt(081)+rrt(087)+rrt(100)+rrt(109)+rrt(120)+rrt(130)+rrt(139)+rrt(140)+  2.d0 * rrt(163)&
             +  2.d0 * rrt(164)+rrt(171)+rrt(172)+rrt(186)+rrt(187)+rrt(200)+rrt(201)+rrt(212)+rrt(213)+rrt(230)+rrt(231)+rrt(241)&
             +rrt(256)+rrt(258)+rrt(270)+rrt(284)+rrt(286)+rrt(288)+rrt(292)+rrt(298)+rrt(299)+rrt(320)+rrt(326)-rrt(328)-rrt(359)&
             +  2.d0 * rrt(374)+rrt(376)+rrt(377)+rrt(378)+rrt(379)+rrt(380)+rrt(381)-rrt(382)+rrt(383)+rrt(385)-rrt(386)-rrt(387)&
             -rrt(388)+rrt(390)+rrt(391)+rrt(393)+rrt(395)+rrt(396)+rrt(397)+rrt(401)-rrt(405)+rrt(407)+rrt(413)+  2.d0 * rrt(414)&
             +rrt(416)+rrt(417)+  2.d0 * rrt(427)-rrt(442)-rrt(462)+rrt(466)+rrt(475)+rrt(477)+rrt(484)-rrt(485)+rrt(489)+rrt(492)&
             -rrt(494)+rrt(495)+rrt(498)+rrt(505) 
  ydot(24) = +rrt(011)-rrt(179)-rrt(181)-rrt(183)-rrt(185)-rrt(187)-rrt(189)-rrt(221)-rrt(223)-rrt(225)-rrt(227)-rrt(229)-rrt(231)&
             -rrt(233)-rrt(235)-rrt(237) 
  ydot(25) = +rrt(002)-rrt(140)-rrt(142)-rrt(144)-rrt(146)-rrt(148)-rrt(150)-rrt(152) 
  ydot(26) = +rrt(098)+rrt(106)+rrt(116)+rrt(126)+rrt(226)+rrt(227)-rrt(252) 
  ydot(27) = +rrt(239)-rrt(255)-rrt(275)-rrt(277)-rrt(338)-rrt(339)-rrt(340)-rrt(341)-rrt(342)-rrt(343)-rrt(344)-rrt(345)-rrt(346)&
             -rrt(347)-rrt(348)-rrt(349)-rrt(350)-rrt(351)-rrt(352)-rrt(353)-rrt(354)-rrt(355)-rrt(356) 
  ydot(28) = +rrt(014)+rrt(016)-rrt(018)-rrt(026)-rrt(027)+rrt(033)+rrt(040)+  2.d0 * rrt(045)+rrt(047)+rrt(052)+rrt(062)+rrt(068)&
             +rrt(080)+rrt(086)+rrt(091)+rrt(108)+rrt(119)+rrt(129)+rrt(136)+rrt(141)+rrt(142)+rrt(161)+rrt(162)+  2.d0 * rrt(173)&
             +  2.d0 * rrt(174)+rrt(184)+rrt(185)+rrt(210)+rrt(211)+rrt(242)+rrt(257)+rrt(259)+rrt(261)+rrt(270)-rrt(278)-rrt(294)&
             -rrt(329)-rrt(341)-rrt(342)-rrt(360)-rrt(361)-rrt(374)-rrt(384)+rrt(386)+rrt(388)-  2.d0 * rrt(389)-rrt(390)-rrt(391)&
             -rrt(392)-rrt(393)-rrt(394)-rrt(395)-rrt(396)-rrt(397)-rrt(398)-rrt(399)-  2.d0 * rrt(400)-rrt(401)+rrt(404)-rrt(406)&
             +rrt(408)-rrt(417)-rrt(446)-rrt(468)-rrt(491)-rrt(498)-rrt(504) 
  ydot(29) = +rrt(031)+rrt(036)+rrt(037)+rrt(041)-rrt(046)-rrt(047)+rrt(059)+rrt(063)+rrt(067)-rrt(088)-rrt(089)-rrt(090)-rrt(091)&
             -rrt(092)+rrt(111)+rrt(123)+rrt(132)+rrt(138)+rrt(157)+rrt(158)+rrt(165)+rrt(166)+rrt(247)+rrt(267)+rrt(271)-rrt(298)&
             +rrt(310)-rrt(313)-rrt(314)+rrt(321)-rrt(324)-rrt(335)-rrt(371)-rrt(372)-rrt(377)-rrt(391)+rrt(405)-rrt(408)-rrt(409)&
             +rrt(419)-rrt(433)+rrt(434)+rrt(437)+  2.d0 * rrt(439)+rrt(440)+rrt(441)-  2.d0 * rrt(448)-rrt(449)-rrt(450)-rrt(451)&
             -rrt(452)-rrt(453)-rrt(454)-rrt(455)-rrt(456)-rrt(457)-rrt(458)+rrt(459)+rrt(461)+rrt(479)+rrt(484)-rrt(499) 
  ydot(30) = +rrt(092)-rrt(254)+rrt(356)-rrt(357)-rrt(358)-rrt(359)-rrt(360)-rrt(361)-rrt(362)-rrt(363)-rrt(364)-rrt(365)-rrt(366)&
             -rrt(367)-rrt(368)-rrt(369)-rrt(370)-rrt(371)-rrt(372)-rrt(373) 
  ydot(31) = -rrt(383)-rrt(415)-rrt(432)+rrt(444)-rrt(447)-rrt(458)-rrt(463)-rrt(467)-rrt(478)+rrt(485)-rrt(486)-rrt(490)-rrt(491)&
             -rrt(492)-rrt(493)+rrt(494)+rrt(497)+rrt(498)+rrt(499)+rrt(500)+rrt(501)+rrt(502)+rrt(503)+rrt(505) 
  ydot(32) = +rrt(010)-rrt(178)-rrt(180)-rrt(182)-rrt(184)-rrt(186)-rrt(188)-rrt(220)-rrt(222)-rrt(224)-rrt(226)-rrt(228)-rrt(230)&
             -rrt(232)-rrt(234)-rrt(236) 
  ydot(33) = +rrt(005)-rrt(165)-rrt(167)-rrt(169)-rrt(171)-rrt(173)-rrt(204)-rrt(206)-rrt(208)-rrt(210)-rrt(212) 
  ydot(34) = +rrt(029)-rrt(035)-rrt(036)-rrt(037)-rrt(038)-rrt(039)-rrt(040)+rrt(053)-rrt(076)-rrt(077)-rrt(078)-rrt(079)-rrt(080)&
             -rrt(081)-rrt(082)+rrt(102)+rrt(113)+rrt(125)+rrt(153)+rrt(154)+rrt(186)+rrt(187)+rrt(234)+rrt(235)+rrt(245)+rrt(264)&
             -rrt(332)-rrt(366)-rrt(367)-rrt(376)-rrt(388)-rrt(390)+rrt(394)+rrt(409)+rrt(410)+rrt(411)+rrt(412)+rrt(415)+rrt(417)&
             -  2.d0 * rrt(418)-rrt(419)-rrt(420)-rrt(421)-rrt(422)-rrt(423)-rrt(424)-rrt(425)-rrt(426)-rrt(427)-rrt(428)-rrt(429)&
             -rrt(430)-  2.d0 * rrt(431)-rrt(432)-rrt(433)+rrt(435)+rrt(436)+rrt(438)+rrt(441)-rrt(444)+rrt(466)+rrt(475)+rrt(490)&
             +  2.d0 * rrt(496)-rrt(502) 
  ydot(35) = +rrt(051)+rrt(057)+rrt(061)+rrt(065)-rrt(067)-rrt(068)-rrt(134)-rrt(135)-rrt(136)-rrt(137)-rrt(138)+rrt(182)+rrt(183)&
             +rrt(253) 
  ydot(36) = -rrt(326) 
  ydot(37) = -rrt(012)+rrt(050)+rrt(055)-rrt(060)-rrt(061)-rrt(062)-rrt(063)-rrt(064)-rrt(115)-rrt(116)-rrt(117)-rrt(118)-rrt(119)&
             -rrt(120)-rrt(121)-rrt(122)-rrt(123)-rrt(124)-rrt(125)+rrt(180)+rrt(181)+rrt(251)+rrt(380)+rrt(395)-rrt(396)+rrt(402)&
             +rrt(411)+rrt(424)-rrt(425)-rrt(438)-rrt(440)+rrt(446)+rrt(447)+rrt(451)-rrt(452)+rrt(453)+rrt(463)+rrt(464)+rrt(469)&
             -rrt(470)+  2.d0 * rrt(471)+rrt(473)+rrt(476)+rrt(478)-rrt(479)-  2.d0 * rrt(480)-rrt(481)-rrt(482)-rrt(483)-rrt(484)&
             -rrt(485)-rrt(486)+rrt(487)+rrt(488)+rrt(492)+rrt(503) 
  ydot(38) = +rrt(009)-rrt(177)-rrt(216)-rrt(219) 
  ydot(39) = +rrt(003)-rrt(153)-rrt(155)-rrt(157)-rrt(159)-rrt(161)-rrt(163)-rrt(190)-rrt(192)-rrt(194)-rrt(196)-rrt(198)-rrt(200)&
             -rrt(202) 
  ydot(40) = +rrt(069)+rrt(190)+rrt(191)-rrt(244)-rrt(264)-rrt(265)-rrt(309)-rrt(310)-rrt(311)+rrt(332)+rrt(345) 
  ydot(41) = +rrt(004)-rrt(154)-rrt(156)-rrt(158)-rrt(160)-rrt(162)-rrt(164)-rrt(191)-rrt(193)-rrt(195)-rrt(197)-rrt(199)-rrt(201)&
             -rrt(203) 
  ydot(42) = -rrt(005)-rrt(006)+rrt(030)+rrt(035)-rrt(041)-rrt(042)-rrt(043)-rrt(044)-rrt(045)+rrt(054)+rrt(058)+rrt(062)-rrt(083)&
             -rrt(084)-rrt(085)-rrt(086)-rrt(087)+rrt(112)+rrt(124)+rrt(133)+rrt(155)+rrt(156)+rrt(188)+rrt(189)+rrt(246)+rrt(265)&
             +rrt(266)-rrt(281)-rrt(286)-rrt(287)-rrt(297)-rrt(309)+rrt(314)+rrt(316)-rrt(317)+rrt(322)-rrt(323)-rrt(333)-rrt(334)&
             -rrt(350)-rrt(351)-rrt(352)-rrt(368)-rrt(369)-rrt(370)+rrt(375)+rrt(377)+rrt(390)+rrt(394)+rrt(409)+rrt(416)+rrt(418)&
             -rrt(419)+rrt(421)+rrt(423)+rrt(428)+rrt(430)+rrt(432)+  2.d0 * rrt(433)-rrt(434)-rrt(435)-rrt(436)-rrt(437)-rrt(438)&
             -rrt(439)-rrt(440)-  2.d0 * rrt(441)-rrt(442)-rrt(443)-rrt(444)-rrt(445)-rrt(446)-rrt(447)+rrt(448)+rrt(449)+rrt(451)&
             +rrt(452)+rrt(454)+rrt(456)+rrt(460)+rrt(477)+rrt(490)+rrt(491)+rrt(499) 
  ydot(43) = -rrt(003)-rrt(004)-rrt(029)-rrt(030)-rrt(031)-rrt(032)-rrt(033)-rrt(034)+rrt(052)-rrt(069)-rrt(070)-rrt(071)-rrt(072)&
             -rrt(073)-rrt(074)-rrt(075)+rrt(103)+rrt(114)+rrt(184)+rrt(185)+rrt(236)+rrt(237)+rrt(244)-rrt(280)-rrt(285)-rrt(296)&
             +rrt(309)-rrt(316)-rrt(321)-rrt(322)-rrt(331)-rrt(345)-rrt(346)-rrt(347)-rrt(348)-rrt(349)-rrt(363)-rrt(364)-rrt(365)&
             +rrt(376)+rrt(382)+rrt(384)+rrt(388)-rrt(402)-rrt(403)-rrt(409)-rrt(410)-rrt(411)-rrt(412)-rrt(413)-rrt(414)-rrt(415)&
             -rrt(416)-rrt(417)+rrt(418)+rrt(419)+rrt(420)+rrt(422)+rrt(424)+rrt(425)+rrt(426)+rrt(429)+rrt(445)+rrt(502) 
  ydot(44) = +rrt(022)+rrt(025)+rrt(027)+rrt(028)+rrt(082)+rrt(087)+rrt(091)+rrt(094)+rrt(114)+rrt(125)+rrt(133)+rrt(138)+rrt(151)&
             +rrt(152)+rrt(212)+rrt(213)+rrt(217)+rrt(218)+rrt(219)-rrt(243)-rrt(305)-rrt(306)-rrt(307)-rrt(308)+rrt(344)+rrt(361)&
             +rrt(362) 
  ydot(45) = -rrt(010)-rrt(011)-rrt(049)-rrt(050)-rrt(051)-rrt(052)-rrt(053)-rrt(054)-rrt(095)-rrt(096)-rrt(097)-rrt(098)-rrt(099)&
             -rrt(100)-rrt(101)-rrt(102)-rrt(103)+rrt(249)+rrt(379)-rrt(393)+rrt(410)-rrt(422)+rrt(423)-rrt(449)+rrt(450)-rrt(464)&
             -rrt(465)-rrt(466)-rrt(467)-rrt(468)+rrt(469)+rrt(470)+rrt(472)+rrt(474)+rrt(500) 
  ydot(46) = +rrt(071)+rrt(077)+rrt(083)+rrt(101)+rrt(109)+rrt(119)+rrt(128)+rrt(194)+rrt(195)+rrt(204)+rrt(205)+rrt(232)+rrt(233)&
             -rrt(246)-rrt(271)-rrt(272)+rrt(285)+rrt(287)+rrt(301)+rrt(305)+rrt(309)+rrt(312)-rrt(313)-rrt(314)-rrt(315)+rrt(322)&
             +rrt(323)+rrt(335)+rrt(347)+rrt(350)+rrt(364)+rrt(366)+rrt(368) 
  ydot(47) = -rrt(256)-rrt(257)-rrt(278)-rrt(279)-rrt(280)-rrt(281)-rrt(282)-rrt(283)+rrt(284)+rrt(290)+rrt(327)+rrt(338) 
  ydot(48) = +rrt(504)-rrt(505) 
  ydot(49) = +rrt(383)+rrt(415)+rrt(431)+rrt(432)+rrt(458)+rrt(467)+rrt(468)+rrt(478)+rrt(486)+rrt(493)-rrt(494)-rrt(495)-rrt(496)&
             -rrt(497)-rrt(498)-rrt(499)-rrt(500)-rrt(501)-rrt(502)-rrt(503)-rrt(504) 
  ydot(50) = +rrt(021)+rrt(024)+rrt(026)+rrt(075)+rrt(081)+rrt(086)+rrt(090)+rrt(103)+rrt(113)+rrt(124)+rrt(132)+rrt(137)+rrt(149)&
             +rrt(150)+rrt(202)+rrt(203)+rrt(210)+rrt(211)+rrt(236)+rrt(237)-rrt(242)-rrt(263)+rrt(279)-rrt(299)-rrt(300)-rrt(301)&
             -rrt(302)-rrt(303)-rrt(304)+rrt(308)+rrt(330)+rrt(342)+rrt(343)+rrt(360) 
  ydot(51) = -rrt(001)-rrt(002)-rrt(013)-rrt(014)-rrt(015)-rrt(019)-rrt(020)-rrt(021)-rrt(022)+rrt(033)+rrt(039)+rrt(054)+rrt(059)&
             +rrt(064)+rrt(075)+rrt(082)+rrt(101)+rrt(110)+rrt(121)+rrt(161)+rrt(162)+rrt(188)+rrt(189)+rrt(202)+rrt(203)+rrt(232)&
             +rrt(233)+rrt(240)+rrt(278)+rrt(279)+rrt(280)+rrt(281)+rrt(282)-rrt(284)+rrt(285)+rrt(287)+rrt(289)-rrt(292)-rrt(293)&
             +rrt(296)+rrt(297)-rrt(299)-rrt(300)-rrt(301)-rrt(302)-rrt(303)-rrt(305)-rrt(306)-rrt(307)-rrt(320)-rrt(326)-rrt(327)&
             -rrt(338)-rrt(339)-rrt(340)-rrt(357)-rrt(358)-rrt(374)-rrt(375)-rrt(376)-rrt(377)-rrt(378)-rrt(379)-rrt(380)-rrt(381)&
             -rrt(382)-rrt(383)-rrt(384)-rrt(385)+rrt(413)+rrt(494) 
  ydot(52) = +rrt(007)-rrt(175)-rrt(214)-rrt(217) 
  if( ldensity_constant ) where( density_constant(:) ) ydot(1:species_max) = 0.0d0
  ydot(53) = 0.0d0
  if( lgas_heating ) then
    ydot(53) = ( ZDPlasKin_cfg(14)/k_B + ydot(53) ) / ( sum(density(1:species_max)) - density(species_electrons) ) &
            + eV_to_K * ZDPlasKin_cfg(11) * density(species_electrons)
    ydot(53) = ydot(53) * ZDPlasKin_cfg(13)
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
  if( lgas_heating ) ZDPlasKin_cfg(1) = y(53)
  density(:) = y(1:species_max)
  call ZDPlasKin_reac_rates(t)
  pd(22,01) = pd(22,01) + rrt(001) * density(51) 
  pd(22,51) = pd(22,51) + rrt(001) * density(01) 
  pd(51,01) = pd(51,01) - rrt(001) * density(51) 
  pd(51,51) = pd(51,51) - rrt(001) * density(01) 
  pd(25,01) = pd(25,01) + rrt(002) * density(51) 
  pd(25,51) = pd(25,51) + rrt(002) * density(01) 
  pd(51,01) = pd(51,01) - rrt(002) * density(51) 
  pd(51,51) = pd(51,51) - rrt(002) * density(01) 
  pd(39,01) = pd(39,01) + rrt(003) * density(43) 
  pd(39,43) = pd(39,43) + rrt(003) * density(01) 
  pd(43,01) = pd(43,01) - rrt(003) * density(43) 
  pd(43,43) = pd(43,43) - rrt(003) * density(01) 
  pd(41,01) = pd(41,01) + rrt(004) * density(43) 
  pd(41,43) = pd(41,43) + rrt(004) * density(01) 
  pd(43,01) = pd(43,01) - rrt(004) * density(43) 
  pd(43,43) = pd(43,43) - rrt(004) * density(01) 
  pd(33,01) = pd(33,01) + rrt(005) * density(42) 
  pd(33,42) = pd(33,42) + rrt(005) * density(01) 
  pd(42,01) = pd(42,01) - rrt(005) * density(42) 
  pd(42,42) = pd(42,42) - rrt(005) * density(01) 
  pd(02,01) = pd(02,01) + rrt(006) * density(42) 
  pd(02,42) = pd(02,42) + rrt(006) * density(01) 
  pd(42,01) = pd(42,01) - rrt(006) * density(42) 
  pd(42,42) = pd(42,42) - rrt(006) * density(01) 
  pd(21,01) = pd(21,01) - rrt(007) * density(21) 
  pd(21,21) = pd(21,21) - rrt(007) * density(01) 
  pd(52,01) = pd(52,01) + rrt(007) * density(21) 
  pd(52,21) = pd(52,21) + rrt(007) * density(01) 
  pd(09,01) = pd(09,01) + rrt(008) * density(21) 
  pd(09,21) = pd(09,21) + rrt(008) * density(01) 
  pd(21,01) = pd(21,01) - rrt(008) * density(21) 
  pd(21,21) = pd(21,21) - rrt(008) * density(01) 
  pd(21,01) = pd(21,01) - rrt(009) * density(21) 
  pd(21,21) = pd(21,21) - rrt(009) * density(01) 
  pd(38,01) = pd(38,01) + rrt(009) * density(21) 
  pd(38,21) = pd(38,21) + rrt(009) * density(01) 
  pd(32,01) = pd(32,01) + rrt(010) * density(45) 
  pd(32,45) = pd(32,45) + rrt(010) * density(01) 
  pd(45,01) = pd(45,01) - rrt(010) * density(45) 
  pd(45,45) = pd(45,45) - rrt(010) * density(01) 
  pd(24,01) = pd(24,01) + rrt(011) * density(45) 
  pd(24,45) = pd(24,45) + rrt(011) * density(01) 
  pd(45,01) = pd(45,01) - rrt(011) * density(45) 
  pd(45,45) = pd(45,45) - rrt(011) * density(01) 
  pd(15,01) = pd(15,01) + rrt(012) * density(37) 
  pd(15,37) = pd(15,37) + rrt(012) * density(01) 
  pd(37,01) = pd(37,01) - rrt(012) * density(37) 
  pd(37,37) = pd(37,37) - rrt(012) * density(01) 
  pd(13,01) = pd(13,01) + rrt(013) * density(51) 
  pd(13,51) = pd(13,51) + rrt(013) * density(01) 
  pd(23,01) = pd(23,01) + rrt(013) * density(51) 
  pd(23,51) = pd(23,51) + rrt(013) * density(01) 
  pd(51,01) = pd(51,01) - rrt(013) * density(51) 
  pd(51,51) = pd(51,51) - rrt(013) * density(01) 
  pd(11,01) = pd(11,01) + rrt(014) * density(51) 
  pd(11,51) = pd(11,51) + rrt(014) * density(01) 
  pd(28,01) = pd(28,01) + rrt(014) * density(51) 
  pd(28,51) = pd(28,51) + rrt(014) * density(01) 
  pd(51,01) = pd(51,01) - rrt(014) * density(51) 
  pd(51,51) = pd(51,51) - rrt(014) * density(01) 
  pd(08,01) = pd(08,01) + rrt(015) * density(51) 
  pd(08,51) = pd(08,51) + rrt(015) * density(01) 
  pd(11,01) = pd(11,01) + rrt(015) * density(51) 
  pd(11,51) = pd(11,51) + rrt(015) * density(01) 
  pd(13,01) = pd(13,01) + rrt(015) * density(51) 
  pd(13,51) = pd(13,51) + rrt(015) * density(01) 
  pd(51,01) = pd(51,01) - rrt(015) * density(51) 
  pd(51,51) = pd(51,51) - rrt(015) * density(01) 
  pd(13,01) = pd(13,01) + rrt(016) * density(23) 
  pd(13,23) = pd(13,23) + rrt(016) * density(01) 
  pd(23,01) = pd(23,01) - rrt(016) * density(23) 
  pd(23,23) = pd(23,23) - rrt(016) * density(01) 
  pd(28,01) = pd(28,01) + rrt(016) * density(23) 
  pd(28,23) = pd(28,23) + rrt(016) * density(01) 
  pd(08,01) = pd(08,01) + rrt(017) * density(23) 
  pd(08,23) = pd(08,23) + rrt(017) * density(01) 
  pd(11,01) = pd(11,01) + rrt(017) * density(23) 
  pd(11,23) = pd(11,23) + rrt(017) * density(01) 
  pd(23,01) = pd(23,01) - rrt(017) * density(23) 
  pd(23,23) = pd(23,23) - rrt(017) * density(01) 
  pd(08,01) = pd(08,01) + rrt(018) * density(28) 
  pd(08,28) = pd(08,28) + rrt(018) * density(01) 
  pd(13,01) = pd(13,01) + rrt(018) * density(28) 
  pd(13,28) = pd(13,28) + rrt(018) * density(01) 
  pd(28,01) = pd(28,01) - rrt(018) * density(28) 
  pd(28,28) = pd(28,28) - rrt(018) * density(01) 
  pd(01,01) = pd(01,01) + rrt(019) * density(51) 
  pd(01,51) = pd(01,51) + rrt(019) * density(01) 
  pd(10,01) = pd(10,01) + rrt(019) * density(51) 
  pd(10,51) = pd(10,51) + rrt(019) * density(01) 
  pd(51,01) = pd(51,01) - rrt(019) * density(51) 
  pd(51,51) = pd(51,51) - rrt(019) * density(01) 
  pd(01,01) = pd(01,01) + rrt(020) * density(51) 
  pd(01,51) = pd(01,51) + rrt(020) * density(01) 
  pd(13,01) = pd(13,01) + rrt(020) * density(51) 
  pd(13,51) = pd(13,51) + rrt(020) * density(01) 
  pd(16,01) = pd(16,01) + rrt(020) * density(51) 
  pd(16,51) = pd(16,51) + rrt(020) * density(01) 
  pd(51,01) = pd(51,01) - rrt(020) * density(51) 
  pd(51,51) = pd(51,51) - rrt(020) * density(01) 
  pd(01,01) = pd(01,01) + rrt(021) * density(51) 
  pd(01,51) = pd(01,51) + rrt(021) * density(01) 
  pd(11,01) = pd(11,01) + rrt(021) * density(51) 
  pd(11,51) = pd(11,51) + rrt(021) * density(01) 
  pd(50,01) = pd(50,01) + rrt(021) * density(51) 
  pd(50,51) = pd(50,51) + rrt(021) * density(01) 
  pd(51,01) = pd(51,01) - rrt(021) * density(51) 
  pd(51,51) = pd(51,51) - rrt(021) * density(01) 
  pd(01,01) = pd(01,01) + rrt(022) * density(51) 
  pd(01,51) = pd(01,51) + rrt(022) * density(01) 
  pd(11,01) = pd(11,01) + rrt(022) * density(51) 
  pd(11,51) = pd(11,51) + rrt(022) * density(01) 
  pd(13,01) = pd(13,01) + rrt(022) * density(51) 
  pd(13,51) = pd(13,51) + rrt(022) * density(01) 
  pd(44,01) = pd(44,01) + rrt(022) * density(51) 
  pd(44,51) = pd(44,51) + rrt(022) * density(01) 
  pd(51,01) = pd(51,01) - rrt(022) * density(51) 
  pd(51,51) = pd(51,51) - rrt(022) * density(01) 
  pd(01,01) = pd(01,01) + rrt(023) * density(23) 
  pd(01,23) = pd(01,23) + rrt(023) * density(01) 
  pd(16,01) = pd(16,01) + rrt(023) * density(23) 
  pd(16,23) = pd(16,23) + rrt(023) * density(01) 
  pd(23,01) = pd(23,01) - rrt(023) * density(23) 
  pd(23,23) = pd(23,23) - rrt(023) * density(01) 
  pd(01,01) = pd(01,01) + rrt(024) * density(23) 
  pd(01,23) = pd(01,23) + rrt(024) * density(01) 
  pd(13,01) = pd(13,01) + rrt(024) * density(23) 
  pd(13,23) = pd(13,23) + rrt(024) * density(01) 
  pd(23,01) = pd(23,01) - rrt(024) * density(23) 
  pd(23,23) = pd(23,23) - rrt(024) * density(01) 
  pd(50,01) = pd(50,01) + rrt(024) * density(23) 
  pd(50,23) = pd(50,23) + rrt(024) * density(01) 
  pd(01,01) = pd(01,01) + rrt(025) * density(23) 
  pd(01,23) = pd(01,23) + rrt(025) * density(01) 
  pd(11,01) = pd(11,01) + rrt(025) * density(23) 
  pd(11,23) = pd(11,23) + rrt(025) * density(01) 
  pd(23,01) = pd(23,01) - rrt(025) * density(23) 
  pd(23,23) = pd(23,23) - rrt(025) * density(01) 
  pd(44,01) = pd(44,01) + rrt(025) * density(23) 
  pd(44,23) = pd(44,23) + rrt(025) * density(01) 
  pd(01,01) = pd(01,01) + rrt(026) * density(28) 
  pd(01,28) = pd(01,28) + rrt(026) * density(01) 
  pd(28,01) = pd(28,01) - rrt(026) * density(28) 
  pd(28,28) = pd(28,28) - rrt(026) * density(01) 
  pd(50,01) = pd(50,01) + rrt(026) * density(28) 
  pd(50,28) = pd(50,28) + rrt(026) * density(01) 
  pd(01,01) = pd(01,01) + rrt(027) * density(28) 
  pd(01,28) = pd(01,28) + rrt(027) * density(01) 
  pd(13,01) = pd(13,01) + rrt(027) * density(28) 
  pd(13,28) = pd(13,28) + rrt(027) * density(01) 
  pd(28,01) = pd(28,01) - rrt(027) * density(28) 
  pd(28,28) = pd(28,28) - rrt(027) * density(01) 
  pd(44,01) = pd(44,01) + rrt(027) * density(28) 
  pd(44,28) = pd(44,28) + rrt(027) * density(01) 
  pd(01,01) = pd(01,01) + rrt(028) * density(08) 
  pd(01,08) = pd(01,08) + rrt(028) * density(01) 
  pd(08,01) = pd(08,01) - rrt(028) * density(08) 
  pd(08,08) = pd(08,08) - rrt(028) * density(01) 
  pd(44,01) = pd(44,01) + rrt(028) * density(08) 
  pd(44,08) = pd(44,08) + rrt(028) * density(01) 
  pd(13,01) = pd(13,01) + rrt(029) * density(43) 
  pd(13,43) = pd(13,43) + rrt(029) * density(01) 
  pd(34,01) = pd(34,01) + rrt(029) * density(43) 
  pd(34,43) = pd(34,43) + rrt(029) * density(01) 
  pd(43,01) = pd(43,01) - rrt(029) * density(43) 
  pd(43,43) = pd(43,43) - rrt(029) * density(01) 
  pd(11,01) = pd(11,01) + rrt(030) * density(43) 
  pd(11,43) = pd(11,43) + rrt(030) * density(01) 
  pd(42,01) = pd(42,01) + rrt(030) * density(43) 
  pd(42,43) = pd(42,43) + rrt(030) * density(01) 
  pd(43,01) = pd(43,01) - rrt(030) * density(43) 
  pd(43,43) = pd(43,43) - rrt(030) * density(01) 
  pd(11,01) = pd(11,01) + rrt(031) * density(43) 
  pd(11,43) = pd(11,43) + rrt(031) * density(01) 
  pd(13,01) = pd(13,01) + rrt(031) * density(43) 
  pd(13,43) = pd(13,43) + rrt(031) * density(01) 
  pd(29,01) = pd(29,01) + rrt(031) * density(43) 
  pd(29,43) = pd(29,43) + rrt(031) * density(01) 
  pd(43,01) = pd(43,01) - rrt(031) * density(43) 
  pd(43,43) = pd(43,43) - rrt(031) * density(01) 
  pd(11,01) = pd(11,01) + rrt(032) * density(43) * 2.0d0
  pd(11,43) = pd(11,43) + rrt(032) * density(01) * 2.0d0
  pd(21,01) = pd(21,01) + rrt(032) * density(43) 
  pd(21,43) = pd(21,43) + rrt(032) * density(01) 
  pd(43,01) = pd(43,01) - rrt(032) * density(43) 
  pd(43,43) = pd(43,43) - rrt(032) * density(01) 
  pd(28,01) = pd(28,01) + rrt(033) * density(43) 
  pd(28,43) = pd(28,43) + rrt(033) * density(01) 
  pd(43,01) = pd(43,01) - rrt(033) * density(43) 
  pd(43,43) = pd(43,43) - rrt(033) * density(01) 
  pd(51,01) = pd(51,01) + rrt(033) * density(43) 
  pd(51,43) = pd(51,43) + rrt(033) * density(01) 
  pd(23,01) = pd(23,01) + rrt(034) * density(43) * 2.0d0
  pd(23,43) = pd(23,43) + rrt(034) * density(01) * 2.0d0
  pd(43,01) = pd(43,01) - rrt(034) * density(43) 
  pd(43,43) = pd(43,43) - rrt(034) * density(01) 
  pd(13,01) = pd(13,01) + rrt(035) * density(34) 
  pd(13,34) = pd(13,34) + rrt(035) * density(01) 
  pd(34,01) = pd(34,01) - rrt(035) * density(34) 
  pd(34,34) = pd(34,34) - rrt(035) * density(01) 
  pd(42,01) = pd(42,01) + rrt(035) * density(34) 
  pd(42,34) = pd(42,34) + rrt(035) * density(01) 
  pd(11,01) = pd(11,01) + rrt(036) * density(34) 
  pd(11,34) = pd(11,34) + rrt(036) * density(01) 
  pd(29,01) = pd(29,01) + rrt(036) * density(34) 
  pd(29,34) = pd(29,34) + rrt(036) * density(01) 
  pd(34,01) = pd(34,01) - rrt(036) * density(34) 
  pd(34,34) = pd(34,34) - rrt(036) * density(01) 
  pd(13,01) = pd(13,01) + rrt(037) * density(34) * 2.0d0
  pd(13,34) = pd(13,34) + rrt(037) * density(01) * 2.0d0
  pd(29,01) = pd(29,01) + rrt(037) * density(34) 
  pd(29,34) = pd(29,34) + rrt(037) * density(01) 
  pd(34,01) = pd(34,01) - rrt(037) * density(34) 
  pd(34,34) = pd(34,34) - rrt(037) * density(01) 
  pd(11,01) = pd(11,01) + rrt(038) * density(34) 
  pd(11,34) = pd(11,34) + rrt(038) * density(01) 
  pd(13,01) = pd(13,01) + rrt(038) * density(34) 
  pd(13,34) = pd(13,34) + rrt(038) * density(01) 
  pd(21,01) = pd(21,01) + rrt(038) * density(34) 
  pd(21,34) = pd(21,34) + rrt(038) * density(01) 
  pd(34,01) = pd(34,01) - rrt(038) * density(34) 
  pd(34,34) = pd(34,34) - rrt(038) * density(01) 
  pd(08,01) = pd(08,01) + rrt(039) * density(34) 
  pd(08,34) = pd(08,34) + rrt(039) * density(01) 
  pd(34,01) = pd(34,01) - rrt(039) * density(34) 
  pd(34,34) = pd(34,34) - rrt(039) * density(01) 
  pd(51,01) = pd(51,01) + rrt(039) * density(34) 
  pd(51,34) = pd(51,34) + rrt(039) * density(01) 
  pd(23,01) = pd(23,01) + rrt(040) * density(34) 
  pd(23,34) = pd(23,34) + rrt(040) * density(01) 
  pd(28,01) = pd(28,01) + rrt(040) * density(34) 
  pd(28,34) = pd(28,34) + rrt(040) * density(01) 
  pd(34,01) = pd(34,01) - rrt(040) * density(34) 
  pd(34,34) = pd(34,34) - rrt(040) * density(01) 
  pd(13,01) = pd(13,01) + rrt(041) * density(42) 
  pd(13,42) = pd(13,42) + rrt(041) * density(01) 
  pd(29,01) = pd(29,01) + rrt(041) * density(42) 
  pd(29,42) = pd(29,42) + rrt(041) * density(01) 
  pd(42,01) = pd(42,01) - rrt(041) * density(42) 
  pd(42,42) = pd(42,42) - rrt(041) * density(01) 
  pd(11,01) = pd(11,01) + rrt(042) * density(42) 
  pd(11,42) = pd(11,42) + rrt(042) * density(01) 
  pd(21,01) = pd(21,01) + rrt(042) * density(42) 
  pd(21,42) = pd(21,42) + rrt(042) * density(01) 
  pd(42,01) = pd(42,01) - rrt(042) * density(42) 
  pd(42,42) = pd(42,42) - rrt(042) * density(01) 
  pd(13,01) = pd(13,01) + rrt(043) * density(42) * 2.0d0
  pd(13,42) = pd(13,42) + rrt(043) * density(01) * 2.0d0
  pd(21,01) = pd(21,01) + rrt(043) * density(42) 
  pd(21,42) = pd(21,42) + rrt(043) * density(01) 
  pd(42,01) = pd(42,01) - rrt(043) * density(42) 
  pd(42,42) = pd(42,42) - rrt(043) * density(01) 
  pd(08,01) = pd(08,01) + rrt(044) * density(42) 
  pd(08,42) = pd(08,42) + rrt(044) * density(01) 
  pd(23,01) = pd(23,01) + rrt(044) * density(42) 
  pd(23,42) = pd(23,42) + rrt(044) * density(01) 
  pd(42,01) = pd(42,01) - rrt(044) * density(42) 
  pd(42,42) = pd(42,42) - rrt(044) * density(01) 
  pd(28,01) = pd(28,01) + rrt(045) * density(42) * 2.0d0
  pd(28,42) = pd(28,42) + rrt(045) * density(01) * 2.0d0
  pd(42,01) = pd(42,01) - rrt(045) * density(42) 
  pd(42,42) = pd(42,42) - rrt(045) * density(01) 
  pd(13,01) = pd(13,01) + rrt(046) * density(29) 
  pd(13,29) = pd(13,29) + rrt(046) * density(01) 
  pd(21,01) = pd(21,01) + rrt(046) * density(29) 
  pd(21,29) = pd(21,29) + rrt(046) * density(01) 
  pd(29,01) = pd(29,01) - rrt(046) * density(29) 
  pd(29,29) = pd(29,29) - rrt(046) * density(01) 
  pd(08,01) = pd(08,01) + rrt(047) * density(29) 
  pd(08,29) = pd(08,29) + rrt(047) * density(01) 
  pd(28,01) = pd(28,01) + rrt(047) * density(29) 
  pd(28,29) = pd(28,29) + rrt(047) * density(01) 
  pd(29,01) = pd(29,01) - rrt(047) * density(29) 
  pd(29,29) = pd(29,29) - rrt(047) * density(01) 
  pd(08,01) = pd(08,01) + rrt(048) * density(21) * 2.0d0
  pd(08,21) = pd(08,21) + rrt(048) * density(01) * 2.0d0
  pd(21,01) = pd(21,01) - rrt(048) * density(21) 
  pd(21,21) = pd(21,21) - rrt(048) * density(01) 
  pd(13,01) = pd(13,01) + rrt(049) * density(45) 
  pd(13,45) = pd(13,45) + rrt(049) * density(01) 
  pd(18,01) = pd(18,01) + rrt(049) * density(45) 
  pd(18,45) = pd(18,45) + rrt(049) * density(01) 
  pd(45,01) = pd(45,01) - rrt(049) * density(45) 
  pd(45,45) = pd(45,45) - rrt(049) * density(01) 
  pd(11,01) = pd(11,01) + rrt(050) * density(45) 
  pd(11,45) = pd(11,45) + rrt(050) * density(01) 
  pd(37,01) = pd(37,01) + rrt(050) * density(45) 
  pd(37,45) = pd(37,45) + rrt(050) * density(01) 
  pd(45,01) = pd(45,01) - rrt(050) * density(45) 
  pd(45,45) = pd(45,45) - rrt(050) * density(01) 
  pd(11,01) = pd(11,01) + rrt(051) * density(45) * 2.0d0
  pd(11,45) = pd(11,45) + rrt(051) * density(01) * 2.0d0
  pd(35,01) = pd(35,01) + rrt(051) * density(45) 
  pd(35,45) = pd(35,45) + rrt(051) * density(01) 
  pd(45,01) = pd(45,01) - rrt(051) * density(45) 
  pd(45,45) = pd(45,45) - rrt(051) * density(01) 
  pd(28,01) = pd(28,01) + rrt(052) * density(45) 
  pd(28,45) = pd(28,45) + rrt(052) * density(01) 
  pd(43,01) = pd(43,01) + rrt(052) * density(45) 
  pd(43,45) = pd(43,45) + rrt(052) * density(01) 
  pd(45,01) = pd(45,01) - rrt(052) * density(45) 
  pd(45,45) = pd(45,45) - rrt(052) * density(01) 
  pd(23,01) = pd(23,01) + rrt(053) * density(45) 
  pd(23,45) = pd(23,45) + rrt(053) * density(01) 
  pd(34,01) = pd(34,01) + rrt(053) * density(45) 
  pd(34,45) = pd(34,45) + rrt(053) * density(01) 
  pd(45,01) = pd(45,01) - rrt(053) * density(45) 
  pd(45,45) = pd(45,45) - rrt(053) * density(01) 
  pd(42,01) = pd(42,01) + rrt(054) * density(45) 
  pd(42,45) = pd(42,45) + rrt(054) * density(01) 
  pd(45,01) = pd(45,01) - rrt(054) * density(45) 
  pd(45,45) = pd(45,45) - rrt(054) * density(01) 
  pd(51,01) = pd(51,01) + rrt(054) * density(45) 
  pd(51,45) = pd(51,45) + rrt(054) * density(01) 
  pd(13,01) = pd(13,01) + rrt(055) * density(18) 
  pd(13,18) = pd(13,18) + rrt(055) * density(01) 
  pd(18,01) = pd(18,01) - rrt(055) * density(18) 
  pd(18,18) = pd(18,18) - rrt(055) * density(01) 
  pd(37,01) = pd(37,01) + rrt(055) * density(18) 
  pd(37,18) = pd(37,18) + rrt(055) * density(01) 
  pd(11,01) = pd(11,01) + rrt(056) * density(18) 
  pd(11,18) = pd(11,18) + rrt(056) * density(01) 
  pd(18,01) = pd(18,01) - rrt(056) * density(18) 
  pd(18,18) = pd(18,18) - rrt(056) * density(01) 
  pd(20,01) = pd(20,01) + rrt(056) * density(18) 
  pd(20,18) = pd(20,18) + rrt(056) * density(01) 
  pd(11,01) = pd(11,01) + rrt(057) * density(18) 
  pd(11,18) = pd(11,18) + rrt(057) * density(01) 
  pd(13,01) = pd(13,01) + rrt(057) * density(18) 
  pd(13,18) = pd(13,18) + rrt(057) * density(01) 
  pd(18,01) = pd(18,01) - rrt(057) * density(18) 
  pd(18,18) = pd(18,18) - rrt(057) * density(01) 
  pd(35,01) = pd(35,01) + rrt(057) * density(18) 
  pd(35,18) = pd(35,18) + rrt(057) * density(01) 
  pd(18,01) = pd(18,01) - rrt(058) * density(18) 
  pd(18,18) = pd(18,18) - rrt(058) * density(01) 
  pd(23,01) = pd(23,01) + rrt(058) * density(18) 
  pd(23,18) = pd(23,18) + rrt(058) * density(01) 
  pd(42,01) = pd(42,01) + rrt(058) * density(18) 
  pd(42,18) = pd(42,18) + rrt(058) * density(01) 
  pd(18,01) = pd(18,01) - rrt(059) * density(18) 
  pd(18,18) = pd(18,18) - rrt(059) * density(01) 
  pd(29,01) = pd(29,01) + rrt(059) * density(18) 
  pd(29,18) = pd(29,18) + rrt(059) * density(01) 
  pd(51,01) = pd(51,01) + rrt(059) * density(18) 
  pd(51,18) = pd(51,18) + rrt(059) * density(01) 
  pd(13,01) = pd(13,01) + rrt(060) * density(37) 
  pd(13,37) = pd(13,37) + rrt(060) * density(01) 
  pd(20,01) = pd(20,01) + rrt(060) * density(37) 
  pd(20,37) = pd(20,37) + rrt(060) * density(01) 
  pd(37,01) = pd(37,01) - rrt(060) * density(37) 
  pd(37,37) = pd(37,37) - rrt(060) * density(01) 
  pd(11,01) = pd(11,01) + rrt(061) * density(37) 
  pd(11,37) = pd(11,37) + rrt(061) * density(01) 
  pd(35,01) = pd(35,01) + rrt(061) * density(37) 
  pd(35,37) = pd(35,37) + rrt(061) * density(01) 
  pd(37,01) = pd(37,01) - rrt(061) * density(37) 
  pd(37,37) = pd(37,37) - rrt(061) * density(01) 
  pd(28,01) = pd(28,01) + rrt(062) * density(37) 
  pd(28,37) = pd(28,37) + rrt(062) * density(01) 
  pd(37,01) = pd(37,01) - rrt(062) * density(37) 
  pd(37,37) = pd(37,37) - rrt(062) * density(01) 
  pd(42,01) = pd(42,01) + rrt(062) * density(37) 
  pd(42,37) = pd(42,37) + rrt(062) * density(01) 
  pd(23,01) = pd(23,01) + rrt(063) * density(37) 
  pd(23,37) = pd(23,37) + rrt(063) * density(01) 
  pd(29,01) = pd(29,01) + rrt(063) * density(37) 
  pd(29,37) = pd(29,37) + rrt(063) * density(01) 
  pd(37,01) = pd(37,01) - rrt(063) * density(37) 
  pd(37,37) = pd(37,37) - rrt(063) * density(01) 
  pd(21,01) = pd(21,01) + rrt(064) * density(37) 
  pd(21,37) = pd(21,37) + rrt(064) * density(01) 
  pd(37,01) = pd(37,01) - rrt(064) * density(37) 
  pd(37,37) = pd(37,37) - rrt(064) * density(01) 
  pd(51,01) = pd(51,01) + rrt(064) * density(37) 
  pd(51,37) = pd(51,37) + rrt(064) * density(01) 
  pd(13,01) = pd(13,01) + rrt(065) * density(20) 
  pd(13,20) = pd(13,20) + rrt(065) * density(01) 
  pd(20,01) = pd(20,01) - rrt(065) * density(20) 
  pd(20,20) = pd(20,20) - rrt(065) * density(01) 
  pd(35,01) = pd(35,01) + rrt(065) * density(20) 
  pd(35,20) = pd(35,20) + rrt(065) * density(01) 
  pd(20,01) = pd(20,01) - rrt(066) * density(20) 
  pd(20,20) = pd(20,20) - rrt(066) * density(01) 
  pd(21,01) = pd(21,01) + rrt(066) * density(20) 
  pd(21,20) = pd(21,20) + rrt(066) * density(01) 
  pd(23,01) = pd(23,01) + rrt(066) * density(20) 
  pd(23,20) = pd(23,20) + rrt(066) * density(01) 
  pd(08,01) = pd(08,01) + rrt(067) * density(35) 
  pd(08,35) = pd(08,35) + rrt(067) * density(01) 
  pd(29,01) = pd(29,01) + rrt(067) * density(35) 
  pd(29,35) = pd(29,35) + rrt(067) * density(01) 
  pd(35,01) = pd(35,01) - rrt(067) * density(35) 
  pd(35,35) = pd(35,35) - rrt(067) * density(01) 
  pd(21,01) = pd(21,01) + rrt(068) * density(35) 
  pd(21,35) = pd(21,35) + rrt(068) * density(01) 
  pd(28,01) = pd(28,01) + rrt(068) * density(35) 
  pd(28,35) = pd(28,35) + rrt(068) * density(01) 
  pd(35,01) = pd(35,01) - rrt(068) * density(35) 
  pd(35,35) = pd(35,35) - rrt(068) * density(01) 
  pd(01,01) = pd(01,01) + rrt(069) * density(43) 
  pd(01,43) = pd(01,43) + rrt(069) * density(01) 
  pd(40,01) = pd(40,01) + rrt(069) * density(43) 
  pd(40,43) = pd(40,43) + rrt(069) * density(01) 
  pd(43,01) = pd(43,01) - rrt(069) * density(43) 
  pd(43,43) = pd(43,43) - rrt(069) * density(01) 
  pd(01,01) = pd(01,01) + rrt(070) * density(43) 
  pd(01,43) = pd(01,43) + rrt(070) * density(01) 
  pd(04,01) = pd(04,01) + rrt(070) * density(43) 
  pd(04,43) = pd(04,43) + rrt(070) * density(01) 
  pd(13,01) = pd(13,01) + rrt(070) * density(43) 
  pd(13,43) = pd(13,43) + rrt(070) * density(01) 
  pd(43,01) = pd(43,01) - rrt(070) * density(43) 
  pd(43,43) = pd(43,43) - rrt(070) * density(01) 
  pd(01,01) = pd(01,01) + rrt(071) * density(43) 
  pd(01,43) = pd(01,43) + rrt(071) * density(01) 
  pd(11,01) = pd(11,01) + rrt(071) * density(43) 
  pd(11,43) = pd(11,43) + rrt(071) * density(01) 
  pd(43,01) = pd(43,01) - rrt(071) * density(43) 
  pd(43,43) = pd(43,43) - rrt(071) * density(01) 
  pd(46,01) = pd(46,01) + rrt(071) * density(43) 
  pd(46,43) = pd(46,43) + rrt(071) * density(01) 
  pd(01,01) = pd(01,01) + rrt(072) * density(43) 
  pd(01,43) = pd(01,43) + rrt(072) * density(01) 
  pd(07,01) = pd(07,01) + rrt(072) * density(43) 
  pd(07,43) = pd(07,43) + rrt(072) * density(01) 
  pd(11,01) = pd(11,01) + rrt(072) * density(43) 
  pd(11,43) = pd(11,43) + rrt(072) * density(01) 
  pd(13,01) = pd(13,01) + rrt(072) * density(43) 
  pd(13,43) = pd(13,43) + rrt(072) * density(01) 
  pd(43,01) = pd(43,01) - rrt(072) * density(43) 
  pd(43,43) = pd(43,43) - rrt(072) * density(01) 
  pd(01,01) = pd(01,01) + rrt(073) * density(43) 
  pd(01,43) = pd(01,43) + rrt(073) * density(01) 
  pd(11,01) = pd(11,01) + rrt(073) * density(43) * 2.0d0
  pd(11,43) = pd(11,43) + rrt(073) * density(01) * 2.0d0
  pd(17,01) = pd(17,01) + rrt(073) * density(43) 
  pd(17,43) = pd(17,43) + rrt(073) * density(01) 
  pd(43,01) = pd(43,01) - rrt(073) * density(43) 
  pd(43,43) = pd(43,43) - rrt(073) * density(01) 
  pd(01,01) = pd(01,01) + rrt(074) * density(43) 
  pd(01,43) = pd(01,43) + rrt(074) * density(01) 
  pd(16,01) = pd(16,01) + rrt(074) * density(43) 
  pd(16,43) = pd(16,43) + rrt(074) * density(01) 
  pd(23,01) = pd(23,01) + rrt(074) * density(43) 
  pd(23,43) = pd(23,43) + rrt(074) * density(01) 
  pd(43,01) = pd(43,01) - rrt(074) * density(43) 
  pd(43,43) = pd(43,43) - rrt(074) * density(01) 
  pd(01,01) = pd(01,01) + rrt(075) * density(43) 
  pd(01,43) = pd(01,43) + rrt(075) * density(01) 
  pd(43,01) = pd(43,01) - rrt(075) * density(43) 
  pd(43,43) = pd(43,43) - rrt(075) * density(01) 
  pd(50,01) = pd(50,01) + rrt(075) * density(43) 
  pd(50,43) = pd(50,43) + rrt(075) * density(01) 
  pd(51,01) = pd(51,01) + rrt(075) * density(43) 
  pd(51,43) = pd(51,43) + rrt(075) * density(01) 
  pd(01,01) = pd(01,01) + rrt(076) * density(34) 
  pd(01,34) = pd(01,34) + rrt(076) * density(01) 
  pd(04,01) = pd(04,01) + rrt(076) * density(34) 
  pd(04,34) = pd(04,34) + rrt(076) * density(01) 
  pd(34,01) = pd(34,01) - rrt(076) * density(34) 
  pd(34,34) = pd(34,34) - rrt(076) * density(01) 
  pd(01,01) = pd(01,01) + rrt(077) * density(34) 
  pd(01,34) = pd(01,34) + rrt(077) * density(01) 
  pd(13,01) = pd(13,01) + rrt(077) * density(34) 
  pd(13,34) = pd(13,34) + rrt(077) * density(01) 
  pd(34,01) = pd(34,01) - rrt(077) * density(34) 
  pd(34,34) = pd(34,34) - rrt(077) * density(01) 
  pd(46,01) = pd(46,01) + rrt(077) * density(34) 
  pd(46,34) = pd(46,34) + rrt(077) * density(01) 
  pd(01,01) = pd(01,01) + rrt(078) * density(34) 
  pd(01,34) = pd(01,34) + rrt(078) * density(01) 
  pd(07,01) = pd(07,01) + rrt(078) * density(34) 
  pd(07,34) = pd(07,34) + rrt(078) * density(01) 
  pd(11,01) = pd(11,01) + rrt(078) * density(34) 
  pd(11,34) = pd(11,34) + rrt(078) * density(01) 
  pd(34,01) = pd(34,01) - rrt(078) * density(34) 
  pd(34,34) = pd(34,34) - rrt(078) * density(01) 
  pd(01,01) = pd(01,01) + rrt(079) * density(34) 
  pd(01,34) = pd(01,34) + rrt(079) * density(01) 
  pd(11,01) = pd(11,01) + rrt(079) * density(34) 
  pd(11,34) = pd(11,34) + rrt(079) * density(01) 
  pd(13,01) = pd(13,01) + rrt(079) * density(34) 
  pd(13,34) = pd(13,34) + rrt(079) * density(01) 
  pd(17,01) = pd(17,01) + rrt(079) * density(34) 
  pd(17,34) = pd(17,34) + rrt(079) * density(01) 
  pd(34,01) = pd(34,01) - rrt(079) * density(34) 
  pd(34,34) = pd(34,34) - rrt(079) * density(01) 
  pd(01,01) = pd(01,01) + rrt(080) * density(34) 
  pd(01,34) = pd(01,34) + rrt(080) * density(01) 
  pd(16,01) = pd(16,01) + rrt(080) * density(34) 
  pd(16,34) = pd(16,34) + rrt(080) * density(01) 
  pd(28,01) = pd(28,01) + rrt(080) * density(34) 
  pd(28,34) = pd(28,34) + rrt(080) * density(01) 
  pd(34,01) = pd(34,01) - rrt(080) * density(34) 
  pd(34,34) = pd(34,34) - rrt(080) * density(01) 
  pd(01,01) = pd(01,01) + rrt(081) * density(34) 
  pd(01,34) = pd(01,34) + rrt(081) * density(01) 
  pd(23,01) = pd(23,01) + rrt(081) * density(34) 
  pd(23,34) = pd(23,34) + rrt(081) * density(01) 
  pd(34,01) = pd(34,01) - rrt(081) * density(34) 
  pd(34,34) = pd(34,34) - rrt(081) * density(01) 
  pd(50,01) = pd(50,01) + rrt(081) * density(34) 
  pd(50,34) = pd(50,34) + rrt(081) * density(01) 
  pd(01,01) = pd(01,01) + rrt(082) * density(34) 
  pd(01,34) = pd(01,34) + rrt(082) * density(01) 
  pd(34,01) = pd(34,01) - rrt(082) * density(34) 
  pd(34,34) = pd(34,34) - rrt(082) * density(01) 
  pd(44,01) = pd(44,01) + rrt(082) * density(34) 
  pd(44,34) = pd(44,34) + rrt(082) * density(01) 
  pd(51,01) = pd(51,01) + rrt(082) * density(34) 
  pd(51,34) = pd(51,34) + rrt(082) * density(01) 
  pd(01,01) = pd(01,01) + rrt(083) * density(42) 
  pd(01,42) = pd(01,42) + rrt(083) * density(01) 
  pd(42,01) = pd(42,01) - rrt(083) * density(42) 
  pd(42,42) = pd(42,42) - rrt(083) * density(01) 
  pd(46,01) = pd(46,01) + rrt(083) * density(42) 
  pd(46,42) = pd(46,42) + rrt(083) * density(01) 
  pd(01,01) = pd(01,01) + rrt(084) * density(42) 
  pd(01,42) = pd(01,42) + rrt(084) * density(01) 
  pd(07,01) = pd(07,01) + rrt(084) * density(42) 
  pd(07,42) = pd(07,42) + rrt(084) * density(01) 
  pd(13,01) = pd(13,01) + rrt(084) * density(42) 
  pd(13,42) = pd(13,42) + rrt(084) * density(01) 
  pd(42,01) = pd(42,01) - rrt(084) * density(42) 
  pd(42,42) = pd(42,42) - rrt(084) * density(01) 
  pd(01,01) = pd(01,01) + rrt(085) * density(42) 
  pd(01,42) = pd(01,42) + rrt(085) * density(01) 
  pd(08,01) = pd(08,01) + rrt(085) * density(42) 
  pd(08,42) = pd(08,42) + rrt(085) * density(01) 
  pd(16,01) = pd(16,01) + rrt(085) * density(42) 
  pd(16,42) = pd(16,42) + rrt(085) * density(01) 
  pd(42,01) = pd(42,01) - rrt(085) * density(42) 
  pd(42,42) = pd(42,42) - rrt(085) * density(01) 
  pd(01,01) = pd(01,01) + rrt(086) * density(42) 
  pd(01,42) = pd(01,42) + rrt(086) * density(01) 
  pd(28,01) = pd(28,01) + rrt(086) * density(42) 
  pd(28,42) = pd(28,42) + rrt(086) * density(01) 
  pd(42,01) = pd(42,01) - rrt(086) * density(42) 
  pd(42,42) = pd(42,42) - rrt(086) * density(01) 
  pd(50,01) = pd(50,01) + rrt(086) * density(42) 
  pd(50,42) = pd(50,42) + rrt(086) * density(01) 
  pd(01,01) = pd(01,01) + rrt(087) * density(42) 
  pd(01,42) = pd(01,42) + rrt(087) * density(01) 
  pd(23,01) = pd(23,01) + rrt(087) * density(42) 
  pd(23,42) = pd(23,42) + rrt(087) * density(01) 
  pd(42,01) = pd(42,01) - rrt(087) * density(42) 
  pd(42,42) = pd(42,42) - rrt(087) * density(01) 
  pd(44,01) = pd(44,01) + rrt(087) * density(42) 
  pd(44,42) = pd(44,42) + rrt(087) * density(01) 
  pd(01,01) = pd(01,01) + rrt(088) * density(29) 
  pd(01,29) = pd(01,29) + rrt(088) * density(01) 
  pd(07,01) = pd(07,01) + rrt(088) * density(29) 
  pd(07,29) = pd(07,29) + rrt(088) * density(01) 
  pd(29,01) = pd(29,01) - rrt(088) * density(29) 
  pd(29,29) = pd(29,29) - rrt(088) * density(01) 
  pd(01,01) = pd(01,01) + rrt(089) * density(29) 
  pd(01,29) = pd(01,29) + rrt(089) * density(01) 
  pd(13,01) = pd(13,01) + rrt(089) * density(29) 
  pd(13,29) = pd(13,29) + rrt(089) * density(01) 
  pd(17,01) = pd(17,01) + rrt(089) * density(29) 
  pd(17,29) = pd(17,29) + rrt(089) * density(01) 
  pd(29,01) = pd(29,01) - rrt(089) * density(29) 
  pd(29,29) = pd(29,29) - rrt(089) * density(01) 
  pd(01,01) = pd(01,01) + rrt(090) * density(29) 
  pd(01,29) = pd(01,29) + rrt(090) * density(01) 
  pd(08,01) = pd(08,01) + rrt(090) * density(29) 
  pd(08,29) = pd(08,29) + rrt(090) * density(01) 
  pd(29,01) = pd(29,01) - rrt(090) * density(29) 
  pd(29,29) = pd(29,29) - rrt(090) * density(01) 
  pd(50,01) = pd(50,01) + rrt(090) * density(29) 
  pd(50,29) = pd(50,29) + rrt(090) * density(01) 
  pd(01,01) = pd(01,01) + rrt(091) * density(29) 
  pd(01,29) = pd(01,29) + rrt(091) * density(01) 
  pd(28,01) = pd(28,01) + rrt(091) * density(29) 
  pd(28,29) = pd(28,29) + rrt(091) * density(01) 
  pd(29,01) = pd(29,01) - rrt(091) * density(29) 
  pd(29,29) = pd(29,29) - rrt(091) * density(01) 
  pd(44,01) = pd(44,01) + rrt(091) * density(29) 
  pd(44,29) = pd(44,29) + rrt(091) * density(01) 
  pd(01,01) = pd(01,01) + rrt(092) * density(29) 
  pd(01,29) = pd(01,29) + rrt(092) * density(01) 
  pd(21,01) = pd(21,01) + rrt(092) * density(29) 
  pd(21,29) = pd(21,29) + rrt(092) * density(01) 
  pd(29,01) = pd(29,01) - rrt(092) * density(29) 
  pd(29,29) = pd(29,29) - rrt(092) * density(01) 
  pd(30,01) = pd(30,01) + rrt(092) * density(29) 
  pd(30,29) = pd(30,29) + rrt(092) * density(01) 
  pd(01,01) = pd(01,01) + rrt(093) * density(21) 
  pd(01,21) = pd(01,21) + rrt(093) * density(01) 
  pd(17,01) = pd(17,01) + rrt(093) * density(21) 
  pd(17,21) = pd(17,21) + rrt(093) * density(01) 
  pd(21,01) = pd(21,01) - rrt(093) * density(21) 
  pd(21,21) = pd(21,21) - rrt(093) * density(01) 
  pd(01,01) = pd(01,01) + rrt(094) * density(21) 
  pd(01,21) = pd(01,21) + rrt(094) * density(01) 
  pd(08,01) = pd(08,01) + rrt(094) * density(21) 
  pd(08,21) = pd(08,21) + rrt(094) * density(01) 
  pd(21,01) = pd(21,01) - rrt(094) * density(21) 
  pd(21,21) = pd(21,21) - rrt(094) * density(01) 
  pd(44,01) = pd(44,01) + rrt(094) * density(21) 
  pd(44,21) = pd(44,21) + rrt(094) * density(01) 
  pd(01,01) = pd(01,01) + rrt(095) * density(45) 
  pd(01,45) = pd(01,45) + rrt(095) * density(01) 
  pd(14,01) = pd(14,01) + rrt(095) * density(45) 
  pd(14,45) = pd(14,45) + rrt(095) * density(01) 
  pd(45,01) = pd(45,01) - rrt(095) * density(45) 
  pd(45,45) = pd(45,45) - rrt(095) * density(01) 
  pd(01,01) = pd(01,01) + rrt(096) * density(45) 
  pd(01,45) = pd(01,45) + rrt(096) * density(01) 
  pd(03,01) = pd(03,01) + rrt(096) * density(45) 
  pd(03,45) = pd(03,45) + rrt(096) * density(01) 
  pd(13,01) = pd(13,01) + rrt(096) * density(45) 
  pd(13,45) = pd(13,45) + rrt(096) * density(01) 
  pd(45,01) = pd(45,01) - rrt(096) * density(45) 
  pd(45,45) = pd(45,45) - rrt(096) * density(01) 
  pd(01,01) = pd(01,01) + rrt(097) * density(45) 
  pd(01,45) = pd(01,45) + rrt(097) * density(01) 
  pd(05,01) = pd(05,01) + rrt(097) * density(45) 
  pd(05,45) = pd(05,45) + rrt(097) * density(01) 
  pd(11,01) = pd(11,01) + rrt(097) * density(45) 
  pd(11,45) = pd(11,45) + rrt(097) * density(01) 
  pd(45,01) = pd(45,01) - rrt(097) * density(45) 
  pd(45,45) = pd(45,45) - rrt(097) * density(01) 
  pd(01,01) = pd(01,01) + rrt(098) * density(45) 
  pd(01,45) = pd(01,45) + rrt(098) * density(01) 
  pd(11,01) = pd(11,01) + rrt(098) * density(45) 
  pd(11,45) = pd(11,45) + rrt(098) * density(01) 
  pd(13,01) = pd(13,01) + rrt(098) * density(45) 
  pd(13,45) = pd(13,45) + rrt(098) * density(01) 
  pd(26,01) = pd(26,01) + rrt(098) * density(45) 
  pd(26,45) = pd(26,45) + rrt(098) * density(01) 
  pd(45,01) = pd(45,01) - rrt(098) * density(45) 
  pd(45,45) = pd(45,45) - rrt(098) * density(01) 
  pd(01,01) = pd(01,01) + rrt(099) * density(45) 
  pd(01,45) = pd(01,45) + rrt(099) * density(01) 
  pd(11,01) = pd(11,01) + rrt(099) * density(45) * 2.0d0
  pd(11,45) = pd(11,45) + rrt(099) * density(01) * 2.0d0
  pd(19,01) = pd(19,01) + rrt(099) * density(45) 
  pd(19,45) = pd(19,45) + rrt(099) * density(01) 
  pd(45,01) = pd(45,01) - rrt(099) * density(45) 
  pd(45,45) = pd(45,45) - rrt(099) * density(01) 
  pd(01,01) = pd(01,01) + rrt(100) * density(45) 
  pd(01,45) = pd(01,45) + rrt(100) * density(01) 
  pd(04,01) = pd(04,01) + rrt(100) * density(45) 
  pd(04,45) = pd(04,45) + rrt(100) * density(01) 
  pd(23,01) = pd(23,01) + rrt(100) * density(45) 
  pd(23,45) = pd(23,45) + rrt(100) * density(01) 
  pd(45,01) = pd(45,01) - rrt(100) * density(45) 
  pd(45,45) = pd(45,45) - rrt(100) * density(01) 
  pd(01,01) = pd(01,01) + rrt(101) * density(45) 
  pd(01,45) = pd(01,45) + rrt(101) * density(01) 
  pd(45,01) = pd(45,01) - rrt(101) * density(45) 
  pd(45,45) = pd(45,45) - rrt(101) * density(01) 
  pd(46,01) = pd(46,01) + rrt(101) * density(45) 
  pd(46,45) = pd(46,45) + rrt(101) * density(01) 
  pd(51,01) = pd(51,01) + rrt(101) * density(45) 
  pd(51,45) = pd(51,45) + rrt(101) * density(01) 
  pd(01,01) = pd(01,01) + rrt(102) * density(45) 
  pd(01,45) = pd(01,45) + rrt(102) * density(01) 
  pd(16,01) = pd(16,01) + rrt(102) * density(45) 
  pd(16,45) = pd(16,45) + rrt(102) * density(01) 
  pd(34,01) = pd(34,01) + rrt(102) * density(45) 
  pd(34,45) = pd(34,45) + rrt(102) * density(01) 
  pd(45,01) = pd(45,01) - rrt(102) * density(45) 
  pd(45,45) = pd(45,45) - rrt(102) * density(01) 
  pd(01,01) = pd(01,01) + rrt(103) * density(45) 
  pd(01,45) = pd(01,45) + rrt(103) * density(01) 
  pd(43,01) = pd(43,01) + rrt(103) * density(45) 
  pd(43,45) = pd(43,45) + rrt(103) * density(01) 
  pd(45,01) = pd(45,01) - rrt(103) * density(45) 
  pd(45,45) = pd(45,45) - rrt(103) * density(01) 
  pd(50,01) = pd(50,01) + rrt(103) * density(45) 
  pd(50,45) = pd(50,45) + rrt(103) * density(01) 
  pd(01,01) = pd(01,01) + rrt(104) * density(18) 
  pd(01,18) = pd(01,18) + rrt(104) * density(01) 
  pd(03,01) = pd(03,01) + rrt(104) * density(18) 
  pd(03,18) = pd(03,18) + rrt(104) * density(01) 
  pd(18,01) = pd(18,01) - rrt(104) * density(18) 
  pd(18,18) = pd(18,18) - rrt(104) * density(01) 
  pd(01,01) = pd(01,01) + rrt(105) * density(18) 
  pd(01,18) = pd(01,18) + rrt(105) * density(01) 
  pd(05,01) = pd(05,01) + rrt(105) * density(18) 
  pd(05,18) = pd(05,18) + rrt(105) * density(01) 
  pd(13,01) = pd(13,01) + rrt(105) * density(18) 
  pd(13,18) = pd(13,18) + rrt(105) * density(01) 
  pd(18,01) = pd(18,01) - rrt(105) * density(18) 
  pd(18,18) = pd(18,18) - rrt(105) * density(01) 
  pd(01,01) = pd(01,01) + rrt(106) * density(18) 
  pd(01,18) = pd(01,18) + rrt(106) * density(01) 
  pd(11,01) = pd(11,01) + rrt(106) * density(18) 
  pd(11,18) = pd(11,18) + rrt(106) * density(01) 
  pd(18,01) = pd(18,01) - rrt(106) * density(18) 
  pd(18,18) = pd(18,18) - rrt(106) * density(01) 
  pd(26,01) = pd(26,01) + rrt(106) * density(18) 
  pd(26,18) = pd(26,18) + rrt(106) * density(01) 
  pd(01,01) = pd(01,01) + rrt(107) * density(18) 
  pd(01,18) = pd(01,18) + rrt(107) * density(01) 
  pd(11,01) = pd(11,01) + rrt(107) * density(18) 
  pd(11,18) = pd(11,18) + rrt(107) * density(01) 
  pd(13,01) = pd(13,01) + rrt(107) * density(18) 
  pd(13,18) = pd(13,18) + rrt(107) * density(01) 
  pd(18,01) = pd(18,01) - rrt(107) * density(18) 
  pd(18,18) = pd(18,18) - rrt(107) * density(01) 
  pd(19,01) = pd(19,01) + rrt(107) * density(18) 
  pd(19,18) = pd(19,18) + rrt(107) * density(01) 
  pd(01,01) = pd(01,01) + rrt(108) * density(18) 
  pd(01,18) = pd(01,18) + rrt(108) * density(01) 
  pd(04,01) = pd(04,01) + rrt(108) * density(18) 
  pd(04,18) = pd(04,18) + rrt(108) * density(01) 
  pd(18,01) = pd(18,01) - rrt(108) * density(18) 
  pd(18,18) = pd(18,18) - rrt(108) * density(01) 
  pd(28,01) = pd(28,01) + rrt(108) * density(18) 
  pd(28,18) = pd(28,18) + rrt(108) * density(01) 
  pd(01,01) = pd(01,01) + rrt(109) * density(18) 
  pd(01,18) = pd(01,18) + rrt(109) * density(01) 
  pd(18,01) = pd(18,01) - rrt(109) * density(18) 
  pd(18,18) = pd(18,18) - rrt(109) * density(01) 
  pd(23,01) = pd(23,01) + rrt(109) * density(18) 
  pd(23,18) = pd(23,18) + rrt(109) * density(01) 
  pd(46,01) = pd(46,01) + rrt(109) * density(18) 
  pd(46,18) = pd(46,18) + rrt(109) * density(01) 
  pd(01,01) = pd(01,01) + rrt(110) * density(18) 
  pd(01,18) = pd(01,18) + rrt(110) * density(01) 
  pd(07,01) = pd(07,01) + rrt(110) * density(18) 
  pd(07,18) = pd(07,18) + rrt(110) * density(01) 
  pd(18,01) = pd(18,01) - rrt(110) * density(18) 
  pd(18,18) = pd(18,18) - rrt(110) * density(01) 
  pd(51,01) = pd(51,01) + rrt(110) * density(18) 
  pd(51,18) = pd(51,18) + rrt(110) * density(01) 
  pd(01,01) = pd(01,01) + rrt(111) * density(18) 
  pd(01,18) = pd(01,18) + rrt(111) * density(01) 
  pd(10,01) = pd(10,01) + rrt(111) * density(18) 
  pd(10,18) = pd(10,18) + rrt(111) * density(01) 
  pd(18,01) = pd(18,01) - rrt(111) * density(18) 
  pd(18,18) = pd(18,18) - rrt(111) * density(01) 
  pd(29,01) = pd(29,01) + rrt(111) * density(18) 
  pd(29,18) = pd(29,18) + rrt(111) * density(01) 
  pd(01,01) = pd(01,01) + rrt(112) * density(18) 
  pd(01,18) = pd(01,18) + rrt(112) * density(01) 
  pd(16,01) = pd(16,01) + rrt(112) * density(18) 
  pd(16,18) = pd(16,18) + rrt(112) * density(01) 
  pd(18,01) = pd(18,01) - rrt(112) * density(18) 
  pd(18,18) = pd(18,18) - rrt(112) * density(01) 
  pd(42,01) = pd(42,01) + rrt(112) * density(18) 
  pd(42,18) = pd(42,18) + rrt(112) * density(01) 
  pd(01,01) = pd(01,01) + rrt(113) * density(18) 
  pd(01,18) = pd(01,18) + rrt(113) * density(01) 
  pd(18,01) = pd(18,01) - rrt(113) * density(18) 
  pd(18,18) = pd(18,18) - rrt(113) * density(01) 
  pd(34,01) = pd(34,01) + rrt(113) * density(18) 
  pd(34,18) = pd(34,18) + rrt(113) * density(01) 
  pd(50,01) = pd(50,01) + rrt(113) * density(18) 
  pd(50,18) = pd(50,18) + rrt(113) * density(01) 
  pd(01,01) = pd(01,01) + rrt(114) * density(18) 
  pd(01,18) = pd(01,18) + rrt(114) * density(01) 
  pd(18,01) = pd(18,01) - rrt(114) * density(18) 
  pd(18,18) = pd(18,18) - rrt(114) * density(01) 
  pd(43,01) = pd(43,01) + rrt(114) * density(18) 
  pd(43,18) = pd(43,18) + rrt(114) * density(01) 
  pd(44,01) = pd(44,01) + rrt(114) * density(18) 
  pd(44,18) = pd(44,18) + rrt(114) * density(01) 
  pd(01,01) = pd(01,01) + rrt(115) * density(37) 
  pd(01,37) = pd(01,37) + rrt(115) * density(01) 
  pd(05,01) = pd(05,01) + rrt(115) * density(37) 
  pd(05,37) = pd(05,37) + rrt(115) * density(01) 
  pd(37,01) = pd(37,01) - rrt(115) * density(37) 
  pd(37,37) = pd(37,37) - rrt(115) * density(01) 
  pd(01,01) = pd(01,01) + rrt(116) * density(37) 
  pd(01,37) = pd(01,37) + rrt(116) * density(01) 
  pd(13,01) = pd(13,01) + rrt(116) * density(37) 
  pd(13,37) = pd(13,37) + rrt(116) * density(01) 
  pd(26,01) = pd(26,01) + rrt(116) * density(37) 
  pd(26,37) = pd(26,37) + rrt(116) * density(01) 
  pd(37,01) = pd(37,01) - rrt(116) * density(37) 
  pd(37,37) = pd(37,37) - rrt(116) * density(01) 
  pd(01,01) = pd(01,01) + rrt(117) * density(37) 
  pd(01,37) = pd(01,37) + rrt(117) * density(01) 
  pd(11,01) = pd(11,01) + rrt(117) * density(37) 
  pd(11,37) = pd(11,37) + rrt(117) * density(01) 
  pd(19,01) = pd(19,01) + rrt(117) * density(37) 
  pd(19,37) = pd(19,37) + rrt(117) * density(01) 
  pd(37,01) = pd(37,01) - rrt(117) * density(37) 
  pd(37,37) = pd(37,37) - rrt(117) * density(01) 
  pd(01,01) = pd(01,01) + rrt(118) * density(37) 
  pd(01,37) = pd(01,37) + rrt(118) * density(01) 
  pd(04,01) = pd(04,01) + rrt(118) * density(37) 
  pd(04,37) = pd(04,37) + rrt(118) * density(01) 
  pd(08,01) = pd(08,01) + rrt(118) * density(37) 
  pd(08,37) = pd(08,37) + rrt(118) * density(01) 
  pd(37,01) = pd(37,01) - rrt(118) * density(37) 
  pd(37,37) = pd(37,37) - rrt(118) * density(01) 
  pd(01,01) = pd(01,01) + rrt(119) * density(37) 
  pd(01,37) = pd(01,37) + rrt(119) * density(01) 
  pd(28,01) = pd(28,01) + rrt(119) * density(37) 
  pd(28,37) = pd(28,37) + rrt(119) * density(01) 
  pd(37,01) = pd(37,01) - rrt(119) * density(37) 
  pd(37,37) = pd(37,37) - rrt(119) * density(01) 
  pd(46,01) = pd(46,01) + rrt(119) * density(37) 
  pd(46,37) = pd(46,37) + rrt(119) * density(01) 
  pd(01,01) = pd(01,01) + rrt(120) * density(37) 
  pd(01,37) = pd(01,37) + rrt(120) * density(01) 
  pd(07,01) = pd(07,01) + rrt(120) * density(37) 
  pd(07,37) = pd(07,37) + rrt(120) * density(01) 
  pd(23,01) = pd(23,01) + rrt(120) * density(37) 
  pd(23,37) = pd(23,37) + rrt(120) * density(01) 
  pd(37,01) = pd(37,01) - rrt(120) * density(37) 
  pd(37,37) = pd(37,37) - rrt(120) * density(01) 
  pd(01,01) = pd(01,01) + rrt(121) * density(37) 
  pd(01,37) = pd(01,37) + rrt(121) * density(01) 
  pd(17,01) = pd(17,01) + rrt(121) * density(37) 
  pd(17,37) = pd(17,37) + rrt(121) * density(01) 
  pd(37,01) = pd(37,01) - rrt(121) * density(37) 
  pd(37,37) = pd(37,37) - rrt(121) * density(01) 
  pd(51,01) = pd(51,01) + rrt(121) * density(37) 
  pd(51,37) = pd(51,37) + rrt(121) * density(01) 
  pd(01,01) = pd(01,01) + rrt(122) * density(37) 
  pd(01,37) = pd(01,37) + rrt(122) * density(01) 
  pd(10,01) = pd(10,01) + rrt(122) * density(37) 
  pd(10,37) = pd(10,37) + rrt(122) * density(01) 
  pd(21,01) = pd(21,01) + rrt(122) * density(37) 
  pd(21,37) = pd(21,37) + rrt(122) * density(01) 
  pd(37,01) = pd(37,01) - rrt(122) * density(37) 
  pd(37,37) = pd(37,37) - rrt(122) * density(01) 
  pd(01,01) = pd(01,01) + rrt(123) * density(37) 
  pd(01,37) = pd(01,37) + rrt(123) * density(01) 
  pd(16,01) = pd(16,01) + rrt(123) * density(37) 
  pd(16,37) = pd(16,37) + rrt(123) * density(01) 
  pd(29,01) = pd(29,01) + rrt(123) * density(37) 
  pd(29,37) = pd(29,37) + rrt(123) * density(01) 
  pd(37,01) = pd(37,01) - rrt(123) * density(37) 
  pd(37,37) = pd(37,37) - rrt(123) * density(01) 
  pd(01,01) = pd(01,01) + rrt(124) * density(37) 
  pd(01,37) = pd(01,37) + rrt(124) * density(01) 
  pd(37,01) = pd(37,01) - rrt(124) * density(37) 
  pd(37,37) = pd(37,37) - rrt(124) * density(01) 
  pd(42,01) = pd(42,01) + rrt(124) * density(37) 
  pd(42,37) = pd(42,37) + rrt(124) * density(01) 
  pd(50,01) = pd(50,01) + rrt(124) * density(37) 
  pd(50,37) = pd(50,37) + rrt(124) * density(01) 
  pd(01,01) = pd(01,01) + rrt(125) * density(37) 
  pd(01,37) = pd(01,37) + rrt(125) * density(01) 
  pd(34,01) = pd(34,01) + rrt(125) * density(37) 
  pd(34,37) = pd(34,37) + rrt(125) * density(01) 
  pd(37,01) = pd(37,01) - rrt(125) * density(37) 
  pd(37,37) = pd(37,37) - rrt(125) * density(01) 
  pd(44,01) = pd(44,01) + rrt(125) * density(37) 
  pd(44,37) = pd(44,37) + rrt(125) * density(01) 
  pd(01,01) = pd(01,01) + rrt(126) * density(20) 
  pd(01,20) = pd(01,20) + rrt(126) * density(01) 
  pd(20,01) = pd(20,01) - rrt(126) * density(20) 
  pd(20,20) = pd(20,20) - rrt(126) * density(01) 
  pd(26,01) = pd(26,01) + rrt(126) * density(20) 
  pd(26,20) = pd(26,20) + rrt(126) * density(01) 
  pd(01,01) = pd(01,01) + rrt(127) * density(20) 
  pd(01,20) = pd(01,20) + rrt(127) * density(01) 
  pd(13,01) = pd(13,01) + rrt(127) * density(20) 
  pd(13,20) = pd(13,20) + rrt(127) * density(01) 
  pd(19,01) = pd(19,01) + rrt(127) * density(20) 
  pd(19,20) = pd(19,20) + rrt(127) * density(01) 
  pd(20,01) = pd(20,01) - rrt(127) * density(20) 
  pd(20,20) = pd(20,20) - rrt(127) * density(01) 
  pd(01,01) = pd(01,01) + rrt(128) * density(20) 
  pd(01,20) = pd(01,20) + rrt(128) * density(01) 
  pd(08,01) = pd(08,01) + rrt(128) * density(20) 
  pd(08,20) = pd(08,20) + rrt(128) * density(01) 
  pd(20,01) = pd(20,01) - rrt(128) * density(20) 
  pd(20,20) = pd(20,20) - rrt(128) * density(01) 
  pd(46,01) = pd(46,01) + rrt(128) * density(20) 
  pd(46,20) = pd(46,20) + rrt(128) * density(01) 
  pd(01,01) = pd(01,01) + rrt(129) * density(20) 
  pd(01,20) = pd(01,20) + rrt(129) * density(01) 
  pd(07,01) = pd(07,01) + rrt(129) * density(20) 
  pd(07,20) = pd(07,20) + rrt(129) * density(01) 
  pd(20,01) = pd(20,01) - rrt(129) * density(20) 
  pd(20,20) = pd(20,20) - rrt(129) * density(01) 
  pd(28,01) = pd(28,01) + rrt(129) * density(20) 
  pd(28,20) = pd(28,20) + rrt(129) * density(01) 
  pd(01,01) = pd(01,01) + rrt(130) * density(20) 
  pd(01,20) = pd(01,20) + rrt(130) * density(01) 
  pd(17,01) = pd(17,01) + rrt(130) * density(20) 
  pd(17,20) = pd(17,20) + rrt(130) * density(01) 
  pd(20,01) = pd(20,01) - rrt(130) * density(20) 
  pd(20,20) = pd(20,20) - rrt(130) * density(01) 
  pd(23,01) = pd(23,01) + rrt(130) * density(20) 
  pd(23,20) = pd(23,20) + rrt(130) * density(01) 
  pd(01,01) = pd(01,01) + rrt(131) * density(20) 
  pd(01,20) = pd(01,20) + rrt(131) * density(01) 
  pd(16,01) = pd(16,01) + rrt(131) * density(20) 
  pd(16,20) = pd(16,20) + rrt(131) * density(01) 
  pd(20,01) = pd(20,01) - rrt(131) * density(20) 
  pd(20,20) = pd(20,20) - rrt(131) * density(01) 
  pd(21,01) = pd(21,01) + rrt(131) * density(20) 
  pd(21,20) = pd(21,20) + rrt(131) * density(01) 
  pd(01,01) = pd(01,01) + rrt(132) * density(20) 
  pd(01,20) = pd(01,20) + rrt(132) * density(01) 
  pd(20,01) = pd(20,01) - rrt(132) * density(20) 
  pd(20,20) = pd(20,20) - rrt(132) * density(01) 
  pd(29,01) = pd(29,01) + rrt(132) * density(20) 
  pd(29,20) = pd(29,20) + rrt(132) * density(01) 
  pd(50,01) = pd(50,01) + rrt(132) * density(20) 
  pd(50,20) = pd(50,20) + rrt(132) * density(01) 
  pd(01,01) = pd(01,01) + rrt(133) * density(20) 
  pd(01,20) = pd(01,20) + rrt(133) * density(01) 
  pd(20,01) = pd(20,01) - rrt(133) * density(20) 
  pd(20,20) = pd(20,20) - rrt(133) * density(01) 
  pd(42,01) = pd(42,01) + rrt(133) * density(20) 
  pd(42,20) = pd(42,20) + rrt(133) * density(01) 
  pd(44,01) = pd(44,01) + rrt(133) * density(20) 
  pd(44,20) = pd(44,20) + rrt(133) * density(01) 
  pd(01,01) = pd(01,01) + rrt(134) * density(35) 
  pd(01,35) = pd(01,35) + rrt(134) * density(01) 
  pd(19,01) = pd(19,01) + rrt(134) * density(35) 
  pd(19,35) = pd(19,35) + rrt(134) * density(01) 
  pd(35,01) = pd(35,01) - rrt(134) * density(35) 
  pd(35,35) = pd(35,35) - rrt(134) * density(01) 
  pd(01,01) = pd(01,01) + rrt(135) * density(35) 
  pd(01,35) = pd(01,35) + rrt(135) * density(01) 
  pd(07,01) = pd(07,01) + rrt(135) * density(35) 
  pd(07,35) = pd(07,35) + rrt(135) * density(01) 
  pd(08,01) = pd(08,01) + rrt(135) * density(35) 
  pd(08,35) = pd(08,35) + rrt(135) * density(01) 
  pd(35,01) = pd(35,01) - rrt(135) * density(35) 
  pd(35,35) = pd(35,35) - rrt(135) * density(01) 
  pd(01,01) = pd(01,01) + rrt(136) * density(35) 
  pd(01,35) = pd(01,35) + rrt(136) * density(01) 
  pd(17,01) = pd(17,01) + rrt(136) * density(35) 
  pd(17,35) = pd(17,35) + rrt(136) * density(01) 
  pd(28,01) = pd(28,01) + rrt(136) * density(35) 
  pd(28,35) = pd(28,35) + rrt(136) * density(01) 
  pd(35,01) = pd(35,01) - rrt(136) * density(35) 
  pd(35,35) = pd(35,35) - rrt(136) * density(01) 
  pd(01,01) = pd(01,01) + rrt(137) * density(35) 
  pd(01,35) = pd(01,35) + rrt(137) * density(01) 
  pd(21,01) = pd(21,01) + rrt(137) * density(35) 
  pd(21,35) = pd(21,35) + rrt(137) * density(01) 
  pd(35,01) = pd(35,01) - rrt(137) * density(35) 
  pd(35,35) = pd(35,35) - rrt(137) * density(01) 
  pd(50,01) = pd(50,01) + rrt(137) * density(35) 
  pd(50,35) = pd(50,35) + rrt(137) * density(01) 
  pd(01,01) = pd(01,01) + rrt(138) * density(35) 
  pd(01,35) = pd(01,35) + rrt(138) * density(01) 
  pd(29,01) = pd(29,01) + rrt(138) * density(35) 
  pd(29,35) = pd(29,35) + rrt(138) * density(01) 
  pd(35,01) = pd(35,01) - rrt(138) * density(35) 
  pd(35,35) = pd(35,35) - rrt(138) * density(01) 
  pd(44,01) = pd(44,01) + rrt(138) * density(35) 
  pd(44,35) = pd(44,35) + rrt(138) * density(01) 
  pd(13,01) = pd(13,01) + rrt(139) * density(22) 
  pd(13,22) = pd(13,22) + rrt(139) * density(01) 
  pd(22,01) = pd(22,01) - rrt(139) * density(22) 
  pd(22,22) = pd(22,22) - rrt(139) * density(01) 
  pd(23,01) = pd(23,01) + rrt(139) * density(22) 
  pd(23,22) = pd(23,22) + rrt(139) * density(01) 
  pd(13,01) = pd(13,01) + rrt(140) * density(25) 
  pd(13,25) = pd(13,25) + rrt(140) * density(01) 
  pd(23,01) = pd(23,01) + rrt(140) * density(25) 
  pd(23,25) = pd(23,25) + rrt(140) * density(01) 
  pd(25,01) = pd(25,01) - rrt(140) * density(25) 
  pd(25,25) = pd(25,25) - rrt(140) * density(01) 
  pd(11,01) = pd(11,01) + rrt(141) * density(22) 
  pd(11,22) = pd(11,22) + rrt(141) * density(01) 
  pd(22,01) = pd(22,01) - rrt(141) * density(22) 
  pd(22,22) = pd(22,22) - rrt(141) * density(01) 
  pd(28,01) = pd(28,01) + rrt(141) * density(22) 
  pd(28,22) = pd(28,22) + rrt(141) * density(01) 
  pd(11,01) = pd(11,01) + rrt(142) * density(25) 
  pd(11,25) = pd(11,25) + rrt(142) * density(01) 
  pd(25,01) = pd(25,01) - rrt(142) * density(25) 
  pd(25,25) = pd(25,25) - rrt(142) * density(01) 
  pd(28,01) = pd(28,01) + rrt(142) * density(25) 
  pd(28,25) = pd(28,25) + rrt(142) * density(01) 
  pd(08,01) = pd(08,01) + rrt(143) * density(22) 
  pd(08,22) = pd(08,22) + rrt(143) * density(01) 
  pd(11,01) = pd(11,01) + rrt(143) * density(22) 
  pd(11,22) = pd(11,22) + rrt(143) * density(01) 
  pd(13,01) = pd(13,01) + rrt(143) * density(22) 
  pd(13,22) = pd(13,22) + rrt(143) * density(01) 
  pd(22,01) = pd(22,01) - rrt(143) * density(22) 
  pd(22,22) = pd(22,22) - rrt(143) * density(01) 
  pd(08,01) = pd(08,01) + rrt(144) * density(25) 
  pd(08,25) = pd(08,25) + rrt(144) * density(01) 
  pd(11,01) = pd(11,01) + rrt(144) * density(25) 
  pd(11,25) = pd(11,25) + rrt(144) * density(01) 
  pd(13,01) = pd(13,01) + rrt(144) * density(25) 
  pd(13,25) = pd(13,25) + rrt(144) * density(01) 
  pd(25,01) = pd(25,01) - rrt(144) * density(25) 
  pd(25,25) = pd(25,25) - rrt(144) * density(01) 
  pd(01,01) = pd(01,01) + rrt(145) * density(22) 
  pd(01,22) = pd(01,22) + rrt(145) * density(01) 
  pd(10,01) = pd(10,01) + rrt(145) * density(22) 
  pd(10,22) = pd(10,22) + rrt(145) * density(01) 
  pd(22,01) = pd(22,01) - rrt(145) * density(22) 
  pd(22,22) = pd(22,22) - rrt(145) * density(01) 
  pd(01,01) = pd(01,01) + rrt(146) * density(25) 
  pd(01,25) = pd(01,25) + rrt(146) * density(01) 
  pd(10,01) = pd(10,01) + rrt(146) * density(25) 
  pd(10,25) = pd(10,25) + rrt(146) * density(01) 
  pd(25,01) = pd(25,01) - rrt(146) * density(25) 
  pd(25,25) = pd(25,25) - rrt(146) * density(01) 
  pd(01,01) = pd(01,01) + rrt(147) * density(22) 
  pd(01,22) = pd(01,22) + rrt(147) * density(01) 
  pd(13,01) = pd(13,01) + rrt(147) * density(22) 
  pd(13,22) = pd(13,22) + rrt(147) * density(01) 
  pd(16,01) = pd(16,01) + rrt(147) * density(22) 
  pd(16,22) = pd(16,22) + rrt(147) * density(01) 
  pd(22,01) = pd(22,01) - rrt(147) * density(22) 
  pd(22,22) = pd(22,22) - rrt(147) * density(01) 
  pd(01,01) = pd(01,01) + rrt(148) * density(25) 
  pd(01,25) = pd(01,25) + rrt(148) * density(01) 
  pd(13,01) = pd(13,01) + rrt(148) * density(25) 
  pd(13,25) = pd(13,25) + rrt(148) * density(01) 
  pd(16,01) = pd(16,01) + rrt(148) * density(25) 
  pd(16,25) = pd(16,25) + rrt(148) * density(01) 
  pd(25,01) = pd(25,01) - rrt(148) * density(25) 
  pd(25,25) = pd(25,25) - rrt(148) * density(01) 
  pd(01,01) = pd(01,01) + rrt(149) * density(22) 
  pd(01,22) = pd(01,22) + rrt(149) * density(01) 
  pd(11,01) = pd(11,01) + rrt(149) * density(22) 
  pd(11,22) = pd(11,22) + rrt(149) * density(01) 
  pd(22,01) = pd(22,01) - rrt(149) * density(22) 
  pd(22,22) = pd(22,22) - rrt(149) * density(01) 
  pd(50,01) = pd(50,01) + rrt(149) * density(22) 
  pd(50,22) = pd(50,22) + rrt(149) * density(01) 
  pd(01,01) = pd(01,01) + rrt(150) * density(25) 
  pd(01,25) = pd(01,25) + rrt(150) * density(01) 
  pd(11,01) = pd(11,01) + rrt(150) * density(25) 
  pd(11,25) = pd(11,25) + rrt(150) * density(01) 
  pd(25,01) = pd(25,01) - rrt(150) * density(25) 
  pd(25,25) = pd(25,25) - rrt(150) * density(01) 
  pd(50,01) = pd(50,01) + rrt(150) * density(25) 
  pd(50,25) = pd(50,25) + rrt(150) * density(01) 
  pd(01,01) = pd(01,01) + rrt(151) * density(22) 
  pd(01,22) = pd(01,22) + rrt(151) * density(01) 
  pd(11,01) = pd(11,01) + rrt(151) * density(22) 
  pd(11,22) = pd(11,22) + rrt(151) * density(01) 
  pd(13,01) = pd(13,01) + rrt(151) * density(22) 
  pd(13,22) = pd(13,22) + rrt(151) * density(01) 
  pd(22,01) = pd(22,01) - rrt(151) * density(22) 
  pd(22,22) = pd(22,22) - rrt(151) * density(01) 
  pd(44,01) = pd(44,01) + rrt(151) * density(22) 
  pd(44,22) = pd(44,22) + rrt(151) * density(01) 
  pd(01,01) = pd(01,01) + rrt(152) * density(25) 
  pd(01,25) = pd(01,25) + rrt(152) * density(01) 
  pd(11,01) = pd(11,01) + rrt(152) * density(25) 
  pd(11,25) = pd(11,25) + rrt(152) * density(01) 
  pd(13,01) = pd(13,01) + rrt(152) * density(25) 
  pd(13,25) = pd(13,25) + rrt(152) * density(01) 
  pd(25,01) = pd(25,01) - rrt(152) * density(25) 
  pd(25,25) = pd(25,25) - rrt(152) * density(01) 
  pd(44,01) = pd(44,01) + rrt(152) * density(25) 
  pd(44,25) = pd(44,25) + rrt(152) * density(01) 
  pd(13,01) = pd(13,01) + rrt(153) * density(39) 
  pd(13,39) = pd(13,39) + rrt(153) * density(01) 
  pd(34,01) = pd(34,01) + rrt(153) * density(39) 
  pd(34,39) = pd(34,39) + rrt(153) * density(01) 
  pd(39,01) = pd(39,01) - rrt(153) * density(39) 
  pd(39,39) = pd(39,39) - rrt(153) * density(01) 
  pd(13,01) = pd(13,01) + rrt(154) * density(41) 
  pd(13,41) = pd(13,41) + rrt(154) * density(01) 
  pd(34,01) = pd(34,01) + rrt(154) * density(41) 
  pd(34,41) = pd(34,41) + rrt(154) * density(01) 
  pd(41,01) = pd(41,01) - rrt(154) * density(41) 
  pd(41,41) = pd(41,41) - rrt(154) * density(01) 
  pd(11,01) = pd(11,01) + rrt(155) * density(39) 
  pd(11,39) = pd(11,39) + rrt(155) * density(01) 
  pd(39,01) = pd(39,01) - rrt(155) * density(39) 
  pd(39,39) = pd(39,39) - rrt(155) * density(01) 
  pd(42,01) = pd(42,01) + rrt(155) * density(39) 
  pd(42,39) = pd(42,39) + rrt(155) * density(01) 
  pd(11,01) = pd(11,01) + rrt(156) * density(41) 
  pd(11,41) = pd(11,41) + rrt(156) * density(01) 
  pd(41,01) = pd(41,01) - rrt(156) * density(41) 
  pd(41,41) = pd(41,41) - rrt(156) * density(01) 
  pd(42,01) = pd(42,01) + rrt(156) * density(41) 
  pd(42,41) = pd(42,41) + rrt(156) * density(01) 
  pd(11,01) = pd(11,01) + rrt(157) * density(39) 
  pd(11,39) = pd(11,39) + rrt(157) * density(01) 
  pd(13,01) = pd(13,01) + rrt(157) * density(39) 
  pd(13,39) = pd(13,39) + rrt(157) * density(01) 
  pd(29,01) = pd(29,01) + rrt(157) * density(39) 
  pd(29,39) = pd(29,39) + rrt(157) * density(01) 
  pd(39,01) = pd(39,01) - rrt(157) * density(39) 
  pd(39,39) = pd(39,39) - rrt(157) * density(01) 
  pd(11,01) = pd(11,01) + rrt(158) * density(41) 
  pd(11,41) = pd(11,41) + rrt(158) * density(01) 
  pd(13,01) = pd(13,01) + rrt(158) * density(41) 
  pd(13,41) = pd(13,41) + rrt(158) * density(01) 
  pd(29,01) = pd(29,01) + rrt(158) * density(41) 
  pd(29,41) = pd(29,41) + rrt(158) * density(01) 
  pd(41,01) = pd(41,01) - rrt(158) * density(41) 
  pd(41,41) = pd(41,41) - rrt(158) * density(01) 
  pd(11,01) = pd(11,01) + rrt(159) * density(39) * 2.0d0
  pd(11,39) = pd(11,39) + rrt(159) * density(01) * 2.0d0
  pd(21,01) = pd(21,01) + rrt(159) * density(39) 
  pd(21,39) = pd(21,39) + rrt(159) * density(01) 
  pd(39,01) = pd(39,01) - rrt(159) * density(39) 
  pd(39,39) = pd(39,39) - rrt(159) * density(01) 
  pd(11,01) = pd(11,01) + rrt(160) * density(41) * 2.0d0
  pd(11,41) = pd(11,41) + rrt(160) * density(01) * 2.0d0
  pd(21,01) = pd(21,01) + rrt(160) * density(41) 
  pd(21,41) = pd(21,41) + rrt(160) * density(01) 
  pd(41,01) = pd(41,01) - rrt(160) * density(41) 
  pd(41,41) = pd(41,41) - rrt(160) * density(01) 
  pd(28,01) = pd(28,01) + rrt(161) * density(39) 
  pd(28,39) = pd(28,39) + rrt(161) * density(01) 
  pd(39,01) = pd(39,01) - rrt(161) * density(39) 
  pd(39,39) = pd(39,39) - rrt(161) * density(01) 
  pd(51,01) = pd(51,01) + rrt(161) * density(39) 
  pd(51,39) = pd(51,39) + rrt(161) * density(01) 
  pd(28,01) = pd(28,01) + rrt(162) * density(41) 
  pd(28,41) = pd(28,41) + rrt(162) * density(01) 
  pd(41,01) = pd(41,01) - rrt(162) * density(41) 
  pd(41,41) = pd(41,41) - rrt(162) * density(01) 
  pd(51,01) = pd(51,01) + rrt(162) * density(41) 
  pd(51,41) = pd(51,41) + rrt(162) * density(01) 
  pd(23,01) = pd(23,01) + rrt(163) * density(39) * 2.0d0
  pd(23,39) = pd(23,39) + rrt(163) * density(01) * 2.0d0
  pd(39,01) = pd(39,01) - rrt(163) * density(39) 
  pd(39,39) = pd(39,39) - rrt(163) * density(01) 
  pd(23,01) = pd(23,01) + rrt(164) * density(41) * 2.0d0
  pd(23,41) = pd(23,41) + rrt(164) * density(01) * 2.0d0
  pd(41,01) = pd(41,01) - rrt(164) * density(41) 
  pd(41,41) = pd(41,41) - rrt(164) * density(01) 
  pd(13,01) = pd(13,01) + rrt(165) * density(33) 
  pd(13,33) = pd(13,33) + rrt(165) * density(01) 
  pd(29,01) = pd(29,01) + rrt(165) * density(33) 
  pd(29,33) = pd(29,33) + rrt(165) * density(01) 
  pd(33,01) = pd(33,01) - rrt(165) * density(33) 
  pd(33,33) = pd(33,33) - rrt(165) * density(01) 
  pd(02,01) = pd(02,01) - rrt(166) * density(02) 
  pd(02,02) = pd(02,02) - rrt(166) * density(01) 
  pd(13,01) = pd(13,01) + rrt(166) * density(02) 
  pd(13,02) = pd(13,02) + rrt(166) * density(01) 
  pd(29,01) = pd(29,01) + rrt(166) * density(02) 
  pd(29,02) = pd(29,02) + rrt(166) * density(01) 
  pd(11,01) = pd(11,01) + rrt(167) * density(33) 
  pd(11,33) = pd(11,33) + rrt(167) * density(01) 
  pd(21,01) = pd(21,01) + rrt(167) * density(33) 
  pd(21,33) = pd(21,33) + rrt(167) * density(01) 
  pd(33,01) = pd(33,01) - rrt(167) * density(33) 
  pd(33,33) = pd(33,33) - rrt(167) * density(01) 
  pd(02,01) = pd(02,01) - rrt(168) * density(02) 
  pd(02,02) = pd(02,02) - rrt(168) * density(01) 
  pd(11,01) = pd(11,01) + rrt(168) * density(02) 
  pd(11,02) = pd(11,02) + rrt(168) * density(01) 
  pd(21,01) = pd(21,01) + rrt(168) * density(02) 
  pd(21,02) = pd(21,02) + rrt(168) * density(01) 
  pd(13,01) = pd(13,01) + rrt(169) * density(33) * 2.0d0
  pd(13,33) = pd(13,33) + rrt(169) * density(01) * 2.0d0
  pd(21,01) = pd(21,01) + rrt(169) * density(33) 
  pd(21,33) = pd(21,33) + rrt(169) * density(01) 
  pd(33,01) = pd(33,01) - rrt(169) * density(33) 
  pd(33,33) = pd(33,33) - rrt(169) * density(01) 
  pd(02,01) = pd(02,01) - rrt(170) * density(02) 
  pd(02,02) = pd(02,02) - rrt(170) * density(01) 
  pd(13,01) = pd(13,01) + rrt(170) * density(02) * 2.0d0
  pd(13,02) = pd(13,02) + rrt(170) * density(01) * 2.0d0
  pd(21,01) = pd(21,01) + rrt(170) * density(02) 
  pd(21,02) = pd(21,02) + rrt(170) * density(01) 
  pd(08,01) = pd(08,01) + rrt(171) * density(33) 
  pd(08,33) = pd(08,33) + rrt(171) * density(01) 
  pd(23,01) = pd(23,01) + rrt(171) * density(33) 
  pd(23,33) = pd(23,33) + rrt(171) * density(01) 
  pd(33,01) = pd(33,01) - rrt(171) * density(33) 
  pd(33,33) = pd(33,33) - rrt(171) * density(01) 
  pd(02,01) = pd(02,01) - rrt(172) * density(02) 
  pd(02,02) = pd(02,02) - rrt(172) * density(01) 
  pd(08,01) = pd(08,01) + rrt(172) * density(02) 
  pd(08,02) = pd(08,02) + rrt(172) * density(01) 
  pd(23,01) = pd(23,01) + rrt(172) * density(02) 
  pd(23,02) = pd(23,02) + rrt(172) * density(01) 
  pd(28,01) = pd(28,01) + rrt(173) * density(33) * 2.0d0
  pd(28,33) = pd(28,33) + rrt(173) * density(01) * 2.0d0
  pd(33,01) = pd(33,01) - rrt(173) * density(33) 
  pd(33,33) = pd(33,33) - rrt(173) * density(01) 
  pd(02,01) = pd(02,01) - rrt(174) * density(02) 
  pd(02,02) = pd(02,02) - rrt(174) * density(01) 
  pd(28,01) = pd(28,01) + rrt(174) * density(02) * 2.0d0
  pd(28,02) = pd(28,02) + rrt(174) * density(01) * 2.0d0
  pd(08,01) = pd(08,01) + rrt(175) * density(52) * 2.0d0
  pd(08,52) = pd(08,52) + rrt(175) * density(01) * 2.0d0
  pd(52,01) = pd(52,01) - rrt(175) * density(52) 
  pd(52,52) = pd(52,52) - rrt(175) * density(01) 
  pd(08,01) = pd(08,01) + rrt(176) * density(09) * 2.0d0
  pd(08,09) = pd(08,09) + rrt(176) * density(01) * 2.0d0
  pd(09,01) = pd(09,01) - rrt(176) * density(09) 
  pd(09,09) = pd(09,09) - rrt(176) * density(01) 
  pd(08,01) = pd(08,01) + rrt(177) * density(38) * 2.0d0
  pd(08,38) = pd(08,38) + rrt(177) * density(01) * 2.0d0
  pd(38,01) = pd(38,01) - rrt(177) * density(38) 
  pd(38,38) = pd(38,38) - rrt(177) * density(01) 
  pd(13,01) = pd(13,01) + rrt(178) * density(32) 
  pd(13,32) = pd(13,32) + rrt(178) * density(01) 
  pd(18,01) = pd(18,01) + rrt(178) * density(32) 
  pd(18,32) = pd(18,32) + rrt(178) * density(01) 
  pd(32,01) = pd(32,01) - rrt(178) * density(32) 
  pd(32,32) = pd(32,32) - rrt(178) * density(01) 
  pd(13,01) = pd(13,01) + rrt(179) * density(24) 
  pd(13,24) = pd(13,24) + rrt(179) * density(01) 
  pd(18,01) = pd(18,01) + rrt(179) * density(24) 
  pd(18,24) = pd(18,24) + rrt(179) * density(01) 
  pd(24,01) = pd(24,01) - rrt(179) * density(24) 
  pd(24,24) = pd(24,24) - rrt(179) * density(01) 
  pd(11,01) = pd(11,01) + rrt(180) * density(32) 
  pd(11,32) = pd(11,32) + rrt(180) * density(01) 
  pd(32,01) = pd(32,01) - rrt(180) * density(32) 
  pd(32,32) = pd(32,32) - rrt(180) * density(01) 
  pd(37,01) = pd(37,01) + rrt(180) * density(32) 
  pd(37,32) = pd(37,32) + rrt(180) * density(01) 
  pd(11,01) = pd(11,01) + rrt(181) * density(24) 
  pd(11,24) = pd(11,24) + rrt(181) * density(01) 
  pd(24,01) = pd(24,01) - rrt(181) * density(24) 
  pd(24,24) = pd(24,24) - rrt(181) * density(01) 
  pd(37,01) = pd(37,01) + rrt(181) * density(24) 
  pd(37,24) = pd(37,24) + rrt(181) * density(01) 
  pd(11,01) = pd(11,01) + rrt(182) * density(32) * 2.0d0
  pd(11,32) = pd(11,32) + rrt(182) * density(01) * 2.0d0
  pd(32,01) = pd(32,01) - rrt(182) * density(32) 
  pd(32,32) = pd(32,32) - rrt(182) * density(01) 
  pd(35,01) = pd(35,01) + rrt(182) * density(32) 
  pd(35,32) = pd(35,32) + rrt(182) * density(01) 
  pd(11,01) = pd(11,01) + rrt(183) * density(24) * 2.0d0
  pd(11,24) = pd(11,24) + rrt(183) * density(01) * 2.0d0
  pd(24,01) = pd(24,01) - rrt(183) * density(24) 
  pd(24,24) = pd(24,24) - rrt(183) * density(01) 
  pd(35,01) = pd(35,01) + rrt(183) * density(24) 
  pd(35,24) = pd(35,24) + rrt(183) * density(01) 
  pd(28,01) = pd(28,01) + rrt(184) * density(32) 
  pd(28,32) = pd(28,32) + rrt(184) * density(01) 
  pd(32,01) = pd(32,01) - rrt(184) * density(32) 
  pd(32,32) = pd(32,32) - rrt(184) * density(01) 
  pd(43,01) = pd(43,01) + rrt(184) * density(32) 
  pd(43,32) = pd(43,32) + rrt(184) * density(01) 
  pd(24,01) = pd(24,01) - rrt(185) * density(24) 
  pd(24,24) = pd(24,24) - rrt(185) * density(01) 
  pd(28,01) = pd(28,01) + rrt(185) * density(24) 
  pd(28,24) = pd(28,24) + rrt(185) * density(01) 
  pd(43,01) = pd(43,01) + rrt(185) * density(24) 
  pd(43,24) = pd(43,24) + rrt(185) * density(01) 
  pd(23,01) = pd(23,01) + rrt(186) * density(32) 
  pd(23,32) = pd(23,32) + rrt(186) * density(01) 
  pd(32,01) = pd(32,01) - rrt(186) * density(32) 
  pd(32,32) = pd(32,32) - rrt(186) * density(01) 
  pd(34,01) = pd(34,01) + rrt(186) * density(32) 
  pd(34,32) = pd(34,32) + rrt(186) * density(01) 
  pd(23,01) = pd(23,01) + rrt(187) * density(24) 
  pd(23,24) = pd(23,24) + rrt(187) * density(01) 
  pd(24,01) = pd(24,01) - rrt(187) * density(24) 
  pd(24,24) = pd(24,24) - rrt(187) * density(01) 
  pd(34,01) = pd(34,01) + rrt(187) * density(24) 
  pd(34,24) = pd(34,24) + rrt(187) * density(01) 
  pd(32,01) = pd(32,01) - rrt(188) * density(32) 
  pd(32,32) = pd(32,32) - rrt(188) * density(01) 
  pd(42,01) = pd(42,01) + rrt(188) * density(32) 
  pd(42,32) = pd(42,32) + rrt(188) * density(01) 
  pd(51,01) = pd(51,01) + rrt(188) * density(32) 
  pd(51,32) = pd(51,32) + rrt(188) * density(01) 
  pd(24,01) = pd(24,01) - rrt(189) * density(24) 
  pd(24,24) = pd(24,24) - rrt(189) * density(01) 
  pd(42,01) = pd(42,01) + rrt(189) * density(24) 
  pd(42,24) = pd(42,24) + rrt(189) * density(01) 
  pd(51,01) = pd(51,01) + rrt(189) * density(24) 
  pd(51,24) = pd(51,24) + rrt(189) * density(01) 
  pd(01,01) = pd(01,01) + rrt(190) * density(39) 
  pd(01,39) = pd(01,39) + rrt(190) * density(01) 
  pd(39,01) = pd(39,01) - rrt(190) * density(39) 
  pd(39,39) = pd(39,39) - rrt(190) * density(01) 
  pd(40,01) = pd(40,01) + rrt(190) * density(39) 
  pd(40,39) = pd(40,39) + rrt(190) * density(01) 
  pd(01,01) = pd(01,01) + rrt(191) * density(41) 
  pd(01,41) = pd(01,41) + rrt(191) * density(01) 
  pd(40,01) = pd(40,01) + rrt(191) * density(41) 
  pd(40,41) = pd(40,41) + rrt(191) * density(01) 
  pd(41,01) = pd(41,01) - rrt(191) * density(41) 
  pd(41,41) = pd(41,41) - rrt(191) * density(01) 
  pd(01,01) = pd(01,01) + rrt(192) * density(39) 
  pd(01,39) = pd(01,39) + rrt(192) * density(01) 
  pd(04,01) = pd(04,01) + rrt(192) * density(39) 
  pd(04,39) = pd(04,39) + rrt(192) * density(01) 
  pd(13,01) = pd(13,01) + rrt(192) * density(39) 
  pd(13,39) = pd(13,39) + rrt(192) * density(01) 
  pd(39,01) = pd(39,01) - rrt(192) * density(39) 
  pd(39,39) = pd(39,39) - rrt(192) * density(01) 
  pd(01,01) = pd(01,01) + rrt(193) * density(41) 
  pd(01,41) = pd(01,41) + rrt(193) * density(01) 
  pd(04,01) = pd(04,01) + rrt(193) * density(41) 
  pd(04,41) = pd(04,41) + rrt(193) * density(01) 
  pd(13,01) = pd(13,01) + rrt(193) * density(41) 
  pd(13,41) = pd(13,41) + rrt(193) * density(01) 
  pd(41,01) = pd(41,01) - rrt(193) * density(41) 
  pd(41,41) = pd(41,41) - rrt(193) * density(01) 
  pd(01,01) = pd(01,01) + rrt(194) * density(39) 
  pd(01,39) = pd(01,39) + rrt(194) * density(01) 
  pd(11,01) = pd(11,01) + rrt(194) * density(39) 
  pd(11,39) = pd(11,39) + rrt(194) * density(01) 
  pd(39,01) = pd(39,01) - rrt(194) * density(39) 
  pd(39,39) = pd(39,39) - rrt(194) * density(01) 
  pd(46,01) = pd(46,01) + rrt(194) * density(39) 
  pd(46,39) = pd(46,39) + rrt(194) * density(01) 
  pd(01,01) = pd(01,01) + rrt(195) * density(41) 
  pd(01,41) = pd(01,41) + rrt(195) * density(01) 
  pd(11,01) = pd(11,01) + rrt(195) * density(41) 
  pd(11,41) = pd(11,41) + rrt(195) * density(01) 
  pd(41,01) = pd(41,01) - rrt(195) * density(41) 
  pd(41,41) = pd(41,41) - rrt(195) * density(01) 
  pd(46,01) = pd(46,01) + rrt(195) * density(41) 
  pd(46,41) = pd(46,41) + rrt(195) * density(01) 
  pd(01,01) = pd(01,01) + rrt(196) * density(39) 
  pd(01,39) = pd(01,39) + rrt(196) * density(01) 
  pd(07,01) = pd(07,01) + rrt(196) * density(39) 
  pd(07,39) = pd(07,39) + rrt(196) * density(01) 
  pd(11,01) = pd(11,01) + rrt(196) * density(39) 
  pd(11,39) = pd(11,39) + rrt(196) * density(01) 
  pd(13,01) = pd(13,01) + rrt(196) * density(39) 
  pd(13,39) = pd(13,39) + rrt(196) * density(01) 
  pd(39,01) = pd(39,01) - rrt(196) * density(39) 
  pd(39,39) = pd(39,39) - rrt(196) * density(01) 
  pd(01,01) = pd(01,01) + rrt(197) * density(41) 
  pd(01,41) = pd(01,41) + rrt(197) * density(01) 
  pd(07,01) = pd(07,01) + rrt(197) * density(41) 
  pd(07,41) = pd(07,41) + rrt(197) * density(01) 
  pd(11,01) = pd(11,01) + rrt(197) * density(41) 
  pd(11,41) = pd(11,41) + rrt(197) * density(01) 
  pd(13,01) = pd(13,01) + rrt(197) * density(41) 
  pd(13,41) = pd(13,41) + rrt(197) * density(01) 
  pd(41,01) = pd(41,01) - rrt(197) * density(41) 
  pd(41,41) = pd(41,41) - rrt(197) * density(01) 
  pd(01,01) = pd(01,01) + rrt(198) * density(39) 
  pd(01,39) = pd(01,39) + rrt(198) * density(01) 
  pd(11,01) = pd(11,01) + rrt(198) * density(39) * 2.0d0
  pd(11,39) = pd(11,39) + rrt(198) * density(01) * 2.0d0
  pd(17,01) = pd(17,01) + rrt(198) * density(39) 
  pd(17,39) = pd(17,39) + rrt(198) * density(01) 
  pd(39,01) = pd(39,01) - rrt(198) * density(39) 
  pd(39,39) = pd(39,39) - rrt(198) * density(01) 
  pd(01,01) = pd(01,01) + rrt(199) * density(41) 
  pd(01,41) = pd(01,41) + rrt(199) * density(01) 
  pd(11,01) = pd(11,01) + rrt(199) * density(41) * 2.0d0
  pd(11,41) = pd(11,41) + rrt(199) * density(01) * 2.0d0
  pd(17,01) = pd(17,01) + rrt(199) * density(41) 
  pd(17,41) = pd(17,41) + rrt(199) * density(01) 
  pd(41,01) = pd(41,01) - rrt(199) * density(41) 
  pd(41,41) = pd(41,41) - rrt(199) * density(01) 
  pd(01,01) = pd(01,01) + rrt(200) * density(39) 
  pd(01,39) = pd(01,39) + rrt(200) * density(01) 
  pd(16,01) = pd(16,01) + rrt(200) * density(39) 
  pd(16,39) = pd(16,39) + rrt(200) * density(01) 
  pd(23,01) = pd(23,01) + rrt(200) * density(39) 
  pd(23,39) = pd(23,39) + rrt(200) * density(01) 
  pd(39,01) = pd(39,01) - rrt(200) * density(39) 
  pd(39,39) = pd(39,39) - rrt(200) * density(01) 
  pd(01,01) = pd(01,01) + rrt(201) * density(41) 
  pd(01,41) = pd(01,41) + rrt(201) * density(01) 
  pd(16,01) = pd(16,01) + rrt(201) * density(41) 
  pd(16,41) = pd(16,41) + rrt(201) * density(01) 
  pd(23,01) = pd(23,01) + rrt(201) * density(41) 
  pd(23,41) = pd(23,41) + rrt(201) * density(01) 
  pd(41,01) = pd(41,01) - rrt(201) * density(41) 
  pd(41,41) = pd(41,41) - rrt(201) * density(01) 
  pd(01,01) = pd(01,01) + rrt(202) * density(39) 
  pd(01,39) = pd(01,39) + rrt(202) * density(01) 
  pd(39,01) = pd(39,01) - rrt(202) * density(39) 
  pd(39,39) = pd(39,39) - rrt(202) * density(01) 
  pd(50,01) = pd(50,01) + rrt(202) * density(39) 
  pd(50,39) = pd(50,39) + rrt(202) * density(01) 
  pd(51,01) = pd(51,01) + rrt(202) * density(39) 
  pd(51,39) = pd(51,39) + rrt(202) * density(01) 
  pd(01,01) = pd(01,01) + rrt(203) * density(41) 
  pd(01,41) = pd(01,41) + rrt(203) * density(01) 
  pd(41,01) = pd(41,01) - rrt(203) * density(41) 
  pd(41,41) = pd(41,41) - rrt(203) * density(01) 
  pd(50,01) = pd(50,01) + rrt(203) * density(41) 
  pd(50,41) = pd(50,41) + rrt(203) * density(01) 
  pd(51,01) = pd(51,01) + rrt(203) * density(41) 
  pd(51,41) = pd(51,41) + rrt(203) * density(01) 
  pd(01,01) = pd(01,01) + rrt(204) * density(33) 
  pd(01,33) = pd(01,33) + rrt(204) * density(01) 
  pd(33,01) = pd(33,01) - rrt(204) * density(33) 
  pd(33,33) = pd(33,33) - rrt(204) * density(01) 
  pd(46,01) = pd(46,01) + rrt(204) * density(33) 
  pd(46,33) = pd(46,33) + rrt(204) * density(01) 
  pd(01,01) = pd(01,01) + rrt(205) * density(02) 
  pd(01,02) = pd(01,02) + rrt(205) * density(01) 
  pd(02,01) = pd(02,01) - rrt(205) * density(02) 
  pd(02,02) = pd(02,02) - rrt(205) * density(01) 
  pd(46,01) = pd(46,01) + rrt(205) * density(02) 
  pd(46,02) = pd(46,02) + rrt(205) * density(01) 
  pd(01,01) = pd(01,01) + rrt(206) * density(33) 
  pd(01,33) = pd(01,33) + rrt(206) * density(01) 
  pd(07,01) = pd(07,01) + rrt(206) * density(33) 
  pd(07,33) = pd(07,33) + rrt(206) * density(01) 
  pd(13,01) = pd(13,01) + rrt(206) * density(33) 
  pd(13,33) = pd(13,33) + rrt(206) * density(01) 
  pd(33,01) = pd(33,01) - rrt(206) * density(33) 
  pd(33,33) = pd(33,33) - rrt(206) * density(01) 
  pd(01,01) = pd(01,01) + rrt(207) * density(02) 
  pd(01,02) = pd(01,02) + rrt(207) * density(01) 
  pd(02,01) = pd(02,01) - rrt(207) * density(02) 
  pd(02,02) = pd(02,02) - rrt(207) * density(01) 
  pd(07,01) = pd(07,01) + rrt(207) * density(02) 
  pd(07,02) = pd(07,02) + rrt(207) * density(01) 
  pd(13,01) = pd(13,01) + rrt(207) * density(02) 
  pd(13,02) = pd(13,02) + rrt(207) * density(01) 
  pd(01,01) = pd(01,01) + rrt(208) * density(33) 
  pd(01,33) = pd(01,33) + rrt(208) * density(01) 
  pd(08,01) = pd(08,01) + rrt(208) * density(33) 
  pd(08,33) = pd(08,33) + rrt(208) * density(01) 
  pd(16,01) = pd(16,01) + rrt(208) * density(33) 
  pd(16,33) = pd(16,33) + rrt(208) * density(01) 
  pd(33,01) = pd(33,01) - rrt(208) * density(33) 
  pd(33,33) = pd(33,33) - rrt(208) * density(01) 
  pd(01,01) = pd(01,01) + rrt(209) * density(02) 
  pd(01,02) = pd(01,02) + rrt(209) * density(01) 
  pd(02,01) = pd(02,01) - rrt(209) * density(02) 
  pd(02,02) = pd(02,02) - rrt(209) * density(01) 
  pd(08,01) = pd(08,01) + rrt(209) * density(02) 
  pd(08,02) = pd(08,02) + rrt(209) * density(01) 
  pd(16,01) = pd(16,01) + rrt(209) * density(02) 
  pd(16,02) = pd(16,02) + rrt(209) * density(01) 
  pd(01,01) = pd(01,01) + rrt(210) * density(33) 
  pd(01,33) = pd(01,33) + rrt(210) * density(01) 
  pd(28,01) = pd(28,01) + rrt(210) * density(33) 
  pd(28,33) = pd(28,33) + rrt(210) * density(01) 
  pd(33,01) = pd(33,01) - rrt(210) * density(33) 
  pd(33,33) = pd(33,33) - rrt(210) * density(01) 
  pd(50,01) = pd(50,01) + rrt(210) * density(33) 
  pd(50,33) = pd(50,33) + rrt(210) * density(01) 
  pd(01,01) = pd(01,01) + rrt(211) * density(02) 
  pd(01,02) = pd(01,02) + rrt(211) * density(01) 
  pd(02,01) = pd(02,01) - rrt(211) * density(02) 
  pd(02,02) = pd(02,02) - rrt(211) * density(01) 
  pd(28,01) = pd(28,01) + rrt(211) * density(02) 
  pd(28,02) = pd(28,02) + rrt(211) * density(01) 
  pd(50,01) = pd(50,01) + rrt(211) * density(02) 
  pd(50,02) = pd(50,02) + rrt(211) * density(01) 
  pd(01,01) = pd(01,01) + rrt(212) * density(33) 
  pd(01,33) = pd(01,33) + rrt(212) * density(01) 
  pd(23,01) = pd(23,01) + rrt(212) * density(33) 
  pd(23,33) = pd(23,33) + rrt(212) * density(01) 
  pd(33,01) = pd(33,01) - rrt(212) * density(33) 
  pd(33,33) = pd(33,33) - rrt(212) * density(01) 
  pd(44,01) = pd(44,01) + rrt(212) * density(33) 
  pd(44,33) = pd(44,33) + rrt(212) * density(01) 
  pd(01,01) = pd(01,01) + rrt(213) * density(02) 
  pd(01,02) = pd(01,02) + rrt(213) * density(01) 
  pd(02,01) = pd(02,01) - rrt(213) * density(02) 
  pd(02,02) = pd(02,02) - rrt(213) * density(01) 
  pd(23,01) = pd(23,01) + rrt(213) * density(02) 
  pd(23,02) = pd(23,02) + rrt(213) * density(01) 
  pd(44,01) = pd(44,01) + rrt(213) * density(02) 
  pd(44,02) = pd(44,02) + rrt(213) * density(01) 
  pd(01,01) = pd(01,01) + rrt(214) * density(52) 
  pd(01,52) = pd(01,52) + rrt(214) * density(01) 
  pd(17,01) = pd(17,01) + rrt(214) * density(52) 
  pd(17,52) = pd(17,52) + rrt(214) * density(01) 
  pd(52,01) = pd(52,01) - rrt(214) * density(52) 
  pd(52,52) = pd(52,52) - rrt(214) * density(01) 
  pd(01,01) = pd(01,01) + rrt(215) * density(09) 
  pd(01,09) = pd(01,09) + rrt(215) * density(01) 
  pd(09,01) = pd(09,01) - rrt(215) * density(09) 
  pd(09,09) = pd(09,09) - rrt(215) * density(01) 
  pd(17,01) = pd(17,01) + rrt(215) * density(09) 
  pd(17,09) = pd(17,09) + rrt(215) * density(01) 
  pd(01,01) = pd(01,01) + rrt(216) * density(38) 
  pd(01,38) = pd(01,38) + rrt(216) * density(01) 
  pd(17,01) = pd(17,01) + rrt(216) * density(38) 
  pd(17,38) = pd(17,38) + rrt(216) * density(01) 
  pd(38,01) = pd(38,01) - rrt(216) * density(38) 
  pd(38,38) = pd(38,38) - rrt(216) * density(01) 
  pd(01,01) = pd(01,01) + rrt(217) * density(52) 
  pd(01,52) = pd(01,52) + rrt(217) * density(01) 
  pd(08,01) = pd(08,01) + rrt(217) * density(52) 
  pd(08,52) = pd(08,52) + rrt(217) * density(01) 
  pd(44,01) = pd(44,01) + rrt(217) * density(52) 
  pd(44,52) = pd(44,52) + rrt(217) * density(01) 
  pd(52,01) = pd(52,01) - rrt(217) * density(52) 
  pd(52,52) = pd(52,52) - rrt(217) * density(01) 
  pd(01,01) = pd(01,01) + rrt(218) * density(09) 
  pd(01,09) = pd(01,09) + rrt(218) * density(01) 
  pd(08,01) = pd(08,01) + rrt(218) * density(09) 
  pd(08,09) = pd(08,09) + rrt(218) * density(01) 
  pd(09,01) = pd(09,01) - rrt(218) * density(09) 
  pd(09,09) = pd(09,09) - rrt(218) * density(01) 
  pd(44,01) = pd(44,01) + rrt(218) * density(09) 
  pd(44,09) = pd(44,09) + rrt(218) * density(01) 
  pd(01,01) = pd(01,01) + rrt(219) * density(38) 
  pd(01,38) = pd(01,38) + rrt(219) * density(01) 
  pd(08,01) = pd(08,01) + rrt(219) * density(38) 
  pd(08,38) = pd(08,38) + rrt(219) * density(01) 
  pd(38,01) = pd(38,01) - rrt(219) * density(38) 
  pd(38,38) = pd(38,38) - rrt(219) * density(01) 
  pd(44,01) = pd(44,01) + rrt(219) * density(38) 
  pd(44,38) = pd(44,38) + rrt(219) * density(01) 
  pd(01,01) = pd(01,01) + rrt(220) * density(32) 
  pd(01,32) = pd(01,32) + rrt(220) * density(01) 
  pd(14,01) = pd(14,01) + rrt(220) * density(32) 
  pd(14,32) = pd(14,32) + rrt(220) * density(01) 
  pd(32,01) = pd(32,01) - rrt(220) * density(32) 
  pd(32,32) = pd(32,32) - rrt(220) * density(01) 
  pd(01,01) = pd(01,01) + rrt(221) * density(24) 
  pd(01,24) = pd(01,24) + rrt(221) * density(01) 
  pd(14,01) = pd(14,01) + rrt(221) * density(24) 
  pd(14,24) = pd(14,24) + rrt(221) * density(01) 
  pd(24,01) = pd(24,01) - rrt(221) * density(24) 
  pd(24,24) = pd(24,24) - rrt(221) * density(01) 
  pd(01,01) = pd(01,01) + rrt(222) * density(32) 
  pd(01,32) = pd(01,32) + rrt(222) * density(01) 
  pd(03,01) = pd(03,01) + rrt(222) * density(32) 
  pd(03,32) = pd(03,32) + rrt(222) * density(01) 
  pd(13,01) = pd(13,01) + rrt(222) * density(32) 
  pd(13,32) = pd(13,32) + rrt(222) * density(01) 
  pd(32,01) = pd(32,01) - rrt(222) * density(32) 
  pd(32,32) = pd(32,32) - rrt(222) * density(01) 
  pd(01,01) = pd(01,01) + rrt(223) * density(24) 
  pd(01,24) = pd(01,24) + rrt(223) * density(01) 
  pd(03,01) = pd(03,01) + rrt(223) * density(24) 
  pd(03,24) = pd(03,24) + rrt(223) * density(01) 
  pd(13,01) = pd(13,01) + rrt(223) * density(24) 
  pd(13,24) = pd(13,24) + rrt(223) * density(01) 
  pd(24,01) = pd(24,01) - rrt(223) * density(24) 
  pd(24,24) = pd(24,24) - rrt(223) * density(01) 
  pd(01,01) = pd(01,01) + rrt(224) * density(32) 
  pd(01,32) = pd(01,32) + rrt(224) * density(01) 
  pd(05,01) = pd(05,01) + rrt(224) * density(32) 
  pd(05,32) = pd(05,32) + rrt(224) * density(01) 
  pd(11,01) = pd(11,01) + rrt(224) * density(32) 
  pd(11,32) = pd(11,32) + rrt(224) * density(01) 
  pd(32,01) = pd(32,01) - rrt(224) * density(32) 
  pd(32,32) = pd(32,32) - rrt(224) * density(01) 
  pd(01,01) = pd(01,01) + rrt(225) * density(24) 
  pd(01,24) = pd(01,24) + rrt(225) * density(01) 
  pd(05,01) = pd(05,01) + rrt(225) * density(24) 
  pd(05,24) = pd(05,24) + rrt(225) * density(01) 
  pd(11,01) = pd(11,01) + rrt(225) * density(24) 
  pd(11,24) = pd(11,24) + rrt(225) * density(01) 
  pd(24,01) = pd(24,01) - rrt(225) * density(24) 
  pd(24,24) = pd(24,24) - rrt(225) * density(01) 
  pd(01,01) = pd(01,01) + rrt(226) * density(32) 
  pd(01,32) = pd(01,32) + rrt(226) * density(01) 
  pd(11,01) = pd(11,01) + rrt(226) * density(32) 
  pd(11,32) = pd(11,32) + rrt(226) * density(01) 
  pd(13,01) = pd(13,01) + rrt(226) * density(32) 
  pd(13,32) = pd(13,32) + rrt(226) * density(01) 
  pd(26,01) = pd(26,01) + rrt(226) * density(32) 
  pd(26,32) = pd(26,32) + rrt(226) * density(01) 
  pd(32,01) = pd(32,01) - rrt(226) * density(32) 
  pd(32,32) = pd(32,32) - rrt(226) * density(01) 
  pd(01,01) = pd(01,01) + rrt(227) * density(24) 
  pd(01,24) = pd(01,24) + rrt(227) * density(01) 
  pd(11,01) = pd(11,01) + rrt(227) * density(24) 
  pd(11,24) = pd(11,24) + rrt(227) * density(01) 
  pd(13,01) = pd(13,01) + rrt(227) * density(24) 
  pd(13,24) = pd(13,24) + rrt(227) * density(01) 
  pd(24,01) = pd(24,01) - rrt(227) * density(24) 
  pd(24,24) = pd(24,24) - rrt(227) * density(01) 
  pd(26,01) = pd(26,01) + rrt(227) * density(24) 
  pd(26,24) = pd(26,24) + rrt(227) * density(01) 
  pd(01,01) = pd(01,01) + rrt(228) * density(32) 
  pd(01,32) = pd(01,32) + rrt(228) * density(01) 
  pd(11,01) = pd(11,01) + rrt(228) * density(32) * 2.0d0
  pd(11,32) = pd(11,32) + rrt(228) * density(01) * 2.0d0
  pd(19,01) = pd(19,01) + rrt(228) * density(32) 
  pd(19,32) = pd(19,32) + rrt(228) * density(01) 
  pd(32,01) = pd(32,01) - rrt(228) * density(32) 
  pd(32,32) = pd(32,32) - rrt(228) * density(01) 
  pd(01,01) = pd(01,01) + rrt(229) * density(24) 
  pd(01,24) = pd(01,24) + rrt(229) * density(01) 
  pd(11,01) = pd(11,01) + rrt(229) * density(24) * 2.0d0
  pd(11,24) = pd(11,24) + rrt(229) * density(01) * 2.0d0
  pd(19,01) = pd(19,01) + rrt(229) * density(24) 
  pd(19,24) = pd(19,24) + rrt(229) * density(01) 
  pd(24,01) = pd(24,01) - rrt(229) * density(24) 
  pd(24,24) = pd(24,24) - rrt(229) * density(01) 
  pd(01,01) = pd(01,01) + rrt(230) * density(32) 
  pd(01,32) = pd(01,32) + rrt(230) * density(01) 
  pd(04,01) = pd(04,01) + rrt(230) * density(32) 
  pd(04,32) = pd(04,32) + rrt(230) * density(01) 
  pd(23,01) = pd(23,01) + rrt(230) * density(32) 
  pd(23,32) = pd(23,32) + rrt(230) * density(01) 
  pd(32,01) = pd(32,01) - rrt(230) * density(32) 
  pd(32,32) = pd(32,32) - rrt(230) * density(01) 
  pd(01,01) = pd(01,01) + rrt(231) * density(24) 
  pd(01,24) = pd(01,24) + rrt(231) * density(01) 
  pd(04,01) = pd(04,01) + rrt(231) * density(24) 
  pd(04,24) = pd(04,24) + rrt(231) * density(01) 
  pd(23,01) = pd(23,01) + rrt(231) * density(24) 
  pd(23,24) = pd(23,24) + rrt(231) * density(01) 
  pd(24,01) = pd(24,01) - rrt(231) * density(24) 
  pd(24,24) = pd(24,24) - rrt(231) * density(01) 
  pd(01,01) = pd(01,01) + rrt(232) * density(32) 
  pd(01,32) = pd(01,32) + rrt(232) * density(01) 
  pd(32,01) = pd(32,01) - rrt(232) * density(32) 
  pd(32,32) = pd(32,32) - rrt(232) * density(01) 
  pd(46,01) = pd(46,01) + rrt(232) * density(32) 
  pd(46,32) = pd(46,32) + rrt(232) * density(01) 
  pd(51,01) = pd(51,01) + rrt(232) * density(32) 
  pd(51,32) = pd(51,32) + rrt(232) * density(01) 
  pd(01,01) = pd(01,01) + rrt(233) * density(24) 
  pd(01,24) = pd(01,24) + rrt(233) * density(01) 
  pd(24,01) = pd(24,01) - rrt(233) * density(24) 
  pd(24,24) = pd(24,24) - rrt(233) * density(01) 
  pd(46,01) = pd(46,01) + rrt(233) * density(24) 
  pd(46,24) = pd(46,24) + rrt(233) * density(01) 
  pd(51,01) = pd(51,01) + rrt(233) * density(24) 
  pd(51,24) = pd(51,24) + rrt(233) * density(01) 
  pd(01,01) = pd(01,01) + rrt(234) * density(32) 
  pd(01,32) = pd(01,32) + rrt(234) * density(01) 
  pd(16,01) = pd(16,01) + rrt(234) * density(32) 
  pd(16,32) = pd(16,32) + rrt(234) * density(01) 
  pd(32,01) = pd(32,01) - rrt(234) * density(32) 
  pd(32,32) = pd(32,32) - rrt(234) * density(01) 
  pd(34,01) = pd(34,01) + rrt(234) * density(32) 
  pd(34,32) = pd(34,32) + rrt(234) * density(01) 
  pd(01,01) = pd(01,01) + rrt(235) * density(24) 
  pd(01,24) = pd(01,24) + rrt(235) * density(01) 
  pd(16,01) = pd(16,01) + rrt(235) * density(24) 
  pd(16,24) = pd(16,24) + rrt(235) * density(01) 
  pd(24,01) = pd(24,01) - rrt(235) * density(24) 
  pd(24,24) = pd(24,24) - rrt(235) * density(01) 
  pd(34,01) = pd(34,01) + rrt(235) * density(24) 
  pd(34,24) = pd(34,24) + rrt(235) * density(01) 
  pd(01,01) = pd(01,01) + rrt(236) * density(32) 
  pd(01,32) = pd(01,32) + rrt(236) * density(01) 
  pd(32,01) = pd(32,01) - rrt(236) * density(32) 
  pd(32,32) = pd(32,32) - rrt(236) * density(01) 
  pd(43,01) = pd(43,01) + rrt(236) * density(32) 
  pd(43,32) = pd(43,32) + rrt(236) * density(01) 
  pd(50,01) = pd(50,01) + rrt(236) * density(32) 
  pd(50,32) = pd(50,32) + rrt(236) * density(01) 
  pd(01,01) = pd(01,01) + rrt(237) * density(24) 
  pd(01,24) = pd(01,24) + rrt(237) * density(01) 
  pd(24,01) = pd(24,01) - rrt(237) * density(24) 
  pd(24,24) = pd(24,24) - rrt(237) * density(01) 
  pd(43,01) = pd(43,01) + rrt(237) * density(24) 
  pd(43,24) = pd(43,24) + rrt(237) * density(01) 
  pd(50,01) = pd(50,01) + rrt(237) * density(24) 
  pd(50,24) = pd(50,24) + rrt(237) * density(01) 
  pd(11,01) = pd(11,01) - rrt(238) * density(11) 
  pd(11,11) = pd(11,11) - rrt(238) * density(01) 
  pd(13,01) = pd(13,01) + rrt(238) * density(11) * 2.0d0
  pd(13,11) = pd(13,11) + rrt(238) * density(01) * 2.0d0
  pd(01,01) = pd(01,01) + rrt(239) * density(11) 
  pd(01,11) = pd(01,11) + rrt(239) * density(01) 
  pd(11,01) = pd(11,01) - rrt(239) * density(11) 
  pd(11,11) = pd(11,11) - rrt(239) * density(01) 
  pd(27,01) = pd(27,01) + rrt(239) * density(11) 
  pd(27,11) = pd(27,11) + rrt(239) * density(01) 
  pd(01,01) = pd(01,01) - rrt(240) * density(01) * density(10) * 2.0d0
  pd(01,10) = pd(01,10) - rrt(240) * density(01)**2 
  pd(10,01) = pd(10,01) - rrt(240) * density(01) * density(10) * 2.0d0
  pd(10,10) = pd(10,10) - rrt(240) * density(01)**2 
  pd(51,01) = pd(51,01) + rrt(240) * density(01) * density(10) * 2.0d0
  pd(51,10) = pd(51,10) + rrt(240) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(241) * density(01) * density(16) * 2.0d0
  pd(01,16) = pd(01,16) - rrt(241) * density(01)**2 
  pd(16,01) = pd(16,01) - rrt(241) * density(01) * density(16) * 2.0d0
  pd(16,16) = pd(16,16) - rrt(241) * density(01)**2 
  pd(23,01) = pd(23,01) + rrt(241) * density(01) * density(16) * 2.0d0
  pd(23,16) = pd(23,16) + rrt(241) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(242) * density(01) * density(50) * 2.0d0
  pd(01,50) = pd(01,50) - rrt(242) * density(01)**2 
  pd(28,01) = pd(28,01) + rrt(242) * density(01) * density(50) * 2.0d0
  pd(28,50) = pd(28,50) + rrt(242) * density(01)**2 
  pd(50,01) = pd(50,01) - rrt(242) * density(01) * density(50) * 2.0d0
  pd(50,50) = pd(50,50) - rrt(242) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(243) * density(01) * density(44) * 2.0d0
  pd(01,44) = pd(01,44) - rrt(243) * density(01)**2 
  pd(08,01) = pd(08,01) + rrt(243) * density(01) * density(44) * 2.0d0
  pd(08,44) = pd(08,44) + rrt(243) * density(01)**2 
  pd(44,01) = pd(44,01) - rrt(243) * density(01) * density(44) * 2.0d0
  pd(44,44) = pd(44,44) - rrt(243) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(244) * density(01) * density(40) * 2.0d0
  pd(01,40) = pd(01,40) - rrt(244) * density(01)**2 
  pd(40,01) = pd(40,01) - rrt(244) * density(01) * density(40) * 2.0d0
  pd(40,40) = pd(40,40) - rrt(244) * density(01)**2 
  pd(43,01) = pd(43,01) + rrt(244) * density(01) * density(40) * 2.0d0
  pd(43,40) = pd(43,40) + rrt(244) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(245) * density(01) * density(04) * 2.0d0
  pd(01,04) = pd(01,04) - rrt(245) * density(01)**2 
  pd(04,01) = pd(04,01) - rrt(245) * density(01) * density(04) * 2.0d0
  pd(04,04) = pd(04,04) - rrt(245) * density(01)**2 
  pd(34,01) = pd(34,01) + rrt(245) * density(01) * density(04) * 2.0d0
  pd(34,04) = pd(34,04) + rrt(245) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(246) * density(01) * density(46) * 2.0d0
  pd(01,46) = pd(01,46) - rrt(246) * density(01)**2 
  pd(42,01) = pd(42,01) + rrt(246) * density(01) * density(46) * 2.0d0
  pd(42,46) = pd(42,46) + rrt(246) * density(01)**2 
  pd(46,01) = pd(46,01) - rrt(246) * density(01) * density(46) * 2.0d0
  pd(46,46) = pd(46,46) - rrt(246) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(247) * density(01) * density(07) * 2.0d0
  pd(01,07) = pd(01,07) - rrt(247) * density(01)**2 
  pd(07,01) = pd(07,01) - rrt(247) * density(01) * density(07) * 2.0d0
  pd(07,07) = pd(07,07) - rrt(247) * density(01)**2 
  pd(29,01) = pd(29,01) + rrt(247) * density(01) * density(07) * 2.0d0
  pd(29,07) = pd(29,07) + rrt(247) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(248) * density(01) * density(17) * 2.0d0
  pd(01,17) = pd(01,17) - rrt(248) * density(01)**2 
  pd(17,01) = pd(17,01) - rrt(248) * density(01) * density(17) * 2.0d0
  pd(17,17) = pd(17,17) - rrt(248) * density(01)**2 
  pd(21,01) = pd(21,01) + rrt(248) * density(01) * density(17) * 2.0d0
  pd(21,17) = pd(21,17) + rrt(248) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(249) * density(01) * density(14) * 2.0d0
  pd(01,14) = pd(01,14) - rrt(249) * density(01)**2 
  pd(14,01) = pd(14,01) - rrt(249) * density(01) * density(14) * 2.0d0
  pd(14,14) = pd(14,14) - rrt(249) * density(01)**2 
  pd(45,01) = pd(45,01) + rrt(249) * density(01) * density(14) * 2.0d0
  pd(45,14) = pd(45,14) + rrt(249) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(250) * density(01) * density(03) * 2.0d0
  pd(01,03) = pd(01,03) - rrt(250) * density(01)**2 
  pd(03,01) = pd(03,01) - rrt(250) * density(01) * density(03) * 2.0d0
  pd(03,03) = pd(03,03) - rrt(250) * density(01)**2 
  pd(18,01) = pd(18,01) + rrt(250) * density(01) * density(03) * 2.0d0
  pd(18,03) = pd(18,03) + rrt(250) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(251) * density(01) * density(05) * 2.0d0
  pd(01,05) = pd(01,05) - rrt(251) * density(01)**2 
  pd(05,01) = pd(05,01) - rrt(251) * density(01) * density(05) * 2.0d0
  pd(05,05) = pd(05,05) - rrt(251) * density(01)**2 
  pd(37,01) = pd(37,01) + rrt(251) * density(01) * density(05) * 2.0d0
  pd(37,05) = pd(37,05) + rrt(251) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(252) * density(01) * density(26) * 2.0d0
  pd(01,26) = pd(01,26) - rrt(252) * density(01)**2 
  pd(20,01) = pd(20,01) + rrt(252) * density(01) * density(26) * 2.0d0
  pd(20,26) = pd(20,26) + rrt(252) * density(01)**2 
  pd(26,01) = pd(26,01) - rrt(252) * density(01) * density(26) * 2.0d0
  pd(26,26) = pd(26,26) - rrt(252) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(253) * density(01) * density(19) * 2.0d0
  pd(01,19) = pd(01,19) - rrt(253) * density(01)**2 
  pd(19,01) = pd(19,01) - rrt(253) * density(01) * density(19) * 2.0d0
  pd(19,19) = pd(19,19) - rrt(253) * density(01)**2 
  pd(35,01) = pd(35,01) + rrt(253) * density(01) * density(19) * 2.0d0
  pd(35,19) = pd(35,19) + rrt(253) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(254) * density(01) * density(30) * 2.0d0
  pd(01,30) = pd(01,30) - rrt(254) * density(01)**2 
  pd(13,01) = pd(13,01) + rrt(254) * density(01) * density(30) * 2.0d0
  pd(13,30) = pd(13,30) + rrt(254) * density(01)**2 
  pd(30,01) = pd(30,01) - rrt(254) * density(01) * density(30) * 2.0d0
  pd(30,30) = pd(30,30) - rrt(254) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(255) * density(01) * density(27) * 2.0d0
  pd(01,27) = pd(01,27) - rrt(255) * density(01)**2 
  pd(11,01) = pd(11,01) + rrt(255) * density(01) * density(27) * 2.0d0
  pd(11,27) = pd(11,27) + rrt(255) * density(01)**2 
  pd(27,01) = pd(27,01) - rrt(255) * density(01) * density(27) * 2.0d0
  pd(27,27) = pd(27,27) - rrt(255) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(256) * density(47) 
  pd(01,47) = pd(01,47) - rrt(256) * density(01) 
  pd(13,01) = pd(13,01) + rrt(256) * density(47) * 2.0d0
  pd(13,47) = pd(13,47) + rrt(256) * density(01) * 2.0d0
  pd(23,01) = pd(23,01) + rrt(256) * density(47) 
  pd(23,47) = pd(23,47) + rrt(256) * density(01) 
  pd(47,01) = pd(47,01) - rrt(256) * density(47) 
  pd(47,47) = pd(47,47) - rrt(256) * density(01) 
  pd(01,01) = pd(01,01) - rrt(257) * density(47) 
  pd(01,47) = pd(01,47) - rrt(257) * density(01) 
  pd(11,01) = pd(11,01) + rrt(257) * density(47) 
  pd(11,47) = pd(11,47) + rrt(257) * density(01) 
  pd(13,01) = pd(13,01) + rrt(257) * density(47) 
  pd(13,47) = pd(13,47) + rrt(257) * density(01) 
  pd(28,01) = pd(28,01) + rrt(257) * density(47) 
  pd(28,47) = pd(28,47) + rrt(257) * density(01) 
  pd(47,01) = pd(47,01) - rrt(257) * density(47) 
  pd(47,47) = pd(47,47) - rrt(257) * density(01) 
  pd(01,01) = pd(01,01) - rrt(258) * density(10) 
  pd(01,10) = pd(01,10) - rrt(258) * density(01) 
  pd(10,01) = pd(10,01) - rrt(258) * density(10) 
  pd(10,10) = pd(10,10) - rrt(258) * density(01) 
  pd(13,01) = pd(13,01) + rrt(258) * density(10) 
  pd(13,10) = pd(13,10) + rrt(258) * density(01) 
  pd(23,01) = pd(23,01) + rrt(258) * density(10) 
  pd(23,10) = pd(23,10) + rrt(258) * density(01) 
  pd(01,01) = pd(01,01) - rrt(259) * density(10) 
  pd(01,10) = pd(01,10) - rrt(259) * density(01) 
  pd(10,01) = pd(10,01) - rrt(259) * density(10) 
  pd(10,10) = pd(10,10) - rrt(259) * density(01) 
  pd(13,01) = pd(13,01) + rrt(259) * density(10) * 2.0d0
  pd(13,10) = pd(13,10) + rrt(259) * density(01) * 2.0d0
  pd(28,01) = pd(28,01) + rrt(259) * density(10) 
  pd(28,10) = pd(28,10) + rrt(259) * density(01) 
  pd(01,01) = pd(01,01) - rrt(260) * density(10) 
  pd(01,10) = pd(01,10) - rrt(260) * density(01) 
  pd(08,01) = pd(08,01) + rrt(260) * density(10) 
  pd(08,10) = pd(08,10) + rrt(260) * density(01) 
  pd(10,01) = pd(10,01) - rrt(260) * density(10) 
  pd(10,10) = pd(10,10) - rrt(260) * density(01) 
  pd(11,01) = pd(11,01) + rrt(260) * density(10) 
  pd(11,10) = pd(11,10) + rrt(260) * density(01) 
  pd(13,01) = pd(13,01) + rrt(260) * density(10) 
  pd(13,10) = pd(13,10) + rrt(260) * density(01) 
  pd(01,01) = pd(01,01) - rrt(261) * density(16) 
  pd(01,16) = pd(01,16) - rrt(261) * density(01) 
  pd(13,01) = pd(13,01) + rrt(261) * density(16) 
  pd(13,16) = pd(13,16) + rrt(261) * density(01) 
  pd(16,01) = pd(16,01) - rrt(261) * density(16) 
  pd(16,16) = pd(16,16) - rrt(261) * density(01) 
  pd(28,01) = pd(28,01) + rrt(261) * density(16) 
  pd(28,16) = pd(28,16) + rrt(261) * density(01) 
  pd(01,01) = pd(01,01) - rrt(262) * density(16) 
  pd(01,16) = pd(01,16) - rrt(262) * density(01) 
  pd(08,01) = pd(08,01) + rrt(262) * density(16) 
  pd(08,16) = pd(08,16) + rrt(262) * density(01) 
  pd(11,01) = pd(11,01) + rrt(262) * density(16) 
  pd(11,16) = pd(11,16) + rrt(262) * density(01) 
  pd(16,01) = pd(16,01) - rrt(262) * density(16) 
  pd(16,16) = pd(16,16) - rrt(262) * density(01) 
  pd(01,01) = pd(01,01) - rrt(263) * density(50) 
  pd(01,50) = pd(01,50) - rrt(263) * density(01) 
  pd(08,01) = pd(08,01) + rrt(263) * density(50) 
  pd(08,50) = pd(08,50) + rrt(263) * density(01) 
  pd(13,01) = pd(13,01) + rrt(263) * density(50) 
  pd(13,50) = pd(13,50) + rrt(263) * density(01) 
  pd(50,01) = pd(50,01) - rrt(263) * density(50) 
  pd(50,50) = pd(50,50) - rrt(263) * density(01) 
  pd(01,01) = pd(01,01) - rrt(264) * density(40) 
  pd(01,40) = pd(01,40) - rrt(264) * density(01) 
  pd(13,01) = pd(13,01) + rrt(264) * density(40) 
  pd(13,40) = pd(13,40) + rrt(264) * density(01) 
  pd(34,01) = pd(34,01) + rrt(264) * density(40) 
  pd(34,40) = pd(34,40) + rrt(264) * density(01) 
  pd(40,01) = pd(40,01) - rrt(264) * density(40) 
  pd(40,40) = pd(40,40) - rrt(264) * density(01) 
  pd(01,01) = pd(01,01) - rrt(265) * density(40) 
  pd(01,40) = pd(01,40) - rrt(265) * density(01) 
  pd(13,01) = pd(13,01) + rrt(265) * density(40) * 2.0d0
  pd(13,40) = pd(13,40) + rrt(265) * density(01) * 2.0d0
  pd(40,01) = pd(40,01) - rrt(265) * density(40) 
  pd(40,40) = pd(40,40) - rrt(265) * density(01) 
  pd(42,01) = pd(42,01) + rrt(265) * density(40) 
  pd(42,40) = pd(42,40) + rrt(265) * density(01) 
  pd(01,01) = pd(01,01) - rrt(266) * density(04) 
  pd(01,04) = pd(01,04) - rrt(266) * density(01) 
  pd(04,01) = pd(04,01) - rrt(266) * density(04) 
  pd(04,04) = pd(04,04) - rrt(266) * density(01) 
  pd(13,01) = pd(13,01) + rrt(266) * density(04) 
  pd(13,04) = pd(13,04) + rrt(266) * density(01) 
  pd(42,01) = pd(42,01) + rrt(266) * density(04) 
  pd(42,04) = pd(42,04) + rrt(266) * density(01) 
  pd(01,01) = pd(01,01) - rrt(267) * density(04) 
  pd(01,04) = pd(01,04) - rrt(267) * density(01) 
  pd(04,01) = pd(04,01) - rrt(267) * density(04) 
  pd(04,04) = pd(04,04) - rrt(267) * density(01) 
  pd(13,01) = pd(13,01) + rrt(267) * density(04) * 2.0d0
  pd(13,04) = pd(13,04) + rrt(267) * density(01) * 2.0d0
  pd(29,01) = pd(29,01) + rrt(267) * density(04) 
  pd(29,04) = pd(29,04) + rrt(267) * density(01) 
  pd(01,01) = pd(01,01) - rrt(268) * density(04) 
  pd(01,04) = pd(01,04) - rrt(268) * density(01) 
  pd(04,01) = pd(04,01) - rrt(268) * density(04) 
  pd(04,04) = pd(04,04) - rrt(268) * density(01) 
  pd(11,01) = pd(11,01) + rrt(268) * density(04) 
  pd(11,04) = pd(11,04) + rrt(268) * density(01) 
  pd(13,01) = pd(13,01) + rrt(268) * density(04) 
  pd(13,04) = pd(13,04) + rrt(268) * density(01) 
  pd(21,01) = pd(21,01) + rrt(268) * density(04) 
  pd(21,04) = pd(21,04) + rrt(268) * density(01) 
  pd(01,01) = pd(01,01) - rrt(269) * density(04) 
  pd(01,04) = pd(01,04) - rrt(269) * density(01) 
  pd(04,01) = pd(04,01) - rrt(269) * density(04) 
  pd(04,04) = pd(04,04) - rrt(269) * density(01) 
  pd(13,01) = pd(13,01) + rrt(269) * density(04) * 3.0d0
  pd(13,04) = pd(13,04) + rrt(269) * density(01) * 3.0d0
  pd(21,01) = pd(21,01) + rrt(269) * density(04) 
  pd(21,04) = pd(21,04) + rrt(269) * density(01) 
  pd(01,01) = pd(01,01) - rrt(270) * density(04) 
  pd(01,04) = pd(01,04) - rrt(270) * density(01) 
  pd(04,01) = pd(04,01) - rrt(270) * density(04) 
  pd(04,04) = pd(04,04) - rrt(270) * density(01) 
  pd(23,01) = pd(23,01) + rrt(270) * density(04) 
  pd(23,04) = pd(23,04) + rrt(270) * density(01) 
  pd(28,01) = pd(28,01) + rrt(270) * density(04) 
  pd(28,04) = pd(28,04) + rrt(270) * density(01) 
  pd(01,01) = pd(01,01) - rrt(271) * density(46) 
  pd(01,46) = pd(01,46) - rrt(271) * density(01) 
  pd(13,01) = pd(13,01) + rrt(271) * density(46) 
  pd(13,46) = pd(13,46) + rrt(271) * density(01) 
  pd(29,01) = pd(29,01) + rrt(271) * density(46) 
  pd(29,46) = pd(29,46) + rrt(271) * density(01) 
  pd(46,01) = pd(46,01) - rrt(271) * density(46) 
  pd(46,46) = pd(46,46) - rrt(271) * density(01) 
  pd(01,01) = pd(01,01) - rrt(272) * density(46) 
  pd(01,46) = pd(01,46) - rrt(272) * density(01) 
  pd(13,01) = pd(13,01) + rrt(272) * density(46) * 2.0d0
  pd(13,46) = pd(13,46) + rrt(272) * density(01) * 2.0d0
  pd(21,01) = pd(21,01) + rrt(272) * density(46) 
  pd(21,46) = pd(21,46) + rrt(272) * density(01) 
  pd(46,01) = pd(46,01) - rrt(272) * density(46) 
  pd(46,46) = pd(46,46) - rrt(272) * density(01) 
  pd(01,01) = pd(01,01) - rrt(273) * density(07) 
  pd(01,07) = pd(01,07) - rrt(273) * density(01) 
  pd(07,01) = pd(07,01) - rrt(273) * density(07) 
  pd(07,07) = pd(07,07) - rrt(273) * density(01) 
  pd(13,01) = pd(13,01) + rrt(273) * density(07) 
  pd(13,07) = pd(13,07) + rrt(273) * density(01) 
  pd(21,01) = pd(21,01) + rrt(273) * density(07) 
  pd(21,07) = pd(21,07) + rrt(273) * density(01) 
  pd(01,01) = pd(01,01) - rrt(274) * density(17) 
  pd(01,17) = pd(01,17) - rrt(274) * density(01) 
  pd(08,01) = pd(08,01) + rrt(274) * density(17) * 2.0d0
  pd(08,17) = pd(08,17) + rrt(274) * density(01) * 2.0d0
  pd(17,01) = pd(17,01) - rrt(274) * density(17) 
  pd(17,17) = pd(17,17) - rrt(274) * density(01) 
  pd(01,01) = pd(01,01) - rrt(275) * density(27) 
  pd(01,27) = pd(01,27) - rrt(275) * density(01) 
  pd(13,01) = pd(13,01) + rrt(275) * density(27) * 2.0d0
  pd(13,27) = pd(13,27) + rrt(275) * density(01) * 2.0d0
  pd(27,01) = pd(27,01) - rrt(275) * density(27) 
  pd(27,27) = pd(27,27) - rrt(275) * density(01) 
  pd(01,01) = pd(01,01) - rrt(276) * density(12) 
  pd(01,12) = pd(01,12) - rrt(276) * density(01) 
  pd(11,01) = pd(11,01) + rrt(276) * density(12) 
  pd(11,12) = pd(11,12) + rrt(276) * density(01) 
  pd(12,01) = pd(12,01) - rrt(276) * density(12) 
  pd(12,12) = pd(12,12) - rrt(276) * density(01) 
  pd(13,01) = pd(13,01) + rrt(276) * density(12) 
  pd(13,12) = pd(13,12) + rrt(276) * density(01) 
  pd(11,11) = pd(11,11) - rrt(277) * density(27) 
  pd(11,27) = pd(11,27) - rrt(277) * density(11) 
  pd(12,11) = pd(12,11) + rrt(277) * density(27) 
  pd(12,27) = pd(12,27) + rrt(277) * density(11) 
  pd(13,11) = pd(13,11) + rrt(277) * density(27) 
  pd(13,27) = pd(13,27) + rrt(277) * density(11) 
  pd(27,11) = pd(27,11) - rrt(277) * density(27) 
  pd(27,27) = pd(27,27) - rrt(277) * density(11) 
  pd(16,28) = pd(16,28) + rrt(278) * density(47) 
  pd(16,47) = pd(16,47) + rrt(278) * density(28) 
  pd(28,28) = pd(28,28) - rrt(278) * density(47) 
  pd(28,47) = pd(28,47) - rrt(278) * density(28) 
  pd(47,28) = pd(47,28) - rrt(278) * density(47) 
  pd(47,47) = pd(47,47) - rrt(278) * density(28) 
  pd(51,28) = pd(51,28) + rrt(278) * density(47) 
  pd(51,47) = pd(51,47) + rrt(278) * density(28) 
  pd(08,08) = pd(08,08) - rrt(279) * density(47) 
  pd(08,47) = pd(08,47) - rrt(279) * density(08) 
  pd(47,08) = pd(47,08) - rrt(279) * density(47) 
  pd(47,47) = pd(47,47) - rrt(279) * density(08) 
  pd(50,08) = pd(50,08) + rrt(279) * density(47) 
  pd(50,47) = pd(50,47) + rrt(279) * density(08) 
  pd(51,08) = pd(51,08) + rrt(279) * density(47) 
  pd(51,47) = pd(51,47) + rrt(279) * density(08) 
  pd(04,43) = pd(04,43) + rrt(280) * density(47) 
  pd(04,47) = pd(04,47) + rrt(280) * density(43) 
  pd(11,43) = pd(11,43) + rrt(280) * density(47) 
  pd(11,47) = pd(11,47) + rrt(280) * density(43) 
  pd(43,43) = pd(43,43) - rrt(280) * density(47) 
  pd(43,47) = pd(43,47) - rrt(280) * density(43) 
  pd(47,43) = pd(47,43) - rrt(280) * density(47) 
  pd(47,47) = pd(47,47) - rrt(280) * density(43) 
  pd(51,43) = pd(51,43) + rrt(280) * density(47) 
  pd(51,47) = pd(51,47) + rrt(280) * density(43) 
  pd(04,42) = pd(04,42) + rrt(281) * density(47) 
  pd(04,47) = pd(04,47) + rrt(281) * density(42) 
  pd(42,42) = pd(42,42) - rrt(281) * density(47) 
  pd(42,47) = pd(42,47) - rrt(281) * density(42) 
  pd(47,42) = pd(47,42) - rrt(281) * density(47) 
  pd(47,47) = pd(47,47) - rrt(281) * density(42) 
  pd(51,42) = pd(51,42) + rrt(281) * density(47) 
  pd(51,47) = pd(51,47) + rrt(281) * density(42) 
  pd(07,21) = pd(07,21) + rrt(282) * density(47) 
  pd(07,47) = pd(07,47) + rrt(282) * density(21) 
  pd(21,21) = pd(21,21) - rrt(282) * density(47) 
  pd(21,47) = pd(21,47) - rrt(282) * density(21) 
  pd(47,21) = pd(47,21) - rrt(282) * density(47) 
  pd(47,47) = pd(47,47) - rrt(282) * density(21) 
  pd(51,21) = pd(51,21) + rrt(282) * density(47) 
  pd(51,47) = pd(51,47) + rrt(282) * density(21) 
  pd(10,13) = pd(10,13) + rrt(283) * density(47) 
  pd(10,47) = pd(10,47) + rrt(283) * density(13) 
  pd(11,13) = pd(11,13) + rrt(283) * density(47) 
  pd(11,47) = pd(11,47) + rrt(283) * density(13) 
  pd(13,13) = pd(13,13) - rrt(283) * density(47) 
  pd(13,47) = pd(13,47) - rrt(283) * density(13) 
  pd(47,13) = pd(47,13) - rrt(283) * density(47) 
  pd(47,47) = pd(47,47) - rrt(283) * density(13) 
  pd(10,10) = pd(10,10) - rrt(284) * density(51) 
  pd(10,51) = pd(10,51) - rrt(284) * density(10) 
  pd(23,10) = pd(23,10) + rrt(284) * density(51) 
  pd(23,51) = pd(23,51) + rrt(284) * density(10) 
  pd(47,10) = pd(47,10) + rrt(284) * density(51) 
  pd(47,51) = pd(47,51) + rrt(284) * density(10) 
  pd(51,10) = pd(51,10) - rrt(284) * density(51) 
  pd(51,51) = pd(51,51) - rrt(284) * density(10) 
  pd(10,10) = pd(10,10) - rrt(285) * density(43) 
  pd(10,43) = pd(10,43) - rrt(285) * density(10) 
  pd(11,10) = pd(11,10) + rrt(285) * density(43) 
  pd(11,43) = pd(11,43) + rrt(285) * density(10) 
  pd(43,10) = pd(43,10) - rrt(285) * density(43) 
  pd(43,43) = pd(43,43) - rrt(285) * density(10) 
  pd(46,10) = pd(46,10) + rrt(285) * density(43) 
  pd(46,43) = pd(46,43) + rrt(285) * density(10) 
  pd(51,10) = pd(51,10) + rrt(285) * density(43) 
  pd(51,43) = pd(51,43) + rrt(285) * density(10) 
  pd(04,10) = pd(04,10) + rrt(286) * density(42) 
  pd(04,42) = pd(04,42) + rrt(286) * density(10) 
  pd(10,10) = pd(10,10) - rrt(286) * density(42) 
  pd(10,42) = pd(10,42) - rrt(286) * density(10) 
  pd(23,10) = pd(23,10) + rrt(286) * density(42) 
  pd(23,42) = pd(23,42) + rrt(286) * density(10) 
  pd(42,10) = pd(42,10) - rrt(286) * density(42) 
  pd(42,42) = pd(42,42) - rrt(286) * density(10) 
  pd(10,10) = pd(10,10) - rrt(287) * density(42) 
  pd(10,42) = pd(10,42) - rrt(287) * density(10) 
  pd(42,10) = pd(42,10) - rrt(287) * density(42) 
  pd(42,42) = pd(42,42) - rrt(287) * density(10) 
  pd(46,10) = pd(46,10) + rrt(287) * density(42) 
  pd(46,42) = pd(46,42) + rrt(287) * density(10) 
  pd(51,10) = pd(51,10) + rrt(287) * density(42) 
  pd(51,42) = pd(51,42) + rrt(287) * density(10) 
  pd(07,10) = pd(07,10) + rrt(288) * density(21) 
  pd(07,21) = pd(07,21) + rrt(288) * density(10) 
  pd(10,10) = pd(10,10) - rrt(288) * density(21) 
  pd(10,21) = pd(10,21) - rrt(288) * density(10) 
  pd(21,10) = pd(21,10) - rrt(288) * density(21) 
  pd(21,21) = pd(21,21) - rrt(288) * density(10) 
  pd(23,10) = pd(23,10) + rrt(288) * density(21) 
  pd(23,21) = pd(23,21) + rrt(288) * density(10) 
  pd(10,10) = pd(10,10) - rrt(289) * density(21) 
  pd(10,21) = pd(10,21) - rrt(289) * density(10) 
  pd(17,10) = pd(17,10) + rrt(289) * density(21) 
  pd(17,21) = pd(17,21) + rrt(289) * density(10) 
  pd(21,10) = pd(21,10) - rrt(289) * density(21) 
  pd(21,21) = pd(21,21) - rrt(289) * density(10) 
  pd(51,10) = pd(51,10) + rrt(289) * density(21) 
  pd(51,21) = pd(51,21) + rrt(289) * density(10) 
  pd(10,10) = pd(10,10) - rrt(290) * density(11) 
  pd(10,11) = pd(10,11) - rrt(290) * density(10) 
  pd(11,10) = pd(11,10) - rrt(290) * density(11) 
  pd(11,11) = pd(11,11) - rrt(290) * density(10) 
  pd(13,10) = pd(13,10) + rrt(290) * density(11) 
  pd(13,11) = pd(13,11) + rrt(290) * density(10) 
  pd(47,10) = pd(47,10) + rrt(290) * density(11) 
  pd(47,11) = pd(47,11) + rrt(290) * density(10) 
  pd(10,10) = pd(10,10) - rrt(291) * density(13) 
  pd(10,13) = pd(10,13) - rrt(291) * density(10) 
  pd(11,10) = pd(11,10) + rrt(291) * density(13) 
  pd(11,13) = pd(11,13) + rrt(291) * density(10) 
  pd(13,10) = pd(13,10) - rrt(291) * density(13) 
  pd(13,13) = pd(13,13) - rrt(291) * density(10) 
  pd(16,10) = pd(16,10) + rrt(291) * density(13) 
  pd(16,13) = pd(16,13) + rrt(291) * density(10) 
  pd(10,16) = pd(10,16) + rrt(292) * density(51) 
  pd(10,51) = pd(10,51) + rrt(292) * density(16) 
  pd(16,16) = pd(16,16) - rrt(292) * density(51) 
  pd(16,51) = pd(16,51) - rrt(292) * density(16) 
  pd(23,16) = pd(23,16) + rrt(292) * density(51) 
  pd(23,51) = pd(23,51) + rrt(292) * density(16) 
  pd(51,16) = pd(51,16) - rrt(292) * density(51) 
  pd(51,51) = pd(51,51) - rrt(292) * density(16) 
  pd(04,16) = pd(04,16) + rrt(293) * density(51) 
  pd(04,51) = pd(04,51) + rrt(293) * density(16) 
  pd(11,16) = pd(11,16) + rrt(293) * density(51) 
  pd(11,51) = pd(11,51) + rrt(293) * density(16) 
  pd(16,16) = pd(16,16) - rrt(293) * density(51) 
  pd(16,51) = pd(16,51) - rrt(293) * density(16) 
  pd(51,16) = pd(51,16) - rrt(293) * density(51) 
  pd(51,51) = pd(51,51) - rrt(293) * density(16) 
  pd(07,16) = pd(07,16) + rrt(294) * density(28) 
  pd(07,28) = pd(07,28) + rrt(294) * density(16) 
  pd(11,16) = pd(11,16) + rrt(294) * density(28) 
  pd(11,28) = pd(11,28) + rrt(294) * density(16) 
  pd(16,16) = pd(16,16) - rrt(294) * density(28) 
  pd(16,28) = pd(16,28) - rrt(294) * density(16) 
  pd(28,16) = pd(28,16) - rrt(294) * density(28) 
  pd(28,28) = pd(28,28) - rrt(294) * density(16) 
  pd(08,08) = pd(08,08) - rrt(295) * density(16) 
  pd(08,16) = pd(08,16) - rrt(295) * density(08) 
  pd(11,08) = pd(11,08) + rrt(295) * density(16) 
  pd(11,16) = pd(11,16) + rrt(295) * density(08) 
  pd(16,08) = pd(16,08) - rrt(295) * density(16) 
  pd(16,16) = pd(16,16) - rrt(295) * density(08) 
  pd(17,08) = pd(17,08) + rrt(295) * density(16) 
  pd(17,16) = pd(17,16) + rrt(295) * density(08) 
  pd(04,16) = pd(04,16) + rrt(296) * density(43) 
  pd(04,43) = pd(04,43) + rrt(296) * density(16) 
  pd(16,16) = pd(16,16) - rrt(296) * density(43) 
  pd(16,43) = pd(16,43) - rrt(296) * density(16) 
  pd(43,16) = pd(43,16) - rrt(296) * density(43) 
  pd(43,43) = pd(43,43) - rrt(296) * density(16) 
  pd(51,16) = pd(51,16) + rrt(296) * density(43) 
  pd(51,43) = pd(51,43) + rrt(296) * density(16) 
  pd(07,16) = pd(07,16) + rrt(297) * density(42) 
  pd(07,42) = pd(07,42) + rrt(297) * density(16) 
  pd(16,16) = pd(16,16) - rrt(297) * density(42) 
  pd(16,42) = pd(16,42) - rrt(297) * density(16) 
  pd(42,16) = pd(42,16) - rrt(297) * density(42) 
  pd(42,42) = pd(42,42) - rrt(297) * density(16) 
  pd(51,16) = pd(51,16) + rrt(297) * density(42) 
  pd(51,42) = pd(51,42) + rrt(297) * density(16) 
  pd(07,16) = pd(07,16) + rrt(298) * density(29) 
  pd(07,29) = pd(07,29) + rrt(298) * density(16) 
  pd(16,16) = pd(16,16) - rrt(298) * density(29) 
  pd(16,29) = pd(16,29) - rrt(298) * density(16) 
  pd(23,16) = pd(23,16) + rrt(298) * density(29) 
  pd(23,29) = pd(23,29) + rrt(298) * density(16) 
  pd(29,16) = pd(29,16) - rrt(298) * density(29) 
  pd(29,29) = pd(29,29) - rrt(298) * density(16) 
  pd(16,50) = pd(16,50) + rrt(299) * density(51) 
  pd(16,51) = pd(16,51) + rrt(299) * density(50) 
  pd(23,50) = pd(23,50) + rrt(299) * density(51) 
  pd(23,51) = pd(23,51) + rrt(299) * density(50) 
  pd(50,50) = pd(50,50) - rrt(299) * density(51) 
  pd(50,51) = pd(50,51) - rrt(299) * density(50) 
  pd(51,50) = pd(51,50) - rrt(299) * density(51) 
  pd(51,51) = pd(51,51) - rrt(299) * density(50) 
  pd(04,50) = pd(04,50) + rrt(300) * density(51) 
  pd(04,51) = pd(04,51) + rrt(300) * density(50) 
  pd(13,50) = pd(13,50) + rrt(300) * density(51) 
  pd(13,51) = pd(13,51) + rrt(300) * density(50) 
  pd(50,50) = pd(50,50) - rrt(300) * density(51) 
  pd(50,51) = pd(50,51) - rrt(300) * density(50) 
  pd(51,50) = pd(51,50) - rrt(300) * density(51) 
  pd(51,51) = pd(51,51) - rrt(300) * density(50) 
  pd(11,50) = pd(11,50) + rrt(301) * density(51) 
  pd(11,51) = pd(11,51) + rrt(301) * density(50) 
  pd(46,50) = pd(46,50) + rrt(301) * density(51) 
  pd(46,51) = pd(46,51) + rrt(301) * density(50) 
  pd(50,50) = pd(50,50) - rrt(301) * density(51) 
  pd(50,51) = pd(50,51) - rrt(301) * density(50) 
  pd(51,50) = pd(51,50) - rrt(301) * density(51) 
  pd(51,51) = pd(51,51) - rrt(301) * density(50) 
  pd(07,50) = pd(07,50) + rrt(302) * density(51) 
  pd(07,51) = pd(07,51) + rrt(302) * density(50) 
  pd(11,50) = pd(11,50) + rrt(302) * density(51) 
  pd(11,51) = pd(11,51) + rrt(302) * density(50) 
  pd(13,50) = pd(13,50) + rrt(302) * density(51) 
  pd(13,51) = pd(13,51) + rrt(302) * density(50) 
  pd(50,50) = pd(50,50) - rrt(302) * density(51) 
  pd(50,51) = pd(50,51) - rrt(302) * density(50) 
  pd(51,50) = pd(51,50) - rrt(302) * density(51) 
  pd(51,51) = pd(51,51) - rrt(302) * density(50) 
  pd(11,50) = pd(11,50) + rrt(303) * density(51) * 2.0d0
  pd(11,51) = pd(11,51) + rrt(303) * density(50) * 2.0d0
  pd(17,50) = pd(17,50) + rrt(303) * density(51) 
  pd(17,51) = pd(17,51) + rrt(303) * density(50) 
  pd(50,50) = pd(50,50) - rrt(303) * density(51) 
  pd(50,51) = pd(50,51) - rrt(303) * density(50) 
  pd(51,50) = pd(51,50) - rrt(303) * density(51) 
  pd(51,51) = pd(51,51) - rrt(303) * density(50) 
  pd(11,11) = pd(11,11) - rrt(304) * density(50) 
  pd(11,50) = pd(11,50) - rrt(304) * density(11) 
  pd(13,11) = pd(13,11) + rrt(304) * density(50) 
  pd(13,50) = pd(13,50) + rrt(304) * density(11) 
  pd(16,11) = pd(16,11) + rrt(304) * density(50) 
  pd(16,50) = pd(16,50) + rrt(304) * density(11) 
  pd(50,11) = pd(50,11) - rrt(304) * density(50) 
  pd(50,50) = pd(50,50) - rrt(304) * density(11) 
  pd(13,44) = pd(13,44) + rrt(305) * density(51) 
  pd(13,51) = pd(13,51) + rrt(305) * density(44) 
  pd(44,44) = pd(44,44) - rrt(305) * density(51) 
  pd(44,51) = pd(44,51) - rrt(305) * density(44) 
  pd(46,44) = pd(46,44) + rrt(305) * density(51) 
  pd(46,51) = pd(46,51) + rrt(305) * density(44) 
  pd(51,44) = pd(51,44) - rrt(305) * density(51) 
  pd(51,51) = pd(51,51) - rrt(305) * density(44) 
  pd(07,44) = pd(07,44) + rrt(306) * density(51) 
  pd(07,51) = pd(07,51) + rrt(306) * density(44) 
  pd(11,44) = pd(11,44) + rrt(306) * density(51) 
  pd(11,51) = pd(11,51) + rrt(306) * density(44) 
  pd(44,44) = pd(44,44) - rrt(306) * density(51) 
  pd(44,51) = pd(44,51) - rrt(306) * density(44) 
  pd(51,44) = pd(51,44) - rrt(306) * density(51) 
  pd(51,51) = pd(51,51) - rrt(306) * density(44) 
  pd(11,44) = pd(11,44) + rrt(307) * density(51) 
  pd(11,51) = pd(11,51) + rrt(307) * density(44) 
  pd(13,44) = pd(13,44) + rrt(307) * density(51) 
  pd(13,51) = pd(13,51) + rrt(307) * density(44) 
  pd(17,44) = pd(17,44) + rrt(307) * density(51) 
  pd(17,51) = pd(17,51) + rrt(307) * density(44) 
  pd(44,44) = pd(44,44) - rrt(307) * density(51) 
  pd(44,51) = pd(44,51) - rrt(307) * density(44) 
  pd(51,44) = pd(51,44) - rrt(307) * density(51) 
  pd(51,51) = pd(51,51) - rrt(307) * density(44) 
  pd(11,11) = pd(11,11) - rrt(308) * density(44) 
  pd(11,44) = pd(11,44) - rrt(308) * density(11) 
  pd(13,11) = pd(13,11) + rrt(308) * density(44) 
  pd(13,44) = pd(13,44) + rrt(308) * density(11) 
  pd(44,11) = pd(44,11) - rrt(308) * density(44) 
  pd(44,44) = pd(44,44) - rrt(308) * density(11) 
  pd(50,11) = pd(50,11) + rrt(308) * density(44) 
  pd(50,44) = pd(50,44) + rrt(308) * density(11) 
  pd(40,40) = pd(40,40) - rrt(309) * density(42) 
  pd(40,42) = pd(40,42) - rrt(309) * density(40) 
  pd(42,40) = pd(42,40) - rrt(309) * density(42) 
  pd(42,42) = pd(42,42) - rrt(309) * density(40) 
  pd(43,40) = pd(43,40) + rrt(309) * density(42) 
  pd(43,42) = pd(43,42) + rrt(309) * density(40) 
  pd(46,40) = pd(46,40) + rrt(309) * density(42) 
  pd(46,42) = pd(46,42) + rrt(309) * density(40) 
  pd(04,21) = pd(04,21) + rrt(310) * density(40) 
  pd(04,40) = pd(04,40) + rrt(310) * density(21) 
  pd(21,21) = pd(21,21) - rrt(310) * density(40) 
  pd(21,40) = pd(21,40) - rrt(310) * density(21) 
  pd(29,21) = pd(29,21) + rrt(310) * density(40) 
  pd(29,40) = pd(29,40) + rrt(310) * density(21) 
  pd(40,21) = pd(40,21) - rrt(310) * density(40) 
  pd(40,40) = pd(40,40) - rrt(310) * density(21) 
  pd(04,13) = pd(04,13) + rrt(311) * density(40) 
  pd(04,40) = pd(04,40) + rrt(311) * density(13) 
  pd(11,13) = pd(11,13) + rrt(311) * density(40) 
  pd(11,40) = pd(11,40) + rrt(311) * density(13) 
  pd(13,13) = pd(13,13) - rrt(311) * density(40) 
  pd(13,40) = pd(13,40) - rrt(311) * density(13) 
  pd(40,13) = pd(40,13) - rrt(311) * density(40) 
  pd(40,40) = pd(40,40) - rrt(311) * density(13) 
  pd(04,04) = pd(04,04) - rrt(312) * density(13) 
  pd(04,13) = pd(04,13) - rrt(312) * density(04) 
  pd(11,04) = pd(11,04) + rrt(312) * density(13) 
  pd(11,13) = pd(11,13) + rrt(312) * density(04) 
  pd(13,04) = pd(13,04) - rrt(312) * density(13) 
  pd(13,13) = pd(13,13) - rrt(312) * density(04) 
  pd(46,04) = pd(46,04) + rrt(312) * density(13) 
  pd(46,13) = pd(46,13) + rrt(312) * density(04) 
  pd(04,29) = pd(04,29) + rrt(313) * density(46) 
  pd(04,46) = pd(04,46) + rrt(313) * density(29) 
  pd(21,29) = pd(21,29) + rrt(313) * density(46) 
  pd(21,46) = pd(21,46) + rrt(313) * density(29) 
  pd(29,29) = pd(29,29) - rrt(313) * density(46) 
  pd(29,46) = pd(29,46) - rrt(313) * density(29) 
  pd(46,29) = pd(46,29) - rrt(313) * density(46) 
  pd(46,46) = pd(46,46) - rrt(313) * density(29) 
  pd(07,29) = pd(07,29) + rrt(314) * density(46) 
  pd(07,46) = pd(07,46) + rrt(314) * density(29) 
  pd(29,29) = pd(29,29) - rrt(314) * density(46) 
  pd(29,46) = pd(29,46) - rrt(314) * density(29) 
  pd(42,29) = pd(42,29) + rrt(314) * density(46) 
  pd(42,46) = pd(42,46) + rrt(314) * density(29) 
  pd(46,29) = pd(46,29) - rrt(314) * density(46) 
  pd(46,46) = pd(46,46) - rrt(314) * density(29) 
  pd(07,13) = pd(07,13) + rrt(315) * density(46) 
  pd(07,46) = pd(07,46) + rrt(315) * density(13) 
  pd(11,13) = pd(11,13) + rrt(315) * density(46) 
  pd(11,46) = pd(11,46) + rrt(315) * density(13) 
  pd(13,13) = pd(13,13) - rrt(315) * density(46) 
  pd(13,46) = pd(13,46) - rrt(315) * density(13) 
  pd(46,13) = pd(46,13) - rrt(315) * density(46) 
  pd(46,46) = pd(46,46) - rrt(315) * density(13) 
  pd(04,07) = pd(04,07) + rrt(316) * density(43) 
  pd(04,43) = pd(04,43) + rrt(316) * density(07) 
  pd(07,07) = pd(07,07) - rrt(316) * density(43) 
  pd(07,43) = pd(07,43) - rrt(316) * density(07) 
  pd(42,07) = pd(42,07) + rrt(316) * density(43) 
  pd(42,43) = pd(42,43) + rrt(316) * density(07) 
  pd(43,07) = pd(43,07) - rrt(316) * density(43) 
  pd(43,43) = pd(43,43) - rrt(316) * density(07) 
  pd(04,07) = pd(04,07) + rrt(317) * density(42) 
  pd(04,42) = pd(04,42) + rrt(317) * density(07) 
  pd(07,07) = pd(07,07) - rrt(317) * density(42) 
  pd(07,42) = pd(07,42) - rrt(317) * density(07) 
  pd(21,07) = pd(21,07) + rrt(317) * density(42) 
  pd(21,42) = pd(21,42) + rrt(317) * density(07) 
  pd(42,07) = pd(42,07) - rrt(317) * density(42) 
  pd(42,42) = pd(42,42) - rrt(317) * density(07) 
  pd(06,06) = pd(06,06) - rrt(318) * density(07) 
  pd(06,07) = pd(06,07) - rrt(318) * density(06) 
  pd(07,06) = pd(07,06) - rrt(318) * density(07) 
  pd(07,07) = pd(07,07) - rrt(318) * density(06) 
  pd(17,06) = pd(17,06) + rrt(318) * density(07) 
  pd(17,07) = pd(17,07) + rrt(318) * density(06) 
  pd(21,06) = pd(21,06) + rrt(318) * density(07) 
  pd(21,07) = pd(21,07) + rrt(318) * density(06) 
  pd(07,07) = pd(07,07) - rrt(319) * density(13) 
  pd(07,13) = pd(07,13) - rrt(319) * density(07) 
  pd(11,07) = pd(11,07) + rrt(319) * density(13) 
  pd(11,13) = pd(11,13) + rrt(319) * density(07) 
  pd(13,07) = pd(13,07) - rrt(319) * density(13) 
  pd(13,13) = pd(13,13) - rrt(319) * density(07) 
  pd(17,07) = pd(17,07) + rrt(319) * density(13) 
  pd(17,13) = pd(17,13) + rrt(319) * density(07) 
  pd(07,17) = pd(07,17) + rrt(320) * density(51) 
  pd(07,51) = pd(07,51) + rrt(320) * density(17) 
  pd(17,17) = pd(17,17) - rrt(320) * density(51) 
  pd(17,51) = pd(17,51) - rrt(320) * density(17) 
  pd(23,17) = pd(23,17) + rrt(320) * density(51) 
  pd(23,51) = pd(23,51) + rrt(320) * density(17) 
  pd(51,17) = pd(51,17) - rrt(320) * density(51) 
  pd(51,51) = pd(51,51) - rrt(320) * density(17) 
  pd(04,17) = pd(04,17) + rrt(321) * density(43) 
  pd(04,43) = pd(04,43) + rrt(321) * density(17) 
  pd(17,17) = pd(17,17) - rrt(321) * density(43) 
  pd(17,43) = pd(17,43) - rrt(321) * density(17) 
  pd(29,17) = pd(29,17) + rrt(321) * density(43) 
  pd(29,43) = pd(29,43) + rrt(321) * density(17) 
  pd(43,17) = pd(43,17) - rrt(321) * density(43) 
  pd(43,43) = pd(43,43) - rrt(321) * density(17) 
  pd(17,17) = pd(17,17) - rrt(322) * density(43) 
  pd(17,43) = pd(17,43) - rrt(322) * density(17) 
  pd(42,17) = pd(42,17) + rrt(322) * density(43) 
  pd(42,43) = pd(42,43) + rrt(322) * density(17) 
  pd(43,17) = pd(43,17) - rrt(322) * density(43) 
  pd(43,43) = pd(43,43) - rrt(322) * density(17) 
  pd(46,17) = pd(46,17) + rrt(322) * density(43) 
  pd(46,43) = pd(46,43) + rrt(322) * density(17) 
  pd(17,17) = pd(17,17) - rrt(323) * density(42) 
  pd(17,42) = pd(17,42) - rrt(323) * density(17) 
  pd(21,17) = pd(21,17) + rrt(323) * density(42) 
  pd(21,42) = pd(21,42) + rrt(323) * density(17) 
  pd(42,17) = pd(42,17) - rrt(323) * density(42) 
  pd(42,42) = pd(42,42) - rrt(323) * density(17) 
  pd(46,17) = pd(46,17) + rrt(323) * density(42) 
  pd(46,42) = pd(46,42) + rrt(323) * density(17) 
  pd(07,17) = pd(07,17) + rrt(324) * density(29) 
  pd(07,29) = pd(07,29) + rrt(324) * density(17) 
  pd(17,17) = pd(17,17) - rrt(324) * density(29) 
  pd(17,29) = pd(17,29) - rrt(324) * density(17) 
  pd(21,17) = pd(21,17) + rrt(324) * density(29) 
  pd(21,29) = pd(21,29) + rrt(324) * density(17) 
  pd(29,17) = pd(29,17) - rrt(324) * density(29) 
  pd(29,29) = pd(29,29) - rrt(324) * density(17) 
  pd(07,11) = pd(07,11) + rrt(325) * density(17) 
  pd(07,17) = pd(07,17) + rrt(325) * density(11) 
  pd(11,11) = pd(11,11) - rrt(325) * density(17) 
  pd(11,17) = pd(11,17) - rrt(325) * density(11) 
  pd(13,11) = pd(13,11) + rrt(325) * density(17) 
  pd(13,17) = pd(13,17) + rrt(325) * density(11) 
  pd(17,11) = pd(17,11) - rrt(325) * density(17) 
  pd(17,17) = pd(17,17) - rrt(325) * density(11) 
  pd(17,36) = pd(17,36) + rrt(326) * density(51) 
  pd(17,51) = pd(17,51) + rrt(326) * density(36) 
  pd(23,36) = pd(23,36) + rrt(326) * density(51) 
  pd(23,51) = pd(23,51) + rrt(326) * density(36) 
  pd(36,36) = pd(36,36) - rrt(326) * density(51) 
  pd(36,51) = pd(36,51) - rrt(326) * density(36) 
  pd(51,36) = pd(51,36) - rrt(326) * density(51) 
  pd(51,51) = pd(51,51) - rrt(326) * density(36) 
  pd(11,12) = pd(11,12) + rrt(327) * density(51) 
  pd(11,51) = pd(11,51) + rrt(327) * density(12) 
  pd(12,12) = pd(12,12) - rrt(327) * density(51) 
  pd(12,51) = pd(12,51) - rrt(327) * density(12) 
  pd(47,12) = pd(47,12) + rrt(327) * density(51) 
  pd(47,51) = pd(47,51) + rrt(327) * density(12) 
  pd(51,12) = pd(51,12) - rrt(327) * density(51) 
  pd(51,51) = pd(51,51) - rrt(327) * density(12) 
  pd(10,12) = pd(10,12) + rrt(328) * density(23) 
  pd(10,23) = pd(10,23) + rrt(328) * density(12) 
  pd(11,12) = pd(11,12) + rrt(328) * density(23) 
  pd(11,23) = pd(11,23) + rrt(328) * density(12) 
  pd(12,12) = pd(12,12) - rrt(328) * density(23) 
  pd(12,23) = pd(12,23) - rrt(328) * density(12) 
  pd(23,12) = pd(23,12) - rrt(328) * density(23) 
  pd(23,23) = pd(23,23) - rrt(328) * density(12) 
  pd(11,12) = pd(11,12) + rrt(329) * density(28) 
  pd(11,28) = pd(11,28) + rrt(329) * density(12) 
  pd(12,12) = pd(12,12) - rrt(329) * density(28) 
  pd(12,28) = pd(12,28) - rrt(329) * density(12) 
  pd(16,12) = pd(16,12) + rrt(329) * density(28) 
  pd(16,28) = pd(16,28) + rrt(329) * density(12) 
  pd(28,12) = pd(28,12) - rrt(329) * density(28) 
  pd(28,28) = pd(28,28) - rrt(329) * density(12) 
  pd(08,08) = pd(08,08) - rrt(330) * density(12) 
  pd(08,12) = pd(08,12) - rrt(330) * density(08) 
  pd(11,08) = pd(11,08) + rrt(330) * density(12) 
  pd(11,12) = pd(11,12) + rrt(330) * density(08) 
  pd(12,08) = pd(12,08) - rrt(330) * density(12) 
  pd(12,12) = pd(12,12) - rrt(330) * density(08) 
  pd(50,08) = pd(50,08) + rrt(330) * density(12) 
  pd(50,12) = pd(50,12) + rrt(330) * density(08) 
  pd(04,12) = pd(04,12) + rrt(331) * density(43) 
  pd(04,43) = pd(04,43) + rrt(331) * density(12) 
  pd(11,12) = pd(11,12) + rrt(331) * density(43) * 2.0d0
  pd(11,43) = pd(11,43) + rrt(331) * density(12) * 2.0d0
  pd(12,12) = pd(12,12) - rrt(331) * density(43) 
  pd(12,43) = pd(12,43) - rrt(331) * density(12) 
  pd(43,12) = pd(43,12) - rrt(331) * density(43) 
  pd(43,43) = pd(43,43) - rrt(331) * density(12) 
  pd(11,12) = pd(11,12) + rrt(332) * density(34) 
  pd(11,34) = pd(11,34) + rrt(332) * density(12) 
  pd(12,12) = pd(12,12) - rrt(332) * density(34) 
  pd(12,34) = pd(12,34) - rrt(332) * density(12) 
  pd(34,12) = pd(34,12) - rrt(332) * density(34) 
  pd(34,34) = pd(34,34) - rrt(332) * density(12) 
  pd(40,12) = pd(40,12) + rrt(332) * density(34) 
  pd(40,34) = pd(40,34) + rrt(332) * density(12) 
  pd(04,12) = pd(04,12) + rrt(333) * density(42) 
  pd(04,42) = pd(04,42) + rrt(333) * density(12) 
  pd(11,12) = pd(11,12) + rrt(333) * density(42) 
  pd(11,42) = pd(11,42) + rrt(333) * density(12) 
  pd(12,12) = pd(12,12) - rrt(333) * density(42) 
  pd(12,42) = pd(12,42) - rrt(333) * density(12) 
  pd(42,12) = pd(42,12) - rrt(333) * density(42) 
  pd(42,42) = pd(42,42) - rrt(333) * density(12) 
  pd(07,12) = pd(07,12) + rrt(334) * density(42) 
  pd(07,42) = pd(07,42) + rrt(334) * density(12) 
  pd(11,12) = pd(11,12) + rrt(334) * density(42) * 2.0d0
  pd(11,42) = pd(11,42) + rrt(334) * density(12) * 2.0d0
  pd(12,12) = pd(12,12) - rrt(334) * density(42) 
  pd(12,42) = pd(12,42) - rrt(334) * density(12) 
  pd(42,12) = pd(42,12) - rrt(334) * density(42) 
  pd(42,42) = pd(42,42) - rrt(334) * density(12) 
  pd(11,12) = pd(11,12) + rrt(335) * density(29) 
  pd(11,29) = pd(11,29) + rrt(335) * density(12) 
  pd(12,12) = pd(12,12) - rrt(335) * density(29) 
  pd(12,29) = pd(12,29) - rrt(335) * density(12) 
  pd(29,12) = pd(29,12) - rrt(335) * density(29) 
  pd(29,29) = pd(29,29) - rrt(335) * density(12) 
  pd(46,12) = pd(46,12) + rrt(335) * density(29) 
  pd(46,29) = pd(46,29) + rrt(335) * density(12) 
  pd(06,06) = pd(06,06) - rrt(336) * density(12) 
  pd(06,12) = pd(06,12) - rrt(336) * density(06) 
  pd(11,06) = pd(11,06) + rrt(336) * density(12) 
  pd(11,12) = pd(11,12) + rrt(336) * density(06) 
  pd(12,06) = pd(12,06) - rrt(336) * density(12) 
  pd(12,12) = pd(12,12) - rrt(336) * density(06) 
  pd(17,06) = pd(17,06) + rrt(336) * density(12) 
  pd(17,12) = pd(17,12) + rrt(336) * density(06) 
  pd(07,12) = pd(07,12) + rrt(337) * density(21) 
  pd(07,21) = pd(07,21) + rrt(337) * density(12) 
  pd(11,12) = pd(11,12) + rrt(337) * density(21) 
  pd(11,21) = pd(11,21) + rrt(337) * density(12) 
  pd(12,12) = pd(12,12) - rrt(337) * density(21) 
  pd(12,21) = pd(12,21) - rrt(337) * density(12) 
  pd(21,12) = pd(21,12) - rrt(337) * density(21) 
  pd(21,21) = pd(21,21) - rrt(337) * density(12) 
  pd(13,27) = pd(13,27) + rrt(338) * density(51) 
  pd(13,51) = pd(13,51) + rrt(338) * density(27) 
  pd(27,27) = pd(27,27) - rrt(338) * density(51) 
  pd(27,51) = pd(27,51) - rrt(338) * density(27) 
  pd(47,27) = pd(47,27) + rrt(338) * density(51) 
  pd(47,51) = pd(47,51) + rrt(338) * density(27) 
  pd(51,27) = pd(51,27) - rrt(338) * density(51) 
  pd(51,51) = pd(51,51) - rrt(338) * density(27) 
  pd(10,27) = pd(10,27) + rrt(339) * density(51) 
  pd(10,51) = pd(10,51) + rrt(339) * density(27) 
  pd(11,27) = pd(11,27) + rrt(339) * density(51) 
  pd(11,51) = pd(11,51) + rrt(339) * density(27) 
  pd(27,27) = pd(27,27) - rrt(339) * density(51) 
  pd(27,51) = pd(27,51) - rrt(339) * density(27) 
  pd(51,27) = pd(51,27) - rrt(339) * density(51) 
  pd(51,51) = pd(51,51) - rrt(339) * density(27) 
  pd(11,27) = pd(11,27) + rrt(340) * density(51) 
  pd(11,51) = pd(11,51) + rrt(340) * density(27) 
  pd(13,27) = pd(13,27) + rrt(340) * density(51) 
  pd(13,51) = pd(13,51) + rrt(340) * density(27) 
  pd(16,27) = pd(16,27) + rrt(340) * density(51) 
  pd(16,51) = pd(16,51) + rrt(340) * density(27) 
  pd(27,27) = pd(27,27) - rrt(340) * density(51) 
  pd(27,51) = pd(27,51) - rrt(340) * density(27) 
  pd(51,27) = pd(51,27) - rrt(340) * density(51) 
  pd(51,51) = pd(51,51) - rrt(340) * density(27) 
  pd(13,27) = pd(13,27) + rrt(341) * density(28) 
  pd(13,28) = pd(13,28) + rrt(341) * density(27) 
  pd(16,27) = pd(16,27) + rrt(341) * density(28) 
  pd(16,28) = pd(16,28) + rrt(341) * density(27) 
  pd(27,27) = pd(27,27) - rrt(341) * density(28) 
  pd(27,28) = pd(27,28) - rrt(341) * density(27) 
  pd(28,27) = pd(28,27) - rrt(341) * density(28) 
  pd(28,28) = pd(28,28) - rrt(341) * density(27) 
  pd(11,27) = pd(11,27) + rrt(342) * density(28) 
  pd(11,28) = pd(11,28) + rrt(342) * density(27) 
  pd(27,27) = pd(27,27) - rrt(342) * density(28) 
  pd(27,28) = pd(27,28) - rrt(342) * density(27) 
  pd(28,27) = pd(28,27) - rrt(342) * density(28) 
  pd(28,28) = pd(28,28) - rrt(342) * density(27) 
  pd(50,27) = pd(50,27) + rrt(342) * density(28) 
  pd(50,28) = pd(50,28) + rrt(342) * density(27) 
  pd(08,08) = pd(08,08) - rrt(343) * density(27) 
  pd(08,27) = pd(08,27) - rrt(343) * density(08) 
  pd(13,08) = pd(13,08) + rrt(343) * density(27) 
  pd(13,27) = pd(13,27) + rrt(343) * density(08) 
  pd(27,08) = pd(27,08) - rrt(343) * density(27) 
  pd(27,27) = pd(27,27) - rrt(343) * density(08) 
  pd(50,08) = pd(50,08) + rrt(343) * density(27) 
  pd(50,27) = pd(50,27) + rrt(343) * density(08) 
  pd(08,08) = pd(08,08) - rrt(344) * density(27) 
  pd(08,27) = pd(08,27) - rrt(344) * density(08) 
  pd(11,08) = pd(11,08) + rrt(344) * density(27) 
  pd(11,27) = pd(11,27) + rrt(344) * density(08) 
  pd(27,08) = pd(27,08) - rrt(344) * density(27) 
  pd(27,27) = pd(27,27) - rrt(344) * density(08) 
  pd(44,08) = pd(44,08) + rrt(344) * density(27) 
  pd(44,27) = pd(44,27) + rrt(344) * density(08) 
  pd(11,27) = pd(11,27) + rrt(345) * density(43) 
  pd(11,43) = pd(11,43) + rrt(345) * density(27) 
  pd(27,27) = pd(27,27) - rrt(345) * density(43) 
  pd(27,43) = pd(27,43) - rrt(345) * density(27) 
  pd(40,27) = pd(40,27) + rrt(345) * density(43) 
  pd(40,43) = pd(40,43) + rrt(345) * density(27) 
  pd(43,27) = pd(43,27) - rrt(345) * density(43) 
  pd(43,43) = pd(43,43) - rrt(345) * density(27) 
  pd(04,27) = pd(04,27) + rrt(346) * density(43) 
  pd(04,43) = pd(04,43) + rrt(346) * density(27) 
  pd(11,27) = pd(11,27) + rrt(346) * density(43) 
  pd(11,43) = pd(11,43) + rrt(346) * density(27) 
  pd(13,27) = pd(13,27) + rrt(346) * density(43) 
  pd(13,43) = pd(13,43) + rrt(346) * density(27) 
  pd(27,27) = pd(27,27) - rrt(346) * density(43) 
  pd(27,43) = pd(27,43) - rrt(346) * density(27) 
  pd(43,27) = pd(43,27) - rrt(346) * density(43) 
  pd(43,43) = pd(43,43) - rrt(346) * density(27) 
  pd(11,27) = pd(11,27) + rrt(347) * density(43) * 2.0d0
  pd(11,43) = pd(11,43) + rrt(347) * density(27) * 2.0d0
  pd(27,27) = pd(27,27) - rrt(347) * density(43) 
  pd(27,43) = pd(27,43) - rrt(347) * density(27) 
  pd(43,27) = pd(43,27) - rrt(347) * density(43) 
  pd(43,43) = pd(43,43) - rrt(347) * density(27) 
  pd(46,27) = pd(46,27) + rrt(347) * density(43) 
  pd(46,43) = pd(46,43) + rrt(347) * density(27) 
  pd(07,27) = pd(07,27) + rrt(348) * density(43) 
  pd(07,43) = pd(07,43) + rrt(348) * density(27) 
  pd(11,27) = pd(11,27) + rrt(348) * density(43) * 2.0d0
  pd(11,43) = pd(11,43) + rrt(348) * density(27) * 2.0d0
  pd(13,27) = pd(13,27) + rrt(348) * density(43) 
  pd(13,43) = pd(13,43) + rrt(348) * density(27) 
  pd(27,27) = pd(27,27) - rrt(348) * density(43) 
  pd(27,43) = pd(27,43) - rrt(348) * density(27) 
  pd(43,27) = pd(43,27) - rrt(348) * density(43) 
  pd(43,43) = pd(43,43) - rrt(348) * density(27) 
  pd(11,27) = pd(11,27) + rrt(349) * density(43) * 3.0d0
  pd(11,43) = pd(11,43) + rrt(349) * density(27) * 3.0d0
  pd(17,27) = pd(17,27) + rrt(349) * density(43) 
  pd(17,43) = pd(17,43) + rrt(349) * density(27) 
  pd(27,27) = pd(27,27) - rrt(349) * density(43) 
  pd(27,43) = pd(27,43) - rrt(349) * density(27) 
  pd(43,27) = pd(43,27) - rrt(349) * density(43) 
  pd(43,43) = pd(43,43) - rrt(349) * density(27) 
  pd(11,27) = pd(11,27) + rrt(350) * density(42) 
  pd(11,42) = pd(11,42) + rrt(350) * density(27) 
  pd(27,27) = pd(27,27) - rrt(350) * density(42) 
  pd(27,42) = pd(27,42) - rrt(350) * density(27) 
  pd(42,27) = pd(42,27) - rrt(350) * density(42) 
  pd(42,42) = pd(42,42) - rrt(350) * density(27) 
  pd(46,27) = pd(46,27) + rrt(350) * density(42) 
  pd(46,42) = pd(46,42) + rrt(350) * density(27) 
  pd(07,27) = pd(07,27) + rrt(351) * density(42) 
  pd(07,42) = pd(07,42) + rrt(351) * density(27) 
  pd(11,27) = pd(11,27) + rrt(351) * density(42) 
  pd(11,42) = pd(11,42) + rrt(351) * density(27) 
  pd(13,27) = pd(13,27) + rrt(351) * density(42) 
  pd(13,42) = pd(13,42) + rrt(351) * density(27) 
  pd(27,27) = pd(27,27) - rrt(351) * density(42) 
  pd(27,42) = pd(27,42) - rrt(351) * density(27) 
  pd(42,27) = pd(42,27) - rrt(351) * density(42) 
  pd(42,42) = pd(42,42) - rrt(351) * density(27) 
  pd(11,27) = pd(11,27) + rrt(352) * density(42) * 2.0d0
  pd(11,42) = pd(11,42) + rrt(352) * density(27) * 2.0d0
  pd(17,27) = pd(17,27) + rrt(352) * density(42) 
  pd(17,42) = pd(17,42) + rrt(352) * density(27) 
  pd(27,27) = pd(27,27) - rrt(352) * density(42) 
  pd(27,42) = pd(27,42) - rrt(352) * density(27) 
  pd(42,27) = pd(42,27) - rrt(352) * density(42) 
  pd(42,42) = pd(42,42) - rrt(352) * density(27) 
  pd(07,21) = pd(07,21) + rrt(353) * density(27) 
  pd(07,27) = pd(07,27) + rrt(353) * density(21) 
  pd(13,21) = pd(13,21) + rrt(353) * density(27) 
  pd(13,27) = pd(13,27) + rrt(353) * density(21) 
  pd(21,21) = pd(21,21) - rrt(353) * density(27) 
  pd(21,27) = pd(21,27) - rrt(353) * density(21) 
  pd(27,21) = pd(27,21) - rrt(353) * density(27) 
  pd(27,27) = pd(27,27) - rrt(353) * density(21) 
  pd(11,21) = pd(11,21) + rrt(354) * density(27) 
  pd(11,27) = pd(11,27) + rrt(354) * density(21) 
  pd(17,21) = pd(17,21) + rrt(354) * density(27) 
  pd(17,27) = pd(17,27) + rrt(354) * density(21) 
  pd(21,21) = pd(21,21) - rrt(354) * density(27) 
  pd(21,27) = pd(21,27) - rrt(354) * density(21) 
  pd(27,21) = pd(27,21) - rrt(354) * density(27) 
  pd(27,27) = pd(27,27) - rrt(354) * density(21) 
  pd(12,13) = pd(12,13) + rrt(355) * density(27) 
  pd(12,27) = pd(12,27) + rrt(355) * density(13) 
  pd(13,13) = pd(13,13) - rrt(355) * density(27) 
  pd(13,27) = pd(13,27) - rrt(355) * density(13) 
  pd(27,13) = pd(27,13) - rrt(355) * density(27) 
  pd(27,27) = pd(27,27) - rrt(355) * density(13) 
  pd(11,13) = pd(11,13) + rrt(356) * density(27) 
  pd(11,27) = pd(11,27) + rrt(356) * density(13) 
  pd(13,13) = pd(13,13) - rrt(356) * density(27) 
  pd(13,27) = pd(13,27) - rrt(356) * density(13) 
  pd(27,13) = pd(27,13) - rrt(356) * density(27) 
  pd(27,27) = pd(27,27) - rrt(356) * density(13) 
  pd(30,13) = pd(30,13) + rrt(356) * density(27) 
  pd(30,27) = pd(30,27) + rrt(356) * density(13) 
  pd(10,30) = pd(10,30) + rrt(357) * density(51) 
  pd(10,51) = pd(10,51) + rrt(357) * density(30) 
  pd(13,30) = pd(13,30) + rrt(357) * density(51) 
  pd(13,51) = pd(13,51) + rrt(357) * density(30) 
  pd(30,30) = pd(30,30) - rrt(357) * density(51) 
  pd(30,51) = pd(30,51) - rrt(357) * density(30) 
  pd(51,30) = pd(51,30) - rrt(357) * density(51) 
  pd(51,51) = pd(51,51) - rrt(357) * density(30) 
  pd(11,30) = pd(11,30) + rrt(358) * density(51) 
  pd(11,51) = pd(11,51) + rrt(358) * density(30) 
  pd(16,30) = pd(16,30) + rrt(358) * density(51) 
  pd(16,51) = pd(16,51) + rrt(358) * density(30) 
  pd(30,30) = pd(30,30) - rrt(358) * density(51) 
  pd(30,51) = pd(30,51) - rrt(358) * density(30) 
  pd(51,30) = pd(51,30) - rrt(358) * density(51) 
  pd(51,51) = pd(51,51) - rrt(358) * density(30) 
  pd(13,23) = pd(13,23) + rrt(359) * density(30) 
  pd(13,30) = pd(13,30) + rrt(359) * density(23) 
  pd(16,23) = pd(16,23) + rrt(359) * density(30) 
  pd(16,30) = pd(16,30) + rrt(359) * density(23) 
  pd(23,23) = pd(23,23) - rrt(359) * density(30) 
  pd(23,30) = pd(23,30) - rrt(359) * density(23) 
  pd(30,23) = pd(30,23) - rrt(359) * density(30) 
  pd(30,30) = pd(30,30) - rrt(359) * density(23) 
  pd(13,28) = pd(13,28) + rrt(360) * density(30) 
  pd(13,30) = pd(13,30) + rrt(360) * density(28) 
  pd(28,28) = pd(28,28) - rrt(360) * density(30) 
  pd(28,30) = pd(28,30) - rrt(360) * density(28) 
  pd(30,28) = pd(30,28) - rrt(360) * density(30) 
  pd(30,30) = pd(30,30) - rrt(360) * density(28) 
  pd(50,28) = pd(50,28) + rrt(360) * density(30) 
  pd(50,30) = pd(50,30) + rrt(360) * density(28) 
  pd(11,28) = pd(11,28) + rrt(361) * density(30) 
  pd(11,30) = pd(11,30) + rrt(361) * density(28) 
  pd(28,28) = pd(28,28) - rrt(361) * density(30) 
  pd(28,30) = pd(28,30) - rrt(361) * density(28) 
  pd(30,28) = pd(30,28) - rrt(361) * density(30) 
  pd(30,30) = pd(30,30) - rrt(361) * density(28) 
  pd(44,28) = pd(44,28) + rrt(361) * density(30) 
  pd(44,30) = pd(44,30) + rrt(361) * density(28) 
  pd(08,08) = pd(08,08) - rrt(362) * density(30) 
  pd(08,30) = pd(08,30) - rrt(362) * density(08) 
  pd(13,08) = pd(13,08) + rrt(362) * density(30) 
  pd(13,30) = pd(13,30) + rrt(362) * density(08) 
  pd(30,08) = pd(30,08) - rrt(362) * density(30) 
  pd(30,30) = pd(30,30) - rrt(362) * density(08) 
  pd(44,08) = pd(44,08) + rrt(362) * density(30) 
  pd(44,30) = pd(44,30) + rrt(362) * density(08) 
  pd(04,30) = pd(04,30) + rrt(363) * density(43) 
  pd(04,43) = pd(04,43) + rrt(363) * density(30) 
  pd(11,30) = pd(11,30) + rrt(363) * density(43) 
  pd(11,43) = pd(11,43) + rrt(363) * density(30) 
  pd(30,30) = pd(30,30) - rrt(363) * density(43) 
  pd(30,43) = pd(30,43) - rrt(363) * density(30) 
  pd(43,30) = pd(43,30) - rrt(363) * density(43) 
  pd(43,43) = pd(43,43) - rrt(363) * density(30) 
  pd(11,30) = pd(11,30) + rrt(364) * density(43) 
  pd(11,43) = pd(11,43) + rrt(364) * density(30) 
  pd(13,30) = pd(13,30) + rrt(364) * density(43) 
  pd(13,43) = pd(13,43) + rrt(364) * density(30) 
  pd(30,30) = pd(30,30) - rrt(364) * density(43) 
  pd(30,43) = pd(30,43) - rrt(364) * density(30) 
  pd(43,30) = pd(43,30) - rrt(364) * density(43) 
  pd(43,43) = pd(43,43) - rrt(364) * density(30) 
  pd(46,30) = pd(46,30) + rrt(364) * density(43) 
  pd(46,43) = pd(46,43) + rrt(364) * density(30) 
  pd(07,30) = pd(07,30) + rrt(365) * density(43) 
  pd(07,43) = pd(07,43) + rrt(365) * density(30) 
  pd(11,30) = pd(11,30) + rrt(365) * density(43) * 2.0d0
  pd(11,43) = pd(11,43) + rrt(365) * density(30) * 2.0d0
  pd(30,30) = pd(30,30) - rrt(365) * density(43) 
  pd(30,43) = pd(30,43) - rrt(365) * density(30) 
  pd(43,30) = pd(43,30) - rrt(365) * density(43) 
  pd(43,43) = pd(43,43) - rrt(365) * density(30) 
  pd(11,30) = pd(11,30) + rrt(366) * density(34) 
  pd(11,34) = pd(11,34) + rrt(366) * density(30) 
  pd(30,30) = pd(30,30) - rrt(366) * density(34) 
  pd(30,34) = pd(30,34) - rrt(366) * density(30) 
  pd(34,30) = pd(34,30) - rrt(366) * density(34) 
  pd(34,34) = pd(34,34) - rrt(366) * density(30) 
  pd(46,30) = pd(46,30) + rrt(366) * density(34) 
  pd(46,34) = pd(46,34) + rrt(366) * density(30) 
  pd(07,30) = pd(07,30) + rrt(367) * density(34) 
  pd(07,34) = pd(07,34) + rrt(367) * density(30) 
  pd(11,30) = pd(11,30) + rrt(367) * density(34) 
  pd(11,34) = pd(11,34) + rrt(367) * density(30) 
  pd(13,30) = pd(13,30) + rrt(367) * density(34) 
  pd(13,34) = pd(13,34) + rrt(367) * density(30) 
  pd(30,30) = pd(30,30) - rrt(367) * density(34) 
  pd(30,34) = pd(30,34) - rrt(367) * density(30) 
  pd(34,30) = pd(34,30) - rrt(367) * density(34) 
  pd(34,34) = pd(34,34) - rrt(367) * density(30) 
  pd(13,30) = pd(13,30) + rrt(368) * density(42) 
  pd(13,42) = pd(13,42) + rrt(368) * density(30) 
  pd(30,30) = pd(30,30) - rrt(368) * density(42) 
  pd(30,42) = pd(30,42) - rrt(368) * density(30) 
  pd(42,30) = pd(42,30) - rrt(368) * density(42) 
  pd(42,42) = pd(42,42) - rrt(368) * density(30) 
  pd(46,30) = pd(46,30) + rrt(368) * density(42) 
  pd(46,42) = pd(46,42) + rrt(368) * density(30) 
  pd(07,30) = pd(07,30) + rrt(369) * density(42) 
  pd(07,42) = pd(07,42) + rrt(369) * density(30) 
  pd(11,30) = pd(11,30) + rrt(369) * density(42) 
  pd(11,42) = pd(11,42) + rrt(369) * density(30) 
  pd(30,30) = pd(30,30) - rrt(369) * density(42) 
  pd(30,42) = pd(30,42) - rrt(369) * density(30) 
  pd(42,30) = pd(42,30) - rrt(369) * density(42) 
  pd(42,42) = pd(42,42) - rrt(369) * density(30) 
  pd(11,30) = pd(11,30) + rrt(370) * density(42) 
  pd(11,42) = pd(11,42) + rrt(370) * density(30) 
  pd(13,30) = pd(13,30) + rrt(370) * density(42) 
  pd(13,42) = pd(13,42) + rrt(370) * density(30) 
  pd(17,30) = pd(17,30) + rrt(370) * density(42) 
  pd(17,42) = pd(17,42) + rrt(370) * density(30) 
  pd(30,30) = pd(30,30) - rrt(370) * density(42) 
  pd(30,42) = pd(30,42) - rrt(370) * density(30) 
  pd(42,30) = pd(42,30) - rrt(370) * density(42) 
  pd(42,42) = pd(42,42) - rrt(370) * density(30) 
  pd(07,29) = pd(07,29) + rrt(371) * density(30) 
  pd(07,30) = pd(07,30) + rrt(371) * density(29) 
  pd(13,29) = pd(13,29) + rrt(371) * density(30) 
  pd(13,30) = pd(13,30) + rrt(371) * density(29) 
  pd(29,29) = pd(29,29) - rrt(371) * density(30) 
  pd(29,30) = pd(29,30) - rrt(371) * density(29) 
  pd(30,29) = pd(30,29) - rrt(371) * density(30) 
  pd(30,30) = pd(30,30) - rrt(371) * density(29) 
  pd(11,29) = pd(11,29) + rrt(372) * density(30) 
  pd(11,30) = pd(11,30) + rrt(372) * density(29) 
  pd(17,29) = pd(17,29) + rrt(372) * density(30) 
  pd(17,30) = pd(17,30) + rrt(372) * density(29) 
  pd(29,29) = pd(29,29) - rrt(372) * density(30) 
  pd(29,30) = pd(29,30) - rrt(372) * density(29) 
  pd(30,29) = pd(30,29) - rrt(372) * density(30) 
  pd(30,30) = pd(30,30) - rrt(372) * density(29) 
  pd(13,21) = pd(13,21) + rrt(373) * density(30) 
  pd(13,30) = pd(13,30) + rrt(373) * density(21) 
  pd(17,21) = pd(17,21) + rrt(373) * density(30) 
  pd(17,30) = pd(17,30) + rrt(373) * density(21) 
  pd(21,21) = pd(21,21) - rrt(373) * density(30) 
  pd(21,30) = pd(21,30) - rrt(373) * density(21) 
  pd(30,21) = pd(30,21) - rrt(373) * density(30) 
  pd(30,30) = pd(30,30) - rrt(373) * density(21) 
  pd(23,28) = pd(23,28) + rrt(374) * density(51) * 2.0d0
  pd(23,51) = pd(23,51) + rrt(374) * density(28) * 2.0d0
  pd(28,28) = pd(28,28) - rrt(374) * density(51) 
  pd(28,51) = pd(28,51) - rrt(374) * density(28) 
  pd(51,28) = pd(51,28) - rrt(374) * density(51) 
  pd(51,51) = pd(51,51) - rrt(374) * density(28) 
  pd(08,08) = pd(08,08) - rrt(375) * density(51) 
  pd(08,51) = pd(08,51) - rrt(375) * density(08) 
  pd(13,08) = pd(13,08) + rrt(375) * density(51) 
  pd(13,51) = pd(13,51) + rrt(375) * density(08) 
  pd(42,08) = pd(42,08) + rrt(375) * density(51) 
  pd(42,51) = pd(42,51) + rrt(375) * density(08) 
  pd(51,08) = pd(51,08) - rrt(375) * density(51) 
  pd(51,51) = pd(51,51) - rrt(375) * density(08) 
  pd(23,34) = pd(23,34) + rrt(376) * density(51) 
  pd(23,51) = pd(23,51) + rrt(376) * density(34) 
  pd(34,34) = pd(34,34) - rrt(376) * density(51) 
  pd(34,51) = pd(34,51) - rrt(376) * density(34) 
  pd(43,34) = pd(43,34) + rrt(376) * density(51) 
  pd(43,51) = pd(43,51) + rrt(376) * density(34) 
  pd(51,34) = pd(51,34) - rrt(376) * density(51) 
  pd(51,51) = pd(51,51) - rrt(376) * density(34) 
  pd(23,29) = pd(23,29) + rrt(377) * density(51) 
  pd(23,51) = pd(23,51) + rrt(377) * density(29) 
  pd(29,29) = pd(29,29) - rrt(377) * density(51) 
  pd(29,51) = pd(29,51) - rrt(377) * density(29) 
  pd(42,29) = pd(42,29) + rrt(377) * density(51) 
  pd(42,51) = pd(42,51) + rrt(377) * density(29) 
  pd(51,29) = pd(51,29) - rrt(377) * density(51) 
  pd(51,51) = pd(51,51) - rrt(377) * density(29) 
  pd(06,06) = pd(06,06) - rrt(378) * density(51) 
  pd(06,51) = pd(06,51) - rrt(378) * density(06) 
  pd(21,06) = pd(21,06) + rrt(378) * density(51) 
  pd(21,51) = pd(21,51) + rrt(378) * density(06) 
  pd(23,06) = pd(23,06) + rrt(378) * density(51) 
  pd(23,51) = pd(23,51) + rrt(378) * density(06) 
  pd(51,06) = pd(51,06) - rrt(378) * density(51) 
  pd(51,51) = pd(51,51) - rrt(378) * density(06) 
  pd(18,18) = pd(18,18) - rrt(379) * density(51) 
  pd(18,51) = pd(18,51) - rrt(379) * density(18) 
  pd(23,18) = pd(23,18) + rrt(379) * density(51) 
  pd(23,51) = pd(23,51) + rrt(379) * density(18) 
  pd(45,18) = pd(45,18) + rrt(379) * density(51) 
  pd(45,51) = pd(45,51) + rrt(379) * density(18) 
  pd(51,18) = pd(51,18) - rrt(379) * density(51) 
  pd(51,51) = pd(51,51) - rrt(379) * density(18) 
  pd(20,20) = pd(20,20) - rrt(380) * density(51) 
  pd(20,51) = pd(20,51) - rrt(380) * density(20) 
  pd(23,20) = pd(23,20) + rrt(380) * density(51) 
  pd(23,51) = pd(23,51) + rrt(380) * density(20) 
  pd(37,20) = pd(37,20) + rrt(380) * density(51) 
  pd(37,51) = pd(37,51) + rrt(380) * density(20) 
  pd(51,20) = pd(51,20) - rrt(380) * density(51) 
  pd(51,51) = pd(51,51) - rrt(380) * density(20) 
  pd(11,13) = pd(11,13) + rrt(381) * density(51) 
  pd(11,51) = pd(11,51) + rrt(381) * density(13) 
  pd(13,13) = pd(13,13) - rrt(381) * density(51) 
  pd(13,51) = pd(13,51) - rrt(381) * density(13) 
  pd(23,13) = pd(23,13) + rrt(381) * density(51) 
  pd(23,51) = pd(23,51) + rrt(381) * density(13) 
  pd(51,13) = pd(51,13) - rrt(381) * density(51) 
  pd(51,51) = pd(51,51) - rrt(381) * density(13) 
  pd(13,23) = pd(13,23) + rrt(382) * density(51) 
  pd(13,51) = pd(13,51) + rrt(382) * density(23) 
  pd(23,23) = pd(23,23) - rrt(382) * density(51) 
  pd(23,51) = pd(23,51) - rrt(382) * density(23) 
  pd(43,23) = pd(43,23) + rrt(382) * density(51) 
  pd(43,51) = pd(43,51) + rrt(382) * density(23) 
  pd(51,23) = pd(51,23) - rrt(382) * density(51) 
  pd(51,51) = pd(51,51) - rrt(382) * density(23) 
  pd(23,31) = pd(23,31) + rrt(383) * density(51) 
  pd(23,51) = pd(23,51) + rrt(383) * density(31) 
  pd(31,31) = pd(31,31) - rrt(383) * density(51) 
  pd(31,51) = pd(31,51) - rrt(383) * density(31) 
  pd(49,31) = pd(49,31) + rrt(383) * density(51) 
  pd(49,51) = pd(49,51) + rrt(383) * density(31) 
  pd(51,31) = pd(51,31) - rrt(383) * density(51) 
  pd(51,51) = pd(51,51) - rrt(383) * density(31) 
  pd(28,28) = pd(28,28) - rrt(384) * density(51) 
  pd(28,51) = pd(28,51) - rrt(384) * density(28) 
  pd(43,28) = pd(43,28) + rrt(384) * density(51) 
  pd(43,51) = pd(43,51) + rrt(384) * density(28) 
  pd(51,28) = pd(51,28) - rrt(384) * density(51) 
  pd(51,51) = pd(51,51) - rrt(384) * density(28) 
  pd(13,51) = pd(13,51) + rrt(385) 
  pd(23,51) = pd(23,51) + rrt(385) 
  pd(51,51) = pd(51,51) - rrt(385) 
  pd(13,23) = pd(13,23) + rrt(386) 
  pd(23,23) = pd(23,23) - rrt(386) 
  pd(28,23) = pd(28,23) + rrt(386) 
  pd(08,23) = pd(08,23) + rrt(387) 
  pd(11,23) = pd(11,23) + rrt(387) 
  pd(23,23) = pd(23,23) - rrt(387) 
  pd(23,23) = pd(23,23) - rrt(388) * density(34) 
  pd(23,34) = pd(23,34) - rrt(388) * density(23) 
  pd(28,23) = pd(28,23) + rrt(388) * density(34) 
  pd(28,34) = pd(28,34) + rrt(388) * density(23) 
  pd(34,23) = pd(34,23) - rrt(388) * density(34) 
  pd(34,34) = pd(34,34) - rrt(388) * density(23) 
  pd(43,23) = pd(43,23) + rrt(388) * density(34) 
  pd(43,34) = pd(43,34) + rrt(388) * density(23) 
  pd(13,28) = pd(13,28) + rrt(389) * density(28) * 4.0d0
  pd(21,28) = pd(21,28) + rrt(389) * density(28) * 2.0d0
  pd(28,28) = pd(28,28) - rrt(389) * density(28) * 4.0d0
  pd(23,28) = pd(23,28) + rrt(390) * density(34) 
  pd(23,34) = pd(23,34) + rrt(390) * density(28) 
  pd(28,28) = pd(28,28) - rrt(390) * density(34) 
  pd(28,34) = pd(28,34) - rrt(390) * density(28) 
  pd(34,28) = pd(34,28) - rrt(390) * density(34) 
  pd(34,34) = pd(34,34) - rrt(390) * density(28) 
  pd(42,28) = pd(42,28) + rrt(390) * density(34) 
  pd(42,34) = pd(42,34) + rrt(390) * density(28) 
  pd(21,28) = pd(21,28) + rrt(391) * density(29) 
  pd(21,29) = pd(21,29) + rrt(391) * density(28) 
  pd(23,28) = pd(23,28) + rrt(391) * density(29) 
  pd(23,29) = pd(23,29) + rrt(391) * density(28) 
  pd(28,28) = pd(28,28) - rrt(391) * density(29) 
  pd(28,29) = pd(28,29) - rrt(391) * density(28) 
  pd(29,28) = pd(29,28) - rrt(391) * density(29) 
  pd(29,29) = pd(29,29) - rrt(391) * density(28) 
  pd(06,06) = pd(06,06) - rrt(392) * density(28) 
  pd(06,28) = pd(06,28) - rrt(392) * density(06) 
  pd(08,06) = pd(08,06) + rrt(392) * density(28) 
  pd(08,28) = pd(08,28) + rrt(392) * density(06) 
  pd(21,06) = pd(21,06) + rrt(392) * density(28) 
  pd(21,28) = pd(21,28) + rrt(392) * density(06) 
  pd(28,06) = pd(28,06) - rrt(392) * density(28) 
  pd(28,28) = pd(28,28) - rrt(392) * density(06) 
  pd(18,28) = pd(18,28) + rrt(393) * density(45) 
  pd(18,45) = pd(18,45) + rrt(393) * density(28) 
  pd(23,28) = pd(23,28) + rrt(393) * density(45) 
  pd(23,45) = pd(23,45) + rrt(393) * density(28) 
  pd(28,28) = pd(28,28) - rrt(393) * density(45) 
  pd(28,45) = pd(28,45) - rrt(393) * density(28) 
  pd(45,28) = pd(45,28) - rrt(393) * density(45) 
  pd(45,45) = pd(45,45) - rrt(393) * density(28) 
  pd(18,18) = pd(18,18) - rrt(394) * density(28) 
  pd(18,28) = pd(18,28) - rrt(394) * density(18) 
  pd(28,18) = pd(28,18) - rrt(394) * density(28) 
  pd(28,28) = pd(28,28) - rrt(394) * density(18) 
  pd(34,18) = pd(34,18) + rrt(394) * density(28) 
  pd(34,28) = pd(34,28) + rrt(394) * density(18) 
  pd(42,18) = pd(42,18) + rrt(394) * density(28) 
  pd(42,28) = pd(42,28) + rrt(394) * density(18) 
  pd(18,18) = pd(18,18) - rrt(395) * density(28) 
  pd(18,28) = pd(18,28) - rrt(395) * density(18) 
  pd(23,18) = pd(23,18) + rrt(395) * density(28) 
  pd(23,28) = pd(23,28) + rrt(395) * density(18) 
  pd(28,18) = pd(28,18) - rrt(395) * density(28) 
  pd(28,28) = pd(28,28) - rrt(395) * density(18) 
  pd(37,18) = pd(37,18) + rrt(395) * density(28) 
  pd(37,28) = pd(37,28) + rrt(395) * density(18) 
  pd(20,28) = pd(20,28) + rrt(396) * density(37) 
  pd(20,37) = pd(20,37) + rrt(396) * density(28) 
  pd(23,28) = pd(23,28) + rrt(396) * density(37) 
  pd(23,37) = pd(23,37) + rrt(396) * density(28) 
  pd(28,28) = pd(28,28) - rrt(396) * density(37) 
  pd(28,37) = pd(28,37) - rrt(396) * density(28) 
  pd(37,28) = pd(37,28) - rrt(396) * density(37) 
  pd(37,37) = pd(37,37) - rrt(396) * density(28) 
  pd(11,11) = pd(11,11) - rrt(397) * density(28) 
  pd(11,28) = pd(11,28) - rrt(397) * density(11) 
  pd(13,11) = pd(13,11) + rrt(397) * density(28) 
  pd(13,28) = pd(13,28) + rrt(397) * density(11) 
  pd(23,11) = pd(23,11) + rrt(397) * density(28) 
  pd(23,28) = pd(23,28) + rrt(397) * density(11) 
  pd(28,11) = pd(28,11) - rrt(397) * density(28) 
  pd(28,28) = pd(28,28) - rrt(397) * density(11) 
  pd(08,13) = pd(08,13) + rrt(398) * density(28) 
  pd(08,28) = pd(08,28) + rrt(398) * density(13) 
  pd(11,13) = pd(11,13) + rrt(398) * density(28) 
  pd(11,28) = pd(11,28) + rrt(398) * density(13) 
  pd(13,13) = pd(13,13) - rrt(398) * density(28) 
  pd(13,28) = pd(13,28) - rrt(398) * density(13) 
  pd(28,13) = pd(28,13) - rrt(398) * density(28) 
  pd(28,28) = pd(28,28) - rrt(398) * density(13) 
  pd(08,28) = pd(08,28) + rrt(399) 
  pd(13,28) = pd(13,28) + rrt(399) 
  pd(28,28) = pd(28,28) - rrt(399) 
  pd(11,28) = pd(11,28) + rrt(400) * density(28) * 2.0d0
  pd(21,28) = pd(21,28) + rrt(400) * density(28) * 2.0d0
  pd(28,28) = pd(28,28) - rrt(400) * density(28) * 4.0d0
  pd(13,13) = pd(13,13) - rrt(401) * density(28) 
  pd(13,28) = pd(13,28) - rrt(401) * density(13) 
  pd(23,13) = pd(23,13) + rrt(401) * density(28) 
  pd(23,28) = pd(23,28) + rrt(401) * density(13) 
  pd(28,13) = pd(28,13) - rrt(401) * density(28) 
  pd(28,28) = pd(28,28) - rrt(401) * density(13) 
  pd(08,08) = pd(08,08) - rrt(402) * density(43) 
  pd(08,43) = pd(08,43) - rrt(402) * density(08) 
  pd(13,08) = pd(13,08) + rrt(402) * density(43) 
  pd(13,43) = pd(13,43) + rrt(402) * density(08) 
  pd(37,08) = pd(37,08) + rrt(402) * density(43) 
  pd(37,43) = pd(37,43) + rrt(402) * density(08) 
  pd(43,08) = pd(43,08) - rrt(402) * density(43) 
  pd(43,43) = pd(43,43) - rrt(402) * density(08) 
  pd(08,08) = pd(08,08) - rrt(403) * density(43) 
  pd(08,43) = pd(08,43) - rrt(403) * density(08) 
  pd(18,08) = pd(18,08) + rrt(403) * density(43) 
  pd(18,43) = pd(18,43) + rrt(403) * density(08) 
  pd(43,08) = pd(43,08) - rrt(403) * density(43) 
  pd(43,43) = pd(43,43) - rrt(403) * density(08) 
  pd(08,08) = pd(08,08) - rrt(404) * density(11) 
  pd(08,11) = pd(08,11) - rrt(404) * density(08) 
  pd(11,08) = pd(11,08) - rrt(404) * density(11) 
  pd(11,11) = pd(11,11) - rrt(404) * density(08) 
  pd(13,08) = pd(13,08) + rrt(404) * density(11) 
  pd(13,11) = pd(13,11) + rrt(404) * density(08) 
  pd(28,08) = pd(28,08) + rrt(404) * density(11) 
  pd(28,11) = pd(28,11) + rrt(404) * density(08) 
  pd(08,08) = pd(08,08) - rrt(405) * density(23) 
  pd(08,23) = pd(08,23) - rrt(405) * density(08) 
  pd(13,08) = pd(13,08) + rrt(405) * density(23) 
  pd(13,23) = pd(13,23) + rrt(405) * density(08) 
  pd(23,08) = pd(23,08) - rrt(405) * density(23) 
  pd(23,23) = pd(23,23) - rrt(405) * density(08) 
  pd(29,08) = pd(29,08) + rrt(405) * density(23) 
  pd(29,23) = pd(29,23) + rrt(405) * density(08) 
  pd(08,08) = pd(08,08) - rrt(406) * density(28) 
  pd(08,28) = pd(08,28) - rrt(406) * density(08) 
  pd(13,08) = pd(13,08) + rrt(406) * density(28) 
  pd(13,28) = pd(13,28) + rrt(406) * density(08) 
  pd(21,08) = pd(21,08) + rrt(406) * density(28) 
  pd(21,28) = pd(21,28) + rrt(406) * density(08) 
  pd(28,08) = pd(28,08) - rrt(406) * density(28) 
  pd(28,28) = pd(28,28) - rrt(406) * density(08) 
  pd(08,08) = pd(08,08) - rrt(407) * density(11) 
  pd(08,11) = pd(08,11) - rrt(407) * density(08) 
  pd(11,08) = pd(11,08) - rrt(407) * density(11) 
  pd(11,11) = pd(11,11) - rrt(407) * density(08) 
  pd(23,08) = pd(23,08) + rrt(407) * density(11) 
  pd(23,11) = pd(23,11) + rrt(407) * density(08) 
  pd(08,08) = pd(08,08) - rrt(408) * density(29) 
  pd(08,29) = pd(08,29) - rrt(408) * density(08) 
  pd(21,08) = pd(21,08) + rrt(408) * density(29) 
  pd(21,29) = pd(21,29) + rrt(408) * density(08) 
  pd(28,08) = pd(28,08) + rrt(408) * density(29) 
  pd(28,29) = pd(28,29) + rrt(408) * density(08) 
  pd(29,08) = pd(29,08) - rrt(408) * density(29) 
  pd(29,29) = pd(29,29) - rrt(408) * density(08) 
  pd(29,29) = pd(29,29) - rrt(409) * density(43) 
  pd(29,43) = pd(29,43) - rrt(409) * density(29) 
  pd(34,29) = pd(34,29) + rrt(409) * density(43) 
  pd(34,43) = pd(34,43) + rrt(409) * density(29) 
  pd(42,29) = pd(42,29) + rrt(409) * density(43) 
  pd(42,43) = pd(42,43) + rrt(409) * density(29) 
  pd(43,29) = pd(43,29) - rrt(409) * density(43) 
  pd(43,43) = pd(43,43) - rrt(409) * density(29) 
  pd(18,18) = pd(18,18) - rrt(410) * density(43) 
  pd(18,43) = pd(18,43) - rrt(410) * density(18) 
  pd(34,18) = pd(34,18) + rrt(410) * density(43) 
  pd(34,43) = pd(34,43) + rrt(410) * density(18) 
  pd(43,18) = pd(43,18) - rrt(410) * density(43) 
  pd(43,43) = pd(43,43) - rrt(410) * density(18) 
  pd(45,18) = pd(45,18) + rrt(410) * density(43) 
  pd(45,43) = pd(45,43) + rrt(410) * density(18) 
  pd(20,20) = pd(20,20) - rrt(411) * density(43) 
  pd(20,43) = pd(20,43) - rrt(411) * density(20) 
  pd(34,20) = pd(34,20) + rrt(411) * density(43) 
  pd(34,43) = pd(34,43) + rrt(411) * density(20) 
  pd(37,20) = pd(37,20) + rrt(411) * density(43) 
  pd(37,43) = pd(37,43) + rrt(411) * density(20) 
  pd(43,20) = pd(43,20) - rrt(411) * density(43) 
  pd(43,43) = pd(43,43) - rrt(411) * density(20) 
  pd(11,13) = pd(11,13) + rrt(412) * density(43) 
  pd(11,43) = pd(11,43) + rrt(412) * density(13) 
  pd(13,13) = pd(13,13) - rrt(412) * density(43) 
  pd(13,43) = pd(13,43) - rrt(412) * density(13) 
  pd(34,13) = pd(34,13) + rrt(412) * density(43) 
  pd(34,43) = pd(34,43) + rrt(412) * density(13) 
  pd(43,13) = pd(43,13) - rrt(412) * density(43) 
  pd(43,43) = pd(43,43) - rrt(412) * density(13) 
  pd(13,13) = pd(13,13) - rrt(413) * density(43) 
  pd(13,43) = pd(13,43) - rrt(413) * density(13) 
  pd(23,13) = pd(23,13) + rrt(413) * density(43) 
  pd(23,43) = pd(23,43) + rrt(413) * density(13) 
  pd(43,13) = pd(43,13) - rrt(413) * density(43) 
  pd(43,43) = pd(43,43) - rrt(413) * density(13) 
  pd(51,13) = pd(51,13) + rrt(413) * density(43) 
  pd(51,43) = pd(51,43) + rrt(413) * density(13) 
  pd(23,43) = pd(23,43) + rrt(414) * 2.0d0
  pd(43,43) = pd(43,43) - rrt(414) 
  pd(31,31) = pd(31,31) - rrt(415) * density(43) 
  pd(31,43) = pd(31,43) - rrt(415) * density(31) 
  pd(34,31) = pd(34,31) + rrt(415) * density(43) 
  pd(34,43) = pd(34,43) + rrt(415) * density(31) 
  pd(43,31) = pd(43,31) - rrt(415) * density(43) 
  pd(43,43) = pd(43,43) - rrt(415) * density(31) 
  pd(49,31) = pd(49,31) + rrt(415) * density(43) 
  pd(49,43) = pd(49,43) + rrt(415) * density(31) 
  pd(08,08) = pd(08,08) - rrt(416) * density(43) 
  pd(08,43) = pd(08,43) - rrt(416) * density(08) 
  pd(23,08) = pd(23,08) + rrt(416) * density(43) 
  pd(23,43) = pd(23,43) + rrt(416) * density(08) 
  pd(42,08) = pd(42,08) + rrt(416) * density(43) 
  pd(42,43) = pd(42,43) + rrt(416) * density(08) 
  pd(43,08) = pd(43,08) - rrt(416) * density(43) 
  pd(43,43) = pd(43,43) - rrt(416) * density(08) 
  pd(23,28) = pd(23,28) + rrt(417) * density(43) 
  pd(23,43) = pd(23,43) + rrt(417) * density(28) 
  pd(28,28) = pd(28,28) - rrt(417) * density(43) 
  pd(28,43) = pd(28,43) - rrt(417) * density(28) 
  pd(34,28) = pd(34,28) + rrt(417) * density(43) 
  pd(34,43) = pd(34,43) + rrt(417) * density(28) 
  pd(43,28) = pd(43,28) - rrt(417) * density(43) 
  pd(43,43) = pd(43,43) - rrt(417) * density(28) 
  pd(34,34) = pd(34,34) - rrt(418) * density(34) * 4.0d0
  pd(42,34) = pd(42,34) + rrt(418) * density(34) * 2.0d0
  pd(43,34) = pd(43,34) + rrt(418) * density(34) * 2.0d0
  pd(29,34) = pd(29,34) + rrt(419) * density(42) 
  pd(29,42) = pd(29,42) + rrt(419) * density(34) 
  pd(34,34) = pd(34,34) - rrt(419) * density(42) 
  pd(34,42) = pd(34,42) - rrt(419) * density(34) 
  pd(42,34) = pd(42,34) - rrt(419) * density(42) 
  pd(42,42) = pd(42,42) - rrt(419) * density(34) 
  pd(43,34) = pd(43,34) + rrt(419) * density(42) 
  pd(43,42) = pd(43,42) + rrt(419) * density(34) 
  pd(06,21) = pd(06,21) + rrt(420) * density(34) 
  pd(06,34) = pd(06,34) + rrt(420) * density(21) 
  pd(21,21) = pd(21,21) - rrt(420) * density(34) 
  pd(21,34) = pd(21,34) - rrt(420) * density(21) 
  pd(34,21) = pd(34,21) - rrt(420) * density(34) 
  pd(34,34) = pd(34,34) - rrt(420) * density(21) 
  pd(43,21) = pd(43,21) + rrt(420) * density(34) 
  pd(43,34) = pd(43,34) + rrt(420) * density(21) 
  pd(06,06) = pd(06,06) - rrt(421) * density(34) 
  pd(06,34) = pd(06,34) - rrt(421) * density(06) 
  pd(21,06) = pd(21,06) + rrt(421) * density(34) 
  pd(21,34) = pd(21,34) + rrt(421) * density(06) 
  pd(34,06) = pd(34,06) - rrt(421) * density(34) 
  pd(34,34) = pd(34,34) - rrt(421) * density(06) 
  pd(42,06) = pd(42,06) + rrt(421) * density(34) 
  pd(42,34) = pd(42,34) + rrt(421) * density(06) 
  pd(18,34) = pd(18,34) + rrt(422) * density(45) 
  pd(18,45) = pd(18,45) + rrt(422) * density(34) 
  pd(34,34) = pd(34,34) - rrt(422) * density(45) 
  pd(34,45) = pd(34,45) - rrt(422) * density(34) 
  pd(43,34) = pd(43,34) + rrt(422) * density(45) 
  pd(43,45) = pd(43,45) + rrt(422) * density(34) 
  pd(45,34) = pd(45,34) - rrt(422) * density(45) 
  pd(45,45) = pd(45,45) - rrt(422) * density(34) 
  pd(18,18) = pd(18,18) - rrt(423) * density(34) 
  pd(18,34) = pd(18,34) - rrt(423) * density(18) 
  pd(34,18) = pd(34,18) - rrt(423) * density(34) 
  pd(34,34) = pd(34,34) - rrt(423) * density(18) 
  pd(42,18) = pd(42,18) + rrt(423) * density(34) 
  pd(42,34) = pd(42,34) + rrt(423) * density(18) 
  pd(45,18) = pd(45,18) + rrt(423) * density(34) 
  pd(45,34) = pd(45,34) + rrt(423) * density(18) 
  pd(18,18) = pd(18,18) - rrt(424) * density(34) 
  pd(18,34) = pd(18,34) - rrt(424) * density(18) 
  pd(34,18) = pd(34,18) - rrt(424) * density(34) 
  pd(34,34) = pd(34,34) - rrt(424) * density(18) 
  pd(37,18) = pd(37,18) + rrt(424) * density(34) 
  pd(37,34) = pd(37,34) + rrt(424) * density(18) 
  pd(43,18) = pd(43,18) + rrt(424) * density(34) 
  pd(43,34) = pd(43,34) + rrt(424) * density(18) 
  pd(20,34) = pd(20,34) + rrt(425) * density(37) 
  pd(20,37) = pd(20,37) + rrt(425) * density(34) 
  pd(34,34) = pd(34,34) - rrt(425) * density(37) 
  pd(34,37) = pd(34,37) - rrt(425) * density(34) 
  pd(37,34) = pd(37,34) - rrt(425) * density(37) 
  pd(37,37) = pd(37,37) - rrt(425) * density(34) 
  pd(43,34) = pd(43,34) + rrt(425) * density(37) 
  pd(43,37) = pd(43,37) + rrt(425) * density(34) 
  pd(11,11) = pd(11,11) - rrt(426) * density(34) 
  pd(11,34) = pd(11,34) - rrt(426) * density(11) 
  pd(13,11) = pd(13,11) + rrt(426) * density(34) 
  pd(13,34) = pd(13,34) + rrt(426) * density(11) 
  pd(34,11) = pd(34,11) - rrt(426) * density(34) 
  pd(34,34) = pd(34,34) - rrt(426) * density(11) 
  pd(43,11) = pd(43,11) + rrt(426) * density(34) 
  pd(43,34) = pd(43,34) + rrt(426) * density(11) 
  pd(13,13) = pd(13,13) - rrt(427) * density(34) 
  pd(13,34) = pd(13,34) - rrt(427) * density(13) 
  pd(23,13) = pd(23,13) + rrt(427) * density(34) * 2.0d0
  pd(23,34) = pd(23,34) + rrt(427) * density(13) * 2.0d0
  pd(34,13) = pd(34,13) - rrt(427) * density(34) 
  pd(34,34) = pd(34,34) - rrt(427) * density(13) 
  pd(11,13) = pd(11,13) + rrt(428) * density(34) 
  pd(11,34) = pd(11,34) + rrt(428) * density(13) 
  pd(13,13) = pd(13,13) - rrt(428) * density(34) 
  pd(13,34) = pd(13,34) - rrt(428) * density(13) 
  pd(34,13) = pd(34,13) - rrt(428) * density(34) 
  pd(34,34) = pd(34,34) - rrt(428) * density(13) 
  pd(42,13) = pd(42,13) + rrt(428) * density(34) 
  pd(42,34) = pd(42,34) + rrt(428) * density(13) 
  pd(13,13) = pd(13,13) - rrt(429) * density(34) 
  pd(13,34) = pd(13,34) - rrt(429) * density(13) 
  pd(34,13) = pd(34,13) - rrt(429) * density(34) 
  pd(34,34) = pd(34,34) - rrt(429) * density(13) 
  pd(43,13) = pd(43,13) + rrt(429) * density(34) 
  pd(43,34) = pd(43,34) + rrt(429) * density(13) 
  pd(13,34) = pd(13,34) + rrt(430) 
  pd(34,34) = pd(34,34) - rrt(430) 
  pd(42,34) = pd(42,34) + rrt(430) 
  pd(34,34) = pd(34,34) - rrt(431) * density(34) * 4.0d0
  pd(49,34) = pd(49,34) + rrt(431) * density(34) * 2.0d0
  pd(31,31) = pd(31,31) - rrt(432) * density(34) 
  pd(31,34) = pd(31,34) - rrt(432) * density(31) 
  pd(34,31) = pd(34,31) - rrt(432) * density(34) 
  pd(34,34) = pd(34,34) - rrt(432) * density(31) 
  pd(42,31) = pd(42,31) + rrt(432) * density(34) 
  pd(42,34) = pd(42,34) + rrt(432) * density(31) 
  pd(49,31) = pd(49,31) + rrt(432) * density(34) 
  pd(49,34) = pd(49,34) + rrt(432) * density(31) 
  pd(29,29) = pd(29,29) - rrt(433) * density(34) 
  pd(29,34) = pd(29,34) - rrt(433) * density(29) 
  pd(34,29) = pd(34,29) - rrt(433) * density(34) 
  pd(34,34) = pd(34,34) - rrt(433) * density(29) 
  pd(42,29) = pd(42,29) + rrt(433) * density(34) * 2.0d0
  pd(42,34) = pd(42,34) + rrt(433) * density(29) * 2.0d0
  pd(11,13) = pd(11,13) + rrt(434) * density(42) 
  pd(11,42) = pd(11,42) + rrt(434) * density(13) 
  pd(13,13) = pd(13,13) - rrt(434) * density(42) 
  pd(13,42) = pd(13,42) - rrt(434) * density(13) 
  pd(29,13) = pd(29,13) + rrt(434) * density(42) 
  pd(29,42) = pd(29,42) + rrt(434) * density(13) 
  pd(42,13) = pd(42,13) - rrt(434) * density(42) 
  pd(42,42) = pd(42,42) - rrt(434) * density(13) 
  pd(13,13) = pd(13,13) - rrt(435) * density(42) 
  pd(13,42) = pd(13,42) - rrt(435) * density(13) 
  pd(34,13) = pd(34,13) + rrt(435) * density(42) 
  pd(34,42) = pd(34,42) + rrt(435) * density(13) 
  pd(42,13) = pd(42,13) - rrt(435) * density(42) 
  pd(42,42) = pd(42,42) - rrt(435) * density(13) 
  pd(11,11) = pd(11,11) - rrt(436) * density(42) 
  pd(11,42) = pd(11,42) - rrt(436) * density(11) 
  pd(13,11) = pd(13,11) + rrt(436) * density(42) 
  pd(13,42) = pd(13,42) + rrt(436) * density(11) 
  pd(34,11) = pd(34,11) + rrt(436) * density(42) 
  pd(34,42) = pd(34,42) + rrt(436) * density(11) 
  pd(42,11) = pd(42,11) - rrt(436) * density(42) 
  pd(42,42) = pd(42,42) - rrt(436) * density(11) 
  pd(13,42) = pd(13,42) + rrt(437) 
  pd(29,42) = pd(29,42) + rrt(437) 
  pd(42,42) = pd(42,42) - rrt(437) 
  pd(20,37) = pd(20,37) + rrt(438) * density(42) 
  pd(20,42) = pd(20,42) + rrt(438) * density(37) 
  pd(34,37) = pd(34,37) + rrt(438) * density(42) 
  pd(34,42) = pd(34,42) + rrt(438) * density(37) 
  pd(37,37) = pd(37,37) - rrt(438) * density(42) 
  pd(37,42) = pd(37,42) - rrt(438) * density(37) 
  pd(42,37) = pd(42,37) - rrt(438) * density(42) 
  pd(42,42) = pd(42,42) - rrt(438) * density(37) 
  pd(21,21) = pd(21,21) - rrt(439) * density(42) 
  pd(21,42) = pd(21,42) - rrt(439) * density(21) 
  pd(29,21) = pd(29,21) + rrt(439) * density(42) * 2.0d0
  pd(29,42) = pd(29,42) + rrt(439) * density(21) * 2.0d0
  pd(42,21) = pd(42,21) - rrt(439) * density(42) 
  pd(42,42) = pd(42,42) - rrt(439) * density(21) 
  pd(18,37) = pd(18,37) + rrt(440) * density(42) 
  pd(18,42) = pd(18,42) + rrt(440) * density(37) 
  pd(29,37) = pd(29,37) + rrt(440) * density(42) 
  pd(29,42) = pd(29,42) + rrt(440) * density(37) 
  pd(37,37) = pd(37,37) - rrt(440) * density(42) 
  pd(37,42) = pd(37,42) - rrt(440) * density(37) 
  pd(42,37) = pd(42,37) - rrt(440) * density(42) 
  pd(42,42) = pd(42,42) - rrt(440) * density(37) 
  pd(29,42) = pd(29,42) + rrt(441) * density(42) * 2.0d0
  pd(34,42) = pd(34,42) + rrt(441) * density(42) * 2.0d0
  pd(42,42) = pd(42,42) - rrt(441) * density(42) * 4.0d0
  pd(18,23) = pd(18,23) + rrt(442) * density(42) 
  pd(18,42) = pd(18,42) + rrt(442) * density(23) 
  pd(23,23) = pd(23,23) - rrt(442) * density(42) 
  pd(23,42) = pd(23,42) - rrt(442) * density(23) 
  pd(42,23) = pd(42,23) - rrt(442) * density(42) 
  pd(42,42) = pd(42,42) - rrt(442) * density(23) 
  pd(11,42) = pd(11,42) + rrt(443) 
  pd(21,42) = pd(21,42) + rrt(443) 
  pd(42,42) = pd(42,42) - rrt(443) 
  pd(31,34) = pd(31,34) + rrt(444) * density(42) 
  pd(31,42) = pd(31,42) + rrt(444) * density(34) 
  pd(34,34) = pd(34,34) - rrt(444) * density(42) 
  pd(34,42) = pd(34,42) - rrt(444) * density(34) 
  pd(42,34) = pd(42,34) - rrt(444) * density(42) 
  pd(42,42) = pd(42,42) - rrt(444) * density(34) 
  pd(11,11) = pd(11,11) - rrt(445) * density(42) 
  pd(11,42) = pd(11,42) - rrt(445) * density(11) 
  pd(42,11) = pd(42,11) - rrt(445) * density(42) 
  pd(42,42) = pd(42,42) - rrt(445) * density(11) 
  pd(43,11) = pd(43,11) + rrt(445) * density(42) 
  pd(43,42) = pd(43,42) + rrt(445) * density(11) 
  pd(28,28) = pd(28,28) - rrt(446) * density(42) 
  pd(28,42) = pd(28,42) - rrt(446) * density(28) 
  pd(37,28) = pd(37,28) + rrt(446) * density(42) 
  pd(37,42) = pd(37,42) + rrt(446) * density(28) 
  pd(42,28) = pd(42,28) - rrt(446) * density(42) 
  pd(42,42) = pd(42,42) - rrt(446) * density(28) 
  pd(18,31) = pd(18,31) + rrt(447) * density(42) 
  pd(18,42) = pd(18,42) + rrt(447) * density(31) 
  pd(31,31) = pd(31,31) - rrt(447) * density(42) 
  pd(31,42) = pd(31,42) - rrt(447) * density(31) 
  pd(37,31) = pd(37,31) + rrt(447) * density(42) 
  pd(37,42) = pd(37,42) + rrt(447) * density(31) 
  pd(42,31) = pd(42,31) - rrt(447) * density(42) 
  pd(42,42) = pd(42,42) - rrt(447) * density(31) 
  pd(21,29) = pd(21,29) + rrt(448) * density(29) * 2.0d0
  pd(29,29) = pd(29,29) - rrt(448) * density(29) * 4.0d0
  pd(42,29) = pd(42,29) + rrt(448) * density(29) * 2.0d0
  pd(18,29) = pd(18,29) + rrt(449) * density(45) 
  pd(18,45) = pd(18,45) + rrt(449) * density(29) 
  pd(29,29) = pd(29,29) - rrt(449) * density(45) 
  pd(29,45) = pd(29,45) - rrt(449) * density(29) 
  pd(42,29) = pd(42,29) + rrt(449) * density(45) 
  pd(42,45) = pd(42,45) + rrt(449) * density(29) 
  pd(45,29) = pd(45,29) - rrt(449) * density(45) 
  pd(45,45) = pd(45,45) - rrt(449) * density(29) 
  pd(18,18) = pd(18,18) - rrt(450) * density(29) 
  pd(18,29) = pd(18,29) - rrt(450) * density(18) 
  pd(21,18) = pd(21,18) + rrt(450) * density(29) 
  pd(21,29) = pd(21,29) + rrt(450) * density(18) 
  pd(29,18) = pd(29,18) - rrt(450) * density(29) 
  pd(29,29) = pd(29,29) - rrt(450) * density(18) 
  pd(45,18) = pd(45,18) + rrt(450) * density(29) 
  pd(45,29) = pd(45,29) + rrt(450) * density(18) 
  pd(18,18) = pd(18,18) - rrt(451) * density(29) 
  pd(18,29) = pd(18,29) - rrt(451) * density(18) 
  pd(29,18) = pd(29,18) - rrt(451) * density(29) 
  pd(29,29) = pd(29,29) - rrt(451) * density(18) 
  pd(37,18) = pd(37,18) + rrt(451) * density(29) 
  pd(37,29) = pd(37,29) + rrt(451) * density(18) 
  pd(42,18) = pd(42,18) + rrt(451) * density(29) 
  pd(42,29) = pd(42,29) + rrt(451) * density(18) 
  pd(20,29) = pd(20,29) + rrt(452) * density(37) 
  pd(20,37) = pd(20,37) + rrt(452) * density(29) 
  pd(29,29) = pd(29,29) - rrt(452) * density(37) 
  pd(29,37) = pd(29,37) - rrt(452) * density(29) 
  pd(37,29) = pd(37,29) - rrt(452) * density(37) 
  pd(37,37) = pd(37,37) - rrt(452) * density(29) 
  pd(42,29) = pd(42,29) + rrt(452) * density(37) 
  pd(42,37) = pd(42,37) + rrt(452) * density(29) 
  pd(20,20) = pd(20,20) - rrt(453) * density(29) 
  pd(20,29) = pd(20,29) - rrt(453) * density(20) 
  pd(21,20) = pd(21,20) + rrt(453) * density(29) 
  pd(21,29) = pd(21,29) + rrt(453) * density(20) 
  pd(29,20) = pd(29,20) - rrt(453) * density(29) 
  pd(29,29) = pd(29,29) - rrt(453) * density(20) 
  pd(37,20) = pd(37,20) + rrt(453) * density(29) 
  pd(37,29) = pd(37,29) + rrt(453) * density(20) 
  pd(11,11) = pd(11,11) - rrt(454) * density(29) 
  pd(11,29) = pd(11,29) - rrt(454) * density(11) 
  pd(13,11) = pd(13,11) + rrt(454) * density(29) 
  pd(13,29) = pd(13,29) + rrt(454) * density(11) 
  pd(29,11) = pd(29,11) - rrt(454) * density(29) 
  pd(29,29) = pd(29,29) - rrt(454) * density(11) 
  pd(42,11) = pd(42,11) + rrt(454) * density(29) 
  pd(42,29) = pd(42,29) + rrt(454) * density(11) 
  pd(11,13) = pd(11,13) + rrt(455) * density(29) 
  pd(11,29) = pd(11,29) + rrt(455) * density(13) 
  pd(13,13) = pd(13,13) - rrt(455) * density(29) 
  pd(13,29) = pd(13,29) - rrt(455) * density(13) 
  pd(21,13) = pd(21,13) + rrt(455) * density(29) 
  pd(21,29) = pd(21,29) + rrt(455) * density(13) 
  pd(29,13) = pd(29,13) - rrt(455) * density(29) 
  pd(29,29) = pd(29,29) - rrt(455) * density(13) 
  pd(13,13) = pd(13,13) - rrt(456) * density(29) 
  pd(13,29) = pd(13,29) - rrt(456) * density(13) 
  pd(29,13) = pd(29,13) - rrt(456) * density(29) 
  pd(29,29) = pd(29,29) - rrt(456) * density(13) 
  pd(42,13) = pd(42,13) + rrt(456) * density(29) 
  pd(42,29) = pd(42,29) + rrt(456) * density(13) 
  pd(13,29) = pd(13,29) + rrt(457) 
  pd(21,29) = pd(21,29) + rrt(457) 
  pd(29,29) = pd(29,29) - rrt(457) 
  pd(21,29) = pd(21,29) + rrt(458) * density(31) 
  pd(21,31) = pd(21,31) + rrt(458) * density(29) 
  pd(29,29) = pd(29,29) - rrt(458) * density(31) 
  pd(29,31) = pd(29,31) - rrt(458) * density(29) 
  pd(31,29) = pd(31,29) - rrt(458) * density(31) 
  pd(31,31) = pd(31,31) - rrt(458) * density(29) 
  pd(49,29) = pd(49,29) + rrt(458) * density(31) 
  pd(49,31) = pd(49,31) + rrt(458) * density(29) 
  pd(13,13) = pd(13,13) - rrt(459) * density(21) 
  pd(13,21) = pd(13,21) - rrt(459) * density(13) 
  pd(21,13) = pd(21,13) - rrt(459) * density(21) 
  pd(21,21) = pd(21,21) - rrt(459) * density(13) 
  pd(29,13) = pd(29,13) + rrt(459) * density(21) 
  pd(29,21) = pd(29,21) + rrt(459) * density(13) 
  pd(11,11) = pd(11,11) - rrt(460) * density(21) 
  pd(11,21) = pd(11,21) - rrt(460) * density(11) 
  pd(21,11) = pd(21,11) - rrt(460) * density(21) 
  pd(21,21) = pd(21,21) - rrt(460) * density(11) 
  pd(42,11) = pd(42,11) + rrt(460) * density(21) 
  pd(42,21) = pd(42,21) + rrt(460) * density(11) 
  pd(11,11) = pd(11,11) - rrt(461) * density(21) 
  pd(11,21) = pd(11,21) - rrt(461) * density(11) 
  pd(13,11) = pd(13,11) + rrt(461) * density(21) 
  pd(13,21) = pd(13,21) + rrt(461) * density(11) 
  pd(21,11) = pd(21,11) - rrt(461) * density(21) 
  pd(21,21) = pd(21,21) - rrt(461) * density(11) 
  pd(29,11) = pd(29,11) + rrt(461) * density(21) 
  pd(29,21) = pd(29,21) + rrt(461) * density(11) 
  pd(20,21) = pd(20,21) + rrt(462) * density(23) 
  pd(20,23) = pd(20,23) + rrt(462) * density(21) 
  pd(21,21) = pd(21,21) - rrt(462) * density(23) 
  pd(21,23) = pd(21,23) - rrt(462) * density(21) 
  pd(23,21) = pd(23,21) - rrt(462) * density(23) 
  pd(23,23) = pd(23,23) - rrt(462) * density(21) 
  pd(20,21) = pd(20,21) + rrt(463) * density(31) 
  pd(20,31) = pd(20,31) + rrt(463) * density(21) 
  pd(21,21) = pd(21,21) - rrt(463) * density(31) 
  pd(21,31) = pd(21,31) - rrt(463) * density(21) 
  pd(31,21) = pd(31,21) - rrt(463) * density(31) 
  pd(31,31) = pd(31,31) - rrt(463) * density(21) 
  pd(37,21) = pd(37,21) + rrt(463) * density(31) 
  pd(37,31) = pd(37,31) + rrt(463) * density(21) 
  pd(18,20) = pd(18,20) + rrt(464) * density(45) 
  pd(18,45) = pd(18,45) + rrt(464) * density(20) 
  pd(20,20) = pd(20,20) - rrt(464) * density(45) 
  pd(20,45) = pd(20,45) - rrt(464) * density(20) 
  pd(37,20) = pd(37,20) + rrt(464) * density(45) 
  pd(37,45) = pd(37,45) + rrt(464) * density(20) 
  pd(45,20) = pd(45,20) - rrt(464) * density(45) 
  pd(45,45) = pd(45,45) - rrt(464) * density(20) 
  pd(11,13) = pd(11,13) + rrt(465) * density(45) 
  pd(11,45) = pd(11,45) + rrt(465) * density(13) 
  pd(13,13) = pd(13,13) - rrt(465) * density(45) 
  pd(13,45) = pd(13,45) - rrt(465) * density(13) 
  pd(18,13) = pd(18,13) + rrt(465) * density(45) 
  pd(18,45) = pd(18,45) + rrt(465) * density(13) 
  pd(45,13) = pd(45,13) - rrt(465) * density(45) 
  pd(45,45) = pd(45,45) - rrt(465) * density(13) 
  pd(23,45) = pd(23,45) + rrt(466) 
  pd(34,45) = pd(34,45) + rrt(466) 
  pd(45,45) = pd(45,45) - rrt(466) 
  pd(18,31) = pd(18,31) + rrt(467) * density(45) 
  pd(18,45) = pd(18,45) + rrt(467) * density(31) 
  pd(31,31) = pd(31,31) - rrt(467) * density(45) 
  pd(31,45) = pd(31,45) - rrt(467) * density(31) 
  pd(45,31) = pd(45,31) - rrt(467) * density(45) 
  pd(45,45) = pd(45,45) - rrt(467) * density(31) 
  pd(49,31) = pd(49,31) + rrt(467) * density(45) 
  pd(49,45) = pd(49,45) + rrt(467) * density(31) 
  pd(28,28) = pd(28,28) - rrt(468) * density(45) 
  pd(28,45) = pd(28,45) - rrt(468) * density(28) 
  pd(45,28) = pd(45,28) - rrt(468) * density(45) 
  pd(45,45) = pd(45,45) - rrt(468) * density(28) 
  pd(49,28) = pd(49,28) + rrt(468) * density(45) 
  pd(49,45) = pd(49,45) + rrt(468) * density(28) 
  pd(18,18) = pd(18,18) - rrt(469) * density(18) * 4.0d0
  pd(37,18) = pd(37,18) + rrt(469) * density(18) * 2.0d0
  pd(45,18) = pd(45,18) + rrt(469) * density(18) * 2.0d0
  pd(18,18) = pd(18,18) - rrt(470) * density(37) 
  pd(18,37) = pd(18,37) - rrt(470) * density(18) 
  pd(20,18) = pd(20,18) + rrt(470) * density(37) 
  pd(20,37) = pd(20,37) + rrt(470) * density(18) 
  pd(37,18) = pd(37,18) - rrt(470) * density(37) 
  pd(37,37) = pd(37,37) - rrt(470) * density(18) 
  pd(45,18) = pd(45,18) + rrt(470) * density(37) 
  pd(45,37) = pd(45,37) + rrt(470) * density(18) 
  pd(18,18) = pd(18,18) - rrt(471) * density(20) 
  pd(18,20) = pd(18,20) - rrt(471) * density(18) 
  pd(20,18) = pd(20,18) - rrt(471) * density(20) 
  pd(20,20) = pd(20,20) - rrt(471) * density(18) 
  pd(37,18) = pd(37,18) + rrt(471) * density(20) * 2.0d0
  pd(37,20) = pd(37,20) + rrt(471) * density(18) * 2.0d0
  pd(11,11) = pd(11,11) - rrt(472) * density(18) 
  pd(11,18) = pd(11,18) - rrt(472) * density(11) 
  pd(13,11) = pd(13,11) + rrt(472) * density(18) 
  pd(13,18) = pd(13,18) + rrt(472) * density(11) 
  pd(18,11) = pd(18,11) - rrt(472) * density(18) 
  pd(18,18) = pd(18,18) - rrt(472) * density(11) 
  pd(45,11) = pd(45,11) + rrt(472) * density(18) 
  pd(45,18) = pd(45,18) + rrt(472) * density(11) 
  pd(11,13) = pd(11,13) + rrt(473) * density(18) 
  pd(11,18) = pd(11,18) + rrt(473) * density(13) 
  pd(13,13) = pd(13,13) - rrt(473) * density(18) 
  pd(13,18) = pd(13,18) - rrt(473) * density(13) 
  pd(18,13) = pd(18,13) - rrt(473) * density(18) 
  pd(18,18) = pd(18,18) - rrt(473) * density(13) 
  pd(37,13) = pd(37,13) + rrt(473) * density(18) 
  pd(37,18) = pd(37,18) + rrt(473) * density(13) 
  pd(13,13) = pd(13,13) - rrt(474) * density(18) 
  pd(13,18) = pd(13,18) - rrt(474) * density(13) 
  pd(18,13) = pd(18,13) - rrt(474) * density(18) 
  pd(18,18) = pd(18,18) - rrt(474) * density(13) 
  pd(45,13) = pd(45,13) + rrt(474) * density(18) 
  pd(45,18) = pd(45,18) + rrt(474) * density(13) 
  pd(13,13) = pd(13,13) - rrt(475) * density(18) 
  pd(13,18) = pd(13,18) - rrt(475) * density(13) 
  pd(18,13) = pd(18,13) - rrt(475) * density(18) 
  pd(18,18) = pd(18,18) - rrt(475) * density(13) 
  pd(23,13) = pd(23,13) + rrt(475) * density(18) 
  pd(23,18) = pd(23,18) + rrt(475) * density(13) 
  pd(34,13) = pd(34,13) + rrt(475) * density(18) 
  pd(34,18) = pd(34,18) + rrt(475) * density(13) 
  pd(13,18) = pd(13,18) + rrt(476) 
  pd(18,18) = pd(18,18) - rrt(476) 
  pd(37,18) = pd(37,18) + rrt(476) 
  pd(18,18) = pd(18,18) - rrt(477) 
  pd(23,18) = pd(23,18) + rrt(477) 
  pd(42,18) = pd(42,18) + rrt(477) 
  pd(18,18) = pd(18,18) - rrt(478) * density(31) 
  pd(18,31) = pd(18,31) - rrt(478) * density(18) 
  pd(31,18) = pd(31,18) - rrt(478) * density(31) 
  pd(31,31) = pd(31,31) - rrt(478) * density(18) 
  pd(37,18) = pd(37,18) + rrt(478) * density(31) 
  pd(37,31) = pd(37,31) + rrt(478) * density(18) 
  pd(49,18) = pd(49,18) + rrt(478) * density(31) 
  pd(49,31) = pd(49,31) + rrt(478) * density(18) 
  pd(20,21) = pd(20,21) + rrt(479) * density(37) 
  pd(20,37) = pd(20,37) + rrt(479) * density(21) 
  pd(21,21) = pd(21,21) - rrt(479) * density(37) 
  pd(21,37) = pd(21,37) - rrt(479) * density(21) 
  pd(29,21) = pd(29,21) + rrt(479) * density(37) 
  pd(29,37) = pd(29,37) + rrt(479) * density(21) 
  pd(37,21) = pd(37,21) - rrt(479) * density(37) 
  pd(37,37) = pd(37,37) - rrt(479) * density(21) 
  pd(18,37) = pd(18,37) + rrt(480) * density(37) * 2.0d0
  pd(20,37) = pd(20,37) + rrt(480) * density(37) * 2.0d0
  pd(37,37) = pd(37,37) - rrt(480) * density(37) * 4.0d0
  pd(13,37) = pd(13,37) + rrt(481) 
  pd(20,37) = pd(20,37) + rrt(481) 
  pd(37,37) = pd(37,37) - rrt(481) 
  pd(11,13) = pd(11,13) + rrt(482) * density(37) 
  pd(11,37) = pd(11,37) + rrt(482) * density(13) 
  pd(13,13) = pd(13,13) - rrt(482) * density(37) 
  pd(13,37) = pd(13,37) - rrt(482) * density(13) 
  pd(20,13) = pd(20,13) + rrt(482) * density(37) 
  pd(20,37) = pd(20,37) + rrt(482) * density(13) 
  pd(37,13) = pd(37,13) - rrt(482) * density(37) 
  pd(37,37) = pd(37,37) - rrt(482) * density(13) 
  pd(13,13) = pd(13,13) - rrt(483) * density(37) 
  pd(13,37) = pd(13,37) - rrt(483) * density(13) 
  pd(18,13) = pd(18,13) + rrt(483) * density(37) 
  pd(18,37) = pd(18,37) + rrt(483) * density(13) 
  pd(37,13) = pd(37,13) - rrt(483) * density(37) 
  pd(37,37) = pd(37,37) - rrt(483) * density(13) 
  pd(23,37) = pd(23,37) + rrt(484) 
  pd(29,37) = pd(29,37) + rrt(484) 
  pd(37,37) = pd(37,37) - rrt(484) 
  pd(23,23) = pd(23,23) - rrt(485) * density(37) 
  pd(23,37) = pd(23,37) - rrt(485) * density(23) 
  pd(31,23) = pd(31,23) + rrt(485) * density(37) 
  pd(31,37) = pd(31,37) + rrt(485) * density(23) 
  pd(37,23) = pd(37,23) - rrt(485) * density(37) 
  pd(37,37) = pd(37,37) - rrt(485) * density(23) 
  pd(20,31) = pd(20,31) + rrt(486) * density(37) 
  pd(20,37) = pd(20,37) + rrt(486) * density(31) 
  pd(31,31) = pd(31,31) - rrt(486) * density(37) 
  pd(31,37) = pd(31,37) - rrt(486) * density(31) 
  pd(37,31) = pd(37,31) - rrt(486) * density(37) 
  pd(37,37) = pd(37,37) - rrt(486) * density(31) 
  pd(49,31) = pd(49,31) + rrt(486) * density(37) 
  pd(49,37) = pd(49,37) + rrt(486) * density(31) 
  pd(11,11) = pd(11,11) - rrt(487) * density(20) 
  pd(11,20) = pd(11,20) - rrt(487) * density(11) 
  pd(13,11) = pd(13,11) + rrt(487) * density(20) 
  pd(13,20) = pd(13,20) + rrt(487) * density(11) 
  pd(20,11) = pd(20,11) - rrt(487) * density(20) 
  pd(20,20) = pd(20,20) - rrt(487) * density(11) 
  pd(37,11) = pd(37,11) + rrt(487) * density(20) 
  pd(37,20) = pd(37,20) + rrt(487) * density(11) 
  pd(13,13) = pd(13,13) - rrt(488) * density(20) 
  pd(13,20) = pd(13,20) - rrt(488) * density(13) 
  pd(20,13) = pd(20,13) - rrt(488) * density(20) 
  pd(20,20) = pd(20,20) - rrt(488) * density(13) 
  pd(37,13) = pd(37,13) + rrt(488) * density(20) 
  pd(37,20) = pd(37,20) + rrt(488) * density(13) 
  pd(20,20) = pd(20,20) - rrt(489) 
  pd(21,20) = pd(21,20) + rrt(489) 
  pd(23,20) = pd(23,20) + rrt(489) 
  pd(31,31) = pd(31,31) - rrt(490) 
  pd(34,31) = pd(34,31) + rrt(490) 
  pd(42,31) = pd(42,31) + rrt(490) 
  pd(18,28) = pd(18,28) + rrt(491) * density(31) 
  pd(18,31) = pd(18,31) + rrt(491) * density(28) 
  pd(28,28) = pd(28,28) - rrt(491) * density(31) 
  pd(28,31) = pd(28,31) - rrt(491) * density(28) 
  pd(31,28) = pd(31,28) - rrt(491) * density(31) 
  pd(31,31) = pd(31,31) - rrt(491) * density(28) 
  pd(42,28) = pd(42,28) + rrt(491) * density(31) 
  pd(42,31) = pd(42,31) + rrt(491) * density(28) 
  pd(23,31) = pd(23,31) + rrt(492) 
  pd(31,31) = pd(31,31) - rrt(492) 
  pd(37,31) = pd(37,31) + rrt(492) 
  pd(11,11) = pd(11,11) - rrt(493) * density(31) 
  pd(11,31) = pd(11,31) - rrt(493) * density(11) 
  pd(13,11) = pd(13,11) + rrt(493) * density(31) 
  pd(13,31) = pd(13,31) + rrt(493) * density(11) 
  pd(31,11) = pd(31,11) - rrt(493) * density(31) 
  pd(31,31) = pd(31,31) - rrt(493) * density(11) 
  pd(49,11) = pd(49,11) + rrt(493) * density(31) 
  pd(49,31) = pd(49,31) + rrt(493) * density(11) 
  pd(23,23) = pd(23,23) - rrt(494) * density(49) 
  pd(23,49) = pd(23,49) - rrt(494) * density(23) 
  pd(31,23) = pd(31,23) + rrt(494) * density(49) 
  pd(31,49) = pd(31,49) + rrt(494) * density(23) 
  pd(49,23) = pd(49,23) - rrt(494) * density(49) 
  pd(49,49) = pd(49,49) - rrt(494) * density(23) 
  pd(51,23) = pd(51,23) + rrt(494) * density(49) 
  pd(51,49) = pd(51,49) + rrt(494) * density(23) 
  pd(18,49) = pd(18,49) + rrt(495) 
  pd(23,49) = pd(23,49) + rrt(495) 
  pd(49,49) = pd(49,49) - rrt(495) 
  pd(34,49) = pd(34,49) + rrt(496) * 2.0d0
  pd(49,49) = pd(49,49) - rrt(496) 
  pd(11,13) = pd(11,13) + rrt(497) * density(49) 
  pd(11,49) = pd(11,49) + rrt(497) * density(13) 
  pd(13,13) = pd(13,13) - rrt(497) * density(49) 
  pd(13,49) = pd(13,49) - rrt(497) * density(13) 
  pd(31,13) = pd(31,13) + rrt(497) * density(49) 
  pd(31,49) = pd(31,49) + rrt(497) * density(13) 
  pd(49,13) = pd(49,13) - rrt(497) * density(49) 
  pd(49,49) = pd(49,49) - rrt(497) * density(13) 
  pd(23,28) = pd(23,28) + rrt(498) * density(49) 
  pd(23,49) = pd(23,49) + rrt(498) * density(28) 
  pd(28,28) = pd(28,28) - rrt(498) * density(49) 
  pd(28,49) = pd(28,49) - rrt(498) * density(28) 
  pd(31,28) = pd(31,28) + rrt(498) * density(49) 
  pd(31,49) = pd(31,49) + rrt(498) * density(28) 
  pd(49,28) = pd(49,28) - rrt(498) * density(49) 
  pd(49,49) = pd(49,49) - rrt(498) * density(28) 
  pd(29,29) = pd(29,29) - rrt(499) * density(49) 
  pd(29,49) = pd(29,49) - rrt(499) * density(29) 
  pd(31,29) = pd(31,29) + rrt(499) * density(49) 
  pd(31,49) = pd(31,49) + rrt(499) * density(29) 
  pd(42,29) = pd(42,29) + rrt(499) * density(49) 
  pd(42,49) = pd(42,49) + rrt(499) * density(29) 
  pd(49,29) = pd(49,29) - rrt(499) * density(49) 
  pd(49,49) = pd(49,49) - rrt(499) * density(29) 
  pd(18,18) = pd(18,18) - rrt(500) * density(49) 
  pd(18,49) = pd(18,49) - rrt(500) * density(18) 
  pd(31,18) = pd(31,18) + rrt(500) * density(49) 
  pd(31,49) = pd(31,49) + rrt(500) * density(18) 
  pd(45,18) = pd(45,18) + rrt(500) * density(49) 
  pd(45,49) = pd(45,49) + rrt(500) * density(18) 
  pd(49,18) = pd(49,18) - rrt(500) * density(49) 
  pd(49,49) = pd(49,49) - rrt(500) * density(18) 
  pd(06,06) = pd(06,06) - rrt(501) * density(49) 
  pd(06,49) = pd(06,49) - rrt(501) * density(06) 
  pd(21,06) = pd(21,06) + rrt(501) * density(49) 
  pd(21,49) = pd(21,49) + rrt(501) * density(06) 
  pd(31,06) = pd(31,06) + rrt(501) * density(49) 
  pd(31,49) = pd(31,49) + rrt(501) * density(06) 
  pd(49,06) = pd(49,06) - rrt(501) * density(49) 
  pd(49,49) = pd(49,49) - rrt(501) * density(06) 
  pd(31,34) = pd(31,34) + rrt(502) * density(49) 
  pd(31,49) = pd(31,49) + rrt(502) * density(34) 
  pd(34,34) = pd(34,34) - rrt(502) * density(49) 
  pd(34,49) = pd(34,49) - rrt(502) * density(34) 
  pd(43,34) = pd(43,34) + rrt(502) * density(49) 
  pd(43,49) = pd(43,49) + rrt(502) * density(34) 
  pd(49,34) = pd(49,34) - rrt(502) * density(49) 
  pd(49,49) = pd(49,49) - rrt(502) * density(34) 
  pd(20,20) = pd(20,20) - rrt(503) * density(49) 
  pd(20,49) = pd(20,49) - rrt(503) * density(20) 
  pd(31,20) = pd(31,20) + rrt(503) * density(49) 
  pd(31,49) = pd(31,49) + rrt(503) * density(20) 
  pd(37,20) = pd(37,20) + rrt(503) * density(49) 
  pd(37,49) = pd(37,49) + rrt(503) * density(20) 
  pd(49,20) = pd(49,20) - rrt(503) * density(49) 
  pd(49,49) = pd(49,49) - rrt(503) * density(20) 
  pd(28,28) = pd(28,28) - rrt(504) * density(49) 
  pd(28,49) = pd(28,49) - rrt(504) * density(28) 
  pd(48,28) = pd(48,28) + rrt(504) * density(49) 
  pd(48,49) = pd(48,49) + rrt(504) * density(28) 
  pd(49,28) = pd(49,28) - rrt(504) * density(49) 
  pd(49,49) = pd(49,49) - rrt(504) * density(28) 
  pd(23,48) = pd(23,48) + rrt(505) 
  pd(31,48) = pd(31,48) + rrt(505) 
  pd(48,48) = pd(48,48) - rrt(505) 
  pd(11,11) = pd(11,11) - rrt(506) 
  pd(13,11) = pd(13,11) + rrt(506) * 2.0d0
  pd(11,13) = pd(11,13) + rrt(507) * density(13) * 2.0d0
  pd(13,13) = pd(13,13) - rrt(507) * density(13) * 4.0d0
  if( ldensity_constant ) then
    do i = 1, species_max
      if( density_constant(i) ) pd(i,:) = 0.0d0
    enddo
  endif
  if( lgas_heating ) then
    pd(53,1) = eV_to_K * ZDPlasKin_cfg(11)
    pd(53,:) = pd(53,:) * ZDPlasKin_cfg(13)
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
  DOUBLE PRECISION, PARAMETER :: F0= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F1= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F2= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F3= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F4= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F5= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F6= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F7= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F8= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F9= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F10= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F11= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F12= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F13= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F14= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F15= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F16= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F17= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F18= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F19= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F20= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F21= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F22= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F23= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F24= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F25= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F26= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F27= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F28= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F29= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F30= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F31= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F32= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F33= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F34= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F35= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F36= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F37= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F38= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F39= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F40= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F41= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F42= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F43= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F44= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F45= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F46= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F47= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F48= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F49= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F50= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F51= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F52= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F53= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F54= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F55= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F56= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F57= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F58= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F59= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F60= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F61= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F62= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F63= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F64= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F65= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F66= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F67= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F68= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F69= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F70= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F71= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F72= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F73= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F74= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F75= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F76= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F77= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F78= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F79= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F80= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F81= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F82= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F83= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F84= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F85= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F86= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F87= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F88= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F89= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F90= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F91= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F92= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F93= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F94= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F95= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F96= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F97= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F98= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F99= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F100= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F101= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F102= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F103= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F104= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F105= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F106= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F107= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F108= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F109= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F110= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F111= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F112= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F113= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F114= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F115= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F116= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F117= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F118= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F119= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F120= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F121= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F122= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F123= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F124= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F125= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F126= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F127= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F128= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F129= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F130= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F131= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F132= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F133= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F134= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F135= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F136= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F137= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F138= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F139= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F140= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F141= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F142= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F143= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F144= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F145= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F146= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F147= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F148= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F149= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F150= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F151= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F152= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F153= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F154= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F155= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F156= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F157= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F158= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F159= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F160= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F161= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F162= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F163= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F164= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F165= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F166= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F167= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F168= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F169= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F170= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F171= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F172= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F173= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F174= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F175= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F176= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F177= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F178= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F179= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F180= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F181= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F182= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F183= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F184= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F185= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F186= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F187= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F188= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F189= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F190= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F191= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F192= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F193= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F194= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F195= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F196= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F197= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F198= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F199= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F200= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F201= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F202= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F203= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F204= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F205= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F206= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F207= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F208= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F209= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F210= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F211= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F212= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F213= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F214= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F215= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F216= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F217= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F218= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F219= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F220= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F221= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F222= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F223= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F224= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F225= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F226= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F227= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F228= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F229= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F230= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F231= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F232= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F233= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F234= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F235= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F236= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F237= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F238= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F239= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F240= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F241= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F242= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F243= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F244= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F245= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F246= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F247= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F248= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F249= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F250= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F251= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F252= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F253= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F254= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F255= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F256= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F257= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F258= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F259= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F260= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F261= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F262= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F263= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F264= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F265= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F266= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F267= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F268= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F269= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F270= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F271= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F272= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F273= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F274= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F275= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F276= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F277= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F278= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F279= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F280= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F281= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F282= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F283= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F284= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F285= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F286= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F287= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F288= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F289= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F290= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F291= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F292= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F293= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F294= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F295= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F296= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F297= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F298= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F299= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F300= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F301= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F302= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F303= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F304= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F305= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F306= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F307= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F308= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F309= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F310= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F311= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F312= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F313= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F314= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F315= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F316= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F317= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F318= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F319= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F320= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F321= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F322= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F323= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F324= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F325= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F326= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F327= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F328= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F329= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F330= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F331= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F332= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F333= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F334= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F335= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F336= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F337= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F338= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F339= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F340= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F341= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F342= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F343= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F344= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F345= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F346= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F347= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F348= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F349= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F350= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F351= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F352= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F353= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F354= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F355= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F356= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F357= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F358= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F359= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F360= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F361= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F362= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F363= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F364= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F365= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F366= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F367= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F368= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F369= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F370= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F371= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F372= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F373= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F374= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F375= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F376= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F377= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F378= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F379= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F380= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F381= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F382= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F383= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F384= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F385= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F386= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F387= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F388= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F389= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F390= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F391= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F392= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F393= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F394= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F395= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F396= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F397= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F398= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F399= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F400= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F401= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F402= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F403= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F404= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F405= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F406= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F407= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F408= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F409= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F410= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F411= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F412= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F413= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F414= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F415= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F416= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F417= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F418= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F419= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F420= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F421= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F422= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F423= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F424= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F425= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F426= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F427= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F428= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F429= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F430= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F431= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F432= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F433= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F434= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F435= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F436= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F437= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F438= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F439= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F440= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F441= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F442= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F443= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F444= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F445= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F446= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F447= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F448= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F449= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F450= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F451= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F452= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F453= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F454= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F455= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F456= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F457= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F458= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F459= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F460= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F461= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F462= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F463= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F464= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F465= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F466= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F467= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F468= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F469= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F470= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F471= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F472= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F473= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F474= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F475= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F476= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F477= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F478= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F479= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F480= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F481= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F482= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F483= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F484= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F485= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F486= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F487= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F488= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F489= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F490= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F491= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F492= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F493= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F494= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F495= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F496= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F497= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F498= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F499= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F500= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F501= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F502= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F503= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F504= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F505= 1.000D-3
  DOUBLE PRECISION, PARAMETER :: F506= 1.000D-3
  call ZDPlasKin_bolsig_rates()
  Tgas = ZDPlasKin_cfg(1)
  Te  = ZDPlasKin_cfg(4)
  rrt(001) = F0*bolsig_rates(bolsig_pointer(1))
  rrt(002) = F1*bolsig_rates(bolsig_pointer(2))
  rrt(003) = F2*bolsig_rates(bolsig_pointer(3))
  rrt(004) = F3*bolsig_rates(bolsig_pointer(4))
  rrt(005) = F4*bolsig_rates(bolsig_pointer(5))
  rrt(006) = F5*bolsig_rates(bolsig_pointer(6))
  rrt(007) = F6*bolsig_rates(bolsig_pointer(7))
  rrt(008) = F7*bolsig_rates(bolsig_pointer(8))
  rrt(009) = F8*bolsig_rates(bolsig_pointer(9))
  rrt(010) = F9*bolsig_rates(bolsig_pointer(10))
  rrt(011) = F10*bolsig_rates(bolsig_pointer(11))
  rrt(012) = F11*bolsig_rates(bolsig_pointer(12))
  rrt(013) = F12*bolsig_rates(bolsig_pointer(13))
  rrt(014) = F13*bolsig_rates(bolsig_pointer(14))
  rrt(015) = F14*bolsig_rates(bolsig_pointer(15))
  rrt(016) = F15*bolsig_rates(bolsig_pointer(16))
  rrt(017) = F16*bolsig_rates(bolsig_pointer(17))
  rrt(018) = F17*bolsig_rates(bolsig_pointer(18))
  rrt(019) = F18*bolsig_rates(bolsig_pointer(19))
  rrt(020) = F19*bolsig_rates(bolsig_pointer(20))
  rrt(021) = F20*bolsig_rates(bolsig_pointer(21))
  rrt(022) = F21*bolsig_rates(bolsig_pointer(22))
  rrt(023) = F22*bolsig_rates(bolsig_pointer(23))
  rrt(024) = F23*bolsig_rates(bolsig_pointer(24))
  rrt(025) = F24*bolsig_rates(bolsig_pointer(25))
  rrt(026) = F25*bolsig_rates(bolsig_pointer(26))
  rrt(027) = F26*bolsig_rates(bolsig_pointer(27))
  rrt(028) = F27*bolsig_rates(bolsig_pointer(28))
  rrt(029) = F28*bolsig_rates(bolsig_pointer(29))
  rrt(030) = F29*bolsig_rates(bolsig_pointer(30))
  rrt(031) = F30*bolsig_rates(bolsig_pointer(31))
  rrt(032) = F31*bolsig_rates(bolsig_pointer(32))
  rrt(033) = F32*bolsig_rates(bolsig_pointer(33))
  rrt(034) = F33*bolsig_rates(bolsig_pointer(34))
  rrt(035) = F34*bolsig_rates(bolsig_pointer(35))
  rrt(036) = F35*bolsig_rates(bolsig_pointer(36))
  rrt(037) = F36*bolsig_rates(bolsig_pointer(37))
  rrt(038) = F37*bolsig_rates(bolsig_pointer(38))
  rrt(039) = F38*bolsig_rates(bolsig_pointer(39))
  rrt(040) = F39*bolsig_rates(bolsig_pointer(40))
  rrt(041) = F40*bolsig_rates(bolsig_pointer(41))
  rrt(042) = F41*bolsig_rates(bolsig_pointer(42))
  rrt(043) = F42*bolsig_rates(bolsig_pointer(43))
  rrt(044) = F43*bolsig_rates(bolsig_pointer(44))
  rrt(045) = F44*bolsig_rates(bolsig_pointer(45))
  rrt(046) = F45*bolsig_rates(bolsig_pointer(46))
  rrt(047) = F46*bolsig_rates(bolsig_pointer(47))
  rrt(048) = F47*bolsig_rates(bolsig_pointer(48))
  rrt(049) = F48*bolsig_rates(bolsig_pointer(49))
  rrt(050) = F49*bolsig_rates(bolsig_pointer(50))
  rrt(051) = F50*bolsig_rates(bolsig_pointer(51))
  rrt(052) = F51*bolsig_rates(bolsig_pointer(52))
  rrt(053) = F52*bolsig_rates(bolsig_pointer(53))
  rrt(054) = F53*bolsig_rates(bolsig_pointer(54))
  rrt(055) = F54*bolsig_rates(bolsig_pointer(55))
  rrt(056) = F55*bolsig_rates(bolsig_pointer(56))
  rrt(057) = F56*bolsig_rates(bolsig_pointer(57))
  rrt(058) = F57*bolsig_rates(bolsig_pointer(58))
  rrt(059) = F58*bolsig_rates(bolsig_pointer(59))
  rrt(060) = F59*bolsig_rates(bolsig_pointer(60))
  rrt(061) = F60*bolsig_rates(bolsig_pointer(61))
  rrt(062) = F61*bolsig_rates(bolsig_pointer(62))
  rrt(063) = F62*bolsig_rates(bolsig_pointer(63))
  rrt(064) = F63*bolsig_rates(bolsig_pointer(64))
  rrt(065) = F64*bolsig_rates(bolsig_pointer(65))
  rrt(066) = F65*bolsig_rates(bolsig_pointer(66))
  rrt(067) = F66*bolsig_rates(bolsig_pointer(67))
  rrt(068) = F67*bolsig_rates(bolsig_pointer(68))
  rrt(069) = F68*bolsig_rates(bolsig_pointer(69))
  rrt(070) = F69*bolsig_rates(bolsig_pointer(70))
  rrt(071) = F70*bolsig_rates(bolsig_pointer(71))
  rrt(072) = F71*bolsig_rates(bolsig_pointer(72))
  rrt(073) = F72*bolsig_rates(bolsig_pointer(73))
  rrt(074) = F73*bolsig_rates(bolsig_pointer(74))
  rrt(075) = F74*bolsig_rates(bolsig_pointer(75))
  rrt(076) = F75*bolsig_rates(bolsig_pointer(76))
  rrt(077) = F76*bolsig_rates(bolsig_pointer(77))
  rrt(078) = F77*bolsig_rates(bolsig_pointer(78))
  rrt(079) = F78*bolsig_rates(bolsig_pointer(79))
  rrt(080) = F79*bolsig_rates(bolsig_pointer(80))
  rrt(081) = F80*bolsig_rates(bolsig_pointer(81))
  rrt(082) = F81*bolsig_rates(bolsig_pointer(82))
  rrt(083) = F82*bolsig_rates(bolsig_pointer(83))
  rrt(084) = F83*bolsig_rates(bolsig_pointer(84))
  rrt(085) = F84*bolsig_rates(bolsig_pointer(85))
  rrt(086) = F85*bolsig_rates(bolsig_pointer(86))
  rrt(087) = F86*bolsig_rates(bolsig_pointer(87))
  rrt(088) = F87*bolsig_rates(bolsig_pointer(88))
  rrt(089) = F88*bolsig_rates(bolsig_pointer(89))
  rrt(090) = F89*bolsig_rates(bolsig_pointer(90))
  rrt(091) = F90*bolsig_rates(bolsig_pointer(91))
  rrt(092) = F91*bolsig_rates(bolsig_pointer(92))
  rrt(093) = F92*bolsig_rates(bolsig_pointer(93))
  rrt(094) = F93*bolsig_rates(bolsig_pointer(94))
  rrt(095) = F94*bolsig_rates(bolsig_pointer(95))
  rrt(096) = F95*bolsig_rates(bolsig_pointer(96))
  rrt(097) = F96*bolsig_rates(bolsig_pointer(97))
  rrt(098) = F97*bolsig_rates(bolsig_pointer(98))
  rrt(099) = F98*bolsig_rates(bolsig_pointer(99))
  rrt(100) = F99*bolsig_rates(bolsig_pointer(100))
  rrt(101) = F100*bolsig_rates(bolsig_pointer(101))
  rrt(102) = F101*bolsig_rates(bolsig_pointer(102))
  rrt(103) = F102*bolsig_rates(bolsig_pointer(103))
  rrt(104) = F103*bolsig_rates(bolsig_pointer(104))
  rrt(105) = F104*bolsig_rates(bolsig_pointer(105))
  rrt(106) = F105*bolsig_rates(bolsig_pointer(106))
  rrt(107) = F106*bolsig_rates(bolsig_pointer(107))
  rrt(108) = F107*bolsig_rates(bolsig_pointer(108))
  rrt(109) = F108*bolsig_rates(bolsig_pointer(109))
  rrt(110) = F109*bolsig_rates(bolsig_pointer(110))
  rrt(111) = F110*bolsig_rates(bolsig_pointer(111))
  rrt(112) = F111*bolsig_rates(bolsig_pointer(112))
  rrt(113) = F112*bolsig_rates(bolsig_pointer(113))
  rrt(114) = F113*bolsig_rates(bolsig_pointer(114))
  rrt(115) = F114*bolsig_rates(bolsig_pointer(115))
  rrt(116) = F115*bolsig_rates(bolsig_pointer(116))
  rrt(117) = F116*bolsig_rates(bolsig_pointer(117))
  rrt(118) = F117*bolsig_rates(bolsig_pointer(118))
  rrt(119) = F118*bolsig_rates(bolsig_pointer(119))
  rrt(120) = F119*bolsig_rates(bolsig_pointer(120))
  rrt(121) = F120*bolsig_rates(bolsig_pointer(121))
  rrt(122) = F121*bolsig_rates(bolsig_pointer(122))
  rrt(123) = F122*bolsig_rates(bolsig_pointer(123))
  rrt(124) = F123*bolsig_rates(bolsig_pointer(124))
  rrt(125) = F124*bolsig_rates(bolsig_pointer(125))
  rrt(126) = F125*bolsig_rates(bolsig_pointer(126))
  rrt(127) = F126*bolsig_rates(bolsig_pointer(127))
  rrt(128) = F127*bolsig_rates(bolsig_pointer(128))
  rrt(129) = F128*bolsig_rates(bolsig_pointer(129))
  rrt(130) = F129*bolsig_rates(bolsig_pointer(130))
  rrt(131) = F130*bolsig_rates(bolsig_pointer(131))
  rrt(132) = F131*bolsig_rates(bolsig_pointer(132))
  rrt(133) = F132*bolsig_rates(bolsig_pointer(133))
  rrt(134) = F133*bolsig_rates(bolsig_pointer(134))
  rrt(135) = F134*bolsig_rates(bolsig_pointer(135))
  rrt(136) = F135*bolsig_rates(bolsig_pointer(136))
  rrt(137) = F136*bolsig_rates(bolsig_pointer(137))
  rrt(138) = F137*bolsig_rates(bolsig_pointer(138))
  rrt(139) = F138*bolsig_rates(bolsig_pointer(139))
  rrt(140) = F139*bolsig_rates(bolsig_pointer(140))
  rrt(141) = F140*bolsig_rates(bolsig_pointer(141))
  rrt(142) = F141*bolsig_rates(bolsig_pointer(142))
  rrt(143) = F142*bolsig_rates(bolsig_pointer(143))
  rrt(144) = F143*bolsig_rates(bolsig_pointer(144))
  rrt(145) = F144*bolsig_rates(bolsig_pointer(145))
  rrt(146) = F145*bolsig_rates(bolsig_pointer(146))
  rrt(147) = F146*bolsig_rates(bolsig_pointer(147))
  rrt(148) = F147*bolsig_rates(bolsig_pointer(148))
  rrt(149) = F148*bolsig_rates(bolsig_pointer(149))
  rrt(150) = F149*bolsig_rates(bolsig_pointer(150))
  rrt(151) = F150*bolsig_rates(bolsig_pointer(151))
  rrt(152) = F151*bolsig_rates(bolsig_pointer(152))
  rrt(153) = F152*bolsig_rates(bolsig_pointer(153))
  rrt(154) = F153*bolsig_rates(bolsig_pointer(154))
  rrt(155) = F154*bolsig_rates(bolsig_pointer(155))
  rrt(156) = F155*bolsig_rates(bolsig_pointer(156))
  rrt(157) = F156*bolsig_rates(bolsig_pointer(157))
  rrt(158) = F157*bolsig_rates(bolsig_pointer(158))
  rrt(159) = F158*bolsig_rates(bolsig_pointer(159))
  rrt(160) = F159*bolsig_rates(bolsig_pointer(160))
  rrt(161) = F160*bolsig_rates(bolsig_pointer(161))
  rrt(162) = F161*bolsig_rates(bolsig_pointer(162))
  rrt(163) = F162*bolsig_rates(bolsig_pointer(163))
  rrt(164) = F163*bolsig_rates(bolsig_pointer(164))
  rrt(165) = F164*bolsig_rates(bolsig_pointer(165))
  rrt(166) = F165*bolsig_rates(bolsig_pointer(166))
  rrt(167) = F166*bolsig_rates(bolsig_pointer(167))
  rrt(168) = F167*bolsig_rates(bolsig_pointer(168))
  rrt(169) = F168*bolsig_rates(bolsig_pointer(169))
  rrt(170) = F169*bolsig_rates(bolsig_pointer(170))
  rrt(171) = F170*bolsig_rates(bolsig_pointer(171))
  rrt(172) = F171*bolsig_rates(bolsig_pointer(172))
  rrt(173) = F172*bolsig_rates(bolsig_pointer(173))
  rrt(174) = F173*bolsig_rates(bolsig_pointer(174))
  rrt(175) = F174*bolsig_rates(bolsig_pointer(175))
  rrt(176) = F175*bolsig_rates(bolsig_pointer(176))
  rrt(177) = F176*bolsig_rates(bolsig_pointer(177))
  rrt(178) = F177*bolsig_rates(bolsig_pointer(178))
  rrt(179) = F178*bolsig_rates(bolsig_pointer(179))
  rrt(180) = F179*bolsig_rates(bolsig_pointer(180))
  rrt(181) = F180*bolsig_rates(bolsig_pointer(181))
  rrt(182) = F181*bolsig_rates(bolsig_pointer(182))
  rrt(183) = F182*bolsig_rates(bolsig_pointer(183))
  rrt(184) = F183*bolsig_rates(bolsig_pointer(184))
  rrt(185) = F184*bolsig_rates(bolsig_pointer(185))
  rrt(186) = F185*bolsig_rates(bolsig_pointer(186))
  rrt(187) = F186*bolsig_rates(bolsig_pointer(187))
  rrt(188) = F187*bolsig_rates(bolsig_pointer(188))
  rrt(189) = F188*bolsig_rates(bolsig_pointer(189))
  rrt(190) = F189*bolsig_rates(bolsig_pointer(190))
  rrt(191) = F190*bolsig_rates(bolsig_pointer(191))
  rrt(192) = F191*bolsig_rates(bolsig_pointer(192))
  rrt(193) = F192*bolsig_rates(bolsig_pointer(193))
  rrt(194) = F193*bolsig_rates(bolsig_pointer(194))
  rrt(195) = F194*bolsig_rates(bolsig_pointer(195))
  rrt(196) = F195*bolsig_rates(bolsig_pointer(196))
  rrt(197) = F196*bolsig_rates(bolsig_pointer(197))
  rrt(198) = F197*bolsig_rates(bolsig_pointer(198))
  rrt(199) = F198*bolsig_rates(bolsig_pointer(199))
  rrt(200) = F199*bolsig_rates(bolsig_pointer(200))
  rrt(201) = F200*bolsig_rates(bolsig_pointer(201))
  rrt(202) = F201*bolsig_rates(bolsig_pointer(202))
  rrt(203) = F202*bolsig_rates(bolsig_pointer(203))
  rrt(204) = F203*bolsig_rates(bolsig_pointer(204))
  rrt(205) = F204*bolsig_rates(bolsig_pointer(205))
  rrt(206) = F205*bolsig_rates(bolsig_pointer(206))
  rrt(207) = F206*bolsig_rates(bolsig_pointer(207))
  rrt(208) = F207*bolsig_rates(bolsig_pointer(208))
  rrt(209) = F208*bolsig_rates(bolsig_pointer(209))
  rrt(210) = F209*bolsig_rates(bolsig_pointer(210))
  rrt(211) = F210*bolsig_rates(bolsig_pointer(211))
  rrt(212) = F211*bolsig_rates(bolsig_pointer(212))
  rrt(213) = F212*bolsig_rates(bolsig_pointer(213))
  rrt(214) = F213*bolsig_rates(bolsig_pointer(214))
  rrt(215) = F214*bolsig_rates(bolsig_pointer(215))
  rrt(216) = F215*bolsig_rates(bolsig_pointer(216))
  rrt(217) = F216*bolsig_rates(bolsig_pointer(217))
  rrt(218) = F217*bolsig_rates(bolsig_pointer(218))
  rrt(219) = F218*bolsig_rates(bolsig_pointer(219))
  rrt(220) = F219*bolsig_rates(bolsig_pointer(220))
  rrt(221) = F220*bolsig_rates(bolsig_pointer(221))
  rrt(222) = F221*bolsig_rates(bolsig_pointer(222))
  rrt(223) = F222*bolsig_rates(bolsig_pointer(223))
  rrt(224) = F223*bolsig_rates(bolsig_pointer(224))
  rrt(225) = F224*bolsig_rates(bolsig_pointer(225))
  rrt(226) = F225*bolsig_rates(bolsig_pointer(226))
  rrt(227) = F226*bolsig_rates(bolsig_pointer(227))
  rrt(228) = F227*bolsig_rates(bolsig_pointer(228))
  rrt(229) = F228*bolsig_rates(bolsig_pointer(229))
  rrt(230) = F229*bolsig_rates(bolsig_pointer(230))
  rrt(231) = F230*bolsig_rates(bolsig_pointer(231))
  rrt(232) = F231*bolsig_rates(bolsig_pointer(232))
  rrt(233) = F232*bolsig_rates(bolsig_pointer(233))
  rrt(234) = F233*bolsig_rates(bolsig_pointer(234))
  rrt(235) = F234*bolsig_rates(bolsig_pointer(235))
  rrt(236) = F235*bolsig_rates(bolsig_pointer(236))
  rrt(237) = F236*bolsig_rates(bolsig_pointer(237))
  rrt(238) = F237*bolsig_rates(bolsig_pointer(238))
  rrt(239) = F238*bolsig_rates(bolsig_pointer(239))
  rrt(240) = F239
  rrt(241) = F240
  rrt(242) = F241
  rrt(243) = F242
  rrt(244) = F243
  rrt(245) = F244
  rrt(246) = F245
  rrt(247) = F246
  rrt(248) = F247
  rrt(249) = F248
  rrt(250) = F249
  rrt(251) = F250
  rrt(252) = F251
  rrt(253) = F252
  rrt(254) = F253
  rrt(255) = F254
  rrt(256) = F255*2.57D-07*(300./TGAS)**0.3
  rrt(257) = F256*6.61D-08*(300./TGAS)**0.3
  rrt(258) = F257*1.18D-08*(300./TGAS)**0.5
  rrt(259) = F258*2.42D-08*(300./TGAS)**0.5
  rrt(260) = F259*1.41D-08*(300./TGAS)**0.5
  rrt(261) = F260*2.25D-08*(300./TGAS)**0.5
  rrt(262) = F261*7.88D-09*(300./TGAS)**0.5
  rrt(263) = F262*1.00D-08*(300./TGAS)**0.5
  rrt(264) = F263*2.19D-08*(300./TGAS)**0.71
  rrt(265) = F264*3.36D-08*(300./TGAS)**0.71
  rrt(266) = F265*7.70D-09*(300./TGAS)**0.71
  rrt(267) = F266*1.92D-08*(300./TGAS)**0.71
  rrt(268) = F267*1.60D-08*(300./TGAS)**0.71
  rrt(269) = F268*8.98D-09*(300./TGAS)**0.71
  rrt(270) = F269*9.62D-09*(300./TGAS)**0.71
  rrt(271) = F270*8.29D-09*(300./TGAS)**0.71
  rrt(272) = F271*3.43D-08*(300./TGAS)**0.71
  rrt(273) = F272*1.34D-08*(300./TGAS)**0.71
  rrt(274) = F273*4.87D-09*(300./TGAS)**0.71
  rrt(275) = F274*3.17D21/(6.022D23*TE**4.5)
  rrt(276) = F275*3.17D21/(6.022D23*TE**4.5)
  rrt(277) = F276*2.11D-09
  rrt(278) = F277*9.60D-10
  rrt(279) = F278*6.90D-10
  rrt(280) = F279*2.25D-10
  rrt(281) = F280*1.50D-9
  rrt(282) = F281*1.60D-9
  rrt(283) = F282*1.50D-10
  rrt(284) = F283*1.50D-9
  rrt(285) = F284*1.91D-9
  rrt(286) = F285*4.23D-10
  rrt(287) = F286*1.38D-9
  rrt(288) = F287*1.23D-9
  rrt(289) = F288*1.13D-9
  rrt(290) = F289*3.30D-11
  rrt(291) = F290*1.00D-11
  rrt(292) = F291*1.36D-10
  rrt(293) = F292*1.20D-9
  rrt(294) = F293*9.90D-10
  rrt(295) = F294*7.10D-10
  rrt(296) = F295*1.48D-9
  rrt(297) = F296*3.50D-10
  rrt(298) = F297*3.00D-10
  rrt(299) = F298*1.38D-10
  rrt(300) = F299*3.60D-10
  rrt(301) = F300*8.40D-10
  rrt(302) = F301*2.31D-10
  rrt(303) = F302*3.97D-10
  rrt(304) = F303*1.60D-9
  rrt(305) = F304*6.50D-11
  rrt(306) = F305*1.09D-9
  rrt(307) = F306*1.43D-10
  rrt(308) = F307*1.20D-9
  rrt(309) = F308*1.15D-9
  rrt(310) = F309*2.47D-10
  rrt(311) = F310*1.00D-10
  rrt(312) = F311*1.00D-11
  rrt(313) = F312*5.00D-10
  rrt(314) = F313*5.00D-10
  rrt(315) = F314*3.00D-10
  rrt(316) = F315*2.91D-10
  rrt(317) = F316*8.90D-10
  rrt(318) = F317*3.30D-10
  rrt(319) = F318*6.80D-11
  rrt(320) = F319*4.10D-9
  rrt(321) = F320*1.31D-10
  rrt(322) = F321*2.48D-10
  rrt(323) = F322*4.14D-10
  rrt(324) = F323*3.30D-10
  rrt(325) = F324*1.00D-11
  rrt(326) = F325*3.74D-10
  rrt(327) = F326*2.40D-9
  rrt(328) = F327*2.10D-9
  rrt(329) = F328*1.70D-9
  rrt(330) = F329*1.20D-9
  rrt(331) = F330*2.40D-9
  rrt(332) = F331*1.40D-9
  rrt(333) = F332*1.15D-9
  rrt(334) = F333*1.15D-9
  rrt(335) = F334*2.00D-9
  rrt(336) = F335*1.70D-9
  rrt(337) = F336*3.50D-9
  rrt(338) = F337*1.14D-10
  rrt(339) = F338*1.40D-9
  rrt(340) = F339*2.30D-9
  rrt(341) = F340*1.00D-9
  rrt(342) = F341*1.00D-9
  rrt(343) = F342*7.10D-10
  rrt(344) = F343*7.10D-10
  rrt(345) = F344*2.94D-10
  rrt(346) = F345*1.37D-9
  rrt(347) = F346*2.35D-9
  rrt(348) = F347*6.86D-10
  rrt(349) = F348*1.96D-10
  rrt(350) = F349*2.21D-9
  rrt(351) = F350*1.81D-9
  rrt(352) = F351*8.82D-10
  rrt(353) = F352*4.80D-10
  rrt(354) = F353*4.82D-9
  rrt(355) = F354*2.10D-9
  rrt(356) = F355*6.39D-10
  rrt(357) = F356*1.50D-9
  rrt(358) = F357*2.30D-9
  rrt(359) = F358*3.40D-9
  rrt(360) = F359*1.40D-9
  rrt(361) = F360*1.40D-9
  rrt(362) = F361*1.90D-9
  rrt(363) = F362*1.30D-9
  rrt(364) = F363*1.40D-9
  rrt(365) = F364*2.80D-9
  rrt(366) = F365*1.65D-9
  rrt(367) = F366*3.06D-9
  rrt(368) = F367*1.00D-9
  rrt(369) = F368*3.00D-9
  rrt(370) = F369*1.00D-9
  rrt(371) = F370*2.00D-9
  rrt(372) = F371*2.00D-9
  rrt(373) = F372*5.40D-10
  rrt(374) = F373*4.08D-18*TGAS**2*EXP(-4163./TGAS)
  rrt(375) = F374*9.97D-11
  rrt(376) = F375*2.51D-15*(TGAS/298.)**4.14*EXP(-52.55/(R*TGAS))
  rrt(377) = F376*4.26D-15*(TGAS/298.)**4.02*EXP(-22.86/(R*TGAS))
  rrt(378) = F377*3.01D-12*EXP(-2.08/(R*TGAS))
  rrt(379) = F378*3.54D-16*(TGAS/298.)**4.02*EXP(-45.48/(R*TGAS))
  rrt(380) = F379*1.71D-14*(TGAS/298.)**3.40*EXP(-97.28/(R*TGAS))
  rrt(381) = F380*9.86D-13*(TGAS/298.)**3.00*EXP(-36.67/(R*TGAS))
  rrt(382) = F381*1.33D-10*EXP(-167.00/(R*TGAS))
  rrt(383) = F382*5.68D-17*(TGAS/298.)**3.72*EXP(-33.42/(R*TGAS))
  rrt(384) = F383*1.90D-12
  rrt(385) = F384*2.40D16*EXP(-52800./TGAS)
  rrt(386) = F385*1.69D-8*EXP(-379./(R*TGAS))
  rrt(387) = F386*6.97D-9*EXP(-345./(R*TGAS))
  rrt(388) = F387*3.00D-44*TGAS**9.10
  rrt(389) = F388*3.32D-10*EXP(-45.98/(R*TGAS))
  rrt(390) = F389*3.01D-11
  rrt(391) = F390*3.01D-11
  rrt(392) = F391*3.01D-11
  rrt(393) = F392*1.61D-15*(TGAS/298.)**3.65*EXP(-29.93/(R*TGAS))
  rrt(394) = F393*3.01D-11
  rrt(395) = F394*3.01D-11
  rrt(396) = F395*1.20D-12*EXP(-25.94/(R*TGAS))
  rrt(397) = F396*8.30D-19*TGAS**2.00*EXP(-3938.65/TGAS)
  rrt(398) = F397*1.00D-11*EXP(7.48/(R*TGAS))
  rrt(399) = F398*6.64D-9*EXP(-348./(R*TGAS))
  rrt(400) = F399*2.66D-9*EXP(-6011.07/TGAS)
  rrt(401) = F400*9.96D-10
  rrt(402) = F401*3.00D-11
  rrt(403) = F402*1.14D-29
  rrt(404) = F403*1.79D-10*EXP(-1565.17/TGAS)
  rrt(405) = F404*4.98D-11
  rrt(406) = F405*6.64D-11
  rrt(407) = F406*3.29D-12*TGAS**0.43*EXP(186.21/TGAS)
  rrt(408) = F407*8.30D-11
  rrt(409) = F408*1.46D-13*(TGAS/298.)**3.30*EXP(-43.90/(R*TGAS))
  rrt(410) = F409*1.19D-15*(TGAS/298.)**3.82*EXP(-37.83/(R*TGAS))
  rrt(411) = F410*5.71D-14*(TGAS/298.)**3.30*EXP(-83.06/(R*TGAS))
  rrt(412) = F411*1.23D-11*(TGAS/298.)**1.50*EXP(-31.01/(R*TGAS))
  rrt(413) = F412*8.97D-20*EXP(-48.64/(R*TGAS))
  rrt(414) = F413*1.80D21*TGAS**(-1.24)*EXP(-45700./TGAS)
  rrt(415) = F414*8.30D-13*EXP(-62.77/(R*TGAS))
  rrt(416) = F415*1.79D-10*EXP(132.36/TGAS)
  rrt(417) = F416*9.00D-33*TGAS**6.43
  rrt(418) = F417*2.41D-12
  rrt(419) = F418*5.83D-14*(TGAS/298.)**3.13*EXP(-75.33/(R*TGAS))
  rrt(420) = F419*4.50D-13*EXP(-98.11/(R*TGAS))
  rrt(421) = F420*3.01D-12
  rrt(422) = F421*1.61D-15*(TGAS/298.)**3.65*EXP(-38.25/(R*TGAS))
  rrt(423) = F422*1.91D-12
  rrt(424) = F423*2.41D-12
  rrt(425) = F424*1.69D-15*(TGAS/298.)**3.50*EXP(-35.34/(R*TGAS))
  rrt(426) = F425*4.12D-15*(TGAS/298.)**3.60*EXP(-27.77/(R*TGAS))
  rrt(427) = F426*5.99D-11
  rrt(428) = F427*3.32D-12
  rrt(429) = F428*8.65D-7*TGAS**(-0.99)*EXP(-795.17/TGAS)
  rrt(430) = F429*4.08D12*(TGAS/298.)**1.04*EXP(-154./(R*TGAS))
  rrt(431) = F430*9.55D-12
  rrt(432) = F431*1.40D-12
  rrt(433) = F432*4.42D-11
  rrt(434) = F433*9.00D-10*EXP(-62.36/(R*TGAS))
  rrt(435) = F434*9.68D-12*(TGAS/298.)**1.28*EXP(-5.40/(R*TGAS))
  rrt(436) = F435*9.55D-12
  rrt(437) = F436*4.30D-7*EXP(-404./(R*TGAS))
  rrt(438) = F437*9.60D-11*EXP(-216./(R*TGAS))
  rrt(439) = F438*4.00D-11*EXP(-286./(R*TGAS))
  rrt(440) = F439*1.00D-10*EXP(-316./(R*TGAS))
  rrt(441) = F440*8.00D-10*EXP(-299./(R*TGAS))
  rrt(442) = F441*4.23D-18*TGAS**1.60*EXP(-2868.65/TGAS)
  rrt(443) = F442*8.00D12*TGAS**0.44*EXP(-44675.39/TGAS)
  rrt(444) = F443*3.00D-14*(TGAS/298.)**2.48*EXP(-25.65/(R*TGAS))
  rrt(445) = F444*4.75D-16*EXP(-180./(R*TGAS))
  rrt(446) = F445*5.30D-12*EXP(-2660./TGAS)
  rrt(447) = F446*5.00D-14*EXP(-25.53/(R*TGAS))
  rrt(448) = F447*3.50D-11
  rrt(449) = F448*1.46D-13*(TGAS/298.)**3.30*EXP(-43.90/(R*TGAS))
  rrt(450) = F449*2.01D-12
  rrt(451) = F450*2.01D-12
  rrt(452) = F451*1.68D-15*(TGAS/298.)**3.50*EXP(-19.62/(R*TGAS))
  rrt(453) = F452*8.00D-12
  rrt(454) = F453*1.61D-13*(TGAS/298.)**2.63*EXP(-35.75/(R*TGAS))
  rrt(455) = F454*1.60D-10
  rrt(456) = F455*1.01D-11*TGAS**0.27*EXP(-140.92/TGAS)
  rrt(457) = F456*2.00D14*EXP(-20000./TGAS)
  rrt(458) = F457*1.40D-12
  rrt(459) = F458*9.30D-12*EXP(-1207.85/TGAS)
  rrt(460) = F459*5.00D-13*EXP(-163./(R*TGAS))
  rrt(461) = F460*4.00D-12*EXP(-272./(R*TGAS))
  rrt(462) = F461*1.58D-5*(TGAS/298.)**-8.58*EXP(-84.81/(R*TGAS))
  rrt(463) = F462*1.20D-12*EXP(-37.66/(R*TGAS))
  rrt(464) = F463*5.71D-14*(TGAS/298.)**3.30*EXP(-83.06/(R*TGAS))
  rrt(465) = F464*2.19D-180*TGAS**2.54*EXP(-3400.10/TGAS)
  rrt(466) = F465*1.10D17*EXP(-42470./TGAS)
  rrt(467) = F466*1.61D-15*(TGAS/298.)**3.50*EXP(-27.77/(R*TGAS))
  rrt(468) = F467*4.42D-12
  rrt(469) = F468*2.81D-12
  rrt(470) = F469*1.69D-15*(TGAS/298.)**3.50*EXP(-27.77/(R*TGAS))
  rrt(471) = F470*2.41D-12*EXP(-0.55/(R*TGAS))
  rrt(472) = F471*3.19D-14*(TGAS/298.)**2.84*EXP(-38.25/(R*TGAS))
  rrt(473) = F472*3.01D-12
  rrt(474) = F473*6.00D-11
  rrt(475) = F474*6.74D-18*TGAS**2.19*EXP(-447.91/TGAS)
  rrt(476) = F475*1.25D17*EXP(-237./(R*TGAS))
  rrt(477) = F476*16.0*(TGAS/298.)**-(10.00)*EXP(-150./(R*TGAS))
  rrt(478) = F477*2.41D-12
  rrt(479) = F478*6.71D-11*EXP(-196./(R*TGAS))
  rrt(480) = F479*4.20D-10*EXP(-231./(R*TGAS))
  rrt(481) = F480*2.50D15*EXP(-363./(R*TGAS))
  rrt(482) = F481*4.40D-13*(TGAS/298.)**2.50*EXP(-10.39/(R*TGAS))
  rrt(483) = F482*1.29D-11*(TGAS/298.)**0.51*EXP(-5.15/(R*TGAS))
  rrt(484) = F483*1.28D13*(TGAS/298.)**(-15.70)*EXP(-502./(R*TGAS))
  rrt(485) = F484*1.27D-14*(TGAS/298.)**2.67*EXP(-28.66/(R*TGAS))
  rrt(486) = F485*1.69D-15*(TGAS/298.)**3.50*EXP(-27.77/(R*TGAS))
  rrt(487) = F486*1.39D-13*(TGAS/298.)**2.38*EXP(-79.49/(R*TGAS))
  rrt(488) = F487*78.80*(TGAS/298.)**(-11.76)*EXP(-98.53/(R*TGAS))
  rrt(489) = F488*1.26D13*EXP(-140./(R*TGAS))
  rrt(490) = F489*1.07D2*(TGAS/298.)**(-11.90)*EXP(-135./(R*TGAS))
  rrt(491) = F490*3.01D-11
  rrt(492) = F491*7.71D13*(TGAS/298.)**0.77*EXP(-128./(R*TGAS))
  rrt(493) = F492*2.52D-14*(TGAS/298.)**2.72*EXP(-40.99/(R*TGAS))
  rrt(494) = F493*8.32D-13*EXP(-56.87/(R*TGAS))
  rrt(495) = F494*8.87D-7*EXP(-180./(R*TGAS))
  rrt(496) = F495*7.84D-6*EXP(-207./(R*TGAS))
  rrt(497) = F496*2.19D-10*EXP(-39.24/(R*TGAS))
  rrt(498) = F497*1.81D-12*EXP(-20.54/(R*TGAS))
  rrt(499) = F498*2.42D-15*(TGAS/298.)**3.65*EXP(-21.62/(R*TGAS))
  rrt(500) = F499*2.47D-15*(TGAS/298.)**3.65*EXP(-38.25/(R*TGAS))
  rrt(501) = F500*1.00D-11
  rrt(502) = F501*2.47D-15*(TGAS/298.)**3.65*EXP(-38.25/(R*TGAS))
  rrt(503) = F502*8.63D-14*(TGAS/298.)**3.30*EXP(-83.06/(R*TGAS))
  rrt(504) = F503*9.61D-13
  rrt(505) = F504*3.16D16*EXP(-331./(R*TGAS))
  rrt(506) = F505*1.88D-8*(TGAS/298.)**(-1.10)*EXP(-437./(R*TGAS))
  rrt(507) = F506*5.52D-30*TGAS**(-1.00)
  where( lreaction_block(:) ) rrt(:) = 0.0d0
  return
end subroutine ZDPlasKin_reac_rates
!-----------------------------------------------------------------------------------------------------------------------------------
!
! END OF FILE
!
!-----------------------------------------------------------------------------------------------------------------------------------
