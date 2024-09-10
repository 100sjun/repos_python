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
! Tue Sep 10 10:52:44 2024
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
  /-1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0,&
    0, 0, 1, 1, 0, 1, 1, 0, 0, 1/
  data species_name(1:species_max) &
  /"E        ","C3H8(V2) ","C2H6(V24)","C2H3^+   ","C2H4(V1) ","CH4      ","C3H8     ","C2H5     ","C2H2^+   ","C3H6(V)  ",&
   "CH3      ","H^+      ","CH^+     ","C2H6     ","C2H4(V2) ","C2H4^+   ","C3H4     ","C3H7^+   ","CH4(V24) ","C3H4^+   ",&
   "C3H6     ","C5H12    ","H2       ","C3H8(V1) ","C4H9H    ","C3H5     ","C2H5^+   ","H3^+     ","C2H6(V13)","C2H4     ",&
   "C4H9     ","C2H2     ","C2H      ","H        ","CH5^+    ","C2H^+    ","C2H2(V2) ","C3H6^+   ","C3H5^+   ","C2H2(V13)",&
   "C3H8^+   ","C2H3     ","CH       ","C3H7     ","H2^+     ","CH2^+    ","CH2      ","CH4^+    ","CH3^+    ","C2H2(V5) ",&
   "CH4(V13) ","C2H6^+   "/
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
  /"C3H8(V2) ","C2H6(V24)","C2H4(V1) ","CH4      ","C3H8     ","C2H5     ","CH3      ","C2H6     ","C2H4(V2) ","C3H4     ",&
   "CH4(V24) ","C3H6     ","H2       ","C3H8(V1) ","C3H5     ","C2H6(V13)","C2H4     ","C2H2     ","C2H2(V2) ","C2H2(V13)",&
   "C2H3     ","CH       ","C3H7     ","CH2      ","C2H2(V5) ","CH4(V13) "/
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
  reac_source_local(06,001) = - reac_rate_local(001) 
  reac_source_local(19,001) = + reac_rate_local(001) 
  reac_source_local(06,002) = - reac_rate_local(002) 
  reac_source_local(51,002) = + reac_rate_local(002) 
  reac_source_local(03,003) = + reac_rate_local(003) 
  reac_source_local(14,003) = - reac_rate_local(003) 
  reac_source_local(14,004) = - reac_rate_local(004) 
  reac_source_local(29,004) = + reac_rate_local(004) 
  reac_source_local(05,005) = + reac_rate_local(005) 
  reac_source_local(30,005) = - reac_rate_local(005) 
  reac_source_local(15,006) = + reac_rate_local(006) 
  reac_source_local(30,006) = - reac_rate_local(006) 
  reac_source_local(32,007) = - reac_rate_local(007) 
  reac_source_local(50,007) = + reac_rate_local(007) 
  reac_source_local(32,008) = - reac_rate_local(008) 
  reac_source_local(37,008) = + reac_rate_local(008) 
  reac_source_local(32,009) = - reac_rate_local(009) 
  reac_source_local(40,009) = + reac_rate_local(009) 
  reac_source_local(07,010) = - reac_rate_local(010) 
  reac_source_local(24,010) = + reac_rate_local(010) 
  reac_source_local(02,011) = + reac_rate_local(011) 
  reac_source_local(07,011) = - reac_rate_local(011) 
  reac_source_local(10,012) = + reac_rate_local(012) 
  reac_source_local(21,012) = - reac_rate_local(012) 
  reac_source_local(06,013) = - reac_rate_local(013) 
  reac_source_local(11,013) = + reac_rate_local(013) 
  reac_source_local(34,013) = + reac_rate_local(013) 
  reac_source_local(06,014) = - reac_rate_local(014) 
  reac_source_local(23,014) = + reac_rate_local(014) 
  reac_source_local(47,014) = + reac_rate_local(014) 
  reac_source_local(06,015) = - reac_rate_local(015) 
  reac_source_local(23,015) = + reac_rate_local(015) 
  reac_source_local(34,015) = + reac_rate_local(015) 
  reac_source_local(43,015) = + reac_rate_local(015) 
  reac_source_local(11,016) = - reac_rate_local(016) 
  reac_source_local(34,016) = + reac_rate_local(016) 
  reac_source_local(47,016) = + reac_rate_local(016) 
  reac_source_local(11,017) = - reac_rate_local(017) 
  reac_source_local(23,017) = + reac_rate_local(017) 
  reac_source_local(43,017) = + reac_rate_local(017) 
  reac_source_local(34,018) = + reac_rate_local(018) 
  reac_source_local(43,018) = + reac_rate_local(018) 
  reac_source_local(47,018) = - reac_rate_local(018) 
  reac_source_local(01,019) = + reac_rate_local(019) 
  reac_source_local(06,019) = - reac_rate_local(019) 
  reac_source_local(48,019) = + reac_rate_local(019) 
  reac_source_local(01,020) = + reac_rate_local(020) 
  reac_source_local(06,020) = - reac_rate_local(020) 
  reac_source_local(34,020) = + reac_rate_local(020) 
  reac_source_local(49,020) = + reac_rate_local(020) 
  reac_source_local(01,021) = + reac_rate_local(021) 
  reac_source_local(06,021) = - reac_rate_local(021) 
  reac_source_local(23,021) = + reac_rate_local(021) 
  reac_source_local(46,021) = + reac_rate_local(021) 
  reac_source_local(01,022) = + reac_rate_local(022) 
  reac_source_local(06,022) = - reac_rate_local(022) 
  reac_source_local(13,022) = + reac_rate_local(022) 
  reac_source_local(23,022) = + reac_rate_local(022) 
  reac_source_local(34,022) = + reac_rate_local(022) 
  reac_source_local(01,023) = + reac_rate_local(023) 
  reac_source_local(11,023) = - reac_rate_local(023) 
  reac_source_local(49,023) = + reac_rate_local(023) 
  reac_source_local(01,024) = + reac_rate_local(024) 
  reac_source_local(11,024) = - reac_rate_local(024) 
  reac_source_local(34,024) = + reac_rate_local(024) 
  reac_source_local(46,024) = + reac_rate_local(024) 
  reac_source_local(01,025) = + reac_rate_local(025) 
  reac_source_local(11,025) = - reac_rate_local(025) 
  reac_source_local(13,025) = + reac_rate_local(025) 
  reac_source_local(23,025) = + reac_rate_local(025) 
  reac_source_local(01,026) = + reac_rate_local(026) 
  reac_source_local(46,026) = + reac_rate_local(026) 
  reac_source_local(47,026) = - reac_rate_local(026) 
  reac_source_local(01,027) = + reac_rate_local(027) 
  reac_source_local(13,027) = + reac_rate_local(027) 
  reac_source_local(34,027) = + reac_rate_local(027) 
  reac_source_local(47,027) = - reac_rate_local(027) 
  reac_source_local(01,028) = + reac_rate_local(028) 
  reac_source_local(13,028) = + reac_rate_local(028) 
  reac_source_local(43,028) = - reac_rate_local(028) 
  reac_source_local(08,029) = + reac_rate_local(029) 
  reac_source_local(14,029) = - reac_rate_local(029) 
  reac_source_local(34,029) = + reac_rate_local(029) 
  reac_source_local(14,030) = - reac_rate_local(030) 
  reac_source_local(23,030) = + reac_rate_local(030) 
  reac_source_local(30,030) = + reac_rate_local(030) 
  reac_source_local(14,031) = - reac_rate_local(031) 
  reac_source_local(23,031) = + reac_rate_local(031) 
  reac_source_local(34,031) = + reac_rate_local(031) 
  reac_source_local(42,031) = + reac_rate_local(031) 
  reac_source_local(14,032) = - reac_rate_local(032) 
  reac_source_local(23,032) = + reac_rate_local(032) * 2.d0
  reac_source_local(32,032) = + reac_rate_local(032) 
  reac_source_local(06,033) = + reac_rate_local(033) 
  reac_source_local(14,033) = - reac_rate_local(033) 
  reac_source_local(47,033) = + reac_rate_local(033) 
  reac_source_local(11,034) = + reac_rate_local(034) * 2.d0
  reac_source_local(14,034) = - reac_rate_local(034) 
  reac_source_local(08,035) = - reac_rate_local(035) 
  reac_source_local(30,035) = + reac_rate_local(035) 
  reac_source_local(34,035) = + reac_rate_local(035) 
  reac_source_local(08,036) = - reac_rate_local(036) 
  reac_source_local(23,036) = + reac_rate_local(036) 
  reac_source_local(42,036) = + reac_rate_local(036) 
  reac_source_local(08,037) = - reac_rate_local(037) 
  reac_source_local(34,037) = + reac_rate_local(037) * 2.d0
  reac_source_local(42,037) = + reac_rate_local(037) 
  reac_source_local(08,038) = - reac_rate_local(038) 
  reac_source_local(23,038) = + reac_rate_local(038) 
  reac_source_local(32,038) = + reac_rate_local(038) 
  reac_source_local(34,038) = + reac_rate_local(038) 
  reac_source_local(06,039) = + reac_rate_local(039) 
  reac_source_local(08,039) = - reac_rate_local(039) 
  reac_source_local(43,039) = + reac_rate_local(039) 
  reac_source_local(08,040) = - reac_rate_local(040) 
  reac_source_local(11,040) = + reac_rate_local(040) 
  reac_source_local(47,040) = + reac_rate_local(040) 
  reac_source_local(30,041) = - reac_rate_local(041) 
  reac_source_local(34,041) = + reac_rate_local(041) 
  reac_source_local(42,041) = + reac_rate_local(041) 
  reac_source_local(23,042) = + reac_rate_local(042) 
  reac_source_local(30,042) = - reac_rate_local(042) 
  reac_source_local(32,042) = + reac_rate_local(042) 
  reac_source_local(30,043) = - reac_rate_local(043) 
  reac_source_local(32,043) = + reac_rate_local(043) 
  reac_source_local(34,043) = + reac_rate_local(043) * 2.d0
  reac_source_local(11,044) = + reac_rate_local(044) 
  reac_source_local(30,044) = - reac_rate_local(044) 
  reac_source_local(43,044) = + reac_rate_local(044) 
  reac_source_local(30,045) = - reac_rate_local(045) 
  reac_source_local(47,045) = + reac_rate_local(045) * 2.d0
  reac_source_local(32,046) = + reac_rate_local(046) 
  reac_source_local(34,046) = + reac_rate_local(046) 
  reac_source_local(42,046) = - reac_rate_local(046) 
  reac_source_local(42,047) = - reac_rate_local(047) 
  reac_source_local(43,047) = + reac_rate_local(047) 
  reac_source_local(47,047) = + reac_rate_local(047) 
  reac_source_local(32,048) = - reac_rate_local(048) 
  reac_source_local(43,048) = + reac_rate_local(048) * 2.d0
  reac_source_local(07,049) = - reac_rate_local(049) 
  reac_source_local(34,049) = + reac_rate_local(049) 
  reac_source_local(44,049) = + reac_rate_local(049) 
  reac_source_local(07,050) = - reac_rate_local(050) 
  reac_source_local(21,050) = + reac_rate_local(050) 
  reac_source_local(23,050) = + reac_rate_local(050) 
  reac_source_local(07,051) = - reac_rate_local(051) 
  reac_source_local(17,051) = + reac_rate_local(051) 
  reac_source_local(23,051) = + reac_rate_local(051) * 2.d0
  reac_source_local(07,052) = - reac_rate_local(052) 
  reac_source_local(14,052) = + reac_rate_local(052) 
  reac_source_local(47,052) = + reac_rate_local(052) 
  reac_source_local(07,053) = - reac_rate_local(053) 
  reac_source_local(08,053) = + reac_rate_local(053) 
  reac_source_local(11,053) = + reac_rate_local(053) 
  reac_source_local(06,054) = + reac_rate_local(054) 
  reac_source_local(07,054) = - reac_rate_local(054) 
  reac_source_local(30,054) = + reac_rate_local(054) 
  reac_source_local(21,055) = + reac_rate_local(055) 
  reac_source_local(34,055) = + reac_rate_local(055) 
  reac_source_local(44,055) = - reac_rate_local(055) 
  reac_source_local(23,056) = + reac_rate_local(056) 
  reac_source_local(26,056) = + reac_rate_local(056) 
  reac_source_local(44,056) = - reac_rate_local(056) 
  reac_source_local(17,057) = + reac_rate_local(057) 
  reac_source_local(23,057) = + reac_rate_local(057) 
  reac_source_local(34,057) = + reac_rate_local(057) 
  reac_source_local(44,057) = - reac_rate_local(057) 
  reac_source_local(11,058) = + reac_rate_local(058) 
  reac_source_local(30,058) = + reac_rate_local(058) 
  reac_source_local(44,058) = - reac_rate_local(058) 
  reac_source_local(06,059) = + reac_rate_local(059) 
  reac_source_local(42,059) = + reac_rate_local(059) 
  reac_source_local(44,059) = - reac_rate_local(059) 
  reac_source_local(21,060) = - reac_rate_local(060) 
  reac_source_local(26,060) = + reac_rate_local(060) 
  reac_source_local(34,060) = + reac_rate_local(060) 
  reac_source_local(17,061) = + reac_rate_local(061) 
  reac_source_local(21,061) = - reac_rate_local(061) 
  reac_source_local(23,061) = + reac_rate_local(061) 
  reac_source_local(21,062) = - reac_rate_local(062) 
  reac_source_local(30,062) = + reac_rate_local(062) 
  reac_source_local(47,062) = + reac_rate_local(062) 
  reac_source_local(11,063) = + reac_rate_local(063) 
  reac_source_local(21,063) = - reac_rate_local(063) 
  reac_source_local(42,063) = + reac_rate_local(063) 
  reac_source_local(06,064) = + reac_rate_local(064) 
  reac_source_local(21,064) = - reac_rate_local(064) 
  reac_source_local(32,064) = + reac_rate_local(064) 
  reac_source_local(17,065) = + reac_rate_local(065) 
  reac_source_local(26,065) = - reac_rate_local(065) 
  reac_source_local(34,065) = + reac_rate_local(065) 
  reac_source_local(11,066) = + reac_rate_local(066) 
  reac_source_local(26,066) = - reac_rate_local(066) 
  reac_source_local(32,066) = + reac_rate_local(066) 
  reac_source_local(17,067) = - reac_rate_local(067) 
  reac_source_local(42,067) = + reac_rate_local(067) 
  reac_source_local(43,067) = + reac_rate_local(067) 
  reac_source_local(17,068) = - reac_rate_local(068) 
  reac_source_local(32,068) = + reac_rate_local(068) 
  reac_source_local(47,068) = + reac_rate_local(068) 
  reac_source_local(01,069) = + reac_rate_local(069) 
  reac_source_local(14,069) = - reac_rate_local(069) 
  reac_source_local(52,069) = + reac_rate_local(069) 
  reac_source_local(01,070) = + reac_rate_local(070) 
  reac_source_local(14,070) = - reac_rate_local(070) 
  reac_source_local(27,070) = + reac_rate_local(070) 
  reac_source_local(34,070) = + reac_rate_local(070) 
  reac_source_local(01,071) = + reac_rate_local(071) 
  reac_source_local(14,071) = - reac_rate_local(071) 
  reac_source_local(16,071) = + reac_rate_local(071) 
  reac_source_local(23,071) = + reac_rate_local(071) 
  reac_source_local(01,072) = + reac_rate_local(072) 
  reac_source_local(04,072) = + reac_rate_local(072) 
  reac_source_local(14,072) = - reac_rate_local(072) 
  reac_source_local(23,072) = + reac_rate_local(072) 
  reac_source_local(34,072) = + reac_rate_local(072) 
  reac_source_local(01,073) = + reac_rate_local(073) 
  reac_source_local(09,073) = + reac_rate_local(073) 
  reac_source_local(14,073) = - reac_rate_local(073) 
  reac_source_local(23,073) = + reac_rate_local(073) * 2.d0
  reac_source_local(01,074) = + reac_rate_local(074) 
  reac_source_local(11,074) = + reac_rate_local(074) 
  reac_source_local(14,074) = - reac_rate_local(074) 
  reac_source_local(49,074) = + reac_rate_local(074) 
  reac_source_local(01,075) = + reac_rate_local(075) 
  reac_source_local(06,075) = + reac_rate_local(075) 
  reac_source_local(14,075) = - reac_rate_local(075) 
  reac_source_local(46,075) = + reac_rate_local(075) 
  reac_source_local(01,076) = + reac_rate_local(076) 
  reac_source_local(08,076) = - reac_rate_local(076) 
  reac_source_local(27,076) = + reac_rate_local(076) 
  reac_source_local(01,077) = + reac_rate_local(077) 
  reac_source_local(08,077) = - reac_rate_local(077) 
  reac_source_local(16,077) = + reac_rate_local(077) 
  reac_source_local(34,077) = + reac_rate_local(077) 
  reac_source_local(01,078) = + reac_rate_local(078) 
  reac_source_local(04,078) = + reac_rate_local(078) 
  reac_source_local(08,078) = - reac_rate_local(078) 
  reac_source_local(23,078) = + reac_rate_local(078) 
  reac_source_local(01,079) = + reac_rate_local(079) 
  reac_source_local(08,079) = - reac_rate_local(079) 
  reac_source_local(09,079) = + reac_rate_local(079) 
  reac_source_local(23,079) = + reac_rate_local(079) 
  reac_source_local(34,079) = + reac_rate_local(079) 
  reac_source_local(01,080) = + reac_rate_local(080) 
  reac_source_local(08,080) = - reac_rate_local(080) 
  reac_source_local(47,080) = + reac_rate_local(080) 
  reac_source_local(49,080) = + reac_rate_local(080) 
  reac_source_local(01,081) = + reac_rate_local(081) 
  reac_source_local(08,081) = - reac_rate_local(081) 
  reac_source_local(11,081) = + reac_rate_local(081) 
  reac_source_local(46,081) = + reac_rate_local(081) 
  reac_source_local(01,082) = + reac_rate_local(082) 
  reac_source_local(06,082) = + reac_rate_local(082) 
  reac_source_local(08,082) = - reac_rate_local(082) 
  reac_source_local(13,082) = + reac_rate_local(082) 
  reac_source_local(01,083) = + reac_rate_local(083) 
  reac_source_local(16,083) = + reac_rate_local(083) 
  reac_source_local(30,083) = - reac_rate_local(083) 
  reac_source_local(01,084) = + reac_rate_local(084) 
  reac_source_local(04,084) = + reac_rate_local(084) 
  reac_source_local(30,084) = - reac_rate_local(084) 
  reac_source_local(34,084) = + reac_rate_local(084) 
  reac_source_local(01,085) = + reac_rate_local(085) 
  reac_source_local(30,085) = - reac_rate_local(085) 
  reac_source_local(43,085) = + reac_rate_local(085) 
  reac_source_local(49,085) = + reac_rate_local(085) 
  reac_source_local(01,086) = + reac_rate_local(086) 
  reac_source_local(30,086) = - reac_rate_local(086) 
  reac_source_local(46,086) = + reac_rate_local(086) 
  reac_source_local(47,086) = + reac_rate_local(086) 
  reac_source_local(01,087) = + reac_rate_local(087) 
  reac_source_local(11,087) = + reac_rate_local(087) 
  reac_source_local(13,087) = + reac_rate_local(087) 
  reac_source_local(30,087) = - reac_rate_local(087) 
  reac_source_local(01,088) = + reac_rate_local(088) 
  reac_source_local(04,088) = + reac_rate_local(088) 
  reac_source_local(42,088) = - reac_rate_local(088) 
  reac_source_local(01,089) = + reac_rate_local(089) 
  reac_source_local(09,089) = + reac_rate_local(089) 
  reac_source_local(34,089) = + reac_rate_local(089) 
  reac_source_local(42,089) = - reac_rate_local(089) 
  reac_source_local(01,090) = + reac_rate_local(090) 
  reac_source_local(42,090) = - reac_rate_local(090) 
  reac_source_local(43,090) = + reac_rate_local(090) 
  reac_source_local(46,090) = + reac_rate_local(090) 
  reac_source_local(01,091) = + reac_rate_local(091) 
  reac_source_local(13,091) = + reac_rate_local(091) 
  reac_source_local(42,091) = - reac_rate_local(091) 
  reac_source_local(47,091) = + reac_rate_local(091) 
  reac_source_local(01,092) = + reac_rate_local(092) 
  reac_source_local(12,092) = + reac_rate_local(092) 
  reac_source_local(32,092) = + reac_rate_local(092) 
  reac_source_local(42,092) = - reac_rate_local(092) 
  reac_source_local(01,093) = + reac_rate_local(093) 
  reac_source_local(09,093) = + reac_rate_local(093) 
  reac_source_local(32,093) = - reac_rate_local(093) 
  reac_source_local(01,094) = + reac_rate_local(094) 
  reac_source_local(13,094) = + reac_rate_local(094) 
  reac_source_local(32,094) = - reac_rate_local(094) 
  reac_source_local(43,094) = + reac_rate_local(094) 
  reac_source_local(01,095) = + reac_rate_local(095) 
  reac_source_local(07,095) = - reac_rate_local(095) 
  reac_source_local(41,095) = + reac_rate_local(095) 
  reac_source_local(01,096) = + reac_rate_local(096) 
  reac_source_local(07,096) = - reac_rate_local(096) 
  reac_source_local(18,096) = + reac_rate_local(096) 
  reac_source_local(34,096) = + reac_rate_local(096) 
  reac_source_local(01,097) = + reac_rate_local(097) 
  reac_source_local(07,097) = - reac_rate_local(097) 
  reac_source_local(23,097) = + reac_rate_local(097) 
  reac_source_local(38,097) = + reac_rate_local(097) 
  reac_source_local(01,098) = + reac_rate_local(098) 
  reac_source_local(07,098) = - reac_rate_local(098) 
  reac_source_local(23,098) = + reac_rate_local(098) 
  reac_source_local(34,098) = + reac_rate_local(098) 
  reac_source_local(39,098) = + reac_rate_local(098) 
  reac_source_local(01,099) = + reac_rate_local(099) 
  reac_source_local(07,099) = - reac_rate_local(099) 
  reac_source_local(20,099) = + reac_rate_local(099) 
  reac_source_local(23,099) = + reac_rate_local(099) * 2.d0
  reac_source_local(01,100) = + reac_rate_local(100) 
  reac_source_local(07,100) = - reac_rate_local(100) 
  reac_source_local(11,100) = + reac_rate_local(100) 
  reac_source_local(27,100) = + reac_rate_local(100) 
  reac_source_local(01,101) = + reac_rate_local(101) 
  reac_source_local(06,101) = + reac_rate_local(101) 
  reac_source_local(07,101) = - reac_rate_local(101) 
  reac_source_local(16,101) = + reac_rate_local(101) 
  reac_source_local(01,102) = + reac_rate_local(102) 
  reac_source_local(07,102) = - reac_rate_local(102) 
  reac_source_local(08,102) = + reac_rate_local(102) 
  reac_source_local(49,102) = + reac_rate_local(102) 
  reac_source_local(01,103) = + reac_rate_local(103) 
  reac_source_local(07,103) = - reac_rate_local(103) 
  reac_source_local(14,103) = + reac_rate_local(103) 
  reac_source_local(46,103) = + reac_rate_local(103) 
  reac_source_local(01,104) = + reac_rate_local(104) 
  reac_source_local(18,104) = + reac_rate_local(104) 
  reac_source_local(44,104) = - reac_rate_local(104) 
  reac_source_local(01,105) = + reac_rate_local(105) 
  reac_source_local(34,105) = + reac_rate_local(105) 
  reac_source_local(38,105) = + reac_rate_local(105) 
  reac_source_local(44,105) = - reac_rate_local(105) 
  reac_source_local(01,106) = + reac_rate_local(106) 
  reac_source_local(23,106) = + reac_rate_local(106) 
  reac_source_local(39,106) = + reac_rate_local(106) 
  reac_source_local(44,106) = - reac_rate_local(106) 
  reac_source_local(01,107) = + reac_rate_local(107) 
  reac_source_local(20,107) = + reac_rate_local(107) 
  reac_source_local(23,107) = + reac_rate_local(107) 
  reac_source_local(34,107) = + reac_rate_local(107) 
  reac_source_local(44,107) = - reac_rate_local(107) 
  reac_source_local(01,108) = + reac_rate_local(108) 
  reac_source_local(27,108) = + reac_rate_local(108) 
  reac_source_local(44,108) = - reac_rate_local(108) 
  reac_source_local(47,108) = + reac_rate_local(108) 
  reac_source_local(01,109) = + reac_rate_local(109) 
  reac_source_local(11,109) = + reac_rate_local(109) 
  reac_source_local(16,109) = + reac_rate_local(109) 
  reac_source_local(44,109) = - reac_rate_local(109) 
  reac_source_local(01,110) = + reac_rate_local(110) 
  reac_source_local(04,110) = + reac_rate_local(110) 
  reac_source_local(06,110) = + reac_rate_local(110) 
  reac_source_local(44,110) = - reac_rate_local(110) 
  reac_source_local(01,111) = + reac_rate_local(111) 
  reac_source_local(42,111) = + reac_rate_local(111) 
  reac_source_local(44,111) = - reac_rate_local(111) 
  reac_source_local(48,111) = + reac_rate_local(111) 
  reac_source_local(01,112) = + reac_rate_local(112) 
  reac_source_local(30,112) = + reac_rate_local(112) 
  reac_source_local(44,112) = - reac_rate_local(112) 
  reac_source_local(49,112) = + reac_rate_local(112) 
  reac_source_local(01,113) = + reac_rate_local(113) 
  reac_source_local(08,113) = + reac_rate_local(113) 
  reac_source_local(44,113) = - reac_rate_local(113) 
  reac_source_local(46,113) = + reac_rate_local(113) 
  reac_source_local(01,114) = + reac_rate_local(114) 
  reac_source_local(13,114) = + reac_rate_local(114) 
  reac_source_local(14,114) = + reac_rate_local(114) 
  reac_source_local(44,114) = - reac_rate_local(114) 
  reac_source_local(01,115) = + reac_rate_local(115) 
  reac_source_local(21,115) = - reac_rate_local(115) 
  reac_source_local(38,115) = + reac_rate_local(115) 
  reac_source_local(01,116) = + reac_rate_local(116) 
  reac_source_local(21,116) = - reac_rate_local(116) 
  reac_source_local(34,116) = + reac_rate_local(116) 
  reac_source_local(39,116) = + reac_rate_local(116) 
  reac_source_local(01,117) = + reac_rate_local(117) 
  reac_source_local(20,117) = + reac_rate_local(117) 
  reac_source_local(21,117) = - reac_rate_local(117) 
  reac_source_local(23,117) = + reac_rate_local(117) 
  reac_source_local(01,118) = + reac_rate_local(118) 
  reac_source_local(21,118) = - reac_rate_local(118) 
  reac_source_local(27,118) = + reac_rate_local(118) 
  reac_source_local(43,118) = + reac_rate_local(118) 
  reac_source_local(01,119) = + reac_rate_local(119) 
  reac_source_local(16,119) = + reac_rate_local(119) 
  reac_source_local(21,119) = - reac_rate_local(119) 
  reac_source_local(47,119) = + reac_rate_local(119) 
  reac_source_local(01,120) = + reac_rate_local(120) 
  reac_source_local(04,120) = + reac_rate_local(120) 
  reac_source_local(11,120) = + reac_rate_local(120) 
  reac_source_local(21,120) = - reac_rate_local(120) 
  reac_source_local(01,121) = + reac_rate_local(121) 
  reac_source_local(06,121) = + reac_rate_local(121) 
  reac_source_local(09,121) = + reac_rate_local(121) 
  reac_source_local(21,121) = - reac_rate_local(121) 
  reac_source_local(01,122) = + reac_rate_local(122) 
  reac_source_local(21,122) = - reac_rate_local(122) 
  reac_source_local(32,122) = + reac_rate_local(122) 
  reac_source_local(48,122) = + reac_rate_local(122) 
  reac_source_local(01,123) = + reac_rate_local(123) 
  reac_source_local(21,123) = - reac_rate_local(123) 
  reac_source_local(42,123) = + reac_rate_local(123) 
  reac_source_local(49,123) = + reac_rate_local(123) 
  reac_source_local(01,124) = + reac_rate_local(124) 
  reac_source_local(21,124) = - reac_rate_local(124) 
  reac_source_local(30,124) = + reac_rate_local(124) 
  reac_source_local(46,124) = + reac_rate_local(124) 
  reac_source_local(01,125) = + reac_rate_local(125) 
  reac_source_local(08,125) = + reac_rate_local(125) 
  reac_source_local(13,125) = + reac_rate_local(125) 
  reac_source_local(21,125) = - reac_rate_local(125) 
  reac_source_local(01,126) = + reac_rate_local(126) 
  reac_source_local(26,126) = - reac_rate_local(126) 
  reac_source_local(39,126) = + reac_rate_local(126) 
  reac_source_local(01,127) = + reac_rate_local(127) 
  reac_source_local(20,127) = + reac_rate_local(127) 
  reac_source_local(26,127) = - reac_rate_local(127) 
  reac_source_local(34,127) = + reac_rate_local(127) 
  reac_source_local(01,128) = + reac_rate_local(128) 
  reac_source_local(16,128) = + reac_rate_local(128) 
  reac_source_local(26,128) = - reac_rate_local(128) 
  reac_source_local(43,128) = + reac_rate_local(128) 
  reac_source_local(01,129) = + reac_rate_local(129) 
  reac_source_local(04,129) = + reac_rate_local(129) 
  reac_source_local(26,129) = - reac_rate_local(129) 
  reac_source_local(47,129) = + reac_rate_local(129) 
  reac_source_local(01,130) = + reac_rate_local(130) 
  reac_source_local(09,130) = + reac_rate_local(130) 
  reac_source_local(11,130) = + reac_rate_local(130) 
  reac_source_local(26,130) = - reac_rate_local(130) 
  reac_source_local(01,131) = + reac_rate_local(131) 
  reac_source_local(26,131) = - reac_rate_local(131) 
  reac_source_local(32,131) = + reac_rate_local(131) 
  reac_source_local(49,131) = + reac_rate_local(131) 
  reac_source_local(01,132) = + reac_rate_local(132) 
  reac_source_local(26,132) = - reac_rate_local(132) 
  reac_source_local(42,132) = + reac_rate_local(132) 
  reac_source_local(46,132) = + reac_rate_local(132) 
  reac_source_local(01,133) = + reac_rate_local(133) 
  reac_source_local(13,133) = + reac_rate_local(133) 
  reac_source_local(26,133) = - reac_rate_local(133) 
  reac_source_local(30,133) = + reac_rate_local(133) 
  reac_source_local(01,134) = + reac_rate_local(134) 
  reac_source_local(17,134) = - reac_rate_local(134) 
  reac_source_local(20,134) = + reac_rate_local(134) 
  reac_source_local(01,135) = + reac_rate_local(135) 
  reac_source_local(04,135) = + reac_rate_local(135) 
  reac_source_local(17,135) = - reac_rate_local(135) 
  reac_source_local(43,135) = + reac_rate_local(135) 
  reac_source_local(01,136) = + reac_rate_local(136) 
  reac_source_local(09,136) = + reac_rate_local(136) 
  reac_source_local(17,136) = - reac_rate_local(136) 
  reac_source_local(47,136) = + reac_rate_local(136) 
  reac_source_local(01,137) = + reac_rate_local(137) 
  reac_source_local(17,137) = - reac_rate_local(137) 
  reac_source_local(32,137) = + reac_rate_local(137) 
  reac_source_local(46,137) = + reac_rate_local(137) 
  reac_source_local(01,138) = + reac_rate_local(138) 
  reac_source_local(13,138) = + reac_rate_local(138) 
  reac_source_local(17,138) = - reac_rate_local(138) 
  reac_source_local(42,138) = + reac_rate_local(138) 
  reac_source_local(11,139) = + reac_rate_local(139) 
  reac_source_local(19,139) = - reac_rate_local(139) 
  reac_source_local(34,139) = + reac_rate_local(139) 
  reac_source_local(11,140) = + reac_rate_local(140) 
  reac_source_local(34,140) = + reac_rate_local(140) 
  reac_source_local(51,140) = - reac_rate_local(140) 
  reac_source_local(19,141) = - reac_rate_local(141) 
  reac_source_local(23,141) = + reac_rate_local(141) 
  reac_source_local(47,141) = + reac_rate_local(141) 
  reac_source_local(23,142) = + reac_rate_local(142) 
  reac_source_local(47,142) = + reac_rate_local(142) 
  reac_source_local(51,142) = - reac_rate_local(142) 
  reac_source_local(19,143) = - reac_rate_local(143) 
  reac_source_local(23,143) = + reac_rate_local(143) 
  reac_source_local(34,143) = + reac_rate_local(143) 
  reac_source_local(43,143) = + reac_rate_local(143) 
  reac_source_local(23,144) = + reac_rate_local(144) 
  reac_source_local(34,144) = + reac_rate_local(144) 
  reac_source_local(43,144) = + reac_rate_local(144) 
  reac_source_local(51,144) = - reac_rate_local(144) 
  reac_source_local(01,145) = + reac_rate_local(145) 
  reac_source_local(19,145) = - reac_rate_local(145) 
  reac_source_local(48,145) = + reac_rate_local(145) 
  reac_source_local(01,146) = + reac_rate_local(146) 
  reac_source_local(48,146) = + reac_rate_local(146) 
  reac_source_local(51,146) = - reac_rate_local(146) 
  reac_source_local(01,147) = + reac_rate_local(147) 
  reac_source_local(19,147) = - reac_rate_local(147) 
  reac_source_local(34,147) = + reac_rate_local(147) 
  reac_source_local(49,147) = + reac_rate_local(147) 
  reac_source_local(01,148) = + reac_rate_local(148) 
  reac_source_local(34,148) = + reac_rate_local(148) 
  reac_source_local(49,148) = + reac_rate_local(148) 
  reac_source_local(51,148) = - reac_rate_local(148) 
  reac_source_local(01,149) = + reac_rate_local(149) 
  reac_source_local(19,149) = - reac_rate_local(149) 
  reac_source_local(23,149) = + reac_rate_local(149) 
  reac_source_local(46,149) = + reac_rate_local(149) 
  reac_source_local(01,150) = + reac_rate_local(150) 
  reac_source_local(23,150) = + reac_rate_local(150) 
  reac_source_local(46,150) = + reac_rate_local(150) 
  reac_source_local(51,150) = - reac_rate_local(150) 
  reac_source_local(01,151) = + reac_rate_local(151) 
  reac_source_local(13,151) = + reac_rate_local(151) 
  reac_source_local(19,151) = - reac_rate_local(151) 
  reac_source_local(23,151) = + reac_rate_local(151) 
  reac_source_local(34,151) = + reac_rate_local(151) 
  reac_source_local(01,152) = + reac_rate_local(152) 
  reac_source_local(13,152) = + reac_rate_local(152) 
  reac_source_local(23,152) = + reac_rate_local(152) 
  reac_source_local(34,152) = + reac_rate_local(152) 
  reac_source_local(51,152) = - reac_rate_local(152) 
  reac_source_local(03,153) = - reac_rate_local(153) 
  reac_source_local(08,153) = + reac_rate_local(153) 
  reac_source_local(34,153) = + reac_rate_local(153) 
  reac_source_local(08,154) = + reac_rate_local(154) 
  reac_source_local(29,154) = - reac_rate_local(154) 
  reac_source_local(34,154) = + reac_rate_local(154) 
  reac_source_local(03,155) = - reac_rate_local(155) 
  reac_source_local(23,155) = + reac_rate_local(155) 
  reac_source_local(30,155) = + reac_rate_local(155) 
  reac_source_local(23,156) = + reac_rate_local(156) 
  reac_source_local(29,156) = - reac_rate_local(156) 
  reac_source_local(30,156) = + reac_rate_local(156) 
  reac_source_local(03,157) = - reac_rate_local(157) 
  reac_source_local(23,157) = + reac_rate_local(157) 
  reac_source_local(34,157) = + reac_rate_local(157) 
  reac_source_local(42,157) = + reac_rate_local(157) 
  reac_source_local(23,158) = + reac_rate_local(158) 
  reac_source_local(29,158) = - reac_rate_local(158) 
  reac_source_local(34,158) = + reac_rate_local(158) 
  reac_source_local(42,158) = + reac_rate_local(158) 
  reac_source_local(03,159) = - reac_rate_local(159) 
  reac_source_local(23,159) = + reac_rate_local(159) * 2.d0
  reac_source_local(32,159) = + reac_rate_local(159) 
  reac_source_local(23,160) = + reac_rate_local(160) * 2.d0
  reac_source_local(29,160) = - reac_rate_local(160) 
  reac_source_local(32,160) = + reac_rate_local(160) 
  reac_source_local(03,161) = - reac_rate_local(161) 
  reac_source_local(06,161) = + reac_rate_local(161) 
  reac_source_local(47,161) = + reac_rate_local(161) 
  reac_source_local(06,162) = + reac_rate_local(162) 
  reac_source_local(29,162) = - reac_rate_local(162) 
  reac_source_local(47,162) = + reac_rate_local(162) 
  reac_source_local(03,163) = - reac_rate_local(163) 
  reac_source_local(11,163) = + reac_rate_local(163) * 2.d0
  reac_source_local(11,164) = + reac_rate_local(164) * 2.d0
  reac_source_local(29,164) = - reac_rate_local(164) 
  reac_source_local(05,165) = - reac_rate_local(165) 
  reac_source_local(34,165) = + reac_rate_local(165) 
  reac_source_local(42,165) = + reac_rate_local(165) 
  reac_source_local(15,166) = - reac_rate_local(166) 
  reac_source_local(34,166) = + reac_rate_local(166) 
  reac_source_local(42,166) = + reac_rate_local(166) 
  reac_source_local(05,167) = - reac_rate_local(167) 
  reac_source_local(23,167) = + reac_rate_local(167) 
  reac_source_local(32,167) = + reac_rate_local(167) 
  reac_source_local(15,168) = - reac_rate_local(168) 
  reac_source_local(23,168) = + reac_rate_local(168) 
  reac_source_local(32,168) = + reac_rate_local(168) 
  reac_source_local(05,169) = - reac_rate_local(169) 
  reac_source_local(32,169) = + reac_rate_local(169) 
  reac_source_local(34,169) = + reac_rate_local(169) * 2.d0
  reac_source_local(15,170) = - reac_rate_local(170) 
  reac_source_local(32,170) = + reac_rate_local(170) 
  reac_source_local(34,170) = + reac_rate_local(170) * 2.d0
  reac_source_local(05,171) = - reac_rate_local(171) 
  reac_source_local(11,171) = + reac_rate_local(171) 
  reac_source_local(43,171) = + reac_rate_local(171) 
  reac_source_local(11,172) = + reac_rate_local(172) 
  reac_source_local(15,172) = - reac_rate_local(172) 
  reac_source_local(43,172) = + reac_rate_local(172) 
  reac_source_local(05,173) = - reac_rate_local(173) 
  reac_source_local(47,173) = + reac_rate_local(173) * 2.d0
  reac_source_local(15,174) = - reac_rate_local(174) 
  reac_source_local(47,174) = + reac_rate_local(174) * 2.d0
  reac_source_local(43,175) = + reac_rate_local(175) * 2.d0
  reac_source_local(50,175) = - reac_rate_local(175) 
  reac_source_local(37,176) = - reac_rate_local(176) 
  reac_source_local(43,176) = + reac_rate_local(176) * 2.d0
  reac_source_local(40,177) = - reac_rate_local(177) 
  reac_source_local(43,177) = + reac_rate_local(177) * 2.d0
  reac_source_local(24,178) = - reac_rate_local(178) 
  reac_source_local(34,178) = + reac_rate_local(178) 
  reac_source_local(44,178) = + reac_rate_local(178) 
  reac_source_local(02,179) = - reac_rate_local(179) 
  reac_source_local(34,179) = + reac_rate_local(179) 
  reac_source_local(44,179) = + reac_rate_local(179) 
  reac_source_local(21,180) = + reac_rate_local(180) 
  reac_source_local(23,180) = + reac_rate_local(180) 
  reac_source_local(24,180) = - reac_rate_local(180) 
  reac_source_local(02,181) = - reac_rate_local(181) 
  reac_source_local(21,181) = + reac_rate_local(181) 
  reac_source_local(23,181) = + reac_rate_local(181) 
  reac_source_local(17,182) = + reac_rate_local(182) 
  reac_source_local(23,182) = + reac_rate_local(182) * 2.d0
  reac_source_local(24,182) = - reac_rate_local(182) 
  reac_source_local(02,183) = - reac_rate_local(183) 
  reac_source_local(17,183) = + reac_rate_local(183) 
  reac_source_local(23,183) = + reac_rate_local(183) * 2.d0
  reac_source_local(14,184) = + reac_rate_local(184) 
  reac_source_local(24,184) = - reac_rate_local(184) 
  reac_source_local(47,184) = + reac_rate_local(184) 
  reac_source_local(02,185) = - reac_rate_local(185) 
  reac_source_local(14,185) = + reac_rate_local(185) 
  reac_source_local(47,185) = + reac_rate_local(185) 
  reac_source_local(08,186) = + reac_rate_local(186) 
  reac_source_local(11,186) = + reac_rate_local(186) 
  reac_source_local(24,186) = - reac_rate_local(186) 
  reac_source_local(02,187) = - reac_rate_local(187) 
  reac_source_local(08,187) = + reac_rate_local(187) 
  reac_source_local(11,187) = + reac_rate_local(187) 
  reac_source_local(06,188) = + reac_rate_local(188) 
  reac_source_local(24,188) = - reac_rate_local(188) 
  reac_source_local(30,188) = + reac_rate_local(188) 
  reac_source_local(02,189) = - reac_rate_local(189) 
  reac_source_local(06,189) = + reac_rate_local(189) 
  reac_source_local(30,189) = + reac_rate_local(189) 
  reac_source_local(01,190) = + reac_rate_local(190) 
  reac_source_local(03,190) = - reac_rate_local(190) 
  reac_source_local(52,190) = + reac_rate_local(190) 
  reac_source_local(01,191) = + reac_rate_local(191) 
  reac_source_local(29,191) = - reac_rate_local(191) 
  reac_source_local(52,191) = + reac_rate_local(191) 
  reac_source_local(01,192) = + reac_rate_local(192) 
  reac_source_local(03,192) = - reac_rate_local(192) 
  reac_source_local(27,192) = + reac_rate_local(192) 
  reac_source_local(34,192) = + reac_rate_local(192) 
  reac_source_local(01,193) = + reac_rate_local(193) 
  reac_source_local(27,193) = + reac_rate_local(193) 
  reac_source_local(29,193) = - reac_rate_local(193) 
  reac_source_local(34,193) = + reac_rate_local(193) 
  reac_source_local(01,194) = + reac_rate_local(194) 
  reac_source_local(03,194) = - reac_rate_local(194) 
  reac_source_local(16,194) = + reac_rate_local(194) 
  reac_source_local(23,194) = + reac_rate_local(194) 
  reac_source_local(01,195) = + reac_rate_local(195) 
  reac_source_local(16,195) = + reac_rate_local(195) 
  reac_source_local(23,195) = + reac_rate_local(195) 
  reac_source_local(29,195) = - reac_rate_local(195) 
  reac_source_local(01,196) = + reac_rate_local(196) 
  reac_source_local(03,196) = - reac_rate_local(196) 
  reac_source_local(04,196) = + reac_rate_local(196) 
  reac_source_local(23,196) = + reac_rate_local(196) 
  reac_source_local(34,196) = + reac_rate_local(196) 
  reac_source_local(01,197) = + reac_rate_local(197) 
  reac_source_local(04,197) = + reac_rate_local(197) 
  reac_source_local(23,197) = + reac_rate_local(197) 
  reac_source_local(29,197) = - reac_rate_local(197) 
  reac_source_local(34,197) = + reac_rate_local(197) 
  reac_source_local(01,198) = + reac_rate_local(198) 
  reac_source_local(03,198) = - reac_rate_local(198) 
  reac_source_local(09,198) = + reac_rate_local(198) 
  reac_source_local(23,198) = + reac_rate_local(198) * 2.d0
  reac_source_local(01,199) = + reac_rate_local(199) 
  reac_source_local(09,199) = + reac_rate_local(199) 
  reac_source_local(23,199) = + reac_rate_local(199) * 2.d0
  reac_source_local(29,199) = - reac_rate_local(199) 
  reac_source_local(01,200) = + reac_rate_local(200) 
  reac_source_local(03,200) = - reac_rate_local(200) 
  reac_source_local(11,200) = + reac_rate_local(200) 
  reac_source_local(49,200) = + reac_rate_local(200) 
  reac_source_local(01,201) = + reac_rate_local(201) 
  reac_source_local(11,201) = + reac_rate_local(201) 
  reac_source_local(29,201) = - reac_rate_local(201) 
  reac_source_local(49,201) = + reac_rate_local(201) 
  reac_source_local(01,202) = + reac_rate_local(202) 
  reac_source_local(03,202) = - reac_rate_local(202) 
  reac_source_local(06,202) = + reac_rate_local(202) 
  reac_source_local(46,202) = + reac_rate_local(202) 
  reac_source_local(01,203) = + reac_rate_local(203) 
  reac_source_local(06,203) = + reac_rate_local(203) 
  reac_source_local(29,203) = - reac_rate_local(203) 
  reac_source_local(46,203) = + reac_rate_local(203) 
  reac_source_local(01,204) = + reac_rate_local(204) 
  reac_source_local(05,204) = - reac_rate_local(204) 
  reac_source_local(16,204) = + reac_rate_local(204) 
  reac_source_local(01,205) = + reac_rate_local(205) 
  reac_source_local(15,205) = - reac_rate_local(205) 
  reac_source_local(16,205) = + reac_rate_local(205) 
  reac_source_local(01,206) = + reac_rate_local(206) 
  reac_source_local(04,206) = + reac_rate_local(206) 
  reac_source_local(05,206) = - reac_rate_local(206) 
  reac_source_local(34,206) = + reac_rate_local(206) 
  reac_source_local(01,207) = + reac_rate_local(207) 
  reac_source_local(04,207) = + reac_rate_local(207) 
  reac_source_local(15,207) = - reac_rate_local(207) 
  reac_source_local(34,207) = + reac_rate_local(207) 
  reac_source_local(01,208) = + reac_rate_local(208) 
  reac_source_local(05,208) = - reac_rate_local(208) 
  reac_source_local(43,208) = + reac_rate_local(208) 
  reac_source_local(49,208) = + reac_rate_local(208) 
  reac_source_local(01,209) = + reac_rate_local(209) 
  reac_source_local(15,209) = - reac_rate_local(209) 
  reac_source_local(43,209) = + reac_rate_local(209) 
  reac_source_local(49,209) = + reac_rate_local(209) 
  reac_source_local(01,210) = + reac_rate_local(210) 
  reac_source_local(05,210) = - reac_rate_local(210) 
  reac_source_local(46,210) = + reac_rate_local(210) 
  reac_source_local(47,210) = + reac_rate_local(210) 
  reac_source_local(01,211) = + reac_rate_local(211) 
  reac_source_local(15,211) = - reac_rate_local(211) 
  reac_source_local(46,211) = + reac_rate_local(211) 
  reac_source_local(47,211) = + reac_rate_local(211) 
  reac_source_local(01,212) = + reac_rate_local(212) 
  reac_source_local(05,212) = - reac_rate_local(212) 
  reac_source_local(11,212) = + reac_rate_local(212) 
  reac_source_local(13,212) = + reac_rate_local(212) 
  reac_source_local(01,213) = + reac_rate_local(213) 
  reac_source_local(11,213) = + reac_rate_local(213) 
  reac_source_local(13,213) = + reac_rate_local(213) 
  reac_source_local(15,213) = - reac_rate_local(213) 
  reac_source_local(01,214) = + reac_rate_local(214) 
  reac_source_local(09,214) = + reac_rate_local(214) 
  reac_source_local(50,214) = - reac_rate_local(214) 
  reac_source_local(01,215) = + reac_rate_local(215) 
  reac_source_local(09,215) = + reac_rate_local(215) 
  reac_source_local(37,215) = - reac_rate_local(215) 
  reac_source_local(01,216) = + reac_rate_local(216) 
  reac_source_local(09,216) = + reac_rate_local(216) 
  reac_source_local(40,216) = - reac_rate_local(216) 
  reac_source_local(01,217) = + reac_rate_local(217) 
  reac_source_local(13,217) = + reac_rate_local(217) 
  reac_source_local(43,217) = + reac_rate_local(217) 
  reac_source_local(50,217) = - reac_rate_local(217) 
  reac_source_local(01,218) = + reac_rate_local(218) 
  reac_source_local(13,218) = + reac_rate_local(218) 
  reac_source_local(37,218) = - reac_rate_local(218) 
  reac_source_local(43,218) = + reac_rate_local(218) 
  reac_source_local(01,219) = + reac_rate_local(219) 
  reac_source_local(13,219) = + reac_rate_local(219) 
  reac_source_local(40,219) = - reac_rate_local(219) 
  reac_source_local(43,219) = + reac_rate_local(219) 
  reac_source_local(01,220) = + reac_rate_local(220) 
  reac_source_local(24,220) = - reac_rate_local(220) 
  reac_source_local(41,220) = + reac_rate_local(220) 
  reac_source_local(01,221) = + reac_rate_local(221) 
  reac_source_local(02,221) = - reac_rate_local(221) 
  reac_source_local(41,221) = + reac_rate_local(221) 
  reac_source_local(01,222) = + reac_rate_local(222) 
  reac_source_local(18,222) = + reac_rate_local(222) 
  reac_source_local(24,222) = - reac_rate_local(222) 
  reac_source_local(34,222) = + reac_rate_local(222) 
  reac_source_local(01,223) = + reac_rate_local(223) 
  reac_source_local(02,223) = - reac_rate_local(223) 
  reac_source_local(18,223) = + reac_rate_local(223) 
  reac_source_local(34,223) = + reac_rate_local(223) 
  reac_source_local(01,224) = + reac_rate_local(224) 
  reac_source_local(23,224) = + reac_rate_local(224) 
  reac_source_local(24,224) = - reac_rate_local(224) 
  reac_source_local(38,224) = + reac_rate_local(224) 
  reac_source_local(01,225) = + reac_rate_local(225) 
  reac_source_local(02,225) = - reac_rate_local(225) 
  reac_source_local(23,225) = + reac_rate_local(225) 
  reac_source_local(38,225) = + reac_rate_local(225) 
  reac_source_local(01,226) = + reac_rate_local(226) 
  reac_source_local(23,226) = + reac_rate_local(226) 
  reac_source_local(24,226) = - reac_rate_local(226) 
  reac_source_local(34,226) = + reac_rate_local(226) 
  reac_source_local(39,226) = + reac_rate_local(226) 
  reac_source_local(01,227) = + reac_rate_local(227) 
  reac_source_local(02,227) = - reac_rate_local(227) 
  reac_source_local(23,227) = + reac_rate_local(227) 
  reac_source_local(34,227) = + reac_rate_local(227) 
  reac_source_local(39,227) = + reac_rate_local(227) 
  reac_source_local(01,228) = + reac_rate_local(228) 
  reac_source_local(20,228) = + reac_rate_local(228) 
  reac_source_local(23,228) = + reac_rate_local(228) * 2.d0
  reac_source_local(24,228) = - reac_rate_local(228) 
  reac_source_local(01,229) = + reac_rate_local(229) 
  reac_source_local(02,229) = - reac_rate_local(229) 
  reac_source_local(20,229) = + reac_rate_local(229) 
  reac_source_local(23,229) = + reac_rate_local(229) * 2.d0
  reac_source_local(01,230) = + reac_rate_local(230) 
  reac_source_local(11,230) = + reac_rate_local(230) 
  reac_source_local(24,230) = - reac_rate_local(230) 
  reac_source_local(27,230) = + reac_rate_local(230) 
  reac_source_local(01,231) = + reac_rate_local(231) 
  reac_source_local(02,231) = - reac_rate_local(231) 
  reac_source_local(11,231) = + reac_rate_local(231) 
  reac_source_local(27,231) = + reac_rate_local(231) 
  reac_source_local(01,232) = + reac_rate_local(232) 
  reac_source_local(06,232) = + reac_rate_local(232) 
  reac_source_local(16,232) = + reac_rate_local(232) 
  reac_source_local(24,232) = - reac_rate_local(232) 
  reac_source_local(01,233) = + reac_rate_local(233) 
  reac_source_local(02,233) = - reac_rate_local(233) 
  reac_source_local(06,233) = + reac_rate_local(233) 
  reac_source_local(16,233) = + reac_rate_local(233) 
  reac_source_local(01,234) = + reac_rate_local(234) 
  reac_source_local(08,234) = + reac_rate_local(234) 
  reac_source_local(24,234) = - reac_rate_local(234) 
  reac_source_local(49,234) = + reac_rate_local(234) 
  reac_source_local(01,235) = + reac_rate_local(235) 
  reac_source_local(02,235) = - reac_rate_local(235) 
  reac_source_local(08,235) = + reac_rate_local(235) 
  reac_source_local(49,235) = + reac_rate_local(235) 
  reac_source_local(01,236) = + reac_rate_local(236) 
  reac_source_local(14,236) = + reac_rate_local(236) 
  reac_source_local(24,236) = - reac_rate_local(236) 
  reac_source_local(46,236) = + reac_rate_local(236) 
  reac_source_local(01,237) = + reac_rate_local(237) 
  reac_source_local(02,237) = - reac_rate_local(237) 
  reac_source_local(14,237) = + reac_rate_local(237) 
  reac_source_local(46,237) = + reac_rate_local(237) 
  reac_source_local(23,238) = - reac_rate_local(238) 
  reac_source_local(34,238) = + reac_rate_local(238) * 2.d0
  reac_source_local(01,239) = + reac_rate_local(239) 
  reac_source_local(23,239) = - reac_rate_local(239) 
  reac_source_local(45,239) = + reac_rate_local(239) 
  reac_source_local(01,240) = - reac_rate_local(240) 
  reac_source_local(06,240) = + reac_rate_local(240) 
  reac_source_local(48,240) = - reac_rate_local(240) 
  reac_source_local(01,241) = - reac_rate_local(241) 
  reac_source_local(11,241) = + reac_rate_local(241) 
  reac_source_local(49,241) = - reac_rate_local(241) 
  reac_source_local(01,242) = - reac_rate_local(242) 
  reac_source_local(46,242) = - reac_rate_local(242) 
  reac_source_local(47,242) = + reac_rate_local(242) 
  reac_source_local(01,243) = - reac_rate_local(243) 
  reac_source_local(13,243) = - reac_rate_local(243) 
  reac_source_local(43,243) = + reac_rate_local(243) 
  reac_source_local(01,244) = - reac_rate_local(244) 
  reac_source_local(14,244) = + reac_rate_local(244) 
  reac_source_local(52,244) = - reac_rate_local(244) 
  reac_source_local(01,245) = - reac_rate_local(245) 
  reac_source_local(08,245) = + reac_rate_local(245) 
  reac_source_local(27,245) = - reac_rate_local(245) 
  reac_source_local(01,246) = - reac_rate_local(246) 
  reac_source_local(16,246) = - reac_rate_local(246) 
  reac_source_local(30,246) = + reac_rate_local(246) 
  reac_source_local(01,247) = - reac_rate_local(247) 
  reac_source_local(04,247) = - reac_rate_local(247) 
  reac_source_local(42,247) = + reac_rate_local(247) 
  reac_source_local(01,248) = - reac_rate_local(248) 
  reac_source_local(09,248) = - reac_rate_local(248) 
  reac_source_local(32,248) = + reac_rate_local(248) 
  reac_source_local(01,249) = - reac_rate_local(249) 
  reac_source_local(07,249) = + reac_rate_local(249) 
  reac_source_local(41,249) = - reac_rate_local(249) 
  reac_source_local(01,250) = - reac_rate_local(250) 
  reac_source_local(18,250) = - reac_rate_local(250) 
  reac_source_local(44,250) = + reac_rate_local(250) 
  reac_source_local(01,251) = - reac_rate_local(251) 
  reac_source_local(21,251) = + reac_rate_local(251) 
  reac_source_local(38,251) = - reac_rate_local(251) 
  reac_source_local(01,252) = - reac_rate_local(252) 
  reac_source_local(26,252) = + reac_rate_local(252) 
  reac_source_local(39,252) = - reac_rate_local(252) 
  reac_source_local(01,253) = - reac_rate_local(253) 
  reac_source_local(17,253) = + reac_rate_local(253) 
  reac_source_local(20,253) = - reac_rate_local(253) 
  reac_source_local(01,254) = - reac_rate_local(254) 
  reac_source_local(12,254) = - reac_rate_local(254) 
  reac_source_local(34,254) = + reac_rate_local(254) 
  reac_source_local(01,255) = - reac_rate_local(255) 
  reac_source_local(23,255) = + reac_rate_local(255) 
  reac_source_local(45,255) = - reac_rate_local(255) 
  reac_source_local(01,256) = - reac_rate_local(256) 
  reac_source_local(11,256) = + reac_rate_local(256) 
  reac_source_local(34,256) = + reac_rate_local(256) * 2.d0
  reac_source_local(35,256) = - reac_rate_local(256) 
  reac_source_local(01,257) = - reac_rate_local(257) 
  reac_source_local(23,257) = + reac_rate_local(257) 
  reac_source_local(34,257) = + reac_rate_local(257) 
  reac_source_local(35,257) = - reac_rate_local(257) 
  reac_source_local(47,257) = + reac_rate_local(257) 
  reac_source_local(01,258) = - reac_rate_local(258) 
  reac_source_local(11,258) = + reac_rate_local(258) 
  reac_source_local(34,258) = + reac_rate_local(258) 
  reac_source_local(48,258) = - reac_rate_local(258) 
  reac_source_local(01,259) = - reac_rate_local(259) 
  reac_source_local(34,259) = + reac_rate_local(259) * 2.d0
  reac_source_local(47,259) = + reac_rate_local(259) 
  reac_source_local(48,259) = - reac_rate_local(259) 
  reac_source_local(01,260) = - reac_rate_local(260) 
  reac_source_local(23,260) = + reac_rate_local(260) 
  reac_source_local(34,260) = + reac_rate_local(260) 
  reac_source_local(43,260) = + reac_rate_local(260) 
  reac_source_local(48,260) = - reac_rate_local(260) 
  reac_source_local(01,261) = - reac_rate_local(261) 
  reac_source_local(34,261) = + reac_rate_local(261) 
  reac_source_local(47,261) = + reac_rate_local(261) 
  reac_source_local(49,261) = - reac_rate_local(261) 
  reac_source_local(01,262) = - reac_rate_local(262) 
  reac_source_local(23,262) = + reac_rate_local(262) 
  reac_source_local(43,262) = + reac_rate_local(262) 
  reac_source_local(49,262) = - reac_rate_local(262) 
  reac_source_local(01,263) = - reac_rate_local(263) 
  reac_source_local(34,263) = + reac_rate_local(263) 
  reac_source_local(43,263) = + reac_rate_local(263) 
  reac_source_local(46,263) = - reac_rate_local(263) 
  reac_source_local(01,264) = - reac_rate_local(264) 
  reac_source_local(08,264) = + reac_rate_local(264) 
  reac_source_local(34,264) = + reac_rate_local(264) 
  reac_source_local(52,264) = - reac_rate_local(264) 
  reac_source_local(01,265) = - reac_rate_local(265) 
  reac_source_local(30,265) = + reac_rate_local(265) 
  reac_source_local(34,265) = + reac_rate_local(265) * 2.d0
  reac_source_local(52,265) = - reac_rate_local(265) 
  reac_source_local(01,266) = - reac_rate_local(266) 
  reac_source_local(27,266) = - reac_rate_local(266) 
  reac_source_local(30,266) = + reac_rate_local(266) 
  reac_source_local(34,266) = + reac_rate_local(266) 
  reac_source_local(01,267) = - reac_rate_local(267) 
  reac_source_local(27,267) = - reac_rate_local(267) 
  reac_source_local(34,267) = + reac_rate_local(267) * 2.d0
  reac_source_local(42,267) = + reac_rate_local(267) 
  reac_source_local(01,268) = - reac_rate_local(268) 
  reac_source_local(23,268) = + reac_rate_local(268) 
  reac_source_local(27,268) = - reac_rate_local(268) 
  reac_source_local(32,268) = + reac_rate_local(268) 
  reac_source_local(34,268) = + reac_rate_local(268) 
  reac_source_local(01,269) = - reac_rate_local(269) 
  reac_source_local(27,269) = - reac_rate_local(269) 
  reac_source_local(32,269) = + reac_rate_local(269) 
  reac_source_local(34,269) = + reac_rate_local(269) * 3.d0
  reac_source_local(01,270) = - reac_rate_local(270) 
  reac_source_local(11,270) = + reac_rate_local(270) 
  reac_source_local(27,270) = - reac_rate_local(270) 
  reac_source_local(47,270) = + reac_rate_local(270) 
  reac_source_local(01,271) = - reac_rate_local(271) 
  reac_source_local(16,271) = - reac_rate_local(271) 
  reac_source_local(34,271) = + reac_rate_local(271) 
  reac_source_local(42,271) = + reac_rate_local(271) 
  reac_source_local(01,272) = - reac_rate_local(272) 
  reac_source_local(16,272) = - reac_rate_local(272) 
  reac_source_local(32,272) = + reac_rate_local(272) 
  reac_source_local(34,272) = + reac_rate_local(272) * 2.d0
  reac_source_local(01,273) = - reac_rate_local(273) 
  reac_source_local(04,273) = - reac_rate_local(273) 
  reac_source_local(32,273) = + reac_rate_local(273) 
  reac_source_local(34,273) = + reac_rate_local(273) 
  reac_source_local(01,274) = - reac_rate_local(274) 
  reac_source_local(09,274) = - reac_rate_local(274) 
  reac_source_local(43,274) = + reac_rate_local(274) * 2.d0
  reac_source_local(01,275) = - reac_rate_local(275) 
  reac_source_local(34,275) = + reac_rate_local(275) * 2.d0
  reac_source_local(45,275) = - reac_rate_local(275) 
  reac_source_local(01,276) = - reac_rate_local(276) 
  reac_source_local(23,276) = + reac_rate_local(276) 
  reac_source_local(28,276) = - reac_rate_local(276) 
  reac_source_local(34,276) = + reac_rate_local(276) 
  reac_source_local(23,277) = - reac_rate_local(277) 
  reac_source_local(28,277) = + reac_rate_local(277) 
  reac_source_local(34,277) = + reac_rate_local(277) 
  reac_source_local(45,277) = - reac_rate_local(277) 
  reac_source_local(06,278) = + reac_rate_local(278) 
  reac_source_local(35,278) = - reac_rate_local(278) 
  reac_source_local(47,278) = - reac_rate_local(278) 
  reac_source_local(49,278) = + reac_rate_local(278) 
  reac_source_local(06,279) = + reac_rate_local(279) 
  reac_source_local(35,279) = - reac_rate_local(279) 
  reac_source_local(43,279) = - reac_rate_local(279) 
  reac_source_local(46,279) = + reac_rate_local(279) 
  reac_source_local(06,280) = + reac_rate_local(280) 
  reac_source_local(14,280) = - reac_rate_local(280) 
  reac_source_local(23,280) = + reac_rate_local(280) 
  reac_source_local(27,280) = + reac_rate_local(280) 
  reac_source_local(35,280) = - reac_rate_local(280) 
  reac_source_local(06,281) = + reac_rate_local(281) 
  reac_source_local(27,281) = + reac_rate_local(281) 
  reac_source_local(30,281) = - reac_rate_local(281) 
  reac_source_local(35,281) = - reac_rate_local(281) 
  reac_source_local(04,282) = + reac_rate_local(282) 
  reac_source_local(06,282) = + reac_rate_local(282) 
  reac_source_local(32,282) = - reac_rate_local(282) 
  reac_source_local(35,282) = - reac_rate_local(282) 
  reac_source_local(23,283) = + reac_rate_local(283) 
  reac_source_local(34,283) = - reac_rate_local(283) 
  reac_source_local(35,283) = - reac_rate_local(283) 
  reac_source_local(48,283) = + reac_rate_local(283) 
  reac_source_local(06,284) = - reac_rate_local(284) 
  reac_source_local(11,284) = + reac_rate_local(284) 
  reac_source_local(35,284) = + reac_rate_local(284) 
  reac_source_local(48,284) = - reac_rate_local(284) 
  reac_source_local(06,285) = + reac_rate_local(285) 
  reac_source_local(14,285) = - reac_rate_local(285) 
  reac_source_local(16,285) = + reac_rate_local(285) 
  reac_source_local(23,285) = + reac_rate_local(285) 
  reac_source_local(48,285) = - reac_rate_local(285) 
  reac_source_local(11,286) = + reac_rate_local(286) 
  reac_source_local(27,286) = + reac_rate_local(286) 
  reac_source_local(30,286) = - reac_rate_local(286) 
  reac_source_local(48,286) = - reac_rate_local(286) 
  reac_source_local(06,287) = + reac_rate_local(287) 
  reac_source_local(16,287) = + reac_rate_local(287) 
  reac_source_local(30,287) = - reac_rate_local(287) 
  reac_source_local(48,287) = - reac_rate_local(287) 
  reac_source_local(04,288) = + reac_rate_local(288) 
  reac_source_local(11,288) = + reac_rate_local(288) 
  reac_source_local(32,288) = - reac_rate_local(288) 
  reac_source_local(48,288) = - reac_rate_local(288) 
  reac_source_local(06,289) = + reac_rate_local(289) 
  reac_source_local(09,289) = + reac_rate_local(289) 
  reac_source_local(32,289) = - reac_rate_local(289) 
  reac_source_local(48,289) = - reac_rate_local(289) 
  reac_source_local(23,290) = - reac_rate_local(290) 
  reac_source_local(34,290) = + reac_rate_local(290) 
  reac_source_local(35,290) = + reac_rate_local(290) 
  reac_source_local(48,290) = - reac_rate_local(290) 
  reac_source_local(23,291) = + reac_rate_local(291) 
  reac_source_local(34,291) = - reac_rate_local(291) 
  reac_source_local(48,291) = - reac_rate_local(291) 
  reac_source_local(49,291) = + reac_rate_local(291) 
  reac_source_local(06,292) = - reac_rate_local(292) 
  reac_source_local(11,292) = + reac_rate_local(292) 
  reac_source_local(48,292) = + reac_rate_local(292) 
  reac_source_local(49,292) = - reac_rate_local(292) 
  reac_source_local(06,293) = - reac_rate_local(293) 
  reac_source_local(23,293) = + reac_rate_local(293) 
  reac_source_local(27,293) = + reac_rate_local(293) 
  reac_source_local(49,293) = - reac_rate_local(293) 
  reac_source_local(04,294) = + reac_rate_local(294) 
  reac_source_local(23,294) = + reac_rate_local(294) 
  reac_source_local(47,294) = - reac_rate_local(294) 
  reac_source_local(49,294) = - reac_rate_local(294) 
  reac_source_local(09,295) = + reac_rate_local(295) 
  reac_source_local(23,295) = + reac_rate_local(295) 
  reac_source_local(43,295) = - reac_rate_local(295) 
  reac_source_local(49,295) = - reac_rate_local(295) 
  reac_source_local(06,296) = + reac_rate_local(296) 
  reac_source_local(14,296) = - reac_rate_local(296) 
  reac_source_local(27,296) = + reac_rate_local(296) 
  reac_source_local(49,296) = - reac_rate_local(296) 
  reac_source_local(04,297) = + reac_rate_local(297) 
  reac_source_local(06,297) = + reac_rate_local(297) 
  reac_source_local(30,297) = - reac_rate_local(297) 
  reac_source_local(49,297) = - reac_rate_local(297) 
  reac_source_local(04,298) = + reac_rate_local(298) 
  reac_source_local(11,298) = + reac_rate_local(298) 
  reac_source_local(42,298) = - reac_rate_local(298) 
  reac_source_local(49,298) = - reac_rate_local(298) 
  reac_source_local(06,299) = - reac_rate_local(299) 
  reac_source_local(11,299) = + reac_rate_local(299) 
  reac_source_local(46,299) = - reac_rate_local(299) 
  reac_source_local(49,299) = + reac_rate_local(299) 
  reac_source_local(06,300) = - reac_rate_local(300) 
  reac_source_local(27,300) = + reac_rate_local(300) 
  reac_source_local(34,300) = + reac_rate_local(300) 
  reac_source_local(46,300) = - reac_rate_local(300) 
  reac_source_local(06,301) = - reac_rate_local(301) 
  reac_source_local(16,301) = + reac_rate_local(301) 
  reac_source_local(23,301) = + reac_rate_local(301) 
  reac_source_local(46,301) = - reac_rate_local(301) 
  reac_source_local(04,302) = + reac_rate_local(302) 
  reac_source_local(06,302) = - reac_rate_local(302) 
  reac_source_local(23,302) = + reac_rate_local(302) 
  reac_source_local(34,302) = + reac_rate_local(302) 
  reac_source_local(46,302) = - reac_rate_local(302) 
  reac_source_local(06,303) = - reac_rate_local(303) 
  reac_source_local(09,303) = + reac_rate_local(303) 
  reac_source_local(23,303) = + reac_rate_local(303) * 2.d0
  reac_source_local(46,303) = - reac_rate_local(303) 
  reac_source_local(23,304) = - reac_rate_local(304) 
  reac_source_local(34,304) = + reac_rate_local(304) 
  reac_source_local(46,304) = - reac_rate_local(304) 
  reac_source_local(49,304) = + reac_rate_local(304) 
  reac_source_local(06,305) = - reac_rate_local(305) 
  reac_source_local(13,305) = - reac_rate_local(305) 
  reac_source_local(16,305) = + reac_rate_local(305) 
  reac_source_local(34,305) = + reac_rate_local(305) 
  reac_source_local(04,306) = + reac_rate_local(306) 
  reac_source_local(06,306) = - reac_rate_local(306) 
  reac_source_local(13,306) = - reac_rate_local(306) 
  reac_source_local(23,306) = + reac_rate_local(306) 
  reac_source_local(06,307) = - reac_rate_local(307) 
  reac_source_local(09,307) = + reac_rate_local(307) 
  reac_source_local(13,307) = - reac_rate_local(307) 
  reac_source_local(23,307) = + reac_rate_local(307) 
  reac_source_local(34,307) = + reac_rate_local(307) 
  reac_source_local(13,308) = - reac_rate_local(308) 
  reac_source_local(23,308) = - reac_rate_local(308) 
  reac_source_local(34,308) = + reac_rate_local(308) 
  reac_source_local(46,308) = + reac_rate_local(308) 
  reac_source_local(14,309) = + reac_rate_local(309) 
  reac_source_local(16,309) = + reac_rate_local(309) 
  reac_source_local(30,309) = - reac_rate_local(309) 
  reac_source_local(52,309) = - reac_rate_local(309) 
  reac_source_local(27,310) = + reac_rate_local(310) 
  reac_source_local(32,310) = - reac_rate_local(310) 
  reac_source_local(42,310) = + reac_rate_local(310) 
  reac_source_local(52,310) = - reac_rate_local(310) 
  reac_source_local(23,311) = + reac_rate_local(311) 
  reac_source_local(27,311) = + reac_rate_local(311) 
  reac_source_local(34,311) = - reac_rate_local(311) 
  reac_source_local(52,311) = - reac_rate_local(311) 
  reac_source_local(16,312) = + reac_rate_local(312) 
  reac_source_local(23,312) = + reac_rate_local(312) 
  reac_source_local(27,312) = - reac_rate_local(312) 
  reac_source_local(34,312) = - reac_rate_local(312) 
  reac_source_local(16,313) = - reac_rate_local(313) 
  reac_source_local(27,313) = + reac_rate_local(313) 
  reac_source_local(32,313) = + reac_rate_local(313) 
  reac_source_local(42,313) = - reac_rate_local(313) 
  reac_source_local(04,314) = + reac_rate_local(314) 
  reac_source_local(16,314) = - reac_rate_local(314) 
  reac_source_local(30,314) = + reac_rate_local(314) 
  reac_source_local(42,314) = - reac_rate_local(314) 
  reac_source_local(04,315) = + reac_rate_local(315) 
  reac_source_local(16,315) = - reac_rate_local(315) 
  reac_source_local(23,315) = + reac_rate_local(315) 
  reac_source_local(34,315) = - reac_rate_local(315) 
  reac_source_local(04,316) = - reac_rate_local(316) 
  reac_source_local(14,316) = - reac_rate_local(316) 
  reac_source_local(27,316) = + reac_rate_local(316) 
  reac_source_local(30,316) = + reac_rate_local(316) 
  reac_source_local(04,317) = - reac_rate_local(317) 
  reac_source_local(27,317) = + reac_rate_local(317) 
  reac_source_local(30,317) = - reac_rate_local(317) 
  reac_source_local(32,317) = + reac_rate_local(317) 
  reac_source_local(04,318) = - reac_rate_local(318) 
  reac_source_local(09,318) = + reac_rate_local(318) 
  reac_source_local(32,318) = + reac_rate_local(318) 
  reac_source_local(33,318) = - reac_rate_local(318) 
  reac_source_local(04,319) = - reac_rate_local(319) 
  reac_source_local(09,319) = + reac_rate_local(319) 
  reac_source_local(23,319) = + reac_rate_local(319) 
  reac_source_local(34,319) = - reac_rate_local(319) 
  reac_source_local(04,320) = + reac_rate_local(320) 
  reac_source_local(06,320) = - reac_rate_local(320) 
  reac_source_local(09,320) = - reac_rate_local(320) 
  reac_source_local(11,320) = + reac_rate_local(320) 
  reac_source_local(09,321) = - reac_rate_local(321) 
  reac_source_local(14,321) = - reac_rate_local(321) 
  reac_source_local(27,321) = + reac_rate_local(321) 
  reac_source_local(42,321) = + reac_rate_local(321) 
  reac_source_local(09,322) = - reac_rate_local(322) 
  reac_source_local(14,322) = - reac_rate_local(322) 
  reac_source_local(16,322) = + reac_rate_local(322) 
  reac_source_local(30,322) = + reac_rate_local(322) 
  reac_source_local(09,323) = - reac_rate_local(323) 
  reac_source_local(16,323) = + reac_rate_local(323) 
  reac_source_local(30,323) = - reac_rate_local(323) 
  reac_source_local(32,323) = + reac_rate_local(323) 
  reac_source_local(04,324) = + reac_rate_local(324) 
  reac_source_local(09,324) = - reac_rate_local(324) 
  reac_source_local(32,324) = + reac_rate_local(324) 
  reac_source_local(42,324) = - reac_rate_local(324) 
  reac_source_local(04,325) = + reac_rate_local(325) 
  reac_source_local(09,325) = - reac_rate_local(325) 
  reac_source_local(23,325) = - reac_rate_local(325) 
  reac_source_local(34,325) = + reac_rate_local(325) 
  reac_source_local(06,326) = - reac_rate_local(326) 
  reac_source_local(09,326) = + reac_rate_local(326) 
  reac_source_local(11,326) = + reac_rate_local(326) 
  reac_source_local(36,326) = - reac_rate_local(326) 
  reac_source_local(06,327) = - reac_rate_local(327) 
  reac_source_local(23,327) = + reac_rate_local(327) 
  reac_source_local(28,327) = - reac_rate_local(327) 
  reac_source_local(35,327) = + reac_rate_local(327) 
  reac_source_local(11,328) = - reac_rate_local(328) 
  reac_source_local(23,328) = + reac_rate_local(328) 
  reac_source_local(28,328) = - reac_rate_local(328) 
  reac_source_local(48,328) = + reac_rate_local(328) 
  reac_source_local(23,329) = + reac_rate_local(329) 
  reac_source_local(28,329) = - reac_rate_local(329) 
  reac_source_local(47,329) = - reac_rate_local(329) 
  reac_source_local(49,329) = + reac_rate_local(329) 
  reac_source_local(23,330) = + reac_rate_local(330) 
  reac_source_local(28,330) = - reac_rate_local(330) 
  reac_source_local(43,330) = - reac_rate_local(330) 
  reac_source_local(46,330) = + reac_rate_local(330) 
  reac_source_local(14,331) = - reac_rate_local(331) 
  reac_source_local(23,331) = + reac_rate_local(331) * 2.d0
  reac_source_local(27,331) = + reac_rate_local(331) 
  reac_source_local(28,331) = - reac_rate_local(331) 
  reac_source_local(08,332) = - reac_rate_local(332) 
  reac_source_local(23,332) = + reac_rate_local(332) 
  reac_source_local(28,332) = - reac_rate_local(332) 
  reac_source_local(52,332) = + reac_rate_local(332) 
  reac_source_local(23,333) = + reac_rate_local(333) 
  reac_source_local(27,333) = + reac_rate_local(333) 
  reac_source_local(28,333) = - reac_rate_local(333) 
  reac_source_local(30,333) = - reac_rate_local(333) 
  reac_source_local(04,334) = + reac_rate_local(334) 
  reac_source_local(23,334) = + reac_rate_local(334) * 2.d0
  reac_source_local(28,334) = - reac_rate_local(334) 
  reac_source_local(30,334) = - reac_rate_local(334) 
  reac_source_local(16,335) = + reac_rate_local(335) 
  reac_source_local(23,335) = + reac_rate_local(335) 
  reac_source_local(28,335) = - reac_rate_local(335) 
  reac_source_local(42,335) = - reac_rate_local(335) 
  reac_source_local(09,336) = + reac_rate_local(336) 
  reac_source_local(23,336) = + reac_rate_local(336) 
  reac_source_local(28,336) = - reac_rate_local(336) 
  reac_source_local(33,336) = - reac_rate_local(336) 
  reac_source_local(04,337) = + reac_rate_local(337) 
  reac_source_local(23,337) = + reac_rate_local(337) 
  reac_source_local(28,337) = - reac_rate_local(337) 
  reac_source_local(32,337) = - reac_rate_local(337) 
  reac_source_local(06,338) = - reac_rate_local(338) 
  reac_source_local(34,338) = + reac_rate_local(338) 
  reac_source_local(35,338) = + reac_rate_local(338) 
  reac_source_local(45,338) = - reac_rate_local(338) 
  reac_source_local(06,339) = - reac_rate_local(339) 
  reac_source_local(23,339) = + reac_rate_local(339) 
  reac_source_local(45,339) = - reac_rate_local(339) 
  reac_source_local(48,339) = + reac_rate_local(339) 
  reac_source_local(06,340) = - reac_rate_local(340) 
  reac_source_local(23,340) = + reac_rate_local(340) 
  reac_source_local(34,340) = + reac_rate_local(340) 
  reac_source_local(45,340) = - reac_rate_local(340) 
  reac_source_local(49,340) = + reac_rate_local(340) 
  reac_source_local(34,341) = + reac_rate_local(341) 
  reac_source_local(45,341) = - reac_rate_local(341) 
  reac_source_local(47,341) = - reac_rate_local(341) 
  reac_source_local(49,341) = + reac_rate_local(341) 
  reac_source_local(23,342) = + reac_rate_local(342) 
  reac_source_local(45,342) = - reac_rate_local(342) 
  reac_source_local(46,342) = + reac_rate_local(342) 
  reac_source_local(47,342) = - reac_rate_local(342) 
  reac_source_local(34,343) = + reac_rate_local(343) 
  reac_source_local(43,343) = - reac_rate_local(343) 
  reac_source_local(45,343) = - reac_rate_local(343) 
  reac_source_local(46,343) = + reac_rate_local(343) 
  reac_source_local(13,344) = + reac_rate_local(344) 
  reac_source_local(23,344) = + reac_rate_local(344) 
  reac_source_local(43,344) = - reac_rate_local(344) 
  reac_source_local(45,344) = - reac_rate_local(344) 
  reac_source_local(14,345) = - reac_rate_local(345) 
  reac_source_local(23,345) = + reac_rate_local(345) 
  reac_source_local(45,345) = - reac_rate_local(345) 
  reac_source_local(52,345) = + reac_rate_local(345) 
  reac_source_local(14,346) = - reac_rate_local(346) 
  reac_source_local(23,346) = + reac_rate_local(346) 
  reac_source_local(27,346) = + reac_rate_local(346) 
  reac_source_local(34,346) = + reac_rate_local(346) 
  reac_source_local(45,346) = - reac_rate_local(346) 
  reac_source_local(14,347) = - reac_rate_local(347) 
  reac_source_local(16,347) = + reac_rate_local(347) 
  reac_source_local(23,347) = + reac_rate_local(347) * 2.d0
  reac_source_local(45,347) = - reac_rate_local(347) 
  reac_source_local(04,348) = + reac_rate_local(348) 
  reac_source_local(14,348) = - reac_rate_local(348) 
  reac_source_local(23,348) = + reac_rate_local(348) * 2.d0
  reac_source_local(34,348) = + reac_rate_local(348) 
  reac_source_local(45,348) = - reac_rate_local(348) 
  reac_source_local(09,349) = + reac_rate_local(349) 
  reac_source_local(14,349) = - reac_rate_local(349) 
  reac_source_local(23,349) = + reac_rate_local(349) * 3.d0
  reac_source_local(45,349) = - reac_rate_local(349) 
  reac_source_local(16,350) = + reac_rate_local(350) 
  reac_source_local(23,350) = + reac_rate_local(350) 
  reac_source_local(30,350) = - reac_rate_local(350) 
  reac_source_local(45,350) = - reac_rate_local(350) 
  reac_source_local(04,351) = + reac_rate_local(351) 
  reac_source_local(23,351) = + reac_rate_local(351) 
  reac_source_local(30,351) = - reac_rate_local(351) 
  reac_source_local(34,351) = + reac_rate_local(351) 
  reac_source_local(45,351) = - reac_rate_local(351) 
  reac_source_local(09,352) = + reac_rate_local(352) 
  reac_source_local(23,352) = + reac_rate_local(352) * 2.d0
  reac_source_local(30,352) = - reac_rate_local(352) 
  reac_source_local(45,352) = - reac_rate_local(352) 
  reac_source_local(04,353) = + reac_rate_local(353) 
  reac_source_local(32,353) = - reac_rate_local(353) 
  reac_source_local(34,353) = + reac_rate_local(353) 
  reac_source_local(45,353) = - reac_rate_local(353) 
  reac_source_local(09,354) = + reac_rate_local(354) 
  reac_source_local(23,354) = + reac_rate_local(354) 
  reac_source_local(32,354) = - reac_rate_local(354) 
  reac_source_local(45,354) = - reac_rate_local(354) 
  reac_source_local(28,355) = + reac_rate_local(355) 
  reac_source_local(34,355) = - reac_rate_local(355) 
  reac_source_local(45,355) = - reac_rate_local(355) 
  reac_source_local(12,356) = + reac_rate_local(356) 
  reac_source_local(23,356) = + reac_rate_local(356) 
  reac_source_local(34,356) = - reac_rate_local(356) 
  reac_source_local(45,356) = - reac_rate_local(356) 
  reac_source_local(06,357) = - reac_rate_local(357) 
  reac_source_local(12,357) = - reac_rate_local(357) 
  reac_source_local(34,357) = + reac_rate_local(357) 
  reac_source_local(48,357) = + reac_rate_local(357) 
  reac_source_local(06,358) = - reac_rate_local(358) 
  reac_source_local(12,358) = - reac_rate_local(358) 
  reac_source_local(23,358) = + reac_rate_local(358) 
  reac_source_local(49,358) = + reac_rate_local(358) 
  reac_source_local(11,359) = - reac_rate_local(359) 
  reac_source_local(12,359) = - reac_rate_local(359) 
  reac_source_local(34,359) = + reac_rate_local(359) 
  reac_source_local(49,359) = + reac_rate_local(359) 
  reac_source_local(12,360) = - reac_rate_local(360) 
  reac_source_local(34,360) = + reac_rate_local(360) 
  reac_source_local(46,360) = + reac_rate_local(360) 
  reac_source_local(47,360) = - reac_rate_local(360) 
  reac_source_local(12,361) = - reac_rate_local(361) 
  reac_source_local(13,361) = + reac_rate_local(361) 
  reac_source_local(23,361) = + reac_rate_local(361) 
  reac_source_local(47,361) = - reac_rate_local(361) 
  reac_source_local(12,362) = - reac_rate_local(362) 
  reac_source_local(13,362) = + reac_rate_local(362) 
  reac_source_local(34,362) = + reac_rate_local(362) 
  reac_source_local(43,362) = - reac_rate_local(362) 
  reac_source_local(12,363) = - reac_rate_local(363) 
  reac_source_local(14,363) = - reac_rate_local(363) 
  reac_source_local(23,363) = + reac_rate_local(363) 
  reac_source_local(27,363) = + reac_rate_local(363) 
  reac_source_local(12,364) = - reac_rate_local(364) 
  reac_source_local(14,364) = - reac_rate_local(364) 
  reac_source_local(16,364) = + reac_rate_local(364) 
  reac_source_local(23,364) = + reac_rate_local(364) 
  reac_source_local(34,364) = + reac_rate_local(364) 
  reac_source_local(04,365) = + reac_rate_local(365) 
  reac_source_local(12,365) = - reac_rate_local(365) 
  reac_source_local(14,365) = - reac_rate_local(365) 
  reac_source_local(23,365) = + reac_rate_local(365) * 2.d0
  reac_source_local(08,366) = - reac_rate_local(366) 
  reac_source_local(12,366) = - reac_rate_local(366) 
  reac_source_local(16,366) = + reac_rate_local(366) 
  reac_source_local(23,366) = + reac_rate_local(366) 
  reac_source_local(04,367) = + reac_rate_local(367) 
  reac_source_local(08,367) = - reac_rate_local(367) 
  reac_source_local(12,367) = - reac_rate_local(367) 
  reac_source_local(23,367) = + reac_rate_local(367) 
  reac_source_local(34,367) = + reac_rate_local(367) 
  reac_source_local(12,368) = - reac_rate_local(368) 
  reac_source_local(16,368) = + reac_rate_local(368) 
  reac_source_local(30,368) = - reac_rate_local(368) 
  reac_source_local(34,368) = + reac_rate_local(368) 
  reac_source_local(04,369) = + reac_rate_local(369) 
  reac_source_local(12,369) = - reac_rate_local(369) 
  reac_source_local(23,369) = + reac_rate_local(369) 
  reac_source_local(30,369) = - reac_rate_local(369) 
  reac_source_local(09,370) = + reac_rate_local(370) 
  reac_source_local(12,370) = - reac_rate_local(370) 
  reac_source_local(23,370) = + reac_rate_local(370) 
  reac_source_local(30,370) = - reac_rate_local(370) 
  reac_source_local(34,370) = + reac_rate_local(370) 
  reac_source_local(04,371) = + reac_rate_local(371) 
  reac_source_local(12,371) = - reac_rate_local(371) 
  reac_source_local(34,371) = + reac_rate_local(371) 
  reac_source_local(42,371) = - reac_rate_local(371) 
  reac_source_local(09,372) = + reac_rate_local(372) 
  reac_source_local(12,372) = - reac_rate_local(372) 
  reac_source_local(23,372) = + reac_rate_local(372) 
  reac_source_local(42,372) = - reac_rate_local(372) 
  reac_source_local(09,373) = + reac_rate_local(373) 
  reac_source_local(12,373) = - reac_rate_local(373) 
  reac_source_local(32,373) = - reac_rate_local(373) 
  reac_source_local(34,373) = + reac_rate_local(373) 
  reac_source_local(06,374) = - reac_rate_local(374) 
  reac_source_local(11,374) = + reac_rate_local(374) * 2.d0
  reac_source_local(47,374) = - reac_rate_local(374) 
  reac_source_local(06,375) = - reac_rate_local(375) 
  reac_source_local(30,375) = + reac_rate_local(375) 
  reac_source_local(34,375) = + reac_rate_local(375) 
  reac_source_local(43,375) = - reac_rate_local(375) 
  reac_source_local(06,376) = - reac_rate_local(376) 
  reac_source_local(08,376) = - reac_rate_local(376) 
  reac_source_local(11,376) = + reac_rate_local(376) 
  reac_source_local(14,376) = + reac_rate_local(376) 
  reac_source_local(06,377) = - reac_rate_local(377) 
  reac_source_local(11,377) = + reac_rate_local(377) 
  reac_source_local(30,377) = + reac_rate_local(377) 
  reac_source_local(42,377) = - reac_rate_local(377) 
  reac_source_local(06,378) = - reac_rate_local(378) 
  reac_source_local(11,378) = + reac_rate_local(378) 
  reac_source_local(32,378) = + reac_rate_local(378) 
  reac_source_local(33,378) = - reac_rate_local(378) 
  reac_source_local(06,379) = - reac_rate_local(379) 
  reac_source_local(07,379) = + reac_rate_local(379) 
  reac_source_local(11,379) = + reac_rate_local(379) 
  reac_source_local(44,379) = - reac_rate_local(379) 
  reac_source_local(06,380) = - reac_rate_local(380) 
  reac_source_local(11,380) = + reac_rate_local(380) 
  reac_source_local(21,380) = + reac_rate_local(380) 
  reac_source_local(26,380) = - reac_rate_local(380) 
  reac_source_local(06,381) = - reac_rate_local(381) 
  reac_source_local(11,381) = + reac_rate_local(381) 
  reac_source_local(23,381) = + reac_rate_local(381) 
  reac_source_local(34,381) = - reac_rate_local(381) 
  reac_source_local(06,382) = - reac_rate_local(382) 
  reac_source_local(11,382) = - reac_rate_local(382) 
  reac_source_local(14,382) = + reac_rate_local(382) 
  reac_source_local(34,382) = + reac_rate_local(382) 
  reac_source_local(06,383) = - reac_rate_local(383) 
  reac_source_local(11,383) = + reac_rate_local(383) 
  reac_source_local(25,383) = + reac_rate_local(383) 
  reac_source_local(31,383) = - reac_rate_local(383) 
  reac_source_local(06,384) = - reac_rate_local(384) 
  reac_source_local(14,384) = + reac_rate_local(384) 
  reac_source_local(47,384) = - reac_rate_local(384) 
  reac_source_local(06,385) = - reac_rate_local(385) 
  reac_source_local(11,385) = + reac_rate_local(385) 
  reac_source_local(34,385) = + reac_rate_local(385) 
  reac_source_local(11,386) = - reac_rate_local(386) 
  reac_source_local(34,386) = + reac_rate_local(386) 
  reac_source_local(47,386) = + reac_rate_local(386) 
  reac_source_local(11,387) = - reac_rate_local(387) 
  reac_source_local(23,387) = + reac_rate_local(387) 
  reac_source_local(43,387) = + reac_rate_local(387) 
  reac_source_local(08,388) = - reac_rate_local(388) 
  reac_source_local(11,388) = - reac_rate_local(388) 
  reac_source_local(14,388) = + reac_rate_local(388) 
  reac_source_local(47,388) = + reac_rate_local(388) 
  reac_source_local(32,389) = + reac_rate_local(389) 
  reac_source_local(34,389) = + reac_rate_local(389) * 2.d0
  reac_source_local(47,389) = - reac_rate_local(389) * 2.d0
  reac_source_local(08,390) = - reac_rate_local(390) 
  reac_source_local(11,390) = + reac_rate_local(390) 
  reac_source_local(30,390) = + reac_rate_local(390) 
  reac_source_local(47,390) = - reac_rate_local(390) 
  reac_source_local(11,391) = + reac_rate_local(391) 
  reac_source_local(32,391) = + reac_rate_local(391) 
  reac_source_local(42,391) = - reac_rate_local(391) 
  reac_source_local(47,391) = - reac_rate_local(391) 
  reac_source_local(32,392) = + reac_rate_local(392) 
  reac_source_local(33,392) = - reac_rate_local(392) 
  reac_source_local(43,392) = + reac_rate_local(392) 
  reac_source_local(47,392) = - reac_rate_local(392) 
  reac_source_local(07,393) = - reac_rate_local(393) 
  reac_source_local(11,393) = + reac_rate_local(393) 
  reac_source_local(44,393) = + reac_rate_local(393) 
  reac_source_local(47,393) = - reac_rate_local(393) 
  reac_source_local(08,394) = + reac_rate_local(394) 
  reac_source_local(30,394) = + reac_rate_local(394) 
  reac_source_local(44,394) = - reac_rate_local(394) 
  reac_source_local(47,394) = - reac_rate_local(394) 
  reac_source_local(11,395) = + reac_rate_local(395) 
  reac_source_local(21,395) = + reac_rate_local(395) 
  reac_source_local(44,395) = - reac_rate_local(395) 
  reac_source_local(47,395) = - reac_rate_local(395) 
  reac_source_local(11,396) = + reac_rate_local(396) 
  reac_source_local(21,396) = - reac_rate_local(396) 
  reac_source_local(26,396) = + reac_rate_local(396) 
  reac_source_local(47,396) = - reac_rate_local(396) 
  reac_source_local(11,397) = + reac_rate_local(397) 
  reac_source_local(23,397) = - reac_rate_local(397) 
  reac_source_local(34,397) = + reac_rate_local(397) 
  reac_source_local(47,397) = - reac_rate_local(397) 
  reac_source_local(23,398) = + reac_rate_local(398) 
  reac_source_local(34,398) = - reac_rate_local(398) 
  reac_source_local(43,398) = + reac_rate_local(398) 
  reac_source_local(47,398) = - reac_rate_local(398) 
  reac_source_local(34,399) = + reac_rate_local(399) 
  reac_source_local(43,399) = + reac_rate_local(399) 
  reac_source_local(47,399) = - reac_rate_local(399) 
  reac_source_local(23,400) = + reac_rate_local(400) 
  reac_source_local(32,400) = + reac_rate_local(400) 
  reac_source_local(47,400) = - reac_rate_local(400) * 2.d0
  reac_source_local(11,401) = + reac_rate_local(401) 
  reac_source_local(34,401) = - reac_rate_local(401) 
  reac_source_local(47,401) = - reac_rate_local(401) 
  reac_source_local(14,402) = - reac_rate_local(402) 
  reac_source_local(21,402) = + reac_rate_local(402) 
  reac_source_local(34,402) = + reac_rate_local(402) 
  reac_source_local(43,402) = - reac_rate_local(402) 
  reac_source_local(14,403) = - reac_rate_local(403) 
  reac_source_local(43,403) = - reac_rate_local(403) 
  reac_source_local(44,403) = + reac_rate_local(403) 
  reac_source_local(23,404) = - reac_rate_local(404) 
  reac_source_local(34,404) = + reac_rate_local(404) 
  reac_source_local(43,404) = - reac_rate_local(404) 
  reac_source_local(47,404) = + reac_rate_local(404) 
  reac_source_local(11,405) = - reac_rate_local(405) 
  reac_source_local(34,405) = + reac_rate_local(405) 
  reac_source_local(42,405) = + reac_rate_local(405) 
  reac_source_local(43,405) = - reac_rate_local(405) 
  reac_source_local(32,406) = + reac_rate_local(406) 
  reac_source_local(34,406) = + reac_rate_local(406) 
  reac_source_local(43,406) = - reac_rate_local(406) 
  reac_source_local(47,406) = - reac_rate_local(406) 
  reac_source_local(11,407) = + reac_rate_local(407) 
  reac_source_local(23,407) = - reac_rate_local(407) 
  reac_source_local(43,407) = - reac_rate_local(407) 
  reac_source_local(32,408) = + reac_rate_local(408) 
  reac_source_local(42,408) = - reac_rate_local(408) 
  reac_source_local(43,408) = - reac_rate_local(408) 
  reac_source_local(47,408) = + reac_rate_local(408) 
  reac_source_local(08,409) = + reac_rate_local(409) 
  reac_source_local(14,409) = - reac_rate_local(409) 
  reac_source_local(30,409) = + reac_rate_local(409) 
  reac_source_local(42,409) = - reac_rate_local(409) 
  reac_source_local(07,410) = + reac_rate_local(410) 
  reac_source_local(08,410) = + reac_rate_local(410) 
  reac_source_local(14,410) = - reac_rate_local(410) 
  reac_source_local(44,410) = - reac_rate_local(410) 
  reac_source_local(08,411) = + reac_rate_local(411) 
  reac_source_local(14,411) = - reac_rate_local(411) 
  reac_source_local(21,411) = + reac_rate_local(411) 
  reac_source_local(26,411) = - reac_rate_local(411) 
  reac_source_local(08,412) = + reac_rate_local(412) 
  reac_source_local(14,412) = - reac_rate_local(412) 
  reac_source_local(23,412) = + reac_rate_local(412) 
  reac_source_local(34,412) = - reac_rate_local(412) 
  reac_source_local(06,413) = + reac_rate_local(413) 
  reac_source_local(11,413) = + reac_rate_local(413) 
  reac_source_local(14,413) = - reac_rate_local(413) 
  reac_source_local(34,413) = - reac_rate_local(413) 
  reac_source_local(11,414) = + reac_rate_local(414) * 2.d0
  reac_source_local(14,414) = - reac_rate_local(414) 
  reac_source_local(08,415) = + reac_rate_local(415) 
  reac_source_local(14,415) = - reac_rate_local(415) 
  reac_source_local(25,415) = + reac_rate_local(415) 
  reac_source_local(31,415) = - reac_rate_local(415) 
  reac_source_local(11,416) = + reac_rate_local(416) 
  reac_source_local(14,416) = - reac_rate_local(416) 
  reac_source_local(30,416) = + reac_rate_local(416) 
  reac_source_local(43,416) = - reac_rate_local(416) 
  reac_source_local(08,417) = + reac_rate_local(417) 
  reac_source_local(11,417) = + reac_rate_local(417) 
  reac_source_local(14,417) = - reac_rate_local(417) 
  reac_source_local(47,417) = - reac_rate_local(417) 
  reac_source_local(08,418) = - reac_rate_local(418) * 2.d0
  reac_source_local(14,418) = + reac_rate_local(418) 
  reac_source_local(30,418) = + reac_rate_local(418) 
  reac_source_local(08,419) = - reac_rate_local(419) 
  reac_source_local(14,419) = + reac_rate_local(419) 
  reac_source_local(30,419) = - reac_rate_local(419) 
  reac_source_local(42,419) = + reac_rate_local(419) 
  reac_source_local(08,420) = - reac_rate_local(420) 
  reac_source_local(14,420) = + reac_rate_local(420) 
  reac_source_local(32,420) = - reac_rate_local(420) 
  reac_source_local(33,420) = + reac_rate_local(420) 
  reac_source_local(08,421) = - reac_rate_local(421) 
  reac_source_local(30,421) = + reac_rate_local(421) 
  reac_source_local(32,421) = + reac_rate_local(421) 
  reac_source_local(33,421) = - reac_rate_local(421) 
  reac_source_local(07,422) = - reac_rate_local(422) 
  reac_source_local(08,422) = - reac_rate_local(422) 
  reac_source_local(14,422) = + reac_rate_local(422) 
  reac_source_local(44,422) = + reac_rate_local(422) 
  reac_source_local(07,423) = + reac_rate_local(423) 
  reac_source_local(08,423) = - reac_rate_local(423) 
  reac_source_local(30,423) = + reac_rate_local(423) 
  reac_source_local(44,423) = - reac_rate_local(423) 
  reac_source_local(08,424) = - reac_rate_local(424) 
  reac_source_local(14,424) = + reac_rate_local(424) 
  reac_source_local(21,424) = + reac_rate_local(424) 
  reac_source_local(44,424) = - reac_rate_local(424) 
  reac_source_local(08,425) = - reac_rate_local(425) 
  reac_source_local(14,425) = + reac_rate_local(425) 
  reac_source_local(21,425) = - reac_rate_local(425) 
  reac_source_local(26,425) = + reac_rate_local(425) 
  reac_source_local(08,426) = - reac_rate_local(426) 
  reac_source_local(14,426) = + reac_rate_local(426) 
  reac_source_local(23,426) = - reac_rate_local(426) 
  reac_source_local(34,426) = + reac_rate_local(426) 
  reac_source_local(08,427) = - reac_rate_local(427) 
  reac_source_local(11,427) = + reac_rate_local(427) * 2.d0
  reac_source_local(34,427) = - reac_rate_local(427) 
  reac_source_local(08,428) = - reac_rate_local(428) 
  reac_source_local(23,428) = + reac_rate_local(428) 
  reac_source_local(30,428) = + reac_rate_local(428) 
  reac_source_local(34,428) = - reac_rate_local(428) 
  reac_source_local(08,429) = - reac_rate_local(429) 
  reac_source_local(14,429) = + reac_rate_local(429) 
  reac_source_local(34,429) = - reac_rate_local(429) 
  reac_source_local(08,430) = - reac_rate_local(430) 
  reac_source_local(30,430) = + reac_rate_local(430) 
  reac_source_local(34,430) = + reac_rate_local(430) 
  reac_source_local(08,431) = - reac_rate_local(431) * 2.d0
  reac_source_local(25,431) = + reac_rate_local(431) 
  reac_source_local(08,432) = - reac_rate_local(432) 
  reac_source_local(25,432) = + reac_rate_local(432) 
  reac_source_local(30,432) = + reac_rate_local(432) 
  reac_source_local(31,432) = - reac_rate_local(432) 
  reac_source_local(08,433) = - reac_rate_local(433) 
  reac_source_local(30,433) = + reac_rate_local(433) * 2.d0
  reac_source_local(42,433) = - reac_rate_local(433) 
  reac_source_local(23,434) = + reac_rate_local(434) 
  reac_source_local(30,434) = - reac_rate_local(434) 
  reac_source_local(34,434) = - reac_rate_local(434) 
  reac_source_local(42,434) = + reac_rate_local(434) 
  reac_source_local(08,435) = + reac_rate_local(435) 
  reac_source_local(30,435) = - reac_rate_local(435) 
  reac_source_local(34,435) = - reac_rate_local(435) 
  reac_source_local(08,436) = + reac_rate_local(436) 
  reac_source_local(23,436) = - reac_rate_local(436) 
  reac_source_local(30,436) = - reac_rate_local(436) 
  reac_source_local(34,436) = + reac_rate_local(436) 
  reac_source_local(30,437) = - reac_rate_local(437) 
  reac_source_local(34,437) = + reac_rate_local(437) 
  reac_source_local(42,437) = + reac_rate_local(437) 
  reac_source_local(08,438) = + reac_rate_local(438) 
  reac_source_local(21,438) = - reac_rate_local(438) 
  reac_source_local(26,438) = + reac_rate_local(438) 
  reac_source_local(30,438) = - reac_rate_local(438) 
  reac_source_local(30,439) = - reac_rate_local(439) 
  reac_source_local(32,439) = - reac_rate_local(439) 
  reac_source_local(42,439) = + reac_rate_local(439) * 2.d0
  reac_source_local(21,440) = - reac_rate_local(440) 
  reac_source_local(30,440) = - reac_rate_local(440) 
  reac_source_local(42,440) = + reac_rate_local(440) 
  reac_source_local(44,440) = + reac_rate_local(440) 
  reac_source_local(08,441) = + reac_rate_local(441) 
  reac_source_local(30,441) = - reac_rate_local(441) * 2.d0
  reac_source_local(42,441) = + reac_rate_local(441) 
  reac_source_local(11,442) = - reac_rate_local(442) 
  reac_source_local(30,442) = - reac_rate_local(442) 
  reac_source_local(44,442) = + reac_rate_local(442) 
  reac_source_local(23,443) = + reac_rate_local(443) 
  reac_source_local(30,443) = - reac_rate_local(443) 
  reac_source_local(32,443) = + reac_rate_local(443) 
  reac_source_local(08,444) = - reac_rate_local(444) 
  reac_source_local(30,444) = - reac_rate_local(444) 
  reac_source_local(31,444) = + reac_rate_local(444) 
  reac_source_local(14,445) = + reac_rate_local(445) 
  reac_source_local(23,445) = - reac_rate_local(445) 
  reac_source_local(30,445) = - reac_rate_local(445) 
  reac_source_local(21,446) = + reac_rate_local(446) 
  reac_source_local(30,446) = - reac_rate_local(446) 
  reac_source_local(47,446) = - reac_rate_local(446) 
  reac_source_local(21,447) = + reac_rate_local(447) 
  reac_source_local(30,447) = - reac_rate_local(447) 
  reac_source_local(31,447) = - reac_rate_local(447) 
  reac_source_local(44,447) = + reac_rate_local(447) 
  reac_source_local(30,448) = + reac_rate_local(448) 
  reac_source_local(32,448) = + reac_rate_local(448) 
  reac_source_local(42,448) = - reac_rate_local(448) * 2.d0
  reac_source_local(07,449) = - reac_rate_local(449) 
  reac_source_local(30,449) = + reac_rate_local(449) 
  reac_source_local(42,449) = - reac_rate_local(449) 
  reac_source_local(44,449) = + reac_rate_local(449) 
  reac_source_local(07,450) = + reac_rate_local(450) 
  reac_source_local(32,450) = + reac_rate_local(450) 
  reac_source_local(42,450) = - reac_rate_local(450) 
  reac_source_local(44,450) = - reac_rate_local(450) 
  reac_source_local(21,451) = + reac_rate_local(451) 
  reac_source_local(30,451) = + reac_rate_local(451) 
  reac_source_local(42,451) = - reac_rate_local(451) 
  reac_source_local(44,451) = - reac_rate_local(451) 
  reac_source_local(21,452) = - reac_rate_local(452) 
  reac_source_local(26,452) = + reac_rate_local(452) 
  reac_source_local(30,452) = + reac_rate_local(452) 
  reac_source_local(42,452) = - reac_rate_local(452) 
  reac_source_local(21,453) = + reac_rate_local(453) 
  reac_source_local(26,453) = - reac_rate_local(453) 
  reac_source_local(32,453) = + reac_rate_local(453) 
  reac_source_local(42,453) = - reac_rate_local(453) 
  reac_source_local(23,454) = - reac_rate_local(454) 
  reac_source_local(30,454) = + reac_rate_local(454) 
  reac_source_local(34,454) = + reac_rate_local(454) 
  reac_source_local(42,454) = - reac_rate_local(454) 
  reac_source_local(23,455) = + reac_rate_local(455) 
  reac_source_local(32,455) = + reac_rate_local(455) 
  reac_source_local(34,455) = - reac_rate_local(455) 
  reac_source_local(42,455) = - reac_rate_local(455) 
  reac_source_local(30,456) = + reac_rate_local(456) 
  reac_source_local(34,456) = - reac_rate_local(456) 
  reac_source_local(42,456) = - reac_rate_local(456) 
  reac_source_local(32,457) = + reac_rate_local(457) 
  reac_source_local(34,457) = + reac_rate_local(457) 
  reac_source_local(42,457) = - reac_rate_local(457) 
  reac_source_local(25,458) = + reac_rate_local(458) 
  reac_source_local(31,458) = - reac_rate_local(458) 
  reac_source_local(32,458) = + reac_rate_local(458) 
  reac_source_local(42,458) = - reac_rate_local(458) 
  reac_source_local(32,459) = - reac_rate_local(459) 
  reac_source_local(34,459) = - reac_rate_local(459) 
  reac_source_local(42,459) = + reac_rate_local(459) 
  reac_source_local(23,460) = - reac_rate_local(460) 
  reac_source_local(30,460) = + reac_rate_local(460) 
  reac_source_local(32,460) = - reac_rate_local(460) 
  reac_source_local(23,461) = - reac_rate_local(461) 
  reac_source_local(32,461) = - reac_rate_local(461) 
  reac_source_local(34,461) = + reac_rate_local(461) 
  reac_source_local(42,461) = + reac_rate_local(461) 
  reac_source_local(11,462) = - reac_rate_local(462) 
  reac_source_local(26,462) = + reac_rate_local(462) 
  reac_source_local(32,462) = - reac_rate_local(462) 
  reac_source_local(21,463) = + reac_rate_local(463) 
  reac_source_local(26,463) = + reac_rate_local(463) 
  reac_source_local(31,463) = - reac_rate_local(463) 
  reac_source_local(32,463) = - reac_rate_local(463) 
  reac_source_local(07,464) = - reac_rate_local(464) 
  reac_source_local(21,464) = + reac_rate_local(464) 
  reac_source_local(26,464) = - reac_rate_local(464) 
  reac_source_local(44,464) = + reac_rate_local(464) 
  reac_source_local(07,465) = - reac_rate_local(465) 
  reac_source_local(23,465) = + reac_rate_local(465) 
  reac_source_local(34,465) = - reac_rate_local(465) 
  reac_source_local(44,465) = + reac_rate_local(465) 
  reac_source_local(07,466) = - reac_rate_local(466) 
  reac_source_local(08,466) = + reac_rate_local(466) 
  reac_source_local(11,466) = + reac_rate_local(466) 
  reac_source_local(07,467) = - reac_rate_local(467) 
  reac_source_local(25,467) = + reac_rate_local(467) 
  reac_source_local(31,467) = - reac_rate_local(467) 
  reac_source_local(44,467) = + reac_rate_local(467) 
  reac_source_local(07,468) = - reac_rate_local(468) 
  reac_source_local(25,468) = + reac_rate_local(468) 
  reac_source_local(47,468) = - reac_rate_local(468) 
  reac_source_local(07,469) = + reac_rate_local(469) 
  reac_source_local(21,469) = + reac_rate_local(469) 
  reac_source_local(44,469) = - reac_rate_local(469) * 2.d0
  reac_source_local(07,470) = + reac_rate_local(470) 
  reac_source_local(21,470) = - reac_rate_local(470) 
  reac_source_local(26,470) = + reac_rate_local(470) 
  reac_source_local(44,470) = - reac_rate_local(470) 
  reac_source_local(21,471) = + reac_rate_local(471) * 2.d0
  reac_source_local(26,471) = - reac_rate_local(471) 
  reac_source_local(44,471) = - reac_rate_local(471) 
  reac_source_local(07,472) = + reac_rate_local(472) 
  reac_source_local(23,472) = - reac_rate_local(472) 
  reac_source_local(34,472) = + reac_rate_local(472) 
  reac_source_local(44,472) = - reac_rate_local(472) 
  reac_source_local(21,473) = + reac_rate_local(473) 
  reac_source_local(23,473) = + reac_rate_local(473) 
  reac_source_local(34,473) = - reac_rate_local(473) 
  reac_source_local(44,473) = - reac_rate_local(473) 
  reac_source_local(07,474) = + reac_rate_local(474) 
  reac_source_local(34,474) = - reac_rate_local(474) 
  reac_source_local(44,474) = - reac_rate_local(474) 
  reac_source_local(08,475) = + reac_rate_local(475) 
  reac_source_local(11,475) = + reac_rate_local(475) 
  reac_source_local(34,475) = - reac_rate_local(475) 
  reac_source_local(44,475) = - reac_rate_local(475) 
  reac_source_local(21,476) = + reac_rate_local(476) 
  reac_source_local(34,476) = + reac_rate_local(476) 
  reac_source_local(44,476) = - reac_rate_local(476) 
  reac_source_local(11,477) = + reac_rate_local(477) 
  reac_source_local(30,477) = + reac_rate_local(477) 
  reac_source_local(44,477) = - reac_rate_local(477) 
  reac_source_local(21,478) = + reac_rate_local(478) 
  reac_source_local(25,478) = + reac_rate_local(478) 
  reac_source_local(31,478) = - reac_rate_local(478) 
  reac_source_local(44,478) = - reac_rate_local(478) 
  reac_source_local(21,479) = - reac_rate_local(479) 
  reac_source_local(26,479) = + reac_rate_local(479) 
  reac_source_local(32,479) = - reac_rate_local(479) 
  reac_source_local(42,479) = + reac_rate_local(479) 
  reac_source_local(21,480) = - reac_rate_local(480) * 2.d0
  reac_source_local(26,480) = + reac_rate_local(480) 
  reac_source_local(44,480) = + reac_rate_local(480) 
  reac_source_local(21,481) = - reac_rate_local(481) 
  reac_source_local(26,481) = + reac_rate_local(481) 
  reac_source_local(34,481) = + reac_rate_local(481) 
  reac_source_local(21,482) = - reac_rate_local(482) 
  reac_source_local(23,482) = + reac_rate_local(482) 
  reac_source_local(26,482) = + reac_rate_local(482) 
  reac_source_local(34,482) = - reac_rate_local(482) 
  reac_source_local(21,483) = - reac_rate_local(483) 
  reac_source_local(34,483) = - reac_rate_local(483) 
  reac_source_local(44,483) = + reac_rate_local(483) 
  reac_source_local(11,484) = + reac_rate_local(484) 
  reac_source_local(21,484) = - reac_rate_local(484) 
  reac_source_local(42,484) = + reac_rate_local(484) 
  reac_source_local(11,485) = - reac_rate_local(485) 
  reac_source_local(21,485) = - reac_rate_local(485) 
  reac_source_local(31,485) = + reac_rate_local(485) 
  reac_source_local(21,486) = - reac_rate_local(486) 
  reac_source_local(25,486) = + reac_rate_local(486) 
  reac_source_local(26,486) = + reac_rate_local(486) 
  reac_source_local(31,486) = - reac_rate_local(486) 
  reac_source_local(21,487) = + reac_rate_local(487) 
  reac_source_local(23,487) = - reac_rate_local(487) 
  reac_source_local(26,487) = - reac_rate_local(487) 
  reac_source_local(34,487) = + reac_rate_local(487) 
  reac_source_local(21,488) = + reac_rate_local(488) 
  reac_source_local(26,488) = - reac_rate_local(488) 
  reac_source_local(34,488) = - reac_rate_local(488) 
  reac_source_local(11,489) = + reac_rate_local(489) 
  reac_source_local(26,489) = - reac_rate_local(489) 
  reac_source_local(32,489) = + reac_rate_local(489) 
  reac_source_local(08,490) = + reac_rate_local(490) 
  reac_source_local(30,490) = + reac_rate_local(490) 
  reac_source_local(31,490) = - reac_rate_local(490) 
  reac_source_local(30,491) = + reac_rate_local(491) 
  reac_source_local(31,491) = - reac_rate_local(491) 
  reac_source_local(44,491) = + reac_rate_local(491) 
  reac_source_local(47,491) = - reac_rate_local(491) 
  reac_source_local(11,492) = + reac_rate_local(492) 
  reac_source_local(21,492) = + reac_rate_local(492) 
  reac_source_local(31,492) = - reac_rate_local(492) 
  reac_source_local(23,493) = - reac_rate_local(493) 
  reac_source_local(25,493) = + reac_rate_local(493) 
  reac_source_local(31,493) = - reac_rate_local(493) 
  reac_source_local(34,493) = + reac_rate_local(493) 
  reac_source_local(06,494) = + reac_rate_local(494) 
  reac_source_local(11,494) = - reac_rate_local(494) 
  reac_source_local(25,494) = - reac_rate_local(494) 
  reac_source_local(31,494) = + reac_rate_local(494) 
  reac_source_local(11,495) = + reac_rate_local(495) 
  reac_source_local(25,495) = - reac_rate_local(495) 
  reac_source_local(44,495) = + reac_rate_local(495) 
  reac_source_local(08,496) = + reac_rate_local(496) * 2.d0
  reac_source_local(25,496) = - reac_rate_local(496) 
  reac_source_local(23,497) = + reac_rate_local(497) 
  reac_source_local(25,497) = - reac_rate_local(497) 
  reac_source_local(31,497) = + reac_rate_local(497) 
  reac_source_local(34,497) = - reac_rate_local(497) 
  reac_source_local(11,498) = + reac_rate_local(498) 
  reac_source_local(25,498) = - reac_rate_local(498) 
  reac_source_local(31,498) = + reac_rate_local(498) 
  reac_source_local(47,498) = - reac_rate_local(498) 
  reac_source_local(25,499) = - reac_rate_local(499) 
  reac_source_local(30,499) = + reac_rate_local(499) 
  reac_source_local(31,499) = + reac_rate_local(499) 
  reac_source_local(42,499) = - reac_rate_local(499) 
  reac_source_local(07,500) = + reac_rate_local(500) 
  reac_source_local(25,500) = - reac_rate_local(500) 
  reac_source_local(31,500) = + reac_rate_local(500) 
  reac_source_local(44,500) = - reac_rate_local(500) 
  reac_source_local(25,501) = - reac_rate_local(501) 
  reac_source_local(31,501) = + reac_rate_local(501) 
  reac_source_local(32,501) = + reac_rate_local(501) 
  reac_source_local(33,501) = - reac_rate_local(501) 
  reac_source_local(08,502) = - reac_rate_local(502) 
  reac_source_local(14,502) = + reac_rate_local(502) 
  reac_source_local(25,502) = - reac_rate_local(502) 
  reac_source_local(31,502) = + reac_rate_local(502) 
  reac_source_local(21,503) = + reac_rate_local(503) 
  reac_source_local(25,503) = - reac_rate_local(503) 
  reac_source_local(26,503) = - reac_rate_local(503) 
  reac_source_local(31,503) = + reac_rate_local(503) 
  reac_source_local(22,504) = + reac_rate_local(504) 
  reac_source_local(25,504) = - reac_rate_local(504) 
  reac_source_local(47,504) = - reac_rate_local(504) 
  reac_source_local(11,505) = + reac_rate_local(505) 
  reac_source_local(22,505) = - reac_rate_local(505) 
  reac_source_local(31,505) = + reac_rate_local(505) 
  reac_source_local(23,506) = - reac_rate_local(506) 
  reac_source_local(34,506) = + reac_rate_local(506) * 2.d0
  reac_source_local(23,507) = + reac_rate_local(507) 
  reac_source_local(34,507) = - reac_rate_local(507) * 2.d0
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
  rrt(001) = rrt(001) * density(01) * density(06) 
  rrt(002) = rrt(002) * density(01) * density(06) 
  rrt(003) = rrt(003) * density(01) * density(14) 
  rrt(004) = rrt(004) * density(01) * density(14) 
  rrt(005) = rrt(005) * density(01) * density(30) 
  rrt(006) = rrt(006) * density(01) * density(30) 
  rrt(007) = rrt(007) * density(01) * density(32) 
  rrt(008) = rrt(008) * density(01) * density(32) 
  rrt(009) = rrt(009) * density(01) * density(32) 
  rrt(010) = rrt(010) * density(01) * density(07) 
  rrt(011) = rrt(011) * density(01) * density(07) 
  rrt(012) = rrt(012) * density(01) * density(21) 
  rrt(013) = rrt(013) * density(01) * density(06) 
  rrt(014) = rrt(014) * density(01) * density(06) 
  rrt(015) = rrt(015) * density(01) * density(06) 
  rrt(016) = rrt(016) * density(01) * density(11) 
  rrt(017) = rrt(017) * density(01) * density(11) 
  rrt(018) = rrt(018) * density(01) * density(47) 
  rrt(019) = rrt(019) * density(01) * density(06) 
  rrt(020) = rrt(020) * density(01) * density(06) 
  rrt(021) = rrt(021) * density(01) * density(06) 
  rrt(022) = rrt(022) * density(01) * density(06) 
  rrt(023) = rrt(023) * density(01) * density(11) 
  rrt(024) = rrt(024) * density(01) * density(11) 
  rrt(025) = rrt(025) * density(01) * density(11) 
  rrt(026) = rrt(026) * density(01) * density(47) 
  rrt(027) = rrt(027) * density(01) * density(47) 
  rrt(028) = rrt(028) * density(01) * density(43) 
  rrt(029) = rrt(029) * density(01) * density(14) 
  rrt(030) = rrt(030) * density(01) * density(14) 
  rrt(031) = rrt(031) * density(01) * density(14) 
  rrt(032) = rrt(032) * density(01) * density(14) 
  rrt(033) = rrt(033) * density(01) * density(14) 
  rrt(034) = rrt(034) * density(01) * density(14) 
  rrt(035) = rrt(035) * density(01) * density(08) 
  rrt(036) = rrt(036) * density(01) * density(08) 
  rrt(037) = rrt(037) * density(01) * density(08) 
  rrt(038) = rrt(038) * density(01) * density(08) 
  rrt(039) = rrt(039) * density(01) * density(08) 
  rrt(040) = rrt(040) * density(01) * density(08) 
  rrt(041) = rrt(041) * density(01) * density(30) 
  rrt(042) = rrt(042) * density(01) * density(30) 
  rrt(043) = rrt(043) * density(01) * density(30) 
  rrt(044) = rrt(044) * density(01) * density(30) 
  rrt(045) = rrt(045) * density(01) * density(30) 
  rrt(046) = rrt(046) * density(01) * density(42) 
  rrt(047) = rrt(047) * density(01) * density(42) 
  rrt(048) = rrt(048) * density(01) * density(32) 
  rrt(049) = rrt(049) * density(01) * density(07) 
  rrt(050) = rrt(050) * density(01) * density(07) 
  rrt(051) = rrt(051) * density(01) * density(07) 
  rrt(052) = rrt(052) * density(01) * density(07) 
  rrt(053) = rrt(053) * density(01) * density(07) 
  rrt(054) = rrt(054) * density(01) * density(07) 
  rrt(055) = rrt(055) * density(01) * density(44) 
  rrt(056) = rrt(056) * density(01) * density(44) 
  rrt(057) = rrt(057) * density(01) * density(44) 
  rrt(058) = rrt(058) * density(01) * density(44) 
  rrt(059) = rrt(059) * density(01) * density(44) 
  rrt(060) = rrt(060) * density(01) * density(21) 
  rrt(061) = rrt(061) * density(01) * density(21) 
  rrt(062) = rrt(062) * density(01) * density(21) 
  rrt(063) = rrt(063) * density(01) * density(21) 
  rrt(064) = rrt(064) * density(01) * density(21) 
  rrt(065) = rrt(065) * density(01) * density(26) 
  rrt(066) = rrt(066) * density(01) * density(26) 
  rrt(067) = rrt(067) * density(01) * density(17) 
  rrt(068) = rrt(068) * density(01) * density(17) 
  rrt(069) = rrt(069) * density(01) * density(14) 
  rrt(070) = rrt(070) * density(01) * density(14) 
  rrt(071) = rrt(071) * density(01) * density(14) 
  rrt(072) = rrt(072) * density(01) * density(14) 
  rrt(073) = rrt(073) * density(01) * density(14) 
  rrt(074) = rrt(074) * density(01) * density(14) 
  rrt(075) = rrt(075) * density(01) * density(14) 
  rrt(076) = rrt(076) * density(01) * density(08) 
  rrt(077) = rrt(077) * density(01) * density(08) 
  rrt(078) = rrt(078) * density(01) * density(08) 
  rrt(079) = rrt(079) * density(01) * density(08) 
  rrt(080) = rrt(080) * density(01) * density(08) 
  rrt(081) = rrt(081) * density(01) * density(08) 
  rrt(082) = rrt(082) * density(01) * density(08) 
  rrt(083) = rrt(083) * density(01) * density(30) 
  rrt(084) = rrt(084) * density(01) * density(30) 
  rrt(085) = rrt(085) * density(01) * density(30) 
  rrt(086) = rrt(086) * density(01) * density(30) 
  rrt(087) = rrt(087) * density(01) * density(30) 
  rrt(088) = rrt(088) * density(01) * density(42) 
  rrt(089) = rrt(089) * density(01) * density(42) 
  rrt(090) = rrt(090) * density(01) * density(42) 
  rrt(091) = rrt(091) * density(01) * density(42) 
  rrt(092) = rrt(092) * density(01) * density(42) 
  rrt(093) = rrt(093) * density(01) * density(32) 
  rrt(094) = rrt(094) * density(01) * density(32) 
  rrt(095) = rrt(095) * density(01) * density(07) 
  rrt(096) = rrt(096) * density(01) * density(07) 
  rrt(097) = rrt(097) * density(01) * density(07) 
  rrt(098) = rrt(098) * density(01) * density(07) 
  rrt(099) = rrt(099) * density(01) * density(07) 
  rrt(100) = rrt(100) * density(01) * density(07) 
  rrt(101) = rrt(101) * density(01) * density(07) 
  rrt(102) = rrt(102) * density(01) * density(07) 
  rrt(103) = rrt(103) * density(01) * density(07) 
  rrt(104) = rrt(104) * density(01) * density(44) 
  rrt(105) = rrt(105) * density(01) * density(44) 
  rrt(106) = rrt(106) * density(01) * density(44) 
  rrt(107) = rrt(107) * density(01) * density(44) 
  rrt(108) = rrt(108) * density(01) * density(44) 
  rrt(109) = rrt(109) * density(01) * density(44) 
  rrt(110) = rrt(110) * density(01) * density(44) 
  rrt(111) = rrt(111) * density(01) * density(44) 
  rrt(112) = rrt(112) * density(01) * density(44) 
  rrt(113) = rrt(113) * density(01) * density(44) 
  rrt(114) = rrt(114) * density(01) * density(44) 
  rrt(115) = rrt(115) * density(01) * density(21) 
  rrt(116) = rrt(116) * density(01) * density(21) 
  rrt(117) = rrt(117) * density(01) * density(21) 
  rrt(118) = rrt(118) * density(01) * density(21) 
  rrt(119) = rrt(119) * density(01) * density(21) 
  rrt(120) = rrt(120) * density(01) * density(21) 
  rrt(121) = rrt(121) * density(01) * density(21) 
  rrt(122) = rrt(122) * density(01) * density(21) 
  rrt(123) = rrt(123) * density(01) * density(21) 
  rrt(124) = rrt(124) * density(01) * density(21) 
  rrt(125) = rrt(125) * density(01) * density(21) 
  rrt(126) = rrt(126) * density(01) * density(26) 
  rrt(127) = rrt(127) * density(01) * density(26) 
  rrt(128) = rrt(128) * density(01) * density(26) 
  rrt(129) = rrt(129) * density(01) * density(26) 
  rrt(130) = rrt(130) * density(01) * density(26) 
  rrt(131) = rrt(131) * density(01) * density(26) 
  rrt(132) = rrt(132) * density(01) * density(26) 
  rrt(133) = rrt(133) * density(01) * density(26) 
  rrt(134) = rrt(134) * density(01) * density(17) 
  rrt(135) = rrt(135) * density(01) * density(17) 
  rrt(136) = rrt(136) * density(01) * density(17) 
  rrt(137) = rrt(137) * density(01) * density(17) 
  rrt(138) = rrt(138) * density(01) * density(17) 
  rrt(139) = rrt(139) * density(01) * density(19) 
  rrt(140) = rrt(140) * density(01) * density(51) 
  rrt(141) = rrt(141) * density(01) * density(19) 
  rrt(142) = rrt(142) * density(01) * density(51) 
  rrt(143) = rrt(143) * density(01) * density(19) 
  rrt(144) = rrt(144) * density(01) * density(51) 
  rrt(145) = rrt(145) * density(01) * density(19) 
  rrt(146) = rrt(146) * density(01) * density(51) 
  rrt(147) = rrt(147) * density(01) * density(19) 
  rrt(148) = rrt(148) * density(01) * density(51) 
  rrt(149) = rrt(149) * density(01) * density(19) 
  rrt(150) = rrt(150) * density(01) * density(51) 
  rrt(151) = rrt(151) * density(01) * density(19) 
  rrt(152) = rrt(152) * density(01) * density(51) 
  rrt(153) = rrt(153) * density(01) * density(03) 
  rrt(154) = rrt(154) * density(01) * density(29) 
  rrt(155) = rrt(155) * density(01) * density(03) 
  rrt(156) = rrt(156) * density(01) * density(29) 
  rrt(157) = rrt(157) * density(01) * density(03) 
  rrt(158) = rrt(158) * density(01) * density(29) 
  rrt(159) = rrt(159) * density(01) * density(03) 
  rrt(160) = rrt(160) * density(01) * density(29) 
  rrt(161) = rrt(161) * density(01) * density(03) 
  rrt(162) = rrt(162) * density(01) * density(29) 
  rrt(163) = rrt(163) * density(01) * density(03) 
  rrt(164) = rrt(164) * density(01) * density(29) 
  rrt(165) = rrt(165) * density(01) * density(05) 
  rrt(166) = rrt(166) * density(01) * density(15) 
  rrt(167) = rrt(167) * density(01) * density(05) 
  rrt(168) = rrt(168) * density(01) * density(15) 
  rrt(169) = rrt(169) * density(01) * density(05) 
  rrt(170) = rrt(170) * density(01) * density(15) 
  rrt(171) = rrt(171) * density(01) * density(05) 
  rrt(172) = rrt(172) * density(01) * density(15) 
  rrt(173) = rrt(173) * density(01) * density(05) 
  rrt(174) = rrt(174) * density(01) * density(15) 
  rrt(175) = rrt(175) * density(01) * density(50) 
  rrt(176) = rrt(176) * density(01) * density(37) 
  rrt(177) = rrt(177) * density(01) * density(40) 
  rrt(178) = rrt(178) * density(01) * density(24) 
  rrt(179) = rrt(179) * density(01) * density(02) 
  rrt(180) = rrt(180) * density(01) * density(24) 
  rrt(181) = rrt(181) * density(01) * density(02) 
  rrt(182) = rrt(182) * density(01) * density(24) 
  rrt(183) = rrt(183) * density(01) * density(02) 
  rrt(184) = rrt(184) * density(01) * density(24) 
  rrt(185) = rrt(185) * density(01) * density(02) 
  rrt(186) = rrt(186) * density(01) * density(24) 
  rrt(187) = rrt(187) * density(01) * density(02) 
  rrt(188) = rrt(188) * density(01) * density(24) 
  rrt(189) = rrt(189) * density(01) * density(02) 
  rrt(190) = rrt(190) * density(01) * density(03) 
  rrt(191) = rrt(191) * density(01) * density(29) 
  rrt(192) = rrt(192) * density(01) * density(03) 
  rrt(193) = rrt(193) * density(01) * density(29) 
  rrt(194) = rrt(194) * density(01) * density(03) 
  rrt(195) = rrt(195) * density(01) * density(29) 
  rrt(196) = rrt(196) * density(01) * density(03) 
  rrt(197) = rrt(197) * density(01) * density(29) 
  rrt(198) = rrt(198) * density(01) * density(03) 
  rrt(199) = rrt(199) * density(01) * density(29) 
  rrt(200) = rrt(200) * density(01) * density(03) 
  rrt(201) = rrt(201) * density(01) * density(29) 
  rrt(202) = rrt(202) * density(01) * density(03) 
  rrt(203) = rrt(203) * density(01) * density(29) 
  rrt(204) = rrt(204) * density(01) * density(05) 
  rrt(205) = rrt(205) * density(01) * density(15) 
  rrt(206) = rrt(206) * density(01) * density(05) 
  rrt(207) = rrt(207) * density(01) * density(15) 
  rrt(208) = rrt(208) * density(01) * density(05) 
  rrt(209) = rrt(209) * density(01) * density(15) 
  rrt(210) = rrt(210) * density(01) * density(05) 
  rrt(211) = rrt(211) * density(01) * density(15) 
  rrt(212) = rrt(212) * density(01) * density(05) 
  rrt(213) = rrt(213) * density(01) * density(15) 
  rrt(214) = rrt(214) * density(01) * density(50) 
  rrt(215) = rrt(215) * density(01) * density(37) 
  rrt(216) = rrt(216) * density(01) * density(40) 
  rrt(217) = rrt(217) * density(01) * density(50) 
  rrt(218) = rrt(218) * density(01) * density(37) 
  rrt(219) = rrt(219) * density(01) * density(40) 
  rrt(220) = rrt(220) * density(01) * density(24) 
  rrt(221) = rrt(221) * density(01) * density(02) 
  rrt(222) = rrt(222) * density(01) * density(24) 
  rrt(223) = rrt(223) * density(01) * density(02) 
  rrt(224) = rrt(224) * density(01) * density(24) 
  rrt(225) = rrt(225) * density(01) * density(02) 
  rrt(226) = rrt(226) * density(01) * density(24) 
  rrt(227) = rrt(227) * density(01) * density(02) 
  rrt(228) = rrt(228) * density(01) * density(24) 
  rrt(229) = rrt(229) * density(01) * density(02) 
  rrt(230) = rrt(230) * density(01) * density(24) 
  rrt(231) = rrt(231) * density(01) * density(02) 
  rrt(232) = rrt(232) * density(01) * density(24) 
  rrt(233) = rrt(233) * density(01) * density(02) 
  rrt(234) = rrt(234) * density(01) * density(24) 
  rrt(235) = rrt(235) * density(01) * density(02) 
  rrt(236) = rrt(236) * density(01) * density(24) 
  rrt(237) = rrt(237) * density(01) * density(02) 
  rrt(238) = rrt(238) * density(01) * density(23) 
  rrt(239) = rrt(239) * density(01) * density(23) 
  rrt(240) = rrt(240) * density(01)**2 * density(48) 
  rrt(241) = rrt(241) * density(01)**2 * density(49) 
  rrt(242) = rrt(242) * density(01)**2 * density(46) 
  rrt(243) = rrt(243) * density(01)**2 * density(13) 
  rrt(244) = rrt(244) * density(01)**2 * density(52) 
  rrt(245) = rrt(245) * density(01)**2 * density(27) 
  rrt(246) = rrt(246) * density(01)**2 * density(16) 
  rrt(247) = rrt(247) * density(01)**2 * density(04) 
  rrt(248) = rrt(248) * density(01)**2 * density(09) 
  rrt(249) = rrt(249) * density(01)**2 * density(41) 
  rrt(250) = rrt(250) * density(01)**2 * density(18) 
  rrt(251) = rrt(251) * density(01)**2 * density(38) 
  rrt(252) = rrt(252) * density(01)**2 * density(39) 
  rrt(253) = rrt(253) * density(01)**2 * density(20) 
  rrt(254) = rrt(254) * density(01)**2 * density(12) 
  rrt(255) = rrt(255) * density(01)**2 * density(45) 
  rrt(256) = rrt(256) * density(01) * density(35) 
  rrt(257) = rrt(257) * density(01) * density(35) 
  rrt(258) = rrt(258) * density(01) * density(48) 
  rrt(259) = rrt(259) * density(01) * density(48) 
  rrt(260) = rrt(260) * density(01) * density(48) 
  rrt(261) = rrt(261) * density(01) * density(49) 
  rrt(262) = rrt(262) * density(01) * density(49) 
  rrt(263) = rrt(263) * density(01) * density(46) 
  rrt(264) = rrt(264) * density(01) * density(52) 
  rrt(265) = rrt(265) * density(01) * density(52) 
  rrt(266) = rrt(266) * density(01) * density(27) 
  rrt(267) = rrt(267) * density(01) * density(27) 
  rrt(268) = rrt(268) * density(01) * density(27) 
  rrt(269) = rrt(269) * density(01) * density(27) 
  rrt(270) = rrt(270) * density(01) * density(27) 
  rrt(271) = rrt(271) * density(01) * density(16) 
  rrt(272) = rrt(272) * density(01) * density(16) 
  rrt(273) = rrt(273) * density(01) * density(04) 
  rrt(274) = rrt(274) * density(01) * density(09) 
  rrt(275) = rrt(275) * density(01) * density(45) 
  rrt(276) = rrt(276) * density(01) * density(28) 
  rrt(277) = rrt(277) * density(23) * density(45) 
  rrt(278) = rrt(278) * density(35) * density(47) 
  rrt(279) = rrt(279) * density(35) * density(43) 
  rrt(280) = rrt(280) * density(14) * density(35) 
  rrt(281) = rrt(281) * density(30) * density(35) 
  rrt(282) = rrt(282) * density(32) * density(35) 
  rrt(283) = rrt(283) * density(34) * density(35) 
  rrt(284) = rrt(284) * density(06) * density(48) 
  rrt(285) = rrt(285) * density(14) * density(48) 
  rrt(286) = rrt(286) * density(30) * density(48) 
  rrt(287) = rrt(287) * density(30) * density(48) 
  rrt(288) = rrt(288) * density(32) * density(48) 
  rrt(289) = rrt(289) * density(32) * density(48) 
  rrt(290) = rrt(290) * density(23) * density(48) 
  rrt(291) = rrt(291) * density(34) * density(48) 
  rrt(292) = rrt(292) * density(06) * density(49) 
  rrt(293) = rrt(293) * density(06) * density(49) 
  rrt(294) = rrt(294) * density(47) * density(49) 
  rrt(295) = rrt(295) * density(43) * density(49) 
  rrt(296) = rrt(296) * density(14) * density(49) 
  rrt(297) = rrt(297) * density(30) * density(49) 
  rrt(298) = rrt(298) * density(42) * density(49) 
  rrt(299) = rrt(299) * density(06) * density(46) 
  rrt(300) = rrt(300) * density(06) * density(46) 
  rrt(301) = rrt(301) * density(06) * density(46) 
  rrt(302) = rrt(302) * density(06) * density(46) 
  rrt(303) = rrt(303) * density(06) * density(46) 
  rrt(304) = rrt(304) * density(23) * density(46) 
  rrt(305) = rrt(305) * density(06) * density(13) 
  rrt(306) = rrt(306) * density(06) * density(13) 
  rrt(307) = rrt(307) * density(06) * density(13) 
  rrt(308) = rrt(308) * density(13) * density(23) 
  rrt(309) = rrt(309) * density(30) * density(52) 
  rrt(310) = rrt(310) * density(32) * density(52) 
  rrt(311) = rrt(311) * density(34) * density(52) 
  rrt(312) = rrt(312) * density(27) * density(34) 
  rrt(313) = rrt(313) * density(16) * density(42) 
  rrt(314) = rrt(314) * density(16) * density(42) 
  rrt(315) = rrt(315) * density(16) * density(34) 
  rrt(316) = rrt(316) * density(04) * density(14) 
  rrt(317) = rrt(317) * density(04) * density(30) 
  rrt(318) = rrt(318) * density(04) * density(33) 
  rrt(319) = rrt(319) * density(04) * density(34) 
  rrt(320) = rrt(320) * density(06) * density(09) 
  rrt(321) = rrt(321) * density(09) * density(14) 
  rrt(322) = rrt(322) * density(09) * density(14) 
  rrt(323) = rrt(323) * density(09) * density(30) 
  rrt(324) = rrt(324) * density(09) * density(42) 
  rrt(325) = rrt(325) * density(09) * density(23) 
  rrt(326) = rrt(326) * density(06) * density(36) 
  rrt(327) = rrt(327) * density(06) * density(28) 
  rrt(328) = rrt(328) * density(11) * density(28) 
  rrt(329) = rrt(329) * density(28) * density(47) 
  rrt(330) = rrt(330) * density(28) * density(43) 
  rrt(331) = rrt(331) * density(14) * density(28) 
  rrt(332) = rrt(332) * density(08) * density(28) 
  rrt(333) = rrt(333) * density(28) * density(30) 
  rrt(334) = rrt(334) * density(28) * density(30) 
  rrt(335) = rrt(335) * density(28) * density(42) 
  rrt(336) = rrt(336) * density(28) * density(33) 
  rrt(337) = rrt(337) * density(28) * density(32) 
  rrt(338) = rrt(338) * density(06) * density(45) 
  rrt(339) = rrt(339) * density(06) * density(45) 
  rrt(340) = rrt(340) * density(06) * density(45) 
  rrt(341) = rrt(341) * density(45) * density(47) 
  rrt(342) = rrt(342) * density(45) * density(47) 
  rrt(343) = rrt(343) * density(43) * density(45) 
  rrt(344) = rrt(344) * density(43) * density(45) 
  rrt(345) = rrt(345) * density(14) * density(45) 
  rrt(346) = rrt(346) * density(14) * density(45) 
  rrt(347) = rrt(347) * density(14) * density(45) 
  rrt(348) = rrt(348) * density(14) * density(45) 
  rrt(349) = rrt(349) * density(14) * density(45) 
  rrt(350) = rrt(350) * density(30) * density(45) 
  rrt(351) = rrt(351) * density(30) * density(45) 
  rrt(352) = rrt(352) * density(30) * density(45) 
  rrt(353) = rrt(353) * density(32) * density(45) 
  rrt(354) = rrt(354) * density(32) * density(45) 
  rrt(355) = rrt(355) * density(34) * density(45) 
  rrt(356) = rrt(356) * density(34) * density(45) 
  rrt(357) = rrt(357) * density(06) * density(12) 
  rrt(358) = rrt(358) * density(06) * density(12) 
  rrt(359) = rrt(359) * density(11) * density(12) 
  rrt(360) = rrt(360) * density(12) * density(47) 
  rrt(361) = rrt(361) * density(12) * density(47) 
  rrt(362) = rrt(362) * density(12) * density(43) 
  rrt(363) = rrt(363) * density(12) * density(14) 
  rrt(364) = rrt(364) * density(12) * density(14) 
  rrt(365) = rrt(365) * density(12) * density(14) 
  rrt(366) = rrt(366) * density(08) * density(12) 
  rrt(367) = rrt(367) * density(08) * density(12) 
  rrt(368) = rrt(368) * density(12) * density(30) 
  rrt(369) = rrt(369) * density(12) * density(30) 
  rrt(370) = rrt(370) * density(12) * density(30) 
  rrt(371) = rrt(371) * density(12) * density(42) 
  rrt(372) = rrt(372) * density(12) * density(42) 
  rrt(373) = rrt(373) * density(12) * density(32) 
  rrt(374) = rrt(374) * density(06) * density(47) 
  rrt(375) = rrt(375) * density(06) * density(43) 
  rrt(376) = rrt(376) * density(06) * density(08) 
  rrt(377) = rrt(377) * density(06) * density(42) 
  rrt(378) = rrt(378) * density(06) * density(33) 
  rrt(379) = rrt(379) * density(06) * density(44) 
  rrt(380) = rrt(380) * density(06) * density(26) 
  rrt(381) = rrt(381) * density(06) * density(34) 
  rrt(382) = rrt(382) * density(06) * density(11) 
  rrt(383) = rrt(383) * density(06) * density(31) 
  rrt(384) = rrt(384) * density(06) * density(47) 
  rrt(385) = rrt(385) * density(06) 
  rrt(386) = rrt(386) * density(11) 
  rrt(387) = rrt(387) * density(11) 
  rrt(388) = rrt(388) * density(08) * density(11) 
  rrt(389) = rrt(389) * density(47)**2 
  rrt(390) = rrt(390) * density(08) * density(47) 
  rrt(391) = rrt(391) * density(42) * density(47) 
  rrt(392) = rrt(392) * density(33) * density(47) 
  rrt(393) = rrt(393) * density(07) * density(47) 
  rrt(394) = rrt(394) * density(44) * density(47) 
  rrt(395) = rrt(395) * density(44) * density(47) 
  rrt(396) = rrt(396) * density(21) * density(47) 
  rrt(397) = rrt(397) * density(23) * density(47) 
  rrt(398) = rrt(398) * density(34) * density(47) 
  rrt(399) = rrt(399) * density(47) 
  rrt(400) = rrt(400) * density(47)**2 
  rrt(401) = rrt(401) * density(34) * density(47) 
  rrt(402) = rrt(402) * density(14) * density(43) 
  rrt(403) = rrt(403) * density(14) * density(43) 
  rrt(404) = rrt(404) * density(23) * density(43) 
  rrt(405) = rrt(405) * density(11) * density(43) 
  rrt(406) = rrt(406) * density(43) * density(47) 
  rrt(407) = rrt(407) * density(23) * density(43) 
  rrt(408) = rrt(408) * density(42) * density(43) 
  rrt(409) = rrt(409) * density(14) * density(42) 
  rrt(410) = rrt(410) * density(14) * density(44) 
  rrt(411) = rrt(411) * density(14) * density(26) 
  rrt(412) = rrt(412) * density(14) * density(34) 
  rrt(413) = rrt(413) * density(14) * density(34) 
  rrt(414) = rrt(414) * density(14) 
  rrt(415) = rrt(415) * density(14) * density(31) 
  rrt(416) = rrt(416) * density(14) * density(43) 
  rrt(417) = rrt(417) * density(14) * density(47) 
  rrt(418) = rrt(418) * density(08)**2 
  rrt(419) = rrt(419) * density(08) * density(30) 
  rrt(420) = rrt(420) * density(08) * density(32) 
  rrt(421) = rrt(421) * density(08) * density(33) 
  rrt(422) = rrt(422) * density(07) * density(08) 
  rrt(423) = rrt(423) * density(08) * density(44) 
  rrt(424) = rrt(424) * density(08) * density(44) 
  rrt(425) = rrt(425) * density(08) * density(21) 
  rrt(426) = rrt(426) * density(08) * density(23) 
  rrt(427) = rrt(427) * density(08) * density(34) 
  rrt(428) = rrt(428) * density(08) * density(34) 
  rrt(429) = rrt(429) * density(08) * density(34) 
  rrt(430) = rrt(430) * density(08) 
  rrt(431) = rrt(431) * density(08)**2 
  rrt(432) = rrt(432) * density(08) * density(31) 
  rrt(433) = rrt(433) * density(08) * density(42) 
  rrt(434) = rrt(434) * density(30) * density(34) 
  rrt(435) = rrt(435) * density(30) * density(34) 
  rrt(436) = rrt(436) * density(23) * density(30) 
  rrt(437) = rrt(437) * density(30) 
  rrt(438) = rrt(438) * density(21) * density(30) 
  rrt(439) = rrt(439) * density(30) * density(32) 
  rrt(440) = rrt(440) * density(21) * density(30) 
  rrt(441) = rrt(441) * density(30)**2 
  rrt(442) = rrt(442) * density(11) * density(30) 
  rrt(443) = rrt(443) * density(30) 
  rrt(444) = rrt(444) * density(08) * density(30) 
  rrt(445) = rrt(445) * density(23) * density(30) 
  rrt(446) = rrt(446) * density(30) * density(47) 
  rrt(447) = rrt(447) * density(30) * density(31) 
  rrt(448) = rrt(448) * density(42)**2 
  rrt(449) = rrt(449) * density(07) * density(42) 
  rrt(450) = rrt(450) * density(42) * density(44) 
  rrt(451) = rrt(451) * density(42) * density(44) 
  rrt(452) = rrt(452) * density(21) * density(42) 
  rrt(453) = rrt(453) * density(26) * density(42) 
  rrt(454) = rrt(454) * density(23) * density(42) 
  rrt(455) = rrt(455) * density(34) * density(42) 
  rrt(456) = rrt(456) * density(34) * density(42) 
  rrt(457) = rrt(457) * density(42) 
  rrt(458) = rrt(458) * density(31) * density(42) 
  rrt(459) = rrt(459) * density(32) * density(34) 
  rrt(460) = rrt(460) * density(23) * density(32) 
  rrt(461) = rrt(461) * density(23) * density(32) 
  rrt(462) = rrt(462) * density(11) * density(32) 
  rrt(463) = rrt(463) * density(31) * density(32) 
  rrt(464) = rrt(464) * density(07) * density(26) 
  rrt(465) = rrt(465) * density(07) * density(34) 
  rrt(466) = rrt(466) * density(07) 
  rrt(467) = rrt(467) * density(07) * density(31) 
  rrt(468) = rrt(468) * density(07) * density(47) 
  rrt(469) = rrt(469) * density(44)**2 
  rrt(470) = rrt(470) * density(21) * density(44) 
  rrt(471) = rrt(471) * density(26) * density(44) 
  rrt(472) = rrt(472) * density(23) * density(44) 
  rrt(473) = rrt(473) * density(34) * density(44) 
  rrt(474) = rrt(474) * density(34) * density(44) 
  rrt(475) = rrt(475) * density(34) * density(44) 
  rrt(476) = rrt(476) * density(44) 
  rrt(477) = rrt(477) * density(44) 
  rrt(478) = rrt(478) * density(31) * density(44) 
  rrt(479) = rrt(479) * density(21) * density(32) 
  rrt(480) = rrt(480) * density(21)**2 
  rrt(481) = rrt(481) * density(21) 
  rrt(482) = rrt(482) * density(21) * density(34) 
  rrt(483) = rrt(483) * density(21) * density(34) 
  rrt(484) = rrt(484) * density(21) 
  rrt(485) = rrt(485) * density(11) * density(21) 
  rrt(486) = rrt(486) * density(21) * density(31) 
  rrt(487) = rrt(487) * density(23) * density(26) 
  rrt(488) = rrt(488) * density(26) * density(34) 
  rrt(489) = rrt(489) * density(26) 
  rrt(490) = rrt(490) * density(31) 
  rrt(491) = rrt(491) * density(31) * density(47) 
  rrt(492) = rrt(492) * density(31) 
  rrt(493) = rrt(493) * density(23) * density(31) 
  rrt(494) = rrt(494) * density(11) * density(25) 
  rrt(495) = rrt(495) * density(25) 
  rrt(496) = rrt(496) * density(25) 
  rrt(497) = rrt(497) * density(25) * density(34) 
  rrt(498) = rrt(498) * density(25) * density(47) 
  rrt(499) = rrt(499) * density(25) * density(42) 
  rrt(500) = rrt(500) * density(25) * density(44) 
  rrt(501) = rrt(501) * density(25) * density(33) 
  rrt(502) = rrt(502) * density(08) * density(25) 
  rrt(503) = rrt(503) * density(25) * density(26) 
  rrt(504) = rrt(504) * density(25) * density(47) 
  rrt(505) = rrt(505) * density(22) 
  rrt(506) = rrt(506) * density(23) 
  rrt(507) = rrt(507) * density(34)**2 
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
  ydot(02) = +rrt(011)-rrt(179)-rrt(181)-rrt(183)-rrt(185)-rrt(187)-rrt(189)-rrt(221)-rrt(223)-rrt(225)-rrt(227)-rrt(229)-rrt(231)&
             -rrt(233)-rrt(235)-rrt(237) 
  ydot(03) = +rrt(003)-rrt(153)-rrt(155)-rrt(157)-rrt(159)-rrt(161)-rrt(163)-rrt(190)-rrt(192)-rrt(194)-rrt(196)-rrt(198)-rrt(200)&
             -rrt(202) 
  ydot(04) = +rrt(072)+rrt(078)+rrt(084)+rrt(088)+rrt(110)+rrt(120)+rrt(129)+rrt(135)+rrt(196)+rrt(197)+rrt(206)+rrt(207)-rrt(247)&
             -rrt(273)+rrt(282)+rrt(288)+rrt(294)+rrt(297)+rrt(298)+rrt(302)+rrt(306)+rrt(314)+rrt(315)-rrt(316)-rrt(317)-rrt(318)&
             -rrt(319)+rrt(320)+rrt(324)+rrt(325)+rrt(334)+rrt(337)+rrt(348)+rrt(351)+rrt(353)+rrt(365)+rrt(367)+rrt(369)+rrt(371) 
  ydot(05) = +rrt(005)-rrt(165)-rrt(167)-rrt(169)-rrt(171)-rrt(173)-rrt(204)-rrt(206)-rrt(208)-rrt(210)-rrt(212) 
  ydot(06) = -rrt(001)-rrt(002)-rrt(013)-rrt(014)-rrt(015)-rrt(019)-rrt(020)-rrt(021)-rrt(022)+rrt(033)+rrt(039)+rrt(054)+rrt(059)&
             +rrt(064)+rrt(075)+rrt(082)+rrt(101)+rrt(110)+rrt(121)+rrt(161)+rrt(162)+rrt(188)+rrt(189)+rrt(202)+rrt(203)+rrt(232)&
             +rrt(233)+rrt(240)+rrt(278)+rrt(279)+rrt(280)+rrt(281)+rrt(282)-rrt(284)+rrt(285)+rrt(287)+rrt(289)-rrt(292)-rrt(293)&
             +rrt(296)+rrt(297)-rrt(299)-rrt(300)-rrt(301)-rrt(302)-rrt(303)-rrt(305)-rrt(306)-rrt(307)-rrt(320)-rrt(326)-rrt(327)&
             -rrt(338)-rrt(339)-rrt(340)-rrt(357)-rrt(358)-rrt(374)-rrt(375)-rrt(376)-rrt(377)-rrt(378)-rrt(379)-rrt(380)-rrt(381)&
             -rrt(382)-rrt(383)-rrt(384)-rrt(385)+rrt(413)+rrt(494) 
  ydot(07) = -rrt(010)-rrt(011)-rrt(049)-rrt(050)-rrt(051)-rrt(052)-rrt(053)-rrt(054)-rrt(095)-rrt(096)-rrt(097)-rrt(098)-rrt(099)&
             -rrt(100)-rrt(101)-rrt(102)-rrt(103)+rrt(249)+rrt(379)-rrt(393)+rrt(410)-rrt(422)+rrt(423)-rrt(449)+rrt(450)-rrt(464)&
             -rrt(465)-rrt(466)-rrt(467)-rrt(468)+rrt(469)+rrt(470)+rrt(472)+rrt(474)+rrt(500) 
  ydot(08) = +rrt(029)-rrt(035)-rrt(036)-rrt(037)-rrt(038)-rrt(039)-rrt(040)+rrt(053)-rrt(076)-rrt(077)-rrt(078)-rrt(079)-rrt(080)&
             -rrt(081)-rrt(082)+rrt(102)+rrt(113)+rrt(125)+rrt(153)+rrt(154)+rrt(186)+rrt(187)+rrt(234)+rrt(235)+rrt(245)+rrt(264)&
             -rrt(332)-rrt(366)-rrt(367)-rrt(376)-rrt(388)-rrt(390)+rrt(394)+rrt(409)+rrt(410)+rrt(411)+rrt(412)+rrt(415)+rrt(417)&
             -  2.d0 * rrt(418)-rrt(419)-rrt(420)-rrt(421)-rrt(422)-rrt(423)-rrt(424)-rrt(425)-rrt(426)-rrt(427)-rrt(428)-rrt(429)&
             -rrt(430)-  2.d0 * rrt(431)-rrt(432)-rrt(433)+rrt(435)+rrt(436)+rrt(438)+rrt(441)-rrt(444)+rrt(466)+rrt(475)+rrt(490)&
             +  2.d0 * rrt(496)-rrt(502) 
  ydot(09) = +rrt(073)+rrt(079)+rrt(089)+rrt(093)+rrt(121)+rrt(130)+rrt(136)+rrt(198)+rrt(199)+rrt(214)+rrt(215)+rrt(216)-rrt(248)&
             -rrt(274)+rrt(289)+rrt(295)+rrt(303)+rrt(307)+rrt(318)+rrt(319)-rrt(320)-rrt(321)-rrt(322)-rrt(323)-rrt(324)-rrt(325)&
             +rrt(326)+rrt(336)+rrt(349)+rrt(352)+rrt(354)+rrt(370)+rrt(372)+rrt(373) 
  ydot(10) = +rrt(012) 
  ydot(11) = +rrt(013)-rrt(016)-rrt(017)-rrt(023)-rrt(024)-rrt(025)+  2.d0 * rrt(034)+rrt(040)+rrt(044)+rrt(053)+rrt(058)+rrt(063)&
             +rrt(066)+rrt(074)+rrt(081)+rrt(087)+rrt(100)+rrt(109)+rrt(120)+rrt(130)+rrt(139)+rrt(140)+  2.d0 * rrt(163)&
             +  2.d0 * rrt(164)+rrt(171)+rrt(172)+rrt(186)+rrt(187)+rrt(200)+rrt(201)+rrt(212)+rrt(213)+rrt(230)+rrt(231)+rrt(241)&
             +rrt(256)+rrt(258)+rrt(270)+rrt(284)+rrt(286)+rrt(288)+rrt(292)+rrt(298)+rrt(299)+rrt(320)+rrt(326)-rrt(328)-rrt(359)&
             +  2.d0 * rrt(374)+rrt(376)+rrt(377)+rrt(378)+rrt(379)+rrt(380)+rrt(381)-rrt(382)+rrt(383)+rrt(385)-rrt(386)-rrt(387)&
             -rrt(388)+rrt(390)+rrt(391)+rrt(393)+rrt(395)+rrt(396)+rrt(397)+rrt(401)-rrt(405)+rrt(407)+rrt(413)+  2.d0 * rrt(414)&
             +rrt(416)+rrt(417)+  2.d0 * rrt(427)-rrt(442)-rrt(462)+rrt(466)+rrt(475)+rrt(477)+rrt(484)-rrt(485)+rrt(489)+rrt(492)&
             -rrt(494)+rrt(495)+rrt(498)+rrt(505) 
  ydot(12) = +rrt(092)-rrt(254)+rrt(356)-rrt(357)-rrt(358)-rrt(359)-rrt(360)-rrt(361)-rrt(362)-rrt(363)-rrt(364)-rrt(365)-rrt(366)&
             -rrt(367)-rrt(368)-rrt(369)-rrt(370)-rrt(371)-rrt(372)-rrt(373) 
  ydot(13) = +rrt(022)+rrt(025)+rrt(027)+rrt(028)+rrt(082)+rrt(087)+rrt(091)+rrt(094)+rrt(114)+rrt(125)+rrt(133)+rrt(138)+rrt(151)&
             +rrt(152)+rrt(212)+rrt(213)+rrt(217)+rrt(218)+rrt(219)-rrt(243)-rrt(305)-rrt(306)-rrt(307)-rrt(308)+rrt(344)+rrt(361)&
             +rrt(362) 
  ydot(14) = -rrt(003)-rrt(004)-rrt(029)-rrt(030)-rrt(031)-rrt(032)-rrt(033)-rrt(034)+rrt(052)-rrt(069)-rrt(070)-rrt(071)-rrt(072)&
             -rrt(073)-rrt(074)-rrt(075)+rrt(103)+rrt(114)+rrt(184)+rrt(185)+rrt(236)+rrt(237)+rrt(244)-rrt(280)-rrt(285)-rrt(296)&
             +rrt(309)-rrt(316)-rrt(321)-rrt(322)-rrt(331)-rrt(345)-rrt(346)-rrt(347)-rrt(348)-rrt(349)-rrt(363)-rrt(364)-rrt(365)&
             +rrt(376)+rrt(382)+rrt(384)+rrt(388)-rrt(402)-rrt(403)-rrt(409)-rrt(410)-rrt(411)-rrt(412)-rrt(413)-rrt(414)-rrt(415)&
             -rrt(416)-rrt(417)+rrt(418)+rrt(419)+rrt(420)+rrt(422)+rrt(424)+rrt(425)+rrt(426)+rrt(429)+rrt(445)+rrt(502) 
  ydot(15) = +rrt(006)-rrt(166)-rrt(168)-rrt(170)-rrt(172)-rrt(174)-rrt(205)-rrt(207)-rrt(209)-rrt(211)-rrt(213) 
  ydot(16) = +rrt(071)+rrt(077)+rrt(083)+rrt(101)+rrt(109)+rrt(119)+rrt(128)+rrt(194)+rrt(195)+rrt(204)+rrt(205)+rrt(232)+rrt(233)&
             -rrt(246)-rrt(271)-rrt(272)+rrt(285)+rrt(287)+rrt(301)+rrt(305)+rrt(309)+rrt(312)-rrt(313)-rrt(314)-rrt(315)+rrt(322)&
             +rrt(323)+rrt(335)+rrt(347)+rrt(350)+rrt(364)+rrt(366)+rrt(368) 
  ydot(17) = +rrt(051)+rrt(057)+rrt(061)+rrt(065)-rrt(067)-rrt(068)-rrt(134)-rrt(135)-rrt(136)-rrt(137)-rrt(138)+rrt(182)+rrt(183)&
             +rrt(253) 
  ydot(18) = +rrt(096)+rrt(104)+rrt(222)+rrt(223)-rrt(250) 
  ydot(19) = +rrt(001)-rrt(139)-rrt(141)-rrt(143)-rrt(145)-rrt(147)-rrt(149)-rrt(151) 
  ydot(20) = +rrt(099)+rrt(107)+rrt(117)+rrt(127)+rrt(134)+rrt(228)+rrt(229)-rrt(253) 
  ydot(21) = -rrt(012)+rrt(050)+rrt(055)-rrt(060)-rrt(061)-rrt(062)-rrt(063)-rrt(064)-rrt(115)-rrt(116)-rrt(117)-rrt(118)-rrt(119)&
             -rrt(120)-rrt(121)-rrt(122)-rrt(123)-rrt(124)-rrt(125)+rrt(180)+rrt(181)+rrt(251)+rrt(380)+rrt(395)-rrt(396)+rrt(402)&
             +rrt(411)+rrt(424)-rrt(425)-rrt(438)-rrt(440)+rrt(446)+rrt(447)+rrt(451)-rrt(452)+rrt(453)+rrt(463)+rrt(464)+rrt(469)&
             -rrt(470)+  2.d0 * rrt(471)+rrt(473)+rrt(476)+rrt(478)-rrt(479)-  2.d0 * rrt(480)-rrt(481)-rrt(482)-rrt(483)-rrt(484)&
             -rrt(485)-rrt(486)+rrt(487)+rrt(488)+rrt(492)+rrt(503) 
  ydot(22) = +rrt(504)-rrt(505) 
  ydot(23) = +rrt(014)+rrt(015)+rrt(017)+rrt(021)+rrt(022)+rrt(025)+rrt(030)+rrt(031)+  2.d0 * rrt(032)+rrt(036)+rrt(038)+rrt(042)&
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
  ydot(24) = +rrt(010)-rrt(178)-rrt(180)-rrt(182)-rrt(184)-rrt(186)-rrt(188)-rrt(220)-rrt(222)-rrt(224)-rrt(226)-rrt(228)-rrt(230)&
             -rrt(232)-rrt(234)-rrt(236) 
  ydot(25) = +rrt(383)+rrt(415)+rrt(431)+rrt(432)+rrt(458)+rrt(467)+rrt(468)+rrt(478)+rrt(486)+rrt(493)-rrt(494)-rrt(495)-rrt(496)&
             -rrt(497)-rrt(498)-rrt(499)-rrt(500)-rrt(501)-rrt(502)-rrt(503)-rrt(504) 
  ydot(26) = +rrt(056)+rrt(060)-rrt(065)-rrt(066)-rrt(126)-rrt(127)-rrt(128)-rrt(129)-rrt(130)-rrt(131)-rrt(132)-rrt(133)+rrt(252)&
             -rrt(380)+rrt(396)-rrt(411)+rrt(425)+rrt(438)+rrt(452)-rrt(453)+rrt(462)+rrt(463)-rrt(464)+rrt(470)-rrt(471)+rrt(479)&
             +rrt(480)+rrt(481)+rrt(482)+rrt(486)-rrt(487)-rrt(488)-rrt(489)-rrt(503) 
  ydot(27) = +rrt(070)+rrt(076)+rrt(100)+rrt(108)+rrt(118)+rrt(192)+rrt(193)+rrt(230)+rrt(231)-rrt(245)-rrt(266)-rrt(267)-rrt(268)&
             -rrt(269)-rrt(270)+rrt(280)+rrt(281)+rrt(286)+rrt(293)+rrt(296)+rrt(300)+rrt(310)+rrt(311)-rrt(312)+rrt(313)+rrt(316)&
             +rrt(317)+rrt(321)+rrt(331)+rrt(333)+rrt(346)+rrt(363) 
  ydot(28) = -rrt(276)+rrt(277)-rrt(327)-rrt(328)-rrt(329)-rrt(330)-rrt(331)-rrt(332)-rrt(333)-rrt(334)-rrt(335)-rrt(336)-rrt(337)&
             +rrt(355) 
  ydot(29) = +rrt(004)-rrt(154)-rrt(156)-rrt(158)-rrt(160)-rrt(162)-rrt(164)-rrt(191)-rrt(193)-rrt(195)-rrt(197)-rrt(199)-rrt(201)&
             -rrt(203) 
  ydot(30) = -rrt(005)-rrt(006)+rrt(030)+rrt(035)-rrt(041)-rrt(042)-rrt(043)-rrt(044)-rrt(045)+rrt(054)+rrt(058)+rrt(062)-rrt(083)&
             -rrt(084)-rrt(085)-rrt(086)-rrt(087)+rrt(112)+rrt(124)+rrt(133)+rrt(155)+rrt(156)+rrt(188)+rrt(189)+rrt(246)+rrt(265)&
             +rrt(266)-rrt(281)-rrt(286)-rrt(287)-rrt(297)-rrt(309)+rrt(314)+rrt(316)-rrt(317)+rrt(322)-rrt(323)-rrt(333)-rrt(334)&
             -rrt(350)-rrt(351)-rrt(352)-rrt(368)-rrt(369)-rrt(370)+rrt(375)+rrt(377)+rrt(390)+rrt(394)+rrt(409)+rrt(416)+rrt(418)&
             -rrt(419)+rrt(421)+rrt(423)+rrt(428)+rrt(430)+rrt(432)+  2.d0 * rrt(433)-rrt(434)-rrt(435)-rrt(436)-rrt(437)-rrt(438)&
             -rrt(439)-rrt(440)-  2.d0 * rrt(441)-rrt(442)-rrt(443)-rrt(444)-rrt(445)-rrt(446)-rrt(447)+rrt(448)+rrt(449)+rrt(451)&
             +rrt(452)+rrt(454)+rrt(456)+rrt(460)+rrt(477)+rrt(490)+rrt(491)+rrt(499) 
  ydot(31) = -rrt(383)-rrt(415)-rrt(432)+rrt(444)-rrt(447)-rrt(458)-rrt(463)-rrt(467)-rrt(478)+rrt(485)-rrt(486)-rrt(490)-rrt(491)&
             -rrt(492)-rrt(493)+rrt(494)+rrt(497)+rrt(498)+rrt(499)+rrt(500)+rrt(501)+rrt(502)+rrt(503)+rrt(505) 
  ydot(32) = -rrt(007)-rrt(008)-rrt(009)+rrt(032)+rrt(038)+rrt(042)+rrt(043)+rrt(046)-rrt(048)+rrt(064)+rrt(066)+rrt(068)+rrt(092)&
             -rrt(093)-rrt(094)+rrt(122)+rrt(131)+rrt(137)+rrt(159)+rrt(160)+rrt(167)+rrt(168)+rrt(169)+rrt(170)+rrt(248)+rrt(268)&
             +rrt(269)+rrt(272)+rrt(273)-rrt(282)-rrt(288)-rrt(289)-rrt(310)+rrt(313)+rrt(317)+rrt(318)+rrt(323)+rrt(324)-rrt(337)&
             -rrt(353)-rrt(354)-rrt(373)+rrt(378)+rrt(389)+rrt(391)+rrt(392)+rrt(400)+rrt(406)+rrt(408)-rrt(420)+rrt(421)-rrt(439)&
             +rrt(443)+rrt(448)+rrt(450)+rrt(453)+rrt(455)+rrt(457)+rrt(458)-rrt(459)-rrt(460)-rrt(461)-rrt(462)-rrt(463)-rrt(479)&
             +rrt(489)+rrt(501) 
  ydot(33) = -rrt(318)-rrt(336)-rrt(378)-rrt(392)+rrt(420)-rrt(421)-rrt(501) 
  ydot(34) = +rrt(013)+rrt(015)+rrt(016)+rrt(018)+rrt(020)+rrt(022)+rrt(024)+rrt(027)+rrt(029)+rrt(031)+rrt(035)+  2.d0 * rrt(037)&
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
  ydot(35) = -rrt(256)-rrt(257)-rrt(278)-rrt(279)-rrt(280)-rrt(281)-rrt(282)-rrt(283)+rrt(284)+rrt(290)+rrt(327)+rrt(338) 
  ydot(36) = -rrt(326) 
  ydot(37) = +rrt(008)-rrt(176)-rrt(215)-rrt(218) 
  ydot(38) = +rrt(097)+rrt(105)+rrt(115)+rrt(224)+rrt(225)-rrt(251) 
  ydot(39) = +rrt(098)+rrt(106)+rrt(116)+rrt(126)+rrt(226)+rrt(227)-rrt(252) 
  ydot(40) = +rrt(009)-rrt(177)-rrt(216)-rrt(219) 
  ydot(41) = +rrt(095)+rrt(220)+rrt(221)-rrt(249) 
  ydot(42) = +rrt(031)+rrt(036)+rrt(037)+rrt(041)-rrt(046)-rrt(047)+rrt(059)+rrt(063)+rrt(067)-rrt(088)-rrt(089)-rrt(090)-rrt(091)&
             -rrt(092)+rrt(111)+rrt(123)+rrt(132)+rrt(138)+rrt(157)+rrt(158)+rrt(165)+rrt(166)+rrt(247)+rrt(267)+rrt(271)-rrt(298)&
             +rrt(310)-rrt(313)-rrt(314)+rrt(321)-rrt(324)-rrt(335)-rrt(371)-rrt(372)-rrt(377)-rrt(391)+rrt(405)-rrt(408)-rrt(409)&
             +rrt(419)-rrt(433)+rrt(434)+rrt(437)+  2.d0 * rrt(439)+rrt(440)+rrt(441)-  2.d0 * rrt(448)-rrt(449)-rrt(450)-rrt(451)&
             -rrt(452)-rrt(453)-rrt(454)-rrt(455)-rrt(456)-rrt(457)-rrt(458)+rrt(459)+rrt(461)+rrt(479)+rrt(484)-rrt(499) 
  ydot(43) = +rrt(015)+rrt(017)+rrt(018)-rrt(028)+rrt(039)+rrt(044)+rrt(047)+  2.d0 * rrt(048)+rrt(067)+rrt(085)+rrt(090)+rrt(094)&
             +rrt(118)+rrt(128)+rrt(135)+rrt(143)+rrt(144)+rrt(171)+rrt(172)+  2.d0 * rrt(175)+  2.d0 * rrt(176)+  2.d0 * rrt(177)&
             +rrt(208)+rrt(209)+rrt(217)+rrt(218)+rrt(219)+rrt(243)+rrt(260)+rrt(262)+rrt(263)+  2.d0 * rrt(274)-rrt(279)-rrt(295)&
             -rrt(330)-rrt(343)-rrt(344)-rrt(362)-rrt(375)+rrt(387)+rrt(392)+rrt(398)+rrt(399)-rrt(402)-rrt(403)-rrt(404)-rrt(405)&
             -rrt(406)-rrt(407)-rrt(408)-rrt(416) 
  ydot(44) = +rrt(049)-rrt(055)-rrt(056)-rrt(057)-rrt(058)-rrt(059)-rrt(104)-rrt(105)-rrt(106)-rrt(107)-rrt(108)-rrt(109)-rrt(110)&
             -rrt(111)-rrt(112)-rrt(113)-rrt(114)+rrt(178)+rrt(179)+rrt(250)-rrt(379)+rrt(393)-rrt(394)-rrt(395)+rrt(403)-rrt(410)&
             +rrt(422)-rrt(423)-rrt(424)+rrt(440)+rrt(442)+rrt(447)+rrt(449)-rrt(450)-rrt(451)+rrt(464)+rrt(465)+rrt(467)&
             -  2.d0 * rrt(469)-rrt(470)-rrt(471)-rrt(472)-rrt(473)-rrt(474)-rrt(475)-rrt(476)-rrt(477)-rrt(478)+rrt(480)+rrt(483)&
             +rrt(491)+rrt(495)-rrt(500) 
  ydot(45) = +rrt(239)-rrt(255)-rrt(275)-rrt(277)-rrt(338)-rrt(339)-rrt(340)-rrt(341)-rrt(342)-rrt(343)-rrt(344)-rrt(345)-rrt(346)&
             -rrt(347)-rrt(348)-rrt(349)-rrt(350)-rrt(351)-rrt(352)-rrt(353)-rrt(354)-rrt(355)-rrt(356) 
  ydot(46) = +rrt(021)+rrt(024)+rrt(026)+rrt(075)+rrt(081)+rrt(086)+rrt(090)+rrt(103)+rrt(113)+rrt(124)+rrt(132)+rrt(137)+rrt(149)&
             +rrt(150)+rrt(202)+rrt(203)+rrt(210)+rrt(211)+rrt(236)+rrt(237)-rrt(242)-rrt(263)+rrt(279)-rrt(299)-rrt(300)-rrt(301)&
             -rrt(302)-rrt(303)-rrt(304)+rrt(308)+rrt(330)+rrt(342)+rrt(343)+rrt(360) 
  ydot(47) = +rrt(014)+rrt(016)-rrt(018)-rrt(026)-rrt(027)+rrt(033)+rrt(040)+  2.d0 * rrt(045)+rrt(047)+rrt(052)+rrt(062)+rrt(068)&
             +rrt(080)+rrt(086)+rrt(091)+rrt(108)+rrt(119)+rrt(129)+rrt(136)+rrt(141)+rrt(142)+rrt(161)+rrt(162)+  2.d0 * rrt(173)&
             +  2.d0 * rrt(174)+rrt(184)+rrt(185)+rrt(210)+rrt(211)+rrt(242)+rrt(257)+rrt(259)+rrt(261)+rrt(270)-rrt(278)-rrt(294)&
             -rrt(329)-rrt(341)-rrt(342)-rrt(360)-rrt(361)-rrt(374)-rrt(384)+rrt(386)+rrt(388)-  2.d0 * rrt(389)-rrt(390)-rrt(391)&
             -rrt(392)-rrt(393)-rrt(394)-rrt(395)-rrt(396)-rrt(397)-rrt(398)-rrt(399)-  2.d0 * rrt(400)-rrt(401)+rrt(404)-rrt(406)&
             +rrt(408)-rrt(417)-rrt(446)-rrt(468)-rrt(491)-rrt(498)-rrt(504) 
  ydot(48) = +rrt(019)+rrt(111)+rrt(122)+rrt(145)+rrt(146)-rrt(240)-rrt(258)-rrt(259)-rrt(260)+rrt(283)-rrt(284)-rrt(285)-rrt(286)&
             -rrt(287)-rrt(288)-rrt(289)-rrt(290)-rrt(291)+rrt(292)+rrt(328)+rrt(339)+rrt(357) 
  ydot(49) = +rrt(020)+rrt(023)+rrt(074)+rrt(080)+rrt(085)+rrt(102)+rrt(112)+rrt(123)+rrt(131)+rrt(147)+rrt(148)+rrt(200)+rrt(201)&
             +rrt(208)+rrt(209)+rrt(234)+rrt(235)-rrt(241)-rrt(261)-rrt(262)+rrt(278)+rrt(291)-rrt(292)-rrt(293)-rrt(294)-rrt(295)&
             -rrt(296)-rrt(297)-rrt(298)+rrt(299)+rrt(304)+rrt(329)+rrt(340)+rrt(341)+rrt(358)+rrt(359) 
  ydot(50) = +rrt(007)-rrt(175)-rrt(214)-rrt(217) 
  ydot(51) = +rrt(002)-rrt(140)-rrt(142)-rrt(144)-rrt(146)-rrt(148)-rrt(150)-rrt(152) 
  ydot(52) = +rrt(069)+rrt(190)+rrt(191)-rrt(244)-rrt(264)-rrt(265)-rrt(309)-rrt(310)-rrt(311)+rrt(332)+rrt(345) 
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
  pd(06,01) = pd(06,01) - rrt(001) * density(06) 
  pd(06,06) = pd(06,06) - rrt(001) * density(01) 
  pd(19,01) = pd(19,01) + rrt(001) * density(06) 
  pd(19,06) = pd(19,06) + rrt(001) * density(01) 
  pd(06,01) = pd(06,01) - rrt(002) * density(06) 
  pd(06,06) = pd(06,06) - rrt(002) * density(01) 
  pd(51,01) = pd(51,01) + rrt(002) * density(06) 
  pd(51,06) = pd(51,06) + rrt(002) * density(01) 
  pd(03,01) = pd(03,01) + rrt(003) * density(14) 
  pd(03,14) = pd(03,14) + rrt(003) * density(01) 
  pd(14,01) = pd(14,01) - rrt(003) * density(14) 
  pd(14,14) = pd(14,14) - rrt(003) * density(01) 
  pd(14,01) = pd(14,01) - rrt(004) * density(14) 
  pd(14,14) = pd(14,14) - rrt(004) * density(01) 
  pd(29,01) = pd(29,01) + rrt(004) * density(14) 
  pd(29,14) = pd(29,14) + rrt(004) * density(01) 
  pd(05,01) = pd(05,01) + rrt(005) * density(30) 
  pd(05,30) = pd(05,30) + rrt(005) * density(01) 
  pd(30,01) = pd(30,01) - rrt(005) * density(30) 
  pd(30,30) = pd(30,30) - rrt(005) * density(01) 
  pd(15,01) = pd(15,01) + rrt(006) * density(30) 
  pd(15,30) = pd(15,30) + rrt(006) * density(01) 
  pd(30,01) = pd(30,01) - rrt(006) * density(30) 
  pd(30,30) = pd(30,30) - rrt(006) * density(01) 
  pd(32,01) = pd(32,01) - rrt(007) * density(32) 
  pd(32,32) = pd(32,32) - rrt(007) * density(01) 
  pd(50,01) = pd(50,01) + rrt(007) * density(32) 
  pd(50,32) = pd(50,32) + rrt(007) * density(01) 
  pd(32,01) = pd(32,01) - rrt(008) * density(32) 
  pd(32,32) = pd(32,32) - rrt(008) * density(01) 
  pd(37,01) = pd(37,01) + rrt(008) * density(32) 
  pd(37,32) = pd(37,32) + rrt(008) * density(01) 
  pd(32,01) = pd(32,01) - rrt(009) * density(32) 
  pd(32,32) = pd(32,32) - rrt(009) * density(01) 
  pd(40,01) = pd(40,01) + rrt(009) * density(32) 
  pd(40,32) = pd(40,32) + rrt(009) * density(01) 
  pd(07,01) = pd(07,01) - rrt(010) * density(07) 
  pd(07,07) = pd(07,07) - rrt(010) * density(01) 
  pd(24,01) = pd(24,01) + rrt(010) * density(07) 
  pd(24,07) = pd(24,07) + rrt(010) * density(01) 
  pd(02,01) = pd(02,01) + rrt(011) * density(07) 
  pd(02,07) = pd(02,07) + rrt(011) * density(01) 
  pd(07,01) = pd(07,01) - rrt(011) * density(07) 
  pd(07,07) = pd(07,07) - rrt(011) * density(01) 
  pd(10,01) = pd(10,01) + rrt(012) * density(21) 
  pd(10,21) = pd(10,21) + rrt(012) * density(01) 
  pd(21,01) = pd(21,01) - rrt(012) * density(21) 
  pd(21,21) = pd(21,21) - rrt(012) * density(01) 
  pd(06,01) = pd(06,01) - rrt(013) * density(06) 
  pd(06,06) = pd(06,06) - rrt(013) * density(01) 
  pd(11,01) = pd(11,01) + rrt(013) * density(06) 
  pd(11,06) = pd(11,06) + rrt(013) * density(01) 
  pd(34,01) = pd(34,01) + rrt(013) * density(06) 
  pd(34,06) = pd(34,06) + rrt(013) * density(01) 
  pd(06,01) = pd(06,01) - rrt(014) * density(06) 
  pd(06,06) = pd(06,06) - rrt(014) * density(01) 
  pd(23,01) = pd(23,01) + rrt(014) * density(06) 
  pd(23,06) = pd(23,06) + rrt(014) * density(01) 
  pd(47,01) = pd(47,01) + rrt(014) * density(06) 
  pd(47,06) = pd(47,06) + rrt(014) * density(01) 
  pd(06,01) = pd(06,01) - rrt(015) * density(06) 
  pd(06,06) = pd(06,06) - rrt(015) * density(01) 
  pd(23,01) = pd(23,01) + rrt(015) * density(06) 
  pd(23,06) = pd(23,06) + rrt(015) * density(01) 
  pd(34,01) = pd(34,01) + rrt(015) * density(06) 
  pd(34,06) = pd(34,06) + rrt(015) * density(01) 
  pd(43,01) = pd(43,01) + rrt(015) * density(06) 
  pd(43,06) = pd(43,06) + rrt(015) * density(01) 
  pd(11,01) = pd(11,01) - rrt(016) * density(11) 
  pd(11,11) = pd(11,11) - rrt(016) * density(01) 
  pd(34,01) = pd(34,01) + rrt(016) * density(11) 
  pd(34,11) = pd(34,11) + rrt(016) * density(01) 
  pd(47,01) = pd(47,01) + rrt(016) * density(11) 
  pd(47,11) = pd(47,11) + rrt(016) * density(01) 
  pd(11,01) = pd(11,01) - rrt(017) * density(11) 
  pd(11,11) = pd(11,11) - rrt(017) * density(01) 
  pd(23,01) = pd(23,01) + rrt(017) * density(11) 
  pd(23,11) = pd(23,11) + rrt(017) * density(01) 
  pd(43,01) = pd(43,01) + rrt(017) * density(11) 
  pd(43,11) = pd(43,11) + rrt(017) * density(01) 
  pd(34,01) = pd(34,01) + rrt(018) * density(47) 
  pd(34,47) = pd(34,47) + rrt(018) * density(01) 
  pd(43,01) = pd(43,01) + rrt(018) * density(47) 
  pd(43,47) = pd(43,47) + rrt(018) * density(01) 
  pd(47,01) = pd(47,01) - rrt(018) * density(47) 
  pd(47,47) = pd(47,47) - rrt(018) * density(01) 
  pd(01,01) = pd(01,01) + rrt(019) * density(06) 
  pd(01,06) = pd(01,06) + rrt(019) * density(01) 
  pd(06,01) = pd(06,01) - rrt(019) * density(06) 
  pd(06,06) = pd(06,06) - rrt(019) * density(01) 
  pd(48,01) = pd(48,01) + rrt(019) * density(06) 
  pd(48,06) = pd(48,06) + rrt(019) * density(01) 
  pd(01,01) = pd(01,01) + rrt(020) * density(06) 
  pd(01,06) = pd(01,06) + rrt(020) * density(01) 
  pd(06,01) = pd(06,01) - rrt(020) * density(06) 
  pd(06,06) = pd(06,06) - rrt(020) * density(01) 
  pd(34,01) = pd(34,01) + rrt(020) * density(06) 
  pd(34,06) = pd(34,06) + rrt(020) * density(01) 
  pd(49,01) = pd(49,01) + rrt(020) * density(06) 
  pd(49,06) = pd(49,06) + rrt(020) * density(01) 
  pd(01,01) = pd(01,01) + rrt(021) * density(06) 
  pd(01,06) = pd(01,06) + rrt(021) * density(01) 
  pd(06,01) = pd(06,01) - rrt(021) * density(06) 
  pd(06,06) = pd(06,06) - rrt(021) * density(01) 
  pd(23,01) = pd(23,01) + rrt(021) * density(06) 
  pd(23,06) = pd(23,06) + rrt(021) * density(01) 
  pd(46,01) = pd(46,01) + rrt(021) * density(06) 
  pd(46,06) = pd(46,06) + rrt(021) * density(01) 
  pd(01,01) = pd(01,01) + rrt(022) * density(06) 
  pd(01,06) = pd(01,06) + rrt(022) * density(01) 
  pd(06,01) = pd(06,01) - rrt(022) * density(06) 
  pd(06,06) = pd(06,06) - rrt(022) * density(01) 
  pd(13,01) = pd(13,01) + rrt(022) * density(06) 
  pd(13,06) = pd(13,06) + rrt(022) * density(01) 
  pd(23,01) = pd(23,01) + rrt(022) * density(06) 
  pd(23,06) = pd(23,06) + rrt(022) * density(01) 
  pd(34,01) = pd(34,01) + rrt(022) * density(06) 
  pd(34,06) = pd(34,06) + rrt(022) * density(01) 
  pd(01,01) = pd(01,01) + rrt(023) * density(11) 
  pd(01,11) = pd(01,11) + rrt(023) * density(01) 
  pd(11,01) = pd(11,01) - rrt(023) * density(11) 
  pd(11,11) = pd(11,11) - rrt(023) * density(01) 
  pd(49,01) = pd(49,01) + rrt(023) * density(11) 
  pd(49,11) = pd(49,11) + rrt(023) * density(01) 
  pd(01,01) = pd(01,01) + rrt(024) * density(11) 
  pd(01,11) = pd(01,11) + rrt(024) * density(01) 
  pd(11,01) = pd(11,01) - rrt(024) * density(11) 
  pd(11,11) = pd(11,11) - rrt(024) * density(01) 
  pd(34,01) = pd(34,01) + rrt(024) * density(11) 
  pd(34,11) = pd(34,11) + rrt(024) * density(01) 
  pd(46,01) = pd(46,01) + rrt(024) * density(11) 
  pd(46,11) = pd(46,11) + rrt(024) * density(01) 
  pd(01,01) = pd(01,01) + rrt(025) * density(11) 
  pd(01,11) = pd(01,11) + rrt(025) * density(01) 
  pd(11,01) = pd(11,01) - rrt(025) * density(11) 
  pd(11,11) = pd(11,11) - rrt(025) * density(01) 
  pd(13,01) = pd(13,01) + rrt(025) * density(11) 
  pd(13,11) = pd(13,11) + rrt(025) * density(01) 
  pd(23,01) = pd(23,01) + rrt(025) * density(11) 
  pd(23,11) = pd(23,11) + rrt(025) * density(01) 
  pd(01,01) = pd(01,01) + rrt(026) * density(47) 
  pd(01,47) = pd(01,47) + rrt(026) * density(01) 
  pd(46,01) = pd(46,01) + rrt(026) * density(47) 
  pd(46,47) = pd(46,47) + rrt(026) * density(01) 
  pd(47,01) = pd(47,01) - rrt(026) * density(47) 
  pd(47,47) = pd(47,47) - rrt(026) * density(01) 
  pd(01,01) = pd(01,01) + rrt(027) * density(47) 
  pd(01,47) = pd(01,47) + rrt(027) * density(01) 
  pd(13,01) = pd(13,01) + rrt(027) * density(47) 
  pd(13,47) = pd(13,47) + rrt(027) * density(01) 
  pd(34,01) = pd(34,01) + rrt(027) * density(47) 
  pd(34,47) = pd(34,47) + rrt(027) * density(01) 
  pd(47,01) = pd(47,01) - rrt(027) * density(47) 
  pd(47,47) = pd(47,47) - rrt(027) * density(01) 
  pd(01,01) = pd(01,01) + rrt(028) * density(43) 
  pd(01,43) = pd(01,43) + rrt(028) * density(01) 
  pd(13,01) = pd(13,01) + rrt(028) * density(43) 
  pd(13,43) = pd(13,43) + rrt(028) * density(01) 
  pd(43,01) = pd(43,01) - rrt(028) * density(43) 
  pd(43,43) = pd(43,43) - rrt(028) * density(01) 
  pd(08,01) = pd(08,01) + rrt(029) * density(14) 
  pd(08,14) = pd(08,14) + rrt(029) * density(01) 
  pd(14,01) = pd(14,01) - rrt(029) * density(14) 
  pd(14,14) = pd(14,14) - rrt(029) * density(01) 
  pd(34,01) = pd(34,01) + rrt(029) * density(14) 
  pd(34,14) = pd(34,14) + rrt(029) * density(01) 
  pd(14,01) = pd(14,01) - rrt(030) * density(14) 
  pd(14,14) = pd(14,14) - rrt(030) * density(01) 
  pd(23,01) = pd(23,01) + rrt(030) * density(14) 
  pd(23,14) = pd(23,14) + rrt(030) * density(01) 
  pd(30,01) = pd(30,01) + rrt(030) * density(14) 
  pd(30,14) = pd(30,14) + rrt(030) * density(01) 
  pd(14,01) = pd(14,01) - rrt(031) * density(14) 
  pd(14,14) = pd(14,14) - rrt(031) * density(01) 
  pd(23,01) = pd(23,01) + rrt(031) * density(14) 
  pd(23,14) = pd(23,14) + rrt(031) * density(01) 
  pd(34,01) = pd(34,01) + rrt(031) * density(14) 
  pd(34,14) = pd(34,14) + rrt(031) * density(01) 
  pd(42,01) = pd(42,01) + rrt(031) * density(14) 
  pd(42,14) = pd(42,14) + rrt(031) * density(01) 
  pd(14,01) = pd(14,01) - rrt(032) * density(14) 
  pd(14,14) = pd(14,14) - rrt(032) * density(01) 
  pd(23,01) = pd(23,01) + rrt(032) * density(14) * 2.0d0
  pd(23,14) = pd(23,14) + rrt(032) * density(01) * 2.0d0
  pd(32,01) = pd(32,01) + rrt(032) * density(14) 
  pd(32,14) = pd(32,14) + rrt(032) * density(01) 
  pd(06,01) = pd(06,01) + rrt(033) * density(14) 
  pd(06,14) = pd(06,14) + rrt(033) * density(01) 
  pd(14,01) = pd(14,01) - rrt(033) * density(14) 
  pd(14,14) = pd(14,14) - rrt(033) * density(01) 
  pd(47,01) = pd(47,01) + rrt(033) * density(14) 
  pd(47,14) = pd(47,14) + rrt(033) * density(01) 
  pd(11,01) = pd(11,01) + rrt(034) * density(14) * 2.0d0
  pd(11,14) = pd(11,14) + rrt(034) * density(01) * 2.0d0
  pd(14,01) = pd(14,01) - rrt(034) * density(14) 
  pd(14,14) = pd(14,14) - rrt(034) * density(01) 
  pd(08,01) = pd(08,01) - rrt(035) * density(08) 
  pd(08,08) = pd(08,08) - rrt(035) * density(01) 
  pd(30,01) = pd(30,01) + rrt(035) * density(08) 
  pd(30,08) = pd(30,08) + rrt(035) * density(01) 
  pd(34,01) = pd(34,01) + rrt(035) * density(08) 
  pd(34,08) = pd(34,08) + rrt(035) * density(01) 
  pd(08,01) = pd(08,01) - rrt(036) * density(08) 
  pd(08,08) = pd(08,08) - rrt(036) * density(01) 
  pd(23,01) = pd(23,01) + rrt(036) * density(08) 
  pd(23,08) = pd(23,08) + rrt(036) * density(01) 
  pd(42,01) = pd(42,01) + rrt(036) * density(08) 
  pd(42,08) = pd(42,08) + rrt(036) * density(01) 
  pd(08,01) = pd(08,01) - rrt(037) * density(08) 
  pd(08,08) = pd(08,08) - rrt(037) * density(01) 
  pd(34,01) = pd(34,01) + rrt(037) * density(08) * 2.0d0
  pd(34,08) = pd(34,08) + rrt(037) * density(01) * 2.0d0
  pd(42,01) = pd(42,01) + rrt(037) * density(08) 
  pd(42,08) = pd(42,08) + rrt(037) * density(01) 
  pd(08,01) = pd(08,01) - rrt(038) * density(08) 
  pd(08,08) = pd(08,08) - rrt(038) * density(01) 
  pd(23,01) = pd(23,01) + rrt(038) * density(08) 
  pd(23,08) = pd(23,08) + rrt(038) * density(01) 
  pd(32,01) = pd(32,01) + rrt(038) * density(08) 
  pd(32,08) = pd(32,08) + rrt(038) * density(01) 
  pd(34,01) = pd(34,01) + rrt(038) * density(08) 
  pd(34,08) = pd(34,08) + rrt(038) * density(01) 
  pd(06,01) = pd(06,01) + rrt(039) * density(08) 
  pd(06,08) = pd(06,08) + rrt(039) * density(01) 
  pd(08,01) = pd(08,01) - rrt(039) * density(08) 
  pd(08,08) = pd(08,08) - rrt(039) * density(01) 
  pd(43,01) = pd(43,01) + rrt(039) * density(08) 
  pd(43,08) = pd(43,08) + rrt(039) * density(01) 
  pd(08,01) = pd(08,01) - rrt(040) * density(08) 
  pd(08,08) = pd(08,08) - rrt(040) * density(01) 
  pd(11,01) = pd(11,01) + rrt(040) * density(08) 
  pd(11,08) = pd(11,08) + rrt(040) * density(01) 
  pd(47,01) = pd(47,01) + rrt(040) * density(08) 
  pd(47,08) = pd(47,08) + rrt(040) * density(01) 
  pd(30,01) = pd(30,01) - rrt(041) * density(30) 
  pd(30,30) = pd(30,30) - rrt(041) * density(01) 
  pd(34,01) = pd(34,01) + rrt(041) * density(30) 
  pd(34,30) = pd(34,30) + rrt(041) * density(01) 
  pd(42,01) = pd(42,01) + rrt(041) * density(30) 
  pd(42,30) = pd(42,30) + rrt(041) * density(01) 
  pd(23,01) = pd(23,01) + rrt(042) * density(30) 
  pd(23,30) = pd(23,30) + rrt(042) * density(01) 
  pd(30,01) = pd(30,01) - rrt(042) * density(30) 
  pd(30,30) = pd(30,30) - rrt(042) * density(01) 
  pd(32,01) = pd(32,01) + rrt(042) * density(30) 
  pd(32,30) = pd(32,30) + rrt(042) * density(01) 
  pd(30,01) = pd(30,01) - rrt(043) * density(30) 
  pd(30,30) = pd(30,30) - rrt(043) * density(01) 
  pd(32,01) = pd(32,01) + rrt(043) * density(30) 
  pd(32,30) = pd(32,30) + rrt(043) * density(01) 
  pd(34,01) = pd(34,01) + rrt(043) * density(30) * 2.0d0
  pd(34,30) = pd(34,30) + rrt(043) * density(01) * 2.0d0
  pd(11,01) = pd(11,01) + rrt(044) * density(30) 
  pd(11,30) = pd(11,30) + rrt(044) * density(01) 
  pd(30,01) = pd(30,01) - rrt(044) * density(30) 
  pd(30,30) = pd(30,30) - rrt(044) * density(01) 
  pd(43,01) = pd(43,01) + rrt(044) * density(30) 
  pd(43,30) = pd(43,30) + rrt(044) * density(01) 
  pd(30,01) = pd(30,01) - rrt(045) * density(30) 
  pd(30,30) = pd(30,30) - rrt(045) * density(01) 
  pd(47,01) = pd(47,01) + rrt(045) * density(30) * 2.0d0
  pd(47,30) = pd(47,30) + rrt(045) * density(01) * 2.0d0
  pd(32,01) = pd(32,01) + rrt(046) * density(42) 
  pd(32,42) = pd(32,42) + rrt(046) * density(01) 
  pd(34,01) = pd(34,01) + rrt(046) * density(42) 
  pd(34,42) = pd(34,42) + rrt(046) * density(01) 
  pd(42,01) = pd(42,01) - rrt(046) * density(42) 
  pd(42,42) = pd(42,42) - rrt(046) * density(01) 
  pd(42,01) = pd(42,01) - rrt(047) * density(42) 
  pd(42,42) = pd(42,42) - rrt(047) * density(01) 
  pd(43,01) = pd(43,01) + rrt(047) * density(42) 
  pd(43,42) = pd(43,42) + rrt(047) * density(01) 
  pd(47,01) = pd(47,01) + rrt(047) * density(42) 
  pd(47,42) = pd(47,42) + rrt(047) * density(01) 
  pd(32,01) = pd(32,01) - rrt(048) * density(32) 
  pd(32,32) = pd(32,32) - rrt(048) * density(01) 
  pd(43,01) = pd(43,01) + rrt(048) * density(32) * 2.0d0
  pd(43,32) = pd(43,32) + rrt(048) * density(01) * 2.0d0
  pd(07,01) = pd(07,01) - rrt(049) * density(07) 
  pd(07,07) = pd(07,07) - rrt(049) * density(01) 
  pd(34,01) = pd(34,01) + rrt(049) * density(07) 
  pd(34,07) = pd(34,07) + rrt(049) * density(01) 
  pd(44,01) = pd(44,01) + rrt(049) * density(07) 
  pd(44,07) = pd(44,07) + rrt(049) * density(01) 
  pd(07,01) = pd(07,01) - rrt(050) * density(07) 
  pd(07,07) = pd(07,07) - rrt(050) * density(01) 
  pd(21,01) = pd(21,01) + rrt(050) * density(07) 
  pd(21,07) = pd(21,07) + rrt(050) * density(01) 
  pd(23,01) = pd(23,01) + rrt(050) * density(07) 
  pd(23,07) = pd(23,07) + rrt(050) * density(01) 
  pd(07,01) = pd(07,01) - rrt(051) * density(07) 
  pd(07,07) = pd(07,07) - rrt(051) * density(01) 
  pd(17,01) = pd(17,01) + rrt(051) * density(07) 
  pd(17,07) = pd(17,07) + rrt(051) * density(01) 
  pd(23,01) = pd(23,01) + rrt(051) * density(07) * 2.0d0
  pd(23,07) = pd(23,07) + rrt(051) * density(01) * 2.0d0
  pd(07,01) = pd(07,01) - rrt(052) * density(07) 
  pd(07,07) = pd(07,07) - rrt(052) * density(01) 
  pd(14,01) = pd(14,01) + rrt(052) * density(07) 
  pd(14,07) = pd(14,07) + rrt(052) * density(01) 
  pd(47,01) = pd(47,01) + rrt(052) * density(07) 
  pd(47,07) = pd(47,07) + rrt(052) * density(01) 
  pd(07,01) = pd(07,01) - rrt(053) * density(07) 
  pd(07,07) = pd(07,07) - rrt(053) * density(01) 
  pd(08,01) = pd(08,01) + rrt(053) * density(07) 
  pd(08,07) = pd(08,07) + rrt(053) * density(01) 
  pd(11,01) = pd(11,01) + rrt(053) * density(07) 
  pd(11,07) = pd(11,07) + rrt(053) * density(01) 
  pd(06,01) = pd(06,01) + rrt(054) * density(07) 
  pd(06,07) = pd(06,07) + rrt(054) * density(01) 
  pd(07,01) = pd(07,01) - rrt(054) * density(07) 
  pd(07,07) = pd(07,07) - rrt(054) * density(01) 
  pd(30,01) = pd(30,01) + rrt(054) * density(07) 
  pd(30,07) = pd(30,07) + rrt(054) * density(01) 
  pd(21,01) = pd(21,01) + rrt(055) * density(44) 
  pd(21,44) = pd(21,44) + rrt(055) * density(01) 
  pd(34,01) = pd(34,01) + rrt(055) * density(44) 
  pd(34,44) = pd(34,44) + rrt(055) * density(01) 
  pd(44,01) = pd(44,01) - rrt(055) * density(44) 
  pd(44,44) = pd(44,44) - rrt(055) * density(01) 
  pd(23,01) = pd(23,01) + rrt(056) * density(44) 
  pd(23,44) = pd(23,44) + rrt(056) * density(01) 
  pd(26,01) = pd(26,01) + rrt(056) * density(44) 
  pd(26,44) = pd(26,44) + rrt(056) * density(01) 
  pd(44,01) = pd(44,01) - rrt(056) * density(44) 
  pd(44,44) = pd(44,44) - rrt(056) * density(01) 
  pd(17,01) = pd(17,01) + rrt(057) * density(44) 
  pd(17,44) = pd(17,44) + rrt(057) * density(01) 
  pd(23,01) = pd(23,01) + rrt(057) * density(44) 
  pd(23,44) = pd(23,44) + rrt(057) * density(01) 
  pd(34,01) = pd(34,01) + rrt(057) * density(44) 
  pd(34,44) = pd(34,44) + rrt(057) * density(01) 
  pd(44,01) = pd(44,01) - rrt(057) * density(44) 
  pd(44,44) = pd(44,44) - rrt(057) * density(01) 
  pd(11,01) = pd(11,01) + rrt(058) * density(44) 
  pd(11,44) = pd(11,44) + rrt(058) * density(01) 
  pd(30,01) = pd(30,01) + rrt(058) * density(44) 
  pd(30,44) = pd(30,44) + rrt(058) * density(01) 
  pd(44,01) = pd(44,01) - rrt(058) * density(44) 
  pd(44,44) = pd(44,44) - rrt(058) * density(01) 
  pd(06,01) = pd(06,01) + rrt(059) * density(44) 
  pd(06,44) = pd(06,44) + rrt(059) * density(01) 
  pd(42,01) = pd(42,01) + rrt(059) * density(44) 
  pd(42,44) = pd(42,44) + rrt(059) * density(01) 
  pd(44,01) = pd(44,01) - rrt(059) * density(44) 
  pd(44,44) = pd(44,44) - rrt(059) * density(01) 
  pd(21,01) = pd(21,01) - rrt(060) * density(21) 
  pd(21,21) = pd(21,21) - rrt(060) * density(01) 
  pd(26,01) = pd(26,01) + rrt(060) * density(21) 
  pd(26,21) = pd(26,21) + rrt(060) * density(01) 
  pd(34,01) = pd(34,01) + rrt(060) * density(21) 
  pd(34,21) = pd(34,21) + rrt(060) * density(01) 
  pd(17,01) = pd(17,01) + rrt(061) * density(21) 
  pd(17,21) = pd(17,21) + rrt(061) * density(01) 
  pd(21,01) = pd(21,01) - rrt(061) * density(21) 
  pd(21,21) = pd(21,21) - rrt(061) * density(01) 
  pd(23,01) = pd(23,01) + rrt(061) * density(21) 
  pd(23,21) = pd(23,21) + rrt(061) * density(01) 
  pd(21,01) = pd(21,01) - rrt(062) * density(21) 
  pd(21,21) = pd(21,21) - rrt(062) * density(01) 
  pd(30,01) = pd(30,01) + rrt(062) * density(21) 
  pd(30,21) = pd(30,21) + rrt(062) * density(01) 
  pd(47,01) = pd(47,01) + rrt(062) * density(21) 
  pd(47,21) = pd(47,21) + rrt(062) * density(01) 
  pd(11,01) = pd(11,01) + rrt(063) * density(21) 
  pd(11,21) = pd(11,21) + rrt(063) * density(01) 
  pd(21,01) = pd(21,01) - rrt(063) * density(21) 
  pd(21,21) = pd(21,21) - rrt(063) * density(01) 
  pd(42,01) = pd(42,01) + rrt(063) * density(21) 
  pd(42,21) = pd(42,21) + rrt(063) * density(01) 
  pd(06,01) = pd(06,01) + rrt(064) * density(21) 
  pd(06,21) = pd(06,21) + rrt(064) * density(01) 
  pd(21,01) = pd(21,01) - rrt(064) * density(21) 
  pd(21,21) = pd(21,21) - rrt(064) * density(01) 
  pd(32,01) = pd(32,01) + rrt(064) * density(21) 
  pd(32,21) = pd(32,21) + rrt(064) * density(01) 
  pd(17,01) = pd(17,01) + rrt(065) * density(26) 
  pd(17,26) = pd(17,26) + rrt(065) * density(01) 
  pd(26,01) = pd(26,01) - rrt(065) * density(26) 
  pd(26,26) = pd(26,26) - rrt(065) * density(01) 
  pd(34,01) = pd(34,01) + rrt(065) * density(26) 
  pd(34,26) = pd(34,26) + rrt(065) * density(01) 
  pd(11,01) = pd(11,01) + rrt(066) * density(26) 
  pd(11,26) = pd(11,26) + rrt(066) * density(01) 
  pd(26,01) = pd(26,01) - rrt(066) * density(26) 
  pd(26,26) = pd(26,26) - rrt(066) * density(01) 
  pd(32,01) = pd(32,01) + rrt(066) * density(26) 
  pd(32,26) = pd(32,26) + rrt(066) * density(01) 
  pd(17,01) = pd(17,01) - rrt(067) * density(17) 
  pd(17,17) = pd(17,17) - rrt(067) * density(01) 
  pd(42,01) = pd(42,01) + rrt(067) * density(17) 
  pd(42,17) = pd(42,17) + rrt(067) * density(01) 
  pd(43,01) = pd(43,01) + rrt(067) * density(17) 
  pd(43,17) = pd(43,17) + rrt(067) * density(01) 
  pd(17,01) = pd(17,01) - rrt(068) * density(17) 
  pd(17,17) = pd(17,17) - rrt(068) * density(01) 
  pd(32,01) = pd(32,01) + rrt(068) * density(17) 
  pd(32,17) = pd(32,17) + rrt(068) * density(01) 
  pd(47,01) = pd(47,01) + rrt(068) * density(17) 
  pd(47,17) = pd(47,17) + rrt(068) * density(01) 
  pd(01,01) = pd(01,01) + rrt(069) * density(14) 
  pd(01,14) = pd(01,14) + rrt(069) * density(01) 
  pd(14,01) = pd(14,01) - rrt(069) * density(14) 
  pd(14,14) = pd(14,14) - rrt(069) * density(01) 
  pd(52,01) = pd(52,01) + rrt(069) * density(14) 
  pd(52,14) = pd(52,14) + rrt(069) * density(01) 
  pd(01,01) = pd(01,01) + rrt(070) * density(14) 
  pd(01,14) = pd(01,14) + rrt(070) * density(01) 
  pd(14,01) = pd(14,01) - rrt(070) * density(14) 
  pd(14,14) = pd(14,14) - rrt(070) * density(01) 
  pd(27,01) = pd(27,01) + rrt(070) * density(14) 
  pd(27,14) = pd(27,14) + rrt(070) * density(01) 
  pd(34,01) = pd(34,01) + rrt(070) * density(14) 
  pd(34,14) = pd(34,14) + rrt(070) * density(01) 
  pd(01,01) = pd(01,01) + rrt(071) * density(14) 
  pd(01,14) = pd(01,14) + rrt(071) * density(01) 
  pd(14,01) = pd(14,01) - rrt(071) * density(14) 
  pd(14,14) = pd(14,14) - rrt(071) * density(01) 
  pd(16,01) = pd(16,01) + rrt(071) * density(14) 
  pd(16,14) = pd(16,14) + rrt(071) * density(01) 
  pd(23,01) = pd(23,01) + rrt(071) * density(14) 
  pd(23,14) = pd(23,14) + rrt(071) * density(01) 
  pd(01,01) = pd(01,01) + rrt(072) * density(14) 
  pd(01,14) = pd(01,14) + rrt(072) * density(01) 
  pd(04,01) = pd(04,01) + rrt(072) * density(14) 
  pd(04,14) = pd(04,14) + rrt(072) * density(01) 
  pd(14,01) = pd(14,01) - rrt(072) * density(14) 
  pd(14,14) = pd(14,14) - rrt(072) * density(01) 
  pd(23,01) = pd(23,01) + rrt(072) * density(14) 
  pd(23,14) = pd(23,14) + rrt(072) * density(01) 
  pd(34,01) = pd(34,01) + rrt(072) * density(14) 
  pd(34,14) = pd(34,14) + rrt(072) * density(01) 
  pd(01,01) = pd(01,01) + rrt(073) * density(14) 
  pd(01,14) = pd(01,14) + rrt(073) * density(01) 
  pd(09,01) = pd(09,01) + rrt(073) * density(14) 
  pd(09,14) = pd(09,14) + rrt(073) * density(01) 
  pd(14,01) = pd(14,01) - rrt(073) * density(14) 
  pd(14,14) = pd(14,14) - rrt(073) * density(01) 
  pd(23,01) = pd(23,01) + rrt(073) * density(14) * 2.0d0
  pd(23,14) = pd(23,14) + rrt(073) * density(01) * 2.0d0
  pd(01,01) = pd(01,01) + rrt(074) * density(14) 
  pd(01,14) = pd(01,14) + rrt(074) * density(01) 
  pd(11,01) = pd(11,01) + rrt(074) * density(14) 
  pd(11,14) = pd(11,14) + rrt(074) * density(01) 
  pd(14,01) = pd(14,01) - rrt(074) * density(14) 
  pd(14,14) = pd(14,14) - rrt(074) * density(01) 
  pd(49,01) = pd(49,01) + rrt(074) * density(14) 
  pd(49,14) = pd(49,14) + rrt(074) * density(01) 
  pd(01,01) = pd(01,01) + rrt(075) * density(14) 
  pd(01,14) = pd(01,14) + rrt(075) * density(01) 
  pd(06,01) = pd(06,01) + rrt(075) * density(14) 
  pd(06,14) = pd(06,14) + rrt(075) * density(01) 
  pd(14,01) = pd(14,01) - rrt(075) * density(14) 
  pd(14,14) = pd(14,14) - rrt(075) * density(01) 
  pd(46,01) = pd(46,01) + rrt(075) * density(14) 
  pd(46,14) = pd(46,14) + rrt(075) * density(01) 
  pd(01,01) = pd(01,01) + rrt(076) * density(08) 
  pd(01,08) = pd(01,08) + rrt(076) * density(01) 
  pd(08,01) = pd(08,01) - rrt(076) * density(08) 
  pd(08,08) = pd(08,08) - rrt(076) * density(01) 
  pd(27,01) = pd(27,01) + rrt(076) * density(08) 
  pd(27,08) = pd(27,08) + rrt(076) * density(01) 
  pd(01,01) = pd(01,01) + rrt(077) * density(08) 
  pd(01,08) = pd(01,08) + rrt(077) * density(01) 
  pd(08,01) = pd(08,01) - rrt(077) * density(08) 
  pd(08,08) = pd(08,08) - rrt(077) * density(01) 
  pd(16,01) = pd(16,01) + rrt(077) * density(08) 
  pd(16,08) = pd(16,08) + rrt(077) * density(01) 
  pd(34,01) = pd(34,01) + rrt(077) * density(08) 
  pd(34,08) = pd(34,08) + rrt(077) * density(01) 
  pd(01,01) = pd(01,01) + rrt(078) * density(08) 
  pd(01,08) = pd(01,08) + rrt(078) * density(01) 
  pd(04,01) = pd(04,01) + rrt(078) * density(08) 
  pd(04,08) = pd(04,08) + rrt(078) * density(01) 
  pd(08,01) = pd(08,01) - rrt(078) * density(08) 
  pd(08,08) = pd(08,08) - rrt(078) * density(01) 
  pd(23,01) = pd(23,01) + rrt(078) * density(08) 
  pd(23,08) = pd(23,08) + rrt(078) * density(01) 
  pd(01,01) = pd(01,01) + rrt(079) * density(08) 
  pd(01,08) = pd(01,08) + rrt(079) * density(01) 
  pd(08,01) = pd(08,01) - rrt(079) * density(08) 
  pd(08,08) = pd(08,08) - rrt(079) * density(01) 
  pd(09,01) = pd(09,01) + rrt(079) * density(08) 
  pd(09,08) = pd(09,08) + rrt(079) * density(01) 
  pd(23,01) = pd(23,01) + rrt(079) * density(08) 
  pd(23,08) = pd(23,08) + rrt(079) * density(01) 
  pd(34,01) = pd(34,01) + rrt(079) * density(08) 
  pd(34,08) = pd(34,08) + rrt(079) * density(01) 
  pd(01,01) = pd(01,01) + rrt(080) * density(08) 
  pd(01,08) = pd(01,08) + rrt(080) * density(01) 
  pd(08,01) = pd(08,01) - rrt(080) * density(08) 
  pd(08,08) = pd(08,08) - rrt(080) * density(01) 
  pd(47,01) = pd(47,01) + rrt(080) * density(08) 
  pd(47,08) = pd(47,08) + rrt(080) * density(01) 
  pd(49,01) = pd(49,01) + rrt(080) * density(08) 
  pd(49,08) = pd(49,08) + rrt(080) * density(01) 
  pd(01,01) = pd(01,01) + rrt(081) * density(08) 
  pd(01,08) = pd(01,08) + rrt(081) * density(01) 
  pd(08,01) = pd(08,01) - rrt(081) * density(08) 
  pd(08,08) = pd(08,08) - rrt(081) * density(01) 
  pd(11,01) = pd(11,01) + rrt(081) * density(08) 
  pd(11,08) = pd(11,08) + rrt(081) * density(01) 
  pd(46,01) = pd(46,01) + rrt(081) * density(08) 
  pd(46,08) = pd(46,08) + rrt(081) * density(01) 
  pd(01,01) = pd(01,01) + rrt(082) * density(08) 
  pd(01,08) = pd(01,08) + rrt(082) * density(01) 
  pd(06,01) = pd(06,01) + rrt(082) * density(08) 
  pd(06,08) = pd(06,08) + rrt(082) * density(01) 
  pd(08,01) = pd(08,01) - rrt(082) * density(08) 
  pd(08,08) = pd(08,08) - rrt(082) * density(01) 
  pd(13,01) = pd(13,01) + rrt(082) * density(08) 
  pd(13,08) = pd(13,08) + rrt(082) * density(01) 
  pd(01,01) = pd(01,01) + rrt(083) * density(30) 
  pd(01,30) = pd(01,30) + rrt(083) * density(01) 
  pd(16,01) = pd(16,01) + rrt(083) * density(30) 
  pd(16,30) = pd(16,30) + rrt(083) * density(01) 
  pd(30,01) = pd(30,01) - rrt(083) * density(30) 
  pd(30,30) = pd(30,30) - rrt(083) * density(01) 
  pd(01,01) = pd(01,01) + rrt(084) * density(30) 
  pd(01,30) = pd(01,30) + rrt(084) * density(01) 
  pd(04,01) = pd(04,01) + rrt(084) * density(30) 
  pd(04,30) = pd(04,30) + rrt(084) * density(01) 
  pd(30,01) = pd(30,01) - rrt(084) * density(30) 
  pd(30,30) = pd(30,30) - rrt(084) * density(01) 
  pd(34,01) = pd(34,01) + rrt(084) * density(30) 
  pd(34,30) = pd(34,30) + rrt(084) * density(01) 
  pd(01,01) = pd(01,01) + rrt(085) * density(30) 
  pd(01,30) = pd(01,30) + rrt(085) * density(01) 
  pd(30,01) = pd(30,01) - rrt(085) * density(30) 
  pd(30,30) = pd(30,30) - rrt(085) * density(01) 
  pd(43,01) = pd(43,01) + rrt(085) * density(30) 
  pd(43,30) = pd(43,30) + rrt(085) * density(01) 
  pd(49,01) = pd(49,01) + rrt(085) * density(30) 
  pd(49,30) = pd(49,30) + rrt(085) * density(01) 
  pd(01,01) = pd(01,01) + rrt(086) * density(30) 
  pd(01,30) = pd(01,30) + rrt(086) * density(01) 
  pd(30,01) = pd(30,01) - rrt(086) * density(30) 
  pd(30,30) = pd(30,30) - rrt(086) * density(01) 
  pd(46,01) = pd(46,01) + rrt(086) * density(30) 
  pd(46,30) = pd(46,30) + rrt(086) * density(01) 
  pd(47,01) = pd(47,01) + rrt(086) * density(30) 
  pd(47,30) = pd(47,30) + rrt(086) * density(01) 
  pd(01,01) = pd(01,01) + rrt(087) * density(30) 
  pd(01,30) = pd(01,30) + rrt(087) * density(01) 
  pd(11,01) = pd(11,01) + rrt(087) * density(30) 
  pd(11,30) = pd(11,30) + rrt(087) * density(01) 
  pd(13,01) = pd(13,01) + rrt(087) * density(30) 
  pd(13,30) = pd(13,30) + rrt(087) * density(01) 
  pd(30,01) = pd(30,01) - rrt(087) * density(30) 
  pd(30,30) = pd(30,30) - rrt(087) * density(01) 
  pd(01,01) = pd(01,01) + rrt(088) * density(42) 
  pd(01,42) = pd(01,42) + rrt(088) * density(01) 
  pd(04,01) = pd(04,01) + rrt(088) * density(42) 
  pd(04,42) = pd(04,42) + rrt(088) * density(01) 
  pd(42,01) = pd(42,01) - rrt(088) * density(42) 
  pd(42,42) = pd(42,42) - rrt(088) * density(01) 
  pd(01,01) = pd(01,01) + rrt(089) * density(42) 
  pd(01,42) = pd(01,42) + rrt(089) * density(01) 
  pd(09,01) = pd(09,01) + rrt(089) * density(42) 
  pd(09,42) = pd(09,42) + rrt(089) * density(01) 
  pd(34,01) = pd(34,01) + rrt(089) * density(42) 
  pd(34,42) = pd(34,42) + rrt(089) * density(01) 
  pd(42,01) = pd(42,01) - rrt(089) * density(42) 
  pd(42,42) = pd(42,42) - rrt(089) * density(01) 
  pd(01,01) = pd(01,01) + rrt(090) * density(42) 
  pd(01,42) = pd(01,42) + rrt(090) * density(01) 
  pd(42,01) = pd(42,01) - rrt(090) * density(42) 
  pd(42,42) = pd(42,42) - rrt(090) * density(01) 
  pd(43,01) = pd(43,01) + rrt(090) * density(42) 
  pd(43,42) = pd(43,42) + rrt(090) * density(01) 
  pd(46,01) = pd(46,01) + rrt(090) * density(42) 
  pd(46,42) = pd(46,42) + rrt(090) * density(01) 
  pd(01,01) = pd(01,01) + rrt(091) * density(42) 
  pd(01,42) = pd(01,42) + rrt(091) * density(01) 
  pd(13,01) = pd(13,01) + rrt(091) * density(42) 
  pd(13,42) = pd(13,42) + rrt(091) * density(01) 
  pd(42,01) = pd(42,01) - rrt(091) * density(42) 
  pd(42,42) = pd(42,42) - rrt(091) * density(01) 
  pd(47,01) = pd(47,01) + rrt(091) * density(42) 
  pd(47,42) = pd(47,42) + rrt(091) * density(01) 
  pd(01,01) = pd(01,01) + rrt(092) * density(42) 
  pd(01,42) = pd(01,42) + rrt(092) * density(01) 
  pd(12,01) = pd(12,01) + rrt(092) * density(42) 
  pd(12,42) = pd(12,42) + rrt(092) * density(01) 
  pd(32,01) = pd(32,01) + rrt(092) * density(42) 
  pd(32,42) = pd(32,42) + rrt(092) * density(01) 
  pd(42,01) = pd(42,01) - rrt(092) * density(42) 
  pd(42,42) = pd(42,42) - rrt(092) * density(01) 
  pd(01,01) = pd(01,01) + rrt(093) * density(32) 
  pd(01,32) = pd(01,32) + rrt(093) * density(01) 
  pd(09,01) = pd(09,01) + rrt(093) * density(32) 
  pd(09,32) = pd(09,32) + rrt(093) * density(01) 
  pd(32,01) = pd(32,01) - rrt(093) * density(32) 
  pd(32,32) = pd(32,32) - rrt(093) * density(01) 
  pd(01,01) = pd(01,01) + rrt(094) * density(32) 
  pd(01,32) = pd(01,32) + rrt(094) * density(01) 
  pd(13,01) = pd(13,01) + rrt(094) * density(32) 
  pd(13,32) = pd(13,32) + rrt(094) * density(01) 
  pd(32,01) = pd(32,01) - rrt(094) * density(32) 
  pd(32,32) = pd(32,32) - rrt(094) * density(01) 
  pd(43,01) = pd(43,01) + rrt(094) * density(32) 
  pd(43,32) = pd(43,32) + rrt(094) * density(01) 
  pd(01,01) = pd(01,01) + rrt(095) * density(07) 
  pd(01,07) = pd(01,07) + rrt(095) * density(01) 
  pd(07,01) = pd(07,01) - rrt(095) * density(07) 
  pd(07,07) = pd(07,07) - rrt(095) * density(01) 
  pd(41,01) = pd(41,01) + rrt(095) * density(07) 
  pd(41,07) = pd(41,07) + rrt(095) * density(01) 
  pd(01,01) = pd(01,01) + rrt(096) * density(07) 
  pd(01,07) = pd(01,07) + rrt(096) * density(01) 
  pd(07,01) = pd(07,01) - rrt(096) * density(07) 
  pd(07,07) = pd(07,07) - rrt(096) * density(01) 
  pd(18,01) = pd(18,01) + rrt(096) * density(07) 
  pd(18,07) = pd(18,07) + rrt(096) * density(01) 
  pd(34,01) = pd(34,01) + rrt(096) * density(07) 
  pd(34,07) = pd(34,07) + rrt(096) * density(01) 
  pd(01,01) = pd(01,01) + rrt(097) * density(07) 
  pd(01,07) = pd(01,07) + rrt(097) * density(01) 
  pd(07,01) = pd(07,01) - rrt(097) * density(07) 
  pd(07,07) = pd(07,07) - rrt(097) * density(01) 
  pd(23,01) = pd(23,01) + rrt(097) * density(07) 
  pd(23,07) = pd(23,07) + rrt(097) * density(01) 
  pd(38,01) = pd(38,01) + rrt(097) * density(07) 
  pd(38,07) = pd(38,07) + rrt(097) * density(01) 
  pd(01,01) = pd(01,01) + rrt(098) * density(07) 
  pd(01,07) = pd(01,07) + rrt(098) * density(01) 
  pd(07,01) = pd(07,01) - rrt(098) * density(07) 
  pd(07,07) = pd(07,07) - rrt(098) * density(01) 
  pd(23,01) = pd(23,01) + rrt(098) * density(07) 
  pd(23,07) = pd(23,07) + rrt(098) * density(01) 
  pd(34,01) = pd(34,01) + rrt(098) * density(07) 
  pd(34,07) = pd(34,07) + rrt(098) * density(01) 
  pd(39,01) = pd(39,01) + rrt(098) * density(07) 
  pd(39,07) = pd(39,07) + rrt(098) * density(01) 
  pd(01,01) = pd(01,01) + rrt(099) * density(07) 
  pd(01,07) = pd(01,07) + rrt(099) * density(01) 
  pd(07,01) = pd(07,01) - rrt(099) * density(07) 
  pd(07,07) = pd(07,07) - rrt(099) * density(01) 
  pd(20,01) = pd(20,01) + rrt(099) * density(07) 
  pd(20,07) = pd(20,07) + rrt(099) * density(01) 
  pd(23,01) = pd(23,01) + rrt(099) * density(07) * 2.0d0
  pd(23,07) = pd(23,07) + rrt(099) * density(01) * 2.0d0
  pd(01,01) = pd(01,01) + rrt(100) * density(07) 
  pd(01,07) = pd(01,07) + rrt(100) * density(01) 
  pd(07,01) = pd(07,01) - rrt(100) * density(07) 
  pd(07,07) = pd(07,07) - rrt(100) * density(01) 
  pd(11,01) = pd(11,01) + rrt(100) * density(07) 
  pd(11,07) = pd(11,07) + rrt(100) * density(01) 
  pd(27,01) = pd(27,01) + rrt(100) * density(07) 
  pd(27,07) = pd(27,07) + rrt(100) * density(01) 
  pd(01,01) = pd(01,01) + rrt(101) * density(07) 
  pd(01,07) = pd(01,07) + rrt(101) * density(01) 
  pd(06,01) = pd(06,01) + rrt(101) * density(07) 
  pd(06,07) = pd(06,07) + rrt(101) * density(01) 
  pd(07,01) = pd(07,01) - rrt(101) * density(07) 
  pd(07,07) = pd(07,07) - rrt(101) * density(01) 
  pd(16,01) = pd(16,01) + rrt(101) * density(07) 
  pd(16,07) = pd(16,07) + rrt(101) * density(01) 
  pd(01,01) = pd(01,01) + rrt(102) * density(07) 
  pd(01,07) = pd(01,07) + rrt(102) * density(01) 
  pd(07,01) = pd(07,01) - rrt(102) * density(07) 
  pd(07,07) = pd(07,07) - rrt(102) * density(01) 
  pd(08,01) = pd(08,01) + rrt(102) * density(07) 
  pd(08,07) = pd(08,07) + rrt(102) * density(01) 
  pd(49,01) = pd(49,01) + rrt(102) * density(07) 
  pd(49,07) = pd(49,07) + rrt(102) * density(01) 
  pd(01,01) = pd(01,01) + rrt(103) * density(07) 
  pd(01,07) = pd(01,07) + rrt(103) * density(01) 
  pd(07,01) = pd(07,01) - rrt(103) * density(07) 
  pd(07,07) = pd(07,07) - rrt(103) * density(01) 
  pd(14,01) = pd(14,01) + rrt(103) * density(07) 
  pd(14,07) = pd(14,07) + rrt(103) * density(01) 
  pd(46,01) = pd(46,01) + rrt(103) * density(07) 
  pd(46,07) = pd(46,07) + rrt(103) * density(01) 
  pd(01,01) = pd(01,01) + rrt(104) * density(44) 
  pd(01,44) = pd(01,44) + rrt(104) * density(01) 
  pd(18,01) = pd(18,01) + rrt(104) * density(44) 
  pd(18,44) = pd(18,44) + rrt(104) * density(01) 
  pd(44,01) = pd(44,01) - rrt(104) * density(44) 
  pd(44,44) = pd(44,44) - rrt(104) * density(01) 
  pd(01,01) = pd(01,01) + rrt(105) * density(44) 
  pd(01,44) = pd(01,44) + rrt(105) * density(01) 
  pd(34,01) = pd(34,01) + rrt(105) * density(44) 
  pd(34,44) = pd(34,44) + rrt(105) * density(01) 
  pd(38,01) = pd(38,01) + rrt(105) * density(44) 
  pd(38,44) = pd(38,44) + rrt(105) * density(01) 
  pd(44,01) = pd(44,01) - rrt(105) * density(44) 
  pd(44,44) = pd(44,44) - rrt(105) * density(01) 
  pd(01,01) = pd(01,01) + rrt(106) * density(44) 
  pd(01,44) = pd(01,44) + rrt(106) * density(01) 
  pd(23,01) = pd(23,01) + rrt(106) * density(44) 
  pd(23,44) = pd(23,44) + rrt(106) * density(01) 
  pd(39,01) = pd(39,01) + rrt(106) * density(44) 
  pd(39,44) = pd(39,44) + rrt(106) * density(01) 
  pd(44,01) = pd(44,01) - rrt(106) * density(44) 
  pd(44,44) = pd(44,44) - rrt(106) * density(01) 
  pd(01,01) = pd(01,01) + rrt(107) * density(44) 
  pd(01,44) = pd(01,44) + rrt(107) * density(01) 
  pd(20,01) = pd(20,01) + rrt(107) * density(44) 
  pd(20,44) = pd(20,44) + rrt(107) * density(01) 
  pd(23,01) = pd(23,01) + rrt(107) * density(44) 
  pd(23,44) = pd(23,44) + rrt(107) * density(01) 
  pd(34,01) = pd(34,01) + rrt(107) * density(44) 
  pd(34,44) = pd(34,44) + rrt(107) * density(01) 
  pd(44,01) = pd(44,01) - rrt(107) * density(44) 
  pd(44,44) = pd(44,44) - rrt(107) * density(01) 
  pd(01,01) = pd(01,01) + rrt(108) * density(44) 
  pd(01,44) = pd(01,44) + rrt(108) * density(01) 
  pd(27,01) = pd(27,01) + rrt(108) * density(44) 
  pd(27,44) = pd(27,44) + rrt(108) * density(01) 
  pd(44,01) = pd(44,01) - rrt(108) * density(44) 
  pd(44,44) = pd(44,44) - rrt(108) * density(01) 
  pd(47,01) = pd(47,01) + rrt(108) * density(44) 
  pd(47,44) = pd(47,44) + rrt(108) * density(01) 
  pd(01,01) = pd(01,01) + rrt(109) * density(44) 
  pd(01,44) = pd(01,44) + rrt(109) * density(01) 
  pd(11,01) = pd(11,01) + rrt(109) * density(44) 
  pd(11,44) = pd(11,44) + rrt(109) * density(01) 
  pd(16,01) = pd(16,01) + rrt(109) * density(44) 
  pd(16,44) = pd(16,44) + rrt(109) * density(01) 
  pd(44,01) = pd(44,01) - rrt(109) * density(44) 
  pd(44,44) = pd(44,44) - rrt(109) * density(01) 
  pd(01,01) = pd(01,01) + rrt(110) * density(44) 
  pd(01,44) = pd(01,44) + rrt(110) * density(01) 
  pd(04,01) = pd(04,01) + rrt(110) * density(44) 
  pd(04,44) = pd(04,44) + rrt(110) * density(01) 
  pd(06,01) = pd(06,01) + rrt(110) * density(44) 
  pd(06,44) = pd(06,44) + rrt(110) * density(01) 
  pd(44,01) = pd(44,01) - rrt(110) * density(44) 
  pd(44,44) = pd(44,44) - rrt(110) * density(01) 
  pd(01,01) = pd(01,01) + rrt(111) * density(44) 
  pd(01,44) = pd(01,44) + rrt(111) * density(01) 
  pd(42,01) = pd(42,01) + rrt(111) * density(44) 
  pd(42,44) = pd(42,44) + rrt(111) * density(01) 
  pd(44,01) = pd(44,01) - rrt(111) * density(44) 
  pd(44,44) = pd(44,44) - rrt(111) * density(01) 
  pd(48,01) = pd(48,01) + rrt(111) * density(44) 
  pd(48,44) = pd(48,44) + rrt(111) * density(01) 
  pd(01,01) = pd(01,01) + rrt(112) * density(44) 
  pd(01,44) = pd(01,44) + rrt(112) * density(01) 
  pd(30,01) = pd(30,01) + rrt(112) * density(44) 
  pd(30,44) = pd(30,44) + rrt(112) * density(01) 
  pd(44,01) = pd(44,01) - rrt(112) * density(44) 
  pd(44,44) = pd(44,44) - rrt(112) * density(01) 
  pd(49,01) = pd(49,01) + rrt(112) * density(44) 
  pd(49,44) = pd(49,44) + rrt(112) * density(01) 
  pd(01,01) = pd(01,01) + rrt(113) * density(44) 
  pd(01,44) = pd(01,44) + rrt(113) * density(01) 
  pd(08,01) = pd(08,01) + rrt(113) * density(44) 
  pd(08,44) = pd(08,44) + rrt(113) * density(01) 
  pd(44,01) = pd(44,01) - rrt(113) * density(44) 
  pd(44,44) = pd(44,44) - rrt(113) * density(01) 
  pd(46,01) = pd(46,01) + rrt(113) * density(44) 
  pd(46,44) = pd(46,44) + rrt(113) * density(01) 
  pd(01,01) = pd(01,01) + rrt(114) * density(44) 
  pd(01,44) = pd(01,44) + rrt(114) * density(01) 
  pd(13,01) = pd(13,01) + rrt(114) * density(44) 
  pd(13,44) = pd(13,44) + rrt(114) * density(01) 
  pd(14,01) = pd(14,01) + rrt(114) * density(44) 
  pd(14,44) = pd(14,44) + rrt(114) * density(01) 
  pd(44,01) = pd(44,01) - rrt(114) * density(44) 
  pd(44,44) = pd(44,44) - rrt(114) * density(01) 
  pd(01,01) = pd(01,01) + rrt(115) * density(21) 
  pd(01,21) = pd(01,21) + rrt(115) * density(01) 
  pd(21,01) = pd(21,01) - rrt(115) * density(21) 
  pd(21,21) = pd(21,21) - rrt(115) * density(01) 
  pd(38,01) = pd(38,01) + rrt(115) * density(21) 
  pd(38,21) = pd(38,21) + rrt(115) * density(01) 
  pd(01,01) = pd(01,01) + rrt(116) * density(21) 
  pd(01,21) = pd(01,21) + rrt(116) * density(01) 
  pd(21,01) = pd(21,01) - rrt(116) * density(21) 
  pd(21,21) = pd(21,21) - rrt(116) * density(01) 
  pd(34,01) = pd(34,01) + rrt(116) * density(21) 
  pd(34,21) = pd(34,21) + rrt(116) * density(01) 
  pd(39,01) = pd(39,01) + rrt(116) * density(21) 
  pd(39,21) = pd(39,21) + rrt(116) * density(01) 
  pd(01,01) = pd(01,01) + rrt(117) * density(21) 
  pd(01,21) = pd(01,21) + rrt(117) * density(01) 
  pd(20,01) = pd(20,01) + rrt(117) * density(21) 
  pd(20,21) = pd(20,21) + rrt(117) * density(01) 
  pd(21,01) = pd(21,01) - rrt(117) * density(21) 
  pd(21,21) = pd(21,21) - rrt(117) * density(01) 
  pd(23,01) = pd(23,01) + rrt(117) * density(21) 
  pd(23,21) = pd(23,21) + rrt(117) * density(01) 
  pd(01,01) = pd(01,01) + rrt(118) * density(21) 
  pd(01,21) = pd(01,21) + rrt(118) * density(01) 
  pd(21,01) = pd(21,01) - rrt(118) * density(21) 
  pd(21,21) = pd(21,21) - rrt(118) * density(01) 
  pd(27,01) = pd(27,01) + rrt(118) * density(21) 
  pd(27,21) = pd(27,21) + rrt(118) * density(01) 
  pd(43,01) = pd(43,01) + rrt(118) * density(21) 
  pd(43,21) = pd(43,21) + rrt(118) * density(01) 
  pd(01,01) = pd(01,01) + rrt(119) * density(21) 
  pd(01,21) = pd(01,21) + rrt(119) * density(01) 
  pd(16,01) = pd(16,01) + rrt(119) * density(21) 
  pd(16,21) = pd(16,21) + rrt(119) * density(01) 
  pd(21,01) = pd(21,01) - rrt(119) * density(21) 
  pd(21,21) = pd(21,21) - rrt(119) * density(01) 
  pd(47,01) = pd(47,01) + rrt(119) * density(21) 
  pd(47,21) = pd(47,21) + rrt(119) * density(01) 
  pd(01,01) = pd(01,01) + rrt(120) * density(21) 
  pd(01,21) = pd(01,21) + rrt(120) * density(01) 
  pd(04,01) = pd(04,01) + rrt(120) * density(21) 
  pd(04,21) = pd(04,21) + rrt(120) * density(01) 
  pd(11,01) = pd(11,01) + rrt(120) * density(21) 
  pd(11,21) = pd(11,21) + rrt(120) * density(01) 
  pd(21,01) = pd(21,01) - rrt(120) * density(21) 
  pd(21,21) = pd(21,21) - rrt(120) * density(01) 
  pd(01,01) = pd(01,01) + rrt(121) * density(21) 
  pd(01,21) = pd(01,21) + rrt(121) * density(01) 
  pd(06,01) = pd(06,01) + rrt(121) * density(21) 
  pd(06,21) = pd(06,21) + rrt(121) * density(01) 
  pd(09,01) = pd(09,01) + rrt(121) * density(21) 
  pd(09,21) = pd(09,21) + rrt(121) * density(01) 
  pd(21,01) = pd(21,01) - rrt(121) * density(21) 
  pd(21,21) = pd(21,21) - rrt(121) * density(01) 
  pd(01,01) = pd(01,01) + rrt(122) * density(21) 
  pd(01,21) = pd(01,21) + rrt(122) * density(01) 
  pd(21,01) = pd(21,01) - rrt(122) * density(21) 
  pd(21,21) = pd(21,21) - rrt(122) * density(01) 
  pd(32,01) = pd(32,01) + rrt(122) * density(21) 
  pd(32,21) = pd(32,21) + rrt(122) * density(01) 
  pd(48,01) = pd(48,01) + rrt(122) * density(21) 
  pd(48,21) = pd(48,21) + rrt(122) * density(01) 
  pd(01,01) = pd(01,01) + rrt(123) * density(21) 
  pd(01,21) = pd(01,21) + rrt(123) * density(01) 
  pd(21,01) = pd(21,01) - rrt(123) * density(21) 
  pd(21,21) = pd(21,21) - rrt(123) * density(01) 
  pd(42,01) = pd(42,01) + rrt(123) * density(21) 
  pd(42,21) = pd(42,21) + rrt(123) * density(01) 
  pd(49,01) = pd(49,01) + rrt(123) * density(21) 
  pd(49,21) = pd(49,21) + rrt(123) * density(01) 
  pd(01,01) = pd(01,01) + rrt(124) * density(21) 
  pd(01,21) = pd(01,21) + rrt(124) * density(01) 
  pd(21,01) = pd(21,01) - rrt(124) * density(21) 
  pd(21,21) = pd(21,21) - rrt(124) * density(01) 
  pd(30,01) = pd(30,01) + rrt(124) * density(21) 
  pd(30,21) = pd(30,21) + rrt(124) * density(01) 
  pd(46,01) = pd(46,01) + rrt(124) * density(21) 
  pd(46,21) = pd(46,21) + rrt(124) * density(01) 
  pd(01,01) = pd(01,01) + rrt(125) * density(21) 
  pd(01,21) = pd(01,21) + rrt(125) * density(01) 
  pd(08,01) = pd(08,01) + rrt(125) * density(21) 
  pd(08,21) = pd(08,21) + rrt(125) * density(01) 
  pd(13,01) = pd(13,01) + rrt(125) * density(21) 
  pd(13,21) = pd(13,21) + rrt(125) * density(01) 
  pd(21,01) = pd(21,01) - rrt(125) * density(21) 
  pd(21,21) = pd(21,21) - rrt(125) * density(01) 
  pd(01,01) = pd(01,01) + rrt(126) * density(26) 
  pd(01,26) = pd(01,26) + rrt(126) * density(01) 
  pd(26,01) = pd(26,01) - rrt(126) * density(26) 
  pd(26,26) = pd(26,26) - rrt(126) * density(01) 
  pd(39,01) = pd(39,01) + rrt(126) * density(26) 
  pd(39,26) = pd(39,26) + rrt(126) * density(01) 
  pd(01,01) = pd(01,01) + rrt(127) * density(26) 
  pd(01,26) = pd(01,26) + rrt(127) * density(01) 
  pd(20,01) = pd(20,01) + rrt(127) * density(26) 
  pd(20,26) = pd(20,26) + rrt(127) * density(01) 
  pd(26,01) = pd(26,01) - rrt(127) * density(26) 
  pd(26,26) = pd(26,26) - rrt(127) * density(01) 
  pd(34,01) = pd(34,01) + rrt(127) * density(26) 
  pd(34,26) = pd(34,26) + rrt(127) * density(01) 
  pd(01,01) = pd(01,01) + rrt(128) * density(26) 
  pd(01,26) = pd(01,26) + rrt(128) * density(01) 
  pd(16,01) = pd(16,01) + rrt(128) * density(26) 
  pd(16,26) = pd(16,26) + rrt(128) * density(01) 
  pd(26,01) = pd(26,01) - rrt(128) * density(26) 
  pd(26,26) = pd(26,26) - rrt(128) * density(01) 
  pd(43,01) = pd(43,01) + rrt(128) * density(26) 
  pd(43,26) = pd(43,26) + rrt(128) * density(01) 
  pd(01,01) = pd(01,01) + rrt(129) * density(26) 
  pd(01,26) = pd(01,26) + rrt(129) * density(01) 
  pd(04,01) = pd(04,01) + rrt(129) * density(26) 
  pd(04,26) = pd(04,26) + rrt(129) * density(01) 
  pd(26,01) = pd(26,01) - rrt(129) * density(26) 
  pd(26,26) = pd(26,26) - rrt(129) * density(01) 
  pd(47,01) = pd(47,01) + rrt(129) * density(26) 
  pd(47,26) = pd(47,26) + rrt(129) * density(01) 
  pd(01,01) = pd(01,01) + rrt(130) * density(26) 
  pd(01,26) = pd(01,26) + rrt(130) * density(01) 
  pd(09,01) = pd(09,01) + rrt(130) * density(26) 
  pd(09,26) = pd(09,26) + rrt(130) * density(01) 
  pd(11,01) = pd(11,01) + rrt(130) * density(26) 
  pd(11,26) = pd(11,26) + rrt(130) * density(01) 
  pd(26,01) = pd(26,01) - rrt(130) * density(26) 
  pd(26,26) = pd(26,26) - rrt(130) * density(01) 
  pd(01,01) = pd(01,01) + rrt(131) * density(26) 
  pd(01,26) = pd(01,26) + rrt(131) * density(01) 
  pd(26,01) = pd(26,01) - rrt(131) * density(26) 
  pd(26,26) = pd(26,26) - rrt(131) * density(01) 
  pd(32,01) = pd(32,01) + rrt(131) * density(26) 
  pd(32,26) = pd(32,26) + rrt(131) * density(01) 
  pd(49,01) = pd(49,01) + rrt(131) * density(26) 
  pd(49,26) = pd(49,26) + rrt(131) * density(01) 
  pd(01,01) = pd(01,01) + rrt(132) * density(26) 
  pd(01,26) = pd(01,26) + rrt(132) * density(01) 
  pd(26,01) = pd(26,01) - rrt(132) * density(26) 
  pd(26,26) = pd(26,26) - rrt(132) * density(01) 
  pd(42,01) = pd(42,01) + rrt(132) * density(26) 
  pd(42,26) = pd(42,26) + rrt(132) * density(01) 
  pd(46,01) = pd(46,01) + rrt(132) * density(26) 
  pd(46,26) = pd(46,26) + rrt(132) * density(01) 
  pd(01,01) = pd(01,01) + rrt(133) * density(26) 
  pd(01,26) = pd(01,26) + rrt(133) * density(01) 
  pd(13,01) = pd(13,01) + rrt(133) * density(26) 
  pd(13,26) = pd(13,26) + rrt(133) * density(01) 
  pd(26,01) = pd(26,01) - rrt(133) * density(26) 
  pd(26,26) = pd(26,26) - rrt(133) * density(01) 
  pd(30,01) = pd(30,01) + rrt(133) * density(26) 
  pd(30,26) = pd(30,26) + rrt(133) * density(01) 
  pd(01,01) = pd(01,01) + rrt(134) * density(17) 
  pd(01,17) = pd(01,17) + rrt(134) * density(01) 
  pd(17,01) = pd(17,01) - rrt(134) * density(17) 
  pd(17,17) = pd(17,17) - rrt(134) * density(01) 
  pd(20,01) = pd(20,01) + rrt(134) * density(17) 
  pd(20,17) = pd(20,17) + rrt(134) * density(01) 
  pd(01,01) = pd(01,01) + rrt(135) * density(17) 
  pd(01,17) = pd(01,17) + rrt(135) * density(01) 
  pd(04,01) = pd(04,01) + rrt(135) * density(17) 
  pd(04,17) = pd(04,17) + rrt(135) * density(01) 
  pd(17,01) = pd(17,01) - rrt(135) * density(17) 
  pd(17,17) = pd(17,17) - rrt(135) * density(01) 
  pd(43,01) = pd(43,01) + rrt(135) * density(17) 
  pd(43,17) = pd(43,17) + rrt(135) * density(01) 
  pd(01,01) = pd(01,01) + rrt(136) * density(17) 
  pd(01,17) = pd(01,17) + rrt(136) * density(01) 
  pd(09,01) = pd(09,01) + rrt(136) * density(17) 
  pd(09,17) = pd(09,17) + rrt(136) * density(01) 
  pd(17,01) = pd(17,01) - rrt(136) * density(17) 
  pd(17,17) = pd(17,17) - rrt(136) * density(01) 
  pd(47,01) = pd(47,01) + rrt(136) * density(17) 
  pd(47,17) = pd(47,17) + rrt(136) * density(01) 
  pd(01,01) = pd(01,01) + rrt(137) * density(17) 
  pd(01,17) = pd(01,17) + rrt(137) * density(01) 
  pd(17,01) = pd(17,01) - rrt(137) * density(17) 
  pd(17,17) = pd(17,17) - rrt(137) * density(01) 
  pd(32,01) = pd(32,01) + rrt(137) * density(17) 
  pd(32,17) = pd(32,17) + rrt(137) * density(01) 
  pd(46,01) = pd(46,01) + rrt(137) * density(17) 
  pd(46,17) = pd(46,17) + rrt(137) * density(01) 
  pd(01,01) = pd(01,01) + rrt(138) * density(17) 
  pd(01,17) = pd(01,17) + rrt(138) * density(01) 
  pd(13,01) = pd(13,01) + rrt(138) * density(17) 
  pd(13,17) = pd(13,17) + rrt(138) * density(01) 
  pd(17,01) = pd(17,01) - rrt(138) * density(17) 
  pd(17,17) = pd(17,17) - rrt(138) * density(01) 
  pd(42,01) = pd(42,01) + rrt(138) * density(17) 
  pd(42,17) = pd(42,17) + rrt(138) * density(01) 
  pd(11,01) = pd(11,01) + rrt(139) * density(19) 
  pd(11,19) = pd(11,19) + rrt(139) * density(01) 
  pd(19,01) = pd(19,01) - rrt(139) * density(19) 
  pd(19,19) = pd(19,19) - rrt(139) * density(01) 
  pd(34,01) = pd(34,01) + rrt(139) * density(19) 
  pd(34,19) = pd(34,19) + rrt(139) * density(01) 
  pd(11,01) = pd(11,01) + rrt(140) * density(51) 
  pd(11,51) = pd(11,51) + rrt(140) * density(01) 
  pd(34,01) = pd(34,01) + rrt(140) * density(51) 
  pd(34,51) = pd(34,51) + rrt(140) * density(01) 
  pd(51,01) = pd(51,01) - rrt(140) * density(51) 
  pd(51,51) = pd(51,51) - rrt(140) * density(01) 
  pd(19,01) = pd(19,01) - rrt(141) * density(19) 
  pd(19,19) = pd(19,19) - rrt(141) * density(01) 
  pd(23,01) = pd(23,01) + rrt(141) * density(19) 
  pd(23,19) = pd(23,19) + rrt(141) * density(01) 
  pd(47,01) = pd(47,01) + rrt(141) * density(19) 
  pd(47,19) = pd(47,19) + rrt(141) * density(01) 
  pd(23,01) = pd(23,01) + rrt(142) * density(51) 
  pd(23,51) = pd(23,51) + rrt(142) * density(01) 
  pd(47,01) = pd(47,01) + rrt(142) * density(51) 
  pd(47,51) = pd(47,51) + rrt(142) * density(01) 
  pd(51,01) = pd(51,01) - rrt(142) * density(51) 
  pd(51,51) = pd(51,51) - rrt(142) * density(01) 
  pd(19,01) = pd(19,01) - rrt(143) * density(19) 
  pd(19,19) = pd(19,19) - rrt(143) * density(01) 
  pd(23,01) = pd(23,01) + rrt(143) * density(19) 
  pd(23,19) = pd(23,19) + rrt(143) * density(01) 
  pd(34,01) = pd(34,01) + rrt(143) * density(19) 
  pd(34,19) = pd(34,19) + rrt(143) * density(01) 
  pd(43,01) = pd(43,01) + rrt(143) * density(19) 
  pd(43,19) = pd(43,19) + rrt(143) * density(01) 
  pd(23,01) = pd(23,01) + rrt(144) * density(51) 
  pd(23,51) = pd(23,51) + rrt(144) * density(01) 
  pd(34,01) = pd(34,01) + rrt(144) * density(51) 
  pd(34,51) = pd(34,51) + rrt(144) * density(01) 
  pd(43,01) = pd(43,01) + rrt(144) * density(51) 
  pd(43,51) = pd(43,51) + rrt(144) * density(01) 
  pd(51,01) = pd(51,01) - rrt(144) * density(51) 
  pd(51,51) = pd(51,51) - rrt(144) * density(01) 
  pd(01,01) = pd(01,01) + rrt(145) * density(19) 
  pd(01,19) = pd(01,19) + rrt(145) * density(01) 
  pd(19,01) = pd(19,01) - rrt(145) * density(19) 
  pd(19,19) = pd(19,19) - rrt(145) * density(01) 
  pd(48,01) = pd(48,01) + rrt(145) * density(19) 
  pd(48,19) = pd(48,19) + rrt(145) * density(01) 
  pd(01,01) = pd(01,01) + rrt(146) * density(51) 
  pd(01,51) = pd(01,51) + rrt(146) * density(01) 
  pd(48,01) = pd(48,01) + rrt(146) * density(51) 
  pd(48,51) = pd(48,51) + rrt(146) * density(01) 
  pd(51,01) = pd(51,01) - rrt(146) * density(51) 
  pd(51,51) = pd(51,51) - rrt(146) * density(01) 
  pd(01,01) = pd(01,01) + rrt(147) * density(19) 
  pd(01,19) = pd(01,19) + rrt(147) * density(01) 
  pd(19,01) = pd(19,01) - rrt(147) * density(19) 
  pd(19,19) = pd(19,19) - rrt(147) * density(01) 
  pd(34,01) = pd(34,01) + rrt(147) * density(19) 
  pd(34,19) = pd(34,19) + rrt(147) * density(01) 
  pd(49,01) = pd(49,01) + rrt(147) * density(19) 
  pd(49,19) = pd(49,19) + rrt(147) * density(01) 
  pd(01,01) = pd(01,01) + rrt(148) * density(51) 
  pd(01,51) = pd(01,51) + rrt(148) * density(01) 
  pd(34,01) = pd(34,01) + rrt(148) * density(51) 
  pd(34,51) = pd(34,51) + rrt(148) * density(01) 
  pd(49,01) = pd(49,01) + rrt(148) * density(51) 
  pd(49,51) = pd(49,51) + rrt(148) * density(01) 
  pd(51,01) = pd(51,01) - rrt(148) * density(51) 
  pd(51,51) = pd(51,51) - rrt(148) * density(01) 
  pd(01,01) = pd(01,01) + rrt(149) * density(19) 
  pd(01,19) = pd(01,19) + rrt(149) * density(01) 
  pd(19,01) = pd(19,01) - rrt(149) * density(19) 
  pd(19,19) = pd(19,19) - rrt(149) * density(01) 
  pd(23,01) = pd(23,01) + rrt(149) * density(19) 
  pd(23,19) = pd(23,19) + rrt(149) * density(01) 
  pd(46,01) = pd(46,01) + rrt(149) * density(19) 
  pd(46,19) = pd(46,19) + rrt(149) * density(01) 
  pd(01,01) = pd(01,01) + rrt(150) * density(51) 
  pd(01,51) = pd(01,51) + rrt(150) * density(01) 
  pd(23,01) = pd(23,01) + rrt(150) * density(51) 
  pd(23,51) = pd(23,51) + rrt(150) * density(01) 
  pd(46,01) = pd(46,01) + rrt(150) * density(51) 
  pd(46,51) = pd(46,51) + rrt(150) * density(01) 
  pd(51,01) = pd(51,01) - rrt(150) * density(51) 
  pd(51,51) = pd(51,51) - rrt(150) * density(01) 
  pd(01,01) = pd(01,01) + rrt(151) * density(19) 
  pd(01,19) = pd(01,19) + rrt(151) * density(01) 
  pd(13,01) = pd(13,01) + rrt(151) * density(19) 
  pd(13,19) = pd(13,19) + rrt(151) * density(01) 
  pd(19,01) = pd(19,01) - rrt(151) * density(19) 
  pd(19,19) = pd(19,19) - rrt(151) * density(01) 
  pd(23,01) = pd(23,01) + rrt(151) * density(19) 
  pd(23,19) = pd(23,19) + rrt(151) * density(01) 
  pd(34,01) = pd(34,01) + rrt(151) * density(19) 
  pd(34,19) = pd(34,19) + rrt(151) * density(01) 
  pd(01,01) = pd(01,01) + rrt(152) * density(51) 
  pd(01,51) = pd(01,51) + rrt(152) * density(01) 
  pd(13,01) = pd(13,01) + rrt(152) * density(51) 
  pd(13,51) = pd(13,51) + rrt(152) * density(01) 
  pd(23,01) = pd(23,01) + rrt(152) * density(51) 
  pd(23,51) = pd(23,51) + rrt(152) * density(01) 
  pd(34,01) = pd(34,01) + rrt(152) * density(51) 
  pd(34,51) = pd(34,51) + rrt(152) * density(01) 
  pd(51,01) = pd(51,01) - rrt(152) * density(51) 
  pd(51,51) = pd(51,51) - rrt(152) * density(01) 
  pd(03,01) = pd(03,01) - rrt(153) * density(03) 
  pd(03,03) = pd(03,03) - rrt(153) * density(01) 
  pd(08,01) = pd(08,01) + rrt(153) * density(03) 
  pd(08,03) = pd(08,03) + rrt(153) * density(01) 
  pd(34,01) = pd(34,01) + rrt(153) * density(03) 
  pd(34,03) = pd(34,03) + rrt(153) * density(01) 
  pd(08,01) = pd(08,01) + rrt(154) * density(29) 
  pd(08,29) = pd(08,29) + rrt(154) * density(01) 
  pd(29,01) = pd(29,01) - rrt(154) * density(29) 
  pd(29,29) = pd(29,29) - rrt(154) * density(01) 
  pd(34,01) = pd(34,01) + rrt(154) * density(29) 
  pd(34,29) = pd(34,29) + rrt(154) * density(01) 
  pd(03,01) = pd(03,01) - rrt(155) * density(03) 
  pd(03,03) = pd(03,03) - rrt(155) * density(01) 
  pd(23,01) = pd(23,01) + rrt(155) * density(03) 
  pd(23,03) = pd(23,03) + rrt(155) * density(01) 
  pd(30,01) = pd(30,01) + rrt(155) * density(03) 
  pd(30,03) = pd(30,03) + rrt(155) * density(01) 
  pd(23,01) = pd(23,01) + rrt(156) * density(29) 
  pd(23,29) = pd(23,29) + rrt(156) * density(01) 
  pd(29,01) = pd(29,01) - rrt(156) * density(29) 
  pd(29,29) = pd(29,29) - rrt(156) * density(01) 
  pd(30,01) = pd(30,01) + rrt(156) * density(29) 
  pd(30,29) = pd(30,29) + rrt(156) * density(01) 
  pd(03,01) = pd(03,01) - rrt(157) * density(03) 
  pd(03,03) = pd(03,03) - rrt(157) * density(01) 
  pd(23,01) = pd(23,01) + rrt(157) * density(03) 
  pd(23,03) = pd(23,03) + rrt(157) * density(01) 
  pd(34,01) = pd(34,01) + rrt(157) * density(03) 
  pd(34,03) = pd(34,03) + rrt(157) * density(01) 
  pd(42,01) = pd(42,01) + rrt(157) * density(03) 
  pd(42,03) = pd(42,03) + rrt(157) * density(01) 
  pd(23,01) = pd(23,01) + rrt(158) * density(29) 
  pd(23,29) = pd(23,29) + rrt(158) * density(01) 
  pd(29,01) = pd(29,01) - rrt(158) * density(29) 
  pd(29,29) = pd(29,29) - rrt(158) * density(01) 
  pd(34,01) = pd(34,01) + rrt(158) * density(29) 
  pd(34,29) = pd(34,29) + rrt(158) * density(01) 
  pd(42,01) = pd(42,01) + rrt(158) * density(29) 
  pd(42,29) = pd(42,29) + rrt(158) * density(01) 
  pd(03,01) = pd(03,01) - rrt(159) * density(03) 
  pd(03,03) = pd(03,03) - rrt(159) * density(01) 
  pd(23,01) = pd(23,01) + rrt(159) * density(03) * 2.0d0
  pd(23,03) = pd(23,03) + rrt(159) * density(01) * 2.0d0
  pd(32,01) = pd(32,01) + rrt(159) * density(03) 
  pd(32,03) = pd(32,03) + rrt(159) * density(01) 
  pd(23,01) = pd(23,01) + rrt(160) * density(29) * 2.0d0
  pd(23,29) = pd(23,29) + rrt(160) * density(01) * 2.0d0
  pd(29,01) = pd(29,01) - rrt(160) * density(29) 
  pd(29,29) = pd(29,29) - rrt(160) * density(01) 
  pd(32,01) = pd(32,01) + rrt(160) * density(29) 
  pd(32,29) = pd(32,29) + rrt(160) * density(01) 
  pd(03,01) = pd(03,01) - rrt(161) * density(03) 
  pd(03,03) = pd(03,03) - rrt(161) * density(01) 
  pd(06,01) = pd(06,01) + rrt(161) * density(03) 
  pd(06,03) = pd(06,03) + rrt(161) * density(01) 
  pd(47,01) = pd(47,01) + rrt(161) * density(03) 
  pd(47,03) = pd(47,03) + rrt(161) * density(01) 
  pd(06,01) = pd(06,01) + rrt(162) * density(29) 
  pd(06,29) = pd(06,29) + rrt(162) * density(01) 
  pd(29,01) = pd(29,01) - rrt(162) * density(29) 
  pd(29,29) = pd(29,29) - rrt(162) * density(01) 
  pd(47,01) = pd(47,01) + rrt(162) * density(29) 
  pd(47,29) = pd(47,29) + rrt(162) * density(01) 
  pd(03,01) = pd(03,01) - rrt(163) * density(03) 
  pd(03,03) = pd(03,03) - rrt(163) * density(01) 
  pd(11,01) = pd(11,01) + rrt(163) * density(03) * 2.0d0
  pd(11,03) = pd(11,03) + rrt(163) * density(01) * 2.0d0
  pd(11,01) = pd(11,01) + rrt(164) * density(29) * 2.0d0
  pd(11,29) = pd(11,29) + rrt(164) * density(01) * 2.0d0
  pd(29,01) = pd(29,01) - rrt(164) * density(29) 
  pd(29,29) = pd(29,29) - rrt(164) * density(01) 
  pd(05,01) = pd(05,01) - rrt(165) * density(05) 
  pd(05,05) = pd(05,05) - rrt(165) * density(01) 
  pd(34,01) = pd(34,01) + rrt(165) * density(05) 
  pd(34,05) = pd(34,05) + rrt(165) * density(01) 
  pd(42,01) = pd(42,01) + rrt(165) * density(05) 
  pd(42,05) = pd(42,05) + rrt(165) * density(01) 
  pd(15,01) = pd(15,01) - rrt(166) * density(15) 
  pd(15,15) = pd(15,15) - rrt(166) * density(01) 
  pd(34,01) = pd(34,01) + rrt(166) * density(15) 
  pd(34,15) = pd(34,15) + rrt(166) * density(01) 
  pd(42,01) = pd(42,01) + rrt(166) * density(15) 
  pd(42,15) = pd(42,15) + rrt(166) * density(01) 
  pd(05,01) = pd(05,01) - rrt(167) * density(05) 
  pd(05,05) = pd(05,05) - rrt(167) * density(01) 
  pd(23,01) = pd(23,01) + rrt(167) * density(05) 
  pd(23,05) = pd(23,05) + rrt(167) * density(01) 
  pd(32,01) = pd(32,01) + rrt(167) * density(05) 
  pd(32,05) = pd(32,05) + rrt(167) * density(01) 
  pd(15,01) = pd(15,01) - rrt(168) * density(15) 
  pd(15,15) = pd(15,15) - rrt(168) * density(01) 
  pd(23,01) = pd(23,01) + rrt(168) * density(15) 
  pd(23,15) = pd(23,15) + rrt(168) * density(01) 
  pd(32,01) = pd(32,01) + rrt(168) * density(15) 
  pd(32,15) = pd(32,15) + rrt(168) * density(01) 
  pd(05,01) = pd(05,01) - rrt(169) * density(05) 
  pd(05,05) = pd(05,05) - rrt(169) * density(01) 
  pd(32,01) = pd(32,01) + rrt(169) * density(05) 
  pd(32,05) = pd(32,05) + rrt(169) * density(01) 
  pd(34,01) = pd(34,01) + rrt(169) * density(05) * 2.0d0
  pd(34,05) = pd(34,05) + rrt(169) * density(01) * 2.0d0
  pd(15,01) = pd(15,01) - rrt(170) * density(15) 
  pd(15,15) = pd(15,15) - rrt(170) * density(01) 
  pd(32,01) = pd(32,01) + rrt(170) * density(15) 
  pd(32,15) = pd(32,15) + rrt(170) * density(01) 
  pd(34,01) = pd(34,01) + rrt(170) * density(15) * 2.0d0
  pd(34,15) = pd(34,15) + rrt(170) * density(01) * 2.0d0
  pd(05,01) = pd(05,01) - rrt(171) * density(05) 
  pd(05,05) = pd(05,05) - rrt(171) * density(01) 
  pd(11,01) = pd(11,01) + rrt(171) * density(05) 
  pd(11,05) = pd(11,05) + rrt(171) * density(01) 
  pd(43,01) = pd(43,01) + rrt(171) * density(05) 
  pd(43,05) = pd(43,05) + rrt(171) * density(01) 
  pd(11,01) = pd(11,01) + rrt(172) * density(15) 
  pd(11,15) = pd(11,15) + rrt(172) * density(01) 
  pd(15,01) = pd(15,01) - rrt(172) * density(15) 
  pd(15,15) = pd(15,15) - rrt(172) * density(01) 
  pd(43,01) = pd(43,01) + rrt(172) * density(15) 
  pd(43,15) = pd(43,15) + rrt(172) * density(01) 
  pd(05,01) = pd(05,01) - rrt(173) * density(05) 
  pd(05,05) = pd(05,05) - rrt(173) * density(01) 
  pd(47,01) = pd(47,01) + rrt(173) * density(05) * 2.0d0
  pd(47,05) = pd(47,05) + rrt(173) * density(01) * 2.0d0
  pd(15,01) = pd(15,01) - rrt(174) * density(15) 
  pd(15,15) = pd(15,15) - rrt(174) * density(01) 
  pd(47,01) = pd(47,01) + rrt(174) * density(15) * 2.0d0
  pd(47,15) = pd(47,15) + rrt(174) * density(01) * 2.0d0
  pd(43,01) = pd(43,01) + rrt(175) * density(50) * 2.0d0
  pd(43,50) = pd(43,50) + rrt(175) * density(01) * 2.0d0
  pd(50,01) = pd(50,01) - rrt(175) * density(50) 
  pd(50,50) = pd(50,50) - rrt(175) * density(01) 
  pd(37,01) = pd(37,01) - rrt(176) * density(37) 
  pd(37,37) = pd(37,37) - rrt(176) * density(01) 
  pd(43,01) = pd(43,01) + rrt(176) * density(37) * 2.0d0
  pd(43,37) = pd(43,37) + rrt(176) * density(01) * 2.0d0
  pd(40,01) = pd(40,01) - rrt(177) * density(40) 
  pd(40,40) = pd(40,40) - rrt(177) * density(01) 
  pd(43,01) = pd(43,01) + rrt(177) * density(40) * 2.0d0
  pd(43,40) = pd(43,40) + rrt(177) * density(01) * 2.0d0
  pd(24,01) = pd(24,01) - rrt(178) * density(24) 
  pd(24,24) = pd(24,24) - rrt(178) * density(01) 
  pd(34,01) = pd(34,01) + rrt(178) * density(24) 
  pd(34,24) = pd(34,24) + rrt(178) * density(01) 
  pd(44,01) = pd(44,01) + rrt(178) * density(24) 
  pd(44,24) = pd(44,24) + rrt(178) * density(01) 
  pd(02,01) = pd(02,01) - rrt(179) * density(02) 
  pd(02,02) = pd(02,02) - rrt(179) * density(01) 
  pd(34,01) = pd(34,01) + rrt(179) * density(02) 
  pd(34,02) = pd(34,02) + rrt(179) * density(01) 
  pd(44,01) = pd(44,01) + rrt(179) * density(02) 
  pd(44,02) = pd(44,02) + rrt(179) * density(01) 
  pd(21,01) = pd(21,01) + rrt(180) * density(24) 
  pd(21,24) = pd(21,24) + rrt(180) * density(01) 
  pd(23,01) = pd(23,01) + rrt(180) * density(24) 
  pd(23,24) = pd(23,24) + rrt(180) * density(01) 
  pd(24,01) = pd(24,01) - rrt(180) * density(24) 
  pd(24,24) = pd(24,24) - rrt(180) * density(01) 
  pd(02,01) = pd(02,01) - rrt(181) * density(02) 
  pd(02,02) = pd(02,02) - rrt(181) * density(01) 
  pd(21,01) = pd(21,01) + rrt(181) * density(02) 
  pd(21,02) = pd(21,02) + rrt(181) * density(01) 
  pd(23,01) = pd(23,01) + rrt(181) * density(02) 
  pd(23,02) = pd(23,02) + rrt(181) * density(01) 
  pd(17,01) = pd(17,01) + rrt(182) * density(24) 
  pd(17,24) = pd(17,24) + rrt(182) * density(01) 
  pd(23,01) = pd(23,01) + rrt(182) * density(24) * 2.0d0
  pd(23,24) = pd(23,24) + rrt(182) * density(01) * 2.0d0
  pd(24,01) = pd(24,01) - rrt(182) * density(24) 
  pd(24,24) = pd(24,24) - rrt(182) * density(01) 
  pd(02,01) = pd(02,01) - rrt(183) * density(02) 
  pd(02,02) = pd(02,02) - rrt(183) * density(01) 
  pd(17,01) = pd(17,01) + rrt(183) * density(02) 
  pd(17,02) = pd(17,02) + rrt(183) * density(01) 
  pd(23,01) = pd(23,01) + rrt(183) * density(02) * 2.0d0
  pd(23,02) = pd(23,02) + rrt(183) * density(01) * 2.0d0
  pd(14,01) = pd(14,01) + rrt(184) * density(24) 
  pd(14,24) = pd(14,24) + rrt(184) * density(01) 
  pd(24,01) = pd(24,01) - rrt(184) * density(24) 
  pd(24,24) = pd(24,24) - rrt(184) * density(01) 
  pd(47,01) = pd(47,01) + rrt(184) * density(24) 
  pd(47,24) = pd(47,24) + rrt(184) * density(01) 
  pd(02,01) = pd(02,01) - rrt(185) * density(02) 
  pd(02,02) = pd(02,02) - rrt(185) * density(01) 
  pd(14,01) = pd(14,01) + rrt(185) * density(02) 
  pd(14,02) = pd(14,02) + rrt(185) * density(01) 
  pd(47,01) = pd(47,01) + rrt(185) * density(02) 
  pd(47,02) = pd(47,02) + rrt(185) * density(01) 
  pd(08,01) = pd(08,01) + rrt(186) * density(24) 
  pd(08,24) = pd(08,24) + rrt(186) * density(01) 
  pd(11,01) = pd(11,01) + rrt(186) * density(24) 
  pd(11,24) = pd(11,24) + rrt(186) * density(01) 
  pd(24,01) = pd(24,01) - rrt(186) * density(24) 
  pd(24,24) = pd(24,24) - rrt(186) * density(01) 
  pd(02,01) = pd(02,01) - rrt(187) * density(02) 
  pd(02,02) = pd(02,02) - rrt(187) * density(01) 
  pd(08,01) = pd(08,01) + rrt(187) * density(02) 
  pd(08,02) = pd(08,02) + rrt(187) * density(01) 
  pd(11,01) = pd(11,01) + rrt(187) * density(02) 
  pd(11,02) = pd(11,02) + rrt(187) * density(01) 
  pd(06,01) = pd(06,01) + rrt(188) * density(24) 
  pd(06,24) = pd(06,24) + rrt(188) * density(01) 
  pd(24,01) = pd(24,01) - rrt(188) * density(24) 
  pd(24,24) = pd(24,24) - rrt(188) * density(01) 
  pd(30,01) = pd(30,01) + rrt(188) * density(24) 
  pd(30,24) = pd(30,24) + rrt(188) * density(01) 
  pd(02,01) = pd(02,01) - rrt(189) * density(02) 
  pd(02,02) = pd(02,02) - rrt(189) * density(01) 
  pd(06,01) = pd(06,01) + rrt(189) * density(02) 
  pd(06,02) = pd(06,02) + rrt(189) * density(01) 
  pd(30,01) = pd(30,01) + rrt(189) * density(02) 
  pd(30,02) = pd(30,02) + rrt(189) * density(01) 
  pd(01,01) = pd(01,01) + rrt(190) * density(03) 
  pd(01,03) = pd(01,03) + rrt(190) * density(01) 
  pd(03,01) = pd(03,01) - rrt(190) * density(03) 
  pd(03,03) = pd(03,03) - rrt(190) * density(01) 
  pd(52,01) = pd(52,01) + rrt(190) * density(03) 
  pd(52,03) = pd(52,03) + rrt(190) * density(01) 
  pd(01,01) = pd(01,01) + rrt(191) * density(29) 
  pd(01,29) = pd(01,29) + rrt(191) * density(01) 
  pd(29,01) = pd(29,01) - rrt(191) * density(29) 
  pd(29,29) = pd(29,29) - rrt(191) * density(01) 
  pd(52,01) = pd(52,01) + rrt(191) * density(29) 
  pd(52,29) = pd(52,29) + rrt(191) * density(01) 
  pd(01,01) = pd(01,01) + rrt(192) * density(03) 
  pd(01,03) = pd(01,03) + rrt(192) * density(01) 
  pd(03,01) = pd(03,01) - rrt(192) * density(03) 
  pd(03,03) = pd(03,03) - rrt(192) * density(01) 
  pd(27,01) = pd(27,01) + rrt(192) * density(03) 
  pd(27,03) = pd(27,03) + rrt(192) * density(01) 
  pd(34,01) = pd(34,01) + rrt(192) * density(03) 
  pd(34,03) = pd(34,03) + rrt(192) * density(01) 
  pd(01,01) = pd(01,01) + rrt(193) * density(29) 
  pd(01,29) = pd(01,29) + rrt(193) * density(01) 
  pd(27,01) = pd(27,01) + rrt(193) * density(29) 
  pd(27,29) = pd(27,29) + rrt(193) * density(01) 
  pd(29,01) = pd(29,01) - rrt(193) * density(29) 
  pd(29,29) = pd(29,29) - rrt(193) * density(01) 
  pd(34,01) = pd(34,01) + rrt(193) * density(29) 
  pd(34,29) = pd(34,29) + rrt(193) * density(01) 
  pd(01,01) = pd(01,01) + rrt(194) * density(03) 
  pd(01,03) = pd(01,03) + rrt(194) * density(01) 
  pd(03,01) = pd(03,01) - rrt(194) * density(03) 
  pd(03,03) = pd(03,03) - rrt(194) * density(01) 
  pd(16,01) = pd(16,01) + rrt(194) * density(03) 
  pd(16,03) = pd(16,03) + rrt(194) * density(01) 
  pd(23,01) = pd(23,01) + rrt(194) * density(03) 
  pd(23,03) = pd(23,03) + rrt(194) * density(01) 
  pd(01,01) = pd(01,01) + rrt(195) * density(29) 
  pd(01,29) = pd(01,29) + rrt(195) * density(01) 
  pd(16,01) = pd(16,01) + rrt(195) * density(29) 
  pd(16,29) = pd(16,29) + rrt(195) * density(01) 
  pd(23,01) = pd(23,01) + rrt(195) * density(29) 
  pd(23,29) = pd(23,29) + rrt(195) * density(01) 
  pd(29,01) = pd(29,01) - rrt(195) * density(29) 
  pd(29,29) = pd(29,29) - rrt(195) * density(01) 
  pd(01,01) = pd(01,01) + rrt(196) * density(03) 
  pd(01,03) = pd(01,03) + rrt(196) * density(01) 
  pd(03,01) = pd(03,01) - rrt(196) * density(03) 
  pd(03,03) = pd(03,03) - rrt(196) * density(01) 
  pd(04,01) = pd(04,01) + rrt(196) * density(03) 
  pd(04,03) = pd(04,03) + rrt(196) * density(01) 
  pd(23,01) = pd(23,01) + rrt(196) * density(03) 
  pd(23,03) = pd(23,03) + rrt(196) * density(01) 
  pd(34,01) = pd(34,01) + rrt(196) * density(03) 
  pd(34,03) = pd(34,03) + rrt(196) * density(01) 
  pd(01,01) = pd(01,01) + rrt(197) * density(29) 
  pd(01,29) = pd(01,29) + rrt(197) * density(01) 
  pd(04,01) = pd(04,01) + rrt(197) * density(29) 
  pd(04,29) = pd(04,29) + rrt(197) * density(01) 
  pd(23,01) = pd(23,01) + rrt(197) * density(29) 
  pd(23,29) = pd(23,29) + rrt(197) * density(01) 
  pd(29,01) = pd(29,01) - rrt(197) * density(29) 
  pd(29,29) = pd(29,29) - rrt(197) * density(01) 
  pd(34,01) = pd(34,01) + rrt(197) * density(29) 
  pd(34,29) = pd(34,29) + rrt(197) * density(01) 
  pd(01,01) = pd(01,01) + rrt(198) * density(03) 
  pd(01,03) = pd(01,03) + rrt(198) * density(01) 
  pd(03,01) = pd(03,01) - rrt(198) * density(03) 
  pd(03,03) = pd(03,03) - rrt(198) * density(01) 
  pd(09,01) = pd(09,01) + rrt(198) * density(03) 
  pd(09,03) = pd(09,03) + rrt(198) * density(01) 
  pd(23,01) = pd(23,01) + rrt(198) * density(03) * 2.0d0
  pd(23,03) = pd(23,03) + rrt(198) * density(01) * 2.0d0
  pd(01,01) = pd(01,01) + rrt(199) * density(29) 
  pd(01,29) = pd(01,29) + rrt(199) * density(01) 
  pd(09,01) = pd(09,01) + rrt(199) * density(29) 
  pd(09,29) = pd(09,29) + rrt(199) * density(01) 
  pd(23,01) = pd(23,01) + rrt(199) * density(29) * 2.0d0
  pd(23,29) = pd(23,29) + rrt(199) * density(01) * 2.0d0
  pd(29,01) = pd(29,01) - rrt(199) * density(29) 
  pd(29,29) = pd(29,29) - rrt(199) * density(01) 
  pd(01,01) = pd(01,01) + rrt(200) * density(03) 
  pd(01,03) = pd(01,03) + rrt(200) * density(01) 
  pd(03,01) = pd(03,01) - rrt(200) * density(03) 
  pd(03,03) = pd(03,03) - rrt(200) * density(01) 
  pd(11,01) = pd(11,01) + rrt(200) * density(03) 
  pd(11,03) = pd(11,03) + rrt(200) * density(01) 
  pd(49,01) = pd(49,01) + rrt(200) * density(03) 
  pd(49,03) = pd(49,03) + rrt(200) * density(01) 
  pd(01,01) = pd(01,01) + rrt(201) * density(29) 
  pd(01,29) = pd(01,29) + rrt(201) * density(01) 
  pd(11,01) = pd(11,01) + rrt(201) * density(29) 
  pd(11,29) = pd(11,29) + rrt(201) * density(01) 
  pd(29,01) = pd(29,01) - rrt(201) * density(29) 
  pd(29,29) = pd(29,29) - rrt(201) * density(01) 
  pd(49,01) = pd(49,01) + rrt(201) * density(29) 
  pd(49,29) = pd(49,29) + rrt(201) * density(01) 
  pd(01,01) = pd(01,01) + rrt(202) * density(03) 
  pd(01,03) = pd(01,03) + rrt(202) * density(01) 
  pd(03,01) = pd(03,01) - rrt(202) * density(03) 
  pd(03,03) = pd(03,03) - rrt(202) * density(01) 
  pd(06,01) = pd(06,01) + rrt(202) * density(03) 
  pd(06,03) = pd(06,03) + rrt(202) * density(01) 
  pd(46,01) = pd(46,01) + rrt(202) * density(03) 
  pd(46,03) = pd(46,03) + rrt(202) * density(01) 
  pd(01,01) = pd(01,01) + rrt(203) * density(29) 
  pd(01,29) = pd(01,29) + rrt(203) * density(01) 
  pd(06,01) = pd(06,01) + rrt(203) * density(29) 
  pd(06,29) = pd(06,29) + rrt(203) * density(01) 
  pd(29,01) = pd(29,01) - rrt(203) * density(29) 
  pd(29,29) = pd(29,29) - rrt(203) * density(01) 
  pd(46,01) = pd(46,01) + rrt(203) * density(29) 
  pd(46,29) = pd(46,29) + rrt(203) * density(01) 
  pd(01,01) = pd(01,01) + rrt(204) * density(05) 
  pd(01,05) = pd(01,05) + rrt(204) * density(01) 
  pd(05,01) = pd(05,01) - rrt(204) * density(05) 
  pd(05,05) = pd(05,05) - rrt(204) * density(01) 
  pd(16,01) = pd(16,01) + rrt(204) * density(05) 
  pd(16,05) = pd(16,05) + rrt(204) * density(01) 
  pd(01,01) = pd(01,01) + rrt(205) * density(15) 
  pd(01,15) = pd(01,15) + rrt(205) * density(01) 
  pd(15,01) = pd(15,01) - rrt(205) * density(15) 
  pd(15,15) = pd(15,15) - rrt(205) * density(01) 
  pd(16,01) = pd(16,01) + rrt(205) * density(15) 
  pd(16,15) = pd(16,15) + rrt(205) * density(01) 
  pd(01,01) = pd(01,01) + rrt(206) * density(05) 
  pd(01,05) = pd(01,05) + rrt(206) * density(01) 
  pd(04,01) = pd(04,01) + rrt(206) * density(05) 
  pd(04,05) = pd(04,05) + rrt(206) * density(01) 
  pd(05,01) = pd(05,01) - rrt(206) * density(05) 
  pd(05,05) = pd(05,05) - rrt(206) * density(01) 
  pd(34,01) = pd(34,01) + rrt(206) * density(05) 
  pd(34,05) = pd(34,05) + rrt(206) * density(01) 
  pd(01,01) = pd(01,01) + rrt(207) * density(15) 
  pd(01,15) = pd(01,15) + rrt(207) * density(01) 
  pd(04,01) = pd(04,01) + rrt(207) * density(15) 
  pd(04,15) = pd(04,15) + rrt(207) * density(01) 
  pd(15,01) = pd(15,01) - rrt(207) * density(15) 
  pd(15,15) = pd(15,15) - rrt(207) * density(01) 
  pd(34,01) = pd(34,01) + rrt(207) * density(15) 
  pd(34,15) = pd(34,15) + rrt(207) * density(01) 
  pd(01,01) = pd(01,01) + rrt(208) * density(05) 
  pd(01,05) = pd(01,05) + rrt(208) * density(01) 
  pd(05,01) = pd(05,01) - rrt(208) * density(05) 
  pd(05,05) = pd(05,05) - rrt(208) * density(01) 
  pd(43,01) = pd(43,01) + rrt(208) * density(05) 
  pd(43,05) = pd(43,05) + rrt(208) * density(01) 
  pd(49,01) = pd(49,01) + rrt(208) * density(05) 
  pd(49,05) = pd(49,05) + rrt(208) * density(01) 
  pd(01,01) = pd(01,01) + rrt(209) * density(15) 
  pd(01,15) = pd(01,15) + rrt(209) * density(01) 
  pd(15,01) = pd(15,01) - rrt(209) * density(15) 
  pd(15,15) = pd(15,15) - rrt(209) * density(01) 
  pd(43,01) = pd(43,01) + rrt(209) * density(15) 
  pd(43,15) = pd(43,15) + rrt(209) * density(01) 
  pd(49,01) = pd(49,01) + rrt(209) * density(15) 
  pd(49,15) = pd(49,15) + rrt(209) * density(01) 
  pd(01,01) = pd(01,01) + rrt(210) * density(05) 
  pd(01,05) = pd(01,05) + rrt(210) * density(01) 
  pd(05,01) = pd(05,01) - rrt(210) * density(05) 
  pd(05,05) = pd(05,05) - rrt(210) * density(01) 
  pd(46,01) = pd(46,01) + rrt(210) * density(05) 
  pd(46,05) = pd(46,05) + rrt(210) * density(01) 
  pd(47,01) = pd(47,01) + rrt(210) * density(05) 
  pd(47,05) = pd(47,05) + rrt(210) * density(01) 
  pd(01,01) = pd(01,01) + rrt(211) * density(15) 
  pd(01,15) = pd(01,15) + rrt(211) * density(01) 
  pd(15,01) = pd(15,01) - rrt(211) * density(15) 
  pd(15,15) = pd(15,15) - rrt(211) * density(01) 
  pd(46,01) = pd(46,01) + rrt(211) * density(15) 
  pd(46,15) = pd(46,15) + rrt(211) * density(01) 
  pd(47,01) = pd(47,01) + rrt(211) * density(15) 
  pd(47,15) = pd(47,15) + rrt(211) * density(01) 
  pd(01,01) = pd(01,01) + rrt(212) * density(05) 
  pd(01,05) = pd(01,05) + rrt(212) * density(01) 
  pd(05,01) = pd(05,01) - rrt(212) * density(05) 
  pd(05,05) = pd(05,05) - rrt(212) * density(01) 
  pd(11,01) = pd(11,01) + rrt(212) * density(05) 
  pd(11,05) = pd(11,05) + rrt(212) * density(01) 
  pd(13,01) = pd(13,01) + rrt(212) * density(05) 
  pd(13,05) = pd(13,05) + rrt(212) * density(01) 
  pd(01,01) = pd(01,01) + rrt(213) * density(15) 
  pd(01,15) = pd(01,15) + rrt(213) * density(01) 
  pd(11,01) = pd(11,01) + rrt(213) * density(15) 
  pd(11,15) = pd(11,15) + rrt(213) * density(01) 
  pd(13,01) = pd(13,01) + rrt(213) * density(15) 
  pd(13,15) = pd(13,15) + rrt(213) * density(01) 
  pd(15,01) = pd(15,01) - rrt(213) * density(15) 
  pd(15,15) = pd(15,15) - rrt(213) * density(01) 
  pd(01,01) = pd(01,01) + rrt(214) * density(50) 
  pd(01,50) = pd(01,50) + rrt(214) * density(01) 
  pd(09,01) = pd(09,01) + rrt(214) * density(50) 
  pd(09,50) = pd(09,50) + rrt(214) * density(01) 
  pd(50,01) = pd(50,01) - rrt(214) * density(50) 
  pd(50,50) = pd(50,50) - rrt(214) * density(01) 
  pd(01,01) = pd(01,01) + rrt(215) * density(37) 
  pd(01,37) = pd(01,37) + rrt(215) * density(01) 
  pd(09,01) = pd(09,01) + rrt(215) * density(37) 
  pd(09,37) = pd(09,37) + rrt(215) * density(01) 
  pd(37,01) = pd(37,01) - rrt(215) * density(37) 
  pd(37,37) = pd(37,37) - rrt(215) * density(01) 
  pd(01,01) = pd(01,01) + rrt(216) * density(40) 
  pd(01,40) = pd(01,40) + rrt(216) * density(01) 
  pd(09,01) = pd(09,01) + rrt(216) * density(40) 
  pd(09,40) = pd(09,40) + rrt(216) * density(01) 
  pd(40,01) = pd(40,01) - rrt(216) * density(40) 
  pd(40,40) = pd(40,40) - rrt(216) * density(01) 
  pd(01,01) = pd(01,01) + rrt(217) * density(50) 
  pd(01,50) = pd(01,50) + rrt(217) * density(01) 
  pd(13,01) = pd(13,01) + rrt(217) * density(50) 
  pd(13,50) = pd(13,50) + rrt(217) * density(01) 
  pd(43,01) = pd(43,01) + rrt(217) * density(50) 
  pd(43,50) = pd(43,50) + rrt(217) * density(01) 
  pd(50,01) = pd(50,01) - rrt(217) * density(50) 
  pd(50,50) = pd(50,50) - rrt(217) * density(01) 
  pd(01,01) = pd(01,01) + rrt(218) * density(37) 
  pd(01,37) = pd(01,37) + rrt(218) * density(01) 
  pd(13,01) = pd(13,01) + rrt(218) * density(37) 
  pd(13,37) = pd(13,37) + rrt(218) * density(01) 
  pd(37,01) = pd(37,01) - rrt(218) * density(37) 
  pd(37,37) = pd(37,37) - rrt(218) * density(01) 
  pd(43,01) = pd(43,01) + rrt(218) * density(37) 
  pd(43,37) = pd(43,37) + rrt(218) * density(01) 
  pd(01,01) = pd(01,01) + rrt(219) * density(40) 
  pd(01,40) = pd(01,40) + rrt(219) * density(01) 
  pd(13,01) = pd(13,01) + rrt(219) * density(40) 
  pd(13,40) = pd(13,40) + rrt(219) * density(01) 
  pd(40,01) = pd(40,01) - rrt(219) * density(40) 
  pd(40,40) = pd(40,40) - rrt(219) * density(01) 
  pd(43,01) = pd(43,01) + rrt(219) * density(40) 
  pd(43,40) = pd(43,40) + rrt(219) * density(01) 
  pd(01,01) = pd(01,01) + rrt(220) * density(24) 
  pd(01,24) = pd(01,24) + rrt(220) * density(01) 
  pd(24,01) = pd(24,01) - rrt(220) * density(24) 
  pd(24,24) = pd(24,24) - rrt(220) * density(01) 
  pd(41,01) = pd(41,01) + rrt(220) * density(24) 
  pd(41,24) = pd(41,24) + rrt(220) * density(01) 
  pd(01,01) = pd(01,01) + rrt(221) * density(02) 
  pd(01,02) = pd(01,02) + rrt(221) * density(01) 
  pd(02,01) = pd(02,01) - rrt(221) * density(02) 
  pd(02,02) = pd(02,02) - rrt(221) * density(01) 
  pd(41,01) = pd(41,01) + rrt(221) * density(02) 
  pd(41,02) = pd(41,02) + rrt(221) * density(01) 
  pd(01,01) = pd(01,01) + rrt(222) * density(24) 
  pd(01,24) = pd(01,24) + rrt(222) * density(01) 
  pd(18,01) = pd(18,01) + rrt(222) * density(24) 
  pd(18,24) = pd(18,24) + rrt(222) * density(01) 
  pd(24,01) = pd(24,01) - rrt(222) * density(24) 
  pd(24,24) = pd(24,24) - rrt(222) * density(01) 
  pd(34,01) = pd(34,01) + rrt(222) * density(24) 
  pd(34,24) = pd(34,24) + rrt(222) * density(01) 
  pd(01,01) = pd(01,01) + rrt(223) * density(02) 
  pd(01,02) = pd(01,02) + rrt(223) * density(01) 
  pd(02,01) = pd(02,01) - rrt(223) * density(02) 
  pd(02,02) = pd(02,02) - rrt(223) * density(01) 
  pd(18,01) = pd(18,01) + rrt(223) * density(02) 
  pd(18,02) = pd(18,02) + rrt(223) * density(01) 
  pd(34,01) = pd(34,01) + rrt(223) * density(02) 
  pd(34,02) = pd(34,02) + rrt(223) * density(01) 
  pd(01,01) = pd(01,01) + rrt(224) * density(24) 
  pd(01,24) = pd(01,24) + rrt(224) * density(01) 
  pd(23,01) = pd(23,01) + rrt(224) * density(24) 
  pd(23,24) = pd(23,24) + rrt(224) * density(01) 
  pd(24,01) = pd(24,01) - rrt(224) * density(24) 
  pd(24,24) = pd(24,24) - rrt(224) * density(01) 
  pd(38,01) = pd(38,01) + rrt(224) * density(24) 
  pd(38,24) = pd(38,24) + rrt(224) * density(01) 
  pd(01,01) = pd(01,01) + rrt(225) * density(02) 
  pd(01,02) = pd(01,02) + rrt(225) * density(01) 
  pd(02,01) = pd(02,01) - rrt(225) * density(02) 
  pd(02,02) = pd(02,02) - rrt(225) * density(01) 
  pd(23,01) = pd(23,01) + rrt(225) * density(02) 
  pd(23,02) = pd(23,02) + rrt(225) * density(01) 
  pd(38,01) = pd(38,01) + rrt(225) * density(02) 
  pd(38,02) = pd(38,02) + rrt(225) * density(01) 
  pd(01,01) = pd(01,01) + rrt(226) * density(24) 
  pd(01,24) = pd(01,24) + rrt(226) * density(01) 
  pd(23,01) = pd(23,01) + rrt(226) * density(24) 
  pd(23,24) = pd(23,24) + rrt(226) * density(01) 
  pd(24,01) = pd(24,01) - rrt(226) * density(24) 
  pd(24,24) = pd(24,24) - rrt(226) * density(01) 
  pd(34,01) = pd(34,01) + rrt(226) * density(24) 
  pd(34,24) = pd(34,24) + rrt(226) * density(01) 
  pd(39,01) = pd(39,01) + rrt(226) * density(24) 
  pd(39,24) = pd(39,24) + rrt(226) * density(01) 
  pd(01,01) = pd(01,01) + rrt(227) * density(02) 
  pd(01,02) = pd(01,02) + rrt(227) * density(01) 
  pd(02,01) = pd(02,01) - rrt(227) * density(02) 
  pd(02,02) = pd(02,02) - rrt(227) * density(01) 
  pd(23,01) = pd(23,01) + rrt(227) * density(02) 
  pd(23,02) = pd(23,02) + rrt(227) * density(01) 
  pd(34,01) = pd(34,01) + rrt(227) * density(02) 
  pd(34,02) = pd(34,02) + rrt(227) * density(01) 
  pd(39,01) = pd(39,01) + rrt(227) * density(02) 
  pd(39,02) = pd(39,02) + rrt(227) * density(01) 
  pd(01,01) = pd(01,01) + rrt(228) * density(24) 
  pd(01,24) = pd(01,24) + rrt(228) * density(01) 
  pd(20,01) = pd(20,01) + rrt(228) * density(24) 
  pd(20,24) = pd(20,24) + rrt(228) * density(01) 
  pd(23,01) = pd(23,01) + rrt(228) * density(24) * 2.0d0
  pd(23,24) = pd(23,24) + rrt(228) * density(01) * 2.0d0
  pd(24,01) = pd(24,01) - rrt(228) * density(24) 
  pd(24,24) = pd(24,24) - rrt(228) * density(01) 
  pd(01,01) = pd(01,01) + rrt(229) * density(02) 
  pd(01,02) = pd(01,02) + rrt(229) * density(01) 
  pd(02,01) = pd(02,01) - rrt(229) * density(02) 
  pd(02,02) = pd(02,02) - rrt(229) * density(01) 
  pd(20,01) = pd(20,01) + rrt(229) * density(02) 
  pd(20,02) = pd(20,02) + rrt(229) * density(01) 
  pd(23,01) = pd(23,01) + rrt(229) * density(02) * 2.0d0
  pd(23,02) = pd(23,02) + rrt(229) * density(01) * 2.0d0
  pd(01,01) = pd(01,01) + rrt(230) * density(24) 
  pd(01,24) = pd(01,24) + rrt(230) * density(01) 
  pd(11,01) = pd(11,01) + rrt(230) * density(24) 
  pd(11,24) = pd(11,24) + rrt(230) * density(01) 
  pd(24,01) = pd(24,01) - rrt(230) * density(24) 
  pd(24,24) = pd(24,24) - rrt(230) * density(01) 
  pd(27,01) = pd(27,01) + rrt(230) * density(24) 
  pd(27,24) = pd(27,24) + rrt(230) * density(01) 
  pd(01,01) = pd(01,01) + rrt(231) * density(02) 
  pd(01,02) = pd(01,02) + rrt(231) * density(01) 
  pd(02,01) = pd(02,01) - rrt(231) * density(02) 
  pd(02,02) = pd(02,02) - rrt(231) * density(01) 
  pd(11,01) = pd(11,01) + rrt(231) * density(02) 
  pd(11,02) = pd(11,02) + rrt(231) * density(01) 
  pd(27,01) = pd(27,01) + rrt(231) * density(02) 
  pd(27,02) = pd(27,02) + rrt(231) * density(01) 
  pd(01,01) = pd(01,01) + rrt(232) * density(24) 
  pd(01,24) = pd(01,24) + rrt(232) * density(01) 
  pd(06,01) = pd(06,01) + rrt(232) * density(24) 
  pd(06,24) = pd(06,24) + rrt(232) * density(01) 
  pd(16,01) = pd(16,01) + rrt(232) * density(24) 
  pd(16,24) = pd(16,24) + rrt(232) * density(01) 
  pd(24,01) = pd(24,01) - rrt(232) * density(24) 
  pd(24,24) = pd(24,24) - rrt(232) * density(01) 
  pd(01,01) = pd(01,01) + rrt(233) * density(02) 
  pd(01,02) = pd(01,02) + rrt(233) * density(01) 
  pd(02,01) = pd(02,01) - rrt(233) * density(02) 
  pd(02,02) = pd(02,02) - rrt(233) * density(01) 
  pd(06,01) = pd(06,01) + rrt(233) * density(02) 
  pd(06,02) = pd(06,02) + rrt(233) * density(01) 
  pd(16,01) = pd(16,01) + rrt(233) * density(02) 
  pd(16,02) = pd(16,02) + rrt(233) * density(01) 
  pd(01,01) = pd(01,01) + rrt(234) * density(24) 
  pd(01,24) = pd(01,24) + rrt(234) * density(01) 
  pd(08,01) = pd(08,01) + rrt(234) * density(24) 
  pd(08,24) = pd(08,24) + rrt(234) * density(01) 
  pd(24,01) = pd(24,01) - rrt(234) * density(24) 
  pd(24,24) = pd(24,24) - rrt(234) * density(01) 
  pd(49,01) = pd(49,01) + rrt(234) * density(24) 
  pd(49,24) = pd(49,24) + rrt(234) * density(01) 
  pd(01,01) = pd(01,01) + rrt(235) * density(02) 
  pd(01,02) = pd(01,02) + rrt(235) * density(01) 
  pd(02,01) = pd(02,01) - rrt(235) * density(02) 
  pd(02,02) = pd(02,02) - rrt(235) * density(01) 
  pd(08,01) = pd(08,01) + rrt(235) * density(02) 
  pd(08,02) = pd(08,02) + rrt(235) * density(01) 
  pd(49,01) = pd(49,01) + rrt(235) * density(02) 
  pd(49,02) = pd(49,02) + rrt(235) * density(01) 
  pd(01,01) = pd(01,01) + rrt(236) * density(24) 
  pd(01,24) = pd(01,24) + rrt(236) * density(01) 
  pd(14,01) = pd(14,01) + rrt(236) * density(24) 
  pd(14,24) = pd(14,24) + rrt(236) * density(01) 
  pd(24,01) = pd(24,01) - rrt(236) * density(24) 
  pd(24,24) = pd(24,24) - rrt(236) * density(01) 
  pd(46,01) = pd(46,01) + rrt(236) * density(24) 
  pd(46,24) = pd(46,24) + rrt(236) * density(01) 
  pd(01,01) = pd(01,01) + rrt(237) * density(02) 
  pd(01,02) = pd(01,02) + rrt(237) * density(01) 
  pd(02,01) = pd(02,01) - rrt(237) * density(02) 
  pd(02,02) = pd(02,02) - rrt(237) * density(01) 
  pd(14,01) = pd(14,01) + rrt(237) * density(02) 
  pd(14,02) = pd(14,02) + rrt(237) * density(01) 
  pd(46,01) = pd(46,01) + rrt(237) * density(02) 
  pd(46,02) = pd(46,02) + rrt(237) * density(01) 
  pd(23,01) = pd(23,01) - rrt(238) * density(23) 
  pd(23,23) = pd(23,23) - rrt(238) * density(01) 
  pd(34,01) = pd(34,01) + rrt(238) * density(23) * 2.0d0
  pd(34,23) = pd(34,23) + rrt(238) * density(01) * 2.0d0
  pd(01,01) = pd(01,01) + rrt(239) * density(23) 
  pd(01,23) = pd(01,23) + rrt(239) * density(01) 
  pd(23,01) = pd(23,01) - rrt(239) * density(23) 
  pd(23,23) = pd(23,23) - rrt(239) * density(01) 
  pd(45,01) = pd(45,01) + rrt(239) * density(23) 
  pd(45,23) = pd(45,23) + rrt(239) * density(01) 
  pd(01,01) = pd(01,01) - rrt(240) * density(01) * density(48) * 2.0d0
  pd(01,48) = pd(01,48) - rrt(240) * density(01)**2 
  pd(06,01) = pd(06,01) + rrt(240) * density(01) * density(48) * 2.0d0
  pd(06,48) = pd(06,48) + rrt(240) * density(01)**2 
  pd(48,01) = pd(48,01) - rrt(240) * density(01) * density(48) * 2.0d0
  pd(48,48) = pd(48,48) - rrt(240) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(241) * density(01) * density(49) * 2.0d0
  pd(01,49) = pd(01,49) - rrt(241) * density(01)**2 
  pd(11,01) = pd(11,01) + rrt(241) * density(01) * density(49) * 2.0d0
  pd(11,49) = pd(11,49) + rrt(241) * density(01)**2 
  pd(49,01) = pd(49,01) - rrt(241) * density(01) * density(49) * 2.0d0
  pd(49,49) = pd(49,49) - rrt(241) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(242) * density(01) * density(46) * 2.0d0
  pd(01,46) = pd(01,46) - rrt(242) * density(01)**2 
  pd(46,01) = pd(46,01) - rrt(242) * density(01) * density(46) * 2.0d0
  pd(46,46) = pd(46,46) - rrt(242) * density(01)**2 
  pd(47,01) = pd(47,01) + rrt(242) * density(01) * density(46) * 2.0d0
  pd(47,46) = pd(47,46) + rrt(242) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(243) * density(01) * density(13) * 2.0d0
  pd(01,13) = pd(01,13) - rrt(243) * density(01)**2 
  pd(13,01) = pd(13,01) - rrt(243) * density(01) * density(13) * 2.0d0
  pd(13,13) = pd(13,13) - rrt(243) * density(01)**2 
  pd(43,01) = pd(43,01) + rrt(243) * density(01) * density(13) * 2.0d0
  pd(43,13) = pd(43,13) + rrt(243) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(244) * density(01) * density(52) * 2.0d0
  pd(01,52) = pd(01,52) - rrt(244) * density(01)**2 
  pd(14,01) = pd(14,01) + rrt(244) * density(01) * density(52) * 2.0d0
  pd(14,52) = pd(14,52) + rrt(244) * density(01)**2 
  pd(52,01) = pd(52,01) - rrt(244) * density(01) * density(52) * 2.0d0
  pd(52,52) = pd(52,52) - rrt(244) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(245) * density(01) * density(27) * 2.0d0
  pd(01,27) = pd(01,27) - rrt(245) * density(01)**2 
  pd(08,01) = pd(08,01) + rrt(245) * density(01) * density(27) * 2.0d0
  pd(08,27) = pd(08,27) + rrt(245) * density(01)**2 
  pd(27,01) = pd(27,01) - rrt(245) * density(01) * density(27) * 2.0d0
  pd(27,27) = pd(27,27) - rrt(245) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(246) * density(01) * density(16) * 2.0d0
  pd(01,16) = pd(01,16) - rrt(246) * density(01)**2 
  pd(16,01) = pd(16,01) - rrt(246) * density(01) * density(16) * 2.0d0
  pd(16,16) = pd(16,16) - rrt(246) * density(01)**2 
  pd(30,01) = pd(30,01) + rrt(246) * density(01) * density(16) * 2.0d0
  pd(30,16) = pd(30,16) + rrt(246) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(247) * density(01) * density(04) * 2.0d0
  pd(01,04) = pd(01,04) - rrt(247) * density(01)**2 
  pd(04,01) = pd(04,01) - rrt(247) * density(01) * density(04) * 2.0d0
  pd(04,04) = pd(04,04) - rrt(247) * density(01)**2 
  pd(42,01) = pd(42,01) + rrt(247) * density(01) * density(04) * 2.0d0
  pd(42,04) = pd(42,04) + rrt(247) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(248) * density(01) * density(09) * 2.0d0
  pd(01,09) = pd(01,09) - rrt(248) * density(01)**2 
  pd(09,01) = pd(09,01) - rrt(248) * density(01) * density(09) * 2.0d0
  pd(09,09) = pd(09,09) - rrt(248) * density(01)**2 
  pd(32,01) = pd(32,01) + rrt(248) * density(01) * density(09) * 2.0d0
  pd(32,09) = pd(32,09) + rrt(248) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(249) * density(01) * density(41) * 2.0d0
  pd(01,41) = pd(01,41) - rrt(249) * density(01)**2 
  pd(07,01) = pd(07,01) + rrt(249) * density(01) * density(41) * 2.0d0
  pd(07,41) = pd(07,41) + rrt(249) * density(01)**2 
  pd(41,01) = pd(41,01) - rrt(249) * density(01) * density(41) * 2.0d0
  pd(41,41) = pd(41,41) - rrt(249) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(250) * density(01) * density(18) * 2.0d0
  pd(01,18) = pd(01,18) - rrt(250) * density(01)**2 
  pd(18,01) = pd(18,01) - rrt(250) * density(01) * density(18) * 2.0d0
  pd(18,18) = pd(18,18) - rrt(250) * density(01)**2 
  pd(44,01) = pd(44,01) + rrt(250) * density(01) * density(18) * 2.0d0
  pd(44,18) = pd(44,18) + rrt(250) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(251) * density(01) * density(38) * 2.0d0
  pd(01,38) = pd(01,38) - rrt(251) * density(01)**2 
  pd(21,01) = pd(21,01) + rrt(251) * density(01) * density(38) * 2.0d0
  pd(21,38) = pd(21,38) + rrt(251) * density(01)**2 
  pd(38,01) = pd(38,01) - rrt(251) * density(01) * density(38) * 2.0d0
  pd(38,38) = pd(38,38) - rrt(251) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(252) * density(01) * density(39) * 2.0d0
  pd(01,39) = pd(01,39) - rrt(252) * density(01)**2 
  pd(26,01) = pd(26,01) + rrt(252) * density(01) * density(39) * 2.0d0
  pd(26,39) = pd(26,39) + rrt(252) * density(01)**2 
  pd(39,01) = pd(39,01) - rrt(252) * density(01) * density(39) * 2.0d0
  pd(39,39) = pd(39,39) - rrt(252) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(253) * density(01) * density(20) * 2.0d0
  pd(01,20) = pd(01,20) - rrt(253) * density(01)**2 
  pd(17,01) = pd(17,01) + rrt(253) * density(01) * density(20) * 2.0d0
  pd(17,20) = pd(17,20) + rrt(253) * density(01)**2 
  pd(20,01) = pd(20,01) - rrt(253) * density(01) * density(20) * 2.0d0
  pd(20,20) = pd(20,20) - rrt(253) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(254) * density(01) * density(12) * 2.0d0
  pd(01,12) = pd(01,12) - rrt(254) * density(01)**2 
  pd(12,01) = pd(12,01) - rrt(254) * density(01) * density(12) * 2.0d0
  pd(12,12) = pd(12,12) - rrt(254) * density(01)**2 
  pd(34,01) = pd(34,01) + rrt(254) * density(01) * density(12) * 2.0d0
  pd(34,12) = pd(34,12) + rrt(254) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(255) * density(01) * density(45) * 2.0d0
  pd(01,45) = pd(01,45) - rrt(255) * density(01)**2 
  pd(23,01) = pd(23,01) + rrt(255) * density(01) * density(45) * 2.0d0
  pd(23,45) = pd(23,45) + rrt(255) * density(01)**2 
  pd(45,01) = pd(45,01) - rrt(255) * density(01) * density(45) * 2.0d0
  pd(45,45) = pd(45,45) - rrt(255) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(256) * density(35) 
  pd(01,35) = pd(01,35) - rrt(256) * density(01) 
  pd(11,01) = pd(11,01) + rrt(256) * density(35) 
  pd(11,35) = pd(11,35) + rrt(256) * density(01) 
  pd(34,01) = pd(34,01) + rrt(256) * density(35) * 2.0d0
  pd(34,35) = pd(34,35) + rrt(256) * density(01) * 2.0d0
  pd(35,01) = pd(35,01) - rrt(256) * density(35) 
  pd(35,35) = pd(35,35) - rrt(256) * density(01) 
  pd(01,01) = pd(01,01) - rrt(257) * density(35) 
  pd(01,35) = pd(01,35) - rrt(257) * density(01) 
  pd(23,01) = pd(23,01) + rrt(257) * density(35) 
  pd(23,35) = pd(23,35) + rrt(257) * density(01) 
  pd(34,01) = pd(34,01) + rrt(257) * density(35) 
  pd(34,35) = pd(34,35) + rrt(257) * density(01) 
  pd(35,01) = pd(35,01) - rrt(257) * density(35) 
  pd(35,35) = pd(35,35) - rrt(257) * density(01) 
  pd(47,01) = pd(47,01) + rrt(257) * density(35) 
  pd(47,35) = pd(47,35) + rrt(257) * density(01) 
  pd(01,01) = pd(01,01) - rrt(258) * density(48) 
  pd(01,48) = pd(01,48) - rrt(258) * density(01) 
  pd(11,01) = pd(11,01) + rrt(258) * density(48) 
  pd(11,48) = pd(11,48) + rrt(258) * density(01) 
  pd(34,01) = pd(34,01) + rrt(258) * density(48) 
  pd(34,48) = pd(34,48) + rrt(258) * density(01) 
  pd(48,01) = pd(48,01) - rrt(258) * density(48) 
  pd(48,48) = pd(48,48) - rrt(258) * density(01) 
  pd(01,01) = pd(01,01) - rrt(259) * density(48) 
  pd(01,48) = pd(01,48) - rrt(259) * density(01) 
  pd(34,01) = pd(34,01) + rrt(259) * density(48) * 2.0d0
  pd(34,48) = pd(34,48) + rrt(259) * density(01) * 2.0d0
  pd(47,01) = pd(47,01) + rrt(259) * density(48) 
  pd(47,48) = pd(47,48) + rrt(259) * density(01) 
  pd(48,01) = pd(48,01) - rrt(259) * density(48) 
  pd(48,48) = pd(48,48) - rrt(259) * density(01) 
  pd(01,01) = pd(01,01) - rrt(260) * density(48) 
  pd(01,48) = pd(01,48) - rrt(260) * density(01) 
  pd(23,01) = pd(23,01) + rrt(260) * density(48) 
  pd(23,48) = pd(23,48) + rrt(260) * density(01) 
  pd(34,01) = pd(34,01) + rrt(260) * density(48) 
  pd(34,48) = pd(34,48) + rrt(260) * density(01) 
  pd(43,01) = pd(43,01) + rrt(260) * density(48) 
  pd(43,48) = pd(43,48) + rrt(260) * density(01) 
  pd(48,01) = pd(48,01) - rrt(260) * density(48) 
  pd(48,48) = pd(48,48) - rrt(260) * density(01) 
  pd(01,01) = pd(01,01) - rrt(261) * density(49) 
  pd(01,49) = pd(01,49) - rrt(261) * density(01) 
  pd(34,01) = pd(34,01) + rrt(261) * density(49) 
  pd(34,49) = pd(34,49) + rrt(261) * density(01) 
  pd(47,01) = pd(47,01) + rrt(261) * density(49) 
  pd(47,49) = pd(47,49) + rrt(261) * density(01) 
  pd(49,01) = pd(49,01) - rrt(261) * density(49) 
  pd(49,49) = pd(49,49) - rrt(261) * density(01) 
  pd(01,01) = pd(01,01) - rrt(262) * density(49) 
  pd(01,49) = pd(01,49) - rrt(262) * density(01) 
  pd(23,01) = pd(23,01) + rrt(262) * density(49) 
  pd(23,49) = pd(23,49) + rrt(262) * density(01) 
  pd(43,01) = pd(43,01) + rrt(262) * density(49) 
  pd(43,49) = pd(43,49) + rrt(262) * density(01) 
  pd(49,01) = pd(49,01) - rrt(262) * density(49) 
  pd(49,49) = pd(49,49) - rrt(262) * density(01) 
  pd(01,01) = pd(01,01) - rrt(263) * density(46) 
  pd(01,46) = pd(01,46) - rrt(263) * density(01) 
  pd(34,01) = pd(34,01) + rrt(263) * density(46) 
  pd(34,46) = pd(34,46) + rrt(263) * density(01) 
  pd(43,01) = pd(43,01) + rrt(263) * density(46) 
  pd(43,46) = pd(43,46) + rrt(263) * density(01) 
  pd(46,01) = pd(46,01) - rrt(263) * density(46) 
  pd(46,46) = pd(46,46) - rrt(263) * density(01) 
  pd(01,01) = pd(01,01) - rrt(264) * density(52) 
  pd(01,52) = pd(01,52) - rrt(264) * density(01) 
  pd(08,01) = pd(08,01) + rrt(264) * density(52) 
  pd(08,52) = pd(08,52) + rrt(264) * density(01) 
  pd(34,01) = pd(34,01) + rrt(264) * density(52) 
  pd(34,52) = pd(34,52) + rrt(264) * density(01) 
  pd(52,01) = pd(52,01) - rrt(264) * density(52) 
  pd(52,52) = pd(52,52) - rrt(264) * density(01) 
  pd(01,01) = pd(01,01) - rrt(265) * density(52) 
  pd(01,52) = pd(01,52) - rrt(265) * density(01) 
  pd(30,01) = pd(30,01) + rrt(265) * density(52) 
  pd(30,52) = pd(30,52) + rrt(265) * density(01) 
  pd(34,01) = pd(34,01) + rrt(265) * density(52) * 2.0d0
  pd(34,52) = pd(34,52) + rrt(265) * density(01) * 2.0d0
  pd(52,01) = pd(52,01) - rrt(265) * density(52) 
  pd(52,52) = pd(52,52) - rrt(265) * density(01) 
  pd(01,01) = pd(01,01) - rrt(266) * density(27) 
  pd(01,27) = pd(01,27) - rrt(266) * density(01) 
  pd(27,01) = pd(27,01) - rrt(266) * density(27) 
  pd(27,27) = pd(27,27) - rrt(266) * density(01) 
  pd(30,01) = pd(30,01) + rrt(266) * density(27) 
  pd(30,27) = pd(30,27) + rrt(266) * density(01) 
  pd(34,01) = pd(34,01) + rrt(266) * density(27) 
  pd(34,27) = pd(34,27) + rrt(266) * density(01) 
  pd(01,01) = pd(01,01) - rrt(267) * density(27) 
  pd(01,27) = pd(01,27) - rrt(267) * density(01) 
  pd(27,01) = pd(27,01) - rrt(267) * density(27) 
  pd(27,27) = pd(27,27) - rrt(267) * density(01) 
  pd(34,01) = pd(34,01) + rrt(267) * density(27) * 2.0d0
  pd(34,27) = pd(34,27) + rrt(267) * density(01) * 2.0d0
  pd(42,01) = pd(42,01) + rrt(267) * density(27) 
  pd(42,27) = pd(42,27) + rrt(267) * density(01) 
  pd(01,01) = pd(01,01) - rrt(268) * density(27) 
  pd(01,27) = pd(01,27) - rrt(268) * density(01) 
  pd(23,01) = pd(23,01) + rrt(268) * density(27) 
  pd(23,27) = pd(23,27) + rrt(268) * density(01) 
  pd(27,01) = pd(27,01) - rrt(268) * density(27) 
  pd(27,27) = pd(27,27) - rrt(268) * density(01) 
  pd(32,01) = pd(32,01) + rrt(268) * density(27) 
  pd(32,27) = pd(32,27) + rrt(268) * density(01) 
  pd(34,01) = pd(34,01) + rrt(268) * density(27) 
  pd(34,27) = pd(34,27) + rrt(268) * density(01) 
  pd(01,01) = pd(01,01) - rrt(269) * density(27) 
  pd(01,27) = pd(01,27) - rrt(269) * density(01) 
  pd(27,01) = pd(27,01) - rrt(269) * density(27) 
  pd(27,27) = pd(27,27) - rrt(269) * density(01) 
  pd(32,01) = pd(32,01) + rrt(269) * density(27) 
  pd(32,27) = pd(32,27) + rrt(269) * density(01) 
  pd(34,01) = pd(34,01) + rrt(269) * density(27) * 3.0d0
  pd(34,27) = pd(34,27) + rrt(269) * density(01) * 3.0d0
  pd(01,01) = pd(01,01) - rrt(270) * density(27) 
  pd(01,27) = pd(01,27) - rrt(270) * density(01) 
  pd(11,01) = pd(11,01) + rrt(270) * density(27) 
  pd(11,27) = pd(11,27) + rrt(270) * density(01) 
  pd(27,01) = pd(27,01) - rrt(270) * density(27) 
  pd(27,27) = pd(27,27) - rrt(270) * density(01) 
  pd(47,01) = pd(47,01) + rrt(270) * density(27) 
  pd(47,27) = pd(47,27) + rrt(270) * density(01) 
  pd(01,01) = pd(01,01) - rrt(271) * density(16) 
  pd(01,16) = pd(01,16) - rrt(271) * density(01) 
  pd(16,01) = pd(16,01) - rrt(271) * density(16) 
  pd(16,16) = pd(16,16) - rrt(271) * density(01) 
  pd(34,01) = pd(34,01) + rrt(271) * density(16) 
  pd(34,16) = pd(34,16) + rrt(271) * density(01) 
  pd(42,01) = pd(42,01) + rrt(271) * density(16) 
  pd(42,16) = pd(42,16) + rrt(271) * density(01) 
  pd(01,01) = pd(01,01) - rrt(272) * density(16) 
  pd(01,16) = pd(01,16) - rrt(272) * density(01) 
  pd(16,01) = pd(16,01) - rrt(272) * density(16) 
  pd(16,16) = pd(16,16) - rrt(272) * density(01) 
  pd(32,01) = pd(32,01) + rrt(272) * density(16) 
  pd(32,16) = pd(32,16) + rrt(272) * density(01) 
  pd(34,01) = pd(34,01) + rrt(272) * density(16) * 2.0d0
  pd(34,16) = pd(34,16) + rrt(272) * density(01) * 2.0d0
  pd(01,01) = pd(01,01) - rrt(273) * density(04) 
  pd(01,04) = pd(01,04) - rrt(273) * density(01) 
  pd(04,01) = pd(04,01) - rrt(273) * density(04) 
  pd(04,04) = pd(04,04) - rrt(273) * density(01) 
  pd(32,01) = pd(32,01) + rrt(273) * density(04) 
  pd(32,04) = pd(32,04) + rrt(273) * density(01) 
  pd(34,01) = pd(34,01) + rrt(273) * density(04) 
  pd(34,04) = pd(34,04) + rrt(273) * density(01) 
  pd(01,01) = pd(01,01) - rrt(274) * density(09) 
  pd(01,09) = pd(01,09) - rrt(274) * density(01) 
  pd(09,01) = pd(09,01) - rrt(274) * density(09) 
  pd(09,09) = pd(09,09) - rrt(274) * density(01) 
  pd(43,01) = pd(43,01) + rrt(274) * density(09) * 2.0d0
  pd(43,09) = pd(43,09) + rrt(274) * density(01) * 2.0d0
  pd(01,01) = pd(01,01) - rrt(275) * density(45) 
  pd(01,45) = pd(01,45) - rrt(275) * density(01) 
  pd(34,01) = pd(34,01) + rrt(275) * density(45) * 2.0d0
  pd(34,45) = pd(34,45) + rrt(275) * density(01) * 2.0d0
  pd(45,01) = pd(45,01) - rrt(275) * density(45) 
  pd(45,45) = pd(45,45) - rrt(275) * density(01) 
  pd(01,01) = pd(01,01) - rrt(276) * density(28) 
  pd(01,28) = pd(01,28) - rrt(276) * density(01) 
  pd(23,01) = pd(23,01) + rrt(276) * density(28) 
  pd(23,28) = pd(23,28) + rrt(276) * density(01) 
  pd(28,01) = pd(28,01) - rrt(276) * density(28) 
  pd(28,28) = pd(28,28) - rrt(276) * density(01) 
  pd(34,01) = pd(34,01) + rrt(276) * density(28) 
  pd(34,28) = pd(34,28) + rrt(276) * density(01) 
  pd(23,23) = pd(23,23) - rrt(277) * density(45) 
  pd(23,45) = pd(23,45) - rrt(277) * density(23) 
  pd(28,23) = pd(28,23) + rrt(277) * density(45) 
  pd(28,45) = pd(28,45) + rrt(277) * density(23) 
  pd(34,23) = pd(34,23) + rrt(277) * density(45) 
  pd(34,45) = pd(34,45) + rrt(277) * density(23) 
  pd(45,23) = pd(45,23) - rrt(277) * density(45) 
  pd(45,45) = pd(45,45) - rrt(277) * density(23) 
  pd(06,35) = pd(06,35) + rrt(278) * density(47) 
  pd(06,47) = pd(06,47) + rrt(278) * density(35) 
  pd(35,35) = pd(35,35) - rrt(278) * density(47) 
  pd(35,47) = pd(35,47) - rrt(278) * density(35) 
  pd(47,35) = pd(47,35) - rrt(278) * density(47) 
  pd(47,47) = pd(47,47) - rrt(278) * density(35) 
  pd(49,35) = pd(49,35) + rrt(278) * density(47) 
  pd(49,47) = pd(49,47) + rrt(278) * density(35) 
  pd(06,35) = pd(06,35) + rrt(279) * density(43) 
  pd(06,43) = pd(06,43) + rrt(279) * density(35) 
  pd(35,35) = pd(35,35) - rrt(279) * density(43) 
  pd(35,43) = pd(35,43) - rrt(279) * density(35) 
  pd(43,35) = pd(43,35) - rrt(279) * density(43) 
  pd(43,43) = pd(43,43) - rrt(279) * density(35) 
  pd(46,35) = pd(46,35) + rrt(279) * density(43) 
  pd(46,43) = pd(46,43) + rrt(279) * density(35) 
  pd(06,14) = pd(06,14) + rrt(280) * density(35) 
  pd(06,35) = pd(06,35) + rrt(280) * density(14) 
  pd(14,14) = pd(14,14) - rrt(280) * density(35) 
  pd(14,35) = pd(14,35) - rrt(280) * density(14) 
  pd(23,14) = pd(23,14) + rrt(280) * density(35) 
  pd(23,35) = pd(23,35) + rrt(280) * density(14) 
  pd(27,14) = pd(27,14) + rrt(280) * density(35) 
  pd(27,35) = pd(27,35) + rrt(280) * density(14) 
  pd(35,14) = pd(35,14) - rrt(280) * density(35) 
  pd(35,35) = pd(35,35) - rrt(280) * density(14) 
  pd(06,30) = pd(06,30) + rrt(281) * density(35) 
  pd(06,35) = pd(06,35) + rrt(281) * density(30) 
  pd(27,30) = pd(27,30) + rrt(281) * density(35) 
  pd(27,35) = pd(27,35) + rrt(281) * density(30) 
  pd(30,30) = pd(30,30) - rrt(281) * density(35) 
  pd(30,35) = pd(30,35) - rrt(281) * density(30) 
  pd(35,30) = pd(35,30) - rrt(281) * density(35) 
  pd(35,35) = pd(35,35) - rrt(281) * density(30) 
  pd(04,32) = pd(04,32) + rrt(282) * density(35) 
  pd(04,35) = pd(04,35) + rrt(282) * density(32) 
  pd(06,32) = pd(06,32) + rrt(282) * density(35) 
  pd(06,35) = pd(06,35) + rrt(282) * density(32) 
  pd(32,32) = pd(32,32) - rrt(282) * density(35) 
  pd(32,35) = pd(32,35) - rrt(282) * density(32) 
  pd(35,32) = pd(35,32) - rrt(282) * density(35) 
  pd(35,35) = pd(35,35) - rrt(282) * density(32) 
  pd(23,34) = pd(23,34) + rrt(283) * density(35) 
  pd(23,35) = pd(23,35) + rrt(283) * density(34) 
  pd(34,34) = pd(34,34) - rrt(283) * density(35) 
  pd(34,35) = pd(34,35) - rrt(283) * density(34) 
  pd(35,34) = pd(35,34) - rrt(283) * density(35) 
  pd(35,35) = pd(35,35) - rrt(283) * density(34) 
  pd(48,34) = pd(48,34) + rrt(283) * density(35) 
  pd(48,35) = pd(48,35) + rrt(283) * density(34) 
  pd(06,06) = pd(06,06) - rrt(284) * density(48) 
  pd(06,48) = pd(06,48) - rrt(284) * density(06) 
  pd(11,06) = pd(11,06) + rrt(284) * density(48) 
  pd(11,48) = pd(11,48) + rrt(284) * density(06) 
  pd(35,06) = pd(35,06) + rrt(284) * density(48) 
  pd(35,48) = pd(35,48) + rrt(284) * density(06) 
  pd(48,06) = pd(48,06) - rrt(284) * density(48) 
  pd(48,48) = pd(48,48) - rrt(284) * density(06) 
  pd(06,14) = pd(06,14) + rrt(285) * density(48) 
  pd(06,48) = pd(06,48) + rrt(285) * density(14) 
  pd(14,14) = pd(14,14) - rrt(285) * density(48) 
  pd(14,48) = pd(14,48) - rrt(285) * density(14) 
  pd(16,14) = pd(16,14) + rrt(285) * density(48) 
  pd(16,48) = pd(16,48) + rrt(285) * density(14) 
  pd(23,14) = pd(23,14) + rrt(285) * density(48) 
  pd(23,48) = pd(23,48) + rrt(285) * density(14) 
  pd(48,14) = pd(48,14) - rrt(285) * density(48) 
  pd(48,48) = pd(48,48) - rrt(285) * density(14) 
  pd(11,30) = pd(11,30) + rrt(286) * density(48) 
  pd(11,48) = pd(11,48) + rrt(286) * density(30) 
  pd(27,30) = pd(27,30) + rrt(286) * density(48) 
  pd(27,48) = pd(27,48) + rrt(286) * density(30) 
  pd(30,30) = pd(30,30) - rrt(286) * density(48) 
  pd(30,48) = pd(30,48) - rrt(286) * density(30) 
  pd(48,30) = pd(48,30) - rrt(286) * density(48) 
  pd(48,48) = pd(48,48) - rrt(286) * density(30) 
  pd(06,30) = pd(06,30) + rrt(287) * density(48) 
  pd(06,48) = pd(06,48) + rrt(287) * density(30) 
  pd(16,30) = pd(16,30) + rrt(287) * density(48) 
  pd(16,48) = pd(16,48) + rrt(287) * density(30) 
  pd(30,30) = pd(30,30) - rrt(287) * density(48) 
  pd(30,48) = pd(30,48) - rrt(287) * density(30) 
  pd(48,30) = pd(48,30) - rrt(287) * density(48) 
  pd(48,48) = pd(48,48) - rrt(287) * density(30) 
  pd(04,32) = pd(04,32) + rrt(288) * density(48) 
  pd(04,48) = pd(04,48) + rrt(288) * density(32) 
  pd(11,32) = pd(11,32) + rrt(288) * density(48) 
  pd(11,48) = pd(11,48) + rrt(288) * density(32) 
  pd(32,32) = pd(32,32) - rrt(288) * density(48) 
  pd(32,48) = pd(32,48) - rrt(288) * density(32) 
  pd(48,32) = pd(48,32) - rrt(288) * density(48) 
  pd(48,48) = pd(48,48) - rrt(288) * density(32) 
  pd(06,32) = pd(06,32) + rrt(289) * density(48) 
  pd(06,48) = pd(06,48) + rrt(289) * density(32) 
  pd(09,32) = pd(09,32) + rrt(289) * density(48) 
  pd(09,48) = pd(09,48) + rrt(289) * density(32) 
  pd(32,32) = pd(32,32) - rrt(289) * density(48) 
  pd(32,48) = pd(32,48) - rrt(289) * density(32) 
  pd(48,32) = pd(48,32) - rrt(289) * density(48) 
  pd(48,48) = pd(48,48) - rrt(289) * density(32) 
  pd(23,23) = pd(23,23) - rrt(290) * density(48) 
  pd(23,48) = pd(23,48) - rrt(290) * density(23) 
  pd(34,23) = pd(34,23) + rrt(290) * density(48) 
  pd(34,48) = pd(34,48) + rrt(290) * density(23) 
  pd(35,23) = pd(35,23) + rrt(290) * density(48) 
  pd(35,48) = pd(35,48) + rrt(290) * density(23) 
  pd(48,23) = pd(48,23) - rrt(290) * density(48) 
  pd(48,48) = pd(48,48) - rrt(290) * density(23) 
  pd(23,34) = pd(23,34) + rrt(291) * density(48) 
  pd(23,48) = pd(23,48) + rrt(291) * density(34) 
  pd(34,34) = pd(34,34) - rrt(291) * density(48) 
  pd(34,48) = pd(34,48) - rrt(291) * density(34) 
  pd(48,34) = pd(48,34) - rrt(291) * density(48) 
  pd(48,48) = pd(48,48) - rrt(291) * density(34) 
  pd(49,34) = pd(49,34) + rrt(291) * density(48) 
  pd(49,48) = pd(49,48) + rrt(291) * density(34) 
  pd(06,06) = pd(06,06) - rrt(292) * density(49) 
  pd(06,49) = pd(06,49) - rrt(292) * density(06) 
  pd(11,06) = pd(11,06) + rrt(292) * density(49) 
  pd(11,49) = pd(11,49) + rrt(292) * density(06) 
  pd(48,06) = pd(48,06) + rrt(292) * density(49) 
  pd(48,49) = pd(48,49) + rrt(292) * density(06) 
  pd(49,06) = pd(49,06) - rrt(292) * density(49) 
  pd(49,49) = pd(49,49) - rrt(292) * density(06) 
  pd(06,06) = pd(06,06) - rrt(293) * density(49) 
  pd(06,49) = pd(06,49) - rrt(293) * density(06) 
  pd(23,06) = pd(23,06) + rrt(293) * density(49) 
  pd(23,49) = pd(23,49) + rrt(293) * density(06) 
  pd(27,06) = pd(27,06) + rrt(293) * density(49) 
  pd(27,49) = pd(27,49) + rrt(293) * density(06) 
  pd(49,06) = pd(49,06) - rrt(293) * density(49) 
  pd(49,49) = pd(49,49) - rrt(293) * density(06) 
  pd(04,47) = pd(04,47) + rrt(294) * density(49) 
  pd(04,49) = pd(04,49) + rrt(294) * density(47) 
  pd(23,47) = pd(23,47) + rrt(294) * density(49) 
  pd(23,49) = pd(23,49) + rrt(294) * density(47) 
  pd(47,47) = pd(47,47) - rrt(294) * density(49) 
  pd(47,49) = pd(47,49) - rrt(294) * density(47) 
  pd(49,47) = pd(49,47) - rrt(294) * density(49) 
  pd(49,49) = pd(49,49) - rrt(294) * density(47) 
  pd(09,43) = pd(09,43) + rrt(295) * density(49) 
  pd(09,49) = pd(09,49) + rrt(295) * density(43) 
  pd(23,43) = pd(23,43) + rrt(295) * density(49) 
  pd(23,49) = pd(23,49) + rrt(295) * density(43) 
  pd(43,43) = pd(43,43) - rrt(295) * density(49) 
  pd(43,49) = pd(43,49) - rrt(295) * density(43) 
  pd(49,43) = pd(49,43) - rrt(295) * density(49) 
  pd(49,49) = pd(49,49) - rrt(295) * density(43) 
  pd(06,14) = pd(06,14) + rrt(296) * density(49) 
  pd(06,49) = pd(06,49) + rrt(296) * density(14) 
  pd(14,14) = pd(14,14) - rrt(296) * density(49) 
  pd(14,49) = pd(14,49) - rrt(296) * density(14) 
  pd(27,14) = pd(27,14) + rrt(296) * density(49) 
  pd(27,49) = pd(27,49) + rrt(296) * density(14) 
  pd(49,14) = pd(49,14) - rrt(296) * density(49) 
  pd(49,49) = pd(49,49) - rrt(296) * density(14) 
  pd(04,30) = pd(04,30) + rrt(297) * density(49) 
  pd(04,49) = pd(04,49) + rrt(297) * density(30) 
  pd(06,30) = pd(06,30) + rrt(297) * density(49) 
  pd(06,49) = pd(06,49) + rrt(297) * density(30) 
  pd(30,30) = pd(30,30) - rrt(297) * density(49) 
  pd(30,49) = pd(30,49) - rrt(297) * density(30) 
  pd(49,30) = pd(49,30) - rrt(297) * density(49) 
  pd(49,49) = pd(49,49) - rrt(297) * density(30) 
  pd(04,42) = pd(04,42) + rrt(298) * density(49) 
  pd(04,49) = pd(04,49) + rrt(298) * density(42) 
  pd(11,42) = pd(11,42) + rrt(298) * density(49) 
  pd(11,49) = pd(11,49) + rrt(298) * density(42) 
  pd(42,42) = pd(42,42) - rrt(298) * density(49) 
  pd(42,49) = pd(42,49) - rrt(298) * density(42) 
  pd(49,42) = pd(49,42) - rrt(298) * density(49) 
  pd(49,49) = pd(49,49) - rrt(298) * density(42) 
  pd(06,06) = pd(06,06) - rrt(299) * density(46) 
  pd(06,46) = pd(06,46) - rrt(299) * density(06) 
  pd(11,06) = pd(11,06) + rrt(299) * density(46) 
  pd(11,46) = pd(11,46) + rrt(299) * density(06) 
  pd(46,06) = pd(46,06) - rrt(299) * density(46) 
  pd(46,46) = pd(46,46) - rrt(299) * density(06) 
  pd(49,06) = pd(49,06) + rrt(299) * density(46) 
  pd(49,46) = pd(49,46) + rrt(299) * density(06) 
  pd(06,06) = pd(06,06) - rrt(300) * density(46) 
  pd(06,46) = pd(06,46) - rrt(300) * density(06) 
  pd(27,06) = pd(27,06) + rrt(300) * density(46) 
  pd(27,46) = pd(27,46) + rrt(300) * density(06) 
  pd(34,06) = pd(34,06) + rrt(300) * density(46) 
  pd(34,46) = pd(34,46) + rrt(300) * density(06) 
  pd(46,06) = pd(46,06) - rrt(300) * density(46) 
  pd(46,46) = pd(46,46) - rrt(300) * density(06) 
  pd(06,06) = pd(06,06) - rrt(301) * density(46) 
  pd(06,46) = pd(06,46) - rrt(301) * density(06) 
  pd(16,06) = pd(16,06) + rrt(301) * density(46) 
  pd(16,46) = pd(16,46) + rrt(301) * density(06) 
  pd(23,06) = pd(23,06) + rrt(301) * density(46) 
  pd(23,46) = pd(23,46) + rrt(301) * density(06) 
  pd(46,06) = pd(46,06) - rrt(301) * density(46) 
  pd(46,46) = pd(46,46) - rrt(301) * density(06) 
  pd(04,06) = pd(04,06) + rrt(302) * density(46) 
  pd(04,46) = pd(04,46) + rrt(302) * density(06) 
  pd(06,06) = pd(06,06) - rrt(302) * density(46) 
  pd(06,46) = pd(06,46) - rrt(302) * density(06) 
  pd(23,06) = pd(23,06) + rrt(302) * density(46) 
  pd(23,46) = pd(23,46) + rrt(302) * density(06) 
  pd(34,06) = pd(34,06) + rrt(302) * density(46) 
  pd(34,46) = pd(34,46) + rrt(302) * density(06) 
  pd(46,06) = pd(46,06) - rrt(302) * density(46) 
  pd(46,46) = pd(46,46) - rrt(302) * density(06) 
  pd(06,06) = pd(06,06) - rrt(303) * density(46) 
  pd(06,46) = pd(06,46) - rrt(303) * density(06) 
  pd(09,06) = pd(09,06) + rrt(303) * density(46) 
  pd(09,46) = pd(09,46) + rrt(303) * density(06) 
  pd(23,06) = pd(23,06) + rrt(303) * density(46) * 2.0d0
  pd(23,46) = pd(23,46) + rrt(303) * density(06) * 2.0d0
  pd(46,06) = pd(46,06) - rrt(303) * density(46) 
  pd(46,46) = pd(46,46) - rrt(303) * density(06) 
  pd(23,23) = pd(23,23) - rrt(304) * density(46) 
  pd(23,46) = pd(23,46) - rrt(304) * density(23) 
  pd(34,23) = pd(34,23) + rrt(304) * density(46) 
  pd(34,46) = pd(34,46) + rrt(304) * density(23) 
  pd(46,23) = pd(46,23) - rrt(304) * density(46) 
  pd(46,46) = pd(46,46) - rrt(304) * density(23) 
  pd(49,23) = pd(49,23) + rrt(304) * density(46) 
  pd(49,46) = pd(49,46) + rrt(304) * density(23) 
  pd(06,06) = pd(06,06) - rrt(305) * density(13) 
  pd(06,13) = pd(06,13) - rrt(305) * density(06) 
  pd(13,06) = pd(13,06) - rrt(305) * density(13) 
  pd(13,13) = pd(13,13) - rrt(305) * density(06) 
  pd(16,06) = pd(16,06) + rrt(305) * density(13) 
  pd(16,13) = pd(16,13) + rrt(305) * density(06) 
  pd(34,06) = pd(34,06) + rrt(305) * density(13) 
  pd(34,13) = pd(34,13) + rrt(305) * density(06) 
  pd(04,06) = pd(04,06) + rrt(306) * density(13) 
  pd(04,13) = pd(04,13) + rrt(306) * density(06) 
  pd(06,06) = pd(06,06) - rrt(306) * density(13) 
  pd(06,13) = pd(06,13) - rrt(306) * density(06) 
  pd(13,06) = pd(13,06) - rrt(306) * density(13) 
  pd(13,13) = pd(13,13) - rrt(306) * density(06) 
  pd(23,06) = pd(23,06) + rrt(306) * density(13) 
  pd(23,13) = pd(23,13) + rrt(306) * density(06) 
  pd(06,06) = pd(06,06) - rrt(307) * density(13) 
  pd(06,13) = pd(06,13) - rrt(307) * density(06) 
  pd(09,06) = pd(09,06) + rrt(307) * density(13) 
  pd(09,13) = pd(09,13) + rrt(307) * density(06) 
  pd(13,06) = pd(13,06) - rrt(307) * density(13) 
  pd(13,13) = pd(13,13) - rrt(307) * density(06) 
  pd(23,06) = pd(23,06) + rrt(307) * density(13) 
  pd(23,13) = pd(23,13) + rrt(307) * density(06) 
  pd(34,06) = pd(34,06) + rrt(307) * density(13) 
  pd(34,13) = pd(34,13) + rrt(307) * density(06) 
  pd(13,13) = pd(13,13) - rrt(308) * density(23) 
  pd(13,23) = pd(13,23) - rrt(308) * density(13) 
  pd(23,13) = pd(23,13) - rrt(308) * density(23) 
  pd(23,23) = pd(23,23) - rrt(308) * density(13) 
  pd(34,13) = pd(34,13) + rrt(308) * density(23) 
  pd(34,23) = pd(34,23) + rrt(308) * density(13) 
  pd(46,13) = pd(46,13) + rrt(308) * density(23) 
  pd(46,23) = pd(46,23) + rrt(308) * density(13) 
  pd(14,30) = pd(14,30) + rrt(309) * density(52) 
  pd(14,52) = pd(14,52) + rrt(309) * density(30) 
  pd(16,30) = pd(16,30) + rrt(309) * density(52) 
  pd(16,52) = pd(16,52) + rrt(309) * density(30) 
  pd(30,30) = pd(30,30) - rrt(309) * density(52) 
  pd(30,52) = pd(30,52) - rrt(309) * density(30) 
  pd(52,30) = pd(52,30) - rrt(309) * density(52) 
  pd(52,52) = pd(52,52) - rrt(309) * density(30) 
  pd(27,32) = pd(27,32) + rrt(310) * density(52) 
  pd(27,52) = pd(27,52) + rrt(310) * density(32) 
  pd(32,32) = pd(32,32) - rrt(310) * density(52) 
  pd(32,52) = pd(32,52) - rrt(310) * density(32) 
  pd(42,32) = pd(42,32) + rrt(310) * density(52) 
  pd(42,52) = pd(42,52) + rrt(310) * density(32) 
  pd(52,32) = pd(52,32) - rrt(310) * density(52) 
  pd(52,52) = pd(52,52) - rrt(310) * density(32) 
  pd(23,34) = pd(23,34) + rrt(311) * density(52) 
  pd(23,52) = pd(23,52) + rrt(311) * density(34) 
  pd(27,34) = pd(27,34) + rrt(311) * density(52) 
  pd(27,52) = pd(27,52) + rrt(311) * density(34) 
  pd(34,34) = pd(34,34) - rrt(311) * density(52) 
  pd(34,52) = pd(34,52) - rrt(311) * density(34) 
  pd(52,34) = pd(52,34) - rrt(311) * density(52) 
  pd(52,52) = pd(52,52) - rrt(311) * density(34) 
  pd(16,27) = pd(16,27) + rrt(312) * density(34) 
  pd(16,34) = pd(16,34) + rrt(312) * density(27) 
  pd(23,27) = pd(23,27) + rrt(312) * density(34) 
  pd(23,34) = pd(23,34) + rrt(312) * density(27) 
  pd(27,27) = pd(27,27) - rrt(312) * density(34) 
  pd(27,34) = pd(27,34) - rrt(312) * density(27) 
  pd(34,27) = pd(34,27) - rrt(312) * density(34) 
  pd(34,34) = pd(34,34) - rrt(312) * density(27) 
  pd(16,16) = pd(16,16) - rrt(313) * density(42) 
  pd(16,42) = pd(16,42) - rrt(313) * density(16) 
  pd(27,16) = pd(27,16) + rrt(313) * density(42) 
  pd(27,42) = pd(27,42) + rrt(313) * density(16) 
  pd(32,16) = pd(32,16) + rrt(313) * density(42) 
  pd(32,42) = pd(32,42) + rrt(313) * density(16) 
  pd(42,16) = pd(42,16) - rrt(313) * density(42) 
  pd(42,42) = pd(42,42) - rrt(313) * density(16) 
  pd(04,16) = pd(04,16) + rrt(314) * density(42) 
  pd(04,42) = pd(04,42) + rrt(314) * density(16) 
  pd(16,16) = pd(16,16) - rrt(314) * density(42) 
  pd(16,42) = pd(16,42) - rrt(314) * density(16) 
  pd(30,16) = pd(30,16) + rrt(314) * density(42) 
  pd(30,42) = pd(30,42) + rrt(314) * density(16) 
  pd(42,16) = pd(42,16) - rrt(314) * density(42) 
  pd(42,42) = pd(42,42) - rrt(314) * density(16) 
  pd(04,16) = pd(04,16) + rrt(315) * density(34) 
  pd(04,34) = pd(04,34) + rrt(315) * density(16) 
  pd(16,16) = pd(16,16) - rrt(315) * density(34) 
  pd(16,34) = pd(16,34) - rrt(315) * density(16) 
  pd(23,16) = pd(23,16) + rrt(315) * density(34) 
  pd(23,34) = pd(23,34) + rrt(315) * density(16) 
  pd(34,16) = pd(34,16) - rrt(315) * density(34) 
  pd(34,34) = pd(34,34) - rrt(315) * density(16) 
  pd(04,04) = pd(04,04) - rrt(316) * density(14) 
  pd(04,14) = pd(04,14) - rrt(316) * density(04) 
  pd(14,04) = pd(14,04) - rrt(316) * density(14) 
  pd(14,14) = pd(14,14) - rrt(316) * density(04) 
  pd(27,04) = pd(27,04) + rrt(316) * density(14) 
  pd(27,14) = pd(27,14) + rrt(316) * density(04) 
  pd(30,04) = pd(30,04) + rrt(316) * density(14) 
  pd(30,14) = pd(30,14) + rrt(316) * density(04) 
  pd(04,04) = pd(04,04) - rrt(317) * density(30) 
  pd(04,30) = pd(04,30) - rrt(317) * density(04) 
  pd(27,04) = pd(27,04) + rrt(317) * density(30) 
  pd(27,30) = pd(27,30) + rrt(317) * density(04) 
  pd(30,04) = pd(30,04) - rrt(317) * density(30) 
  pd(30,30) = pd(30,30) - rrt(317) * density(04) 
  pd(32,04) = pd(32,04) + rrt(317) * density(30) 
  pd(32,30) = pd(32,30) + rrt(317) * density(04) 
  pd(04,04) = pd(04,04) - rrt(318) * density(33) 
  pd(04,33) = pd(04,33) - rrt(318) * density(04) 
  pd(09,04) = pd(09,04) + rrt(318) * density(33) 
  pd(09,33) = pd(09,33) + rrt(318) * density(04) 
  pd(32,04) = pd(32,04) + rrt(318) * density(33) 
  pd(32,33) = pd(32,33) + rrt(318) * density(04) 
  pd(33,04) = pd(33,04) - rrt(318) * density(33) 
  pd(33,33) = pd(33,33) - rrt(318) * density(04) 
  pd(04,04) = pd(04,04) - rrt(319) * density(34) 
  pd(04,34) = pd(04,34) - rrt(319) * density(04) 
  pd(09,04) = pd(09,04) + rrt(319) * density(34) 
  pd(09,34) = pd(09,34) + rrt(319) * density(04) 
  pd(23,04) = pd(23,04) + rrt(319) * density(34) 
  pd(23,34) = pd(23,34) + rrt(319) * density(04) 
  pd(34,04) = pd(34,04) - rrt(319) * density(34) 
  pd(34,34) = pd(34,34) - rrt(319) * density(04) 
  pd(04,06) = pd(04,06) + rrt(320) * density(09) 
  pd(04,09) = pd(04,09) + rrt(320) * density(06) 
  pd(06,06) = pd(06,06) - rrt(320) * density(09) 
  pd(06,09) = pd(06,09) - rrt(320) * density(06) 
  pd(09,06) = pd(09,06) - rrt(320) * density(09) 
  pd(09,09) = pd(09,09) - rrt(320) * density(06) 
  pd(11,06) = pd(11,06) + rrt(320) * density(09) 
  pd(11,09) = pd(11,09) + rrt(320) * density(06) 
  pd(09,09) = pd(09,09) - rrt(321) * density(14) 
  pd(09,14) = pd(09,14) - rrt(321) * density(09) 
  pd(14,09) = pd(14,09) - rrt(321) * density(14) 
  pd(14,14) = pd(14,14) - rrt(321) * density(09) 
  pd(27,09) = pd(27,09) + rrt(321) * density(14) 
  pd(27,14) = pd(27,14) + rrt(321) * density(09) 
  pd(42,09) = pd(42,09) + rrt(321) * density(14) 
  pd(42,14) = pd(42,14) + rrt(321) * density(09) 
  pd(09,09) = pd(09,09) - rrt(322) * density(14) 
  pd(09,14) = pd(09,14) - rrt(322) * density(09) 
  pd(14,09) = pd(14,09) - rrt(322) * density(14) 
  pd(14,14) = pd(14,14) - rrt(322) * density(09) 
  pd(16,09) = pd(16,09) + rrt(322) * density(14) 
  pd(16,14) = pd(16,14) + rrt(322) * density(09) 
  pd(30,09) = pd(30,09) + rrt(322) * density(14) 
  pd(30,14) = pd(30,14) + rrt(322) * density(09) 
  pd(09,09) = pd(09,09) - rrt(323) * density(30) 
  pd(09,30) = pd(09,30) - rrt(323) * density(09) 
  pd(16,09) = pd(16,09) + rrt(323) * density(30) 
  pd(16,30) = pd(16,30) + rrt(323) * density(09) 
  pd(30,09) = pd(30,09) - rrt(323) * density(30) 
  pd(30,30) = pd(30,30) - rrt(323) * density(09) 
  pd(32,09) = pd(32,09) + rrt(323) * density(30) 
  pd(32,30) = pd(32,30) + rrt(323) * density(09) 
  pd(04,09) = pd(04,09) + rrt(324) * density(42) 
  pd(04,42) = pd(04,42) + rrt(324) * density(09) 
  pd(09,09) = pd(09,09) - rrt(324) * density(42) 
  pd(09,42) = pd(09,42) - rrt(324) * density(09) 
  pd(32,09) = pd(32,09) + rrt(324) * density(42) 
  pd(32,42) = pd(32,42) + rrt(324) * density(09) 
  pd(42,09) = pd(42,09) - rrt(324) * density(42) 
  pd(42,42) = pd(42,42) - rrt(324) * density(09) 
  pd(04,09) = pd(04,09) + rrt(325) * density(23) 
  pd(04,23) = pd(04,23) + rrt(325) * density(09) 
  pd(09,09) = pd(09,09) - rrt(325) * density(23) 
  pd(09,23) = pd(09,23) - rrt(325) * density(09) 
  pd(23,09) = pd(23,09) - rrt(325) * density(23) 
  pd(23,23) = pd(23,23) - rrt(325) * density(09) 
  pd(34,09) = pd(34,09) + rrt(325) * density(23) 
  pd(34,23) = pd(34,23) + rrt(325) * density(09) 
  pd(06,06) = pd(06,06) - rrt(326) * density(36) 
  pd(06,36) = pd(06,36) - rrt(326) * density(06) 
  pd(09,06) = pd(09,06) + rrt(326) * density(36) 
  pd(09,36) = pd(09,36) + rrt(326) * density(06) 
  pd(11,06) = pd(11,06) + rrt(326) * density(36) 
  pd(11,36) = pd(11,36) + rrt(326) * density(06) 
  pd(36,06) = pd(36,06) - rrt(326) * density(36) 
  pd(36,36) = pd(36,36) - rrt(326) * density(06) 
  pd(06,06) = pd(06,06) - rrt(327) * density(28) 
  pd(06,28) = pd(06,28) - rrt(327) * density(06) 
  pd(23,06) = pd(23,06) + rrt(327) * density(28) 
  pd(23,28) = pd(23,28) + rrt(327) * density(06) 
  pd(28,06) = pd(28,06) - rrt(327) * density(28) 
  pd(28,28) = pd(28,28) - rrt(327) * density(06) 
  pd(35,06) = pd(35,06) + rrt(327) * density(28) 
  pd(35,28) = pd(35,28) + rrt(327) * density(06) 
  pd(11,11) = pd(11,11) - rrt(328) * density(28) 
  pd(11,28) = pd(11,28) - rrt(328) * density(11) 
  pd(23,11) = pd(23,11) + rrt(328) * density(28) 
  pd(23,28) = pd(23,28) + rrt(328) * density(11) 
  pd(28,11) = pd(28,11) - rrt(328) * density(28) 
  pd(28,28) = pd(28,28) - rrt(328) * density(11) 
  pd(48,11) = pd(48,11) + rrt(328) * density(28) 
  pd(48,28) = pd(48,28) + rrt(328) * density(11) 
  pd(23,28) = pd(23,28) + rrt(329) * density(47) 
  pd(23,47) = pd(23,47) + rrt(329) * density(28) 
  pd(28,28) = pd(28,28) - rrt(329) * density(47) 
  pd(28,47) = pd(28,47) - rrt(329) * density(28) 
  pd(47,28) = pd(47,28) - rrt(329) * density(47) 
  pd(47,47) = pd(47,47) - rrt(329) * density(28) 
  pd(49,28) = pd(49,28) + rrt(329) * density(47) 
  pd(49,47) = pd(49,47) + rrt(329) * density(28) 
  pd(23,28) = pd(23,28) + rrt(330) * density(43) 
  pd(23,43) = pd(23,43) + rrt(330) * density(28) 
  pd(28,28) = pd(28,28) - rrt(330) * density(43) 
  pd(28,43) = pd(28,43) - rrt(330) * density(28) 
  pd(43,28) = pd(43,28) - rrt(330) * density(43) 
  pd(43,43) = pd(43,43) - rrt(330) * density(28) 
  pd(46,28) = pd(46,28) + rrt(330) * density(43) 
  pd(46,43) = pd(46,43) + rrt(330) * density(28) 
  pd(14,14) = pd(14,14) - rrt(331) * density(28) 
  pd(14,28) = pd(14,28) - rrt(331) * density(14) 
  pd(23,14) = pd(23,14) + rrt(331) * density(28) * 2.0d0
  pd(23,28) = pd(23,28) + rrt(331) * density(14) * 2.0d0
  pd(27,14) = pd(27,14) + rrt(331) * density(28) 
  pd(27,28) = pd(27,28) + rrt(331) * density(14) 
  pd(28,14) = pd(28,14) - rrt(331) * density(28) 
  pd(28,28) = pd(28,28) - rrt(331) * density(14) 
  pd(08,08) = pd(08,08) - rrt(332) * density(28) 
  pd(08,28) = pd(08,28) - rrt(332) * density(08) 
  pd(23,08) = pd(23,08) + rrt(332) * density(28) 
  pd(23,28) = pd(23,28) + rrt(332) * density(08) 
  pd(28,08) = pd(28,08) - rrt(332) * density(28) 
  pd(28,28) = pd(28,28) - rrt(332) * density(08) 
  pd(52,08) = pd(52,08) + rrt(332) * density(28) 
  pd(52,28) = pd(52,28) + rrt(332) * density(08) 
  pd(23,28) = pd(23,28) + rrt(333) * density(30) 
  pd(23,30) = pd(23,30) + rrt(333) * density(28) 
  pd(27,28) = pd(27,28) + rrt(333) * density(30) 
  pd(27,30) = pd(27,30) + rrt(333) * density(28) 
  pd(28,28) = pd(28,28) - rrt(333) * density(30) 
  pd(28,30) = pd(28,30) - rrt(333) * density(28) 
  pd(30,28) = pd(30,28) - rrt(333) * density(30) 
  pd(30,30) = pd(30,30) - rrt(333) * density(28) 
  pd(04,28) = pd(04,28) + rrt(334) * density(30) 
  pd(04,30) = pd(04,30) + rrt(334) * density(28) 
  pd(23,28) = pd(23,28) + rrt(334) * density(30) * 2.0d0
  pd(23,30) = pd(23,30) + rrt(334) * density(28) * 2.0d0
  pd(28,28) = pd(28,28) - rrt(334) * density(30) 
  pd(28,30) = pd(28,30) - rrt(334) * density(28) 
  pd(30,28) = pd(30,28) - rrt(334) * density(30) 
  pd(30,30) = pd(30,30) - rrt(334) * density(28) 
  pd(16,28) = pd(16,28) + rrt(335) * density(42) 
  pd(16,42) = pd(16,42) + rrt(335) * density(28) 
  pd(23,28) = pd(23,28) + rrt(335) * density(42) 
  pd(23,42) = pd(23,42) + rrt(335) * density(28) 
  pd(28,28) = pd(28,28) - rrt(335) * density(42) 
  pd(28,42) = pd(28,42) - rrt(335) * density(28) 
  pd(42,28) = pd(42,28) - rrt(335) * density(42) 
  pd(42,42) = pd(42,42) - rrt(335) * density(28) 
  pd(09,28) = pd(09,28) + rrt(336) * density(33) 
  pd(09,33) = pd(09,33) + rrt(336) * density(28) 
  pd(23,28) = pd(23,28) + rrt(336) * density(33) 
  pd(23,33) = pd(23,33) + rrt(336) * density(28) 
  pd(28,28) = pd(28,28) - rrt(336) * density(33) 
  pd(28,33) = pd(28,33) - rrt(336) * density(28) 
  pd(33,28) = pd(33,28) - rrt(336) * density(33) 
  pd(33,33) = pd(33,33) - rrt(336) * density(28) 
  pd(04,28) = pd(04,28) + rrt(337) * density(32) 
  pd(04,32) = pd(04,32) + rrt(337) * density(28) 
  pd(23,28) = pd(23,28) + rrt(337) * density(32) 
  pd(23,32) = pd(23,32) + rrt(337) * density(28) 
  pd(28,28) = pd(28,28) - rrt(337) * density(32) 
  pd(28,32) = pd(28,32) - rrt(337) * density(28) 
  pd(32,28) = pd(32,28) - rrt(337) * density(32) 
  pd(32,32) = pd(32,32) - rrt(337) * density(28) 
  pd(06,06) = pd(06,06) - rrt(338) * density(45) 
  pd(06,45) = pd(06,45) - rrt(338) * density(06) 
  pd(34,06) = pd(34,06) + rrt(338) * density(45) 
  pd(34,45) = pd(34,45) + rrt(338) * density(06) 
  pd(35,06) = pd(35,06) + rrt(338) * density(45) 
  pd(35,45) = pd(35,45) + rrt(338) * density(06) 
  pd(45,06) = pd(45,06) - rrt(338) * density(45) 
  pd(45,45) = pd(45,45) - rrt(338) * density(06) 
  pd(06,06) = pd(06,06) - rrt(339) * density(45) 
  pd(06,45) = pd(06,45) - rrt(339) * density(06) 
  pd(23,06) = pd(23,06) + rrt(339) * density(45) 
  pd(23,45) = pd(23,45) + rrt(339) * density(06) 
  pd(45,06) = pd(45,06) - rrt(339) * density(45) 
  pd(45,45) = pd(45,45) - rrt(339) * density(06) 
  pd(48,06) = pd(48,06) + rrt(339) * density(45) 
  pd(48,45) = pd(48,45) + rrt(339) * density(06) 
  pd(06,06) = pd(06,06) - rrt(340) * density(45) 
  pd(06,45) = pd(06,45) - rrt(340) * density(06) 
  pd(23,06) = pd(23,06) + rrt(340) * density(45) 
  pd(23,45) = pd(23,45) + rrt(340) * density(06) 
  pd(34,06) = pd(34,06) + rrt(340) * density(45) 
  pd(34,45) = pd(34,45) + rrt(340) * density(06) 
  pd(45,06) = pd(45,06) - rrt(340) * density(45) 
  pd(45,45) = pd(45,45) - rrt(340) * density(06) 
  pd(49,06) = pd(49,06) + rrt(340) * density(45) 
  pd(49,45) = pd(49,45) + rrt(340) * density(06) 
  pd(34,45) = pd(34,45) + rrt(341) * density(47) 
  pd(34,47) = pd(34,47) + rrt(341) * density(45) 
  pd(45,45) = pd(45,45) - rrt(341) * density(47) 
  pd(45,47) = pd(45,47) - rrt(341) * density(45) 
  pd(47,45) = pd(47,45) - rrt(341) * density(47) 
  pd(47,47) = pd(47,47) - rrt(341) * density(45) 
  pd(49,45) = pd(49,45) + rrt(341) * density(47) 
  pd(49,47) = pd(49,47) + rrt(341) * density(45) 
  pd(23,45) = pd(23,45) + rrt(342) * density(47) 
  pd(23,47) = pd(23,47) + rrt(342) * density(45) 
  pd(45,45) = pd(45,45) - rrt(342) * density(47) 
  pd(45,47) = pd(45,47) - rrt(342) * density(45) 
  pd(46,45) = pd(46,45) + rrt(342) * density(47) 
  pd(46,47) = pd(46,47) + rrt(342) * density(45) 
  pd(47,45) = pd(47,45) - rrt(342) * density(47) 
  pd(47,47) = pd(47,47) - rrt(342) * density(45) 
  pd(34,43) = pd(34,43) + rrt(343) * density(45) 
  pd(34,45) = pd(34,45) + rrt(343) * density(43) 
  pd(43,43) = pd(43,43) - rrt(343) * density(45) 
  pd(43,45) = pd(43,45) - rrt(343) * density(43) 
  pd(45,43) = pd(45,43) - rrt(343) * density(45) 
  pd(45,45) = pd(45,45) - rrt(343) * density(43) 
  pd(46,43) = pd(46,43) + rrt(343) * density(45) 
  pd(46,45) = pd(46,45) + rrt(343) * density(43) 
  pd(13,43) = pd(13,43) + rrt(344) * density(45) 
  pd(13,45) = pd(13,45) + rrt(344) * density(43) 
  pd(23,43) = pd(23,43) + rrt(344) * density(45) 
  pd(23,45) = pd(23,45) + rrt(344) * density(43) 
  pd(43,43) = pd(43,43) - rrt(344) * density(45) 
  pd(43,45) = pd(43,45) - rrt(344) * density(43) 
  pd(45,43) = pd(45,43) - rrt(344) * density(45) 
  pd(45,45) = pd(45,45) - rrt(344) * density(43) 
  pd(14,14) = pd(14,14) - rrt(345) * density(45) 
  pd(14,45) = pd(14,45) - rrt(345) * density(14) 
  pd(23,14) = pd(23,14) + rrt(345) * density(45) 
  pd(23,45) = pd(23,45) + rrt(345) * density(14) 
  pd(45,14) = pd(45,14) - rrt(345) * density(45) 
  pd(45,45) = pd(45,45) - rrt(345) * density(14) 
  pd(52,14) = pd(52,14) + rrt(345) * density(45) 
  pd(52,45) = pd(52,45) + rrt(345) * density(14) 
  pd(14,14) = pd(14,14) - rrt(346) * density(45) 
  pd(14,45) = pd(14,45) - rrt(346) * density(14) 
  pd(23,14) = pd(23,14) + rrt(346) * density(45) 
  pd(23,45) = pd(23,45) + rrt(346) * density(14) 
  pd(27,14) = pd(27,14) + rrt(346) * density(45) 
  pd(27,45) = pd(27,45) + rrt(346) * density(14) 
  pd(34,14) = pd(34,14) + rrt(346) * density(45) 
  pd(34,45) = pd(34,45) + rrt(346) * density(14) 
  pd(45,14) = pd(45,14) - rrt(346) * density(45) 
  pd(45,45) = pd(45,45) - rrt(346) * density(14) 
  pd(14,14) = pd(14,14) - rrt(347) * density(45) 
  pd(14,45) = pd(14,45) - rrt(347) * density(14) 
  pd(16,14) = pd(16,14) + rrt(347) * density(45) 
  pd(16,45) = pd(16,45) + rrt(347) * density(14) 
  pd(23,14) = pd(23,14) + rrt(347) * density(45) * 2.0d0
  pd(23,45) = pd(23,45) + rrt(347) * density(14) * 2.0d0
  pd(45,14) = pd(45,14) - rrt(347) * density(45) 
  pd(45,45) = pd(45,45) - rrt(347) * density(14) 
  pd(04,14) = pd(04,14) + rrt(348) * density(45) 
  pd(04,45) = pd(04,45) + rrt(348) * density(14) 
  pd(14,14) = pd(14,14) - rrt(348) * density(45) 
  pd(14,45) = pd(14,45) - rrt(348) * density(14) 
  pd(23,14) = pd(23,14) + rrt(348) * density(45) * 2.0d0
  pd(23,45) = pd(23,45) + rrt(348) * density(14) * 2.0d0
  pd(34,14) = pd(34,14) + rrt(348) * density(45) 
  pd(34,45) = pd(34,45) + rrt(348) * density(14) 
  pd(45,14) = pd(45,14) - rrt(348) * density(45) 
  pd(45,45) = pd(45,45) - rrt(348) * density(14) 
  pd(09,14) = pd(09,14) + rrt(349) * density(45) 
  pd(09,45) = pd(09,45) + rrt(349) * density(14) 
  pd(14,14) = pd(14,14) - rrt(349) * density(45) 
  pd(14,45) = pd(14,45) - rrt(349) * density(14) 
  pd(23,14) = pd(23,14) + rrt(349) * density(45) * 3.0d0
  pd(23,45) = pd(23,45) + rrt(349) * density(14) * 3.0d0
  pd(45,14) = pd(45,14) - rrt(349) * density(45) 
  pd(45,45) = pd(45,45) - rrt(349) * density(14) 
  pd(16,30) = pd(16,30) + rrt(350) * density(45) 
  pd(16,45) = pd(16,45) + rrt(350) * density(30) 
  pd(23,30) = pd(23,30) + rrt(350) * density(45) 
  pd(23,45) = pd(23,45) + rrt(350) * density(30) 
  pd(30,30) = pd(30,30) - rrt(350) * density(45) 
  pd(30,45) = pd(30,45) - rrt(350) * density(30) 
  pd(45,30) = pd(45,30) - rrt(350) * density(45) 
  pd(45,45) = pd(45,45) - rrt(350) * density(30) 
  pd(04,30) = pd(04,30) + rrt(351) * density(45) 
  pd(04,45) = pd(04,45) + rrt(351) * density(30) 
  pd(23,30) = pd(23,30) + rrt(351) * density(45) 
  pd(23,45) = pd(23,45) + rrt(351) * density(30) 
  pd(30,30) = pd(30,30) - rrt(351) * density(45) 
  pd(30,45) = pd(30,45) - rrt(351) * density(30) 
  pd(34,30) = pd(34,30) + rrt(351) * density(45) 
  pd(34,45) = pd(34,45) + rrt(351) * density(30) 
  pd(45,30) = pd(45,30) - rrt(351) * density(45) 
  pd(45,45) = pd(45,45) - rrt(351) * density(30) 
  pd(09,30) = pd(09,30) + rrt(352) * density(45) 
  pd(09,45) = pd(09,45) + rrt(352) * density(30) 
  pd(23,30) = pd(23,30) + rrt(352) * density(45) * 2.0d0
  pd(23,45) = pd(23,45) + rrt(352) * density(30) * 2.0d0
  pd(30,30) = pd(30,30) - rrt(352) * density(45) 
  pd(30,45) = pd(30,45) - rrt(352) * density(30) 
  pd(45,30) = pd(45,30) - rrt(352) * density(45) 
  pd(45,45) = pd(45,45) - rrt(352) * density(30) 
  pd(04,32) = pd(04,32) + rrt(353) * density(45) 
  pd(04,45) = pd(04,45) + rrt(353) * density(32) 
  pd(32,32) = pd(32,32) - rrt(353) * density(45) 
  pd(32,45) = pd(32,45) - rrt(353) * density(32) 
  pd(34,32) = pd(34,32) + rrt(353) * density(45) 
  pd(34,45) = pd(34,45) + rrt(353) * density(32) 
  pd(45,32) = pd(45,32) - rrt(353) * density(45) 
  pd(45,45) = pd(45,45) - rrt(353) * density(32) 
  pd(09,32) = pd(09,32) + rrt(354) * density(45) 
  pd(09,45) = pd(09,45) + rrt(354) * density(32) 
  pd(23,32) = pd(23,32) + rrt(354) * density(45) 
  pd(23,45) = pd(23,45) + rrt(354) * density(32) 
  pd(32,32) = pd(32,32) - rrt(354) * density(45) 
  pd(32,45) = pd(32,45) - rrt(354) * density(32) 
  pd(45,32) = pd(45,32) - rrt(354) * density(45) 
  pd(45,45) = pd(45,45) - rrt(354) * density(32) 
  pd(28,34) = pd(28,34) + rrt(355) * density(45) 
  pd(28,45) = pd(28,45) + rrt(355) * density(34) 
  pd(34,34) = pd(34,34) - rrt(355) * density(45) 
  pd(34,45) = pd(34,45) - rrt(355) * density(34) 
  pd(45,34) = pd(45,34) - rrt(355) * density(45) 
  pd(45,45) = pd(45,45) - rrt(355) * density(34) 
  pd(12,34) = pd(12,34) + rrt(356) * density(45) 
  pd(12,45) = pd(12,45) + rrt(356) * density(34) 
  pd(23,34) = pd(23,34) + rrt(356) * density(45) 
  pd(23,45) = pd(23,45) + rrt(356) * density(34) 
  pd(34,34) = pd(34,34) - rrt(356) * density(45) 
  pd(34,45) = pd(34,45) - rrt(356) * density(34) 
  pd(45,34) = pd(45,34) - rrt(356) * density(45) 
  pd(45,45) = pd(45,45) - rrt(356) * density(34) 
  pd(06,06) = pd(06,06) - rrt(357) * density(12) 
  pd(06,12) = pd(06,12) - rrt(357) * density(06) 
  pd(12,06) = pd(12,06) - rrt(357) * density(12) 
  pd(12,12) = pd(12,12) - rrt(357) * density(06) 
  pd(34,06) = pd(34,06) + rrt(357) * density(12) 
  pd(34,12) = pd(34,12) + rrt(357) * density(06) 
  pd(48,06) = pd(48,06) + rrt(357) * density(12) 
  pd(48,12) = pd(48,12) + rrt(357) * density(06) 
  pd(06,06) = pd(06,06) - rrt(358) * density(12) 
  pd(06,12) = pd(06,12) - rrt(358) * density(06) 
  pd(12,06) = pd(12,06) - rrt(358) * density(12) 
  pd(12,12) = pd(12,12) - rrt(358) * density(06) 
  pd(23,06) = pd(23,06) + rrt(358) * density(12) 
  pd(23,12) = pd(23,12) + rrt(358) * density(06) 
  pd(49,06) = pd(49,06) + rrt(358) * density(12) 
  pd(49,12) = pd(49,12) + rrt(358) * density(06) 
  pd(11,11) = pd(11,11) - rrt(359) * density(12) 
  pd(11,12) = pd(11,12) - rrt(359) * density(11) 
  pd(12,11) = pd(12,11) - rrt(359) * density(12) 
  pd(12,12) = pd(12,12) - rrt(359) * density(11) 
  pd(34,11) = pd(34,11) + rrt(359) * density(12) 
  pd(34,12) = pd(34,12) + rrt(359) * density(11) 
  pd(49,11) = pd(49,11) + rrt(359) * density(12) 
  pd(49,12) = pd(49,12) + rrt(359) * density(11) 
  pd(12,12) = pd(12,12) - rrt(360) * density(47) 
  pd(12,47) = pd(12,47) - rrt(360) * density(12) 
  pd(34,12) = pd(34,12) + rrt(360) * density(47) 
  pd(34,47) = pd(34,47) + rrt(360) * density(12) 
  pd(46,12) = pd(46,12) + rrt(360) * density(47) 
  pd(46,47) = pd(46,47) + rrt(360) * density(12) 
  pd(47,12) = pd(47,12) - rrt(360) * density(47) 
  pd(47,47) = pd(47,47) - rrt(360) * density(12) 
  pd(12,12) = pd(12,12) - rrt(361) * density(47) 
  pd(12,47) = pd(12,47) - rrt(361) * density(12) 
  pd(13,12) = pd(13,12) + rrt(361) * density(47) 
  pd(13,47) = pd(13,47) + rrt(361) * density(12) 
  pd(23,12) = pd(23,12) + rrt(361) * density(47) 
  pd(23,47) = pd(23,47) + rrt(361) * density(12) 
  pd(47,12) = pd(47,12) - rrt(361) * density(47) 
  pd(47,47) = pd(47,47) - rrt(361) * density(12) 
  pd(12,12) = pd(12,12) - rrt(362) * density(43) 
  pd(12,43) = pd(12,43) - rrt(362) * density(12) 
  pd(13,12) = pd(13,12) + rrt(362) * density(43) 
  pd(13,43) = pd(13,43) + rrt(362) * density(12) 
  pd(34,12) = pd(34,12) + rrt(362) * density(43) 
  pd(34,43) = pd(34,43) + rrt(362) * density(12) 
  pd(43,12) = pd(43,12) - rrt(362) * density(43) 
  pd(43,43) = pd(43,43) - rrt(362) * density(12) 
  pd(12,12) = pd(12,12) - rrt(363) * density(14) 
  pd(12,14) = pd(12,14) - rrt(363) * density(12) 
  pd(14,12) = pd(14,12) - rrt(363) * density(14) 
  pd(14,14) = pd(14,14) - rrt(363) * density(12) 
  pd(23,12) = pd(23,12) + rrt(363) * density(14) 
  pd(23,14) = pd(23,14) + rrt(363) * density(12) 
  pd(27,12) = pd(27,12) + rrt(363) * density(14) 
  pd(27,14) = pd(27,14) + rrt(363) * density(12) 
  pd(12,12) = pd(12,12) - rrt(364) * density(14) 
  pd(12,14) = pd(12,14) - rrt(364) * density(12) 
  pd(14,12) = pd(14,12) - rrt(364) * density(14) 
  pd(14,14) = pd(14,14) - rrt(364) * density(12) 
  pd(16,12) = pd(16,12) + rrt(364) * density(14) 
  pd(16,14) = pd(16,14) + rrt(364) * density(12) 
  pd(23,12) = pd(23,12) + rrt(364) * density(14) 
  pd(23,14) = pd(23,14) + rrt(364) * density(12) 
  pd(34,12) = pd(34,12) + rrt(364) * density(14) 
  pd(34,14) = pd(34,14) + rrt(364) * density(12) 
  pd(04,12) = pd(04,12) + rrt(365) * density(14) 
  pd(04,14) = pd(04,14) + rrt(365) * density(12) 
  pd(12,12) = pd(12,12) - rrt(365) * density(14) 
  pd(12,14) = pd(12,14) - rrt(365) * density(12) 
  pd(14,12) = pd(14,12) - rrt(365) * density(14) 
  pd(14,14) = pd(14,14) - rrt(365) * density(12) 
  pd(23,12) = pd(23,12) + rrt(365) * density(14) * 2.0d0
  pd(23,14) = pd(23,14) + rrt(365) * density(12) * 2.0d0
  pd(08,08) = pd(08,08) - rrt(366) * density(12) 
  pd(08,12) = pd(08,12) - rrt(366) * density(08) 
  pd(12,08) = pd(12,08) - rrt(366) * density(12) 
  pd(12,12) = pd(12,12) - rrt(366) * density(08) 
  pd(16,08) = pd(16,08) + rrt(366) * density(12) 
  pd(16,12) = pd(16,12) + rrt(366) * density(08) 
  pd(23,08) = pd(23,08) + rrt(366) * density(12) 
  pd(23,12) = pd(23,12) + rrt(366) * density(08) 
  pd(04,08) = pd(04,08) + rrt(367) * density(12) 
  pd(04,12) = pd(04,12) + rrt(367) * density(08) 
  pd(08,08) = pd(08,08) - rrt(367) * density(12) 
  pd(08,12) = pd(08,12) - rrt(367) * density(08) 
  pd(12,08) = pd(12,08) - rrt(367) * density(12) 
  pd(12,12) = pd(12,12) - rrt(367) * density(08) 
  pd(23,08) = pd(23,08) + rrt(367) * density(12) 
  pd(23,12) = pd(23,12) + rrt(367) * density(08) 
  pd(34,08) = pd(34,08) + rrt(367) * density(12) 
  pd(34,12) = pd(34,12) + rrt(367) * density(08) 
  pd(12,12) = pd(12,12) - rrt(368) * density(30) 
  pd(12,30) = pd(12,30) - rrt(368) * density(12) 
  pd(16,12) = pd(16,12) + rrt(368) * density(30) 
  pd(16,30) = pd(16,30) + rrt(368) * density(12) 
  pd(30,12) = pd(30,12) - rrt(368) * density(30) 
  pd(30,30) = pd(30,30) - rrt(368) * density(12) 
  pd(34,12) = pd(34,12) + rrt(368) * density(30) 
  pd(34,30) = pd(34,30) + rrt(368) * density(12) 
  pd(04,12) = pd(04,12) + rrt(369) * density(30) 
  pd(04,30) = pd(04,30) + rrt(369) * density(12) 
  pd(12,12) = pd(12,12) - rrt(369) * density(30) 
  pd(12,30) = pd(12,30) - rrt(369) * density(12) 
  pd(23,12) = pd(23,12) + rrt(369) * density(30) 
  pd(23,30) = pd(23,30) + rrt(369) * density(12) 
  pd(30,12) = pd(30,12) - rrt(369) * density(30) 
  pd(30,30) = pd(30,30) - rrt(369) * density(12) 
  pd(09,12) = pd(09,12) + rrt(370) * density(30) 
  pd(09,30) = pd(09,30) + rrt(370) * density(12) 
  pd(12,12) = pd(12,12) - rrt(370) * density(30) 
  pd(12,30) = pd(12,30) - rrt(370) * density(12) 
  pd(23,12) = pd(23,12) + rrt(370) * density(30) 
  pd(23,30) = pd(23,30) + rrt(370) * density(12) 
  pd(30,12) = pd(30,12) - rrt(370) * density(30) 
  pd(30,30) = pd(30,30) - rrt(370) * density(12) 
  pd(34,12) = pd(34,12) + rrt(370) * density(30) 
  pd(34,30) = pd(34,30) + rrt(370) * density(12) 
  pd(04,12) = pd(04,12) + rrt(371) * density(42) 
  pd(04,42) = pd(04,42) + rrt(371) * density(12) 
  pd(12,12) = pd(12,12) - rrt(371) * density(42) 
  pd(12,42) = pd(12,42) - rrt(371) * density(12) 
  pd(34,12) = pd(34,12) + rrt(371) * density(42) 
  pd(34,42) = pd(34,42) + rrt(371) * density(12) 
  pd(42,12) = pd(42,12) - rrt(371) * density(42) 
  pd(42,42) = pd(42,42) - rrt(371) * density(12) 
  pd(09,12) = pd(09,12) + rrt(372) * density(42) 
  pd(09,42) = pd(09,42) + rrt(372) * density(12) 
  pd(12,12) = pd(12,12) - rrt(372) * density(42) 
  pd(12,42) = pd(12,42) - rrt(372) * density(12) 
  pd(23,12) = pd(23,12) + rrt(372) * density(42) 
  pd(23,42) = pd(23,42) + rrt(372) * density(12) 
  pd(42,12) = pd(42,12) - rrt(372) * density(42) 
  pd(42,42) = pd(42,42) - rrt(372) * density(12) 
  pd(09,12) = pd(09,12) + rrt(373) * density(32) 
  pd(09,32) = pd(09,32) + rrt(373) * density(12) 
  pd(12,12) = pd(12,12) - rrt(373) * density(32) 
  pd(12,32) = pd(12,32) - rrt(373) * density(12) 
  pd(32,12) = pd(32,12) - rrt(373) * density(32) 
  pd(32,32) = pd(32,32) - rrt(373) * density(12) 
  pd(34,12) = pd(34,12) + rrt(373) * density(32) 
  pd(34,32) = pd(34,32) + rrt(373) * density(12) 
  pd(06,06) = pd(06,06) - rrt(374) * density(47) 
  pd(06,47) = pd(06,47) - rrt(374) * density(06) 
  pd(11,06) = pd(11,06) + rrt(374) * density(47) * 2.0d0
  pd(11,47) = pd(11,47) + rrt(374) * density(06) * 2.0d0
  pd(47,06) = pd(47,06) - rrt(374) * density(47) 
  pd(47,47) = pd(47,47) - rrt(374) * density(06) 
  pd(06,06) = pd(06,06) - rrt(375) * density(43) 
  pd(06,43) = pd(06,43) - rrt(375) * density(06) 
  pd(30,06) = pd(30,06) + rrt(375) * density(43) 
  pd(30,43) = pd(30,43) + rrt(375) * density(06) 
  pd(34,06) = pd(34,06) + rrt(375) * density(43) 
  pd(34,43) = pd(34,43) + rrt(375) * density(06) 
  pd(43,06) = pd(43,06) - rrt(375) * density(43) 
  pd(43,43) = pd(43,43) - rrt(375) * density(06) 
  pd(06,06) = pd(06,06) - rrt(376) * density(08) 
  pd(06,08) = pd(06,08) - rrt(376) * density(06) 
  pd(08,06) = pd(08,06) - rrt(376) * density(08) 
  pd(08,08) = pd(08,08) - rrt(376) * density(06) 
  pd(11,06) = pd(11,06) + rrt(376) * density(08) 
  pd(11,08) = pd(11,08) + rrt(376) * density(06) 
  pd(14,06) = pd(14,06) + rrt(376) * density(08) 
  pd(14,08) = pd(14,08) + rrt(376) * density(06) 
  pd(06,06) = pd(06,06) - rrt(377) * density(42) 
  pd(06,42) = pd(06,42) - rrt(377) * density(06) 
  pd(11,06) = pd(11,06) + rrt(377) * density(42) 
  pd(11,42) = pd(11,42) + rrt(377) * density(06) 
  pd(30,06) = pd(30,06) + rrt(377) * density(42) 
  pd(30,42) = pd(30,42) + rrt(377) * density(06) 
  pd(42,06) = pd(42,06) - rrt(377) * density(42) 
  pd(42,42) = pd(42,42) - rrt(377) * density(06) 
  pd(06,06) = pd(06,06) - rrt(378) * density(33) 
  pd(06,33) = pd(06,33) - rrt(378) * density(06) 
  pd(11,06) = pd(11,06) + rrt(378) * density(33) 
  pd(11,33) = pd(11,33) + rrt(378) * density(06) 
  pd(32,06) = pd(32,06) + rrt(378) * density(33) 
  pd(32,33) = pd(32,33) + rrt(378) * density(06) 
  pd(33,06) = pd(33,06) - rrt(378) * density(33) 
  pd(33,33) = pd(33,33) - rrt(378) * density(06) 
  pd(06,06) = pd(06,06) - rrt(379) * density(44) 
  pd(06,44) = pd(06,44) - rrt(379) * density(06) 
  pd(07,06) = pd(07,06) + rrt(379) * density(44) 
  pd(07,44) = pd(07,44) + rrt(379) * density(06) 
  pd(11,06) = pd(11,06) + rrt(379) * density(44) 
  pd(11,44) = pd(11,44) + rrt(379) * density(06) 
  pd(44,06) = pd(44,06) - rrt(379) * density(44) 
  pd(44,44) = pd(44,44) - rrt(379) * density(06) 
  pd(06,06) = pd(06,06) - rrt(380) * density(26) 
  pd(06,26) = pd(06,26) - rrt(380) * density(06) 
  pd(11,06) = pd(11,06) + rrt(380) * density(26) 
  pd(11,26) = pd(11,26) + rrt(380) * density(06) 
  pd(21,06) = pd(21,06) + rrt(380) * density(26) 
  pd(21,26) = pd(21,26) + rrt(380) * density(06) 
  pd(26,06) = pd(26,06) - rrt(380) * density(26) 
  pd(26,26) = pd(26,26) - rrt(380) * density(06) 
  pd(06,06) = pd(06,06) - rrt(381) * density(34) 
  pd(06,34) = pd(06,34) - rrt(381) * density(06) 
  pd(11,06) = pd(11,06) + rrt(381) * density(34) 
  pd(11,34) = pd(11,34) + rrt(381) * density(06) 
  pd(23,06) = pd(23,06) + rrt(381) * density(34) 
  pd(23,34) = pd(23,34) + rrt(381) * density(06) 
  pd(34,06) = pd(34,06) - rrt(381) * density(34) 
  pd(34,34) = pd(34,34) - rrt(381) * density(06) 
  pd(06,06) = pd(06,06) - rrt(382) * density(11) 
  pd(06,11) = pd(06,11) - rrt(382) * density(06) 
  pd(11,06) = pd(11,06) - rrt(382) * density(11) 
  pd(11,11) = pd(11,11) - rrt(382) * density(06) 
  pd(14,06) = pd(14,06) + rrt(382) * density(11) 
  pd(14,11) = pd(14,11) + rrt(382) * density(06) 
  pd(34,06) = pd(34,06) + rrt(382) * density(11) 
  pd(34,11) = pd(34,11) + rrt(382) * density(06) 
  pd(06,06) = pd(06,06) - rrt(383) * density(31) 
  pd(06,31) = pd(06,31) - rrt(383) * density(06) 
  pd(11,06) = pd(11,06) + rrt(383) * density(31) 
  pd(11,31) = pd(11,31) + rrt(383) * density(06) 
  pd(25,06) = pd(25,06) + rrt(383) * density(31) 
  pd(25,31) = pd(25,31) + rrt(383) * density(06) 
  pd(31,06) = pd(31,06) - rrt(383) * density(31) 
  pd(31,31) = pd(31,31) - rrt(383) * density(06) 
  pd(06,06) = pd(06,06) - rrt(384) * density(47) 
  pd(06,47) = pd(06,47) - rrt(384) * density(06) 
  pd(14,06) = pd(14,06) + rrt(384) * density(47) 
  pd(14,47) = pd(14,47) + rrt(384) * density(06) 
  pd(47,06) = pd(47,06) - rrt(384) * density(47) 
  pd(47,47) = pd(47,47) - rrt(384) * density(06) 
  pd(06,06) = pd(06,06) - rrt(385) 
  pd(11,06) = pd(11,06) + rrt(385) 
  pd(34,06) = pd(34,06) + rrt(385) 
  pd(11,11) = pd(11,11) - rrt(386) 
  pd(34,11) = pd(34,11) + rrt(386) 
  pd(47,11) = pd(47,11) + rrt(386) 
  pd(11,11) = pd(11,11) - rrt(387) 
  pd(23,11) = pd(23,11) + rrt(387) 
  pd(43,11) = pd(43,11) + rrt(387) 
  pd(08,08) = pd(08,08) - rrt(388) * density(11) 
  pd(08,11) = pd(08,11) - rrt(388) * density(08) 
  pd(11,08) = pd(11,08) - rrt(388) * density(11) 
  pd(11,11) = pd(11,11) - rrt(388) * density(08) 
  pd(14,08) = pd(14,08) + rrt(388) * density(11) 
  pd(14,11) = pd(14,11) + rrt(388) * density(08) 
  pd(47,08) = pd(47,08) + rrt(388) * density(11) 
  pd(47,11) = pd(47,11) + rrt(388) * density(08) 
  pd(32,47) = pd(32,47) + rrt(389) * density(47) * 2.0d0
  pd(34,47) = pd(34,47) + rrt(389) * density(47) * 4.0d0
  pd(47,47) = pd(47,47) - rrt(389) * density(47) * 4.0d0
  pd(08,08) = pd(08,08) - rrt(390) * density(47) 
  pd(08,47) = pd(08,47) - rrt(390) * density(08) 
  pd(11,08) = pd(11,08) + rrt(390) * density(47) 
  pd(11,47) = pd(11,47) + rrt(390) * density(08) 
  pd(30,08) = pd(30,08) + rrt(390) * density(47) 
  pd(30,47) = pd(30,47) + rrt(390) * density(08) 
  pd(47,08) = pd(47,08) - rrt(390) * density(47) 
  pd(47,47) = pd(47,47) - rrt(390) * density(08) 
  pd(11,42) = pd(11,42) + rrt(391) * density(47) 
  pd(11,47) = pd(11,47) + rrt(391) * density(42) 
  pd(32,42) = pd(32,42) + rrt(391) * density(47) 
  pd(32,47) = pd(32,47) + rrt(391) * density(42) 
  pd(42,42) = pd(42,42) - rrt(391) * density(47) 
  pd(42,47) = pd(42,47) - rrt(391) * density(42) 
  pd(47,42) = pd(47,42) - rrt(391) * density(47) 
  pd(47,47) = pd(47,47) - rrt(391) * density(42) 
  pd(32,33) = pd(32,33) + rrt(392) * density(47) 
  pd(32,47) = pd(32,47) + rrt(392) * density(33) 
  pd(33,33) = pd(33,33) - rrt(392) * density(47) 
  pd(33,47) = pd(33,47) - rrt(392) * density(33) 
  pd(43,33) = pd(43,33) + rrt(392) * density(47) 
  pd(43,47) = pd(43,47) + rrt(392) * density(33) 
  pd(47,33) = pd(47,33) - rrt(392) * density(47) 
  pd(47,47) = pd(47,47) - rrt(392) * density(33) 
  pd(07,07) = pd(07,07) - rrt(393) * density(47) 
  pd(07,47) = pd(07,47) - rrt(393) * density(07) 
  pd(11,07) = pd(11,07) + rrt(393) * density(47) 
  pd(11,47) = pd(11,47) + rrt(393) * density(07) 
  pd(44,07) = pd(44,07) + rrt(393) * density(47) 
  pd(44,47) = pd(44,47) + rrt(393) * density(07) 
  pd(47,07) = pd(47,07) - rrt(393) * density(47) 
  pd(47,47) = pd(47,47) - rrt(393) * density(07) 
  pd(08,44) = pd(08,44) + rrt(394) * density(47) 
  pd(08,47) = pd(08,47) + rrt(394) * density(44) 
  pd(30,44) = pd(30,44) + rrt(394) * density(47) 
  pd(30,47) = pd(30,47) + rrt(394) * density(44) 
  pd(44,44) = pd(44,44) - rrt(394) * density(47) 
  pd(44,47) = pd(44,47) - rrt(394) * density(44) 
  pd(47,44) = pd(47,44) - rrt(394) * density(47) 
  pd(47,47) = pd(47,47) - rrt(394) * density(44) 
  pd(11,44) = pd(11,44) + rrt(395) * density(47) 
  pd(11,47) = pd(11,47) + rrt(395) * density(44) 
  pd(21,44) = pd(21,44) + rrt(395) * density(47) 
  pd(21,47) = pd(21,47) + rrt(395) * density(44) 
  pd(44,44) = pd(44,44) - rrt(395) * density(47) 
  pd(44,47) = pd(44,47) - rrt(395) * density(44) 
  pd(47,44) = pd(47,44) - rrt(395) * density(47) 
  pd(47,47) = pd(47,47) - rrt(395) * density(44) 
  pd(11,21) = pd(11,21) + rrt(396) * density(47) 
  pd(11,47) = pd(11,47) + rrt(396) * density(21) 
  pd(21,21) = pd(21,21) - rrt(396) * density(47) 
  pd(21,47) = pd(21,47) - rrt(396) * density(21) 
  pd(26,21) = pd(26,21) + rrt(396) * density(47) 
  pd(26,47) = pd(26,47) + rrt(396) * density(21) 
  pd(47,21) = pd(47,21) - rrt(396) * density(47) 
  pd(47,47) = pd(47,47) - rrt(396) * density(21) 
  pd(11,23) = pd(11,23) + rrt(397) * density(47) 
  pd(11,47) = pd(11,47) + rrt(397) * density(23) 
  pd(23,23) = pd(23,23) - rrt(397) * density(47) 
  pd(23,47) = pd(23,47) - rrt(397) * density(23) 
  pd(34,23) = pd(34,23) + rrt(397) * density(47) 
  pd(34,47) = pd(34,47) + rrt(397) * density(23) 
  pd(47,23) = pd(47,23) - rrt(397) * density(47) 
  pd(47,47) = pd(47,47) - rrt(397) * density(23) 
  pd(23,34) = pd(23,34) + rrt(398) * density(47) 
  pd(23,47) = pd(23,47) + rrt(398) * density(34) 
  pd(34,34) = pd(34,34) - rrt(398) * density(47) 
  pd(34,47) = pd(34,47) - rrt(398) * density(34) 
  pd(43,34) = pd(43,34) + rrt(398) * density(47) 
  pd(43,47) = pd(43,47) + rrt(398) * density(34) 
  pd(47,34) = pd(47,34) - rrt(398) * density(47) 
  pd(47,47) = pd(47,47) - rrt(398) * density(34) 
  pd(34,47) = pd(34,47) + rrt(399) 
  pd(43,47) = pd(43,47) + rrt(399) 
  pd(47,47) = pd(47,47) - rrt(399) 
  pd(23,47) = pd(23,47) + rrt(400) * density(47) * 2.0d0
  pd(32,47) = pd(32,47) + rrt(400) * density(47) * 2.0d0
  pd(47,47) = pd(47,47) - rrt(400) * density(47) * 4.0d0
  pd(11,34) = pd(11,34) + rrt(401) * density(47) 
  pd(11,47) = pd(11,47) + rrt(401) * density(34) 
  pd(34,34) = pd(34,34) - rrt(401) * density(47) 
  pd(34,47) = pd(34,47) - rrt(401) * density(34) 
  pd(47,34) = pd(47,34) - rrt(401) * density(47) 
  pd(47,47) = pd(47,47) - rrt(401) * density(34) 
  pd(14,14) = pd(14,14) - rrt(402) * density(43) 
  pd(14,43) = pd(14,43) - rrt(402) * density(14) 
  pd(21,14) = pd(21,14) + rrt(402) * density(43) 
  pd(21,43) = pd(21,43) + rrt(402) * density(14) 
  pd(34,14) = pd(34,14) + rrt(402) * density(43) 
  pd(34,43) = pd(34,43) + rrt(402) * density(14) 
  pd(43,14) = pd(43,14) - rrt(402) * density(43) 
  pd(43,43) = pd(43,43) - rrt(402) * density(14) 
  pd(14,14) = pd(14,14) - rrt(403) * density(43) 
  pd(14,43) = pd(14,43) - rrt(403) * density(14) 
  pd(43,14) = pd(43,14) - rrt(403) * density(43) 
  pd(43,43) = pd(43,43) - rrt(403) * density(14) 
  pd(44,14) = pd(44,14) + rrt(403) * density(43) 
  pd(44,43) = pd(44,43) + rrt(403) * density(14) 
  pd(23,23) = pd(23,23) - rrt(404) * density(43) 
  pd(23,43) = pd(23,43) - rrt(404) * density(23) 
  pd(34,23) = pd(34,23) + rrt(404) * density(43) 
  pd(34,43) = pd(34,43) + rrt(404) * density(23) 
  pd(43,23) = pd(43,23) - rrt(404) * density(43) 
  pd(43,43) = pd(43,43) - rrt(404) * density(23) 
  pd(47,23) = pd(47,23) + rrt(404) * density(43) 
  pd(47,43) = pd(47,43) + rrt(404) * density(23) 
  pd(11,11) = pd(11,11) - rrt(405) * density(43) 
  pd(11,43) = pd(11,43) - rrt(405) * density(11) 
  pd(34,11) = pd(34,11) + rrt(405) * density(43) 
  pd(34,43) = pd(34,43) + rrt(405) * density(11) 
  pd(42,11) = pd(42,11) + rrt(405) * density(43) 
  pd(42,43) = pd(42,43) + rrt(405) * density(11) 
  pd(43,11) = pd(43,11) - rrt(405) * density(43) 
  pd(43,43) = pd(43,43) - rrt(405) * density(11) 
  pd(32,43) = pd(32,43) + rrt(406) * density(47) 
  pd(32,47) = pd(32,47) + rrt(406) * density(43) 
  pd(34,43) = pd(34,43) + rrt(406) * density(47) 
  pd(34,47) = pd(34,47) + rrt(406) * density(43) 
  pd(43,43) = pd(43,43) - rrt(406) * density(47) 
  pd(43,47) = pd(43,47) - rrt(406) * density(43) 
  pd(47,43) = pd(47,43) - rrt(406) * density(47) 
  pd(47,47) = pd(47,47) - rrt(406) * density(43) 
  pd(11,23) = pd(11,23) + rrt(407) * density(43) 
  pd(11,43) = pd(11,43) + rrt(407) * density(23) 
  pd(23,23) = pd(23,23) - rrt(407) * density(43) 
  pd(23,43) = pd(23,43) - rrt(407) * density(23) 
  pd(43,23) = pd(43,23) - rrt(407) * density(43) 
  pd(43,43) = pd(43,43) - rrt(407) * density(23) 
  pd(32,42) = pd(32,42) + rrt(408) * density(43) 
  pd(32,43) = pd(32,43) + rrt(408) * density(42) 
  pd(42,42) = pd(42,42) - rrt(408) * density(43) 
  pd(42,43) = pd(42,43) - rrt(408) * density(42) 
  pd(43,42) = pd(43,42) - rrt(408) * density(43) 
  pd(43,43) = pd(43,43) - rrt(408) * density(42) 
  pd(47,42) = pd(47,42) + rrt(408) * density(43) 
  pd(47,43) = pd(47,43) + rrt(408) * density(42) 
  pd(08,14) = pd(08,14) + rrt(409) * density(42) 
  pd(08,42) = pd(08,42) + rrt(409) * density(14) 
  pd(14,14) = pd(14,14) - rrt(409) * density(42) 
  pd(14,42) = pd(14,42) - rrt(409) * density(14) 
  pd(30,14) = pd(30,14) + rrt(409) * density(42) 
  pd(30,42) = pd(30,42) + rrt(409) * density(14) 
  pd(42,14) = pd(42,14) - rrt(409) * density(42) 
  pd(42,42) = pd(42,42) - rrt(409) * density(14) 
  pd(07,14) = pd(07,14) + rrt(410) * density(44) 
  pd(07,44) = pd(07,44) + rrt(410) * density(14) 
  pd(08,14) = pd(08,14) + rrt(410) * density(44) 
  pd(08,44) = pd(08,44) + rrt(410) * density(14) 
  pd(14,14) = pd(14,14) - rrt(410) * density(44) 
  pd(14,44) = pd(14,44) - rrt(410) * density(14) 
  pd(44,14) = pd(44,14) - rrt(410) * density(44) 
  pd(44,44) = pd(44,44) - rrt(410) * density(14) 
  pd(08,14) = pd(08,14) + rrt(411) * density(26) 
  pd(08,26) = pd(08,26) + rrt(411) * density(14) 
  pd(14,14) = pd(14,14) - rrt(411) * density(26) 
  pd(14,26) = pd(14,26) - rrt(411) * density(14) 
  pd(21,14) = pd(21,14) + rrt(411) * density(26) 
  pd(21,26) = pd(21,26) + rrt(411) * density(14) 
  pd(26,14) = pd(26,14) - rrt(411) * density(26) 
  pd(26,26) = pd(26,26) - rrt(411) * density(14) 
  pd(08,14) = pd(08,14) + rrt(412) * density(34) 
  pd(08,34) = pd(08,34) + rrt(412) * density(14) 
  pd(14,14) = pd(14,14) - rrt(412) * density(34) 
  pd(14,34) = pd(14,34) - rrt(412) * density(14) 
  pd(23,14) = pd(23,14) + rrt(412) * density(34) 
  pd(23,34) = pd(23,34) + rrt(412) * density(14) 
  pd(34,14) = pd(34,14) - rrt(412) * density(34) 
  pd(34,34) = pd(34,34) - rrt(412) * density(14) 
  pd(06,14) = pd(06,14) + rrt(413) * density(34) 
  pd(06,34) = pd(06,34) + rrt(413) * density(14) 
  pd(11,14) = pd(11,14) + rrt(413) * density(34) 
  pd(11,34) = pd(11,34) + rrt(413) * density(14) 
  pd(14,14) = pd(14,14) - rrt(413) * density(34) 
  pd(14,34) = pd(14,34) - rrt(413) * density(14) 
  pd(34,14) = pd(34,14) - rrt(413) * density(34) 
  pd(34,34) = pd(34,34) - rrt(413) * density(14) 
  pd(11,14) = pd(11,14) + rrt(414) * 2.0d0
  pd(14,14) = pd(14,14) - rrt(414) 
  pd(08,14) = pd(08,14) + rrt(415) * density(31) 
  pd(08,31) = pd(08,31) + rrt(415) * density(14) 
  pd(14,14) = pd(14,14) - rrt(415) * density(31) 
  pd(14,31) = pd(14,31) - rrt(415) * density(14) 
  pd(25,14) = pd(25,14) + rrt(415) * density(31) 
  pd(25,31) = pd(25,31) + rrt(415) * density(14) 
  pd(31,14) = pd(31,14) - rrt(415) * density(31) 
  pd(31,31) = pd(31,31) - rrt(415) * density(14) 
  pd(11,14) = pd(11,14) + rrt(416) * density(43) 
  pd(11,43) = pd(11,43) + rrt(416) * density(14) 
  pd(14,14) = pd(14,14) - rrt(416) * density(43) 
  pd(14,43) = pd(14,43) - rrt(416) * density(14) 
  pd(30,14) = pd(30,14) + rrt(416) * density(43) 
  pd(30,43) = pd(30,43) + rrt(416) * density(14) 
  pd(43,14) = pd(43,14) - rrt(416) * density(43) 
  pd(43,43) = pd(43,43) - rrt(416) * density(14) 
  pd(08,14) = pd(08,14) + rrt(417) * density(47) 
  pd(08,47) = pd(08,47) + rrt(417) * density(14) 
  pd(11,14) = pd(11,14) + rrt(417) * density(47) 
  pd(11,47) = pd(11,47) + rrt(417) * density(14) 
  pd(14,14) = pd(14,14) - rrt(417) * density(47) 
  pd(14,47) = pd(14,47) - rrt(417) * density(14) 
  pd(47,14) = pd(47,14) - rrt(417) * density(47) 
  pd(47,47) = pd(47,47) - rrt(417) * density(14) 
  pd(08,08) = pd(08,08) - rrt(418) * density(08) * 4.0d0
  pd(14,08) = pd(14,08) + rrt(418) * density(08) * 2.0d0
  pd(30,08) = pd(30,08) + rrt(418) * density(08) * 2.0d0
  pd(08,08) = pd(08,08) - rrt(419) * density(30) 
  pd(08,30) = pd(08,30) - rrt(419) * density(08) 
  pd(14,08) = pd(14,08) + rrt(419) * density(30) 
  pd(14,30) = pd(14,30) + rrt(419) * density(08) 
  pd(30,08) = pd(30,08) - rrt(419) * density(30) 
  pd(30,30) = pd(30,30) - rrt(419) * density(08) 
  pd(42,08) = pd(42,08) + rrt(419) * density(30) 
  pd(42,30) = pd(42,30) + rrt(419) * density(08) 
  pd(08,08) = pd(08,08) - rrt(420) * density(32) 
  pd(08,32) = pd(08,32) - rrt(420) * density(08) 
  pd(14,08) = pd(14,08) + rrt(420) * density(32) 
  pd(14,32) = pd(14,32) + rrt(420) * density(08) 
  pd(32,08) = pd(32,08) - rrt(420) * density(32) 
  pd(32,32) = pd(32,32) - rrt(420) * density(08) 
  pd(33,08) = pd(33,08) + rrt(420) * density(32) 
  pd(33,32) = pd(33,32) + rrt(420) * density(08) 
  pd(08,08) = pd(08,08) - rrt(421) * density(33) 
  pd(08,33) = pd(08,33) - rrt(421) * density(08) 
  pd(30,08) = pd(30,08) + rrt(421) * density(33) 
  pd(30,33) = pd(30,33) + rrt(421) * density(08) 
  pd(32,08) = pd(32,08) + rrt(421) * density(33) 
  pd(32,33) = pd(32,33) + rrt(421) * density(08) 
  pd(33,08) = pd(33,08) - rrt(421) * density(33) 
  pd(33,33) = pd(33,33) - rrt(421) * density(08) 
  pd(07,07) = pd(07,07) - rrt(422) * density(08) 
  pd(07,08) = pd(07,08) - rrt(422) * density(07) 
  pd(08,07) = pd(08,07) - rrt(422) * density(08) 
  pd(08,08) = pd(08,08) - rrt(422) * density(07) 
  pd(14,07) = pd(14,07) + rrt(422) * density(08) 
  pd(14,08) = pd(14,08) + rrt(422) * density(07) 
  pd(44,07) = pd(44,07) + rrt(422) * density(08) 
  pd(44,08) = pd(44,08) + rrt(422) * density(07) 
  pd(07,08) = pd(07,08) + rrt(423) * density(44) 
  pd(07,44) = pd(07,44) + rrt(423) * density(08) 
  pd(08,08) = pd(08,08) - rrt(423) * density(44) 
  pd(08,44) = pd(08,44) - rrt(423) * density(08) 
  pd(30,08) = pd(30,08) + rrt(423) * density(44) 
  pd(30,44) = pd(30,44) + rrt(423) * density(08) 
  pd(44,08) = pd(44,08) - rrt(423) * density(44) 
  pd(44,44) = pd(44,44) - rrt(423) * density(08) 
  pd(08,08) = pd(08,08) - rrt(424) * density(44) 
  pd(08,44) = pd(08,44) - rrt(424) * density(08) 
  pd(14,08) = pd(14,08) + rrt(424) * density(44) 
  pd(14,44) = pd(14,44) + rrt(424) * density(08) 
  pd(21,08) = pd(21,08) + rrt(424) * density(44) 
  pd(21,44) = pd(21,44) + rrt(424) * density(08) 
  pd(44,08) = pd(44,08) - rrt(424) * density(44) 
  pd(44,44) = pd(44,44) - rrt(424) * density(08) 
  pd(08,08) = pd(08,08) - rrt(425) * density(21) 
  pd(08,21) = pd(08,21) - rrt(425) * density(08) 
  pd(14,08) = pd(14,08) + rrt(425) * density(21) 
  pd(14,21) = pd(14,21) + rrt(425) * density(08) 
  pd(21,08) = pd(21,08) - rrt(425) * density(21) 
  pd(21,21) = pd(21,21) - rrt(425) * density(08) 
  pd(26,08) = pd(26,08) + rrt(425) * density(21) 
  pd(26,21) = pd(26,21) + rrt(425) * density(08) 
  pd(08,08) = pd(08,08) - rrt(426) * density(23) 
  pd(08,23) = pd(08,23) - rrt(426) * density(08) 
  pd(14,08) = pd(14,08) + rrt(426) * density(23) 
  pd(14,23) = pd(14,23) + rrt(426) * density(08) 
  pd(23,08) = pd(23,08) - rrt(426) * density(23) 
  pd(23,23) = pd(23,23) - rrt(426) * density(08) 
  pd(34,08) = pd(34,08) + rrt(426) * density(23) 
  pd(34,23) = pd(34,23) + rrt(426) * density(08) 
  pd(08,08) = pd(08,08) - rrt(427) * density(34) 
  pd(08,34) = pd(08,34) - rrt(427) * density(08) 
  pd(11,08) = pd(11,08) + rrt(427) * density(34) * 2.0d0
  pd(11,34) = pd(11,34) + rrt(427) * density(08) * 2.0d0
  pd(34,08) = pd(34,08) - rrt(427) * density(34) 
  pd(34,34) = pd(34,34) - rrt(427) * density(08) 
  pd(08,08) = pd(08,08) - rrt(428) * density(34) 
  pd(08,34) = pd(08,34) - rrt(428) * density(08) 
  pd(23,08) = pd(23,08) + rrt(428) * density(34) 
  pd(23,34) = pd(23,34) + rrt(428) * density(08) 
  pd(30,08) = pd(30,08) + rrt(428) * density(34) 
  pd(30,34) = pd(30,34) + rrt(428) * density(08) 
  pd(34,08) = pd(34,08) - rrt(428) * density(34) 
  pd(34,34) = pd(34,34) - rrt(428) * density(08) 
  pd(08,08) = pd(08,08) - rrt(429) * density(34) 
  pd(08,34) = pd(08,34) - rrt(429) * density(08) 
  pd(14,08) = pd(14,08) + rrt(429) * density(34) 
  pd(14,34) = pd(14,34) + rrt(429) * density(08) 
  pd(34,08) = pd(34,08) - rrt(429) * density(34) 
  pd(34,34) = pd(34,34) - rrt(429) * density(08) 
  pd(08,08) = pd(08,08) - rrt(430) 
  pd(30,08) = pd(30,08) + rrt(430) 
  pd(34,08) = pd(34,08) + rrt(430) 
  pd(08,08) = pd(08,08) - rrt(431) * density(08) * 4.0d0
  pd(25,08) = pd(25,08) + rrt(431) * density(08) * 2.0d0
  pd(08,08) = pd(08,08) - rrt(432) * density(31) 
  pd(08,31) = pd(08,31) - rrt(432) * density(08) 
  pd(25,08) = pd(25,08) + rrt(432) * density(31) 
  pd(25,31) = pd(25,31) + rrt(432) * density(08) 
  pd(30,08) = pd(30,08) + rrt(432) * density(31) 
  pd(30,31) = pd(30,31) + rrt(432) * density(08) 
  pd(31,08) = pd(31,08) - rrt(432) * density(31) 
  pd(31,31) = pd(31,31) - rrt(432) * density(08) 
  pd(08,08) = pd(08,08) - rrt(433) * density(42) 
  pd(08,42) = pd(08,42) - rrt(433) * density(08) 
  pd(30,08) = pd(30,08) + rrt(433) * density(42) * 2.0d0
  pd(30,42) = pd(30,42) + rrt(433) * density(08) * 2.0d0
  pd(42,08) = pd(42,08) - rrt(433) * density(42) 
  pd(42,42) = pd(42,42) - rrt(433) * density(08) 
  pd(23,30) = pd(23,30) + rrt(434) * density(34) 
  pd(23,34) = pd(23,34) + rrt(434) * density(30) 
  pd(30,30) = pd(30,30) - rrt(434) * density(34) 
  pd(30,34) = pd(30,34) - rrt(434) * density(30) 
  pd(34,30) = pd(34,30) - rrt(434) * density(34) 
  pd(34,34) = pd(34,34) - rrt(434) * density(30) 
  pd(42,30) = pd(42,30) + rrt(434) * density(34) 
  pd(42,34) = pd(42,34) + rrt(434) * density(30) 
  pd(08,30) = pd(08,30) + rrt(435) * density(34) 
  pd(08,34) = pd(08,34) + rrt(435) * density(30) 
  pd(30,30) = pd(30,30) - rrt(435) * density(34) 
  pd(30,34) = pd(30,34) - rrt(435) * density(30) 
  pd(34,30) = pd(34,30) - rrt(435) * density(34) 
  pd(34,34) = pd(34,34) - rrt(435) * density(30) 
  pd(08,23) = pd(08,23) + rrt(436) * density(30) 
  pd(08,30) = pd(08,30) + rrt(436) * density(23) 
  pd(23,23) = pd(23,23) - rrt(436) * density(30) 
  pd(23,30) = pd(23,30) - rrt(436) * density(23) 
  pd(30,23) = pd(30,23) - rrt(436) * density(30) 
  pd(30,30) = pd(30,30) - rrt(436) * density(23) 
  pd(34,23) = pd(34,23) + rrt(436) * density(30) 
  pd(34,30) = pd(34,30) + rrt(436) * density(23) 
  pd(30,30) = pd(30,30) - rrt(437) 
  pd(34,30) = pd(34,30) + rrt(437) 
  pd(42,30) = pd(42,30) + rrt(437) 
  pd(08,21) = pd(08,21) + rrt(438) * density(30) 
  pd(08,30) = pd(08,30) + rrt(438) * density(21) 
  pd(21,21) = pd(21,21) - rrt(438) * density(30) 
  pd(21,30) = pd(21,30) - rrt(438) * density(21) 
  pd(26,21) = pd(26,21) + rrt(438) * density(30) 
  pd(26,30) = pd(26,30) + rrt(438) * density(21) 
  pd(30,21) = pd(30,21) - rrt(438) * density(30) 
  pd(30,30) = pd(30,30) - rrt(438) * density(21) 
  pd(30,30) = pd(30,30) - rrt(439) * density(32) 
  pd(30,32) = pd(30,32) - rrt(439) * density(30) 
  pd(32,30) = pd(32,30) - rrt(439) * density(32) 
  pd(32,32) = pd(32,32) - rrt(439) * density(30) 
  pd(42,30) = pd(42,30) + rrt(439) * density(32) * 2.0d0
  pd(42,32) = pd(42,32) + rrt(439) * density(30) * 2.0d0
  pd(21,21) = pd(21,21) - rrt(440) * density(30) 
  pd(21,30) = pd(21,30) - rrt(440) * density(21) 
  pd(30,21) = pd(30,21) - rrt(440) * density(30) 
  pd(30,30) = pd(30,30) - rrt(440) * density(21) 
  pd(42,21) = pd(42,21) + rrt(440) * density(30) 
  pd(42,30) = pd(42,30) + rrt(440) * density(21) 
  pd(44,21) = pd(44,21) + rrt(440) * density(30) 
  pd(44,30) = pd(44,30) + rrt(440) * density(21) 
  pd(08,30) = pd(08,30) + rrt(441) * density(30) * 2.0d0
  pd(30,30) = pd(30,30) - rrt(441) * density(30) * 4.0d0
  pd(42,30) = pd(42,30) + rrt(441) * density(30) * 2.0d0
  pd(11,11) = pd(11,11) - rrt(442) * density(30) 
  pd(11,30) = pd(11,30) - rrt(442) * density(11) 
  pd(30,11) = pd(30,11) - rrt(442) * density(30) 
  pd(30,30) = pd(30,30) - rrt(442) * density(11) 
  pd(44,11) = pd(44,11) + rrt(442) * density(30) 
  pd(44,30) = pd(44,30) + rrt(442) * density(11) 
  pd(23,30) = pd(23,30) + rrt(443) 
  pd(30,30) = pd(30,30) - rrt(443) 
  pd(32,30) = pd(32,30) + rrt(443) 
  pd(08,08) = pd(08,08) - rrt(444) * density(30) 
  pd(08,30) = pd(08,30) - rrt(444) * density(08) 
  pd(30,08) = pd(30,08) - rrt(444) * density(30) 
  pd(30,30) = pd(30,30) - rrt(444) * density(08) 
  pd(31,08) = pd(31,08) + rrt(444) * density(30) 
  pd(31,30) = pd(31,30) + rrt(444) * density(08) 
  pd(14,23) = pd(14,23) + rrt(445) * density(30) 
  pd(14,30) = pd(14,30) + rrt(445) * density(23) 
  pd(23,23) = pd(23,23) - rrt(445) * density(30) 
  pd(23,30) = pd(23,30) - rrt(445) * density(23) 
  pd(30,23) = pd(30,23) - rrt(445) * density(30) 
  pd(30,30) = pd(30,30) - rrt(445) * density(23) 
  pd(21,30) = pd(21,30) + rrt(446) * density(47) 
  pd(21,47) = pd(21,47) + rrt(446) * density(30) 
  pd(30,30) = pd(30,30) - rrt(446) * density(47) 
  pd(30,47) = pd(30,47) - rrt(446) * density(30) 
  pd(47,30) = pd(47,30) - rrt(446) * density(47) 
  pd(47,47) = pd(47,47) - rrt(446) * density(30) 
  pd(21,30) = pd(21,30) + rrt(447) * density(31) 
  pd(21,31) = pd(21,31) + rrt(447) * density(30) 
  pd(30,30) = pd(30,30) - rrt(447) * density(31) 
  pd(30,31) = pd(30,31) - rrt(447) * density(30) 
  pd(31,30) = pd(31,30) - rrt(447) * density(31) 
  pd(31,31) = pd(31,31) - rrt(447) * density(30) 
  pd(44,30) = pd(44,30) + rrt(447) * density(31) 
  pd(44,31) = pd(44,31) + rrt(447) * density(30) 
  pd(30,42) = pd(30,42) + rrt(448) * density(42) * 2.0d0
  pd(32,42) = pd(32,42) + rrt(448) * density(42) * 2.0d0
  pd(42,42) = pd(42,42) - rrt(448) * density(42) * 4.0d0
  pd(07,07) = pd(07,07) - rrt(449) * density(42) 
  pd(07,42) = pd(07,42) - rrt(449) * density(07) 
  pd(30,07) = pd(30,07) + rrt(449) * density(42) 
  pd(30,42) = pd(30,42) + rrt(449) * density(07) 
  pd(42,07) = pd(42,07) - rrt(449) * density(42) 
  pd(42,42) = pd(42,42) - rrt(449) * density(07) 
  pd(44,07) = pd(44,07) + rrt(449) * density(42) 
  pd(44,42) = pd(44,42) + rrt(449) * density(07) 
  pd(07,42) = pd(07,42) + rrt(450) * density(44) 
  pd(07,44) = pd(07,44) + rrt(450) * density(42) 
  pd(32,42) = pd(32,42) + rrt(450) * density(44) 
  pd(32,44) = pd(32,44) + rrt(450) * density(42) 
  pd(42,42) = pd(42,42) - rrt(450) * density(44) 
  pd(42,44) = pd(42,44) - rrt(450) * density(42) 
  pd(44,42) = pd(44,42) - rrt(450) * density(44) 
  pd(44,44) = pd(44,44) - rrt(450) * density(42) 
  pd(21,42) = pd(21,42) + rrt(451) * density(44) 
  pd(21,44) = pd(21,44) + rrt(451) * density(42) 
  pd(30,42) = pd(30,42) + rrt(451) * density(44) 
  pd(30,44) = pd(30,44) + rrt(451) * density(42) 
  pd(42,42) = pd(42,42) - rrt(451) * density(44) 
  pd(42,44) = pd(42,44) - rrt(451) * density(42) 
  pd(44,42) = pd(44,42) - rrt(451) * density(44) 
  pd(44,44) = pd(44,44) - rrt(451) * density(42) 
  pd(21,21) = pd(21,21) - rrt(452) * density(42) 
  pd(21,42) = pd(21,42) - rrt(452) * density(21) 
  pd(26,21) = pd(26,21) + rrt(452) * density(42) 
  pd(26,42) = pd(26,42) + rrt(452) * density(21) 
  pd(30,21) = pd(30,21) + rrt(452) * density(42) 
  pd(30,42) = pd(30,42) + rrt(452) * density(21) 
  pd(42,21) = pd(42,21) - rrt(452) * density(42) 
  pd(42,42) = pd(42,42) - rrt(452) * density(21) 
  pd(21,26) = pd(21,26) + rrt(453) * density(42) 
  pd(21,42) = pd(21,42) + rrt(453) * density(26) 
  pd(26,26) = pd(26,26) - rrt(453) * density(42) 
  pd(26,42) = pd(26,42) - rrt(453) * density(26) 
  pd(32,26) = pd(32,26) + rrt(453) * density(42) 
  pd(32,42) = pd(32,42) + rrt(453) * density(26) 
  pd(42,26) = pd(42,26) - rrt(453) * density(42) 
  pd(42,42) = pd(42,42) - rrt(453) * density(26) 
  pd(23,23) = pd(23,23) - rrt(454) * density(42) 
  pd(23,42) = pd(23,42) - rrt(454) * density(23) 
  pd(30,23) = pd(30,23) + rrt(454) * density(42) 
  pd(30,42) = pd(30,42) + rrt(454) * density(23) 
  pd(34,23) = pd(34,23) + rrt(454) * density(42) 
  pd(34,42) = pd(34,42) + rrt(454) * density(23) 
  pd(42,23) = pd(42,23) - rrt(454) * density(42) 
  pd(42,42) = pd(42,42) - rrt(454) * density(23) 
  pd(23,34) = pd(23,34) + rrt(455) * density(42) 
  pd(23,42) = pd(23,42) + rrt(455) * density(34) 
  pd(32,34) = pd(32,34) + rrt(455) * density(42) 
  pd(32,42) = pd(32,42) + rrt(455) * density(34) 
  pd(34,34) = pd(34,34) - rrt(455) * density(42) 
  pd(34,42) = pd(34,42) - rrt(455) * density(34) 
  pd(42,34) = pd(42,34) - rrt(455) * density(42) 
  pd(42,42) = pd(42,42) - rrt(455) * density(34) 
  pd(30,34) = pd(30,34) + rrt(456) * density(42) 
  pd(30,42) = pd(30,42) + rrt(456) * density(34) 
  pd(34,34) = pd(34,34) - rrt(456) * density(42) 
  pd(34,42) = pd(34,42) - rrt(456) * density(34) 
  pd(42,34) = pd(42,34) - rrt(456) * density(42) 
  pd(42,42) = pd(42,42) - rrt(456) * density(34) 
  pd(32,42) = pd(32,42) + rrt(457) 
  pd(34,42) = pd(34,42) + rrt(457) 
  pd(42,42) = pd(42,42) - rrt(457) 
  pd(25,31) = pd(25,31) + rrt(458) * density(42) 
  pd(25,42) = pd(25,42) + rrt(458) * density(31) 
  pd(31,31) = pd(31,31) - rrt(458) * density(42) 
  pd(31,42) = pd(31,42) - rrt(458) * density(31) 
  pd(32,31) = pd(32,31) + rrt(458) * density(42) 
  pd(32,42) = pd(32,42) + rrt(458) * density(31) 
  pd(42,31) = pd(42,31) - rrt(458) * density(42) 
  pd(42,42) = pd(42,42) - rrt(458) * density(31) 
  pd(32,32) = pd(32,32) - rrt(459) * density(34) 
  pd(32,34) = pd(32,34) - rrt(459) * density(32) 
  pd(34,32) = pd(34,32) - rrt(459) * density(34) 
  pd(34,34) = pd(34,34) - rrt(459) * density(32) 
  pd(42,32) = pd(42,32) + rrt(459) * density(34) 
  pd(42,34) = pd(42,34) + rrt(459) * density(32) 
  pd(23,23) = pd(23,23) - rrt(460) * density(32) 
  pd(23,32) = pd(23,32) - rrt(460) * density(23) 
  pd(30,23) = pd(30,23) + rrt(460) * density(32) 
  pd(30,32) = pd(30,32) + rrt(460) * density(23) 
  pd(32,23) = pd(32,23) - rrt(460) * density(32) 
  pd(32,32) = pd(32,32) - rrt(460) * density(23) 
  pd(23,23) = pd(23,23) - rrt(461) * density(32) 
  pd(23,32) = pd(23,32) - rrt(461) * density(23) 
  pd(32,23) = pd(32,23) - rrt(461) * density(32) 
  pd(32,32) = pd(32,32) - rrt(461) * density(23) 
  pd(34,23) = pd(34,23) + rrt(461) * density(32) 
  pd(34,32) = pd(34,32) + rrt(461) * density(23) 
  pd(42,23) = pd(42,23) + rrt(461) * density(32) 
  pd(42,32) = pd(42,32) + rrt(461) * density(23) 
  pd(11,11) = pd(11,11) - rrt(462) * density(32) 
  pd(11,32) = pd(11,32) - rrt(462) * density(11) 
  pd(26,11) = pd(26,11) + rrt(462) * density(32) 
  pd(26,32) = pd(26,32) + rrt(462) * density(11) 
  pd(32,11) = pd(32,11) - rrt(462) * density(32) 
  pd(32,32) = pd(32,32) - rrt(462) * density(11) 
  pd(21,31) = pd(21,31) + rrt(463) * density(32) 
  pd(21,32) = pd(21,32) + rrt(463) * density(31) 
  pd(26,31) = pd(26,31) + rrt(463) * density(32) 
  pd(26,32) = pd(26,32) + rrt(463) * density(31) 
  pd(31,31) = pd(31,31) - rrt(463) * density(32) 
  pd(31,32) = pd(31,32) - rrt(463) * density(31) 
  pd(32,31) = pd(32,31) - rrt(463) * density(32) 
  pd(32,32) = pd(32,32) - rrt(463) * density(31) 
  pd(07,07) = pd(07,07) - rrt(464) * density(26) 
  pd(07,26) = pd(07,26) - rrt(464) * density(07) 
  pd(21,07) = pd(21,07) + rrt(464) * density(26) 
  pd(21,26) = pd(21,26) + rrt(464) * density(07) 
  pd(26,07) = pd(26,07) - rrt(464) * density(26) 
  pd(26,26) = pd(26,26) - rrt(464) * density(07) 
  pd(44,07) = pd(44,07) + rrt(464) * density(26) 
  pd(44,26) = pd(44,26) + rrt(464) * density(07) 
  pd(07,07) = pd(07,07) - rrt(465) * density(34) 
  pd(07,34) = pd(07,34) - rrt(465) * density(07) 
  pd(23,07) = pd(23,07) + rrt(465) * density(34) 
  pd(23,34) = pd(23,34) + rrt(465) * density(07) 
  pd(34,07) = pd(34,07) - rrt(465) * density(34) 
  pd(34,34) = pd(34,34) - rrt(465) * density(07) 
  pd(44,07) = pd(44,07) + rrt(465) * density(34) 
  pd(44,34) = pd(44,34) + rrt(465) * density(07) 
  pd(07,07) = pd(07,07) - rrt(466) 
  pd(08,07) = pd(08,07) + rrt(466) 
  pd(11,07) = pd(11,07) + rrt(466) 
  pd(07,07) = pd(07,07) - rrt(467) * density(31) 
  pd(07,31) = pd(07,31) - rrt(467) * density(07) 
  pd(25,07) = pd(25,07) + rrt(467) * density(31) 
  pd(25,31) = pd(25,31) + rrt(467) * density(07) 
  pd(31,07) = pd(31,07) - rrt(467) * density(31) 
  pd(31,31) = pd(31,31) - rrt(467) * density(07) 
  pd(44,07) = pd(44,07) + rrt(467) * density(31) 
  pd(44,31) = pd(44,31) + rrt(467) * density(07) 
  pd(07,07) = pd(07,07) - rrt(468) * density(47) 
  pd(07,47) = pd(07,47) - rrt(468) * density(07) 
  pd(25,07) = pd(25,07) + rrt(468) * density(47) 
  pd(25,47) = pd(25,47) + rrt(468) * density(07) 
  pd(47,07) = pd(47,07) - rrt(468) * density(47) 
  pd(47,47) = pd(47,47) - rrt(468) * density(07) 
  pd(07,44) = pd(07,44) + rrt(469) * density(44) * 2.0d0
  pd(21,44) = pd(21,44) + rrt(469) * density(44) * 2.0d0
  pd(44,44) = pd(44,44) - rrt(469) * density(44) * 4.0d0
  pd(07,21) = pd(07,21) + rrt(470) * density(44) 
  pd(07,44) = pd(07,44) + rrt(470) * density(21) 
  pd(21,21) = pd(21,21) - rrt(470) * density(44) 
  pd(21,44) = pd(21,44) - rrt(470) * density(21) 
  pd(26,21) = pd(26,21) + rrt(470) * density(44) 
  pd(26,44) = pd(26,44) + rrt(470) * density(21) 
  pd(44,21) = pd(44,21) - rrt(470) * density(44) 
  pd(44,44) = pd(44,44) - rrt(470) * density(21) 
  pd(21,26) = pd(21,26) + rrt(471) * density(44) * 2.0d0
  pd(21,44) = pd(21,44) + rrt(471) * density(26) * 2.0d0
  pd(26,26) = pd(26,26) - rrt(471) * density(44) 
  pd(26,44) = pd(26,44) - rrt(471) * density(26) 
  pd(44,26) = pd(44,26) - rrt(471) * density(44) 
  pd(44,44) = pd(44,44) - rrt(471) * density(26) 
  pd(07,23) = pd(07,23) + rrt(472) * density(44) 
  pd(07,44) = pd(07,44) + rrt(472) * density(23) 
  pd(23,23) = pd(23,23) - rrt(472) * density(44) 
  pd(23,44) = pd(23,44) - rrt(472) * density(23) 
  pd(34,23) = pd(34,23) + rrt(472) * density(44) 
  pd(34,44) = pd(34,44) + rrt(472) * density(23) 
  pd(44,23) = pd(44,23) - rrt(472) * density(44) 
  pd(44,44) = pd(44,44) - rrt(472) * density(23) 
  pd(21,34) = pd(21,34) + rrt(473) * density(44) 
  pd(21,44) = pd(21,44) + rrt(473) * density(34) 
  pd(23,34) = pd(23,34) + rrt(473) * density(44) 
  pd(23,44) = pd(23,44) + rrt(473) * density(34) 
  pd(34,34) = pd(34,34) - rrt(473) * density(44) 
  pd(34,44) = pd(34,44) - rrt(473) * density(34) 
  pd(44,34) = pd(44,34) - rrt(473) * density(44) 
  pd(44,44) = pd(44,44) - rrt(473) * density(34) 
  pd(07,34) = pd(07,34) + rrt(474) * density(44) 
  pd(07,44) = pd(07,44) + rrt(474) * density(34) 
  pd(34,34) = pd(34,34) - rrt(474) * density(44) 
  pd(34,44) = pd(34,44) - rrt(474) * density(34) 
  pd(44,34) = pd(44,34) - rrt(474) * density(44) 
  pd(44,44) = pd(44,44) - rrt(474) * density(34) 
  pd(08,34) = pd(08,34) + rrt(475) * density(44) 
  pd(08,44) = pd(08,44) + rrt(475) * density(34) 
  pd(11,34) = pd(11,34) + rrt(475) * density(44) 
  pd(11,44) = pd(11,44) + rrt(475) * density(34) 
  pd(34,34) = pd(34,34) - rrt(475) * density(44) 
  pd(34,44) = pd(34,44) - rrt(475) * density(34) 
  pd(44,34) = pd(44,34) - rrt(475) * density(44) 
  pd(44,44) = pd(44,44) - rrt(475) * density(34) 
  pd(21,44) = pd(21,44) + rrt(476) 
  pd(34,44) = pd(34,44) + rrt(476) 
  pd(44,44) = pd(44,44) - rrt(476) 
  pd(11,44) = pd(11,44) + rrt(477) 
  pd(30,44) = pd(30,44) + rrt(477) 
  pd(44,44) = pd(44,44) - rrt(477) 
  pd(21,31) = pd(21,31) + rrt(478) * density(44) 
  pd(21,44) = pd(21,44) + rrt(478) * density(31) 
  pd(25,31) = pd(25,31) + rrt(478) * density(44) 
  pd(25,44) = pd(25,44) + rrt(478) * density(31) 
  pd(31,31) = pd(31,31) - rrt(478) * density(44) 
  pd(31,44) = pd(31,44) - rrt(478) * density(31) 
  pd(44,31) = pd(44,31) - rrt(478) * density(44) 
  pd(44,44) = pd(44,44) - rrt(478) * density(31) 
  pd(21,21) = pd(21,21) - rrt(479) * density(32) 
  pd(21,32) = pd(21,32) - rrt(479) * density(21) 
  pd(26,21) = pd(26,21) + rrt(479) * density(32) 
  pd(26,32) = pd(26,32) + rrt(479) * density(21) 
  pd(32,21) = pd(32,21) - rrt(479) * density(32) 
  pd(32,32) = pd(32,32) - rrt(479) * density(21) 
  pd(42,21) = pd(42,21) + rrt(479) * density(32) 
  pd(42,32) = pd(42,32) + rrt(479) * density(21) 
  pd(21,21) = pd(21,21) - rrt(480) * density(21) * 4.0d0
  pd(26,21) = pd(26,21) + rrt(480) * density(21) * 2.0d0
  pd(44,21) = pd(44,21) + rrt(480) * density(21) * 2.0d0
  pd(21,21) = pd(21,21) - rrt(481) 
  pd(26,21) = pd(26,21) + rrt(481) 
  pd(34,21) = pd(34,21) + rrt(481) 
  pd(21,21) = pd(21,21) - rrt(482) * density(34) 
  pd(21,34) = pd(21,34) - rrt(482) * density(21) 
  pd(23,21) = pd(23,21) + rrt(482) * density(34) 
  pd(23,34) = pd(23,34) + rrt(482) * density(21) 
  pd(26,21) = pd(26,21) + rrt(482) * density(34) 
  pd(26,34) = pd(26,34) + rrt(482) * density(21) 
  pd(34,21) = pd(34,21) - rrt(482) * density(34) 
  pd(34,34) = pd(34,34) - rrt(482) * density(21) 
  pd(21,21) = pd(21,21) - rrt(483) * density(34) 
  pd(21,34) = pd(21,34) - rrt(483) * density(21) 
  pd(34,21) = pd(34,21) - rrt(483) * density(34) 
  pd(34,34) = pd(34,34) - rrt(483) * density(21) 
  pd(44,21) = pd(44,21) + rrt(483) * density(34) 
  pd(44,34) = pd(44,34) + rrt(483) * density(21) 
  pd(11,21) = pd(11,21) + rrt(484) 
  pd(21,21) = pd(21,21) - rrt(484) 
  pd(42,21) = pd(42,21) + rrt(484) 
  pd(11,11) = pd(11,11) - rrt(485) * density(21) 
  pd(11,21) = pd(11,21) - rrt(485) * density(11) 
  pd(21,11) = pd(21,11) - rrt(485) * density(21) 
  pd(21,21) = pd(21,21) - rrt(485) * density(11) 
  pd(31,11) = pd(31,11) + rrt(485) * density(21) 
  pd(31,21) = pd(31,21) + rrt(485) * density(11) 
  pd(21,21) = pd(21,21) - rrt(486) * density(31) 
  pd(21,31) = pd(21,31) - rrt(486) * density(21) 
  pd(25,21) = pd(25,21) + rrt(486) * density(31) 
  pd(25,31) = pd(25,31) + rrt(486) * density(21) 
  pd(26,21) = pd(26,21) + rrt(486) * density(31) 
  pd(26,31) = pd(26,31) + rrt(486) * density(21) 
  pd(31,21) = pd(31,21) - rrt(486) * density(31) 
  pd(31,31) = pd(31,31) - rrt(486) * density(21) 
  pd(21,23) = pd(21,23) + rrt(487) * density(26) 
  pd(21,26) = pd(21,26) + rrt(487) * density(23) 
  pd(23,23) = pd(23,23) - rrt(487) * density(26) 
  pd(23,26) = pd(23,26) - rrt(487) * density(23) 
  pd(26,23) = pd(26,23) - rrt(487) * density(26) 
  pd(26,26) = pd(26,26) - rrt(487) * density(23) 
  pd(34,23) = pd(34,23) + rrt(487) * density(26) 
  pd(34,26) = pd(34,26) + rrt(487) * density(23) 
  pd(21,26) = pd(21,26) + rrt(488) * density(34) 
  pd(21,34) = pd(21,34) + rrt(488) * density(26) 
  pd(26,26) = pd(26,26) - rrt(488) * density(34) 
  pd(26,34) = pd(26,34) - rrt(488) * density(26) 
  pd(34,26) = pd(34,26) - rrt(488) * density(34) 
  pd(34,34) = pd(34,34) - rrt(488) * density(26) 
  pd(11,26) = pd(11,26) + rrt(489) 
  pd(26,26) = pd(26,26) - rrt(489) 
  pd(32,26) = pd(32,26) + rrt(489) 
  pd(08,31) = pd(08,31) + rrt(490) 
  pd(30,31) = pd(30,31) + rrt(490) 
  pd(31,31) = pd(31,31) - rrt(490) 
  pd(30,31) = pd(30,31) + rrt(491) * density(47) 
  pd(30,47) = pd(30,47) + rrt(491) * density(31) 
  pd(31,31) = pd(31,31) - rrt(491) * density(47) 
  pd(31,47) = pd(31,47) - rrt(491) * density(31) 
  pd(44,31) = pd(44,31) + rrt(491) * density(47) 
  pd(44,47) = pd(44,47) + rrt(491) * density(31) 
  pd(47,31) = pd(47,31) - rrt(491) * density(47) 
  pd(47,47) = pd(47,47) - rrt(491) * density(31) 
  pd(11,31) = pd(11,31) + rrt(492) 
  pd(21,31) = pd(21,31) + rrt(492) 
  pd(31,31) = pd(31,31) - rrt(492) 
  pd(23,23) = pd(23,23) - rrt(493) * density(31) 
  pd(23,31) = pd(23,31) - rrt(493) * density(23) 
  pd(25,23) = pd(25,23) + rrt(493) * density(31) 
  pd(25,31) = pd(25,31) + rrt(493) * density(23) 
  pd(31,23) = pd(31,23) - rrt(493) * density(31) 
  pd(31,31) = pd(31,31) - rrt(493) * density(23) 
  pd(34,23) = pd(34,23) + rrt(493) * density(31) 
  pd(34,31) = pd(34,31) + rrt(493) * density(23) 
  pd(06,11) = pd(06,11) + rrt(494) * density(25) 
  pd(06,25) = pd(06,25) + rrt(494) * density(11) 
  pd(11,11) = pd(11,11) - rrt(494) * density(25) 
  pd(11,25) = pd(11,25) - rrt(494) * density(11) 
  pd(25,11) = pd(25,11) - rrt(494) * density(25) 
  pd(25,25) = pd(25,25) - rrt(494) * density(11) 
  pd(31,11) = pd(31,11) + rrt(494) * density(25) 
  pd(31,25) = pd(31,25) + rrt(494) * density(11) 
  pd(11,25) = pd(11,25) + rrt(495) 
  pd(25,25) = pd(25,25) - rrt(495) 
  pd(44,25) = pd(44,25) + rrt(495) 
  pd(08,25) = pd(08,25) + rrt(496) * 2.0d0
  pd(25,25) = pd(25,25) - rrt(496) 
  pd(23,25) = pd(23,25) + rrt(497) * density(34) 
  pd(23,34) = pd(23,34) + rrt(497) * density(25) 
  pd(25,25) = pd(25,25) - rrt(497) * density(34) 
  pd(25,34) = pd(25,34) - rrt(497) * density(25) 
  pd(31,25) = pd(31,25) + rrt(497) * density(34) 
  pd(31,34) = pd(31,34) + rrt(497) * density(25) 
  pd(34,25) = pd(34,25) - rrt(497) * density(34) 
  pd(34,34) = pd(34,34) - rrt(497) * density(25) 
  pd(11,25) = pd(11,25) + rrt(498) * density(47) 
  pd(11,47) = pd(11,47) + rrt(498) * density(25) 
  pd(25,25) = pd(25,25) - rrt(498) * density(47) 
  pd(25,47) = pd(25,47) - rrt(498) * density(25) 
  pd(31,25) = pd(31,25) + rrt(498) * density(47) 
  pd(31,47) = pd(31,47) + rrt(498) * density(25) 
  pd(47,25) = pd(47,25) - rrt(498) * density(47) 
  pd(47,47) = pd(47,47) - rrt(498) * density(25) 
  pd(25,25) = pd(25,25) - rrt(499) * density(42) 
  pd(25,42) = pd(25,42) - rrt(499) * density(25) 
  pd(30,25) = pd(30,25) + rrt(499) * density(42) 
  pd(30,42) = pd(30,42) + rrt(499) * density(25) 
  pd(31,25) = pd(31,25) + rrt(499) * density(42) 
  pd(31,42) = pd(31,42) + rrt(499) * density(25) 
  pd(42,25) = pd(42,25) - rrt(499) * density(42) 
  pd(42,42) = pd(42,42) - rrt(499) * density(25) 
  pd(07,25) = pd(07,25) + rrt(500) * density(44) 
  pd(07,44) = pd(07,44) + rrt(500) * density(25) 
  pd(25,25) = pd(25,25) - rrt(500) * density(44) 
  pd(25,44) = pd(25,44) - rrt(500) * density(25) 
  pd(31,25) = pd(31,25) + rrt(500) * density(44) 
  pd(31,44) = pd(31,44) + rrt(500) * density(25) 
  pd(44,25) = pd(44,25) - rrt(500) * density(44) 
  pd(44,44) = pd(44,44) - rrt(500) * density(25) 
  pd(25,25) = pd(25,25) - rrt(501) * density(33) 
  pd(25,33) = pd(25,33) - rrt(501) * density(25) 
  pd(31,25) = pd(31,25) + rrt(501) * density(33) 
  pd(31,33) = pd(31,33) + rrt(501) * density(25) 
  pd(32,25) = pd(32,25) + rrt(501) * density(33) 
  pd(32,33) = pd(32,33) + rrt(501) * density(25) 
  pd(33,25) = pd(33,25) - rrt(501) * density(33) 
  pd(33,33) = pd(33,33) - rrt(501) * density(25) 
  pd(08,08) = pd(08,08) - rrt(502) * density(25) 
  pd(08,25) = pd(08,25) - rrt(502) * density(08) 
  pd(14,08) = pd(14,08) + rrt(502) * density(25) 
  pd(14,25) = pd(14,25) + rrt(502) * density(08) 
  pd(25,08) = pd(25,08) - rrt(502) * density(25) 
  pd(25,25) = pd(25,25) - rrt(502) * density(08) 
  pd(31,08) = pd(31,08) + rrt(502) * density(25) 
  pd(31,25) = pd(31,25) + rrt(502) * density(08) 
  pd(21,25) = pd(21,25) + rrt(503) * density(26) 
  pd(21,26) = pd(21,26) + rrt(503) * density(25) 
  pd(25,25) = pd(25,25) - rrt(503) * density(26) 
  pd(25,26) = pd(25,26) - rrt(503) * density(25) 
  pd(26,25) = pd(26,25) - rrt(503) * density(26) 
  pd(26,26) = pd(26,26) - rrt(503) * density(25) 
  pd(31,25) = pd(31,25) + rrt(503) * density(26) 
  pd(31,26) = pd(31,26) + rrt(503) * density(25) 
  pd(22,25) = pd(22,25) + rrt(504) * density(47) 
  pd(22,47) = pd(22,47) + rrt(504) * density(25) 
  pd(25,25) = pd(25,25) - rrt(504) * density(47) 
  pd(25,47) = pd(25,47) - rrt(504) * density(25) 
  pd(47,25) = pd(47,25) - rrt(504) * density(47) 
  pd(47,47) = pd(47,47) - rrt(504) * density(25) 
  pd(11,22) = pd(11,22) + rrt(505) 
  pd(22,22) = pd(22,22) - rrt(505) 
  pd(31,22) = pd(31,22) + rrt(505) 
  pd(23,23) = pd(23,23) - rrt(506) 
  pd(34,23) = pd(34,23) + rrt(506) * 2.0d0
  pd(23,34) = pd(23,34) + rrt(507) * density(34) * 2.0d0
  pd(34,34) = pd(34,34) - rrt(507) * density(34) * 4.0d0
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
  DOUBLE PRECISION, PARAMETER :: FEC = 1.000D0
  DOUBLE PRECISION, PARAMETER :: FIR = 1.000D0
  DOUBLE PRECISION, PARAMETER :: FID = 1.000D0
  DOUBLE PRECISION, PARAMETER :: FIN = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F0= 5.084D-3
  DOUBLE PRECISION, PARAMETER :: F1= 1.865D-3
  DOUBLE PRECISION, PARAMETER :: F2= 1.707D-3
  DOUBLE PRECISION, PARAMETER :: F3= 1.383D0
  DOUBLE PRECISION, PARAMETER :: F4= 1.229D-9
  DOUBLE PRECISION, PARAMETER :: F5= 1.027D-7
  DOUBLE PRECISION, PARAMETER :: F6= 7.604D-8
  DOUBLE PRECISION, PARAMETER :: F7= 3.930D-3
  DOUBLE PRECISION, PARAMETER :: F8= 2.155D-6
  DOUBLE PRECISION, PARAMETER :: F9= 9.574D-8
  DOUBLE PRECISION, PARAMETER :: F10= 4.808D-8
  DOUBLE PRECISION, PARAMETER :: F11= 2.175D-8
  DOUBLE PRECISION, PARAMETER :: F12= 1.234D-9
  DOUBLE PRECISION, PARAMETER :: F13= 4.040D-2
  DOUBLE PRECISION, PARAMETER :: F14= 3.939D-1
  DOUBLE PRECISION, PARAMETER :: F15= 8.343D-10
  DOUBLE PRECISION, PARAMETER :: F16= 1.743D-6
  DOUBLE PRECISION, PARAMETER :: F17= 2.374D-10
  DOUBLE PRECISION, PARAMETER :: F18= 5.917D-5
  DOUBLE PRECISION, PARAMETER :: F19= 9.485D-10
  DOUBLE PRECISION, PARAMETER :: F20= 4.882D-9
  DOUBLE PRECISION, PARAMETER :: F21= 2.781D-1
  DOUBLE PRECISION, PARAMETER :: F22= 2.316D-5
  DOUBLE PRECISION, PARAMETER :: F23= 7.587D-10
  DOUBLE PRECISION, PARAMETER :: F24= 3.471D-8
  DOUBLE PRECISION, PARAMETER :: F25= 2.444D-6
  DOUBLE PRECISION, PARAMETER :: F26= 3.181D-7
  DOUBLE PRECISION, PARAMETER :: F27= 1.497D-5
  DOUBLE PRECISION, PARAMETER :: F28= 9.834D-10
  DOUBLE PRECISION, PARAMETER :: F29= 4.028D-7
  DOUBLE PRECISION, PARAMETER :: F30= 1.930D-5
  DOUBLE PRECISION, PARAMETER :: F31= 6.000D-2
  DOUBLE PRECISION, PARAMETER :: F32= 1.880D-3
  DOUBLE PRECISION, PARAMETER :: F33= 6.340D-7
  DOUBLE PRECISION, PARAMETER :: F34= 8.628D-4
  DOUBLE PRECISION, PARAMETER :: F35= 1.187D-5
  DOUBLE PRECISION, PARAMETER :: F36= 5.169D-2
  DOUBLE PRECISION, PARAMETER :: F37= 3.012D-5
  DOUBLE PRECISION, PARAMETER :: F38= 1.165D-5
  DOUBLE PRECISION, PARAMETER :: F39= 4.931D-5
  DOUBLE PRECISION, PARAMETER :: F40= 1.390D-7
  DOUBLE PRECISION, PARAMETER :: F41= 5.836D-4
  DOUBLE PRECISION, PARAMETER :: F42= 2.884D-7
  DOUBLE PRECISION, PARAMETER :: F43= 8.384D-9
  DOUBLE PRECISION, PARAMETER :: F44= 1.411D-5
  DOUBLE PRECISION, PARAMETER :: F45= 1.259D1
  DOUBLE PRECISION, PARAMETER :: F46= 4.026D-10
  DOUBLE PRECISION, PARAMETER :: F47= 2.626D-9
  DOUBLE PRECISION, PARAMETER :: F48= 5.202D-1
  DOUBLE PRECISION, PARAMETER :: F49= 3.160D-5
  DOUBLE PRECISION, PARAMETER :: F50= 3.950D-4
  DOUBLE PRECISION, PARAMETER :: F51= 2.573D1
  DOUBLE PRECISION, PARAMETER :: F52= 1.396D0
  DOUBLE PRECISION, PARAMETER :: F53= 1.812D-5
  DOUBLE PRECISION, PARAMETER :: F54= 1.357D-8
  DOUBLE PRECISION, PARAMETER :: F55= 2.471D-7
  DOUBLE PRECISION, PARAMETER :: F56= 6.963D-9
  DOUBLE PRECISION, PARAMETER :: F57= 4.017D-9
  DOUBLE PRECISION, PARAMETER :: F58= 3.924D-10
  DOUBLE PRECISION, PARAMETER :: F59= 2.337D0
  DOUBLE PRECISION, PARAMETER :: F60= 5.724D-10
  DOUBLE PRECISION, PARAMETER :: F61= 2.485D1
  DOUBLE PRECISION, PARAMETER :: F62= 2.749D-9
  DOUBLE PRECISION, PARAMETER :: F63= 1.998D-8
  DOUBLE PRECISION, PARAMETER :: F64= 1.496D-4
  DOUBLE PRECISION, PARAMETER :: F65= 2.782D0
  DOUBLE PRECISION, PARAMETER :: F66= 9.322D-2
  DOUBLE PRECISION, PARAMETER :: F67= 8.397D-1
  DOUBLE PRECISION, PARAMETER :: F68= 1.642D1
  DOUBLE PRECISION, PARAMETER :: F69= 3.732D-4
  DOUBLE PRECISION, PARAMETER :: F70= 1.096D-10
  DOUBLE PRECISION, PARAMETER :: F71= 1.214D-10
  DOUBLE PRECISION, PARAMETER :: F72= 1.026D1
  DOUBLE PRECISION, PARAMETER :: F73= 1.306D1
  DOUBLE PRECISION, PARAMETER :: F74= 6.544D-2
  DOUBLE PRECISION, PARAMETER :: F75= 8.913D-3
  DOUBLE PRECISION, PARAMETER :: F76= 6.871D-7
  DOUBLE PRECISION, PARAMETER :: F77= 3.115D-9
  DOUBLE PRECISION, PARAMETER :: F78= 3.084D0
  DOUBLE PRECISION, PARAMETER :: F79= 3.594D-5
  DOUBLE PRECISION, PARAMETER :: F80= 3.414D-9
  DOUBLE PRECISION, PARAMETER :: F81= 4.779D-4
  DOUBLE PRECISION, PARAMETER :: F82= 1.188D-7
  DOUBLE PRECISION, PARAMETER :: F83= 7.497D-10
  DOUBLE PRECISION, PARAMETER :: F84= 2.914D-7
  DOUBLE PRECISION, PARAMETER :: F85= 9.159D-9
  DOUBLE PRECISION, PARAMETER :: F86= 2.681D-6
  DOUBLE PRECISION, PARAMETER :: F87= 1.626D-6
  DOUBLE PRECISION, PARAMETER :: F88= 5.791D-8
  DOUBLE PRECISION, PARAMETER :: F89= 2.562D-3
  DOUBLE PRECISION, PARAMETER :: F90= 1.774D-6
  DOUBLE PRECISION, PARAMETER :: F91= 1.702D-3
  DOUBLE PRECISION, PARAMETER :: F92= 6.724D-9
  DOUBLE PRECISION, PARAMETER :: F93= 3.002D1
  DOUBLE PRECISION, PARAMETER :: F94= 1.388D-4
  DOUBLE PRECISION, PARAMETER :: F95= 2.202D-3
  DOUBLE PRECISION, PARAMETER :: F96= 7.974D-1
  DOUBLE PRECISION, PARAMETER :: F97= 9.477D-6
  DOUBLE PRECISION, PARAMETER :: F98= 1.203D-10
  DOUBLE PRECISION, PARAMETER :: F99= 8.281D-8
  DOUBLE PRECISION, PARAMETER :: F100= 7.622D-6
  DOUBLE PRECISION, PARAMETER :: F101= 2.900D-10
  DOUBLE PRECISION, PARAMETER :: F102= 1.535D-6
  DOUBLE PRECISION, PARAMETER :: F103= 2.390D1
  DOUBLE PRECISION, PARAMETER :: F104= 7.393D-9
  DOUBLE PRECISION, PARAMETER :: F105= 3.518D0
  DOUBLE PRECISION, PARAMETER :: F106= 1.552D-7
  DOUBLE PRECISION, PARAMETER :: F107= 3.574D1
  DOUBLE PRECISION, PARAMETER :: F108= 2.943D-1
  DOUBLE PRECISION, PARAMETER :: F109= 1.576D-4
  DOUBLE PRECISION, PARAMETER :: F110= 8.019D-7
  DOUBLE PRECISION, PARAMETER :: F111= 6.506D-3
  DOUBLE PRECISION, PARAMETER :: F112= 1.276D-5
  DOUBLE PRECISION, PARAMETER :: F113= 2.701D-4
  DOUBLE PRECISION, PARAMETER :: F114= 4.207D-2
  DOUBLE PRECISION, PARAMETER :: F115= 5.217D-9
  DOUBLE PRECISION, PARAMETER :: F116= 2.476D-7
  DOUBLE PRECISION, PARAMETER :: F117= 1.385D-3
  DOUBLE PRECISION, PARAMETER :: F118= 1.995D-7
  DOUBLE PRECISION, PARAMETER :: F119= 1.586D-1
  DOUBLE PRECISION, PARAMETER :: F120= 6.573D-2
  DOUBLE PRECISION, PARAMETER :: F121= 2.432D-10
  DOUBLE PRECISION, PARAMETER :: F122= 2.876D-4
  DOUBLE PRECISION, PARAMETER :: F123= 5.244D-3
  DOUBLE PRECISION, PARAMETER :: F124= 5.284D-8
  DOUBLE PRECISION, PARAMETER :: F125= 2.259D-3
  DOUBLE PRECISION, PARAMETER :: F126= 1.011D-9
  DOUBLE PRECISION, PARAMETER :: F127= 1.342D-2
  DOUBLE PRECISION, PARAMETER :: F128= 3.279D0
  DOUBLE PRECISION, PARAMETER :: F129= 5.744D-8
  DOUBLE PRECISION, PARAMETER :: F130= 1.297D-6
  DOUBLE PRECISION, PARAMETER :: F131= 3.029D0
  DOUBLE PRECISION, PARAMETER :: F132= 1.068D-1
  DOUBLE PRECISION, PARAMETER :: F133= 1.295D-5
  call ZDPlasKin_bolsig_rates()
  Tgas = ZDPlasKin_cfg(1)
  Te  = ZDPlasKin_cfg(4)
  rrt(001) = FEC*bolsig_rates(bolsig_pointer(1))
  rrt(002) = FEC*bolsig_rates(bolsig_pointer(2))
  rrt(003) = FEC*bolsig_rates(bolsig_pointer(3))
  rrt(004) = FEC*bolsig_rates(bolsig_pointer(4))
  rrt(005) = FEC*bolsig_rates(bolsig_pointer(5))
  rrt(006) = FEC*bolsig_rates(bolsig_pointer(6))
  rrt(007) = FEC*bolsig_rates(bolsig_pointer(7))
  rrt(008) = FEC*bolsig_rates(bolsig_pointer(8))
  rrt(009) = FEC*bolsig_rates(bolsig_pointer(9))
  rrt(010) = FEC*bolsig_rates(bolsig_pointer(10))
  rrt(011) = FEC*bolsig_rates(bolsig_pointer(11))
  rrt(012) = FEC*bolsig_rates(bolsig_pointer(12))
  rrt(013) = FEC*bolsig_rates(bolsig_pointer(13))
  rrt(014) = FEC*bolsig_rates(bolsig_pointer(14))
  rrt(015) = FEC*bolsig_rates(bolsig_pointer(15))
  rrt(016) = FEC*bolsig_rates(bolsig_pointer(16))
  rrt(017) = FEC*bolsig_rates(bolsig_pointer(17))
  rrt(018) = FEC*bolsig_rates(bolsig_pointer(18))
  rrt(019) = FEC*bolsig_rates(bolsig_pointer(19))
  rrt(020) = FEC*bolsig_rates(bolsig_pointer(20))
  rrt(021) = FEC*bolsig_rates(bolsig_pointer(21))
  rrt(022) = FEC*bolsig_rates(bolsig_pointer(22))
  rrt(023) = FEC*bolsig_rates(bolsig_pointer(23))
  rrt(024) = FEC*bolsig_rates(bolsig_pointer(24))
  rrt(025) = FEC*bolsig_rates(bolsig_pointer(25))
  rrt(026) = FEC*bolsig_rates(bolsig_pointer(26))
  rrt(027) = FEC*bolsig_rates(bolsig_pointer(27))
  rrt(028) = FEC*bolsig_rates(bolsig_pointer(28))
  rrt(029) = FEC*bolsig_rates(bolsig_pointer(29))
  rrt(030) = FEC*bolsig_rates(bolsig_pointer(30))
  rrt(031) = FEC*bolsig_rates(bolsig_pointer(31))
  rrt(032) = FEC*bolsig_rates(bolsig_pointer(32))
  rrt(033) = FEC*bolsig_rates(bolsig_pointer(33))
  rrt(034) = FEC*bolsig_rates(bolsig_pointer(34))
  rrt(035) = FEC*bolsig_rates(bolsig_pointer(35))
  rrt(036) = FEC*bolsig_rates(bolsig_pointer(36))
  rrt(037) = FEC*bolsig_rates(bolsig_pointer(37))
  rrt(038) = FEC*bolsig_rates(bolsig_pointer(38))
  rrt(039) = FEC*bolsig_rates(bolsig_pointer(39))
  rrt(040) = FEC*bolsig_rates(bolsig_pointer(40))
  rrt(041) = FEC*bolsig_rates(bolsig_pointer(41))
  rrt(042) = FEC*bolsig_rates(bolsig_pointer(42))
  rrt(043) = FEC*bolsig_rates(bolsig_pointer(43))
  rrt(044) = FEC*bolsig_rates(bolsig_pointer(44))
  rrt(045) = FEC*bolsig_rates(bolsig_pointer(45))
  rrt(046) = FEC*bolsig_rates(bolsig_pointer(46))
  rrt(047) = FEC*bolsig_rates(bolsig_pointer(47))
  rrt(048) = FEC*bolsig_rates(bolsig_pointer(48))
  rrt(049) = FEC*bolsig_rates(bolsig_pointer(49))
  rrt(050) = FEC*bolsig_rates(bolsig_pointer(50))
  rrt(051) = FEC*bolsig_rates(bolsig_pointer(51))
  rrt(052) = FEC*bolsig_rates(bolsig_pointer(52))
  rrt(053) = FEC*bolsig_rates(bolsig_pointer(53))
  rrt(054) = FEC*bolsig_rates(bolsig_pointer(54))
  rrt(055) = FEC*bolsig_rates(bolsig_pointer(55))
  rrt(056) = FEC*bolsig_rates(bolsig_pointer(56))
  rrt(057) = FEC*bolsig_rates(bolsig_pointer(57))
  rrt(058) = FEC*bolsig_rates(bolsig_pointer(58))
  rrt(059) = FEC*bolsig_rates(bolsig_pointer(59))
  rrt(060) = FEC*bolsig_rates(bolsig_pointer(60))
  rrt(061) = FEC*bolsig_rates(bolsig_pointer(61))
  rrt(062) = FEC*bolsig_rates(bolsig_pointer(62))
  rrt(063) = FEC*bolsig_rates(bolsig_pointer(63))
  rrt(064) = FEC*bolsig_rates(bolsig_pointer(64))
  rrt(065) = FEC*bolsig_rates(bolsig_pointer(65))
  rrt(066) = FEC*bolsig_rates(bolsig_pointer(66))
  rrt(067) = FEC*bolsig_rates(bolsig_pointer(67))
  rrt(068) = FEC*bolsig_rates(bolsig_pointer(68))
  rrt(069) = FEC*bolsig_rates(bolsig_pointer(69))
  rrt(070) = FEC*bolsig_rates(bolsig_pointer(70))
  rrt(071) = FEC*bolsig_rates(bolsig_pointer(71))
  rrt(072) = FEC*bolsig_rates(bolsig_pointer(72))
  rrt(073) = FEC*bolsig_rates(bolsig_pointer(73))
  rrt(074) = FEC*bolsig_rates(bolsig_pointer(74))
  rrt(075) = FEC*bolsig_rates(bolsig_pointer(75))
  rrt(076) = FEC*bolsig_rates(bolsig_pointer(76))
  rrt(077) = FEC*bolsig_rates(bolsig_pointer(77))
  rrt(078) = FEC*bolsig_rates(bolsig_pointer(78))
  rrt(079) = FEC*bolsig_rates(bolsig_pointer(79))
  rrt(080) = FEC*bolsig_rates(bolsig_pointer(80))
  rrt(081) = FEC*bolsig_rates(bolsig_pointer(81))
  rrt(082) = FEC*bolsig_rates(bolsig_pointer(82))
  rrt(083) = FEC*bolsig_rates(bolsig_pointer(83))
  rrt(084) = FEC*bolsig_rates(bolsig_pointer(84))
  rrt(085) = FEC*bolsig_rates(bolsig_pointer(85))
  rrt(086) = FEC*bolsig_rates(bolsig_pointer(86))
  rrt(087) = FEC*bolsig_rates(bolsig_pointer(87))
  rrt(088) = FEC*bolsig_rates(bolsig_pointer(88))
  rrt(089) = FEC*bolsig_rates(bolsig_pointer(89))
  rrt(090) = FEC*bolsig_rates(bolsig_pointer(90))
  rrt(091) = FEC*bolsig_rates(bolsig_pointer(91))
  rrt(092) = FEC*bolsig_rates(bolsig_pointer(92))
  rrt(093) = FEC*bolsig_rates(bolsig_pointer(93))
  rrt(094) = FEC*bolsig_rates(bolsig_pointer(94))
  rrt(095) = FEC*bolsig_rates(bolsig_pointer(95))
  rrt(096) = FEC*bolsig_rates(bolsig_pointer(96))
  rrt(097) = FEC*bolsig_rates(bolsig_pointer(97))
  rrt(098) = FEC*bolsig_rates(bolsig_pointer(98))
  rrt(099) = FEC*bolsig_rates(bolsig_pointer(99))
  rrt(100) = FEC*bolsig_rates(bolsig_pointer(100))
  rrt(101) = FEC*bolsig_rates(bolsig_pointer(101))
  rrt(102) = FEC*bolsig_rates(bolsig_pointer(102))
  rrt(103) = FEC*bolsig_rates(bolsig_pointer(103))
  rrt(104) = FEC*bolsig_rates(bolsig_pointer(104))
  rrt(105) = FEC*bolsig_rates(bolsig_pointer(105))
  rrt(106) = FEC*bolsig_rates(bolsig_pointer(106))
  rrt(107) = FEC*bolsig_rates(bolsig_pointer(107))
  rrt(108) = FEC*bolsig_rates(bolsig_pointer(108))
  rrt(109) = FEC*bolsig_rates(bolsig_pointer(109))
  rrt(110) = FEC*bolsig_rates(bolsig_pointer(110))
  rrt(111) = FEC*bolsig_rates(bolsig_pointer(111))
  rrt(112) = FEC*bolsig_rates(bolsig_pointer(112))
  rrt(113) = FEC*bolsig_rates(bolsig_pointer(113))
  rrt(114) = FEC*bolsig_rates(bolsig_pointer(114))
  rrt(115) = FEC*bolsig_rates(bolsig_pointer(115))
  rrt(116) = FEC*bolsig_rates(bolsig_pointer(116))
  rrt(117) = FEC*bolsig_rates(bolsig_pointer(117))
  rrt(118) = FEC*bolsig_rates(bolsig_pointer(118))
  rrt(119) = FEC*bolsig_rates(bolsig_pointer(119))
  rrt(120) = FEC*bolsig_rates(bolsig_pointer(120))
  rrt(121) = FEC*bolsig_rates(bolsig_pointer(121))
  rrt(122) = FEC*bolsig_rates(bolsig_pointer(122))
  rrt(123) = FEC*bolsig_rates(bolsig_pointer(123))
  rrt(124) = FEC*bolsig_rates(bolsig_pointer(124))
  rrt(125) = FEC*bolsig_rates(bolsig_pointer(125))
  rrt(126) = FEC*bolsig_rates(bolsig_pointer(126))
  rrt(127) = FEC*bolsig_rates(bolsig_pointer(127))
  rrt(128) = FEC*bolsig_rates(bolsig_pointer(128))
  rrt(129) = FEC*bolsig_rates(bolsig_pointer(129))
  rrt(130) = FEC*bolsig_rates(bolsig_pointer(130))
  rrt(131) = FEC*bolsig_rates(bolsig_pointer(131))
  rrt(132) = FEC*bolsig_rates(bolsig_pointer(132))
  rrt(133) = FEC*bolsig_rates(bolsig_pointer(133))
  rrt(134) = FEC*bolsig_rates(bolsig_pointer(134))
  rrt(135) = FEC*bolsig_rates(bolsig_pointer(135))
  rrt(136) = FEC*bolsig_rates(bolsig_pointer(136))
  rrt(137) = FEC*bolsig_rates(bolsig_pointer(137))
  rrt(138) = FEC*bolsig_rates(bolsig_pointer(138))
  rrt(139) = FEC*bolsig_rates(bolsig_pointer(139))
  rrt(140) = FEC*bolsig_rates(bolsig_pointer(140))
  rrt(141) = FEC*bolsig_rates(bolsig_pointer(141))
  rrt(142) = FEC*bolsig_rates(bolsig_pointer(142))
  rrt(143) = FEC*bolsig_rates(bolsig_pointer(143))
  rrt(144) = FEC*bolsig_rates(bolsig_pointer(144))
  rrt(145) = FEC*bolsig_rates(bolsig_pointer(145))
  rrt(146) = FEC*bolsig_rates(bolsig_pointer(146))
  rrt(147) = FEC*bolsig_rates(bolsig_pointer(147))
  rrt(148) = FEC*bolsig_rates(bolsig_pointer(148))
  rrt(149) = FEC*bolsig_rates(bolsig_pointer(149))
  rrt(150) = FEC*bolsig_rates(bolsig_pointer(150))
  rrt(151) = FEC*bolsig_rates(bolsig_pointer(151))
  rrt(152) = FEC*bolsig_rates(bolsig_pointer(152))
  rrt(153) = FEC*bolsig_rates(bolsig_pointer(153))
  rrt(154) = FEC*bolsig_rates(bolsig_pointer(154))
  rrt(155) = FEC*bolsig_rates(bolsig_pointer(155))
  rrt(156) = FEC*bolsig_rates(bolsig_pointer(156))
  rrt(157) = FEC*bolsig_rates(bolsig_pointer(157))
  rrt(158) = FEC*bolsig_rates(bolsig_pointer(158))
  rrt(159) = FEC*bolsig_rates(bolsig_pointer(159))
  rrt(160) = FEC*bolsig_rates(bolsig_pointer(160))
  rrt(161) = FEC*bolsig_rates(bolsig_pointer(161))
  rrt(162) = FEC*bolsig_rates(bolsig_pointer(162))
  rrt(163) = FEC*bolsig_rates(bolsig_pointer(163))
  rrt(164) = FEC*bolsig_rates(bolsig_pointer(164))
  rrt(165) = FEC*bolsig_rates(bolsig_pointer(165))
  rrt(166) = FEC*bolsig_rates(bolsig_pointer(166))
  rrt(167) = FEC*bolsig_rates(bolsig_pointer(167))
  rrt(168) = FEC*bolsig_rates(bolsig_pointer(168))
  rrt(169) = FEC*bolsig_rates(bolsig_pointer(169))
  rrt(170) = FEC*bolsig_rates(bolsig_pointer(170))
  rrt(171) = FEC*bolsig_rates(bolsig_pointer(171))
  rrt(172) = FEC*bolsig_rates(bolsig_pointer(172))
  rrt(173) = FEC*bolsig_rates(bolsig_pointer(173))
  rrt(174) = FEC*bolsig_rates(bolsig_pointer(174))
  rrt(175) = FEC*bolsig_rates(bolsig_pointer(175))
  rrt(176) = FEC*bolsig_rates(bolsig_pointer(176))
  rrt(177) = FEC*bolsig_rates(bolsig_pointer(177))
  rrt(178) = FEC*bolsig_rates(bolsig_pointer(178))
  rrt(179) = FEC*bolsig_rates(bolsig_pointer(179))
  rrt(180) = FEC*bolsig_rates(bolsig_pointer(180))
  rrt(181) = FEC*bolsig_rates(bolsig_pointer(181))
  rrt(182) = FEC*bolsig_rates(bolsig_pointer(182))
  rrt(183) = FEC*bolsig_rates(bolsig_pointer(183))
  rrt(184) = FEC*bolsig_rates(bolsig_pointer(184))
  rrt(185) = FEC*bolsig_rates(bolsig_pointer(185))
  rrt(186) = FEC*bolsig_rates(bolsig_pointer(186))
  rrt(187) = FEC*bolsig_rates(bolsig_pointer(187))
  rrt(188) = FEC*bolsig_rates(bolsig_pointer(188))
  rrt(189) = FEC*bolsig_rates(bolsig_pointer(189))
  rrt(190) = FEC*bolsig_rates(bolsig_pointer(190))
  rrt(191) = FEC*bolsig_rates(bolsig_pointer(191))
  rrt(192) = FEC*bolsig_rates(bolsig_pointer(192))
  rrt(193) = FEC*bolsig_rates(bolsig_pointer(193))
  rrt(194) = FEC*bolsig_rates(bolsig_pointer(194))
  rrt(195) = FEC*bolsig_rates(bolsig_pointer(195))
  rrt(196) = FEC*bolsig_rates(bolsig_pointer(196))
  rrt(197) = FEC*bolsig_rates(bolsig_pointer(197))
  rrt(198) = FEC*bolsig_rates(bolsig_pointer(198))
  rrt(199) = FEC*bolsig_rates(bolsig_pointer(199))
  rrt(200) = FEC*bolsig_rates(bolsig_pointer(200))
  rrt(201) = FEC*bolsig_rates(bolsig_pointer(201))
  rrt(202) = FEC*bolsig_rates(bolsig_pointer(202))
  rrt(203) = FEC*bolsig_rates(bolsig_pointer(203))
  rrt(204) = FEC*bolsig_rates(bolsig_pointer(204))
  rrt(205) = FEC*bolsig_rates(bolsig_pointer(205))
  rrt(206) = FEC*bolsig_rates(bolsig_pointer(206))
  rrt(207) = FEC*bolsig_rates(bolsig_pointer(207))
  rrt(208) = FEC*bolsig_rates(bolsig_pointer(208))
  rrt(209) = FEC*bolsig_rates(bolsig_pointer(209))
  rrt(210) = FEC*bolsig_rates(bolsig_pointer(210))
  rrt(211) = FEC*bolsig_rates(bolsig_pointer(211))
  rrt(212) = FEC*bolsig_rates(bolsig_pointer(212))
  rrt(213) = FEC*bolsig_rates(bolsig_pointer(213))
  rrt(214) = FEC*bolsig_rates(bolsig_pointer(214))
  rrt(215) = FEC*bolsig_rates(bolsig_pointer(215))
  rrt(216) = FEC*bolsig_rates(bolsig_pointer(216))
  rrt(217) = FEC*bolsig_rates(bolsig_pointer(217))
  rrt(218) = FEC*bolsig_rates(bolsig_pointer(218))
  rrt(219) = FEC*bolsig_rates(bolsig_pointer(219))
  rrt(220) = FEC*bolsig_rates(bolsig_pointer(220))
  rrt(221) = FEC*bolsig_rates(bolsig_pointer(221))
  rrt(222) = FEC*bolsig_rates(bolsig_pointer(222))
  rrt(223) = FEC*bolsig_rates(bolsig_pointer(223))
  rrt(224) = FEC*bolsig_rates(bolsig_pointer(224))
  rrt(225) = FEC*bolsig_rates(bolsig_pointer(225))
  rrt(226) = FEC*bolsig_rates(bolsig_pointer(226))
  rrt(227) = FEC*bolsig_rates(bolsig_pointer(227))
  rrt(228) = FEC*bolsig_rates(bolsig_pointer(228))
  rrt(229) = FEC*bolsig_rates(bolsig_pointer(229))
  rrt(230) = FEC*bolsig_rates(bolsig_pointer(230))
  rrt(231) = FEC*bolsig_rates(bolsig_pointer(231))
  rrt(232) = FEC*bolsig_rates(bolsig_pointer(232))
  rrt(233) = FEC*bolsig_rates(bolsig_pointer(233))
  rrt(234) = FEC*bolsig_rates(bolsig_pointer(234))
  rrt(235) = FEC*bolsig_rates(bolsig_pointer(235))
  rrt(236) = FEC*bolsig_rates(bolsig_pointer(236))
  rrt(237) = FEC*bolsig_rates(bolsig_pointer(237))
  rrt(238) = FEC*bolsig_rates(bolsig_pointer(238))
  rrt(239) = FEC*bolsig_rates(bolsig_pointer(239))
  rrt(240) = FIR
  rrt(241) = FIR
  rrt(242) = FIR
  rrt(243) = FIR
  rrt(244) = FIR
  rrt(245) = FIR
  rrt(246) = FIR
  rrt(247) = FIR
  rrt(248) = FIR
  rrt(249) = FIR
  rrt(250) = FIR
  rrt(251) = FIR
  rrt(252) = FIR
  rrt(253) = FIR
  rrt(254) = FIR
  rrt(255) = FIR
  rrt(256) = FID*2.57D-07*(300./TGAS)**0.3
  rrt(257) = FID*6.61D-08*(300./TGAS)**0.3
  rrt(258) = FID*1.18D-08*(300./TGAS)**0.5
  rrt(259) = FID*2.42D-08*(300./TGAS)**0.5
  rrt(260) = FID*1.41D-08*(300./TGAS)**0.5
  rrt(261) = FID*2.25D-08*(300./TGAS)**0.5
  rrt(262) = FID*7.88D-09*(300./TGAS)**0.5
  rrt(263) = FID*1.00D-08*(300./TGAS)**0.5
  rrt(264) = FID*2.19D-08*(300./TGAS)**0.71
  rrt(265) = FID*3.36D-08*(300./TGAS)**0.71
  rrt(266) = FID*7.70D-09*(300./TGAS)**0.71
  rrt(267) = FID*1.92D-08*(300./TGAS)**0.71
  rrt(268) = FID*1.60D-08*(300./TGAS)**0.71
  rrt(269) = FID*8.98D-09*(300./TGAS)**0.71
  rrt(270) = FID*9.62D-09*(300./TGAS)**0.71
  rrt(271) = FID*8.29D-09*(300./TGAS)**0.71
  rrt(272) = FID*3.43D-08*(300./TGAS)**0.71
  rrt(273) = FID*1.34D-08*(300./TGAS)**0.71
  rrt(274) = FID*4.87D-09*(300./TGAS)**0.71
  rrt(275) = FID*3.17D21/(6.022D23*TE**4.5)
  rrt(276) = rrt(275)
  rrt(277) = FIN*2.11D-09
  rrt(278) = FIN*9.60D-10
  rrt(279) = FIN*6.90D-10
  rrt(280) = FIN*2.25D-10
  rrt(281) = FIN*1.50D-9
  rrt(282) = FIN*1.60D-9
  rrt(283) = FIN*1.50D-10
  rrt(284) = rrt(281)
  rrt(285) = FIN*1.91D-9
  rrt(286) = FIN*4.23D-10
  rrt(287) = FIN*1.38D-9
  rrt(288) = FIN*1.23D-9
  rrt(289) = FIN*1.13D-9
  rrt(290) = FIN*3.30D-11
  rrt(291) = FIN*1.00D-11
  rrt(292) = FIN*1.36D-10
  rrt(293) = FIN*1.20D-9
  rrt(294) = FIN*9.90D-10
  rrt(295) = FIN*7.10D-10
  rrt(296) = FIN*1.48D-9
  rrt(297) = FIN*3.50D-10
  rrt(298) = FIN*3.00D-10
  rrt(299) = FIN*1.38D-10
  rrt(300) = FIN*3.60D-10
  rrt(301) = FIN*8.40D-10
  rrt(302) = FIN*2.31D-10
  rrt(303) = FIN*3.97D-10
  rrt(304) = rrt(282)
  rrt(305) = FIN*6.50D-11
  rrt(306) = FIN*1.09D-9
  rrt(307) = FIN*1.43D-10
  rrt(308) = rrt(293)
  rrt(309) = FIN*1.15D-9
  rrt(310) = FIN*2.47D-10
  rrt(311) = FIN*1.00D-10
  rrt(312) = rrt(291)
  rrt(313) = FIN*5.00D-10
  rrt(314) = rrt(313)
  rrt(315) = rrt(298)
  rrt(316) = FIN*2.91D-10
  rrt(317) = FIN*8.90D-10
  rrt(318) = FIN*3.30D-10
  rrt(319) = FIN*6.80D-11
  rrt(320) = FIN*4.10D-9
  rrt(321) = FIN*1.31D-10
  rrt(322) = FIN*2.48D-10
  rrt(323) = FIN*4.14D-10
  rrt(324) = rrt(318)
  rrt(325) = rrt(291)
  rrt(326) = FIN*3.74D-10
  rrt(327) = FIN*2.40D-9
  rrt(328) = FIN*2.10D-9
  rrt(329) = FIN*1.70D-9
  rrt(330) = rrt(293)
  rrt(331) = rrt(327)
  rrt(332) = FIN*1.40D-9
  rrt(333) = rrt(309)
  rrt(334) = rrt(309)
  rrt(335) = FIN*2.00D-9
  rrt(336) = rrt(329)
  rrt(337) = FIN*3.50D-9
  rrt(338) = FIN*1.14D-10
  rrt(339) = rrt(332)
  rrt(340) = FIN*2.30D-9
  rrt(341) = FIN*1.00D-9
  rrt(342) = rrt(341)
  rrt(343) = rrt(295)
  rrt(344) = rrt(295)
  rrt(345) = FIN*2.94D-10
  rrt(346) = FIN*1.37D-9
  rrt(347) = FIN*2.35D-9
  rrt(348) = FIN*6.86D-10
  rrt(349) = FIN*1.96D-10
  rrt(350) = FIN*2.21D-9
  rrt(351) = FIN*1.81D-9
  rrt(352) = FIN*8.82D-10
  rrt(353) = FIN*4.80D-10
  rrt(354) = FIN*4.82D-9
  rrt(355) = rrt(328)
  rrt(356) = FIN*6.39D-10
  rrt(357) = rrt(281)
  rrt(358) = rrt(340)
  rrt(359) = FIN*3.40D-9
  rrt(360) = rrt(332)
  rrt(361) = rrt(332)
  rrt(362) = FIN*1.90D-9
  rrt(363) = FIN*1.30D-9
  rrt(364) = rrt(332)
  rrt(365) = FIN*2.80D-9
  rrt(366) = FIN*1.65D-9
  rrt(367) = FIN*3.06D-9
  rrt(368) = rrt(341)
  rrt(369) = FIN*3.00D-9
  rrt(370) = rrt(341)
  rrt(371) = rrt(335)
  rrt(372) = rrt(335)
  rrt(373) = FIN*5.40D-10
  rrt(374) = F0*4.08D-18*TGAS**2*EXP(-4163./TGAS)
  rrt(375) = F1*9.97D-11
  rrt(376) = F2*2.51D-15*(TGAS/298.)**4.14*EXP(-52.55/(R*TGAS))
  rrt(377) = F3*4.26D-15*(TGAS/298.)**4.02*EXP(-22.86/(R*TGAS))
  rrt(378) = F4*3.01D-12*EXP(-2.08/(R*TGAS))
  rrt(379) = F5*3.54D-16*(TGAS/298.)**4.02*EXP(-45.48/(R*TGAS))
  rrt(380) = F6*1.71D-14*(TGAS/298.)**3.40*EXP(-97.28/(R*TGAS))
  rrt(381) = F7*9.86D-13*(TGAS/298.)**3.00*EXP(-36.67/(R*TGAS))
  rrt(382) = F8*1.33D-10*EXP(-167.00/(R*TGAS))
  rrt(383) = F9*5.68D-17*(TGAS/298.)**3.72*EXP(-33.42/(R*TGAS))
  rrt(384) = F10*1.90D-12
  rrt(385) = F11*2.40D16*EXP(-52800./TGAS)
  rrt(386) = F12*1.69D-8*EXP(-379./(R*TGAS))
  rrt(387) = F13*6.97D-9*EXP(-345./(R*TGAS))
  rrt(388) = F14*3.00D-44*TGAS**9.10
  rrt(389) = F15*3.32D-10*EXP(-45.98/(R*TGAS))
  rrt(390) = F16*3.01D-11
  rrt(391) = F17*3.01D-11
  rrt(392) = F18*3.01D-11
  rrt(393) = F19*1.61D-15*(TGAS/298.)**3.65*EXP(-29.93/(R*TGAS))
  rrt(394) = F20*3.01D-11
  rrt(395) = F21*3.01D-11
  rrt(396) = F22*1.20D-12*EXP(-25.94/(R*TGAS))
  rrt(397) = F23*8.30D-19*TGAS**2.00*EXP(-3938.65/TGAS)
  rrt(398) = F24*1.00D-11*EXP(7.48/(R*TGAS))
  rrt(399) = F25*6.64D-9*EXP(-348./(R*TGAS))
  rrt(400) = F26*2.66D-9*EXP(-6011.07/TGAS)
  rrt(401) = F27*9.96D-10
  rrt(402) = F28*3.00D-11
  rrt(403) = F29*1.14D-29
  rrt(404) = F30*1.79D-10*EXP(-1565.17/TGAS)
  rrt(405) = F31*4.98D-11
  rrt(406) = F32*6.64D-11
  rrt(407) = F33*3.29D-12*TGAS**0.43*EXP(186.21/TGAS)
  rrt(408) = F34*8.30D-11
  rrt(409) = F35*1.46D-13*(TGAS/298.)**3.30*EXP(-43.90/(R*TGAS))
  rrt(410) = F36*1.19D-15*(TGAS/298.)**3.82*EXP(-37.83/(R*TGAS))
  rrt(411) = F37*5.71D-14*(TGAS/298.)**3.30*EXP(-83.06/(R*TGAS))
  rrt(412) = F38*1.23D-11*(TGAS/298.)**1.50*EXP(-31.01/(R*TGAS))
  rrt(413) = F39*8.97D-20*EXP(-48.64/(R*TGAS))
  rrt(414) = F40*1.80D21*TGAS**(-1.24)*EXP(-45700./TGAS)
  rrt(415) = F41*8.30D-13*EXP(-62.77/(R*TGAS))
  rrt(416) = F42*1.79D-10*EXP(132.36/TGAS)
  rrt(417) = F43*9.00D-33*TGAS**6.43
  rrt(418) = F44*2.41D-12
  rrt(419) = F45*5.83D-14*(TGAS/298.)**3.13*EXP(-75.33/(R*TGAS))
  rrt(420) = F46*4.50D-13*EXP(-98.11/(R*TGAS))
  rrt(421) = F47*3.01D-12
  rrt(422) = F48*1.61D-15*(TGAS/298.)**3.65*EXP(-38.25/(R*TGAS))
  rrt(423) = F49*1.91D-12
  rrt(424) = F50*2.41D-12
  rrt(425) = F51*1.69D-15*(TGAS/298.)**3.50*EXP(-35.34/(R*TGAS))
  rrt(426) = F52*4.12D-15*(TGAS/298.)**3.60*EXP(-27.77/(R*TGAS))
  rrt(427) = F53*5.99D-11
  rrt(428) = F54*3.32D-12
  rrt(429) = F55*8.65D-7*TGAS**(-0.99)*EXP(-795.17/TGAS)
  rrt(430) = F56*4.08D12*(TGAS/298.)**1.04*EXP(-154./(R*TGAS))
  rrt(431) = F57*9.55D-12
  rrt(432) = F58*1.40D-12
  rrt(433) = F59*4.42D-11
  rrt(434) = F60*9.00D-10*EXP(-62.36/(R*TGAS))
  rrt(435) = F61*9.68D-12*(TGAS/298.)**1.28*EXP(-5.40/(R*TGAS))
  rrt(436) = F62*9.55D-12
  rrt(437) = F63*4.30D-7*EXP(-404./(R*TGAS))
  rrt(438) = F64*9.60D-11*EXP(-216./(R*TGAS))
  rrt(439) = F65*4.00D-11*EXP(-286./(R*TGAS))
  rrt(440) = F66*1.00D-10*EXP(-316./(R*TGAS))
  rrt(441) = F67*8.00D-10*EXP(-299./(R*TGAS))
  rrt(442) = F68*4.23D-18*TGAS**1.60*EXP(-2868.65/TGAS)
  rrt(443) = F69*8.00D12*TGAS**0.44*EXP(-44675.39/TGAS)
  rrt(444) = F70*3.00D-14*(TGAS/298.)**2.48*EXP(-25.65/(R*TGAS))
  rrt(445) = F71*4.75D-16*EXP(-180./(R*TGAS))
  rrt(446) = F72*5.30D-12*EXP(-2660./TGAS)
  rrt(447) = F73*5.00D-14*EXP(-25.53/(R*TGAS))
  rrt(448) = F74*3.50D-11
  rrt(449) = F75*1.46D-13*(TGAS/298.)**3.30*EXP(-43.90/(R*TGAS))
  rrt(450) = F76*2.01D-12
  rrt(451) = F77*2.01D-12
  rrt(452) = F78*1.68D-15*(TGAS/298.)**3.50*EXP(-19.62/(R*TGAS))
  rrt(453) = F79*8.00D-12
  rrt(454) = F80*1.61D-13*(TGAS/298.)**2.63*EXP(-35.75/(R*TGAS))
  rrt(455) = F81*1.60D-10
  rrt(456) = F82*1.01D-11*TGAS**0.27*EXP(-140.92/TGAS)
  rrt(457) = F83*2.00D14*EXP(-20000./TGAS)
  rrt(458) = F84*1.40D-12
  rrt(459) = F85*9.30D-12*EXP(-1207.85/TGAS)
  rrt(460) = F86*5.00D-13*EXP(-163./(R*TGAS))
  rrt(461) = F87*4.00D-12*EXP(-272./(R*TGAS))
  rrt(462) = F88*1.58D-5*(TGAS/298.)**-8.58*EXP(-84.81/(R*TGAS))
  rrt(463) = F89*1.20D-12*EXP(-37.66/(R*TGAS))
  rrt(464) = F90*5.71D-14*(TGAS/298.)**3.30*EXP(-83.06/(R*TGAS))
  rrt(465) = F91*2.19D-180*TGAS**2.54*EXP(-3400.10/TGAS)
  rrt(466) = F92*1.10D17*EXP(-42470./TGAS)
  rrt(467) = F93*1.61D-15*(TGAS/298.)**3.50*EXP(-27.77/(R*TGAS))
  rrt(468) = F94*4.42D-12
  rrt(469) = F95*2.81D-12
  rrt(470) = F96*1.69D-15*(TGAS/298.)**3.50*EXP(-27.77/(R*TGAS))
  rrt(471) = F97*2.41D-12*EXP(-0.55/(R*TGAS))
  rrt(472) = F98*3.19D-14*(TGAS/298.)**2.84*EXP(-38.25/(R*TGAS))
  rrt(473) = F99*3.01D-12
  rrt(474) = F100*6.00D-11
  rrt(475) = F101*6.74D-18*TGAS**2.19*EXP(-447.91/TGAS)
  rrt(476) = F102*1.25D17*EXP(-237./(R*TGAS))
  rrt(477) = F103*16.0*(TGAS/298.)**-(10.00)*EXP(-150./(R*TGAS))
  rrt(478) = F104*2.41D-12
  rrt(479) = F105*6.71D-11*EXP(-196./(R*TGAS))
  rrt(480) = F106*4.20D-10*EXP(-231./(R*TGAS))
  rrt(481) = F107*2.50D15*EXP(-363./(R*TGAS))
  rrt(482) = F108*4.40D-13*(TGAS/298.)**2.50*EXP(-10.39/(R*TGAS))
  rrt(483) = F109*1.29D-11*(TGAS/298.)**0.51*EXP(-5.15/(R*TGAS))
  rrt(484) = F110*1.28D13*(TGAS/298.)**(-15.70)*EXP(-502./(R*TGAS))
  rrt(485) = F111*1.27D-14*(TGAS/298.)**2.67*EXP(-28.66/(R*TGAS))
  rrt(486) = F112*1.69D-15*(TGAS/298.)**3.50*EXP(-27.77/(R*TGAS))
  rrt(487) = F113*1.39D-13*(TGAS/298.)**2.38*EXP(-79.49/(R*TGAS))
  rrt(488) = F114*78.80*(TGAS/298.)**(-11.76)*EXP(-98.53/(R*TGAS))
  rrt(489) = F115*1.26D13*EXP(-140./(R*TGAS))
  rrt(490) = F116*1.07D2*(TGAS/298.)**(-11.90)*EXP(-135./(R*TGAS))
  rrt(491) = F117*3.01D-11
  rrt(492) = F118*7.71D13*(TGAS/298.)**0.77*EXP(-128./(R*TGAS))
  rrt(493) = F119*2.52D-14*(TGAS/298.)**2.72*EXP(-40.99/(R*TGAS))
  rrt(494) = F120*8.32D-13*EXP(-56.87/(R*TGAS))
  rrt(495) = F121*8.87D-7*EXP(-180./(R*TGAS))
  rrt(496) = F122*7.84D-6*EXP(-207./(R*TGAS))
  rrt(497) = F123*2.19D-10*EXP(-39.24/(R*TGAS))
  rrt(498) = F124*1.81D-12*EXP(-20.54/(R*TGAS))
  rrt(499) = F125*2.42D-15*(TGAS/298.)**3.65*EXP(-21.62/(R*TGAS))
  rrt(500) = F126*2.47D-15*(TGAS/298.)**3.65*EXP(-38.25/(R*TGAS))
  rrt(501) = F127*1.00D-11
  rrt(502) = F128*2.47D-15*(TGAS/298.)**3.65*EXP(-38.25/(R*TGAS))
  rrt(503) = F129*8.63D-14*(TGAS/298.)**3.30*EXP(-83.06/(R*TGAS))
  rrt(504) = F130*9.61D-13
  rrt(505) = F131*3.16D16*EXP(-331./(R*TGAS))
  rrt(506) = F132*1.88D-8*(TGAS/298.)**(-1.10)*EXP(-437./(R*TGAS))
  rrt(507) = F133*5.52D-30*TGAS**(-1.00)
  where( lreaction_block(:) ) rrt(:) = 0.0d0
  return
end subroutine ZDPlasKin_reac_rates
!-----------------------------------------------------------------------------------------------------------------------------------
!
! END OF FILE
!
!-----------------------------------------------------------------------------------------------------------------------------------
