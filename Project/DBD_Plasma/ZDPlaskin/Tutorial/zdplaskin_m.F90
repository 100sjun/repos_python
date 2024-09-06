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
! Wed Sep  4 13:40:39 2024
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
  integer, parameter :: species_max = 50, species_electrons = 1, species_length = 9, reactions_max = 485, reactions_length = 28
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
  /-1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,&
    0, 1, 1, 0, 1, 1, 1, 0/
  data species_name(1:species_max) &
  /"E        ","H2       ","C3H5     ","C2H4(V2) ","CH4      ","C2H3^+   ","H^+      ","C3H5^+   ","C3H6     ","CH2      ",&
   "C5H12    ","C3H8     ","H        ","C2H2^+   ","C2H6^+   ","C2H4(V1) ","C3H4     ","C2H      ","C2H2     ","C2H6     ",&
   "C2H2(V5) ","C3H8^+   ","C3H7     ","CH3      ","CH4^+    ","H3^+     ","CH2^+    ","CH^+     ","C2H2(V13)","C3H7^+   ",&
   "C2H2(V2) ","C2H4^+   ","C2H6(V24)","C2H6(V13)","C2H4     ","CH4(V13) ","C3H8(V1) ","C2H5     ","C4H9     ","C3H8(V2) ",&
   "CH5^+    ","C3H6(V)  ","CH4(V24) ","C2H5^+   ","H2^+     ","CH       ","C3H4^+   ","CH3^+    ","C3H6^+   ","C2H3     "/
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
   "E+E+C3H4^+=>E+C3H4          ","E+E+H^+=>E+H                ","E+E+H2^+=>E+H2              ","E+E+CH5^+=>E+CH3+H+H        ",&
   "E+E+CH5^+=>E+CH2+H2+H       ","E+E+CH4^+=>E+CH3+H          ","E+E+CH4^+=>E+CH2+H+H        ","E+E+CH4^+=>E+CH+H2+H        ",&
   "E+E+CH3^+=>E+CH2+H          ","E+E+CH3^+=>E+CH+H2          ","E+E+CH2^+=>E+CH+H           ","E+E+C2H6^+=>E+C2H5+H        ",&
   "E+E+C2H6^+=>E+C2H4+H+H      ","E+E+C2H5^+=>E+C2H4+H        ","E+E+C2H5^+=>E+C2H3+H+H      ","E+E+C2H5^+=>E+C2H2+H2+H     ",&
   "E+E+C2H5^+=>E+C2H2+H+H+H    ","E+E+C2H5^+=>E+CH3+CH2       ","E+E+C2H4^+=>E+C2H3+H        ","E+E+C2H4^+=>E+C2H2+H+H      ",&
   "E+E+C2H3^+=>E+C2H2+H        ","E+E+C2H2^+=>E+CH+CH         ","E+E+H2^+=>E+H+H             ","E+E+H3^+=>E+H+H2            ",&
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
   "H2+C2H2^+=>H+C2H3^+         ","CH4+H3^+=>H2+CH5^+          ","CH3+H3^+=>H2+CH4^+          ","CH2+H3^+=>H2+CH3^+          ",&
   "CH+H3^+=>H2+CH2^+           ","C2H6+H3^+=>H2+H2+C2H5^+     ","C2H5+H3^+=>H2+C2H6^+        ","C2H4+H3^+=>H2+C2H5^+        ",&
   "C2H4+H3^+=>H2+H2+C2H3^+     ","C2H3+H3^+=>H2+C2H4^+        ","C2H+H3^+=>H2+C2H2^+         ","C2H2+H3^+=>H2+C2H3^+        ",&
   "CH4+H2^+=>H+CH5^+           ","CH4+H2^+=>H2+CH4^+          ","CH4+H2^+=>H2+H+CH3^+        ","CH2+H2^+=>H+CH3^+           ",&
   "CH2+H2^+=>H2+CH2^+          ","CH+H2^+=>H+CH2^+            ","CH+H2^+=>H2+CH^+            ","C2H6+H2^+=>H2+C2H6^+        ",&
   "C2H6+H2^+=>H2+H+C2H5^+      ","C2H6+H2^+=>H2+H2+C2H4^+     ","C2H6+H2^+=>H2+H2+H+C2H3^+   ","C2H6+H2^+=>H2+H2+H2+C2H2^+  ",&
   "C2H4+H2^+=>H2+C2H4^+        ","C2H4+H2^+=>H2+H+C2H3^+      ","C2H4+H2^+=>H2+H2+C2H2^+     ","C2H2+H2^+=>H+C2H3^+         ",&
   "C2H2+H2^+=>H2+C2H2^+        ","H+H2^+=>H3^+                ","H+H2^+=>H2+H^+              ","CH4+H^+=>H+CH4^+            ",&
   "CH4+H^+=>H2+CH3^+           ","CH3+H^+=>H+CH3^+            ","CH2+H^+=>H+CH2^+            ","CH2+H^+=>H2+CH^+            "/
  data reaction_sign(361:432) &
  /"CH+H^+=>H+CH^+              ","C2H6+H^+=>H2+C2H5^+         ","C2H6+H^+=>H2+H+C2H4^+       ","C2H6+H^+=>H2+H2+C2H3^+      ",&
   "C2H5+H^+=>H2+C2H4^+         ","C2H5+H^+=>H2+H+C2H3^+       ","C2H4+H^+=>H+C2H4^+          ","C2H4+H^+=>H2+C2H3^+         ",&
   "C2H4+H^+=>H2+H+C2H2^+       ","C2H3+H^+=>H+C2H3^+          ","C2H3+H^+=>H2+C2H2^+         ","C2H2+H^+=>H+C2H2^+          ",&
   "CH4+CH2=>CH3+CH3            ","CH4+CH=>C2H4+H              ","CH4+C2H5=>C2H6+CH3          ","CH4+C2H3=>C2H4+CH3          ",&
   "CH4+C2H=>C2H2+CH3           ","CH4+C3H7=>C3H8+CH3          ","CH4+C3H5=>C3H6+CH3          ","CH4+H=>CH3+H2               ",&
   "CH4+CH3=>H+C2H6             ","CH4+CH2=>C2H6               ","CH4=>CH3+H                  ","CH3=>CH2+H                  ",&
   "CH3=>CH+H2                  ","CH3+C2H5=>C2H6+CH2          ","CH2+CH2=>C2H2+H+H           ","CH2+C2H5=>C2H4+CH3          ",&
   "CH2+C2H3=>C2H2+CH3          ","CH2+C2H=>C2H2+CH            ","CH2+C3H8=>C3H7+CH3          ","CH2+C3H7=>C2H4+C2H5         ",&
   "CH2+C3H7=>C3H6+CH3          ","CH2+C3H6=>C3H5+CH3          ","CH2+H2=>CH3+H               ","CH2+H=>CH+H2                ",&
   "CH2=>CH+H                   ","CH2+CH2=>C2H2+H2            ","CH2+H=>CH3                  ","CH+C2H6=>C3H6+H             ",&
   "CH+C2H6=>C3H7               ","CH+H2=>CH2+H                ","CH+CH3=>C2H3+H              ","CH+CH2=>C2H2+H              ",&
   "CH+H2=>CH3                  ","CH+C2H3=>CH2+C2H2           ","C2H6+C2H3=>C2H5+C2H4        ","C2H6+C3H7=>C3H8+C2H5        ",&
   "C2H6+C3H5=>C3H6+C2H5        ","C2H6+H=>C2H5+H2             ","C2H6+H=>CH4+CH3             ","C2H6=>CH3+CH3               ",&
   "C2H6+CH=>C2H4+CH3           ","C2H6+CH2=>C2H5+CH3          ","C2H5+C2H5=>C2H6+C2H4        ","C2H5+C2H4=>C2H6+C2H3        ",&
   "C2H5+C2H2=>C2H6+C2H         ","C2H5+C2H=>C2H4+C2H2         ","C2H5+C3H8=>C2H6+C3H7        ","C2H5+C3H7=>C3H8+C2H4        ",&
   "C2H5+C3H7=>C3H6+C2H6        ","C2H5+C3H6=>C3H5+C2H6        ","C2H5+H2=>C2H6+H             ","C2H5+H=>CH3+CH3             ",&
   "C2H5+H=>C2H4+H2             ","C2H5+H=>C2H6                ","C2H5=>C2H4+H                ","C2H5+C2H3=>C2H4+C2H4        ",&
   "C2H4+H=>C2H3+H2             ","C2H4+H=>C2H5                ","C2H4+H2=>C2H5+H             ","C2H4=>C2H3+H                "/
  data reaction_sign(433:485) &
  /"C2H4+C3H6=>C3H5+C2H5        ","C2H4+C2H2=>C2H3+C2H3        ","C2H4+C3H6=>C2H3+C3H7        ","C2H4+C2H4=>C2H5+C2H3        ",&
   "C2H4+CH3=>C3H7              ","C2H4=>C2H2+H2               ","C2H4+C2H5=>C4H9             ","C2H4+H2=>C2H6               ",&
   "C2H4+CH2=>C3H6              ","C2H4+C4H9=>C3H6+C3H7        ","C2H3+C2H3=>C2H4+C2H2        ","C2H3+C3H8=>C2H4+C3H7        ",&
   "C2H3+C3H7=>C3H8+C2H2        ","C2H3+C3H7=>C3H6+C2H4        ","C2H3+C3H6=>C3H5+C2H4        ","C2H3+C3H5=>C3H6+C2H2        ",&
   "C2H3+H2=>C2H4+H             ","C2H3+H=>C2H2+H2             ","C2H3+H=>C2H4                ","C2H3=>C2H2+H                ",&
   "C2H2+H=>C2H3                ","C2H2+H2=>C2H4               ","C2H2+H2=>C2H3+H             ","C2H2+CH3=>C3H5              ",&
   "C2H2+C4H9=>C3H6+C3H5        ","C3H8+C3H5=>C3H6+C3H7        ","C3H8+H=>C3H7+H2             ","C3H8=>C2H5+CH3              ",&
   "C3H7+C3H7=>C3H6+C3H8        ","C3H7+C3H6=>C3H5+C3H8        ","C3H7+C3H5=>C3H6+C3H6        ","C3H7+H2=>C3H8+H             ",&
   "C3H7+H=>C3H6+H2             ","C3H7+H=>C3H8                ","C3H7+H=>CH3+C2H5            ","C3H7=>C3H6+H                ",&
   "C3H7=>C2H4+CH3              ","C3H6+C2H2=>C2H3+C3H5        ","C3H6+C3H6=>C3H7+C3H5        ","C3H6=>C3H5+H                ",&
   "C3H6+H=>C3H5+H2             ","C3H6+H=>C3H7                ","C3H6=>CH3+C2H3              ","C3H6+CH3=>C4H9              ",&
   "C3H5+H2=>C3H6+H             ","C3H5+H=>C3H6                ","C3H5=>C2H2+CH3              ","C4H9=>C2H4+C2H5             ",&
   "C4H9+CH2=>C2H4+C3H7         ","C4H9=>C3H6+CH3              ","C5H12=>CH3+C4H9             ","H2=>H+H                     ",&
   "H+H=>H2                     "/
  data bolsig_species(1:bolsig_species_max) &
  /"H2       ","C3H5     ","C2H4(V2) ","CH4      ","C3H6     ","CH2      ","C3H8     ","C2H4(V1) ","C3H4     ","C2H2     ",&
   "C2H6     ","C2H2(V5) ","C3H7     ","CH3      ","C2H2(V13)","C2H2(V2) ","C2H6(V24)","C2H5     ","C2H4     ","CH4(V13) ",&
   "C2H6(V13)","C3H8(V1) ","C3H8(V2) ","CH4(V24) ","CH       ","C2H3     "/
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
311 format(331x,50(1x,i9))
312 format(A3,1x,A28,1x,50(1x,A9))
313 format(i3,1x,A28,1x,50(1x,1pd9.2))
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
      write(ifile_unit,"(485(i3))",err=200) int(mrtm(i,:))
    enddo
    lerror = .false.
200 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file " // &
                                     "<qt_matrix.txt> (subroutine writer_save_qtplaskin)")
    close(ifile_unit)
    open(ifile_unit,file="qt_densities.txt",action="write",err=300)
    write(ifile_unit,"(1x,A14,50(111x,i2.2))",err=300) "Time_s", ( i, i = 1, species_max )
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
    write(ifile_unit,"(1x,A12,485(101x,i3.3))",err=500) "Time_s", ( i, i = 1, reactions_max )
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
      write(ifile_unit,"(1pe15.6,50(1pe13.4))") densav(0,2), densav(1:,2)
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
      write(ifile_unit,"(486(1pe13.4))") densav(0,2), rrt_loc(:)
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
  reac_source_local(05,001) = - reac_rate_local(001) 
  reac_source_local(43,001) = + reac_rate_local(001) 
  reac_source_local(05,002) = - reac_rate_local(002) 
  reac_source_local(36,002) = + reac_rate_local(002) 
  reac_source_local(20,003) = - reac_rate_local(003) 
  reac_source_local(33,003) = + reac_rate_local(003) 
  reac_source_local(20,004) = - reac_rate_local(004) 
  reac_source_local(34,004) = + reac_rate_local(004) 
  reac_source_local(16,005) = + reac_rate_local(005) 
  reac_source_local(35,005) = - reac_rate_local(005) 
  reac_source_local(04,006) = + reac_rate_local(006) 
  reac_source_local(35,006) = - reac_rate_local(006) 
  reac_source_local(19,007) = - reac_rate_local(007) 
  reac_source_local(21,007) = + reac_rate_local(007) 
  reac_source_local(19,008) = - reac_rate_local(008) 
  reac_source_local(31,008) = + reac_rate_local(008) 
  reac_source_local(19,009) = - reac_rate_local(009) 
  reac_source_local(29,009) = + reac_rate_local(009) 
  reac_source_local(12,010) = - reac_rate_local(010) 
  reac_source_local(37,010) = + reac_rate_local(010) 
  reac_source_local(12,011) = - reac_rate_local(011) 
  reac_source_local(40,011) = + reac_rate_local(011) 
  reac_source_local(09,012) = - reac_rate_local(012) 
  reac_source_local(42,012) = + reac_rate_local(012) 
  reac_source_local(05,013) = - reac_rate_local(013) 
  reac_source_local(13,013) = + reac_rate_local(013) 
  reac_source_local(24,013) = + reac_rate_local(013) 
  reac_source_local(02,014) = + reac_rate_local(014) 
  reac_source_local(05,014) = - reac_rate_local(014) 
  reac_source_local(10,014) = + reac_rate_local(014) 
  reac_source_local(02,015) = + reac_rate_local(015) 
  reac_source_local(05,015) = - reac_rate_local(015) 
  reac_source_local(13,015) = + reac_rate_local(015) 
  reac_source_local(46,015) = + reac_rate_local(015) 
  reac_source_local(10,016) = + reac_rate_local(016) 
  reac_source_local(13,016) = + reac_rate_local(016) 
  reac_source_local(24,016) = - reac_rate_local(016) 
  reac_source_local(02,017) = + reac_rate_local(017) 
  reac_source_local(24,017) = - reac_rate_local(017) 
  reac_source_local(46,017) = + reac_rate_local(017) 
  reac_source_local(10,018) = - reac_rate_local(018) 
  reac_source_local(13,018) = + reac_rate_local(018) 
  reac_source_local(46,018) = + reac_rate_local(018) 
  reac_source_local(01,019) = + reac_rate_local(019) 
  reac_source_local(05,019) = - reac_rate_local(019) 
  reac_source_local(25,019) = + reac_rate_local(019) 
  reac_source_local(01,020) = + reac_rate_local(020) 
  reac_source_local(05,020) = - reac_rate_local(020) 
  reac_source_local(13,020) = + reac_rate_local(020) 
  reac_source_local(48,020) = + reac_rate_local(020) 
  reac_source_local(01,021) = + reac_rate_local(021) 
  reac_source_local(02,021) = + reac_rate_local(021) 
  reac_source_local(05,021) = - reac_rate_local(021) 
  reac_source_local(27,021) = + reac_rate_local(021) 
  reac_source_local(01,022) = + reac_rate_local(022) 
  reac_source_local(02,022) = + reac_rate_local(022) 
  reac_source_local(05,022) = - reac_rate_local(022) 
  reac_source_local(13,022) = + reac_rate_local(022) 
  reac_source_local(28,022) = + reac_rate_local(022) 
  reac_source_local(01,023) = + reac_rate_local(023) 
  reac_source_local(24,023) = - reac_rate_local(023) 
  reac_source_local(48,023) = + reac_rate_local(023) 
  reac_source_local(01,024) = + reac_rate_local(024) 
  reac_source_local(13,024) = + reac_rate_local(024) 
  reac_source_local(24,024) = - reac_rate_local(024) 
  reac_source_local(27,024) = + reac_rate_local(024) 
  reac_source_local(01,025) = + reac_rate_local(025) 
  reac_source_local(02,025) = + reac_rate_local(025) 
  reac_source_local(24,025) = - reac_rate_local(025) 
  reac_source_local(28,025) = + reac_rate_local(025) 
  reac_source_local(01,026) = + reac_rate_local(026) 
  reac_source_local(10,026) = - reac_rate_local(026) 
  reac_source_local(27,026) = + reac_rate_local(026) 
  reac_source_local(01,027) = + reac_rate_local(027) 
  reac_source_local(10,027) = - reac_rate_local(027) 
  reac_source_local(13,027) = + reac_rate_local(027) 
  reac_source_local(28,027) = + reac_rate_local(027) 
  reac_source_local(01,028) = + reac_rate_local(028) 
  reac_source_local(28,028) = + reac_rate_local(028) 
  reac_source_local(46,028) = - reac_rate_local(028) 
  reac_source_local(13,029) = + reac_rate_local(029) 
  reac_source_local(20,029) = - reac_rate_local(029) 
  reac_source_local(38,029) = + reac_rate_local(029) 
  reac_source_local(02,030) = + reac_rate_local(030) 
  reac_source_local(20,030) = - reac_rate_local(030) 
  reac_source_local(35,030) = + reac_rate_local(030) 
  reac_source_local(02,031) = + reac_rate_local(031) 
  reac_source_local(13,031) = + reac_rate_local(031) 
  reac_source_local(20,031) = - reac_rate_local(031) 
  reac_source_local(50,031) = + reac_rate_local(031) 
  reac_source_local(02,032) = + reac_rate_local(032) * 2.d0
  reac_source_local(19,032) = + reac_rate_local(032) 
  reac_source_local(20,032) = - reac_rate_local(032) 
  reac_source_local(05,033) = + reac_rate_local(033) 
  reac_source_local(10,033) = + reac_rate_local(033) 
  reac_source_local(20,033) = - reac_rate_local(033) 
  reac_source_local(20,034) = - reac_rate_local(034) 
  reac_source_local(24,034) = + reac_rate_local(034) * 2.d0
  reac_source_local(13,035) = + reac_rate_local(035) 
  reac_source_local(35,035) = + reac_rate_local(035) 
  reac_source_local(38,035) = - reac_rate_local(035) 
  reac_source_local(02,036) = + reac_rate_local(036) 
  reac_source_local(38,036) = - reac_rate_local(036) 
  reac_source_local(50,036) = + reac_rate_local(036) 
  reac_source_local(13,037) = + reac_rate_local(037) * 2.d0
  reac_source_local(38,037) = - reac_rate_local(037) 
  reac_source_local(50,037) = + reac_rate_local(037) 
  reac_source_local(02,038) = + reac_rate_local(038) 
  reac_source_local(13,038) = + reac_rate_local(038) 
  reac_source_local(19,038) = + reac_rate_local(038) 
  reac_source_local(38,038) = - reac_rate_local(038) 
  reac_source_local(05,039) = + reac_rate_local(039) 
  reac_source_local(38,039) = - reac_rate_local(039) 
  reac_source_local(46,039) = + reac_rate_local(039) 
  reac_source_local(10,040) = + reac_rate_local(040) 
  reac_source_local(24,040) = + reac_rate_local(040) 
  reac_source_local(38,040) = - reac_rate_local(040) 
  reac_source_local(13,041) = + reac_rate_local(041) 
  reac_source_local(35,041) = - reac_rate_local(041) 
  reac_source_local(50,041) = + reac_rate_local(041) 
  reac_source_local(02,042) = + reac_rate_local(042) 
  reac_source_local(19,042) = + reac_rate_local(042) 
  reac_source_local(35,042) = - reac_rate_local(042) 
  reac_source_local(13,043) = + reac_rate_local(043) * 2.d0
  reac_source_local(19,043) = + reac_rate_local(043) 
  reac_source_local(35,043) = - reac_rate_local(043) 
  reac_source_local(24,044) = + reac_rate_local(044) 
  reac_source_local(35,044) = - reac_rate_local(044) 
  reac_source_local(46,044) = + reac_rate_local(044) 
  reac_source_local(10,045) = + reac_rate_local(045) * 2.d0
  reac_source_local(35,045) = - reac_rate_local(045) 
  reac_source_local(13,046) = + reac_rate_local(046) 
  reac_source_local(19,046) = + reac_rate_local(046) 
  reac_source_local(50,046) = - reac_rate_local(046) 
  reac_source_local(10,047) = + reac_rate_local(047) 
  reac_source_local(46,047) = + reac_rate_local(047) 
  reac_source_local(50,047) = - reac_rate_local(047) 
  reac_source_local(19,048) = - reac_rate_local(048) 
  reac_source_local(46,048) = + reac_rate_local(048) * 2.d0
  reac_source_local(12,049) = - reac_rate_local(049) 
  reac_source_local(13,049) = + reac_rate_local(049) 
  reac_source_local(23,049) = + reac_rate_local(049) 
  reac_source_local(02,050) = + reac_rate_local(050) 
  reac_source_local(09,050) = + reac_rate_local(050) 
  reac_source_local(12,050) = - reac_rate_local(050) 
  reac_source_local(02,051) = + reac_rate_local(051) * 2.d0
  reac_source_local(12,051) = - reac_rate_local(051) 
  reac_source_local(17,051) = + reac_rate_local(051) 
  reac_source_local(10,052) = + reac_rate_local(052) 
  reac_source_local(12,052) = - reac_rate_local(052) 
  reac_source_local(20,052) = + reac_rate_local(052) 
  reac_source_local(12,053) = - reac_rate_local(053) 
  reac_source_local(24,053) = + reac_rate_local(053) 
  reac_source_local(38,053) = + reac_rate_local(053) 
  reac_source_local(05,054) = + reac_rate_local(054) 
  reac_source_local(12,054) = - reac_rate_local(054) 
  reac_source_local(35,054) = + reac_rate_local(054) 
  reac_source_local(09,055) = + reac_rate_local(055) 
  reac_source_local(13,055) = + reac_rate_local(055) 
  reac_source_local(23,055) = - reac_rate_local(055) 
  reac_source_local(02,056) = + reac_rate_local(056) 
  reac_source_local(03,056) = + reac_rate_local(056) 
  reac_source_local(23,056) = - reac_rate_local(056) 
  reac_source_local(02,057) = + reac_rate_local(057) 
  reac_source_local(13,057) = + reac_rate_local(057) 
  reac_source_local(17,057) = + reac_rate_local(057) 
  reac_source_local(23,057) = - reac_rate_local(057) 
  reac_source_local(23,058) = - reac_rate_local(058) 
  reac_source_local(24,058) = + reac_rate_local(058) 
  reac_source_local(35,058) = + reac_rate_local(058) 
  reac_source_local(05,059) = + reac_rate_local(059) 
  reac_source_local(23,059) = - reac_rate_local(059) 
  reac_source_local(50,059) = + reac_rate_local(059) 
  reac_source_local(03,060) = + reac_rate_local(060) 
  reac_source_local(09,060) = - reac_rate_local(060) 
  reac_source_local(13,060) = + reac_rate_local(060) 
  reac_source_local(02,061) = + reac_rate_local(061) 
  reac_source_local(09,061) = - reac_rate_local(061) 
  reac_source_local(17,061) = + reac_rate_local(061) 
  reac_source_local(09,062) = - reac_rate_local(062) 
  reac_source_local(10,062) = + reac_rate_local(062) 
  reac_source_local(35,062) = + reac_rate_local(062) 
  reac_source_local(09,063) = - reac_rate_local(063) 
  reac_source_local(24,063) = + reac_rate_local(063) 
  reac_source_local(50,063) = + reac_rate_local(063) 
  reac_source_local(05,064) = + reac_rate_local(064) 
  reac_source_local(09,064) = - reac_rate_local(064) 
  reac_source_local(19,064) = + reac_rate_local(064) 
  reac_source_local(03,065) = - reac_rate_local(065) 
  reac_source_local(13,065) = + reac_rate_local(065) 
  reac_source_local(17,065) = + reac_rate_local(065) 
  reac_source_local(03,066) = - reac_rate_local(066) 
  reac_source_local(19,066) = + reac_rate_local(066) 
  reac_source_local(24,066) = + reac_rate_local(066) 
  reac_source_local(17,067) = - reac_rate_local(067) 
  reac_source_local(46,067) = + reac_rate_local(067) 
  reac_source_local(50,067) = + reac_rate_local(067) 
  reac_source_local(10,068) = + reac_rate_local(068) 
  reac_source_local(17,068) = - reac_rate_local(068) 
  reac_source_local(19,068) = + reac_rate_local(068) 
  reac_source_local(01,069) = + reac_rate_local(069) 
  reac_source_local(15,069) = + reac_rate_local(069) 
  reac_source_local(20,069) = - reac_rate_local(069) 
  reac_source_local(01,070) = + reac_rate_local(070) 
  reac_source_local(13,070) = + reac_rate_local(070) 
  reac_source_local(20,070) = - reac_rate_local(070) 
  reac_source_local(44,070) = + reac_rate_local(070) 
  reac_source_local(01,071) = + reac_rate_local(071) 
  reac_source_local(02,071) = + reac_rate_local(071) 
  reac_source_local(20,071) = - reac_rate_local(071) 
  reac_source_local(32,071) = + reac_rate_local(071) 
  reac_source_local(01,072) = + reac_rate_local(072) 
  reac_source_local(02,072) = + reac_rate_local(072) 
  reac_source_local(06,072) = + reac_rate_local(072) 
  reac_source_local(13,072) = + reac_rate_local(072) 
  reac_source_local(20,072) = - reac_rate_local(072) 
  reac_source_local(01,073) = + reac_rate_local(073) 
  reac_source_local(02,073) = + reac_rate_local(073) * 2.d0
  reac_source_local(14,073) = + reac_rate_local(073) 
  reac_source_local(20,073) = - reac_rate_local(073) 
  reac_source_local(01,074) = + reac_rate_local(074) 
  reac_source_local(20,074) = - reac_rate_local(074) 
  reac_source_local(24,074) = + reac_rate_local(074) 
  reac_source_local(48,074) = + reac_rate_local(074) 
  reac_source_local(01,075) = + reac_rate_local(075) 
  reac_source_local(05,075) = + reac_rate_local(075) 
  reac_source_local(20,075) = - reac_rate_local(075) 
  reac_source_local(27,075) = + reac_rate_local(075) 
  reac_source_local(01,076) = + reac_rate_local(076) 
  reac_source_local(38,076) = - reac_rate_local(076) 
  reac_source_local(44,076) = + reac_rate_local(076) 
  reac_source_local(01,077) = + reac_rate_local(077) 
  reac_source_local(13,077) = + reac_rate_local(077) 
  reac_source_local(32,077) = + reac_rate_local(077) 
  reac_source_local(38,077) = - reac_rate_local(077) 
  reac_source_local(01,078) = + reac_rate_local(078) 
  reac_source_local(02,078) = + reac_rate_local(078) 
  reac_source_local(06,078) = + reac_rate_local(078) 
  reac_source_local(38,078) = - reac_rate_local(078) 
  reac_source_local(01,079) = + reac_rate_local(079) 
  reac_source_local(02,079) = + reac_rate_local(079) 
  reac_source_local(13,079) = + reac_rate_local(079) 
  reac_source_local(14,079) = + reac_rate_local(079) 
  reac_source_local(38,079) = - reac_rate_local(079) 
  reac_source_local(01,080) = + reac_rate_local(080) 
  reac_source_local(10,080) = + reac_rate_local(080) 
  reac_source_local(38,080) = - reac_rate_local(080) 
  reac_source_local(48,080) = + reac_rate_local(080) 
  reac_source_local(01,081) = + reac_rate_local(081) 
  reac_source_local(24,081) = + reac_rate_local(081) 
  reac_source_local(27,081) = + reac_rate_local(081) 
  reac_source_local(38,081) = - reac_rate_local(081) 
  reac_source_local(01,082) = + reac_rate_local(082) 
  reac_source_local(05,082) = + reac_rate_local(082) 
  reac_source_local(28,082) = + reac_rate_local(082) 
  reac_source_local(38,082) = - reac_rate_local(082) 
  reac_source_local(01,083) = + reac_rate_local(083) 
  reac_source_local(32,083) = + reac_rate_local(083) 
  reac_source_local(35,083) = - reac_rate_local(083) 
  reac_source_local(01,084) = + reac_rate_local(084) 
  reac_source_local(06,084) = + reac_rate_local(084) 
  reac_source_local(13,084) = + reac_rate_local(084) 
  reac_source_local(35,084) = - reac_rate_local(084) 
  reac_source_local(01,085) = + reac_rate_local(085) 
  reac_source_local(35,085) = - reac_rate_local(085) 
  reac_source_local(46,085) = + reac_rate_local(085) 
  reac_source_local(48,085) = + reac_rate_local(085) 
  reac_source_local(01,086) = + reac_rate_local(086) 
  reac_source_local(10,086) = + reac_rate_local(086) 
  reac_source_local(27,086) = + reac_rate_local(086) 
  reac_source_local(35,086) = - reac_rate_local(086) 
  reac_source_local(01,087) = + reac_rate_local(087) 
  reac_source_local(24,087) = + reac_rate_local(087) 
  reac_source_local(28,087) = + reac_rate_local(087) 
  reac_source_local(35,087) = - reac_rate_local(087) 
  reac_source_local(01,088) = + reac_rate_local(088) 
  reac_source_local(06,088) = + reac_rate_local(088) 
  reac_source_local(50,088) = - reac_rate_local(088) 
  reac_source_local(01,089) = + reac_rate_local(089) 
  reac_source_local(13,089) = + reac_rate_local(089) 
  reac_source_local(14,089) = + reac_rate_local(089) 
  reac_source_local(50,089) = - reac_rate_local(089) 
  reac_source_local(01,090) = + reac_rate_local(090) 
  reac_source_local(27,090) = + reac_rate_local(090) 
  reac_source_local(46,090) = + reac_rate_local(090) 
  reac_source_local(50,090) = - reac_rate_local(090) 
  reac_source_local(01,091) = + reac_rate_local(091) 
  reac_source_local(10,091) = + reac_rate_local(091) 
  reac_source_local(28,091) = + reac_rate_local(091) 
  reac_source_local(50,091) = - reac_rate_local(091) 
  reac_source_local(01,092) = + reac_rate_local(092) 
  reac_source_local(07,092) = + reac_rate_local(092) 
  reac_source_local(19,092) = + reac_rate_local(092) 
  reac_source_local(50,092) = - reac_rate_local(092) 
  reac_source_local(01,093) = + reac_rate_local(093) 
  reac_source_local(14,093) = + reac_rate_local(093) 
  reac_source_local(19,093) = - reac_rate_local(093) 
  reac_source_local(01,094) = + reac_rate_local(094) 
  reac_source_local(19,094) = - reac_rate_local(094) 
  reac_source_local(28,094) = + reac_rate_local(094) 
  reac_source_local(46,094) = + reac_rate_local(094) 
  reac_source_local(01,095) = + reac_rate_local(095) 
  reac_source_local(12,095) = - reac_rate_local(095) 
  reac_source_local(22,095) = + reac_rate_local(095) 
  reac_source_local(01,096) = + reac_rate_local(096) 
  reac_source_local(12,096) = - reac_rate_local(096) 
  reac_source_local(13,096) = + reac_rate_local(096) 
  reac_source_local(30,096) = + reac_rate_local(096) 
  reac_source_local(01,097) = + reac_rate_local(097) 
  reac_source_local(02,097) = + reac_rate_local(097) 
  reac_source_local(12,097) = - reac_rate_local(097) 
  reac_source_local(49,097) = + reac_rate_local(097) 
  reac_source_local(01,098) = + reac_rate_local(098) 
  reac_source_local(02,098) = + reac_rate_local(098) 
  reac_source_local(08,098) = + reac_rate_local(098) 
  reac_source_local(12,098) = - reac_rate_local(098) 
  reac_source_local(13,098) = + reac_rate_local(098) 
  reac_source_local(01,099) = + reac_rate_local(099) 
  reac_source_local(02,099) = + reac_rate_local(099) * 2.d0
  reac_source_local(12,099) = - reac_rate_local(099) 
  reac_source_local(47,099) = + reac_rate_local(099) 
  reac_source_local(01,100) = + reac_rate_local(100) 
  reac_source_local(12,100) = - reac_rate_local(100) 
  reac_source_local(24,100) = + reac_rate_local(100) 
  reac_source_local(44,100) = + reac_rate_local(100) 
  reac_source_local(01,101) = + reac_rate_local(101) 
  reac_source_local(05,101) = + reac_rate_local(101) 
  reac_source_local(12,101) = - reac_rate_local(101) 
  reac_source_local(32,101) = + reac_rate_local(101) 
  reac_source_local(01,102) = + reac_rate_local(102) 
  reac_source_local(12,102) = - reac_rate_local(102) 
  reac_source_local(38,102) = + reac_rate_local(102) 
  reac_source_local(48,102) = + reac_rate_local(102) 
  reac_source_local(01,103) = + reac_rate_local(103) 
  reac_source_local(12,103) = - reac_rate_local(103) 
  reac_source_local(20,103) = + reac_rate_local(103) 
  reac_source_local(27,103) = + reac_rate_local(103) 
  reac_source_local(01,104) = + reac_rate_local(104) 
  reac_source_local(23,104) = - reac_rate_local(104) 
  reac_source_local(30,104) = + reac_rate_local(104) 
  reac_source_local(01,105) = + reac_rate_local(105) 
  reac_source_local(13,105) = + reac_rate_local(105) 
  reac_source_local(23,105) = - reac_rate_local(105) 
  reac_source_local(49,105) = + reac_rate_local(105) 
  reac_source_local(01,106) = + reac_rate_local(106) 
  reac_source_local(02,106) = + reac_rate_local(106) 
  reac_source_local(08,106) = + reac_rate_local(106) 
  reac_source_local(23,106) = - reac_rate_local(106) 
  reac_source_local(01,107) = + reac_rate_local(107) 
  reac_source_local(02,107) = + reac_rate_local(107) 
  reac_source_local(13,107) = + reac_rate_local(107) 
  reac_source_local(23,107) = - reac_rate_local(107) 
  reac_source_local(47,107) = + reac_rate_local(107) 
  reac_source_local(01,108) = + reac_rate_local(108) 
  reac_source_local(10,108) = + reac_rate_local(108) 
  reac_source_local(23,108) = - reac_rate_local(108) 
  reac_source_local(44,108) = + reac_rate_local(108) 
  reac_source_local(01,109) = + reac_rate_local(109) 
  reac_source_local(23,109) = - reac_rate_local(109) 
  reac_source_local(24,109) = + reac_rate_local(109) 
  reac_source_local(32,109) = + reac_rate_local(109) 
  reac_source_local(01,110) = + reac_rate_local(110) 
  reac_source_local(05,110) = + reac_rate_local(110) 
  reac_source_local(06,110) = + reac_rate_local(110) 
  reac_source_local(23,110) = - reac_rate_local(110) 
  reac_source_local(01,111) = + reac_rate_local(111) 
  reac_source_local(23,111) = - reac_rate_local(111) 
  reac_source_local(25,111) = + reac_rate_local(111) 
  reac_source_local(50,111) = + reac_rate_local(111) 
  reac_source_local(01,112) = + reac_rate_local(112) 
  reac_source_local(23,112) = - reac_rate_local(112) 
  reac_source_local(35,112) = + reac_rate_local(112) 
  reac_source_local(48,112) = + reac_rate_local(112) 
  reac_source_local(01,113) = + reac_rate_local(113) 
  reac_source_local(23,113) = - reac_rate_local(113) 
  reac_source_local(27,113) = + reac_rate_local(113) 
  reac_source_local(38,113) = + reac_rate_local(113) 
  reac_source_local(01,114) = + reac_rate_local(114) 
  reac_source_local(20,114) = + reac_rate_local(114) 
  reac_source_local(23,114) = - reac_rate_local(114) 
  reac_source_local(28,114) = + reac_rate_local(114) 
  reac_source_local(01,115) = + reac_rate_local(115) 
  reac_source_local(09,115) = - reac_rate_local(115) 
  reac_source_local(49,115) = + reac_rate_local(115) 
  reac_source_local(01,116) = + reac_rate_local(116) 
  reac_source_local(08,116) = + reac_rate_local(116) 
  reac_source_local(09,116) = - reac_rate_local(116) 
  reac_source_local(13,116) = + reac_rate_local(116) 
  reac_source_local(01,117) = + reac_rate_local(117) 
  reac_source_local(02,117) = + reac_rate_local(117) 
  reac_source_local(09,117) = - reac_rate_local(117) 
  reac_source_local(47,117) = + reac_rate_local(117) 
  reac_source_local(01,118) = + reac_rate_local(118) 
  reac_source_local(09,118) = - reac_rate_local(118) 
  reac_source_local(44,118) = + reac_rate_local(118) 
  reac_source_local(46,118) = + reac_rate_local(118) 
  reac_source_local(01,119) = + reac_rate_local(119) 
  reac_source_local(09,119) = - reac_rate_local(119) 
  reac_source_local(10,119) = + reac_rate_local(119) 
  reac_source_local(32,119) = + reac_rate_local(119) 
  reac_source_local(01,120) = + reac_rate_local(120) 
  reac_source_local(06,120) = + reac_rate_local(120) 
  reac_source_local(09,120) = - reac_rate_local(120) 
  reac_source_local(24,120) = + reac_rate_local(120) 
  reac_source_local(01,121) = + reac_rate_local(121) 
  reac_source_local(05,121) = + reac_rate_local(121) 
  reac_source_local(09,121) = - reac_rate_local(121) 
  reac_source_local(14,121) = + reac_rate_local(121) 
  reac_source_local(01,122) = + reac_rate_local(122) 
  reac_source_local(09,122) = - reac_rate_local(122) 
  reac_source_local(19,122) = + reac_rate_local(122) 
  reac_source_local(25,122) = + reac_rate_local(122) 
  reac_source_local(01,123) = + reac_rate_local(123) 
  reac_source_local(09,123) = - reac_rate_local(123) 
  reac_source_local(48,123) = + reac_rate_local(123) 
  reac_source_local(50,123) = + reac_rate_local(123) 
  reac_source_local(01,124) = + reac_rate_local(124) 
  reac_source_local(09,124) = - reac_rate_local(124) 
  reac_source_local(27,124) = + reac_rate_local(124) 
  reac_source_local(35,124) = + reac_rate_local(124) 
  reac_source_local(01,125) = + reac_rate_local(125) 
  reac_source_local(09,125) = - reac_rate_local(125) 
  reac_source_local(28,125) = + reac_rate_local(125) 
  reac_source_local(38,125) = + reac_rate_local(125) 
  reac_source_local(01,126) = + reac_rate_local(126) 
  reac_source_local(03,126) = - reac_rate_local(126) 
  reac_source_local(08,126) = + reac_rate_local(126) 
  reac_source_local(01,127) = + reac_rate_local(127) 
  reac_source_local(03,127) = - reac_rate_local(127) 
  reac_source_local(13,127) = + reac_rate_local(127) 
  reac_source_local(47,127) = + reac_rate_local(127) 
  reac_source_local(01,128) = + reac_rate_local(128) 
  reac_source_local(03,128) = - reac_rate_local(128) 
  reac_source_local(32,128) = + reac_rate_local(128) 
  reac_source_local(46,128) = + reac_rate_local(128) 
  reac_source_local(01,129) = + reac_rate_local(129) 
  reac_source_local(03,129) = - reac_rate_local(129) 
  reac_source_local(06,129) = + reac_rate_local(129) 
  reac_source_local(10,129) = + reac_rate_local(129) 
  reac_source_local(01,130) = + reac_rate_local(130) 
  reac_source_local(03,130) = - reac_rate_local(130) 
  reac_source_local(14,130) = + reac_rate_local(130) 
  reac_source_local(24,130) = + reac_rate_local(130) 
  reac_source_local(01,131) = + reac_rate_local(131) 
  reac_source_local(03,131) = - reac_rate_local(131) 
  reac_source_local(19,131) = + reac_rate_local(131) 
  reac_source_local(48,131) = + reac_rate_local(131) 
  reac_source_local(01,132) = + reac_rate_local(132) 
  reac_source_local(03,132) = - reac_rate_local(132) 
  reac_source_local(27,132) = + reac_rate_local(132) 
  reac_source_local(50,132) = + reac_rate_local(132) 
  reac_source_local(01,133) = + reac_rate_local(133) 
  reac_source_local(03,133) = - reac_rate_local(133) 
  reac_source_local(28,133) = + reac_rate_local(133) 
  reac_source_local(35,133) = + reac_rate_local(133) 
  reac_source_local(01,134) = + reac_rate_local(134) 
  reac_source_local(17,134) = - reac_rate_local(134) 
  reac_source_local(47,134) = + reac_rate_local(134) 
  reac_source_local(01,135) = + reac_rate_local(135) 
  reac_source_local(06,135) = + reac_rate_local(135) 
  reac_source_local(17,135) = - reac_rate_local(135) 
  reac_source_local(46,135) = + reac_rate_local(135) 
  reac_source_local(01,136) = + reac_rate_local(136) 
  reac_source_local(10,136) = + reac_rate_local(136) 
  reac_source_local(14,136) = + reac_rate_local(136) 
  reac_source_local(17,136) = - reac_rate_local(136) 
  reac_source_local(01,137) = + reac_rate_local(137) 
  reac_source_local(17,137) = - reac_rate_local(137) 
  reac_source_local(19,137) = + reac_rate_local(137) 
  reac_source_local(27,137) = + reac_rate_local(137) 
  reac_source_local(01,138) = + reac_rate_local(138) 
  reac_source_local(17,138) = - reac_rate_local(138) 
  reac_source_local(28,138) = + reac_rate_local(138) 
  reac_source_local(50,138) = + reac_rate_local(138) 
  reac_source_local(13,139) = + reac_rate_local(139) 
  reac_source_local(24,139) = + reac_rate_local(139) 
  reac_source_local(43,139) = - reac_rate_local(139) 
  reac_source_local(13,140) = + reac_rate_local(140) 
  reac_source_local(24,140) = + reac_rate_local(140) 
  reac_source_local(36,140) = - reac_rate_local(140) 
  reac_source_local(02,141) = + reac_rate_local(141) 
  reac_source_local(10,141) = + reac_rate_local(141) 
  reac_source_local(43,141) = - reac_rate_local(141) 
  reac_source_local(02,142) = + reac_rate_local(142) 
  reac_source_local(10,142) = + reac_rate_local(142) 
  reac_source_local(36,142) = - reac_rate_local(142) 
  reac_source_local(02,143) = + reac_rate_local(143) 
  reac_source_local(13,143) = + reac_rate_local(143) 
  reac_source_local(43,143) = - reac_rate_local(143) 
  reac_source_local(46,143) = + reac_rate_local(143) 
  reac_source_local(02,144) = + reac_rate_local(144) 
  reac_source_local(13,144) = + reac_rate_local(144) 
  reac_source_local(36,144) = - reac_rate_local(144) 
  reac_source_local(46,144) = + reac_rate_local(144) 
  reac_source_local(01,145) = + reac_rate_local(145) 
  reac_source_local(25,145) = + reac_rate_local(145) 
  reac_source_local(43,145) = - reac_rate_local(145) 
  reac_source_local(01,146) = + reac_rate_local(146) 
  reac_source_local(25,146) = + reac_rate_local(146) 
  reac_source_local(36,146) = - reac_rate_local(146) 
  reac_source_local(01,147) = + reac_rate_local(147) 
  reac_source_local(13,147) = + reac_rate_local(147) 
  reac_source_local(43,147) = - reac_rate_local(147) 
  reac_source_local(48,147) = + reac_rate_local(147) 
  reac_source_local(01,148) = + reac_rate_local(148) 
  reac_source_local(13,148) = + reac_rate_local(148) 
  reac_source_local(36,148) = - reac_rate_local(148) 
  reac_source_local(48,148) = + reac_rate_local(148) 
  reac_source_local(01,149) = + reac_rate_local(149) 
  reac_source_local(02,149) = + reac_rate_local(149) 
  reac_source_local(27,149) = + reac_rate_local(149) 
  reac_source_local(43,149) = - reac_rate_local(149) 
  reac_source_local(01,150) = + reac_rate_local(150) 
  reac_source_local(02,150) = + reac_rate_local(150) 
  reac_source_local(27,150) = + reac_rate_local(150) 
  reac_source_local(36,150) = - reac_rate_local(150) 
  reac_source_local(01,151) = + reac_rate_local(151) 
  reac_source_local(02,151) = + reac_rate_local(151) 
  reac_source_local(13,151) = + reac_rate_local(151) 
  reac_source_local(28,151) = + reac_rate_local(151) 
  reac_source_local(43,151) = - reac_rate_local(151) 
  reac_source_local(01,152) = + reac_rate_local(152) 
  reac_source_local(02,152) = + reac_rate_local(152) 
  reac_source_local(13,152) = + reac_rate_local(152) 
  reac_source_local(28,152) = + reac_rate_local(152) 
  reac_source_local(36,152) = - reac_rate_local(152) 
  reac_source_local(13,153) = + reac_rate_local(153) 
  reac_source_local(33,153) = - reac_rate_local(153) 
  reac_source_local(38,153) = + reac_rate_local(153) 
  reac_source_local(13,154) = + reac_rate_local(154) 
  reac_source_local(34,154) = - reac_rate_local(154) 
  reac_source_local(38,154) = + reac_rate_local(154) 
  reac_source_local(02,155) = + reac_rate_local(155) 
  reac_source_local(33,155) = - reac_rate_local(155) 
  reac_source_local(35,155) = + reac_rate_local(155) 
  reac_source_local(02,156) = + reac_rate_local(156) 
  reac_source_local(34,156) = - reac_rate_local(156) 
  reac_source_local(35,156) = + reac_rate_local(156) 
  reac_source_local(02,157) = + reac_rate_local(157) 
  reac_source_local(13,157) = + reac_rate_local(157) 
  reac_source_local(33,157) = - reac_rate_local(157) 
  reac_source_local(50,157) = + reac_rate_local(157) 
  reac_source_local(02,158) = + reac_rate_local(158) 
  reac_source_local(13,158) = + reac_rate_local(158) 
  reac_source_local(34,158) = - reac_rate_local(158) 
  reac_source_local(50,158) = + reac_rate_local(158) 
  reac_source_local(02,159) = + reac_rate_local(159) * 2.d0
  reac_source_local(19,159) = + reac_rate_local(159) 
  reac_source_local(33,159) = - reac_rate_local(159) 
  reac_source_local(02,160) = + reac_rate_local(160) * 2.d0
  reac_source_local(19,160) = + reac_rate_local(160) 
  reac_source_local(34,160) = - reac_rate_local(160) 
  reac_source_local(05,161) = + reac_rate_local(161) 
  reac_source_local(10,161) = + reac_rate_local(161) 
  reac_source_local(33,161) = - reac_rate_local(161) 
  reac_source_local(05,162) = + reac_rate_local(162) 
  reac_source_local(10,162) = + reac_rate_local(162) 
  reac_source_local(34,162) = - reac_rate_local(162) 
  reac_source_local(24,163) = + reac_rate_local(163) * 2.d0
  reac_source_local(33,163) = - reac_rate_local(163) 
  reac_source_local(24,164) = + reac_rate_local(164) * 2.d0
  reac_source_local(34,164) = - reac_rate_local(164) 
  reac_source_local(13,165) = + reac_rate_local(165) 
  reac_source_local(16,165) = - reac_rate_local(165) 
  reac_source_local(50,165) = + reac_rate_local(165) 
  reac_source_local(04,166) = - reac_rate_local(166) 
  reac_source_local(13,166) = + reac_rate_local(166) 
  reac_source_local(50,166) = + reac_rate_local(166) 
  reac_source_local(02,167) = + reac_rate_local(167) 
  reac_source_local(16,167) = - reac_rate_local(167) 
  reac_source_local(19,167) = + reac_rate_local(167) 
  reac_source_local(02,168) = + reac_rate_local(168) 
  reac_source_local(04,168) = - reac_rate_local(168) 
  reac_source_local(19,168) = + reac_rate_local(168) 
  reac_source_local(13,169) = + reac_rate_local(169) * 2.d0
  reac_source_local(16,169) = - reac_rate_local(169) 
  reac_source_local(19,169) = + reac_rate_local(169) 
  reac_source_local(04,170) = - reac_rate_local(170) 
  reac_source_local(13,170) = + reac_rate_local(170) * 2.d0
  reac_source_local(19,170) = + reac_rate_local(170) 
  reac_source_local(16,171) = - reac_rate_local(171) 
  reac_source_local(24,171) = + reac_rate_local(171) 
  reac_source_local(46,171) = + reac_rate_local(171) 
  reac_source_local(04,172) = - reac_rate_local(172) 
  reac_source_local(24,172) = + reac_rate_local(172) 
  reac_source_local(46,172) = + reac_rate_local(172) 
  reac_source_local(10,173) = + reac_rate_local(173) * 2.d0
  reac_source_local(16,173) = - reac_rate_local(173) 
  reac_source_local(04,174) = - reac_rate_local(174) 
  reac_source_local(10,174) = + reac_rate_local(174) * 2.d0
  reac_source_local(21,175) = - reac_rate_local(175) 
  reac_source_local(46,175) = + reac_rate_local(175) * 2.d0
  reac_source_local(31,176) = - reac_rate_local(176) 
  reac_source_local(46,176) = + reac_rate_local(176) * 2.d0
  reac_source_local(29,177) = - reac_rate_local(177) 
  reac_source_local(46,177) = + reac_rate_local(177) * 2.d0
  reac_source_local(13,178) = + reac_rate_local(178) 
  reac_source_local(23,178) = + reac_rate_local(178) 
  reac_source_local(37,178) = - reac_rate_local(178) 
  reac_source_local(13,179) = + reac_rate_local(179) 
  reac_source_local(23,179) = + reac_rate_local(179) 
  reac_source_local(40,179) = - reac_rate_local(179) 
  reac_source_local(02,180) = + reac_rate_local(180) 
  reac_source_local(09,180) = + reac_rate_local(180) 
  reac_source_local(37,180) = - reac_rate_local(180) 
  reac_source_local(02,181) = + reac_rate_local(181) 
  reac_source_local(09,181) = + reac_rate_local(181) 
  reac_source_local(40,181) = - reac_rate_local(181) 
  reac_source_local(02,182) = + reac_rate_local(182) * 2.d0
  reac_source_local(17,182) = + reac_rate_local(182) 
  reac_source_local(37,182) = - reac_rate_local(182) 
  reac_source_local(02,183) = + reac_rate_local(183) * 2.d0
  reac_source_local(17,183) = + reac_rate_local(183) 
  reac_source_local(40,183) = - reac_rate_local(183) 
  reac_source_local(10,184) = + reac_rate_local(184) 
  reac_source_local(20,184) = + reac_rate_local(184) 
  reac_source_local(37,184) = - reac_rate_local(184) 
  reac_source_local(10,185) = + reac_rate_local(185) 
  reac_source_local(20,185) = + reac_rate_local(185) 
  reac_source_local(40,185) = - reac_rate_local(185) 
  reac_source_local(24,186) = + reac_rate_local(186) 
  reac_source_local(37,186) = - reac_rate_local(186) 
  reac_source_local(38,186) = + reac_rate_local(186) 
  reac_source_local(24,187) = + reac_rate_local(187) 
  reac_source_local(38,187) = + reac_rate_local(187) 
  reac_source_local(40,187) = - reac_rate_local(187) 
  reac_source_local(05,188) = + reac_rate_local(188) 
  reac_source_local(35,188) = + reac_rate_local(188) 
  reac_source_local(37,188) = - reac_rate_local(188) 
  reac_source_local(05,189) = + reac_rate_local(189) 
  reac_source_local(35,189) = + reac_rate_local(189) 
  reac_source_local(40,189) = - reac_rate_local(189) 
  reac_source_local(01,190) = + reac_rate_local(190) 
  reac_source_local(15,190) = + reac_rate_local(190) 
  reac_source_local(33,190) = - reac_rate_local(190) 
  reac_source_local(01,191) = + reac_rate_local(191) 
  reac_source_local(15,191) = + reac_rate_local(191) 
  reac_source_local(34,191) = - reac_rate_local(191) 
  reac_source_local(01,192) = + reac_rate_local(192) 
  reac_source_local(13,192) = + reac_rate_local(192) 
  reac_source_local(33,192) = - reac_rate_local(192) 
  reac_source_local(44,192) = + reac_rate_local(192) 
  reac_source_local(01,193) = + reac_rate_local(193) 
  reac_source_local(13,193) = + reac_rate_local(193) 
  reac_source_local(34,193) = - reac_rate_local(193) 
  reac_source_local(44,193) = + reac_rate_local(193) 
  reac_source_local(01,194) = + reac_rate_local(194) 
  reac_source_local(02,194) = + reac_rate_local(194) 
  reac_source_local(32,194) = + reac_rate_local(194) 
  reac_source_local(33,194) = - reac_rate_local(194) 
  reac_source_local(01,195) = + reac_rate_local(195) 
  reac_source_local(02,195) = + reac_rate_local(195) 
  reac_source_local(32,195) = + reac_rate_local(195) 
  reac_source_local(34,195) = - reac_rate_local(195) 
  reac_source_local(01,196) = + reac_rate_local(196) 
  reac_source_local(02,196) = + reac_rate_local(196) 
  reac_source_local(06,196) = + reac_rate_local(196) 
  reac_source_local(13,196) = + reac_rate_local(196) 
  reac_source_local(33,196) = - reac_rate_local(196) 
  reac_source_local(01,197) = + reac_rate_local(197) 
  reac_source_local(02,197) = + reac_rate_local(197) 
  reac_source_local(06,197) = + reac_rate_local(197) 
  reac_source_local(13,197) = + reac_rate_local(197) 
  reac_source_local(34,197) = - reac_rate_local(197) 
  reac_source_local(01,198) = + reac_rate_local(198) 
  reac_source_local(02,198) = + reac_rate_local(198) * 2.d0
  reac_source_local(14,198) = + reac_rate_local(198) 
  reac_source_local(33,198) = - reac_rate_local(198) 
  reac_source_local(01,199) = + reac_rate_local(199) 
  reac_source_local(02,199) = + reac_rate_local(199) * 2.d0
  reac_source_local(14,199) = + reac_rate_local(199) 
  reac_source_local(34,199) = - reac_rate_local(199) 
  reac_source_local(01,200) = + reac_rate_local(200) 
  reac_source_local(24,200) = + reac_rate_local(200) 
  reac_source_local(33,200) = - reac_rate_local(200) 
  reac_source_local(48,200) = + reac_rate_local(200) 
  reac_source_local(01,201) = + reac_rate_local(201) 
  reac_source_local(24,201) = + reac_rate_local(201) 
  reac_source_local(34,201) = - reac_rate_local(201) 
  reac_source_local(48,201) = + reac_rate_local(201) 
  reac_source_local(01,202) = + reac_rate_local(202) 
  reac_source_local(05,202) = + reac_rate_local(202) 
  reac_source_local(27,202) = + reac_rate_local(202) 
  reac_source_local(33,202) = - reac_rate_local(202) 
  reac_source_local(01,203) = + reac_rate_local(203) 
  reac_source_local(05,203) = + reac_rate_local(203) 
  reac_source_local(27,203) = + reac_rate_local(203) 
  reac_source_local(34,203) = - reac_rate_local(203) 
  reac_source_local(01,204) = + reac_rate_local(204) 
  reac_source_local(16,204) = - reac_rate_local(204) 
  reac_source_local(32,204) = + reac_rate_local(204) 
  reac_source_local(01,205) = + reac_rate_local(205) 
  reac_source_local(04,205) = - reac_rate_local(205) 
  reac_source_local(32,205) = + reac_rate_local(205) 
  reac_source_local(01,206) = + reac_rate_local(206) 
  reac_source_local(06,206) = + reac_rate_local(206) 
  reac_source_local(13,206) = + reac_rate_local(206) 
  reac_source_local(16,206) = - reac_rate_local(206) 
  reac_source_local(01,207) = + reac_rate_local(207) 
  reac_source_local(04,207) = - reac_rate_local(207) 
  reac_source_local(06,207) = + reac_rate_local(207) 
  reac_source_local(13,207) = + reac_rate_local(207) 
  reac_source_local(01,208) = + reac_rate_local(208) 
  reac_source_local(16,208) = - reac_rate_local(208) 
  reac_source_local(46,208) = + reac_rate_local(208) 
  reac_source_local(48,208) = + reac_rate_local(208) 
  reac_source_local(01,209) = + reac_rate_local(209) 
  reac_source_local(04,209) = - reac_rate_local(209) 
  reac_source_local(46,209) = + reac_rate_local(209) 
  reac_source_local(48,209) = + reac_rate_local(209) 
  reac_source_local(01,210) = + reac_rate_local(210) 
  reac_source_local(10,210) = + reac_rate_local(210) 
  reac_source_local(16,210) = - reac_rate_local(210) 
  reac_source_local(27,210) = + reac_rate_local(210) 
  reac_source_local(01,211) = + reac_rate_local(211) 
  reac_source_local(04,211) = - reac_rate_local(211) 
  reac_source_local(10,211) = + reac_rate_local(211) 
  reac_source_local(27,211) = + reac_rate_local(211) 
  reac_source_local(01,212) = + reac_rate_local(212) 
  reac_source_local(16,212) = - reac_rate_local(212) 
  reac_source_local(24,212) = + reac_rate_local(212) 
  reac_source_local(28,212) = + reac_rate_local(212) 
  reac_source_local(01,213) = + reac_rate_local(213) 
  reac_source_local(04,213) = - reac_rate_local(213) 
  reac_source_local(24,213) = + reac_rate_local(213) 
  reac_source_local(28,213) = + reac_rate_local(213) 
  reac_source_local(01,214) = + reac_rate_local(214) 
  reac_source_local(14,214) = + reac_rate_local(214) 
  reac_source_local(21,214) = - reac_rate_local(214) 
  reac_source_local(01,215) = + reac_rate_local(215) 
  reac_source_local(14,215) = + reac_rate_local(215) 
  reac_source_local(31,215) = - reac_rate_local(215) 
  reac_source_local(01,216) = + reac_rate_local(216) 
  reac_source_local(14,216) = + reac_rate_local(216) 
  reac_source_local(29,216) = - reac_rate_local(216) 
  reac_source_local(01,217) = + reac_rate_local(217) 
  reac_source_local(21,217) = - reac_rate_local(217) 
  reac_source_local(28,217) = + reac_rate_local(217) 
  reac_source_local(46,217) = + reac_rate_local(217) 
  reac_source_local(01,218) = + reac_rate_local(218) 
  reac_source_local(28,218) = + reac_rate_local(218) 
  reac_source_local(31,218) = - reac_rate_local(218) 
  reac_source_local(46,218) = + reac_rate_local(218) 
  reac_source_local(01,219) = + reac_rate_local(219) 
  reac_source_local(28,219) = + reac_rate_local(219) 
  reac_source_local(29,219) = - reac_rate_local(219) 
  reac_source_local(46,219) = + reac_rate_local(219) 
  reac_source_local(01,220) = + reac_rate_local(220) 
  reac_source_local(22,220) = + reac_rate_local(220) 
  reac_source_local(37,220) = - reac_rate_local(220) 
  reac_source_local(01,221) = + reac_rate_local(221) 
  reac_source_local(22,221) = + reac_rate_local(221) 
  reac_source_local(40,221) = - reac_rate_local(221) 
  reac_source_local(01,222) = + reac_rate_local(222) 
  reac_source_local(13,222) = + reac_rate_local(222) 
  reac_source_local(30,222) = + reac_rate_local(222) 
  reac_source_local(37,222) = - reac_rate_local(222) 
  reac_source_local(01,223) = + reac_rate_local(223) 
  reac_source_local(13,223) = + reac_rate_local(223) 
  reac_source_local(30,223) = + reac_rate_local(223) 
  reac_source_local(40,223) = - reac_rate_local(223) 
  reac_source_local(01,224) = + reac_rate_local(224) 
  reac_source_local(02,224) = + reac_rate_local(224) 
  reac_source_local(37,224) = - reac_rate_local(224) 
  reac_source_local(49,224) = + reac_rate_local(224) 
  reac_source_local(01,225) = + reac_rate_local(225) 
  reac_source_local(02,225) = + reac_rate_local(225) 
  reac_source_local(40,225) = - reac_rate_local(225) 
  reac_source_local(49,225) = + reac_rate_local(225) 
  reac_source_local(01,226) = + reac_rate_local(226) 
  reac_source_local(02,226) = + reac_rate_local(226) 
  reac_source_local(08,226) = + reac_rate_local(226) 
  reac_source_local(13,226) = + reac_rate_local(226) 
  reac_source_local(37,226) = - reac_rate_local(226) 
  reac_source_local(01,227) = + reac_rate_local(227) 
  reac_source_local(02,227) = + reac_rate_local(227) 
  reac_source_local(08,227) = + reac_rate_local(227) 
  reac_source_local(13,227) = + reac_rate_local(227) 
  reac_source_local(40,227) = - reac_rate_local(227) 
  reac_source_local(01,228) = + reac_rate_local(228) 
  reac_source_local(02,228) = + reac_rate_local(228) * 2.d0
  reac_source_local(37,228) = - reac_rate_local(228) 
  reac_source_local(47,228) = + reac_rate_local(228) 
  reac_source_local(01,229) = + reac_rate_local(229) 
  reac_source_local(02,229) = + reac_rate_local(229) * 2.d0
  reac_source_local(40,229) = - reac_rate_local(229) 
  reac_source_local(47,229) = + reac_rate_local(229) 
  reac_source_local(01,230) = + reac_rate_local(230) 
  reac_source_local(24,230) = + reac_rate_local(230) 
  reac_source_local(37,230) = - reac_rate_local(230) 
  reac_source_local(44,230) = + reac_rate_local(230) 
  reac_source_local(01,231) = + reac_rate_local(231) 
  reac_source_local(24,231) = + reac_rate_local(231) 
  reac_source_local(40,231) = - reac_rate_local(231) 
  reac_source_local(44,231) = + reac_rate_local(231) 
  reac_source_local(01,232) = + reac_rate_local(232) 
  reac_source_local(05,232) = + reac_rate_local(232) 
  reac_source_local(32,232) = + reac_rate_local(232) 
  reac_source_local(37,232) = - reac_rate_local(232) 
  reac_source_local(01,233) = + reac_rate_local(233) 
  reac_source_local(05,233) = + reac_rate_local(233) 
  reac_source_local(32,233) = + reac_rate_local(233) 
  reac_source_local(40,233) = - reac_rate_local(233) 
  reac_source_local(01,234) = + reac_rate_local(234) 
  reac_source_local(37,234) = - reac_rate_local(234) 
  reac_source_local(38,234) = + reac_rate_local(234) 
  reac_source_local(48,234) = + reac_rate_local(234) 
  reac_source_local(01,235) = + reac_rate_local(235) 
  reac_source_local(38,235) = + reac_rate_local(235) 
  reac_source_local(40,235) = - reac_rate_local(235) 
  reac_source_local(48,235) = + reac_rate_local(235) 
  reac_source_local(01,236) = + reac_rate_local(236) 
  reac_source_local(20,236) = + reac_rate_local(236) 
  reac_source_local(27,236) = + reac_rate_local(236) 
  reac_source_local(37,236) = - reac_rate_local(236) 
  reac_source_local(01,237) = + reac_rate_local(237) 
  reac_source_local(20,237) = + reac_rate_local(237) 
  reac_source_local(27,237) = + reac_rate_local(237) 
  reac_source_local(40,237) = - reac_rate_local(237) 
  reac_source_local(02,238) = - reac_rate_local(238) 
  reac_source_local(13,238) = + reac_rate_local(238) * 2.d0
  reac_source_local(01,239) = + reac_rate_local(239) 
  reac_source_local(02,239) = - reac_rate_local(239) 
  reac_source_local(45,239) = + reac_rate_local(239) 
  reac_source_local(01,240) = - reac_rate_local(240) 
  reac_source_local(05,240) = + reac_rate_local(240) 
  reac_source_local(25,240) = - reac_rate_local(240) 
  reac_source_local(01,241) = - reac_rate_local(241) 
  reac_source_local(24,241) = + reac_rate_local(241) 
  reac_source_local(48,241) = - reac_rate_local(241) 
  reac_source_local(01,242) = - reac_rate_local(242) 
  reac_source_local(10,242) = + reac_rate_local(242) 
  reac_source_local(27,242) = - reac_rate_local(242) 
  reac_source_local(01,243) = - reac_rate_local(243) 
  reac_source_local(28,243) = - reac_rate_local(243) 
  reac_source_local(46,243) = + reac_rate_local(243) 
  reac_source_local(01,244) = - reac_rate_local(244) 
  reac_source_local(15,244) = - reac_rate_local(244) 
  reac_source_local(20,244) = + reac_rate_local(244) 
  reac_source_local(01,245) = - reac_rate_local(245) 
  reac_source_local(38,245) = + reac_rate_local(245) 
  reac_source_local(44,245) = - reac_rate_local(245) 
  reac_source_local(01,246) = - reac_rate_local(246) 
  reac_source_local(32,246) = - reac_rate_local(246) 
  reac_source_local(35,246) = + reac_rate_local(246) 
  reac_source_local(01,247) = - reac_rate_local(247) 
  reac_source_local(06,247) = - reac_rate_local(247) 
  reac_source_local(50,247) = + reac_rate_local(247) 
  reac_source_local(01,248) = - reac_rate_local(248) 
  reac_source_local(14,248) = - reac_rate_local(248) 
  reac_source_local(19,248) = + reac_rate_local(248) 
  reac_source_local(01,249) = - reac_rate_local(249) 
  reac_source_local(12,249) = + reac_rate_local(249) 
  reac_source_local(22,249) = - reac_rate_local(249) 
  reac_source_local(01,250) = - reac_rate_local(250) 
  reac_source_local(23,250) = + reac_rate_local(250) 
  reac_source_local(30,250) = - reac_rate_local(250) 
  reac_source_local(01,251) = - reac_rate_local(251) 
  reac_source_local(09,251) = + reac_rate_local(251) 
  reac_source_local(49,251) = - reac_rate_local(251) 
  reac_source_local(01,252) = - reac_rate_local(252) 
  reac_source_local(03,252) = + reac_rate_local(252) 
  reac_source_local(08,252) = - reac_rate_local(252) 
  reac_source_local(01,253) = - reac_rate_local(253) 
  reac_source_local(17,253) = + reac_rate_local(253) 
  reac_source_local(47,253) = - reac_rate_local(253) 
  reac_source_local(01,254) = - reac_rate_local(254) 
  reac_source_local(07,254) = - reac_rate_local(254) 
  reac_source_local(13,254) = + reac_rate_local(254) 
  reac_source_local(01,255) = - reac_rate_local(255) 
  reac_source_local(02,255) = + reac_rate_local(255) 
  reac_source_local(45,255) = - reac_rate_local(255) 
  reac_source_local(01,256) = - reac_rate_local(256) 
  reac_source_local(13,256) = + reac_rate_local(256) * 2.d0
  reac_source_local(24,256) = + reac_rate_local(256) 
  reac_source_local(41,256) = - reac_rate_local(256) 
  reac_source_local(01,257) = - reac_rate_local(257) 
  reac_source_local(02,257) = + reac_rate_local(257) 
  reac_source_local(10,257) = + reac_rate_local(257) 
  reac_source_local(13,257) = + reac_rate_local(257) 
  reac_source_local(41,257) = - reac_rate_local(257) 
  reac_source_local(01,258) = - reac_rate_local(258) 
  reac_source_local(13,258) = + reac_rate_local(258) 
  reac_source_local(24,258) = + reac_rate_local(258) 
  reac_source_local(25,258) = - reac_rate_local(258) 
  reac_source_local(01,259) = - reac_rate_local(259) 
  reac_source_local(10,259) = + reac_rate_local(259) 
  reac_source_local(13,259) = + reac_rate_local(259) * 2.d0
  reac_source_local(25,259) = - reac_rate_local(259) 
  reac_source_local(01,260) = - reac_rate_local(260) 
  reac_source_local(02,260) = + reac_rate_local(260) 
  reac_source_local(13,260) = + reac_rate_local(260) 
  reac_source_local(25,260) = - reac_rate_local(260) 
  reac_source_local(46,260) = + reac_rate_local(260) 
  reac_source_local(01,261) = - reac_rate_local(261) 
  reac_source_local(10,261) = + reac_rate_local(261) 
  reac_source_local(13,261) = + reac_rate_local(261) 
  reac_source_local(48,261) = - reac_rate_local(261) 
  reac_source_local(01,262) = - reac_rate_local(262) 
  reac_source_local(02,262) = + reac_rate_local(262) 
  reac_source_local(46,262) = + reac_rate_local(262) 
  reac_source_local(48,262) = - reac_rate_local(262) 
  reac_source_local(01,263) = - reac_rate_local(263) 
  reac_source_local(13,263) = + reac_rate_local(263) 
  reac_source_local(27,263) = - reac_rate_local(263) 
  reac_source_local(46,263) = + reac_rate_local(263) 
  reac_source_local(01,264) = - reac_rate_local(264) 
  reac_source_local(13,264) = + reac_rate_local(264) 
  reac_source_local(15,264) = - reac_rate_local(264) 
  reac_source_local(38,264) = + reac_rate_local(264) 
  reac_source_local(01,265) = - reac_rate_local(265) 
  reac_source_local(13,265) = + reac_rate_local(265) * 2.d0
  reac_source_local(15,265) = - reac_rate_local(265) 
  reac_source_local(35,265) = + reac_rate_local(265) 
  reac_source_local(01,266) = - reac_rate_local(266) 
  reac_source_local(13,266) = + reac_rate_local(266) 
  reac_source_local(35,266) = + reac_rate_local(266) 
  reac_source_local(44,266) = - reac_rate_local(266) 
  reac_source_local(01,267) = - reac_rate_local(267) 
  reac_source_local(13,267) = + reac_rate_local(267) * 2.d0
  reac_source_local(44,267) = - reac_rate_local(267) 
  reac_source_local(50,267) = + reac_rate_local(267) 
  reac_source_local(01,268) = - reac_rate_local(268) 
  reac_source_local(02,268) = + reac_rate_local(268) 
  reac_source_local(13,268) = + reac_rate_local(268) 
  reac_source_local(19,268) = + reac_rate_local(268) 
  reac_source_local(44,268) = - reac_rate_local(268) 
  reac_source_local(01,269) = - reac_rate_local(269) 
  reac_source_local(13,269) = + reac_rate_local(269) * 3.d0
  reac_source_local(19,269) = + reac_rate_local(269) 
  reac_source_local(44,269) = - reac_rate_local(269) 
  reac_source_local(01,270) = - reac_rate_local(270) 
  reac_source_local(10,270) = + reac_rate_local(270) 
  reac_source_local(24,270) = + reac_rate_local(270) 
  reac_source_local(44,270) = - reac_rate_local(270) 
  reac_source_local(01,271) = - reac_rate_local(271) 
  reac_source_local(13,271) = + reac_rate_local(271) 
  reac_source_local(32,271) = - reac_rate_local(271) 
  reac_source_local(50,271) = + reac_rate_local(271) 
  reac_source_local(01,272) = - reac_rate_local(272) 
  reac_source_local(13,272) = + reac_rate_local(272) * 2.d0
  reac_source_local(19,272) = + reac_rate_local(272) 
  reac_source_local(32,272) = - reac_rate_local(272) 
  reac_source_local(01,273) = - reac_rate_local(273) 
  reac_source_local(06,273) = - reac_rate_local(273) 
  reac_source_local(13,273) = + reac_rate_local(273) 
  reac_source_local(19,273) = + reac_rate_local(273) 
  reac_source_local(01,274) = - reac_rate_local(274) 
  reac_source_local(14,274) = - reac_rate_local(274) 
  reac_source_local(46,274) = + reac_rate_local(274) * 2.d0
  reac_source_local(01,275) = - reac_rate_local(275) 
  reac_source_local(13,275) = + reac_rate_local(275) * 2.d0
  reac_source_local(45,275) = - reac_rate_local(275) 
  reac_source_local(01,276) = - reac_rate_local(276) 
  reac_source_local(02,276) = + reac_rate_local(276) 
  reac_source_local(13,276) = + reac_rate_local(276) 
  reac_source_local(26,276) = - reac_rate_local(276) 
  reac_source_local(02,277) = - reac_rate_local(277) 
  reac_source_local(13,277) = + reac_rate_local(277) 
  reac_source_local(26,277) = + reac_rate_local(277) 
  reac_source_local(45,277) = - reac_rate_local(277) 
  reac_source_local(05,278) = + reac_rate_local(278) 
  reac_source_local(10,278) = - reac_rate_local(278) 
  reac_source_local(41,278) = - reac_rate_local(278) 
  reac_source_local(48,278) = + reac_rate_local(278) 
  reac_source_local(05,279) = + reac_rate_local(279) 
  reac_source_local(27,279) = + reac_rate_local(279) 
  reac_source_local(41,279) = - reac_rate_local(279) 
  reac_source_local(46,279) = - reac_rate_local(279) 
  reac_source_local(02,280) = + reac_rate_local(280) 
  reac_source_local(05,280) = + reac_rate_local(280) 
  reac_source_local(20,280) = - reac_rate_local(280) 
  reac_source_local(41,280) = - reac_rate_local(280) 
  reac_source_local(44,280) = + reac_rate_local(280) 
  reac_source_local(05,281) = + reac_rate_local(281) 
  reac_source_local(35,281) = - reac_rate_local(281) 
  reac_source_local(41,281) = - reac_rate_local(281) 
  reac_source_local(44,281) = + reac_rate_local(281) 
  reac_source_local(05,282) = + reac_rate_local(282) 
  reac_source_local(06,282) = + reac_rate_local(282) 
  reac_source_local(19,282) = - reac_rate_local(282) 
  reac_source_local(41,282) = - reac_rate_local(282) 
  reac_source_local(02,283) = + reac_rate_local(283) 
  reac_source_local(13,283) = - reac_rate_local(283) 
  reac_source_local(25,283) = + reac_rate_local(283) 
  reac_source_local(41,283) = - reac_rate_local(283) 
  reac_source_local(05,284) = - reac_rate_local(284) 
  reac_source_local(24,284) = + reac_rate_local(284) 
  reac_source_local(25,284) = - reac_rate_local(284) 
  reac_source_local(41,284) = + reac_rate_local(284) 
  reac_source_local(02,285) = + reac_rate_local(285) 
  reac_source_local(05,285) = + reac_rate_local(285) 
  reac_source_local(20,285) = - reac_rate_local(285) 
  reac_source_local(25,285) = - reac_rate_local(285) 
  reac_source_local(32,285) = + reac_rate_local(285) 
  reac_source_local(24,286) = + reac_rate_local(286) 
  reac_source_local(25,286) = - reac_rate_local(286) 
  reac_source_local(35,286) = - reac_rate_local(286) 
  reac_source_local(44,286) = + reac_rate_local(286) 
  reac_source_local(05,287) = + reac_rate_local(287) 
  reac_source_local(25,287) = - reac_rate_local(287) 
  reac_source_local(32,287) = + reac_rate_local(287) 
  reac_source_local(35,287) = - reac_rate_local(287) 
  reac_source_local(06,288) = + reac_rate_local(288) 
  reac_source_local(19,288) = - reac_rate_local(288) 
  reac_source_local(24,288) = + reac_rate_local(288) 
  reac_source_local(25,288) = - reac_rate_local(288) 
  reac_source_local(05,289) = + reac_rate_local(289) 
  reac_source_local(14,289) = + reac_rate_local(289) 
  reac_source_local(19,289) = - reac_rate_local(289) 
  reac_source_local(25,289) = - reac_rate_local(289) 
  reac_source_local(02,290) = - reac_rate_local(290) 
  reac_source_local(13,290) = + reac_rate_local(290) 
  reac_source_local(25,290) = - reac_rate_local(290) 
  reac_source_local(41,290) = + reac_rate_local(290) 
  reac_source_local(02,291) = + reac_rate_local(291) 
  reac_source_local(13,291) = - reac_rate_local(291) 
  reac_source_local(25,291) = - reac_rate_local(291) 
  reac_source_local(48,291) = + reac_rate_local(291) 
  reac_source_local(05,292) = - reac_rate_local(292) 
  reac_source_local(24,292) = + reac_rate_local(292) 
  reac_source_local(25,292) = + reac_rate_local(292) 
  reac_source_local(48,292) = - reac_rate_local(292) 
  reac_source_local(02,293) = + reac_rate_local(293) 
  reac_source_local(05,293) = - reac_rate_local(293) 
  reac_source_local(44,293) = + reac_rate_local(293) 
  reac_source_local(48,293) = - reac_rate_local(293) 
  reac_source_local(02,294) = + reac_rate_local(294) 
  reac_source_local(06,294) = + reac_rate_local(294) 
  reac_source_local(10,294) = - reac_rate_local(294) 
  reac_source_local(48,294) = - reac_rate_local(294) 
  reac_source_local(02,295) = + reac_rate_local(295) 
  reac_source_local(14,295) = + reac_rate_local(295) 
  reac_source_local(46,295) = - reac_rate_local(295) 
  reac_source_local(48,295) = - reac_rate_local(295) 
  reac_source_local(05,296) = + reac_rate_local(296) 
  reac_source_local(20,296) = - reac_rate_local(296) 
  reac_source_local(44,296) = + reac_rate_local(296) 
  reac_source_local(48,296) = - reac_rate_local(296) 
  reac_source_local(05,297) = + reac_rate_local(297) 
  reac_source_local(06,297) = + reac_rate_local(297) 
  reac_source_local(35,297) = - reac_rate_local(297) 
  reac_source_local(48,297) = - reac_rate_local(297) 
  reac_source_local(06,298) = + reac_rate_local(298) 
  reac_source_local(24,298) = + reac_rate_local(298) 
  reac_source_local(48,298) = - reac_rate_local(298) 
  reac_source_local(50,298) = - reac_rate_local(298) 
  reac_source_local(05,299) = - reac_rate_local(299) 
  reac_source_local(24,299) = + reac_rate_local(299) 
  reac_source_local(27,299) = - reac_rate_local(299) 
  reac_source_local(48,299) = + reac_rate_local(299) 
  reac_source_local(05,300) = - reac_rate_local(300) 
  reac_source_local(13,300) = + reac_rate_local(300) 
  reac_source_local(27,300) = - reac_rate_local(300) 
  reac_source_local(44,300) = + reac_rate_local(300) 
  reac_source_local(02,301) = + reac_rate_local(301) 
  reac_source_local(05,301) = - reac_rate_local(301) 
  reac_source_local(27,301) = - reac_rate_local(301) 
  reac_source_local(32,301) = + reac_rate_local(301) 
  reac_source_local(02,302) = + reac_rate_local(302) 
  reac_source_local(05,302) = - reac_rate_local(302) 
  reac_source_local(06,302) = + reac_rate_local(302) 
  reac_source_local(13,302) = + reac_rate_local(302) 
  reac_source_local(27,302) = - reac_rate_local(302) 
  reac_source_local(02,303) = + reac_rate_local(303) * 2.d0
  reac_source_local(05,303) = - reac_rate_local(303) 
  reac_source_local(14,303) = + reac_rate_local(303) 
  reac_source_local(27,303) = - reac_rate_local(303) 
  reac_source_local(02,304) = - reac_rate_local(304) 
  reac_source_local(13,304) = + reac_rate_local(304) 
  reac_source_local(27,304) = - reac_rate_local(304) 
  reac_source_local(48,304) = + reac_rate_local(304) 
  reac_source_local(05,305) = - reac_rate_local(305) 
  reac_source_local(13,305) = + reac_rate_local(305) 
  reac_source_local(28,305) = - reac_rate_local(305) 
  reac_source_local(32,305) = + reac_rate_local(305) 
  reac_source_local(02,306) = + reac_rate_local(306) 
  reac_source_local(05,306) = - reac_rate_local(306) 
  reac_source_local(06,306) = + reac_rate_local(306) 
  reac_source_local(28,306) = - reac_rate_local(306) 
  reac_source_local(02,307) = + reac_rate_local(307) 
  reac_source_local(05,307) = - reac_rate_local(307) 
  reac_source_local(13,307) = + reac_rate_local(307) 
  reac_source_local(14,307) = + reac_rate_local(307) 
  reac_source_local(28,307) = - reac_rate_local(307) 
  reac_source_local(02,308) = - reac_rate_local(308) 
  reac_source_local(13,308) = + reac_rate_local(308) 
  reac_source_local(27,308) = + reac_rate_local(308) 
  reac_source_local(28,308) = - reac_rate_local(308) 
  reac_source_local(15,309) = - reac_rate_local(309) 
  reac_source_local(20,309) = + reac_rate_local(309) 
  reac_source_local(32,309) = + reac_rate_local(309) 
  reac_source_local(35,309) = - reac_rate_local(309) 
  reac_source_local(15,310) = - reac_rate_local(310) 
  reac_source_local(19,310) = - reac_rate_local(310) 
  reac_source_local(44,310) = + reac_rate_local(310) 
  reac_source_local(50,310) = + reac_rate_local(310) 
  reac_source_local(02,311) = + reac_rate_local(311) 
  reac_source_local(13,311) = - reac_rate_local(311) 
  reac_source_local(15,311) = - reac_rate_local(311) 
  reac_source_local(44,311) = + reac_rate_local(311) 
  reac_source_local(02,312) = + reac_rate_local(312) 
  reac_source_local(13,312) = - reac_rate_local(312) 
  reac_source_local(32,312) = + reac_rate_local(312) 
  reac_source_local(44,312) = - reac_rate_local(312) 
  reac_source_local(19,313) = + reac_rate_local(313) 
  reac_source_local(32,313) = - reac_rate_local(313) 
  reac_source_local(44,313) = + reac_rate_local(313) 
  reac_source_local(50,313) = - reac_rate_local(313) 
  reac_source_local(06,314) = + reac_rate_local(314) 
  reac_source_local(32,314) = - reac_rate_local(314) 
  reac_source_local(35,314) = + reac_rate_local(314) 
  reac_source_local(50,314) = - reac_rate_local(314) 
  reac_source_local(02,315) = + reac_rate_local(315) 
  reac_source_local(06,315) = + reac_rate_local(315) 
  reac_source_local(13,315) = - reac_rate_local(315) 
  reac_source_local(32,315) = - reac_rate_local(315) 
  reac_source_local(06,316) = - reac_rate_local(316) 
  reac_source_local(20,316) = - reac_rate_local(316) 
  reac_source_local(35,316) = + reac_rate_local(316) 
  reac_source_local(44,316) = + reac_rate_local(316) 
  reac_source_local(06,317) = - reac_rate_local(317) 
  reac_source_local(19,317) = + reac_rate_local(317) 
  reac_source_local(35,317) = - reac_rate_local(317) 
  reac_source_local(44,317) = + reac_rate_local(317) 
  reac_source_local(06,318) = - reac_rate_local(318) 
  reac_source_local(14,318) = + reac_rate_local(318) 
  reac_source_local(18,318) = - reac_rate_local(318) 
  reac_source_local(19,318) = + reac_rate_local(318) 
  reac_source_local(02,319) = + reac_rate_local(319) 
  reac_source_local(06,319) = - reac_rate_local(319) 
  reac_source_local(13,319) = - reac_rate_local(319) 
  reac_source_local(14,319) = + reac_rate_local(319) 
  reac_source_local(05,320) = - reac_rate_local(320) 
  reac_source_local(06,320) = + reac_rate_local(320) 
  reac_source_local(14,320) = - reac_rate_local(320) 
  reac_source_local(24,320) = + reac_rate_local(320) 
  reac_source_local(14,321) = - reac_rate_local(321) 
  reac_source_local(20,321) = - reac_rate_local(321) 
  reac_source_local(44,321) = + reac_rate_local(321) 
  reac_source_local(50,321) = + reac_rate_local(321) 
  reac_source_local(14,322) = - reac_rate_local(322) 
  reac_source_local(20,322) = - reac_rate_local(322) 
  reac_source_local(32,322) = + reac_rate_local(322) 
  reac_source_local(35,322) = + reac_rate_local(322) 
  reac_source_local(14,323) = - reac_rate_local(323) 
  reac_source_local(19,323) = + reac_rate_local(323) 
  reac_source_local(32,323) = + reac_rate_local(323) 
  reac_source_local(35,323) = - reac_rate_local(323) 
  reac_source_local(06,324) = + reac_rate_local(324) 
  reac_source_local(14,324) = - reac_rate_local(324) 
  reac_source_local(19,324) = + reac_rate_local(324) 
  reac_source_local(50,324) = - reac_rate_local(324) 
  reac_source_local(02,325) = - reac_rate_local(325) 
  reac_source_local(06,325) = + reac_rate_local(325) 
  reac_source_local(13,325) = + reac_rate_local(325) 
  reac_source_local(14,325) = - reac_rate_local(325) 
  reac_source_local(02,326) = + reac_rate_local(326) 
  reac_source_local(05,326) = - reac_rate_local(326) 
  reac_source_local(26,326) = - reac_rate_local(326) 
  reac_source_local(41,326) = + reac_rate_local(326) 
  reac_source_local(02,327) = + reac_rate_local(327) 
  reac_source_local(24,327) = - reac_rate_local(327) 
  reac_source_local(25,327) = + reac_rate_local(327) 
  reac_source_local(26,327) = - reac_rate_local(327) 
  reac_source_local(02,328) = + reac_rate_local(328) 
  reac_source_local(10,328) = - reac_rate_local(328) 
  reac_source_local(26,328) = - reac_rate_local(328) 
  reac_source_local(48,328) = + reac_rate_local(328) 
  reac_source_local(02,329) = + reac_rate_local(329) 
  reac_source_local(26,329) = - reac_rate_local(329) 
  reac_source_local(27,329) = + reac_rate_local(329) 
  reac_source_local(46,329) = - reac_rate_local(329) 
  reac_source_local(02,330) = + reac_rate_local(330) * 2.d0
  reac_source_local(20,330) = - reac_rate_local(330) 
  reac_source_local(26,330) = - reac_rate_local(330) 
  reac_source_local(44,330) = + reac_rate_local(330) 
  reac_source_local(02,331) = + reac_rate_local(331) 
  reac_source_local(15,331) = + reac_rate_local(331) 
  reac_source_local(26,331) = - reac_rate_local(331) 
  reac_source_local(38,331) = - reac_rate_local(331) 
  reac_source_local(02,332) = + reac_rate_local(332) 
  reac_source_local(26,332) = - reac_rate_local(332) 
  reac_source_local(35,332) = - reac_rate_local(332) 
  reac_source_local(44,332) = + reac_rate_local(332) 
  reac_source_local(02,333) = + reac_rate_local(333) * 2.d0
  reac_source_local(06,333) = + reac_rate_local(333) 
  reac_source_local(26,333) = - reac_rate_local(333) 
  reac_source_local(35,333) = - reac_rate_local(333) 
  reac_source_local(02,334) = + reac_rate_local(334) 
  reac_source_local(26,334) = - reac_rate_local(334) 
  reac_source_local(32,334) = + reac_rate_local(334) 
  reac_source_local(50,334) = - reac_rate_local(334) 
  reac_source_local(02,335) = + reac_rate_local(335) 
  reac_source_local(14,335) = + reac_rate_local(335) 
  reac_source_local(18,335) = - reac_rate_local(335) 
  reac_source_local(26,335) = - reac_rate_local(335) 
  reac_source_local(02,336) = + reac_rate_local(336) 
  reac_source_local(06,336) = + reac_rate_local(336) 
  reac_source_local(19,336) = - reac_rate_local(336) 
  reac_source_local(26,336) = - reac_rate_local(336) 
  reac_source_local(05,337) = - reac_rate_local(337) 
  reac_source_local(13,337) = + reac_rate_local(337) 
  reac_source_local(41,337) = + reac_rate_local(337) 
  reac_source_local(45,337) = - reac_rate_local(337) 
  reac_source_local(02,338) = + reac_rate_local(338) 
  reac_source_local(05,338) = - reac_rate_local(338) 
  reac_source_local(25,338) = + reac_rate_local(338) 
  reac_source_local(45,338) = - reac_rate_local(338) 
  reac_source_local(02,339) = + reac_rate_local(339) 
  reac_source_local(05,339) = - reac_rate_local(339) 
  reac_source_local(13,339) = + reac_rate_local(339) 
  reac_source_local(45,339) = - reac_rate_local(339) 
  reac_source_local(48,339) = + reac_rate_local(339) 
  reac_source_local(10,340) = - reac_rate_local(340) 
  reac_source_local(13,340) = + reac_rate_local(340) 
  reac_source_local(45,340) = - reac_rate_local(340) 
  reac_source_local(48,340) = + reac_rate_local(340) 
  reac_source_local(02,341) = + reac_rate_local(341) 
  reac_source_local(10,341) = - reac_rate_local(341) 
  reac_source_local(27,341) = + reac_rate_local(341) 
  reac_source_local(45,341) = - reac_rate_local(341) 
  reac_source_local(13,342) = + reac_rate_local(342) 
  reac_source_local(27,342) = + reac_rate_local(342) 
  reac_source_local(45,342) = - reac_rate_local(342) 
  reac_source_local(46,342) = - reac_rate_local(342) 
  reac_source_local(02,343) = + reac_rate_local(343) 
  reac_source_local(28,343) = + reac_rate_local(343) 
  reac_source_local(45,343) = - reac_rate_local(343) 
  reac_source_local(46,343) = - reac_rate_local(343) 
  reac_source_local(02,344) = + reac_rate_local(344) 
  reac_source_local(15,344) = + reac_rate_local(344) 
  reac_source_local(20,344) = - reac_rate_local(344) 
  reac_source_local(45,344) = - reac_rate_local(344) 
  reac_source_local(02,345) = + reac_rate_local(345) 
  reac_source_local(13,345) = + reac_rate_local(345) 
  reac_source_local(20,345) = - reac_rate_local(345) 
  reac_source_local(44,345) = + reac_rate_local(345) 
  reac_source_local(45,345) = - reac_rate_local(345) 
  reac_source_local(02,346) = + reac_rate_local(346) * 2.d0
  reac_source_local(20,346) = - reac_rate_local(346) 
  reac_source_local(32,346) = + reac_rate_local(346) 
  reac_source_local(45,346) = - reac_rate_local(346) 
  reac_source_local(02,347) = + reac_rate_local(347) * 2.d0
  reac_source_local(06,347) = + reac_rate_local(347) 
  reac_source_local(13,347) = + reac_rate_local(347) 
  reac_source_local(20,347) = - reac_rate_local(347) 
  reac_source_local(45,347) = - reac_rate_local(347) 
  reac_source_local(02,348) = + reac_rate_local(348) * 3.d0
  reac_source_local(14,348) = + reac_rate_local(348) 
  reac_source_local(20,348) = - reac_rate_local(348) 
  reac_source_local(45,348) = - reac_rate_local(348) 
  reac_source_local(02,349) = + reac_rate_local(349) 
  reac_source_local(32,349) = + reac_rate_local(349) 
  reac_source_local(35,349) = - reac_rate_local(349) 
  reac_source_local(45,349) = - reac_rate_local(349) 
  reac_source_local(02,350) = + reac_rate_local(350) 
  reac_source_local(06,350) = + reac_rate_local(350) 
  reac_source_local(13,350) = + reac_rate_local(350) 
  reac_source_local(35,350) = - reac_rate_local(350) 
  reac_source_local(45,350) = - reac_rate_local(350) 
  reac_source_local(02,351) = + reac_rate_local(351) * 2.d0
  reac_source_local(14,351) = + reac_rate_local(351) 
  reac_source_local(35,351) = - reac_rate_local(351) 
  reac_source_local(45,351) = - reac_rate_local(351) 
  reac_source_local(06,352) = + reac_rate_local(352) 
  reac_source_local(13,352) = + reac_rate_local(352) 
  reac_source_local(19,352) = - reac_rate_local(352) 
  reac_source_local(45,352) = - reac_rate_local(352) 
  reac_source_local(02,353) = + reac_rate_local(353) 
  reac_source_local(14,353) = + reac_rate_local(353) 
  reac_source_local(19,353) = - reac_rate_local(353) 
  reac_source_local(45,353) = - reac_rate_local(353) 
  reac_source_local(13,354) = - reac_rate_local(354) 
  reac_source_local(26,354) = + reac_rate_local(354) 
  reac_source_local(45,354) = - reac_rate_local(354) 
  reac_source_local(02,355) = + reac_rate_local(355) 
  reac_source_local(07,355) = + reac_rate_local(355) 
  reac_source_local(13,355) = - reac_rate_local(355) 
  reac_source_local(45,355) = - reac_rate_local(355) 
  reac_source_local(05,356) = - reac_rate_local(356) 
  reac_source_local(07,356) = - reac_rate_local(356) 
  reac_source_local(13,356) = + reac_rate_local(356) 
  reac_source_local(25,356) = + reac_rate_local(356) 
  reac_source_local(02,357) = + reac_rate_local(357) 
  reac_source_local(05,357) = - reac_rate_local(357) 
  reac_source_local(07,357) = - reac_rate_local(357) 
  reac_source_local(48,357) = + reac_rate_local(357) 
  reac_source_local(07,358) = - reac_rate_local(358) 
  reac_source_local(13,358) = + reac_rate_local(358) 
  reac_source_local(24,358) = - reac_rate_local(358) 
  reac_source_local(48,358) = + reac_rate_local(358) 
  reac_source_local(07,359) = - reac_rate_local(359) 
  reac_source_local(10,359) = - reac_rate_local(359) 
  reac_source_local(13,359) = + reac_rate_local(359) 
  reac_source_local(27,359) = + reac_rate_local(359) 
  reac_source_local(02,360) = + reac_rate_local(360) 
  reac_source_local(07,360) = - reac_rate_local(360) 
  reac_source_local(10,360) = - reac_rate_local(360) 
  reac_source_local(28,360) = + reac_rate_local(360) 
  reac_source_local(07,361) = - reac_rate_local(361) 
  reac_source_local(13,361) = + reac_rate_local(361) 
  reac_source_local(28,361) = + reac_rate_local(361) 
  reac_source_local(46,361) = - reac_rate_local(361) 
  reac_source_local(02,362) = + reac_rate_local(362) 
  reac_source_local(07,362) = - reac_rate_local(362) 
  reac_source_local(20,362) = - reac_rate_local(362) 
  reac_source_local(44,362) = + reac_rate_local(362) 
  reac_source_local(02,363) = + reac_rate_local(363) 
  reac_source_local(07,363) = - reac_rate_local(363) 
  reac_source_local(13,363) = + reac_rate_local(363) 
  reac_source_local(20,363) = - reac_rate_local(363) 
  reac_source_local(32,363) = + reac_rate_local(363) 
  reac_source_local(02,364) = + reac_rate_local(364) * 2.d0
  reac_source_local(06,364) = + reac_rate_local(364) 
  reac_source_local(07,364) = - reac_rate_local(364) 
  reac_source_local(20,364) = - reac_rate_local(364) 
  reac_source_local(02,365) = + reac_rate_local(365) 
  reac_source_local(07,365) = - reac_rate_local(365) 
  reac_source_local(32,365) = + reac_rate_local(365) 
  reac_source_local(38,365) = - reac_rate_local(365) 
  reac_source_local(02,366) = + reac_rate_local(366) 
  reac_source_local(06,366) = + reac_rate_local(366) 
  reac_source_local(07,366) = - reac_rate_local(366) 
  reac_source_local(13,366) = + reac_rate_local(366) 
  reac_source_local(38,366) = - reac_rate_local(366) 
  reac_source_local(07,367) = - reac_rate_local(367) 
  reac_source_local(13,367) = + reac_rate_local(367) 
  reac_source_local(32,367) = + reac_rate_local(367) 
  reac_source_local(35,367) = - reac_rate_local(367) 
  reac_source_local(02,368) = + reac_rate_local(368) 
  reac_source_local(06,368) = + reac_rate_local(368) 
  reac_source_local(07,368) = - reac_rate_local(368) 
  reac_source_local(35,368) = - reac_rate_local(368) 
  reac_source_local(02,369) = + reac_rate_local(369) 
  reac_source_local(07,369) = - reac_rate_local(369) 
  reac_source_local(13,369) = + reac_rate_local(369) 
  reac_source_local(14,369) = + reac_rate_local(369) 
  reac_source_local(35,369) = - reac_rate_local(369) 
  reac_source_local(06,370) = + reac_rate_local(370) 
  reac_source_local(07,370) = - reac_rate_local(370) 
  reac_source_local(13,370) = + reac_rate_local(370) 
  reac_source_local(50,370) = - reac_rate_local(370) 
  reac_source_local(02,371) = + reac_rate_local(371) 
  reac_source_local(07,371) = - reac_rate_local(371) 
  reac_source_local(14,371) = + reac_rate_local(371) 
  reac_source_local(50,371) = - reac_rate_local(371) 
  reac_source_local(07,372) = - reac_rate_local(372) 
  reac_source_local(13,372) = + reac_rate_local(372) 
  reac_source_local(14,372) = + reac_rate_local(372) 
  reac_source_local(19,372) = - reac_rate_local(372) 
  reac_source_local(05,373) = - reac_rate_local(373) 
  reac_source_local(10,373) = - reac_rate_local(373) 
  reac_source_local(24,373) = + reac_rate_local(373) * 2.d0
  reac_source_local(05,374) = - reac_rate_local(374) 
  reac_source_local(13,374) = + reac_rate_local(374) 
  reac_source_local(35,374) = + reac_rate_local(374) 
  reac_source_local(46,374) = - reac_rate_local(374) 
  reac_source_local(05,375) = - reac_rate_local(375) 
  reac_source_local(20,375) = + reac_rate_local(375) 
  reac_source_local(24,375) = + reac_rate_local(375) 
  reac_source_local(38,375) = - reac_rate_local(375) 
  reac_source_local(05,376) = - reac_rate_local(376) 
  reac_source_local(24,376) = + reac_rate_local(376) 
  reac_source_local(35,376) = + reac_rate_local(376) 
  reac_source_local(50,376) = - reac_rate_local(376) 
  reac_source_local(05,377) = - reac_rate_local(377) 
  reac_source_local(18,377) = - reac_rate_local(377) 
  reac_source_local(19,377) = + reac_rate_local(377) 
  reac_source_local(24,377) = + reac_rate_local(377) 
  reac_source_local(05,378) = - reac_rate_local(378) 
  reac_source_local(12,378) = + reac_rate_local(378) 
  reac_source_local(23,378) = - reac_rate_local(378) 
  reac_source_local(24,378) = + reac_rate_local(378) 
  reac_source_local(03,379) = - reac_rate_local(379) 
  reac_source_local(05,379) = - reac_rate_local(379) 
  reac_source_local(09,379) = + reac_rate_local(379) 
  reac_source_local(24,379) = + reac_rate_local(379) 
  reac_source_local(02,380) = + reac_rate_local(380) 
  reac_source_local(05,380) = - reac_rate_local(380) 
  reac_source_local(13,380) = - reac_rate_local(380) 
  reac_source_local(24,380) = + reac_rate_local(380) 
  reac_source_local(05,381) = - reac_rate_local(381) 
  reac_source_local(13,381) = + reac_rate_local(381) 
  reac_source_local(20,381) = + reac_rate_local(381) 
  reac_source_local(24,381) = - reac_rate_local(381) 
  reac_source_local(05,382) = - reac_rate_local(382) 
  reac_source_local(10,382) = - reac_rate_local(382) 
  reac_source_local(20,382) = + reac_rate_local(382) 
  reac_source_local(05,383) = - reac_rate_local(383) 
  reac_source_local(13,383) = + reac_rate_local(383) 
  reac_source_local(24,383) = + reac_rate_local(383) 
  reac_source_local(10,384) = + reac_rate_local(384) 
  reac_source_local(13,384) = + reac_rate_local(384) 
  reac_source_local(24,384) = - reac_rate_local(384) 
  reac_source_local(02,385) = + reac_rate_local(385) 
  reac_source_local(24,385) = - reac_rate_local(385) 
  reac_source_local(46,385) = + reac_rate_local(385) 
  reac_source_local(10,386) = + reac_rate_local(386) 
  reac_source_local(20,386) = + reac_rate_local(386) 
  reac_source_local(24,386) = - reac_rate_local(386) 
  reac_source_local(38,386) = - reac_rate_local(386) 
  reac_source_local(10,387) = - reac_rate_local(387) * 2.d0
  reac_source_local(13,387) = + reac_rate_local(387) * 2.d0
  reac_source_local(19,387) = + reac_rate_local(387) 
  reac_source_local(10,388) = - reac_rate_local(388) 
  reac_source_local(24,388) = + reac_rate_local(388) 
  reac_source_local(35,388) = + reac_rate_local(388) 
  reac_source_local(38,388) = - reac_rate_local(388) 
  reac_source_local(10,389) = - reac_rate_local(389) 
  reac_source_local(19,389) = + reac_rate_local(389) 
  reac_source_local(24,389) = + reac_rate_local(389) 
  reac_source_local(50,389) = - reac_rate_local(389) 
  reac_source_local(10,390) = - reac_rate_local(390) 
  reac_source_local(18,390) = - reac_rate_local(390) 
  reac_source_local(19,390) = + reac_rate_local(390) 
  reac_source_local(46,390) = + reac_rate_local(390) 
  reac_source_local(10,391) = - reac_rate_local(391) 
  reac_source_local(12,391) = - reac_rate_local(391) 
  reac_source_local(23,391) = + reac_rate_local(391) 
  reac_source_local(24,391) = + reac_rate_local(391) 
  reac_source_local(10,392) = - reac_rate_local(392) 
  reac_source_local(23,392) = - reac_rate_local(392) 
  reac_source_local(35,392) = + reac_rate_local(392) 
  reac_source_local(38,392) = + reac_rate_local(392) 
  reac_source_local(09,393) = + reac_rate_local(393) 
  reac_source_local(10,393) = - reac_rate_local(393) 
  reac_source_local(23,393) = - reac_rate_local(393) 
  reac_source_local(24,393) = + reac_rate_local(393) 
  reac_source_local(03,394) = + reac_rate_local(394) 
  reac_source_local(09,394) = - reac_rate_local(394) 
  reac_source_local(10,394) = - reac_rate_local(394) 
  reac_source_local(24,394) = + reac_rate_local(394) 
  reac_source_local(02,395) = - reac_rate_local(395) 
  reac_source_local(10,395) = - reac_rate_local(395) 
  reac_source_local(13,395) = + reac_rate_local(395) 
  reac_source_local(24,395) = + reac_rate_local(395) 
  reac_source_local(02,396) = + reac_rate_local(396) 
  reac_source_local(10,396) = - reac_rate_local(396) 
  reac_source_local(13,396) = - reac_rate_local(396) 
  reac_source_local(46,396) = + reac_rate_local(396) 
  reac_source_local(10,397) = - reac_rate_local(397) 
  reac_source_local(13,397) = + reac_rate_local(397) 
  reac_source_local(46,397) = + reac_rate_local(397) 
  reac_source_local(02,398) = + reac_rate_local(398) 
  reac_source_local(10,398) = - reac_rate_local(398) * 2.d0
  reac_source_local(19,398) = + reac_rate_local(398) 
  reac_source_local(10,399) = - reac_rate_local(399) 
  reac_source_local(13,399) = - reac_rate_local(399) 
  reac_source_local(24,399) = + reac_rate_local(399) 
  reac_source_local(09,400) = + reac_rate_local(400) 
  reac_source_local(13,400) = + reac_rate_local(400) 
  reac_source_local(20,400) = - reac_rate_local(400) 
  reac_source_local(46,400) = - reac_rate_local(400) 
  reac_source_local(20,401) = - reac_rate_local(401) 
  reac_source_local(23,401) = + reac_rate_local(401) 
  reac_source_local(46,401) = - reac_rate_local(401) 
  reac_source_local(02,402) = - reac_rate_local(402) 
  reac_source_local(10,402) = + reac_rate_local(402) 
  reac_source_local(13,402) = + reac_rate_local(402) 
  reac_source_local(46,402) = - reac_rate_local(402) 
  reac_source_local(13,403) = + reac_rate_local(403) 
  reac_source_local(24,403) = - reac_rate_local(403) 
  reac_source_local(46,403) = - reac_rate_local(403) 
  reac_source_local(50,403) = + reac_rate_local(403) 
  reac_source_local(10,404) = - reac_rate_local(404) 
  reac_source_local(13,404) = + reac_rate_local(404) 
  reac_source_local(19,404) = + reac_rate_local(404) 
  reac_source_local(46,404) = - reac_rate_local(404) 
  reac_source_local(02,405) = - reac_rate_local(405) 
  reac_source_local(24,405) = + reac_rate_local(405) 
  reac_source_local(46,405) = - reac_rate_local(405) 
  reac_source_local(10,406) = + reac_rate_local(406) 
  reac_source_local(19,406) = + reac_rate_local(406) 
  reac_source_local(46,406) = - reac_rate_local(406) 
  reac_source_local(50,406) = - reac_rate_local(406) 
  reac_source_local(20,407) = - reac_rate_local(407) 
  reac_source_local(35,407) = + reac_rate_local(407) 
  reac_source_local(38,407) = + reac_rate_local(407) 
  reac_source_local(50,407) = - reac_rate_local(407) 
  reac_source_local(12,408) = + reac_rate_local(408) 
  reac_source_local(20,408) = - reac_rate_local(408) 
  reac_source_local(23,408) = - reac_rate_local(408) 
  reac_source_local(38,408) = + reac_rate_local(408) 
  reac_source_local(03,409) = - reac_rate_local(409) 
  reac_source_local(09,409) = + reac_rate_local(409) 
  reac_source_local(20,409) = - reac_rate_local(409) 
  reac_source_local(38,409) = + reac_rate_local(409) 
  reac_source_local(02,410) = + reac_rate_local(410) 
  reac_source_local(13,410) = - reac_rate_local(410) 
  reac_source_local(20,410) = - reac_rate_local(410) 
  reac_source_local(38,410) = + reac_rate_local(410) 
  reac_source_local(05,411) = + reac_rate_local(411) 
  reac_source_local(13,411) = - reac_rate_local(411) 
  reac_source_local(20,411) = - reac_rate_local(411) 
  reac_source_local(24,411) = + reac_rate_local(411) 
  reac_source_local(20,412) = - reac_rate_local(412) 
  reac_source_local(24,412) = + reac_rate_local(412) * 2.d0
  reac_source_local(20,413) = - reac_rate_local(413) 
  reac_source_local(24,413) = + reac_rate_local(413) 
  reac_source_local(35,413) = + reac_rate_local(413) 
  reac_source_local(46,413) = - reac_rate_local(413) 
  reac_source_local(10,414) = - reac_rate_local(414) 
  reac_source_local(20,414) = - reac_rate_local(414) 
  reac_source_local(24,414) = + reac_rate_local(414) 
  reac_source_local(38,414) = + reac_rate_local(414) 
  reac_source_local(20,415) = + reac_rate_local(415) 
  reac_source_local(35,415) = + reac_rate_local(415) 
  reac_source_local(38,415) = - reac_rate_local(415) * 2.d0
  reac_source_local(20,416) = + reac_rate_local(416) 
  reac_source_local(35,416) = - reac_rate_local(416) 
  reac_source_local(38,416) = - reac_rate_local(416) 
  reac_source_local(50,416) = + reac_rate_local(416) 
  reac_source_local(18,417) = + reac_rate_local(417) 
  reac_source_local(19,417) = - reac_rate_local(417) 
  reac_source_local(20,417) = + reac_rate_local(417) 
  reac_source_local(38,417) = - reac_rate_local(417) 
  reac_source_local(18,418) = - reac_rate_local(418) 
  reac_source_local(19,418) = + reac_rate_local(418) 
  reac_source_local(35,418) = + reac_rate_local(418) 
  reac_source_local(38,418) = - reac_rate_local(418) 
  reac_source_local(12,419) = - reac_rate_local(419) 
  reac_source_local(20,419) = + reac_rate_local(419) 
  reac_source_local(23,419) = + reac_rate_local(419) 
  reac_source_local(38,419) = - reac_rate_local(419) 
  reac_source_local(12,420) = + reac_rate_local(420) 
  reac_source_local(23,420) = - reac_rate_local(420) 
  reac_source_local(35,420) = + reac_rate_local(420) 
  reac_source_local(38,420) = - reac_rate_local(420) 
  reac_source_local(09,421) = + reac_rate_local(421) 
  reac_source_local(20,421) = + reac_rate_local(421) 
  reac_source_local(23,421) = - reac_rate_local(421) 
  reac_source_local(38,421) = - reac_rate_local(421) 
  reac_source_local(03,422) = + reac_rate_local(422) 
  reac_source_local(09,422) = - reac_rate_local(422) 
  reac_source_local(20,422) = + reac_rate_local(422) 
  reac_source_local(38,422) = - reac_rate_local(422) 
  reac_source_local(02,423) = - reac_rate_local(423) 
  reac_source_local(13,423) = + reac_rate_local(423) 
  reac_source_local(20,423) = + reac_rate_local(423) 
  reac_source_local(38,423) = - reac_rate_local(423) 
  reac_source_local(13,424) = - reac_rate_local(424) 
  reac_source_local(24,424) = + reac_rate_local(424) * 2.d0
  reac_source_local(38,424) = - reac_rate_local(424) 
  reac_source_local(02,425) = + reac_rate_local(425) 
  reac_source_local(13,425) = - reac_rate_local(425) 
  reac_source_local(35,425) = + reac_rate_local(425) 
  reac_source_local(38,425) = - reac_rate_local(425) 
  reac_source_local(13,426) = - reac_rate_local(426) 
  reac_source_local(20,426) = + reac_rate_local(426) 
  reac_source_local(38,426) = - reac_rate_local(426) 
  reac_source_local(13,427) = + reac_rate_local(427) 
  reac_source_local(35,427) = + reac_rate_local(427) 
  reac_source_local(38,427) = - reac_rate_local(427) 
  reac_source_local(35,428) = + reac_rate_local(428) * 2.d0
  reac_source_local(38,428) = - reac_rate_local(428) 
  reac_source_local(50,428) = - reac_rate_local(428) 
  reac_source_local(02,429) = + reac_rate_local(429) 
  reac_source_local(13,429) = - reac_rate_local(429) 
  reac_source_local(35,429) = - reac_rate_local(429) 
  reac_source_local(50,429) = + reac_rate_local(429) 
  reac_source_local(13,430) = - reac_rate_local(430) 
  reac_source_local(35,430) = - reac_rate_local(430) 
  reac_source_local(38,430) = + reac_rate_local(430) 
  reac_source_local(02,431) = - reac_rate_local(431) 
  reac_source_local(13,431) = + reac_rate_local(431) 
  reac_source_local(35,431) = - reac_rate_local(431) 
  reac_source_local(38,431) = + reac_rate_local(431) 
  reac_source_local(13,432) = + reac_rate_local(432) 
  reac_source_local(35,432) = - reac_rate_local(432) 
  reac_source_local(50,432) = + reac_rate_local(432) 
  reac_source_local(03,433) = + reac_rate_local(433) 
  reac_source_local(09,433) = - reac_rate_local(433) 
  reac_source_local(35,433) = - reac_rate_local(433) 
  reac_source_local(38,433) = + reac_rate_local(433) 
  reac_source_local(19,434) = - reac_rate_local(434) 
  reac_source_local(35,434) = - reac_rate_local(434) 
  reac_source_local(50,434) = + reac_rate_local(434) * 2.d0
  reac_source_local(09,435) = - reac_rate_local(435) 
  reac_source_local(23,435) = + reac_rate_local(435) 
  reac_source_local(35,435) = - reac_rate_local(435) 
  reac_source_local(50,435) = + reac_rate_local(435) 
  reac_source_local(35,436) = - reac_rate_local(436) * 2.d0
  reac_source_local(38,436) = + reac_rate_local(436) 
  reac_source_local(50,436) = + reac_rate_local(436) 
  reac_source_local(23,437) = + reac_rate_local(437) 
  reac_source_local(24,437) = - reac_rate_local(437) 
  reac_source_local(35,437) = - reac_rate_local(437) 
  reac_source_local(02,438) = + reac_rate_local(438) 
  reac_source_local(19,438) = + reac_rate_local(438) 
  reac_source_local(35,438) = - reac_rate_local(438) 
  reac_source_local(35,439) = - reac_rate_local(439) 
  reac_source_local(38,439) = - reac_rate_local(439) 
  reac_source_local(39,439) = + reac_rate_local(439) 
  reac_source_local(02,440) = - reac_rate_local(440) 
  reac_source_local(20,440) = + reac_rate_local(440) 
  reac_source_local(35,440) = - reac_rate_local(440) 
  reac_source_local(09,441) = + reac_rate_local(441) 
  reac_source_local(10,441) = - reac_rate_local(441) 
  reac_source_local(35,441) = - reac_rate_local(441) 
  reac_source_local(09,442) = + reac_rate_local(442) 
  reac_source_local(23,442) = + reac_rate_local(442) 
  reac_source_local(35,442) = - reac_rate_local(442) 
  reac_source_local(39,442) = - reac_rate_local(442) 
  reac_source_local(19,443) = + reac_rate_local(443) 
  reac_source_local(35,443) = + reac_rate_local(443) 
  reac_source_local(50,443) = - reac_rate_local(443) * 2.d0
  reac_source_local(12,444) = - reac_rate_local(444) 
  reac_source_local(23,444) = + reac_rate_local(444) 
  reac_source_local(35,444) = + reac_rate_local(444) 
  reac_source_local(50,444) = - reac_rate_local(444) 
  reac_source_local(12,445) = + reac_rate_local(445) 
  reac_source_local(19,445) = + reac_rate_local(445) 
  reac_source_local(23,445) = - reac_rate_local(445) 
  reac_source_local(50,445) = - reac_rate_local(445) 
  reac_source_local(09,446) = + reac_rate_local(446) 
  reac_source_local(23,446) = - reac_rate_local(446) 
  reac_source_local(35,446) = + reac_rate_local(446) 
  reac_source_local(50,446) = - reac_rate_local(446) 
  reac_source_local(03,447) = + reac_rate_local(447) 
  reac_source_local(09,447) = - reac_rate_local(447) 
  reac_source_local(35,447) = + reac_rate_local(447) 
  reac_source_local(50,447) = - reac_rate_local(447) 
  reac_source_local(03,448) = - reac_rate_local(448) 
  reac_source_local(09,448) = + reac_rate_local(448) 
  reac_source_local(19,448) = + reac_rate_local(448) 
  reac_source_local(50,448) = - reac_rate_local(448) 
  reac_source_local(02,449) = - reac_rate_local(449) 
  reac_source_local(13,449) = + reac_rate_local(449) 
  reac_source_local(35,449) = + reac_rate_local(449) 
  reac_source_local(50,449) = - reac_rate_local(449) 
  reac_source_local(02,450) = + reac_rate_local(450) 
  reac_source_local(13,450) = - reac_rate_local(450) 
  reac_source_local(19,450) = + reac_rate_local(450) 
  reac_source_local(50,450) = - reac_rate_local(450) 
  reac_source_local(13,451) = - reac_rate_local(451) 
  reac_source_local(35,451) = + reac_rate_local(451) 
  reac_source_local(50,451) = - reac_rate_local(451) 
  reac_source_local(13,452) = + reac_rate_local(452) 
  reac_source_local(19,452) = + reac_rate_local(452) 
  reac_source_local(50,452) = - reac_rate_local(452) 
  reac_source_local(13,453) = - reac_rate_local(453) 
  reac_source_local(19,453) = - reac_rate_local(453) 
  reac_source_local(50,453) = + reac_rate_local(453) 
  reac_source_local(02,454) = - reac_rate_local(454) 
  reac_source_local(19,454) = - reac_rate_local(454) 
  reac_source_local(35,454) = + reac_rate_local(454) 
  reac_source_local(02,455) = - reac_rate_local(455) 
  reac_source_local(13,455) = + reac_rate_local(455) 
  reac_source_local(19,455) = - reac_rate_local(455) 
  reac_source_local(50,455) = + reac_rate_local(455) 
  reac_source_local(03,456) = + reac_rate_local(456) 
  reac_source_local(19,456) = - reac_rate_local(456) 
  reac_source_local(24,456) = - reac_rate_local(456) 
  reac_source_local(03,457) = + reac_rate_local(457) 
  reac_source_local(09,457) = + reac_rate_local(457) 
  reac_source_local(19,457) = - reac_rate_local(457) 
  reac_source_local(39,457) = - reac_rate_local(457) 
  reac_source_local(03,458) = - reac_rate_local(458) 
  reac_source_local(09,458) = + reac_rate_local(458) 
  reac_source_local(12,458) = - reac_rate_local(458) 
  reac_source_local(23,458) = + reac_rate_local(458) 
  reac_source_local(02,459) = + reac_rate_local(459) 
  reac_source_local(12,459) = - reac_rate_local(459) 
  reac_source_local(13,459) = - reac_rate_local(459) 
  reac_source_local(23,459) = + reac_rate_local(459) 
  reac_source_local(12,460) = - reac_rate_local(460) 
  reac_source_local(24,460) = + reac_rate_local(460) 
  reac_source_local(38,460) = + reac_rate_local(460) 
  reac_source_local(09,461) = + reac_rate_local(461) 
  reac_source_local(12,461) = + reac_rate_local(461) 
  reac_source_local(23,461) = - reac_rate_local(461) * 2.d0
  reac_source_local(03,462) = + reac_rate_local(462) 
  reac_source_local(09,462) = - reac_rate_local(462) 
  reac_source_local(12,462) = + reac_rate_local(462) 
  reac_source_local(23,462) = - reac_rate_local(462) 
  reac_source_local(03,463) = - reac_rate_local(463) 
  reac_source_local(09,463) = + reac_rate_local(463) * 2.d0
  reac_source_local(23,463) = - reac_rate_local(463) 
  reac_source_local(02,464) = - reac_rate_local(464) 
  reac_source_local(12,464) = + reac_rate_local(464) 
  reac_source_local(13,464) = + reac_rate_local(464) 
  reac_source_local(23,464) = - reac_rate_local(464) 
  reac_source_local(02,465) = + reac_rate_local(465) 
  reac_source_local(09,465) = + reac_rate_local(465) 
  reac_source_local(13,465) = - reac_rate_local(465) 
  reac_source_local(23,465) = - reac_rate_local(465) 
  reac_source_local(12,466) = + reac_rate_local(466) 
  reac_source_local(13,466) = - reac_rate_local(466) 
  reac_source_local(23,466) = - reac_rate_local(466) 
  reac_source_local(13,467) = - reac_rate_local(467) 
  reac_source_local(23,467) = - reac_rate_local(467) 
  reac_source_local(24,467) = + reac_rate_local(467) 
  reac_source_local(38,467) = + reac_rate_local(467) 
  reac_source_local(09,468) = + reac_rate_local(468) 
  reac_source_local(13,468) = + reac_rate_local(468) 
  reac_source_local(23,468) = - reac_rate_local(468) 
  reac_source_local(23,469) = - reac_rate_local(469) 
  reac_source_local(24,469) = + reac_rate_local(469) 
  reac_source_local(35,469) = + reac_rate_local(469) 
  reac_source_local(03,470) = + reac_rate_local(470) 
  reac_source_local(09,470) = - reac_rate_local(470) 
  reac_source_local(19,470) = - reac_rate_local(470) 
  reac_source_local(50,470) = + reac_rate_local(470) 
  reac_source_local(03,471) = + reac_rate_local(471) 
  reac_source_local(09,471) = - reac_rate_local(471) * 2.d0
  reac_source_local(23,471) = + reac_rate_local(471) 
  reac_source_local(03,472) = + reac_rate_local(472) 
  reac_source_local(09,472) = - reac_rate_local(472) 
  reac_source_local(13,472) = + reac_rate_local(472) 
  reac_source_local(02,473) = + reac_rate_local(473) 
  reac_source_local(03,473) = + reac_rate_local(473) 
  reac_source_local(09,473) = - reac_rate_local(473) 
  reac_source_local(13,473) = - reac_rate_local(473) 
  reac_source_local(09,474) = - reac_rate_local(474) 
  reac_source_local(13,474) = - reac_rate_local(474) 
  reac_source_local(23,474) = + reac_rate_local(474) 
  reac_source_local(09,475) = - reac_rate_local(475) 
  reac_source_local(24,475) = + reac_rate_local(475) 
  reac_source_local(50,475) = + reac_rate_local(475) 
  reac_source_local(09,476) = - reac_rate_local(476) 
  reac_source_local(24,476) = - reac_rate_local(476) 
  reac_source_local(39,476) = + reac_rate_local(476) 
  reac_source_local(02,477) = - reac_rate_local(477) 
  reac_source_local(03,477) = - reac_rate_local(477) 
  reac_source_local(09,477) = + reac_rate_local(477) 
  reac_source_local(13,477) = + reac_rate_local(477) 
  reac_source_local(03,478) = - reac_rate_local(478) 
  reac_source_local(09,478) = + reac_rate_local(478) 
  reac_source_local(13,478) = - reac_rate_local(478) 
  reac_source_local(03,479) = - reac_rate_local(479) 
  reac_source_local(19,479) = + reac_rate_local(479) 
  reac_source_local(24,479) = + reac_rate_local(479) 
  reac_source_local(35,480) = + reac_rate_local(480) 
  reac_source_local(38,480) = + reac_rate_local(480) 
  reac_source_local(39,480) = - reac_rate_local(480) 
  reac_source_local(10,481) = - reac_rate_local(481) 
  reac_source_local(23,481) = + reac_rate_local(481) 
  reac_source_local(35,481) = + reac_rate_local(481) 
  reac_source_local(39,481) = - reac_rate_local(481) 
  reac_source_local(09,482) = + reac_rate_local(482) 
  reac_source_local(24,482) = + reac_rate_local(482) 
  reac_source_local(39,482) = - reac_rate_local(482) 
  reac_source_local(11,483) = - reac_rate_local(483) 
  reac_source_local(24,483) = + reac_rate_local(483) 
  reac_source_local(39,483) = + reac_rate_local(483) 
  reac_source_local(02,484) = - reac_rate_local(484) 
  reac_source_local(13,484) = + reac_rate_local(484) * 2.d0
  reac_source_local(02,485) = + reac_rate_local(485) 
  reac_source_local(13,485) = - reac_rate_local(485) * 2.d0
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
  if( lgas_heating ) ZDPlasKin_cfg(1) = y(51)
  density(:) = y(1:species_max)
  call ZDPlasKin_reac_rates(t)
  rrt(001) = rrt(001) * density(01) * density(05) 
  rrt(002) = rrt(002) * density(01) * density(05) 
  rrt(003) = rrt(003) * density(01) * density(20) 
  rrt(004) = rrt(004) * density(01) * density(20) 
  rrt(005) = rrt(005) * density(01) * density(35) 
  rrt(006) = rrt(006) * density(01) * density(35) 
  rrt(007) = rrt(007) * density(01) * density(19) 
  rrt(008) = rrt(008) * density(01) * density(19) 
  rrt(009) = rrt(009) * density(01) * density(19) 
  rrt(010) = rrt(010) * density(01) * density(12) 
  rrt(011) = rrt(011) * density(01) * density(12) 
  rrt(012) = rrt(012) * density(01) * density(09) 
  rrt(013) = rrt(013) * density(01) * density(05) 
  rrt(014) = rrt(014) * density(01) * density(05) 
  rrt(015) = rrt(015) * density(01) * density(05) 
  rrt(016) = rrt(016) * density(01) * density(24) 
  rrt(017) = rrt(017) * density(01) * density(24) 
  rrt(018) = rrt(018) * density(01) * density(10) 
  rrt(019) = rrt(019) * density(01) * density(05) 
  rrt(020) = rrt(020) * density(01) * density(05) 
  rrt(021) = rrt(021) * density(01) * density(05) 
  rrt(022) = rrt(022) * density(01) * density(05) 
  rrt(023) = rrt(023) * density(01) * density(24) 
  rrt(024) = rrt(024) * density(01) * density(24) 
  rrt(025) = rrt(025) * density(01) * density(24) 
  rrt(026) = rrt(026) * density(01) * density(10) 
  rrt(027) = rrt(027) * density(01) * density(10) 
  rrt(028) = rrt(028) * density(01) * density(46) 
  rrt(029) = rrt(029) * density(01) * density(20) 
  rrt(030) = rrt(030) * density(01) * density(20) 
  rrt(031) = rrt(031) * density(01) * density(20) 
  rrt(032) = rrt(032) * density(01) * density(20) 
  rrt(033) = rrt(033) * density(01) * density(20) 
  rrt(034) = rrt(034) * density(01) * density(20) 
  rrt(035) = rrt(035) * density(01) * density(38) 
  rrt(036) = rrt(036) * density(01) * density(38) 
  rrt(037) = rrt(037) * density(01) * density(38) 
  rrt(038) = rrt(038) * density(01) * density(38) 
  rrt(039) = rrt(039) * density(01) * density(38) 
  rrt(040) = rrt(040) * density(01) * density(38) 
  rrt(041) = rrt(041) * density(01) * density(35) 
  rrt(042) = rrt(042) * density(01) * density(35) 
  rrt(043) = rrt(043) * density(01) * density(35) 
  rrt(044) = rrt(044) * density(01) * density(35) 
  rrt(045) = rrt(045) * density(01) * density(35) 
  rrt(046) = rrt(046) * density(01) * density(50) 
  rrt(047) = rrt(047) * density(01) * density(50) 
  rrt(048) = rrt(048) * density(01) * density(19) 
  rrt(049) = rrt(049) * density(01) * density(12) 
  rrt(050) = rrt(050) * density(01) * density(12) 
  rrt(051) = rrt(051) * density(01) * density(12) 
  rrt(052) = rrt(052) * density(01) * density(12) 
  rrt(053) = rrt(053) * density(01) * density(12) 
  rrt(054) = rrt(054) * density(01) * density(12) 
  rrt(055) = rrt(055) * density(01) * density(23) 
  rrt(056) = rrt(056) * density(01) * density(23) 
  rrt(057) = rrt(057) * density(01) * density(23) 
  rrt(058) = rrt(058) * density(01) * density(23) 
  rrt(059) = rrt(059) * density(01) * density(23) 
  rrt(060) = rrt(060) * density(01) * density(09) 
  rrt(061) = rrt(061) * density(01) * density(09) 
  rrt(062) = rrt(062) * density(01) * density(09) 
  rrt(063) = rrt(063) * density(01) * density(09) 
  rrt(064) = rrt(064) * density(01) * density(09) 
  rrt(065) = rrt(065) * density(01) * density(03) 
  rrt(066) = rrt(066) * density(01) * density(03) 
  rrt(067) = rrt(067) * density(01) * density(17) 
  rrt(068) = rrt(068) * density(01) * density(17) 
  rrt(069) = rrt(069) * density(01) * density(20) 
  rrt(070) = rrt(070) * density(01) * density(20) 
  rrt(071) = rrt(071) * density(01) * density(20) 
  rrt(072) = rrt(072) * density(01) * density(20) 
  rrt(073) = rrt(073) * density(01) * density(20) 
  rrt(074) = rrt(074) * density(01) * density(20) 
  rrt(075) = rrt(075) * density(01) * density(20) 
  rrt(076) = rrt(076) * density(01) * density(38) 
  rrt(077) = rrt(077) * density(01) * density(38) 
  rrt(078) = rrt(078) * density(01) * density(38) 
  rrt(079) = rrt(079) * density(01) * density(38) 
  rrt(080) = rrt(080) * density(01) * density(38) 
  rrt(081) = rrt(081) * density(01) * density(38) 
  rrt(082) = rrt(082) * density(01) * density(38) 
  rrt(083) = rrt(083) * density(01) * density(35) 
  rrt(084) = rrt(084) * density(01) * density(35) 
  rrt(085) = rrt(085) * density(01) * density(35) 
  rrt(086) = rrt(086) * density(01) * density(35) 
  rrt(087) = rrt(087) * density(01) * density(35) 
  rrt(088) = rrt(088) * density(01) * density(50) 
  rrt(089) = rrt(089) * density(01) * density(50) 
  rrt(090) = rrt(090) * density(01) * density(50) 
  rrt(091) = rrt(091) * density(01) * density(50) 
  rrt(092) = rrt(092) * density(01) * density(50) 
  rrt(093) = rrt(093) * density(01) * density(19) 
  rrt(094) = rrt(094) * density(01) * density(19) 
  rrt(095) = rrt(095) * density(01) * density(12) 
  rrt(096) = rrt(096) * density(01) * density(12) 
  rrt(097) = rrt(097) * density(01) * density(12) 
  rrt(098) = rrt(098) * density(01) * density(12) 
  rrt(099) = rrt(099) * density(01) * density(12) 
  rrt(100) = rrt(100) * density(01) * density(12) 
  rrt(101) = rrt(101) * density(01) * density(12) 
  rrt(102) = rrt(102) * density(01) * density(12) 
  rrt(103) = rrt(103) * density(01) * density(12) 
  rrt(104) = rrt(104) * density(01) * density(23) 
  rrt(105) = rrt(105) * density(01) * density(23) 
  rrt(106) = rrt(106) * density(01) * density(23) 
  rrt(107) = rrt(107) * density(01) * density(23) 
  rrt(108) = rrt(108) * density(01) * density(23) 
  rrt(109) = rrt(109) * density(01) * density(23) 
  rrt(110) = rrt(110) * density(01) * density(23) 
  rrt(111) = rrt(111) * density(01) * density(23) 
  rrt(112) = rrt(112) * density(01) * density(23) 
  rrt(113) = rrt(113) * density(01) * density(23) 
  rrt(114) = rrt(114) * density(01) * density(23) 
  rrt(115) = rrt(115) * density(01) * density(09) 
  rrt(116) = rrt(116) * density(01) * density(09) 
  rrt(117) = rrt(117) * density(01) * density(09) 
  rrt(118) = rrt(118) * density(01) * density(09) 
  rrt(119) = rrt(119) * density(01) * density(09) 
  rrt(120) = rrt(120) * density(01) * density(09) 
  rrt(121) = rrt(121) * density(01) * density(09) 
  rrt(122) = rrt(122) * density(01) * density(09) 
  rrt(123) = rrt(123) * density(01) * density(09) 
  rrt(124) = rrt(124) * density(01) * density(09) 
  rrt(125) = rrt(125) * density(01) * density(09) 
  rrt(126) = rrt(126) * density(01) * density(03) 
  rrt(127) = rrt(127) * density(01) * density(03) 
  rrt(128) = rrt(128) * density(01) * density(03) 
  rrt(129) = rrt(129) * density(01) * density(03) 
  rrt(130) = rrt(130) * density(01) * density(03) 
  rrt(131) = rrt(131) * density(01) * density(03) 
  rrt(132) = rrt(132) * density(01) * density(03) 
  rrt(133) = rrt(133) * density(01) * density(03) 
  rrt(134) = rrt(134) * density(01) * density(17) 
  rrt(135) = rrt(135) * density(01) * density(17) 
  rrt(136) = rrt(136) * density(01) * density(17) 
  rrt(137) = rrt(137) * density(01) * density(17) 
  rrt(138) = rrt(138) * density(01) * density(17) 
  rrt(139) = rrt(139) * density(01) * density(43) 
  rrt(140) = rrt(140) * density(01) * density(36) 
  rrt(141) = rrt(141) * density(01) * density(43) 
  rrt(142) = rrt(142) * density(01) * density(36) 
  rrt(143) = rrt(143) * density(01) * density(43) 
  rrt(144) = rrt(144) * density(01) * density(36) 
  rrt(145) = rrt(145) * density(01) * density(43) 
  rrt(146) = rrt(146) * density(01) * density(36) 
  rrt(147) = rrt(147) * density(01) * density(43) 
  rrt(148) = rrt(148) * density(01) * density(36) 
  rrt(149) = rrt(149) * density(01) * density(43) 
  rrt(150) = rrt(150) * density(01) * density(36) 
  rrt(151) = rrt(151) * density(01) * density(43) 
  rrt(152) = rrt(152) * density(01) * density(36) 
  rrt(153) = rrt(153) * density(01) * density(33) 
  rrt(154) = rrt(154) * density(01) * density(34) 
  rrt(155) = rrt(155) * density(01) * density(33) 
  rrt(156) = rrt(156) * density(01) * density(34) 
  rrt(157) = rrt(157) * density(01) * density(33) 
  rrt(158) = rrt(158) * density(01) * density(34) 
  rrt(159) = rrt(159) * density(01) * density(33) 
  rrt(160) = rrt(160) * density(01) * density(34) 
  rrt(161) = rrt(161) * density(01) * density(33) 
  rrt(162) = rrt(162) * density(01) * density(34) 
  rrt(163) = rrt(163) * density(01) * density(33) 
  rrt(164) = rrt(164) * density(01) * density(34) 
  rrt(165) = rrt(165) * density(01) * density(16) 
  rrt(166) = rrt(166) * density(01) * density(04) 
  rrt(167) = rrt(167) * density(01) * density(16) 
  rrt(168) = rrt(168) * density(01) * density(04) 
  rrt(169) = rrt(169) * density(01) * density(16) 
  rrt(170) = rrt(170) * density(01) * density(04) 
  rrt(171) = rrt(171) * density(01) * density(16) 
  rrt(172) = rrt(172) * density(01) * density(04) 
  rrt(173) = rrt(173) * density(01) * density(16) 
  rrt(174) = rrt(174) * density(01) * density(04) 
  rrt(175) = rrt(175) * density(01) * density(21) 
  rrt(176) = rrt(176) * density(01) * density(31) 
  rrt(177) = rrt(177) * density(01) * density(29) 
  rrt(178) = rrt(178) * density(01) * density(37) 
  rrt(179) = rrt(179) * density(01) * density(40) 
  rrt(180) = rrt(180) * density(01) * density(37) 
  rrt(181) = rrt(181) * density(01) * density(40) 
  rrt(182) = rrt(182) * density(01) * density(37) 
  rrt(183) = rrt(183) * density(01) * density(40) 
  rrt(184) = rrt(184) * density(01) * density(37) 
  rrt(185) = rrt(185) * density(01) * density(40) 
  rrt(186) = rrt(186) * density(01) * density(37) 
  rrt(187) = rrt(187) * density(01) * density(40) 
  rrt(188) = rrt(188) * density(01) * density(37) 
  rrt(189) = rrt(189) * density(01) * density(40) 
  rrt(190) = rrt(190) * density(01) * density(33) 
  rrt(191) = rrt(191) * density(01) * density(34) 
  rrt(192) = rrt(192) * density(01) * density(33) 
  rrt(193) = rrt(193) * density(01) * density(34) 
  rrt(194) = rrt(194) * density(01) * density(33) 
  rrt(195) = rrt(195) * density(01) * density(34) 
  rrt(196) = rrt(196) * density(01) * density(33) 
  rrt(197) = rrt(197) * density(01) * density(34) 
  rrt(198) = rrt(198) * density(01) * density(33) 
  rrt(199) = rrt(199) * density(01) * density(34) 
  rrt(200) = rrt(200) * density(01) * density(33) 
  rrt(201) = rrt(201) * density(01) * density(34) 
  rrt(202) = rrt(202) * density(01) * density(33) 
  rrt(203) = rrt(203) * density(01) * density(34) 
  rrt(204) = rrt(204) * density(01) * density(16) 
  rrt(205) = rrt(205) * density(01) * density(04) 
  rrt(206) = rrt(206) * density(01) * density(16) 
  rrt(207) = rrt(207) * density(01) * density(04) 
  rrt(208) = rrt(208) * density(01) * density(16) 
  rrt(209) = rrt(209) * density(01) * density(04) 
  rrt(210) = rrt(210) * density(01) * density(16) 
  rrt(211) = rrt(211) * density(01) * density(04) 
  rrt(212) = rrt(212) * density(01) * density(16) 
  rrt(213) = rrt(213) * density(01) * density(04) 
  rrt(214) = rrt(214) * density(01) * density(21) 
  rrt(215) = rrt(215) * density(01) * density(31) 
  rrt(216) = rrt(216) * density(01) * density(29) 
  rrt(217) = rrt(217) * density(01) * density(21) 
  rrt(218) = rrt(218) * density(01) * density(31) 
  rrt(219) = rrt(219) * density(01) * density(29) 
  rrt(220) = rrt(220) * density(01) * density(37) 
  rrt(221) = rrt(221) * density(01) * density(40) 
  rrt(222) = rrt(222) * density(01) * density(37) 
  rrt(223) = rrt(223) * density(01) * density(40) 
  rrt(224) = rrt(224) * density(01) * density(37) 
  rrt(225) = rrt(225) * density(01) * density(40) 
  rrt(226) = rrt(226) * density(01) * density(37) 
  rrt(227) = rrt(227) * density(01) * density(40) 
  rrt(228) = rrt(228) * density(01) * density(37) 
  rrt(229) = rrt(229) * density(01) * density(40) 
  rrt(230) = rrt(230) * density(01) * density(37) 
  rrt(231) = rrt(231) * density(01) * density(40) 
  rrt(232) = rrt(232) * density(01) * density(37) 
  rrt(233) = rrt(233) * density(01) * density(40) 
  rrt(234) = rrt(234) * density(01) * density(37) 
  rrt(235) = rrt(235) * density(01) * density(40) 
  rrt(236) = rrt(236) * density(01) * density(37) 
  rrt(237) = rrt(237) * density(01) * density(40) 
  rrt(238) = rrt(238) * density(01) * density(02) 
  rrt(239) = rrt(239) * density(01) * density(02) 
  rrt(240) = rrt(240) * density(01)**2 * density(25) 
  rrt(241) = rrt(241) * density(01)**2 * density(48) 
  rrt(242) = rrt(242) * density(01)**2 * density(27) 
  rrt(243) = rrt(243) * density(01)**2 * density(28) 
  rrt(244) = rrt(244) * density(01)**2 * density(15) 
  rrt(245) = rrt(245) * density(01)**2 * density(44) 
  rrt(246) = rrt(246) * density(01)**2 * density(32) 
  rrt(247) = rrt(247) * density(01)**2 * density(06) 
  rrt(248) = rrt(248) * density(01)**2 * density(14) 
  rrt(249) = rrt(249) * density(01)**2 * density(22) 
  rrt(250) = rrt(250) * density(01)**2 * density(30) 
  rrt(251) = rrt(251) * density(01)**2 * density(49) 
  rrt(252) = rrt(252) * density(01)**2 * density(08) 
  rrt(253) = rrt(253) * density(01)**2 * density(47) 
  rrt(254) = rrt(254) * density(01)**2 * density(07) 
  rrt(255) = rrt(255) * density(01)**2 * density(45) 
  rrt(256) = rrt(256) * density(01)**2 * density(41) 
  rrt(257) = rrt(257) * density(01)**2 * density(41) 
  rrt(258) = rrt(258) * density(01)**2 * density(25) 
  rrt(259) = rrt(259) * density(01)**2 * density(25) 
  rrt(260) = rrt(260) * density(01)**2 * density(25) 
  rrt(261) = rrt(261) * density(01)**2 * density(48) 
  rrt(262) = rrt(262) * density(01)**2 * density(48) 
  rrt(263) = rrt(263) * density(01)**2 * density(27) 
  rrt(264) = rrt(264) * density(01)**2 * density(15) 
  rrt(265) = rrt(265) * density(01)**2 * density(15) 
  rrt(266) = rrt(266) * density(01)**2 * density(44) 
  rrt(267) = rrt(267) * density(01)**2 * density(44) 
  rrt(268) = rrt(268) * density(01)**2 * density(44) 
  rrt(269) = rrt(269) * density(01)**2 * density(44) 
  rrt(270) = rrt(270) * density(01)**2 * density(44) 
  rrt(271) = rrt(271) * density(01)**2 * density(32) 
  rrt(272) = rrt(272) * density(01)**2 * density(32) 
  rrt(273) = rrt(273) * density(01)**2 * density(06) 
  rrt(274) = rrt(274) * density(01)**2 * density(14) 
  rrt(275) = rrt(275) * density(01)**2 * density(45) 
  rrt(276) = rrt(276) * density(01)**2 * density(26) 
  rrt(277) = rrt(277) * density(02) * density(45) 
  rrt(278) = rrt(278) * density(10) * density(41) 
  rrt(279) = rrt(279) * density(41) * density(46) 
  rrt(280) = rrt(280) * density(20) * density(41) 
  rrt(281) = rrt(281) * density(35) * density(41) 
  rrt(282) = rrt(282) * density(19) * density(41) 
  rrt(283) = rrt(283) * density(13) * density(41) 
  rrt(284) = rrt(284) * density(05) * density(25) 
  rrt(285) = rrt(285) * density(20) * density(25) 
  rrt(286) = rrt(286) * density(25) * density(35) 
  rrt(287) = rrt(287) * density(25) * density(35) 
  rrt(288) = rrt(288) * density(19) * density(25) 
  rrt(289) = rrt(289) * density(19) * density(25) 
  rrt(290) = rrt(290) * density(02) * density(25) 
  rrt(291) = rrt(291) * density(13) * density(25) 
  rrt(292) = rrt(292) * density(05) * density(48) 
  rrt(293) = rrt(293) * density(05) * density(48) 
  rrt(294) = rrt(294) * density(10) * density(48) 
  rrt(295) = rrt(295) * density(46) * density(48) 
  rrt(296) = rrt(296) * density(20) * density(48) 
  rrt(297) = rrt(297) * density(35) * density(48) 
  rrt(298) = rrt(298) * density(48) * density(50) 
  rrt(299) = rrt(299) * density(05) * density(27) 
  rrt(300) = rrt(300) * density(05) * density(27) 
  rrt(301) = rrt(301) * density(05) * density(27) 
  rrt(302) = rrt(302) * density(05) * density(27) 
  rrt(303) = rrt(303) * density(05) * density(27) 
  rrt(304) = rrt(304) * density(02) * density(27) 
  rrt(305) = rrt(305) * density(05) * density(28) 
  rrt(306) = rrt(306) * density(05) * density(28) 
  rrt(307) = rrt(307) * density(05) * density(28) 
  rrt(308) = rrt(308) * density(02) * density(28) 
  rrt(309) = rrt(309) * density(15) * density(35) 
  rrt(310) = rrt(310) * density(15) * density(19) 
  rrt(311) = rrt(311) * density(13) * density(15) 
  rrt(312) = rrt(312) * density(13) * density(44) 
  rrt(313) = rrt(313) * density(32) * density(50) 
  rrt(314) = rrt(314) * density(32) * density(50) 
  rrt(315) = rrt(315) * density(13) * density(32) 
  rrt(316) = rrt(316) * density(06) * density(20) 
  rrt(317) = rrt(317) * density(06) * density(35) 
  rrt(318) = rrt(318) * density(06) * density(18) 
  rrt(319) = rrt(319) * density(06) * density(13) 
  rrt(320) = rrt(320) * density(05) * density(14) 
  rrt(321) = rrt(321) * density(14) * density(20) 
  rrt(322) = rrt(322) * density(14) * density(20) 
  rrt(323) = rrt(323) * density(14) * density(35) 
  rrt(324) = rrt(324) * density(14) * density(50) 
  rrt(325) = rrt(325) * density(02) * density(14) 
  rrt(326) = rrt(326) * density(05) * density(26) 
  rrt(327) = rrt(327) * density(24) * density(26) 
  rrt(328) = rrt(328) * density(10) * density(26) 
  rrt(329) = rrt(329) * density(26) * density(46) 
  rrt(330) = rrt(330) * density(20) * density(26) 
  rrt(331) = rrt(331) * density(26) * density(38) 
  rrt(332) = rrt(332) * density(26) * density(35) 
  rrt(333) = rrt(333) * density(26) * density(35) 
  rrt(334) = rrt(334) * density(26) * density(50) 
  rrt(335) = rrt(335) * density(18) * density(26) 
  rrt(336) = rrt(336) * density(19) * density(26) 
  rrt(337) = rrt(337) * density(05) * density(45) 
  rrt(338) = rrt(338) * density(05) * density(45) 
  rrt(339) = rrt(339) * density(05) * density(45) 
  rrt(340) = rrt(340) * density(10) * density(45) 
  rrt(341) = rrt(341) * density(10) * density(45) 
  rrt(342) = rrt(342) * density(45) * density(46) 
  rrt(343) = rrt(343) * density(45) * density(46) 
  rrt(344) = rrt(344) * density(20) * density(45) 
  rrt(345) = rrt(345) * density(20) * density(45) 
  rrt(346) = rrt(346) * density(20) * density(45) 
  rrt(347) = rrt(347) * density(20) * density(45) 
  rrt(348) = rrt(348) * density(20) * density(45) 
  rrt(349) = rrt(349) * density(35) * density(45) 
  rrt(350) = rrt(350) * density(35) * density(45) 
  rrt(351) = rrt(351) * density(35) * density(45) 
  rrt(352) = rrt(352) * density(19) * density(45) 
  rrt(353) = rrt(353) * density(19) * density(45) 
  rrt(354) = rrt(354) * density(13) * density(45) 
  rrt(355) = rrt(355) * density(13) * density(45) 
  rrt(356) = rrt(356) * density(05) * density(07) 
  rrt(357) = rrt(357) * density(05) * density(07) 
  rrt(358) = rrt(358) * density(07) * density(24) 
  rrt(359) = rrt(359) * density(07) * density(10) 
  rrt(360) = rrt(360) * density(07) * density(10) 
  rrt(361) = rrt(361) * density(07) * density(46) 
  rrt(362) = rrt(362) * density(07) * density(20) 
  rrt(363) = rrt(363) * density(07) * density(20) 
  rrt(364) = rrt(364) * density(07) * density(20) 
  rrt(365) = rrt(365) * density(07) * density(38) 
  rrt(366) = rrt(366) * density(07) * density(38) 
  rrt(367) = rrt(367) * density(07) * density(35) 
  rrt(368) = rrt(368) * density(07) * density(35) 
  rrt(369) = rrt(369) * density(07) * density(35) 
  rrt(370) = rrt(370) * density(07) * density(50) 
  rrt(371) = rrt(371) * density(07) * density(50) 
  rrt(372) = rrt(372) * density(07) * density(19) 
  rrt(373) = rrt(373) * density(05) * density(10) 
  rrt(374) = rrt(374) * density(05) * density(46) 
  rrt(375) = rrt(375) * density(05) * density(38) 
  rrt(376) = rrt(376) * density(05) * density(50) 
  rrt(377) = rrt(377) * density(05) * density(18) 
  rrt(378) = rrt(378) * density(05) * density(23) 
  rrt(379) = rrt(379) * density(03) * density(05) 
  rrt(380) = rrt(380) * density(05) * density(13) 
  rrt(381) = rrt(381) * density(05) * density(24) 
  rrt(382) = rrt(382) * density(05) * density(10) 
  rrt(383) = rrt(383) * density(05) 
  rrt(384) = rrt(384) * density(24) 
  rrt(385) = rrt(385) * density(24) 
  rrt(386) = rrt(386) * density(24) * density(38) 
  rrt(387) = rrt(387) * density(10)**2 
  rrt(388) = rrt(388) * density(10) * density(38) 
  rrt(389) = rrt(389) * density(10) * density(50) 
  rrt(390) = rrt(390) * density(10) * density(18) 
  rrt(391) = rrt(391) * density(10) * density(12) 
  rrt(392) = rrt(392) * density(10) * density(23) 
  rrt(393) = rrt(393) * density(10) * density(23) 
  rrt(394) = rrt(394) * density(09) * density(10) 
  rrt(395) = rrt(395) * density(02) * density(10) 
  rrt(396) = rrt(396) * density(10) * density(13) 
  rrt(397) = rrt(397) * density(10) 
  rrt(398) = rrt(398) * density(10)**2 
  rrt(399) = rrt(399) * density(10) * density(13) 
  rrt(400) = rrt(400) * density(20) * density(46) 
  rrt(401) = rrt(401) * density(20) * density(46) 
  rrt(402) = rrt(402) * density(02) * density(46) 
  rrt(403) = rrt(403) * density(24) * density(46) 
  rrt(404) = rrt(404) * density(10) * density(46) 
  rrt(405) = rrt(405) * density(02) * density(46) 
  rrt(406) = rrt(406) * density(46) * density(50) 
  rrt(407) = rrt(407) * density(20) * density(50) 
  rrt(408) = rrt(408) * density(20) * density(23) 
  rrt(409) = rrt(409) * density(03) * density(20) 
  rrt(410) = rrt(410) * density(13) * density(20) 
  rrt(411) = rrt(411) * density(13) * density(20) 
  rrt(412) = rrt(412) * density(20) 
  rrt(413) = rrt(413) * density(20) * density(46) 
  rrt(414) = rrt(414) * density(10) * density(20) 
  rrt(415) = rrt(415) * density(38)**2 
  rrt(416) = rrt(416) * density(35) * density(38) 
  rrt(417) = rrt(417) * density(19) * density(38) 
  rrt(418) = rrt(418) * density(18) * density(38) 
  rrt(419) = rrt(419) * density(12) * density(38) 
  rrt(420) = rrt(420) * density(23) * density(38) 
  rrt(421) = rrt(421) * density(23) * density(38) 
  rrt(422) = rrt(422) * density(09) * density(38) 
  rrt(423) = rrt(423) * density(02) * density(38) 
  rrt(424) = rrt(424) * density(13) * density(38) 
  rrt(425) = rrt(425) * density(13) * density(38) 
  rrt(426) = rrt(426) * density(13) * density(38) 
  rrt(427) = rrt(427) * density(38) 
  rrt(428) = rrt(428) * density(38) * density(50) 
  rrt(429) = rrt(429) * density(13) * density(35) 
  rrt(430) = rrt(430) * density(13) * density(35) 
  rrt(431) = rrt(431) * density(02) * density(35) 
  rrt(432) = rrt(432) * density(35) 
  rrt(433) = rrt(433) * density(09) * density(35) 
  rrt(434) = rrt(434) * density(19) * density(35) 
  rrt(435) = rrt(435) * density(09) * density(35) 
  rrt(436) = rrt(436) * density(35)**2 
  rrt(437) = rrt(437) * density(24) * density(35) 
  rrt(438) = rrt(438) * density(35) 
  rrt(439) = rrt(439) * density(35) * density(38) 
  rrt(440) = rrt(440) * density(02) * density(35) 
  rrt(441) = rrt(441) * density(10) * density(35) 
  rrt(442) = rrt(442) * density(35) * density(39) 
  rrt(443) = rrt(443) * density(50)**2 
  rrt(444) = rrt(444) * density(12) * density(50) 
  rrt(445) = rrt(445) * density(23) * density(50) 
  rrt(446) = rrt(446) * density(23) * density(50) 
  rrt(447) = rrt(447) * density(09) * density(50) 
  rrt(448) = rrt(448) * density(03) * density(50) 
  rrt(449) = rrt(449) * density(02) * density(50) 
  rrt(450) = rrt(450) * density(13) * density(50) 
  rrt(451) = rrt(451) * density(13) * density(50) 
  rrt(452) = rrt(452) * density(50) 
  rrt(453) = rrt(453) * density(13) * density(19) 
  rrt(454) = rrt(454) * density(02) * density(19) 
  rrt(455) = rrt(455) * density(02) * density(19) 
  rrt(456) = rrt(456) * density(19) * density(24) 
  rrt(457) = rrt(457) * density(19) * density(39) 
  rrt(458) = rrt(458) * density(03) * density(12) 
  rrt(459) = rrt(459) * density(12) * density(13) 
  rrt(460) = rrt(460) * density(12) 
  rrt(461) = rrt(461) * density(23)**2 
  rrt(462) = rrt(462) * density(09) * density(23) 
  rrt(463) = rrt(463) * density(03) * density(23) 
  rrt(464) = rrt(464) * density(02) * density(23) 
  rrt(465) = rrt(465) * density(13) * density(23) 
  rrt(466) = rrt(466) * density(13) * density(23) 
  rrt(467) = rrt(467) * density(13) * density(23) 
  rrt(468) = rrt(468) * density(23) 
  rrt(469) = rrt(469) * density(23) 
  rrt(470) = rrt(470) * density(09) * density(19) 
  rrt(471) = rrt(471) * density(09)**2 
  rrt(472) = rrt(472) * density(09) 
  rrt(473) = rrt(473) * density(09) * density(13) 
  rrt(474) = rrt(474) * density(09) * density(13) 
  rrt(475) = rrt(475) * density(09) 
  rrt(476) = rrt(476) * density(09) * density(24) 
  rrt(477) = rrt(477) * density(02) * density(03) 
  rrt(478) = rrt(478) * density(03) * density(13) 
  rrt(479) = rrt(479) * density(03) 
  rrt(480) = rrt(480) * density(39) 
  rrt(481) = rrt(481) * density(10) * density(39) 
  rrt(482) = rrt(482) * density(39) 
  rrt(483) = rrt(483) * density(11) 
  rrt(484) = rrt(484) * density(02) 
  rrt(485) = rrt(485) * density(13)**2 
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
  ydot(02) = +rrt(014)+rrt(015)+rrt(017)+rrt(021)+rrt(022)+rrt(025)+rrt(030)+rrt(031)+  2.d0 * rrt(032)+rrt(036)+rrt(038)+rrt(042)&
             +rrt(050)+  2.d0 * rrt(051)+rrt(056)+rrt(057)+rrt(061)+rrt(071)+rrt(072)+  2.d0 * rrt(073)+rrt(078)+rrt(079)+rrt(097)&
             +rrt(098)+  2.d0 * rrt(099)+rrt(106)+rrt(107)+rrt(117)+rrt(141)+rrt(142)+rrt(143)+rrt(144)+rrt(149)+rrt(150)+rrt(151)&
             +rrt(152)+rrt(155)+rrt(156)+rrt(157)+rrt(158)+  2.d0 * rrt(159)+  2.d0 * rrt(160)+rrt(167)+rrt(168)+rrt(180)+rrt(181)&
             +  2.d0 * rrt(182)+  2.d0 * rrt(183)+rrt(194)+rrt(195)+rrt(196)+rrt(197)+  2.d0 * rrt(198)+  2.d0 * rrt(199)+rrt(224)&
             +rrt(225)+rrt(226)+rrt(227)+  2.d0 * rrt(228)+  2.d0 * rrt(229)-rrt(238)-rrt(239)+rrt(255)+rrt(257)+rrt(260)+rrt(262)&
             +rrt(268)+rrt(276)-rrt(277)+rrt(280)+rrt(283)+rrt(285)-rrt(290)+rrt(291)+rrt(293)+rrt(294)+rrt(295)+rrt(301)+rrt(302)&
             +  2.d0 * rrt(303)-rrt(304)+rrt(306)+rrt(307)-rrt(308)+rrt(311)+rrt(312)+rrt(315)+rrt(319)-rrt(325)+rrt(326)+rrt(327)&
             +rrt(328)+rrt(329)+  2.d0 * rrt(330)+rrt(331)+rrt(332)+  2.d0 * rrt(333)+rrt(334)+rrt(335)+rrt(336)+rrt(338)+rrt(339)&
             +rrt(341)+rrt(343)+rrt(344)+rrt(345)+  2.d0 * rrt(346)+  2.d0 * rrt(347)+  3.d0 * rrt(348)+rrt(349)+rrt(350)&
             +  2.d0 * rrt(351)+rrt(353)+rrt(355)+rrt(357)+rrt(360)+rrt(362)+rrt(363)+  2.d0 * rrt(364)+rrt(365)+rrt(366)+rrt(368)&
             +rrt(369)+rrt(371)+rrt(380)+rrt(385)-rrt(395)+rrt(396)+rrt(398)-rrt(402)-rrt(405)+rrt(410)-rrt(423)+rrt(425)+rrt(429)&
             -rrt(431)+rrt(438)-rrt(440)-rrt(449)+rrt(450)-rrt(454)-rrt(455)+rrt(459)-rrt(464)+rrt(465)+rrt(473)-rrt(477)-rrt(484)&
             +rrt(485) 
  ydot(03) = +rrt(056)+rrt(060)-rrt(065)-rrt(066)-rrt(126)-rrt(127)-rrt(128)-rrt(129)-rrt(130)-rrt(131)-rrt(132)-rrt(133)+rrt(252)&
             -rrt(379)+rrt(394)-rrt(409)+rrt(422)+rrt(433)+rrt(447)-rrt(448)+rrt(456)+rrt(457)-rrt(458)+rrt(462)-rrt(463)+rrt(470)&
             +rrt(471)+rrt(472)+rrt(473)-rrt(477)-rrt(478)-rrt(479) 
  ydot(04) = +rrt(006)-rrt(166)-rrt(168)-rrt(170)-rrt(172)-rrt(174)-rrt(205)-rrt(207)-rrt(209)-rrt(211)-rrt(213) 
  ydot(05) = -rrt(001)-rrt(002)-rrt(013)-rrt(014)-rrt(015)-rrt(019)-rrt(020)-rrt(021)-rrt(022)+rrt(033)+rrt(039)+rrt(054)+rrt(059)&
             +rrt(064)+rrt(075)+rrt(082)+rrt(101)+rrt(110)+rrt(121)+rrt(161)+rrt(162)+rrt(188)+rrt(189)+rrt(202)+rrt(203)+rrt(232)&
             +rrt(233)+rrt(240)+rrt(278)+rrt(279)+rrt(280)+rrt(281)+rrt(282)-rrt(284)+rrt(285)+rrt(287)+rrt(289)-rrt(292)-rrt(293)&
             +rrt(296)+rrt(297)-rrt(299)-rrt(300)-rrt(301)-rrt(302)-rrt(303)-rrt(305)-rrt(306)-rrt(307)-rrt(320)-rrt(326)-rrt(337)&
             -rrt(338)-rrt(339)-rrt(356)-rrt(357)-rrt(373)-rrt(374)-rrt(375)-rrt(376)-rrt(377)-rrt(378)-rrt(379)-rrt(380)-rrt(381)&
             -rrt(382)-rrt(383)+rrt(411) 
  ydot(06) = +rrt(072)+rrt(078)+rrt(084)+rrt(088)+rrt(110)+rrt(120)+rrt(129)+rrt(135)+rrt(196)+rrt(197)+rrt(206)+rrt(207)-rrt(247)&
             -rrt(273)+rrt(282)+rrt(288)+rrt(294)+rrt(297)+rrt(298)+rrt(302)+rrt(306)+rrt(314)+rrt(315)-rrt(316)-rrt(317)-rrt(318)&
             -rrt(319)+rrt(320)+rrt(324)+rrt(325)+rrt(333)+rrt(336)+rrt(347)+rrt(350)+rrt(352)+rrt(364)+rrt(366)+rrt(368)+rrt(370) 
  ydot(07) = +rrt(092)-rrt(254)+rrt(355)-rrt(356)-rrt(357)-rrt(358)-rrt(359)-rrt(360)-rrt(361)-rrt(362)-rrt(363)-rrt(364)-rrt(365)&
             -rrt(366)-rrt(367)-rrt(368)-rrt(369)-rrt(370)-rrt(371)-rrt(372) 
  ydot(08) = +rrt(098)+rrt(106)+rrt(116)+rrt(126)+rrt(226)+rrt(227)-rrt(252) 
  ydot(09) = -rrt(012)+rrt(050)+rrt(055)-rrt(060)-rrt(061)-rrt(062)-rrt(063)-rrt(064)-rrt(115)-rrt(116)-rrt(117)-rrt(118)-rrt(119)&
             -rrt(120)-rrt(121)-rrt(122)-rrt(123)-rrt(124)-rrt(125)+rrt(180)+rrt(181)+rrt(251)+rrt(379)+rrt(393)-rrt(394)+rrt(400)&
             +rrt(409)+rrt(421)-rrt(422)-rrt(433)-rrt(435)+rrt(441)+rrt(442)+rrt(446)-rrt(447)+rrt(448)+rrt(457)+rrt(458)+rrt(461)&
             -rrt(462)+  2.d0 * rrt(463)+rrt(465)+rrt(468)-rrt(470)-  2.d0 * rrt(471)-rrt(472)-rrt(473)-rrt(474)-rrt(475)-rrt(476)&
             +rrt(477)+rrt(478)+rrt(482) 
  ydot(10) = +rrt(014)+rrt(016)-rrt(018)-rrt(026)-rrt(027)+rrt(033)+rrt(040)+  2.d0 * rrt(045)+rrt(047)+rrt(052)+rrt(062)+rrt(068)&
             +rrt(080)+rrt(086)+rrt(091)+rrt(108)+rrt(119)+rrt(129)+rrt(136)+rrt(141)+rrt(142)+rrt(161)+rrt(162)+  2.d0 * rrt(173)&
             +  2.d0 * rrt(174)+rrt(184)+rrt(185)+rrt(210)+rrt(211)+rrt(242)+rrt(257)+rrt(259)+rrt(261)+rrt(270)-rrt(278)-rrt(294)&
             -rrt(328)-rrt(340)-rrt(341)-rrt(359)-rrt(360)-rrt(373)-rrt(382)+rrt(384)+rrt(386)-  2.d0 * rrt(387)-rrt(388)-rrt(389)&
             -rrt(390)-rrt(391)-rrt(392)-rrt(393)-rrt(394)-rrt(395)-rrt(396)-rrt(397)-  2.d0 * rrt(398)-rrt(399)+rrt(402)-rrt(404)&
             +rrt(406)-rrt(414)-rrt(441)-rrt(481) 
  ydot(11) = -rrt(483) 
  ydot(12) = -rrt(010)-rrt(011)-rrt(049)-rrt(050)-rrt(051)-rrt(052)-rrt(053)-rrt(054)-rrt(095)-rrt(096)-rrt(097)-rrt(098)-rrt(099)&
             -rrt(100)-rrt(101)-rrt(102)-rrt(103)+rrt(249)+rrt(378)-rrt(391)+rrt(408)-rrt(419)+rrt(420)-rrt(444)+rrt(445)-rrt(458)&
             -rrt(459)-rrt(460)+rrt(461)+rrt(462)+rrt(464)+rrt(466) 
  ydot(13) = +rrt(013)+rrt(015)+rrt(016)+rrt(018)+rrt(020)+rrt(022)+rrt(024)+rrt(027)+rrt(029)+rrt(031)+rrt(035)+  2.d0 * rrt(037)&
             +rrt(038)+rrt(041)+  2.d0 * rrt(043)+rrt(046)+rrt(049)+rrt(055)+rrt(057)+rrt(060)+rrt(065)+rrt(070)+rrt(072)+rrt(077)&
             +rrt(079)+rrt(084)+rrt(089)+rrt(096)+rrt(098)+rrt(105)+rrt(107)+rrt(116)+rrt(127)+rrt(139)+rrt(140)+rrt(143)+rrt(144)&
             +rrt(147)+rrt(148)+rrt(151)+rrt(152)+rrt(153)+rrt(154)+rrt(157)+rrt(158)+rrt(165)+rrt(166)+  2.d0 * rrt(169)&
             +  2.d0 * rrt(170)+rrt(178)+rrt(179)+rrt(192)+rrt(193)+rrt(196)+rrt(197)+rrt(206)+rrt(207)+rrt(222)+rrt(223)+rrt(226)&
             +rrt(227)+  2.d0 * rrt(238)+rrt(254)+  2.d0 * rrt(256)+rrt(257)+rrt(258)+  2.d0 * rrt(259)+rrt(260)+rrt(261)+rrt(263)&
             +rrt(264)+  2.d0 * rrt(265)+rrt(266)+  2.d0 * rrt(267)+rrt(268)+  3.d0 * rrt(269)+rrt(271)+  2.d0 * rrt(272)+rrt(273)&
             +  2.d0 * rrt(275)+rrt(276)+rrt(277)-rrt(283)+rrt(290)-rrt(291)+rrt(300)+rrt(302)+rrt(304)+rrt(305)+rrt(307)+rrt(308)&
             -rrt(311)-rrt(312)-rrt(315)-rrt(319)+rrt(325)+rrt(337)+rrt(339)+rrt(340)+rrt(342)+rrt(345)+rrt(347)+rrt(350)+rrt(352)&
             -rrt(354)-rrt(355)+rrt(356)+rrt(358)+rrt(359)+rrt(361)+rrt(363)+rrt(366)+rrt(367)+rrt(369)+rrt(370)+rrt(372)+rrt(374)&
             -rrt(380)+rrt(381)+rrt(383)+rrt(384)+  2.d0 * rrt(387)+rrt(395)-rrt(396)+rrt(397)-rrt(399)+rrt(400)+rrt(402)+rrt(403)&
             +rrt(404)-rrt(410)-rrt(411)+rrt(423)-rrt(424)-rrt(425)-rrt(426)+rrt(427)-rrt(429)-rrt(430)+rrt(431)+rrt(432)+rrt(449)&
             -rrt(450)-rrt(451)+rrt(452)-rrt(453)+rrt(455)-rrt(459)+rrt(464)-rrt(465)-rrt(466)-rrt(467)+rrt(468)+rrt(472)-rrt(473)&
             -rrt(474)+rrt(477)-rrt(478)+  2.d0 * rrt(484)-  2.d0 * rrt(485) 
  ydot(14) = +rrt(073)+rrt(079)+rrt(089)+rrt(093)+rrt(121)+rrt(130)+rrt(136)+rrt(198)+rrt(199)+rrt(214)+rrt(215)+rrt(216)-rrt(248)&
             -rrt(274)+rrt(289)+rrt(295)+rrt(303)+rrt(307)+rrt(318)+rrt(319)-rrt(320)-rrt(321)-rrt(322)-rrt(323)-rrt(324)-rrt(325)&
             +rrt(335)+rrt(348)+rrt(351)+rrt(353)+rrt(369)+rrt(371)+rrt(372) 
  ydot(15) = +rrt(069)+rrt(190)+rrt(191)-rrt(244)-rrt(264)-rrt(265)-rrt(309)-rrt(310)-rrt(311)+rrt(331)+rrt(344) 
  ydot(16) = +rrt(005)-rrt(165)-rrt(167)-rrt(169)-rrt(171)-rrt(173)-rrt(204)-rrt(206)-rrt(208)-rrt(210)-rrt(212) 
  ydot(17) = +rrt(051)+rrt(057)+rrt(061)+rrt(065)-rrt(067)-rrt(068)-rrt(134)-rrt(135)-rrt(136)-rrt(137)-rrt(138)+rrt(182)+rrt(183)&
             +rrt(253) 
  ydot(18) = -rrt(318)-rrt(335)-rrt(377)-rrt(390)+rrt(417)-rrt(418) 
  ydot(19) = -rrt(007)-rrt(008)-rrt(009)+rrt(032)+rrt(038)+rrt(042)+rrt(043)+rrt(046)-rrt(048)+rrt(064)+rrt(066)+rrt(068)+rrt(092)&
             -rrt(093)-rrt(094)+rrt(122)+rrt(131)+rrt(137)+rrt(159)+rrt(160)+rrt(167)+rrt(168)+rrt(169)+rrt(170)+rrt(248)+rrt(268)&
             +rrt(269)+rrt(272)+rrt(273)-rrt(282)-rrt(288)-rrt(289)-rrt(310)+rrt(313)+rrt(317)+rrt(318)+rrt(323)+rrt(324)-rrt(336)&
             -rrt(352)-rrt(353)-rrt(372)+rrt(377)+rrt(387)+rrt(389)+rrt(390)+rrt(398)+rrt(404)+rrt(406)-rrt(417)+rrt(418)-rrt(434)&
             +rrt(438)+rrt(443)+rrt(445)+rrt(448)+rrt(450)+rrt(452)-rrt(453)-rrt(454)-rrt(455)-rrt(456)-rrt(457)-rrt(470)+rrt(479) 
  ydot(20) = -rrt(003)-rrt(004)-rrt(029)-rrt(030)-rrt(031)-rrt(032)-rrt(033)-rrt(034)+rrt(052)-rrt(069)-rrt(070)-rrt(071)-rrt(072)&
             -rrt(073)-rrt(074)-rrt(075)+rrt(103)+rrt(114)+rrt(184)+rrt(185)+rrt(236)+rrt(237)+rrt(244)-rrt(280)-rrt(285)-rrt(296)&
             +rrt(309)-rrt(316)-rrt(321)-rrt(322)-rrt(330)-rrt(344)-rrt(345)-rrt(346)-rrt(347)-rrt(348)-rrt(362)-rrt(363)-rrt(364)&
             +rrt(375)+rrt(381)+rrt(382)+rrt(386)-rrt(400)-rrt(401)-rrt(407)-rrt(408)-rrt(409)-rrt(410)-rrt(411)-rrt(412)-rrt(413)&
             -rrt(414)+rrt(415)+rrt(416)+rrt(417)+rrt(419)+rrt(421)+rrt(422)+rrt(423)+rrt(426)+rrt(440) 
  ydot(21) = +rrt(007)-rrt(175)-rrt(214)-rrt(217) 
  ydot(22) = +rrt(095)+rrt(220)+rrt(221)-rrt(249) 
  ydot(23) = +rrt(049)-rrt(055)-rrt(056)-rrt(057)-rrt(058)-rrt(059)-rrt(104)-rrt(105)-rrt(106)-rrt(107)-rrt(108)-rrt(109)-rrt(110)&
             -rrt(111)-rrt(112)-rrt(113)-rrt(114)+rrt(178)+rrt(179)+rrt(250)-rrt(378)+rrt(391)-rrt(392)-rrt(393)+rrt(401)-rrt(408)&
             +rrt(419)-rrt(420)-rrt(421)+rrt(435)+rrt(437)+rrt(442)+rrt(444)-rrt(445)-rrt(446)+rrt(458)+rrt(459)-  2.d0 * rrt(461)&
             -rrt(462)-rrt(463)-rrt(464)-rrt(465)-rrt(466)-rrt(467)-rrt(468)-rrt(469)+rrt(471)+rrt(474)+rrt(481) 
  ydot(24) = +rrt(013)-rrt(016)-rrt(017)-rrt(023)-rrt(024)-rrt(025)+  2.d0 * rrt(034)+rrt(040)+rrt(044)+rrt(053)+rrt(058)+rrt(063)&
             +rrt(066)+rrt(074)+rrt(081)+rrt(087)+rrt(100)+rrt(109)+rrt(120)+rrt(130)+rrt(139)+rrt(140)+  2.d0 * rrt(163)&
             +  2.d0 * rrt(164)+rrt(171)+rrt(172)+rrt(186)+rrt(187)+rrt(200)+rrt(201)+rrt(212)+rrt(213)+rrt(230)+rrt(231)+rrt(241)&
             +rrt(256)+rrt(258)+rrt(270)+rrt(284)+rrt(286)+rrt(288)+rrt(292)+rrt(298)+rrt(299)+rrt(320)-rrt(327)-rrt(358)&
             +  2.d0 * rrt(373)+rrt(375)+rrt(376)+rrt(377)+rrt(378)+rrt(379)+rrt(380)-rrt(381)+rrt(383)-rrt(384)-rrt(385)-rrt(386)&
             +rrt(388)+rrt(389)+rrt(391)+rrt(393)+rrt(394)+rrt(395)+rrt(399)-rrt(403)+rrt(405)+rrt(411)+  2.d0 * rrt(412)+rrt(413)&
             +rrt(414)+  2.d0 * rrt(424)-rrt(437)-rrt(456)+rrt(460)+rrt(467)+rrt(469)+rrt(475)-rrt(476)+rrt(479)+rrt(482)+rrt(483) 
  ydot(25) = +rrt(019)+rrt(111)+rrt(122)+rrt(145)+rrt(146)-rrt(240)-rrt(258)-rrt(259)-rrt(260)+rrt(283)-rrt(284)-rrt(285)-rrt(286)&
             -rrt(287)-rrt(288)-rrt(289)-rrt(290)-rrt(291)+rrt(292)+rrt(327)+rrt(338)+rrt(356) 
  ydot(26) = -rrt(276)+rrt(277)-rrt(326)-rrt(327)-rrt(328)-rrt(329)-rrt(330)-rrt(331)-rrt(332)-rrt(333)-rrt(334)-rrt(335)-rrt(336)&
             +rrt(354) 
  ydot(27) = +rrt(021)+rrt(024)+rrt(026)+rrt(075)+rrt(081)+rrt(086)+rrt(090)+rrt(103)+rrt(113)+rrt(124)+rrt(132)+rrt(137)+rrt(149)&
             +rrt(150)+rrt(202)+rrt(203)+rrt(210)+rrt(211)+rrt(236)+rrt(237)-rrt(242)-rrt(263)+rrt(279)-rrt(299)-rrt(300)-rrt(301)&
             -rrt(302)-rrt(303)-rrt(304)+rrt(308)+rrt(329)+rrt(341)+rrt(342)+rrt(359) 
  ydot(28) = +rrt(022)+rrt(025)+rrt(027)+rrt(028)+rrt(082)+rrt(087)+rrt(091)+rrt(094)+rrt(114)+rrt(125)+rrt(133)+rrt(138)+rrt(151)&
             +rrt(152)+rrt(212)+rrt(213)+rrt(217)+rrt(218)+rrt(219)-rrt(243)-rrt(305)-rrt(306)-rrt(307)-rrt(308)+rrt(343)+rrt(360)&
             +rrt(361) 
  ydot(29) = +rrt(009)-rrt(177)-rrt(216)-rrt(219) 
  ydot(30) = +rrt(096)+rrt(104)+rrt(222)+rrt(223)-rrt(250) 
  ydot(31) = +rrt(008)-rrt(176)-rrt(215)-rrt(218) 
  ydot(32) = +rrt(071)+rrt(077)+rrt(083)+rrt(101)+rrt(109)+rrt(119)+rrt(128)+rrt(194)+rrt(195)+rrt(204)+rrt(205)+rrt(232)+rrt(233)&
             -rrt(246)-rrt(271)-rrt(272)+rrt(285)+rrt(287)+rrt(301)+rrt(305)+rrt(309)+rrt(312)-rrt(313)-rrt(314)-rrt(315)+rrt(322)&
             +rrt(323)+rrt(334)+rrt(346)+rrt(349)+rrt(363)+rrt(365)+rrt(367) 
  ydot(33) = +rrt(003)-rrt(153)-rrt(155)-rrt(157)-rrt(159)-rrt(161)-rrt(163)-rrt(190)-rrt(192)-rrt(194)-rrt(196)-rrt(198)-rrt(200)&
             -rrt(202) 
  ydot(34) = +rrt(004)-rrt(154)-rrt(156)-rrt(158)-rrt(160)-rrt(162)-rrt(164)-rrt(191)-rrt(193)-rrt(195)-rrt(197)-rrt(199)-rrt(201)&
             -rrt(203) 
  ydot(35) = -rrt(005)-rrt(006)+rrt(030)+rrt(035)-rrt(041)-rrt(042)-rrt(043)-rrt(044)-rrt(045)+rrt(054)+rrt(058)+rrt(062)-rrt(083)&
             -rrt(084)-rrt(085)-rrt(086)-rrt(087)+rrt(112)+rrt(124)+rrt(133)+rrt(155)+rrt(156)+rrt(188)+rrt(189)+rrt(246)+rrt(265)&
             +rrt(266)-rrt(281)-rrt(286)-rrt(287)-rrt(297)-rrt(309)+rrt(314)+rrt(316)-rrt(317)+rrt(322)-rrt(323)-rrt(332)-rrt(333)&
             -rrt(349)-rrt(350)-rrt(351)-rrt(367)-rrt(368)-rrt(369)+rrt(374)+rrt(376)+rrt(388)+rrt(392)+rrt(407)+rrt(413)+rrt(415)&
             -rrt(416)+rrt(418)+rrt(420)+rrt(425)+rrt(427)+  2.d0 * rrt(428)-rrt(429)-rrt(430)-rrt(431)-rrt(432)-rrt(433)-rrt(434)&
             -rrt(435)-  2.d0 * rrt(436)-rrt(437)-rrt(438)-rrt(439)-rrt(440)-rrt(441)-rrt(442)+rrt(443)+rrt(444)+rrt(446)+rrt(447)&
             +rrt(449)+rrt(451)+rrt(454)+rrt(469)+rrt(480)+rrt(481) 
  ydot(36) = +rrt(002)-rrt(140)-rrt(142)-rrt(144)-rrt(146)-rrt(148)-rrt(150)-rrt(152) 
  ydot(37) = +rrt(010)-rrt(178)-rrt(180)-rrt(182)-rrt(184)-rrt(186)-rrt(188)-rrt(220)-rrt(222)-rrt(224)-rrt(226)-rrt(228)-rrt(230)&
             -rrt(232)-rrt(234)-rrt(236) 
  ydot(38) = +rrt(029)-rrt(035)-rrt(036)-rrt(037)-rrt(038)-rrt(039)-rrt(040)+rrt(053)-rrt(076)-rrt(077)-rrt(078)-rrt(079)-rrt(080)&
             -rrt(081)-rrt(082)+rrt(102)+rrt(113)+rrt(125)+rrt(153)+rrt(154)+rrt(186)+rrt(187)+rrt(234)+rrt(235)+rrt(245)+rrt(264)&
             -rrt(331)-rrt(365)-rrt(366)-rrt(375)-rrt(386)-rrt(388)+rrt(392)+rrt(407)+rrt(408)+rrt(409)+rrt(410)+rrt(414)&
             -  2.d0 * rrt(415)-rrt(416)-rrt(417)-rrt(418)-rrt(419)-rrt(420)-rrt(421)-rrt(422)-rrt(423)-rrt(424)-rrt(425)-rrt(426)&
             -rrt(427)-rrt(428)+rrt(430)+rrt(431)+rrt(433)+rrt(436)-rrt(439)+rrt(460)+rrt(467)+rrt(480) 
  ydot(39) = +rrt(439)-rrt(442)-rrt(457)+rrt(476)-rrt(480)-rrt(481)-rrt(482)+rrt(483) 
  ydot(40) = +rrt(011)-rrt(179)-rrt(181)-rrt(183)-rrt(185)-rrt(187)-rrt(189)-rrt(221)-rrt(223)-rrt(225)-rrt(227)-rrt(229)-rrt(231)&
             -rrt(233)-rrt(235)-rrt(237) 
  ydot(41) = -rrt(256)-rrt(257)-rrt(278)-rrt(279)-rrt(280)-rrt(281)-rrt(282)-rrt(283)+rrt(284)+rrt(290)+rrt(326)+rrt(337) 
  ydot(42) = +rrt(012) 
  ydot(43) = +rrt(001)-rrt(139)-rrt(141)-rrt(143)-rrt(145)-rrt(147)-rrt(149)-rrt(151) 
  ydot(44) = +rrt(070)+rrt(076)+rrt(100)+rrt(108)+rrt(118)+rrt(192)+rrt(193)+rrt(230)+rrt(231)-rrt(245)-rrt(266)-rrt(267)-rrt(268)&
             -rrt(269)-rrt(270)+rrt(280)+rrt(281)+rrt(286)+rrt(293)+rrt(296)+rrt(300)+rrt(310)+rrt(311)-rrt(312)+rrt(313)+rrt(316)&
             +rrt(317)+rrt(321)+rrt(330)+rrt(332)+rrt(345)+rrt(362) 
  ydot(45) = +rrt(239)-rrt(255)-rrt(275)-rrt(277)-rrt(337)-rrt(338)-rrt(339)-rrt(340)-rrt(341)-rrt(342)-rrt(343)-rrt(344)-rrt(345)&
             -rrt(346)-rrt(347)-rrt(348)-rrt(349)-rrt(350)-rrt(351)-rrt(352)-rrt(353)-rrt(354)-rrt(355) 
  ydot(46) = +rrt(015)+rrt(017)+rrt(018)-rrt(028)+rrt(039)+rrt(044)+rrt(047)+  2.d0 * rrt(048)+rrt(067)+rrt(085)+rrt(090)+rrt(094)&
             +rrt(118)+rrt(128)+rrt(135)+rrt(143)+rrt(144)+rrt(171)+rrt(172)+  2.d0 * rrt(175)+  2.d0 * rrt(176)+  2.d0 * rrt(177)&
             +rrt(208)+rrt(209)+rrt(217)+rrt(218)+rrt(219)+rrt(243)+rrt(260)+rrt(262)+rrt(263)+  2.d0 * rrt(274)-rrt(279)-rrt(295)&
             -rrt(329)-rrt(342)-rrt(343)-rrt(361)-rrt(374)+rrt(385)+rrt(390)+rrt(396)+rrt(397)-rrt(400)-rrt(401)-rrt(402)-rrt(403)&
             -rrt(404)-rrt(405)-rrt(406)-rrt(413) 
  ydot(47) = +rrt(099)+rrt(107)+rrt(117)+rrt(127)+rrt(134)+rrt(228)+rrt(229)-rrt(253) 
  ydot(48) = +rrt(020)+rrt(023)+rrt(074)+rrt(080)+rrt(085)+rrt(102)+rrt(112)+rrt(123)+rrt(131)+rrt(147)+rrt(148)+rrt(200)+rrt(201)&
             +rrt(208)+rrt(209)+rrt(234)+rrt(235)-rrt(241)-rrt(261)-rrt(262)+rrt(278)+rrt(291)-rrt(292)-rrt(293)-rrt(294)-rrt(295)&
             -rrt(296)-rrt(297)-rrt(298)+rrt(299)+rrt(304)+rrt(328)+rrt(339)+rrt(340)+rrt(357)+rrt(358) 
  ydot(49) = +rrt(097)+rrt(105)+rrt(115)+rrt(224)+rrt(225)-rrt(251) 
  ydot(50) = +rrt(031)+rrt(036)+rrt(037)+rrt(041)-rrt(046)-rrt(047)+rrt(059)+rrt(063)+rrt(067)-rrt(088)-rrt(089)-rrt(090)-rrt(091)&
             -rrt(092)+rrt(111)+rrt(123)+rrt(132)+rrt(138)+rrt(157)+rrt(158)+rrt(165)+rrt(166)+rrt(247)+rrt(267)+rrt(271)-rrt(298)&
             +rrt(310)-rrt(313)-rrt(314)+rrt(321)-rrt(324)-rrt(334)-rrt(370)-rrt(371)-rrt(376)-rrt(389)+rrt(403)-rrt(406)-rrt(407)&
             +rrt(416)-rrt(428)+rrt(429)+rrt(432)+  2.d0 * rrt(434)+rrt(435)+rrt(436)-  2.d0 * rrt(443)-rrt(444)-rrt(445)-rrt(446)&
             -rrt(447)-rrt(448)-rrt(449)-rrt(450)-rrt(451)-rrt(452)+rrt(453)+rrt(455)+rrt(470)+rrt(475) 
  if( ldensity_constant ) where( density_constant(:) ) ydot(1:species_max) = 0.0d0
  ydot(51) = 0.0d0
  if( lgas_heating ) then
    ydot(51) = ( ZDPlasKin_cfg(14)/k_B + ydot(51) ) / ( sum(density(1:species_max)) - density(species_electrons) ) &
            + eV_to_K * ZDPlasKin_cfg(11) * density(species_electrons)
    ydot(51) = ydot(51) * ZDPlasKin_cfg(13)
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
  if( lgas_heating ) ZDPlasKin_cfg(1) = y(51)
  density(:) = y(1:species_max)
  call ZDPlasKin_reac_rates(t)
  pd(05,01) = pd(05,01) - rrt(001) * density(05) 
  pd(05,05) = pd(05,05) - rrt(001) * density(01) 
  pd(43,01) = pd(43,01) + rrt(001) * density(05) 
  pd(43,05) = pd(43,05) + rrt(001) * density(01) 
  pd(05,01) = pd(05,01) - rrt(002) * density(05) 
  pd(05,05) = pd(05,05) - rrt(002) * density(01) 
  pd(36,01) = pd(36,01) + rrt(002) * density(05) 
  pd(36,05) = pd(36,05) + rrt(002) * density(01) 
  pd(20,01) = pd(20,01) - rrt(003) * density(20) 
  pd(20,20) = pd(20,20) - rrt(003) * density(01) 
  pd(33,01) = pd(33,01) + rrt(003) * density(20) 
  pd(33,20) = pd(33,20) + rrt(003) * density(01) 
  pd(20,01) = pd(20,01) - rrt(004) * density(20) 
  pd(20,20) = pd(20,20) - rrt(004) * density(01) 
  pd(34,01) = pd(34,01) + rrt(004) * density(20) 
  pd(34,20) = pd(34,20) + rrt(004) * density(01) 
  pd(16,01) = pd(16,01) + rrt(005) * density(35) 
  pd(16,35) = pd(16,35) + rrt(005) * density(01) 
  pd(35,01) = pd(35,01) - rrt(005) * density(35) 
  pd(35,35) = pd(35,35) - rrt(005) * density(01) 
  pd(04,01) = pd(04,01) + rrt(006) * density(35) 
  pd(04,35) = pd(04,35) + rrt(006) * density(01) 
  pd(35,01) = pd(35,01) - rrt(006) * density(35) 
  pd(35,35) = pd(35,35) - rrt(006) * density(01) 
  pd(19,01) = pd(19,01) - rrt(007) * density(19) 
  pd(19,19) = pd(19,19) - rrt(007) * density(01) 
  pd(21,01) = pd(21,01) + rrt(007) * density(19) 
  pd(21,19) = pd(21,19) + rrt(007) * density(01) 
  pd(19,01) = pd(19,01) - rrt(008) * density(19) 
  pd(19,19) = pd(19,19) - rrt(008) * density(01) 
  pd(31,01) = pd(31,01) + rrt(008) * density(19) 
  pd(31,19) = pd(31,19) + rrt(008) * density(01) 
  pd(19,01) = pd(19,01) - rrt(009) * density(19) 
  pd(19,19) = pd(19,19) - rrt(009) * density(01) 
  pd(29,01) = pd(29,01) + rrt(009) * density(19) 
  pd(29,19) = pd(29,19) + rrt(009) * density(01) 
  pd(12,01) = pd(12,01) - rrt(010) * density(12) 
  pd(12,12) = pd(12,12) - rrt(010) * density(01) 
  pd(37,01) = pd(37,01) + rrt(010) * density(12) 
  pd(37,12) = pd(37,12) + rrt(010) * density(01) 
  pd(12,01) = pd(12,01) - rrt(011) * density(12) 
  pd(12,12) = pd(12,12) - rrt(011) * density(01) 
  pd(40,01) = pd(40,01) + rrt(011) * density(12) 
  pd(40,12) = pd(40,12) + rrt(011) * density(01) 
  pd(09,01) = pd(09,01) - rrt(012) * density(09) 
  pd(09,09) = pd(09,09) - rrt(012) * density(01) 
  pd(42,01) = pd(42,01) + rrt(012) * density(09) 
  pd(42,09) = pd(42,09) + rrt(012) * density(01) 
  pd(05,01) = pd(05,01) - rrt(013) * density(05) 
  pd(05,05) = pd(05,05) - rrt(013) * density(01) 
  pd(13,01) = pd(13,01) + rrt(013) * density(05) 
  pd(13,05) = pd(13,05) + rrt(013) * density(01) 
  pd(24,01) = pd(24,01) + rrt(013) * density(05) 
  pd(24,05) = pd(24,05) + rrt(013) * density(01) 
  pd(02,01) = pd(02,01) + rrt(014) * density(05) 
  pd(02,05) = pd(02,05) + rrt(014) * density(01) 
  pd(05,01) = pd(05,01) - rrt(014) * density(05) 
  pd(05,05) = pd(05,05) - rrt(014) * density(01) 
  pd(10,01) = pd(10,01) + rrt(014) * density(05) 
  pd(10,05) = pd(10,05) + rrt(014) * density(01) 
  pd(02,01) = pd(02,01) + rrt(015) * density(05) 
  pd(02,05) = pd(02,05) + rrt(015) * density(01) 
  pd(05,01) = pd(05,01) - rrt(015) * density(05) 
  pd(05,05) = pd(05,05) - rrt(015) * density(01) 
  pd(13,01) = pd(13,01) + rrt(015) * density(05) 
  pd(13,05) = pd(13,05) + rrt(015) * density(01) 
  pd(46,01) = pd(46,01) + rrt(015) * density(05) 
  pd(46,05) = pd(46,05) + rrt(015) * density(01) 
  pd(10,01) = pd(10,01) + rrt(016) * density(24) 
  pd(10,24) = pd(10,24) + rrt(016) * density(01) 
  pd(13,01) = pd(13,01) + rrt(016) * density(24) 
  pd(13,24) = pd(13,24) + rrt(016) * density(01) 
  pd(24,01) = pd(24,01) - rrt(016) * density(24) 
  pd(24,24) = pd(24,24) - rrt(016) * density(01) 
  pd(02,01) = pd(02,01) + rrt(017) * density(24) 
  pd(02,24) = pd(02,24) + rrt(017) * density(01) 
  pd(24,01) = pd(24,01) - rrt(017) * density(24) 
  pd(24,24) = pd(24,24) - rrt(017) * density(01) 
  pd(46,01) = pd(46,01) + rrt(017) * density(24) 
  pd(46,24) = pd(46,24) + rrt(017) * density(01) 
  pd(10,01) = pd(10,01) - rrt(018) * density(10) 
  pd(10,10) = pd(10,10) - rrt(018) * density(01) 
  pd(13,01) = pd(13,01) + rrt(018) * density(10) 
  pd(13,10) = pd(13,10) + rrt(018) * density(01) 
  pd(46,01) = pd(46,01) + rrt(018) * density(10) 
  pd(46,10) = pd(46,10) + rrt(018) * density(01) 
  pd(01,01) = pd(01,01) + rrt(019) * density(05) 
  pd(01,05) = pd(01,05) + rrt(019) * density(01) 
  pd(05,01) = pd(05,01) - rrt(019) * density(05) 
  pd(05,05) = pd(05,05) - rrt(019) * density(01) 
  pd(25,01) = pd(25,01) + rrt(019) * density(05) 
  pd(25,05) = pd(25,05) + rrt(019) * density(01) 
  pd(01,01) = pd(01,01) + rrt(020) * density(05) 
  pd(01,05) = pd(01,05) + rrt(020) * density(01) 
  pd(05,01) = pd(05,01) - rrt(020) * density(05) 
  pd(05,05) = pd(05,05) - rrt(020) * density(01) 
  pd(13,01) = pd(13,01) + rrt(020) * density(05) 
  pd(13,05) = pd(13,05) + rrt(020) * density(01) 
  pd(48,01) = pd(48,01) + rrt(020) * density(05) 
  pd(48,05) = pd(48,05) + rrt(020) * density(01) 
  pd(01,01) = pd(01,01) + rrt(021) * density(05) 
  pd(01,05) = pd(01,05) + rrt(021) * density(01) 
  pd(02,01) = pd(02,01) + rrt(021) * density(05) 
  pd(02,05) = pd(02,05) + rrt(021) * density(01) 
  pd(05,01) = pd(05,01) - rrt(021) * density(05) 
  pd(05,05) = pd(05,05) - rrt(021) * density(01) 
  pd(27,01) = pd(27,01) + rrt(021) * density(05) 
  pd(27,05) = pd(27,05) + rrt(021) * density(01) 
  pd(01,01) = pd(01,01) + rrt(022) * density(05) 
  pd(01,05) = pd(01,05) + rrt(022) * density(01) 
  pd(02,01) = pd(02,01) + rrt(022) * density(05) 
  pd(02,05) = pd(02,05) + rrt(022) * density(01) 
  pd(05,01) = pd(05,01) - rrt(022) * density(05) 
  pd(05,05) = pd(05,05) - rrt(022) * density(01) 
  pd(13,01) = pd(13,01) + rrt(022) * density(05) 
  pd(13,05) = pd(13,05) + rrt(022) * density(01) 
  pd(28,01) = pd(28,01) + rrt(022) * density(05) 
  pd(28,05) = pd(28,05) + rrt(022) * density(01) 
  pd(01,01) = pd(01,01) + rrt(023) * density(24) 
  pd(01,24) = pd(01,24) + rrt(023) * density(01) 
  pd(24,01) = pd(24,01) - rrt(023) * density(24) 
  pd(24,24) = pd(24,24) - rrt(023) * density(01) 
  pd(48,01) = pd(48,01) + rrt(023) * density(24) 
  pd(48,24) = pd(48,24) + rrt(023) * density(01) 
  pd(01,01) = pd(01,01) + rrt(024) * density(24) 
  pd(01,24) = pd(01,24) + rrt(024) * density(01) 
  pd(13,01) = pd(13,01) + rrt(024) * density(24) 
  pd(13,24) = pd(13,24) + rrt(024) * density(01) 
  pd(24,01) = pd(24,01) - rrt(024) * density(24) 
  pd(24,24) = pd(24,24) - rrt(024) * density(01) 
  pd(27,01) = pd(27,01) + rrt(024) * density(24) 
  pd(27,24) = pd(27,24) + rrt(024) * density(01) 
  pd(01,01) = pd(01,01) + rrt(025) * density(24) 
  pd(01,24) = pd(01,24) + rrt(025) * density(01) 
  pd(02,01) = pd(02,01) + rrt(025) * density(24) 
  pd(02,24) = pd(02,24) + rrt(025) * density(01) 
  pd(24,01) = pd(24,01) - rrt(025) * density(24) 
  pd(24,24) = pd(24,24) - rrt(025) * density(01) 
  pd(28,01) = pd(28,01) + rrt(025) * density(24) 
  pd(28,24) = pd(28,24) + rrt(025) * density(01) 
  pd(01,01) = pd(01,01) + rrt(026) * density(10) 
  pd(01,10) = pd(01,10) + rrt(026) * density(01) 
  pd(10,01) = pd(10,01) - rrt(026) * density(10) 
  pd(10,10) = pd(10,10) - rrt(026) * density(01) 
  pd(27,01) = pd(27,01) + rrt(026) * density(10) 
  pd(27,10) = pd(27,10) + rrt(026) * density(01) 
  pd(01,01) = pd(01,01) + rrt(027) * density(10) 
  pd(01,10) = pd(01,10) + rrt(027) * density(01) 
  pd(10,01) = pd(10,01) - rrt(027) * density(10) 
  pd(10,10) = pd(10,10) - rrt(027) * density(01) 
  pd(13,01) = pd(13,01) + rrt(027) * density(10) 
  pd(13,10) = pd(13,10) + rrt(027) * density(01) 
  pd(28,01) = pd(28,01) + rrt(027) * density(10) 
  pd(28,10) = pd(28,10) + rrt(027) * density(01) 
  pd(01,01) = pd(01,01) + rrt(028) * density(46) 
  pd(01,46) = pd(01,46) + rrt(028) * density(01) 
  pd(28,01) = pd(28,01) + rrt(028) * density(46) 
  pd(28,46) = pd(28,46) + rrt(028) * density(01) 
  pd(46,01) = pd(46,01) - rrt(028) * density(46) 
  pd(46,46) = pd(46,46) - rrt(028) * density(01) 
  pd(13,01) = pd(13,01) + rrt(029) * density(20) 
  pd(13,20) = pd(13,20) + rrt(029) * density(01) 
  pd(20,01) = pd(20,01) - rrt(029) * density(20) 
  pd(20,20) = pd(20,20) - rrt(029) * density(01) 
  pd(38,01) = pd(38,01) + rrt(029) * density(20) 
  pd(38,20) = pd(38,20) + rrt(029) * density(01) 
  pd(02,01) = pd(02,01) + rrt(030) * density(20) 
  pd(02,20) = pd(02,20) + rrt(030) * density(01) 
  pd(20,01) = pd(20,01) - rrt(030) * density(20) 
  pd(20,20) = pd(20,20) - rrt(030) * density(01) 
  pd(35,01) = pd(35,01) + rrt(030) * density(20) 
  pd(35,20) = pd(35,20) + rrt(030) * density(01) 
  pd(02,01) = pd(02,01) + rrt(031) * density(20) 
  pd(02,20) = pd(02,20) + rrt(031) * density(01) 
  pd(13,01) = pd(13,01) + rrt(031) * density(20) 
  pd(13,20) = pd(13,20) + rrt(031) * density(01) 
  pd(20,01) = pd(20,01) - rrt(031) * density(20) 
  pd(20,20) = pd(20,20) - rrt(031) * density(01) 
  pd(50,01) = pd(50,01) + rrt(031) * density(20) 
  pd(50,20) = pd(50,20) + rrt(031) * density(01) 
  pd(02,01) = pd(02,01) + rrt(032) * density(20) * 2.0d0
  pd(02,20) = pd(02,20) + rrt(032) * density(01) * 2.0d0
  pd(19,01) = pd(19,01) + rrt(032) * density(20) 
  pd(19,20) = pd(19,20) + rrt(032) * density(01) 
  pd(20,01) = pd(20,01) - rrt(032) * density(20) 
  pd(20,20) = pd(20,20) - rrt(032) * density(01) 
  pd(05,01) = pd(05,01) + rrt(033) * density(20) 
  pd(05,20) = pd(05,20) + rrt(033) * density(01) 
  pd(10,01) = pd(10,01) + rrt(033) * density(20) 
  pd(10,20) = pd(10,20) + rrt(033) * density(01) 
  pd(20,01) = pd(20,01) - rrt(033) * density(20) 
  pd(20,20) = pd(20,20) - rrt(033) * density(01) 
  pd(20,01) = pd(20,01) - rrt(034) * density(20) 
  pd(20,20) = pd(20,20) - rrt(034) * density(01) 
  pd(24,01) = pd(24,01) + rrt(034) * density(20) * 2.0d0
  pd(24,20) = pd(24,20) + rrt(034) * density(01) * 2.0d0
  pd(13,01) = pd(13,01) + rrt(035) * density(38) 
  pd(13,38) = pd(13,38) + rrt(035) * density(01) 
  pd(35,01) = pd(35,01) + rrt(035) * density(38) 
  pd(35,38) = pd(35,38) + rrt(035) * density(01) 
  pd(38,01) = pd(38,01) - rrt(035) * density(38) 
  pd(38,38) = pd(38,38) - rrt(035) * density(01) 
  pd(02,01) = pd(02,01) + rrt(036) * density(38) 
  pd(02,38) = pd(02,38) + rrt(036) * density(01) 
  pd(38,01) = pd(38,01) - rrt(036) * density(38) 
  pd(38,38) = pd(38,38) - rrt(036) * density(01) 
  pd(50,01) = pd(50,01) + rrt(036) * density(38) 
  pd(50,38) = pd(50,38) + rrt(036) * density(01) 
  pd(13,01) = pd(13,01) + rrt(037) * density(38) * 2.0d0
  pd(13,38) = pd(13,38) + rrt(037) * density(01) * 2.0d0
  pd(38,01) = pd(38,01) - rrt(037) * density(38) 
  pd(38,38) = pd(38,38) - rrt(037) * density(01) 
  pd(50,01) = pd(50,01) + rrt(037) * density(38) 
  pd(50,38) = pd(50,38) + rrt(037) * density(01) 
  pd(02,01) = pd(02,01) + rrt(038) * density(38) 
  pd(02,38) = pd(02,38) + rrt(038) * density(01) 
  pd(13,01) = pd(13,01) + rrt(038) * density(38) 
  pd(13,38) = pd(13,38) + rrt(038) * density(01) 
  pd(19,01) = pd(19,01) + rrt(038) * density(38) 
  pd(19,38) = pd(19,38) + rrt(038) * density(01) 
  pd(38,01) = pd(38,01) - rrt(038) * density(38) 
  pd(38,38) = pd(38,38) - rrt(038) * density(01) 
  pd(05,01) = pd(05,01) + rrt(039) * density(38) 
  pd(05,38) = pd(05,38) + rrt(039) * density(01) 
  pd(38,01) = pd(38,01) - rrt(039) * density(38) 
  pd(38,38) = pd(38,38) - rrt(039) * density(01) 
  pd(46,01) = pd(46,01) + rrt(039) * density(38) 
  pd(46,38) = pd(46,38) + rrt(039) * density(01) 
  pd(10,01) = pd(10,01) + rrt(040) * density(38) 
  pd(10,38) = pd(10,38) + rrt(040) * density(01) 
  pd(24,01) = pd(24,01) + rrt(040) * density(38) 
  pd(24,38) = pd(24,38) + rrt(040) * density(01) 
  pd(38,01) = pd(38,01) - rrt(040) * density(38) 
  pd(38,38) = pd(38,38) - rrt(040) * density(01) 
  pd(13,01) = pd(13,01) + rrt(041) * density(35) 
  pd(13,35) = pd(13,35) + rrt(041) * density(01) 
  pd(35,01) = pd(35,01) - rrt(041) * density(35) 
  pd(35,35) = pd(35,35) - rrt(041) * density(01) 
  pd(50,01) = pd(50,01) + rrt(041) * density(35) 
  pd(50,35) = pd(50,35) + rrt(041) * density(01) 
  pd(02,01) = pd(02,01) + rrt(042) * density(35) 
  pd(02,35) = pd(02,35) + rrt(042) * density(01) 
  pd(19,01) = pd(19,01) + rrt(042) * density(35) 
  pd(19,35) = pd(19,35) + rrt(042) * density(01) 
  pd(35,01) = pd(35,01) - rrt(042) * density(35) 
  pd(35,35) = pd(35,35) - rrt(042) * density(01) 
  pd(13,01) = pd(13,01) + rrt(043) * density(35) * 2.0d0
  pd(13,35) = pd(13,35) + rrt(043) * density(01) * 2.0d0
  pd(19,01) = pd(19,01) + rrt(043) * density(35) 
  pd(19,35) = pd(19,35) + rrt(043) * density(01) 
  pd(35,01) = pd(35,01) - rrt(043) * density(35) 
  pd(35,35) = pd(35,35) - rrt(043) * density(01) 
  pd(24,01) = pd(24,01) + rrt(044) * density(35) 
  pd(24,35) = pd(24,35) + rrt(044) * density(01) 
  pd(35,01) = pd(35,01) - rrt(044) * density(35) 
  pd(35,35) = pd(35,35) - rrt(044) * density(01) 
  pd(46,01) = pd(46,01) + rrt(044) * density(35) 
  pd(46,35) = pd(46,35) + rrt(044) * density(01) 
  pd(10,01) = pd(10,01) + rrt(045) * density(35) * 2.0d0
  pd(10,35) = pd(10,35) + rrt(045) * density(01) * 2.0d0
  pd(35,01) = pd(35,01) - rrt(045) * density(35) 
  pd(35,35) = pd(35,35) - rrt(045) * density(01) 
  pd(13,01) = pd(13,01) + rrt(046) * density(50) 
  pd(13,50) = pd(13,50) + rrt(046) * density(01) 
  pd(19,01) = pd(19,01) + rrt(046) * density(50) 
  pd(19,50) = pd(19,50) + rrt(046) * density(01) 
  pd(50,01) = pd(50,01) - rrt(046) * density(50) 
  pd(50,50) = pd(50,50) - rrt(046) * density(01) 
  pd(10,01) = pd(10,01) + rrt(047) * density(50) 
  pd(10,50) = pd(10,50) + rrt(047) * density(01) 
  pd(46,01) = pd(46,01) + rrt(047) * density(50) 
  pd(46,50) = pd(46,50) + rrt(047) * density(01) 
  pd(50,01) = pd(50,01) - rrt(047) * density(50) 
  pd(50,50) = pd(50,50) - rrt(047) * density(01) 
  pd(19,01) = pd(19,01) - rrt(048) * density(19) 
  pd(19,19) = pd(19,19) - rrt(048) * density(01) 
  pd(46,01) = pd(46,01) + rrt(048) * density(19) * 2.0d0
  pd(46,19) = pd(46,19) + rrt(048) * density(01) * 2.0d0
  pd(12,01) = pd(12,01) - rrt(049) * density(12) 
  pd(12,12) = pd(12,12) - rrt(049) * density(01) 
  pd(13,01) = pd(13,01) + rrt(049) * density(12) 
  pd(13,12) = pd(13,12) + rrt(049) * density(01) 
  pd(23,01) = pd(23,01) + rrt(049) * density(12) 
  pd(23,12) = pd(23,12) + rrt(049) * density(01) 
  pd(02,01) = pd(02,01) + rrt(050) * density(12) 
  pd(02,12) = pd(02,12) + rrt(050) * density(01) 
  pd(09,01) = pd(09,01) + rrt(050) * density(12) 
  pd(09,12) = pd(09,12) + rrt(050) * density(01) 
  pd(12,01) = pd(12,01) - rrt(050) * density(12) 
  pd(12,12) = pd(12,12) - rrt(050) * density(01) 
  pd(02,01) = pd(02,01) + rrt(051) * density(12) * 2.0d0
  pd(02,12) = pd(02,12) + rrt(051) * density(01) * 2.0d0
  pd(12,01) = pd(12,01) - rrt(051) * density(12) 
  pd(12,12) = pd(12,12) - rrt(051) * density(01) 
  pd(17,01) = pd(17,01) + rrt(051) * density(12) 
  pd(17,12) = pd(17,12) + rrt(051) * density(01) 
  pd(10,01) = pd(10,01) + rrt(052) * density(12) 
  pd(10,12) = pd(10,12) + rrt(052) * density(01) 
  pd(12,01) = pd(12,01) - rrt(052) * density(12) 
  pd(12,12) = pd(12,12) - rrt(052) * density(01) 
  pd(20,01) = pd(20,01) + rrt(052) * density(12) 
  pd(20,12) = pd(20,12) + rrt(052) * density(01) 
  pd(12,01) = pd(12,01) - rrt(053) * density(12) 
  pd(12,12) = pd(12,12) - rrt(053) * density(01) 
  pd(24,01) = pd(24,01) + rrt(053) * density(12) 
  pd(24,12) = pd(24,12) + rrt(053) * density(01) 
  pd(38,01) = pd(38,01) + rrt(053) * density(12) 
  pd(38,12) = pd(38,12) + rrt(053) * density(01) 
  pd(05,01) = pd(05,01) + rrt(054) * density(12) 
  pd(05,12) = pd(05,12) + rrt(054) * density(01) 
  pd(12,01) = pd(12,01) - rrt(054) * density(12) 
  pd(12,12) = pd(12,12) - rrt(054) * density(01) 
  pd(35,01) = pd(35,01) + rrt(054) * density(12) 
  pd(35,12) = pd(35,12) + rrt(054) * density(01) 
  pd(09,01) = pd(09,01) + rrt(055) * density(23) 
  pd(09,23) = pd(09,23) + rrt(055) * density(01) 
  pd(13,01) = pd(13,01) + rrt(055) * density(23) 
  pd(13,23) = pd(13,23) + rrt(055) * density(01) 
  pd(23,01) = pd(23,01) - rrt(055) * density(23) 
  pd(23,23) = pd(23,23) - rrt(055) * density(01) 
  pd(02,01) = pd(02,01) + rrt(056) * density(23) 
  pd(02,23) = pd(02,23) + rrt(056) * density(01) 
  pd(03,01) = pd(03,01) + rrt(056) * density(23) 
  pd(03,23) = pd(03,23) + rrt(056) * density(01) 
  pd(23,01) = pd(23,01) - rrt(056) * density(23) 
  pd(23,23) = pd(23,23) - rrt(056) * density(01) 
  pd(02,01) = pd(02,01) + rrt(057) * density(23) 
  pd(02,23) = pd(02,23) + rrt(057) * density(01) 
  pd(13,01) = pd(13,01) + rrt(057) * density(23) 
  pd(13,23) = pd(13,23) + rrt(057) * density(01) 
  pd(17,01) = pd(17,01) + rrt(057) * density(23) 
  pd(17,23) = pd(17,23) + rrt(057) * density(01) 
  pd(23,01) = pd(23,01) - rrt(057) * density(23) 
  pd(23,23) = pd(23,23) - rrt(057) * density(01) 
  pd(23,01) = pd(23,01) - rrt(058) * density(23) 
  pd(23,23) = pd(23,23) - rrt(058) * density(01) 
  pd(24,01) = pd(24,01) + rrt(058) * density(23) 
  pd(24,23) = pd(24,23) + rrt(058) * density(01) 
  pd(35,01) = pd(35,01) + rrt(058) * density(23) 
  pd(35,23) = pd(35,23) + rrt(058) * density(01) 
  pd(05,01) = pd(05,01) + rrt(059) * density(23) 
  pd(05,23) = pd(05,23) + rrt(059) * density(01) 
  pd(23,01) = pd(23,01) - rrt(059) * density(23) 
  pd(23,23) = pd(23,23) - rrt(059) * density(01) 
  pd(50,01) = pd(50,01) + rrt(059) * density(23) 
  pd(50,23) = pd(50,23) + rrt(059) * density(01) 
  pd(03,01) = pd(03,01) + rrt(060) * density(09) 
  pd(03,09) = pd(03,09) + rrt(060) * density(01) 
  pd(09,01) = pd(09,01) - rrt(060) * density(09) 
  pd(09,09) = pd(09,09) - rrt(060) * density(01) 
  pd(13,01) = pd(13,01) + rrt(060) * density(09) 
  pd(13,09) = pd(13,09) + rrt(060) * density(01) 
  pd(02,01) = pd(02,01) + rrt(061) * density(09) 
  pd(02,09) = pd(02,09) + rrt(061) * density(01) 
  pd(09,01) = pd(09,01) - rrt(061) * density(09) 
  pd(09,09) = pd(09,09) - rrt(061) * density(01) 
  pd(17,01) = pd(17,01) + rrt(061) * density(09) 
  pd(17,09) = pd(17,09) + rrt(061) * density(01) 
  pd(09,01) = pd(09,01) - rrt(062) * density(09) 
  pd(09,09) = pd(09,09) - rrt(062) * density(01) 
  pd(10,01) = pd(10,01) + rrt(062) * density(09) 
  pd(10,09) = pd(10,09) + rrt(062) * density(01) 
  pd(35,01) = pd(35,01) + rrt(062) * density(09) 
  pd(35,09) = pd(35,09) + rrt(062) * density(01) 
  pd(09,01) = pd(09,01) - rrt(063) * density(09) 
  pd(09,09) = pd(09,09) - rrt(063) * density(01) 
  pd(24,01) = pd(24,01) + rrt(063) * density(09) 
  pd(24,09) = pd(24,09) + rrt(063) * density(01) 
  pd(50,01) = pd(50,01) + rrt(063) * density(09) 
  pd(50,09) = pd(50,09) + rrt(063) * density(01) 
  pd(05,01) = pd(05,01) + rrt(064) * density(09) 
  pd(05,09) = pd(05,09) + rrt(064) * density(01) 
  pd(09,01) = pd(09,01) - rrt(064) * density(09) 
  pd(09,09) = pd(09,09) - rrt(064) * density(01) 
  pd(19,01) = pd(19,01) + rrt(064) * density(09) 
  pd(19,09) = pd(19,09) + rrt(064) * density(01) 
  pd(03,01) = pd(03,01) - rrt(065) * density(03) 
  pd(03,03) = pd(03,03) - rrt(065) * density(01) 
  pd(13,01) = pd(13,01) + rrt(065) * density(03) 
  pd(13,03) = pd(13,03) + rrt(065) * density(01) 
  pd(17,01) = pd(17,01) + rrt(065) * density(03) 
  pd(17,03) = pd(17,03) + rrt(065) * density(01) 
  pd(03,01) = pd(03,01) - rrt(066) * density(03) 
  pd(03,03) = pd(03,03) - rrt(066) * density(01) 
  pd(19,01) = pd(19,01) + rrt(066) * density(03) 
  pd(19,03) = pd(19,03) + rrt(066) * density(01) 
  pd(24,01) = pd(24,01) + rrt(066) * density(03) 
  pd(24,03) = pd(24,03) + rrt(066) * density(01) 
  pd(17,01) = pd(17,01) - rrt(067) * density(17) 
  pd(17,17) = pd(17,17) - rrt(067) * density(01) 
  pd(46,01) = pd(46,01) + rrt(067) * density(17) 
  pd(46,17) = pd(46,17) + rrt(067) * density(01) 
  pd(50,01) = pd(50,01) + rrt(067) * density(17) 
  pd(50,17) = pd(50,17) + rrt(067) * density(01) 
  pd(10,01) = pd(10,01) + rrt(068) * density(17) 
  pd(10,17) = pd(10,17) + rrt(068) * density(01) 
  pd(17,01) = pd(17,01) - rrt(068) * density(17) 
  pd(17,17) = pd(17,17) - rrt(068) * density(01) 
  pd(19,01) = pd(19,01) + rrt(068) * density(17) 
  pd(19,17) = pd(19,17) + rrt(068) * density(01) 
  pd(01,01) = pd(01,01) + rrt(069) * density(20) 
  pd(01,20) = pd(01,20) + rrt(069) * density(01) 
  pd(15,01) = pd(15,01) + rrt(069) * density(20) 
  pd(15,20) = pd(15,20) + rrt(069) * density(01) 
  pd(20,01) = pd(20,01) - rrt(069) * density(20) 
  pd(20,20) = pd(20,20) - rrt(069) * density(01) 
  pd(01,01) = pd(01,01) + rrt(070) * density(20) 
  pd(01,20) = pd(01,20) + rrt(070) * density(01) 
  pd(13,01) = pd(13,01) + rrt(070) * density(20) 
  pd(13,20) = pd(13,20) + rrt(070) * density(01) 
  pd(20,01) = pd(20,01) - rrt(070) * density(20) 
  pd(20,20) = pd(20,20) - rrt(070) * density(01) 
  pd(44,01) = pd(44,01) + rrt(070) * density(20) 
  pd(44,20) = pd(44,20) + rrt(070) * density(01) 
  pd(01,01) = pd(01,01) + rrt(071) * density(20) 
  pd(01,20) = pd(01,20) + rrt(071) * density(01) 
  pd(02,01) = pd(02,01) + rrt(071) * density(20) 
  pd(02,20) = pd(02,20) + rrt(071) * density(01) 
  pd(20,01) = pd(20,01) - rrt(071) * density(20) 
  pd(20,20) = pd(20,20) - rrt(071) * density(01) 
  pd(32,01) = pd(32,01) + rrt(071) * density(20) 
  pd(32,20) = pd(32,20) + rrt(071) * density(01) 
  pd(01,01) = pd(01,01) + rrt(072) * density(20) 
  pd(01,20) = pd(01,20) + rrt(072) * density(01) 
  pd(02,01) = pd(02,01) + rrt(072) * density(20) 
  pd(02,20) = pd(02,20) + rrt(072) * density(01) 
  pd(06,01) = pd(06,01) + rrt(072) * density(20) 
  pd(06,20) = pd(06,20) + rrt(072) * density(01) 
  pd(13,01) = pd(13,01) + rrt(072) * density(20) 
  pd(13,20) = pd(13,20) + rrt(072) * density(01) 
  pd(20,01) = pd(20,01) - rrt(072) * density(20) 
  pd(20,20) = pd(20,20) - rrt(072) * density(01) 
  pd(01,01) = pd(01,01) + rrt(073) * density(20) 
  pd(01,20) = pd(01,20) + rrt(073) * density(01) 
  pd(02,01) = pd(02,01) + rrt(073) * density(20) * 2.0d0
  pd(02,20) = pd(02,20) + rrt(073) * density(01) * 2.0d0
  pd(14,01) = pd(14,01) + rrt(073) * density(20) 
  pd(14,20) = pd(14,20) + rrt(073) * density(01) 
  pd(20,01) = pd(20,01) - rrt(073) * density(20) 
  pd(20,20) = pd(20,20) - rrt(073) * density(01) 
  pd(01,01) = pd(01,01) + rrt(074) * density(20) 
  pd(01,20) = pd(01,20) + rrt(074) * density(01) 
  pd(20,01) = pd(20,01) - rrt(074) * density(20) 
  pd(20,20) = pd(20,20) - rrt(074) * density(01) 
  pd(24,01) = pd(24,01) + rrt(074) * density(20) 
  pd(24,20) = pd(24,20) + rrt(074) * density(01) 
  pd(48,01) = pd(48,01) + rrt(074) * density(20) 
  pd(48,20) = pd(48,20) + rrt(074) * density(01) 
  pd(01,01) = pd(01,01) + rrt(075) * density(20) 
  pd(01,20) = pd(01,20) + rrt(075) * density(01) 
  pd(05,01) = pd(05,01) + rrt(075) * density(20) 
  pd(05,20) = pd(05,20) + rrt(075) * density(01) 
  pd(20,01) = pd(20,01) - rrt(075) * density(20) 
  pd(20,20) = pd(20,20) - rrt(075) * density(01) 
  pd(27,01) = pd(27,01) + rrt(075) * density(20) 
  pd(27,20) = pd(27,20) + rrt(075) * density(01) 
  pd(01,01) = pd(01,01) + rrt(076) * density(38) 
  pd(01,38) = pd(01,38) + rrt(076) * density(01) 
  pd(38,01) = pd(38,01) - rrt(076) * density(38) 
  pd(38,38) = pd(38,38) - rrt(076) * density(01) 
  pd(44,01) = pd(44,01) + rrt(076) * density(38) 
  pd(44,38) = pd(44,38) + rrt(076) * density(01) 
  pd(01,01) = pd(01,01) + rrt(077) * density(38) 
  pd(01,38) = pd(01,38) + rrt(077) * density(01) 
  pd(13,01) = pd(13,01) + rrt(077) * density(38) 
  pd(13,38) = pd(13,38) + rrt(077) * density(01) 
  pd(32,01) = pd(32,01) + rrt(077) * density(38) 
  pd(32,38) = pd(32,38) + rrt(077) * density(01) 
  pd(38,01) = pd(38,01) - rrt(077) * density(38) 
  pd(38,38) = pd(38,38) - rrt(077) * density(01) 
  pd(01,01) = pd(01,01) + rrt(078) * density(38) 
  pd(01,38) = pd(01,38) + rrt(078) * density(01) 
  pd(02,01) = pd(02,01) + rrt(078) * density(38) 
  pd(02,38) = pd(02,38) + rrt(078) * density(01) 
  pd(06,01) = pd(06,01) + rrt(078) * density(38) 
  pd(06,38) = pd(06,38) + rrt(078) * density(01) 
  pd(38,01) = pd(38,01) - rrt(078) * density(38) 
  pd(38,38) = pd(38,38) - rrt(078) * density(01) 
  pd(01,01) = pd(01,01) + rrt(079) * density(38) 
  pd(01,38) = pd(01,38) + rrt(079) * density(01) 
  pd(02,01) = pd(02,01) + rrt(079) * density(38) 
  pd(02,38) = pd(02,38) + rrt(079) * density(01) 
  pd(13,01) = pd(13,01) + rrt(079) * density(38) 
  pd(13,38) = pd(13,38) + rrt(079) * density(01) 
  pd(14,01) = pd(14,01) + rrt(079) * density(38) 
  pd(14,38) = pd(14,38) + rrt(079) * density(01) 
  pd(38,01) = pd(38,01) - rrt(079) * density(38) 
  pd(38,38) = pd(38,38) - rrt(079) * density(01) 
  pd(01,01) = pd(01,01) + rrt(080) * density(38) 
  pd(01,38) = pd(01,38) + rrt(080) * density(01) 
  pd(10,01) = pd(10,01) + rrt(080) * density(38) 
  pd(10,38) = pd(10,38) + rrt(080) * density(01) 
  pd(38,01) = pd(38,01) - rrt(080) * density(38) 
  pd(38,38) = pd(38,38) - rrt(080) * density(01) 
  pd(48,01) = pd(48,01) + rrt(080) * density(38) 
  pd(48,38) = pd(48,38) + rrt(080) * density(01) 
  pd(01,01) = pd(01,01) + rrt(081) * density(38) 
  pd(01,38) = pd(01,38) + rrt(081) * density(01) 
  pd(24,01) = pd(24,01) + rrt(081) * density(38) 
  pd(24,38) = pd(24,38) + rrt(081) * density(01) 
  pd(27,01) = pd(27,01) + rrt(081) * density(38) 
  pd(27,38) = pd(27,38) + rrt(081) * density(01) 
  pd(38,01) = pd(38,01) - rrt(081) * density(38) 
  pd(38,38) = pd(38,38) - rrt(081) * density(01) 
  pd(01,01) = pd(01,01) + rrt(082) * density(38) 
  pd(01,38) = pd(01,38) + rrt(082) * density(01) 
  pd(05,01) = pd(05,01) + rrt(082) * density(38) 
  pd(05,38) = pd(05,38) + rrt(082) * density(01) 
  pd(28,01) = pd(28,01) + rrt(082) * density(38) 
  pd(28,38) = pd(28,38) + rrt(082) * density(01) 
  pd(38,01) = pd(38,01) - rrt(082) * density(38) 
  pd(38,38) = pd(38,38) - rrt(082) * density(01) 
  pd(01,01) = pd(01,01) + rrt(083) * density(35) 
  pd(01,35) = pd(01,35) + rrt(083) * density(01) 
  pd(32,01) = pd(32,01) + rrt(083) * density(35) 
  pd(32,35) = pd(32,35) + rrt(083) * density(01) 
  pd(35,01) = pd(35,01) - rrt(083) * density(35) 
  pd(35,35) = pd(35,35) - rrt(083) * density(01) 
  pd(01,01) = pd(01,01) + rrt(084) * density(35) 
  pd(01,35) = pd(01,35) + rrt(084) * density(01) 
  pd(06,01) = pd(06,01) + rrt(084) * density(35) 
  pd(06,35) = pd(06,35) + rrt(084) * density(01) 
  pd(13,01) = pd(13,01) + rrt(084) * density(35) 
  pd(13,35) = pd(13,35) + rrt(084) * density(01) 
  pd(35,01) = pd(35,01) - rrt(084) * density(35) 
  pd(35,35) = pd(35,35) - rrt(084) * density(01) 
  pd(01,01) = pd(01,01) + rrt(085) * density(35) 
  pd(01,35) = pd(01,35) + rrt(085) * density(01) 
  pd(35,01) = pd(35,01) - rrt(085) * density(35) 
  pd(35,35) = pd(35,35) - rrt(085) * density(01) 
  pd(46,01) = pd(46,01) + rrt(085) * density(35) 
  pd(46,35) = pd(46,35) + rrt(085) * density(01) 
  pd(48,01) = pd(48,01) + rrt(085) * density(35) 
  pd(48,35) = pd(48,35) + rrt(085) * density(01) 
  pd(01,01) = pd(01,01) + rrt(086) * density(35) 
  pd(01,35) = pd(01,35) + rrt(086) * density(01) 
  pd(10,01) = pd(10,01) + rrt(086) * density(35) 
  pd(10,35) = pd(10,35) + rrt(086) * density(01) 
  pd(27,01) = pd(27,01) + rrt(086) * density(35) 
  pd(27,35) = pd(27,35) + rrt(086) * density(01) 
  pd(35,01) = pd(35,01) - rrt(086) * density(35) 
  pd(35,35) = pd(35,35) - rrt(086) * density(01) 
  pd(01,01) = pd(01,01) + rrt(087) * density(35) 
  pd(01,35) = pd(01,35) + rrt(087) * density(01) 
  pd(24,01) = pd(24,01) + rrt(087) * density(35) 
  pd(24,35) = pd(24,35) + rrt(087) * density(01) 
  pd(28,01) = pd(28,01) + rrt(087) * density(35) 
  pd(28,35) = pd(28,35) + rrt(087) * density(01) 
  pd(35,01) = pd(35,01) - rrt(087) * density(35) 
  pd(35,35) = pd(35,35) - rrt(087) * density(01) 
  pd(01,01) = pd(01,01) + rrt(088) * density(50) 
  pd(01,50) = pd(01,50) + rrt(088) * density(01) 
  pd(06,01) = pd(06,01) + rrt(088) * density(50) 
  pd(06,50) = pd(06,50) + rrt(088) * density(01) 
  pd(50,01) = pd(50,01) - rrt(088) * density(50) 
  pd(50,50) = pd(50,50) - rrt(088) * density(01) 
  pd(01,01) = pd(01,01) + rrt(089) * density(50) 
  pd(01,50) = pd(01,50) + rrt(089) * density(01) 
  pd(13,01) = pd(13,01) + rrt(089) * density(50) 
  pd(13,50) = pd(13,50) + rrt(089) * density(01) 
  pd(14,01) = pd(14,01) + rrt(089) * density(50) 
  pd(14,50) = pd(14,50) + rrt(089) * density(01) 
  pd(50,01) = pd(50,01) - rrt(089) * density(50) 
  pd(50,50) = pd(50,50) - rrt(089) * density(01) 
  pd(01,01) = pd(01,01) + rrt(090) * density(50) 
  pd(01,50) = pd(01,50) + rrt(090) * density(01) 
  pd(27,01) = pd(27,01) + rrt(090) * density(50) 
  pd(27,50) = pd(27,50) + rrt(090) * density(01) 
  pd(46,01) = pd(46,01) + rrt(090) * density(50) 
  pd(46,50) = pd(46,50) + rrt(090) * density(01) 
  pd(50,01) = pd(50,01) - rrt(090) * density(50) 
  pd(50,50) = pd(50,50) - rrt(090) * density(01) 
  pd(01,01) = pd(01,01) + rrt(091) * density(50) 
  pd(01,50) = pd(01,50) + rrt(091) * density(01) 
  pd(10,01) = pd(10,01) + rrt(091) * density(50) 
  pd(10,50) = pd(10,50) + rrt(091) * density(01) 
  pd(28,01) = pd(28,01) + rrt(091) * density(50) 
  pd(28,50) = pd(28,50) + rrt(091) * density(01) 
  pd(50,01) = pd(50,01) - rrt(091) * density(50) 
  pd(50,50) = pd(50,50) - rrt(091) * density(01) 
  pd(01,01) = pd(01,01) + rrt(092) * density(50) 
  pd(01,50) = pd(01,50) + rrt(092) * density(01) 
  pd(07,01) = pd(07,01) + rrt(092) * density(50) 
  pd(07,50) = pd(07,50) + rrt(092) * density(01) 
  pd(19,01) = pd(19,01) + rrt(092) * density(50) 
  pd(19,50) = pd(19,50) + rrt(092) * density(01) 
  pd(50,01) = pd(50,01) - rrt(092) * density(50) 
  pd(50,50) = pd(50,50) - rrt(092) * density(01) 
  pd(01,01) = pd(01,01) + rrt(093) * density(19) 
  pd(01,19) = pd(01,19) + rrt(093) * density(01) 
  pd(14,01) = pd(14,01) + rrt(093) * density(19) 
  pd(14,19) = pd(14,19) + rrt(093) * density(01) 
  pd(19,01) = pd(19,01) - rrt(093) * density(19) 
  pd(19,19) = pd(19,19) - rrt(093) * density(01) 
  pd(01,01) = pd(01,01) + rrt(094) * density(19) 
  pd(01,19) = pd(01,19) + rrt(094) * density(01) 
  pd(19,01) = pd(19,01) - rrt(094) * density(19) 
  pd(19,19) = pd(19,19) - rrt(094) * density(01) 
  pd(28,01) = pd(28,01) + rrt(094) * density(19) 
  pd(28,19) = pd(28,19) + rrt(094) * density(01) 
  pd(46,01) = pd(46,01) + rrt(094) * density(19) 
  pd(46,19) = pd(46,19) + rrt(094) * density(01) 
  pd(01,01) = pd(01,01) + rrt(095) * density(12) 
  pd(01,12) = pd(01,12) + rrt(095) * density(01) 
  pd(12,01) = pd(12,01) - rrt(095) * density(12) 
  pd(12,12) = pd(12,12) - rrt(095) * density(01) 
  pd(22,01) = pd(22,01) + rrt(095) * density(12) 
  pd(22,12) = pd(22,12) + rrt(095) * density(01) 
  pd(01,01) = pd(01,01) + rrt(096) * density(12) 
  pd(01,12) = pd(01,12) + rrt(096) * density(01) 
  pd(12,01) = pd(12,01) - rrt(096) * density(12) 
  pd(12,12) = pd(12,12) - rrt(096) * density(01) 
  pd(13,01) = pd(13,01) + rrt(096) * density(12) 
  pd(13,12) = pd(13,12) + rrt(096) * density(01) 
  pd(30,01) = pd(30,01) + rrt(096) * density(12) 
  pd(30,12) = pd(30,12) + rrt(096) * density(01) 
  pd(01,01) = pd(01,01) + rrt(097) * density(12) 
  pd(01,12) = pd(01,12) + rrt(097) * density(01) 
  pd(02,01) = pd(02,01) + rrt(097) * density(12) 
  pd(02,12) = pd(02,12) + rrt(097) * density(01) 
  pd(12,01) = pd(12,01) - rrt(097) * density(12) 
  pd(12,12) = pd(12,12) - rrt(097) * density(01) 
  pd(49,01) = pd(49,01) + rrt(097) * density(12) 
  pd(49,12) = pd(49,12) + rrt(097) * density(01) 
  pd(01,01) = pd(01,01) + rrt(098) * density(12) 
  pd(01,12) = pd(01,12) + rrt(098) * density(01) 
  pd(02,01) = pd(02,01) + rrt(098) * density(12) 
  pd(02,12) = pd(02,12) + rrt(098) * density(01) 
  pd(08,01) = pd(08,01) + rrt(098) * density(12) 
  pd(08,12) = pd(08,12) + rrt(098) * density(01) 
  pd(12,01) = pd(12,01) - rrt(098) * density(12) 
  pd(12,12) = pd(12,12) - rrt(098) * density(01) 
  pd(13,01) = pd(13,01) + rrt(098) * density(12) 
  pd(13,12) = pd(13,12) + rrt(098) * density(01) 
  pd(01,01) = pd(01,01) + rrt(099) * density(12) 
  pd(01,12) = pd(01,12) + rrt(099) * density(01) 
  pd(02,01) = pd(02,01) + rrt(099) * density(12) * 2.0d0
  pd(02,12) = pd(02,12) + rrt(099) * density(01) * 2.0d0
  pd(12,01) = pd(12,01) - rrt(099) * density(12) 
  pd(12,12) = pd(12,12) - rrt(099) * density(01) 
  pd(47,01) = pd(47,01) + rrt(099) * density(12) 
  pd(47,12) = pd(47,12) + rrt(099) * density(01) 
  pd(01,01) = pd(01,01) + rrt(100) * density(12) 
  pd(01,12) = pd(01,12) + rrt(100) * density(01) 
  pd(12,01) = pd(12,01) - rrt(100) * density(12) 
  pd(12,12) = pd(12,12) - rrt(100) * density(01) 
  pd(24,01) = pd(24,01) + rrt(100) * density(12) 
  pd(24,12) = pd(24,12) + rrt(100) * density(01) 
  pd(44,01) = pd(44,01) + rrt(100) * density(12) 
  pd(44,12) = pd(44,12) + rrt(100) * density(01) 
  pd(01,01) = pd(01,01) + rrt(101) * density(12) 
  pd(01,12) = pd(01,12) + rrt(101) * density(01) 
  pd(05,01) = pd(05,01) + rrt(101) * density(12) 
  pd(05,12) = pd(05,12) + rrt(101) * density(01) 
  pd(12,01) = pd(12,01) - rrt(101) * density(12) 
  pd(12,12) = pd(12,12) - rrt(101) * density(01) 
  pd(32,01) = pd(32,01) + rrt(101) * density(12) 
  pd(32,12) = pd(32,12) + rrt(101) * density(01) 
  pd(01,01) = pd(01,01) + rrt(102) * density(12) 
  pd(01,12) = pd(01,12) + rrt(102) * density(01) 
  pd(12,01) = pd(12,01) - rrt(102) * density(12) 
  pd(12,12) = pd(12,12) - rrt(102) * density(01) 
  pd(38,01) = pd(38,01) + rrt(102) * density(12) 
  pd(38,12) = pd(38,12) + rrt(102) * density(01) 
  pd(48,01) = pd(48,01) + rrt(102) * density(12) 
  pd(48,12) = pd(48,12) + rrt(102) * density(01) 
  pd(01,01) = pd(01,01) + rrt(103) * density(12) 
  pd(01,12) = pd(01,12) + rrt(103) * density(01) 
  pd(12,01) = pd(12,01) - rrt(103) * density(12) 
  pd(12,12) = pd(12,12) - rrt(103) * density(01) 
  pd(20,01) = pd(20,01) + rrt(103) * density(12) 
  pd(20,12) = pd(20,12) + rrt(103) * density(01) 
  pd(27,01) = pd(27,01) + rrt(103) * density(12) 
  pd(27,12) = pd(27,12) + rrt(103) * density(01) 
  pd(01,01) = pd(01,01) + rrt(104) * density(23) 
  pd(01,23) = pd(01,23) + rrt(104) * density(01) 
  pd(23,01) = pd(23,01) - rrt(104) * density(23) 
  pd(23,23) = pd(23,23) - rrt(104) * density(01) 
  pd(30,01) = pd(30,01) + rrt(104) * density(23) 
  pd(30,23) = pd(30,23) + rrt(104) * density(01) 
  pd(01,01) = pd(01,01) + rrt(105) * density(23) 
  pd(01,23) = pd(01,23) + rrt(105) * density(01) 
  pd(13,01) = pd(13,01) + rrt(105) * density(23) 
  pd(13,23) = pd(13,23) + rrt(105) * density(01) 
  pd(23,01) = pd(23,01) - rrt(105) * density(23) 
  pd(23,23) = pd(23,23) - rrt(105) * density(01) 
  pd(49,01) = pd(49,01) + rrt(105) * density(23) 
  pd(49,23) = pd(49,23) + rrt(105) * density(01) 
  pd(01,01) = pd(01,01) + rrt(106) * density(23) 
  pd(01,23) = pd(01,23) + rrt(106) * density(01) 
  pd(02,01) = pd(02,01) + rrt(106) * density(23) 
  pd(02,23) = pd(02,23) + rrt(106) * density(01) 
  pd(08,01) = pd(08,01) + rrt(106) * density(23) 
  pd(08,23) = pd(08,23) + rrt(106) * density(01) 
  pd(23,01) = pd(23,01) - rrt(106) * density(23) 
  pd(23,23) = pd(23,23) - rrt(106) * density(01) 
  pd(01,01) = pd(01,01) + rrt(107) * density(23) 
  pd(01,23) = pd(01,23) + rrt(107) * density(01) 
  pd(02,01) = pd(02,01) + rrt(107) * density(23) 
  pd(02,23) = pd(02,23) + rrt(107) * density(01) 
  pd(13,01) = pd(13,01) + rrt(107) * density(23) 
  pd(13,23) = pd(13,23) + rrt(107) * density(01) 
  pd(23,01) = pd(23,01) - rrt(107) * density(23) 
  pd(23,23) = pd(23,23) - rrt(107) * density(01) 
  pd(47,01) = pd(47,01) + rrt(107) * density(23) 
  pd(47,23) = pd(47,23) + rrt(107) * density(01) 
  pd(01,01) = pd(01,01) + rrt(108) * density(23) 
  pd(01,23) = pd(01,23) + rrt(108) * density(01) 
  pd(10,01) = pd(10,01) + rrt(108) * density(23) 
  pd(10,23) = pd(10,23) + rrt(108) * density(01) 
  pd(23,01) = pd(23,01) - rrt(108) * density(23) 
  pd(23,23) = pd(23,23) - rrt(108) * density(01) 
  pd(44,01) = pd(44,01) + rrt(108) * density(23) 
  pd(44,23) = pd(44,23) + rrt(108) * density(01) 
  pd(01,01) = pd(01,01) + rrt(109) * density(23) 
  pd(01,23) = pd(01,23) + rrt(109) * density(01) 
  pd(23,01) = pd(23,01) - rrt(109) * density(23) 
  pd(23,23) = pd(23,23) - rrt(109) * density(01) 
  pd(24,01) = pd(24,01) + rrt(109) * density(23) 
  pd(24,23) = pd(24,23) + rrt(109) * density(01) 
  pd(32,01) = pd(32,01) + rrt(109) * density(23) 
  pd(32,23) = pd(32,23) + rrt(109) * density(01) 
  pd(01,01) = pd(01,01) + rrt(110) * density(23) 
  pd(01,23) = pd(01,23) + rrt(110) * density(01) 
  pd(05,01) = pd(05,01) + rrt(110) * density(23) 
  pd(05,23) = pd(05,23) + rrt(110) * density(01) 
  pd(06,01) = pd(06,01) + rrt(110) * density(23) 
  pd(06,23) = pd(06,23) + rrt(110) * density(01) 
  pd(23,01) = pd(23,01) - rrt(110) * density(23) 
  pd(23,23) = pd(23,23) - rrt(110) * density(01) 
  pd(01,01) = pd(01,01) + rrt(111) * density(23) 
  pd(01,23) = pd(01,23) + rrt(111) * density(01) 
  pd(23,01) = pd(23,01) - rrt(111) * density(23) 
  pd(23,23) = pd(23,23) - rrt(111) * density(01) 
  pd(25,01) = pd(25,01) + rrt(111) * density(23) 
  pd(25,23) = pd(25,23) + rrt(111) * density(01) 
  pd(50,01) = pd(50,01) + rrt(111) * density(23) 
  pd(50,23) = pd(50,23) + rrt(111) * density(01) 
  pd(01,01) = pd(01,01) + rrt(112) * density(23) 
  pd(01,23) = pd(01,23) + rrt(112) * density(01) 
  pd(23,01) = pd(23,01) - rrt(112) * density(23) 
  pd(23,23) = pd(23,23) - rrt(112) * density(01) 
  pd(35,01) = pd(35,01) + rrt(112) * density(23) 
  pd(35,23) = pd(35,23) + rrt(112) * density(01) 
  pd(48,01) = pd(48,01) + rrt(112) * density(23) 
  pd(48,23) = pd(48,23) + rrt(112) * density(01) 
  pd(01,01) = pd(01,01) + rrt(113) * density(23) 
  pd(01,23) = pd(01,23) + rrt(113) * density(01) 
  pd(23,01) = pd(23,01) - rrt(113) * density(23) 
  pd(23,23) = pd(23,23) - rrt(113) * density(01) 
  pd(27,01) = pd(27,01) + rrt(113) * density(23) 
  pd(27,23) = pd(27,23) + rrt(113) * density(01) 
  pd(38,01) = pd(38,01) + rrt(113) * density(23) 
  pd(38,23) = pd(38,23) + rrt(113) * density(01) 
  pd(01,01) = pd(01,01) + rrt(114) * density(23) 
  pd(01,23) = pd(01,23) + rrt(114) * density(01) 
  pd(20,01) = pd(20,01) + rrt(114) * density(23) 
  pd(20,23) = pd(20,23) + rrt(114) * density(01) 
  pd(23,01) = pd(23,01) - rrt(114) * density(23) 
  pd(23,23) = pd(23,23) - rrt(114) * density(01) 
  pd(28,01) = pd(28,01) + rrt(114) * density(23) 
  pd(28,23) = pd(28,23) + rrt(114) * density(01) 
  pd(01,01) = pd(01,01) + rrt(115) * density(09) 
  pd(01,09) = pd(01,09) + rrt(115) * density(01) 
  pd(09,01) = pd(09,01) - rrt(115) * density(09) 
  pd(09,09) = pd(09,09) - rrt(115) * density(01) 
  pd(49,01) = pd(49,01) + rrt(115) * density(09) 
  pd(49,09) = pd(49,09) + rrt(115) * density(01) 
  pd(01,01) = pd(01,01) + rrt(116) * density(09) 
  pd(01,09) = pd(01,09) + rrt(116) * density(01) 
  pd(08,01) = pd(08,01) + rrt(116) * density(09) 
  pd(08,09) = pd(08,09) + rrt(116) * density(01) 
  pd(09,01) = pd(09,01) - rrt(116) * density(09) 
  pd(09,09) = pd(09,09) - rrt(116) * density(01) 
  pd(13,01) = pd(13,01) + rrt(116) * density(09) 
  pd(13,09) = pd(13,09) + rrt(116) * density(01) 
  pd(01,01) = pd(01,01) + rrt(117) * density(09) 
  pd(01,09) = pd(01,09) + rrt(117) * density(01) 
  pd(02,01) = pd(02,01) + rrt(117) * density(09) 
  pd(02,09) = pd(02,09) + rrt(117) * density(01) 
  pd(09,01) = pd(09,01) - rrt(117) * density(09) 
  pd(09,09) = pd(09,09) - rrt(117) * density(01) 
  pd(47,01) = pd(47,01) + rrt(117) * density(09) 
  pd(47,09) = pd(47,09) + rrt(117) * density(01) 
  pd(01,01) = pd(01,01) + rrt(118) * density(09) 
  pd(01,09) = pd(01,09) + rrt(118) * density(01) 
  pd(09,01) = pd(09,01) - rrt(118) * density(09) 
  pd(09,09) = pd(09,09) - rrt(118) * density(01) 
  pd(44,01) = pd(44,01) + rrt(118) * density(09) 
  pd(44,09) = pd(44,09) + rrt(118) * density(01) 
  pd(46,01) = pd(46,01) + rrt(118) * density(09) 
  pd(46,09) = pd(46,09) + rrt(118) * density(01) 
  pd(01,01) = pd(01,01) + rrt(119) * density(09) 
  pd(01,09) = pd(01,09) + rrt(119) * density(01) 
  pd(09,01) = pd(09,01) - rrt(119) * density(09) 
  pd(09,09) = pd(09,09) - rrt(119) * density(01) 
  pd(10,01) = pd(10,01) + rrt(119) * density(09) 
  pd(10,09) = pd(10,09) + rrt(119) * density(01) 
  pd(32,01) = pd(32,01) + rrt(119) * density(09) 
  pd(32,09) = pd(32,09) + rrt(119) * density(01) 
  pd(01,01) = pd(01,01) + rrt(120) * density(09) 
  pd(01,09) = pd(01,09) + rrt(120) * density(01) 
  pd(06,01) = pd(06,01) + rrt(120) * density(09) 
  pd(06,09) = pd(06,09) + rrt(120) * density(01) 
  pd(09,01) = pd(09,01) - rrt(120) * density(09) 
  pd(09,09) = pd(09,09) - rrt(120) * density(01) 
  pd(24,01) = pd(24,01) + rrt(120) * density(09) 
  pd(24,09) = pd(24,09) + rrt(120) * density(01) 
  pd(01,01) = pd(01,01) + rrt(121) * density(09) 
  pd(01,09) = pd(01,09) + rrt(121) * density(01) 
  pd(05,01) = pd(05,01) + rrt(121) * density(09) 
  pd(05,09) = pd(05,09) + rrt(121) * density(01) 
  pd(09,01) = pd(09,01) - rrt(121) * density(09) 
  pd(09,09) = pd(09,09) - rrt(121) * density(01) 
  pd(14,01) = pd(14,01) + rrt(121) * density(09) 
  pd(14,09) = pd(14,09) + rrt(121) * density(01) 
  pd(01,01) = pd(01,01) + rrt(122) * density(09) 
  pd(01,09) = pd(01,09) + rrt(122) * density(01) 
  pd(09,01) = pd(09,01) - rrt(122) * density(09) 
  pd(09,09) = pd(09,09) - rrt(122) * density(01) 
  pd(19,01) = pd(19,01) + rrt(122) * density(09) 
  pd(19,09) = pd(19,09) + rrt(122) * density(01) 
  pd(25,01) = pd(25,01) + rrt(122) * density(09) 
  pd(25,09) = pd(25,09) + rrt(122) * density(01) 
  pd(01,01) = pd(01,01) + rrt(123) * density(09) 
  pd(01,09) = pd(01,09) + rrt(123) * density(01) 
  pd(09,01) = pd(09,01) - rrt(123) * density(09) 
  pd(09,09) = pd(09,09) - rrt(123) * density(01) 
  pd(48,01) = pd(48,01) + rrt(123) * density(09) 
  pd(48,09) = pd(48,09) + rrt(123) * density(01) 
  pd(50,01) = pd(50,01) + rrt(123) * density(09) 
  pd(50,09) = pd(50,09) + rrt(123) * density(01) 
  pd(01,01) = pd(01,01) + rrt(124) * density(09) 
  pd(01,09) = pd(01,09) + rrt(124) * density(01) 
  pd(09,01) = pd(09,01) - rrt(124) * density(09) 
  pd(09,09) = pd(09,09) - rrt(124) * density(01) 
  pd(27,01) = pd(27,01) + rrt(124) * density(09) 
  pd(27,09) = pd(27,09) + rrt(124) * density(01) 
  pd(35,01) = pd(35,01) + rrt(124) * density(09) 
  pd(35,09) = pd(35,09) + rrt(124) * density(01) 
  pd(01,01) = pd(01,01) + rrt(125) * density(09) 
  pd(01,09) = pd(01,09) + rrt(125) * density(01) 
  pd(09,01) = pd(09,01) - rrt(125) * density(09) 
  pd(09,09) = pd(09,09) - rrt(125) * density(01) 
  pd(28,01) = pd(28,01) + rrt(125) * density(09) 
  pd(28,09) = pd(28,09) + rrt(125) * density(01) 
  pd(38,01) = pd(38,01) + rrt(125) * density(09) 
  pd(38,09) = pd(38,09) + rrt(125) * density(01) 
  pd(01,01) = pd(01,01) + rrt(126) * density(03) 
  pd(01,03) = pd(01,03) + rrt(126) * density(01) 
  pd(03,01) = pd(03,01) - rrt(126) * density(03) 
  pd(03,03) = pd(03,03) - rrt(126) * density(01) 
  pd(08,01) = pd(08,01) + rrt(126) * density(03) 
  pd(08,03) = pd(08,03) + rrt(126) * density(01) 
  pd(01,01) = pd(01,01) + rrt(127) * density(03) 
  pd(01,03) = pd(01,03) + rrt(127) * density(01) 
  pd(03,01) = pd(03,01) - rrt(127) * density(03) 
  pd(03,03) = pd(03,03) - rrt(127) * density(01) 
  pd(13,01) = pd(13,01) + rrt(127) * density(03) 
  pd(13,03) = pd(13,03) + rrt(127) * density(01) 
  pd(47,01) = pd(47,01) + rrt(127) * density(03) 
  pd(47,03) = pd(47,03) + rrt(127) * density(01) 
  pd(01,01) = pd(01,01) + rrt(128) * density(03) 
  pd(01,03) = pd(01,03) + rrt(128) * density(01) 
  pd(03,01) = pd(03,01) - rrt(128) * density(03) 
  pd(03,03) = pd(03,03) - rrt(128) * density(01) 
  pd(32,01) = pd(32,01) + rrt(128) * density(03) 
  pd(32,03) = pd(32,03) + rrt(128) * density(01) 
  pd(46,01) = pd(46,01) + rrt(128) * density(03) 
  pd(46,03) = pd(46,03) + rrt(128) * density(01) 
  pd(01,01) = pd(01,01) + rrt(129) * density(03) 
  pd(01,03) = pd(01,03) + rrt(129) * density(01) 
  pd(03,01) = pd(03,01) - rrt(129) * density(03) 
  pd(03,03) = pd(03,03) - rrt(129) * density(01) 
  pd(06,01) = pd(06,01) + rrt(129) * density(03) 
  pd(06,03) = pd(06,03) + rrt(129) * density(01) 
  pd(10,01) = pd(10,01) + rrt(129) * density(03) 
  pd(10,03) = pd(10,03) + rrt(129) * density(01) 
  pd(01,01) = pd(01,01) + rrt(130) * density(03) 
  pd(01,03) = pd(01,03) + rrt(130) * density(01) 
  pd(03,01) = pd(03,01) - rrt(130) * density(03) 
  pd(03,03) = pd(03,03) - rrt(130) * density(01) 
  pd(14,01) = pd(14,01) + rrt(130) * density(03) 
  pd(14,03) = pd(14,03) + rrt(130) * density(01) 
  pd(24,01) = pd(24,01) + rrt(130) * density(03) 
  pd(24,03) = pd(24,03) + rrt(130) * density(01) 
  pd(01,01) = pd(01,01) + rrt(131) * density(03) 
  pd(01,03) = pd(01,03) + rrt(131) * density(01) 
  pd(03,01) = pd(03,01) - rrt(131) * density(03) 
  pd(03,03) = pd(03,03) - rrt(131) * density(01) 
  pd(19,01) = pd(19,01) + rrt(131) * density(03) 
  pd(19,03) = pd(19,03) + rrt(131) * density(01) 
  pd(48,01) = pd(48,01) + rrt(131) * density(03) 
  pd(48,03) = pd(48,03) + rrt(131) * density(01) 
  pd(01,01) = pd(01,01) + rrt(132) * density(03) 
  pd(01,03) = pd(01,03) + rrt(132) * density(01) 
  pd(03,01) = pd(03,01) - rrt(132) * density(03) 
  pd(03,03) = pd(03,03) - rrt(132) * density(01) 
  pd(27,01) = pd(27,01) + rrt(132) * density(03) 
  pd(27,03) = pd(27,03) + rrt(132) * density(01) 
  pd(50,01) = pd(50,01) + rrt(132) * density(03) 
  pd(50,03) = pd(50,03) + rrt(132) * density(01) 
  pd(01,01) = pd(01,01) + rrt(133) * density(03) 
  pd(01,03) = pd(01,03) + rrt(133) * density(01) 
  pd(03,01) = pd(03,01) - rrt(133) * density(03) 
  pd(03,03) = pd(03,03) - rrt(133) * density(01) 
  pd(28,01) = pd(28,01) + rrt(133) * density(03) 
  pd(28,03) = pd(28,03) + rrt(133) * density(01) 
  pd(35,01) = pd(35,01) + rrt(133) * density(03) 
  pd(35,03) = pd(35,03) + rrt(133) * density(01) 
  pd(01,01) = pd(01,01) + rrt(134) * density(17) 
  pd(01,17) = pd(01,17) + rrt(134) * density(01) 
  pd(17,01) = pd(17,01) - rrt(134) * density(17) 
  pd(17,17) = pd(17,17) - rrt(134) * density(01) 
  pd(47,01) = pd(47,01) + rrt(134) * density(17) 
  pd(47,17) = pd(47,17) + rrt(134) * density(01) 
  pd(01,01) = pd(01,01) + rrt(135) * density(17) 
  pd(01,17) = pd(01,17) + rrt(135) * density(01) 
  pd(06,01) = pd(06,01) + rrt(135) * density(17) 
  pd(06,17) = pd(06,17) + rrt(135) * density(01) 
  pd(17,01) = pd(17,01) - rrt(135) * density(17) 
  pd(17,17) = pd(17,17) - rrt(135) * density(01) 
  pd(46,01) = pd(46,01) + rrt(135) * density(17) 
  pd(46,17) = pd(46,17) + rrt(135) * density(01) 
  pd(01,01) = pd(01,01) + rrt(136) * density(17) 
  pd(01,17) = pd(01,17) + rrt(136) * density(01) 
  pd(10,01) = pd(10,01) + rrt(136) * density(17) 
  pd(10,17) = pd(10,17) + rrt(136) * density(01) 
  pd(14,01) = pd(14,01) + rrt(136) * density(17) 
  pd(14,17) = pd(14,17) + rrt(136) * density(01) 
  pd(17,01) = pd(17,01) - rrt(136) * density(17) 
  pd(17,17) = pd(17,17) - rrt(136) * density(01) 
  pd(01,01) = pd(01,01) + rrt(137) * density(17) 
  pd(01,17) = pd(01,17) + rrt(137) * density(01) 
  pd(17,01) = pd(17,01) - rrt(137) * density(17) 
  pd(17,17) = pd(17,17) - rrt(137) * density(01) 
  pd(19,01) = pd(19,01) + rrt(137) * density(17) 
  pd(19,17) = pd(19,17) + rrt(137) * density(01) 
  pd(27,01) = pd(27,01) + rrt(137) * density(17) 
  pd(27,17) = pd(27,17) + rrt(137) * density(01) 
  pd(01,01) = pd(01,01) + rrt(138) * density(17) 
  pd(01,17) = pd(01,17) + rrt(138) * density(01) 
  pd(17,01) = pd(17,01) - rrt(138) * density(17) 
  pd(17,17) = pd(17,17) - rrt(138) * density(01) 
  pd(28,01) = pd(28,01) + rrt(138) * density(17) 
  pd(28,17) = pd(28,17) + rrt(138) * density(01) 
  pd(50,01) = pd(50,01) + rrt(138) * density(17) 
  pd(50,17) = pd(50,17) + rrt(138) * density(01) 
  pd(13,01) = pd(13,01) + rrt(139) * density(43) 
  pd(13,43) = pd(13,43) + rrt(139) * density(01) 
  pd(24,01) = pd(24,01) + rrt(139) * density(43) 
  pd(24,43) = pd(24,43) + rrt(139) * density(01) 
  pd(43,01) = pd(43,01) - rrt(139) * density(43) 
  pd(43,43) = pd(43,43) - rrt(139) * density(01) 
  pd(13,01) = pd(13,01) + rrt(140) * density(36) 
  pd(13,36) = pd(13,36) + rrt(140) * density(01) 
  pd(24,01) = pd(24,01) + rrt(140) * density(36) 
  pd(24,36) = pd(24,36) + rrt(140) * density(01) 
  pd(36,01) = pd(36,01) - rrt(140) * density(36) 
  pd(36,36) = pd(36,36) - rrt(140) * density(01) 
  pd(02,01) = pd(02,01) + rrt(141) * density(43) 
  pd(02,43) = pd(02,43) + rrt(141) * density(01) 
  pd(10,01) = pd(10,01) + rrt(141) * density(43) 
  pd(10,43) = pd(10,43) + rrt(141) * density(01) 
  pd(43,01) = pd(43,01) - rrt(141) * density(43) 
  pd(43,43) = pd(43,43) - rrt(141) * density(01) 
  pd(02,01) = pd(02,01) + rrt(142) * density(36) 
  pd(02,36) = pd(02,36) + rrt(142) * density(01) 
  pd(10,01) = pd(10,01) + rrt(142) * density(36) 
  pd(10,36) = pd(10,36) + rrt(142) * density(01) 
  pd(36,01) = pd(36,01) - rrt(142) * density(36) 
  pd(36,36) = pd(36,36) - rrt(142) * density(01) 
  pd(02,01) = pd(02,01) + rrt(143) * density(43) 
  pd(02,43) = pd(02,43) + rrt(143) * density(01) 
  pd(13,01) = pd(13,01) + rrt(143) * density(43) 
  pd(13,43) = pd(13,43) + rrt(143) * density(01) 
  pd(43,01) = pd(43,01) - rrt(143) * density(43) 
  pd(43,43) = pd(43,43) - rrt(143) * density(01) 
  pd(46,01) = pd(46,01) + rrt(143) * density(43) 
  pd(46,43) = pd(46,43) + rrt(143) * density(01) 
  pd(02,01) = pd(02,01) + rrt(144) * density(36) 
  pd(02,36) = pd(02,36) + rrt(144) * density(01) 
  pd(13,01) = pd(13,01) + rrt(144) * density(36) 
  pd(13,36) = pd(13,36) + rrt(144) * density(01) 
  pd(36,01) = pd(36,01) - rrt(144) * density(36) 
  pd(36,36) = pd(36,36) - rrt(144) * density(01) 
  pd(46,01) = pd(46,01) + rrt(144) * density(36) 
  pd(46,36) = pd(46,36) + rrt(144) * density(01) 
  pd(01,01) = pd(01,01) + rrt(145) * density(43) 
  pd(01,43) = pd(01,43) + rrt(145) * density(01) 
  pd(25,01) = pd(25,01) + rrt(145) * density(43) 
  pd(25,43) = pd(25,43) + rrt(145) * density(01) 
  pd(43,01) = pd(43,01) - rrt(145) * density(43) 
  pd(43,43) = pd(43,43) - rrt(145) * density(01) 
  pd(01,01) = pd(01,01) + rrt(146) * density(36) 
  pd(01,36) = pd(01,36) + rrt(146) * density(01) 
  pd(25,01) = pd(25,01) + rrt(146) * density(36) 
  pd(25,36) = pd(25,36) + rrt(146) * density(01) 
  pd(36,01) = pd(36,01) - rrt(146) * density(36) 
  pd(36,36) = pd(36,36) - rrt(146) * density(01) 
  pd(01,01) = pd(01,01) + rrt(147) * density(43) 
  pd(01,43) = pd(01,43) + rrt(147) * density(01) 
  pd(13,01) = pd(13,01) + rrt(147) * density(43) 
  pd(13,43) = pd(13,43) + rrt(147) * density(01) 
  pd(43,01) = pd(43,01) - rrt(147) * density(43) 
  pd(43,43) = pd(43,43) - rrt(147) * density(01) 
  pd(48,01) = pd(48,01) + rrt(147) * density(43) 
  pd(48,43) = pd(48,43) + rrt(147) * density(01) 
  pd(01,01) = pd(01,01) + rrt(148) * density(36) 
  pd(01,36) = pd(01,36) + rrt(148) * density(01) 
  pd(13,01) = pd(13,01) + rrt(148) * density(36) 
  pd(13,36) = pd(13,36) + rrt(148) * density(01) 
  pd(36,01) = pd(36,01) - rrt(148) * density(36) 
  pd(36,36) = pd(36,36) - rrt(148) * density(01) 
  pd(48,01) = pd(48,01) + rrt(148) * density(36) 
  pd(48,36) = pd(48,36) + rrt(148) * density(01) 
  pd(01,01) = pd(01,01) + rrt(149) * density(43) 
  pd(01,43) = pd(01,43) + rrt(149) * density(01) 
  pd(02,01) = pd(02,01) + rrt(149) * density(43) 
  pd(02,43) = pd(02,43) + rrt(149) * density(01) 
  pd(27,01) = pd(27,01) + rrt(149) * density(43) 
  pd(27,43) = pd(27,43) + rrt(149) * density(01) 
  pd(43,01) = pd(43,01) - rrt(149) * density(43) 
  pd(43,43) = pd(43,43) - rrt(149) * density(01) 
  pd(01,01) = pd(01,01) + rrt(150) * density(36) 
  pd(01,36) = pd(01,36) + rrt(150) * density(01) 
  pd(02,01) = pd(02,01) + rrt(150) * density(36) 
  pd(02,36) = pd(02,36) + rrt(150) * density(01) 
  pd(27,01) = pd(27,01) + rrt(150) * density(36) 
  pd(27,36) = pd(27,36) + rrt(150) * density(01) 
  pd(36,01) = pd(36,01) - rrt(150) * density(36) 
  pd(36,36) = pd(36,36) - rrt(150) * density(01) 
  pd(01,01) = pd(01,01) + rrt(151) * density(43) 
  pd(01,43) = pd(01,43) + rrt(151) * density(01) 
  pd(02,01) = pd(02,01) + rrt(151) * density(43) 
  pd(02,43) = pd(02,43) + rrt(151) * density(01) 
  pd(13,01) = pd(13,01) + rrt(151) * density(43) 
  pd(13,43) = pd(13,43) + rrt(151) * density(01) 
  pd(28,01) = pd(28,01) + rrt(151) * density(43) 
  pd(28,43) = pd(28,43) + rrt(151) * density(01) 
  pd(43,01) = pd(43,01) - rrt(151) * density(43) 
  pd(43,43) = pd(43,43) - rrt(151) * density(01) 
  pd(01,01) = pd(01,01) + rrt(152) * density(36) 
  pd(01,36) = pd(01,36) + rrt(152) * density(01) 
  pd(02,01) = pd(02,01) + rrt(152) * density(36) 
  pd(02,36) = pd(02,36) + rrt(152) * density(01) 
  pd(13,01) = pd(13,01) + rrt(152) * density(36) 
  pd(13,36) = pd(13,36) + rrt(152) * density(01) 
  pd(28,01) = pd(28,01) + rrt(152) * density(36) 
  pd(28,36) = pd(28,36) + rrt(152) * density(01) 
  pd(36,01) = pd(36,01) - rrt(152) * density(36) 
  pd(36,36) = pd(36,36) - rrt(152) * density(01) 
  pd(13,01) = pd(13,01) + rrt(153) * density(33) 
  pd(13,33) = pd(13,33) + rrt(153) * density(01) 
  pd(33,01) = pd(33,01) - rrt(153) * density(33) 
  pd(33,33) = pd(33,33) - rrt(153) * density(01) 
  pd(38,01) = pd(38,01) + rrt(153) * density(33) 
  pd(38,33) = pd(38,33) + rrt(153) * density(01) 
  pd(13,01) = pd(13,01) + rrt(154) * density(34) 
  pd(13,34) = pd(13,34) + rrt(154) * density(01) 
  pd(34,01) = pd(34,01) - rrt(154) * density(34) 
  pd(34,34) = pd(34,34) - rrt(154) * density(01) 
  pd(38,01) = pd(38,01) + rrt(154) * density(34) 
  pd(38,34) = pd(38,34) + rrt(154) * density(01) 
  pd(02,01) = pd(02,01) + rrt(155) * density(33) 
  pd(02,33) = pd(02,33) + rrt(155) * density(01) 
  pd(33,01) = pd(33,01) - rrt(155) * density(33) 
  pd(33,33) = pd(33,33) - rrt(155) * density(01) 
  pd(35,01) = pd(35,01) + rrt(155) * density(33) 
  pd(35,33) = pd(35,33) + rrt(155) * density(01) 
  pd(02,01) = pd(02,01) + rrt(156) * density(34) 
  pd(02,34) = pd(02,34) + rrt(156) * density(01) 
  pd(34,01) = pd(34,01) - rrt(156) * density(34) 
  pd(34,34) = pd(34,34) - rrt(156) * density(01) 
  pd(35,01) = pd(35,01) + rrt(156) * density(34) 
  pd(35,34) = pd(35,34) + rrt(156) * density(01) 
  pd(02,01) = pd(02,01) + rrt(157) * density(33) 
  pd(02,33) = pd(02,33) + rrt(157) * density(01) 
  pd(13,01) = pd(13,01) + rrt(157) * density(33) 
  pd(13,33) = pd(13,33) + rrt(157) * density(01) 
  pd(33,01) = pd(33,01) - rrt(157) * density(33) 
  pd(33,33) = pd(33,33) - rrt(157) * density(01) 
  pd(50,01) = pd(50,01) + rrt(157) * density(33) 
  pd(50,33) = pd(50,33) + rrt(157) * density(01) 
  pd(02,01) = pd(02,01) + rrt(158) * density(34) 
  pd(02,34) = pd(02,34) + rrt(158) * density(01) 
  pd(13,01) = pd(13,01) + rrt(158) * density(34) 
  pd(13,34) = pd(13,34) + rrt(158) * density(01) 
  pd(34,01) = pd(34,01) - rrt(158) * density(34) 
  pd(34,34) = pd(34,34) - rrt(158) * density(01) 
  pd(50,01) = pd(50,01) + rrt(158) * density(34) 
  pd(50,34) = pd(50,34) + rrt(158) * density(01) 
  pd(02,01) = pd(02,01) + rrt(159) * density(33) * 2.0d0
  pd(02,33) = pd(02,33) + rrt(159) * density(01) * 2.0d0
  pd(19,01) = pd(19,01) + rrt(159) * density(33) 
  pd(19,33) = pd(19,33) + rrt(159) * density(01) 
  pd(33,01) = pd(33,01) - rrt(159) * density(33) 
  pd(33,33) = pd(33,33) - rrt(159) * density(01) 
  pd(02,01) = pd(02,01) + rrt(160) * density(34) * 2.0d0
  pd(02,34) = pd(02,34) + rrt(160) * density(01) * 2.0d0
  pd(19,01) = pd(19,01) + rrt(160) * density(34) 
  pd(19,34) = pd(19,34) + rrt(160) * density(01) 
  pd(34,01) = pd(34,01) - rrt(160) * density(34) 
  pd(34,34) = pd(34,34) - rrt(160) * density(01) 
  pd(05,01) = pd(05,01) + rrt(161) * density(33) 
  pd(05,33) = pd(05,33) + rrt(161) * density(01) 
  pd(10,01) = pd(10,01) + rrt(161) * density(33) 
  pd(10,33) = pd(10,33) + rrt(161) * density(01) 
  pd(33,01) = pd(33,01) - rrt(161) * density(33) 
  pd(33,33) = pd(33,33) - rrt(161) * density(01) 
  pd(05,01) = pd(05,01) + rrt(162) * density(34) 
  pd(05,34) = pd(05,34) + rrt(162) * density(01) 
  pd(10,01) = pd(10,01) + rrt(162) * density(34) 
  pd(10,34) = pd(10,34) + rrt(162) * density(01) 
  pd(34,01) = pd(34,01) - rrt(162) * density(34) 
  pd(34,34) = pd(34,34) - rrt(162) * density(01) 
  pd(24,01) = pd(24,01) + rrt(163) * density(33) * 2.0d0
  pd(24,33) = pd(24,33) + rrt(163) * density(01) * 2.0d0
  pd(33,01) = pd(33,01) - rrt(163) * density(33) 
  pd(33,33) = pd(33,33) - rrt(163) * density(01) 
  pd(24,01) = pd(24,01) + rrt(164) * density(34) * 2.0d0
  pd(24,34) = pd(24,34) + rrt(164) * density(01) * 2.0d0
  pd(34,01) = pd(34,01) - rrt(164) * density(34) 
  pd(34,34) = pd(34,34) - rrt(164) * density(01) 
  pd(13,01) = pd(13,01) + rrt(165) * density(16) 
  pd(13,16) = pd(13,16) + rrt(165) * density(01) 
  pd(16,01) = pd(16,01) - rrt(165) * density(16) 
  pd(16,16) = pd(16,16) - rrt(165) * density(01) 
  pd(50,01) = pd(50,01) + rrt(165) * density(16) 
  pd(50,16) = pd(50,16) + rrt(165) * density(01) 
  pd(04,01) = pd(04,01) - rrt(166) * density(04) 
  pd(04,04) = pd(04,04) - rrt(166) * density(01) 
  pd(13,01) = pd(13,01) + rrt(166) * density(04) 
  pd(13,04) = pd(13,04) + rrt(166) * density(01) 
  pd(50,01) = pd(50,01) + rrt(166) * density(04) 
  pd(50,04) = pd(50,04) + rrt(166) * density(01) 
  pd(02,01) = pd(02,01) + rrt(167) * density(16) 
  pd(02,16) = pd(02,16) + rrt(167) * density(01) 
  pd(16,01) = pd(16,01) - rrt(167) * density(16) 
  pd(16,16) = pd(16,16) - rrt(167) * density(01) 
  pd(19,01) = pd(19,01) + rrt(167) * density(16) 
  pd(19,16) = pd(19,16) + rrt(167) * density(01) 
  pd(02,01) = pd(02,01) + rrt(168) * density(04) 
  pd(02,04) = pd(02,04) + rrt(168) * density(01) 
  pd(04,01) = pd(04,01) - rrt(168) * density(04) 
  pd(04,04) = pd(04,04) - rrt(168) * density(01) 
  pd(19,01) = pd(19,01) + rrt(168) * density(04) 
  pd(19,04) = pd(19,04) + rrt(168) * density(01) 
  pd(13,01) = pd(13,01) + rrt(169) * density(16) * 2.0d0
  pd(13,16) = pd(13,16) + rrt(169) * density(01) * 2.0d0
  pd(16,01) = pd(16,01) - rrt(169) * density(16) 
  pd(16,16) = pd(16,16) - rrt(169) * density(01) 
  pd(19,01) = pd(19,01) + rrt(169) * density(16) 
  pd(19,16) = pd(19,16) + rrt(169) * density(01) 
  pd(04,01) = pd(04,01) - rrt(170) * density(04) 
  pd(04,04) = pd(04,04) - rrt(170) * density(01) 
  pd(13,01) = pd(13,01) + rrt(170) * density(04) * 2.0d0
  pd(13,04) = pd(13,04) + rrt(170) * density(01) * 2.0d0
  pd(19,01) = pd(19,01) + rrt(170) * density(04) 
  pd(19,04) = pd(19,04) + rrt(170) * density(01) 
  pd(16,01) = pd(16,01) - rrt(171) * density(16) 
  pd(16,16) = pd(16,16) - rrt(171) * density(01) 
  pd(24,01) = pd(24,01) + rrt(171) * density(16) 
  pd(24,16) = pd(24,16) + rrt(171) * density(01) 
  pd(46,01) = pd(46,01) + rrt(171) * density(16) 
  pd(46,16) = pd(46,16) + rrt(171) * density(01) 
  pd(04,01) = pd(04,01) - rrt(172) * density(04) 
  pd(04,04) = pd(04,04) - rrt(172) * density(01) 
  pd(24,01) = pd(24,01) + rrt(172) * density(04) 
  pd(24,04) = pd(24,04) + rrt(172) * density(01) 
  pd(46,01) = pd(46,01) + rrt(172) * density(04) 
  pd(46,04) = pd(46,04) + rrt(172) * density(01) 
  pd(10,01) = pd(10,01) + rrt(173) * density(16) * 2.0d0
  pd(10,16) = pd(10,16) + rrt(173) * density(01) * 2.0d0
  pd(16,01) = pd(16,01) - rrt(173) * density(16) 
  pd(16,16) = pd(16,16) - rrt(173) * density(01) 
  pd(04,01) = pd(04,01) - rrt(174) * density(04) 
  pd(04,04) = pd(04,04) - rrt(174) * density(01) 
  pd(10,01) = pd(10,01) + rrt(174) * density(04) * 2.0d0
  pd(10,04) = pd(10,04) + rrt(174) * density(01) * 2.0d0
  pd(21,01) = pd(21,01) - rrt(175) * density(21) 
  pd(21,21) = pd(21,21) - rrt(175) * density(01) 
  pd(46,01) = pd(46,01) + rrt(175) * density(21) * 2.0d0
  pd(46,21) = pd(46,21) + rrt(175) * density(01) * 2.0d0
  pd(31,01) = pd(31,01) - rrt(176) * density(31) 
  pd(31,31) = pd(31,31) - rrt(176) * density(01) 
  pd(46,01) = pd(46,01) + rrt(176) * density(31) * 2.0d0
  pd(46,31) = pd(46,31) + rrt(176) * density(01) * 2.0d0
  pd(29,01) = pd(29,01) - rrt(177) * density(29) 
  pd(29,29) = pd(29,29) - rrt(177) * density(01) 
  pd(46,01) = pd(46,01) + rrt(177) * density(29) * 2.0d0
  pd(46,29) = pd(46,29) + rrt(177) * density(01) * 2.0d0
  pd(13,01) = pd(13,01) + rrt(178) * density(37) 
  pd(13,37) = pd(13,37) + rrt(178) * density(01) 
  pd(23,01) = pd(23,01) + rrt(178) * density(37) 
  pd(23,37) = pd(23,37) + rrt(178) * density(01) 
  pd(37,01) = pd(37,01) - rrt(178) * density(37) 
  pd(37,37) = pd(37,37) - rrt(178) * density(01) 
  pd(13,01) = pd(13,01) + rrt(179) * density(40) 
  pd(13,40) = pd(13,40) + rrt(179) * density(01) 
  pd(23,01) = pd(23,01) + rrt(179) * density(40) 
  pd(23,40) = pd(23,40) + rrt(179) * density(01) 
  pd(40,01) = pd(40,01) - rrt(179) * density(40) 
  pd(40,40) = pd(40,40) - rrt(179) * density(01) 
  pd(02,01) = pd(02,01) + rrt(180) * density(37) 
  pd(02,37) = pd(02,37) + rrt(180) * density(01) 
  pd(09,01) = pd(09,01) + rrt(180) * density(37) 
  pd(09,37) = pd(09,37) + rrt(180) * density(01) 
  pd(37,01) = pd(37,01) - rrt(180) * density(37) 
  pd(37,37) = pd(37,37) - rrt(180) * density(01) 
  pd(02,01) = pd(02,01) + rrt(181) * density(40) 
  pd(02,40) = pd(02,40) + rrt(181) * density(01) 
  pd(09,01) = pd(09,01) + rrt(181) * density(40) 
  pd(09,40) = pd(09,40) + rrt(181) * density(01) 
  pd(40,01) = pd(40,01) - rrt(181) * density(40) 
  pd(40,40) = pd(40,40) - rrt(181) * density(01) 
  pd(02,01) = pd(02,01) + rrt(182) * density(37) * 2.0d0
  pd(02,37) = pd(02,37) + rrt(182) * density(01) * 2.0d0
  pd(17,01) = pd(17,01) + rrt(182) * density(37) 
  pd(17,37) = pd(17,37) + rrt(182) * density(01) 
  pd(37,01) = pd(37,01) - rrt(182) * density(37) 
  pd(37,37) = pd(37,37) - rrt(182) * density(01) 
  pd(02,01) = pd(02,01) + rrt(183) * density(40) * 2.0d0
  pd(02,40) = pd(02,40) + rrt(183) * density(01) * 2.0d0
  pd(17,01) = pd(17,01) + rrt(183) * density(40) 
  pd(17,40) = pd(17,40) + rrt(183) * density(01) 
  pd(40,01) = pd(40,01) - rrt(183) * density(40) 
  pd(40,40) = pd(40,40) - rrt(183) * density(01) 
  pd(10,01) = pd(10,01) + rrt(184) * density(37) 
  pd(10,37) = pd(10,37) + rrt(184) * density(01) 
  pd(20,01) = pd(20,01) + rrt(184) * density(37) 
  pd(20,37) = pd(20,37) + rrt(184) * density(01) 
  pd(37,01) = pd(37,01) - rrt(184) * density(37) 
  pd(37,37) = pd(37,37) - rrt(184) * density(01) 
  pd(10,01) = pd(10,01) + rrt(185) * density(40) 
  pd(10,40) = pd(10,40) + rrt(185) * density(01) 
  pd(20,01) = pd(20,01) + rrt(185) * density(40) 
  pd(20,40) = pd(20,40) + rrt(185) * density(01) 
  pd(40,01) = pd(40,01) - rrt(185) * density(40) 
  pd(40,40) = pd(40,40) - rrt(185) * density(01) 
  pd(24,01) = pd(24,01) + rrt(186) * density(37) 
  pd(24,37) = pd(24,37) + rrt(186) * density(01) 
  pd(37,01) = pd(37,01) - rrt(186) * density(37) 
  pd(37,37) = pd(37,37) - rrt(186) * density(01) 
  pd(38,01) = pd(38,01) + rrt(186) * density(37) 
  pd(38,37) = pd(38,37) + rrt(186) * density(01) 
  pd(24,01) = pd(24,01) + rrt(187) * density(40) 
  pd(24,40) = pd(24,40) + rrt(187) * density(01) 
  pd(38,01) = pd(38,01) + rrt(187) * density(40) 
  pd(38,40) = pd(38,40) + rrt(187) * density(01) 
  pd(40,01) = pd(40,01) - rrt(187) * density(40) 
  pd(40,40) = pd(40,40) - rrt(187) * density(01) 
  pd(05,01) = pd(05,01) + rrt(188) * density(37) 
  pd(05,37) = pd(05,37) + rrt(188) * density(01) 
  pd(35,01) = pd(35,01) + rrt(188) * density(37) 
  pd(35,37) = pd(35,37) + rrt(188) * density(01) 
  pd(37,01) = pd(37,01) - rrt(188) * density(37) 
  pd(37,37) = pd(37,37) - rrt(188) * density(01) 
  pd(05,01) = pd(05,01) + rrt(189) * density(40) 
  pd(05,40) = pd(05,40) + rrt(189) * density(01) 
  pd(35,01) = pd(35,01) + rrt(189) * density(40) 
  pd(35,40) = pd(35,40) + rrt(189) * density(01) 
  pd(40,01) = pd(40,01) - rrt(189) * density(40) 
  pd(40,40) = pd(40,40) - rrt(189) * density(01) 
  pd(01,01) = pd(01,01) + rrt(190) * density(33) 
  pd(01,33) = pd(01,33) + rrt(190) * density(01) 
  pd(15,01) = pd(15,01) + rrt(190) * density(33) 
  pd(15,33) = pd(15,33) + rrt(190) * density(01) 
  pd(33,01) = pd(33,01) - rrt(190) * density(33) 
  pd(33,33) = pd(33,33) - rrt(190) * density(01) 
  pd(01,01) = pd(01,01) + rrt(191) * density(34) 
  pd(01,34) = pd(01,34) + rrt(191) * density(01) 
  pd(15,01) = pd(15,01) + rrt(191) * density(34) 
  pd(15,34) = pd(15,34) + rrt(191) * density(01) 
  pd(34,01) = pd(34,01) - rrt(191) * density(34) 
  pd(34,34) = pd(34,34) - rrt(191) * density(01) 
  pd(01,01) = pd(01,01) + rrt(192) * density(33) 
  pd(01,33) = pd(01,33) + rrt(192) * density(01) 
  pd(13,01) = pd(13,01) + rrt(192) * density(33) 
  pd(13,33) = pd(13,33) + rrt(192) * density(01) 
  pd(33,01) = pd(33,01) - rrt(192) * density(33) 
  pd(33,33) = pd(33,33) - rrt(192) * density(01) 
  pd(44,01) = pd(44,01) + rrt(192) * density(33) 
  pd(44,33) = pd(44,33) + rrt(192) * density(01) 
  pd(01,01) = pd(01,01) + rrt(193) * density(34) 
  pd(01,34) = pd(01,34) + rrt(193) * density(01) 
  pd(13,01) = pd(13,01) + rrt(193) * density(34) 
  pd(13,34) = pd(13,34) + rrt(193) * density(01) 
  pd(34,01) = pd(34,01) - rrt(193) * density(34) 
  pd(34,34) = pd(34,34) - rrt(193) * density(01) 
  pd(44,01) = pd(44,01) + rrt(193) * density(34) 
  pd(44,34) = pd(44,34) + rrt(193) * density(01) 
  pd(01,01) = pd(01,01) + rrt(194) * density(33) 
  pd(01,33) = pd(01,33) + rrt(194) * density(01) 
  pd(02,01) = pd(02,01) + rrt(194) * density(33) 
  pd(02,33) = pd(02,33) + rrt(194) * density(01) 
  pd(32,01) = pd(32,01) + rrt(194) * density(33) 
  pd(32,33) = pd(32,33) + rrt(194) * density(01) 
  pd(33,01) = pd(33,01) - rrt(194) * density(33) 
  pd(33,33) = pd(33,33) - rrt(194) * density(01) 
  pd(01,01) = pd(01,01) + rrt(195) * density(34) 
  pd(01,34) = pd(01,34) + rrt(195) * density(01) 
  pd(02,01) = pd(02,01) + rrt(195) * density(34) 
  pd(02,34) = pd(02,34) + rrt(195) * density(01) 
  pd(32,01) = pd(32,01) + rrt(195) * density(34) 
  pd(32,34) = pd(32,34) + rrt(195) * density(01) 
  pd(34,01) = pd(34,01) - rrt(195) * density(34) 
  pd(34,34) = pd(34,34) - rrt(195) * density(01) 
  pd(01,01) = pd(01,01) + rrt(196) * density(33) 
  pd(01,33) = pd(01,33) + rrt(196) * density(01) 
  pd(02,01) = pd(02,01) + rrt(196) * density(33) 
  pd(02,33) = pd(02,33) + rrt(196) * density(01) 
  pd(06,01) = pd(06,01) + rrt(196) * density(33) 
  pd(06,33) = pd(06,33) + rrt(196) * density(01) 
  pd(13,01) = pd(13,01) + rrt(196) * density(33) 
  pd(13,33) = pd(13,33) + rrt(196) * density(01) 
  pd(33,01) = pd(33,01) - rrt(196) * density(33) 
  pd(33,33) = pd(33,33) - rrt(196) * density(01) 
  pd(01,01) = pd(01,01) + rrt(197) * density(34) 
  pd(01,34) = pd(01,34) + rrt(197) * density(01) 
  pd(02,01) = pd(02,01) + rrt(197) * density(34) 
  pd(02,34) = pd(02,34) + rrt(197) * density(01) 
  pd(06,01) = pd(06,01) + rrt(197) * density(34) 
  pd(06,34) = pd(06,34) + rrt(197) * density(01) 
  pd(13,01) = pd(13,01) + rrt(197) * density(34) 
  pd(13,34) = pd(13,34) + rrt(197) * density(01) 
  pd(34,01) = pd(34,01) - rrt(197) * density(34) 
  pd(34,34) = pd(34,34) - rrt(197) * density(01) 
  pd(01,01) = pd(01,01) + rrt(198) * density(33) 
  pd(01,33) = pd(01,33) + rrt(198) * density(01) 
  pd(02,01) = pd(02,01) + rrt(198) * density(33) * 2.0d0
  pd(02,33) = pd(02,33) + rrt(198) * density(01) * 2.0d0
  pd(14,01) = pd(14,01) + rrt(198) * density(33) 
  pd(14,33) = pd(14,33) + rrt(198) * density(01) 
  pd(33,01) = pd(33,01) - rrt(198) * density(33) 
  pd(33,33) = pd(33,33) - rrt(198) * density(01) 
  pd(01,01) = pd(01,01) + rrt(199) * density(34) 
  pd(01,34) = pd(01,34) + rrt(199) * density(01) 
  pd(02,01) = pd(02,01) + rrt(199) * density(34) * 2.0d0
  pd(02,34) = pd(02,34) + rrt(199) * density(01) * 2.0d0
  pd(14,01) = pd(14,01) + rrt(199) * density(34) 
  pd(14,34) = pd(14,34) + rrt(199) * density(01) 
  pd(34,01) = pd(34,01) - rrt(199) * density(34) 
  pd(34,34) = pd(34,34) - rrt(199) * density(01) 
  pd(01,01) = pd(01,01) + rrt(200) * density(33) 
  pd(01,33) = pd(01,33) + rrt(200) * density(01) 
  pd(24,01) = pd(24,01) + rrt(200) * density(33) 
  pd(24,33) = pd(24,33) + rrt(200) * density(01) 
  pd(33,01) = pd(33,01) - rrt(200) * density(33) 
  pd(33,33) = pd(33,33) - rrt(200) * density(01) 
  pd(48,01) = pd(48,01) + rrt(200) * density(33) 
  pd(48,33) = pd(48,33) + rrt(200) * density(01) 
  pd(01,01) = pd(01,01) + rrt(201) * density(34) 
  pd(01,34) = pd(01,34) + rrt(201) * density(01) 
  pd(24,01) = pd(24,01) + rrt(201) * density(34) 
  pd(24,34) = pd(24,34) + rrt(201) * density(01) 
  pd(34,01) = pd(34,01) - rrt(201) * density(34) 
  pd(34,34) = pd(34,34) - rrt(201) * density(01) 
  pd(48,01) = pd(48,01) + rrt(201) * density(34) 
  pd(48,34) = pd(48,34) + rrt(201) * density(01) 
  pd(01,01) = pd(01,01) + rrt(202) * density(33) 
  pd(01,33) = pd(01,33) + rrt(202) * density(01) 
  pd(05,01) = pd(05,01) + rrt(202) * density(33) 
  pd(05,33) = pd(05,33) + rrt(202) * density(01) 
  pd(27,01) = pd(27,01) + rrt(202) * density(33) 
  pd(27,33) = pd(27,33) + rrt(202) * density(01) 
  pd(33,01) = pd(33,01) - rrt(202) * density(33) 
  pd(33,33) = pd(33,33) - rrt(202) * density(01) 
  pd(01,01) = pd(01,01) + rrt(203) * density(34) 
  pd(01,34) = pd(01,34) + rrt(203) * density(01) 
  pd(05,01) = pd(05,01) + rrt(203) * density(34) 
  pd(05,34) = pd(05,34) + rrt(203) * density(01) 
  pd(27,01) = pd(27,01) + rrt(203) * density(34) 
  pd(27,34) = pd(27,34) + rrt(203) * density(01) 
  pd(34,01) = pd(34,01) - rrt(203) * density(34) 
  pd(34,34) = pd(34,34) - rrt(203) * density(01) 
  pd(01,01) = pd(01,01) + rrt(204) * density(16) 
  pd(01,16) = pd(01,16) + rrt(204) * density(01) 
  pd(16,01) = pd(16,01) - rrt(204) * density(16) 
  pd(16,16) = pd(16,16) - rrt(204) * density(01) 
  pd(32,01) = pd(32,01) + rrt(204) * density(16) 
  pd(32,16) = pd(32,16) + rrt(204) * density(01) 
  pd(01,01) = pd(01,01) + rrt(205) * density(04) 
  pd(01,04) = pd(01,04) + rrt(205) * density(01) 
  pd(04,01) = pd(04,01) - rrt(205) * density(04) 
  pd(04,04) = pd(04,04) - rrt(205) * density(01) 
  pd(32,01) = pd(32,01) + rrt(205) * density(04) 
  pd(32,04) = pd(32,04) + rrt(205) * density(01) 
  pd(01,01) = pd(01,01) + rrt(206) * density(16) 
  pd(01,16) = pd(01,16) + rrt(206) * density(01) 
  pd(06,01) = pd(06,01) + rrt(206) * density(16) 
  pd(06,16) = pd(06,16) + rrt(206) * density(01) 
  pd(13,01) = pd(13,01) + rrt(206) * density(16) 
  pd(13,16) = pd(13,16) + rrt(206) * density(01) 
  pd(16,01) = pd(16,01) - rrt(206) * density(16) 
  pd(16,16) = pd(16,16) - rrt(206) * density(01) 
  pd(01,01) = pd(01,01) + rrt(207) * density(04) 
  pd(01,04) = pd(01,04) + rrt(207) * density(01) 
  pd(04,01) = pd(04,01) - rrt(207) * density(04) 
  pd(04,04) = pd(04,04) - rrt(207) * density(01) 
  pd(06,01) = pd(06,01) + rrt(207) * density(04) 
  pd(06,04) = pd(06,04) + rrt(207) * density(01) 
  pd(13,01) = pd(13,01) + rrt(207) * density(04) 
  pd(13,04) = pd(13,04) + rrt(207) * density(01) 
  pd(01,01) = pd(01,01) + rrt(208) * density(16) 
  pd(01,16) = pd(01,16) + rrt(208) * density(01) 
  pd(16,01) = pd(16,01) - rrt(208) * density(16) 
  pd(16,16) = pd(16,16) - rrt(208) * density(01) 
  pd(46,01) = pd(46,01) + rrt(208) * density(16) 
  pd(46,16) = pd(46,16) + rrt(208) * density(01) 
  pd(48,01) = pd(48,01) + rrt(208) * density(16) 
  pd(48,16) = pd(48,16) + rrt(208) * density(01) 
  pd(01,01) = pd(01,01) + rrt(209) * density(04) 
  pd(01,04) = pd(01,04) + rrt(209) * density(01) 
  pd(04,01) = pd(04,01) - rrt(209) * density(04) 
  pd(04,04) = pd(04,04) - rrt(209) * density(01) 
  pd(46,01) = pd(46,01) + rrt(209) * density(04) 
  pd(46,04) = pd(46,04) + rrt(209) * density(01) 
  pd(48,01) = pd(48,01) + rrt(209) * density(04) 
  pd(48,04) = pd(48,04) + rrt(209) * density(01) 
  pd(01,01) = pd(01,01) + rrt(210) * density(16) 
  pd(01,16) = pd(01,16) + rrt(210) * density(01) 
  pd(10,01) = pd(10,01) + rrt(210) * density(16) 
  pd(10,16) = pd(10,16) + rrt(210) * density(01) 
  pd(16,01) = pd(16,01) - rrt(210) * density(16) 
  pd(16,16) = pd(16,16) - rrt(210) * density(01) 
  pd(27,01) = pd(27,01) + rrt(210) * density(16) 
  pd(27,16) = pd(27,16) + rrt(210) * density(01) 
  pd(01,01) = pd(01,01) + rrt(211) * density(04) 
  pd(01,04) = pd(01,04) + rrt(211) * density(01) 
  pd(04,01) = pd(04,01) - rrt(211) * density(04) 
  pd(04,04) = pd(04,04) - rrt(211) * density(01) 
  pd(10,01) = pd(10,01) + rrt(211) * density(04) 
  pd(10,04) = pd(10,04) + rrt(211) * density(01) 
  pd(27,01) = pd(27,01) + rrt(211) * density(04) 
  pd(27,04) = pd(27,04) + rrt(211) * density(01) 
  pd(01,01) = pd(01,01) + rrt(212) * density(16) 
  pd(01,16) = pd(01,16) + rrt(212) * density(01) 
  pd(16,01) = pd(16,01) - rrt(212) * density(16) 
  pd(16,16) = pd(16,16) - rrt(212) * density(01) 
  pd(24,01) = pd(24,01) + rrt(212) * density(16) 
  pd(24,16) = pd(24,16) + rrt(212) * density(01) 
  pd(28,01) = pd(28,01) + rrt(212) * density(16) 
  pd(28,16) = pd(28,16) + rrt(212) * density(01) 
  pd(01,01) = pd(01,01) + rrt(213) * density(04) 
  pd(01,04) = pd(01,04) + rrt(213) * density(01) 
  pd(04,01) = pd(04,01) - rrt(213) * density(04) 
  pd(04,04) = pd(04,04) - rrt(213) * density(01) 
  pd(24,01) = pd(24,01) + rrt(213) * density(04) 
  pd(24,04) = pd(24,04) + rrt(213) * density(01) 
  pd(28,01) = pd(28,01) + rrt(213) * density(04) 
  pd(28,04) = pd(28,04) + rrt(213) * density(01) 
  pd(01,01) = pd(01,01) + rrt(214) * density(21) 
  pd(01,21) = pd(01,21) + rrt(214) * density(01) 
  pd(14,01) = pd(14,01) + rrt(214) * density(21) 
  pd(14,21) = pd(14,21) + rrt(214) * density(01) 
  pd(21,01) = pd(21,01) - rrt(214) * density(21) 
  pd(21,21) = pd(21,21) - rrt(214) * density(01) 
  pd(01,01) = pd(01,01) + rrt(215) * density(31) 
  pd(01,31) = pd(01,31) + rrt(215) * density(01) 
  pd(14,01) = pd(14,01) + rrt(215) * density(31) 
  pd(14,31) = pd(14,31) + rrt(215) * density(01) 
  pd(31,01) = pd(31,01) - rrt(215) * density(31) 
  pd(31,31) = pd(31,31) - rrt(215) * density(01) 
  pd(01,01) = pd(01,01) + rrt(216) * density(29) 
  pd(01,29) = pd(01,29) + rrt(216) * density(01) 
  pd(14,01) = pd(14,01) + rrt(216) * density(29) 
  pd(14,29) = pd(14,29) + rrt(216) * density(01) 
  pd(29,01) = pd(29,01) - rrt(216) * density(29) 
  pd(29,29) = pd(29,29) - rrt(216) * density(01) 
  pd(01,01) = pd(01,01) + rrt(217) * density(21) 
  pd(01,21) = pd(01,21) + rrt(217) * density(01) 
  pd(21,01) = pd(21,01) - rrt(217) * density(21) 
  pd(21,21) = pd(21,21) - rrt(217) * density(01) 
  pd(28,01) = pd(28,01) + rrt(217) * density(21) 
  pd(28,21) = pd(28,21) + rrt(217) * density(01) 
  pd(46,01) = pd(46,01) + rrt(217) * density(21) 
  pd(46,21) = pd(46,21) + rrt(217) * density(01) 
  pd(01,01) = pd(01,01) + rrt(218) * density(31) 
  pd(01,31) = pd(01,31) + rrt(218) * density(01) 
  pd(28,01) = pd(28,01) + rrt(218) * density(31) 
  pd(28,31) = pd(28,31) + rrt(218) * density(01) 
  pd(31,01) = pd(31,01) - rrt(218) * density(31) 
  pd(31,31) = pd(31,31) - rrt(218) * density(01) 
  pd(46,01) = pd(46,01) + rrt(218) * density(31) 
  pd(46,31) = pd(46,31) + rrt(218) * density(01) 
  pd(01,01) = pd(01,01) + rrt(219) * density(29) 
  pd(01,29) = pd(01,29) + rrt(219) * density(01) 
  pd(28,01) = pd(28,01) + rrt(219) * density(29) 
  pd(28,29) = pd(28,29) + rrt(219) * density(01) 
  pd(29,01) = pd(29,01) - rrt(219) * density(29) 
  pd(29,29) = pd(29,29) - rrt(219) * density(01) 
  pd(46,01) = pd(46,01) + rrt(219) * density(29) 
  pd(46,29) = pd(46,29) + rrt(219) * density(01) 
  pd(01,01) = pd(01,01) + rrt(220) * density(37) 
  pd(01,37) = pd(01,37) + rrt(220) * density(01) 
  pd(22,01) = pd(22,01) + rrt(220) * density(37) 
  pd(22,37) = pd(22,37) + rrt(220) * density(01) 
  pd(37,01) = pd(37,01) - rrt(220) * density(37) 
  pd(37,37) = pd(37,37) - rrt(220) * density(01) 
  pd(01,01) = pd(01,01) + rrt(221) * density(40) 
  pd(01,40) = pd(01,40) + rrt(221) * density(01) 
  pd(22,01) = pd(22,01) + rrt(221) * density(40) 
  pd(22,40) = pd(22,40) + rrt(221) * density(01) 
  pd(40,01) = pd(40,01) - rrt(221) * density(40) 
  pd(40,40) = pd(40,40) - rrt(221) * density(01) 
  pd(01,01) = pd(01,01) + rrt(222) * density(37) 
  pd(01,37) = pd(01,37) + rrt(222) * density(01) 
  pd(13,01) = pd(13,01) + rrt(222) * density(37) 
  pd(13,37) = pd(13,37) + rrt(222) * density(01) 
  pd(30,01) = pd(30,01) + rrt(222) * density(37) 
  pd(30,37) = pd(30,37) + rrt(222) * density(01) 
  pd(37,01) = pd(37,01) - rrt(222) * density(37) 
  pd(37,37) = pd(37,37) - rrt(222) * density(01) 
  pd(01,01) = pd(01,01) + rrt(223) * density(40) 
  pd(01,40) = pd(01,40) + rrt(223) * density(01) 
  pd(13,01) = pd(13,01) + rrt(223) * density(40) 
  pd(13,40) = pd(13,40) + rrt(223) * density(01) 
  pd(30,01) = pd(30,01) + rrt(223) * density(40) 
  pd(30,40) = pd(30,40) + rrt(223) * density(01) 
  pd(40,01) = pd(40,01) - rrt(223) * density(40) 
  pd(40,40) = pd(40,40) - rrt(223) * density(01) 
  pd(01,01) = pd(01,01) + rrt(224) * density(37) 
  pd(01,37) = pd(01,37) + rrt(224) * density(01) 
  pd(02,01) = pd(02,01) + rrt(224) * density(37) 
  pd(02,37) = pd(02,37) + rrt(224) * density(01) 
  pd(37,01) = pd(37,01) - rrt(224) * density(37) 
  pd(37,37) = pd(37,37) - rrt(224) * density(01) 
  pd(49,01) = pd(49,01) + rrt(224) * density(37) 
  pd(49,37) = pd(49,37) + rrt(224) * density(01) 
  pd(01,01) = pd(01,01) + rrt(225) * density(40) 
  pd(01,40) = pd(01,40) + rrt(225) * density(01) 
  pd(02,01) = pd(02,01) + rrt(225) * density(40) 
  pd(02,40) = pd(02,40) + rrt(225) * density(01) 
  pd(40,01) = pd(40,01) - rrt(225) * density(40) 
  pd(40,40) = pd(40,40) - rrt(225) * density(01) 
  pd(49,01) = pd(49,01) + rrt(225) * density(40) 
  pd(49,40) = pd(49,40) + rrt(225) * density(01) 
  pd(01,01) = pd(01,01) + rrt(226) * density(37) 
  pd(01,37) = pd(01,37) + rrt(226) * density(01) 
  pd(02,01) = pd(02,01) + rrt(226) * density(37) 
  pd(02,37) = pd(02,37) + rrt(226) * density(01) 
  pd(08,01) = pd(08,01) + rrt(226) * density(37) 
  pd(08,37) = pd(08,37) + rrt(226) * density(01) 
  pd(13,01) = pd(13,01) + rrt(226) * density(37) 
  pd(13,37) = pd(13,37) + rrt(226) * density(01) 
  pd(37,01) = pd(37,01) - rrt(226) * density(37) 
  pd(37,37) = pd(37,37) - rrt(226) * density(01) 
  pd(01,01) = pd(01,01) + rrt(227) * density(40) 
  pd(01,40) = pd(01,40) + rrt(227) * density(01) 
  pd(02,01) = pd(02,01) + rrt(227) * density(40) 
  pd(02,40) = pd(02,40) + rrt(227) * density(01) 
  pd(08,01) = pd(08,01) + rrt(227) * density(40) 
  pd(08,40) = pd(08,40) + rrt(227) * density(01) 
  pd(13,01) = pd(13,01) + rrt(227) * density(40) 
  pd(13,40) = pd(13,40) + rrt(227) * density(01) 
  pd(40,01) = pd(40,01) - rrt(227) * density(40) 
  pd(40,40) = pd(40,40) - rrt(227) * density(01) 
  pd(01,01) = pd(01,01) + rrt(228) * density(37) 
  pd(01,37) = pd(01,37) + rrt(228) * density(01) 
  pd(02,01) = pd(02,01) + rrt(228) * density(37) * 2.0d0
  pd(02,37) = pd(02,37) + rrt(228) * density(01) * 2.0d0
  pd(37,01) = pd(37,01) - rrt(228) * density(37) 
  pd(37,37) = pd(37,37) - rrt(228) * density(01) 
  pd(47,01) = pd(47,01) + rrt(228) * density(37) 
  pd(47,37) = pd(47,37) + rrt(228) * density(01) 
  pd(01,01) = pd(01,01) + rrt(229) * density(40) 
  pd(01,40) = pd(01,40) + rrt(229) * density(01) 
  pd(02,01) = pd(02,01) + rrt(229) * density(40) * 2.0d0
  pd(02,40) = pd(02,40) + rrt(229) * density(01) * 2.0d0
  pd(40,01) = pd(40,01) - rrt(229) * density(40) 
  pd(40,40) = pd(40,40) - rrt(229) * density(01) 
  pd(47,01) = pd(47,01) + rrt(229) * density(40) 
  pd(47,40) = pd(47,40) + rrt(229) * density(01) 
  pd(01,01) = pd(01,01) + rrt(230) * density(37) 
  pd(01,37) = pd(01,37) + rrt(230) * density(01) 
  pd(24,01) = pd(24,01) + rrt(230) * density(37) 
  pd(24,37) = pd(24,37) + rrt(230) * density(01) 
  pd(37,01) = pd(37,01) - rrt(230) * density(37) 
  pd(37,37) = pd(37,37) - rrt(230) * density(01) 
  pd(44,01) = pd(44,01) + rrt(230) * density(37) 
  pd(44,37) = pd(44,37) + rrt(230) * density(01) 
  pd(01,01) = pd(01,01) + rrt(231) * density(40) 
  pd(01,40) = pd(01,40) + rrt(231) * density(01) 
  pd(24,01) = pd(24,01) + rrt(231) * density(40) 
  pd(24,40) = pd(24,40) + rrt(231) * density(01) 
  pd(40,01) = pd(40,01) - rrt(231) * density(40) 
  pd(40,40) = pd(40,40) - rrt(231) * density(01) 
  pd(44,01) = pd(44,01) + rrt(231) * density(40) 
  pd(44,40) = pd(44,40) + rrt(231) * density(01) 
  pd(01,01) = pd(01,01) + rrt(232) * density(37) 
  pd(01,37) = pd(01,37) + rrt(232) * density(01) 
  pd(05,01) = pd(05,01) + rrt(232) * density(37) 
  pd(05,37) = pd(05,37) + rrt(232) * density(01) 
  pd(32,01) = pd(32,01) + rrt(232) * density(37) 
  pd(32,37) = pd(32,37) + rrt(232) * density(01) 
  pd(37,01) = pd(37,01) - rrt(232) * density(37) 
  pd(37,37) = pd(37,37) - rrt(232) * density(01) 
  pd(01,01) = pd(01,01) + rrt(233) * density(40) 
  pd(01,40) = pd(01,40) + rrt(233) * density(01) 
  pd(05,01) = pd(05,01) + rrt(233) * density(40) 
  pd(05,40) = pd(05,40) + rrt(233) * density(01) 
  pd(32,01) = pd(32,01) + rrt(233) * density(40) 
  pd(32,40) = pd(32,40) + rrt(233) * density(01) 
  pd(40,01) = pd(40,01) - rrt(233) * density(40) 
  pd(40,40) = pd(40,40) - rrt(233) * density(01) 
  pd(01,01) = pd(01,01) + rrt(234) * density(37) 
  pd(01,37) = pd(01,37) + rrt(234) * density(01) 
  pd(37,01) = pd(37,01) - rrt(234) * density(37) 
  pd(37,37) = pd(37,37) - rrt(234) * density(01) 
  pd(38,01) = pd(38,01) + rrt(234) * density(37) 
  pd(38,37) = pd(38,37) + rrt(234) * density(01) 
  pd(48,01) = pd(48,01) + rrt(234) * density(37) 
  pd(48,37) = pd(48,37) + rrt(234) * density(01) 
  pd(01,01) = pd(01,01) + rrt(235) * density(40) 
  pd(01,40) = pd(01,40) + rrt(235) * density(01) 
  pd(38,01) = pd(38,01) + rrt(235) * density(40) 
  pd(38,40) = pd(38,40) + rrt(235) * density(01) 
  pd(40,01) = pd(40,01) - rrt(235) * density(40) 
  pd(40,40) = pd(40,40) - rrt(235) * density(01) 
  pd(48,01) = pd(48,01) + rrt(235) * density(40) 
  pd(48,40) = pd(48,40) + rrt(235) * density(01) 
  pd(01,01) = pd(01,01) + rrt(236) * density(37) 
  pd(01,37) = pd(01,37) + rrt(236) * density(01) 
  pd(20,01) = pd(20,01) + rrt(236) * density(37) 
  pd(20,37) = pd(20,37) + rrt(236) * density(01) 
  pd(27,01) = pd(27,01) + rrt(236) * density(37) 
  pd(27,37) = pd(27,37) + rrt(236) * density(01) 
  pd(37,01) = pd(37,01) - rrt(236) * density(37) 
  pd(37,37) = pd(37,37) - rrt(236) * density(01) 
  pd(01,01) = pd(01,01) + rrt(237) * density(40) 
  pd(01,40) = pd(01,40) + rrt(237) * density(01) 
  pd(20,01) = pd(20,01) + rrt(237) * density(40) 
  pd(20,40) = pd(20,40) + rrt(237) * density(01) 
  pd(27,01) = pd(27,01) + rrt(237) * density(40) 
  pd(27,40) = pd(27,40) + rrt(237) * density(01) 
  pd(40,01) = pd(40,01) - rrt(237) * density(40) 
  pd(40,40) = pd(40,40) - rrt(237) * density(01) 
  pd(02,01) = pd(02,01) - rrt(238) * density(02) 
  pd(02,02) = pd(02,02) - rrt(238) * density(01) 
  pd(13,01) = pd(13,01) + rrt(238) * density(02) * 2.0d0
  pd(13,02) = pd(13,02) + rrt(238) * density(01) * 2.0d0
  pd(01,01) = pd(01,01) + rrt(239) * density(02) 
  pd(01,02) = pd(01,02) + rrt(239) * density(01) 
  pd(02,01) = pd(02,01) - rrt(239) * density(02) 
  pd(02,02) = pd(02,02) - rrt(239) * density(01) 
  pd(45,01) = pd(45,01) + rrt(239) * density(02) 
  pd(45,02) = pd(45,02) + rrt(239) * density(01) 
  pd(01,01) = pd(01,01) - rrt(240) * density(01) * density(25) * 2.0d0
  pd(01,25) = pd(01,25) - rrt(240) * density(01)**2 
  pd(05,01) = pd(05,01) + rrt(240) * density(01) * density(25) * 2.0d0
  pd(05,25) = pd(05,25) + rrt(240) * density(01)**2 
  pd(25,01) = pd(25,01) - rrt(240) * density(01) * density(25) * 2.0d0
  pd(25,25) = pd(25,25) - rrt(240) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(241) * density(01) * density(48) * 2.0d0
  pd(01,48) = pd(01,48) - rrt(241) * density(01)**2 
  pd(24,01) = pd(24,01) + rrt(241) * density(01) * density(48) * 2.0d0
  pd(24,48) = pd(24,48) + rrt(241) * density(01)**2 
  pd(48,01) = pd(48,01) - rrt(241) * density(01) * density(48) * 2.0d0
  pd(48,48) = pd(48,48) - rrt(241) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(242) * density(01) * density(27) * 2.0d0
  pd(01,27) = pd(01,27) - rrt(242) * density(01)**2 
  pd(10,01) = pd(10,01) + rrt(242) * density(01) * density(27) * 2.0d0
  pd(10,27) = pd(10,27) + rrt(242) * density(01)**2 
  pd(27,01) = pd(27,01) - rrt(242) * density(01) * density(27) * 2.0d0
  pd(27,27) = pd(27,27) - rrt(242) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(243) * density(01) * density(28) * 2.0d0
  pd(01,28) = pd(01,28) - rrt(243) * density(01)**2 
  pd(28,01) = pd(28,01) - rrt(243) * density(01) * density(28) * 2.0d0
  pd(28,28) = pd(28,28) - rrt(243) * density(01)**2 
  pd(46,01) = pd(46,01) + rrt(243) * density(01) * density(28) * 2.0d0
  pd(46,28) = pd(46,28) + rrt(243) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(244) * density(01) * density(15) * 2.0d0
  pd(01,15) = pd(01,15) - rrt(244) * density(01)**2 
  pd(15,01) = pd(15,01) - rrt(244) * density(01) * density(15) * 2.0d0
  pd(15,15) = pd(15,15) - rrt(244) * density(01)**2 
  pd(20,01) = pd(20,01) + rrt(244) * density(01) * density(15) * 2.0d0
  pd(20,15) = pd(20,15) + rrt(244) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(245) * density(01) * density(44) * 2.0d0
  pd(01,44) = pd(01,44) - rrt(245) * density(01)**2 
  pd(38,01) = pd(38,01) + rrt(245) * density(01) * density(44) * 2.0d0
  pd(38,44) = pd(38,44) + rrt(245) * density(01)**2 
  pd(44,01) = pd(44,01) - rrt(245) * density(01) * density(44) * 2.0d0
  pd(44,44) = pd(44,44) - rrt(245) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(246) * density(01) * density(32) * 2.0d0
  pd(01,32) = pd(01,32) - rrt(246) * density(01)**2 
  pd(32,01) = pd(32,01) - rrt(246) * density(01) * density(32) * 2.0d0
  pd(32,32) = pd(32,32) - rrt(246) * density(01)**2 
  pd(35,01) = pd(35,01) + rrt(246) * density(01) * density(32) * 2.0d0
  pd(35,32) = pd(35,32) + rrt(246) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(247) * density(01) * density(06) * 2.0d0
  pd(01,06) = pd(01,06) - rrt(247) * density(01)**2 
  pd(06,01) = pd(06,01) - rrt(247) * density(01) * density(06) * 2.0d0
  pd(06,06) = pd(06,06) - rrt(247) * density(01)**2 
  pd(50,01) = pd(50,01) + rrt(247) * density(01) * density(06) * 2.0d0
  pd(50,06) = pd(50,06) + rrt(247) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(248) * density(01) * density(14) * 2.0d0
  pd(01,14) = pd(01,14) - rrt(248) * density(01)**2 
  pd(14,01) = pd(14,01) - rrt(248) * density(01) * density(14) * 2.0d0
  pd(14,14) = pd(14,14) - rrt(248) * density(01)**2 
  pd(19,01) = pd(19,01) + rrt(248) * density(01) * density(14) * 2.0d0
  pd(19,14) = pd(19,14) + rrt(248) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(249) * density(01) * density(22) * 2.0d0
  pd(01,22) = pd(01,22) - rrt(249) * density(01)**2 
  pd(12,01) = pd(12,01) + rrt(249) * density(01) * density(22) * 2.0d0
  pd(12,22) = pd(12,22) + rrt(249) * density(01)**2 
  pd(22,01) = pd(22,01) - rrt(249) * density(01) * density(22) * 2.0d0
  pd(22,22) = pd(22,22) - rrt(249) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(250) * density(01) * density(30) * 2.0d0
  pd(01,30) = pd(01,30) - rrt(250) * density(01)**2 
  pd(23,01) = pd(23,01) + rrt(250) * density(01) * density(30) * 2.0d0
  pd(23,30) = pd(23,30) + rrt(250) * density(01)**2 
  pd(30,01) = pd(30,01) - rrt(250) * density(01) * density(30) * 2.0d0
  pd(30,30) = pd(30,30) - rrt(250) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(251) * density(01) * density(49) * 2.0d0
  pd(01,49) = pd(01,49) - rrt(251) * density(01)**2 
  pd(09,01) = pd(09,01) + rrt(251) * density(01) * density(49) * 2.0d0
  pd(09,49) = pd(09,49) + rrt(251) * density(01)**2 
  pd(49,01) = pd(49,01) - rrt(251) * density(01) * density(49) * 2.0d0
  pd(49,49) = pd(49,49) - rrt(251) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(252) * density(01) * density(08) * 2.0d0
  pd(01,08) = pd(01,08) - rrt(252) * density(01)**2 
  pd(03,01) = pd(03,01) + rrt(252) * density(01) * density(08) * 2.0d0
  pd(03,08) = pd(03,08) + rrt(252) * density(01)**2 
  pd(08,01) = pd(08,01) - rrt(252) * density(01) * density(08) * 2.0d0
  pd(08,08) = pd(08,08) - rrt(252) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(253) * density(01) * density(47) * 2.0d0
  pd(01,47) = pd(01,47) - rrt(253) * density(01)**2 
  pd(17,01) = pd(17,01) + rrt(253) * density(01) * density(47) * 2.0d0
  pd(17,47) = pd(17,47) + rrt(253) * density(01)**2 
  pd(47,01) = pd(47,01) - rrt(253) * density(01) * density(47) * 2.0d0
  pd(47,47) = pd(47,47) - rrt(253) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(254) * density(01) * density(07) * 2.0d0
  pd(01,07) = pd(01,07) - rrt(254) * density(01)**2 
  pd(07,01) = pd(07,01) - rrt(254) * density(01) * density(07) * 2.0d0
  pd(07,07) = pd(07,07) - rrt(254) * density(01)**2 
  pd(13,01) = pd(13,01) + rrt(254) * density(01) * density(07) * 2.0d0
  pd(13,07) = pd(13,07) + rrt(254) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(255) * density(01) * density(45) * 2.0d0
  pd(01,45) = pd(01,45) - rrt(255) * density(01)**2 
  pd(02,01) = pd(02,01) + rrt(255) * density(01) * density(45) * 2.0d0
  pd(02,45) = pd(02,45) + rrt(255) * density(01)**2 
  pd(45,01) = pd(45,01) - rrt(255) * density(01) * density(45) * 2.0d0
  pd(45,45) = pd(45,45) - rrt(255) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(256) * density(01) * density(41) * 2.0d0
  pd(01,41) = pd(01,41) - rrt(256) * density(01)**2 
  pd(13,01) = pd(13,01) + rrt(256) * density(01) * density(41) * 4.0d0
  pd(13,41) = pd(13,41) + rrt(256) * density(01)**2 * 2.0d0
  pd(24,01) = pd(24,01) + rrt(256) * density(01) * density(41) * 2.0d0
  pd(24,41) = pd(24,41) + rrt(256) * density(01)**2 
  pd(41,01) = pd(41,01) - rrt(256) * density(01) * density(41) * 2.0d0
  pd(41,41) = pd(41,41) - rrt(256) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(257) * density(01) * density(41) * 2.0d0
  pd(01,41) = pd(01,41) - rrt(257) * density(01)**2 
  pd(02,01) = pd(02,01) + rrt(257) * density(01) * density(41) * 2.0d0
  pd(02,41) = pd(02,41) + rrt(257) * density(01)**2 
  pd(10,01) = pd(10,01) + rrt(257) * density(01) * density(41) * 2.0d0
  pd(10,41) = pd(10,41) + rrt(257) * density(01)**2 
  pd(13,01) = pd(13,01) + rrt(257) * density(01) * density(41) * 2.0d0
  pd(13,41) = pd(13,41) + rrt(257) * density(01)**2 
  pd(41,01) = pd(41,01) - rrt(257) * density(01) * density(41) * 2.0d0
  pd(41,41) = pd(41,41) - rrt(257) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(258) * density(01) * density(25) * 2.0d0
  pd(01,25) = pd(01,25) - rrt(258) * density(01)**2 
  pd(13,01) = pd(13,01) + rrt(258) * density(01) * density(25) * 2.0d0
  pd(13,25) = pd(13,25) + rrt(258) * density(01)**2 
  pd(24,01) = pd(24,01) + rrt(258) * density(01) * density(25) * 2.0d0
  pd(24,25) = pd(24,25) + rrt(258) * density(01)**2 
  pd(25,01) = pd(25,01) - rrt(258) * density(01) * density(25) * 2.0d0
  pd(25,25) = pd(25,25) - rrt(258) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(259) * density(01) * density(25) * 2.0d0
  pd(01,25) = pd(01,25) - rrt(259) * density(01)**2 
  pd(10,01) = pd(10,01) + rrt(259) * density(01) * density(25) * 2.0d0
  pd(10,25) = pd(10,25) + rrt(259) * density(01)**2 
  pd(13,01) = pd(13,01) + rrt(259) * density(01) * density(25) * 4.0d0
  pd(13,25) = pd(13,25) + rrt(259) * density(01)**2 * 2.0d0
  pd(25,01) = pd(25,01) - rrt(259) * density(01) * density(25) * 2.0d0
  pd(25,25) = pd(25,25) - rrt(259) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(260) * density(01) * density(25) * 2.0d0
  pd(01,25) = pd(01,25) - rrt(260) * density(01)**2 
  pd(02,01) = pd(02,01) + rrt(260) * density(01) * density(25) * 2.0d0
  pd(02,25) = pd(02,25) + rrt(260) * density(01)**2 
  pd(13,01) = pd(13,01) + rrt(260) * density(01) * density(25) * 2.0d0
  pd(13,25) = pd(13,25) + rrt(260) * density(01)**2 
  pd(25,01) = pd(25,01) - rrt(260) * density(01) * density(25) * 2.0d0
  pd(25,25) = pd(25,25) - rrt(260) * density(01)**2 
  pd(46,01) = pd(46,01) + rrt(260) * density(01) * density(25) * 2.0d0
  pd(46,25) = pd(46,25) + rrt(260) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(261) * density(01) * density(48) * 2.0d0
  pd(01,48) = pd(01,48) - rrt(261) * density(01)**2 
  pd(10,01) = pd(10,01) + rrt(261) * density(01) * density(48) * 2.0d0
  pd(10,48) = pd(10,48) + rrt(261) * density(01)**2 
  pd(13,01) = pd(13,01) + rrt(261) * density(01) * density(48) * 2.0d0
  pd(13,48) = pd(13,48) + rrt(261) * density(01)**2 
  pd(48,01) = pd(48,01) - rrt(261) * density(01) * density(48) * 2.0d0
  pd(48,48) = pd(48,48) - rrt(261) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(262) * density(01) * density(48) * 2.0d0
  pd(01,48) = pd(01,48) - rrt(262) * density(01)**2 
  pd(02,01) = pd(02,01) + rrt(262) * density(01) * density(48) * 2.0d0
  pd(02,48) = pd(02,48) + rrt(262) * density(01)**2 
  pd(46,01) = pd(46,01) + rrt(262) * density(01) * density(48) * 2.0d0
  pd(46,48) = pd(46,48) + rrt(262) * density(01)**2 
  pd(48,01) = pd(48,01) - rrt(262) * density(01) * density(48) * 2.0d0
  pd(48,48) = pd(48,48) - rrt(262) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(263) * density(01) * density(27) * 2.0d0
  pd(01,27) = pd(01,27) - rrt(263) * density(01)**2 
  pd(13,01) = pd(13,01) + rrt(263) * density(01) * density(27) * 2.0d0
  pd(13,27) = pd(13,27) + rrt(263) * density(01)**2 
  pd(27,01) = pd(27,01) - rrt(263) * density(01) * density(27) * 2.0d0
  pd(27,27) = pd(27,27) - rrt(263) * density(01)**2 
  pd(46,01) = pd(46,01) + rrt(263) * density(01) * density(27) * 2.0d0
  pd(46,27) = pd(46,27) + rrt(263) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(264) * density(01) * density(15) * 2.0d0
  pd(01,15) = pd(01,15) - rrt(264) * density(01)**2 
  pd(13,01) = pd(13,01) + rrt(264) * density(01) * density(15) * 2.0d0
  pd(13,15) = pd(13,15) + rrt(264) * density(01)**2 
  pd(15,01) = pd(15,01) - rrt(264) * density(01) * density(15) * 2.0d0
  pd(15,15) = pd(15,15) - rrt(264) * density(01)**2 
  pd(38,01) = pd(38,01) + rrt(264) * density(01) * density(15) * 2.0d0
  pd(38,15) = pd(38,15) + rrt(264) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(265) * density(01) * density(15) * 2.0d0
  pd(01,15) = pd(01,15) - rrt(265) * density(01)**2 
  pd(13,01) = pd(13,01) + rrt(265) * density(01) * density(15) * 4.0d0
  pd(13,15) = pd(13,15) + rrt(265) * density(01)**2 * 2.0d0
  pd(15,01) = pd(15,01) - rrt(265) * density(01) * density(15) * 2.0d0
  pd(15,15) = pd(15,15) - rrt(265) * density(01)**2 
  pd(35,01) = pd(35,01) + rrt(265) * density(01) * density(15) * 2.0d0
  pd(35,15) = pd(35,15) + rrt(265) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(266) * density(01) * density(44) * 2.0d0
  pd(01,44) = pd(01,44) - rrt(266) * density(01)**2 
  pd(13,01) = pd(13,01) + rrt(266) * density(01) * density(44) * 2.0d0
  pd(13,44) = pd(13,44) + rrt(266) * density(01)**2 
  pd(35,01) = pd(35,01) + rrt(266) * density(01) * density(44) * 2.0d0
  pd(35,44) = pd(35,44) + rrt(266) * density(01)**2 
  pd(44,01) = pd(44,01) - rrt(266) * density(01) * density(44) * 2.0d0
  pd(44,44) = pd(44,44) - rrt(266) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(267) * density(01) * density(44) * 2.0d0
  pd(01,44) = pd(01,44) - rrt(267) * density(01)**2 
  pd(13,01) = pd(13,01) + rrt(267) * density(01) * density(44) * 4.0d0
  pd(13,44) = pd(13,44) + rrt(267) * density(01)**2 * 2.0d0
  pd(44,01) = pd(44,01) - rrt(267) * density(01) * density(44) * 2.0d0
  pd(44,44) = pd(44,44) - rrt(267) * density(01)**2 
  pd(50,01) = pd(50,01) + rrt(267) * density(01) * density(44) * 2.0d0
  pd(50,44) = pd(50,44) + rrt(267) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(268) * density(01) * density(44) * 2.0d0
  pd(01,44) = pd(01,44) - rrt(268) * density(01)**2 
  pd(02,01) = pd(02,01) + rrt(268) * density(01) * density(44) * 2.0d0
  pd(02,44) = pd(02,44) + rrt(268) * density(01)**2 
  pd(13,01) = pd(13,01) + rrt(268) * density(01) * density(44) * 2.0d0
  pd(13,44) = pd(13,44) + rrt(268) * density(01)**2 
  pd(19,01) = pd(19,01) + rrt(268) * density(01) * density(44) * 2.0d0
  pd(19,44) = pd(19,44) + rrt(268) * density(01)**2 
  pd(44,01) = pd(44,01) - rrt(268) * density(01) * density(44) * 2.0d0
  pd(44,44) = pd(44,44) - rrt(268) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(269) * density(01) * density(44) * 2.0d0
  pd(01,44) = pd(01,44) - rrt(269) * density(01)**2 
  pd(13,01) = pd(13,01) + rrt(269) * density(01) * density(44) * 6.0d0
  pd(13,44) = pd(13,44) + rrt(269) * density(01)**2 * 3.0d0
  pd(19,01) = pd(19,01) + rrt(269) * density(01) * density(44) * 2.0d0
  pd(19,44) = pd(19,44) + rrt(269) * density(01)**2 
  pd(44,01) = pd(44,01) - rrt(269) * density(01) * density(44) * 2.0d0
  pd(44,44) = pd(44,44) - rrt(269) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(270) * density(01) * density(44) * 2.0d0
  pd(01,44) = pd(01,44) - rrt(270) * density(01)**2 
  pd(10,01) = pd(10,01) + rrt(270) * density(01) * density(44) * 2.0d0
  pd(10,44) = pd(10,44) + rrt(270) * density(01)**2 
  pd(24,01) = pd(24,01) + rrt(270) * density(01) * density(44) * 2.0d0
  pd(24,44) = pd(24,44) + rrt(270) * density(01)**2 
  pd(44,01) = pd(44,01) - rrt(270) * density(01) * density(44) * 2.0d0
  pd(44,44) = pd(44,44) - rrt(270) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(271) * density(01) * density(32) * 2.0d0
  pd(01,32) = pd(01,32) - rrt(271) * density(01)**2 
  pd(13,01) = pd(13,01) + rrt(271) * density(01) * density(32) * 2.0d0
  pd(13,32) = pd(13,32) + rrt(271) * density(01)**2 
  pd(32,01) = pd(32,01) - rrt(271) * density(01) * density(32) * 2.0d0
  pd(32,32) = pd(32,32) - rrt(271) * density(01)**2 
  pd(50,01) = pd(50,01) + rrt(271) * density(01) * density(32) * 2.0d0
  pd(50,32) = pd(50,32) + rrt(271) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(272) * density(01) * density(32) * 2.0d0
  pd(01,32) = pd(01,32) - rrt(272) * density(01)**2 
  pd(13,01) = pd(13,01) + rrt(272) * density(01) * density(32) * 4.0d0
  pd(13,32) = pd(13,32) + rrt(272) * density(01)**2 * 2.0d0
  pd(19,01) = pd(19,01) + rrt(272) * density(01) * density(32) * 2.0d0
  pd(19,32) = pd(19,32) + rrt(272) * density(01)**2 
  pd(32,01) = pd(32,01) - rrt(272) * density(01) * density(32) * 2.0d0
  pd(32,32) = pd(32,32) - rrt(272) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(273) * density(01) * density(06) * 2.0d0
  pd(01,06) = pd(01,06) - rrt(273) * density(01)**2 
  pd(06,01) = pd(06,01) - rrt(273) * density(01) * density(06) * 2.0d0
  pd(06,06) = pd(06,06) - rrt(273) * density(01)**2 
  pd(13,01) = pd(13,01) + rrt(273) * density(01) * density(06) * 2.0d0
  pd(13,06) = pd(13,06) + rrt(273) * density(01)**2 
  pd(19,01) = pd(19,01) + rrt(273) * density(01) * density(06) * 2.0d0
  pd(19,06) = pd(19,06) + rrt(273) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(274) * density(01) * density(14) * 2.0d0
  pd(01,14) = pd(01,14) - rrt(274) * density(01)**2 
  pd(14,01) = pd(14,01) - rrt(274) * density(01) * density(14) * 2.0d0
  pd(14,14) = pd(14,14) - rrt(274) * density(01)**2 
  pd(46,01) = pd(46,01) + rrt(274) * density(01) * density(14) * 4.0d0
  pd(46,14) = pd(46,14) + rrt(274) * density(01)**2 * 2.0d0
  pd(01,01) = pd(01,01) - rrt(275) * density(01) * density(45) * 2.0d0
  pd(01,45) = pd(01,45) - rrt(275) * density(01)**2 
  pd(13,01) = pd(13,01) + rrt(275) * density(01) * density(45) * 4.0d0
  pd(13,45) = pd(13,45) + rrt(275) * density(01)**2 * 2.0d0
  pd(45,01) = pd(45,01) - rrt(275) * density(01) * density(45) * 2.0d0
  pd(45,45) = pd(45,45) - rrt(275) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(276) * density(01) * density(26) * 2.0d0
  pd(01,26) = pd(01,26) - rrt(276) * density(01)**2 
  pd(02,01) = pd(02,01) + rrt(276) * density(01) * density(26) * 2.0d0
  pd(02,26) = pd(02,26) + rrt(276) * density(01)**2 
  pd(13,01) = pd(13,01) + rrt(276) * density(01) * density(26) * 2.0d0
  pd(13,26) = pd(13,26) + rrt(276) * density(01)**2 
  pd(26,01) = pd(26,01) - rrt(276) * density(01) * density(26) * 2.0d0
  pd(26,26) = pd(26,26) - rrt(276) * density(01)**2 
  pd(02,02) = pd(02,02) - rrt(277) * density(45) 
  pd(02,45) = pd(02,45) - rrt(277) * density(02) 
  pd(13,02) = pd(13,02) + rrt(277) * density(45) 
  pd(13,45) = pd(13,45) + rrt(277) * density(02) 
  pd(26,02) = pd(26,02) + rrt(277) * density(45) 
  pd(26,45) = pd(26,45) + rrt(277) * density(02) 
  pd(45,02) = pd(45,02) - rrt(277) * density(45) 
  pd(45,45) = pd(45,45) - rrt(277) * density(02) 
  pd(05,10) = pd(05,10) + rrt(278) * density(41) 
  pd(05,41) = pd(05,41) + rrt(278) * density(10) 
  pd(10,10) = pd(10,10) - rrt(278) * density(41) 
  pd(10,41) = pd(10,41) - rrt(278) * density(10) 
  pd(41,10) = pd(41,10) - rrt(278) * density(41) 
  pd(41,41) = pd(41,41) - rrt(278) * density(10) 
  pd(48,10) = pd(48,10) + rrt(278) * density(41) 
  pd(48,41) = pd(48,41) + rrt(278) * density(10) 
  pd(05,41) = pd(05,41) + rrt(279) * density(46) 
  pd(05,46) = pd(05,46) + rrt(279) * density(41) 
  pd(27,41) = pd(27,41) + rrt(279) * density(46) 
  pd(27,46) = pd(27,46) + rrt(279) * density(41) 
  pd(41,41) = pd(41,41) - rrt(279) * density(46) 
  pd(41,46) = pd(41,46) - rrt(279) * density(41) 
  pd(46,41) = pd(46,41) - rrt(279) * density(46) 
  pd(46,46) = pd(46,46) - rrt(279) * density(41) 
  pd(02,20) = pd(02,20) + rrt(280) * density(41) 
  pd(02,41) = pd(02,41) + rrt(280) * density(20) 
  pd(05,20) = pd(05,20) + rrt(280) * density(41) 
  pd(05,41) = pd(05,41) + rrt(280) * density(20) 
  pd(20,20) = pd(20,20) - rrt(280) * density(41) 
  pd(20,41) = pd(20,41) - rrt(280) * density(20) 
  pd(41,20) = pd(41,20) - rrt(280) * density(41) 
  pd(41,41) = pd(41,41) - rrt(280) * density(20) 
  pd(44,20) = pd(44,20) + rrt(280) * density(41) 
  pd(44,41) = pd(44,41) + rrt(280) * density(20) 
  pd(05,35) = pd(05,35) + rrt(281) * density(41) 
  pd(05,41) = pd(05,41) + rrt(281) * density(35) 
  pd(35,35) = pd(35,35) - rrt(281) * density(41) 
  pd(35,41) = pd(35,41) - rrt(281) * density(35) 
  pd(41,35) = pd(41,35) - rrt(281) * density(41) 
  pd(41,41) = pd(41,41) - rrt(281) * density(35) 
  pd(44,35) = pd(44,35) + rrt(281) * density(41) 
  pd(44,41) = pd(44,41) + rrt(281) * density(35) 
  pd(05,19) = pd(05,19) + rrt(282) * density(41) 
  pd(05,41) = pd(05,41) + rrt(282) * density(19) 
  pd(06,19) = pd(06,19) + rrt(282) * density(41) 
  pd(06,41) = pd(06,41) + rrt(282) * density(19) 
  pd(19,19) = pd(19,19) - rrt(282) * density(41) 
  pd(19,41) = pd(19,41) - rrt(282) * density(19) 
  pd(41,19) = pd(41,19) - rrt(282) * density(41) 
  pd(41,41) = pd(41,41) - rrt(282) * density(19) 
  pd(02,13) = pd(02,13) + rrt(283) * density(41) 
  pd(02,41) = pd(02,41) + rrt(283) * density(13) 
  pd(13,13) = pd(13,13) - rrt(283) * density(41) 
  pd(13,41) = pd(13,41) - rrt(283) * density(13) 
  pd(25,13) = pd(25,13) + rrt(283) * density(41) 
  pd(25,41) = pd(25,41) + rrt(283) * density(13) 
  pd(41,13) = pd(41,13) - rrt(283) * density(41) 
  pd(41,41) = pd(41,41) - rrt(283) * density(13) 
  pd(05,05) = pd(05,05) - rrt(284) * density(25) 
  pd(05,25) = pd(05,25) - rrt(284) * density(05) 
  pd(24,05) = pd(24,05) + rrt(284) * density(25) 
  pd(24,25) = pd(24,25) + rrt(284) * density(05) 
  pd(25,05) = pd(25,05) - rrt(284) * density(25) 
  pd(25,25) = pd(25,25) - rrt(284) * density(05) 
  pd(41,05) = pd(41,05) + rrt(284) * density(25) 
  pd(41,25) = pd(41,25) + rrt(284) * density(05) 
  pd(02,20) = pd(02,20) + rrt(285) * density(25) 
  pd(02,25) = pd(02,25) + rrt(285) * density(20) 
  pd(05,20) = pd(05,20) + rrt(285) * density(25) 
  pd(05,25) = pd(05,25) + rrt(285) * density(20) 
  pd(20,20) = pd(20,20) - rrt(285) * density(25) 
  pd(20,25) = pd(20,25) - rrt(285) * density(20) 
  pd(25,20) = pd(25,20) - rrt(285) * density(25) 
  pd(25,25) = pd(25,25) - rrt(285) * density(20) 
  pd(32,20) = pd(32,20) + rrt(285) * density(25) 
  pd(32,25) = pd(32,25) + rrt(285) * density(20) 
  pd(24,25) = pd(24,25) + rrt(286) * density(35) 
  pd(24,35) = pd(24,35) + rrt(286) * density(25) 
  pd(25,25) = pd(25,25) - rrt(286) * density(35) 
  pd(25,35) = pd(25,35) - rrt(286) * density(25) 
  pd(35,25) = pd(35,25) - rrt(286) * density(35) 
  pd(35,35) = pd(35,35) - rrt(286) * density(25) 
  pd(44,25) = pd(44,25) + rrt(286) * density(35) 
  pd(44,35) = pd(44,35) + rrt(286) * density(25) 
  pd(05,25) = pd(05,25) + rrt(287) * density(35) 
  pd(05,35) = pd(05,35) + rrt(287) * density(25) 
  pd(25,25) = pd(25,25) - rrt(287) * density(35) 
  pd(25,35) = pd(25,35) - rrt(287) * density(25) 
  pd(32,25) = pd(32,25) + rrt(287) * density(35) 
  pd(32,35) = pd(32,35) + rrt(287) * density(25) 
  pd(35,25) = pd(35,25) - rrt(287) * density(35) 
  pd(35,35) = pd(35,35) - rrt(287) * density(25) 
  pd(06,19) = pd(06,19) + rrt(288) * density(25) 
  pd(06,25) = pd(06,25) + rrt(288) * density(19) 
  pd(19,19) = pd(19,19) - rrt(288) * density(25) 
  pd(19,25) = pd(19,25) - rrt(288) * density(19) 
  pd(24,19) = pd(24,19) + rrt(288) * density(25) 
  pd(24,25) = pd(24,25) + rrt(288) * density(19) 
  pd(25,19) = pd(25,19) - rrt(288) * density(25) 
  pd(25,25) = pd(25,25) - rrt(288) * density(19) 
  pd(05,19) = pd(05,19) + rrt(289) * density(25) 
  pd(05,25) = pd(05,25) + rrt(289) * density(19) 
  pd(14,19) = pd(14,19) + rrt(289) * density(25) 
  pd(14,25) = pd(14,25) + rrt(289) * density(19) 
  pd(19,19) = pd(19,19) - rrt(289) * density(25) 
  pd(19,25) = pd(19,25) - rrt(289) * density(19) 
  pd(25,19) = pd(25,19) - rrt(289) * density(25) 
  pd(25,25) = pd(25,25) - rrt(289) * density(19) 
  pd(02,02) = pd(02,02) - rrt(290) * density(25) 
  pd(02,25) = pd(02,25) - rrt(290) * density(02) 
  pd(13,02) = pd(13,02) + rrt(290) * density(25) 
  pd(13,25) = pd(13,25) + rrt(290) * density(02) 
  pd(25,02) = pd(25,02) - rrt(290) * density(25) 
  pd(25,25) = pd(25,25) - rrt(290) * density(02) 
  pd(41,02) = pd(41,02) + rrt(290) * density(25) 
  pd(41,25) = pd(41,25) + rrt(290) * density(02) 
  pd(02,13) = pd(02,13) + rrt(291) * density(25) 
  pd(02,25) = pd(02,25) + rrt(291) * density(13) 
  pd(13,13) = pd(13,13) - rrt(291) * density(25) 
  pd(13,25) = pd(13,25) - rrt(291) * density(13) 
  pd(25,13) = pd(25,13) - rrt(291) * density(25) 
  pd(25,25) = pd(25,25) - rrt(291) * density(13) 
  pd(48,13) = pd(48,13) + rrt(291) * density(25) 
  pd(48,25) = pd(48,25) + rrt(291) * density(13) 
  pd(05,05) = pd(05,05) - rrt(292) * density(48) 
  pd(05,48) = pd(05,48) - rrt(292) * density(05) 
  pd(24,05) = pd(24,05) + rrt(292) * density(48) 
  pd(24,48) = pd(24,48) + rrt(292) * density(05) 
  pd(25,05) = pd(25,05) + rrt(292) * density(48) 
  pd(25,48) = pd(25,48) + rrt(292) * density(05) 
  pd(48,05) = pd(48,05) - rrt(292) * density(48) 
  pd(48,48) = pd(48,48) - rrt(292) * density(05) 
  pd(02,05) = pd(02,05) + rrt(293) * density(48) 
  pd(02,48) = pd(02,48) + rrt(293) * density(05) 
  pd(05,05) = pd(05,05) - rrt(293) * density(48) 
  pd(05,48) = pd(05,48) - rrt(293) * density(05) 
  pd(44,05) = pd(44,05) + rrt(293) * density(48) 
  pd(44,48) = pd(44,48) + rrt(293) * density(05) 
  pd(48,05) = pd(48,05) - rrt(293) * density(48) 
  pd(48,48) = pd(48,48) - rrt(293) * density(05) 
  pd(02,10) = pd(02,10) + rrt(294) * density(48) 
  pd(02,48) = pd(02,48) + rrt(294) * density(10) 
  pd(06,10) = pd(06,10) + rrt(294) * density(48) 
  pd(06,48) = pd(06,48) + rrt(294) * density(10) 
  pd(10,10) = pd(10,10) - rrt(294) * density(48) 
  pd(10,48) = pd(10,48) - rrt(294) * density(10) 
  pd(48,10) = pd(48,10) - rrt(294) * density(48) 
  pd(48,48) = pd(48,48) - rrt(294) * density(10) 
  pd(02,46) = pd(02,46) + rrt(295) * density(48) 
  pd(02,48) = pd(02,48) + rrt(295) * density(46) 
  pd(14,46) = pd(14,46) + rrt(295) * density(48) 
  pd(14,48) = pd(14,48) + rrt(295) * density(46) 
  pd(46,46) = pd(46,46) - rrt(295) * density(48) 
  pd(46,48) = pd(46,48) - rrt(295) * density(46) 
  pd(48,46) = pd(48,46) - rrt(295) * density(48) 
  pd(48,48) = pd(48,48) - rrt(295) * density(46) 
  pd(05,20) = pd(05,20) + rrt(296) * density(48) 
  pd(05,48) = pd(05,48) + rrt(296) * density(20) 
  pd(20,20) = pd(20,20) - rrt(296) * density(48) 
  pd(20,48) = pd(20,48) - rrt(296) * density(20) 
  pd(44,20) = pd(44,20) + rrt(296) * density(48) 
  pd(44,48) = pd(44,48) + rrt(296) * density(20) 
  pd(48,20) = pd(48,20) - rrt(296) * density(48) 
  pd(48,48) = pd(48,48) - rrt(296) * density(20) 
  pd(05,35) = pd(05,35) + rrt(297) * density(48) 
  pd(05,48) = pd(05,48) + rrt(297) * density(35) 
  pd(06,35) = pd(06,35) + rrt(297) * density(48) 
  pd(06,48) = pd(06,48) + rrt(297) * density(35) 
  pd(35,35) = pd(35,35) - rrt(297) * density(48) 
  pd(35,48) = pd(35,48) - rrt(297) * density(35) 
  pd(48,35) = pd(48,35) - rrt(297) * density(48) 
  pd(48,48) = pd(48,48) - rrt(297) * density(35) 
  pd(06,48) = pd(06,48) + rrt(298) * density(50) 
  pd(06,50) = pd(06,50) + rrt(298) * density(48) 
  pd(24,48) = pd(24,48) + rrt(298) * density(50) 
  pd(24,50) = pd(24,50) + rrt(298) * density(48) 
  pd(48,48) = pd(48,48) - rrt(298) * density(50) 
  pd(48,50) = pd(48,50) - rrt(298) * density(48) 
  pd(50,48) = pd(50,48) - rrt(298) * density(50) 
  pd(50,50) = pd(50,50) - rrt(298) * density(48) 
  pd(05,05) = pd(05,05) - rrt(299) * density(27) 
  pd(05,27) = pd(05,27) - rrt(299) * density(05) 
  pd(24,05) = pd(24,05) + rrt(299) * density(27) 
  pd(24,27) = pd(24,27) + rrt(299) * density(05) 
  pd(27,05) = pd(27,05) - rrt(299) * density(27) 
  pd(27,27) = pd(27,27) - rrt(299) * density(05) 
  pd(48,05) = pd(48,05) + rrt(299) * density(27) 
  pd(48,27) = pd(48,27) + rrt(299) * density(05) 
  pd(05,05) = pd(05,05) - rrt(300) * density(27) 
  pd(05,27) = pd(05,27) - rrt(300) * density(05) 
  pd(13,05) = pd(13,05) + rrt(300) * density(27) 
  pd(13,27) = pd(13,27) + rrt(300) * density(05) 
  pd(27,05) = pd(27,05) - rrt(300) * density(27) 
  pd(27,27) = pd(27,27) - rrt(300) * density(05) 
  pd(44,05) = pd(44,05) + rrt(300) * density(27) 
  pd(44,27) = pd(44,27) + rrt(300) * density(05) 
  pd(02,05) = pd(02,05) + rrt(301) * density(27) 
  pd(02,27) = pd(02,27) + rrt(301) * density(05) 
  pd(05,05) = pd(05,05) - rrt(301) * density(27) 
  pd(05,27) = pd(05,27) - rrt(301) * density(05) 
  pd(27,05) = pd(27,05) - rrt(301) * density(27) 
  pd(27,27) = pd(27,27) - rrt(301) * density(05) 
  pd(32,05) = pd(32,05) + rrt(301) * density(27) 
  pd(32,27) = pd(32,27) + rrt(301) * density(05) 
  pd(02,05) = pd(02,05) + rrt(302) * density(27) 
  pd(02,27) = pd(02,27) + rrt(302) * density(05) 
  pd(05,05) = pd(05,05) - rrt(302) * density(27) 
  pd(05,27) = pd(05,27) - rrt(302) * density(05) 
  pd(06,05) = pd(06,05) + rrt(302) * density(27) 
  pd(06,27) = pd(06,27) + rrt(302) * density(05) 
  pd(13,05) = pd(13,05) + rrt(302) * density(27) 
  pd(13,27) = pd(13,27) + rrt(302) * density(05) 
  pd(27,05) = pd(27,05) - rrt(302) * density(27) 
  pd(27,27) = pd(27,27) - rrt(302) * density(05) 
  pd(02,05) = pd(02,05) + rrt(303) * density(27) * 2.0d0
  pd(02,27) = pd(02,27) + rrt(303) * density(05) * 2.0d0
  pd(05,05) = pd(05,05) - rrt(303) * density(27) 
  pd(05,27) = pd(05,27) - rrt(303) * density(05) 
  pd(14,05) = pd(14,05) + rrt(303) * density(27) 
  pd(14,27) = pd(14,27) + rrt(303) * density(05) 
  pd(27,05) = pd(27,05) - rrt(303) * density(27) 
  pd(27,27) = pd(27,27) - rrt(303) * density(05) 
  pd(02,02) = pd(02,02) - rrt(304) * density(27) 
  pd(02,27) = pd(02,27) - rrt(304) * density(02) 
  pd(13,02) = pd(13,02) + rrt(304) * density(27) 
  pd(13,27) = pd(13,27) + rrt(304) * density(02) 
  pd(27,02) = pd(27,02) - rrt(304) * density(27) 
  pd(27,27) = pd(27,27) - rrt(304) * density(02) 
  pd(48,02) = pd(48,02) + rrt(304) * density(27) 
  pd(48,27) = pd(48,27) + rrt(304) * density(02) 
  pd(05,05) = pd(05,05) - rrt(305) * density(28) 
  pd(05,28) = pd(05,28) - rrt(305) * density(05) 
  pd(13,05) = pd(13,05) + rrt(305) * density(28) 
  pd(13,28) = pd(13,28) + rrt(305) * density(05) 
  pd(28,05) = pd(28,05) - rrt(305) * density(28) 
  pd(28,28) = pd(28,28) - rrt(305) * density(05) 
  pd(32,05) = pd(32,05) + rrt(305) * density(28) 
  pd(32,28) = pd(32,28) + rrt(305) * density(05) 
  pd(02,05) = pd(02,05) + rrt(306) * density(28) 
  pd(02,28) = pd(02,28) + rrt(306) * density(05) 
  pd(05,05) = pd(05,05) - rrt(306) * density(28) 
  pd(05,28) = pd(05,28) - rrt(306) * density(05) 
  pd(06,05) = pd(06,05) + rrt(306) * density(28) 
  pd(06,28) = pd(06,28) + rrt(306) * density(05) 
  pd(28,05) = pd(28,05) - rrt(306) * density(28) 
  pd(28,28) = pd(28,28) - rrt(306) * density(05) 
  pd(02,05) = pd(02,05) + rrt(307) * density(28) 
  pd(02,28) = pd(02,28) + rrt(307) * density(05) 
  pd(05,05) = pd(05,05) - rrt(307) * density(28) 
  pd(05,28) = pd(05,28) - rrt(307) * density(05) 
  pd(13,05) = pd(13,05) + rrt(307) * density(28) 
  pd(13,28) = pd(13,28) + rrt(307) * density(05) 
  pd(14,05) = pd(14,05) + rrt(307) * density(28) 
  pd(14,28) = pd(14,28) + rrt(307) * density(05) 
  pd(28,05) = pd(28,05) - rrt(307) * density(28) 
  pd(28,28) = pd(28,28) - rrt(307) * density(05) 
  pd(02,02) = pd(02,02) - rrt(308) * density(28) 
  pd(02,28) = pd(02,28) - rrt(308) * density(02) 
  pd(13,02) = pd(13,02) + rrt(308) * density(28) 
  pd(13,28) = pd(13,28) + rrt(308) * density(02) 
  pd(27,02) = pd(27,02) + rrt(308) * density(28) 
  pd(27,28) = pd(27,28) + rrt(308) * density(02) 
  pd(28,02) = pd(28,02) - rrt(308) * density(28) 
  pd(28,28) = pd(28,28) - rrt(308) * density(02) 
  pd(15,15) = pd(15,15) - rrt(309) * density(35) 
  pd(15,35) = pd(15,35) - rrt(309) * density(15) 
  pd(20,15) = pd(20,15) + rrt(309) * density(35) 
  pd(20,35) = pd(20,35) + rrt(309) * density(15) 
  pd(32,15) = pd(32,15) + rrt(309) * density(35) 
  pd(32,35) = pd(32,35) + rrt(309) * density(15) 
  pd(35,15) = pd(35,15) - rrt(309) * density(35) 
  pd(35,35) = pd(35,35) - rrt(309) * density(15) 
  pd(15,15) = pd(15,15) - rrt(310) * density(19) 
  pd(15,19) = pd(15,19) - rrt(310) * density(15) 
  pd(19,15) = pd(19,15) - rrt(310) * density(19) 
  pd(19,19) = pd(19,19) - rrt(310) * density(15) 
  pd(44,15) = pd(44,15) + rrt(310) * density(19) 
  pd(44,19) = pd(44,19) + rrt(310) * density(15) 
  pd(50,15) = pd(50,15) + rrt(310) * density(19) 
  pd(50,19) = pd(50,19) + rrt(310) * density(15) 
  pd(02,13) = pd(02,13) + rrt(311) * density(15) 
  pd(02,15) = pd(02,15) + rrt(311) * density(13) 
  pd(13,13) = pd(13,13) - rrt(311) * density(15) 
  pd(13,15) = pd(13,15) - rrt(311) * density(13) 
  pd(15,13) = pd(15,13) - rrt(311) * density(15) 
  pd(15,15) = pd(15,15) - rrt(311) * density(13) 
  pd(44,13) = pd(44,13) + rrt(311) * density(15) 
  pd(44,15) = pd(44,15) + rrt(311) * density(13) 
  pd(02,13) = pd(02,13) + rrt(312) * density(44) 
  pd(02,44) = pd(02,44) + rrt(312) * density(13) 
  pd(13,13) = pd(13,13) - rrt(312) * density(44) 
  pd(13,44) = pd(13,44) - rrt(312) * density(13) 
  pd(32,13) = pd(32,13) + rrt(312) * density(44) 
  pd(32,44) = pd(32,44) + rrt(312) * density(13) 
  pd(44,13) = pd(44,13) - rrt(312) * density(44) 
  pd(44,44) = pd(44,44) - rrt(312) * density(13) 
  pd(19,32) = pd(19,32) + rrt(313) * density(50) 
  pd(19,50) = pd(19,50) + rrt(313) * density(32) 
  pd(32,32) = pd(32,32) - rrt(313) * density(50) 
  pd(32,50) = pd(32,50) - rrt(313) * density(32) 
  pd(44,32) = pd(44,32) + rrt(313) * density(50) 
  pd(44,50) = pd(44,50) + rrt(313) * density(32) 
  pd(50,32) = pd(50,32) - rrt(313) * density(50) 
  pd(50,50) = pd(50,50) - rrt(313) * density(32) 
  pd(06,32) = pd(06,32) + rrt(314) * density(50) 
  pd(06,50) = pd(06,50) + rrt(314) * density(32) 
  pd(32,32) = pd(32,32) - rrt(314) * density(50) 
  pd(32,50) = pd(32,50) - rrt(314) * density(32) 
  pd(35,32) = pd(35,32) + rrt(314) * density(50) 
  pd(35,50) = pd(35,50) + rrt(314) * density(32) 
  pd(50,32) = pd(50,32) - rrt(314) * density(50) 
  pd(50,50) = pd(50,50) - rrt(314) * density(32) 
  pd(02,13) = pd(02,13) + rrt(315) * density(32) 
  pd(02,32) = pd(02,32) + rrt(315) * density(13) 
  pd(06,13) = pd(06,13) + rrt(315) * density(32) 
  pd(06,32) = pd(06,32) + rrt(315) * density(13) 
  pd(13,13) = pd(13,13) - rrt(315) * density(32) 
  pd(13,32) = pd(13,32) - rrt(315) * density(13) 
  pd(32,13) = pd(32,13) - rrt(315) * density(32) 
  pd(32,32) = pd(32,32) - rrt(315) * density(13) 
  pd(06,06) = pd(06,06) - rrt(316) * density(20) 
  pd(06,20) = pd(06,20) - rrt(316) * density(06) 
  pd(20,06) = pd(20,06) - rrt(316) * density(20) 
  pd(20,20) = pd(20,20) - rrt(316) * density(06) 
  pd(35,06) = pd(35,06) + rrt(316) * density(20) 
  pd(35,20) = pd(35,20) + rrt(316) * density(06) 
  pd(44,06) = pd(44,06) + rrt(316) * density(20) 
  pd(44,20) = pd(44,20) + rrt(316) * density(06) 
  pd(06,06) = pd(06,06) - rrt(317) * density(35) 
  pd(06,35) = pd(06,35) - rrt(317) * density(06) 
  pd(19,06) = pd(19,06) + rrt(317) * density(35) 
  pd(19,35) = pd(19,35) + rrt(317) * density(06) 
  pd(35,06) = pd(35,06) - rrt(317) * density(35) 
  pd(35,35) = pd(35,35) - rrt(317) * density(06) 
  pd(44,06) = pd(44,06) + rrt(317) * density(35) 
  pd(44,35) = pd(44,35) + rrt(317) * density(06) 
  pd(06,06) = pd(06,06) - rrt(318) * density(18) 
  pd(06,18) = pd(06,18) - rrt(318) * density(06) 
  pd(14,06) = pd(14,06) + rrt(318) * density(18) 
  pd(14,18) = pd(14,18) + rrt(318) * density(06) 
  pd(18,06) = pd(18,06) - rrt(318) * density(18) 
  pd(18,18) = pd(18,18) - rrt(318) * density(06) 
  pd(19,06) = pd(19,06) + rrt(318) * density(18) 
  pd(19,18) = pd(19,18) + rrt(318) * density(06) 
  pd(02,06) = pd(02,06) + rrt(319) * density(13) 
  pd(02,13) = pd(02,13) + rrt(319) * density(06) 
  pd(06,06) = pd(06,06) - rrt(319) * density(13) 
  pd(06,13) = pd(06,13) - rrt(319) * density(06) 
  pd(13,06) = pd(13,06) - rrt(319) * density(13) 
  pd(13,13) = pd(13,13) - rrt(319) * density(06) 
  pd(14,06) = pd(14,06) + rrt(319) * density(13) 
  pd(14,13) = pd(14,13) + rrt(319) * density(06) 
  pd(05,05) = pd(05,05) - rrt(320) * density(14) 
  pd(05,14) = pd(05,14) - rrt(320) * density(05) 
  pd(06,05) = pd(06,05) + rrt(320) * density(14) 
  pd(06,14) = pd(06,14) + rrt(320) * density(05) 
  pd(14,05) = pd(14,05) - rrt(320) * density(14) 
  pd(14,14) = pd(14,14) - rrt(320) * density(05) 
  pd(24,05) = pd(24,05) + rrt(320) * density(14) 
  pd(24,14) = pd(24,14) + rrt(320) * density(05) 
  pd(14,14) = pd(14,14) - rrt(321) * density(20) 
  pd(14,20) = pd(14,20) - rrt(321) * density(14) 
  pd(20,14) = pd(20,14) - rrt(321) * density(20) 
  pd(20,20) = pd(20,20) - rrt(321) * density(14) 
  pd(44,14) = pd(44,14) + rrt(321) * density(20) 
  pd(44,20) = pd(44,20) + rrt(321) * density(14) 
  pd(50,14) = pd(50,14) + rrt(321) * density(20) 
  pd(50,20) = pd(50,20) + rrt(321) * density(14) 
  pd(14,14) = pd(14,14) - rrt(322) * density(20) 
  pd(14,20) = pd(14,20) - rrt(322) * density(14) 
  pd(20,14) = pd(20,14) - rrt(322) * density(20) 
  pd(20,20) = pd(20,20) - rrt(322) * density(14) 
  pd(32,14) = pd(32,14) + rrt(322) * density(20) 
  pd(32,20) = pd(32,20) + rrt(322) * density(14) 
  pd(35,14) = pd(35,14) + rrt(322) * density(20) 
  pd(35,20) = pd(35,20) + rrt(322) * density(14) 
  pd(14,14) = pd(14,14) - rrt(323) * density(35) 
  pd(14,35) = pd(14,35) - rrt(323) * density(14) 
  pd(19,14) = pd(19,14) + rrt(323) * density(35) 
  pd(19,35) = pd(19,35) + rrt(323) * density(14) 
  pd(32,14) = pd(32,14) + rrt(323) * density(35) 
  pd(32,35) = pd(32,35) + rrt(323) * density(14) 
  pd(35,14) = pd(35,14) - rrt(323) * density(35) 
  pd(35,35) = pd(35,35) - rrt(323) * density(14) 
  pd(06,14) = pd(06,14) + rrt(324) * density(50) 
  pd(06,50) = pd(06,50) + rrt(324) * density(14) 
  pd(14,14) = pd(14,14) - rrt(324) * density(50) 
  pd(14,50) = pd(14,50) - rrt(324) * density(14) 
  pd(19,14) = pd(19,14) + rrt(324) * density(50) 
  pd(19,50) = pd(19,50) + rrt(324) * density(14) 
  pd(50,14) = pd(50,14) - rrt(324) * density(50) 
  pd(50,50) = pd(50,50) - rrt(324) * density(14) 
  pd(02,02) = pd(02,02) - rrt(325) * density(14) 
  pd(02,14) = pd(02,14) - rrt(325) * density(02) 
  pd(06,02) = pd(06,02) + rrt(325) * density(14) 
  pd(06,14) = pd(06,14) + rrt(325) * density(02) 
  pd(13,02) = pd(13,02) + rrt(325) * density(14) 
  pd(13,14) = pd(13,14) + rrt(325) * density(02) 
  pd(14,02) = pd(14,02) - rrt(325) * density(14) 
  pd(14,14) = pd(14,14) - rrt(325) * density(02) 
  pd(02,05) = pd(02,05) + rrt(326) * density(26) 
  pd(02,26) = pd(02,26) + rrt(326) * density(05) 
  pd(05,05) = pd(05,05) - rrt(326) * density(26) 
  pd(05,26) = pd(05,26) - rrt(326) * density(05) 
  pd(26,05) = pd(26,05) - rrt(326) * density(26) 
  pd(26,26) = pd(26,26) - rrt(326) * density(05) 
  pd(41,05) = pd(41,05) + rrt(326) * density(26) 
  pd(41,26) = pd(41,26) + rrt(326) * density(05) 
  pd(02,24) = pd(02,24) + rrt(327) * density(26) 
  pd(02,26) = pd(02,26) + rrt(327) * density(24) 
  pd(24,24) = pd(24,24) - rrt(327) * density(26) 
  pd(24,26) = pd(24,26) - rrt(327) * density(24) 
  pd(25,24) = pd(25,24) + rrt(327) * density(26) 
  pd(25,26) = pd(25,26) + rrt(327) * density(24) 
  pd(26,24) = pd(26,24) - rrt(327) * density(26) 
  pd(26,26) = pd(26,26) - rrt(327) * density(24) 
  pd(02,10) = pd(02,10) + rrt(328) * density(26) 
  pd(02,26) = pd(02,26) + rrt(328) * density(10) 
  pd(10,10) = pd(10,10) - rrt(328) * density(26) 
  pd(10,26) = pd(10,26) - rrt(328) * density(10) 
  pd(26,10) = pd(26,10) - rrt(328) * density(26) 
  pd(26,26) = pd(26,26) - rrt(328) * density(10) 
  pd(48,10) = pd(48,10) + rrt(328) * density(26) 
  pd(48,26) = pd(48,26) + rrt(328) * density(10) 
  pd(02,26) = pd(02,26) + rrt(329) * density(46) 
  pd(02,46) = pd(02,46) + rrt(329) * density(26) 
  pd(26,26) = pd(26,26) - rrt(329) * density(46) 
  pd(26,46) = pd(26,46) - rrt(329) * density(26) 
  pd(27,26) = pd(27,26) + rrt(329) * density(46) 
  pd(27,46) = pd(27,46) + rrt(329) * density(26) 
  pd(46,26) = pd(46,26) - rrt(329) * density(46) 
  pd(46,46) = pd(46,46) - rrt(329) * density(26) 
  pd(02,20) = pd(02,20) + rrt(330) * density(26) * 2.0d0
  pd(02,26) = pd(02,26) + rrt(330) * density(20) * 2.0d0
  pd(20,20) = pd(20,20) - rrt(330) * density(26) 
  pd(20,26) = pd(20,26) - rrt(330) * density(20) 
  pd(26,20) = pd(26,20) - rrt(330) * density(26) 
  pd(26,26) = pd(26,26) - rrt(330) * density(20) 
  pd(44,20) = pd(44,20) + rrt(330) * density(26) 
  pd(44,26) = pd(44,26) + rrt(330) * density(20) 
  pd(02,26) = pd(02,26) + rrt(331) * density(38) 
  pd(02,38) = pd(02,38) + rrt(331) * density(26) 
  pd(15,26) = pd(15,26) + rrt(331) * density(38) 
  pd(15,38) = pd(15,38) + rrt(331) * density(26) 
  pd(26,26) = pd(26,26) - rrt(331) * density(38) 
  pd(26,38) = pd(26,38) - rrt(331) * density(26) 
  pd(38,26) = pd(38,26) - rrt(331) * density(38) 
  pd(38,38) = pd(38,38) - rrt(331) * density(26) 
  pd(02,26) = pd(02,26) + rrt(332) * density(35) 
  pd(02,35) = pd(02,35) + rrt(332) * density(26) 
  pd(26,26) = pd(26,26) - rrt(332) * density(35) 
  pd(26,35) = pd(26,35) - rrt(332) * density(26) 
  pd(35,26) = pd(35,26) - rrt(332) * density(35) 
  pd(35,35) = pd(35,35) - rrt(332) * density(26) 
  pd(44,26) = pd(44,26) + rrt(332) * density(35) 
  pd(44,35) = pd(44,35) + rrt(332) * density(26) 
  pd(02,26) = pd(02,26) + rrt(333) * density(35) * 2.0d0
  pd(02,35) = pd(02,35) + rrt(333) * density(26) * 2.0d0
  pd(06,26) = pd(06,26) + rrt(333) * density(35) 
  pd(06,35) = pd(06,35) + rrt(333) * density(26) 
  pd(26,26) = pd(26,26) - rrt(333) * density(35) 
  pd(26,35) = pd(26,35) - rrt(333) * density(26) 
  pd(35,26) = pd(35,26) - rrt(333) * density(35) 
  pd(35,35) = pd(35,35) - rrt(333) * density(26) 
  pd(02,26) = pd(02,26) + rrt(334) * density(50) 
  pd(02,50) = pd(02,50) + rrt(334) * density(26) 
  pd(26,26) = pd(26,26) - rrt(334) * density(50) 
  pd(26,50) = pd(26,50) - rrt(334) * density(26) 
  pd(32,26) = pd(32,26) + rrt(334) * density(50) 
  pd(32,50) = pd(32,50) + rrt(334) * density(26) 
  pd(50,26) = pd(50,26) - rrt(334) * density(50) 
  pd(50,50) = pd(50,50) - rrt(334) * density(26) 
  pd(02,18) = pd(02,18) + rrt(335) * density(26) 
  pd(02,26) = pd(02,26) + rrt(335) * density(18) 
  pd(14,18) = pd(14,18) + rrt(335) * density(26) 
  pd(14,26) = pd(14,26) + rrt(335) * density(18) 
  pd(18,18) = pd(18,18) - rrt(335) * density(26) 
  pd(18,26) = pd(18,26) - rrt(335) * density(18) 
  pd(26,18) = pd(26,18) - rrt(335) * density(26) 
  pd(26,26) = pd(26,26) - rrt(335) * density(18) 
  pd(02,19) = pd(02,19) + rrt(336) * density(26) 
  pd(02,26) = pd(02,26) + rrt(336) * density(19) 
  pd(06,19) = pd(06,19) + rrt(336) * density(26) 
  pd(06,26) = pd(06,26) + rrt(336) * density(19) 
  pd(19,19) = pd(19,19) - rrt(336) * density(26) 
  pd(19,26) = pd(19,26) - rrt(336) * density(19) 
  pd(26,19) = pd(26,19) - rrt(336) * density(26) 
  pd(26,26) = pd(26,26) - rrt(336) * density(19) 
  pd(05,05) = pd(05,05) - rrt(337) * density(45) 
  pd(05,45) = pd(05,45) - rrt(337) * density(05) 
  pd(13,05) = pd(13,05) + rrt(337) * density(45) 
  pd(13,45) = pd(13,45) + rrt(337) * density(05) 
  pd(41,05) = pd(41,05) + rrt(337) * density(45) 
  pd(41,45) = pd(41,45) + rrt(337) * density(05) 
  pd(45,05) = pd(45,05) - rrt(337) * density(45) 
  pd(45,45) = pd(45,45) - rrt(337) * density(05) 
  pd(02,05) = pd(02,05) + rrt(338) * density(45) 
  pd(02,45) = pd(02,45) + rrt(338) * density(05) 
  pd(05,05) = pd(05,05) - rrt(338) * density(45) 
  pd(05,45) = pd(05,45) - rrt(338) * density(05) 
  pd(25,05) = pd(25,05) + rrt(338) * density(45) 
  pd(25,45) = pd(25,45) + rrt(338) * density(05) 
  pd(45,05) = pd(45,05) - rrt(338) * density(45) 
  pd(45,45) = pd(45,45) - rrt(338) * density(05) 
  pd(02,05) = pd(02,05) + rrt(339) * density(45) 
  pd(02,45) = pd(02,45) + rrt(339) * density(05) 
  pd(05,05) = pd(05,05) - rrt(339) * density(45) 
  pd(05,45) = pd(05,45) - rrt(339) * density(05) 
  pd(13,05) = pd(13,05) + rrt(339) * density(45) 
  pd(13,45) = pd(13,45) + rrt(339) * density(05) 
  pd(45,05) = pd(45,05) - rrt(339) * density(45) 
  pd(45,45) = pd(45,45) - rrt(339) * density(05) 
  pd(48,05) = pd(48,05) + rrt(339) * density(45) 
  pd(48,45) = pd(48,45) + rrt(339) * density(05) 
  pd(10,10) = pd(10,10) - rrt(340) * density(45) 
  pd(10,45) = pd(10,45) - rrt(340) * density(10) 
  pd(13,10) = pd(13,10) + rrt(340) * density(45) 
  pd(13,45) = pd(13,45) + rrt(340) * density(10) 
  pd(45,10) = pd(45,10) - rrt(340) * density(45) 
  pd(45,45) = pd(45,45) - rrt(340) * density(10) 
  pd(48,10) = pd(48,10) + rrt(340) * density(45) 
  pd(48,45) = pd(48,45) + rrt(340) * density(10) 
  pd(02,10) = pd(02,10) + rrt(341) * density(45) 
  pd(02,45) = pd(02,45) + rrt(341) * density(10) 
  pd(10,10) = pd(10,10) - rrt(341) * density(45) 
  pd(10,45) = pd(10,45) - rrt(341) * density(10) 
  pd(27,10) = pd(27,10) + rrt(341) * density(45) 
  pd(27,45) = pd(27,45) + rrt(341) * density(10) 
  pd(45,10) = pd(45,10) - rrt(341) * density(45) 
  pd(45,45) = pd(45,45) - rrt(341) * density(10) 
  pd(13,45) = pd(13,45) + rrt(342) * density(46) 
  pd(13,46) = pd(13,46) + rrt(342) * density(45) 
  pd(27,45) = pd(27,45) + rrt(342) * density(46) 
  pd(27,46) = pd(27,46) + rrt(342) * density(45) 
  pd(45,45) = pd(45,45) - rrt(342) * density(46) 
  pd(45,46) = pd(45,46) - rrt(342) * density(45) 
  pd(46,45) = pd(46,45) - rrt(342) * density(46) 
  pd(46,46) = pd(46,46) - rrt(342) * density(45) 
  pd(02,45) = pd(02,45) + rrt(343) * density(46) 
  pd(02,46) = pd(02,46) + rrt(343) * density(45) 
  pd(28,45) = pd(28,45) + rrt(343) * density(46) 
  pd(28,46) = pd(28,46) + rrt(343) * density(45) 
  pd(45,45) = pd(45,45) - rrt(343) * density(46) 
  pd(45,46) = pd(45,46) - rrt(343) * density(45) 
  pd(46,45) = pd(46,45) - rrt(343) * density(46) 
  pd(46,46) = pd(46,46) - rrt(343) * density(45) 
  pd(02,20) = pd(02,20) + rrt(344) * density(45) 
  pd(02,45) = pd(02,45) + rrt(344) * density(20) 
  pd(15,20) = pd(15,20) + rrt(344) * density(45) 
  pd(15,45) = pd(15,45) + rrt(344) * density(20) 
  pd(20,20) = pd(20,20) - rrt(344) * density(45) 
  pd(20,45) = pd(20,45) - rrt(344) * density(20) 
  pd(45,20) = pd(45,20) - rrt(344) * density(45) 
  pd(45,45) = pd(45,45) - rrt(344) * density(20) 
  pd(02,20) = pd(02,20) + rrt(345) * density(45) 
  pd(02,45) = pd(02,45) + rrt(345) * density(20) 
  pd(13,20) = pd(13,20) + rrt(345) * density(45) 
  pd(13,45) = pd(13,45) + rrt(345) * density(20) 
  pd(20,20) = pd(20,20) - rrt(345) * density(45) 
  pd(20,45) = pd(20,45) - rrt(345) * density(20) 
  pd(44,20) = pd(44,20) + rrt(345) * density(45) 
  pd(44,45) = pd(44,45) + rrt(345) * density(20) 
  pd(45,20) = pd(45,20) - rrt(345) * density(45) 
  pd(45,45) = pd(45,45) - rrt(345) * density(20) 
  pd(02,20) = pd(02,20) + rrt(346) * density(45) * 2.0d0
  pd(02,45) = pd(02,45) + rrt(346) * density(20) * 2.0d0
  pd(20,20) = pd(20,20) - rrt(346) * density(45) 
  pd(20,45) = pd(20,45) - rrt(346) * density(20) 
  pd(32,20) = pd(32,20) + rrt(346) * density(45) 
  pd(32,45) = pd(32,45) + rrt(346) * density(20) 
  pd(45,20) = pd(45,20) - rrt(346) * density(45) 
  pd(45,45) = pd(45,45) - rrt(346) * density(20) 
  pd(02,20) = pd(02,20) + rrt(347) * density(45) * 2.0d0
  pd(02,45) = pd(02,45) + rrt(347) * density(20) * 2.0d0
  pd(06,20) = pd(06,20) + rrt(347) * density(45) 
  pd(06,45) = pd(06,45) + rrt(347) * density(20) 
  pd(13,20) = pd(13,20) + rrt(347) * density(45) 
  pd(13,45) = pd(13,45) + rrt(347) * density(20) 
  pd(20,20) = pd(20,20) - rrt(347) * density(45) 
  pd(20,45) = pd(20,45) - rrt(347) * density(20) 
  pd(45,20) = pd(45,20) - rrt(347) * density(45) 
  pd(45,45) = pd(45,45) - rrt(347) * density(20) 
  pd(02,20) = pd(02,20) + rrt(348) * density(45) * 3.0d0
  pd(02,45) = pd(02,45) + rrt(348) * density(20) * 3.0d0
  pd(14,20) = pd(14,20) + rrt(348) * density(45) 
  pd(14,45) = pd(14,45) + rrt(348) * density(20) 
  pd(20,20) = pd(20,20) - rrt(348) * density(45) 
  pd(20,45) = pd(20,45) - rrt(348) * density(20) 
  pd(45,20) = pd(45,20) - rrt(348) * density(45) 
  pd(45,45) = pd(45,45) - rrt(348) * density(20) 
  pd(02,35) = pd(02,35) + rrt(349) * density(45) 
  pd(02,45) = pd(02,45) + rrt(349) * density(35) 
  pd(32,35) = pd(32,35) + rrt(349) * density(45) 
  pd(32,45) = pd(32,45) + rrt(349) * density(35) 
  pd(35,35) = pd(35,35) - rrt(349) * density(45) 
  pd(35,45) = pd(35,45) - rrt(349) * density(35) 
  pd(45,35) = pd(45,35) - rrt(349) * density(45) 
  pd(45,45) = pd(45,45) - rrt(349) * density(35) 
  pd(02,35) = pd(02,35) + rrt(350) * density(45) 
  pd(02,45) = pd(02,45) + rrt(350) * density(35) 
  pd(06,35) = pd(06,35) + rrt(350) * density(45) 
  pd(06,45) = pd(06,45) + rrt(350) * density(35) 
  pd(13,35) = pd(13,35) + rrt(350) * density(45) 
  pd(13,45) = pd(13,45) + rrt(350) * density(35) 
  pd(35,35) = pd(35,35) - rrt(350) * density(45) 
  pd(35,45) = pd(35,45) - rrt(350) * density(35) 
  pd(45,35) = pd(45,35) - rrt(350) * density(45) 
  pd(45,45) = pd(45,45) - rrt(350) * density(35) 
  pd(02,35) = pd(02,35) + rrt(351) * density(45) * 2.0d0
  pd(02,45) = pd(02,45) + rrt(351) * density(35) * 2.0d0
  pd(14,35) = pd(14,35) + rrt(351) * density(45) 
  pd(14,45) = pd(14,45) + rrt(351) * density(35) 
  pd(35,35) = pd(35,35) - rrt(351) * density(45) 
  pd(35,45) = pd(35,45) - rrt(351) * density(35) 
  pd(45,35) = pd(45,35) - rrt(351) * density(45) 
  pd(45,45) = pd(45,45) - rrt(351) * density(35) 
  pd(06,19) = pd(06,19) + rrt(352) * density(45) 
  pd(06,45) = pd(06,45) + rrt(352) * density(19) 
  pd(13,19) = pd(13,19) + rrt(352) * density(45) 
  pd(13,45) = pd(13,45) + rrt(352) * density(19) 
  pd(19,19) = pd(19,19) - rrt(352) * density(45) 
  pd(19,45) = pd(19,45) - rrt(352) * density(19) 
  pd(45,19) = pd(45,19) - rrt(352) * density(45) 
  pd(45,45) = pd(45,45) - rrt(352) * density(19) 
  pd(02,19) = pd(02,19) + rrt(353) * density(45) 
  pd(02,45) = pd(02,45) + rrt(353) * density(19) 
  pd(14,19) = pd(14,19) + rrt(353) * density(45) 
  pd(14,45) = pd(14,45) + rrt(353) * density(19) 
  pd(19,19) = pd(19,19) - rrt(353) * density(45) 
  pd(19,45) = pd(19,45) - rrt(353) * density(19) 
  pd(45,19) = pd(45,19) - rrt(353) * density(45) 
  pd(45,45) = pd(45,45) - rrt(353) * density(19) 
  pd(13,13) = pd(13,13) - rrt(354) * density(45) 
  pd(13,45) = pd(13,45) - rrt(354) * density(13) 
  pd(26,13) = pd(26,13) + rrt(354) * density(45) 
  pd(26,45) = pd(26,45) + rrt(354) * density(13) 
  pd(45,13) = pd(45,13) - rrt(354) * density(45) 
  pd(45,45) = pd(45,45) - rrt(354) * density(13) 
  pd(02,13) = pd(02,13) + rrt(355) * density(45) 
  pd(02,45) = pd(02,45) + rrt(355) * density(13) 
  pd(07,13) = pd(07,13) + rrt(355) * density(45) 
  pd(07,45) = pd(07,45) + rrt(355) * density(13) 
  pd(13,13) = pd(13,13) - rrt(355) * density(45) 
  pd(13,45) = pd(13,45) - rrt(355) * density(13) 
  pd(45,13) = pd(45,13) - rrt(355) * density(45) 
  pd(45,45) = pd(45,45) - rrt(355) * density(13) 
  pd(05,05) = pd(05,05) - rrt(356) * density(07) 
  pd(05,07) = pd(05,07) - rrt(356) * density(05) 
  pd(07,05) = pd(07,05) - rrt(356) * density(07) 
  pd(07,07) = pd(07,07) - rrt(356) * density(05) 
  pd(13,05) = pd(13,05) + rrt(356) * density(07) 
  pd(13,07) = pd(13,07) + rrt(356) * density(05) 
  pd(25,05) = pd(25,05) + rrt(356) * density(07) 
  pd(25,07) = pd(25,07) + rrt(356) * density(05) 
  pd(02,05) = pd(02,05) + rrt(357) * density(07) 
  pd(02,07) = pd(02,07) + rrt(357) * density(05) 
  pd(05,05) = pd(05,05) - rrt(357) * density(07) 
  pd(05,07) = pd(05,07) - rrt(357) * density(05) 
  pd(07,05) = pd(07,05) - rrt(357) * density(07) 
  pd(07,07) = pd(07,07) - rrt(357) * density(05) 
  pd(48,05) = pd(48,05) + rrt(357) * density(07) 
  pd(48,07) = pd(48,07) + rrt(357) * density(05) 
  pd(07,07) = pd(07,07) - rrt(358) * density(24) 
  pd(07,24) = pd(07,24) - rrt(358) * density(07) 
  pd(13,07) = pd(13,07) + rrt(358) * density(24) 
  pd(13,24) = pd(13,24) + rrt(358) * density(07) 
  pd(24,07) = pd(24,07) - rrt(358) * density(24) 
  pd(24,24) = pd(24,24) - rrt(358) * density(07) 
  pd(48,07) = pd(48,07) + rrt(358) * density(24) 
  pd(48,24) = pd(48,24) + rrt(358) * density(07) 
  pd(07,07) = pd(07,07) - rrt(359) * density(10) 
  pd(07,10) = pd(07,10) - rrt(359) * density(07) 
  pd(10,07) = pd(10,07) - rrt(359) * density(10) 
  pd(10,10) = pd(10,10) - rrt(359) * density(07) 
  pd(13,07) = pd(13,07) + rrt(359) * density(10) 
  pd(13,10) = pd(13,10) + rrt(359) * density(07) 
  pd(27,07) = pd(27,07) + rrt(359) * density(10) 
  pd(27,10) = pd(27,10) + rrt(359) * density(07) 
  pd(02,07) = pd(02,07) + rrt(360) * density(10) 
  pd(02,10) = pd(02,10) + rrt(360) * density(07) 
  pd(07,07) = pd(07,07) - rrt(360) * density(10) 
  pd(07,10) = pd(07,10) - rrt(360) * density(07) 
  pd(10,07) = pd(10,07) - rrt(360) * density(10) 
  pd(10,10) = pd(10,10) - rrt(360) * density(07) 
  pd(28,07) = pd(28,07) + rrt(360) * density(10) 
  pd(28,10) = pd(28,10) + rrt(360) * density(07) 
  pd(07,07) = pd(07,07) - rrt(361) * density(46) 
  pd(07,46) = pd(07,46) - rrt(361) * density(07) 
  pd(13,07) = pd(13,07) + rrt(361) * density(46) 
  pd(13,46) = pd(13,46) + rrt(361) * density(07) 
  pd(28,07) = pd(28,07) + rrt(361) * density(46) 
  pd(28,46) = pd(28,46) + rrt(361) * density(07) 
  pd(46,07) = pd(46,07) - rrt(361) * density(46) 
  pd(46,46) = pd(46,46) - rrt(361) * density(07) 
  pd(02,07) = pd(02,07) + rrt(362) * density(20) 
  pd(02,20) = pd(02,20) + rrt(362) * density(07) 
  pd(07,07) = pd(07,07) - rrt(362) * density(20) 
  pd(07,20) = pd(07,20) - rrt(362) * density(07) 
  pd(20,07) = pd(20,07) - rrt(362) * density(20) 
  pd(20,20) = pd(20,20) - rrt(362) * density(07) 
  pd(44,07) = pd(44,07) + rrt(362) * density(20) 
  pd(44,20) = pd(44,20) + rrt(362) * density(07) 
  pd(02,07) = pd(02,07) + rrt(363) * density(20) 
  pd(02,20) = pd(02,20) + rrt(363) * density(07) 
  pd(07,07) = pd(07,07) - rrt(363) * density(20) 
  pd(07,20) = pd(07,20) - rrt(363) * density(07) 
  pd(13,07) = pd(13,07) + rrt(363) * density(20) 
  pd(13,20) = pd(13,20) + rrt(363) * density(07) 
  pd(20,07) = pd(20,07) - rrt(363) * density(20) 
  pd(20,20) = pd(20,20) - rrt(363) * density(07) 
  pd(32,07) = pd(32,07) + rrt(363) * density(20) 
  pd(32,20) = pd(32,20) + rrt(363) * density(07) 
  pd(02,07) = pd(02,07) + rrt(364) * density(20) * 2.0d0
  pd(02,20) = pd(02,20) + rrt(364) * density(07) * 2.0d0
  pd(06,07) = pd(06,07) + rrt(364) * density(20) 
  pd(06,20) = pd(06,20) + rrt(364) * density(07) 
  pd(07,07) = pd(07,07) - rrt(364) * density(20) 
  pd(07,20) = pd(07,20) - rrt(364) * density(07) 
  pd(20,07) = pd(20,07) - rrt(364) * density(20) 
  pd(20,20) = pd(20,20) - rrt(364) * density(07) 
  pd(02,07) = pd(02,07) + rrt(365) * density(38) 
  pd(02,38) = pd(02,38) + rrt(365) * density(07) 
  pd(07,07) = pd(07,07) - rrt(365) * density(38) 
  pd(07,38) = pd(07,38) - rrt(365) * density(07) 
  pd(32,07) = pd(32,07) + rrt(365) * density(38) 
  pd(32,38) = pd(32,38) + rrt(365) * density(07) 
  pd(38,07) = pd(38,07) - rrt(365) * density(38) 
  pd(38,38) = pd(38,38) - rrt(365) * density(07) 
  pd(02,07) = pd(02,07) + rrt(366) * density(38) 
  pd(02,38) = pd(02,38) + rrt(366) * density(07) 
  pd(06,07) = pd(06,07) + rrt(366) * density(38) 
  pd(06,38) = pd(06,38) + rrt(366) * density(07) 
  pd(07,07) = pd(07,07) - rrt(366) * density(38) 
  pd(07,38) = pd(07,38) - rrt(366) * density(07) 
  pd(13,07) = pd(13,07) + rrt(366) * density(38) 
  pd(13,38) = pd(13,38) + rrt(366) * density(07) 
  pd(38,07) = pd(38,07) - rrt(366) * density(38) 
  pd(38,38) = pd(38,38) - rrt(366) * density(07) 
  pd(07,07) = pd(07,07) - rrt(367) * density(35) 
  pd(07,35) = pd(07,35) - rrt(367) * density(07) 
  pd(13,07) = pd(13,07) + rrt(367) * density(35) 
  pd(13,35) = pd(13,35) + rrt(367) * density(07) 
  pd(32,07) = pd(32,07) + rrt(367) * density(35) 
  pd(32,35) = pd(32,35) + rrt(367) * density(07) 
  pd(35,07) = pd(35,07) - rrt(367) * density(35) 
  pd(35,35) = pd(35,35) - rrt(367) * density(07) 
  pd(02,07) = pd(02,07) + rrt(368) * density(35) 
  pd(02,35) = pd(02,35) + rrt(368) * density(07) 
  pd(06,07) = pd(06,07) + rrt(368) * density(35) 
  pd(06,35) = pd(06,35) + rrt(368) * density(07) 
  pd(07,07) = pd(07,07) - rrt(368) * density(35) 
  pd(07,35) = pd(07,35) - rrt(368) * density(07) 
  pd(35,07) = pd(35,07) - rrt(368) * density(35) 
  pd(35,35) = pd(35,35) - rrt(368) * density(07) 
  pd(02,07) = pd(02,07) + rrt(369) * density(35) 
  pd(02,35) = pd(02,35) + rrt(369) * density(07) 
  pd(07,07) = pd(07,07) - rrt(369) * density(35) 
  pd(07,35) = pd(07,35) - rrt(369) * density(07) 
  pd(13,07) = pd(13,07) + rrt(369) * density(35) 
  pd(13,35) = pd(13,35) + rrt(369) * density(07) 
  pd(14,07) = pd(14,07) + rrt(369) * density(35) 
  pd(14,35) = pd(14,35) + rrt(369) * density(07) 
  pd(35,07) = pd(35,07) - rrt(369) * density(35) 
  pd(35,35) = pd(35,35) - rrt(369) * density(07) 
  pd(06,07) = pd(06,07) + rrt(370) * density(50) 
  pd(06,50) = pd(06,50) + rrt(370) * density(07) 
  pd(07,07) = pd(07,07) - rrt(370) * density(50) 
  pd(07,50) = pd(07,50) - rrt(370) * density(07) 
  pd(13,07) = pd(13,07) + rrt(370) * density(50) 
  pd(13,50) = pd(13,50) + rrt(370) * density(07) 
  pd(50,07) = pd(50,07) - rrt(370) * density(50) 
  pd(50,50) = pd(50,50) - rrt(370) * density(07) 
  pd(02,07) = pd(02,07) + rrt(371) * density(50) 
  pd(02,50) = pd(02,50) + rrt(371) * density(07) 
  pd(07,07) = pd(07,07) - rrt(371) * density(50) 
  pd(07,50) = pd(07,50) - rrt(371) * density(07) 
  pd(14,07) = pd(14,07) + rrt(371) * density(50) 
  pd(14,50) = pd(14,50) + rrt(371) * density(07) 
  pd(50,07) = pd(50,07) - rrt(371) * density(50) 
  pd(50,50) = pd(50,50) - rrt(371) * density(07) 
  pd(07,07) = pd(07,07) - rrt(372) * density(19) 
  pd(07,19) = pd(07,19) - rrt(372) * density(07) 
  pd(13,07) = pd(13,07) + rrt(372) * density(19) 
  pd(13,19) = pd(13,19) + rrt(372) * density(07) 
  pd(14,07) = pd(14,07) + rrt(372) * density(19) 
  pd(14,19) = pd(14,19) + rrt(372) * density(07) 
  pd(19,07) = pd(19,07) - rrt(372) * density(19) 
  pd(19,19) = pd(19,19) - rrt(372) * density(07) 
  pd(05,05) = pd(05,05) - rrt(373) * density(10) 
  pd(05,10) = pd(05,10) - rrt(373) * density(05) 
  pd(10,05) = pd(10,05) - rrt(373) * density(10) 
  pd(10,10) = pd(10,10) - rrt(373) * density(05) 
  pd(24,05) = pd(24,05) + rrt(373) * density(10) * 2.0d0
  pd(24,10) = pd(24,10) + rrt(373) * density(05) * 2.0d0
  pd(05,05) = pd(05,05) - rrt(374) * density(46) 
  pd(05,46) = pd(05,46) - rrt(374) * density(05) 
  pd(13,05) = pd(13,05) + rrt(374) * density(46) 
  pd(13,46) = pd(13,46) + rrt(374) * density(05) 
  pd(35,05) = pd(35,05) + rrt(374) * density(46) 
  pd(35,46) = pd(35,46) + rrt(374) * density(05) 
  pd(46,05) = pd(46,05) - rrt(374) * density(46) 
  pd(46,46) = pd(46,46) - rrt(374) * density(05) 
  pd(05,05) = pd(05,05) - rrt(375) * density(38) 
  pd(05,38) = pd(05,38) - rrt(375) * density(05) 
  pd(20,05) = pd(20,05) + rrt(375) * density(38) 
  pd(20,38) = pd(20,38) + rrt(375) * density(05) 
  pd(24,05) = pd(24,05) + rrt(375) * density(38) 
  pd(24,38) = pd(24,38) + rrt(375) * density(05) 
  pd(38,05) = pd(38,05) - rrt(375) * density(38) 
  pd(38,38) = pd(38,38) - rrt(375) * density(05) 
  pd(05,05) = pd(05,05) - rrt(376) * density(50) 
  pd(05,50) = pd(05,50) - rrt(376) * density(05) 
  pd(24,05) = pd(24,05) + rrt(376) * density(50) 
  pd(24,50) = pd(24,50) + rrt(376) * density(05) 
  pd(35,05) = pd(35,05) + rrt(376) * density(50) 
  pd(35,50) = pd(35,50) + rrt(376) * density(05) 
  pd(50,05) = pd(50,05) - rrt(376) * density(50) 
  pd(50,50) = pd(50,50) - rrt(376) * density(05) 
  pd(05,05) = pd(05,05) - rrt(377) * density(18) 
  pd(05,18) = pd(05,18) - rrt(377) * density(05) 
  pd(18,05) = pd(18,05) - rrt(377) * density(18) 
  pd(18,18) = pd(18,18) - rrt(377) * density(05) 
  pd(19,05) = pd(19,05) + rrt(377) * density(18) 
  pd(19,18) = pd(19,18) + rrt(377) * density(05) 
  pd(24,05) = pd(24,05) + rrt(377) * density(18) 
  pd(24,18) = pd(24,18) + rrt(377) * density(05) 
  pd(05,05) = pd(05,05) - rrt(378) * density(23) 
  pd(05,23) = pd(05,23) - rrt(378) * density(05) 
  pd(12,05) = pd(12,05) + rrt(378) * density(23) 
  pd(12,23) = pd(12,23) + rrt(378) * density(05) 
  pd(23,05) = pd(23,05) - rrt(378) * density(23) 
  pd(23,23) = pd(23,23) - rrt(378) * density(05) 
  pd(24,05) = pd(24,05) + rrt(378) * density(23) 
  pd(24,23) = pd(24,23) + rrt(378) * density(05) 
  pd(03,03) = pd(03,03) - rrt(379) * density(05) 
  pd(03,05) = pd(03,05) - rrt(379) * density(03) 
  pd(05,03) = pd(05,03) - rrt(379) * density(05) 
  pd(05,05) = pd(05,05) - rrt(379) * density(03) 
  pd(09,03) = pd(09,03) + rrt(379) * density(05) 
  pd(09,05) = pd(09,05) + rrt(379) * density(03) 
  pd(24,03) = pd(24,03) + rrt(379) * density(05) 
  pd(24,05) = pd(24,05) + rrt(379) * density(03) 
  pd(02,05) = pd(02,05) + rrt(380) * density(13) 
  pd(02,13) = pd(02,13) + rrt(380) * density(05) 
  pd(05,05) = pd(05,05) - rrt(380) * density(13) 
  pd(05,13) = pd(05,13) - rrt(380) * density(05) 
  pd(13,05) = pd(13,05) - rrt(380) * density(13) 
  pd(13,13) = pd(13,13) - rrt(380) * density(05) 
  pd(24,05) = pd(24,05) + rrt(380) * density(13) 
  pd(24,13) = pd(24,13) + rrt(380) * density(05) 
  pd(05,05) = pd(05,05) - rrt(381) * density(24) 
  pd(05,24) = pd(05,24) - rrt(381) * density(05) 
  pd(13,05) = pd(13,05) + rrt(381) * density(24) 
  pd(13,24) = pd(13,24) + rrt(381) * density(05) 
  pd(20,05) = pd(20,05) + rrt(381) * density(24) 
  pd(20,24) = pd(20,24) + rrt(381) * density(05) 
  pd(24,05) = pd(24,05) - rrt(381) * density(24) 
  pd(24,24) = pd(24,24) - rrt(381) * density(05) 
  pd(05,05) = pd(05,05) - rrt(382) * density(10) 
  pd(05,10) = pd(05,10) - rrt(382) * density(05) 
  pd(10,05) = pd(10,05) - rrt(382) * density(10) 
  pd(10,10) = pd(10,10) - rrt(382) * density(05) 
  pd(20,05) = pd(20,05) + rrt(382) * density(10) 
  pd(20,10) = pd(20,10) + rrt(382) * density(05) 
  pd(05,05) = pd(05,05) - rrt(383) 
  pd(13,05) = pd(13,05) + rrt(383) 
  pd(24,05) = pd(24,05) + rrt(383) 
  pd(10,24) = pd(10,24) + rrt(384) 
  pd(13,24) = pd(13,24) + rrt(384) 
  pd(24,24) = pd(24,24) - rrt(384) 
  pd(02,24) = pd(02,24) + rrt(385) 
  pd(24,24) = pd(24,24) - rrt(385) 
  pd(46,24) = pd(46,24) + rrt(385) 
  pd(10,24) = pd(10,24) + rrt(386) * density(38) 
  pd(10,38) = pd(10,38) + rrt(386) * density(24) 
  pd(20,24) = pd(20,24) + rrt(386) * density(38) 
  pd(20,38) = pd(20,38) + rrt(386) * density(24) 
  pd(24,24) = pd(24,24) - rrt(386) * density(38) 
  pd(24,38) = pd(24,38) - rrt(386) * density(24) 
  pd(38,24) = pd(38,24) - rrt(386) * density(38) 
  pd(38,38) = pd(38,38) - rrt(386) * density(24) 
  pd(10,10) = pd(10,10) - rrt(387) * density(10) * 4.0d0
  pd(13,10) = pd(13,10) + rrt(387) * density(10) * 4.0d0
  pd(19,10) = pd(19,10) + rrt(387) * density(10) * 2.0d0
  pd(10,10) = pd(10,10) - rrt(388) * density(38) 
  pd(10,38) = pd(10,38) - rrt(388) * density(10) 
  pd(24,10) = pd(24,10) + rrt(388) * density(38) 
  pd(24,38) = pd(24,38) + rrt(388) * density(10) 
  pd(35,10) = pd(35,10) + rrt(388) * density(38) 
  pd(35,38) = pd(35,38) + rrt(388) * density(10) 
  pd(38,10) = pd(38,10) - rrt(388) * density(38) 
  pd(38,38) = pd(38,38) - rrt(388) * density(10) 
  pd(10,10) = pd(10,10) - rrt(389) * density(50) 
  pd(10,50) = pd(10,50) - rrt(389) * density(10) 
  pd(19,10) = pd(19,10) + rrt(389) * density(50) 
  pd(19,50) = pd(19,50) + rrt(389) * density(10) 
  pd(24,10) = pd(24,10) + rrt(389) * density(50) 
  pd(24,50) = pd(24,50) + rrt(389) * density(10) 
  pd(50,10) = pd(50,10) - rrt(389) * density(50) 
  pd(50,50) = pd(50,50) - rrt(389) * density(10) 
  pd(10,10) = pd(10,10) - rrt(390) * density(18) 
  pd(10,18) = pd(10,18) - rrt(390) * density(10) 
  pd(18,10) = pd(18,10) - rrt(390) * density(18) 
  pd(18,18) = pd(18,18) - rrt(390) * density(10) 
  pd(19,10) = pd(19,10) + rrt(390) * density(18) 
  pd(19,18) = pd(19,18) + rrt(390) * density(10) 
  pd(46,10) = pd(46,10) + rrt(390) * density(18) 
  pd(46,18) = pd(46,18) + rrt(390) * density(10) 
  pd(10,10) = pd(10,10) - rrt(391) * density(12) 
  pd(10,12) = pd(10,12) - rrt(391) * density(10) 
  pd(12,10) = pd(12,10) - rrt(391) * density(12) 
  pd(12,12) = pd(12,12) - rrt(391) * density(10) 
  pd(23,10) = pd(23,10) + rrt(391) * density(12) 
  pd(23,12) = pd(23,12) + rrt(391) * density(10) 
  pd(24,10) = pd(24,10) + rrt(391) * density(12) 
  pd(24,12) = pd(24,12) + rrt(391) * density(10) 
  pd(10,10) = pd(10,10) - rrt(392) * density(23) 
  pd(10,23) = pd(10,23) - rrt(392) * density(10) 
  pd(23,10) = pd(23,10) - rrt(392) * density(23) 
  pd(23,23) = pd(23,23) - rrt(392) * density(10) 
  pd(35,10) = pd(35,10) + rrt(392) * density(23) 
  pd(35,23) = pd(35,23) + rrt(392) * density(10) 
  pd(38,10) = pd(38,10) + rrt(392) * density(23) 
  pd(38,23) = pd(38,23) + rrt(392) * density(10) 
  pd(09,10) = pd(09,10) + rrt(393) * density(23) 
  pd(09,23) = pd(09,23) + rrt(393) * density(10) 
  pd(10,10) = pd(10,10) - rrt(393) * density(23) 
  pd(10,23) = pd(10,23) - rrt(393) * density(10) 
  pd(23,10) = pd(23,10) - rrt(393) * density(23) 
  pd(23,23) = pd(23,23) - rrt(393) * density(10) 
  pd(24,10) = pd(24,10) + rrt(393) * density(23) 
  pd(24,23) = pd(24,23) + rrt(393) * density(10) 
  pd(03,09) = pd(03,09) + rrt(394) * density(10) 
  pd(03,10) = pd(03,10) + rrt(394) * density(09) 
  pd(09,09) = pd(09,09) - rrt(394) * density(10) 
  pd(09,10) = pd(09,10) - rrt(394) * density(09) 
  pd(10,09) = pd(10,09) - rrt(394) * density(10) 
  pd(10,10) = pd(10,10) - rrt(394) * density(09) 
  pd(24,09) = pd(24,09) + rrt(394) * density(10) 
  pd(24,10) = pd(24,10) + rrt(394) * density(09) 
  pd(02,02) = pd(02,02) - rrt(395) * density(10) 
  pd(02,10) = pd(02,10) - rrt(395) * density(02) 
  pd(10,02) = pd(10,02) - rrt(395) * density(10) 
  pd(10,10) = pd(10,10) - rrt(395) * density(02) 
  pd(13,02) = pd(13,02) + rrt(395) * density(10) 
  pd(13,10) = pd(13,10) + rrt(395) * density(02) 
  pd(24,02) = pd(24,02) + rrt(395) * density(10) 
  pd(24,10) = pd(24,10) + rrt(395) * density(02) 
  pd(02,10) = pd(02,10) + rrt(396) * density(13) 
  pd(02,13) = pd(02,13) + rrt(396) * density(10) 
  pd(10,10) = pd(10,10) - rrt(396) * density(13) 
  pd(10,13) = pd(10,13) - rrt(396) * density(10) 
  pd(13,10) = pd(13,10) - rrt(396) * density(13) 
  pd(13,13) = pd(13,13) - rrt(396) * density(10) 
  pd(46,10) = pd(46,10) + rrt(396) * density(13) 
  pd(46,13) = pd(46,13) + rrt(396) * density(10) 
  pd(10,10) = pd(10,10) - rrt(397) 
  pd(13,10) = pd(13,10) + rrt(397) 
  pd(46,10) = pd(46,10) + rrt(397) 
  pd(02,10) = pd(02,10) + rrt(398) * density(10) * 2.0d0
  pd(10,10) = pd(10,10) - rrt(398) * density(10) * 4.0d0
  pd(19,10) = pd(19,10) + rrt(398) * density(10) * 2.0d0
  pd(10,10) = pd(10,10) - rrt(399) * density(13) 
  pd(10,13) = pd(10,13) - rrt(399) * density(10) 
  pd(13,10) = pd(13,10) - rrt(399) * density(13) 
  pd(13,13) = pd(13,13) - rrt(399) * density(10) 
  pd(24,10) = pd(24,10) + rrt(399) * density(13) 
  pd(24,13) = pd(24,13) + rrt(399) * density(10) 
  pd(09,20) = pd(09,20) + rrt(400) * density(46) 
  pd(09,46) = pd(09,46) + rrt(400) * density(20) 
  pd(13,20) = pd(13,20) + rrt(400) * density(46) 
  pd(13,46) = pd(13,46) + rrt(400) * density(20) 
  pd(20,20) = pd(20,20) - rrt(400) * density(46) 
  pd(20,46) = pd(20,46) - rrt(400) * density(20) 
  pd(46,20) = pd(46,20) - rrt(400) * density(46) 
  pd(46,46) = pd(46,46) - rrt(400) * density(20) 
  pd(20,20) = pd(20,20) - rrt(401) * density(46) 
  pd(20,46) = pd(20,46) - rrt(401) * density(20) 
  pd(23,20) = pd(23,20) + rrt(401) * density(46) 
  pd(23,46) = pd(23,46) + rrt(401) * density(20) 
  pd(46,20) = pd(46,20) - rrt(401) * density(46) 
  pd(46,46) = pd(46,46) - rrt(401) * density(20) 
  pd(02,02) = pd(02,02) - rrt(402) * density(46) 
  pd(02,46) = pd(02,46) - rrt(402) * density(02) 
  pd(10,02) = pd(10,02) + rrt(402) * density(46) 
  pd(10,46) = pd(10,46) + rrt(402) * density(02) 
  pd(13,02) = pd(13,02) + rrt(402) * density(46) 
  pd(13,46) = pd(13,46) + rrt(402) * density(02) 
  pd(46,02) = pd(46,02) - rrt(402) * density(46) 
  pd(46,46) = pd(46,46) - rrt(402) * density(02) 
  pd(13,24) = pd(13,24) + rrt(403) * density(46) 
  pd(13,46) = pd(13,46) + rrt(403) * density(24) 
  pd(24,24) = pd(24,24) - rrt(403) * density(46) 
  pd(24,46) = pd(24,46) - rrt(403) * density(24) 
  pd(46,24) = pd(46,24) - rrt(403) * density(46) 
  pd(46,46) = pd(46,46) - rrt(403) * density(24) 
  pd(50,24) = pd(50,24) + rrt(403) * density(46) 
  pd(50,46) = pd(50,46) + rrt(403) * density(24) 
  pd(10,10) = pd(10,10) - rrt(404) * density(46) 
  pd(10,46) = pd(10,46) - rrt(404) * density(10) 
  pd(13,10) = pd(13,10) + rrt(404) * density(46) 
  pd(13,46) = pd(13,46) + rrt(404) * density(10) 
  pd(19,10) = pd(19,10) + rrt(404) * density(46) 
  pd(19,46) = pd(19,46) + rrt(404) * density(10) 
  pd(46,10) = pd(46,10) - rrt(404) * density(46) 
  pd(46,46) = pd(46,46) - rrt(404) * density(10) 
  pd(02,02) = pd(02,02) - rrt(405) * density(46) 
  pd(02,46) = pd(02,46) - rrt(405) * density(02) 
  pd(24,02) = pd(24,02) + rrt(405) * density(46) 
  pd(24,46) = pd(24,46) + rrt(405) * density(02) 
  pd(46,02) = pd(46,02) - rrt(405) * density(46) 
  pd(46,46) = pd(46,46) - rrt(405) * density(02) 
  pd(10,46) = pd(10,46) + rrt(406) * density(50) 
  pd(10,50) = pd(10,50) + rrt(406) * density(46) 
  pd(19,46) = pd(19,46) + rrt(406) * density(50) 
  pd(19,50) = pd(19,50) + rrt(406) * density(46) 
  pd(46,46) = pd(46,46) - rrt(406) * density(50) 
  pd(46,50) = pd(46,50) - rrt(406) * density(46) 
  pd(50,46) = pd(50,46) - rrt(406) * density(50) 
  pd(50,50) = pd(50,50) - rrt(406) * density(46) 
  pd(20,20) = pd(20,20) - rrt(407) * density(50) 
  pd(20,50) = pd(20,50) - rrt(407) * density(20) 
  pd(35,20) = pd(35,20) + rrt(407) * density(50) 
  pd(35,50) = pd(35,50) + rrt(407) * density(20) 
  pd(38,20) = pd(38,20) + rrt(407) * density(50) 
  pd(38,50) = pd(38,50) + rrt(407) * density(20) 
  pd(50,20) = pd(50,20) - rrt(407) * density(50) 
  pd(50,50) = pd(50,50) - rrt(407) * density(20) 
  pd(12,20) = pd(12,20) + rrt(408) * density(23) 
  pd(12,23) = pd(12,23) + rrt(408) * density(20) 
  pd(20,20) = pd(20,20) - rrt(408) * density(23) 
  pd(20,23) = pd(20,23) - rrt(408) * density(20) 
  pd(23,20) = pd(23,20) - rrt(408) * density(23) 
  pd(23,23) = pd(23,23) - rrt(408) * density(20) 
  pd(38,20) = pd(38,20) + rrt(408) * density(23) 
  pd(38,23) = pd(38,23) + rrt(408) * density(20) 
  pd(03,03) = pd(03,03) - rrt(409) * density(20) 
  pd(03,20) = pd(03,20) - rrt(409) * density(03) 
  pd(09,03) = pd(09,03) + rrt(409) * density(20) 
  pd(09,20) = pd(09,20) + rrt(409) * density(03) 
  pd(20,03) = pd(20,03) - rrt(409) * density(20) 
  pd(20,20) = pd(20,20) - rrt(409) * density(03) 
  pd(38,03) = pd(38,03) + rrt(409) * density(20) 
  pd(38,20) = pd(38,20) + rrt(409) * density(03) 
  pd(02,13) = pd(02,13) + rrt(410) * density(20) 
  pd(02,20) = pd(02,20) + rrt(410) * density(13) 
  pd(13,13) = pd(13,13) - rrt(410) * density(20) 
  pd(13,20) = pd(13,20) - rrt(410) * density(13) 
  pd(20,13) = pd(20,13) - rrt(410) * density(20) 
  pd(20,20) = pd(20,20) - rrt(410) * density(13) 
  pd(38,13) = pd(38,13) + rrt(410) * density(20) 
  pd(38,20) = pd(38,20) + rrt(410) * density(13) 
  pd(05,13) = pd(05,13) + rrt(411) * density(20) 
  pd(05,20) = pd(05,20) + rrt(411) * density(13) 
  pd(13,13) = pd(13,13) - rrt(411) * density(20) 
  pd(13,20) = pd(13,20) - rrt(411) * density(13) 
  pd(20,13) = pd(20,13) - rrt(411) * density(20) 
  pd(20,20) = pd(20,20) - rrt(411) * density(13) 
  pd(24,13) = pd(24,13) + rrt(411) * density(20) 
  pd(24,20) = pd(24,20) + rrt(411) * density(13) 
  pd(20,20) = pd(20,20) - rrt(412) 
  pd(24,20) = pd(24,20) + rrt(412) * 2.0d0
  pd(20,20) = pd(20,20) - rrt(413) * density(46) 
  pd(20,46) = pd(20,46) - rrt(413) * density(20) 
  pd(24,20) = pd(24,20) + rrt(413) * density(46) 
  pd(24,46) = pd(24,46) + rrt(413) * density(20) 
  pd(35,20) = pd(35,20) + rrt(413) * density(46) 
  pd(35,46) = pd(35,46) + rrt(413) * density(20) 
  pd(46,20) = pd(46,20) - rrt(413) * density(46) 
  pd(46,46) = pd(46,46) - rrt(413) * density(20) 
  pd(10,10) = pd(10,10) - rrt(414) * density(20) 
  pd(10,20) = pd(10,20) - rrt(414) * density(10) 
  pd(20,10) = pd(20,10) - rrt(414) * density(20) 
  pd(20,20) = pd(20,20) - rrt(414) * density(10) 
  pd(24,10) = pd(24,10) + rrt(414) * density(20) 
  pd(24,20) = pd(24,20) + rrt(414) * density(10) 
  pd(38,10) = pd(38,10) + rrt(414) * density(20) 
  pd(38,20) = pd(38,20) + rrt(414) * density(10) 
  pd(20,38) = pd(20,38) + rrt(415) * density(38) * 2.0d0
  pd(35,38) = pd(35,38) + rrt(415) * density(38) * 2.0d0
  pd(38,38) = pd(38,38) - rrt(415) * density(38) * 4.0d0
  pd(20,35) = pd(20,35) + rrt(416) * density(38) 
  pd(20,38) = pd(20,38) + rrt(416) * density(35) 
  pd(35,35) = pd(35,35) - rrt(416) * density(38) 
  pd(35,38) = pd(35,38) - rrt(416) * density(35) 
  pd(38,35) = pd(38,35) - rrt(416) * density(38) 
  pd(38,38) = pd(38,38) - rrt(416) * density(35) 
  pd(50,35) = pd(50,35) + rrt(416) * density(38) 
  pd(50,38) = pd(50,38) + rrt(416) * density(35) 
  pd(18,19) = pd(18,19) + rrt(417) * density(38) 
  pd(18,38) = pd(18,38) + rrt(417) * density(19) 
  pd(19,19) = pd(19,19) - rrt(417) * density(38) 
  pd(19,38) = pd(19,38) - rrt(417) * density(19) 
  pd(20,19) = pd(20,19) + rrt(417) * density(38) 
  pd(20,38) = pd(20,38) + rrt(417) * density(19) 
  pd(38,19) = pd(38,19) - rrt(417) * density(38) 
  pd(38,38) = pd(38,38) - rrt(417) * density(19) 
  pd(18,18) = pd(18,18) - rrt(418) * density(38) 
  pd(18,38) = pd(18,38) - rrt(418) * density(18) 
  pd(19,18) = pd(19,18) + rrt(418) * density(38) 
  pd(19,38) = pd(19,38) + rrt(418) * density(18) 
  pd(35,18) = pd(35,18) + rrt(418) * density(38) 
  pd(35,38) = pd(35,38) + rrt(418) * density(18) 
  pd(38,18) = pd(38,18) - rrt(418) * density(38) 
  pd(38,38) = pd(38,38) - rrt(418) * density(18) 
  pd(12,12) = pd(12,12) - rrt(419) * density(38) 
  pd(12,38) = pd(12,38) - rrt(419) * density(12) 
  pd(20,12) = pd(20,12) + rrt(419) * density(38) 
  pd(20,38) = pd(20,38) + rrt(419) * density(12) 
  pd(23,12) = pd(23,12) + rrt(419) * density(38) 
  pd(23,38) = pd(23,38) + rrt(419) * density(12) 
  pd(38,12) = pd(38,12) - rrt(419) * density(38) 
  pd(38,38) = pd(38,38) - rrt(419) * density(12) 
  pd(12,23) = pd(12,23) + rrt(420) * density(38) 
  pd(12,38) = pd(12,38) + rrt(420) * density(23) 
  pd(23,23) = pd(23,23) - rrt(420) * density(38) 
  pd(23,38) = pd(23,38) - rrt(420) * density(23) 
  pd(35,23) = pd(35,23) + rrt(420) * density(38) 
  pd(35,38) = pd(35,38) + rrt(420) * density(23) 
  pd(38,23) = pd(38,23) - rrt(420) * density(38) 
  pd(38,38) = pd(38,38) - rrt(420) * density(23) 
  pd(09,23) = pd(09,23) + rrt(421) * density(38) 
  pd(09,38) = pd(09,38) + rrt(421) * density(23) 
  pd(20,23) = pd(20,23) + rrt(421) * density(38) 
  pd(20,38) = pd(20,38) + rrt(421) * density(23) 
  pd(23,23) = pd(23,23) - rrt(421) * density(38) 
  pd(23,38) = pd(23,38) - rrt(421) * density(23) 
  pd(38,23) = pd(38,23) - rrt(421) * density(38) 
  pd(38,38) = pd(38,38) - rrt(421) * density(23) 
  pd(03,09) = pd(03,09) + rrt(422) * density(38) 
  pd(03,38) = pd(03,38) + rrt(422) * density(09) 
  pd(09,09) = pd(09,09) - rrt(422) * density(38) 
  pd(09,38) = pd(09,38) - rrt(422) * density(09) 
  pd(20,09) = pd(20,09) + rrt(422) * density(38) 
  pd(20,38) = pd(20,38) + rrt(422) * density(09) 
  pd(38,09) = pd(38,09) - rrt(422) * density(38) 
  pd(38,38) = pd(38,38) - rrt(422) * density(09) 
  pd(02,02) = pd(02,02) - rrt(423) * density(38) 
  pd(02,38) = pd(02,38) - rrt(423) * density(02) 
  pd(13,02) = pd(13,02) + rrt(423) * density(38) 
  pd(13,38) = pd(13,38) + rrt(423) * density(02) 
  pd(20,02) = pd(20,02) + rrt(423) * density(38) 
  pd(20,38) = pd(20,38) + rrt(423) * density(02) 
  pd(38,02) = pd(38,02) - rrt(423) * density(38) 
  pd(38,38) = pd(38,38) - rrt(423) * density(02) 
  pd(13,13) = pd(13,13) - rrt(424) * density(38) 
  pd(13,38) = pd(13,38) - rrt(424) * density(13) 
  pd(24,13) = pd(24,13) + rrt(424) * density(38) * 2.0d0
  pd(24,38) = pd(24,38) + rrt(424) * density(13) * 2.0d0
  pd(38,13) = pd(38,13) - rrt(424) * density(38) 
  pd(38,38) = pd(38,38) - rrt(424) * density(13) 
  pd(02,13) = pd(02,13) + rrt(425) * density(38) 
  pd(02,38) = pd(02,38) + rrt(425) * density(13) 
  pd(13,13) = pd(13,13) - rrt(425) * density(38) 
  pd(13,38) = pd(13,38) - rrt(425) * density(13) 
  pd(35,13) = pd(35,13) + rrt(425) * density(38) 
  pd(35,38) = pd(35,38) + rrt(425) * density(13) 
  pd(38,13) = pd(38,13) - rrt(425) * density(38) 
  pd(38,38) = pd(38,38) - rrt(425) * density(13) 
  pd(13,13) = pd(13,13) - rrt(426) * density(38) 
  pd(13,38) = pd(13,38) - rrt(426) * density(13) 
  pd(20,13) = pd(20,13) + rrt(426) * density(38) 
  pd(20,38) = pd(20,38) + rrt(426) * density(13) 
  pd(38,13) = pd(38,13) - rrt(426) * density(38) 
  pd(38,38) = pd(38,38) - rrt(426) * density(13) 
  pd(13,38) = pd(13,38) + rrt(427) 
  pd(35,38) = pd(35,38) + rrt(427) 
  pd(38,38) = pd(38,38) - rrt(427) 
  pd(35,38) = pd(35,38) + rrt(428) * density(50) * 2.0d0
  pd(35,50) = pd(35,50) + rrt(428) * density(38) * 2.0d0
  pd(38,38) = pd(38,38) - rrt(428) * density(50) 
  pd(38,50) = pd(38,50) - rrt(428) * density(38) 
  pd(50,38) = pd(50,38) - rrt(428) * density(50) 
  pd(50,50) = pd(50,50) - rrt(428) * density(38) 
  pd(02,13) = pd(02,13) + rrt(429) * density(35) 
  pd(02,35) = pd(02,35) + rrt(429) * density(13) 
  pd(13,13) = pd(13,13) - rrt(429) * density(35) 
  pd(13,35) = pd(13,35) - rrt(429) * density(13) 
  pd(35,13) = pd(35,13) - rrt(429) * density(35) 
  pd(35,35) = pd(35,35) - rrt(429) * density(13) 
  pd(50,13) = pd(50,13) + rrt(429) * density(35) 
  pd(50,35) = pd(50,35) + rrt(429) * density(13) 
  pd(13,13) = pd(13,13) - rrt(430) * density(35) 
  pd(13,35) = pd(13,35) - rrt(430) * density(13) 
  pd(35,13) = pd(35,13) - rrt(430) * density(35) 
  pd(35,35) = pd(35,35) - rrt(430) * density(13) 
  pd(38,13) = pd(38,13) + rrt(430) * density(35) 
  pd(38,35) = pd(38,35) + rrt(430) * density(13) 
  pd(02,02) = pd(02,02) - rrt(431) * density(35) 
  pd(02,35) = pd(02,35) - rrt(431) * density(02) 
  pd(13,02) = pd(13,02) + rrt(431) * density(35) 
  pd(13,35) = pd(13,35) + rrt(431) * density(02) 
  pd(35,02) = pd(35,02) - rrt(431) * density(35) 
  pd(35,35) = pd(35,35) - rrt(431) * density(02) 
  pd(38,02) = pd(38,02) + rrt(431) * density(35) 
  pd(38,35) = pd(38,35) + rrt(431) * density(02) 
  pd(13,35) = pd(13,35) + rrt(432) 
  pd(35,35) = pd(35,35) - rrt(432) 
  pd(50,35) = pd(50,35) + rrt(432) 
  pd(03,09) = pd(03,09) + rrt(433) * density(35) 
  pd(03,35) = pd(03,35) + rrt(433) * density(09) 
  pd(09,09) = pd(09,09) - rrt(433) * density(35) 
  pd(09,35) = pd(09,35) - rrt(433) * density(09) 
  pd(35,09) = pd(35,09) - rrt(433) * density(35) 
  pd(35,35) = pd(35,35) - rrt(433) * density(09) 
  pd(38,09) = pd(38,09) + rrt(433) * density(35) 
  pd(38,35) = pd(38,35) + rrt(433) * density(09) 
  pd(19,19) = pd(19,19) - rrt(434) * density(35) 
  pd(19,35) = pd(19,35) - rrt(434) * density(19) 
  pd(35,19) = pd(35,19) - rrt(434) * density(35) 
  pd(35,35) = pd(35,35) - rrt(434) * density(19) 
  pd(50,19) = pd(50,19) + rrt(434) * density(35) * 2.0d0
  pd(50,35) = pd(50,35) + rrt(434) * density(19) * 2.0d0
  pd(09,09) = pd(09,09) - rrt(435) * density(35) 
  pd(09,35) = pd(09,35) - rrt(435) * density(09) 
  pd(23,09) = pd(23,09) + rrt(435) * density(35) 
  pd(23,35) = pd(23,35) + rrt(435) * density(09) 
  pd(35,09) = pd(35,09) - rrt(435) * density(35) 
  pd(35,35) = pd(35,35) - rrt(435) * density(09) 
  pd(50,09) = pd(50,09) + rrt(435) * density(35) 
  pd(50,35) = pd(50,35) + rrt(435) * density(09) 
  pd(35,35) = pd(35,35) - rrt(436) * density(35) * 4.0d0
  pd(38,35) = pd(38,35) + rrt(436) * density(35) * 2.0d0
  pd(50,35) = pd(50,35) + rrt(436) * density(35) * 2.0d0
  pd(23,24) = pd(23,24) + rrt(437) * density(35) 
  pd(23,35) = pd(23,35) + rrt(437) * density(24) 
  pd(24,24) = pd(24,24) - rrt(437) * density(35) 
  pd(24,35) = pd(24,35) - rrt(437) * density(24) 
  pd(35,24) = pd(35,24) - rrt(437) * density(35) 
  pd(35,35) = pd(35,35) - rrt(437) * density(24) 
  pd(02,35) = pd(02,35) + rrt(438) 
  pd(19,35) = pd(19,35) + rrt(438) 
  pd(35,35) = pd(35,35) - rrt(438) 
  pd(35,35) = pd(35,35) - rrt(439) * density(38) 
  pd(35,38) = pd(35,38) - rrt(439) * density(35) 
  pd(38,35) = pd(38,35) - rrt(439) * density(38) 
  pd(38,38) = pd(38,38) - rrt(439) * density(35) 
  pd(39,35) = pd(39,35) + rrt(439) * density(38) 
  pd(39,38) = pd(39,38) + rrt(439) * density(35) 
  pd(02,02) = pd(02,02) - rrt(440) * density(35) 
  pd(02,35) = pd(02,35) - rrt(440) * density(02) 
  pd(20,02) = pd(20,02) + rrt(440) * density(35) 
  pd(20,35) = pd(20,35) + rrt(440) * density(02) 
  pd(35,02) = pd(35,02) - rrt(440) * density(35) 
  pd(35,35) = pd(35,35) - rrt(440) * density(02) 
  pd(09,10) = pd(09,10) + rrt(441) * density(35) 
  pd(09,35) = pd(09,35) + rrt(441) * density(10) 
  pd(10,10) = pd(10,10) - rrt(441) * density(35) 
  pd(10,35) = pd(10,35) - rrt(441) * density(10) 
  pd(35,10) = pd(35,10) - rrt(441) * density(35) 
  pd(35,35) = pd(35,35) - rrt(441) * density(10) 
  pd(09,35) = pd(09,35) + rrt(442) * density(39) 
  pd(09,39) = pd(09,39) + rrt(442) * density(35) 
  pd(23,35) = pd(23,35) + rrt(442) * density(39) 
  pd(23,39) = pd(23,39) + rrt(442) * density(35) 
  pd(35,35) = pd(35,35) - rrt(442) * density(39) 
  pd(35,39) = pd(35,39) - rrt(442) * density(35) 
  pd(39,35) = pd(39,35) - rrt(442) * density(39) 
  pd(39,39) = pd(39,39) - rrt(442) * density(35) 
  pd(19,50) = pd(19,50) + rrt(443) * density(50) * 2.0d0
  pd(35,50) = pd(35,50) + rrt(443) * density(50) * 2.0d0
  pd(50,50) = pd(50,50) - rrt(443) * density(50) * 4.0d0
  pd(12,12) = pd(12,12) - rrt(444) * density(50) 
  pd(12,50) = pd(12,50) - rrt(444) * density(12) 
  pd(23,12) = pd(23,12) + rrt(444) * density(50) 
  pd(23,50) = pd(23,50) + rrt(444) * density(12) 
  pd(35,12) = pd(35,12) + rrt(444) * density(50) 
  pd(35,50) = pd(35,50) + rrt(444) * density(12) 
  pd(50,12) = pd(50,12) - rrt(444) * density(50) 
  pd(50,50) = pd(50,50) - rrt(444) * density(12) 
  pd(12,23) = pd(12,23) + rrt(445) * density(50) 
  pd(12,50) = pd(12,50) + rrt(445) * density(23) 
  pd(19,23) = pd(19,23) + rrt(445) * density(50) 
  pd(19,50) = pd(19,50) + rrt(445) * density(23) 
  pd(23,23) = pd(23,23) - rrt(445) * density(50) 
  pd(23,50) = pd(23,50) - rrt(445) * density(23) 
  pd(50,23) = pd(50,23) - rrt(445) * density(50) 
  pd(50,50) = pd(50,50) - rrt(445) * density(23) 
  pd(09,23) = pd(09,23) + rrt(446) * density(50) 
  pd(09,50) = pd(09,50) + rrt(446) * density(23) 
  pd(23,23) = pd(23,23) - rrt(446) * density(50) 
  pd(23,50) = pd(23,50) - rrt(446) * density(23) 
  pd(35,23) = pd(35,23) + rrt(446) * density(50) 
  pd(35,50) = pd(35,50) + rrt(446) * density(23) 
  pd(50,23) = pd(50,23) - rrt(446) * density(50) 
  pd(50,50) = pd(50,50) - rrt(446) * density(23) 
  pd(03,09) = pd(03,09) + rrt(447) * density(50) 
  pd(03,50) = pd(03,50) + rrt(447) * density(09) 
  pd(09,09) = pd(09,09) - rrt(447) * density(50) 
  pd(09,50) = pd(09,50) - rrt(447) * density(09) 
  pd(35,09) = pd(35,09) + rrt(447) * density(50) 
  pd(35,50) = pd(35,50) + rrt(447) * density(09) 
  pd(50,09) = pd(50,09) - rrt(447) * density(50) 
  pd(50,50) = pd(50,50) - rrt(447) * density(09) 
  pd(03,03) = pd(03,03) - rrt(448) * density(50) 
  pd(03,50) = pd(03,50) - rrt(448) * density(03) 
  pd(09,03) = pd(09,03) + rrt(448) * density(50) 
  pd(09,50) = pd(09,50) + rrt(448) * density(03) 
  pd(19,03) = pd(19,03) + rrt(448) * density(50) 
  pd(19,50) = pd(19,50) + rrt(448) * density(03) 
  pd(50,03) = pd(50,03) - rrt(448) * density(50) 
  pd(50,50) = pd(50,50) - rrt(448) * density(03) 
  pd(02,02) = pd(02,02) - rrt(449) * density(50) 
  pd(02,50) = pd(02,50) - rrt(449) * density(02) 
  pd(13,02) = pd(13,02) + rrt(449) * density(50) 
  pd(13,50) = pd(13,50) + rrt(449) * density(02) 
  pd(35,02) = pd(35,02) + rrt(449) * density(50) 
  pd(35,50) = pd(35,50) + rrt(449) * density(02) 
  pd(50,02) = pd(50,02) - rrt(449) * density(50) 
  pd(50,50) = pd(50,50) - rrt(449) * density(02) 
  pd(02,13) = pd(02,13) + rrt(450) * density(50) 
  pd(02,50) = pd(02,50) + rrt(450) * density(13) 
  pd(13,13) = pd(13,13) - rrt(450) * density(50) 
  pd(13,50) = pd(13,50) - rrt(450) * density(13) 
  pd(19,13) = pd(19,13) + rrt(450) * density(50) 
  pd(19,50) = pd(19,50) + rrt(450) * density(13) 
  pd(50,13) = pd(50,13) - rrt(450) * density(50) 
  pd(50,50) = pd(50,50) - rrt(450) * density(13) 
  pd(13,13) = pd(13,13) - rrt(451) * density(50) 
  pd(13,50) = pd(13,50) - rrt(451) * density(13) 
  pd(35,13) = pd(35,13) + rrt(451) * density(50) 
  pd(35,50) = pd(35,50) + rrt(451) * density(13) 
  pd(50,13) = pd(50,13) - rrt(451) * density(50) 
  pd(50,50) = pd(50,50) - rrt(451) * density(13) 
  pd(13,50) = pd(13,50) + rrt(452) 
  pd(19,50) = pd(19,50) + rrt(452) 
  pd(50,50) = pd(50,50) - rrt(452) 
  pd(13,13) = pd(13,13) - rrt(453) * density(19) 
  pd(13,19) = pd(13,19) - rrt(453) * density(13) 
  pd(19,13) = pd(19,13) - rrt(453) * density(19) 
  pd(19,19) = pd(19,19) - rrt(453) * density(13) 
  pd(50,13) = pd(50,13) + rrt(453) * density(19) 
  pd(50,19) = pd(50,19) + rrt(453) * density(13) 
  pd(02,02) = pd(02,02) - rrt(454) * density(19) 
  pd(02,19) = pd(02,19) - rrt(454) * density(02) 
  pd(19,02) = pd(19,02) - rrt(454) * density(19) 
  pd(19,19) = pd(19,19) - rrt(454) * density(02) 
  pd(35,02) = pd(35,02) + rrt(454) * density(19) 
  pd(35,19) = pd(35,19) + rrt(454) * density(02) 
  pd(02,02) = pd(02,02) - rrt(455) * density(19) 
  pd(02,19) = pd(02,19) - rrt(455) * density(02) 
  pd(13,02) = pd(13,02) + rrt(455) * density(19) 
  pd(13,19) = pd(13,19) + rrt(455) * density(02) 
  pd(19,02) = pd(19,02) - rrt(455) * density(19) 
  pd(19,19) = pd(19,19) - rrt(455) * density(02) 
  pd(50,02) = pd(50,02) + rrt(455) * density(19) 
  pd(50,19) = pd(50,19) + rrt(455) * density(02) 
  pd(03,19) = pd(03,19) + rrt(456) * density(24) 
  pd(03,24) = pd(03,24) + rrt(456) * density(19) 
  pd(19,19) = pd(19,19) - rrt(456) * density(24) 
  pd(19,24) = pd(19,24) - rrt(456) * density(19) 
  pd(24,19) = pd(24,19) - rrt(456) * density(24) 
  pd(24,24) = pd(24,24) - rrt(456) * density(19) 
  pd(03,19) = pd(03,19) + rrt(457) * density(39) 
  pd(03,39) = pd(03,39) + rrt(457) * density(19) 
  pd(09,19) = pd(09,19) + rrt(457) * density(39) 
  pd(09,39) = pd(09,39) + rrt(457) * density(19) 
  pd(19,19) = pd(19,19) - rrt(457) * density(39) 
  pd(19,39) = pd(19,39) - rrt(457) * density(19) 
  pd(39,19) = pd(39,19) - rrt(457) * density(39) 
  pd(39,39) = pd(39,39) - rrt(457) * density(19) 
  pd(03,03) = pd(03,03) - rrt(458) * density(12) 
  pd(03,12) = pd(03,12) - rrt(458) * density(03) 
  pd(09,03) = pd(09,03) + rrt(458) * density(12) 
  pd(09,12) = pd(09,12) + rrt(458) * density(03) 
  pd(12,03) = pd(12,03) - rrt(458) * density(12) 
  pd(12,12) = pd(12,12) - rrt(458) * density(03) 
  pd(23,03) = pd(23,03) + rrt(458) * density(12) 
  pd(23,12) = pd(23,12) + rrt(458) * density(03) 
  pd(02,12) = pd(02,12) + rrt(459) * density(13) 
  pd(02,13) = pd(02,13) + rrt(459) * density(12) 
  pd(12,12) = pd(12,12) - rrt(459) * density(13) 
  pd(12,13) = pd(12,13) - rrt(459) * density(12) 
  pd(13,12) = pd(13,12) - rrt(459) * density(13) 
  pd(13,13) = pd(13,13) - rrt(459) * density(12) 
  pd(23,12) = pd(23,12) + rrt(459) * density(13) 
  pd(23,13) = pd(23,13) + rrt(459) * density(12) 
  pd(12,12) = pd(12,12) - rrt(460) 
  pd(24,12) = pd(24,12) + rrt(460) 
  pd(38,12) = pd(38,12) + rrt(460) 
  pd(09,23) = pd(09,23) + rrt(461) * density(23) * 2.0d0
  pd(12,23) = pd(12,23) + rrt(461) * density(23) * 2.0d0
  pd(23,23) = pd(23,23) - rrt(461) * density(23) * 4.0d0
  pd(03,09) = pd(03,09) + rrt(462) * density(23) 
  pd(03,23) = pd(03,23) + rrt(462) * density(09) 
  pd(09,09) = pd(09,09) - rrt(462) * density(23) 
  pd(09,23) = pd(09,23) - rrt(462) * density(09) 
  pd(12,09) = pd(12,09) + rrt(462) * density(23) 
  pd(12,23) = pd(12,23) + rrt(462) * density(09) 
  pd(23,09) = pd(23,09) - rrt(462) * density(23) 
  pd(23,23) = pd(23,23) - rrt(462) * density(09) 
  pd(03,03) = pd(03,03) - rrt(463) * density(23) 
  pd(03,23) = pd(03,23) - rrt(463) * density(03) 
  pd(09,03) = pd(09,03) + rrt(463) * density(23) * 2.0d0
  pd(09,23) = pd(09,23) + rrt(463) * density(03) * 2.0d0
  pd(23,03) = pd(23,03) - rrt(463) * density(23) 
  pd(23,23) = pd(23,23) - rrt(463) * density(03) 
  pd(02,02) = pd(02,02) - rrt(464) * density(23) 
  pd(02,23) = pd(02,23) - rrt(464) * density(02) 
  pd(12,02) = pd(12,02) + rrt(464) * density(23) 
  pd(12,23) = pd(12,23) + rrt(464) * density(02) 
  pd(13,02) = pd(13,02) + rrt(464) * density(23) 
  pd(13,23) = pd(13,23) + rrt(464) * density(02) 
  pd(23,02) = pd(23,02) - rrt(464) * density(23) 
  pd(23,23) = pd(23,23) - rrt(464) * density(02) 
  pd(02,13) = pd(02,13) + rrt(465) * density(23) 
  pd(02,23) = pd(02,23) + rrt(465) * density(13) 
  pd(09,13) = pd(09,13) + rrt(465) * density(23) 
  pd(09,23) = pd(09,23) + rrt(465) * density(13) 
  pd(13,13) = pd(13,13) - rrt(465) * density(23) 
  pd(13,23) = pd(13,23) - rrt(465) * density(13) 
  pd(23,13) = pd(23,13) - rrt(465) * density(23) 
  pd(23,23) = pd(23,23) - rrt(465) * density(13) 
  pd(12,13) = pd(12,13) + rrt(466) * density(23) 
  pd(12,23) = pd(12,23) + rrt(466) * density(13) 
  pd(13,13) = pd(13,13) - rrt(466) * density(23) 
  pd(13,23) = pd(13,23) - rrt(466) * density(13) 
  pd(23,13) = pd(23,13) - rrt(466) * density(23) 
  pd(23,23) = pd(23,23) - rrt(466) * density(13) 
  pd(13,13) = pd(13,13) - rrt(467) * density(23) 
  pd(13,23) = pd(13,23) - rrt(467) * density(13) 
  pd(23,13) = pd(23,13) - rrt(467) * density(23) 
  pd(23,23) = pd(23,23) - rrt(467) * density(13) 
  pd(24,13) = pd(24,13) + rrt(467) * density(23) 
  pd(24,23) = pd(24,23) + rrt(467) * density(13) 
  pd(38,13) = pd(38,13) + rrt(467) * density(23) 
  pd(38,23) = pd(38,23) + rrt(467) * density(13) 
  pd(09,23) = pd(09,23) + rrt(468) 
  pd(13,23) = pd(13,23) + rrt(468) 
  pd(23,23) = pd(23,23) - rrt(468) 
  pd(23,23) = pd(23,23) - rrt(469) 
  pd(24,23) = pd(24,23) + rrt(469) 
  pd(35,23) = pd(35,23) + rrt(469) 
  pd(03,09) = pd(03,09) + rrt(470) * density(19) 
  pd(03,19) = pd(03,19) + rrt(470) * density(09) 
  pd(09,09) = pd(09,09) - rrt(470) * density(19) 
  pd(09,19) = pd(09,19) - rrt(470) * density(09) 
  pd(19,09) = pd(19,09) - rrt(470) * density(19) 
  pd(19,19) = pd(19,19) - rrt(470) * density(09) 
  pd(50,09) = pd(50,09) + rrt(470) * density(19) 
  pd(50,19) = pd(50,19) + rrt(470) * density(09) 
  pd(03,09) = pd(03,09) + rrt(471) * density(09) * 2.0d0
  pd(09,09) = pd(09,09) - rrt(471) * density(09) * 4.0d0
  pd(23,09) = pd(23,09) + rrt(471) * density(09) * 2.0d0
  pd(03,09) = pd(03,09) + rrt(472) 
  pd(09,09) = pd(09,09) - rrt(472) 
  pd(13,09) = pd(13,09) + rrt(472) 
  pd(02,09) = pd(02,09) + rrt(473) * density(13) 
  pd(02,13) = pd(02,13) + rrt(473) * density(09) 
  pd(03,09) = pd(03,09) + rrt(473) * density(13) 
  pd(03,13) = pd(03,13) + rrt(473) * density(09) 
  pd(09,09) = pd(09,09) - rrt(473) * density(13) 
  pd(09,13) = pd(09,13) - rrt(473) * density(09) 
  pd(13,09) = pd(13,09) - rrt(473) * density(13) 
  pd(13,13) = pd(13,13) - rrt(473) * density(09) 
  pd(09,09) = pd(09,09) - rrt(474) * density(13) 
  pd(09,13) = pd(09,13) - rrt(474) * density(09) 
  pd(13,09) = pd(13,09) - rrt(474) * density(13) 
  pd(13,13) = pd(13,13) - rrt(474) * density(09) 
  pd(23,09) = pd(23,09) + rrt(474) * density(13) 
  pd(23,13) = pd(23,13) + rrt(474) * density(09) 
  pd(09,09) = pd(09,09) - rrt(475) 
  pd(24,09) = pd(24,09) + rrt(475) 
  pd(50,09) = pd(50,09) + rrt(475) 
  pd(09,09) = pd(09,09) - rrt(476) * density(24) 
  pd(09,24) = pd(09,24) - rrt(476) * density(09) 
  pd(24,09) = pd(24,09) - rrt(476) * density(24) 
  pd(24,24) = pd(24,24) - rrt(476) * density(09) 
  pd(39,09) = pd(39,09) + rrt(476) * density(24) 
  pd(39,24) = pd(39,24) + rrt(476) * density(09) 
  pd(02,02) = pd(02,02) - rrt(477) * density(03) 
  pd(02,03) = pd(02,03) - rrt(477) * density(02) 
  pd(03,02) = pd(03,02) - rrt(477) * density(03) 
  pd(03,03) = pd(03,03) - rrt(477) * density(02) 
  pd(09,02) = pd(09,02) + rrt(477) * density(03) 
  pd(09,03) = pd(09,03) + rrt(477) * density(02) 
  pd(13,02) = pd(13,02) + rrt(477) * density(03) 
  pd(13,03) = pd(13,03) + rrt(477) * density(02) 
  pd(03,03) = pd(03,03) - rrt(478) * density(13) 
  pd(03,13) = pd(03,13) - rrt(478) * density(03) 
  pd(09,03) = pd(09,03) + rrt(478) * density(13) 
  pd(09,13) = pd(09,13) + rrt(478) * density(03) 
  pd(13,03) = pd(13,03) - rrt(478) * density(13) 
  pd(13,13) = pd(13,13) - rrt(478) * density(03) 
  pd(03,03) = pd(03,03) - rrt(479) 
  pd(19,03) = pd(19,03) + rrt(479) 
  pd(24,03) = pd(24,03) + rrt(479) 
  pd(35,39) = pd(35,39) + rrt(480) 
  pd(38,39) = pd(38,39) + rrt(480) 
  pd(39,39) = pd(39,39) - rrt(480) 
  pd(10,10) = pd(10,10) - rrt(481) * density(39) 
  pd(10,39) = pd(10,39) - rrt(481) * density(10) 
  pd(23,10) = pd(23,10) + rrt(481) * density(39) 
  pd(23,39) = pd(23,39) + rrt(481) * density(10) 
  pd(35,10) = pd(35,10) + rrt(481) * density(39) 
  pd(35,39) = pd(35,39) + rrt(481) * density(10) 
  pd(39,10) = pd(39,10) - rrt(481) * density(39) 
  pd(39,39) = pd(39,39) - rrt(481) * density(10) 
  pd(09,39) = pd(09,39) + rrt(482) 
  pd(24,39) = pd(24,39) + rrt(482) 
  pd(39,39) = pd(39,39) - rrt(482) 
  pd(11,11) = pd(11,11) - rrt(483) 
  pd(24,11) = pd(24,11) + rrt(483) 
  pd(39,11) = pd(39,11) + rrt(483) 
  pd(02,02) = pd(02,02) - rrt(484) 
  pd(13,02) = pd(13,02) + rrt(484) * 2.0d0
  pd(02,13) = pd(02,13) + rrt(485) * density(13) * 2.0d0
  pd(13,13) = pd(13,13) - rrt(485) * density(13) * 4.0d0
  if( ldensity_constant ) then
    do i = 1, species_max
      if( density_constant(i) ) pd(i,:) = 0.0d0
    enddo
  endif
  if( lgas_heating ) then
    pd(51,1) = eV_to_K * ZDPlasKin_cfg(11)
    pd(51,:) = pd(51,:) * ZDPlasKin_cfg(13)
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
  DOUBLE PRECISION, PARAMETER :: R = 8.314D-3, FB = 7.5D-3, FI = 7.5D-3, FG = 1.0D0, FRI = 7.5D-8, FCH3 = 1.0D6, FC2H4 = 1.0D-8
  call ZDPlasKin_bolsig_rates()
  Tgas = ZDPlasKin_cfg(1)
  Te  = ZDPlasKin_cfg(4)
  rrt(001) = FB*bolsig_rates(bolsig_pointer(1))
  rrt(002) = FB*bolsig_rates(bolsig_pointer(2))
  rrt(003) = FB*bolsig_rates(bolsig_pointer(3))
  rrt(004) = FB*bolsig_rates(bolsig_pointer(4))
  rrt(005) = FB*bolsig_rates(bolsig_pointer(5))
  rrt(006) = FB*bolsig_rates(bolsig_pointer(6))
  rrt(007) = FB*bolsig_rates(bolsig_pointer(7))
  rrt(008) = FB*bolsig_rates(bolsig_pointer(8))
  rrt(009) = FB*bolsig_rates(bolsig_pointer(9))
  rrt(010) = FB*bolsig_rates(bolsig_pointer(10))
  rrt(011) = FB*bolsig_rates(bolsig_pointer(11))
  rrt(012) = FB*bolsig_rates(bolsig_pointer(12))
  rrt(013) = FB*bolsig_rates(bolsig_pointer(13))
  rrt(014) = FB*bolsig_rates(bolsig_pointer(14))
  rrt(015) = FB*bolsig_rates(bolsig_pointer(15))
  rrt(016) = FB*bolsig_rates(bolsig_pointer(16))
  rrt(017) = FB*bolsig_rates(bolsig_pointer(17))
  rrt(018) = FB*bolsig_rates(bolsig_pointer(18))
  rrt(019) = FB*bolsig_rates(bolsig_pointer(19))
  rrt(020) = FB*bolsig_rates(bolsig_pointer(20))
  rrt(021) = FB*bolsig_rates(bolsig_pointer(21))
  rrt(022) = FB*bolsig_rates(bolsig_pointer(22))
  rrt(023) = FB*bolsig_rates(bolsig_pointer(23))
  rrt(024) = FB*bolsig_rates(bolsig_pointer(24))
  rrt(025) = FB*bolsig_rates(bolsig_pointer(25))
  rrt(026) = FB*bolsig_rates(bolsig_pointer(26))
  rrt(027) = FB*bolsig_rates(bolsig_pointer(27))
  rrt(028) = FB*bolsig_rates(bolsig_pointer(28))
  rrt(029) = FB*bolsig_rates(bolsig_pointer(29))
  rrt(030) = FB*bolsig_rates(bolsig_pointer(30))
  rrt(031) = FB*bolsig_rates(bolsig_pointer(31))
  rrt(032) = FB*bolsig_rates(bolsig_pointer(32))
  rrt(033) = FB*bolsig_rates(bolsig_pointer(33))
  rrt(034) = FB*bolsig_rates(bolsig_pointer(34))
  rrt(035) = FB*bolsig_rates(bolsig_pointer(35))
  rrt(036) = FB*bolsig_rates(bolsig_pointer(36))
  rrt(037) = FB*bolsig_rates(bolsig_pointer(37))
  rrt(038) = FB*bolsig_rates(bolsig_pointer(38))
  rrt(039) = FB*bolsig_rates(bolsig_pointer(39))
  rrt(040) = FB*bolsig_rates(bolsig_pointer(40))
  rrt(041) = FB*bolsig_rates(bolsig_pointer(41))
  rrt(042) = FB*bolsig_rates(bolsig_pointer(42))
  rrt(043) = FB*bolsig_rates(bolsig_pointer(43))
  rrt(044) = FB*bolsig_rates(bolsig_pointer(44))
  rrt(045) = FB*bolsig_rates(bolsig_pointer(45))
  rrt(046) = FB*bolsig_rates(bolsig_pointer(46))
  rrt(047) = FB*bolsig_rates(bolsig_pointer(47))
  rrt(048) = FB*bolsig_rates(bolsig_pointer(48))
  rrt(049) = FB*bolsig_rates(bolsig_pointer(49))
  rrt(050) = FB*bolsig_rates(bolsig_pointer(50))
  rrt(051) = FB*bolsig_rates(bolsig_pointer(51))
  rrt(052) = FB*bolsig_rates(bolsig_pointer(52))
  rrt(053) = FB*bolsig_rates(bolsig_pointer(53))
  rrt(054) = FB*bolsig_rates(bolsig_pointer(54))
  rrt(055) = FB*bolsig_rates(bolsig_pointer(55))
  rrt(056) = FB*bolsig_rates(bolsig_pointer(56))
  rrt(057) = FB*bolsig_rates(bolsig_pointer(57))
  rrt(058) = FB*bolsig_rates(bolsig_pointer(58))
  rrt(059) = FB*bolsig_rates(bolsig_pointer(59))
  rrt(060) = FB*bolsig_rates(bolsig_pointer(60))
  rrt(061) = FB*bolsig_rates(bolsig_pointer(61))
  rrt(062) = FB*bolsig_rates(bolsig_pointer(62))
  rrt(063) = FB*bolsig_rates(bolsig_pointer(63))
  rrt(064) = FB*bolsig_rates(bolsig_pointer(64))
  rrt(065) = FB*bolsig_rates(bolsig_pointer(65))
  rrt(066) = FB*bolsig_rates(bolsig_pointer(66))
  rrt(067) = FB*bolsig_rates(bolsig_pointer(67))
  rrt(068) = FB*bolsig_rates(bolsig_pointer(68))
  rrt(069) = FB*bolsig_rates(bolsig_pointer(69))
  rrt(070) = FB*bolsig_rates(bolsig_pointer(70))
  rrt(071) = FB*bolsig_rates(bolsig_pointer(71))
  rrt(072) = FB*bolsig_rates(bolsig_pointer(72))
  rrt(073) = FB*bolsig_rates(bolsig_pointer(73))
  rrt(074) = FB*bolsig_rates(bolsig_pointer(74))
  rrt(075) = FB*bolsig_rates(bolsig_pointer(75))
  rrt(076) = FB*bolsig_rates(bolsig_pointer(76))
  rrt(077) = FB*bolsig_rates(bolsig_pointer(77))
  rrt(078) = FB*bolsig_rates(bolsig_pointer(78))
  rrt(079) = FB*bolsig_rates(bolsig_pointer(79))
  rrt(080) = FB*bolsig_rates(bolsig_pointer(80))
  rrt(081) = FB*bolsig_rates(bolsig_pointer(81))
  rrt(082) = FB*bolsig_rates(bolsig_pointer(82))
  rrt(083) = FB*bolsig_rates(bolsig_pointer(83))
  rrt(084) = FB*bolsig_rates(bolsig_pointer(84))
  rrt(085) = FB*bolsig_rates(bolsig_pointer(85))
  rrt(086) = FB*bolsig_rates(bolsig_pointer(86))
  rrt(087) = FB*bolsig_rates(bolsig_pointer(87))
  rrt(088) = FB*bolsig_rates(bolsig_pointer(88))
  rrt(089) = FB*bolsig_rates(bolsig_pointer(89))
  rrt(090) = FB*bolsig_rates(bolsig_pointer(90))
  rrt(091) = FB*bolsig_rates(bolsig_pointer(91))
  rrt(092) = FB*bolsig_rates(bolsig_pointer(92))
  rrt(093) = FB*bolsig_rates(bolsig_pointer(93))
  rrt(094) = FB*bolsig_rates(bolsig_pointer(94))
  rrt(095) = FB*bolsig_rates(bolsig_pointer(95))
  rrt(096) = FB*bolsig_rates(bolsig_pointer(96))
  rrt(097) = FB*bolsig_rates(bolsig_pointer(97))
  rrt(098) = FB*bolsig_rates(bolsig_pointer(98))
  rrt(099) = FB*bolsig_rates(bolsig_pointer(99))
  rrt(100) = FB*bolsig_rates(bolsig_pointer(100))
  rrt(101) = FB*bolsig_rates(bolsig_pointer(101))
  rrt(102) = FB*bolsig_rates(bolsig_pointer(102))
  rrt(103) = FB*bolsig_rates(bolsig_pointer(103))
  rrt(104) = FB*bolsig_rates(bolsig_pointer(104))
  rrt(105) = FB*bolsig_rates(bolsig_pointer(105))
  rrt(106) = FB*bolsig_rates(bolsig_pointer(106))
  rrt(107) = FB*bolsig_rates(bolsig_pointer(107))
  rrt(108) = FB*bolsig_rates(bolsig_pointer(108))
  rrt(109) = FB*bolsig_rates(bolsig_pointer(109))
  rrt(110) = FB*bolsig_rates(bolsig_pointer(110))
  rrt(111) = FB*bolsig_rates(bolsig_pointer(111))
  rrt(112) = FB*bolsig_rates(bolsig_pointer(112))
  rrt(113) = FB*bolsig_rates(bolsig_pointer(113))
  rrt(114) = FB*bolsig_rates(bolsig_pointer(114))
  rrt(115) = FB*bolsig_rates(bolsig_pointer(115))
  rrt(116) = FB*bolsig_rates(bolsig_pointer(116))
  rrt(117) = FB*bolsig_rates(bolsig_pointer(117))
  rrt(118) = FB*bolsig_rates(bolsig_pointer(118))
  rrt(119) = FB*bolsig_rates(bolsig_pointer(119))
  rrt(120) = FB*bolsig_rates(bolsig_pointer(120))
  rrt(121) = FB*bolsig_rates(bolsig_pointer(121))
  rrt(122) = FB*bolsig_rates(bolsig_pointer(122))
  rrt(123) = FB*bolsig_rates(bolsig_pointer(123))
  rrt(124) = FB*bolsig_rates(bolsig_pointer(124))
  rrt(125) = FB*bolsig_rates(bolsig_pointer(125))
  rrt(126) = FB*bolsig_rates(bolsig_pointer(126))
  rrt(127) = FB*bolsig_rates(bolsig_pointer(127))
  rrt(128) = FB*bolsig_rates(bolsig_pointer(128))
  rrt(129) = FB*bolsig_rates(bolsig_pointer(129))
  rrt(130) = FB*bolsig_rates(bolsig_pointer(130))
  rrt(131) = FB*bolsig_rates(bolsig_pointer(131))
  rrt(132) = FB*bolsig_rates(bolsig_pointer(132))
  rrt(133) = FB*bolsig_rates(bolsig_pointer(133))
  rrt(134) = FB*bolsig_rates(bolsig_pointer(134))
  rrt(135) = FB*bolsig_rates(bolsig_pointer(135))
  rrt(136) = FB*bolsig_rates(bolsig_pointer(136))
  rrt(137) = FB*bolsig_rates(bolsig_pointer(137))
  rrt(138) = FB*bolsig_rates(bolsig_pointer(138))
  rrt(139) = FB*bolsig_rates(bolsig_pointer(139))
  rrt(140) = FB*bolsig_rates(bolsig_pointer(140))
  rrt(141) = FB*bolsig_rates(bolsig_pointer(141))
  rrt(142) = FB*bolsig_rates(bolsig_pointer(142))
  rrt(143) = FB*bolsig_rates(bolsig_pointer(143))
  rrt(144) = FB*bolsig_rates(bolsig_pointer(144))
  rrt(145) = FB*bolsig_rates(bolsig_pointer(145))
  rrt(146) = FB*bolsig_rates(bolsig_pointer(146))
  rrt(147) = FB*bolsig_rates(bolsig_pointer(147))
  rrt(148) = FB*bolsig_rates(bolsig_pointer(148))
  rrt(149) = FB*bolsig_rates(bolsig_pointer(149))
  rrt(150) = FB*bolsig_rates(bolsig_pointer(150))
  rrt(151) = FB*bolsig_rates(bolsig_pointer(151))
  rrt(152) = FB*bolsig_rates(bolsig_pointer(152))
  rrt(153) = FB*bolsig_rates(bolsig_pointer(153))
  rrt(154) = FB*bolsig_rates(bolsig_pointer(154))
  rrt(155) = FB*bolsig_rates(bolsig_pointer(155))
  rrt(156) = FB*bolsig_rates(bolsig_pointer(156))
  rrt(157) = FB*bolsig_rates(bolsig_pointer(157))
  rrt(158) = FB*bolsig_rates(bolsig_pointer(158))
  rrt(159) = FB*bolsig_rates(bolsig_pointer(159))
  rrt(160) = FB*bolsig_rates(bolsig_pointer(160))
  rrt(161) = FB*bolsig_rates(bolsig_pointer(161))
  rrt(162) = FB*bolsig_rates(bolsig_pointer(162))
  rrt(163) = FB*bolsig_rates(bolsig_pointer(163))
  rrt(164) = FB*bolsig_rates(bolsig_pointer(164))
  rrt(165) = FB*bolsig_rates(bolsig_pointer(165))
  rrt(166) = FB*bolsig_rates(bolsig_pointer(166))
  rrt(167) = FB*bolsig_rates(bolsig_pointer(167))
  rrt(168) = FB*bolsig_rates(bolsig_pointer(168))
  rrt(169) = FB*bolsig_rates(bolsig_pointer(169))
  rrt(170) = FB*bolsig_rates(bolsig_pointer(170))
  rrt(171) = FB*bolsig_rates(bolsig_pointer(171))
  rrt(172) = FB*bolsig_rates(bolsig_pointer(172))
  rrt(173) = FB*bolsig_rates(bolsig_pointer(173))
  rrt(174) = FB*bolsig_rates(bolsig_pointer(174))
  rrt(175) = FB*bolsig_rates(bolsig_pointer(175))
  rrt(176) = FB*bolsig_rates(bolsig_pointer(176))
  rrt(177) = FB*bolsig_rates(bolsig_pointer(177))
  rrt(178) = FB*bolsig_rates(bolsig_pointer(178))
  rrt(179) = FB*bolsig_rates(bolsig_pointer(179))
  rrt(180) = FB*bolsig_rates(bolsig_pointer(180))
  rrt(181) = FB*bolsig_rates(bolsig_pointer(181))
  rrt(182) = FB*bolsig_rates(bolsig_pointer(182))
  rrt(183) = FB*bolsig_rates(bolsig_pointer(183))
  rrt(184) = FB*bolsig_rates(bolsig_pointer(184))
  rrt(185) = FB*bolsig_rates(bolsig_pointer(185))
  rrt(186) = FB*bolsig_rates(bolsig_pointer(186))
  rrt(187) = FB*bolsig_rates(bolsig_pointer(187))
  rrt(188) = FB*bolsig_rates(bolsig_pointer(188))
  rrt(189) = FB*bolsig_rates(bolsig_pointer(189))
  rrt(190) = FB*bolsig_rates(bolsig_pointer(190))
  rrt(191) = FB*bolsig_rates(bolsig_pointer(191))
  rrt(192) = FB*bolsig_rates(bolsig_pointer(192))
  rrt(193) = FB*bolsig_rates(bolsig_pointer(193))
  rrt(194) = FB*bolsig_rates(bolsig_pointer(194))
  rrt(195) = FB*bolsig_rates(bolsig_pointer(195))
  rrt(196) = FB*bolsig_rates(bolsig_pointer(196))
  rrt(197) = FB*bolsig_rates(bolsig_pointer(197))
  rrt(198) = FB*bolsig_rates(bolsig_pointer(198))
  rrt(199) = FB*bolsig_rates(bolsig_pointer(199))
  rrt(200) = FB*bolsig_rates(bolsig_pointer(200))
  rrt(201) = FB*bolsig_rates(bolsig_pointer(201))
  rrt(202) = FB*bolsig_rates(bolsig_pointer(202))
  rrt(203) = FB*bolsig_rates(bolsig_pointer(203))
  rrt(204) = FB*bolsig_rates(bolsig_pointer(204))
  rrt(205) = FB*bolsig_rates(bolsig_pointer(205))
  rrt(206) = FB*bolsig_rates(bolsig_pointer(206))
  rrt(207) = FB*bolsig_rates(bolsig_pointer(207))
  rrt(208) = FB*bolsig_rates(bolsig_pointer(208))
  rrt(209) = FB*bolsig_rates(bolsig_pointer(209))
  rrt(210) = FB*bolsig_rates(bolsig_pointer(210))
  rrt(211) = FB*bolsig_rates(bolsig_pointer(211))
  rrt(212) = FB*bolsig_rates(bolsig_pointer(212))
  rrt(213) = FB*bolsig_rates(bolsig_pointer(213))
  rrt(214) = FB*bolsig_rates(bolsig_pointer(214))
  rrt(215) = FB*bolsig_rates(bolsig_pointer(215))
  rrt(216) = FB*bolsig_rates(bolsig_pointer(216))
  rrt(217) = FB*bolsig_rates(bolsig_pointer(217))
  rrt(218) = FB*bolsig_rates(bolsig_pointer(218))
  rrt(219) = FB*bolsig_rates(bolsig_pointer(219))
  rrt(220) = FB*bolsig_rates(bolsig_pointer(220))
  rrt(221) = FB*bolsig_rates(bolsig_pointer(221))
  rrt(222) = FB*bolsig_rates(bolsig_pointer(222))
  rrt(223) = FB*bolsig_rates(bolsig_pointer(223))
  rrt(224) = FB*bolsig_rates(bolsig_pointer(224))
  rrt(225) = FB*bolsig_rates(bolsig_pointer(225))
  rrt(226) = FB*bolsig_rates(bolsig_pointer(226))
  rrt(227) = FB*bolsig_rates(bolsig_pointer(227))
  rrt(228) = FB*bolsig_rates(bolsig_pointer(228))
  rrt(229) = FB*bolsig_rates(bolsig_pointer(229))
  rrt(230) = FB*bolsig_rates(bolsig_pointer(230))
  rrt(231) = FB*bolsig_rates(bolsig_pointer(231))
  rrt(232) = FB*bolsig_rates(bolsig_pointer(232))
  rrt(233) = FB*bolsig_rates(bolsig_pointer(233))
  rrt(234) = FB*bolsig_rates(bolsig_pointer(234))
  rrt(235) = FB*bolsig_rates(bolsig_pointer(235))
  rrt(236) = FB*bolsig_rates(bolsig_pointer(236))
  rrt(237) = FB*bolsig_rates(bolsig_pointer(237))
  rrt(238) = FB*bolsig_rates(bolsig_pointer(238))
  rrt(239) = FB*bolsig_rates(bolsig_pointer(239))
  rrt(240) = FRI
  rrt(241) = FRI
  rrt(242) = FRI
  rrt(243) = FRI
  rrt(244) = FRI
  rrt(245) = FRI
  rrt(246) = FRI
  rrt(247) = FRI
  rrt(248) = FRI
  rrt(249) = FRI
  rrt(250) = FRI
  rrt(251) = FRI
  rrt(252) = FRI
  rrt(253) = FRI
  rrt(254) = FRI
  rrt(255) = FRI
  rrt(256) = FI*2.57D-07*(300./TGAS)**0.3
  rrt(257) = FI*6.61D-08*(300./TGAS)**0.3
  rrt(258) = FI*1.18D-08*(300./TGAS)**0.5
  rrt(259) = FI*2.42D-08*(300./TGAS)**0.5
  rrt(260) = FI*1.41D-08*(300./TGAS)**0.5
  rrt(261) = FI*2.25D-08*(300./TGAS)**0.5
  rrt(262) = FI*7.88D-09*(300./TGAS)**0.5
  rrt(263) = FI*1.00D-08*(300./TGAS)**0.5
  rrt(264) = FI*2.19D-08*(300./TGAS)**0.71
  rrt(265) = FI*3.36D-08*(300./TGAS)**0.71
  rrt(266) = FI*7.70D-09*(300./TGAS)**0.71
  rrt(267) = FI*1.92D-08*(300./TGAS)**0.71
  rrt(268) = FI*1.60D-08*(300./TGAS)**0.71
  rrt(269) = FI*8.98D-09*(300./TGAS)**0.71
  rrt(270) = FI*9.62D-09*(300./TGAS)**0.71
  rrt(271) = FI*8.29D-09*(300./TGAS)**0.71
  rrt(272) = FI*3.43D-08*(300./TGAS)**0.71
  rrt(273) = FI*1.34D-08*(300./TGAS)**0.71
  rrt(274) = FI*4.87D-09*(300./TGAS)**0.71
  rrt(275) = FI*3.17D21/(6.022D23*TE**4.5)
  rrt(276) = rrt(275)
  rrt(277) = FG*2.11D-09
  rrt(278) = FG*9.60D-10
  rrt(279) = FG*6.90D-10
  rrt(280) = FG*2.25D-10
  rrt(281) = FC2H4*FG*1.50D-9
  rrt(282) = FC2H4*FG*1.60D-9
  rrt(283) = FG*1.50D-10
  rrt(284) = FG*1.50D-9
  rrt(285) = FG*1.91D-9
  rrt(286) = FC2H4*FG*4.23D-10
  rrt(287) = FC2H4*FG*1.38D-9
  rrt(288) = FC2H4*FG*1.23D-9
  rrt(289) = FC2H4*FG*1.13D-9
  rrt(290) = FG*3.30D-11
  rrt(291) = FG*1.00D-11
  rrt(292) = FG*1.36D-10
  rrt(293) = FG*1.20D-9
  rrt(294) = FG*9.90D-10
  rrt(295) = FG*7.10D-10
  rrt(296) = FG*1.48D-9
  rrt(297) = FC2H4*FG*3.50D-10
  rrt(298) = FG*3.00D-10
  rrt(299) = FG*1.38D-10
  rrt(300) = FG*3.60D-10
  rrt(301) = FG*8.40D-10
  rrt(302) = FG*2.31D-10
  rrt(303) = FG*3.97D-10
  rrt(304) = FG*1.60D-9
  rrt(305) = FG*6.50D-11
  rrt(306) = FG*1.09D-9
  rrt(307) = FG*1.43D-10
  rrt(308) = rrt(293)
  rrt(309) = FC2H4*FG*1.15D-9
  rrt(310) = FC2H4*FG*2.47D-10
  rrt(311) = FG*1.00D-10
  rrt(312) = rrt(291)
  rrt(313) = FG*5.00D-10
  rrt(314) = rrt(313)
  rrt(315) = rrt(298)
  rrt(316) = FG*2.91D-10
  rrt(317) = FC2H4*FG*8.90D-10
  rrt(318) = FG*3.30D-10
  rrt(319) = FG*6.80D-11
  rrt(320) = FG*4.10D-9
  rrt(321) = FG*1.31D-10
  rrt(322) = FG*2.48D-10
  rrt(323) = FC2H4*FG*4.14D-10
  rrt(324) = rrt(318)
  rrt(325) = rrt(291)
  rrt(326) = FG*2.40D-9
  rrt(327) = FCH3*FG*2.10D-9
  rrt(328) = FG*1.70D-9
  rrt(329) = rrt(293)
  rrt(330) = rrt(326)
  rrt(331) = FG*1.40D-9
  rrt(332) = rrt(309)
  rrt(333) = rrt(309)
  rrt(334) = FG*2.00D-9
  rrt(335) = rrt(328)
  rrt(336) = FC2H4*FG*3.50D-9
  rrt(337) = FG*1.14D-10
  rrt(338) = rrt(331)
  rrt(339) = FG*2.30D-9
  rrt(340) = FG*1.00D-9
  rrt(341) = rrt(340)
  rrt(342) = rrt(295)
  rrt(343) = rrt(295)
  rrt(344) = FG*2.94D-10
  rrt(345) = FG*1.37D-9
  rrt(346) = FG*2.35D-9
  rrt(347) = FG*6.86D-10
  rrt(348) = FG*1.96D-10
  rrt(349) = FC2H4*FG*2.21D-9
  rrt(350) = FC2H4*FG*1.81D-9
  rrt(351) = FC2H4*FG*8.82D-10
  rrt(352) = FC2H4*FG*4.80D-10
  rrt(353) = FC2H4*FG*4.82D-9
  rrt(354) = FG*2.10D-9
  rrt(355) = FG*6.39D-10
  rrt(356) = rrt(284)
  rrt(357) = rrt(339)
  rrt(358) = FCH3*FG*3.40D-9
  rrt(359) = rrt(331)
  rrt(360) = rrt(331)
  rrt(361) = FG*1.90D-9
  rrt(362) = FG*1.30D-9
  rrt(363) = rrt(331)
  rrt(364) = FG*2.80D-9
  rrt(365) = FG*1.65D-9
  rrt(366) = FG*3.06D-9
  rrt(367) = FC2H4*FG*1.00D-9
  rrt(368) = FC2H4*FG*3.00D-9
  rrt(369) = rrt(367)
  rrt(370) = rrt(334)
  rrt(371) = rrt(334)
  rrt(372) = FC2H4*FG*5.40D-10
  rrt(373) = FG*4.08D-18*TGAS**2*EXP(-4163./TGAS)
  rrt(374) = FG*9.97D-11
  rrt(375) = FG*2.51D-15*(TGAS/298.)**4.14*EXP(-52.55/(R*TGAS))
  rrt(376) = FG*4.26D-15*(TGAS/298.)**4.02*EXP(-22.86/(R*TGAS))
  rrt(377) = FG*3.01D-12*EXP(-2.08/(R*TGAS))
  rrt(378) = FG*3.54D-16*(TGAS/298.)**4.02*EXP(-45.48/(R*TGAS))
  rrt(379) = FG*1.71D-14*(TGAS/298.)**3.40*EXP(-97.28/(R*TGAS))
  rrt(380) = FG*9.86D-13*(TGAS/298.)**3.00*EXP(-36.67/(R*TGAS))
  rrt(381) = FCH3*FG*1.33D-10*EXP(-167.00/(R*TGAS))
  rrt(382) = FG*1.90D-12
  rrt(383) = FG*2.40D16*EXP(-52800./TGAS)
  rrt(384) = FCH3*FG*1.69D-8*EXP(-379./(R*TGAS))
  rrt(385) = FCH3*FG*6.97D-9*EXP(-345./(R*TGAS))
  rrt(386) = FCH3*FG*3.00D-44*TGAS**9.10
  rrt(387) = FG*3.32D-10*EXP(-45.98/(R*TGAS))
  rrt(388) = FG*3.01D-11
  rrt(389) = rrt(388)
  rrt(390) = rrt(388)
  rrt(391) = FG*1.61D-15*(TGAS/298.)**3.65*EXP(-29.93/(R*TGAS))
  rrt(392) = rrt(388)
  rrt(393) = rrt(388)
  rrt(394) = FG*1.20D-12*EXP(-25.94/(R*TGAS))
  rrt(395) = FG*8.30D-19*TGAS**2.00*EXP(-3938.65/TGAS)
  rrt(396) = FG*1.00D-11*EXP(7.48/(R*TGAS))
  rrt(397) = FG*6.64D-9*EXP(-348./(R*TGAS))
  rrt(398) = FG*2.66D-9*EXP(-6011.07/TGAS)
  rrt(399) = FG*9.96D-10
  rrt(400) = FG*3.00D-11
  rrt(401) = FG*1.14D-29
  rrt(402) = FG*1.79D-10*EXP(-1565.17/TGAS)
  rrt(403) = FCH3*FG*4.98D-11
  rrt(404) = FG*6.64D-11
  rrt(405) = FG*3.29D-12*TGAS**0.43*EXP(186.21/TGAS)
  rrt(406) = FG*8.30D-11
  rrt(407) = FG*1.46D-13*(TGAS/298.)**3.30*EXP(-43.90/(R*TGAS))
  rrt(408) = FG*1.19D-15*(TGAS/298.)**3.82*EXP(-37.83/(R*TGAS))
  rrt(409) = FG*5.71D-14*(TGAS/298.)**3.30*EXP(-83.06/(R*TGAS))
  rrt(410) = FG*1.23D-11*(TGAS/298.)**1.50*EXP(-31.01/(R*TGAS))
  rrt(411) = FG*8.97D-20*EXP(-48.64/(R*TGAS))
  rrt(412) = FG*1.80D21*TGAS**(-1.24)*EXP(-45700./TGAS)
  rrt(413) = FG*1.79D-10*EXP(132.36/TGAS)
  rrt(414) = FG*9.00D-33*TGAS**6.43
  rrt(415) = FG*2.41D-12
  rrt(416) = FC2H4*FG*5.83D-14*(TGAS/298.)**3.13*EXP(-75.33/(R*TGAS))
  rrt(417) = FC2H4*FG*4.50D-13*EXP(-98.11/(R*TGAS))
  rrt(418) = FG*3.01D-12
  rrt(419) = FG*1.61D-15*(TGAS/298.)**3.65*EXP(-38.25/(R*TGAS))
  rrt(420) = FG*1.91D-12
  rrt(421) = rrt(415)
  rrt(422) = FG*1.69D-15*(TGAS/298.)**3.50*EXP(-35.34/(R*TGAS))
  rrt(423) = FG*4.12D-15*(TGAS/298.)**3.60*EXP(-27.77/(R*TGAS))
  rrt(424) = FG*5.99D-11
  rrt(425) = FG*3.32D-12
  rrt(426) = FG*8.65D-7*TGAS**(-0.99)*EXP(-795.17/TGAS)
  rrt(427) = FG*4.08D12*(TGAS/298.)**1.04*EXP(-154./(R*TGAS))
  rrt(428) = FG*4.42D-11
  rrt(429) = FC2H4*FG*9.00D-10*EXP(-62.36/(R*TGAS))
  rrt(430) = FC2H4*FG*9.68D-12*(TGAS/298.)**1.28*EXP(-5.40/(R*TGAS))
  rrt(431) = FC2H4*FG*9.55D-12
  rrt(432) = FC2H4*FG*4.30D-7*EXP(-404./(R*TGAS))
  rrt(433) = FC2H4*FG*9.60D-11*EXP(-216./(R*TGAS))
  rrt(434) = FC2H4*FG*4.00D-11*EXP(-286./(R*TGAS))
  rrt(435) = FC2H4*FG*1.00D-10*EXP(-316./(R*TGAS))
  rrt(436) = FC2H4*FG*8.00D-10*EXP(-299./(R*TGAS))
  rrt(437) = FG*4.23D-18*TGAS**1.60*EXP(-2868.65/TGAS)
  rrt(438) = FC2H4*FG*8.00D12*TGAS**0.44*EXP(-44675.39/TGAS)
  rrt(439) = FC2H4*FG*3.00D-14*(TGAS/298.)**2.48*EXP(-25.65/(R*TGAS))
  rrt(440) = FC2H4*FG*4.75D-16*EXP(-180./(R*TGAS))
  rrt(441) = FC2H4*FG*5.30D-12*EXP(-2660./TGAS)
  rrt(442) = FC2H4*FG*5.00D-14*EXP(-25.53/(R*TGAS))
  rrt(443) = FG*3.50D-11
  rrt(444) = rrt(407)
  rrt(445) = FG*2.01D-12
  rrt(446) = rrt(445)
  rrt(447) = FG*1.68D-15*(TGAS/298.)**3.50*EXP(-19.62/(R*TGAS))
  rrt(448) = FG*8.00D-12
  rrt(449) = FG*1.61D-13*(TGAS/298.)**2.63*EXP(-35.75/(R*TGAS))
  rrt(450) = FG*1.60D-10
  rrt(451) = FG*1.01D-11*TGAS**0.27*EXP(-140.92/TGAS)
  rrt(452) = FG*2.00D14*EXP(-20000./TGAS)
  rrt(453) = FC2H4*FG*9.30D-12*EXP(-1207.85/TGAS)
  rrt(454) = FC2H4*FG*5.00D-13*EXP(-163./(R*TGAS))
  rrt(455) = FC2H4*FG*4.00D-12*EXP(-272./(R*TGAS))
  rrt(456) = FG*1.58D-5*(TGAS/298.)**(-8.58)*EXP(-84.81/(R*TGAS))
  rrt(457) = FC2H4*FG*1.20D-12*EXP(-37.66/(R*TGAS))
  rrt(458) = rrt(409)
  rrt(459) = FG*2.19D-180*TGAS**2.54*EXP(-3400.10/TGAS)
  rrt(460) = FG*1.10D17*EXP(-42470./TGAS)
  rrt(461) = FG*2.81D-12
  rrt(462) = FG*1.69D-15*(TGAS/298.)**3.50*EXP(-27.77/(R*TGAS))
  rrt(463) = FG*2.41D-12*EXP(-0.55/(R*TGAS))
  rrt(464) = FG*3.19D-14*(TGAS/298.)**2.84*EXP(-38.25/(R*TGAS))
  rrt(465) = rrt(418)
  rrt(466) = FG*6.00D-11
  rrt(467) = FG*6.74D-18*TGAS**2.19*EXP(-447.91/TGAS)
  rrt(468) = FG*1.25D17*EXP(-237./(R*TGAS))
  rrt(469) = FG*16.0*(TGAS/298.)**(-10.00)*EXP(-150./(R*TGAS))
  rrt(470) = FC2H4*FG*6.71D-11*EXP(-196./(R*TGAS))
  rrt(471) = FG*4.20D-10*EXP(-231./(R*TGAS))
  rrt(472) = FG*2.50D15*EXP(-363./(R*TGAS))
  rrt(473) = FG*4.40D-13*(TGAS/298.)**2.50*EXP(-10.39/(R*TGAS))
  rrt(474) = FG*1.29D-11*(TGAS/298.)**0.51*EXP(-5.15/(R*TGAS))
  rrt(475) = FG*1.28D13*(TGAS/298.)**(-15.70)*EXP(-502./(R*TGAS))
  rrt(476) = FCH3*FG*1.27D-14*(TGAS/298.)**2.67*EXP(-28.66/(R*TGAS))
  rrt(477) = FG*1.39D-13*(TGAS/298.)**2.38*EXP(-79.49/(R*TGAS))
  rrt(478) = FG*78.80*(TGAS/298.)**(-11.76)*EXP(-98.53/(R*TGAS))
  rrt(479) = FC2H4*FG*1.26D13*EXP(-140./(R*TGAS))
  rrt(480) = FC2H4*FG*1.07D2*(TGAS/298.)**(-11.90)*EXP(-135./(R*TGAS))
  rrt(481) = rrt(388)
  rrt(482) = FG*7.71D13*(TGAS/298.)**0.77*EXP(-128./(R*TGAS))
  rrt(483) = FG*3.16D16*EXP(-331./(R*TGAS))
  rrt(484) = FG*1.88D-8*(TGAS/298.)**(-1.10)*EXP(-437./(R*TGAS))
  rrt(485) = FG*5.52D-30*TGAS**(-1.00)
  where( lreaction_block(:) ) rrt(:) = 0.0d0
  return
end subroutine ZDPlasKin_reac_rates
!-----------------------------------------------------------------------------------------------------------------------------------
!
! END OF FILE
!
!-----------------------------------------------------------------------------------------------------------------------------------
