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
! Wed Nov 13 17:53:43 2024
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
  integer, parameter :: species_max = 47, species_electrons = 1, species_length = 9, reactions_max = 577, reactions_length = 31
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
  integer, parameter, private               :: bolsig_species_max = 24, bolsig_species_length = 9, bolsig_rates_max = 174 
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
  /-1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0,&
    0, 0, 1, 1, 1/
  data species_name(1:species_max) &
  /"E        ","C        ","CH       ","CH^+     ","CH2      ","CH2^+    ","CH3      ","CH3^+    ","CH4      ","CH4(V13) ",&
   "CH4(V24) ","CH4^+    ","C2H2     ","C2H2(V13)","C2H2(V2) ","C2H2(V5) ","C2H2^+   ","C2H3     ","C2H3^+   ","C2H4     ",&
   "C2H4(V1) ","C2H4(V2) ","C2H4^+   ","C2H5     ","C2H5^+   ","C2H6     ","C2H6(V13)","C2H6(V24)","C2H6^+   ","C3H5     ",&
   "C3H6     ","C3H6(V)  ","C3H6^+   ","C3H7     ","C3H7^+   ","C3H8     ","C3H8(V1) ","C3H8(V2) ","C3H8^+   ","C4H9     ",&
   "C4H9H    ","C5H12    ","H        ","H2       ","H2^+     ","H3^+     ","H^+      "/
  data reaction_sign(1:54) &
  /"bolsig:CH4->CH4(V24)           ","bolsig:CH4->CH4(V13)           ","bolsig:C2H6->C2H6(V24)         ",&
   "bolsig:C2H6->C2H6(V13)         ","bolsig:C2H4->C2H4(V1)          ","bolsig:C2H4->C2H4(V2)          ",&
   "bolsig:C2H2->C2H2(V5)          ","bolsig:C2H2->C2H2(V2)          ","bolsig:C2H2->C2H2(V13)         ",&
   "bolsig:C3H8->C3H8(V1)          ","bolsig:C3H8->C3H8(V2)          ","bolsig:C3H6->C3H6(V)           ",&
   "bolsig:CH4->CH3H               ","bolsig:CH4(V24)->CH3H          ","bolsig:CH4(V13)->CH3H          ",&
   "bolsig:CH4->CH2H2              ","bolsig:CH4(V24)->CH2H2         ","bolsig:CH4(V13)->CH2H2         ",&
   "bolsig:CH4->CHH2H              ","bolsig:CH4(V24)->CHH2H         ","bolsig:CH4(V13)->CHH2H         ",&
   "bolsig:CH4->CH2H2              ","bolsig:CH4(V24)->CH2H2         ","bolsig:CH4(V13)->CH2H2         ",&
   "bolsig:CH3->CH2H               ","bolsig:CH3->CHH2               ","bolsig:CH3->CH2H               ",&
   "bolsig:CH2->CHH                ","bolsig:CH2->CH2                ","bolsig:CH2->CHH                ",&
   "bolsig:CH->CH                  ","bolsig:CH4->CH4^+              ","bolsig:CH4(V24)->CH4^+         ",&
   "bolsig:CH4(V13)->CH4^+         ","bolsig:CH4->HCH3^+             ","bolsig:CH4(V24)->HCH3^+        ",&
   "bolsig:CH4(V13)->HCH3^+        ","bolsig:CH4->H2CH2^+            ","bolsig:CH4(V24)->H2CH2^+       ",&
   "bolsig:CH4(V13)->H2CH2^+       ","bolsig:CH4->H2HCH^+            ","bolsig:CH4(V24)->H2HCH^+       ",&
   "bolsig:CH4(V13)->H2HCH^+       ","bolsig:CH3->CH3^+              ","bolsig:CH3->HCH2^+             ",&
   "bolsig:CH3->H2CH^+             ","bolsig:CH2->CH2^+              ","bolsig:CH2->HCH^+              ",&
   "bolsig:CH->CH^+                ","bolsig:C2H6->C2H5H             ","bolsig:C2H6(V24)->C2H5H        ",&
   "bolsig:C2H6(V13)->C2H5H        ","bolsig:C2H6->C2H4H2            ","bolsig:C2H6(V24)->C2H4H2       "/
  data reaction_sign(55:108) &
  /"bolsig:C2H6(V13)->C2H4H2       ","bolsig:C2H6->C2H3H2H           ","bolsig:C2H6(V24)->C2H3H2H      ",&
   "bolsig:C2H6(V13)->C2H3H2H      ","bolsig:C2H6->C2H2H2H2          ","bolsig:C2H6(V24)->C2H2H2H2     ",&
   "bolsig:C2H6(V13)->C2H2H2H2     ","bolsig:C2H6->CH4CH2            ","bolsig:C2H6(V24)->CH4CH2       ",&
   "bolsig:C2H6(V13)->CH4CH2       ","bolsig:C2H6->CH3CH3            ","bolsig:C2H6(V24)->CH3CH3       ",&
   "bolsig:C2H6(V13)->CH3CH3       ","bolsig:C2H5->C2H4H             ","bolsig:C2H5->C2H3H2            ",&
   "bolsig:C2H5->C2H3HH            ","bolsig:C2H5->C2H2H2H           ","bolsig:C2H5->CH4CH             ",&
   "bolsig:C2H5->CH3CH2            ","bolsig:C2H4->C2H3H             ","bolsig:C2H4(V1)->C2H3H         ",&
   "bolsig:C2H4(V2)->C2H3H         ","bolsig:C2H4->C2H2H2            ","bolsig:C2H4(V1)->C2H2H2        ",&
   "bolsig:C2H4(V2)->C2H2H2        ","bolsig:C2H4->C2H2HH            ","bolsig:C2H4(V1)->C2H2HH        ",&
   "bolsig:C2H4(V2)->C2H2HH        ","bolsig:C2H4->CH3CH             ","bolsig:C2H4(V1)->CH3CH         ",&
   "bolsig:C2H4(V2)->CH3CH         ","bolsig:C2H4->CH2CH2            ","bolsig:C2H4(V1)->CH2CH2        ",&
   "bolsig:C2H4(V2)->CH2CH2        ","bolsig:C2H3->C2H2H             ","bolsig:C2H3->CH2CH             ",&
   "bolsig:C2H2->CHCH              ","bolsig:C2H2(V5)->CHCH          ","bolsig:C2H2(V2)->CHCH          ",&
   "bolsig:C2H2(V13)->CHCH         ","bolsig:C2H6->C2H6^+            ","bolsig:C2H6(V24)->C2H6^+       ",&
   "bolsig:C2H6(V13)->C2H6^+       ","bolsig:C2H6->HC2H5^+           ","bolsig:C2H6(V24)->HC2H5^+      ",&
   "bolsig:C2H6(V13)->HC2H5^+      ","bolsig:C2H6->H2C2H4^+          ","bolsig:C2H6(V24)->H2C2H4^+     ",&
   "bolsig:C2H6(V13)->H2C2H4^+     ","bolsig:C2H6->H2HC2H3^+         ","bolsig:C2H6(V24)->H2HC2H3^+    ",&
   "bolsig:C2H6(V13)->H2HC2H3^+    ","bolsig:C2H6->H2H2C2H2^+        ","bolsig:C2H6(V24)->H2H2C2H2^+   "/
  data reaction_sign(109:162) &
  /"bolsig:C2H6(V13)->H2H2C2H2^+   ","bolsig:C2H6->CH3CH3^+          ","bolsig:C2H6(V24)->CH3CH3^+     ",&
   "bolsig:C2H6(V13)->CH3CH3^+     ","bolsig:C2H6->CH4CH2^+          ","bolsig:C2H6(V24)->CH4CH2^+     ",&
   "bolsig:C2H6(V13)->CH4CH2^+     ","bolsig:C2H5->C2H5^+            ","bolsig:C2H5->HC2H4^+           ",&
   "bolsig:C2H5->H2C2H3^+          ","bolsig:C2H5->H2HC2H2^+         ","bolsig:C2H5->CH2CH3^+          ",&
   "bolsig:C2H5->CH3CH2^+          ","bolsig:C2H5->CH4CH^+           ","bolsig:C2H4->C2H4^+            ",&
   "bolsig:C2H4(V1)->C2H4^+        ","bolsig:C2H4(V2)->C2H4^+        ","bolsig:C2H4->HC2H3^+           ",&
   "bolsig:C2H4(V1)->HC2H3^+       ","bolsig:C2H4(V2)->HC2H3^+       ","bolsig:C2H4->CHCH3^+           ",&
   "bolsig:C2H4(V1)->CHCH3^+       ","bolsig:C2H4(V2)->CHCH3^+       ","bolsig:C2H4->CH2CH2^+          ",&
   "bolsig:C2H4(V1)->CH2CH2^+      ","bolsig:C2H4(V2)->CH2CH2^+      ","bolsig:C2H4->CH3CH^+           ",&
   "bolsig:C2H4(V1)->CH3CH^+       ","bolsig:C2H4(V2)->CH3CH^+       ","bolsig:C2H3->C2H3^+            ",&
   "bolsig:C2H3->HC2H2^+           ","bolsig:C2H3->CHCH2^+           ","bolsig:C2H3->CH2CH^+           ",&
   "bolsig:C2H3->C2H2H^+           ","bolsig:C2H2->C2H2^+            ","bolsig:C2H2(V5)->C2H2^+        ",&
   "bolsig:C2H2(V2)->C2H2^+        ","bolsig:C2H2(V13)->C2H2^+       ","bolsig:C2H2->CHCH^+            ",&
   "bolsig:C2H2(V5)->CHCH^+        ","bolsig:C2H2(V2)->CHCH^+        ","bolsig:C2H2(V13)->CHCH^+       ",&
   "bolsig:C3H8->C3H7H             ","bolsig:C3H8(V1)->C3H7H         ","bolsig:C3H8(V2)->C3H7H         ",&
   "bolsig:C3H8->C3H6H2            ","bolsig:C3H8(V1)->C3H6H2        ","bolsig:C3H8(V2)->C3H6H2        ",&
   "bolsig:C3H7->C3H6H             ","bolsig:C3H8->C3H8^+            ","bolsig:C3H8(V1)->C3H8^+        ",&
   "bolsig:C3H8(V2)->C3H8^+        ","bolsig:C3H8->HC3H7^+           ","bolsig:C3H8(V1)->HC3H7^+       "/
  data reaction_sign(163:216) &
  /"bolsig:C3H8(V2)->HC3H7^+       ","bolsig:C3H8->H2C3H6^+          ","bolsig:C3H8(V1)->H2C3H6^+      ",&
   "bolsig:C3H8(V2)->H2C3H6^+      ","bolsig:C3H7->C3H7^+            ","bolsig:C3H7->HC3H6^+           ",&
   "bolsig:C3H6->C3H6^+            ","bolsig:C3H6(V)->C3H6^+         ","bolsig:C3H6->C2H2CH4           ",&
   "bolsig:C3H6(V)->C2H2CH4        ","bolsig:C3H6->C2H2CH4^+         ","bolsig:C3H6(V)->C2H2CH4^+      ",&
   "E+E+CH^+=>E+CH                 ","E+E+CH2^+=>E+CH2               ","E+E+CH3^+=>E+CH3               ",&
   "E+E+CH4^+=>E+CH4               ","E+E+C2H2^+=>E+C2H2             ","E+E+C2H3^+=>E+C2H3             ",&
   "E+E+C2H4^+=>E+C2H4             ","E+E+C2H5^+=>E+C2H5             ","E+E+C2H6^+=>E+C2H6             ",&
   "E+E+C3H6^+=>E+C3H6             ","E+E+C3H7^+=>E+C3H7             ","E+E+C3H8^+=>E+C3H8             ",&
   "E+CH4^+=>CH3+H                 ","E+CH4^+=>CH2+H+H               ","E+CH4^+=>CH+H2+H               ",&
   "E+CH3^+=>CH2+H                 ","E+CH3^+=>CH+H2                 ","E+CH3^+=>CH+H+H                ",&
   "E+CH3^+=>C+H2+H                ","E+CH2^+=>CH+H                  ","E+CH2^+=>C+H2                  ",&
   "E+CH2^+=>C+H+H                 ","E+CH^+=>C+H                    ","E+C2H6^+=>C2H5+H               ",&
   "E+C2H6^+=>C2H4+H+H             ","E+C2H5^+=>C2H4+H               ","E+C2H5^+=>C2H3+H+H             ",&
   "E+C2H5^+=>C2H2+H2+H            ","E+C2H5^+=>C2H2+H+H+H           ","E+C2H5^+=>CH3+CH2              ",&
   "E+C2H4^+=>C2H3+H               ","E+C2H4^+=>C2H2+H+H             ","E+C2H3^+=>C2H2+H               ",&
   "E+C2H2^+=>CH+CH                ","E+H2^+=>H+H                    ","E+H3^+=>H+H2                   ",&
   "H2+H2^+=>H+H3^+                ","C2H6+CH4^+=>CH4+H2+C2H4^+      ","C2H6(V24)+CH4^+=>CH4+H2+C2H4^+ ",&
   "C2H6(V13)+CH4^+=>CH4+H2+C2H4^+ ","C2H4+CH4^+=>CH3+C2H5^+         ","C2H4(V1)+CH4^+=>CH3+C2H5^+     "/
  data reaction_sign(217:270) &
  /"C2H4(V2)+CH4^+=>CH3+C2H5^+     ","C2H4+CH4^+=>CH4+C2H4^+         ","C2H4(V1)+CH4^+=>CH4+C2H4^+     ",&
   "C2H4(V2)+CH4^+=>CH4+C2H4^+     ","C2H2+CH4^+=>CH3+C2H3^+         ","C2H2(V5)+CH4^+=>CH3+C2H3^+     ",&
   "C2H2(V2)+CH4^+=>CH3+C2H3^+     ","C2H2(V13)+CH4^+=>CH3+C2H3^+    ","C2H2+CH4^+=>CH4+C2H2^+         ",&
   "C2H2(V5)+CH4^+=>CH4+C2H2^+     ","C2H2(V2)+CH4^+=>CH4+C2H2^+     ","C2H2(V13)+CH4^+=>CH4+C2H2^+    ",&
   "H+CH4^+=>H2+CH3^+              ","CH4+CH3^+=>CH3+CH4^+           ","CH4(V24)+CH3^+=>CH3+CH4^+      ",&
   "CH4(V13)+CH3^+=>CH3+CH4^+      ","CH4+CH3^+=>H2+C2H5^+           ","CH4(V24)+CH3^+=>H2+C2H5^+      ",&
   "CH4(V13)+CH3^+=>H2+C2H5^+      ","CH2+CH3^+=>H2+C2H3^+           ","CH+CH3^+=>H2+C2H2^+            ",&
   "C2H6+CH3^+=>CH4+C2H5^+         ","C2H6(V24)+CH3^+=>CH4+C2H5^+    ","C2H6(V13)+CH3^+=>CH4+C2H5^+    ",&
   "C2H4+CH3^+=>CH4+C2H3^+         ","C2H4(V1)+CH3^+=>CH4+C2H3^+     ","C2H4(V2)+CH3^+=>CH4+C2H3^+     ",&
   "C2H3+CH3^+=>CH3+C2H3^+         ","CH4+CH2^+=>CH3+CH3^+           ","CH4(V24)+CH2^+=>CH3+CH3^+      ",&
   "CH4(V13)+CH2^+=>CH3+CH3^+      ","CH4+CH2^+=>H+C2H5^+            ","CH4(V24)+CH2^+=>H+C2H5^+       ",&
   "CH4(V13)+CH2^+=>H+C2H5^+       ","CH4+CH2^+=>H2+C2H4^+           ","CH4(V24)+CH2^+=>H2+C2H4^+      ",&
   "CH4(V13)+CH2^+=>H2+C2H4^+      ","CH4+CH2^+=>H+H2+C2H3^+         ","CH4(V24)+CH2^+=>H+H2+C2H3^+    ",&
   "CH4(V13)+CH2^+=>H+H2+C2H3^+    ","CH4+CH2^+=>H2+H2+C2H2^+        ","CH4(V24)+CH2^+=>H2+H2+C2H2^+   ",&
   "CH4(V13)+CH2^+=>H2+H2+C2H2^+   ","H2+CH2^+=>H+CH3^+              ","CH4+CH^+=>H+C2H4^+             ",&
   "CH4(V24)+CH^+=>H+C2H4^+        ","CH4(V13)+CH^+=>H+C2H4^+        ","CH4+CH^+=>H2+C2H3^+            ",&
   "CH4(V24)+CH^+=>H2+C2H3^+       ","CH4(V13)+CH^+=>H2+C2H3^+       ","CH4+CH^+=>H2+H+C2H2^+          ",&
   "CH4(V24)+CH^+=>H2+H+C2H2^+     ","CH4(V13)+CH^+=>H2+H+C2H2^+     ","H2+CH^+=>H+CH2^+               "/
  data reaction_sign(271:324) &
  /"C2H4+C2H6^+=>C2H6+C2H4^+       ","C2H4(V1)+C2H6^+=>C2H6+C2H4^+   ","C2H4(V2)+C2H6^+=>C2H6+C2H4^+   ",&
   "C2H2+C2H6^+=>C2H3+C2H5^+       ","C2H2(V5)+C2H6^+=>C2H3+C2H5^+   ","C2H2(V2)+C2H6^+=>C2H3+C2H5^+   ",&
   "C2H2(V13)+C2H6^+=>C2H3+C2H5^+  ","H+C2H6^+=>H2+C2H5^+            ","H+C2H5^+=>H2+C2H4^+            ",&
   "C2H3+C2H4^+=>C2H2+C2H5^+       ","C2H3+C2H4^+=>C2H4+C2H3^+       ","H+C2H4^+=>H2+C2H3^+            ",&
   "C2H6+C2H3^+=>C2H4+C2H5^+       ","C2H6(V24)+C2H3^+=>C2H4+C2H5^+  ","C2H6(V13)+C2H3^+=>C2H4+C2H5^+  ",&
   "C2H4+C2H3^+=>C2H2+C2H5^+       ","C2H4(V1)+C2H3^+=>C2H2+C2H5^+   ","C2H4(V2)+C2H3^+=>C2H2+C2H5^+   ",&
   "H+C2H3^+=>H2+C2H2^+            ","CH4+C2H2^+=>CH3+C2H3^+         ","CH4(V24)+C2H2^+=>CH3+C2H3^+    ",&
   "CH4(V13)+C2H2^+=>CH3+C2H3^+    ","C2H6+C2H2^+=>C2H3+C2H5^+       ","C2H6(V24)+C2H2^+=>C2H3+C2H5^+  ",&
   "C2H6(V13)+C2H2^+=>C2H3+C2H5^+  ","C2H6+C2H2^+=>C2H4+C2H4^+       ","C2H6(V24)+C2H2^+=>C2H4+C2H4^+  ",&
   "C2H6(V13)+C2H2^+=>C2H4+C2H4^+  ","C2H4+C2H2^+=>C2H2+C2H4^+       ","C2H4(V1)+C2H2^+=>C2H2+C2H4^+   ",&
   "C2H4(V2)+C2H2^+=>C2H2+C2H4^+   ","C2H3+C2H2^+=>C2H2+C2H3^+       ","H2+C2H2^+=>H+C2H3^+            ",&
   "CH3+H3^+=>H2+CH4^+             ","CH2+H3^+=>H2+CH3^+             ","CH+H3^+=>H2+CH2^+              ",&
   "C2H6+H3^+=>H2+H2+C2H5^+        ","C2H6(V24)+H3^+=>H2+H2+C2H5^+   ","C2H6(V13)+H3^+=>H2+H2+C2H5^+   ",&
   "C2H5+H3^+=>H2+C2H6^+           ","C2H4+H3^+=>H2+C2H5^+           ","C2H4(V1)+H3^+=>H2+C2H5^+       ",&
   "C2H4(V2)+H3^+=>H2+C2H5^+       ","C2H4+H3^+=>H2+H2+C2H3^+        ","C2H4(V1)+H3^+=>H2+H2+C2H3^+    ",&
   "C2H4(V2)+H3^+=>H2+H2+C2H3^+    ","C2H3+H3^+=>H2+C2H4^+           ","C2H2+H3^+=>H2+C2H3^+           ",&
   "C2H2(V5)+H3^+=>H2+C2H3^+       ","C2H2(V2)+H3^+=>H2+C2H3^+       ","C2H2(V13)+H3^+=>H2+C2H3^+      ",&
   "CH4+H2^+=>H2+CH4^+             ","CH4(V24)+H2^+=>H2+CH4^+        ","CH4(V13)+H2^+=>H2+CH4^+        "/
  data reaction_sign(325:378) &
  /"CH4+H2^+=>H2+H+CH3^+           ","CH4(V24)+H2^+=>H2+H+CH3^+      ","CH4(V13)+H2^+=>H2+H+CH3^+      ",&
   "CH2+H2^+=>H+CH3^+              ","CH2+H2^+=>H2+CH2^+             ","CH+H2^+=>H+CH2^+               ",&
   "CH+H2^+=>H2+CH^+               ","C2H6+H2^+=>H2+C2H6^+           ","C2H6(V24)+H2^+=>H2+C2H6^+      ",&
   "C2H6(V13)+H2^+=>H2+C2H6^+      ","C2H6+H2^+=>H2+H+C2H5^+         ","C2H6(V24)+H2^+=>H2+H+C2H5^+    ",&
   "C2H6(V13)+H2^+=>H2+H+C2H5^+    ","C2H6+H2^+=>H2+H2+C2H4^+        ","C2H6(V24)+H2^+=>H2+H2+C2H4^+   ",&
   "C2H6(V13)+H2^+=>H2+H2+C2H4^+   ","C2H6+H2^+=>H2+H2+H+C2H3^+      ","C2H6(V24)+H2^+=>H2+H2+H+C2H3^+ ",&
   "C2H6(V13)+H2^+=>H2+H2+H+C2H3^+ ","C2H6+H2^+=>H2+H2+H2+C2H2^+     ","C2H6(V24)+H2^+=>H2+H2+H2+C2H2^+",&
   "C2H6(V13)+H2^+=>H2+H2+H2+C2H2^+","C2H4+H2^+=>H2+C2H4^+           ","C2H4(V1)+H2^+=>H2+C2H4^+       ",&
   "C2H4(V2)+H2^+=>H2+C2H4^+       ","C2H4+H2^+=>H2+H+C2H3^+         ","C2H4(V1)+H2^+=>H2+H+C2H3^+     ",&
   "C2H4(V2)+H2^+=>H2+H+C2H3^+     ","C2H4+H2^+=>H2+H2+C2H2^+        ","C2H4(V1)+H2^+=>H2+H2+C2H2^+    ",&
   "C2H4(V2)+H2^+=>H2+H2+C2H2^+    ","C2H2+H2^+=>H+C2H3^+            ","C2H2(V5)+H2^+=>H+C2H3^+        ",&
   "C2H2(V2)+H2^+=>H+C2H3^+        ","C2H2(V13)+H2^+=>H+C2H3^+       ","C2H2+H2^+=>H2+C2H2^+           ",&
   "C2H2(V5)+H2^+=>H2+C2H2^+       ","C2H2(V2)+H2^+=>H2+C2H2^+       ","C2H2(V13)+H2^+=>H2+C2H2^+      ",&
   "H+H2^+=>H3^+                   ","H+H2^+=>H2+H^+                 ","CH4+H^+=>H+CH4^+               ",&
   "CH4(V24)+H^+=>H+CH4^+          ","CH4(V13)+H^+=>H+CH4^+          ","CH4+H^+=>H2+CH3^+              ",&
   "CH4(V24)+H^+=>H2+CH3^+         ","CH4(V13)+H^+=>H2+CH3^+         ","CH3+H^+=>H+CH3^+               ",&
   "CH2+H^+=>H+CH2^+               ","CH2+H^+=>H2+CH^+               ","CH+H^+=>H+CH^+                 ",&
   "C2H6+H^+=>H2+C2H5^+            ","C2H6(V24)+H^+=>H2+C2H5^+       ","C2H6(V13)+H^+=>H2+C2H5^+       "/
  data reaction_sign(379:432) &
  /"C2H6+H^+=>H2+H+C2H4^+          ","C2H6(V24)+H^+=>H2+H+C2H4^+     ","C2H6(V13)+H^+=>H2+H+C2H4^+     ",&
   "C2H6+H^+=>H2+H2+C2H3^+         ","C2H6(V24)+H^+=>H2+H2+C2H3^+    ","C2H6(V13)+H^+=>H2+H2+C2H3^+    ",&
   "C2H5+H^+=>H2+C2H4^+            ","C2H5+H^+=>H2+H+C2H3^+          ","C2H4+H^+=>H+C2H4^+             ",&
   "C2H4(V1)+H^+=>H+C2H4^+         ","C2H4(V2)+H^+=>H+C2H4^+         ","C2H4+H^+=>H2+C2H3^+            ",&
   "C2H4(V1)+H^+=>H2+C2H3^+        ","C2H4(V2)+H^+=>H2+C2H3^+        ","C2H4+H^+=>H2+H+C2H2^+          ",&
   "C2H3+H^+=>H+C2H3^+             ","C2H3+H^+=>H2+C2H2^+            ","C2H2+H^+=>H+C2H2^+             ",&
   "C2H2(V5)+H^+=>H+C2H2^+         ","C2H2(V2)+H^+=>H+C2H2^+         ","C2H2(V13)+H^+=>H+C2H2^+        ",&
   "CH4+CH2=>CH3+CH3               ","CH4(V24)+CH2=>CH3+CH3          ","CH4(V13)+CH2=>CH3+CH3          ",&
   "CH4+CH=>C2H4+H                 ","CH4(V24)+CH=>C2H4+H            ","CH4(V13)+CH=>C2H4+H            ",&
   "CH4+C2H5=>C2H6+CH3             ","CH4(V24)+C2H5=>C2H6+CH3        ","CH4(V13)+C2H5=>C2H6+CH3        ",&
   "CH4+C2H3=>C2H4+CH3             ","CH4(V24)+C2H3=>C2H4+CH3        ","CH4(V13)+C2H3=>C2H4+CH3        ",&
   "CH4+C3H7=>C3H8+CH3             ","CH4(V24)+C3H7=>C3H8+CH3        ","CH4(V13)+C3H7=>C3H8+CH3        ",&
   "CH4+H=>CH3+H2                  ","CH4(V24)+H=>CH3+H2             ","CH4(V13)+H=>CH3+H2             ",&
   "CH4+CH3=>H+C2H6                ","CH4(V24)+CH3=>H+C2H6           ","CH4(V13)+CH3=>H+C2H6           ",&
   "CH4+CH2=>C2H6                  ","CH4(V24)+CH2=>C2H6             ","CH4(V13)+CH2=>C2H6             ",&
   "CH4=>CH3+H                     ","CH4(V24)=>CH3+H                ","CH4(V13)=>CH3+H                ",&
   "CH3=>CH2+H                     ","CH3=>CH+H2                     ","CH3+C2H5=>C2H6+CH2             ",&
   "CH2+CH2=>C2H2+H2               ","CH2+C2H5=>C2H4+CH3             ","CH2+C2H3=>C2H2+CH3             "/
  data reaction_sign(433:486) &
  /"CH2+C3H8=>C3H7+CH3             ","CH2+C3H8(V1)=>C3H7+CH3         ","CH2+C3H8(V2)=>C3H7+CH3         ",&
   "CH2+C3H7=>C2H4+C2H5            ","CH2+C3H7=>C3H6+CH3             ","CH2+H2=>CH3+H                  ",&
   "CH2+H=>CH+H2                   ","CH2=>CH+H                      ","CH2+H=>CH3                     ",&
   "CH+C2H6=>C3H6+H                ","CH+C2H6(V24)=>C3H6+H           ","CH+C2H6(V13)=>C3H6+H           ",&
   "CH+C2H6=>C3H7                  ","CH+C2H6(V24)=>C3H7             ","CH+C2H6(V13)=>C3H7             ",&
   "CH+H2=>CH2+H                   ","CH+CH3=>C2H3+H                 ","CH+CH2=>C2H2+H                 ",&
   "CH+H2=>CH3                     ","CH+C2H3=>CH2+C2H2              ","C2H6+C2H3=>C2H5+C2H4           ",&
   "C2H6(V24)+C2H3=>C2H5+C2H4      ","C2H6(V13)+C2H3=>C2H5+C2H4      ","C2H6+C3H7=>C3H8+C2H5           ",&
   "C2H6(V24)+C3H7=>C3H8+C2H5      ","C2H6(V13)+C3H7=>C3H8+C2H5      ","C2H6+H=>C2H5+H2                ",&
   "C2H6(V24)+H=>C2H5+H2           ","C2H6(V13)+H=>C2H5+H2           ","C2H6+H=>CH4+CH3                ",&
   "C2H6(V24)+H=>CH4+CH3           ","C2H6(V13)+H=>CH4+CH3           ","C2H6=>CH3+CH3                  ",&
   "C2H6(V24)=>CH3+CH3             ","C2H6(V13)=>CH3+CH3             ","C2H6+CH=>C2H4+CH3              ",&
   "C2H6(V24)+CH=>C2H4+CH3         ","C2H6(V13)+CH=>C2H4+CH3         ","C2H6+CH2=>C2H5+CH3             ",&
   "C2H6(V24)+CH2=>C2H5+CH3        ","C2H6(V13)+CH2=>C2H5+CH3        ","C2H5+C2H5=>C2H6+C2H4           ",&
   "C2H5+C2H4=>C2H6+C2H3           ","C2H5+C2H4(V1)=>C2H6+C2H3       ","C2H5+C2H4(V2)=>C2H6+C2H3       ",&
   "C2H5+C3H8=>C2H6+C3H7           ","C2H5+C3H8(V1)=>C2H6+C3H7       ","C2H5+C3H8(V2)=>C2H6+C3H7       ",&
   "C2H5+C3H7=>C3H8+C2H4           ","C2H5+C3H7=>C3H6+C2H6           ","C2H5+C3H6=>C3H5+C2H6           ",&
   "C2H5+C3H6(V)=>C3H5+C2H6        ","C2H5+H2=>C2H6+H                ","C2H5+H=>CH3+CH3                "/
  data reaction_sign(487:540) &
  /"C2H5+H=>C2H4+H2                ","C2H5+H=>C2H6                   ","C2H5=>C2H4+H                   ",&
   "C2H5+C2H5=>C4H9H               ","C2H5+C2H3=>C2H4+C2H4           ","C2H4+H=>C2H3+H2                ",&
   "C2H4(V1)+H=>C2H3+H2            ","C2H4(V2)+H=>C2H3+H2            ","C2H4+H=>C2H5                   ",&
   "C2H4(V1)+H=>C2H5               ","C2H4(V2)+H=>C2H5               ","C2H4+H2=>C2H5+H                ",&
   "C2H4=>C2H3+H                   ","C2H4(V1)=>C2H3+H               ","C2H4(V2)=>C2H3+H               ",&
   "C2H4+C2H2=>C2H3+C2H3           ","C2H4(V1)+C2H2=>C2H3+C2H3       ","C2H4(V2)+C2H2=>C2H3+C2H3       ",&
   "C2H4+C2H2(V5)=>C2H3+C2H3       ","C2H4+C2H2(V2)=>C2H3+C2H3       ","C2H4+C2H2(V13)=>C2H3+C2H3      ",&
   "C2H4+C3H6=>C2H3+C3H7           ","C2H4(V1)+C3H6=>C2H3+C3H7       ","C2H4(V2)+C3H6=>C2H3+C3H7       ",&
   "C2H4+C3H6(V)=>C2H3+C3H7        ","C2H4+C2H4=>C2H5+C2H3           ","C2H4(V1)+C2H4(V1)=>C2H5+C2H3   ",&
   "C2H4(V2)+C2H4(V2)=>C2H5+C2H3   ","C2H4+CH3=>C3H7                 ","C2H4(V1)+CH3=>C3H7             ",&
   "C2H4(V2)+CH3=>C3H7             ","C2H4=>C2H2+H2                  ","C2H4(V1)=>C2H2+H2              ",&
   "C2H4(V2)=>C2H2+H2              ","C2H4+H2=>C2H6                  ","C2H4(V1)+H2=>C2H6              ",&
   "C2H4(V2)+H2=>C2H6              ","C2H4+CH2=>C3H6                 ","C2H4(V1)+CH2=>C3H6             ",&
   "C2H4(V2)+CH2=>C3H6             ","C2H3+C2H3=>C2H4+C2H2           ","C2H3+C3H8=>C2H4+C3H7           ",&
   "C2H3+C3H8(V1)=>C2H4+C3H7       ","C2H3+C3H8(V2)=>C2H4+C3H7       ","C2H3+C3H7=>C3H8+C2H2           ",&
   "C2H3+C3H7=>C3H6+C2H4           ","C2H3+H2=>C2H4+H                ","C2H3+H=>C2H2+H2                ",&
   "C2H3+H=>C2H4                   ","C2H3=>C2H2+H                   ","C2H2+H=>C2H3                   ",&
   "C2H2(V5)+H=>C2H3               ","C2H2(V2)+H=>C2H3               ","C2H2(V13)+H=>C2H3              "/
  data reaction_sign(541:577) &
  /"C2H2+H2=>C2H4                  ","C2H2(V5)+H2=>C2H4              ","C2H2(V2)+H2=>C2H4              ",&
   "C2H2(V13)+H2=>C2H4             ","C2H2+H2=>C2H3+H                ","C2H2(V5)+H2=>C2H3+H            ",&
   "C2H2(V2)+H2=>C2H3+H            ","C2H2(V13)+H2=>C2H3+H           ","C2H2+CH3=>C3H5                 ",&
   "C2H2(V5)+CH3=>C3H5             ","C2H2(V2)+CH3=>C3H5             ","C2H2(V13)+CH3=>C3H5            ",&
   "C3H8+H=>C3H7+H2                ","C3H8(V1)+H=>C3H7+H2            ","C3H8(V2)+H=>C3H7+H2            ",&
   "C3H8=>C2H5+CH3                 ","C3H8(V1)=>C2H5+CH3             ","C3H8(V2)=>C2H5+CH3             ",&
   "C3H8+CH2=>C4H9H                ","C3H8(V1)+CH2=>C4H9H            ","C3H8(V2)+CH2=>C4H9H            ",&
   "C3H7+C3H7=>C3H6+C3H8           ","C3H7+H2=>C3H8+H                ","C3H7+H=>C3H6+H2                ",&
   "C3H7+H=>C3H8                   ","C3H7+H=>CH3+C2H5               ","C3H7=>C3H6+H                   ",&
   "C3H7=>C2H4+CH3                 ","C3H6+H=>C3H7                   ","C3H6(V)+H=>C3H7                ",&
   "C3H6=>CH3+C2H3                 ","C4H9H+CH3=>CH4+C4H9            ","C4H9H=>C3H7+CH3                ",&
   "C4H9H=>C2H5+C2H5               ","C4H9H+CH2=>C5H12               ","H2=>H+H                        ",&
   "H+H=>H2                        "/
  data bolsig_species(1:bolsig_species_max) &
  /"CH       ","CH2      ","CH3      ","CH4      ","CH4(V13) ","CH4(V24) ","C2H2     ","C2H2(V13)","C2H2(V2) ","C2H2(V5) ",&
   "C2H3     ","C2H4     ","C2H4(V1) ","C2H4(V2) ","C2H5     ","C2H6     ","C2H6(V13)","C2H6(V24)","C3H6     ","C3H6(V)  ",&
   "C3H7     ","C3H8     ","C3H8(V1) ","C3H8(V2) "/
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
211 format(i3,1x,A31)
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
311 format(361x,47(1x,i9))
312 format(A3,1x,A31,1x,47(1x,A9))
313 format(i3,1x,A31,1x,47(1x,1pd9.2))
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
      write(ifile_unit,"(577(i3))",err=200) int(mrtm(i,:))
    enddo
    lerror = .false.
200 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file " // &
                                     "<qt_matrix.txt> (subroutine writer_save_qtplaskin)")
    close(ifile_unit)
    open(ifile_unit,file="qt_densities.txt",action="write",err=300)
    write(ifile_unit,"(1x,A14,47(111x,i2.2))",err=300) "Time_s", ( i, i = 1, species_max )
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
    write(ifile_unit,"(1x,A12,577(101x,i3.3))",err=500) "Time_s", ( i, i = 1, reactions_max )
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
      write(ifile_unit,"(1pe15.6,47(1pe13.4))") densav(0,2), densav(1:,2)
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
      write(ifile_unit,"(578(1pe13.4))") densav(0,2), rrt_loc(:)
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
  reac_source_local(09,001) = - reac_rate_local(001) 
  reac_source_local(11,001) = + reac_rate_local(001) 
  reac_source_local(09,002) = - reac_rate_local(002) 
  reac_source_local(10,002) = + reac_rate_local(002) 
  reac_source_local(26,003) = - reac_rate_local(003) 
  reac_source_local(28,003) = + reac_rate_local(003) 
  reac_source_local(26,004) = - reac_rate_local(004) 
  reac_source_local(27,004) = + reac_rate_local(004) 
  reac_source_local(20,005) = - reac_rate_local(005) 
  reac_source_local(21,005) = + reac_rate_local(005) 
  reac_source_local(20,006) = - reac_rate_local(006) 
  reac_source_local(22,006) = + reac_rate_local(006) 
  reac_source_local(13,007) = - reac_rate_local(007) 
  reac_source_local(16,007) = + reac_rate_local(007) 
  reac_source_local(13,008) = - reac_rate_local(008) 
  reac_source_local(15,008) = + reac_rate_local(008) 
  reac_source_local(13,009) = - reac_rate_local(009) 
  reac_source_local(14,009) = + reac_rate_local(009) 
  reac_source_local(36,010) = - reac_rate_local(010) 
  reac_source_local(37,010) = + reac_rate_local(010) 
  reac_source_local(36,011) = - reac_rate_local(011) 
  reac_source_local(38,011) = + reac_rate_local(011) 
  reac_source_local(31,012) = - reac_rate_local(012) 
  reac_source_local(32,012) = + reac_rate_local(012) 
  reac_source_local(07,013) = + reac_rate_local(013) 
  reac_source_local(09,013) = - reac_rate_local(013) 
  reac_source_local(43,013) = + reac_rate_local(013) 
  reac_source_local(07,014) = + reac_rate_local(014) 
  reac_source_local(11,014) = - reac_rate_local(014) 
  reac_source_local(43,014) = + reac_rate_local(014) 
  reac_source_local(07,015) = + reac_rate_local(015) 
  reac_source_local(10,015) = - reac_rate_local(015) 
  reac_source_local(43,015) = + reac_rate_local(015) 
  reac_source_local(05,016) = + reac_rate_local(016) 
  reac_source_local(09,016) = - reac_rate_local(016) 
  reac_source_local(44,016) = + reac_rate_local(016) 
  reac_source_local(05,017) = + reac_rate_local(017) 
  reac_source_local(11,017) = - reac_rate_local(017) 
  reac_source_local(44,017) = + reac_rate_local(017) 
  reac_source_local(05,018) = + reac_rate_local(018) 
  reac_source_local(10,018) = - reac_rate_local(018) 
  reac_source_local(44,018) = + reac_rate_local(018) 
  reac_source_local(03,019) = + reac_rate_local(019) 
  reac_source_local(09,019) = - reac_rate_local(019) 
  reac_source_local(43,019) = + reac_rate_local(019) 
  reac_source_local(44,019) = + reac_rate_local(019) 
  reac_source_local(03,020) = + reac_rate_local(020) 
  reac_source_local(11,020) = - reac_rate_local(020) 
  reac_source_local(43,020) = + reac_rate_local(020) 
  reac_source_local(44,020) = + reac_rate_local(020) 
  reac_source_local(03,021) = + reac_rate_local(021) 
  reac_source_local(10,021) = - reac_rate_local(021) 
  reac_source_local(43,021) = + reac_rate_local(021) 
  reac_source_local(44,021) = + reac_rate_local(021) 
  reac_source_local(02,022) = + reac_rate_local(022) 
  reac_source_local(09,022) = - reac_rate_local(022) 
  reac_source_local(44,022) = + reac_rate_local(022) * 2.d0
  reac_source_local(02,023) = + reac_rate_local(023) 
  reac_source_local(11,023) = - reac_rate_local(023) 
  reac_source_local(44,023) = + reac_rate_local(023) * 2.d0
  reac_source_local(02,024) = + reac_rate_local(024) 
  reac_source_local(10,024) = - reac_rate_local(024) 
  reac_source_local(44,024) = + reac_rate_local(024) * 2.d0
  reac_source_local(05,025) = + reac_rate_local(025) 
  reac_source_local(07,025) = - reac_rate_local(025) 
  reac_source_local(43,025) = + reac_rate_local(025) 
  reac_source_local(03,026) = + reac_rate_local(026) 
  reac_source_local(07,026) = - reac_rate_local(026) 
  reac_source_local(44,026) = + reac_rate_local(026) 
  reac_source_local(02,027) = + reac_rate_local(027) 
  reac_source_local(07,027) = - reac_rate_local(027) 
  reac_source_local(43,027) = + reac_rate_local(027) 
  reac_source_local(44,027) = + reac_rate_local(027) 
  reac_source_local(03,028) = + reac_rate_local(028) 
  reac_source_local(05,028) = - reac_rate_local(028) 
  reac_source_local(43,028) = + reac_rate_local(028) 
  reac_source_local(02,029) = + reac_rate_local(029) 
  reac_source_local(05,029) = - reac_rate_local(029) 
  reac_source_local(44,029) = + reac_rate_local(029) 
  reac_source_local(02,030) = + reac_rate_local(030) 
  reac_source_local(05,030) = - reac_rate_local(030) 
  reac_source_local(43,030) = + reac_rate_local(030) * 2.d0
  reac_source_local(02,031) = + reac_rate_local(031) 
  reac_source_local(03,031) = - reac_rate_local(031) 
  reac_source_local(43,031) = + reac_rate_local(031) 
  reac_source_local(01,032) = + reac_rate_local(032) 
  reac_source_local(09,032) = - reac_rate_local(032) 
  reac_source_local(12,032) = + reac_rate_local(032) 
  reac_source_local(01,033) = + reac_rate_local(033) 
  reac_source_local(11,033) = - reac_rate_local(033) 
  reac_source_local(12,033) = + reac_rate_local(033) 
  reac_source_local(01,034) = + reac_rate_local(034) 
  reac_source_local(10,034) = - reac_rate_local(034) 
  reac_source_local(12,034) = + reac_rate_local(034) 
  reac_source_local(01,035) = + reac_rate_local(035) 
  reac_source_local(08,035) = + reac_rate_local(035) 
  reac_source_local(09,035) = - reac_rate_local(035) 
  reac_source_local(43,035) = + reac_rate_local(035) 
  reac_source_local(01,036) = + reac_rate_local(036) 
  reac_source_local(08,036) = + reac_rate_local(036) 
  reac_source_local(11,036) = - reac_rate_local(036) 
  reac_source_local(43,036) = + reac_rate_local(036) 
  reac_source_local(01,037) = + reac_rate_local(037) 
  reac_source_local(08,037) = + reac_rate_local(037) 
  reac_source_local(10,037) = - reac_rate_local(037) 
  reac_source_local(43,037) = + reac_rate_local(037) 
  reac_source_local(01,038) = + reac_rate_local(038) 
  reac_source_local(06,038) = + reac_rate_local(038) 
  reac_source_local(09,038) = - reac_rate_local(038) 
  reac_source_local(44,038) = + reac_rate_local(038) 
  reac_source_local(01,039) = + reac_rate_local(039) 
  reac_source_local(06,039) = + reac_rate_local(039) 
  reac_source_local(11,039) = - reac_rate_local(039) 
  reac_source_local(44,039) = + reac_rate_local(039) 
  reac_source_local(01,040) = + reac_rate_local(040) 
  reac_source_local(06,040) = + reac_rate_local(040) 
  reac_source_local(10,040) = - reac_rate_local(040) 
  reac_source_local(44,040) = + reac_rate_local(040) 
  reac_source_local(01,041) = + reac_rate_local(041) 
  reac_source_local(04,041) = + reac_rate_local(041) 
  reac_source_local(09,041) = - reac_rate_local(041) 
  reac_source_local(43,041) = + reac_rate_local(041) 
  reac_source_local(44,041) = + reac_rate_local(041) 
  reac_source_local(01,042) = + reac_rate_local(042) 
  reac_source_local(04,042) = + reac_rate_local(042) 
  reac_source_local(11,042) = - reac_rate_local(042) 
  reac_source_local(43,042) = + reac_rate_local(042) 
  reac_source_local(44,042) = + reac_rate_local(042) 
  reac_source_local(01,043) = + reac_rate_local(043) 
  reac_source_local(04,043) = + reac_rate_local(043) 
  reac_source_local(10,043) = - reac_rate_local(043) 
  reac_source_local(43,043) = + reac_rate_local(043) 
  reac_source_local(44,043) = + reac_rate_local(043) 
  reac_source_local(01,044) = + reac_rate_local(044) 
  reac_source_local(07,044) = - reac_rate_local(044) 
  reac_source_local(08,044) = + reac_rate_local(044) 
  reac_source_local(01,045) = + reac_rate_local(045) 
  reac_source_local(06,045) = + reac_rate_local(045) 
  reac_source_local(07,045) = - reac_rate_local(045) 
  reac_source_local(43,045) = + reac_rate_local(045) 
  reac_source_local(01,046) = + reac_rate_local(046) 
  reac_source_local(04,046) = + reac_rate_local(046) 
  reac_source_local(07,046) = - reac_rate_local(046) 
  reac_source_local(44,046) = + reac_rate_local(046) 
  reac_source_local(01,047) = + reac_rate_local(047) 
  reac_source_local(05,047) = - reac_rate_local(047) 
  reac_source_local(06,047) = + reac_rate_local(047) 
  reac_source_local(01,048) = + reac_rate_local(048) 
  reac_source_local(04,048) = + reac_rate_local(048) 
  reac_source_local(05,048) = - reac_rate_local(048) 
  reac_source_local(43,048) = + reac_rate_local(048) 
  reac_source_local(01,049) = + reac_rate_local(049) 
  reac_source_local(03,049) = - reac_rate_local(049) 
  reac_source_local(04,049) = + reac_rate_local(049) 
  reac_source_local(24,050) = + reac_rate_local(050) 
  reac_source_local(26,050) = - reac_rate_local(050) 
  reac_source_local(43,050) = + reac_rate_local(050) 
  reac_source_local(24,051) = + reac_rate_local(051) 
  reac_source_local(28,051) = - reac_rate_local(051) 
  reac_source_local(43,051) = + reac_rate_local(051) 
  reac_source_local(24,052) = + reac_rate_local(052) 
  reac_source_local(27,052) = - reac_rate_local(052) 
  reac_source_local(43,052) = + reac_rate_local(052) 
  reac_source_local(20,053) = + reac_rate_local(053) 
  reac_source_local(26,053) = - reac_rate_local(053) 
  reac_source_local(44,053) = + reac_rate_local(053) 
  reac_source_local(20,054) = + reac_rate_local(054) 
  reac_source_local(28,054) = - reac_rate_local(054) 
  reac_source_local(44,054) = + reac_rate_local(054) 
  reac_source_local(20,055) = + reac_rate_local(055) 
  reac_source_local(27,055) = - reac_rate_local(055) 
  reac_source_local(44,055) = + reac_rate_local(055) 
  reac_source_local(18,056) = + reac_rate_local(056) 
  reac_source_local(26,056) = - reac_rate_local(056) 
  reac_source_local(43,056) = + reac_rate_local(056) 
  reac_source_local(44,056) = + reac_rate_local(056) 
  reac_source_local(18,057) = + reac_rate_local(057) 
  reac_source_local(28,057) = - reac_rate_local(057) 
  reac_source_local(43,057) = + reac_rate_local(057) 
  reac_source_local(44,057) = + reac_rate_local(057) 
  reac_source_local(18,058) = + reac_rate_local(058) 
  reac_source_local(27,058) = - reac_rate_local(058) 
  reac_source_local(43,058) = + reac_rate_local(058) 
  reac_source_local(44,058) = + reac_rate_local(058) 
  reac_source_local(13,059) = + reac_rate_local(059) 
  reac_source_local(26,059) = - reac_rate_local(059) 
  reac_source_local(44,059) = + reac_rate_local(059) * 2.d0
  reac_source_local(13,060) = + reac_rate_local(060) 
  reac_source_local(28,060) = - reac_rate_local(060) 
  reac_source_local(44,060) = + reac_rate_local(060) * 2.d0
  reac_source_local(13,061) = + reac_rate_local(061) 
  reac_source_local(27,061) = - reac_rate_local(061) 
  reac_source_local(44,061) = + reac_rate_local(061) * 2.d0
  reac_source_local(05,062) = + reac_rate_local(062) 
  reac_source_local(09,062) = + reac_rate_local(062) 
  reac_source_local(26,062) = - reac_rate_local(062) 
  reac_source_local(05,063) = + reac_rate_local(063) 
  reac_source_local(09,063) = + reac_rate_local(063) 
  reac_source_local(28,063) = - reac_rate_local(063) 
  reac_source_local(05,064) = + reac_rate_local(064) 
  reac_source_local(09,064) = + reac_rate_local(064) 
  reac_source_local(27,064) = - reac_rate_local(064) 
  reac_source_local(07,065) = + reac_rate_local(065) * 2.d0
  reac_source_local(26,065) = - reac_rate_local(065) 
  reac_source_local(07,066) = + reac_rate_local(066) * 2.d0
  reac_source_local(28,066) = - reac_rate_local(066) 
  reac_source_local(07,067) = + reac_rate_local(067) * 2.d0
  reac_source_local(27,067) = - reac_rate_local(067) 
  reac_source_local(20,068) = + reac_rate_local(068) 
  reac_source_local(24,068) = - reac_rate_local(068) 
  reac_source_local(43,068) = + reac_rate_local(068) 
  reac_source_local(18,069) = + reac_rate_local(069) 
  reac_source_local(24,069) = - reac_rate_local(069) 
  reac_source_local(44,069) = + reac_rate_local(069) 
  reac_source_local(18,070) = + reac_rate_local(070) 
  reac_source_local(24,070) = - reac_rate_local(070) 
  reac_source_local(43,070) = + reac_rate_local(070) * 2.d0
  reac_source_local(13,071) = + reac_rate_local(071) 
  reac_source_local(24,071) = - reac_rate_local(071) 
  reac_source_local(43,071) = + reac_rate_local(071) 
  reac_source_local(44,071) = + reac_rate_local(071) 
  reac_source_local(03,072) = + reac_rate_local(072) 
  reac_source_local(09,072) = + reac_rate_local(072) 
  reac_source_local(24,072) = - reac_rate_local(072) 
  reac_source_local(05,073) = + reac_rate_local(073) 
  reac_source_local(07,073) = + reac_rate_local(073) 
  reac_source_local(24,073) = - reac_rate_local(073) 
  reac_source_local(18,074) = + reac_rate_local(074) 
  reac_source_local(20,074) = - reac_rate_local(074) 
  reac_source_local(43,074) = + reac_rate_local(074) 
  reac_source_local(18,075) = + reac_rate_local(075) 
  reac_source_local(21,075) = - reac_rate_local(075) 
  reac_source_local(43,075) = + reac_rate_local(075) 
  reac_source_local(18,076) = + reac_rate_local(076) 
  reac_source_local(22,076) = - reac_rate_local(076) 
  reac_source_local(43,076) = + reac_rate_local(076) 
  reac_source_local(13,077) = + reac_rate_local(077) 
  reac_source_local(20,077) = - reac_rate_local(077) 
  reac_source_local(44,077) = + reac_rate_local(077) 
  reac_source_local(13,078) = + reac_rate_local(078) 
  reac_source_local(21,078) = - reac_rate_local(078) 
  reac_source_local(44,078) = + reac_rate_local(078) 
  reac_source_local(13,079) = + reac_rate_local(079) 
  reac_source_local(22,079) = - reac_rate_local(079) 
  reac_source_local(44,079) = + reac_rate_local(079) 
  reac_source_local(13,080) = + reac_rate_local(080) 
  reac_source_local(20,080) = - reac_rate_local(080) 
  reac_source_local(43,080) = + reac_rate_local(080) * 2.d0
  reac_source_local(13,081) = + reac_rate_local(081) 
  reac_source_local(21,081) = - reac_rate_local(081) 
  reac_source_local(43,081) = + reac_rate_local(081) * 2.d0
  reac_source_local(13,082) = + reac_rate_local(082) 
  reac_source_local(22,082) = - reac_rate_local(082) 
  reac_source_local(43,082) = + reac_rate_local(082) * 2.d0
  reac_source_local(03,083) = + reac_rate_local(083) 
  reac_source_local(07,083) = + reac_rate_local(083) 
  reac_source_local(20,083) = - reac_rate_local(083) 
  reac_source_local(03,084) = + reac_rate_local(084) 
  reac_source_local(07,084) = + reac_rate_local(084) 
  reac_source_local(21,084) = - reac_rate_local(084) 
  reac_source_local(03,085) = + reac_rate_local(085) 
  reac_source_local(07,085) = + reac_rate_local(085) 
  reac_source_local(22,085) = - reac_rate_local(085) 
  reac_source_local(05,086) = + reac_rate_local(086) * 2.d0
  reac_source_local(20,086) = - reac_rate_local(086) 
  reac_source_local(05,087) = + reac_rate_local(087) * 2.d0
  reac_source_local(21,087) = - reac_rate_local(087) 
  reac_source_local(05,088) = + reac_rate_local(088) * 2.d0
  reac_source_local(22,088) = - reac_rate_local(088) 
  reac_source_local(13,089) = + reac_rate_local(089) 
  reac_source_local(18,089) = - reac_rate_local(089) 
  reac_source_local(43,089) = + reac_rate_local(089) 
  reac_source_local(03,090) = + reac_rate_local(090) 
  reac_source_local(05,090) = + reac_rate_local(090) 
  reac_source_local(18,090) = - reac_rate_local(090) 
  reac_source_local(03,091) = + reac_rate_local(091) * 2.d0
  reac_source_local(13,091) = - reac_rate_local(091) 
  reac_source_local(03,092) = + reac_rate_local(092) * 2.d0
  reac_source_local(16,092) = - reac_rate_local(092) 
  reac_source_local(03,093) = + reac_rate_local(093) * 2.d0
  reac_source_local(15,093) = - reac_rate_local(093) 
  reac_source_local(03,094) = + reac_rate_local(094) * 2.d0
  reac_source_local(14,094) = - reac_rate_local(094) 
  reac_source_local(01,095) = + reac_rate_local(095) 
  reac_source_local(26,095) = - reac_rate_local(095) 
  reac_source_local(29,095) = + reac_rate_local(095) 
  reac_source_local(01,096) = + reac_rate_local(096) 
  reac_source_local(28,096) = - reac_rate_local(096) 
  reac_source_local(29,096) = + reac_rate_local(096) 
  reac_source_local(01,097) = + reac_rate_local(097) 
  reac_source_local(27,097) = - reac_rate_local(097) 
  reac_source_local(29,097) = + reac_rate_local(097) 
  reac_source_local(01,098) = + reac_rate_local(098) 
  reac_source_local(25,098) = + reac_rate_local(098) 
  reac_source_local(26,098) = - reac_rate_local(098) 
  reac_source_local(43,098) = + reac_rate_local(098) 
  reac_source_local(01,099) = + reac_rate_local(099) 
  reac_source_local(25,099) = + reac_rate_local(099) 
  reac_source_local(28,099) = - reac_rate_local(099) 
  reac_source_local(43,099) = + reac_rate_local(099) 
  reac_source_local(01,100) = + reac_rate_local(100) 
  reac_source_local(25,100) = + reac_rate_local(100) 
  reac_source_local(27,100) = - reac_rate_local(100) 
  reac_source_local(43,100) = + reac_rate_local(100) 
  reac_source_local(01,101) = + reac_rate_local(101) 
  reac_source_local(23,101) = + reac_rate_local(101) 
  reac_source_local(26,101) = - reac_rate_local(101) 
  reac_source_local(44,101) = + reac_rate_local(101) 
  reac_source_local(01,102) = + reac_rate_local(102) 
  reac_source_local(23,102) = + reac_rate_local(102) 
  reac_source_local(28,102) = - reac_rate_local(102) 
  reac_source_local(44,102) = + reac_rate_local(102) 
  reac_source_local(01,103) = + reac_rate_local(103) 
  reac_source_local(23,103) = + reac_rate_local(103) 
  reac_source_local(27,103) = - reac_rate_local(103) 
  reac_source_local(44,103) = + reac_rate_local(103) 
  reac_source_local(01,104) = + reac_rate_local(104) 
  reac_source_local(19,104) = + reac_rate_local(104) 
  reac_source_local(26,104) = - reac_rate_local(104) 
  reac_source_local(43,104) = + reac_rate_local(104) 
  reac_source_local(44,104) = + reac_rate_local(104) 
  reac_source_local(01,105) = + reac_rate_local(105) 
  reac_source_local(19,105) = + reac_rate_local(105) 
  reac_source_local(28,105) = - reac_rate_local(105) 
  reac_source_local(43,105) = + reac_rate_local(105) 
  reac_source_local(44,105) = + reac_rate_local(105) 
  reac_source_local(01,106) = + reac_rate_local(106) 
  reac_source_local(19,106) = + reac_rate_local(106) 
  reac_source_local(27,106) = - reac_rate_local(106) 
  reac_source_local(43,106) = + reac_rate_local(106) 
  reac_source_local(44,106) = + reac_rate_local(106) 
  reac_source_local(01,107) = + reac_rate_local(107) 
  reac_source_local(17,107) = + reac_rate_local(107) 
  reac_source_local(26,107) = - reac_rate_local(107) 
  reac_source_local(44,107) = + reac_rate_local(107) * 2.d0
  reac_source_local(01,108) = + reac_rate_local(108) 
  reac_source_local(17,108) = + reac_rate_local(108) 
  reac_source_local(28,108) = - reac_rate_local(108) 
  reac_source_local(44,108) = + reac_rate_local(108) * 2.d0
  reac_source_local(01,109) = + reac_rate_local(109) 
  reac_source_local(17,109) = + reac_rate_local(109) 
  reac_source_local(27,109) = - reac_rate_local(109) 
  reac_source_local(44,109) = + reac_rate_local(109) * 2.d0
  reac_source_local(01,110) = + reac_rate_local(110) 
  reac_source_local(07,110) = + reac_rate_local(110) 
  reac_source_local(08,110) = + reac_rate_local(110) 
  reac_source_local(26,110) = - reac_rate_local(110) 
  reac_source_local(01,111) = + reac_rate_local(111) 
  reac_source_local(07,111) = + reac_rate_local(111) 
  reac_source_local(08,111) = + reac_rate_local(111) 
  reac_source_local(28,111) = - reac_rate_local(111) 
  reac_source_local(01,112) = + reac_rate_local(112) 
  reac_source_local(07,112) = + reac_rate_local(112) 
  reac_source_local(08,112) = + reac_rate_local(112) 
  reac_source_local(27,112) = - reac_rate_local(112) 
  reac_source_local(01,113) = + reac_rate_local(113) 
  reac_source_local(06,113) = + reac_rate_local(113) 
  reac_source_local(09,113) = + reac_rate_local(113) 
  reac_source_local(26,113) = - reac_rate_local(113) 
  reac_source_local(01,114) = + reac_rate_local(114) 
  reac_source_local(06,114) = + reac_rate_local(114) 
  reac_source_local(09,114) = + reac_rate_local(114) 
  reac_source_local(28,114) = - reac_rate_local(114) 
  reac_source_local(01,115) = + reac_rate_local(115) 
  reac_source_local(06,115) = + reac_rate_local(115) 
  reac_source_local(09,115) = + reac_rate_local(115) 
  reac_source_local(27,115) = - reac_rate_local(115) 
  reac_source_local(01,116) = + reac_rate_local(116) 
  reac_source_local(24,116) = - reac_rate_local(116) 
  reac_source_local(25,116) = + reac_rate_local(116) 
  reac_source_local(01,117) = + reac_rate_local(117) 
  reac_source_local(23,117) = + reac_rate_local(117) 
  reac_source_local(24,117) = - reac_rate_local(117) 
  reac_source_local(43,117) = + reac_rate_local(117) 
  reac_source_local(01,118) = + reac_rate_local(118) 
  reac_source_local(19,118) = + reac_rate_local(118) 
  reac_source_local(24,118) = - reac_rate_local(118) 
  reac_source_local(44,118) = + reac_rate_local(118) 
  reac_source_local(01,119) = + reac_rate_local(119) 
  reac_source_local(17,119) = + reac_rate_local(119) 
  reac_source_local(24,119) = - reac_rate_local(119) 
  reac_source_local(43,119) = + reac_rate_local(119) 
  reac_source_local(44,119) = + reac_rate_local(119) 
  reac_source_local(01,120) = + reac_rate_local(120) 
  reac_source_local(05,120) = + reac_rate_local(120) 
  reac_source_local(08,120) = + reac_rate_local(120) 
  reac_source_local(24,120) = - reac_rate_local(120) 
  reac_source_local(01,121) = + reac_rate_local(121) 
  reac_source_local(06,121) = + reac_rate_local(121) 
  reac_source_local(07,121) = + reac_rate_local(121) 
  reac_source_local(24,121) = - reac_rate_local(121) 
  reac_source_local(01,122) = + reac_rate_local(122) 
  reac_source_local(04,122) = + reac_rate_local(122) 
  reac_source_local(09,122) = + reac_rate_local(122) 
  reac_source_local(24,122) = - reac_rate_local(122) 
  reac_source_local(01,123) = + reac_rate_local(123) 
  reac_source_local(20,123) = - reac_rate_local(123) 
  reac_source_local(23,123) = + reac_rate_local(123) 
  reac_source_local(01,124) = + reac_rate_local(124) 
  reac_source_local(21,124) = - reac_rate_local(124) 
  reac_source_local(23,124) = + reac_rate_local(124) 
  reac_source_local(01,125) = + reac_rate_local(125) 
  reac_source_local(22,125) = - reac_rate_local(125) 
  reac_source_local(23,125) = + reac_rate_local(125) 
  reac_source_local(01,126) = + reac_rate_local(126) 
  reac_source_local(19,126) = + reac_rate_local(126) 
  reac_source_local(20,126) = - reac_rate_local(126) 
  reac_source_local(43,126) = + reac_rate_local(126) 
  reac_source_local(01,127) = + reac_rate_local(127) 
  reac_source_local(19,127) = + reac_rate_local(127) 
  reac_source_local(21,127) = - reac_rate_local(127) 
  reac_source_local(43,127) = + reac_rate_local(127) 
  reac_source_local(01,128) = + reac_rate_local(128) 
  reac_source_local(19,128) = + reac_rate_local(128) 
  reac_source_local(22,128) = - reac_rate_local(128) 
  reac_source_local(43,128) = + reac_rate_local(128) 
  reac_source_local(01,129) = + reac_rate_local(129) 
  reac_source_local(03,129) = + reac_rate_local(129) 
  reac_source_local(08,129) = + reac_rate_local(129) 
  reac_source_local(20,129) = - reac_rate_local(129) 
  reac_source_local(01,130) = + reac_rate_local(130) 
  reac_source_local(03,130) = + reac_rate_local(130) 
  reac_source_local(08,130) = + reac_rate_local(130) 
  reac_source_local(21,130) = - reac_rate_local(130) 
  reac_source_local(01,131) = + reac_rate_local(131) 
  reac_source_local(03,131) = + reac_rate_local(131) 
  reac_source_local(08,131) = + reac_rate_local(131) 
  reac_source_local(22,131) = - reac_rate_local(131) 
  reac_source_local(01,132) = + reac_rate_local(132) 
  reac_source_local(05,132) = + reac_rate_local(132) 
  reac_source_local(06,132) = + reac_rate_local(132) 
  reac_source_local(20,132) = - reac_rate_local(132) 
  reac_source_local(01,133) = + reac_rate_local(133) 
  reac_source_local(05,133) = + reac_rate_local(133) 
  reac_source_local(06,133) = + reac_rate_local(133) 
  reac_source_local(21,133) = - reac_rate_local(133) 
  reac_source_local(01,134) = + reac_rate_local(134) 
  reac_source_local(05,134) = + reac_rate_local(134) 
  reac_source_local(06,134) = + reac_rate_local(134) 
  reac_source_local(22,134) = - reac_rate_local(134) 
  reac_source_local(01,135) = + reac_rate_local(135) 
  reac_source_local(04,135) = + reac_rate_local(135) 
  reac_source_local(07,135) = + reac_rate_local(135) 
  reac_source_local(20,135) = - reac_rate_local(135) 
  reac_source_local(01,136) = + reac_rate_local(136) 
  reac_source_local(04,136) = + reac_rate_local(136) 
  reac_source_local(07,136) = + reac_rate_local(136) 
  reac_source_local(21,136) = - reac_rate_local(136) 
  reac_source_local(01,137) = + reac_rate_local(137) 
  reac_source_local(04,137) = + reac_rate_local(137) 
  reac_source_local(07,137) = + reac_rate_local(137) 
  reac_source_local(22,137) = - reac_rate_local(137) 
  reac_source_local(01,138) = + reac_rate_local(138) 
  reac_source_local(18,138) = - reac_rate_local(138) 
  reac_source_local(19,138) = + reac_rate_local(138) 
  reac_source_local(01,139) = + reac_rate_local(139) 
  reac_source_local(17,139) = + reac_rate_local(139) 
  reac_source_local(18,139) = - reac_rate_local(139) 
  reac_source_local(43,139) = + reac_rate_local(139) 
  reac_source_local(01,140) = + reac_rate_local(140) 
  reac_source_local(03,140) = + reac_rate_local(140) 
  reac_source_local(06,140) = + reac_rate_local(140) 
  reac_source_local(18,140) = - reac_rate_local(140) 
  reac_source_local(01,141) = + reac_rate_local(141) 
  reac_source_local(04,141) = + reac_rate_local(141) 
  reac_source_local(05,141) = + reac_rate_local(141) 
  reac_source_local(18,141) = - reac_rate_local(141) 
  reac_source_local(01,142) = + reac_rate_local(142) 
  reac_source_local(13,142) = + reac_rate_local(142) 
  reac_source_local(18,142) = - reac_rate_local(142) 
  reac_source_local(47,142) = + reac_rate_local(142) 
  reac_source_local(01,143) = + reac_rate_local(143) 
  reac_source_local(13,143) = - reac_rate_local(143) 
  reac_source_local(17,143) = + reac_rate_local(143) 
  reac_source_local(01,144) = + reac_rate_local(144) 
  reac_source_local(16,144) = - reac_rate_local(144) 
  reac_source_local(17,144) = + reac_rate_local(144) 
  reac_source_local(01,145) = + reac_rate_local(145) 
  reac_source_local(15,145) = - reac_rate_local(145) 
  reac_source_local(17,145) = + reac_rate_local(145) 
  reac_source_local(01,146) = + reac_rate_local(146) 
  reac_source_local(14,146) = - reac_rate_local(146) 
  reac_source_local(17,146) = + reac_rate_local(146) 
  reac_source_local(01,147) = + reac_rate_local(147) 
  reac_source_local(03,147) = + reac_rate_local(147) 
  reac_source_local(04,147) = + reac_rate_local(147) 
  reac_source_local(13,147) = - reac_rate_local(147) 
  reac_source_local(01,148) = + reac_rate_local(148) 
  reac_source_local(03,148) = + reac_rate_local(148) 
  reac_source_local(04,148) = + reac_rate_local(148) 
  reac_source_local(16,148) = - reac_rate_local(148) 
  reac_source_local(01,149) = + reac_rate_local(149) 
  reac_source_local(03,149) = + reac_rate_local(149) 
  reac_source_local(04,149) = + reac_rate_local(149) 
  reac_source_local(15,149) = - reac_rate_local(149) 
  reac_source_local(01,150) = + reac_rate_local(150) 
  reac_source_local(03,150) = + reac_rate_local(150) 
  reac_source_local(04,150) = + reac_rate_local(150) 
  reac_source_local(14,150) = - reac_rate_local(150) 
  reac_source_local(34,151) = + reac_rate_local(151) 
  reac_source_local(36,151) = - reac_rate_local(151) 
  reac_source_local(43,151) = + reac_rate_local(151) 
  reac_source_local(34,152) = + reac_rate_local(152) 
  reac_source_local(37,152) = - reac_rate_local(152) 
  reac_source_local(43,152) = + reac_rate_local(152) 
  reac_source_local(34,153) = + reac_rate_local(153) 
  reac_source_local(38,153) = - reac_rate_local(153) 
  reac_source_local(43,153) = + reac_rate_local(153) 
  reac_source_local(31,154) = + reac_rate_local(154) 
  reac_source_local(36,154) = - reac_rate_local(154) 
  reac_source_local(44,154) = + reac_rate_local(154) 
  reac_source_local(31,155) = + reac_rate_local(155) 
  reac_source_local(37,155) = - reac_rate_local(155) 
  reac_source_local(44,155) = + reac_rate_local(155) 
  reac_source_local(31,156) = + reac_rate_local(156) 
  reac_source_local(38,156) = - reac_rate_local(156) 
  reac_source_local(44,156) = + reac_rate_local(156) 
  reac_source_local(31,157) = + reac_rate_local(157) 
  reac_source_local(34,157) = - reac_rate_local(157) 
  reac_source_local(43,157) = + reac_rate_local(157) 
  reac_source_local(01,158) = + reac_rate_local(158) 
  reac_source_local(36,158) = - reac_rate_local(158) 
  reac_source_local(39,158) = + reac_rate_local(158) 
  reac_source_local(01,159) = + reac_rate_local(159) 
  reac_source_local(37,159) = - reac_rate_local(159) 
  reac_source_local(39,159) = + reac_rate_local(159) 
  reac_source_local(01,160) = + reac_rate_local(160) 
  reac_source_local(38,160) = - reac_rate_local(160) 
  reac_source_local(39,160) = + reac_rate_local(160) 
  reac_source_local(01,161) = + reac_rate_local(161) 
  reac_source_local(35,161) = + reac_rate_local(161) 
  reac_source_local(36,161) = - reac_rate_local(161) 
  reac_source_local(43,161) = + reac_rate_local(161) 
  reac_source_local(01,162) = + reac_rate_local(162) 
  reac_source_local(35,162) = + reac_rate_local(162) 
  reac_source_local(37,162) = - reac_rate_local(162) 
  reac_source_local(43,162) = + reac_rate_local(162) 
  reac_source_local(01,163) = + reac_rate_local(163) 
  reac_source_local(35,163) = + reac_rate_local(163) 
  reac_source_local(38,163) = - reac_rate_local(163) 
  reac_source_local(43,163) = + reac_rate_local(163) 
  reac_source_local(01,164) = + reac_rate_local(164) 
  reac_source_local(33,164) = + reac_rate_local(164) 
  reac_source_local(36,164) = - reac_rate_local(164) 
  reac_source_local(44,164) = + reac_rate_local(164) 
  reac_source_local(01,165) = + reac_rate_local(165) 
  reac_source_local(33,165) = + reac_rate_local(165) 
  reac_source_local(37,165) = - reac_rate_local(165) 
  reac_source_local(44,165) = + reac_rate_local(165) 
  reac_source_local(01,166) = + reac_rate_local(166) 
  reac_source_local(33,166) = + reac_rate_local(166) 
  reac_source_local(38,166) = - reac_rate_local(166) 
  reac_source_local(44,166) = + reac_rate_local(166) 
  reac_source_local(01,167) = + reac_rate_local(167) 
  reac_source_local(34,167) = - reac_rate_local(167) 
  reac_source_local(35,167) = + reac_rate_local(167) 
  reac_source_local(01,168) = + reac_rate_local(168) 
  reac_source_local(33,168) = + reac_rate_local(168) 
  reac_source_local(34,168) = - reac_rate_local(168) 
  reac_source_local(43,168) = + reac_rate_local(168) 
  reac_source_local(01,169) = + reac_rate_local(169) 
  reac_source_local(31,169) = - reac_rate_local(169) 
  reac_source_local(33,169) = + reac_rate_local(169) 
  reac_source_local(01,170) = + reac_rate_local(170) 
  reac_source_local(32,170) = - reac_rate_local(170) 
  reac_source_local(33,170) = + reac_rate_local(170) 
  reac_source_local(09,171) = + reac_rate_local(171) 
  reac_source_local(13,171) = + reac_rate_local(171) 
  reac_source_local(31,171) = - reac_rate_local(171) 
  reac_source_local(09,172) = + reac_rate_local(172) 
  reac_source_local(13,172) = + reac_rate_local(172) 
  reac_source_local(32,172) = - reac_rate_local(172) 
  reac_source_local(01,173) = + reac_rate_local(173) 
  reac_source_local(12,173) = + reac_rate_local(173) 
  reac_source_local(13,173) = + reac_rate_local(173) 
  reac_source_local(31,173) = - reac_rate_local(173) 
  reac_source_local(01,174) = + reac_rate_local(174) 
  reac_source_local(12,174) = + reac_rate_local(174) 
  reac_source_local(13,174) = + reac_rate_local(174) 
  reac_source_local(32,174) = - reac_rate_local(174) 
  reac_source_local(01,175) = - reac_rate_local(175) 
  reac_source_local(03,175) = + reac_rate_local(175) 
  reac_source_local(04,175) = - reac_rate_local(175) 
  reac_source_local(01,176) = - reac_rate_local(176) 
  reac_source_local(05,176) = + reac_rate_local(176) 
  reac_source_local(06,176) = - reac_rate_local(176) 
  reac_source_local(01,177) = - reac_rate_local(177) 
  reac_source_local(07,177) = + reac_rate_local(177) 
  reac_source_local(08,177) = - reac_rate_local(177) 
  reac_source_local(01,178) = - reac_rate_local(178) 
  reac_source_local(09,178) = + reac_rate_local(178) 
  reac_source_local(12,178) = - reac_rate_local(178) 
  reac_source_local(01,179) = - reac_rate_local(179) 
  reac_source_local(13,179) = + reac_rate_local(179) 
  reac_source_local(17,179) = - reac_rate_local(179) 
  reac_source_local(01,180) = - reac_rate_local(180) 
  reac_source_local(18,180) = + reac_rate_local(180) 
  reac_source_local(19,180) = - reac_rate_local(180) 
  reac_source_local(01,181) = - reac_rate_local(181) 
  reac_source_local(20,181) = + reac_rate_local(181) 
  reac_source_local(23,181) = - reac_rate_local(181) 
  reac_source_local(01,182) = - reac_rate_local(182) 
  reac_source_local(24,182) = + reac_rate_local(182) 
  reac_source_local(25,182) = - reac_rate_local(182) 
  reac_source_local(01,183) = - reac_rate_local(183) 
  reac_source_local(26,183) = + reac_rate_local(183) 
  reac_source_local(29,183) = - reac_rate_local(183) 
  reac_source_local(01,184) = - reac_rate_local(184) 
  reac_source_local(31,184) = + reac_rate_local(184) 
  reac_source_local(33,184) = - reac_rate_local(184) 
  reac_source_local(01,185) = - reac_rate_local(185) 
  reac_source_local(34,185) = + reac_rate_local(185) 
  reac_source_local(35,185) = - reac_rate_local(185) 
  reac_source_local(01,186) = - reac_rate_local(186) 
  reac_source_local(36,186) = + reac_rate_local(186) 
  reac_source_local(39,186) = - reac_rate_local(186) 
  reac_source_local(01,187) = - reac_rate_local(187) 
  reac_source_local(07,187) = + reac_rate_local(187) 
  reac_source_local(12,187) = - reac_rate_local(187) 
  reac_source_local(43,187) = + reac_rate_local(187) 
  reac_source_local(01,188) = - reac_rate_local(188) 
  reac_source_local(05,188) = + reac_rate_local(188) 
  reac_source_local(12,188) = - reac_rate_local(188) 
  reac_source_local(43,188) = + reac_rate_local(188) * 2.d0
  reac_source_local(01,189) = - reac_rate_local(189) 
  reac_source_local(03,189) = + reac_rate_local(189) 
  reac_source_local(12,189) = - reac_rate_local(189) 
  reac_source_local(43,189) = + reac_rate_local(189) 
  reac_source_local(44,189) = + reac_rate_local(189) 
  reac_source_local(01,190) = - reac_rate_local(190) 
  reac_source_local(05,190) = + reac_rate_local(190) 
  reac_source_local(08,190) = - reac_rate_local(190) 
  reac_source_local(43,190) = + reac_rate_local(190) 
  reac_source_local(01,191) = - reac_rate_local(191) 
  reac_source_local(03,191) = + reac_rate_local(191) 
  reac_source_local(08,191) = - reac_rate_local(191) 
  reac_source_local(44,191) = + reac_rate_local(191) 
  reac_source_local(01,192) = - reac_rate_local(192) 
  reac_source_local(03,192) = + reac_rate_local(192) 
  reac_source_local(08,192) = - reac_rate_local(192) 
  reac_source_local(43,192) = + reac_rate_local(192) * 2.d0
  reac_source_local(01,193) = - reac_rate_local(193) 
  reac_source_local(02,193) = + reac_rate_local(193) 
  reac_source_local(08,193) = - reac_rate_local(193) 
  reac_source_local(43,193) = + reac_rate_local(193) 
  reac_source_local(44,193) = + reac_rate_local(193) 
  reac_source_local(01,194) = - reac_rate_local(194) 
  reac_source_local(03,194) = + reac_rate_local(194) 
  reac_source_local(06,194) = - reac_rate_local(194) 
  reac_source_local(43,194) = + reac_rate_local(194) 
  reac_source_local(01,195) = - reac_rate_local(195) 
  reac_source_local(02,195) = + reac_rate_local(195) 
  reac_source_local(06,195) = - reac_rate_local(195) 
  reac_source_local(44,195) = + reac_rate_local(195) 
  reac_source_local(01,196) = - reac_rate_local(196) 
  reac_source_local(02,196) = + reac_rate_local(196) 
  reac_source_local(06,196) = - reac_rate_local(196) 
  reac_source_local(43,196) = + reac_rate_local(196) * 2.d0
  reac_source_local(01,197) = - reac_rate_local(197) 
  reac_source_local(02,197) = + reac_rate_local(197) 
  reac_source_local(04,197) = - reac_rate_local(197) 
  reac_source_local(43,197) = + reac_rate_local(197) 
  reac_source_local(01,198) = - reac_rate_local(198) 
  reac_source_local(24,198) = + reac_rate_local(198) 
  reac_source_local(29,198) = - reac_rate_local(198) 
  reac_source_local(43,198) = + reac_rate_local(198) 
  reac_source_local(01,199) = - reac_rate_local(199) 
  reac_source_local(20,199) = + reac_rate_local(199) 
  reac_source_local(29,199) = - reac_rate_local(199) 
  reac_source_local(43,199) = + reac_rate_local(199) * 2.d0
  reac_source_local(01,200) = - reac_rate_local(200) 
  reac_source_local(20,200) = + reac_rate_local(200) 
  reac_source_local(25,200) = - reac_rate_local(200) 
  reac_source_local(43,200) = + reac_rate_local(200) 
  reac_source_local(01,201) = - reac_rate_local(201) 
  reac_source_local(18,201) = + reac_rate_local(201) 
  reac_source_local(25,201) = - reac_rate_local(201) 
  reac_source_local(43,201) = + reac_rate_local(201) * 2.d0
  reac_source_local(01,202) = - reac_rate_local(202) 
  reac_source_local(13,202) = + reac_rate_local(202) 
  reac_source_local(25,202) = - reac_rate_local(202) 
  reac_source_local(43,202) = + reac_rate_local(202) 
  reac_source_local(44,202) = + reac_rate_local(202) 
  reac_source_local(01,203) = - reac_rate_local(203) 
  reac_source_local(13,203) = + reac_rate_local(203) 
  reac_source_local(25,203) = - reac_rate_local(203) 
  reac_source_local(43,203) = + reac_rate_local(203) * 3.d0
  reac_source_local(01,204) = - reac_rate_local(204) 
  reac_source_local(05,204) = + reac_rate_local(204) 
  reac_source_local(07,204) = + reac_rate_local(204) 
  reac_source_local(25,204) = - reac_rate_local(204) 
  reac_source_local(01,205) = - reac_rate_local(205) 
  reac_source_local(18,205) = + reac_rate_local(205) 
  reac_source_local(23,205) = - reac_rate_local(205) 
  reac_source_local(43,205) = + reac_rate_local(205) 
  reac_source_local(01,206) = - reac_rate_local(206) 
  reac_source_local(13,206) = + reac_rate_local(206) 
  reac_source_local(23,206) = - reac_rate_local(206) 
  reac_source_local(43,206) = + reac_rate_local(206) * 2.d0
  reac_source_local(01,207) = - reac_rate_local(207) 
  reac_source_local(13,207) = + reac_rate_local(207) 
  reac_source_local(19,207) = - reac_rate_local(207) 
  reac_source_local(43,207) = + reac_rate_local(207) 
  reac_source_local(01,208) = - reac_rate_local(208) 
  reac_source_local(03,208) = + reac_rate_local(208) * 2.d0
  reac_source_local(17,208) = - reac_rate_local(208) 
  reac_source_local(01,209) = - reac_rate_local(209) 
  reac_source_local(43,209) = + reac_rate_local(209) * 2.d0
  reac_source_local(45,209) = - reac_rate_local(209) 
  reac_source_local(01,210) = - reac_rate_local(210) 
  reac_source_local(43,210) = + reac_rate_local(210) 
  reac_source_local(44,210) = + reac_rate_local(210) 
  reac_source_local(46,210) = - reac_rate_local(210) 
  reac_source_local(43,211) = + reac_rate_local(211) 
  reac_source_local(44,211) = - reac_rate_local(211) 
  reac_source_local(45,211) = - reac_rate_local(211) 
  reac_source_local(46,211) = + reac_rate_local(211) 
  reac_source_local(09,212) = + reac_rate_local(212) 
  reac_source_local(12,212) = - reac_rate_local(212) 
  reac_source_local(23,212) = + reac_rate_local(212) 
  reac_source_local(26,212) = - reac_rate_local(212) 
  reac_source_local(44,212) = + reac_rate_local(212) 
  reac_source_local(09,213) = + reac_rate_local(213) 
  reac_source_local(12,213) = - reac_rate_local(213) 
  reac_source_local(23,213) = + reac_rate_local(213) 
  reac_source_local(28,213) = - reac_rate_local(213) 
  reac_source_local(44,213) = + reac_rate_local(213) 
  reac_source_local(09,214) = + reac_rate_local(214) 
  reac_source_local(12,214) = - reac_rate_local(214) 
  reac_source_local(23,214) = + reac_rate_local(214) 
  reac_source_local(27,214) = - reac_rate_local(214) 
  reac_source_local(44,214) = + reac_rate_local(214) 
  reac_source_local(07,215) = + reac_rate_local(215) 
  reac_source_local(12,215) = - reac_rate_local(215) 
  reac_source_local(20,215) = - reac_rate_local(215) 
  reac_source_local(25,215) = + reac_rate_local(215) 
  reac_source_local(07,216) = + reac_rate_local(216) 
  reac_source_local(12,216) = - reac_rate_local(216) 
  reac_source_local(21,216) = - reac_rate_local(216) 
  reac_source_local(25,216) = + reac_rate_local(216) 
  reac_source_local(07,217) = + reac_rate_local(217) 
  reac_source_local(12,217) = - reac_rate_local(217) 
  reac_source_local(22,217) = - reac_rate_local(217) 
  reac_source_local(25,217) = + reac_rate_local(217) 
  reac_source_local(09,218) = + reac_rate_local(218) 
  reac_source_local(12,218) = - reac_rate_local(218) 
  reac_source_local(20,218) = - reac_rate_local(218) 
  reac_source_local(23,218) = + reac_rate_local(218) 
  reac_source_local(09,219) = + reac_rate_local(219) 
  reac_source_local(12,219) = - reac_rate_local(219) 
  reac_source_local(21,219) = - reac_rate_local(219) 
  reac_source_local(23,219) = + reac_rate_local(219) 
  reac_source_local(09,220) = + reac_rate_local(220) 
  reac_source_local(12,220) = - reac_rate_local(220) 
  reac_source_local(22,220) = - reac_rate_local(220) 
  reac_source_local(23,220) = + reac_rate_local(220) 
  reac_source_local(07,221) = + reac_rate_local(221) 
  reac_source_local(12,221) = - reac_rate_local(221) 
  reac_source_local(13,221) = - reac_rate_local(221) 
  reac_source_local(19,221) = + reac_rate_local(221) 
  reac_source_local(07,222) = + reac_rate_local(222) 
  reac_source_local(12,222) = - reac_rate_local(222) 
  reac_source_local(16,222) = - reac_rate_local(222) 
  reac_source_local(19,222) = + reac_rate_local(222) 
  reac_source_local(07,223) = + reac_rate_local(223) 
  reac_source_local(12,223) = - reac_rate_local(223) 
  reac_source_local(15,223) = - reac_rate_local(223) 
  reac_source_local(19,223) = + reac_rate_local(223) 
  reac_source_local(07,224) = + reac_rate_local(224) 
  reac_source_local(12,224) = - reac_rate_local(224) 
  reac_source_local(14,224) = - reac_rate_local(224) 
  reac_source_local(19,224) = + reac_rate_local(224) 
  reac_source_local(09,225) = + reac_rate_local(225) 
  reac_source_local(12,225) = - reac_rate_local(225) 
  reac_source_local(13,225) = - reac_rate_local(225) 
  reac_source_local(17,225) = + reac_rate_local(225) 
  reac_source_local(09,226) = + reac_rate_local(226) 
  reac_source_local(12,226) = - reac_rate_local(226) 
  reac_source_local(16,226) = - reac_rate_local(226) 
  reac_source_local(17,226) = + reac_rate_local(226) 
  reac_source_local(09,227) = + reac_rate_local(227) 
  reac_source_local(12,227) = - reac_rate_local(227) 
  reac_source_local(15,227) = - reac_rate_local(227) 
  reac_source_local(17,227) = + reac_rate_local(227) 
  reac_source_local(09,228) = + reac_rate_local(228) 
  reac_source_local(12,228) = - reac_rate_local(228) 
  reac_source_local(14,228) = - reac_rate_local(228) 
  reac_source_local(17,228) = + reac_rate_local(228) 
  reac_source_local(08,229) = + reac_rate_local(229) 
  reac_source_local(12,229) = - reac_rate_local(229) 
  reac_source_local(43,229) = - reac_rate_local(229) 
  reac_source_local(44,229) = + reac_rate_local(229) 
  reac_source_local(07,230) = + reac_rate_local(230) 
  reac_source_local(08,230) = - reac_rate_local(230) 
  reac_source_local(09,230) = - reac_rate_local(230) 
  reac_source_local(12,230) = + reac_rate_local(230) 
  reac_source_local(07,231) = + reac_rate_local(231) 
  reac_source_local(08,231) = - reac_rate_local(231) 
  reac_source_local(11,231) = - reac_rate_local(231) 
  reac_source_local(12,231) = + reac_rate_local(231) 
  reac_source_local(07,232) = + reac_rate_local(232) 
  reac_source_local(08,232) = - reac_rate_local(232) 
  reac_source_local(10,232) = - reac_rate_local(232) 
  reac_source_local(12,232) = + reac_rate_local(232) 
  reac_source_local(08,233) = - reac_rate_local(233) 
  reac_source_local(09,233) = - reac_rate_local(233) 
  reac_source_local(25,233) = + reac_rate_local(233) 
  reac_source_local(44,233) = + reac_rate_local(233) 
  reac_source_local(08,234) = - reac_rate_local(234) 
  reac_source_local(11,234) = - reac_rate_local(234) 
  reac_source_local(25,234) = + reac_rate_local(234) 
  reac_source_local(44,234) = + reac_rate_local(234) 
  reac_source_local(08,235) = - reac_rate_local(235) 
  reac_source_local(10,235) = - reac_rate_local(235) 
  reac_source_local(25,235) = + reac_rate_local(235) 
  reac_source_local(44,235) = + reac_rate_local(235) 
  reac_source_local(05,236) = - reac_rate_local(236) 
  reac_source_local(08,236) = - reac_rate_local(236) 
  reac_source_local(19,236) = + reac_rate_local(236) 
  reac_source_local(44,236) = + reac_rate_local(236) 
  reac_source_local(03,237) = - reac_rate_local(237) 
  reac_source_local(08,237) = - reac_rate_local(237) 
  reac_source_local(17,237) = + reac_rate_local(237) 
  reac_source_local(44,237) = + reac_rate_local(237) 
  reac_source_local(08,238) = - reac_rate_local(238) 
  reac_source_local(09,238) = + reac_rate_local(238) 
  reac_source_local(25,238) = + reac_rate_local(238) 
  reac_source_local(26,238) = - reac_rate_local(238) 
  reac_source_local(08,239) = - reac_rate_local(239) 
  reac_source_local(09,239) = + reac_rate_local(239) 
  reac_source_local(25,239) = + reac_rate_local(239) 
  reac_source_local(28,239) = - reac_rate_local(239) 
  reac_source_local(08,240) = - reac_rate_local(240) 
  reac_source_local(09,240) = + reac_rate_local(240) 
  reac_source_local(25,240) = + reac_rate_local(240) 
  reac_source_local(27,240) = - reac_rate_local(240) 
  reac_source_local(08,241) = - reac_rate_local(241) 
  reac_source_local(09,241) = + reac_rate_local(241) 
  reac_source_local(19,241) = + reac_rate_local(241) 
  reac_source_local(20,241) = - reac_rate_local(241) 
  reac_source_local(08,242) = - reac_rate_local(242) 
  reac_source_local(09,242) = + reac_rate_local(242) 
  reac_source_local(19,242) = + reac_rate_local(242) 
  reac_source_local(21,242) = - reac_rate_local(242) 
  reac_source_local(08,243) = - reac_rate_local(243) 
  reac_source_local(09,243) = + reac_rate_local(243) 
  reac_source_local(19,243) = + reac_rate_local(243) 
  reac_source_local(22,243) = - reac_rate_local(243) 
  reac_source_local(07,244) = + reac_rate_local(244) 
  reac_source_local(08,244) = - reac_rate_local(244) 
  reac_source_local(18,244) = - reac_rate_local(244) 
  reac_source_local(19,244) = + reac_rate_local(244) 
  reac_source_local(06,245) = - reac_rate_local(245) 
  reac_source_local(07,245) = + reac_rate_local(245) 
  reac_source_local(08,245) = + reac_rate_local(245) 
  reac_source_local(09,245) = - reac_rate_local(245) 
  reac_source_local(06,246) = - reac_rate_local(246) 
  reac_source_local(07,246) = + reac_rate_local(246) 
  reac_source_local(08,246) = + reac_rate_local(246) 
  reac_source_local(11,246) = - reac_rate_local(246) 
  reac_source_local(06,247) = - reac_rate_local(247) 
  reac_source_local(07,247) = + reac_rate_local(247) 
  reac_source_local(08,247) = + reac_rate_local(247) 
  reac_source_local(10,247) = - reac_rate_local(247) 
  reac_source_local(06,248) = - reac_rate_local(248) 
  reac_source_local(09,248) = - reac_rate_local(248) 
  reac_source_local(25,248) = + reac_rate_local(248) 
  reac_source_local(43,248) = + reac_rate_local(248) 
  reac_source_local(06,249) = - reac_rate_local(249) 
  reac_source_local(11,249) = - reac_rate_local(249) 
  reac_source_local(25,249) = + reac_rate_local(249) 
  reac_source_local(43,249) = + reac_rate_local(249) 
  reac_source_local(06,250) = - reac_rate_local(250) 
  reac_source_local(10,250) = - reac_rate_local(250) 
  reac_source_local(25,250) = + reac_rate_local(250) 
  reac_source_local(43,250) = + reac_rate_local(250) 
  reac_source_local(06,251) = - reac_rate_local(251) 
  reac_source_local(09,251) = - reac_rate_local(251) 
  reac_source_local(23,251) = + reac_rate_local(251) 
  reac_source_local(44,251) = + reac_rate_local(251) 
  reac_source_local(06,252) = - reac_rate_local(252) 
  reac_source_local(11,252) = - reac_rate_local(252) 
  reac_source_local(23,252) = + reac_rate_local(252) 
  reac_source_local(44,252) = + reac_rate_local(252) 
  reac_source_local(06,253) = - reac_rate_local(253) 
  reac_source_local(10,253) = - reac_rate_local(253) 
  reac_source_local(23,253) = + reac_rate_local(253) 
  reac_source_local(44,253) = + reac_rate_local(253) 
  reac_source_local(06,254) = - reac_rate_local(254) 
  reac_source_local(09,254) = - reac_rate_local(254) 
  reac_source_local(19,254) = + reac_rate_local(254) 
  reac_source_local(43,254) = + reac_rate_local(254) 
  reac_source_local(44,254) = + reac_rate_local(254) 
  reac_source_local(06,255) = - reac_rate_local(255) 
  reac_source_local(11,255) = - reac_rate_local(255) 
  reac_source_local(19,255) = + reac_rate_local(255) 
  reac_source_local(43,255) = + reac_rate_local(255) 
  reac_source_local(44,255) = + reac_rate_local(255) 
  reac_source_local(06,256) = - reac_rate_local(256) 
  reac_source_local(10,256) = - reac_rate_local(256) 
  reac_source_local(19,256) = + reac_rate_local(256) 
  reac_source_local(43,256) = + reac_rate_local(256) 
  reac_source_local(44,256) = + reac_rate_local(256) 
  reac_source_local(06,257) = - reac_rate_local(257) 
  reac_source_local(09,257) = - reac_rate_local(257) 
  reac_source_local(17,257) = + reac_rate_local(257) 
  reac_source_local(44,257) = + reac_rate_local(257) * 2.d0
  reac_source_local(06,258) = - reac_rate_local(258) 
  reac_source_local(11,258) = - reac_rate_local(258) 
  reac_source_local(17,258) = + reac_rate_local(258) 
  reac_source_local(44,258) = + reac_rate_local(258) * 2.d0
  reac_source_local(06,259) = - reac_rate_local(259) 
  reac_source_local(10,259) = - reac_rate_local(259) 
  reac_source_local(17,259) = + reac_rate_local(259) 
  reac_source_local(44,259) = + reac_rate_local(259) * 2.d0
  reac_source_local(06,260) = - reac_rate_local(260) 
  reac_source_local(08,260) = + reac_rate_local(260) 
  reac_source_local(43,260) = + reac_rate_local(260) 
  reac_source_local(44,260) = - reac_rate_local(260) 
  reac_source_local(04,261) = - reac_rate_local(261) 
  reac_source_local(09,261) = - reac_rate_local(261) 
  reac_source_local(23,261) = + reac_rate_local(261) 
  reac_source_local(43,261) = + reac_rate_local(261) 
  reac_source_local(04,262) = - reac_rate_local(262) 
  reac_source_local(11,262) = - reac_rate_local(262) 
  reac_source_local(23,262) = + reac_rate_local(262) 
  reac_source_local(43,262) = + reac_rate_local(262) 
  reac_source_local(04,263) = - reac_rate_local(263) 
  reac_source_local(10,263) = - reac_rate_local(263) 
  reac_source_local(23,263) = + reac_rate_local(263) 
  reac_source_local(43,263) = + reac_rate_local(263) 
  reac_source_local(04,264) = - reac_rate_local(264) 
  reac_source_local(09,264) = - reac_rate_local(264) 
  reac_source_local(19,264) = + reac_rate_local(264) 
  reac_source_local(44,264) = + reac_rate_local(264) 
  reac_source_local(04,265) = - reac_rate_local(265) 
  reac_source_local(11,265) = - reac_rate_local(265) 
  reac_source_local(19,265) = + reac_rate_local(265) 
  reac_source_local(44,265) = + reac_rate_local(265) 
  reac_source_local(04,266) = - reac_rate_local(266) 
  reac_source_local(10,266) = - reac_rate_local(266) 
  reac_source_local(19,266) = + reac_rate_local(266) 
  reac_source_local(44,266) = + reac_rate_local(266) 
  reac_source_local(04,267) = - reac_rate_local(267) 
  reac_source_local(09,267) = - reac_rate_local(267) 
  reac_source_local(17,267) = + reac_rate_local(267) 
  reac_source_local(43,267) = + reac_rate_local(267) 
  reac_source_local(44,267) = + reac_rate_local(267) 
  reac_source_local(04,268) = - reac_rate_local(268) 
  reac_source_local(11,268) = - reac_rate_local(268) 
  reac_source_local(17,268) = + reac_rate_local(268) 
  reac_source_local(43,268) = + reac_rate_local(268) 
  reac_source_local(44,268) = + reac_rate_local(268) 
  reac_source_local(04,269) = - reac_rate_local(269) 
  reac_source_local(10,269) = - reac_rate_local(269) 
  reac_source_local(17,269) = + reac_rate_local(269) 
  reac_source_local(43,269) = + reac_rate_local(269) 
  reac_source_local(44,269) = + reac_rate_local(269) 
  reac_source_local(04,270) = - reac_rate_local(270) 
  reac_source_local(06,270) = + reac_rate_local(270) 
  reac_source_local(43,270) = + reac_rate_local(270) 
  reac_source_local(44,270) = - reac_rate_local(270) 
  reac_source_local(20,271) = - reac_rate_local(271) 
  reac_source_local(23,271) = + reac_rate_local(271) 
  reac_source_local(26,271) = + reac_rate_local(271) 
  reac_source_local(29,271) = - reac_rate_local(271) 
  reac_source_local(21,272) = - reac_rate_local(272) 
  reac_source_local(23,272) = + reac_rate_local(272) 
  reac_source_local(26,272) = + reac_rate_local(272) 
  reac_source_local(29,272) = - reac_rate_local(272) 
  reac_source_local(22,273) = - reac_rate_local(273) 
  reac_source_local(23,273) = + reac_rate_local(273) 
  reac_source_local(26,273) = + reac_rate_local(273) 
  reac_source_local(29,273) = - reac_rate_local(273) 
  reac_source_local(13,274) = - reac_rate_local(274) 
  reac_source_local(18,274) = + reac_rate_local(274) 
  reac_source_local(25,274) = + reac_rate_local(274) 
  reac_source_local(29,274) = - reac_rate_local(274) 
  reac_source_local(16,275) = - reac_rate_local(275) 
  reac_source_local(18,275) = + reac_rate_local(275) 
  reac_source_local(25,275) = + reac_rate_local(275) 
  reac_source_local(29,275) = - reac_rate_local(275) 
  reac_source_local(15,276) = - reac_rate_local(276) 
  reac_source_local(18,276) = + reac_rate_local(276) 
  reac_source_local(25,276) = + reac_rate_local(276) 
  reac_source_local(29,276) = - reac_rate_local(276) 
  reac_source_local(14,277) = - reac_rate_local(277) 
  reac_source_local(18,277) = + reac_rate_local(277) 
  reac_source_local(25,277) = + reac_rate_local(277) 
  reac_source_local(29,277) = - reac_rate_local(277) 
  reac_source_local(25,278) = + reac_rate_local(278) 
  reac_source_local(29,278) = - reac_rate_local(278) 
  reac_source_local(43,278) = - reac_rate_local(278) 
  reac_source_local(44,278) = + reac_rate_local(278) 
  reac_source_local(23,279) = + reac_rate_local(279) 
  reac_source_local(25,279) = - reac_rate_local(279) 
  reac_source_local(43,279) = - reac_rate_local(279) 
  reac_source_local(44,279) = + reac_rate_local(279) 
  reac_source_local(13,280) = + reac_rate_local(280) 
  reac_source_local(18,280) = - reac_rate_local(280) 
  reac_source_local(23,280) = - reac_rate_local(280) 
  reac_source_local(25,280) = + reac_rate_local(280) 
  reac_source_local(18,281) = - reac_rate_local(281) 
  reac_source_local(19,281) = + reac_rate_local(281) 
  reac_source_local(20,281) = + reac_rate_local(281) 
  reac_source_local(23,281) = - reac_rate_local(281) 
  reac_source_local(19,282) = + reac_rate_local(282) 
  reac_source_local(23,282) = - reac_rate_local(282) 
  reac_source_local(43,282) = - reac_rate_local(282) 
  reac_source_local(44,282) = + reac_rate_local(282) 
  reac_source_local(19,283) = - reac_rate_local(283) 
  reac_source_local(20,283) = + reac_rate_local(283) 
  reac_source_local(25,283) = + reac_rate_local(283) 
  reac_source_local(26,283) = - reac_rate_local(283) 
  reac_source_local(19,284) = - reac_rate_local(284) 
  reac_source_local(20,284) = + reac_rate_local(284) 
  reac_source_local(25,284) = + reac_rate_local(284) 
  reac_source_local(28,284) = - reac_rate_local(284) 
  reac_source_local(19,285) = - reac_rate_local(285) 
  reac_source_local(20,285) = + reac_rate_local(285) 
  reac_source_local(25,285) = + reac_rate_local(285) 
  reac_source_local(27,285) = - reac_rate_local(285) 
  reac_source_local(13,286) = + reac_rate_local(286) 
  reac_source_local(19,286) = - reac_rate_local(286) 
  reac_source_local(20,286) = - reac_rate_local(286) 
  reac_source_local(25,286) = + reac_rate_local(286) 
  reac_source_local(13,287) = + reac_rate_local(287) 
  reac_source_local(19,287) = - reac_rate_local(287) 
  reac_source_local(21,287) = - reac_rate_local(287) 
  reac_source_local(25,287) = + reac_rate_local(287) 
  reac_source_local(13,288) = + reac_rate_local(288) 
  reac_source_local(19,288) = - reac_rate_local(288) 
  reac_source_local(22,288) = - reac_rate_local(288) 
  reac_source_local(25,288) = + reac_rate_local(288) 
  reac_source_local(17,289) = + reac_rate_local(289) 
  reac_source_local(19,289) = - reac_rate_local(289) 
  reac_source_local(43,289) = - reac_rate_local(289) 
  reac_source_local(44,289) = + reac_rate_local(289) 
  reac_source_local(07,290) = + reac_rate_local(290) 
  reac_source_local(09,290) = - reac_rate_local(290) 
  reac_source_local(17,290) = - reac_rate_local(290) 
  reac_source_local(19,290) = + reac_rate_local(290) 
  reac_source_local(07,291) = + reac_rate_local(291) 
  reac_source_local(11,291) = - reac_rate_local(291) 
  reac_source_local(17,291) = - reac_rate_local(291) 
  reac_source_local(19,291) = + reac_rate_local(291) 
  reac_source_local(07,292) = + reac_rate_local(292) 
  reac_source_local(10,292) = - reac_rate_local(292) 
  reac_source_local(17,292) = - reac_rate_local(292) 
  reac_source_local(19,292) = + reac_rate_local(292) 
  reac_source_local(17,293) = - reac_rate_local(293) 
  reac_source_local(18,293) = + reac_rate_local(293) 
  reac_source_local(25,293) = + reac_rate_local(293) 
  reac_source_local(26,293) = - reac_rate_local(293) 
  reac_source_local(17,294) = - reac_rate_local(294) 
  reac_source_local(18,294) = + reac_rate_local(294) 
  reac_source_local(25,294) = + reac_rate_local(294) 
  reac_source_local(28,294) = - reac_rate_local(294) 
  reac_source_local(17,295) = - reac_rate_local(295) 
  reac_source_local(18,295) = + reac_rate_local(295) 
  reac_source_local(25,295) = + reac_rate_local(295) 
  reac_source_local(27,295) = - reac_rate_local(295) 
  reac_source_local(17,296) = - reac_rate_local(296) 
  reac_source_local(20,296) = + reac_rate_local(296) 
  reac_source_local(23,296) = + reac_rate_local(296) 
  reac_source_local(26,296) = - reac_rate_local(296) 
  reac_source_local(17,297) = - reac_rate_local(297) 
  reac_source_local(20,297) = + reac_rate_local(297) 
  reac_source_local(23,297) = + reac_rate_local(297) 
  reac_source_local(28,297) = - reac_rate_local(297) 
  reac_source_local(17,298) = - reac_rate_local(298) 
  reac_source_local(20,298) = + reac_rate_local(298) 
  reac_source_local(23,298) = + reac_rate_local(298) 
  reac_source_local(27,298) = - reac_rate_local(298) 
  reac_source_local(13,299) = + reac_rate_local(299) 
  reac_source_local(17,299) = - reac_rate_local(299) 
  reac_source_local(20,299) = - reac_rate_local(299) 
  reac_source_local(23,299) = + reac_rate_local(299) 
  reac_source_local(13,300) = + reac_rate_local(300) 
  reac_source_local(17,300) = - reac_rate_local(300) 
  reac_source_local(21,300) = - reac_rate_local(300) 
  reac_source_local(23,300) = + reac_rate_local(300) 
  reac_source_local(13,301) = + reac_rate_local(301) 
  reac_source_local(17,301) = - reac_rate_local(301) 
  reac_source_local(22,301) = - reac_rate_local(301) 
  reac_source_local(23,301) = + reac_rate_local(301) 
  reac_source_local(13,302) = + reac_rate_local(302) 
  reac_source_local(17,302) = - reac_rate_local(302) 
  reac_source_local(18,302) = - reac_rate_local(302) 
  reac_source_local(19,302) = + reac_rate_local(302) 
  reac_source_local(17,303) = - reac_rate_local(303) 
  reac_source_local(19,303) = + reac_rate_local(303) 
  reac_source_local(43,303) = + reac_rate_local(303) 
  reac_source_local(44,303) = - reac_rate_local(303) 
  reac_source_local(07,304) = - reac_rate_local(304) 
  reac_source_local(12,304) = + reac_rate_local(304) 
  reac_source_local(44,304) = + reac_rate_local(304) 
  reac_source_local(46,304) = - reac_rate_local(304) 
  reac_source_local(05,305) = - reac_rate_local(305) 
  reac_source_local(08,305) = + reac_rate_local(305) 
  reac_source_local(44,305) = + reac_rate_local(305) 
  reac_source_local(46,305) = - reac_rate_local(305) 
  reac_source_local(03,306) = - reac_rate_local(306) 
  reac_source_local(06,306) = + reac_rate_local(306) 
  reac_source_local(44,306) = + reac_rate_local(306) 
  reac_source_local(46,306) = - reac_rate_local(306) 
  reac_source_local(25,307) = + reac_rate_local(307) 
  reac_source_local(26,307) = - reac_rate_local(307) 
  reac_source_local(44,307) = + reac_rate_local(307) * 2.d0
  reac_source_local(46,307) = - reac_rate_local(307) 
  reac_source_local(25,308) = + reac_rate_local(308) 
  reac_source_local(28,308) = - reac_rate_local(308) 
  reac_source_local(44,308) = + reac_rate_local(308) * 2.d0
  reac_source_local(46,308) = - reac_rate_local(308) 
  reac_source_local(25,309) = + reac_rate_local(309) 
  reac_source_local(27,309) = - reac_rate_local(309) 
  reac_source_local(44,309) = + reac_rate_local(309) * 2.d0
  reac_source_local(46,309) = - reac_rate_local(309) 
  reac_source_local(24,310) = - reac_rate_local(310) 
  reac_source_local(29,310) = + reac_rate_local(310) 
  reac_source_local(44,310) = + reac_rate_local(310) 
  reac_source_local(46,310) = - reac_rate_local(310) 
  reac_source_local(20,311) = - reac_rate_local(311) 
  reac_source_local(25,311) = + reac_rate_local(311) 
  reac_source_local(44,311) = + reac_rate_local(311) 
  reac_source_local(46,311) = - reac_rate_local(311) 
  reac_source_local(21,312) = - reac_rate_local(312) 
  reac_source_local(25,312) = + reac_rate_local(312) 
  reac_source_local(44,312) = + reac_rate_local(312) 
  reac_source_local(46,312) = - reac_rate_local(312) 
  reac_source_local(22,313) = - reac_rate_local(313) 
  reac_source_local(25,313) = + reac_rate_local(313) 
  reac_source_local(44,313) = + reac_rate_local(313) 
  reac_source_local(46,313) = - reac_rate_local(313) 
  reac_source_local(19,314) = + reac_rate_local(314) 
  reac_source_local(20,314) = - reac_rate_local(314) 
  reac_source_local(44,314) = + reac_rate_local(314) * 2.d0
  reac_source_local(46,314) = - reac_rate_local(314) 
  reac_source_local(19,315) = + reac_rate_local(315) 
  reac_source_local(21,315) = - reac_rate_local(315) 
  reac_source_local(44,315) = + reac_rate_local(315) * 2.d0
  reac_source_local(46,315) = - reac_rate_local(315) 
  reac_source_local(19,316) = + reac_rate_local(316) 
  reac_source_local(22,316) = - reac_rate_local(316) 
  reac_source_local(44,316) = + reac_rate_local(316) * 2.d0
  reac_source_local(46,316) = - reac_rate_local(316) 
  reac_source_local(18,317) = - reac_rate_local(317) 
  reac_source_local(23,317) = + reac_rate_local(317) 
  reac_source_local(44,317) = + reac_rate_local(317) 
  reac_source_local(46,317) = - reac_rate_local(317) 
  reac_source_local(13,318) = - reac_rate_local(318) 
  reac_source_local(19,318) = + reac_rate_local(318) 
  reac_source_local(44,318) = + reac_rate_local(318) 
  reac_source_local(46,318) = - reac_rate_local(318) 
  reac_source_local(16,319) = - reac_rate_local(319) 
  reac_source_local(19,319) = + reac_rate_local(319) 
  reac_source_local(44,319) = + reac_rate_local(319) 
  reac_source_local(46,319) = - reac_rate_local(319) 
  reac_source_local(15,320) = - reac_rate_local(320) 
  reac_source_local(19,320) = + reac_rate_local(320) 
  reac_source_local(44,320) = + reac_rate_local(320) 
  reac_source_local(46,320) = - reac_rate_local(320) 
  reac_source_local(14,321) = - reac_rate_local(321) 
  reac_source_local(19,321) = + reac_rate_local(321) 
  reac_source_local(44,321) = + reac_rate_local(321) 
  reac_source_local(46,321) = - reac_rate_local(321) 
  reac_source_local(09,322) = - reac_rate_local(322) 
  reac_source_local(12,322) = + reac_rate_local(322) 
  reac_source_local(44,322) = + reac_rate_local(322) 
  reac_source_local(45,322) = - reac_rate_local(322) 
  reac_source_local(11,323) = - reac_rate_local(323) 
  reac_source_local(12,323) = + reac_rate_local(323) 
  reac_source_local(44,323) = + reac_rate_local(323) 
  reac_source_local(45,323) = - reac_rate_local(323) 
  reac_source_local(10,324) = - reac_rate_local(324) 
  reac_source_local(12,324) = + reac_rate_local(324) 
  reac_source_local(44,324) = + reac_rate_local(324) 
  reac_source_local(45,324) = - reac_rate_local(324) 
  reac_source_local(08,325) = + reac_rate_local(325) 
  reac_source_local(09,325) = - reac_rate_local(325) 
  reac_source_local(43,325) = + reac_rate_local(325) 
  reac_source_local(44,325) = + reac_rate_local(325) 
  reac_source_local(45,325) = - reac_rate_local(325) 
  reac_source_local(08,326) = + reac_rate_local(326) 
  reac_source_local(11,326) = - reac_rate_local(326) 
  reac_source_local(43,326) = + reac_rate_local(326) 
  reac_source_local(44,326) = + reac_rate_local(326) 
  reac_source_local(45,326) = - reac_rate_local(326) 
  reac_source_local(08,327) = + reac_rate_local(327) 
  reac_source_local(10,327) = - reac_rate_local(327) 
  reac_source_local(43,327) = + reac_rate_local(327) 
  reac_source_local(44,327) = + reac_rate_local(327) 
  reac_source_local(45,327) = - reac_rate_local(327) 
  reac_source_local(05,328) = - reac_rate_local(328) 
  reac_source_local(08,328) = + reac_rate_local(328) 
  reac_source_local(43,328) = + reac_rate_local(328) 
  reac_source_local(45,328) = - reac_rate_local(328) 
  reac_source_local(05,329) = - reac_rate_local(329) 
  reac_source_local(06,329) = + reac_rate_local(329) 
  reac_source_local(44,329) = + reac_rate_local(329) 
  reac_source_local(45,329) = - reac_rate_local(329) 
  reac_source_local(03,330) = - reac_rate_local(330) 
  reac_source_local(06,330) = + reac_rate_local(330) 
  reac_source_local(43,330) = + reac_rate_local(330) 
  reac_source_local(45,330) = - reac_rate_local(330) 
  reac_source_local(03,331) = - reac_rate_local(331) 
  reac_source_local(04,331) = + reac_rate_local(331) 
  reac_source_local(44,331) = + reac_rate_local(331) 
  reac_source_local(45,331) = - reac_rate_local(331) 
  reac_source_local(26,332) = - reac_rate_local(332) 
  reac_source_local(29,332) = + reac_rate_local(332) 
  reac_source_local(44,332) = + reac_rate_local(332) 
  reac_source_local(45,332) = - reac_rate_local(332) 
  reac_source_local(28,333) = - reac_rate_local(333) 
  reac_source_local(29,333) = + reac_rate_local(333) 
  reac_source_local(44,333) = + reac_rate_local(333) 
  reac_source_local(45,333) = - reac_rate_local(333) 
  reac_source_local(27,334) = - reac_rate_local(334) 
  reac_source_local(29,334) = + reac_rate_local(334) 
  reac_source_local(44,334) = + reac_rate_local(334) 
  reac_source_local(45,334) = - reac_rate_local(334) 
  reac_source_local(25,335) = + reac_rate_local(335) 
  reac_source_local(26,335) = - reac_rate_local(335) 
  reac_source_local(43,335) = + reac_rate_local(335) 
  reac_source_local(44,335) = + reac_rate_local(335) 
  reac_source_local(45,335) = - reac_rate_local(335) 
  reac_source_local(25,336) = + reac_rate_local(336) 
  reac_source_local(28,336) = - reac_rate_local(336) 
  reac_source_local(43,336) = + reac_rate_local(336) 
  reac_source_local(44,336) = + reac_rate_local(336) 
  reac_source_local(45,336) = - reac_rate_local(336) 
  reac_source_local(25,337) = + reac_rate_local(337) 
  reac_source_local(27,337) = - reac_rate_local(337) 
  reac_source_local(43,337) = + reac_rate_local(337) 
  reac_source_local(44,337) = + reac_rate_local(337) 
  reac_source_local(45,337) = - reac_rate_local(337) 
  reac_source_local(23,338) = + reac_rate_local(338) 
  reac_source_local(26,338) = - reac_rate_local(338) 
  reac_source_local(44,338) = + reac_rate_local(338) * 2.d0
  reac_source_local(45,338) = - reac_rate_local(338) 
  reac_source_local(23,339) = + reac_rate_local(339) 
  reac_source_local(28,339) = - reac_rate_local(339) 
  reac_source_local(44,339) = + reac_rate_local(339) * 2.d0
  reac_source_local(45,339) = - reac_rate_local(339) 
  reac_source_local(23,340) = + reac_rate_local(340) 
  reac_source_local(27,340) = - reac_rate_local(340) 
  reac_source_local(44,340) = + reac_rate_local(340) * 2.d0
  reac_source_local(45,340) = - reac_rate_local(340) 
  reac_source_local(19,341) = + reac_rate_local(341) 
  reac_source_local(26,341) = - reac_rate_local(341) 
  reac_source_local(43,341) = + reac_rate_local(341) 
  reac_source_local(44,341) = + reac_rate_local(341) * 2.d0
  reac_source_local(45,341) = - reac_rate_local(341) 
  reac_source_local(19,342) = + reac_rate_local(342) 
  reac_source_local(28,342) = - reac_rate_local(342) 
  reac_source_local(43,342) = + reac_rate_local(342) 
  reac_source_local(44,342) = + reac_rate_local(342) * 2.d0
  reac_source_local(45,342) = - reac_rate_local(342) 
  reac_source_local(19,343) = + reac_rate_local(343) 
  reac_source_local(27,343) = - reac_rate_local(343) 
  reac_source_local(43,343) = + reac_rate_local(343) 
  reac_source_local(44,343) = + reac_rate_local(343) * 2.d0
  reac_source_local(45,343) = - reac_rate_local(343) 
  reac_source_local(17,344) = + reac_rate_local(344) 
  reac_source_local(26,344) = - reac_rate_local(344) 
  reac_source_local(44,344) = + reac_rate_local(344) * 3.d0
  reac_source_local(45,344) = - reac_rate_local(344) 
  reac_source_local(17,345) = + reac_rate_local(345) 
  reac_source_local(28,345) = - reac_rate_local(345) 
  reac_source_local(44,345) = + reac_rate_local(345) * 3.d0
  reac_source_local(45,345) = - reac_rate_local(345) 
  reac_source_local(17,346) = + reac_rate_local(346) 
  reac_source_local(27,346) = - reac_rate_local(346) 
  reac_source_local(44,346) = + reac_rate_local(346) * 3.d0
  reac_source_local(45,346) = - reac_rate_local(346) 
  reac_source_local(20,347) = - reac_rate_local(347) 
  reac_source_local(23,347) = + reac_rate_local(347) 
  reac_source_local(44,347) = + reac_rate_local(347) 
  reac_source_local(45,347) = - reac_rate_local(347) 
  reac_source_local(21,348) = - reac_rate_local(348) 
  reac_source_local(23,348) = + reac_rate_local(348) 
  reac_source_local(44,348) = + reac_rate_local(348) 
  reac_source_local(45,348) = - reac_rate_local(348) 
  reac_source_local(22,349) = - reac_rate_local(349) 
  reac_source_local(23,349) = + reac_rate_local(349) 
  reac_source_local(44,349) = + reac_rate_local(349) 
  reac_source_local(45,349) = - reac_rate_local(349) 
  reac_source_local(19,350) = + reac_rate_local(350) 
  reac_source_local(20,350) = - reac_rate_local(350) 
  reac_source_local(43,350) = + reac_rate_local(350) 
  reac_source_local(44,350) = + reac_rate_local(350) 
  reac_source_local(45,350) = - reac_rate_local(350) 
  reac_source_local(19,351) = + reac_rate_local(351) 
  reac_source_local(21,351) = - reac_rate_local(351) 
  reac_source_local(43,351) = + reac_rate_local(351) 
  reac_source_local(44,351) = + reac_rate_local(351) 
  reac_source_local(45,351) = - reac_rate_local(351) 
  reac_source_local(19,352) = + reac_rate_local(352) 
  reac_source_local(22,352) = - reac_rate_local(352) 
  reac_source_local(43,352) = + reac_rate_local(352) 
  reac_source_local(44,352) = + reac_rate_local(352) 
  reac_source_local(45,352) = - reac_rate_local(352) 
  reac_source_local(17,353) = + reac_rate_local(353) 
  reac_source_local(20,353) = - reac_rate_local(353) 
  reac_source_local(44,353) = + reac_rate_local(353) * 2.d0
  reac_source_local(45,353) = - reac_rate_local(353) 
  reac_source_local(17,354) = + reac_rate_local(354) 
  reac_source_local(21,354) = - reac_rate_local(354) 
  reac_source_local(44,354) = + reac_rate_local(354) * 2.d0
  reac_source_local(45,354) = - reac_rate_local(354) 
  reac_source_local(17,355) = + reac_rate_local(355) 
  reac_source_local(22,355) = - reac_rate_local(355) 
  reac_source_local(44,355) = + reac_rate_local(355) * 2.d0
  reac_source_local(45,355) = - reac_rate_local(355) 
  reac_source_local(13,356) = - reac_rate_local(356) 
  reac_source_local(19,356) = + reac_rate_local(356) 
  reac_source_local(43,356) = + reac_rate_local(356) 
  reac_source_local(45,356) = - reac_rate_local(356) 
  reac_source_local(16,357) = - reac_rate_local(357) 
  reac_source_local(19,357) = + reac_rate_local(357) 
  reac_source_local(43,357) = + reac_rate_local(357) 
  reac_source_local(45,357) = - reac_rate_local(357) 
  reac_source_local(15,358) = - reac_rate_local(358) 
  reac_source_local(19,358) = + reac_rate_local(358) 
  reac_source_local(43,358) = + reac_rate_local(358) 
  reac_source_local(45,358) = - reac_rate_local(358) 
  reac_source_local(14,359) = - reac_rate_local(359) 
  reac_source_local(19,359) = + reac_rate_local(359) 
  reac_source_local(43,359) = + reac_rate_local(359) 
  reac_source_local(45,359) = - reac_rate_local(359) 
  reac_source_local(13,360) = - reac_rate_local(360) 
  reac_source_local(17,360) = + reac_rate_local(360) 
  reac_source_local(44,360) = + reac_rate_local(360) 
  reac_source_local(45,360) = - reac_rate_local(360) 
  reac_source_local(16,361) = - reac_rate_local(361) 
  reac_source_local(17,361) = + reac_rate_local(361) 
  reac_source_local(44,361) = + reac_rate_local(361) 
  reac_source_local(45,361) = - reac_rate_local(361) 
  reac_source_local(15,362) = - reac_rate_local(362) 
  reac_source_local(17,362) = + reac_rate_local(362) 
  reac_source_local(44,362) = + reac_rate_local(362) 
  reac_source_local(45,362) = - reac_rate_local(362) 
  reac_source_local(14,363) = - reac_rate_local(363) 
  reac_source_local(17,363) = + reac_rate_local(363) 
  reac_source_local(44,363) = + reac_rate_local(363) 
  reac_source_local(45,363) = - reac_rate_local(363) 
  reac_source_local(43,364) = - reac_rate_local(364) 
  reac_source_local(45,364) = - reac_rate_local(364) 
  reac_source_local(46,364) = + reac_rate_local(364) 
  reac_source_local(43,365) = - reac_rate_local(365) 
  reac_source_local(44,365) = + reac_rate_local(365) 
  reac_source_local(45,365) = - reac_rate_local(365) 
  reac_source_local(47,365) = + reac_rate_local(365) 
  reac_source_local(09,366) = - reac_rate_local(366) 
  reac_source_local(12,366) = + reac_rate_local(366) 
  reac_source_local(43,366) = + reac_rate_local(366) 
  reac_source_local(47,366) = - reac_rate_local(366) 
  reac_source_local(11,367) = - reac_rate_local(367) 
  reac_source_local(12,367) = + reac_rate_local(367) 
  reac_source_local(43,367) = + reac_rate_local(367) 
  reac_source_local(47,367) = - reac_rate_local(367) 
  reac_source_local(10,368) = - reac_rate_local(368) 
  reac_source_local(12,368) = + reac_rate_local(368) 
  reac_source_local(43,368) = + reac_rate_local(368) 
  reac_source_local(47,368) = - reac_rate_local(368) 
  reac_source_local(08,369) = + reac_rate_local(369) 
  reac_source_local(09,369) = - reac_rate_local(369) 
  reac_source_local(44,369) = + reac_rate_local(369) 
  reac_source_local(47,369) = - reac_rate_local(369) 
  reac_source_local(08,370) = + reac_rate_local(370) 
  reac_source_local(11,370) = - reac_rate_local(370) 
  reac_source_local(44,370) = + reac_rate_local(370) 
  reac_source_local(47,370) = - reac_rate_local(370) 
  reac_source_local(08,371) = + reac_rate_local(371) 
  reac_source_local(10,371) = - reac_rate_local(371) 
  reac_source_local(44,371) = + reac_rate_local(371) 
  reac_source_local(47,371) = - reac_rate_local(371) 
  reac_source_local(07,372) = - reac_rate_local(372) 
  reac_source_local(08,372) = + reac_rate_local(372) 
  reac_source_local(43,372) = + reac_rate_local(372) 
  reac_source_local(47,372) = - reac_rate_local(372) 
  reac_source_local(05,373) = - reac_rate_local(373) 
  reac_source_local(06,373) = + reac_rate_local(373) 
  reac_source_local(43,373) = + reac_rate_local(373) 
  reac_source_local(47,373) = - reac_rate_local(373) 
  reac_source_local(04,374) = + reac_rate_local(374) 
  reac_source_local(05,374) = - reac_rate_local(374) 
  reac_source_local(44,374) = + reac_rate_local(374) 
  reac_source_local(47,374) = - reac_rate_local(374) 
  reac_source_local(03,375) = - reac_rate_local(375) 
  reac_source_local(04,375) = + reac_rate_local(375) 
  reac_source_local(43,375) = + reac_rate_local(375) 
  reac_source_local(47,375) = - reac_rate_local(375) 
  reac_source_local(25,376) = + reac_rate_local(376) 
  reac_source_local(26,376) = - reac_rate_local(376) 
  reac_source_local(44,376) = + reac_rate_local(376) 
  reac_source_local(47,376) = - reac_rate_local(376) 
  reac_source_local(25,377) = + reac_rate_local(377) 
  reac_source_local(28,377) = - reac_rate_local(377) 
  reac_source_local(44,377) = + reac_rate_local(377) 
  reac_source_local(47,377) = - reac_rate_local(377) 
  reac_source_local(25,378) = + reac_rate_local(378) 
  reac_source_local(27,378) = - reac_rate_local(378) 
  reac_source_local(44,378) = + reac_rate_local(378) 
  reac_source_local(47,378) = - reac_rate_local(378) 
  reac_source_local(23,379) = + reac_rate_local(379) 
  reac_source_local(26,379) = - reac_rate_local(379) 
  reac_source_local(43,379) = + reac_rate_local(379) 
  reac_source_local(44,379) = + reac_rate_local(379) 
  reac_source_local(47,379) = - reac_rate_local(379) 
  reac_source_local(23,380) = + reac_rate_local(380) 
  reac_source_local(28,380) = - reac_rate_local(380) 
  reac_source_local(43,380) = + reac_rate_local(380) 
  reac_source_local(44,380) = + reac_rate_local(380) 
  reac_source_local(47,380) = - reac_rate_local(380) 
  reac_source_local(23,381) = + reac_rate_local(381) 
  reac_source_local(27,381) = - reac_rate_local(381) 
  reac_source_local(43,381) = + reac_rate_local(381) 
  reac_source_local(44,381) = + reac_rate_local(381) 
  reac_source_local(47,381) = - reac_rate_local(381) 
  reac_source_local(19,382) = + reac_rate_local(382) 
  reac_source_local(26,382) = - reac_rate_local(382) 
  reac_source_local(44,382) = + reac_rate_local(382) * 2.d0
  reac_source_local(47,382) = - reac_rate_local(382) 
  reac_source_local(19,383) = + reac_rate_local(383) 
  reac_source_local(28,383) = - reac_rate_local(383) 
  reac_source_local(44,383) = + reac_rate_local(383) * 2.d0
  reac_source_local(47,383) = - reac_rate_local(383) 
  reac_source_local(19,384) = + reac_rate_local(384) 
  reac_source_local(27,384) = - reac_rate_local(384) 
  reac_source_local(44,384) = + reac_rate_local(384) * 2.d0
  reac_source_local(47,384) = - reac_rate_local(384) 
  reac_source_local(23,385) = + reac_rate_local(385) 
  reac_source_local(24,385) = - reac_rate_local(385) 
  reac_source_local(44,385) = + reac_rate_local(385) 
  reac_source_local(47,385) = - reac_rate_local(385) 
  reac_source_local(19,386) = + reac_rate_local(386) 
  reac_source_local(24,386) = - reac_rate_local(386) 
  reac_source_local(43,386) = + reac_rate_local(386) 
  reac_source_local(44,386) = + reac_rate_local(386) 
  reac_source_local(47,386) = - reac_rate_local(386) 
  reac_source_local(20,387) = - reac_rate_local(387) 
  reac_source_local(23,387) = + reac_rate_local(387) 
  reac_source_local(43,387) = + reac_rate_local(387) 
  reac_source_local(47,387) = - reac_rate_local(387) 
  reac_source_local(21,388) = - reac_rate_local(388) 
  reac_source_local(23,388) = + reac_rate_local(388) 
  reac_source_local(43,388) = + reac_rate_local(388) 
  reac_source_local(47,388) = - reac_rate_local(388) 
  reac_source_local(22,389) = - reac_rate_local(389) 
  reac_source_local(23,389) = + reac_rate_local(389) 
  reac_source_local(43,389) = + reac_rate_local(389) 
  reac_source_local(47,389) = - reac_rate_local(389) 
  reac_source_local(19,390) = + reac_rate_local(390) 
  reac_source_local(20,390) = - reac_rate_local(390) 
  reac_source_local(44,390) = + reac_rate_local(390) 
  reac_source_local(47,390) = - reac_rate_local(390) 
  reac_source_local(19,391) = + reac_rate_local(391) 
  reac_source_local(21,391) = - reac_rate_local(391) 
  reac_source_local(44,391) = + reac_rate_local(391) 
  reac_source_local(47,391) = - reac_rate_local(391) 
  reac_source_local(19,392) = + reac_rate_local(392) 
  reac_source_local(22,392) = - reac_rate_local(392) 
  reac_source_local(44,392) = + reac_rate_local(392) 
  reac_source_local(47,392) = - reac_rate_local(392) 
  reac_source_local(17,393) = + reac_rate_local(393) 
  reac_source_local(20,393) = - reac_rate_local(393) 
  reac_source_local(43,393) = + reac_rate_local(393) 
  reac_source_local(44,393) = + reac_rate_local(393) 
  reac_source_local(47,393) = - reac_rate_local(393) 
  reac_source_local(18,394) = - reac_rate_local(394) 
  reac_source_local(19,394) = + reac_rate_local(394) 
  reac_source_local(43,394) = + reac_rate_local(394) 
  reac_source_local(47,394) = - reac_rate_local(394) 
  reac_source_local(17,395) = + reac_rate_local(395) 
  reac_source_local(18,395) = - reac_rate_local(395) 
  reac_source_local(44,395) = + reac_rate_local(395) 
  reac_source_local(47,395) = - reac_rate_local(395) 
  reac_source_local(13,396) = - reac_rate_local(396) 
  reac_source_local(17,396) = + reac_rate_local(396) 
  reac_source_local(43,396) = + reac_rate_local(396) 
  reac_source_local(47,396) = - reac_rate_local(396) 
  reac_source_local(16,397) = - reac_rate_local(397) 
  reac_source_local(17,397) = + reac_rate_local(397) 
  reac_source_local(43,397) = + reac_rate_local(397) 
  reac_source_local(47,397) = - reac_rate_local(397) 
  reac_source_local(15,398) = - reac_rate_local(398) 
  reac_source_local(17,398) = + reac_rate_local(398) 
  reac_source_local(43,398) = + reac_rate_local(398) 
  reac_source_local(47,398) = - reac_rate_local(398) 
  reac_source_local(14,399) = - reac_rate_local(399) 
  reac_source_local(17,399) = + reac_rate_local(399) 
  reac_source_local(43,399) = + reac_rate_local(399) 
  reac_source_local(47,399) = - reac_rate_local(399) 
  reac_source_local(05,400) = - reac_rate_local(400) 
  reac_source_local(07,400) = + reac_rate_local(400) * 2.d0
  reac_source_local(09,400) = - reac_rate_local(400) 
  reac_source_local(05,401) = - reac_rate_local(401) 
  reac_source_local(07,401) = + reac_rate_local(401) * 2.d0
  reac_source_local(11,401) = - reac_rate_local(401) 
  reac_source_local(05,402) = - reac_rate_local(402) 
  reac_source_local(07,402) = + reac_rate_local(402) * 2.d0
  reac_source_local(10,402) = - reac_rate_local(402) 
  reac_source_local(03,403) = - reac_rate_local(403) 
  reac_source_local(09,403) = - reac_rate_local(403) 
  reac_source_local(20,403) = + reac_rate_local(403) 
  reac_source_local(43,403) = + reac_rate_local(403) 
  reac_source_local(03,404) = - reac_rate_local(404) 
  reac_source_local(11,404) = - reac_rate_local(404) 
  reac_source_local(20,404) = + reac_rate_local(404) 
  reac_source_local(43,404) = + reac_rate_local(404) 
  reac_source_local(03,405) = - reac_rate_local(405) 
  reac_source_local(10,405) = - reac_rate_local(405) 
  reac_source_local(20,405) = + reac_rate_local(405) 
  reac_source_local(43,405) = + reac_rate_local(405) 
  reac_source_local(07,406) = + reac_rate_local(406) 
  reac_source_local(09,406) = - reac_rate_local(406) 
  reac_source_local(24,406) = - reac_rate_local(406) 
  reac_source_local(26,406) = + reac_rate_local(406) 
  reac_source_local(07,407) = + reac_rate_local(407) 
  reac_source_local(11,407) = - reac_rate_local(407) 
  reac_source_local(24,407) = - reac_rate_local(407) 
  reac_source_local(26,407) = + reac_rate_local(407) 
  reac_source_local(07,408) = + reac_rate_local(408) 
  reac_source_local(10,408) = - reac_rate_local(408) 
  reac_source_local(24,408) = - reac_rate_local(408) 
  reac_source_local(26,408) = + reac_rate_local(408) 
  reac_source_local(07,409) = + reac_rate_local(409) 
  reac_source_local(09,409) = - reac_rate_local(409) 
  reac_source_local(18,409) = - reac_rate_local(409) 
  reac_source_local(20,409) = + reac_rate_local(409) 
  reac_source_local(07,410) = + reac_rate_local(410) 
  reac_source_local(11,410) = - reac_rate_local(410) 
  reac_source_local(18,410) = - reac_rate_local(410) 
  reac_source_local(20,410) = + reac_rate_local(410) 
  reac_source_local(07,411) = + reac_rate_local(411) 
  reac_source_local(10,411) = - reac_rate_local(411) 
  reac_source_local(18,411) = - reac_rate_local(411) 
  reac_source_local(20,411) = + reac_rate_local(411) 
  reac_source_local(07,412) = + reac_rate_local(412) 
  reac_source_local(09,412) = - reac_rate_local(412) 
  reac_source_local(34,412) = - reac_rate_local(412) 
  reac_source_local(36,412) = + reac_rate_local(412) 
  reac_source_local(07,413) = + reac_rate_local(413) 
  reac_source_local(11,413) = - reac_rate_local(413) 
  reac_source_local(34,413) = - reac_rate_local(413) 
  reac_source_local(36,413) = + reac_rate_local(413) 
  reac_source_local(07,414) = + reac_rate_local(414) 
  reac_source_local(10,414) = - reac_rate_local(414) 
  reac_source_local(34,414) = - reac_rate_local(414) 
  reac_source_local(36,414) = + reac_rate_local(414) 
  reac_source_local(07,415) = + reac_rate_local(415) 
  reac_source_local(09,415) = - reac_rate_local(415) 
  reac_source_local(43,415) = - reac_rate_local(415) 
  reac_source_local(44,415) = + reac_rate_local(415) 
  reac_source_local(07,416) = + reac_rate_local(416) 
  reac_source_local(11,416) = - reac_rate_local(416) 
  reac_source_local(43,416) = - reac_rate_local(416) 
  reac_source_local(44,416) = + reac_rate_local(416) 
  reac_source_local(07,417) = + reac_rate_local(417) 
  reac_source_local(10,417) = - reac_rate_local(417) 
  reac_source_local(43,417) = - reac_rate_local(417) 
  reac_source_local(44,417) = + reac_rate_local(417) 
  reac_source_local(07,418) = - reac_rate_local(418) 
  reac_source_local(09,418) = - reac_rate_local(418) 
  reac_source_local(26,418) = + reac_rate_local(418) 
  reac_source_local(43,418) = + reac_rate_local(418) 
  reac_source_local(07,419) = - reac_rate_local(419) 
  reac_source_local(11,419) = - reac_rate_local(419) 
  reac_source_local(26,419) = + reac_rate_local(419) 
  reac_source_local(43,419) = + reac_rate_local(419) 
  reac_source_local(07,420) = - reac_rate_local(420) 
  reac_source_local(10,420) = - reac_rate_local(420) 
  reac_source_local(26,420) = + reac_rate_local(420) 
  reac_source_local(43,420) = + reac_rate_local(420) 
  reac_source_local(05,421) = - reac_rate_local(421) 
  reac_source_local(09,421) = - reac_rate_local(421) 
  reac_source_local(26,421) = + reac_rate_local(421) 
  reac_source_local(05,422) = - reac_rate_local(422) 
  reac_source_local(11,422) = - reac_rate_local(422) 
  reac_source_local(26,422) = + reac_rate_local(422) 
  reac_source_local(05,423) = - reac_rate_local(423) 
  reac_source_local(10,423) = - reac_rate_local(423) 
  reac_source_local(26,423) = + reac_rate_local(423) 
  reac_source_local(07,424) = + reac_rate_local(424) 
  reac_source_local(09,424) = - reac_rate_local(424) 
  reac_source_local(43,424) = + reac_rate_local(424) 
  reac_source_local(07,425) = + reac_rate_local(425) 
  reac_source_local(11,425) = - reac_rate_local(425) 
  reac_source_local(43,425) = + reac_rate_local(425) 
  reac_source_local(07,426) = + reac_rate_local(426) 
  reac_source_local(10,426) = - reac_rate_local(426) 
  reac_source_local(43,426) = + reac_rate_local(426) 
  reac_source_local(05,427) = + reac_rate_local(427) 
  reac_source_local(07,427) = - reac_rate_local(427) 
  reac_source_local(43,427) = + reac_rate_local(427) 
  reac_source_local(03,428) = + reac_rate_local(428) 
  reac_source_local(07,428) = - reac_rate_local(428) 
  reac_source_local(44,428) = + reac_rate_local(428) 
  reac_source_local(05,429) = + reac_rate_local(429) 
  reac_source_local(07,429) = - reac_rate_local(429) 
  reac_source_local(24,429) = - reac_rate_local(429) 
  reac_source_local(26,429) = + reac_rate_local(429) 
  reac_source_local(05,430) = - reac_rate_local(430) * 2.d0
  reac_source_local(13,430) = + reac_rate_local(430) 
  reac_source_local(44,430) = + reac_rate_local(430) 
  reac_source_local(05,431) = - reac_rate_local(431) 
  reac_source_local(07,431) = + reac_rate_local(431) 
  reac_source_local(20,431) = + reac_rate_local(431) 
  reac_source_local(24,431) = - reac_rate_local(431) 
  reac_source_local(05,432) = - reac_rate_local(432) 
  reac_source_local(07,432) = + reac_rate_local(432) 
  reac_source_local(13,432) = + reac_rate_local(432) 
  reac_source_local(18,432) = - reac_rate_local(432) 
  reac_source_local(05,433) = - reac_rate_local(433) 
  reac_source_local(07,433) = + reac_rate_local(433) 
  reac_source_local(34,433) = + reac_rate_local(433) 
  reac_source_local(36,433) = - reac_rate_local(433) 
  reac_source_local(05,434) = - reac_rate_local(434) 
  reac_source_local(07,434) = + reac_rate_local(434) 
  reac_source_local(34,434) = + reac_rate_local(434) 
  reac_source_local(37,434) = - reac_rate_local(434) 
  reac_source_local(05,435) = - reac_rate_local(435) 
  reac_source_local(07,435) = + reac_rate_local(435) 
  reac_source_local(34,435) = + reac_rate_local(435) 
  reac_source_local(38,435) = - reac_rate_local(435) 
  reac_source_local(05,436) = - reac_rate_local(436) 
  reac_source_local(20,436) = + reac_rate_local(436) 
  reac_source_local(24,436) = + reac_rate_local(436) 
  reac_source_local(34,436) = - reac_rate_local(436) 
  reac_source_local(05,437) = - reac_rate_local(437) 
  reac_source_local(07,437) = + reac_rate_local(437) 
  reac_source_local(31,437) = + reac_rate_local(437) 
  reac_source_local(34,437) = - reac_rate_local(437) 
  reac_source_local(05,438) = - reac_rate_local(438) 
  reac_source_local(07,438) = + reac_rate_local(438) 
  reac_source_local(43,438) = + reac_rate_local(438) 
  reac_source_local(44,438) = - reac_rate_local(438) 
  reac_source_local(03,439) = + reac_rate_local(439) 
  reac_source_local(05,439) = - reac_rate_local(439) 
  reac_source_local(43,439) = - reac_rate_local(439) 
  reac_source_local(44,439) = + reac_rate_local(439) 
  reac_source_local(03,440) = + reac_rate_local(440) 
  reac_source_local(05,440) = - reac_rate_local(440) 
  reac_source_local(43,440) = + reac_rate_local(440) 
  reac_source_local(05,441) = - reac_rate_local(441) 
  reac_source_local(07,441) = + reac_rate_local(441) 
  reac_source_local(43,441) = - reac_rate_local(441) 
  reac_source_local(03,442) = - reac_rate_local(442) 
  reac_source_local(26,442) = - reac_rate_local(442) 
  reac_source_local(31,442) = + reac_rate_local(442) 
  reac_source_local(43,442) = + reac_rate_local(442) 
  reac_source_local(03,443) = - reac_rate_local(443) 
  reac_source_local(28,443) = - reac_rate_local(443) 
  reac_source_local(31,443) = + reac_rate_local(443) 
  reac_source_local(43,443) = + reac_rate_local(443) 
  reac_source_local(03,444) = - reac_rate_local(444) 
  reac_source_local(27,444) = - reac_rate_local(444) 
  reac_source_local(31,444) = + reac_rate_local(444) 
  reac_source_local(43,444) = + reac_rate_local(444) 
  reac_source_local(03,445) = - reac_rate_local(445) 
  reac_source_local(26,445) = - reac_rate_local(445) 
  reac_source_local(34,445) = + reac_rate_local(445) 
  reac_source_local(03,446) = - reac_rate_local(446) 
  reac_source_local(28,446) = - reac_rate_local(446) 
  reac_source_local(34,446) = + reac_rate_local(446) 
  reac_source_local(03,447) = - reac_rate_local(447) 
  reac_source_local(27,447) = - reac_rate_local(447) 
  reac_source_local(34,447) = + reac_rate_local(447) 
  reac_source_local(03,448) = - reac_rate_local(448) 
  reac_source_local(05,448) = + reac_rate_local(448) 
  reac_source_local(43,448) = + reac_rate_local(448) 
  reac_source_local(44,448) = - reac_rate_local(448) 
  reac_source_local(03,449) = - reac_rate_local(449) 
  reac_source_local(07,449) = - reac_rate_local(449) 
  reac_source_local(18,449) = + reac_rate_local(449) 
  reac_source_local(43,449) = + reac_rate_local(449) 
  reac_source_local(03,450) = - reac_rate_local(450) 
  reac_source_local(05,450) = - reac_rate_local(450) 
  reac_source_local(13,450) = + reac_rate_local(450) 
  reac_source_local(43,450) = + reac_rate_local(450) 
  reac_source_local(03,451) = - reac_rate_local(451) 
  reac_source_local(07,451) = + reac_rate_local(451) 
  reac_source_local(44,451) = - reac_rate_local(451) 
  reac_source_local(03,452) = - reac_rate_local(452) 
  reac_source_local(05,452) = + reac_rate_local(452) 
  reac_source_local(13,452) = + reac_rate_local(452) 
  reac_source_local(18,452) = - reac_rate_local(452) 
  reac_source_local(18,453) = - reac_rate_local(453) 
  reac_source_local(20,453) = + reac_rate_local(453) 
  reac_source_local(24,453) = + reac_rate_local(453) 
  reac_source_local(26,453) = - reac_rate_local(453) 
  reac_source_local(18,454) = - reac_rate_local(454) 
  reac_source_local(20,454) = + reac_rate_local(454) 
  reac_source_local(24,454) = + reac_rate_local(454) 
  reac_source_local(28,454) = - reac_rate_local(454) 
  reac_source_local(18,455) = - reac_rate_local(455) 
  reac_source_local(20,455) = + reac_rate_local(455) 
  reac_source_local(24,455) = + reac_rate_local(455) 
  reac_source_local(27,455) = - reac_rate_local(455) 
  reac_source_local(24,456) = + reac_rate_local(456) 
  reac_source_local(26,456) = - reac_rate_local(456) 
  reac_source_local(34,456) = - reac_rate_local(456) 
  reac_source_local(36,456) = + reac_rate_local(456) 
  reac_source_local(24,457) = + reac_rate_local(457) 
  reac_source_local(28,457) = - reac_rate_local(457) 
  reac_source_local(34,457) = - reac_rate_local(457) 
  reac_source_local(36,457) = + reac_rate_local(457) 
  reac_source_local(24,458) = + reac_rate_local(458) 
  reac_source_local(27,458) = - reac_rate_local(458) 
  reac_source_local(34,458) = - reac_rate_local(458) 
  reac_source_local(36,458) = + reac_rate_local(458) 
  reac_source_local(24,459) = + reac_rate_local(459) 
  reac_source_local(26,459) = - reac_rate_local(459) 
  reac_source_local(43,459) = - reac_rate_local(459) 
  reac_source_local(44,459) = + reac_rate_local(459) 
  reac_source_local(24,460) = + reac_rate_local(460) 
  reac_source_local(28,460) = - reac_rate_local(460) 
  reac_source_local(43,460) = - reac_rate_local(460) 
  reac_source_local(44,460) = + reac_rate_local(460) 
  reac_source_local(24,461) = + reac_rate_local(461) 
  reac_source_local(27,461) = - reac_rate_local(461) 
  reac_source_local(43,461) = - reac_rate_local(461) 
  reac_source_local(44,461) = + reac_rate_local(461) 
  reac_source_local(07,462) = + reac_rate_local(462) 
  reac_source_local(09,462) = + reac_rate_local(462) 
  reac_source_local(26,462) = - reac_rate_local(462) 
  reac_source_local(43,462) = - reac_rate_local(462) 
  reac_source_local(07,463) = + reac_rate_local(463) 
  reac_source_local(09,463) = + reac_rate_local(463) 
  reac_source_local(28,463) = - reac_rate_local(463) 
  reac_source_local(43,463) = - reac_rate_local(463) 
  reac_source_local(07,464) = + reac_rate_local(464) 
  reac_source_local(09,464) = + reac_rate_local(464) 
  reac_source_local(27,464) = - reac_rate_local(464) 
  reac_source_local(43,464) = - reac_rate_local(464) 
  reac_source_local(07,465) = + reac_rate_local(465) * 2.d0
  reac_source_local(26,465) = - reac_rate_local(465) 
  reac_source_local(07,466) = + reac_rate_local(466) * 2.d0
  reac_source_local(28,466) = - reac_rate_local(466) 
  reac_source_local(07,467) = + reac_rate_local(467) * 2.d0
  reac_source_local(27,467) = - reac_rate_local(467) 
  reac_source_local(03,468) = - reac_rate_local(468) 
  reac_source_local(07,468) = + reac_rate_local(468) 
  reac_source_local(20,468) = + reac_rate_local(468) 
  reac_source_local(26,468) = - reac_rate_local(468) 
  reac_source_local(03,469) = - reac_rate_local(469) 
  reac_source_local(07,469) = + reac_rate_local(469) 
  reac_source_local(20,469) = + reac_rate_local(469) 
  reac_source_local(28,469) = - reac_rate_local(469) 
  reac_source_local(03,470) = - reac_rate_local(470) 
  reac_source_local(07,470) = + reac_rate_local(470) 
  reac_source_local(20,470) = + reac_rate_local(470) 
  reac_source_local(27,470) = - reac_rate_local(470) 
  reac_source_local(05,471) = - reac_rate_local(471) 
  reac_source_local(07,471) = + reac_rate_local(471) 
  reac_source_local(24,471) = + reac_rate_local(471) 
  reac_source_local(26,471) = - reac_rate_local(471) 
  reac_source_local(05,472) = - reac_rate_local(472) 
  reac_source_local(07,472) = + reac_rate_local(472) 
  reac_source_local(24,472) = + reac_rate_local(472) 
  reac_source_local(28,472) = - reac_rate_local(472) 
  reac_source_local(05,473) = - reac_rate_local(473) 
  reac_source_local(07,473) = + reac_rate_local(473) 
  reac_source_local(24,473) = + reac_rate_local(473) 
  reac_source_local(27,473) = - reac_rate_local(473) 
  reac_source_local(20,474) = + reac_rate_local(474) 
  reac_source_local(24,474) = - reac_rate_local(474) * 2.d0
  reac_source_local(26,474) = + reac_rate_local(474) 
  reac_source_local(18,475) = + reac_rate_local(475) 
  reac_source_local(20,475) = - reac_rate_local(475) 
  reac_source_local(24,475) = - reac_rate_local(475) 
  reac_source_local(26,475) = + reac_rate_local(475) 
  reac_source_local(18,476) = + reac_rate_local(476) 
  reac_source_local(21,476) = - reac_rate_local(476) 
  reac_source_local(24,476) = - reac_rate_local(476) 
  reac_source_local(26,476) = + reac_rate_local(476) 
  reac_source_local(18,477) = + reac_rate_local(477) 
  reac_source_local(22,477) = - reac_rate_local(477) 
  reac_source_local(24,477) = - reac_rate_local(477) 
  reac_source_local(26,477) = + reac_rate_local(477) 
  reac_source_local(24,478) = - reac_rate_local(478) 
  reac_source_local(26,478) = + reac_rate_local(478) 
  reac_source_local(34,478) = + reac_rate_local(478) 
  reac_source_local(36,478) = - reac_rate_local(478) 
  reac_source_local(24,479) = - reac_rate_local(479) 
  reac_source_local(26,479) = + reac_rate_local(479) 
  reac_source_local(34,479) = + reac_rate_local(479) 
  reac_source_local(37,479) = - reac_rate_local(479) 
  reac_source_local(24,480) = - reac_rate_local(480) 
  reac_source_local(26,480) = + reac_rate_local(480) 
  reac_source_local(34,480) = + reac_rate_local(480) 
  reac_source_local(38,480) = - reac_rate_local(480) 
  reac_source_local(20,481) = + reac_rate_local(481) 
  reac_source_local(24,481) = - reac_rate_local(481) 
  reac_source_local(34,481) = - reac_rate_local(481) 
  reac_source_local(36,481) = + reac_rate_local(481) 
  reac_source_local(24,482) = - reac_rate_local(482) 
  reac_source_local(26,482) = + reac_rate_local(482) 
  reac_source_local(31,482) = + reac_rate_local(482) 
  reac_source_local(34,482) = - reac_rate_local(482) 
  reac_source_local(24,483) = - reac_rate_local(483) 
  reac_source_local(26,483) = + reac_rate_local(483) 
  reac_source_local(30,483) = + reac_rate_local(483) 
  reac_source_local(31,483) = - reac_rate_local(483) 
  reac_source_local(24,484) = - reac_rate_local(484) 
  reac_source_local(26,484) = + reac_rate_local(484) 
  reac_source_local(30,484) = + reac_rate_local(484) 
  reac_source_local(32,484) = - reac_rate_local(484) 
  reac_source_local(24,485) = - reac_rate_local(485) 
  reac_source_local(26,485) = + reac_rate_local(485) 
  reac_source_local(43,485) = + reac_rate_local(485) 
  reac_source_local(44,485) = - reac_rate_local(485) 
  reac_source_local(07,486) = + reac_rate_local(486) * 2.d0
  reac_source_local(24,486) = - reac_rate_local(486) 
  reac_source_local(43,486) = - reac_rate_local(486) 
  reac_source_local(20,487) = + reac_rate_local(487) 
  reac_source_local(24,487) = - reac_rate_local(487) 
  reac_source_local(43,487) = - reac_rate_local(487) 
  reac_source_local(44,487) = + reac_rate_local(487) 
  reac_source_local(24,488) = - reac_rate_local(488) 
  reac_source_local(26,488) = + reac_rate_local(488) 
  reac_source_local(43,488) = - reac_rate_local(488) 
  reac_source_local(20,489) = + reac_rate_local(489) 
  reac_source_local(24,489) = - reac_rate_local(489) 
  reac_source_local(43,489) = + reac_rate_local(489) 
  reac_source_local(24,490) = - reac_rate_local(490) * 2.d0
  reac_source_local(41,490) = + reac_rate_local(490) 
  reac_source_local(18,491) = - reac_rate_local(491) 
  reac_source_local(20,491) = + reac_rate_local(491) * 2.d0
  reac_source_local(24,491) = - reac_rate_local(491) 
  reac_source_local(18,492) = + reac_rate_local(492) 
  reac_source_local(20,492) = - reac_rate_local(492) 
  reac_source_local(43,492) = - reac_rate_local(492) 
  reac_source_local(44,492) = + reac_rate_local(492) 
  reac_source_local(18,493) = + reac_rate_local(493) 
  reac_source_local(21,493) = - reac_rate_local(493) 
  reac_source_local(43,493) = - reac_rate_local(493) 
  reac_source_local(44,493) = + reac_rate_local(493) 
  reac_source_local(18,494) = + reac_rate_local(494) 
  reac_source_local(22,494) = - reac_rate_local(494) 
  reac_source_local(43,494) = - reac_rate_local(494) 
  reac_source_local(44,494) = + reac_rate_local(494) 
  reac_source_local(20,495) = - reac_rate_local(495) 
  reac_source_local(24,495) = + reac_rate_local(495) 
  reac_source_local(43,495) = - reac_rate_local(495) 
  reac_source_local(21,496) = - reac_rate_local(496) 
  reac_source_local(24,496) = + reac_rate_local(496) 
  reac_source_local(43,496) = - reac_rate_local(496) 
  reac_source_local(22,497) = - reac_rate_local(497) 
  reac_source_local(24,497) = + reac_rate_local(497) 
  reac_source_local(43,497) = - reac_rate_local(497) 
  reac_source_local(20,498) = - reac_rate_local(498) 
  reac_source_local(24,498) = + reac_rate_local(498) 
  reac_source_local(43,498) = + reac_rate_local(498) 
  reac_source_local(44,498) = - reac_rate_local(498) 
  reac_source_local(18,499) = + reac_rate_local(499) 
  reac_source_local(20,499) = - reac_rate_local(499) 
  reac_source_local(43,499) = + reac_rate_local(499) 
  reac_source_local(18,500) = + reac_rate_local(500) 
  reac_source_local(21,500) = - reac_rate_local(500) 
  reac_source_local(43,500) = + reac_rate_local(500) 
  reac_source_local(18,501) = + reac_rate_local(501) 
  reac_source_local(22,501) = - reac_rate_local(501) 
  reac_source_local(43,501) = + reac_rate_local(501) 
  reac_source_local(13,502) = - reac_rate_local(502) 
  reac_source_local(18,502) = + reac_rate_local(502) * 2.d0
  reac_source_local(20,502) = - reac_rate_local(502) 
  reac_source_local(13,503) = - reac_rate_local(503) 
  reac_source_local(18,503) = + reac_rate_local(503) * 2.d0
  reac_source_local(21,503) = - reac_rate_local(503) 
  reac_source_local(13,504) = - reac_rate_local(504) 
  reac_source_local(18,504) = + reac_rate_local(504) * 2.d0
  reac_source_local(22,504) = - reac_rate_local(504) 
  reac_source_local(16,505) = - reac_rate_local(505) 
  reac_source_local(18,505) = + reac_rate_local(505) * 2.d0
  reac_source_local(20,505) = - reac_rate_local(505) 
  reac_source_local(15,506) = - reac_rate_local(506) 
  reac_source_local(18,506) = + reac_rate_local(506) * 2.d0
  reac_source_local(20,506) = - reac_rate_local(506) 
  reac_source_local(14,507) = - reac_rate_local(507) 
  reac_source_local(18,507) = + reac_rate_local(507) * 2.d0
  reac_source_local(20,507) = - reac_rate_local(507) 
  reac_source_local(18,508) = + reac_rate_local(508) 
  reac_source_local(20,508) = - reac_rate_local(508) 
  reac_source_local(31,508) = - reac_rate_local(508) 
  reac_source_local(34,508) = + reac_rate_local(508) 
  reac_source_local(18,509) = + reac_rate_local(509) 
  reac_source_local(21,509) = - reac_rate_local(509) 
  reac_source_local(31,509) = - reac_rate_local(509) 
  reac_source_local(34,509) = + reac_rate_local(509) 
  reac_source_local(18,510) = + reac_rate_local(510) 
  reac_source_local(22,510) = - reac_rate_local(510) 
  reac_source_local(31,510) = - reac_rate_local(510) 
  reac_source_local(34,510) = + reac_rate_local(510) 
  reac_source_local(18,511) = + reac_rate_local(511) 
  reac_source_local(20,511) = - reac_rate_local(511) 
  reac_source_local(32,511) = - reac_rate_local(511) 
  reac_source_local(34,511) = + reac_rate_local(511) 
  reac_source_local(18,512) = + reac_rate_local(512) 
  reac_source_local(20,512) = - reac_rate_local(512) * 2.d0
  reac_source_local(24,512) = + reac_rate_local(512) 
  reac_source_local(18,513) = + reac_rate_local(513) 
  reac_source_local(21,513) = - reac_rate_local(513) * 2.d0
  reac_source_local(24,513) = + reac_rate_local(513) 
  reac_source_local(18,514) = + reac_rate_local(514) 
  reac_source_local(22,514) = - reac_rate_local(514) * 2.d0
  reac_source_local(24,514) = + reac_rate_local(514) 
  reac_source_local(07,515) = - reac_rate_local(515) 
  reac_source_local(20,515) = - reac_rate_local(515) 
  reac_source_local(34,515) = + reac_rate_local(515) 
  reac_source_local(07,516) = - reac_rate_local(516) 
  reac_source_local(21,516) = - reac_rate_local(516) 
  reac_source_local(34,516) = + reac_rate_local(516) 
  reac_source_local(07,517) = - reac_rate_local(517) 
  reac_source_local(22,517) = - reac_rate_local(517) 
  reac_source_local(34,517) = + reac_rate_local(517) 
  reac_source_local(13,518) = + reac_rate_local(518) 
  reac_source_local(20,518) = - reac_rate_local(518) 
  reac_source_local(44,518) = + reac_rate_local(518) 
  reac_source_local(13,519) = + reac_rate_local(519) 
  reac_source_local(21,519) = - reac_rate_local(519) 
  reac_source_local(44,519) = + reac_rate_local(519) 
  reac_source_local(13,520) = + reac_rate_local(520) 
  reac_source_local(22,520) = - reac_rate_local(520) 
  reac_source_local(44,520) = + reac_rate_local(520) 
  reac_source_local(20,521) = - reac_rate_local(521) 
  reac_source_local(26,521) = + reac_rate_local(521) 
  reac_source_local(44,521) = - reac_rate_local(521) 
  reac_source_local(21,522) = - reac_rate_local(522) 
  reac_source_local(26,522) = + reac_rate_local(522) 
  reac_source_local(44,522) = - reac_rate_local(522) 
  reac_source_local(22,523) = - reac_rate_local(523) 
  reac_source_local(26,523) = + reac_rate_local(523) 
  reac_source_local(44,523) = - reac_rate_local(523) 
  reac_source_local(05,524) = - reac_rate_local(524) 
  reac_source_local(20,524) = - reac_rate_local(524) 
  reac_source_local(31,524) = + reac_rate_local(524) 
  reac_source_local(05,525) = - reac_rate_local(525) 
  reac_source_local(21,525) = - reac_rate_local(525) 
  reac_source_local(31,525) = + reac_rate_local(525) 
  reac_source_local(05,526) = - reac_rate_local(526) 
  reac_source_local(22,526) = - reac_rate_local(526) 
  reac_source_local(31,526) = + reac_rate_local(526) 
  reac_source_local(13,527) = + reac_rate_local(527) 
  reac_source_local(18,527) = - reac_rate_local(527) * 2.d0
  reac_source_local(20,527) = + reac_rate_local(527) 
  reac_source_local(18,528) = - reac_rate_local(528) 
  reac_source_local(20,528) = + reac_rate_local(528) 
  reac_source_local(34,528) = + reac_rate_local(528) 
  reac_source_local(36,528) = - reac_rate_local(528) 
  reac_source_local(18,529) = - reac_rate_local(529) 
  reac_source_local(20,529) = + reac_rate_local(529) 
  reac_source_local(34,529) = + reac_rate_local(529) 
  reac_source_local(37,529) = - reac_rate_local(529) 
  reac_source_local(18,530) = - reac_rate_local(530) 
  reac_source_local(20,530) = + reac_rate_local(530) 
  reac_source_local(34,530) = + reac_rate_local(530) 
  reac_source_local(38,530) = - reac_rate_local(530) 
  reac_source_local(13,531) = + reac_rate_local(531) 
  reac_source_local(18,531) = - reac_rate_local(531) 
  reac_source_local(34,531) = - reac_rate_local(531) 
  reac_source_local(36,531) = + reac_rate_local(531) 
  reac_source_local(18,532) = - reac_rate_local(532) 
  reac_source_local(20,532) = + reac_rate_local(532) 
  reac_source_local(31,532) = + reac_rate_local(532) 
  reac_source_local(34,532) = - reac_rate_local(532) 
  reac_source_local(18,533) = - reac_rate_local(533) 
  reac_source_local(20,533) = + reac_rate_local(533) 
  reac_source_local(43,533) = + reac_rate_local(533) 
  reac_source_local(44,533) = - reac_rate_local(533) 
  reac_source_local(13,534) = + reac_rate_local(534) 
  reac_source_local(18,534) = - reac_rate_local(534) 
  reac_source_local(43,534) = - reac_rate_local(534) 
  reac_source_local(44,534) = + reac_rate_local(534) 
  reac_source_local(18,535) = - reac_rate_local(535) 
  reac_source_local(20,535) = + reac_rate_local(535) 
  reac_source_local(43,535) = - reac_rate_local(535) 
  reac_source_local(13,536) = + reac_rate_local(536) 
  reac_source_local(18,536) = - reac_rate_local(536) 
  reac_source_local(43,536) = + reac_rate_local(536) 
  reac_source_local(13,537) = - reac_rate_local(537) 
  reac_source_local(18,537) = + reac_rate_local(537) 
  reac_source_local(43,537) = - reac_rate_local(537) 
  reac_source_local(16,538) = - reac_rate_local(538) 
  reac_source_local(18,538) = + reac_rate_local(538) 
  reac_source_local(43,538) = - reac_rate_local(538) 
  reac_source_local(15,539) = - reac_rate_local(539) 
  reac_source_local(18,539) = + reac_rate_local(539) 
  reac_source_local(43,539) = - reac_rate_local(539) 
  reac_source_local(14,540) = - reac_rate_local(540) 
  reac_source_local(18,540) = + reac_rate_local(540) 
  reac_source_local(43,540) = - reac_rate_local(540) 
  reac_source_local(13,541) = - reac_rate_local(541) 
  reac_source_local(20,541) = + reac_rate_local(541) 
  reac_source_local(44,541) = - reac_rate_local(541) 
  reac_source_local(16,542) = - reac_rate_local(542) 
  reac_source_local(20,542) = + reac_rate_local(542) 
  reac_source_local(44,542) = - reac_rate_local(542) 
  reac_source_local(15,543) = - reac_rate_local(543) 
  reac_source_local(20,543) = + reac_rate_local(543) 
  reac_source_local(44,543) = - reac_rate_local(543) 
  reac_source_local(14,544) = - reac_rate_local(544) 
  reac_source_local(20,544) = + reac_rate_local(544) 
  reac_source_local(44,544) = - reac_rate_local(544) 
  reac_source_local(13,545) = - reac_rate_local(545) 
  reac_source_local(18,545) = + reac_rate_local(545) 
  reac_source_local(43,545) = + reac_rate_local(545) 
  reac_source_local(44,545) = - reac_rate_local(545) 
  reac_source_local(16,546) = - reac_rate_local(546) 
  reac_source_local(18,546) = + reac_rate_local(546) 
  reac_source_local(43,546) = + reac_rate_local(546) 
  reac_source_local(44,546) = - reac_rate_local(546) 
  reac_source_local(15,547) = - reac_rate_local(547) 
  reac_source_local(18,547) = + reac_rate_local(547) 
  reac_source_local(43,547) = + reac_rate_local(547) 
  reac_source_local(44,547) = - reac_rate_local(547) 
  reac_source_local(14,548) = - reac_rate_local(548) 
  reac_source_local(18,548) = + reac_rate_local(548) 
  reac_source_local(43,548) = + reac_rate_local(548) 
  reac_source_local(44,548) = - reac_rate_local(548) 
  reac_source_local(07,549) = - reac_rate_local(549) 
  reac_source_local(13,549) = - reac_rate_local(549) 
  reac_source_local(30,549) = + reac_rate_local(549) 
  reac_source_local(07,550) = - reac_rate_local(550) 
  reac_source_local(16,550) = - reac_rate_local(550) 
  reac_source_local(30,550) = + reac_rate_local(550) 
  reac_source_local(07,551) = - reac_rate_local(551) 
  reac_source_local(15,551) = - reac_rate_local(551) 
  reac_source_local(30,551) = + reac_rate_local(551) 
  reac_source_local(07,552) = - reac_rate_local(552) 
  reac_source_local(14,552) = - reac_rate_local(552) 
  reac_source_local(30,552) = + reac_rate_local(552) 
  reac_source_local(34,553) = + reac_rate_local(553) 
  reac_source_local(36,553) = - reac_rate_local(553) 
  reac_source_local(43,553) = - reac_rate_local(553) 
  reac_source_local(44,553) = + reac_rate_local(553) 
  reac_source_local(34,554) = + reac_rate_local(554) 
  reac_source_local(37,554) = - reac_rate_local(554) 
  reac_source_local(43,554) = - reac_rate_local(554) 
  reac_source_local(44,554) = + reac_rate_local(554) 
  reac_source_local(34,555) = + reac_rate_local(555) 
  reac_source_local(38,555) = - reac_rate_local(555) 
  reac_source_local(43,555) = - reac_rate_local(555) 
  reac_source_local(44,555) = + reac_rate_local(555) 
  reac_source_local(07,556) = + reac_rate_local(556) 
  reac_source_local(24,556) = + reac_rate_local(556) 
  reac_source_local(36,556) = - reac_rate_local(556) 
  reac_source_local(07,557) = + reac_rate_local(557) 
  reac_source_local(24,557) = + reac_rate_local(557) 
  reac_source_local(37,557) = - reac_rate_local(557) 
  reac_source_local(07,558) = + reac_rate_local(558) 
  reac_source_local(24,558) = + reac_rate_local(558) 
  reac_source_local(38,558) = - reac_rate_local(558) 
  reac_source_local(05,559) = - reac_rate_local(559) 
  reac_source_local(36,559) = - reac_rate_local(559) 
  reac_source_local(41,559) = + reac_rate_local(559) 
  reac_source_local(05,560) = - reac_rate_local(560) 
  reac_source_local(37,560) = - reac_rate_local(560) 
  reac_source_local(41,560) = + reac_rate_local(560) 
  reac_source_local(05,561) = - reac_rate_local(561) 
  reac_source_local(38,561) = - reac_rate_local(561) 
  reac_source_local(41,561) = + reac_rate_local(561) 
  reac_source_local(31,562) = + reac_rate_local(562) 
  reac_source_local(34,562) = - reac_rate_local(562) * 2.d0
  reac_source_local(36,562) = + reac_rate_local(562) 
  reac_source_local(34,563) = - reac_rate_local(563) 
  reac_source_local(36,563) = + reac_rate_local(563) 
  reac_source_local(43,563) = + reac_rate_local(563) 
  reac_source_local(44,563) = - reac_rate_local(563) 
  reac_source_local(31,564) = + reac_rate_local(564) 
  reac_source_local(34,564) = - reac_rate_local(564) 
  reac_source_local(43,564) = - reac_rate_local(564) 
  reac_source_local(44,564) = + reac_rate_local(564) 
  reac_source_local(34,565) = - reac_rate_local(565) 
  reac_source_local(36,565) = + reac_rate_local(565) 
  reac_source_local(43,565) = - reac_rate_local(565) 
  reac_source_local(07,566) = + reac_rate_local(566) 
  reac_source_local(24,566) = + reac_rate_local(566) 
  reac_source_local(34,566) = - reac_rate_local(566) 
  reac_source_local(43,566) = - reac_rate_local(566) 
  reac_source_local(31,567) = + reac_rate_local(567) 
  reac_source_local(34,567) = - reac_rate_local(567) 
  reac_source_local(43,567) = + reac_rate_local(567) 
  reac_source_local(07,568) = + reac_rate_local(568) 
  reac_source_local(20,568) = + reac_rate_local(568) 
  reac_source_local(34,568) = - reac_rate_local(568) 
  reac_source_local(31,569) = - reac_rate_local(569) 
  reac_source_local(34,569) = + reac_rate_local(569) 
  reac_source_local(43,569) = - reac_rate_local(569) 
  reac_source_local(32,570) = - reac_rate_local(570) 
  reac_source_local(34,570) = + reac_rate_local(570) 
  reac_source_local(43,570) = - reac_rate_local(570) 
  reac_source_local(07,571) = + reac_rate_local(571) 
  reac_source_local(18,571) = + reac_rate_local(571) 
  reac_source_local(31,571) = - reac_rate_local(571) 
  reac_source_local(07,572) = - reac_rate_local(572) 
  reac_source_local(09,572) = + reac_rate_local(572) 
  reac_source_local(40,572) = + reac_rate_local(572) 
  reac_source_local(41,572) = - reac_rate_local(572) 
  reac_source_local(07,573) = + reac_rate_local(573) 
  reac_source_local(34,573) = + reac_rate_local(573) 
  reac_source_local(41,573) = - reac_rate_local(573) 
  reac_source_local(24,574) = + reac_rate_local(574) * 2.d0
  reac_source_local(41,574) = - reac_rate_local(574) 
  reac_source_local(05,575) = - reac_rate_local(575) 
  reac_source_local(41,575) = - reac_rate_local(575) 
  reac_source_local(42,575) = + reac_rate_local(575) 
  reac_source_local(43,576) = + reac_rate_local(576) * 2.d0
  reac_source_local(44,576) = - reac_rate_local(576) 
  reac_source_local(43,577) = - reac_rate_local(577) * 2.d0
  reac_source_local(44,577) = + reac_rate_local(577) 
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
  if( lgas_heating ) ZDPlasKin_cfg(1) = y(48)
  density(:) = y(1:species_max)
  call ZDPlasKin_reac_rates(t)
  rrt(001) = rrt(001) * density(01) * density(09) 
  rrt(002) = rrt(002) * density(01) * density(09) 
  rrt(003) = rrt(003) * density(01) * density(26) 
  rrt(004) = rrt(004) * density(01) * density(26) 
  rrt(005) = rrt(005) * density(01) * density(20) 
  rrt(006) = rrt(006) * density(01) * density(20) 
  rrt(007) = rrt(007) * density(01) * density(13) 
  rrt(008) = rrt(008) * density(01) * density(13) 
  rrt(009) = rrt(009) * density(01) * density(13) 
  rrt(010) = rrt(010) * density(01) * density(36) 
  rrt(011) = rrt(011) * density(01) * density(36) 
  rrt(012) = rrt(012) * density(01) * density(31) 
  rrt(013) = rrt(013) * density(01) * density(09) 
  rrt(014) = rrt(014) * density(01) * density(11) 
  rrt(015) = rrt(015) * density(01) * density(10) 
  rrt(016) = rrt(016) * density(01) * density(09) 
  rrt(017) = rrt(017) * density(01) * density(11) 
  rrt(018) = rrt(018) * density(01) * density(10) 
  rrt(019) = rrt(019) * density(01) * density(09) 
  rrt(020) = rrt(020) * density(01) * density(11) 
  rrt(021) = rrt(021) * density(01) * density(10) 
  rrt(022) = rrt(022) * density(01) * density(09) 
  rrt(023) = rrt(023) * density(01) * density(11) 
  rrt(024) = rrt(024) * density(01) * density(10) 
  rrt(025) = rrt(025) * density(01) * density(07) 
  rrt(026) = rrt(026) * density(01) * density(07) 
  rrt(027) = rrt(027) * density(01) * density(07) 
  rrt(028) = rrt(028) * density(01) * density(05) 
  rrt(029) = rrt(029) * density(01) * density(05) 
  rrt(030) = rrt(030) * density(01) * density(05) 
  rrt(031) = rrt(031) * density(01) * density(03) 
  rrt(032) = rrt(032) * density(01) * density(09) 
  rrt(033) = rrt(033) * density(01) * density(11) 
  rrt(034) = rrt(034) * density(01) * density(10) 
  rrt(035) = rrt(035) * density(01) * density(09) 
  rrt(036) = rrt(036) * density(01) * density(11) 
  rrt(037) = rrt(037) * density(01) * density(10) 
  rrt(038) = rrt(038) * density(01) * density(09) 
  rrt(039) = rrt(039) * density(01) * density(11) 
  rrt(040) = rrt(040) * density(01) * density(10) 
  rrt(041) = rrt(041) * density(01) * density(09) 
  rrt(042) = rrt(042) * density(01) * density(11) 
  rrt(043) = rrt(043) * density(01) * density(10) 
  rrt(044) = rrt(044) * density(01) * density(07) 
  rrt(045) = rrt(045) * density(01) * density(07) 
  rrt(046) = rrt(046) * density(01) * density(07) 
  rrt(047) = rrt(047) * density(01) * density(05) 
  rrt(048) = rrt(048) * density(01) * density(05) 
  rrt(049) = rrt(049) * density(01) * density(03) 
  rrt(050) = rrt(050) * density(01) * density(26) 
  rrt(051) = rrt(051) * density(01) * density(28) 
  rrt(052) = rrt(052) * density(01) * density(27) 
  rrt(053) = rrt(053) * density(01) * density(26) 
  rrt(054) = rrt(054) * density(01) * density(28) 
  rrt(055) = rrt(055) * density(01) * density(27) 
  rrt(056) = rrt(056) * density(01) * density(26) 
  rrt(057) = rrt(057) * density(01) * density(28) 
  rrt(058) = rrt(058) * density(01) * density(27) 
  rrt(059) = rrt(059) * density(01) * density(26) 
  rrt(060) = rrt(060) * density(01) * density(28) 
  rrt(061) = rrt(061) * density(01) * density(27) 
  rrt(062) = rrt(062) * density(01) * density(26) 
  rrt(063) = rrt(063) * density(01) * density(28) 
  rrt(064) = rrt(064) * density(01) * density(27) 
  rrt(065) = rrt(065) * density(01) * density(26) 
  rrt(066) = rrt(066) * density(01) * density(28) 
  rrt(067) = rrt(067) * density(01) * density(27) 
  rrt(068) = rrt(068) * density(01) * density(24) 
  rrt(069) = rrt(069) * density(01) * density(24) 
  rrt(070) = rrt(070) * density(01) * density(24) 
  rrt(071) = rrt(071) * density(01) * density(24) 
  rrt(072) = rrt(072) * density(01) * density(24) 
  rrt(073) = rrt(073) * density(01) * density(24) 
  rrt(074) = rrt(074) * density(01) * density(20) 
  rrt(075) = rrt(075) * density(01) * density(21) 
  rrt(076) = rrt(076) * density(01) * density(22) 
  rrt(077) = rrt(077) * density(01) * density(20) 
  rrt(078) = rrt(078) * density(01) * density(21) 
  rrt(079) = rrt(079) * density(01) * density(22) 
  rrt(080) = rrt(080) * density(01) * density(20) 
  rrt(081) = rrt(081) * density(01) * density(21) 
  rrt(082) = rrt(082) * density(01) * density(22) 
  rrt(083) = rrt(083) * density(01) * density(20) 
  rrt(084) = rrt(084) * density(01) * density(21) 
  rrt(085) = rrt(085) * density(01) * density(22) 
  rrt(086) = rrt(086) * density(01) * density(20) 
  rrt(087) = rrt(087) * density(01) * density(21) 
  rrt(088) = rrt(088) * density(01) * density(22) 
  rrt(089) = rrt(089) * density(01) * density(18) 
  rrt(090) = rrt(090) * density(01) * density(18) 
  rrt(091) = rrt(091) * density(01) * density(13) 
  rrt(092) = rrt(092) * density(01) * density(16) 
  rrt(093) = rrt(093) * density(01) * density(15) 
  rrt(094) = rrt(094) * density(01) * density(14) 
  rrt(095) = rrt(095) * density(01) * density(26) 
  rrt(096) = rrt(096) * density(01) * density(28) 
  rrt(097) = rrt(097) * density(01) * density(27) 
  rrt(098) = rrt(098) * density(01) * density(26) 
  rrt(099) = rrt(099) * density(01) * density(28) 
  rrt(100) = rrt(100) * density(01) * density(27) 
  rrt(101) = rrt(101) * density(01) * density(26) 
  rrt(102) = rrt(102) * density(01) * density(28) 
  rrt(103) = rrt(103) * density(01) * density(27) 
  rrt(104) = rrt(104) * density(01) * density(26) 
  rrt(105) = rrt(105) * density(01) * density(28) 
  rrt(106) = rrt(106) * density(01) * density(27) 
  rrt(107) = rrt(107) * density(01) * density(26) 
  rrt(108) = rrt(108) * density(01) * density(28) 
  rrt(109) = rrt(109) * density(01) * density(27) 
  rrt(110) = rrt(110) * density(01) * density(26) 
  rrt(111) = rrt(111) * density(01) * density(28) 
  rrt(112) = rrt(112) * density(01) * density(27) 
  rrt(113) = rrt(113) * density(01) * density(26) 
  rrt(114) = rrt(114) * density(01) * density(28) 
  rrt(115) = rrt(115) * density(01) * density(27) 
  rrt(116) = rrt(116) * density(01) * density(24) 
  rrt(117) = rrt(117) * density(01) * density(24) 
  rrt(118) = rrt(118) * density(01) * density(24) 
  rrt(119) = rrt(119) * density(01) * density(24) 
  rrt(120) = rrt(120) * density(01) * density(24) 
  rrt(121) = rrt(121) * density(01) * density(24) 
  rrt(122) = rrt(122) * density(01) * density(24) 
  rrt(123) = rrt(123) * density(01) * density(20) 
  rrt(124) = rrt(124) * density(01) * density(21) 
  rrt(125) = rrt(125) * density(01) * density(22) 
  rrt(126) = rrt(126) * density(01) * density(20) 
  rrt(127) = rrt(127) * density(01) * density(21) 
  rrt(128) = rrt(128) * density(01) * density(22) 
  rrt(129) = rrt(129) * density(01) * density(20) 
  rrt(130) = rrt(130) * density(01) * density(21) 
  rrt(131) = rrt(131) * density(01) * density(22) 
  rrt(132) = rrt(132) * density(01) * density(20) 
  rrt(133) = rrt(133) * density(01) * density(21) 
  rrt(134) = rrt(134) * density(01) * density(22) 
  rrt(135) = rrt(135) * density(01) * density(20) 
  rrt(136) = rrt(136) * density(01) * density(21) 
  rrt(137) = rrt(137) * density(01) * density(22) 
  rrt(138) = rrt(138) * density(01) * density(18) 
  rrt(139) = rrt(139) * density(01) * density(18) 
  rrt(140) = rrt(140) * density(01) * density(18) 
  rrt(141) = rrt(141) * density(01) * density(18) 
  rrt(142) = rrt(142) * density(01) * density(18) 
  rrt(143) = rrt(143) * density(01) * density(13) 
  rrt(144) = rrt(144) * density(01) * density(16) 
  rrt(145) = rrt(145) * density(01) * density(15) 
  rrt(146) = rrt(146) * density(01) * density(14) 
  rrt(147) = rrt(147) * density(01) * density(13) 
  rrt(148) = rrt(148) * density(01) * density(16) 
  rrt(149) = rrt(149) * density(01) * density(15) 
  rrt(150) = rrt(150) * density(01) * density(14) 
  rrt(151) = rrt(151) * density(01) * density(36) 
  rrt(152) = rrt(152) * density(01) * density(37) 
  rrt(153) = rrt(153) * density(01) * density(38) 
  rrt(154) = rrt(154) * density(01) * density(36) 
  rrt(155) = rrt(155) * density(01) * density(37) 
  rrt(156) = rrt(156) * density(01) * density(38) 
  rrt(157) = rrt(157) * density(01) * density(34) 
  rrt(158) = rrt(158) * density(01) * density(36) 
  rrt(159) = rrt(159) * density(01) * density(37) 
  rrt(160) = rrt(160) * density(01) * density(38) 
  rrt(161) = rrt(161) * density(01) * density(36) 
  rrt(162) = rrt(162) * density(01) * density(37) 
  rrt(163) = rrt(163) * density(01) * density(38) 
  rrt(164) = rrt(164) * density(01) * density(36) 
  rrt(165) = rrt(165) * density(01) * density(37) 
  rrt(166) = rrt(166) * density(01) * density(38) 
  rrt(167) = rrt(167) * density(01) * density(34) 
  rrt(168) = rrt(168) * density(01) * density(34) 
  rrt(169) = rrt(169) * density(01) * density(31) 
  rrt(170) = rrt(170) * density(01) * density(32) 
  rrt(171) = rrt(171) * density(01) * density(31) 
  rrt(172) = rrt(172) * density(01) * density(32) 
  rrt(173) = rrt(173) * density(01) * density(31) 
  rrt(174) = rrt(174) * density(01) * density(32) 
  rrt(175) = rrt(175) * density(01)**2 * density(04) 
  rrt(176) = rrt(176) * density(01)**2 * density(06) 
  rrt(177) = rrt(177) * density(01)**2 * density(08) 
  rrt(178) = rrt(178) * density(01)**2 * density(12) 
  rrt(179) = rrt(179) * density(01)**2 * density(17) 
  rrt(180) = rrt(180) * density(01)**2 * density(19) 
  rrt(181) = rrt(181) * density(01)**2 * density(23) 
  rrt(182) = rrt(182) * density(01)**2 * density(25) 
  rrt(183) = rrt(183) * density(01)**2 * density(29) 
  rrt(184) = rrt(184) * density(01)**2 * density(33) 
  rrt(185) = rrt(185) * density(01)**2 * density(35) 
  rrt(186) = rrt(186) * density(01)**2 * density(39) 
  rrt(187) = rrt(187) * density(01) * density(12) 
  rrt(188) = rrt(188) * density(01) * density(12) 
  rrt(189) = rrt(189) * density(01) * density(12) 
  rrt(190) = rrt(190) * density(01) * density(08) 
  rrt(191) = rrt(191) * density(01) * density(08) 
  rrt(192) = rrt(192) * density(01) * density(08) 
  rrt(193) = rrt(193) * density(01) * density(08) 
  rrt(194) = rrt(194) * density(01) * density(06) 
  rrt(195) = rrt(195) * density(01) * density(06) 
  rrt(196) = rrt(196) * density(01) * density(06) 
  rrt(197) = rrt(197) * density(01) * density(04) 
  rrt(198) = rrt(198) * density(01) * density(29) 
  rrt(199) = rrt(199) * density(01) * density(29) 
  rrt(200) = rrt(200) * density(01) * density(25) 
  rrt(201) = rrt(201) * density(01) * density(25) 
  rrt(202) = rrt(202) * density(01) * density(25) 
  rrt(203) = rrt(203) * density(01) * density(25) 
  rrt(204) = rrt(204) * density(01) * density(25) 
  rrt(205) = rrt(205) * density(01) * density(23) 
  rrt(206) = rrt(206) * density(01) * density(23) 
  rrt(207) = rrt(207) * density(01) * density(19) 
  rrt(208) = rrt(208) * density(01) * density(17) 
  rrt(209) = rrt(209) * density(01) * density(45) 
  rrt(210) = rrt(210) * density(01) * density(46) 
  rrt(211) = rrt(211) * density(44) * density(45) 
  rrt(212) = rrt(212) * density(12) * density(26) 
  rrt(213) = rrt(213) * density(12) * density(28) 
  rrt(214) = rrt(214) * density(12) * density(27) 
  rrt(215) = rrt(215) * density(12) * density(20) 
  rrt(216) = rrt(216) * density(12) * density(21) 
  rrt(217) = rrt(217) * density(12) * density(22) 
  rrt(218) = rrt(218) * density(12) * density(20) 
  rrt(219) = rrt(219) * density(12) * density(21) 
  rrt(220) = rrt(220) * density(12) * density(22) 
  rrt(221) = rrt(221) * density(12) * density(13) 
  rrt(222) = rrt(222) * density(12) * density(16) 
  rrt(223) = rrt(223) * density(12) * density(15) 
  rrt(224) = rrt(224) * density(12) * density(14) 
  rrt(225) = rrt(225) * density(12) * density(13) 
  rrt(226) = rrt(226) * density(12) * density(16) 
  rrt(227) = rrt(227) * density(12) * density(15) 
  rrt(228) = rrt(228) * density(12) * density(14) 
  rrt(229) = rrt(229) * density(12) * density(43) 
  rrt(230) = rrt(230) * density(08) * density(09) 
  rrt(231) = rrt(231) * density(08) * density(11) 
  rrt(232) = rrt(232) * density(08) * density(10) 
  rrt(233) = rrt(233) * density(08) * density(09) 
  rrt(234) = rrt(234) * density(08) * density(11) 
  rrt(235) = rrt(235) * density(08) * density(10) 
  rrt(236) = rrt(236) * density(05) * density(08) 
  rrt(237) = rrt(237) * density(03) * density(08) 
  rrt(238) = rrt(238) * density(08) * density(26) 
  rrt(239) = rrt(239) * density(08) * density(28) 
  rrt(240) = rrt(240) * density(08) * density(27) 
  rrt(241) = rrt(241) * density(08) * density(20) 
  rrt(242) = rrt(242) * density(08) * density(21) 
  rrt(243) = rrt(243) * density(08) * density(22) 
  rrt(244) = rrt(244) * density(08) * density(18) 
  rrt(245) = rrt(245) * density(06) * density(09) 
  rrt(246) = rrt(246) * density(06) * density(11) 
  rrt(247) = rrt(247) * density(06) * density(10) 
  rrt(248) = rrt(248) * density(06) * density(09) 
  rrt(249) = rrt(249) * density(06) * density(11) 
  rrt(250) = rrt(250) * density(06) * density(10) 
  rrt(251) = rrt(251) * density(06) * density(09) 
  rrt(252) = rrt(252) * density(06) * density(11) 
  rrt(253) = rrt(253) * density(06) * density(10) 
  rrt(254) = rrt(254) * density(06) * density(09) 
  rrt(255) = rrt(255) * density(06) * density(11) 
  rrt(256) = rrt(256) * density(06) * density(10) 
  rrt(257) = rrt(257) * density(06) * density(09) 
  rrt(258) = rrt(258) * density(06) * density(11) 
  rrt(259) = rrt(259) * density(06) * density(10) 
  rrt(260) = rrt(260) * density(06) * density(44) 
  rrt(261) = rrt(261) * density(04) * density(09) 
  rrt(262) = rrt(262) * density(04) * density(11) 
  rrt(263) = rrt(263) * density(04) * density(10) 
  rrt(264) = rrt(264) * density(04) * density(09) 
  rrt(265) = rrt(265) * density(04) * density(11) 
  rrt(266) = rrt(266) * density(04) * density(10) 
  rrt(267) = rrt(267) * density(04) * density(09) 
  rrt(268) = rrt(268) * density(04) * density(11) 
  rrt(269) = rrt(269) * density(04) * density(10) 
  rrt(270) = rrt(270) * density(04) * density(44) 
  rrt(271) = rrt(271) * density(20) * density(29) 
  rrt(272) = rrt(272) * density(21) * density(29) 
  rrt(273) = rrt(273) * density(22) * density(29) 
  rrt(274) = rrt(274) * density(13) * density(29) 
  rrt(275) = rrt(275) * density(16) * density(29) 
  rrt(276) = rrt(276) * density(15) * density(29) 
  rrt(277) = rrt(277) * density(14) * density(29) 
  rrt(278) = rrt(278) * density(29) * density(43) 
  rrt(279) = rrt(279) * density(25) * density(43) 
  rrt(280) = rrt(280) * density(18) * density(23) 
  rrt(281) = rrt(281) * density(18) * density(23) 
  rrt(282) = rrt(282) * density(23) * density(43) 
  rrt(283) = rrt(283) * density(19) * density(26) 
  rrt(284) = rrt(284) * density(19) * density(28) 
  rrt(285) = rrt(285) * density(19) * density(27) 
  rrt(286) = rrt(286) * density(19) * density(20) 
  rrt(287) = rrt(287) * density(19) * density(21) 
  rrt(288) = rrt(288) * density(19) * density(22) 
  rrt(289) = rrt(289) * density(19) * density(43) 
  rrt(290) = rrt(290) * density(09) * density(17) 
  rrt(291) = rrt(291) * density(11) * density(17) 
  rrt(292) = rrt(292) * density(10) * density(17) 
  rrt(293) = rrt(293) * density(17) * density(26) 
  rrt(294) = rrt(294) * density(17) * density(28) 
  rrt(295) = rrt(295) * density(17) * density(27) 
  rrt(296) = rrt(296) * density(17) * density(26) 
  rrt(297) = rrt(297) * density(17) * density(28) 
  rrt(298) = rrt(298) * density(17) * density(27) 
  rrt(299) = rrt(299) * density(17) * density(20) 
  rrt(300) = rrt(300) * density(17) * density(21) 
  rrt(301) = rrt(301) * density(17) * density(22) 
  rrt(302) = rrt(302) * density(17) * density(18) 
  rrt(303) = rrt(303) * density(17) * density(44) 
  rrt(304) = rrt(304) * density(07) * density(46) 
  rrt(305) = rrt(305) * density(05) * density(46) 
  rrt(306) = rrt(306) * density(03) * density(46) 
  rrt(307) = rrt(307) * density(26) * density(46) 
  rrt(308) = rrt(308) * density(28) * density(46) 
  rrt(309) = rrt(309) * density(27) * density(46) 
  rrt(310) = rrt(310) * density(24) * density(46) 
  rrt(311) = rrt(311) * density(20) * density(46) 
  rrt(312) = rrt(312) * density(21) * density(46) 
  rrt(313) = rrt(313) * density(22) * density(46) 
  rrt(314) = rrt(314) * density(20) * density(46) 
  rrt(315) = rrt(315) * density(21) * density(46) 
  rrt(316) = rrt(316) * density(22) * density(46) 
  rrt(317) = rrt(317) * density(18) * density(46) 
  rrt(318) = rrt(318) * density(13) * density(46) 
  rrt(319) = rrt(319) * density(16) * density(46) 
  rrt(320) = rrt(320) * density(15) * density(46) 
  rrt(321) = rrt(321) * density(14) * density(46) 
  rrt(322) = rrt(322) * density(09) * density(45) 
  rrt(323) = rrt(323) * density(11) * density(45) 
  rrt(324) = rrt(324) * density(10) * density(45) 
  rrt(325) = rrt(325) * density(09) * density(45) 
  rrt(326) = rrt(326) * density(11) * density(45) 
  rrt(327) = rrt(327) * density(10) * density(45) 
  rrt(328) = rrt(328) * density(05) * density(45) 
  rrt(329) = rrt(329) * density(05) * density(45) 
  rrt(330) = rrt(330) * density(03) * density(45) 
  rrt(331) = rrt(331) * density(03) * density(45) 
  rrt(332) = rrt(332) * density(26) * density(45) 
  rrt(333) = rrt(333) * density(28) * density(45) 
  rrt(334) = rrt(334) * density(27) * density(45) 
  rrt(335) = rrt(335) * density(26) * density(45) 
  rrt(336) = rrt(336) * density(28) * density(45) 
  rrt(337) = rrt(337) * density(27) * density(45) 
  rrt(338) = rrt(338) * density(26) * density(45) 
  rrt(339) = rrt(339) * density(28) * density(45) 
  rrt(340) = rrt(340) * density(27) * density(45) 
  rrt(341) = rrt(341) * density(26) * density(45) 
  rrt(342) = rrt(342) * density(28) * density(45) 
  rrt(343) = rrt(343) * density(27) * density(45) 
  rrt(344) = rrt(344) * density(26) * density(45) 
  rrt(345) = rrt(345) * density(28) * density(45) 
  rrt(346) = rrt(346) * density(27) * density(45) 
  rrt(347) = rrt(347) * density(20) * density(45) 
  rrt(348) = rrt(348) * density(21) * density(45) 
  rrt(349) = rrt(349) * density(22) * density(45) 
  rrt(350) = rrt(350) * density(20) * density(45) 
  rrt(351) = rrt(351) * density(21) * density(45) 
  rrt(352) = rrt(352) * density(22) * density(45) 
  rrt(353) = rrt(353) * density(20) * density(45) 
  rrt(354) = rrt(354) * density(21) * density(45) 
  rrt(355) = rrt(355) * density(22) * density(45) 
  rrt(356) = rrt(356) * density(13) * density(45) 
  rrt(357) = rrt(357) * density(16) * density(45) 
  rrt(358) = rrt(358) * density(15) * density(45) 
  rrt(359) = rrt(359) * density(14) * density(45) 
  rrt(360) = rrt(360) * density(13) * density(45) 
  rrt(361) = rrt(361) * density(16) * density(45) 
  rrt(362) = rrt(362) * density(15) * density(45) 
  rrt(363) = rrt(363) * density(14) * density(45) 
  rrt(364) = rrt(364) * density(43) * density(45) 
  rrt(365) = rrt(365) * density(43) * density(45) 
  rrt(366) = rrt(366) * density(09) * density(47) 
  rrt(367) = rrt(367) * density(11) * density(47) 
  rrt(368) = rrt(368) * density(10) * density(47) 
  rrt(369) = rrt(369) * density(09) * density(47) 
  rrt(370) = rrt(370) * density(11) * density(47) 
  rrt(371) = rrt(371) * density(10) * density(47) 
  rrt(372) = rrt(372) * density(07) * density(47) 
  rrt(373) = rrt(373) * density(05) * density(47) 
  rrt(374) = rrt(374) * density(05) * density(47) 
  rrt(375) = rrt(375) * density(03) * density(47) 
  rrt(376) = rrt(376) * density(26) * density(47) 
  rrt(377) = rrt(377) * density(28) * density(47) 
  rrt(378) = rrt(378) * density(27) * density(47) 
  rrt(379) = rrt(379) * density(26) * density(47) 
  rrt(380) = rrt(380) * density(28) * density(47) 
  rrt(381) = rrt(381) * density(27) * density(47) 
  rrt(382) = rrt(382) * density(26) * density(47) 
  rrt(383) = rrt(383) * density(28) * density(47) 
  rrt(384) = rrt(384) * density(27) * density(47) 
  rrt(385) = rrt(385) * density(24) * density(47) 
  rrt(386) = rrt(386) * density(24) * density(47) 
  rrt(387) = rrt(387) * density(20) * density(47) 
  rrt(388) = rrt(388) * density(21) * density(47) 
  rrt(389) = rrt(389) * density(22) * density(47) 
  rrt(390) = rrt(390) * density(20) * density(47) 
  rrt(391) = rrt(391) * density(21) * density(47) 
  rrt(392) = rrt(392) * density(22) * density(47) 
  rrt(393) = rrt(393) * density(20) * density(47) 
  rrt(394) = rrt(394) * density(18) * density(47) 
  rrt(395) = rrt(395) * density(18) * density(47) 
  rrt(396) = rrt(396) * density(13) * density(47) 
  rrt(397) = rrt(397) * density(16) * density(47) 
  rrt(398) = rrt(398) * density(15) * density(47) 
  rrt(399) = rrt(399) * density(14) * density(47) 
  rrt(400) = rrt(400) * density(05) * density(09) 
  rrt(401) = rrt(401) * density(05) * density(11) 
  rrt(402) = rrt(402) * density(05) * density(10) 
  rrt(403) = rrt(403) * density(03) * density(09) 
  rrt(404) = rrt(404) * density(03) * density(11) 
  rrt(405) = rrt(405) * density(03) * density(10) 
  rrt(406) = rrt(406) * density(09) * density(24) 
  rrt(407) = rrt(407) * density(11) * density(24) 
  rrt(408) = rrt(408) * density(10) * density(24) 
  rrt(409) = rrt(409) * density(09) * density(18) 
  rrt(410) = rrt(410) * density(11) * density(18) 
  rrt(411) = rrt(411) * density(10) * density(18) 
  rrt(412) = rrt(412) * density(09) * density(34) 
  rrt(413) = rrt(413) * density(11) * density(34) 
  rrt(414) = rrt(414) * density(10) * density(34) 
  rrt(415) = rrt(415) * density(09) * density(43) 
  rrt(416) = rrt(416) * density(11) * density(43) 
  rrt(417) = rrt(417) * density(10) * density(43) 
  rrt(418) = rrt(418) * density(07) * density(09) 
  rrt(419) = rrt(419) * density(07) * density(11) 
  rrt(420) = rrt(420) * density(07) * density(10) 
  rrt(421) = rrt(421) * density(05) * density(09) 
  rrt(422) = rrt(422) * density(05) * density(11) 
  rrt(423) = rrt(423) * density(05) * density(10) 
  rrt(424) = rrt(424) * density(09) 
  rrt(425) = rrt(425) * density(11) 
  rrt(426) = rrt(426) * density(10) 
  rrt(427) = rrt(427) * density(07) 
  rrt(428) = rrt(428) * density(07) 
  rrt(429) = rrt(429) * density(07) * density(24) 
  rrt(430) = rrt(430) * density(05)**2 
  rrt(431) = rrt(431) * density(05) * density(24) 
  rrt(432) = rrt(432) * density(05) * density(18) 
  rrt(433) = rrt(433) * density(05) * density(36) 
  rrt(434) = rrt(434) * density(05) * density(37) 
  rrt(435) = rrt(435) * density(05) * density(38) 
  rrt(436) = rrt(436) * density(05) * density(34) 
  rrt(437) = rrt(437) * density(05) * density(34) 
  rrt(438) = rrt(438) * density(05) * density(44) 
  rrt(439) = rrt(439) * density(05) * density(43) 
  rrt(440) = rrt(440) * density(05) 
  rrt(441) = rrt(441) * density(05) * density(43) 
  rrt(442) = rrt(442) * density(03) * density(26) 
  rrt(443) = rrt(443) * density(03) * density(28) 
  rrt(444) = rrt(444) * density(03) * density(27) 
  rrt(445) = rrt(445) * density(03) * density(26) 
  rrt(446) = rrt(446) * density(03) * density(28) 
  rrt(447) = rrt(447) * density(03) * density(27) 
  rrt(448) = rrt(448) * density(03) * density(44) 
  rrt(449) = rrt(449) * density(03) * density(07) 
  rrt(450) = rrt(450) * density(03) * density(05) 
  rrt(451) = rrt(451) * density(03) * density(44) 
  rrt(452) = rrt(452) * density(03) * density(18) 
  rrt(453) = rrt(453) * density(18) * density(26) 
  rrt(454) = rrt(454) * density(18) * density(28) 
  rrt(455) = rrt(455) * density(18) * density(27) 
  rrt(456) = rrt(456) * density(26) * density(34) 
  rrt(457) = rrt(457) * density(28) * density(34) 
  rrt(458) = rrt(458) * density(27) * density(34) 
  rrt(459) = rrt(459) * density(26) * density(43) 
  rrt(460) = rrt(460) * density(28) * density(43) 
  rrt(461) = rrt(461) * density(27) * density(43) 
  rrt(462) = rrt(462) * density(26) * density(43) 
  rrt(463) = rrt(463) * density(28) * density(43) 
  rrt(464) = rrt(464) * density(27) * density(43) 
  rrt(465) = rrt(465) * density(26) 
  rrt(466) = rrt(466) * density(28) 
  rrt(467) = rrt(467) * density(27) 
  rrt(468) = rrt(468) * density(03) * density(26) 
  rrt(469) = rrt(469) * density(03) * density(28) 
  rrt(470) = rrt(470) * density(03) * density(27) 
  rrt(471) = rrt(471) * density(05) * density(26) 
  rrt(472) = rrt(472) * density(05) * density(28) 
  rrt(473) = rrt(473) * density(05) * density(27) 
  rrt(474) = rrt(474) * density(24)**2 
  rrt(475) = rrt(475) * density(20) * density(24) 
  rrt(476) = rrt(476) * density(21) * density(24) 
  rrt(477) = rrt(477) * density(22) * density(24) 
  rrt(478) = rrt(478) * density(24) * density(36) 
  rrt(479) = rrt(479) * density(24) * density(37) 
  rrt(480) = rrt(480) * density(24) * density(38) 
  rrt(481) = rrt(481) * density(24) * density(34) 
  rrt(482) = rrt(482) * density(24) * density(34) 
  rrt(483) = rrt(483) * density(24) * density(31) 
  rrt(484) = rrt(484) * density(24) * density(32) 
  rrt(485) = rrt(485) * density(24) * density(44) 
  rrt(486) = rrt(486) * density(24) * density(43) 
  rrt(487) = rrt(487) * density(24) * density(43) 
  rrt(488) = rrt(488) * density(24) * density(43) 
  rrt(489) = rrt(489) * density(24) 
  rrt(490) = rrt(490) * density(24)**2 
  rrt(491) = rrt(491) * density(18) * density(24) 
  rrt(492) = rrt(492) * density(20) * density(43) 
  rrt(493) = rrt(493) * density(21) * density(43) 
  rrt(494) = rrt(494) * density(22) * density(43) 
  rrt(495) = rrt(495) * density(20) * density(43) 
  rrt(496) = rrt(496) * density(21) * density(43) 
  rrt(497) = rrt(497) * density(22) * density(43) 
  rrt(498) = rrt(498) * density(20) * density(44) 
  rrt(499) = rrt(499) * density(20) 
  rrt(500) = rrt(500) * density(21) 
  rrt(501) = rrt(501) * density(22) 
  rrt(502) = rrt(502) * density(13) * density(20) 
  rrt(503) = rrt(503) * density(13) * density(21) 
  rrt(504) = rrt(504) * density(13) * density(22) 
  rrt(505) = rrt(505) * density(16) * density(20) 
  rrt(506) = rrt(506) * density(15) * density(20) 
  rrt(507) = rrt(507) * density(14) * density(20) 
  rrt(508) = rrt(508) * density(20) * density(31) 
  rrt(509) = rrt(509) * density(21) * density(31) 
  rrt(510) = rrt(510) * density(22) * density(31) 
  rrt(511) = rrt(511) * density(20) * density(32) 
  rrt(512) = rrt(512) * density(20)**2 
  rrt(513) = rrt(513) * density(21)**2 
  rrt(514) = rrt(514) * density(22)**2 
  rrt(515) = rrt(515) * density(07) * density(20) 
  rrt(516) = rrt(516) * density(07) * density(21) 
  rrt(517) = rrt(517) * density(07) * density(22) 
  rrt(518) = rrt(518) * density(20) 
  rrt(519) = rrt(519) * density(21) 
  rrt(520) = rrt(520) * density(22) 
  rrt(521) = rrt(521) * density(20) * density(44) 
  rrt(522) = rrt(522) * density(21) * density(44) 
  rrt(523) = rrt(523) * density(22) * density(44) 
  rrt(524) = rrt(524) * density(05) * density(20) 
  rrt(525) = rrt(525) * density(05) * density(21) 
  rrt(526) = rrt(526) * density(05) * density(22) 
  rrt(527) = rrt(527) * density(18)**2 
  rrt(528) = rrt(528) * density(18) * density(36) 
  rrt(529) = rrt(529) * density(18) * density(37) 
  rrt(530) = rrt(530) * density(18) * density(38) 
  rrt(531) = rrt(531) * density(18) * density(34) 
  rrt(532) = rrt(532) * density(18) * density(34) 
  rrt(533) = rrt(533) * density(18) * density(44) 
  rrt(534) = rrt(534) * density(18) * density(43) 
  rrt(535) = rrt(535) * density(18) * density(43) 
  rrt(536) = rrt(536) * density(18) 
  rrt(537) = rrt(537) * density(13) * density(43) 
  rrt(538) = rrt(538) * density(16) * density(43) 
  rrt(539) = rrt(539) * density(15) * density(43) 
  rrt(540) = rrt(540) * density(14) * density(43) 
  rrt(541) = rrt(541) * density(13) * density(44) 
  rrt(542) = rrt(542) * density(16) * density(44) 
  rrt(543) = rrt(543) * density(15) * density(44) 
  rrt(544) = rrt(544) * density(14) * density(44) 
  rrt(545) = rrt(545) * density(13) * density(44) 
  rrt(546) = rrt(546) * density(16) * density(44) 
  rrt(547) = rrt(547) * density(15) * density(44) 
  rrt(548) = rrt(548) * density(14) * density(44) 
  rrt(549) = rrt(549) * density(07) * density(13) 
  rrt(550) = rrt(550) * density(07) * density(16) 
  rrt(551) = rrt(551) * density(07) * density(15) 
  rrt(552) = rrt(552) * density(07) * density(14) 
  rrt(553) = rrt(553) * density(36) * density(43) 
  rrt(554) = rrt(554) * density(37) * density(43) 
  rrt(555) = rrt(555) * density(38) * density(43) 
  rrt(556) = rrt(556) * density(36) 
  rrt(557) = rrt(557) * density(37) 
  rrt(558) = rrt(558) * density(38) 
  rrt(559) = rrt(559) * density(05) * density(36) 
  rrt(560) = rrt(560) * density(05) * density(37) 
  rrt(561) = rrt(561) * density(05) * density(38) 
  rrt(562) = rrt(562) * density(34)**2 
  rrt(563) = rrt(563) * density(34) * density(44) 
  rrt(564) = rrt(564) * density(34) * density(43) 
  rrt(565) = rrt(565) * density(34) * density(43) 
  rrt(566) = rrt(566) * density(34) * density(43) 
  rrt(567) = rrt(567) * density(34) 
  rrt(568) = rrt(568) * density(34) 
  rrt(569) = rrt(569) * density(31) * density(43) 
  rrt(570) = rrt(570) * density(32) * density(43) 
  rrt(571) = rrt(571) * density(31) 
  rrt(572) = rrt(572) * density(07) * density(41) 
  rrt(573) = rrt(573) * density(41) 
  rrt(574) = rrt(574) * density(41) 
  rrt(575) = rrt(575) * density(05) * density(41) 
  rrt(576) = rrt(576) * density(44) 
  rrt(577) = rrt(577) * density(43)**2 
  ydot(01) = +rrt(032)+rrt(033)+rrt(034)+rrt(035)+rrt(036)+rrt(037)+rrt(038)+rrt(039)+rrt(040)+rrt(041)+rrt(042)+rrt(043)+rrt(044)&
             +rrt(045)+rrt(046)+rrt(047)+rrt(048)+rrt(049)+rrt(095)+rrt(096)+rrt(097)+rrt(098)+rrt(099)+rrt(100)+rrt(101)+rrt(102)&
             +rrt(103)+rrt(104)+rrt(105)+rrt(106)+rrt(107)+rrt(108)+rrt(109)+rrt(110)+rrt(111)+rrt(112)+rrt(113)+rrt(114)+rrt(115)&
             +rrt(116)+rrt(117)+rrt(118)+rrt(119)+rrt(120)+rrt(121)+rrt(122)+rrt(123)+rrt(124)+rrt(125)+rrt(126)+rrt(127)+rrt(128)&
             +rrt(129)+rrt(130)+rrt(131)+rrt(132)+rrt(133)+rrt(134)+rrt(135)+rrt(136)+rrt(137)+rrt(138)+rrt(139)+rrt(140)+rrt(141)&
             +rrt(142)+rrt(143)+rrt(144)+rrt(145)+rrt(146)+rrt(147)+rrt(148)+rrt(149)+rrt(150)+rrt(158)+rrt(159)+rrt(160)+rrt(161)&
             +rrt(162)+rrt(163)+rrt(164)+rrt(165)+rrt(166)+rrt(167)+rrt(168)+rrt(169)+rrt(170)+rrt(173)+rrt(174)-rrt(175)-rrt(176)&
             -rrt(177)-rrt(178)-rrt(179)-rrt(180)-rrt(181)-rrt(182)-rrt(183)-rrt(184)-rrt(185)-rrt(186)-rrt(187)-rrt(188)-rrt(189)&
             -rrt(190)-rrt(191)-rrt(192)-rrt(193)-rrt(194)-rrt(195)-rrt(196)-rrt(197)-rrt(198)-rrt(199)-rrt(200)-rrt(201)-rrt(202)&
             -rrt(203)-rrt(204)-rrt(205)-rrt(206)-rrt(207)-rrt(208)-rrt(209)-rrt(210) 
  ydot(02) = +rrt(022)+rrt(023)+rrt(024)+rrt(027)+rrt(029)+rrt(030)+rrt(031)+rrt(193)+rrt(195)+rrt(196)+rrt(197) 
  ydot(03) = +rrt(019)+rrt(020)+rrt(021)+rrt(026)+rrt(028)-rrt(031)-rrt(049)+rrt(072)+rrt(083)+rrt(084)+rrt(085)+rrt(090)&
             +  2.d0 * rrt(091)+  2.d0 * rrt(092)+  2.d0 * rrt(093)+  2.d0 * rrt(094)+rrt(129)+rrt(130)+rrt(131)+rrt(140)+rrt(147)&
             +rrt(148)+rrt(149)+rrt(150)+rrt(175)+rrt(189)+rrt(191)+rrt(192)+rrt(194)+  2.d0 * rrt(208)-rrt(237)-rrt(306)-rrt(330)&
             -rrt(331)-rrt(375)-rrt(403)-rrt(404)-rrt(405)+rrt(428)+rrt(439)+rrt(440)-rrt(442)-rrt(443)-rrt(444)-rrt(445)-rrt(446)&
             -rrt(447)-rrt(448)-rrt(449)-rrt(450)-rrt(451)-rrt(452)-rrt(468)-rrt(469)-rrt(470) 
  ydot(04) = +rrt(041)+rrt(042)+rrt(043)+rrt(046)+rrt(048)+rrt(049)+rrt(122)+rrt(135)+rrt(136)+rrt(137)+rrt(141)+rrt(147)+rrt(148)&
             +rrt(149)+rrt(150)-rrt(175)-rrt(197)-rrt(261)-rrt(262)-rrt(263)-rrt(264)-rrt(265)-rrt(266)-rrt(267)-rrt(268)-rrt(269)&
             -rrt(270)+rrt(331)+rrt(374)+rrt(375) 
  ydot(05) = +rrt(016)+rrt(017)+rrt(018)+rrt(025)-rrt(028)-rrt(029)-rrt(030)-rrt(047)-rrt(048)+rrt(062)+rrt(063)+rrt(064)+rrt(073)&
             +  2.d0 * rrt(086)+  2.d0 * rrt(087)+  2.d0 * rrt(088)+rrt(090)+rrt(120)+rrt(132)+rrt(133)+rrt(134)+rrt(141)+rrt(176)&
             +rrt(188)+rrt(190)+rrt(204)-rrt(236)-rrt(305)-rrt(328)-rrt(329)-rrt(373)-rrt(374)-rrt(400)-rrt(401)-rrt(402)-rrt(421)&
             -rrt(422)-rrt(423)+rrt(427)+rrt(429)-  2.d0 * rrt(430)-rrt(431)-rrt(432)-rrt(433)-rrt(434)-rrt(435)-rrt(436)-rrt(437)&
             -rrt(438)-rrt(439)-rrt(440)-rrt(441)+rrt(448)-rrt(450)+rrt(452)-rrt(471)-rrt(472)-rrt(473)-rrt(524)-rrt(525)-rrt(526)&
             -rrt(559)-rrt(560)-rrt(561)-rrt(575) 
  ydot(06) = +rrt(038)+rrt(039)+rrt(040)+rrt(045)+rrt(047)+rrt(113)+rrt(114)+rrt(115)+rrt(121)+rrt(132)+rrt(133)+rrt(134)+rrt(140)&
             -rrt(176)-rrt(194)-rrt(195)-rrt(196)-rrt(245)-rrt(246)-rrt(247)-rrt(248)-rrt(249)-rrt(250)-rrt(251)-rrt(252)-rrt(253)&
             -rrt(254)-rrt(255)-rrt(256)-rrt(257)-rrt(258)-rrt(259)-rrt(260)+rrt(270)+rrt(306)+rrt(329)+rrt(330)+rrt(373) 
  ydot(07) = +rrt(013)+rrt(014)+rrt(015)-rrt(025)-rrt(026)-rrt(027)-rrt(044)-rrt(045)-rrt(046)+  2.d0 * rrt(065)+  2.d0 * rrt(066)&
             +  2.d0 * rrt(067)+rrt(073)+rrt(083)+rrt(084)+rrt(085)+rrt(110)+rrt(111)+rrt(112)+rrt(121)+rrt(135)+rrt(136)+rrt(137)&
             +rrt(177)+rrt(187)+rrt(204)+rrt(215)+rrt(216)+rrt(217)+rrt(221)+rrt(222)+rrt(223)+rrt(224)+rrt(230)+rrt(231)+rrt(232)&
             +rrt(244)+rrt(245)+rrt(246)+rrt(247)+rrt(290)+rrt(291)+rrt(292)-rrt(304)-rrt(372)+  2.d0 * rrt(400)+  2.d0 * rrt(401)&
             +  2.d0 * rrt(402)+rrt(406)+rrt(407)+rrt(408)+rrt(409)+rrt(410)+rrt(411)+rrt(412)+rrt(413)+rrt(414)+rrt(415)+rrt(416)&
             +rrt(417)-rrt(418)-rrt(419)-rrt(420)+rrt(424)+rrt(425)+rrt(426)-rrt(427)-rrt(428)-rrt(429)+rrt(431)+rrt(432)+rrt(433)&
             +rrt(434)+rrt(435)+rrt(437)+rrt(438)+rrt(441)-rrt(449)+rrt(451)+rrt(462)+rrt(463)+rrt(464)+  2.d0 * rrt(465)&
             +  2.d0 * rrt(466)+  2.d0 * rrt(467)+rrt(468)+rrt(469)+rrt(470)+rrt(471)+rrt(472)+rrt(473)+  2.d0 * rrt(486)-rrt(515)&
             -rrt(516)-rrt(517)-rrt(549)-rrt(550)-rrt(551)-rrt(552)+rrt(556)+rrt(557)+rrt(558)+rrt(566)+rrt(568)+rrt(571)-rrt(572)&
             +rrt(573) 
  ydot(08) = +rrt(035)+rrt(036)+rrt(037)+rrt(044)+rrt(110)+rrt(111)+rrt(112)+rrt(120)+rrt(129)+rrt(130)+rrt(131)-rrt(177)-rrt(190)&
             -rrt(191)-rrt(192)-rrt(193)+rrt(229)-rrt(230)-rrt(231)-rrt(232)-rrt(233)-rrt(234)-rrt(235)-rrt(236)-rrt(237)-rrt(238)&
             -rrt(239)-rrt(240)-rrt(241)-rrt(242)-rrt(243)-rrt(244)+rrt(245)+rrt(246)+rrt(247)+rrt(260)+rrt(305)+rrt(325)+rrt(326)&
             +rrt(327)+rrt(328)+rrt(369)+rrt(370)+rrt(371)+rrt(372) 
  ydot(09) = -rrt(001)-rrt(002)-rrt(013)-rrt(016)-rrt(019)-rrt(022)-rrt(032)-rrt(035)-rrt(038)-rrt(041)+rrt(062)+rrt(063)+rrt(064)&
             +rrt(072)+rrt(113)+rrt(114)+rrt(115)+rrt(122)+rrt(171)+rrt(172)+rrt(178)+rrt(212)+rrt(213)+rrt(214)+rrt(218)+rrt(219)&
             +rrt(220)+rrt(225)+rrt(226)+rrt(227)+rrt(228)-rrt(230)-rrt(233)+rrt(238)+rrt(239)+rrt(240)+rrt(241)+rrt(242)+rrt(243)&
             -rrt(245)-rrt(248)-rrt(251)-rrt(254)-rrt(257)-rrt(261)-rrt(264)-rrt(267)-rrt(290)-rrt(322)-rrt(325)-rrt(366)-rrt(369)&
             -rrt(400)-rrt(403)-rrt(406)-rrt(409)-rrt(412)-rrt(415)-rrt(418)-rrt(421)-rrt(424)+rrt(462)+rrt(463)+rrt(464)+rrt(572) 
  ydot(10) = +rrt(002)-rrt(015)-rrt(018)-rrt(021)-rrt(024)-rrt(034)-rrt(037)-rrt(040)-rrt(043)-rrt(232)-rrt(235)-rrt(247)-rrt(250)&
             -rrt(253)-rrt(256)-rrt(259)-rrt(263)-rrt(266)-rrt(269)-rrt(292)-rrt(324)-rrt(327)-rrt(368)-rrt(371)-rrt(402)-rrt(405)&
             -rrt(408)-rrt(411)-rrt(414)-rrt(417)-rrt(420)-rrt(423)-rrt(426) 
  ydot(11) = +rrt(001)-rrt(014)-rrt(017)-rrt(020)-rrt(023)-rrt(033)-rrt(036)-rrt(039)-rrt(042)-rrt(231)-rrt(234)-rrt(246)-rrt(249)&
             -rrt(252)-rrt(255)-rrt(258)-rrt(262)-rrt(265)-rrt(268)-rrt(291)-rrt(323)-rrt(326)-rrt(367)-rrt(370)-rrt(401)-rrt(404)&
             -rrt(407)-rrt(410)-rrt(413)-rrt(416)-rrt(419)-rrt(422)-rrt(425) 
  ydot(12) = +rrt(032)+rrt(033)+rrt(034)+rrt(173)+rrt(174)-rrt(178)-rrt(187)-rrt(188)-rrt(189)-rrt(212)-rrt(213)-rrt(214)-rrt(215)&
             -rrt(216)-rrt(217)-rrt(218)-rrt(219)-rrt(220)-rrt(221)-rrt(222)-rrt(223)-rrt(224)-rrt(225)-rrt(226)-rrt(227)-rrt(228)&
             -rrt(229)+rrt(230)+rrt(231)+rrt(232)+rrt(304)+rrt(322)+rrt(323)+rrt(324)+rrt(366)+rrt(367)+rrt(368) 
  ydot(13) = -rrt(007)-rrt(008)-rrt(009)+rrt(059)+rrt(060)+rrt(061)+rrt(071)+rrt(077)+rrt(078)+rrt(079)+rrt(080)+rrt(081)+rrt(082)&
             +rrt(089)-rrt(091)+rrt(142)-rrt(143)-rrt(147)+rrt(171)+rrt(172)+rrt(173)+rrt(174)+rrt(179)+rrt(202)+rrt(203)+rrt(206)&
             +rrt(207)-rrt(221)-rrt(225)-rrt(274)+rrt(280)+rrt(286)+rrt(287)+rrt(288)+rrt(299)+rrt(300)+rrt(301)+rrt(302)-rrt(318)&
             -rrt(356)-rrt(360)-rrt(396)+rrt(430)+rrt(432)+rrt(450)+rrt(452)-rrt(502)-rrt(503)-rrt(504)+rrt(518)+rrt(519)+rrt(520)&
             +rrt(527)+rrt(531)+rrt(534)+rrt(536)-rrt(537)-rrt(541)-rrt(545)-rrt(549) 
  ydot(14) = +rrt(009)-rrt(094)-rrt(146)-rrt(150)-rrt(224)-rrt(228)-rrt(277)-rrt(321)-rrt(359)-rrt(363)-rrt(399)-rrt(507)-rrt(540)&
             -rrt(544)-rrt(548)-rrt(552) 
  ydot(15) = +rrt(008)-rrt(093)-rrt(145)-rrt(149)-rrt(223)-rrt(227)-rrt(276)-rrt(320)-rrt(358)-rrt(362)-rrt(398)-rrt(506)-rrt(539)&
             -rrt(543)-rrt(547)-rrt(551) 
  ydot(16) = +rrt(007)-rrt(092)-rrt(144)-rrt(148)-rrt(222)-rrt(226)-rrt(275)-rrt(319)-rrt(357)-rrt(361)-rrt(397)-rrt(505)-rrt(538)&
             -rrt(542)-rrt(546)-rrt(550) 
  ydot(17) = +rrt(107)+rrt(108)+rrt(109)+rrt(119)+rrt(139)+rrt(143)+rrt(144)+rrt(145)+rrt(146)-rrt(179)-rrt(208)+rrt(225)+rrt(226)&
             +rrt(227)+rrt(228)+rrt(237)+rrt(257)+rrt(258)+rrt(259)+rrt(267)+rrt(268)+rrt(269)+rrt(289)-rrt(290)-rrt(291)-rrt(292)&
             -rrt(293)-rrt(294)-rrt(295)-rrt(296)-rrt(297)-rrt(298)-rrt(299)-rrt(300)-rrt(301)-rrt(302)-rrt(303)+rrt(344)+rrt(345)&
             +rrt(346)+rrt(353)+rrt(354)+rrt(355)+rrt(360)+rrt(361)+rrt(362)+rrt(363)+rrt(393)+rrt(395)+rrt(396)+rrt(397)+rrt(398)&
             +rrt(399) 
  ydot(18) = +rrt(056)+rrt(057)+rrt(058)+rrt(069)+rrt(070)+rrt(074)+rrt(075)+rrt(076)-rrt(089)-rrt(090)-rrt(138)-rrt(139)-rrt(140)&
             -rrt(141)-rrt(142)+rrt(180)+rrt(201)+rrt(205)-rrt(244)+rrt(274)+rrt(275)+rrt(276)+rrt(277)-rrt(280)-rrt(281)+rrt(293)&
             +rrt(294)+rrt(295)-rrt(302)-rrt(317)-rrt(394)-rrt(395)-rrt(409)-rrt(410)-rrt(411)-rrt(432)+rrt(449)-rrt(452)-rrt(453)&
             -rrt(454)-rrt(455)+rrt(475)+rrt(476)+rrt(477)-rrt(491)+rrt(492)+rrt(493)+rrt(494)+rrt(499)+rrt(500)+rrt(501)&
             +  2.d0 * rrt(502)+  2.d0 * rrt(503)+  2.d0 * rrt(504)+  2.d0 * rrt(505)+  2.d0 * rrt(506)+  2.d0 * rrt(507)+rrt(508)&
             +rrt(509)+rrt(510)+rrt(511)+rrt(512)+rrt(513)+rrt(514)-  2.d0 * rrt(527)-rrt(528)-rrt(529)-rrt(530)-rrt(531)-rrt(532)&
             -rrt(533)-rrt(534)-rrt(535)-rrt(536)+rrt(537)+rrt(538)+rrt(539)+rrt(540)+rrt(545)+rrt(546)+rrt(547)+rrt(548)+rrt(571) 
  ydot(19) = +rrt(104)+rrt(105)+rrt(106)+rrt(118)+rrt(126)+rrt(127)+rrt(128)+rrt(138)-rrt(180)-rrt(207)+rrt(221)+rrt(222)+rrt(223)&
             +rrt(224)+rrt(236)+rrt(241)+rrt(242)+rrt(243)+rrt(244)+rrt(254)+rrt(255)+rrt(256)+rrt(264)+rrt(265)+rrt(266)+rrt(281)&
             +rrt(282)-rrt(283)-rrt(284)-rrt(285)-rrt(286)-rrt(287)-rrt(288)-rrt(289)+rrt(290)+rrt(291)+rrt(292)+rrt(302)+rrt(303)&
             +rrt(314)+rrt(315)+rrt(316)+rrt(318)+rrt(319)+rrt(320)+rrt(321)+rrt(341)+rrt(342)+rrt(343)+rrt(350)+rrt(351)+rrt(352)&
             +rrt(356)+rrt(357)+rrt(358)+rrt(359)+rrt(382)+rrt(383)+rrt(384)+rrt(386)+rrt(390)+rrt(391)+rrt(392)+rrt(394) 
  ydot(20) = -rrt(005)-rrt(006)+rrt(053)+rrt(054)+rrt(055)+rrt(068)-rrt(074)-rrt(077)-rrt(080)-rrt(083)-rrt(086)-rrt(123)-rrt(126)&
             -rrt(129)-rrt(132)-rrt(135)+rrt(181)+rrt(199)+rrt(200)-rrt(215)-rrt(218)-rrt(241)-rrt(271)+rrt(281)+rrt(283)+rrt(284)&
             +rrt(285)-rrt(286)+rrt(296)+rrt(297)+rrt(298)-rrt(299)-rrt(311)-rrt(314)-rrt(347)-rrt(350)-rrt(353)-rrt(387)-rrt(390)&
             -rrt(393)+rrt(403)+rrt(404)+rrt(405)+rrt(409)+rrt(410)+rrt(411)+rrt(431)+rrt(436)+rrt(453)+rrt(454)+rrt(455)+rrt(468)&
             +rrt(469)+rrt(470)+rrt(474)-rrt(475)+rrt(481)+rrt(487)+rrt(489)+  2.d0 * rrt(491)-rrt(492)-rrt(495)-rrt(498)-rrt(499)&
             -rrt(502)-rrt(505)-rrt(506)-rrt(507)-rrt(508)-rrt(511)-  2.d0 * rrt(512)-rrt(515)-rrt(518)-rrt(521)-rrt(524)+rrt(527)&
             +rrt(528)+rrt(529)+rrt(530)+rrt(532)+rrt(533)+rrt(535)+rrt(541)+rrt(542)+rrt(543)+rrt(544)+rrt(568) 
  ydot(21) = +rrt(005)-rrt(075)-rrt(078)-rrt(081)-rrt(084)-rrt(087)-rrt(124)-rrt(127)-rrt(130)-rrt(133)-rrt(136)-rrt(216)-rrt(219)&
             -rrt(242)-rrt(272)-rrt(287)-rrt(300)-rrt(312)-rrt(315)-rrt(348)-rrt(351)-rrt(354)-rrt(388)-rrt(391)-rrt(476)-rrt(493)&
             -rrt(496)-rrt(500)-rrt(503)-rrt(509)-  2.d0 * rrt(513)-rrt(516)-rrt(519)-rrt(522)-rrt(525) 
  ydot(22) = +rrt(006)-rrt(076)-rrt(079)-rrt(082)-rrt(085)-rrt(088)-rrt(125)-rrt(128)-rrt(131)-rrt(134)-rrt(137)-rrt(217)-rrt(220)&
             -rrt(243)-rrt(273)-rrt(288)-rrt(301)-rrt(313)-rrt(316)-rrt(349)-rrt(352)-rrt(355)-rrt(389)-rrt(392)-rrt(477)-rrt(494)&
             -rrt(497)-rrt(501)-rrt(504)-rrt(510)-  2.d0 * rrt(514)-rrt(517)-rrt(520)-rrt(523)-rrt(526) 
  ydot(23) = +rrt(101)+rrt(102)+rrt(103)+rrt(117)+rrt(123)+rrt(124)+rrt(125)-rrt(181)-rrt(205)-rrt(206)+rrt(212)+rrt(213)+rrt(214)&
             +rrt(218)+rrt(219)+rrt(220)+rrt(251)+rrt(252)+rrt(253)+rrt(261)+rrt(262)+rrt(263)+rrt(271)+rrt(272)+rrt(273)+rrt(279)&
             -rrt(280)-rrt(281)-rrt(282)+rrt(296)+rrt(297)+rrt(298)+rrt(299)+rrt(300)+rrt(301)+rrt(317)+rrt(338)+rrt(339)+rrt(340)&
             +rrt(347)+rrt(348)+rrt(349)+rrt(379)+rrt(380)+rrt(381)+rrt(385)+rrt(387)+rrt(388)+rrt(389) 
  ydot(24) = +rrt(050)+rrt(051)+rrt(052)-rrt(068)-rrt(069)-rrt(070)-rrt(071)-rrt(072)-rrt(073)-rrt(116)-rrt(117)-rrt(118)-rrt(119)&
             -rrt(120)-rrt(121)-rrt(122)+rrt(182)+rrt(198)-rrt(310)-rrt(385)-rrt(386)-rrt(406)-rrt(407)-rrt(408)-rrt(429)-rrt(431)&
             +rrt(436)+rrt(453)+rrt(454)+rrt(455)+rrt(456)+rrt(457)+rrt(458)+rrt(459)+rrt(460)+rrt(461)+rrt(471)+rrt(472)+rrt(473)&
             -  2.d0 * rrt(474)-rrt(475)-rrt(476)-rrt(477)-rrt(478)-rrt(479)-rrt(480)-rrt(481)-rrt(482)-rrt(483)-rrt(484)-rrt(485)&
             -rrt(486)-rrt(487)-rrt(488)-rrt(489)-  2.d0 * rrt(490)-rrt(491)+rrt(495)+rrt(496)+rrt(497)+rrt(498)+rrt(512)+rrt(513)&
             +rrt(514)+rrt(556)+rrt(557)+rrt(558)+rrt(566)+  2.d0 * rrt(574) 
  ydot(25) = +rrt(098)+rrt(099)+rrt(100)+rrt(116)-rrt(182)-rrt(200)-rrt(201)-rrt(202)-rrt(203)-rrt(204)+rrt(215)+rrt(216)+rrt(217)&
             +rrt(233)+rrt(234)+rrt(235)+rrt(238)+rrt(239)+rrt(240)+rrt(248)+rrt(249)+rrt(250)+rrt(274)+rrt(275)+rrt(276)+rrt(277)&
             +rrt(278)-rrt(279)+rrt(280)+rrt(283)+rrt(284)+rrt(285)+rrt(286)+rrt(287)+rrt(288)+rrt(293)+rrt(294)+rrt(295)+rrt(307)&
             +rrt(308)+rrt(309)+rrt(311)+rrt(312)+rrt(313)+rrt(335)+rrt(336)+rrt(337)+rrt(376)+rrt(377)+rrt(378) 
  ydot(26) = -rrt(003)-rrt(004)-rrt(050)-rrt(053)-rrt(056)-rrt(059)-rrt(062)-rrt(065)-rrt(095)-rrt(098)-rrt(101)-rrt(104)-rrt(107)&
             -rrt(110)-rrt(113)+rrt(183)-rrt(212)-rrt(238)+rrt(271)+rrt(272)+rrt(273)-rrt(283)-rrt(293)-rrt(296)-rrt(307)-rrt(332)&
             -rrt(335)-rrt(338)-rrt(341)-rrt(344)-rrt(376)-rrt(379)-rrt(382)+rrt(406)+rrt(407)+rrt(408)+rrt(418)+rrt(419)+rrt(420)&
             +rrt(421)+rrt(422)+rrt(423)+rrt(429)-rrt(442)-rrt(445)-rrt(453)-rrt(456)-rrt(459)-rrt(462)-rrt(465)-rrt(468)-rrt(471)&
             +rrt(474)+rrt(475)+rrt(476)+rrt(477)+rrt(478)+rrt(479)+rrt(480)+rrt(482)+rrt(483)+rrt(484)+rrt(485)+rrt(488)+rrt(521)&
             +rrt(522)+rrt(523) 
  ydot(27) = +rrt(004)-rrt(052)-rrt(055)-rrt(058)-rrt(061)-rrt(064)-rrt(067)-rrt(097)-rrt(100)-rrt(103)-rrt(106)-rrt(109)-rrt(112)&
             -rrt(115)-rrt(214)-rrt(240)-rrt(285)-rrt(295)-rrt(298)-rrt(309)-rrt(334)-rrt(337)-rrt(340)-rrt(343)-rrt(346)-rrt(378)&
             -rrt(381)-rrt(384)-rrt(444)-rrt(447)-rrt(455)-rrt(458)-rrt(461)-rrt(464)-rrt(467)-rrt(470)-rrt(473) 
  ydot(28) = +rrt(003)-rrt(051)-rrt(054)-rrt(057)-rrt(060)-rrt(063)-rrt(066)-rrt(096)-rrt(099)-rrt(102)-rrt(105)-rrt(108)-rrt(111)&
             -rrt(114)-rrt(213)-rrt(239)-rrt(284)-rrt(294)-rrt(297)-rrt(308)-rrt(333)-rrt(336)-rrt(339)-rrt(342)-rrt(345)-rrt(377)&
             -rrt(380)-rrt(383)-rrt(443)-rrt(446)-rrt(454)-rrt(457)-rrt(460)-rrt(463)-rrt(466)-rrt(469)-rrt(472) 
  ydot(29) = +rrt(095)+rrt(096)+rrt(097)-rrt(183)-rrt(198)-rrt(199)-rrt(271)-rrt(272)-rrt(273)-rrt(274)-rrt(275)-rrt(276)-rrt(277)&
             -rrt(278)+rrt(310)+rrt(332)+rrt(333)+rrt(334) 
  ydot(30) = +rrt(483)+rrt(484)+rrt(549)+rrt(550)+rrt(551)+rrt(552) 
  ydot(31) = -rrt(012)+rrt(154)+rrt(155)+rrt(156)+rrt(157)-rrt(169)-rrt(171)-rrt(173)+rrt(184)+rrt(437)+rrt(442)+rrt(443)+rrt(444)&
             +rrt(482)-rrt(483)-rrt(508)-rrt(509)-rrt(510)+rrt(524)+rrt(525)+rrt(526)+rrt(532)+rrt(562)+rrt(564)+rrt(567)-rrt(569)&
             -rrt(571) 
  ydot(32) = +rrt(012)-rrt(170)-rrt(172)-rrt(174)-rrt(484)-rrt(511)-rrt(570) 
  ydot(33) = +rrt(164)+rrt(165)+rrt(166)+rrt(168)+rrt(169)+rrt(170)-rrt(184) 
  ydot(34) = +rrt(151)+rrt(152)+rrt(153)-rrt(157)-rrt(167)-rrt(168)+rrt(185)-rrt(412)-rrt(413)-rrt(414)+rrt(433)+rrt(434)+rrt(435)&
             -rrt(436)-rrt(437)+rrt(445)+rrt(446)+rrt(447)-rrt(456)-rrt(457)-rrt(458)+rrt(478)+rrt(479)+rrt(480)-rrt(481)-rrt(482)&
             +rrt(508)+rrt(509)+rrt(510)+rrt(511)+rrt(515)+rrt(516)+rrt(517)+rrt(528)+rrt(529)+rrt(530)-rrt(531)-rrt(532)+rrt(553)&
             +rrt(554)+rrt(555)-  2.d0 * rrt(562)-rrt(563)-rrt(564)-rrt(565)-rrt(566)-rrt(567)-rrt(568)+rrt(569)+rrt(570)+rrt(573) 
  ydot(35) = +rrt(161)+rrt(162)+rrt(163)+rrt(167)-rrt(185) 
  ydot(36) = -rrt(010)-rrt(011)-rrt(151)-rrt(154)-rrt(158)-rrt(161)-rrt(164)+rrt(186)+rrt(412)+rrt(413)+rrt(414)-rrt(433)+rrt(456)&
             +rrt(457)+rrt(458)-rrt(478)+rrt(481)-rrt(528)+rrt(531)-rrt(553)-rrt(556)-rrt(559)+rrt(562)+rrt(563)+rrt(565) 
  ydot(37) = +rrt(010)-rrt(152)-rrt(155)-rrt(159)-rrt(162)-rrt(165)-rrt(434)-rrt(479)-rrt(529)-rrt(554)-rrt(557)-rrt(560) 
  ydot(38) = +rrt(011)-rrt(153)-rrt(156)-rrt(160)-rrt(163)-rrt(166)-rrt(435)-rrt(480)-rrt(530)-rrt(555)-rrt(558)-rrt(561) 
  ydot(39) = +rrt(158)+rrt(159)+rrt(160)-rrt(186) 
  ydot(40) = +rrt(572) 
  ydot(41) = +rrt(490)+rrt(559)+rrt(560)+rrt(561)-rrt(572)-rrt(573)-rrt(574)-rrt(575) 
  ydot(42) = +rrt(575) 
  ydot(43) = +rrt(013)+rrt(014)+rrt(015)+rrt(019)+rrt(020)+rrt(021)+rrt(025)+rrt(027)+rrt(028)+  2.d0 * rrt(030)+rrt(031)+rrt(035)&
             +rrt(036)+rrt(037)+rrt(041)+rrt(042)+rrt(043)+rrt(045)+rrt(048)+rrt(050)+rrt(051)+rrt(052)+rrt(056)+rrt(057)+rrt(058)&
             +rrt(068)+  2.d0 * rrt(070)+rrt(071)+rrt(074)+rrt(075)+rrt(076)+  2.d0 * rrt(080)+  2.d0 * rrt(081)+  2.d0 * rrt(082)&
             +rrt(089)+rrt(098)+rrt(099)+rrt(100)+rrt(104)+rrt(105)+rrt(106)+rrt(117)+rrt(119)+rrt(126)+rrt(127)+rrt(128)+rrt(139)&
             +rrt(151)+rrt(152)+rrt(153)+rrt(157)+rrt(161)+rrt(162)+rrt(163)+rrt(168)+rrt(187)+  2.d0 * rrt(188)+rrt(189)+rrt(190)&
             +  2.d0 * rrt(192)+rrt(193)+rrt(194)+  2.d0 * rrt(196)+rrt(197)+rrt(198)+  2.d0 * rrt(199)+rrt(200)+  2.d0 * rrt(201)&
             +rrt(202)+  3.d0 * rrt(203)+rrt(205)+  2.d0 * rrt(206)+rrt(207)+  2.d0 * rrt(209)+rrt(210)+rrt(211)-rrt(229)+rrt(248)&
             +rrt(249)+rrt(250)+rrt(254)+rrt(255)+rrt(256)+rrt(260)+rrt(261)+rrt(262)+rrt(263)+rrt(267)+rrt(268)+rrt(269)+rrt(270)&
             -rrt(278)-rrt(279)-rrt(282)-rrt(289)+rrt(303)+rrt(325)+rrt(326)+rrt(327)+rrt(328)+rrt(330)+rrt(335)+rrt(336)+rrt(337)&
             +rrt(341)+rrt(342)+rrt(343)+rrt(350)+rrt(351)+rrt(352)+rrt(356)+rrt(357)+rrt(358)+rrt(359)-rrt(364)-rrt(365)+rrt(366)&
             +rrt(367)+rrt(368)+rrt(372)+rrt(373)+rrt(375)+rrt(379)+rrt(380)+rrt(381)+rrt(386)+rrt(387)+rrt(388)+rrt(389)+rrt(393)&
             +rrt(394)+rrt(396)+rrt(397)+rrt(398)+rrt(399)+rrt(403)+rrt(404)+rrt(405)-rrt(415)-rrt(416)-rrt(417)+rrt(418)+rrt(419)&
             +rrt(420)+rrt(424)+rrt(425)+rrt(426)+rrt(427)+rrt(438)-rrt(439)+rrt(440)-rrt(441)+rrt(442)+rrt(443)+rrt(444)+rrt(448)&
             +rrt(449)+rrt(450)-rrt(459)-rrt(460)-rrt(461)-rrt(462)-rrt(463)-rrt(464)+rrt(485)-rrt(486)-rrt(487)-rrt(488)+rrt(489)&
             -rrt(492)-rrt(493)-rrt(494)-rrt(495)-rrt(496)-rrt(497)+rrt(498)+rrt(499)+rrt(500)+rrt(501)+rrt(533)-rrt(534)-rrt(535)&
             +rrt(536)-rrt(537)-rrt(538)-rrt(539)-rrt(540)+rrt(545)+rrt(546)+rrt(547)+rrt(548)-rrt(553)-rrt(554)-rrt(555)+rrt(563)&
             -rrt(564)-rrt(565)-rrt(566)+rrt(567)-rrt(569)-rrt(570)+  2.d0 * rrt(576)-  2.d0 * rrt(577) 
  ydot(44) = +rrt(016)+rrt(017)+rrt(018)+rrt(019)+rrt(020)+rrt(021)+  2.d0 * rrt(022)+  2.d0 * rrt(023)+  2.d0 * rrt(024)+rrt(026)&
             +rrt(027)+rrt(029)+rrt(038)+rrt(039)+rrt(040)+rrt(041)+rrt(042)+rrt(043)+rrt(046)+rrt(053)+rrt(054)+rrt(055)+rrt(056)&
             +rrt(057)+rrt(058)+  2.d0 * rrt(059)+  2.d0 * rrt(060)+  2.d0 * rrt(061)+rrt(069)+rrt(071)+rrt(077)+rrt(078)+rrt(079)&
             +rrt(101)+rrt(102)+rrt(103)+rrt(104)+rrt(105)+rrt(106)+  2.d0 * rrt(107)+  2.d0 * rrt(108)+  2.d0 * rrt(109)+rrt(118)&
             +rrt(119)+rrt(154)+rrt(155)+rrt(156)+rrt(164)+rrt(165)+rrt(166)+rrt(189)+rrt(191)+rrt(193)+rrt(195)+rrt(202)+rrt(210)&
             -rrt(211)+rrt(212)+rrt(213)+rrt(214)+rrt(229)+rrt(233)+rrt(234)+rrt(235)+rrt(236)+rrt(237)+rrt(251)+rrt(252)+rrt(253)&
             +rrt(254)+rrt(255)+rrt(256)+  2.d0 * rrt(257)+  2.d0 * rrt(258)+  2.d0 * rrt(259)-rrt(260)+rrt(264)+rrt(265)+rrt(266)&
             +rrt(267)+rrt(268)+rrt(269)-rrt(270)+rrt(278)+rrt(279)+rrt(282)+rrt(289)-rrt(303)+rrt(304)+rrt(305)+rrt(306)&
             +  2.d0 * rrt(307)+  2.d0 * rrt(308)+  2.d0 * rrt(309)+rrt(310)+rrt(311)+rrt(312)+rrt(313)+  2.d0 * rrt(314)&
             +  2.d0 * rrt(315)+  2.d0 * rrt(316)+rrt(317)+rrt(318)+rrt(319)+rrt(320)+rrt(321)+rrt(322)+rrt(323)+rrt(324)+rrt(325)&
             +rrt(326)+rrt(327)+rrt(329)+rrt(331)+rrt(332)+rrt(333)+rrt(334)+rrt(335)+rrt(336)+rrt(337)+  2.d0 * rrt(338)&
             +  2.d0 * rrt(339)+  2.d0 * rrt(340)+  2.d0 * rrt(341)+  2.d0 * rrt(342)+  2.d0 * rrt(343)+  3.d0 * rrt(344)&
             +  3.d0 * rrt(345)+  3.d0 * rrt(346)+rrt(347)+rrt(348)+rrt(349)+rrt(350)+rrt(351)+rrt(352)+  2.d0 * rrt(353)&
             +  2.d0 * rrt(354)+  2.d0 * rrt(355)+rrt(360)+rrt(361)+rrt(362)+rrt(363)+rrt(365)+rrt(369)+rrt(370)+rrt(371)+rrt(374)&
             +rrt(376)+rrt(377)+rrt(378)+rrt(379)+rrt(380)+rrt(381)+  2.d0 * rrt(382)+  2.d0 * rrt(383)+  2.d0 * rrt(384)+rrt(385)&
             +rrt(386)+rrt(390)+rrt(391)+rrt(392)+rrt(393)+rrt(395)+rrt(415)+rrt(416)+rrt(417)+rrt(428)+rrt(430)-rrt(438)+rrt(439)&
             -rrt(448)-rrt(451)+rrt(459)+rrt(460)+rrt(461)-rrt(485)+rrt(487)+rrt(492)+rrt(493)+rrt(494)-rrt(498)+rrt(518)+rrt(519)&
             +rrt(520)-rrt(521)-rrt(522)-rrt(523)-rrt(533)+rrt(534)-rrt(541)-rrt(542)-rrt(543)-rrt(544)-rrt(545)-rrt(546)-rrt(547)&
             -rrt(548)+rrt(553)+rrt(554)+rrt(555)-rrt(563)+rrt(564)-rrt(576)+rrt(577) 
  ydot(45) = -rrt(209)-rrt(211)-rrt(322)-rrt(323)-rrt(324)-rrt(325)-rrt(326)-rrt(327)-rrt(328)-rrt(329)-rrt(330)-rrt(331)-rrt(332)&
             -rrt(333)-rrt(334)-rrt(335)-rrt(336)-rrt(337)-rrt(338)-rrt(339)-rrt(340)-rrt(341)-rrt(342)-rrt(343)-rrt(344)-rrt(345)&
             -rrt(346)-rrt(347)-rrt(348)-rrt(349)-rrt(350)-rrt(351)-rrt(352)-rrt(353)-rrt(354)-rrt(355)-rrt(356)-rrt(357)-rrt(358)&
             -rrt(359)-rrt(360)-rrt(361)-rrt(362)-rrt(363)-rrt(364)-rrt(365) 
  ydot(46) = -rrt(210)+rrt(211)-rrt(304)-rrt(305)-rrt(306)-rrt(307)-rrt(308)-rrt(309)-rrt(310)-rrt(311)-rrt(312)-rrt(313)-rrt(314)&
             -rrt(315)-rrt(316)-rrt(317)-rrt(318)-rrt(319)-rrt(320)-rrt(321)+rrt(364) 
  ydot(47) = +rrt(142)+rrt(365)-rrt(366)-rrt(367)-rrt(368)-rrt(369)-rrt(370)-rrt(371)-rrt(372)-rrt(373)-rrt(374)-rrt(375)-rrt(376)&
             -rrt(377)-rrt(378)-rrt(379)-rrt(380)-rrt(381)-rrt(382)-rrt(383)-rrt(384)-rrt(385)-rrt(386)-rrt(387)-rrt(388)-rrt(389)&
             -rrt(390)-rrt(391)-rrt(392)-rrt(393)-rrt(394)-rrt(395)-rrt(396)-rrt(397)-rrt(398)-rrt(399) 
  if( ldensity_constant ) where( density_constant(:) ) ydot(1:species_max) = 0.0d0
  ydot(48) = 0.0d0
  if( lgas_heating ) then
    ydot(48) = ( ZDPlasKin_cfg(14)/k_B + ydot(48) ) / ( sum(density(1:species_max)) - density(species_electrons) ) &
            + eV_to_K * ZDPlasKin_cfg(11) * density(species_electrons)
    ydot(48) = ydot(48) * ZDPlasKin_cfg(13)
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
  if( lgas_heating ) ZDPlasKin_cfg(1) = y(48)
  density(:) = y(1:species_max)
  call ZDPlasKin_reac_rates(t)
  pd(09,01) = pd(09,01) - rrt(001) * density(09) 
  pd(09,09) = pd(09,09) - rrt(001) * density(01) 
  pd(11,01) = pd(11,01) + rrt(001) * density(09) 
  pd(11,09) = pd(11,09) + rrt(001) * density(01) 
  pd(09,01) = pd(09,01) - rrt(002) * density(09) 
  pd(09,09) = pd(09,09) - rrt(002) * density(01) 
  pd(10,01) = pd(10,01) + rrt(002) * density(09) 
  pd(10,09) = pd(10,09) + rrt(002) * density(01) 
  pd(26,01) = pd(26,01) - rrt(003) * density(26) 
  pd(26,26) = pd(26,26) - rrt(003) * density(01) 
  pd(28,01) = pd(28,01) + rrt(003) * density(26) 
  pd(28,26) = pd(28,26) + rrt(003) * density(01) 
  pd(26,01) = pd(26,01) - rrt(004) * density(26) 
  pd(26,26) = pd(26,26) - rrt(004) * density(01) 
  pd(27,01) = pd(27,01) + rrt(004) * density(26) 
  pd(27,26) = pd(27,26) + rrt(004) * density(01) 
  pd(20,01) = pd(20,01) - rrt(005) * density(20) 
  pd(20,20) = pd(20,20) - rrt(005) * density(01) 
  pd(21,01) = pd(21,01) + rrt(005) * density(20) 
  pd(21,20) = pd(21,20) + rrt(005) * density(01) 
  pd(20,01) = pd(20,01) - rrt(006) * density(20) 
  pd(20,20) = pd(20,20) - rrt(006) * density(01) 
  pd(22,01) = pd(22,01) + rrt(006) * density(20) 
  pd(22,20) = pd(22,20) + rrt(006) * density(01) 
  pd(13,01) = pd(13,01) - rrt(007) * density(13) 
  pd(13,13) = pd(13,13) - rrt(007) * density(01) 
  pd(16,01) = pd(16,01) + rrt(007) * density(13) 
  pd(16,13) = pd(16,13) + rrt(007) * density(01) 
  pd(13,01) = pd(13,01) - rrt(008) * density(13) 
  pd(13,13) = pd(13,13) - rrt(008) * density(01) 
  pd(15,01) = pd(15,01) + rrt(008) * density(13) 
  pd(15,13) = pd(15,13) + rrt(008) * density(01) 
  pd(13,01) = pd(13,01) - rrt(009) * density(13) 
  pd(13,13) = pd(13,13) - rrt(009) * density(01) 
  pd(14,01) = pd(14,01) + rrt(009) * density(13) 
  pd(14,13) = pd(14,13) + rrt(009) * density(01) 
  pd(36,01) = pd(36,01) - rrt(010) * density(36) 
  pd(36,36) = pd(36,36) - rrt(010) * density(01) 
  pd(37,01) = pd(37,01) + rrt(010) * density(36) 
  pd(37,36) = pd(37,36) + rrt(010) * density(01) 
  pd(36,01) = pd(36,01) - rrt(011) * density(36) 
  pd(36,36) = pd(36,36) - rrt(011) * density(01) 
  pd(38,01) = pd(38,01) + rrt(011) * density(36) 
  pd(38,36) = pd(38,36) + rrt(011) * density(01) 
  pd(31,01) = pd(31,01) - rrt(012) * density(31) 
  pd(31,31) = pd(31,31) - rrt(012) * density(01) 
  pd(32,01) = pd(32,01) + rrt(012) * density(31) 
  pd(32,31) = pd(32,31) + rrt(012) * density(01) 
  pd(07,01) = pd(07,01) + rrt(013) * density(09) 
  pd(07,09) = pd(07,09) + rrt(013) * density(01) 
  pd(09,01) = pd(09,01) - rrt(013) * density(09) 
  pd(09,09) = pd(09,09) - rrt(013) * density(01) 
  pd(43,01) = pd(43,01) + rrt(013) * density(09) 
  pd(43,09) = pd(43,09) + rrt(013) * density(01) 
  pd(07,01) = pd(07,01) + rrt(014) * density(11) 
  pd(07,11) = pd(07,11) + rrt(014) * density(01) 
  pd(11,01) = pd(11,01) - rrt(014) * density(11) 
  pd(11,11) = pd(11,11) - rrt(014) * density(01) 
  pd(43,01) = pd(43,01) + rrt(014) * density(11) 
  pd(43,11) = pd(43,11) + rrt(014) * density(01) 
  pd(07,01) = pd(07,01) + rrt(015) * density(10) 
  pd(07,10) = pd(07,10) + rrt(015) * density(01) 
  pd(10,01) = pd(10,01) - rrt(015) * density(10) 
  pd(10,10) = pd(10,10) - rrt(015) * density(01) 
  pd(43,01) = pd(43,01) + rrt(015) * density(10) 
  pd(43,10) = pd(43,10) + rrt(015) * density(01) 
  pd(05,01) = pd(05,01) + rrt(016) * density(09) 
  pd(05,09) = pd(05,09) + rrt(016) * density(01) 
  pd(09,01) = pd(09,01) - rrt(016) * density(09) 
  pd(09,09) = pd(09,09) - rrt(016) * density(01) 
  pd(44,01) = pd(44,01) + rrt(016) * density(09) 
  pd(44,09) = pd(44,09) + rrt(016) * density(01) 
  pd(05,01) = pd(05,01) + rrt(017) * density(11) 
  pd(05,11) = pd(05,11) + rrt(017) * density(01) 
  pd(11,01) = pd(11,01) - rrt(017) * density(11) 
  pd(11,11) = pd(11,11) - rrt(017) * density(01) 
  pd(44,01) = pd(44,01) + rrt(017) * density(11) 
  pd(44,11) = pd(44,11) + rrt(017) * density(01) 
  pd(05,01) = pd(05,01) + rrt(018) * density(10) 
  pd(05,10) = pd(05,10) + rrt(018) * density(01) 
  pd(10,01) = pd(10,01) - rrt(018) * density(10) 
  pd(10,10) = pd(10,10) - rrt(018) * density(01) 
  pd(44,01) = pd(44,01) + rrt(018) * density(10) 
  pd(44,10) = pd(44,10) + rrt(018) * density(01) 
  pd(03,01) = pd(03,01) + rrt(019) * density(09) 
  pd(03,09) = pd(03,09) + rrt(019) * density(01) 
  pd(09,01) = pd(09,01) - rrt(019) * density(09) 
  pd(09,09) = pd(09,09) - rrt(019) * density(01) 
  pd(43,01) = pd(43,01) + rrt(019) * density(09) 
  pd(43,09) = pd(43,09) + rrt(019) * density(01) 
  pd(44,01) = pd(44,01) + rrt(019) * density(09) 
  pd(44,09) = pd(44,09) + rrt(019) * density(01) 
  pd(03,01) = pd(03,01) + rrt(020) * density(11) 
  pd(03,11) = pd(03,11) + rrt(020) * density(01) 
  pd(11,01) = pd(11,01) - rrt(020) * density(11) 
  pd(11,11) = pd(11,11) - rrt(020) * density(01) 
  pd(43,01) = pd(43,01) + rrt(020) * density(11) 
  pd(43,11) = pd(43,11) + rrt(020) * density(01) 
  pd(44,01) = pd(44,01) + rrt(020) * density(11) 
  pd(44,11) = pd(44,11) + rrt(020) * density(01) 
  pd(03,01) = pd(03,01) + rrt(021) * density(10) 
  pd(03,10) = pd(03,10) + rrt(021) * density(01) 
  pd(10,01) = pd(10,01) - rrt(021) * density(10) 
  pd(10,10) = pd(10,10) - rrt(021) * density(01) 
  pd(43,01) = pd(43,01) + rrt(021) * density(10) 
  pd(43,10) = pd(43,10) + rrt(021) * density(01) 
  pd(44,01) = pd(44,01) + rrt(021) * density(10) 
  pd(44,10) = pd(44,10) + rrt(021) * density(01) 
  pd(02,01) = pd(02,01) + rrt(022) * density(09) 
  pd(02,09) = pd(02,09) + rrt(022) * density(01) 
  pd(09,01) = pd(09,01) - rrt(022) * density(09) 
  pd(09,09) = pd(09,09) - rrt(022) * density(01) 
  pd(44,01) = pd(44,01) + rrt(022) * density(09) * 2.0d0
  pd(44,09) = pd(44,09) + rrt(022) * density(01) * 2.0d0
  pd(02,01) = pd(02,01) + rrt(023) * density(11) 
  pd(02,11) = pd(02,11) + rrt(023) * density(01) 
  pd(11,01) = pd(11,01) - rrt(023) * density(11) 
  pd(11,11) = pd(11,11) - rrt(023) * density(01) 
  pd(44,01) = pd(44,01) + rrt(023) * density(11) * 2.0d0
  pd(44,11) = pd(44,11) + rrt(023) * density(01) * 2.0d0
  pd(02,01) = pd(02,01) + rrt(024) * density(10) 
  pd(02,10) = pd(02,10) + rrt(024) * density(01) 
  pd(10,01) = pd(10,01) - rrt(024) * density(10) 
  pd(10,10) = pd(10,10) - rrt(024) * density(01) 
  pd(44,01) = pd(44,01) + rrt(024) * density(10) * 2.0d0
  pd(44,10) = pd(44,10) + rrt(024) * density(01) * 2.0d0
  pd(05,01) = pd(05,01) + rrt(025) * density(07) 
  pd(05,07) = pd(05,07) + rrt(025) * density(01) 
  pd(07,01) = pd(07,01) - rrt(025) * density(07) 
  pd(07,07) = pd(07,07) - rrt(025) * density(01) 
  pd(43,01) = pd(43,01) + rrt(025) * density(07) 
  pd(43,07) = pd(43,07) + rrt(025) * density(01) 
  pd(03,01) = pd(03,01) + rrt(026) * density(07) 
  pd(03,07) = pd(03,07) + rrt(026) * density(01) 
  pd(07,01) = pd(07,01) - rrt(026) * density(07) 
  pd(07,07) = pd(07,07) - rrt(026) * density(01) 
  pd(44,01) = pd(44,01) + rrt(026) * density(07) 
  pd(44,07) = pd(44,07) + rrt(026) * density(01) 
  pd(02,01) = pd(02,01) + rrt(027) * density(07) 
  pd(02,07) = pd(02,07) + rrt(027) * density(01) 
  pd(07,01) = pd(07,01) - rrt(027) * density(07) 
  pd(07,07) = pd(07,07) - rrt(027) * density(01) 
  pd(43,01) = pd(43,01) + rrt(027) * density(07) 
  pd(43,07) = pd(43,07) + rrt(027) * density(01) 
  pd(44,01) = pd(44,01) + rrt(027) * density(07) 
  pd(44,07) = pd(44,07) + rrt(027) * density(01) 
  pd(03,01) = pd(03,01) + rrt(028) * density(05) 
  pd(03,05) = pd(03,05) + rrt(028) * density(01) 
  pd(05,01) = pd(05,01) - rrt(028) * density(05) 
  pd(05,05) = pd(05,05) - rrt(028) * density(01) 
  pd(43,01) = pd(43,01) + rrt(028) * density(05) 
  pd(43,05) = pd(43,05) + rrt(028) * density(01) 
  pd(02,01) = pd(02,01) + rrt(029) * density(05) 
  pd(02,05) = pd(02,05) + rrt(029) * density(01) 
  pd(05,01) = pd(05,01) - rrt(029) * density(05) 
  pd(05,05) = pd(05,05) - rrt(029) * density(01) 
  pd(44,01) = pd(44,01) + rrt(029) * density(05) 
  pd(44,05) = pd(44,05) + rrt(029) * density(01) 
  pd(02,01) = pd(02,01) + rrt(030) * density(05) 
  pd(02,05) = pd(02,05) + rrt(030) * density(01) 
  pd(05,01) = pd(05,01) - rrt(030) * density(05) 
  pd(05,05) = pd(05,05) - rrt(030) * density(01) 
  pd(43,01) = pd(43,01) + rrt(030) * density(05) * 2.0d0
  pd(43,05) = pd(43,05) + rrt(030) * density(01) * 2.0d0
  pd(02,01) = pd(02,01) + rrt(031) * density(03) 
  pd(02,03) = pd(02,03) + rrt(031) * density(01) 
  pd(03,01) = pd(03,01) - rrt(031) * density(03) 
  pd(03,03) = pd(03,03) - rrt(031) * density(01) 
  pd(43,01) = pd(43,01) + rrt(031) * density(03) 
  pd(43,03) = pd(43,03) + rrt(031) * density(01) 
  pd(01,01) = pd(01,01) + rrt(032) * density(09) 
  pd(01,09) = pd(01,09) + rrt(032) * density(01) 
  pd(09,01) = pd(09,01) - rrt(032) * density(09) 
  pd(09,09) = pd(09,09) - rrt(032) * density(01) 
  pd(12,01) = pd(12,01) + rrt(032) * density(09) 
  pd(12,09) = pd(12,09) + rrt(032) * density(01) 
  pd(01,01) = pd(01,01) + rrt(033) * density(11) 
  pd(01,11) = pd(01,11) + rrt(033) * density(01) 
  pd(11,01) = pd(11,01) - rrt(033) * density(11) 
  pd(11,11) = pd(11,11) - rrt(033) * density(01) 
  pd(12,01) = pd(12,01) + rrt(033) * density(11) 
  pd(12,11) = pd(12,11) + rrt(033) * density(01) 
  pd(01,01) = pd(01,01) + rrt(034) * density(10) 
  pd(01,10) = pd(01,10) + rrt(034) * density(01) 
  pd(10,01) = pd(10,01) - rrt(034) * density(10) 
  pd(10,10) = pd(10,10) - rrt(034) * density(01) 
  pd(12,01) = pd(12,01) + rrt(034) * density(10) 
  pd(12,10) = pd(12,10) + rrt(034) * density(01) 
  pd(01,01) = pd(01,01) + rrt(035) * density(09) 
  pd(01,09) = pd(01,09) + rrt(035) * density(01) 
  pd(08,01) = pd(08,01) + rrt(035) * density(09) 
  pd(08,09) = pd(08,09) + rrt(035) * density(01) 
  pd(09,01) = pd(09,01) - rrt(035) * density(09) 
  pd(09,09) = pd(09,09) - rrt(035) * density(01) 
  pd(43,01) = pd(43,01) + rrt(035) * density(09) 
  pd(43,09) = pd(43,09) + rrt(035) * density(01) 
  pd(01,01) = pd(01,01) + rrt(036) * density(11) 
  pd(01,11) = pd(01,11) + rrt(036) * density(01) 
  pd(08,01) = pd(08,01) + rrt(036) * density(11) 
  pd(08,11) = pd(08,11) + rrt(036) * density(01) 
  pd(11,01) = pd(11,01) - rrt(036) * density(11) 
  pd(11,11) = pd(11,11) - rrt(036) * density(01) 
  pd(43,01) = pd(43,01) + rrt(036) * density(11) 
  pd(43,11) = pd(43,11) + rrt(036) * density(01) 
  pd(01,01) = pd(01,01) + rrt(037) * density(10) 
  pd(01,10) = pd(01,10) + rrt(037) * density(01) 
  pd(08,01) = pd(08,01) + rrt(037) * density(10) 
  pd(08,10) = pd(08,10) + rrt(037) * density(01) 
  pd(10,01) = pd(10,01) - rrt(037) * density(10) 
  pd(10,10) = pd(10,10) - rrt(037) * density(01) 
  pd(43,01) = pd(43,01) + rrt(037) * density(10) 
  pd(43,10) = pd(43,10) + rrt(037) * density(01) 
  pd(01,01) = pd(01,01) + rrt(038) * density(09) 
  pd(01,09) = pd(01,09) + rrt(038) * density(01) 
  pd(06,01) = pd(06,01) + rrt(038) * density(09) 
  pd(06,09) = pd(06,09) + rrt(038) * density(01) 
  pd(09,01) = pd(09,01) - rrt(038) * density(09) 
  pd(09,09) = pd(09,09) - rrt(038) * density(01) 
  pd(44,01) = pd(44,01) + rrt(038) * density(09) 
  pd(44,09) = pd(44,09) + rrt(038) * density(01) 
  pd(01,01) = pd(01,01) + rrt(039) * density(11) 
  pd(01,11) = pd(01,11) + rrt(039) * density(01) 
  pd(06,01) = pd(06,01) + rrt(039) * density(11) 
  pd(06,11) = pd(06,11) + rrt(039) * density(01) 
  pd(11,01) = pd(11,01) - rrt(039) * density(11) 
  pd(11,11) = pd(11,11) - rrt(039) * density(01) 
  pd(44,01) = pd(44,01) + rrt(039) * density(11) 
  pd(44,11) = pd(44,11) + rrt(039) * density(01) 
  pd(01,01) = pd(01,01) + rrt(040) * density(10) 
  pd(01,10) = pd(01,10) + rrt(040) * density(01) 
  pd(06,01) = pd(06,01) + rrt(040) * density(10) 
  pd(06,10) = pd(06,10) + rrt(040) * density(01) 
  pd(10,01) = pd(10,01) - rrt(040) * density(10) 
  pd(10,10) = pd(10,10) - rrt(040) * density(01) 
  pd(44,01) = pd(44,01) + rrt(040) * density(10) 
  pd(44,10) = pd(44,10) + rrt(040) * density(01) 
  pd(01,01) = pd(01,01) + rrt(041) * density(09) 
  pd(01,09) = pd(01,09) + rrt(041) * density(01) 
  pd(04,01) = pd(04,01) + rrt(041) * density(09) 
  pd(04,09) = pd(04,09) + rrt(041) * density(01) 
  pd(09,01) = pd(09,01) - rrt(041) * density(09) 
  pd(09,09) = pd(09,09) - rrt(041) * density(01) 
  pd(43,01) = pd(43,01) + rrt(041) * density(09) 
  pd(43,09) = pd(43,09) + rrt(041) * density(01) 
  pd(44,01) = pd(44,01) + rrt(041) * density(09) 
  pd(44,09) = pd(44,09) + rrt(041) * density(01) 
  pd(01,01) = pd(01,01) + rrt(042) * density(11) 
  pd(01,11) = pd(01,11) + rrt(042) * density(01) 
  pd(04,01) = pd(04,01) + rrt(042) * density(11) 
  pd(04,11) = pd(04,11) + rrt(042) * density(01) 
  pd(11,01) = pd(11,01) - rrt(042) * density(11) 
  pd(11,11) = pd(11,11) - rrt(042) * density(01) 
  pd(43,01) = pd(43,01) + rrt(042) * density(11) 
  pd(43,11) = pd(43,11) + rrt(042) * density(01) 
  pd(44,01) = pd(44,01) + rrt(042) * density(11) 
  pd(44,11) = pd(44,11) + rrt(042) * density(01) 
  pd(01,01) = pd(01,01) + rrt(043) * density(10) 
  pd(01,10) = pd(01,10) + rrt(043) * density(01) 
  pd(04,01) = pd(04,01) + rrt(043) * density(10) 
  pd(04,10) = pd(04,10) + rrt(043) * density(01) 
  pd(10,01) = pd(10,01) - rrt(043) * density(10) 
  pd(10,10) = pd(10,10) - rrt(043) * density(01) 
  pd(43,01) = pd(43,01) + rrt(043) * density(10) 
  pd(43,10) = pd(43,10) + rrt(043) * density(01) 
  pd(44,01) = pd(44,01) + rrt(043) * density(10) 
  pd(44,10) = pd(44,10) + rrt(043) * density(01) 
  pd(01,01) = pd(01,01) + rrt(044) * density(07) 
  pd(01,07) = pd(01,07) + rrt(044) * density(01) 
  pd(07,01) = pd(07,01) - rrt(044) * density(07) 
  pd(07,07) = pd(07,07) - rrt(044) * density(01) 
  pd(08,01) = pd(08,01) + rrt(044) * density(07) 
  pd(08,07) = pd(08,07) + rrt(044) * density(01) 
  pd(01,01) = pd(01,01) + rrt(045) * density(07) 
  pd(01,07) = pd(01,07) + rrt(045) * density(01) 
  pd(06,01) = pd(06,01) + rrt(045) * density(07) 
  pd(06,07) = pd(06,07) + rrt(045) * density(01) 
  pd(07,01) = pd(07,01) - rrt(045) * density(07) 
  pd(07,07) = pd(07,07) - rrt(045) * density(01) 
  pd(43,01) = pd(43,01) + rrt(045) * density(07) 
  pd(43,07) = pd(43,07) + rrt(045) * density(01) 
  pd(01,01) = pd(01,01) + rrt(046) * density(07) 
  pd(01,07) = pd(01,07) + rrt(046) * density(01) 
  pd(04,01) = pd(04,01) + rrt(046) * density(07) 
  pd(04,07) = pd(04,07) + rrt(046) * density(01) 
  pd(07,01) = pd(07,01) - rrt(046) * density(07) 
  pd(07,07) = pd(07,07) - rrt(046) * density(01) 
  pd(44,01) = pd(44,01) + rrt(046) * density(07) 
  pd(44,07) = pd(44,07) + rrt(046) * density(01) 
  pd(01,01) = pd(01,01) + rrt(047) * density(05) 
  pd(01,05) = pd(01,05) + rrt(047) * density(01) 
  pd(05,01) = pd(05,01) - rrt(047) * density(05) 
  pd(05,05) = pd(05,05) - rrt(047) * density(01) 
  pd(06,01) = pd(06,01) + rrt(047) * density(05) 
  pd(06,05) = pd(06,05) + rrt(047) * density(01) 
  pd(01,01) = pd(01,01) + rrt(048) * density(05) 
  pd(01,05) = pd(01,05) + rrt(048) * density(01) 
  pd(04,01) = pd(04,01) + rrt(048) * density(05) 
  pd(04,05) = pd(04,05) + rrt(048) * density(01) 
  pd(05,01) = pd(05,01) - rrt(048) * density(05) 
  pd(05,05) = pd(05,05) - rrt(048) * density(01) 
  pd(43,01) = pd(43,01) + rrt(048) * density(05) 
  pd(43,05) = pd(43,05) + rrt(048) * density(01) 
  pd(01,01) = pd(01,01) + rrt(049) * density(03) 
  pd(01,03) = pd(01,03) + rrt(049) * density(01) 
  pd(03,01) = pd(03,01) - rrt(049) * density(03) 
  pd(03,03) = pd(03,03) - rrt(049) * density(01) 
  pd(04,01) = pd(04,01) + rrt(049) * density(03) 
  pd(04,03) = pd(04,03) + rrt(049) * density(01) 
  pd(24,01) = pd(24,01) + rrt(050) * density(26) 
  pd(24,26) = pd(24,26) + rrt(050) * density(01) 
  pd(26,01) = pd(26,01) - rrt(050) * density(26) 
  pd(26,26) = pd(26,26) - rrt(050) * density(01) 
  pd(43,01) = pd(43,01) + rrt(050) * density(26) 
  pd(43,26) = pd(43,26) + rrt(050) * density(01) 
  pd(24,01) = pd(24,01) + rrt(051) * density(28) 
  pd(24,28) = pd(24,28) + rrt(051) * density(01) 
  pd(28,01) = pd(28,01) - rrt(051) * density(28) 
  pd(28,28) = pd(28,28) - rrt(051) * density(01) 
  pd(43,01) = pd(43,01) + rrt(051) * density(28) 
  pd(43,28) = pd(43,28) + rrt(051) * density(01) 
  pd(24,01) = pd(24,01) + rrt(052) * density(27) 
  pd(24,27) = pd(24,27) + rrt(052) * density(01) 
  pd(27,01) = pd(27,01) - rrt(052) * density(27) 
  pd(27,27) = pd(27,27) - rrt(052) * density(01) 
  pd(43,01) = pd(43,01) + rrt(052) * density(27) 
  pd(43,27) = pd(43,27) + rrt(052) * density(01) 
  pd(20,01) = pd(20,01) + rrt(053) * density(26) 
  pd(20,26) = pd(20,26) + rrt(053) * density(01) 
  pd(26,01) = pd(26,01) - rrt(053) * density(26) 
  pd(26,26) = pd(26,26) - rrt(053) * density(01) 
  pd(44,01) = pd(44,01) + rrt(053) * density(26) 
  pd(44,26) = pd(44,26) + rrt(053) * density(01) 
  pd(20,01) = pd(20,01) + rrt(054) * density(28) 
  pd(20,28) = pd(20,28) + rrt(054) * density(01) 
  pd(28,01) = pd(28,01) - rrt(054) * density(28) 
  pd(28,28) = pd(28,28) - rrt(054) * density(01) 
  pd(44,01) = pd(44,01) + rrt(054) * density(28) 
  pd(44,28) = pd(44,28) + rrt(054) * density(01) 
  pd(20,01) = pd(20,01) + rrt(055) * density(27) 
  pd(20,27) = pd(20,27) + rrt(055) * density(01) 
  pd(27,01) = pd(27,01) - rrt(055) * density(27) 
  pd(27,27) = pd(27,27) - rrt(055) * density(01) 
  pd(44,01) = pd(44,01) + rrt(055) * density(27) 
  pd(44,27) = pd(44,27) + rrt(055) * density(01) 
  pd(18,01) = pd(18,01) + rrt(056) * density(26) 
  pd(18,26) = pd(18,26) + rrt(056) * density(01) 
  pd(26,01) = pd(26,01) - rrt(056) * density(26) 
  pd(26,26) = pd(26,26) - rrt(056) * density(01) 
  pd(43,01) = pd(43,01) + rrt(056) * density(26) 
  pd(43,26) = pd(43,26) + rrt(056) * density(01) 
  pd(44,01) = pd(44,01) + rrt(056) * density(26) 
  pd(44,26) = pd(44,26) + rrt(056) * density(01) 
  pd(18,01) = pd(18,01) + rrt(057) * density(28) 
  pd(18,28) = pd(18,28) + rrt(057) * density(01) 
  pd(28,01) = pd(28,01) - rrt(057) * density(28) 
  pd(28,28) = pd(28,28) - rrt(057) * density(01) 
  pd(43,01) = pd(43,01) + rrt(057) * density(28) 
  pd(43,28) = pd(43,28) + rrt(057) * density(01) 
  pd(44,01) = pd(44,01) + rrt(057) * density(28) 
  pd(44,28) = pd(44,28) + rrt(057) * density(01) 
  pd(18,01) = pd(18,01) + rrt(058) * density(27) 
  pd(18,27) = pd(18,27) + rrt(058) * density(01) 
  pd(27,01) = pd(27,01) - rrt(058) * density(27) 
  pd(27,27) = pd(27,27) - rrt(058) * density(01) 
  pd(43,01) = pd(43,01) + rrt(058) * density(27) 
  pd(43,27) = pd(43,27) + rrt(058) * density(01) 
  pd(44,01) = pd(44,01) + rrt(058) * density(27) 
  pd(44,27) = pd(44,27) + rrt(058) * density(01) 
  pd(13,01) = pd(13,01) + rrt(059) * density(26) 
  pd(13,26) = pd(13,26) + rrt(059) * density(01) 
  pd(26,01) = pd(26,01) - rrt(059) * density(26) 
  pd(26,26) = pd(26,26) - rrt(059) * density(01) 
  pd(44,01) = pd(44,01) + rrt(059) * density(26) * 2.0d0
  pd(44,26) = pd(44,26) + rrt(059) * density(01) * 2.0d0
  pd(13,01) = pd(13,01) + rrt(060) * density(28) 
  pd(13,28) = pd(13,28) + rrt(060) * density(01) 
  pd(28,01) = pd(28,01) - rrt(060) * density(28) 
  pd(28,28) = pd(28,28) - rrt(060) * density(01) 
  pd(44,01) = pd(44,01) + rrt(060) * density(28) * 2.0d0
  pd(44,28) = pd(44,28) + rrt(060) * density(01) * 2.0d0
  pd(13,01) = pd(13,01) + rrt(061) * density(27) 
  pd(13,27) = pd(13,27) + rrt(061) * density(01) 
  pd(27,01) = pd(27,01) - rrt(061) * density(27) 
  pd(27,27) = pd(27,27) - rrt(061) * density(01) 
  pd(44,01) = pd(44,01) + rrt(061) * density(27) * 2.0d0
  pd(44,27) = pd(44,27) + rrt(061) * density(01) * 2.0d0
  pd(05,01) = pd(05,01) + rrt(062) * density(26) 
  pd(05,26) = pd(05,26) + rrt(062) * density(01) 
  pd(09,01) = pd(09,01) + rrt(062) * density(26) 
  pd(09,26) = pd(09,26) + rrt(062) * density(01) 
  pd(26,01) = pd(26,01) - rrt(062) * density(26) 
  pd(26,26) = pd(26,26) - rrt(062) * density(01) 
  pd(05,01) = pd(05,01) + rrt(063) * density(28) 
  pd(05,28) = pd(05,28) + rrt(063) * density(01) 
  pd(09,01) = pd(09,01) + rrt(063) * density(28) 
  pd(09,28) = pd(09,28) + rrt(063) * density(01) 
  pd(28,01) = pd(28,01) - rrt(063) * density(28) 
  pd(28,28) = pd(28,28) - rrt(063) * density(01) 
  pd(05,01) = pd(05,01) + rrt(064) * density(27) 
  pd(05,27) = pd(05,27) + rrt(064) * density(01) 
  pd(09,01) = pd(09,01) + rrt(064) * density(27) 
  pd(09,27) = pd(09,27) + rrt(064) * density(01) 
  pd(27,01) = pd(27,01) - rrt(064) * density(27) 
  pd(27,27) = pd(27,27) - rrt(064) * density(01) 
  pd(07,01) = pd(07,01) + rrt(065) * density(26) * 2.0d0
  pd(07,26) = pd(07,26) + rrt(065) * density(01) * 2.0d0
  pd(26,01) = pd(26,01) - rrt(065) * density(26) 
  pd(26,26) = pd(26,26) - rrt(065) * density(01) 
  pd(07,01) = pd(07,01) + rrt(066) * density(28) * 2.0d0
  pd(07,28) = pd(07,28) + rrt(066) * density(01) * 2.0d0
  pd(28,01) = pd(28,01) - rrt(066) * density(28) 
  pd(28,28) = pd(28,28) - rrt(066) * density(01) 
  pd(07,01) = pd(07,01) + rrt(067) * density(27) * 2.0d0
  pd(07,27) = pd(07,27) + rrt(067) * density(01) * 2.0d0
  pd(27,01) = pd(27,01) - rrt(067) * density(27) 
  pd(27,27) = pd(27,27) - rrt(067) * density(01) 
  pd(20,01) = pd(20,01) + rrt(068) * density(24) 
  pd(20,24) = pd(20,24) + rrt(068) * density(01) 
  pd(24,01) = pd(24,01) - rrt(068) * density(24) 
  pd(24,24) = pd(24,24) - rrt(068) * density(01) 
  pd(43,01) = pd(43,01) + rrt(068) * density(24) 
  pd(43,24) = pd(43,24) + rrt(068) * density(01) 
  pd(18,01) = pd(18,01) + rrt(069) * density(24) 
  pd(18,24) = pd(18,24) + rrt(069) * density(01) 
  pd(24,01) = pd(24,01) - rrt(069) * density(24) 
  pd(24,24) = pd(24,24) - rrt(069) * density(01) 
  pd(44,01) = pd(44,01) + rrt(069) * density(24) 
  pd(44,24) = pd(44,24) + rrt(069) * density(01) 
  pd(18,01) = pd(18,01) + rrt(070) * density(24) 
  pd(18,24) = pd(18,24) + rrt(070) * density(01) 
  pd(24,01) = pd(24,01) - rrt(070) * density(24) 
  pd(24,24) = pd(24,24) - rrt(070) * density(01) 
  pd(43,01) = pd(43,01) + rrt(070) * density(24) * 2.0d0
  pd(43,24) = pd(43,24) + rrt(070) * density(01) * 2.0d0
  pd(13,01) = pd(13,01) + rrt(071) * density(24) 
  pd(13,24) = pd(13,24) + rrt(071) * density(01) 
  pd(24,01) = pd(24,01) - rrt(071) * density(24) 
  pd(24,24) = pd(24,24) - rrt(071) * density(01) 
  pd(43,01) = pd(43,01) + rrt(071) * density(24) 
  pd(43,24) = pd(43,24) + rrt(071) * density(01) 
  pd(44,01) = pd(44,01) + rrt(071) * density(24) 
  pd(44,24) = pd(44,24) + rrt(071) * density(01) 
  pd(03,01) = pd(03,01) + rrt(072) * density(24) 
  pd(03,24) = pd(03,24) + rrt(072) * density(01) 
  pd(09,01) = pd(09,01) + rrt(072) * density(24) 
  pd(09,24) = pd(09,24) + rrt(072) * density(01) 
  pd(24,01) = pd(24,01) - rrt(072) * density(24) 
  pd(24,24) = pd(24,24) - rrt(072) * density(01) 
  pd(05,01) = pd(05,01) + rrt(073) * density(24) 
  pd(05,24) = pd(05,24) + rrt(073) * density(01) 
  pd(07,01) = pd(07,01) + rrt(073) * density(24) 
  pd(07,24) = pd(07,24) + rrt(073) * density(01) 
  pd(24,01) = pd(24,01) - rrt(073) * density(24) 
  pd(24,24) = pd(24,24) - rrt(073) * density(01) 
  pd(18,01) = pd(18,01) + rrt(074) * density(20) 
  pd(18,20) = pd(18,20) + rrt(074) * density(01) 
  pd(20,01) = pd(20,01) - rrt(074) * density(20) 
  pd(20,20) = pd(20,20) - rrt(074) * density(01) 
  pd(43,01) = pd(43,01) + rrt(074) * density(20) 
  pd(43,20) = pd(43,20) + rrt(074) * density(01) 
  pd(18,01) = pd(18,01) + rrt(075) * density(21) 
  pd(18,21) = pd(18,21) + rrt(075) * density(01) 
  pd(21,01) = pd(21,01) - rrt(075) * density(21) 
  pd(21,21) = pd(21,21) - rrt(075) * density(01) 
  pd(43,01) = pd(43,01) + rrt(075) * density(21) 
  pd(43,21) = pd(43,21) + rrt(075) * density(01) 
  pd(18,01) = pd(18,01) + rrt(076) * density(22) 
  pd(18,22) = pd(18,22) + rrt(076) * density(01) 
  pd(22,01) = pd(22,01) - rrt(076) * density(22) 
  pd(22,22) = pd(22,22) - rrt(076) * density(01) 
  pd(43,01) = pd(43,01) + rrt(076) * density(22) 
  pd(43,22) = pd(43,22) + rrt(076) * density(01) 
  pd(13,01) = pd(13,01) + rrt(077) * density(20) 
  pd(13,20) = pd(13,20) + rrt(077) * density(01) 
  pd(20,01) = pd(20,01) - rrt(077) * density(20) 
  pd(20,20) = pd(20,20) - rrt(077) * density(01) 
  pd(44,01) = pd(44,01) + rrt(077) * density(20) 
  pd(44,20) = pd(44,20) + rrt(077) * density(01) 
  pd(13,01) = pd(13,01) + rrt(078) * density(21) 
  pd(13,21) = pd(13,21) + rrt(078) * density(01) 
  pd(21,01) = pd(21,01) - rrt(078) * density(21) 
  pd(21,21) = pd(21,21) - rrt(078) * density(01) 
  pd(44,01) = pd(44,01) + rrt(078) * density(21) 
  pd(44,21) = pd(44,21) + rrt(078) * density(01) 
  pd(13,01) = pd(13,01) + rrt(079) * density(22) 
  pd(13,22) = pd(13,22) + rrt(079) * density(01) 
  pd(22,01) = pd(22,01) - rrt(079) * density(22) 
  pd(22,22) = pd(22,22) - rrt(079) * density(01) 
  pd(44,01) = pd(44,01) + rrt(079) * density(22) 
  pd(44,22) = pd(44,22) + rrt(079) * density(01) 
  pd(13,01) = pd(13,01) + rrt(080) * density(20) 
  pd(13,20) = pd(13,20) + rrt(080) * density(01) 
  pd(20,01) = pd(20,01) - rrt(080) * density(20) 
  pd(20,20) = pd(20,20) - rrt(080) * density(01) 
  pd(43,01) = pd(43,01) + rrt(080) * density(20) * 2.0d0
  pd(43,20) = pd(43,20) + rrt(080) * density(01) * 2.0d0
  pd(13,01) = pd(13,01) + rrt(081) * density(21) 
  pd(13,21) = pd(13,21) + rrt(081) * density(01) 
  pd(21,01) = pd(21,01) - rrt(081) * density(21) 
  pd(21,21) = pd(21,21) - rrt(081) * density(01) 
  pd(43,01) = pd(43,01) + rrt(081) * density(21) * 2.0d0
  pd(43,21) = pd(43,21) + rrt(081) * density(01) * 2.0d0
  pd(13,01) = pd(13,01) + rrt(082) * density(22) 
  pd(13,22) = pd(13,22) + rrt(082) * density(01) 
  pd(22,01) = pd(22,01) - rrt(082) * density(22) 
  pd(22,22) = pd(22,22) - rrt(082) * density(01) 
  pd(43,01) = pd(43,01) + rrt(082) * density(22) * 2.0d0
  pd(43,22) = pd(43,22) + rrt(082) * density(01) * 2.0d0
  pd(03,01) = pd(03,01) + rrt(083) * density(20) 
  pd(03,20) = pd(03,20) + rrt(083) * density(01) 
  pd(07,01) = pd(07,01) + rrt(083) * density(20) 
  pd(07,20) = pd(07,20) + rrt(083) * density(01) 
  pd(20,01) = pd(20,01) - rrt(083) * density(20) 
  pd(20,20) = pd(20,20) - rrt(083) * density(01) 
  pd(03,01) = pd(03,01) + rrt(084) * density(21) 
  pd(03,21) = pd(03,21) + rrt(084) * density(01) 
  pd(07,01) = pd(07,01) + rrt(084) * density(21) 
  pd(07,21) = pd(07,21) + rrt(084) * density(01) 
  pd(21,01) = pd(21,01) - rrt(084) * density(21) 
  pd(21,21) = pd(21,21) - rrt(084) * density(01) 
  pd(03,01) = pd(03,01) + rrt(085) * density(22) 
  pd(03,22) = pd(03,22) + rrt(085) * density(01) 
  pd(07,01) = pd(07,01) + rrt(085) * density(22) 
  pd(07,22) = pd(07,22) + rrt(085) * density(01) 
  pd(22,01) = pd(22,01) - rrt(085) * density(22) 
  pd(22,22) = pd(22,22) - rrt(085) * density(01) 
  pd(05,01) = pd(05,01) + rrt(086) * density(20) * 2.0d0
  pd(05,20) = pd(05,20) + rrt(086) * density(01) * 2.0d0
  pd(20,01) = pd(20,01) - rrt(086) * density(20) 
  pd(20,20) = pd(20,20) - rrt(086) * density(01) 
  pd(05,01) = pd(05,01) + rrt(087) * density(21) * 2.0d0
  pd(05,21) = pd(05,21) + rrt(087) * density(01) * 2.0d0
  pd(21,01) = pd(21,01) - rrt(087) * density(21) 
  pd(21,21) = pd(21,21) - rrt(087) * density(01) 
  pd(05,01) = pd(05,01) + rrt(088) * density(22) * 2.0d0
  pd(05,22) = pd(05,22) + rrt(088) * density(01) * 2.0d0
  pd(22,01) = pd(22,01) - rrt(088) * density(22) 
  pd(22,22) = pd(22,22) - rrt(088) * density(01) 
  pd(13,01) = pd(13,01) + rrt(089) * density(18) 
  pd(13,18) = pd(13,18) + rrt(089) * density(01) 
  pd(18,01) = pd(18,01) - rrt(089) * density(18) 
  pd(18,18) = pd(18,18) - rrt(089) * density(01) 
  pd(43,01) = pd(43,01) + rrt(089) * density(18) 
  pd(43,18) = pd(43,18) + rrt(089) * density(01) 
  pd(03,01) = pd(03,01) + rrt(090) * density(18) 
  pd(03,18) = pd(03,18) + rrt(090) * density(01) 
  pd(05,01) = pd(05,01) + rrt(090) * density(18) 
  pd(05,18) = pd(05,18) + rrt(090) * density(01) 
  pd(18,01) = pd(18,01) - rrt(090) * density(18) 
  pd(18,18) = pd(18,18) - rrt(090) * density(01) 
  pd(03,01) = pd(03,01) + rrt(091) * density(13) * 2.0d0
  pd(03,13) = pd(03,13) + rrt(091) * density(01) * 2.0d0
  pd(13,01) = pd(13,01) - rrt(091) * density(13) 
  pd(13,13) = pd(13,13) - rrt(091) * density(01) 
  pd(03,01) = pd(03,01) + rrt(092) * density(16) * 2.0d0
  pd(03,16) = pd(03,16) + rrt(092) * density(01) * 2.0d0
  pd(16,01) = pd(16,01) - rrt(092) * density(16) 
  pd(16,16) = pd(16,16) - rrt(092) * density(01) 
  pd(03,01) = pd(03,01) + rrt(093) * density(15) * 2.0d0
  pd(03,15) = pd(03,15) + rrt(093) * density(01) * 2.0d0
  pd(15,01) = pd(15,01) - rrt(093) * density(15) 
  pd(15,15) = pd(15,15) - rrt(093) * density(01) 
  pd(03,01) = pd(03,01) + rrt(094) * density(14) * 2.0d0
  pd(03,14) = pd(03,14) + rrt(094) * density(01) * 2.0d0
  pd(14,01) = pd(14,01) - rrt(094) * density(14) 
  pd(14,14) = pd(14,14) - rrt(094) * density(01) 
  pd(01,01) = pd(01,01) + rrt(095) * density(26) 
  pd(01,26) = pd(01,26) + rrt(095) * density(01) 
  pd(26,01) = pd(26,01) - rrt(095) * density(26) 
  pd(26,26) = pd(26,26) - rrt(095) * density(01) 
  pd(29,01) = pd(29,01) + rrt(095) * density(26) 
  pd(29,26) = pd(29,26) + rrt(095) * density(01) 
  pd(01,01) = pd(01,01) + rrt(096) * density(28) 
  pd(01,28) = pd(01,28) + rrt(096) * density(01) 
  pd(28,01) = pd(28,01) - rrt(096) * density(28) 
  pd(28,28) = pd(28,28) - rrt(096) * density(01) 
  pd(29,01) = pd(29,01) + rrt(096) * density(28) 
  pd(29,28) = pd(29,28) + rrt(096) * density(01) 
  pd(01,01) = pd(01,01) + rrt(097) * density(27) 
  pd(01,27) = pd(01,27) + rrt(097) * density(01) 
  pd(27,01) = pd(27,01) - rrt(097) * density(27) 
  pd(27,27) = pd(27,27) - rrt(097) * density(01) 
  pd(29,01) = pd(29,01) + rrt(097) * density(27) 
  pd(29,27) = pd(29,27) + rrt(097) * density(01) 
  pd(01,01) = pd(01,01) + rrt(098) * density(26) 
  pd(01,26) = pd(01,26) + rrt(098) * density(01) 
  pd(25,01) = pd(25,01) + rrt(098) * density(26) 
  pd(25,26) = pd(25,26) + rrt(098) * density(01) 
  pd(26,01) = pd(26,01) - rrt(098) * density(26) 
  pd(26,26) = pd(26,26) - rrt(098) * density(01) 
  pd(43,01) = pd(43,01) + rrt(098) * density(26) 
  pd(43,26) = pd(43,26) + rrt(098) * density(01) 
  pd(01,01) = pd(01,01) + rrt(099) * density(28) 
  pd(01,28) = pd(01,28) + rrt(099) * density(01) 
  pd(25,01) = pd(25,01) + rrt(099) * density(28) 
  pd(25,28) = pd(25,28) + rrt(099) * density(01) 
  pd(28,01) = pd(28,01) - rrt(099) * density(28) 
  pd(28,28) = pd(28,28) - rrt(099) * density(01) 
  pd(43,01) = pd(43,01) + rrt(099) * density(28) 
  pd(43,28) = pd(43,28) + rrt(099) * density(01) 
  pd(01,01) = pd(01,01) + rrt(100) * density(27) 
  pd(01,27) = pd(01,27) + rrt(100) * density(01) 
  pd(25,01) = pd(25,01) + rrt(100) * density(27) 
  pd(25,27) = pd(25,27) + rrt(100) * density(01) 
  pd(27,01) = pd(27,01) - rrt(100) * density(27) 
  pd(27,27) = pd(27,27) - rrt(100) * density(01) 
  pd(43,01) = pd(43,01) + rrt(100) * density(27) 
  pd(43,27) = pd(43,27) + rrt(100) * density(01) 
  pd(01,01) = pd(01,01) + rrt(101) * density(26) 
  pd(01,26) = pd(01,26) + rrt(101) * density(01) 
  pd(23,01) = pd(23,01) + rrt(101) * density(26) 
  pd(23,26) = pd(23,26) + rrt(101) * density(01) 
  pd(26,01) = pd(26,01) - rrt(101) * density(26) 
  pd(26,26) = pd(26,26) - rrt(101) * density(01) 
  pd(44,01) = pd(44,01) + rrt(101) * density(26) 
  pd(44,26) = pd(44,26) + rrt(101) * density(01) 
  pd(01,01) = pd(01,01) + rrt(102) * density(28) 
  pd(01,28) = pd(01,28) + rrt(102) * density(01) 
  pd(23,01) = pd(23,01) + rrt(102) * density(28) 
  pd(23,28) = pd(23,28) + rrt(102) * density(01) 
  pd(28,01) = pd(28,01) - rrt(102) * density(28) 
  pd(28,28) = pd(28,28) - rrt(102) * density(01) 
  pd(44,01) = pd(44,01) + rrt(102) * density(28) 
  pd(44,28) = pd(44,28) + rrt(102) * density(01) 
  pd(01,01) = pd(01,01) + rrt(103) * density(27) 
  pd(01,27) = pd(01,27) + rrt(103) * density(01) 
  pd(23,01) = pd(23,01) + rrt(103) * density(27) 
  pd(23,27) = pd(23,27) + rrt(103) * density(01) 
  pd(27,01) = pd(27,01) - rrt(103) * density(27) 
  pd(27,27) = pd(27,27) - rrt(103) * density(01) 
  pd(44,01) = pd(44,01) + rrt(103) * density(27) 
  pd(44,27) = pd(44,27) + rrt(103) * density(01) 
  pd(01,01) = pd(01,01) + rrt(104) * density(26) 
  pd(01,26) = pd(01,26) + rrt(104) * density(01) 
  pd(19,01) = pd(19,01) + rrt(104) * density(26) 
  pd(19,26) = pd(19,26) + rrt(104) * density(01) 
  pd(26,01) = pd(26,01) - rrt(104) * density(26) 
  pd(26,26) = pd(26,26) - rrt(104) * density(01) 
  pd(43,01) = pd(43,01) + rrt(104) * density(26) 
  pd(43,26) = pd(43,26) + rrt(104) * density(01) 
  pd(44,01) = pd(44,01) + rrt(104) * density(26) 
  pd(44,26) = pd(44,26) + rrt(104) * density(01) 
  pd(01,01) = pd(01,01) + rrt(105) * density(28) 
  pd(01,28) = pd(01,28) + rrt(105) * density(01) 
  pd(19,01) = pd(19,01) + rrt(105) * density(28) 
  pd(19,28) = pd(19,28) + rrt(105) * density(01) 
  pd(28,01) = pd(28,01) - rrt(105) * density(28) 
  pd(28,28) = pd(28,28) - rrt(105) * density(01) 
  pd(43,01) = pd(43,01) + rrt(105) * density(28) 
  pd(43,28) = pd(43,28) + rrt(105) * density(01) 
  pd(44,01) = pd(44,01) + rrt(105) * density(28) 
  pd(44,28) = pd(44,28) + rrt(105) * density(01) 
  pd(01,01) = pd(01,01) + rrt(106) * density(27) 
  pd(01,27) = pd(01,27) + rrt(106) * density(01) 
  pd(19,01) = pd(19,01) + rrt(106) * density(27) 
  pd(19,27) = pd(19,27) + rrt(106) * density(01) 
  pd(27,01) = pd(27,01) - rrt(106) * density(27) 
  pd(27,27) = pd(27,27) - rrt(106) * density(01) 
  pd(43,01) = pd(43,01) + rrt(106) * density(27) 
  pd(43,27) = pd(43,27) + rrt(106) * density(01) 
  pd(44,01) = pd(44,01) + rrt(106) * density(27) 
  pd(44,27) = pd(44,27) + rrt(106) * density(01) 
  pd(01,01) = pd(01,01) + rrt(107) * density(26) 
  pd(01,26) = pd(01,26) + rrt(107) * density(01) 
  pd(17,01) = pd(17,01) + rrt(107) * density(26) 
  pd(17,26) = pd(17,26) + rrt(107) * density(01) 
  pd(26,01) = pd(26,01) - rrt(107) * density(26) 
  pd(26,26) = pd(26,26) - rrt(107) * density(01) 
  pd(44,01) = pd(44,01) + rrt(107) * density(26) * 2.0d0
  pd(44,26) = pd(44,26) + rrt(107) * density(01) * 2.0d0
  pd(01,01) = pd(01,01) + rrt(108) * density(28) 
  pd(01,28) = pd(01,28) + rrt(108) * density(01) 
  pd(17,01) = pd(17,01) + rrt(108) * density(28) 
  pd(17,28) = pd(17,28) + rrt(108) * density(01) 
  pd(28,01) = pd(28,01) - rrt(108) * density(28) 
  pd(28,28) = pd(28,28) - rrt(108) * density(01) 
  pd(44,01) = pd(44,01) + rrt(108) * density(28) * 2.0d0
  pd(44,28) = pd(44,28) + rrt(108) * density(01) * 2.0d0
  pd(01,01) = pd(01,01) + rrt(109) * density(27) 
  pd(01,27) = pd(01,27) + rrt(109) * density(01) 
  pd(17,01) = pd(17,01) + rrt(109) * density(27) 
  pd(17,27) = pd(17,27) + rrt(109) * density(01) 
  pd(27,01) = pd(27,01) - rrt(109) * density(27) 
  pd(27,27) = pd(27,27) - rrt(109) * density(01) 
  pd(44,01) = pd(44,01) + rrt(109) * density(27) * 2.0d0
  pd(44,27) = pd(44,27) + rrt(109) * density(01) * 2.0d0
  pd(01,01) = pd(01,01) + rrt(110) * density(26) 
  pd(01,26) = pd(01,26) + rrt(110) * density(01) 
  pd(07,01) = pd(07,01) + rrt(110) * density(26) 
  pd(07,26) = pd(07,26) + rrt(110) * density(01) 
  pd(08,01) = pd(08,01) + rrt(110) * density(26) 
  pd(08,26) = pd(08,26) + rrt(110) * density(01) 
  pd(26,01) = pd(26,01) - rrt(110) * density(26) 
  pd(26,26) = pd(26,26) - rrt(110) * density(01) 
  pd(01,01) = pd(01,01) + rrt(111) * density(28) 
  pd(01,28) = pd(01,28) + rrt(111) * density(01) 
  pd(07,01) = pd(07,01) + rrt(111) * density(28) 
  pd(07,28) = pd(07,28) + rrt(111) * density(01) 
  pd(08,01) = pd(08,01) + rrt(111) * density(28) 
  pd(08,28) = pd(08,28) + rrt(111) * density(01) 
  pd(28,01) = pd(28,01) - rrt(111) * density(28) 
  pd(28,28) = pd(28,28) - rrt(111) * density(01) 
  pd(01,01) = pd(01,01) + rrt(112) * density(27) 
  pd(01,27) = pd(01,27) + rrt(112) * density(01) 
  pd(07,01) = pd(07,01) + rrt(112) * density(27) 
  pd(07,27) = pd(07,27) + rrt(112) * density(01) 
  pd(08,01) = pd(08,01) + rrt(112) * density(27) 
  pd(08,27) = pd(08,27) + rrt(112) * density(01) 
  pd(27,01) = pd(27,01) - rrt(112) * density(27) 
  pd(27,27) = pd(27,27) - rrt(112) * density(01) 
  pd(01,01) = pd(01,01) + rrt(113) * density(26) 
  pd(01,26) = pd(01,26) + rrt(113) * density(01) 
  pd(06,01) = pd(06,01) + rrt(113) * density(26) 
  pd(06,26) = pd(06,26) + rrt(113) * density(01) 
  pd(09,01) = pd(09,01) + rrt(113) * density(26) 
  pd(09,26) = pd(09,26) + rrt(113) * density(01) 
  pd(26,01) = pd(26,01) - rrt(113) * density(26) 
  pd(26,26) = pd(26,26) - rrt(113) * density(01) 
  pd(01,01) = pd(01,01) + rrt(114) * density(28) 
  pd(01,28) = pd(01,28) + rrt(114) * density(01) 
  pd(06,01) = pd(06,01) + rrt(114) * density(28) 
  pd(06,28) = pd(06,28) + rrt(114) * density(01) 
  pd(09,01) = pd(09,01) + rrt(114) * density(28) 
  pd(09,28) = pd(09,28) + rrt(114) * density(01) 
  pd(28,01) = pd(28,01) - rrt(114) * density(28) 
  pd(28,28) = pd(28,28) - rrt(114) * density(01) 
  pd(01,01) = pd(01,01) + rrt(115) * density(27) 
  pd(01,27) = pd(01,27) + rrt(115) * density(01) 
  pd(06,01) = pd(06,01) + rrt(115) * density(27) 
  pd(06,27) = pd(06,27) + rrt(115) * density(01) 
  pd(09,01) = pd(09,01) + rrt(115) * density(27) 
  pd(09,27) = pd(09,27) + rrt(115) * density(01) 
  pd(27,01) = pd(27,01) - rrt(115) * density(27) 
  pd(27,27) = pd(27,27) - rrt(115) * density(01) 
  pd(01,01) = pd(01,01) + rrt(116) * density(24) 
  pd(01,24) = pd(01,24) + rrt(116) * density(01) 
  pd(24,01) = pd(24,01) - rrt(116) * density(24) 
  pd(24,24) = pd(24,24) - rrt(116) * density(01) 
  pd(25,01) = pd(25,01) + rrt(116) * density(24) 
  pd(25,24) = pd(25,24) + rrt(116) * density(01) 
  pd(01,01) = pd(01,01) + rrt(117) * density(24) 
  pd(01,24) = pd(01,24) + rrt(117) * density(01) 
  pd(23,01) = pd(23,01) + rrt(117) * density(24) 
  pd(23,24) = pd(23,24) + rrt(117) * density(01) 
  pd(24,01) = pd(24,01) - rrt(117) * density(24) 
  pd(24,24) = pd(24,24) - rrt(117) * density(01) 
  pd(43,01) = pd(43,01) + rrt(117) * density(24) 
  pd(43,24) = pd(43,24) + rrt(117) * density(01) 
  pd(01,01) = pd(01,01) + rrt(118) * density(24) 
  pd(01,24) = pd(01,24) + rrt(118) * density(01) 
  pd(19,01) = pd(19,01) + rrt(118) * density(24) 
  pd(19,24) = pd(19,24) + rrt(118) * density(01) 
  pd(24,01) = pd(24,01) - rrt(118) * density(24) 
  pd(24,24) = pd(24,24) - rrt(118) * density(01) 
  pd(44,01) = pd(44,01) + rrt(118) * density(24) 
  pd(44,24) = pd(44,24) + rrt(118) * density(01) 
  pd(01,01) = pd(01,01) + rrt(119) * density(24) 
  pd(01,24) = pd(01,24) + rrt(119) * density(01) 
  pd(17,01) = pd(17,01) + rrt(119) * density(24) 
  pd(17,24) = pd(17,24) + rrt(119) * density(01) 
  pd(24,01) = pd(24,01) - rrt(119) * density(24) 
  pd(24,24) = pd(24,24) - rrt(119) * density(01) 
  pd(43,01) = pd(43,01) + rrt(119) * density(24) 
  pd(43,24) = pd(43,24) + rrt(119) * density(01) 
  pd(44,01) = pd(44,01) + rrt(119) * density(24) 
  pd(44,24) = pd(44,24) + rrt(119) * density(01) 
  pd(01,01) = pd(01,01) + rrt(120) * density(24) 
  pd(01,24) = pd(01,24) + rrt(120) * density(01) 
  pd(05,01) = pd(05,01) + rrt(120) * density(24) 
  pd(05,24) = pd(05,24) + rrt(120) * density(01) 
  pd(08,01) = pd(08,01) + rrt(120) * density(24) 
  pd(08,24) = pd(08,24) + rrt(120) * density(01) 
  pd(24,01) = pd(24,01) - rrt(120) * density(24) 
  pd(24,24) = pd(24,24) - rrt(120) * density(01) 
  pd(01,01) = pd(01,01) + rrt(121) * density(24) 
  pd(01,24) = pd(01,24) + rrt(121) * density(01) 
  pd(06,01) = pd(06,01) + rrt(121) * density(24) 
  pd(06,24) = pd(06,24) + rrt(121) * density(01) 
  pd(07,01) = pd(07,01) + rrt(121) * density(24) 
  pd(07,24) = pd(07,24) + rrt(121) * density(01) 
  pd(24,01) = pd(24,01) - rrt(121) * density(24) 
  pd(24,24) = pd(24,24) - rrt(121) * density(01) 
  pd(01,01) = pd(01,01) + rrt(122) * density(24) 
  pd(01,24) = pd(01,24) + rrt(122) * density(01) 
  pd(04,01) = pd(04,01) + rrt(122) * density(24) 
  pd(04,24) = pd(04,24) + rrt(122) * density(01) 
  pd(09,01) = pd(09,01) + rrt(122) * density(24) 
  pd(09,24) = pd(09,24) + rrt(122) * density(01) 
  pd(24,01) = pd(24,01) - rrt(122) * density(24) 
  pd(24,24) = pd(24,24) - rrt(122) * density(01) 
  pd(01,01) = pd(01,01) + rrt(123) * density(20) 
  pd(01,20) = pd(01,20) + rrt(123) * density(01) 
  pd(20,01) = pd(20,01) - rrt(123) * density(20) 
  pd(20,20) = pd(20,20) - rrt(123) * density(01) 
  pd(23,01) = pd(23,01) + rrt(123) * density(20) 
  pd(23,20) = pd(23,20) + rrt(123) * density(01) 
  pd(01,01) = pd(01,01) + rrt(124) * density(21) 
  pd(01,21) = pd(01,21) + rrt(124) * density(01) 
  pd(21,01) = pd(21,01) - rrt(124) * density(21) 
  pd(21,21) = pd(21,21) - rrt(124) * density(01) 
  pd(23,01) = pd(23,01) + rrt(124) * density(21) 
  pd(23,21) = pd(23,21) + rrt(124) * density(01) 
  pd(01,01) = pd(01,01) + rrt(125) * density(22) 
  pd(01,22) = pd(01,22) + rrt(125) * density(01) 
  pd(22,01) = pd(22,01) - rrt(125) * density(22) 
  pd(22,22) = pd(22,22) - rrt(125) * density(01) 
  pd(23,01) = pd(23,01) + rrt(125) * density(22) 
  pd(23,22) = pd(23,22) + rrt(125) * density(01) 
  pd(01,01) = pd(01,01) + rrt(126) * density(20) 
  pd(01,20) = pd(01,20) + rrt(126) * density(01) 
  pd(19,01) = pd(19,01) + rrt(126) * density(20) 
  pd(19,20) = pd(19,20) + rrt(126) * density(01) 
  pd(20,01) = pd(20,01) - rrt(126) * density(20) 
  pd(20,20) = pd(20,20) - rrt(126) * density(01) 
  pd(43,01) = pd(43,01) + rrt(126) * density(20) 
  pd(43,20) = pd(43,20) + rrt(126) * density(01) 
  pd(01,01) = pd(01,01) + rrt(127) * density(21) 
  pd(01,21) = pd(01,21) + rrt(127) * density(01) 
  pd(19,01) = pd(19,01) + rrt(127) * density(21) 
  pd(19,21) = pd(19,21) + rrt(127) * density(01) 
  pd(21,01) = pd(21,01) - rrt(127) * density(21) 
  pd(21,21) = pd(21,21) - rrt(127) * density(01) 
  pd(43,01) = pd(43,01) + rrt(127) * density(21) 
  pd(43,21) = pd(43,21) + rrt(127) * density(01) 
  pd(01,01) = pd(01,01) + rrt(128) * density(22) 
  pd(01,22) = pd(01,22) + rrt(128) * density(01) 
  pd(19,01) = pd(19,01) + rrt(128) * density(22) 
  pd(19,22) = pd(19,22) + rrt(128) * density(01) 
  pd(22,01) = pd(22,01) - rrt(128) * density(22) 
  pd(22,22) = pd(22,22) - rrt(128) * density(01) 
  pd(43,01) = pd(43,01) + rrt(128) * density(22) 
  pd(43,22) = pd(43,22) + rrt(128) * density(01) 
  pd(01,01) = pd(01,01) + rrt(129) * density(20) 
  pd(01,20) = pd(01,20) + rrt(129) * density(01) 
  pd(03,01) = pd(03,01) + rrt(129) * density(20) 
  pd(03,20) = pd(03,20) + rrt(129) * density(01) 
  pd(08,01) = pd(08,01) + rrt(129) * density(20) 
  pd(08,20) = pd(08,20) + rrt(129) * density(01) 
  pd(20,01) = pd(20,01) - rrt(129) * density(20) 
  pd(20,20) = pd(20,20) - rrt(129) * density(01) 
  pd(01,01) = pd(01,01) + rrt(130) * density(21) 
  pd(01,21) = pd(01,21) + rrt(130) * density(01) 
  pd(03,01) = pd(03,01) + rrt(130) * density(21) 
  pd(03,21) = pd(03,21) + rrt(130) * density(01) 
  pd(08,01) = pd(08,01) + rrt(130) * density(21) 
  pd(08,21) = pd(08,21) + rrt(130) * density(01) 
  pd(21,01) = pd(21,01) - rrt(130) * density(21) 
  pd(21,21) = pd(21,21) - rrt(130) * density(01) 
  pd(01,01) = pd(01,01) + rrt(131) * density(22) 
  pd(01,22) = pd(01,22) + rrt(131) * density(01) 
  pd(03,01) = pd(03,01) + rrt(131) * density(22) 
  pd(03,22) = pd(03,22) + rrt(131) * density(01) 
  pd(08,01) = pd(08,01) + rrt(131) * density(22) 
  pd(08,22) = pd(08,22) + rrt(131) * density(01) 
  pd(22,01) = pd(22,01) - rrt(131) * density(22) 
  pd(22,22) = pd(22,22) - rrt(131) * density(01) 
  pd(01,01) = pd(01,01) + rrt(132) * density(20) 
  pd(01,20) = pd(01,20) + rrt(132) * density(01) 
  pd(05,01) = pd(05,01) + rrt(132) * density(20) 
  pd(05,20) = pd(05,20) + rrt(132) * density(01) 
  pd(06,01) = pd(06,01) + rrt(132) * density(20) 
  pd(06,20) = pd(06,20) + rrt(132) * density(01) 
  pd(20,01) = pd(20,01) - rrt(132) * density(20) 
  pd(20,20) = pd(20,20) - rrt(132) * density(01) 
  pd(01,01) = pd(01,01) + rrt(133) * density(21) 
  pd(01,21) = pd(01,21) + rrt(133) * density(01) 
  pd(05,01) = pd(05,01) + rrt(133) * density(21) 
  pd(05,21) = pd(05,21) + rrt(133) * density(01) 
  pd(06,01) = pd(06,01) + rrt(133) * density(21) 
  pd(06,21) = pd(06,21) + rrt(133) * density(01) 
  pd(21,01) = pd(21,01) - rrt(133) * density(21) 
  pd(21,21) = pd(21,21) - rrt(133) * density(01) 
  pd(01,01) = pd(01,01) + rrt(134) * density(22) 
  pd(01,22) = pd(01,22) + rrt(134) * density(01) 
  pd(05,01) = pd(05,01) + rrt(134) * density(22) 
  pd(05,22) = pd(05,22) + rrt(134) * density(01) 
  pd(06,01) = pd(06,01) + rrt(134) * density(22) 
  pd(06,22) = pd(06,22) + rrt(134) * density(01) 
  pd(22,01) = pd(22,01) - rrt(134) * density(22) 
  pd(22,22) = pd(22,22) - rrt(134) * density(01) 
  pd(01,01) = pd(01,01) + rrt(135) * density(20) 
  pd(01,20) = pd(01,20) + rrt(135) * density(01) 
  pd(04,01) = pd(04,01) + rrt(135) * density(20) 
  pd(04,20) = pd(04,20) + rrt(135) * density(01) 
  pd(07,01) = pd(07,01) + rrt(135) * density(20) 
  pd(07,20) = pd(07,20) + rrt(135) * density(01) 
  pd(20,01) = pd(20,01) - rrt(135) * density(20) 
  pd(20,20) = pd(20,20) - rrt(135) * density(01) 
  pd(01,01) = pd(01,01) + rrt(136) * density(21) 
  pd(01,21) = pd(01,21) + rrt(136) * density(01) 
  pd(04,01) = pd(04,01) + rrt(136) * density(21) 
  pd(04,21) = pd(04,21) + rrt(136) * density(01) 
  pd(07,01) = pd(07,01) + rrt(136) * density(21) 
  pd(07,21) = pd(07,21) + rrt(136) * density(01) 
  pd(21,01) = pd(21,01) - rrt(136) * density(21) 
  pd(21,21) = pd(21,21) - rrt(136) * density(01) 
  pd(01,01) = pd(01,01) + rrt(137) * density(22) 
  pd(01,22) = pd(01,22) + rrt(137) * density(01) 
  pd(04,01) = pd(04,01) + rrt(137) * density(22) 
  pd(04,22) = pd(04,22) + rrt(137) * density(01) 
  pd(07,01) = pd(07,01) + rrt(137) * density(22) 
  pd(07,22) = pd(07,22) + rrt(137) * density(01) 
  pd(22,01) = pd(22,01) - rrt(137) * density(22) 
  pd(22,22) = pd(22,22) - rrt(137) * density(01) 
  pd(01,01) = pd(01,01) + rrt(138) * density(18) 
  pd(01,18) = pd(01,18) + rrt(138) * density(01) 
  pd(18,01) = pd(18,01) - rrt(138) * density(18) 
  pd(18,18) = pd(18,18) - rrt(138) * density(01) 
  pd(19,01) = pd(19,01) + rrt(138) * density(18) 
  pd(19,18) = pd(19,18) + rrt(138) * density(01) 
  pd(01,01) = pd(01,01) + rrt(139) * density(18) 
  pd(01,18) = pd(01,18) + rrt(139) * density(01) 
  pd(17,01) = pd(17,01) + rrt(139) * density(18) 
  pd(17,18) = pd(17,18) + rrt(139) * density(01) 
  pd(18,01) = pd(18,01) - rrt(139) * density(18) 
  pd(18,18) = pd(18,18) - rrt(139) * density(01) 
  pd(43,01) = pd(43,01) + rrt(139) * density(18) 
  pd(43,18) = pd(43,18) + rrt(139) * density(01) 
  pd(01,01) = pd(01,01) + rrt(140) * density(18) 
  pd(01,18) = pd(01,18) + rrt(140) * density(01) 
  pd(03,01) = pd(03,01) + rrt(140) * density(18) 
  pd(03,18) = pd(03,18) + rrt(140) * density(01) 
  pd(06,01) = pd(06,01) + rrt(140) * density(18) 
  pd(06,18) = pd(06,18) + rrt(140) * density(01) 
  pd(18,01) = pd(18,01) - rrt(140) * density(18) 
  pd(18,18) = pd(18,18) - rrt(140) * density(01) 
  pd(01,01) = pd(01,01) + rrt(141) * density(18) 
  pd(01,18) = pd(01,18) + rrt(141) * density(01) 
  pd(04,01) = pd(04,01) + rrt(141) * density(18) 
  pd(04,18) = pd(04,18) + rrt(141) * density(01) 
  pd(05,01) = pd(05,01) + rrt(141) * density(18) 
  pd(05,18) = pd(05,18) + rrt(141) * density(01) 
  pd(18,01) = pd(18,01) - rrt(141) * density(18) 
  pd(18,18) = pd(18,18) - rrt(141) * density(01) 
  pd(01,01) = pd(01,01) + rrt(142) * density(18) 
  pd(01,18) = pd(01,18) + rrt(142) * density(01) 
  pd(13,01) = pd(13,01) + rrt(142) * density(18) 
  pd(13,18) = pd(13,18) + rrt(142) * density(01) 
  pd(18,01) = pd(18,01) - rrt(142) * density(18) 
  pd(18,18) = pd(18,18) - rrt(142) * density(01) 
  pd(47,01) = pd(47,01) + rrt(142) * density(18) 
  pd(47,18) = pd(47,18) + rrt(142) * density(01) 
  pd(01,01) = pd(01,01) + rrt(143) * density(13) 
  pd(01,13) = pd(01,13) + rrt(143) * density(01) 
  pd(13,01) = pd(13,01) - rrt(143) * density(13) 
  pd(13,13) = pd(13,13) - rrt(143) * density(01) 
  pd(17,01) = pd(17,01) + rrt(143) * density(13) 
  pd(17,13) = pd(17,13) + rrt(143) * density(01) 
  pd(01,01) = pd(01,01) + rrt(144) * density(16) 
  pd(01,16) = pd(01,16) + rrt(144) * density(01) 
  pd(16,01) = pd(16,01) - rrt(144) * density(16) 
  pd(16,16) = pd(16,16) - rrt(144) * density(01) 
  pd(17,01) = pd(17,01) + rrt(144) * density(16) 
  pd(17,16) = pd(17,16) + rrt(144) * density(01) 
  pd(01,01) = pd(01,01) + rrt(145) * density(15) 
  pd(01,15) = pd(01,15) + rrt(145) * density(01) 
  pd(15,01) = pd(15,01) - rrt(145) * density(15) 
  pd(15,15) = pd(15,15) - rrt(145) * density(01) 
  pd(17,01) = pd(17,01) + rrt(145) * density(15) 
  pd(17,15) = pd(17,15) + rrt(145) * density(01) 
  pd(01,01) = pd(01,01) + rrt(146) * density(14) 
  pd(01,14) = pd(01,14) + rrt(146) * density(01) 
  pd(14,01) = pd(14,01) - rrt(146) * density(14) 
  pd(14,14) = pd(14,14) - rrt(146) * density(01) 
  pd(17,01) = pd(17,01) + rrt(146) * density(14) 
  pd(17,14) = pd(17,14) + rrt(146) * density(01) 
  pd(01,01) = pd(01,01) + rrt(147) * density(13) 
  pd(01,13) = pd(01,13) + rrt(147) * density(01) 
  pd(03,01) = pd(03,01) + rrt(147) * density(13) 
  pd(03,13) = pd(03,13) + rrt(147) * density(01) 
  pd(04,01) = pd(04,01) + rrt(147) * density(13) 
  pd(04,13) = pd(04,13) + rrt(147) * density(01) 
  pd(13,01) = pd(13,01) - rrt(147) * density(13) 
  pd(13,13) = pd(13,13) - rrt(147) * density(01) 
  pd(01,01) = pd(01,01) + rrt(148) * density(16) 
  pd(01,16) = pd(01,16) + rrt(148) * density(01) 
  pd(03,01) = pd(03,01) + rrt(148) * density(16) 
  pd(03,16) = pd(03,16) + rrt(148) * density(01) 
  pd(04,01) = pd(04,01) + rrt(148) * density(16) 
  pd(04,16) = pd(04,16) + rrt(148) * density(01) 
  pd(16,01) = pd(16,01) - rrt(148) * density(16) 
  pd(16,16) = pd(16,16) - rrt(148) * density(01) 
  pd(01,01) = pd(01,01) + rrt(149) * density(15) 
  pd(01,15) = pd(01,15) + rrt(149) * density(01) 
  pd(03,01) = pd(03,01) + rrt(149) * density(15) 
  pd(03,15) = pd(03,15) + rrt(149) * density(01) 
  pd(04,01) = pd(04,01) + rrt(149) * density(15) 
  pd(04,15) = pd(04,15) + rrt(149) * density(01) 
  pd(15,01) = pd(15,01) - rrt(149) * density(15) 
  pd(15,15) = pd(15,15) - rrt(149) * density(01) 
  pd(01,01) = pd(01,01) + rrt(150) * density(14) 
  pd(01,14) = pd(01,14) + rrt(150) * density(01) 
  pd(03,01) = pd(03,01) + rrt(150) * density(14) 
  pd(03,14) = pd(03,14) + rrt(150) * density(01) 
  pd(04,01) = pd(04,01) + rrt(150) * density(14) 
  pd(04,14) = pd(04,14) + rrt(150) * density(01) 
  pd(14,01) = pd(14,01) - rrt(150) * density(14) 
  pd(14,14) = pd(14,14) - rrt(150) * density(01) 
  pd(34,01) = pd(34,01) + rrt(151) * density(36) 
  pd(34,36) = pd(34,36) + rrt(151) * density(01) 
  pd(36,01) = pd(36,01) - rrt(151) * density(36) 
  pd(36,36) = pd(36,36) - rrt(151) * density(01) 
  pd(43,01) = pd(43,01) + rrt(151) * density(36) 
  pd(43,36) = pd(43,36) + rrt(151) * density(01) 
  pd(34,01) = pd(34,01) + rrt(152) * density(37) 
  pd(34,37) = pd(34,37) + rrt(152) * density(01) 
  pd(37,01) = pd(37,01) - rrt(152) * density(37) 
  pd(37,37) = pd(37,37) - rrt(152) * density(01) 
  pd(43,01) = pd(43,01) + rrt(152) * density(37) 
  pd(43,37) = pd(43,37) + rrt(152) * density(01) 
  pd(34,01) = pd(34,01) + rrt(153) * density(38) 
  pd(34,38) = pd(34,38) + rrt(153) * density(01) 
  pd(38,01) = pd(38,01) - rrt(153) * density(38) 
  pd(38,38) = pd(38,38) - rrt(153) * density(01) 
  pd(43,01) = pd(43,01) + rrt(153) * density(38) 
  pd(43,38) = pd(43,38) + rrt(153) * density(01) 
  pd(31,01) = pd(31,01) + rrt(154) * density(36) 
  pd(31,36) = pd(31,36) + rrt(154) * density(01) 
  pd(36,01) = pd(36,01) - rrt(154) * density(36) 
  pd(36,36) = pd(36,36) - rrt(154) * density(01) 
  pd(44,01) = pd(44,01) + rrt(154) * density(36) 
  pd(44,36) = pd(44,36) + rrt(154) * density(01) 
  pd(31,01) = pd(31,01) + rrt(155) * density(37) 
  pd(31,37) = pd(31,37) + rrt(155) * density(01) 
  pd(37,01) = pd(37,01) - rrt(155) * density(37) 
  pd(37,37) = pd(37,37) - rrt(155) * density(01) 
  pd(44,01) = pd(44,01) + rrt(155) * density(37) 
  pd(44,37) = pd(44,37) + rrt(155) * density(01) 
  pd(31,01) = pd(31,01) + rrt(156) * density(38) 
  pd(31,38) = pd(31,38) + rrt(156) * density(01) 
  pd(38,01) = pd(38,01) - rrt(156) * density(38) 
  pd(38,38) = pd(38,38) - rrt(156) * density(01) 
  pd(44,01) = pd(44,01) + rrt(156) * density(38) 
  pd(44,38) = pd(44,38) + rrt(156) * density(01) 
  pd(31,01) = pd(31,01) + rrt(157) * density(34) 
  pd(31,34) = pd(31,34) + rrt(157) * density(01) 
  pd(34,01) = pd(34,01) - rrt(157) * density(34) 
  pd(34,34) = pd(34,34) - rrt(157) * density(01) 
  pd(43,01) = pd(43,01) + rrt(157) * density(34) 
  pd(43,34) = pd(43,34) + rrt(157) * density(01) 
  pd(01,01) = pd(01,01) + rrt(158) * density(36) 
  pd(01,36) = pd(01,36) + rrt(158) * density(01) 
  pd(36,01) = pd(36,01) - rrt(158) * density(36) 
  pd(36,36) = pd(36,36) - rrt(158) * density(01) 
  pd(39,01) = pd(39,01) + rrt(158) * density(36) 
  pd(39,36) = pd(39,36) + rrt(158) * density(01) 
  pd(01,01) = pd(01,01) + rrt(159) * density(37) 
  pd(01,37) = pd(01,37) + rrt(159) * density(01) 
  pd(37,01) = pd(37,01) - rrt(159) * density(37) 
  pd(37,37) = pd(37,37) - rrt(159) * density(01) 
  pd(39,01) = pd(39,01) + rrt(159) * density(37) 
  pd(39,37) = pd(39,37) + rrt(159) * density(01) 
  pd(01,01) = pd(01,01) + rrt(160) * density(38) 
  pd(01,38) = pd(01,38) + rrt(160) * density(01) 
  pd(38,01) = pd(38,01) - rrt(160) * density(38) 
  pd(38,38) = pd(38,38) - rrt(160) * density(01) 
  pd(39,01) = pd(39,01) + rrt(160) * density(38) 
  pd(39,38) = pd(39,38) + rrt(160) * density(01) 
  pd(01,01) = pd(01,01) + rrt(161) * density(36) 
  pd(01,36) = pd(01,36) + rrt(161) * density(01) 
  pd(35,01) = pd(35,01) + rrt(161) * density(36) 
  pd(35,36) = pd(35,36) + rrt(161) * density(01) 
  pd(36,01) = pd(36,01) - rrt(161) * density(36) 
  pd(36,36) = pd(36,36) - rrt(161) * density(01) 
  pd(43,01) = pd(43,01) + rrt(161) * density(36) 
  pd(43,36) = pd(43,36) + rrt(161) * density(01) 
  pd(01,01) = pd(01,01) + rrt(162) * density(37) 
  pd(01,37) = pd(01,37) + rrt(162) * density(01) 
  pd(35,01) = pd(35,01) + rrt(162) * density(37) 
  pd(35,37) = pd(35,37) + rrt(162) * density(01) 
  pd(37,01) = pd(37,01) - rrt(162) * density(37) 
  pd(37,37) = pd(37,37) - rrt(162) * density(01) 
  pd(43,01) = pd(43,01) + rrt(162) * density(37) 
  pd(43,37) = pd(43,37) + rrt(162) * density(01) 
  pd(01,01) = pd(01,01) + rrt(163) * density(38) 
  pd(01,38) = pd(01,38) + rrt(163) * density(01) 
  pd(35,01) = pd(35,01) + rrt(163) * density(38) 
  pd(35,38) = pd(35,38) + rrt(163) * density(01) 
  pd(38,01) = pd(38,01) - rrt(163) * density(38) 
  pd(38,38) = pd(38,38) - rrt(163) * density(01) 
  pd(43,01) = pd(43,01) + rrt(163) * density(38) 
  pd(43,38) = pd(43,38) + rrt(163) * density(01) 
  pd(01,01) = pd(01,01) + rrt(164) * density(36) 
  pd(01,36) = pd(01,36) + rrt(164) * density(01) 
  pd(33,01) = pd(33,01) + rrt(164) * density(36) 
  pd(33,36) = pd(33,36) + rrt(164) * density(01) 
  pd(36,01) = pd(36,01) - rrt(164) * density(36) 
  pd(36,36) = pd(36,36) - rrt(164) * density(01) 
  pd(44,01) = pd(44,01) + rrt(164) * density(36) 
  pd(44,36) = pd(44,36) + rrt(164) * density(01) 
  pd(01,01) = pd(01,01) + rrt(165) * density(37) 
  pd(01,37) = pd(01,37) + rrt(165) * density(01) 
  pd(33,01) = pd(33,01) + rrt(165) * density(37) 
  pd(33,37) = pd(33,37) + rrt(165) * density(01) 
  pd(37,01) = pd(37,01) - rrt(165) * density(37) 
  pd(37,37) = pd(37,37) - rrt(165) * density(01) 
  pd(44,01) = pd(44,01) + rrt(165) * density(37) 
  pd(44,37) = pd(44,37) + rrt(165) * density(01) 
  pd(01,01) = pd(01,01) + rrt(166) * density(38) 
  pd(01,38) = pd(01,38) + rrt(166) * density(01) 
  pd(33,01) = pd(33,01) + rrt(166) * density(38) 
  pd(33,38) = pd(33,38) + rrt(166) * density(01) 
  pd(38,01) = pd(38,01) - rrt(166) * density(38) 
  pd(38,38) = pd(38,38) - rrt(166) * density(01) 
  pd(44,01) = pd(44,01) + rrt(166) * density(38) 
  pd(44,38) = pd(44,38) + rrt(166) * density(01) 
  pd(01,01) = pd(01,01) + rrt(167) * density(34) 
  pd(01,34) = pd(01,34) + rrt(167) * density(01) 
  pd(34,01) = pd(34,01) - rrt(167) * density(34) 
  pd(34,34) = pd(34,34) - rrt(167) * density(01) 
  pd(35,01) = pd(35,01) + rrt(167) * density(34) 
  pd(35,34) = pd(35,34) + rrt(167) * density(01) 
  pd(01,01) = pd(01,01) + rrt(168) * density(34) 
  pd(01,34) = pd(01,34) + rrt(168) * density(01) 
  pd(33,01) = pd(33,01) + rrt(168) * density(34) 
  pd(33,34) = pd(33,34) + rrt(168) * density(01) 
  pd(34,01) = pd(34,01) - rrt(168) * density(34) 
  pd(34,34) = pd(34,34) - rrt(168) * density(01) 
  pd(43,01) = pd(43,01) + rrt(168) * density(34) 
  pd(43,34) = pd(43,34) + rrt(168) * density(01) 
  pd(01,01) = pd(01,01) + rrt(169) * density(31) 
  pd(01,31) = pd(01,31) + rrt(169) * density(01) 
  pd(31,01) = pd(31,01) - rrt(169) * density(31) 
  pd(31,31) = pd(31,31) - rrt(169) * density(01) 
  pd(33,01) = pd(33,01) + rrt(169) * density(31) 
  pd(33,31) = pd(33,31) + rrt(169) * density(01) 
  pd(01,01) = pd(01,01) + rrt(170) * density(32) 
  pd(01,32) = pd(01,32) + rrt(170) * density(01) 
  pd(32,01) = pd(32,01) - rrt(170) * density(32) 
  pd(32,32) = pd(32,32) - rrt(170) * density(01) 
  pd(33,01) = pd(33,01) + rrt(170) * density(32) 
  pd(33,32) = pd(33,32) + rrt(170) * density(01) 
  pd(09,01) = pd(09,01) + rrt(171) * density(31) 
  pd(09,31) = pd(09,31) + rrt(171) * density(01) 
  pd(13,01) = pd(13,01) + rrt(171) * density(31) 
  pd(13,31) = pd(13,31) + rrt(171) * density(01) 
  pd(31,01) = pd(31,01) - rrt(171) * density(31) 
  pd(31,31) = pd(31,31) - rrt(171) * density(01) 
  pd(09,01) = pd(09,01) + rrt(172) * density(32) 
  pd(09,32) = pd(09,32) + rrt(172) * density(01) 
  pd(13,01) = pd(13,01) + rrt(172) * density(32) 
  pd(13,32) = pd(13,32) + rrt(172) * density(01) 
  pd(32,01) = pd(32,01) - rrt(172) * density(32) 
  pd(32,32) = pd(32,32) - rrt(172) * density(01) 
  pd(01,01) = pd(01,01) + rrt(173) * density(31) 
  pd(01,31) = pd(01,31) + rrt(173) * density(01) 
  pd(12,01) = pd(12,01) + rrt(173) * density(31) 
  pd(12,31) = pd(12,31) + rrt(173) * density(01) 
  pd(13,01) = pd(13,01) + rrt(173) * density(31) 
  pd(13,31) = pd(13,31) + rrt(173) * density(01) 
  pd(31,01) = pd(31,01) - rrt(173) * density(31) 
  pd(31,31) = pd(31,31) - rrt(173) * density(01) 
  pd(01,01) = pd(01,01) + rrt(174) * density(32) 
  pd(01,32) = pd(01,32) + rrt(174) * density(01) 
  pd(12,01) = pd(12,01) + rrt(174) * density(32) 
  pd(12,32) = pd(12,32) + rrt(174) * density(01) 
  pd(13,01) = pd(13,01) + rrt(174) * density(32) 
  pd(13,32) = pd(13,32) + rrt(174) * density(01) 
  pd(32,01) = pd(32,01) - rrt(174) * density(32) 
  pd(32,32) = pd(32,32) - rrt(174) * density(01) 
  pd(01,01) = pd(01,01) - rrt(175) * density(01) * density(04) * 2.0d0
  pd(01,04) = pd(01,04) - rrt(175) * density(01)**2 
  pd(03,01) = pd(03,01) + rrt(175) * density(01) * density(04) * 2.0d0
  pd(03,04) = pd(03,04) + rrt(175) * density(01)**2 
  pd(04,01) = pd(04,01) - rrt(175) * density(01) * density(04) * 2.0d0
  pd(04,04) = pd(04,04) - rrt(175) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(176) * density(01) * density(06) * 2.0d0
  pd(01,06) = pd(01,06) - rrt(176) * density(01)**2 
  pd(05,01) = pd(05,01) + rrt(176) * density(01) * density(06) * 2.0d0
  pd(05,06) = pd(05,06) + rrt(176) * density(01)**2 
  pd(06,01) = pd(06,01) - rrt(176) * density(01) * density(06) * 2.0d0
  pd(06,06) = pd(06,06) - rrt(176) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(177) * density(01) * density(08) * 2.0d0
  pd(01,08) = pd(01,08) - rrt(177) * density(01)**2 
  pd(07,01) = pd(07,01) + rrt(177) * density(01) * density(08) * 2.0d0
  pd(07,08) = pd(07,08) + rrt(177) * density(01)**2 
  pd(08,01) = pd(08,01) - rrt(177) * density(01) * density(08) * 2.0d0
  pd(08,08) = pd(08,08) - rrt(177) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(178) * density(01) * density(12) * 2.0d0
  pd(01,12) = pd(01,12) - rrt(178) * density(01)**2 
  pd(09,01) = pd(09,01) + rrt(178) * density(01) * density(12) * 2.0d0
  pd(09,12) = pd(09,12) + rrt(178) * density(01)**2 
  pd(12,01) = pd(12,01) - rrt(178) * density(01) * density(12) * 2.0d0
  pd(12,12) = pd(12,12) - rrt(178) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(179) * density(01) * density(17) * 2.0d0
  pd(01,17) = pd(01,17) - rrt(179) * density(01)**2 
  pd(13,01) = pd(13,01) + rrt(179) * density(01) * density(17) * 2.0d0
  pd(13,17) = pd(13,17) + rrt(179) * density(01)**2 
  pd(17,01) = pd(17,01) - rrt(179) * density(01) * density(17) * 2.0d0
  pd(17,17) = pd(17,17) - rrt(179) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(180) * density(01) * density(19) * 2.0d0
  pd(01,19) = pd(01,19) - rrt(180) * density(01)**2 
  pd(18,01) = pd(18,01) + rrt(180) * density(01) * density(19) * 2.0d0
  pd(18,19) = pd(18,19) + rrt(180) * density(01)**2 
  pd(19,01) = pd(19,01) - rrt(180) * density(01) * density(19) * 2.0d0
  pd(19,19) = pd(19,19) - rrt(180) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(181) * density(01) * density(23) * 2.0d0
  pd(01,23) = pd(01,23) - rrt(181) * density(01)**2 
  pd(20,01) = pd(20,01) + rrt(181) * density(01) * density(23) * 2.0d0
  pd(20,23) = pd(20,23) + rrt(181) * density(01)**2 
  pd(23,01) = pd(23,01) - rrt(181) * density(01) * density(23) * 2.0d0
  pd(23,23) = pd(23,23) - rrt(181) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(182) * density(01) * density(25) * 2.0d0
  pd(01,25) = pd(01,25) - rrt(182) * density(01)**2 
  pd(24,01) = pd(24,01) + rrt(182) * density(01) * density(25) * 2.0d0
  pd(24,25) = pd(24,25) + rrt(182) * density(01)**2 
  pd(25,01) = pd(25,01) - rrt(182) * density(01) * density(25) * 2.0d0
  pd(25,25) = pd(25,25) - rrt(182) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(183) * density(01) * density(29) * 2.0d0
  pd(01,29) = pd(01,29) - rrt(183) * density(01)**2 
  pd(26,01) = pd(26,01) + rrt(183) * density(01) * density(29) * 2.0d0
  pd(26,29) = pd(26,29) + rrt(183) * density(01)**2 
  pd(29,01) = pd(29,01) - rrt(183) * density(01) * density(29) * 2.0d0
  pd(29,29) = pd(29,29) - rrt(183) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(184) * density(01) * density(33) * 2.0d0
  pd(01,33) = pd(01,33) - rrt(184) * density(01)**2 
  pd(31,01) = pd(31,01) + rrt(184) * density(01) * density(33) * 2.0d0
  pd(31,33) = pd(31,33) + rrt(184) * density(01)**2 
  pd(33,01) = pd(33,01) - rrt(184) * density(01) * density(33) * 2.0d0
  pd(33,33) = pd(33,33) - rrt(184) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(185) * density(01) * density(35) * 2.0d0
  pd(01,35) = pd(01,35) - rrt(185) * density(01)**2 
  pd(34,01) = pd(34,01) + rrt(185) * density(01) * density(35) * 2.0d0
  pd(34,35) = pd(34,35) + rrt(185) * density(01)**2 
  pd(35,01) = pd(35,01) - rrt(185) * density(01) * density(35) * 2.0d0
  pd(35,35) = pd(35,35) - rrt(185) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(186) * density(01) * density(39) * 2.0d0
  pd(01,39) = pd(01,39) - rrt(186) * density(01)**2 
  pd(36,01) = pd(36,01) + rrt(186) * density(01) * density(39) * 2.0d0
  pd(36,39) = pd(36,39) + rrt(186) * density(01)**2 
  pd(39,01) = pd(39,01) - rrt(186) * density(01) * density(39) * 2.0d0
  pd(39,39) = pd(39,39) - rrt(186) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(187) * density(12) 
  pd(01,12) = pd(01,12) - rrt(187) * density(01) 
  pd(07,01) = pd(07,01) + rrt(187) * density(12) 
  pd(07,12) = pd(07,12) + rrt(187) * density(01) 
  pd(12,01) = pd(12,01) - rrt(187) * density(12) 
  pd(12,12) = pd(12,12) - rrt(187) * density(01) 
  pd(43,01) = pd(43,01) + rrt(187) * density(12) 
  pd(43,12) = pd(43,12) + rrt(187) * density(01) 
  pd(01,01) = pd(01,01) - rrt(188) * density(12) 
  pd(01,12) = pd(01,12) - rrt(188) * density(01) 
  pd(05,01) = pd(05,01) + rrt(188) * density(12) 
  pd(05,12) = pd(05,12) + rrt(188) * density(01) 
  pd(12,01) = pd(12,01) - rrt(188) * density(12) 
  pd(12,12) = pd(12,12) - rrt(188) * density(01) 
  pd(43,01) = pd(43,01) + rrt(188) * density(12) * 2.0d0
  pd(43,12) = pd(43,12) + rrt(188) * density(01) * 2.0d0
  pd(01,01) = pd(01,01) - rrt(189) * density(12) 
  pd(01,12) = pd(01,12) - rrt(189) * density(01) 
  pd(03,01) = pd(03,01) + rrt(189) * density(12) 
  pd(03,12) = pd(03,12) + rrt(189) * density(01) 
  pd(12,01) = pd(12,01) - rrt(189) * density(12) 
  pd(12,12) = pd(12,12) - rrt(189) * density(01) 
  pd(43,01) = pd(43,01) + rrt(189) * density(12) 
  pd(43,12) = pd(43,12) + rrt(189) * density(01) 
  pd(44,01) = pd(44,01) + rrt(189) * density(12) 
  pd(44,12) = pd(44,12) + rrt(189) * density(01) 
  pd(01,01) = pd(01,01) - rrt(190) * density(08) 
  pd(01,08) = pd(01,08) - rrt(190) * density(01) 
  pd(05,01) = pd(05,01) + rrt(190) * density(08) 
  pd(05,08) = pd(05,08) + rrt(190) * density(01) 
  pd(08,01) = pd(08,01) - rrt(190) * density(08) 
  pd(08,08) = pd(08,08) - rrt(190) * density(01) 
  pd(43,01) = pd(43,01) + rrt(190) * density(08) 
  pd(43,08) = pd(43,08) + rrt(190) * density(01) 
  pd(01,01) = pd(01,01) - rrt(191) * density(08) 
  pd(01,08) = pd(01,08) - rrt(191) * density(01) 
  pd(03,01) = pd(03,01) + rrt(191) * density(08) 
  pd(03,08) = pd(03,08) + rrt(191) * density(01) 
  pd(08,01) = pd(08,01) - rrt(191) * density(08) 
  pd(08,08) = pd(08,08) - rrt(191) * density(01) 
  pd(44,01) = pd(44,01) + rrt(191) * density(08) 
  pd(44,08) = pd(44,08) + rrt(191) * density(01) 
  pd(01,01) = pd(01,01) - rrt(192) * density(08) 
  pd(01,08) = pd(01,08) - rrt(192) * density(01) 
  pd(03,01) = pd(03,01) + rrt(192) * density(08) 
  pd(03,08) = pd(03,08) + rrt(192) * density(01) 
  pd(08,01) = pd(08,01) - rrt(192) * density(08) 
  pd(08,08) = pd(08,08) - rrt(192) * density(01) 
  pd(43,01) = pd(43,01) + rrt(192) * density(08) * 2.0d0
  pd(43,08) = pd(43,08) + rrt(192) * density(01) * 2.0d0
  pd(01,01) = pd(01,01) - rrt(193) * density(08) 
  pd(01,08) = pd(01,08) - rrt(193) * density(01) 
  pd(02,01) = pd(02,01) + rrt(193) * density(08) 
  pd(02,08) = pd(02,08) + rrt(193) * density(01) 
  pd(08,01) = pd(08,01) - rrt(193) * density(08) 
  pd(08,08) = pd(08,08) - rrt(193) * density(01) 
  pd(43,01) = pd(43,01) + rrt(193) * density(08) 
  pd(43,08) = pd(43,08) + rrt(193) * density(01) 
  pd(44,01) = pd(44,01) + rrt(193) * density(08) 
  pd(44,08) = pd(44,08) + rrt(193) * density(01) 
  pd(01,01) = pd(01,01) - rrt(194) * density(06) 
  pd(01,06) = pd(01,06) - rrt(194) * density(01) 
  pd(03,01) = pd(03,01) + rrt(194) * density(06) 
  pd(03,06) = pd(03,06) + rrt(194) * density(01) 
  pd(06,01) = pd(06,01) - rrt(194) * density(06) 
  pd(06,06) = pd(06,06) - rrt(194) * density(01) 
  pd(43,01) = pd(43,01) + rrt(194) * density(06) 
  pd(43,06) = pd(43,06) + rrt(194) * density(01) 
  pd(01,01) = pd(01,01) - rrt(195) * density(06) 
  pd(01,06) = pd(01,06) - rrt(195) * density(01) 
  pd(02,01) = pd(02,01) + rrt(195) * density(06) 
  pd(02,06) = pd(02,06) + rrt(195) * density(01) 
  pd(06,01) = pd(06,01) - rrt(195) * density(06) 
  pd(06,06) = pd(06,06) - rrt(195) * density(01) 
  pd(44,01) = pd(44,01) + rrt(195) * density(06) 
  pd(44,06) = pd(44,06) + rrt(195) * density(01) 
  pd(01,01) = pd(01,01) - rrt(196) * density(06) 
  pd(01,06) = pd(01,06) - rrt(196) * density(01) 
  pd(02,01) = pd(02,01) + rrt(196) * density(06) 
  pd(02,06) = pd(02,06) + rrt(196) * density(01) 
  pd(06,01) = pd(06,01) - rrt(196) * density(06) 
  pd(06,06) = pd(06,06) - rrt(196) * density(01) 
  pd(43,01) = pd(43,01) + rrt(196) * density(06) * 2.0d0
  pd(43,06) = pd(43,06) + rrt(196) * density(01) * 2.0d0
  pd(01,01) = pd(01,01) - rrt(197) * density(04) 
  pd(01,04) = pd(01,04) - rrt(197) * density(01) 
  pd(02,01) = pd(02,01) + rrt(197) * density(04) 
  pd(02,04) = pd(02,04) + rrt(197) * density(01) 
  pd(04,01) = pd(04,01) - rrt(197) * density(04) 
  pd(04,04) = pd(04,04) - rrt(197) * density(01) 
  pd(43,01) = pd(43,01) + rrt(197) * density(04) 
  pd(43,04) = pd(43,04) + rrt(197) * density(01) 
  pd(01,01) = pd(01,01) - rrt(198) * density(29) 
  pd(01,29) = pd(01,29) - rrt(198) * density(01) 
  pd(24,01) = pd(24,01) + rrt(198) * density(29) 
  pd(24,29) = pd(24,29) + rrt(198) * density(01) 
  pd(29,01) = pd(29,01) - rrt(198) * density(29) 
  pd(29,29) = pd(29,29) - rrt(198) * density(01) 
  pd(43,01) = pd(43,01) + rrt(198) * density(29) 
  pd(43,29) = pd(43,29) + rrt(198) * density(01) 
  pd(01,01) = pd(01,01) - rrt(199) * density(29) 
  pd(01,29) = pd(01,29) - rrt(199) * density(01) 
  pd(20,01) = pd(20,01) + rrt(199) * density(29) 
  pd(20,29) = pd(20,29) + rrt(199) * density(01) 
  pd(29,01) = pd(29,01) - rrt(199) * density(29) 
  pd(29,29) = pd(29,29) - rrt(199) * density(01) 
  pd(43,01) = pd(43,01) + rrt(199) * density(29) * 2.0d0
  pd(43,29) = pd(43,29) + rrt(199) * density(01) * 2.0d0
  pd(01,01) = pd(01,01) - rrt(200) * density(25) 
  pd(01,25) = pd(01,25) - rrt(200) * density(01) 
  pd(20,01) = pd(20,01) + rrt(200) * density(25) 
  pd(20,25) = pd(20,25) + rrt(200) * density(01) 
  pd(25,01) = pd(25,01) - rrt(200) * density(25) 
  pd(25,25) = pd(25,25) - rrt(200) * density(01) 
  pd(43,01) = pd(43,01) + rrt(200) * density(25) 
  pd(43,25) = pd(43,25) + rrt(200) * density(01) 
  pd(01,01) = pd(01,01) - rrt(201) * density(25) 
  pd(01,25) = pd(01,25) - rrt(201) * density(01) 
  pd(18,01) = pd(18,01) + rrt(201) * density(25) 
  pd(18,25) = pd(18,25) + rrt(201) * density(01) 
  pd(25,01) = pd(25,01) - rrt(201) * density(25) 
  pd(25,25) = pd(25,25) - rrt(201) * density(01) 
  pd(43,01) = pd(43,01) + rrt(201) * density(25) * 2.0d0
  pd(43,25) = pd(43,25) + rrt(201) * density(01) * 2.0d0
  pd(01,01) = pd(01,01) - rrt(202) * density(25) 
  pd(01,25) = pd(01,25) - rrt(202) * density(01) 
  pd(13,01) = pd(13,01) + rrt(202) * density(25) 
  pd(13,25) = pd(13,25) + rrt(202) * density(01) 
  pd(25,01) = pd(25,01) - rrt(202) * density(25) 
  pd(25,25) = pd(25,25) - rrt(202) * density(01) 
  pd(43,01) = pd(43,01) + rrt(202) * density(25) 
  pd(43,25) = pd(43,25) + rrt(202) * density(01) 
  pd(44,01) = pd(44,01) + rrt(202) * density(25) 
  pd(44,25) = pd(44,25) + rrt(202) * density(01) 
  pd(01,01) = pd(01,01) - rrt(203) * density(25) 
  pd(01,25) = pd(01,25) - rrt(203) * density(01) 
  pd(13,01) = pd(13,01) + rrt(203) * density(25) 
  pd(13,25) = pd(13,25) + rrt(203) * density(01) 
  pd(25,01) = pd(25,01) - rrt(203) * density(25) 
  pd(25,25) = pd(25,25) - rrt(203) * density(01) 
  pd(43,01) = pd(43,01) + rrt(203) * density(25) * 3.0d0
  pd(43,25) = pd(43,25) + rrt(203) * density(01) * 3.0d0
  pd(01,01) = pd(01,01) - rrt(204) * density(25) 
  pd(01,25) = pd(01,25) - rrt(204) * density(01) 
  pd(05,01) = pd(05,01) + rrt(204) * density(25) 
  pd(05,25) = pd(05,25) + rrt(204) * density(01) 
  pd(07,01) = pd(07,01) + rrt(204) * density(25) 
  pd(07,25) = pd(07,25) + rrt(204) * density(01) 
  pd(25,01) = pd(25,01) - rrt(204) * density(25) 
  pd(25,25) = pd(25,25) - rrt(204) * density(01) 
  pd(01,01) = pd(01,01) - rrt(205) * density(23) 
  pd(01,23) = pd(01,23) - rrt(205) * density(01) 
  pd(18,01) = pd(18,01) + rrt(205) * density(23) 
  pd(18,23) = pd(18,23) + rrt(205) * density(01) 
  pd(23,01) = pd(23,01) - rrt(205) * density(23) 
  pd(23,23) = pd(23,23) - rrt(205) * density(01) 
  pd(43,01) = pd(43,01) + rrt(205) * density(23) 
  pd(43,23) = pd(43,23) + rrt(205) * density(01) 
  pd(01,01) = pd(01,01) - rrt(206) * density(23) 
  pd(01,23) = pd(01,23) - rrt(206) * density(01) 
  pd(13,01) = pd(13,01) + rrt(206) * density(23) 
  pd(13,23) = pd(13,23) + rrt(206) * density(01) 
  pd(23,01) = pd(23,01) - rrt(206) * density(23) 
  pd(23,23) = pd(23,23) - rrt(206) * density(01) 
  pd(43,01) = pd(43,01) + rrt(206) * density(23) * 2.0d0
  pd(43,23) = pd(43,23) + rrt(206) * density(01) * 2.0d0
  pd(01,01) = pd(01,01) - rrt(207) * density(19) 
  pd(01,19) = pd(01,19) - rrt(207) * density(01) 
  pd(13,01) = pd(13,01) + rrt(207) * density(19) 
  pd(13,19) = pd(13,19) + rrt(207) * density(01) 
  pd(19,01) = pd(19,01) - rrt(207) * density(19) 
  pd(19,19) = pd(19,19) - rrt(207) * density(01) 
  pd(43,01) = pd(43,01) + rrt(207) * density(19) 
  pd(43,19) = pd(43,19) + rrt(207) * density(01) 
  pd(01,01) = pd(01,01) - rrt(208) * density(17) 
  pd(01,17) = pd(01,17) - rrt(208) * density(01) 
  pd(03,01) = pd(03,01) + rrt(208) * density(17) * 2.0d0
  pd(03,17) = pd(03,17) + rrt(208) * density(01) * 2.0d0
  pd(17,01) = pd(17,01) - rrt(208) * density(17) 
  pd(17,17) = pd(17,17) - rrt(208) * density(01) 
  pd(01,01) = pd(01,01) - rrt(209) * density(45) 
  pd(01,45) = pd(01,45) - rrt(209) * density(01) 
  pd(43,01) = pd(43,01) + rrt(209) * density(45) * 2.0d0
  pd(43,45) = pd(43,45) + rrt(209) * density(01) * 2.0d0
  pd(45,01) = pd(45,01) - rrt(209) * density(45) 
  pd(45,45) = pd(45,45) - rrt(209) * density(01) 
  pd(01,01) = pd(01,01) - rrt(210) * density(46) 
  pd(01,46) = pd(01,46) - rrt(210) * density(01) 
  pd(43,01) = pd(43,01) + rrt(210) * density(46) 
  pd(43,46) = pd(43,46) + rrt(210) * density(01) 
  pd(44,01) = pd(44,01) + rrt(210) * density(46) 
  pd(44,46) = pd(44,46) + rrt(210) * density(01) 
  pd(46,01) = pd(46,01) - rrt(210) * density(46) 
  pd(46,46) = pd(46,46) - rrt(210) * density(01) 
  pd(43,44) = pd(43,44) + rrt(211) * density(45) 
  pd(43,45) = pd(43,45) + rrt(211) * density(44) 
  pd(44,44) = pd(44,44) - rrt(211) * density(45) 
  pd(44,45) = pd(44,45) - rrt(211) * density(44) 
  pd(45,44) = pd(45,44) - rrt(211) * density(45) 
  pd(45,45) = pd(45,45) - rrt(211) * density(44) 
  pd(46,44) = pd(46,44) + rrt(211) * density(45) 
  pd(46,45) = pd(46,45) + rrt(211) * density(44) 
  pd(09,12) = pd(09,12) + rrt(212) * density(26) 
  pd(09,26) = pd(09,26) + rrt(212) * density(12) 
  pd(12,12) = pd(12,12) - rrt(212) * density(26) 
  pd(12,26) = pd(12,26) - rrt(212) * density(12) 
  pd(23,12) = pd(23,12) + rrt(212) * density(26) 
  pd(23,26) = pd(23,26) + rrt(212) * density(12) 
  pd(26,12) = pd(26,12) - rrt(212) * density(26) 
  pd(26,26) = pd(26,26) - rrt(212) * density(12) 
  pd(44,12) = pd(44,12) + rrt(212) * density(26) 
  pd(44,26) = pd(44,26) + rrt(212) * density(12) 
  pd(09,12) = pd(09,12) + rrt(213) * density(28) 
  pd(09,28) = pd(09,28) + rrt(213) * density(12) 
  pd(12,12) = pd(12,12) - rrt(213) * density(28) 
  pd(12,28) = pd(12,28) - rrt(213) * density(12) 
  pd(23,12) = pd(23,12) + rrt(213) * density(28) 
  pd(23,28) = pd(23,28) + rrt(213) * density(12) 
  pd(28,12) = pd(28,12) - rrt(213) * density(28) 
  pd(28,28) = pd(28,28) - rrt(213) * density(12) 
  pd(44,12) = pd(44,12) + rrt(213) * density(28) 
  pd(44,28) = pd(44,28) + rrt(213) * density(12) 
  pd(09,12) = pd(09,12) + rrt(214) * density(27) 
  pd(09,27) = pd(09,27) + rrt(214) * density(12) 
  pd(12,12) = pd(12,12) - rrt(214) * density(27) 
  pd(12,27) = pd(12,27) - rrt(214) * density(12) 
  pd(23,12) = pd(23,12) + rrt(214) * density(27) 
  pd(23,27) = pd(23,27) + rrt(214) * density(12) 
  pd(27,12) = pd(27,12) - rrt(214) * density(27) 
  pd(27,27) = pd(27,27) - rrt(214) * density(12) 
  pd(44,12) = pd(44,12) + rrt(214) * density(27) 
  pd(44,27) = pd(44,27) + rrt(214) * density(12) 
  pd(07,12) = pd(07,12) + rrt(215) * density(20) 
  pd(07,20) = pd(07,20) + rrt(215) * density(12) 
  pd(12,12) = pd(12,12) - rrt(215) * density(20) 
  pd(12,20) = pd(12,20) - rrt(215) * density(12) 
  pd(20,12) = pd(20,12) - rrt(215) * density(20) 
  pd(20,20) = pd(20,20) - rrt(215) * density(12) 
  pd(25,12) = pd(25,12) + rrt(215) * density(20) 
  pd(25,20) = pd(25,20) + rrt(215) * density(12) 
  pd(07,12) = pd(07,12) + rrt(216) * density(21) 
  pd(07,21) = pd(07,21) + rrt(216) * density(12) 
  pd(12,12) = pd(12,12) - rrt(216) * density(21) 
  pd(12,21) = pd(12,21) - rrt(216) * density(12) 
  pd(21,12) = pd(21,12) - rrt(216) * density(21) 
  pd(21,21) = pd(21,21) - rrt(216) * density(12) 
  pd(25,12) = pd(25,12) + rrt(216) * density(21) 
  pd(25,21) = pd(25,21) + rrt(216) * density(12) 
  pd(07,12) = pd(07,12) + rrt(217) * density(22) 
  pd(07,22) = pd(07,22) + rrt(217) * density(12) 
  pd(12,12) = pd(12,12) - rrt(217) * density(22) 
  pd(12,22) = pd(12,22) - rrt(217) * density(12) 
  pd(22,12) = pd(22,12) - rrt(217) * density(22) 
  pd(22,22) = pd(22,22) - rrt(217) * density(12) 
  pd(25,12) = pd(25,12) + rrt(217) * density(22) 
  pd(25,22) = pd(25,22) + rrt(217) * density(12) 
  pd(09,12) = pd(09,12) + rrt(218) * density(20) 
  pd(09,20) = pd(09,20) + rrt(218) * density(12) 
  pd(12,12) = pd(12,12) - rrt(218) * density(20) 
  pd(12,20) = pd(12,20) - rrt(218) * density(12) 
  pd(20,12) = pd(20,12) - rrt(218) * density(20) 
  pd(20,20) = pd(20,20) - rrt(218) * density(12) 
  pd(23,12) = pd(23,12) + rrt(218) * density(20) 
  pd(23,20) = pd(23,20) + rrt(218) * density(12) 
  pd(09,12) = pd(09,12) + rrt(219) * density(21) 
  pd(09,21) = pd(09,21) + rrt(219) * density(12) 
  pd(12,12) = pd(12,12) - rrt(219) * density(21) 
  pd(12,21) = pd(12,21) - rrt(219) * density(12) 
  pd(21,12) = pd(21,12) - rrt(219) * density(21) 
  pd(21,21) = pd(21,21) - rrt(219) * density(12) 
  pd(23,12) = pd(23,12) + rrt(219) * density(21) 
  pd(23,21) = pd(23,21) + rrt(219) * density(12) 
  pd(09,12) = pd(09,12) + rrt(220) * density(22) 
  pd(09,22) = pd(09,22) + rrt(220) * density(12) 
  pd(12,12) = pd(12,12) - rrt(220) * density(22) 
  pd(12,22) = pd(12,22) - rrt(220) * density(12) 
  pd(22,12) = pd(22,12) - rrt(220) * density(22) 
  pd(22,22) = pd(22,22) - rrt(220) * density(12) 
  pd(23,12) = pd(23,12) + rrt(220) * density(22) 
  pd(23,22) = pd(23,22) + rrt(220) * density(12) 
  pd(07,12) = pd(07,12) + rrt(221) * density(13) 
  pd(07,13) = pd(07,13) + rrt(221) * density(12) 
  pd(12,12) = pd(12,12) - rrt(221) * density(13) 
  pd(12,13) = pd(12,13) - rrt(221) * density(12) 
  pd(13,12) = pd(13,12) - rrt(221) * density(13) 
  pd(13,13) = pd(13,13) - rrt(221) * density(12) 
  pd(19,12) = pd(19,12) + rrt(221) * density(13) 
  pd(19,13) = pd(19,13) + rrt(221) * density(12) 
  pd(07,12) = pd(07,12) + rrt(222) * density(16) 
  pd(07,16) = pd(07,16) + rrt(222) * density(12) 
  pd(12,12) = pd(12,12) - rrt(222) * density(16) 
  pd(12,16) = pd(12,16) - rrt(222) * density(12) 
  pd(16,12) = pd(16,12) - rrt(222) * density(16) 
  pd(16,16) = pd(16,16) - rrt(222) * density(12) 
  pd(19,12) = pd(19,12) + rrt(222) * density(16) 
  pd(19,16) = pd(19,16) + rrt(222) * density(12) 
  pd(07,12) = pd(07,12) + rrt(223) * density(15) 
  pd(07,15) = pd(07,15) + rrt(223) * density(12) 
  pd(12,12) = pd(12,12) - rrt(223) * density(15) 
  pd(12,15) = pd(12,15) - rrt(223) * density(12) 
  pd(15,12) = pd(15,12) - rrt(223) * density(15) 
  pd(15,15) = pd(15,15) - rrt(223) * density(12) 
  pd(19,12) = pd(19,12) + rrt(223) * density(15) 
  pd(19,15) = pd(19,15) + rrt(223) * density(12) 
  pd(07,12) = pd(07,12) + rrt(224) * density(14) 
  pd(07,14) = pd(07,14) + rrt(224) * density(12) 
  pd(12,12) = pd(12,12) - rrt(224) * density(14) 
  pd(12,14) = pd(12,14) - rrt(224) * density(12) 
  pd(14,12) = pd(14,12) - rrt(224) * density(14) 
  pd(14,14) = pd(14,14) - rrt(224) * density(12) 
  pd(19,12) = pd(19,12) + rrt(224) * density(14) 
  pd(19,14) = pd(19,14) + rrt(224) * density(12) 
  pd(09,12) = pd(09,12) + rrt(225) * density(13) 
  pd(09,13) = pd(09,13) + rrt(225) * density(12) 
  pd(12,12) = pd(12,12) - rrt(225) * density(13) 
  pd(12,13) = pd(12,13) - rrt(225) * density(12) 
  pd(13,12) = pd(13,12) - rrt(225) * density(13) 
  pd(13,13) = pd(13,13) - rrt(225) * density(12) 
  pd(17,12) = pd(17,12) + rrt(225) * density(13) 
  pd(17,13) = pd(17,13) + rrt(225) * density(12) 
  pd(09,12) = pd(09,12) + rrt(226) * density(16) 
  pd(09,16) = pd(09,16) + rrt(226) * density(12) 
  pd(12,12) = pd(12,12) - rrt(226) * density(16) 
  pd(12,16) = pd(12,16) - rrt(226) * density(12) 
  pd(16,12) = pd(16,12) - rrt(226) * density(16) 
  pd(16,16) = pd(16,16) - rrt(226) * density(12) 
  pd(17,12) = pd(17,12) + rrt(226) * density(16) 
  pd(17,16) = pd(17,16) + rrt(226) * density(12) 
  pd(09,12) = pd(09,12) + rrt(227) * density(15) 
  pd(09,15) = pd(09,15) + rrt(227) * density(12) 
  pd(12,12) = pd(12,12) - rrt(227) * density(15) 
  pd(12,15) = pd(12,15) - rrt(227) * density(12) 
  pd(15,12) = pd(15,12) - rrt(227) * density(15) 
  pd(15,15) = pd(15,15) - rrt(227) * density(12) 
  pd(17,12) = pd(17,12) + rrt(227) * density(15) 
  pd(17,15) = pd(17,15) + rrt(227) * density(12) 
  pd(09,12) = pd(09,12) + rrt(228) * density(14) 
  pd(09,14) = pd(09,14) + rrt(228) * density(12) 
  pd(12,12) = pd(12,12) - rrt(228) * density(14) 
  pd(12,14) = pd(12,14) - rrt(228) * density(12) 
  pd(14,12) = pd(14,12) - rrt(228) * density(14) 
  pd(14,14) = pd(14,14) - rrt(228) * density(12) 
  pd(17,12) = pd(17,12) + rrt(228) * density(14) 
  pd(17,14) = pd(17,14) + rrt(228) * density(12) 
  pd(08,12) = pd(08,12) + rrt(229) * density(43) 
  pd(08,43) = pd(08,43) + rrt(229) * density(12) 
  pd(12,12) = pd(12,12) - rrt(229) * density(43) 
  pd(12,43) = pd(12,43) - rrt(229) * density(12) 
  pd(43,12) = pd(43,12) - rrt(229) * density(43) 
  pd(43,43) = pd(43,43) - rrt(229) * density(12) 
  pd(44,12) = pd(44,12) + rrt(229) * density(43) 
  pd(44,43) = pd(44,43) + rrt(229) * density(12) 
  pd(07,08) = pd(07,08) + rrt(230) * density(09) 
  pd(07,09) = pd(07,09) + rrt(230) * density(08) 
  pd(08,08) = pd(08,08) - rrt(230) * density(09) 
  pd(08,09) = pd(08,09) - rrt(230) * density(08) 
  pd(09,08) = pd(09,08) - rrt(230) * density(09) 
  pd(09,09) = pd(09,09) - rrt(230) * density(08) 
  pd(12,08) = pd(12,08) + rrt(230) * density(09) 
  pd(12,09) = pd(12,09) + rrt(230) * density(08) 
  pd(07,08) = pd(07,08) + rrt(231) * density(11) 
  pd(07,11) = pd(07,11) + rrt(231) * density(08) 
  pd(08,08) = pd(08,08) - rrt(231) * density(11) 
  pd(08,11) = pd(08,11) - rrt(231) * density(08) 
  pd(11,08) = pd(11,08) - rrt(231) * density(11) 
  pd(11,11) = pd(11,11) - rrt(231) * density(08) 
  pd(12,08) = pd(12,08) + rrt(231) * density(11) 
  pd(12,11) = pd(12,11) + rrt(231) * density(08) 
  pd(07,08) = pd(07,08) + rrt(232) * density(10) 
  pd(07,10) = pd(07,10) + rrt(232) * density(08) 
  pd(08,08) = pd(08,08) - rrt(232) * density(10) 
  pd(08,10) = pd(08,10) - rrt(232) * density(08) 
  pd(10,08) = pd(10,08) - rrt(232) * density(10) 
  pd(10,10) = pd(10,10) - rrt(232) * density(08) 
  pd(12,08) = pd(12,08) + rrt(232) * density(10) 
  pd(12,10) = pd(12,10) + rrt(232) * density(08) 
  pd(08,08) = pd(08,08) - rrt(233) * density(09) 
  pd(08,09) = pd(08,09) - rrt(233) * density(08) 
  pd(09,08) = pd(09,08) - rrt(233) * density(09) 
  pd(09,09) = pd(09,09) - rrt(233) * density(08) 
  pd(25,08) = pd(25,08) + rrt(233) * density(09) 
  pd(25,09) = pd(25,09) + rrt(233) * density(08) 
  pd(44,08) = pd(44,08) + rrt(233) * density(09) 
  pd(44,09) = pd(44,09) + rrt(233) * density(08) 
  pd(08,08) = pd(08,08) - rrt(234) * density(11) 
  pd(08,11) = pd(08,11) - rrt(234) * density(08) 
  pd(11,08) = pd(11,08) - rrt(234) * density(11) 
  pd(11,11) = pd(11,11) - rrt(234) * density(08) 
  pd(25,08) = pd(25,08) + rrt(234) * density(11) 
  pd(25,11) = pd(25,11) + rrt(234) * density(08) 
  pd(44,08) = pd(44,08) + rrt(234) * density(11) 
  pd(44,11) = pd(44,11) + rrt(234) * density(08) 
  pd(08,08) = pd(08,08) - rrt(235) * density(10) 
  pd(08,10) = pd(08,10) - rrt(235) * density(08) 
  pd(10,08) = pd(10,08) - rrt(235) * density(10) 
  pd(10,10) = pd(10,10) - rrt(235) * density(08) 
  pd(25,08) = pd(25,08) + rrt(235) * density(10) 
  pd(25,10) = pd(25,10) + rrt(235) * density(08) 
  pd(44,08) = pd(44,08) + rrt(235) * density(10) 
  pd(44,10) = pd(44,10) + rrt(235) * density(08) 
  pd(05,05) = pd(05,05) - rrt(236) * density(08) 
  pd(05,08) = pd(05,08) - rrt(236) * density(05) 
  pd(08,05) = pd(08,05) - rrt(236) * density(08) 
  pd(08,08) = pd(08,08) - rrt(236) * density(05) 
  pd(19,05) = pd(19,05) + rrt(236) * density(08) 
  pd(19,08) = pd(19,08) + rrt(236) * density(05) 
  pd(44,05) = pd(44,05) + rrt(236) * density(08) 
  pd(44,08) = pd(44,08) + rrt(236) * density(05) 
  pd(03,03) = pd(03,03) - rrt(237) * density(08) 
  pd(03,08) = pd(03,08) - rrt(237) * density(03) 
  pd(08,03) = pd(08,03) - rrt(237) * density(08) 
  pd(08,08) = pd(08,08) - rrt(237) * density(03) 
  pd(17,03) = pd(17,03) + rrt(237) * density(08) 
  pd(17,08) = pd(17,08) + rrt(237) * density(03) 
  pd(44,03) = pd(44,03) + rrt(237) * density(08) 
  pd(44,08) = pd(44,08) + rrt(237) * density(03) 
  pd(08,08) = pd(08,08) - rrt(238) * density(26) 
  pd(08,26) = pd(08,26) - rrt(238) * density(08) 
  pd(09,08) = pd(09,08) + rrt(238) * density(26) 
  pd(09,26) = pd(09,26) + rrt(238) * density(08) 
  pd(25,08) = pd(25,08) + rrt(238) * density(26) 
  pd(25,26) = pd(25,26) + rrt(238) * density(08) 
  pd(26,08) = pd(26,08) - rrt(238) * density(26) 
  pd(26,26) = pd(26,26) - rrt(238) * density(08) 
  pd(08,08) = pd(08,08) - rrt(239) * density(28) 
  pd(08,28) = pd(08,28) - rrt(239) * density(08) 
  pd(09,08) = pd(09,08) + rrt(239) * density(28) 
  pd(09,28) = pd(09,28) + rrt(239) * density(08) 
  pd(25,08) = pd(25,08) + rrt(239) * density(28) 
  pd(25,28) = pd(25,28) + rrt(239) * density(08) 
  pd(28,08) = pd(28,08) - rrt(239) * density(28) 
  pd(28,28) = pd(28,28) - rrt(239) * density(08) 
  pd(08,08) = pd(08,08) - rrt(240) * density(27) 
  pd(08,27) = pd(08,27) - rrt(240) * density(08) 
  pd(09,08) = pd(09,08) + rrt(240) * density(27) 
  pd(09,27) = pd(09,27) + rrt(240) * density(08) 
  pd(25,08) = pd(25,08) + rrt(240) * density(27) 
  pd(25,27) = pd(25,27) + rrt(240) * density(08) 
  pd(27,08) = pd(27,08) - rrt(240) * density(27) 
  pd(27,27) = pd(27,27) - rrt(240) * density(08) 
  pd(08,08) = pd(08,08) - rrt(241) * density(20) 
  pd(08,20) = pd(08,20) - rrt(241) * density(08) 
  pd(09,08) = pd(09,08) + rrt(241) * density(20) 
  pd(09,20) = pd(09,20) + rrt(241) * density(08) 
  pd(19,08) = pd(19,08) + rrt(241) * density(20) 
  pd(19,20) = pd(19,20) + rrt(241) * density(08) 
  pd(20,08) = pd(20,08) - rrt(241) * density(20) 
  pd(20,20) = pd(20,20) - rrt(241) * density(08) 
  pd(08,08) = pd(08,08) - rrt(242) * density(21) 
  pd(08,21) = pd(08,21) - rrt(242) * density(08) 
  pd(09,08) = pd(09,08) + rrt(242) * density(21) 
  pd(09,21) = pd(09,21) + rrt(242) * density(08) 
  pd(19,08) = pd(19,08) + rrt(242) * density(21) 
  pd(19,21) = pd(19,21) + rrt(242) * density(08) 
  pd(21,08) = pd(21,08) - rrt(242) * density(21) 
  pd(21,21) = pd(21,21) - rrt(242) * density(08) 
  pd(08,08) = pd(08,08) - rrt(243) * density(22) 
  pd(08,22) = pd(08,22) - rrt(243) * density(08) 
  pd(09,08) = pd(09,08) + rrt(243) * density(22) 
  pd(09,22) = pd(09,22) + rrt(243) * density(08) 
  pd(19,08) = pd(19,08) + rrt(243) * density(22) 
  pd(19,22) = pd(19,22) + rrt(243) * density(08) 
  pd(22,08) = pd(22,08) - rrt(243) * density(22) 
  pd(22,22) = pd(22,22) - rrt(243) * density(08) 
  pd(07,08) = pd(07,08) + rrt(244) * density(18) 
  pd(07,18) = pd(07,18) + rrt(244) * density(08) 
  pd(08,08) = pd(08,08) - rrt(244) * density(18) 
  pd(08,18) = pd(08,18) - rrt(244) * density(08) 
  pd(18,08) = pd(18,08) - rrt(244) * density(18) 
  pd(18,18) = pd(18,18) - rrt(244) * density(08) 
  pd(19,08) = pd(19,08) + rrt(244) * density(18) 
  pd(19,18) = pd(19,18) + rrt(244) * density(08) 
  pd(06,06) = pd(06,06) - rrt(245) * density(09) 
  pd(06,09) = pd(06,09) - rrt(245) * density(06) 
  pd(07,06) = pd(07,06) + rrt(245) * density(09) 
  pd(07,09) = pd(07,09) + rrt(245) * density(06) 
  pd(08,06) = pd(08,06) + rrt(245) * density(09) 
  pd(08,09) = pd(08,09) + rrt(245) * density(06) 
  pd(09,06) = pd(09,06) - rrt(245) * density(09) 
  pd(09,09) = pd(09,09) - rrt(245) * density(06) 
  pd(06,06) = pd(06,06) - rrt(246) * density(11) 
  pd(06,11) = pd(06,11) - rrt(246) * density(06) 
  pd(07,06) = pd(07,06) + rrt(246) * density(11) 
  pd(07,11) = pd(07,11) + rrt(246) * density(06) 
  pd(08,06) = pd(08,06) + rrt(246) * density(11) 
  pd(08,11) = pd(08,11) + rrt(246) * density(06) 
  pd(11,06) = pd(11,06) - rrt(246) * density(11) 
  pd(11,11) = pd(11,11) - rrt(246) * density(06) 
  pd(06,06) = pd(06,06) - rrt(247) * density(10) 
  pd(06,10) = pd(06,10) - rrt(247) * density(06) 
  pd(07,06) = pd(07,06) + rrt(247) * density(10) 
  pd(07,10) = pd(07,10) + rrt(247) * density(06) 
  pd(08,06) = pd(08,06) + rrt(247) * density(10) 
  pd(08,10) = pd(08,10) + rrt(247) * density(06) 
  pd(10,06) = pd(10,06) - rrt(247) * density(10) 
  pd(10,10) = pd(10,10) - rrt(247) * density(06) 
  pd(06,06) = pd(06,06) - rrt(248) * density(09) 
  pd(06,09) = pd(06,09) - rrt(248) * density(06) 
  pd(09,06) = pd(09,06) - rrt(248) * density(09) 
  pd(09,09) = pd(09,09) - rrt(248) * density(06) 
  pd(25,06) = pd(25,06) + rrt(248) * density(09) 
  pd(25,09) = pd(25,09) + rrt(248) * density(06) 
  pd(43,06) = pd(43,06) + rrt(248) * density(09) 
  pd(43,09) = pd(43,09) + rrt(248) * density(06) 
  pd(06,06) = pd(06,06) - rrt(249) * density(11) 
  pd(06,11) = pd(06,11) - rrt(249) * density(06) 
  pd(11,06) = pd(11,06) - rrt(249) * density(11) 
  pd(11,11) = pd(11,11) - rrt(249) * density(06) 
  pd(25,06) = pd(25,06) + rrt(249) * density(11) 
  pd(25,11) = pd(25,11) + rrt(249) * density(06) 
  pd(43,06) = pd(43,06) + rrt(249) * density(11) 
  pd(43,11) = pd(43,11) + rrt(249) * density(06) 
  pd(06,06) = pd(06,06) - rrt(250) * density(10) 
  pd(06,10) = pd(06,10) - rrt(250) * density(06) 
  pd(10,06) = pd(10,06) - rrt(250) * density(10) 
  pd(10,10) = pd(10,10) - rrt(250) * density(06) 
  pd(25,06) = pd(25,06) + rrt(250) * density(10) 
  pd(25,10) = pd(25,10) + rrt(250) * density(06) 
  pd(43,06) = pd(43,06) + rrt(250) * density(10) 
  pd(43,10) = pd(43,10) + rrt(250) * density(06) 
  pd(06,06) = pd(06,06) - rrt(251) * density(09) 
  pd(06,09) = pd(06,09) - rrt(251) * density(06) 
  pd(09,06) = pd(09,06) - rrt(251) * density(09) 
  pd(09,09) = pd(09,09) - rrt(251) * density(06) 
  pd(23,06) = pd(23,06) + rrt(251) * density(09) 
  pd(23,09) = pd(23,09) + rrt(251) * density(06) 
  pd(44,06) = pd(44,06) + rrt(251) * density(09) 
  pd(44,09) = pd(44,09) + rrt(251) * density(06) 
  pd(06,06) = pd(06,06) - rrt(252) * density(11) 
  pd(06,11) = pd(06,11) - rrt(252) * density(06) 
  pd(11,06) = pd(11,06) - rrt(252) * density(11) 
  pd(11,11) = pd(11,11) - rrt(252) * density(06) 
  pd(23,06) = pd(23,06) + rrt(252) * density(11) 
  pd(23,11) = pd(23,11) + rrt(252) * density(06) 
  pd(44,06) = pd(44,06) + rrt(252) * density(11) 
  pd(44,11) = pd(44,11) + rrt(252) * density(06) 
  pd(06,06) = pd(06,06) - rrt(253) * density(10) 
  pd(06,10) = pd(06,10) - rrt(253) * density(06) 
  pd(10,06) = pd(10,06) - rrt(253) * density(10) 
  pd(10,10) = pd(10,10) - rrt(253) * density(06) 
  pd(23,06) = pd(23,06) + rrt(253) * density(10) 
  pd(23,10) = pd(23,10) + rrt(253) * density(06) 
  pd(44,06) = pd(44,06) + rrt(253) * density(10) 
  pd(44,10) = pd(44,10) + rrt(253) * density(06) 
  pd(06,06) = pd(06,06) - rrt(254) * density(09) 
  pd(06,09) = pd(06,09) - rrt(254) * density(06) 
  pd(09,06) = pd(09,06) - rrt(254) * density(09) 
  pd(09,09) = pd(09,09) - rrt(254) * density(06) 
  pd(19,06) = pd(19,06) + rrt(254) * density(09) 
  pd(19,09) = pd(19,09) + rrt(254) * density(06) 
  pd(43,06) = pd(43,06) + rrt(254) * density(09) 
  pd(43,09) = pd(43,09) + rrt(254) * density(06) 
  pd(44,06) = pd(44,06) + rrt(254) * density(09) 
  pd(44,09) = pd(44,09) + rrt(254) * density(06) 
  pd(06,06) = pd(06,06) - rrt(255) * density(11) 
  pd(06,11) = pd(06,11) - rrt(255) * density(06) 
  pd(11,06) = pd(11,06) - rrt(255) * density(11) 
  pd(11,11) = pd(11,11) - rrt(255) * density(06) 
  pd(19,06) = pd(19,06) + rrt(255) * density(11) 
  pd(19,11) = pd(19,11) + rrt(255) * density(06) 
  pd(43,06) = pd(43,06) + rrt(255) * density(11) 
  pd(43,11) = pd(43,11) + rrt(255) * density(06) 
  pd(44,06) = pd(44,06) + rrt(255) * density(11) 
  pd(44,11) = pd(44,11) + rrt(255) * density(06) 
  pd(06,06) = pd(06,06) - rrt(256) * density(10) 
  pd(06,10) = pd(06,10) - rrt(256) * density(06) 
  pd(10,06) = pd(10,06) - rrt(256) * density(10) 
  pd(10,10) = pd(10,10) - rrt(256) * density(06) 
  pd(19,06) = pd(19,06) + rrt(256) * density(10) 
  pd(19,10) = pd(19,10) + rrt(256) * density(06) 
  pd(43,06) = pd(43,06) + rrt(256) * density(10) 
  pd(43,10) = pd(43,10) + rrt(256) * density(06) 
  pd(44,06) = pd(44,06) + rrt(256) * density(10) 
  pd(44,10) = pd(44,10) + rrt(256) * density(06) 
  pd(06,06) = pd(06,06) - rrt(257) * density(09) 
  pd(06,09) = pd(06,09) - rrt(257) * density(06) 
  pd(09,06) = pd(09,06) - rrt(257) * density(09) 
  pd(09,09) = pd(09,09) - rrt(257) * density(06) 
  pd(17,06) = pd(17,06) + rrt(257) * density(09) 
  pd(17,09) = pd(17,09) + rrt(257) * density(06) 
  pd(44,06) = pd(44,06) + rrt(257) * density(09) * 2.0d0
  pd(44,09) = pd(44,09) + rrt(257) * density(06) * 2.0d0
  pd(06,06) = pd(06,06) - rrt(258) * density(11) 
  pd(06,11) = pd(06,11) - rrt(258) * density(06) 
  pd(11,06) = pd(11,06) - rrt(258) * density(11) 
  pd(11,11) = pd(11,11) - rrt(258) * density(06) 
  pd(17,06) = pd(17,06) + rrt(258) * density(11) 
  pd(17,11) = pd(17,11) + rrt(258) * density(06) 
  pd(44,06) = pd(44,06) + rrt(258) * density(11) * 2.0d0
  pd(44,11) = pd(44,11) + rrt(258) * density(06) * 2.0d0
  pd(06,06) = pd(06,06) - rrt(259) * density(10) 
  pd(06,10) = pd(06,10) - rrt(259) * density(06) 
  pd(10,06) = pd(10,06) - rrt(259) * density(10) 
  pd(10,10) = pd(10,10) - rrt(259) * density(06) 
  pd(17,06) = pd(17,06) + rrt(259) * density(10) 
  pd(17,10) = pd(17,10) + rrt(259) * density(06) 
  pd(44,06) = pd(44,06) + rrt(259) * density(10) * 2.0d0
  pd(44,10) = pd(44,10) + rrt(259) * density(06) * 2.0d0
  pd(06,06) = pd(06,06) - rrt(260) * density(44) 
  pd(06,44) = pd(06,44) - rrt(260) * density(06) 
  pd(08,06) = pd(08,06) + rrt(260) * density(44) 
  pd(08,44) = pd(08,44) + rrt(260) * density(06) 
  pd(43,06) = pd(43,06) + rrt(260) * density(44) 
  pd(43,44) = pd(43,44) + rrt(260) * density(06) 
  pd(44,06) = pd(44,06) - rrt(260) * density(44) 
  pd(44,44) = pd(44,44) - rrt(260) * density(06) 
  pd(04,04) = pd(04,04) - rrt(261) * density(09) 
  pd(04,09) = pd(04,09) - rrt(261) * density(04) 
  pd(09,04) = pd(09,04) - rrt(261) * density(09) 
  pd(09,09) = pd(09,09) - rrt(261) * density(04) 
  pd(23,04) = pd(23,04) + rrt(261) * density(09) 
  pd(23,09) = pd(23,09) + rrt(261) * density(04) 
  pd(43,04) = pd(43,04) + rrt(261) * density(09) 
  pd(43,09) = pd(43,09) + rrt(261) * density(04) 
  pd(04,04) = pd(04,04) - rrt(262) * density(11) 
  pd(04,11) = pd(04,11) - rrt(262) * density(04) 
  pd(11,04) = pd(11,04) - rrt(262) * density(11) 
  pd(11,11) = pd(11,11) - rrt(262) * density(04) 
  pd(23,04) = pd(23,04) + rrt(262) * density(11) 
  pd(23,11) = pd(23,11) + rrt(262) * density(04) 
  pd(43,04) = pd(43,04) + rrt(262) * density(11) 
  pd(43,11) = pd(43,11) + rrt(262) * density(04) 
  pd(04,04) = pd(04,04) - rrt(263) * density(10) 
  pd(04,10) = pd(04,10) - rrt(263) * density(04) 
  pd(10,04) = pd(10,04) - rrt(263) * density(10) 
  pd(10,10) = pd(10,10) - rrt(263) * density(04) 
  pd(23,04) = pd(23,04) + rrt(263) * density(10) 
  pd(23,10) = pd(23,10) + rrt(263) * density(04) 
  pd(43,04) = pd(43,04) + rrt(263) * density(10) 
  pd(43,10) = pd(43,10) + rrt(263) * density(04) 
  pd(04,04) = pd(04,04) - rrt(264) * density(09) 
  pd(04,09) = pd(04,09) - rrt(264) * density(04) 
  pd(09,04) = pd(09,04) - rrt(264) * density(09) 
  pd(09,09) = pd(09,09) - rrt(264) * density(04) 
  pd(19,04) = pd(19,04) + rrt(264) * density(09) 
  pd(19,09) = pd(19,09) + rrt(264) * density(04) 
  pd(44,04) = pd(44,04) + rrt(264) * density(09) 
  pd(44,09) = pd(44,09) + rrt(264) * density(04) 
  pd(04,04) = pd(04,04) - rrt(265) * density(11) 
  pd(04,11) = pd(04,11) - rrt(265) * density(04) 
  pd(11,04) = pd(11,04) - rrt(265) * density(11) 
  pd(11,11) = pd(11,11) - rrt(265) * density(04) 
  pd(19,04) = pd(19,04) + rrt(265) * density(11) 
  pd(19,11) = pd(19,11) + rrt(265) * density(04) 
  pd(44,04) = pd(44,04) + rrt(265) * density(11) 
  pd(44,11) = pd(44,11) + rrt(265) * density(04) 
  pd(04,04) = pd(04,04) - rrt(266) * density(10) 
  pd(04,10) = pd(04,10) - rrt(266) * density(04) 
  pd(10,04) = pd(10,04) - rrt(266) * density(10) 
  pd(10,10) = pd(10,10) - rrt(266) * density(04) 
  pd(19,04) = pd(19,04) + rrt(266) * density(10) 
  pd(19,10) = pd(19,10) + rrt(266) * density(04) 
  pd(44,04) = pd(44,04) + rrt(266) * density(10) 
  pd(44,10) = pd(44,10) + rrt(266) * density(04) 
  pd(04,04) = pd(04,04) - rrt(267) * density(09) 
  pd(04,09) = pd(04,09) - rrt(267) * density(04) 
  pd(09,04) = pd(09,04) - rrt(267) * density(09) 
  pd(09,09) = pd(09,09) - rrt(267) * density(04) 
  pd(17,04) = pd(17,04) + rrt(267) * density(09) 
  pd(17,09) = pd(17,09) + rrt(267) * density(04) 
  pd(43,04) = pd(43,04) + rrt(267) * density(09) 
  pd(43,09) = pd(43,09) + rrt(267) * density(04) 
  pd(44,04) = pd(44,04) + rrt(267) * density(09) 
  pd(44,09) = pd(44,09) + rrt(267) * density(04) 
  pd(04,04) = pd(04,04) - rrt(268) * density(11) 
  pd(04,11) = pd(04,11) - rrt(268) * density(04) 
  pd(11,04) = pd(11,04) - rrt(268) * density(11) 
  pd(11,11) = pd(11,11) - rrt(268) * density(04) 
  pd(17,04) = pd(17,04) + rrt(268) * density(11) 
  pd(17,11) = pd(17,11) + rrt(268) * density(04) 
  pd(43,04) = pd(43,04) + rrt(268) * density(11) 
  pd(43,11) = pd(43,11) + rrt(268) * density(04) 
  pd(44,04) = pd(44,04) + rrt(268) * density(11) 
  pd(44,11) = pd(44,11) + rrt(268) * density(04) 
  pd(04,04) = pd(04,04) - rrt(269) * density(10) 
  pd(04,10) = pd(04,10) - rrt(269) * density(04) 
  pd(10,04) = pd(10,04) - rrt(269) * density(10) 
  pd(10,10) = pd(10,10) - rrt(269) * density(04) 
  pd(17,04) = pd(17,04) + rrt(269) * density(10) 
  pd(17,10) = pd(17,10) + rrt(269) * density(04) 
  pd(43,04) = pd(43,04) + rrt(269) * density(10) 
  pd(43,10) = pd(43,10) + rrt(269) * density(04) 
  pd(44,04) = pd(44,04) + rrt(269) * density(10) 
  pd(44,10) = pd(44,10) + rrt(269) * density(04) 
  pd(04,04) = pd(04,04) - rrt(270) * density(44) 
  pd(04,44) = pd(04,44) - rrt(270) * density(04) 
  pd(06,04) = pd(06,04) + rrt(270) * density(44) 
  pd(06,44) = pd(06,44) + rrt(270) * density(04) 
  pd(43,04) = pd(43,04) + rrt(270) * density(44) 
  pd(43,44) = pd(43,44) + rrt(270) * density(04) 
  pd(44,04) = pd(44,04) - rrt(270) * density(44) 
  pd(44,44) = pd(44,44) - rrt(270) * density(04) 
  pd(20,20) = pd(20,20) - rrt(271) * density(29) 
  pd(20,29) = pd(20,29) - rrt(271) * density(20) 
  pd(23,20) = pd(23,20) + rrt(271) * density(29) 
  pd(23,29) = pd(23,29) + rrt(271) * density(20) 
  pd(26,20) = pd(26,20) + rrt(271) * density(29) 
  pd(26,29) = pd(26,29) + rrt(271) * density(20) 
  pd(29,20) = pd(29,20) - rrt(271) * density(29) 
  pd(29,29) = pd(29,29) - rrt(271) * density(20) 
  pd(21,21) = pd(21,21) - rrt(272) * density(29) 
  pd(21,29) = pd(21,29) - rrt(272) * density(21) 
  pd(23,21) = pd(23,21) + rrt(272) * density(29) 
  pd(23,29) = pd(23,29) + rrt(272) * density(21) 
  pd(26,21) = pd(26,21) + rrt(272) * density(29) 
  pd(26,29) = pd(26,29) + rrt(272) * density(21) 
  pd(29,21) = pd(29,21) - rrt(272) * density(29) 
  pd(29,29) = pd(29,29) - rrt(272) * density(21) 
  pd(22,22) = pd(22,22) - rrt(273) * density(29) 
  pd(22,29) = pd(22,29) - rrt(273) * density(22) 
  pd(23,22) = pd(23,22) + rrt(273) * density(29) 
  pd(23,29) = pd(23,29) + rrt(273) * density(22) 
  pd(26,22) = pd(26,22) + rrt(273) * density(29) 
  pd(26,29) = pd(26,29) + rrt(273) * density(22) 
  pd(29,22) = pd(29,22) - rrt(273) * density(29) 
  pd(29,29) = pd(29,29) - rrt(273) * density(22) 
  pd(13,13) = pd(13,13) - rrt(274) * density(29) 
  pd(13,29) = pd(13,29) - rrt(274) * density(13) 
  pd(18,13) = pd(18,13) + rrt(274) * density(29) 
  pd(18,29) = pd(18,29) + rrt(274) * density(13) 
  pd(25,13) = pd(25,13) + rrt(274) * density(29) 
  pd(25,29) = pd(25,29) + rrt(274) * density(13) 
  pd(29,13) = pd(29,13) - rrt(274) * density(29) 
  pd(29,29) = pd(29,29) - rrt(274) * density(13) 
  pd(16,16) = pd(16,16) - rrt(275) * density(29) 
  pd(16,29) = pd(16,29) - rrt(275) * density(16) 
  pd(18,16) = pd(18,16) + rrt(275) * density(29) 
  pd(18,29) = pd(18,29) + rrt(275) * density(16) 
  pd(25,16) = pd(25,16) + rrt(275) * density(29) 
  pd(25,29) = pd(25,29) + rrt(275) * density(16) 
  pd(29,16) = pd(29,16) - rrt(275) * density(29) 
  pd(29,29) = pd(29,29) - rrt(275) * density(16) 
  pd(15,15) = pd(15,15) - rrt(276) * density(29) 
  pd(15,29) = pd(15,29) - rrt(276) * density(15) 
  pd(18,15) = pd(18,15) + rrt(276) * density(29) 
  pd(18,29) = pd(18,29) + rrt(276) * density(15) 
  pd(25,15) = pd(25,15) + rrt(276) * density(29) 
  pd(25,29) = pd(25,29) + rrt(276) * density(15) 
  pd(29,15) = pd(29,15) - rrt(276) * density(29) 
  pd(29,29) = pd(29,29) - rrt(276) * density(15) 
  pd(14,14) = pd(14,14) - rrt(277) * density(29) 
  pd(14,29) = pd(14,29) - rrt(277) * density(14) 
  pd(18,14) = pd(18,14) + rrt(277) * density(29) 
  pd(18,29) = pd(18,29) + rrt(277) * density(14) 
  pd(25,14) = pd(25,14) + rrt(277) * density(29) 
  pd(25,29) = pd(25,29) + rrt(277) * density(14) 
  pd(29,14) = pd(29,14) - rrt(277) * density(29) 
  pd(29,29) = pd(29,29) - rrt(277) * density(14) 
  pd(25,29) = pd(25,29) + rrt(278) * density(43) 
  pd(25,43) = pd(25,43) + rrt(278) * density(29) 
  pd(29,29) = pd(29,29) - rrt(278) * density(43) 
  pd(29,43) = pd(29,43) - rrt(278) * density(29) 
  pd(43,29) = pd(43,29) - rrt(278) * density(43) 
  pd(43,43) = pd(43,43) - rrt(278) * density(29) 
  pd(44,29) = pd(44,29) + rrt(278) * density(43) 
  pd(44,43) = pd(44,43) + rrt(278) * density(29) 
  pd(23,25) = pd(23,25) + rrt(279) * density(43) 
  pd(23,43) = pd(23,43) + rrt(279) * density(25) 
  pd(25,25) = pd(25,25) - rrt(279) * density(43) 
  pd(25,43) = pd(25,43) - rrt(279) * density(25) 
  pd(43,25) = pd(43,25) - rrt(279) * density(43) 
  pd(43,43) = pd(43,43) - rrt(279) * density(25) 
  pd(44,25) = pd(44,25) + rrt(279) * density(43) 
  pd(44,43) = pd(44,43) + rrt(279) * density(25) 
  pd(13,18) = pd(13,18) + rrt(280) * density(23) 
  pd(13,23) = pd(13,23) + rrt(280) * density(18) 
  pd(18,18) = pd(18,18) - rrt(280) * density(23) 
  pd(18,23) = pd(18,23) - rrt(280) * density(18) 
  pd(23,18) = pd(23,18) - rrt(280) * density(23) 
  pd(23,23) = pd(23,23) - rrt(280) * density(18) 
  pd(25,18) = pd(25,18) + rrt(280) * density(23) 
  pd(25,23) = pd(25,23) + rrt(280) * density(18) 
  pd(18,18) = pd(18,18) - rrt(281) * density(23) 
  pd(18,23) = pd(18,23) - rrt(281) * density(18) 
  pd(19,18) = pd(19,18) + rrt(281) * density(23) 
  pd(19,23) = pd(19,23) + rrt(281) * density(18) 
  pd(20,18) = pd(20,18) + rrt(281) * density(23) 
  pd(20,23) = pd(20,23) + rrt(281) * density(18) 
  pd(23,18) = pd(23,18) - rrt(281) * density(23) 
  pd(23,23) = pd(23,23) - rrt(281) * density(18) 
  pd(19,23) = pd(19,23) + rrt(282) * density(43) 
  pd(19,43) = pd(19,43) + rrt(282) * density(23) 
  pd(23,23) = pd(23,23) - rrt(282) * density(43) 
  pd(23,43) = pd(23,43) - rrt(282) * density(23) 
  pd(43,23) = pd(43,23) - rrt(282) * density(43) 
  pd(43,43) = pd(43,43) - rrt(282) * density(23) 
  pd(44,23) = pd(44,23) + rrt(282) * density(43) 
  pd(44,43) = pd(44,43) + rrt(282) * density(23) 
  pd(19,19) = pd(19,19) - rrt(283) * density(26) 
  pd(19,26) = pd(19,26) - rrt(283) * density(19) 
  pd(20,19) = pd(20,19) + rrt(283) * density(26) 
  pd(20,26) = pd(20,26) + rrt(283) * density(19) 
  pd(25,19) = pd(25,19) + rrt(283) * density(26) 
  pd(25,26) = pd(25,26) + rrt(283) * density(19) 
  pd(26,19) = pd(26,19) - rrt(283) * density(26) 
  pd(26,26) = pd(26,26) - rrt(283) * density(19) 
  pd(19,19) = pd(19,19) - rrt(284) * density(28) 
  pd(19,28) = pd(19,28) - rrt(284) * density(19) 
  pd(20,19) = pd(20,19) + rrt(284) * density(28) 
  pd(20,28) = pd(20,28) + rrt(284) * density(19) 
  pd(25,19) = pd(25,19) + rrt(284) * density(28) 
  pd(25,28) = pd(25,28) + rrt(284) * density(19) 
  pd(28,19) = pd(28,19) - rrt(284) * density(28) 
  pd(28,28) = pd(28,28) - rrt(284) * density(19) 
  pd(19,19) = pd(19,19) - rrt(285) * density(27) 
  pd(19,27) = pd(19,27) - rrt(285) * density(19) 
  pd(20,19) = pd(20,19) + rrt(285) * density(27) 
  pd(20,27) = pd(20,27) + rrt(285) * density(19) 
  pd(25,19) = pd(25,19) + rrt(285) * density(27) 
  pd(25,27) = pd(25,27) + rrt(285) * density(19) 
  pd(27,19) = pd(27,19) - rrt(285) * density(27) 
  pd(27,27) = pd(27,27) - rrt(285) * density(19) 
  pd(13,19) = pd(13,19) + rrt(286) * density(20) 
  pd(13,20) = pd(13,20) + rrt(286) * density(19) 
  pd(19,19) = pd(19,19) - rrt(286) * density(20) 
  pd(19,20) = pd(19,20) - rrt(286) * density(19) 
  pd(20,19) = pd(20,19) - rrt(286) * density(20) 
  pd(20,20) = pd(20,20) - rrt(286) * density(19) 
  pd(25,19) = pd(25,19) + rrt(286) * density(20) 
  pd(25,20) = pd(25,20) + rrt(286) * density(19) 
  pd(13,19) = pd(13,19) + rrt(287) * density(21) 
  pd(13,21) = pd(13,21) + rrt(287) * density(19) 
  pd(19,19) = pd(19,19) - rrt(287) * density(21) 
  pd(19,21) = pd(19,21) - rrt(287) * density(19) 
  pd(21,19) = pd(21,19) - rrt(287) * density(21) 
  pd(21,21) = pd(21,21) - rrt(287) * density(19) 
  pd(25,19) = pd(25,19) + rrt(287) * density(21) 
  pd(25,21) = pd(25,21) + rrt(287) * density(19) 
  pd(13,19) = pd(13,19) + rrt(288) * density(22) 
  pd(13,22) = pd(13,22) + rrt(288) * density(19) 
  pd(19,19) = pd(19,19) - rrt(288) * density(22) 
  pd(19,22) = pd(19,22) - rrt(288) * density(19) 
  pd(22,19) = pd(22,19) - rrt(288) * density(22) 
  pd(22,22) = pd(22,22) - rrt(288) * density(19) 
  pd(25,19) = pd(25,19) + rrt(288) * density(22) 
  pd(25,22) = pd(25,22) + rrt(288) * density(19) 
  pd(17,19) = pd(17,19) + rrt(289) * density(43) 
  pd(17,43) = pd(17,43) + rrt(289) * density(19) 
  pd(19,19) = pd(19,19) - rrt(289) * density(43) 
  pd(19,43) = pd(19,43) - rrt(289) * density(19) 
  pd(43,19) = pd(43,19) - rrt(289) * density(43) 
  pd(43,43) = pd(43,43) - rrt(289) * density(19) 
  pd(44,19) = pd(44,19) + rrt(289) * density(43) 
  pd(44,43) = pd(44,43) + rrt(289) * density(19) 
  pd(07,09) = pd(07,09) + rrt(290) * density(17) 
  pd(07,17) = pd(07,17) + rrt(290) * density(09) 
  pd(09,09) = pd(09,09) - rrt(290) * density(17) 
  pd(09,17) = pd(09,17) - rrt(290) * density(09) 
  pd(17,09) = pd(17,09) - rrt(290) * density(17) 
  pd(17,17) = pd(17,17) - rrt(290) * density(09) 
  pd(19,09) = pd(19,09) + rrt(290) * density(17) 
  pd(19,17) = pd(19,17) + rrt(290) * density(09) 
  pd(07,11) = pd(07,11) + rrt(291) * density(17) 
  pd(07,17) = pd(07,17) + rrt(291) * density(11) 
  pd(11,11) = pd(11,11) - rrt(291) * density(17) 
  pd(11,17) = pd(11,17) - rrt(291) * density(11) 
  pd(17,11) = pd(17,11) - rrt(291) * density(17) 
  pd(17,17) = pd(17,17) - rrt(291) * density(11) 
  pd(19,11) = pd(19,11) + rrt(291) * density(17) 
  pd(19,17) = pd(19,17) + rrt(291) * density(11) 
  pd(07,10) = pd(07,10) + rrt(292) * density(17) 
  pd(07,17) = pd(07,17) + rrt(292) * density(10) 
  pd(10,10) = pd(10,10) - rrt(292) * density(17) 
  pd(10,17) = pd(10,17) - rrt(292) * density(10) 
  pd(17,10) = pd(17,10) - rrt(292) * density(17) 
  pd(17,17) = pd(17,17) - rrt(292) * density(10) 
  pd(19,10) = pd(19,10) + rrt(292) * density(17) 
  pd(19,17) = pd(19,17) + rrt(292) * density(10) 
  pd(17,17) = pd(17,17) - rrt(293) * density(26) 
  pd(17,26) = pd(17,26) - rrt(293) * density(17) 
  pd(18,17) = pd(18,17) + rrt(293) * density(26) 
  pd(18,26) = pd(18,26) + rrt(293) * density(17) 
  pd(25,17) = pd(25,17) + rrt(293) * density(26) 
  pd(25,26) = pd(25,26) + rrt(293) * density(17) 
  pd(26,17) = pd(26,17) - rrt(293) * density(26) 
  pd(26,26) = pd(26,26) - rrt(293) * density(17) 
  pd(17,17) = pd(17,17) - rrt(294) * density(28) 
  pd(17,28) = pd(17,28) - rrt(294) * density(17) 
  pd(18,17) = pd(18,17) + rrt(294) * density(28) 
  pd(18,28) = pd(18,28) + rrt(294) * density(17) 
  pd(25,17) = pd(25,17) + rrt(294) * density(28) 
  pd(25,28) = pd(25,28) + rrt(294) * density(17) 
  pd(28,17) = pd(28,17) - rrt(294) * density(28) 
  pd(28,28) = pd(28,28) - rrt(294) * density(17) 
  pd(17,17) = pd(17,17) - rrt(295) * density(27) 
  pd(17,27) = pd(17,27) - rrt(295) * density(17) 
  pd(18,17) = pd(18,17) + rrt(295) * density(27) 
  pd(18,27) = pd(18,27) + rrt(295) * density(17) 
  pd(25,17) = pd(25,17) + rrt(295) * density(27) 
  pd(25,27) = pd(25,27) + rrt(295) * density(17) 
  pd(27,17) = pd(27,17) - rrt(295) * density(27) 
  pd(27,27) = pd(27,27) - rrt(295) * density(17) 
  pd(17,17) = pd(17,17) - rrt(296) * density(26) 
  pd(17,26) = pd(17,26) - rrt(296) * density(17) 
  pd(20,17) = pd(20,17) + rrt(296) * density(26) 
  pd(20,26) = pd(20,26) + rrt(296) * density(17) 
  pd(23,17) = pd(23,17) + rrt(296) * density(26) 
  pd(23,26) = pd(23,26) + rrt(296) * density(17) 
  pd(26,17) = pd(26,17) - rrt(296) * density(26) 
  pd(26,26) = pd(26,26) - rrt(296) * density(17) 
  pd(17,17) = pd(17,17) - rrt(297) * density(28) 
  pd(17,28) = pd(17,28) - rrt(297) * density(17) 
  pd(20,17) = pd(20,17) + rrt(297) * density(28) 
  pd(20,28) = pd(20,28) + rrt(297) * density(17) 
  pd(23,17) = pd(23,17) + rrt(297) * density(28) 
  pd(23,28) = pd(23,28) + rrt(297) * density(17) 
  pd(28,17) = pd(28,17) - rrt(297) * density(28) 
  pd(28,28) = pd(28,28) - rrt(297) * density(17) 
  pd(17,17) = pd(17,17) - rrt(298) * density(27) 
  pd(17,27) = pd(17,27) - rrt(298) * density(17) 
  pd(20,17) = pd(20,17) + rrt(298) * density(27) 
  pd(20,27) = pd(20,27) + rrt(298) * density(17) 
  pd(23,17) = pd(23,17) + rrt(298) * density(27) 
  pd(23,27) = pd(23,27) + rrt(298) * density(17) 
  pd(27,17) = pd(27,17) - rrt(298) * density(27) 
  pd(27,27) = pd(27,27) - rrt(298) * density(17) 
  pd(13,17) = pd(13,17) + rrt(299) * density(20) 
  pd(13,20) = pd(13,20) + rrt(299) * density(17) 
  pd(17,17) = pd(17,17) - rrt(299) * density(20) 
  pd(17,20) = pd(17,20) - rrt(299) * density(17) 
  pd(20,17) = pd(20,17) - rrt(299) * density(20) 
  pd(20,20) = pd(20,20) - rrt(299) * density(17) 
  pd(23,17) = pd(23,17) + rrt(299) * density(20) 
  pd(23,20) = pd(23,20) + rrt(299) * density(17) 
  pd(13,17) = pd(13,17) + rrt(300) * density(21) 
  pd(13,21) = pd(13,21) + rrt(300) * density(17) 
  pd(17,17) = pd(17,17) - rrt(300) * density(21) 
  pd(17,21) = pd(17,21) - rrt(300) * density(17) 
  pd(21,17) = pd(21,17) - rrt(300) * density(21) 
  pd(21,21) = pd(21,21) - rrt(300) * density(17) 
  pd(23,17) = pd(23,17) + rrt(300) * density(21) 
  pd(23,21) = pd(23,21) + rrt(300) * density(17) 
  pd(13,17) = pd(13,17) + rrt(301) * density(22) 
  pd(13,22) = pd(13,22) + rrt(301) * density(17) 
  pd(17,17) = pd(17,17) - rrt(301) * density(22) 
  pd(17,22) = pd(17,22) - rrt(301) * density(17) 
  pd(22,17) = pd(22,17) - rrt(301) * density(22) 
  pd(22,22) = pd(22,22) - rrt(301) * density(17) 
  pd(23,17) = pd(23,17) + rrt(301) * density(22) 
  pd(23,22) = pd(23,22) + rrt(301) * density(17) 
  pd(13,17) = pd(13,17) + rrt(302) * density(18) 
  pd(13,18) = pd(13,18) + rrt(302) * density(17) 
  pd(17,17) = pd(17,17) - rrt(302) * density(18) 
  pd(17,18) = pd(17,18) - rrt(302) * density(17) 
  pd(18,17) = pd(18,17) - rrt(302) * density(18) 
  pd(18,18) = pd(18,18) - rrt(302) * density(17) 
  pd(19,17) = pd(19,17) + rrt(302) * density(18) 
  pd(19,18) = pd(19,18) + rrt(302) * density(17) 
  pd(17,17) = pd(17,17) - rrt(303) * density(44) 
  pd(17,44) = pd(17,44) - rrt(303) * density(17) 
  pd(19,17) = pd(19,17) + rrt(303) * density(44) 
  pd(19,44) = pd(19,44) + rrt(303) * density(17) 
  pd(43,17) = pd(43,17) + rrt(303) * density(44) 
  pd(43,44) = pd(43,44) + rrt(303) * density(17) 
  pd(44,17) = pd(44,17) - rrt(303) * density(44) 
  pd(44,44) = pd(44,44) - rrt(303) * density(17) 
  pd(07,07) = pd(07,07) - rrt(304) * density(46) 
  pd(07,46) = pd(07,46) - rrt(304) * density(07) 
  pd(12,07) = pd(12,07) + rrt(304) * density(46) 
  pd(12,46) = pd(12,46) + rrt(304) * density(07) 
  pd(44,07) = pd(44,07) + rrt(304) * density(46) 
  pd(44,46) = pd(44,46) + rrt(304) * density(07) 
  pd(46,07) = pd(46,07) - rrt(304) * density(46) 
  pd(46,46) = pd(46,46) - rrt(304) * density(07) 
  pd(05,05) = pd(05,05) - rrt(305) * density(46) 
  pd(05,46) = pd(05,46) - rrt(305) * density(05) 
  pd(08,05) = pd(08,05) + rrt(305) * density(46) 
  pd(08,46) = pd(08,46) + rrt(305) * density(05) 
  pd(44,05) = pd(44,05) + rrt(305) * density(46) 
  pd(44,46) = pd(44,46) + rrt(305) * density(05) 
  pd(46,05) = pd(46,05) - rrt(305) * density(46) 
  pd(46,46) = pd(46,46) - rrt(305) * density(05) 
  pd(03,03) = pd(03,03) - rrt(306) * density(46) 
  pd(03,46) = pd(03,46) - rrt(306) * density(03) 
  pd(06,03) = pd(06,03) + rrt(306) * density(46) 
  pd(06,46) = pd(06,46) + rrt(306) * density(03) 
  pd(44,03) = pd(44,03) + rrt(306) * density(46) 
  pd(44,46) = pd(44,46) + rrt(306) * density(03) 
  pd(46,03) = pd(46,03) - rrt(306) * density(46) 
  pd(46,46) = pd(46,46) - rrt(306) * density(03) 
  pd(25,26) = pd(25,26) + rrt(307) * density(46) 
  pd(25,46) = pd(25,46) + rrt(307) * density(26) 
  pd(26,26) = pd(26,26) - rrt(307) * density(46) 
  pd(26,46) = pd(26,46) - rrt(307) * density(26) 
  pd(44,26) = pd(44,26) + rrt(307) * density(46) * 2.0d0
  pd(44,46) = pd(44,46) + rrt(307) * density(26) * 2.0d0
  pd(46,26) = pd(46,26) - rrt(307) * density(46) 
  pd(46,46) = pd(46,46) - rrt(307) * density(26) 
  pd(25,28) = pd(25,28) + rrt(308) * density(46) 
  pd(25,46) = pd(25,46) + rrt(308) * density(28) 
  pd(28,28) = pd(28,28) - rrt(308) * density(46) 
  pd(28,46) = pd(28,46) - rrt(308) * density(28) 
  pd(44,28) = pd(44,28) + rrt(308) * density(46) * 2.0d0
  pd(44,46) = pd(44,46) + rrt(308) * density(28) * 2.0d0
  pd(46,28) = pd(46,28) - rrt(308) * density(46) 
  pd(46,46) = pd(46,46) - rrt(308) * density(28) 
  pd(25,27) = pd(25,27) + rrt(309) * density(46) 
  pd(25,46) = pd(25,46) + rrt(309) * density(27) 
  pd(27,27) = pd(27,27) - rrt(309) * density(46) 
  pd(27,46) = pd(27,46) - rrt(309) * density(27) 
  pd(44,27) = pd(44,27) + rrt(309) * density(46) * 2.0d0
  pd(44,46) = pd(44,46) + rrt(309) * density(27) * 2.0d0
  pd(46,27) = pd(46,27) - rrt(309) * density(46) 
  pd(46,46) = pd(46,46) - rrt(309) * density(27) 
  pd(24,24) = pd(24,24) - rrt(310) * density(46) 
  pd(24,46) = pd(24,46) - rrt(310) * density(24) 
  pd(29,24) = pd(29,24) + rrt(310) * density(46) 
  pd(29,46) = pd(29,46) + rrt(310) * density(24) 
  pd(44,24) = pd(44,24) + rrt(310) * density(46) 
  pd(44,46) = pd(44,46) + rrt(310) * density(24) 
  pd(46,24) = pd(46,24) - rrt(310) * density(46) 
  pd(46,46) = pd(46,46) - rrt(310) * density(24) 
  pd(20,20) = pd(20,20) - rrt(311) * density(46) 
  pd(20,46) = pd(20,46) - rrt(311) * density(20) 
  pd(25,20) = pd(25,20) + rrt(311) * density(46) 
  pd(25,46) = pd(25,46) + rrt(311) * density(20) 
  pd(44,20) = pd(44,20) + rrt(311) * density(46) 
  pd(44,46) = pd(44,46) + rrt(311) * density(20) 
  pd(46,20) = pd(46,20) - rrt(311) * density(46) 
  pd(46,46) = pd(46,46) - rrt(311) * density(20) 
  pd(21,21) = pd(21,21) - rrt(312) * density(46) 
  pd(21,46) = pd(21,46) - rrt(312) * density(21) 
  pd(25,21) = pd(25,21) + rrt(312) * density(46) 
  pd(25,46) = pd(25,46) + rrt(312) * density(21) 
  pd(44,21) = pd(44,21) + rrt(312) * density(46) 
  pd(44,46) = pd(44,46) + rrt(312) * density(21) 
  pd(46,21) = pd(46,21) - rrt(312) * density(46) 
  pd(46,46) = pd(46,46) - rrt(312) * density(21) 
  pd(22,22) = pd(22,22) - rrt(313) * density(46) 
  pd(22,46) = pd(22,46) - rrt(313) * density(22) 
  pd(25,22) = pd(25,22) + rrt(313) * density(46) 
  pd(25,46) = pd(25,46) + rrt(313) * density(22) 
  pd(44,22) = pd(44,22) + rrt(313) * density(46) 
  pd(44,46) = pd(44,46) + rrt(313) * density(22) 
  pd(46,22) = pd(46,22) - rrt(313) * density(46) 
  pd(46,46) = pd(46,46) - rrt(313) * density(22) 
  pd(19,20) = pd(19,20) + rrt(314) * density(46) 
  pd(19,46) = pd(19,46) + rrt(314) * density(20) 
  pd(20,20) = pd(20,20) - rrt(314) * density(46) 
  pd(20,46) = pd(20,46) - rrt(314) * density(20) 
  pd(44,20) = pd(44,20) + rrt(314) * density(46) * 2.0d0
  pd(44,46) = pd(44,46) + rrt(314) * density(20) * 2.0d0
  pd(46,20) = pd(46,20) - rrt(314) * density(46) 
  pd(46,46) = pd(46,46) - rrt(314) * density(20) 
  pd(19,21) = pd(19,21) + rrt(315) * density(46) 
  pd(19,46) = pd(19,46) + rrt(315) * density(21) 
  pd(21,21) = pd(21,21) - rrt(315) * density(46) 
  pd(21,46) = pd(21,46) - rrt(315) * density(21) 
  pd(44,21) = pd(44,21) + rrt(315) * density(46) * 2.0d0
  pd(44,46) = pd(44,46) + rrt(315) * density(21) * 2.0d0
  pd(46,21) = pd(46,21) - rrt(315) * density(46) 
  pd(46,46) = pd(46,46) - rrt(315) * density(21) 
  pd(19,22) = pd(19,22) + rrt(316) * density(46) 
  pd(19,46) = pd(19,46) + rrt(316) * density(22) 
  pd(22,22) = pd(22,22) - rrt(316) * density(46) 
  pd(22,46) = pd(22,46) - rrt(316) * density(22) 
  pd(44,22) = pd(44,22) + rrt(316) * density(46) * 2.0d0
  pd(44,46) = pd(44,46) + rrt(316) * density(22) * 2.0d0
  pd(46,22) = pd(46,22) - rrt(316) * density(46) 
  pd(46,46) = pd(46,46) - rrt(316) * density(22) 
  pd(18,18) = pd(18,18) - rrt(317) * density(46) 
  pd(18,46) = pd(18,46) - rrt(317) * density(18) 
  pd(23,18) = pd(23,18) + rrt(317) * density(46) 
  pd(23,46) = pd(23,46) + rrt(317) * density(18) 
  pd(44,18) = pd(44,18) + rrt(317) * density(46) 
  pd(44,46) = pd(44,46) + rrt(317) * density(18) 
  pd(46,18) = pd(46,18) - rrt(317) * density(46) 
  pd(46,46) = pd(46,46) - rrt(317) * density(18) 
  pd(13,13) = pd(13,13) - rrt(318) * density(46) 
  pd(13,46) = pd(13,46) - rrt(318) * density(13) 
  pd(19,13) = pd(19,13) + rrt(318) * density(46) 
  pd(19,46) = pd(19,46) + rrt(318) * density(13) 
  pd(44,13) = pd(44,13) + rrt(318) * density(46) 
  pd(44,46) = pd(44,46) + rrt(318) * density(13) 
  pd(46,13) = pd(46,13) - rrt(318) * density(46) 
  pd(46,46) = pd(46,46) - rrt(318) * density(13) 
  pd(16,16) = pd(16,16) - rrt(319) * density(46) 
  pd(16,46) = pd(16,46) - rrt(319) * density(16) 
  pd(19,16) = pd(19,16) + rrt(319) * density(46) 
  pd(19,46) = pd(19,46) + rrt(319) * density(16) 
  pd(44,16) = pd(44,16) + rrt(319) * density(46) 
  pd(44,46) = pd(44,46) + rrt(319) * density(16) 
  pd(46,16) = pd(46,16) - rrt(319) * density(46) 
  pd(46,46) = pd(46,46) - rrt(319) * density(16) 
  pd(15,15) = pd(15,15) - rrt(320) * density(46) 
  pd(15,46) = pd(15,46) - rrt(320) * density(15) 
  pd(19,15) = pd(19,15) + rrt(320) * density(46) 
  pd(19,46) = pd(19,46) + rrt(320) * density(15) 
  pd(44,15) = pd(44,15) + rrt(320) * density(46) 
  pd(44,46) = pd(44,46) + rrt(320) * density(15) 
  pd(46,15) = pd(46,15) - rrt(320) * density(46) 
  pd(46,46) = pd(46,46) - rrt(320) * density(15) 
  pd(14,14) = pd(14,14) - rrt(321) * density(46) 
  pd(14,46) = pd(14,46) - rrt(321) * density(14) 
  pd(19,14) = pd(19,14) + rrt(321) * density(46) 
  pd(19,46) = pd(19,46) + rrt(321) * density(14) 
  pd(44,14) = pd(44,14) + rrt(321) * density(46) 
  pd(44,46) = pd(44,46) + rrt(321) * density(14) 
  pd(46,14) = pd(46,14) - rrt(321) * density(46) 
  pd(46,46) = pd(46,46) - rrt(321) * density(14) 
  pd(09,09) = pd(09,09) - rrt(322) * density(45) 
  pd(09,45) = pd(09,45) - rrt(322) * density(09) 
  pd(12,09) = pd(12,09) + rrt(322) * density(45) 
  pd(12,45) = pd(12,45) + rrt(322) * density(09) 
  pd(44,09) = pd(44,09) + rrt(322) * density(45) 
  pd(44,45) = pd(44,45) + rrt(322) * density(09) 
  pd(45,09) = pd(45,09) - rrt(322) * density(45) 
  pd(45,45) = pd(45,45) - rrt(322) * density(09) 
  pd(11,11) = pd(11,11) - rrt(323) * density(45) 
  pd(11,45) = pd(11,45) - rrt(323) * density(11) 
  pd(12,11) = pd(12,11) + rrt(323) * density(45) 
  pd(12,45) = pd(12,45) + rrt(323) * density(11) 
  pd(44,11) = pd(44,11) + rrt(323) * density(45) 
  pd(44,45) = pd(44,45) + rrt(323) * density(11) 
  pd(45,11) = pd(45,11) - rrt(323) * density(45) 
  pd(45,45) = pd(45,45) - rrt(323) * density(11) 
  pd(10,10) = pd(10,10) - rrt(324) * density(45) 
  pd(10,45) = pd(10,45) - rrt(324) * density(10) 
  pd(12,10) = pd(12,10) + rrt(324) * density(45) 
  pd(12,45) = pd(12,45) + rrt(324) * density(10) 
  pd(44,10) = pd(44,10) + rrt(324) * density(45) 
  pd(44,45) = pd(44,45) + rrt(324) * density(10) 
  pd(45,10) = pd(45,10) - rrt(324) * density(45) 
  pd(45,45) = pd(45,45) - rrt(324) * density(10) 
  pd(08,09) = pd(08,09) + rrt(325) * density(45) 
  pd(08,45) = pd(08,45) + rrt(325) * density(09) 
  pd(09,09) = pd(09,09) - rrt(325) * density(45) 
  pd(09,45) = pd(09,45) - rrt(325) * density(09) 
  pd(43,09) = pd(43,09) + rrt(325) * density(45) 
  pd(43,45) = pd(43,45) + rrt(325) * density(09) 
  pd(44,09) = pd(44,09) + rrt(325) * density(45) 
  pd(44,45) = pd(44,45) + rrt(325) * density(09) 
  pd(45,09) = pd(45,09) - rrt(325) * density(45) 
  pd(45,45) = pd(45,45) - rrt(325) * density(09) 
  pd(08,11) = pd(08,11) + rrt(326) * density(45) 
  pd(08,45) = pd(08,45) + rrt(326) * density(11) 
  pd(11,11) = pd(11,11) - rrt(326) * density(45) 
  pd(11,45) = pd(11,45) - rrt(326) * density(11) 
  pd(43,11) = pd(43,11) + rrt(326) * density(45) 
  pd(43,45) = pd(43,45) + rrt(326) * density(11) 
  pd(44,11) = pd(44,11) + rrt(326) * density(45) 
  pd(44,45) = pd(44,45) + rrt(326) * density(11) 
  pd(45,11) = pd(45,11) - rrt(326) * density(45) 
  pd(45,45) = pd(45,45) - rrt(326) * density(11) 
  pd(08,10) = pd(08,10) + rrt(327) * density(45) 
  pd(08,45) = pd(08,45) + rrt(327) * density(10) 
  pd(10,10) = pd(10,10) - rrt(327) * density(45) 
  pd(10,45) = pd(10,45) - rrt(327) * density(10) 
  pd(43,10) = pd(43,10) + rrt(327) * density(45) 
  pd(43,45) = pd(43,45) + rrt(327) * density(10) 
  pd(44,10) = pd(44,10) + rrt(327) * density(45) 
  pd(44,45) = pd(44,45) + rrt(327) * density(10) 
  pd(45,10) = pd(45,10) - rrt(327) * density(45) 
  pd(45,45) = pd(45,45) - rrt(327) * density(10) 
  pd(05,05) = pd(05,05) - rrt(328) * density(45) 
  pd(05,45) = pd(05,45) - rrt(328) * density(05) 
  pd(08,05) = pd(08,05) + rrt(328) * density(45) 
  pd(08,45) = pd(08,45) + rrt(328) * density(05) 
  pd(43,05) = pd(43,05) + rrt(328) * density(45) 
  pd(43,45) = pd(43,45) + rrt(328) * density(05) 
  pd(45,05) = pd(45,05) - rrt(328) * density(45) 
  pd(45,45) = pd(45,45) - rrt(328) * density(05) 
  pd(05,05) = pd(05,05) - rrt(329) * density(45) 
  pd(05,45) = pd(05,45) - rrt(329) * density(05) 
  pd(06,05) = pd(06,05) + rrt(329) * density(45) 
  pd(06,45) = pd(06,45) + rrt(329) * density(05) 
  pd(44,05) = pd(44,05) + rrt(329) * density(45) 
  pd(44,45) = pd(44,45) + rrt(329) * density(05) 
  pd(45,05) = pd(45,05) - rrt(329) * density(45) 
  pd(45,45) = pd(45,45) - rrt(329) * density(05) 
  pd(03,03) = pd(03,03) - rrt(330) * density(45) 
  pd(03,45) = pd(03,45) - rrt(330) * density(03) 
  pd(06,03) = pd(06,03) + rrt(330) * density(45) 
  pd(06,45) = pd(06,45) + rrt(330) * density(03) 
  pd(43,03) = pd(43,03) + rrt(330) * density(45) 
  pd(43,45) = pd(43,45) + rrt(330) * density(03) 
  pd(45,03) = pd(45,03) - rrt(330) * density(45) 
  pd(45,45) = pd(45,45) - rrt(330) * density(03) 
  pd(03,03) = pd(03,03) - rrt(331) * density(45) 
  pd(03,45) = pd(03,45) - rrt(331) * density(03) 
  pd(04,03) = pd(04,03) + rrt(331) * density(45) 
  pd(04,45) = pd(04,45) + rrt(331) * density(03) 
  pd(44,03) = pd(44,03) + rrt(331) * density(45) 
  pd(44,45) = pd(44,45) + rrt(331) * density(03) 
  pd(45,03) = pd(45,03) - rrt(331) * density(45) 
  pd(45,45) = pd(45,45) - rrt(331) * density(03) 
  pd(26,26) = pd(26,26) - rrt(332) * density(45) 
  pd(26,45) = pd(26,45) - rrt(332) * density(26) 
  pd(29,26) = pd(29,26) + rrt(332) * density(45) 
  pd(29,45) = pd(29,45) + rrt(332) * density(26) 
  pd(44,26) = pd(44,26) + rrt(332) * density(45) 
  pd(44,45) = pd(44,45) + rrt(332) * density(26) 
  pd(45,26) = pd(45,26) - rrt(332) * density(45) 
  pd(45,45) = pd(45,45) - rrt(332) * density(26) 
  pd(28,28) = pd(28,28) - rrt(333) * density(45) 
  pd(28,45) = pd(28,45) - rrt(333) * density(28) 
  pd(29,28) = pd(29,28) + rrt(333) * density(45) 
  pd(29,45) = pd(29,45) + rrt(333) * density(28) 
  pd(44,28) = pd(44,28) + rrt(333) * density(45) 
  pd(44,45) = pd(44,45) + rrt(333) * density(28) 
  pd(45,28) = pd(45,28) - rrt(333) * density(45) 
  pd(45,45) = pd(45,45) - rrt(333) * density(28) 
  pd(27,27) = pd(27,27) - rrt(334) * density(45) 
  pd(27,45) = pd(27,45) - rrt(334) * density(27) 
  pd(29,27) = pd(29,27) + rrt(334) * density(45) 
  pd(29,45) = pd(29,45) + rrt(334) * density(27) 
  pd(44,27) = pd(44,27) + rrt(334) * density(45) 
  pd(44,45) = pd(44,45) + rrt(334) * density(27) 
  pd(45,27) = pd(45,27) - rrt(334) * density(45) 
  pd(45,45) = pd(45,45) - rrt(334) * density(27) 
  pd(25,26) = pd(25,26) + rrt(335) * density(45) 
  pd(25,45) = pd(25,45) + rrt(335) * density(26) 
  pd(26,26) = pd(26,26) - rrt(335) * density(45) 
  pd(26,45) = pd(26,45) - rrt(335) * density(26) 
  pd(43,26) = pd(43,26) + rrt(335) * density(45) 
  pd(43,45) = pd(43,45) + rrt(335) * density(26) 
  pd(44,26) = pd(44,26) + rrt(335) * density(45) 
  pd(44,45) = pd(44,45) + rrt(335) * density(26) 
  pd(45,26) = pd(45,26) - rrt(335) * density(45) 
  pd(45,45) = pd(45,45) - rrt(335) * density(26) 
  pd(25,28) = pd(25,28) + rrt(336) * density(45) 
  pd(25,45) = pd(25,45) + rrt(336) * density(28) 
  pd(28,28) = pd(28,28) - rrt(336) * density(45) 
  pd(28,45) = pd(28,45) - rrt(336) * density(28) 
  pd(43,28) = pd(43,28) + rrt(336) * density(45) 
  pd(43,45) = pd(43,45) + rrt(336) * density(28) 
  pd(44,28) = pd(44,28) + rrt(336) * density(45) 
  pd(44,45) = pd(44,45) + rrt(336) * density(28) 
  pd(45,28) = pd(45,28) - rrt(336) * density(45) 
  pd(45,45) = pd(45,45) - rrt(336) * density(28) 
  pd(25,27) = pd(25,27) + rrt(337) * density(45) 
  pd(25,45) = pd(25,45) + rrt(337) * density(27) 
  pd(27,27) = pd(27,27) - rrt(337) * density(45) 
  pd(27,45) = pd(27,45) - rrt(337) * density(27) 
  pd(43,27) = pd(43,27) + rrt(337) * density(45) 
  pd(43,45) = pd(43,45) + rrt(337) * density(27) 
  pd(44,27) = pd(44,27) + rrt(337) * density(45) 
  pd(44,45) = pd(44,45) + rrt(337) * density(27) 
  pd(45,27) = pd(45,27) - rrt(337) * density(45) 
  pd(45,45) = pd(45,45) - rrt(337) * density(27) 
  pd(23,26) = pd(23,26) + rrt(338) * density(45) 
  pd(23,45) = pd(23,45) + rrt(338) * density(26) 
  pd(26,26) = pd(26,26) - rrt(338) * density(45) 
  pd(26,45) = pd(26,45) - rrt(338) * density(26) 
  pd(44,26) = pd(44,26) + rrt(338) * density(45) * 2.0d0
  pd(44,45) = pd(44,45) + rrt(338) * density(26) * 2.0d0
  pd(45,26) = pd(45,26) - rrt(338) * density(45) 
  pd(45,45) = pd(45,45) - rrt(338) * density(26) 
  pd(23,28) = pd(23,28) + rrt(339) * density(45) 
  pd(23,45) = pd(23,45) + rrt(339) * density(28) 
  pd(28,28) = pd(28,28) - rrt(339) * density(45) 
  pd(28,45) = pd(28,45) - rrt(339) * density(28) 
  pd(44,28) = pd(44,28) + rrt(339) * density(45) * 2.0d0
  pd(44,45) = pd(44,45) + rrt(339) * density(28) * 2.0d0
  pd(45,28) = pd(45,28) - rrt(339) * density(45) 
  pd(45,45) = pd(45,45) - rrt(339) * density(28) 
  pd(23,27) = pd(23,27) + rrt(340) * density(45) 
  pd(23,45) = pd(23,45) + rrt(340) * density(27) 
  pd(27,27) = pd(27,27) - rrt(340) * density(45) 
  pd(27,45) = pd(27,45) - rrt(340) * density(27) 
  pd(44,27) = pd(44,27) + rrt(340) * density(45) * 2.0d0
  pd(44,45) = pd(44,45) + rrt(340) * density(27) * 2.0d0
  pd(45,27) = pd(45,27) - rrt(340) * density(45) 
  pd(45,45) = pd(45,45) - rrt(340) * density(27) 
  pd(19,26) = pd(19,26) + rrt(341) * density(45) 
  pd(19,45) = pd(19,45) + rrt(341) * density(26) 
  pd(26,26) = pd(26,26) - rrt(341) * density(45) 
  pd(26,45) = pd(26,45) - rrt(341) * density(26) 
  pd(43,26) = pd(43,26) + rrt(341) * density(45) 
  pd(43,45) = pd(43,45) + rrt(341) * density(26) 
  pd(44,26) = pd(44,26) + rrt(341) * density(45) * 2.0d0
  pd(44,45) = pd(44,45) + rrt(341) * density(26) * 2.0d0
  pd(45,26) = pd(45,26) - rrt(341) * density(45) 
  pd(45,45) = pd(45,45) - rrt(341) * density(26) 
  pd(19,28) = pd(19,28) + rrt(342) * density(45) 
  pd(19,45) = pd(19,45) + rrt(342) * density(28) 
  pd(28,28) = pd(28,28) - rrt(342) * density(45) 
  pd(28,45) = pd(28,45) - rrt(342) * density(28) 
  pd(43,28) = pd(43,28) + rrt(342) * density(45) 
  pd(43,45) = pd(43,45) + rrt(342) * density(28) 
  pd(44,28) = pd(44,28) + rrt(342) * density(45) * 2.0d0
  pd(44,45) = pd(44,45) + rrt(342) * density(28) * 2.0d0
  pd(45,28) = pd(45,28) - rrt(342) * density(45) 
  pd(45,45) = pd(45,45) - rrt(342) * density(28) 
  pd(19,27) = pd(19,27) + rrt(343) * density(45) 
  pd(19,45) = pd(19,45) + rrt(343) * density(27) 
  pd(27,27) = pd(27,27) - rrt(343) * density(45) 
  pd(27,45) = pd(27,45) - rrt(343) * density(27) 
  pd(43,27) = pd(43,27) + rrt(343) * density(45) 
  pd(43,45) = pd(43,45) + rrt(343) * density(27) 
  pd(44,27) = pd(44,27) + rrt(343) * density(45) * 2.0d0
  pd(44,45) = pd(44,45) + rrt(343) * density(27) * 2.0d0
  pd(45,27) = pd(45,27) - rrt(343) * density(45) 
  pd(45,45) = pd(45,45) - rrt(343) * density(27) 
  pd(17,26) = pd(17,26) + rrt(344) * density(45) 
  pd(17,45) = pd(17,45) + rrt(344) * density(26) 
  pd(26,26) = pd(26,26) - rrt(344) * density(45) 
  pd(26,45) = pd(26,45) - rrt(344) * density(26) 
  pd(44,26) = pd(44,26) + rrt(344) * density(45) * 3.0d0
  pd(44,45) = pd(44,45) + rrt(344) * density(26) * 3.0d0
  pd(45,26) = pd(45,26) - rrt(344) * density(45) 
  pd(45,45) = pd(45,45) - rrt(344) * density(26) 
  pd(17,28) = pd(17,28) + rrt(345) * density(45) 
  pd(17,45) = pd(17,45) + rrt(345) * density(28) 
  pd(28,28) = pd(28,28) - rrt(345) * density(45) 
  pd(28,45) = pd(28,45) - rrt(345) * density(28) 
  pd(44,28) = pd(44,28) + rrt(345) * density(45) * 3.0d0
  pd(44,45) = pd(44,45) + rrt(345) * density(28) * 3.0d0
  pd(45,28) = pd(45,28) - rrt(345) * density(45) 
  pd(45,45) = pd(45,45) - rrt(345) * density(28) 
  pd(17,27) = pd(17,27) + rrt(346) * density(45) 
  pd(17,45) = pd(17,45) + rrt(346) * density(27) 
  pd(27,27) = pd(27,27) - rrt(346) * density(45) 
  pd(27,45) = pd(27,45) - rrt(346) * density(27) 
  pd(44,27) = pd(44,27) + rrt(346) * density(45) * 3.0d0
  pd(44,45) = pd(44,45) + rrt(346) * density(27) * 3.0d0
  pd(45,27) = pd(45,27) - rrt(346) * density(45) 
  pd(45,45) = pd(45,45) - rrt(346) * density(27) 
  pd(20,20) = pd(20,20) - rrt(347) * density(45) 
  pd(20,45) = pd(20,45) - rrt(347) * density(20) 
  pd(23,20) = pd(23,20) + rrt(347) * density(45) 
  pd(23,45) = pd(23,45) + rrt(347) * density(20) 
  pd(44,20) = pd(44,20) + rrt(347) * density(45) 
  pd(44,45) = pd(44,45) + rrt(347) * density(20) 
  pd(45,20) = pd(45,20) - rrt(347) * density(45) 
  pd(45,45) = pd(45,45) - rrt(347) * density(20) 
  pd(21,21) = pd(21,21) - rrt(348) * density(45) 
  pd(21,45) = pd(21,45) - rrt(348) * density(21) 
  pd(23,21) = pd(23,21) + rrt(348) * density(45) 
  pd(23,45) = pd(23,45) + rrt(348) * density(21) 
  pd(44,21) = pd(44,21) + rrt(348) * density(45) 
  pd(44,45) = pd(44,45) + rrt(348) * density(21) 
  pd(45,21) = pd(45,21) - rrt(348) * density(45) 
  pd(45,45) = pd(45,45) - rrt(348) * density(21) 
  pd(22,22) = pd(22,22) - rrt(349) * density(45) 
  pd(22,45) = pd(22,45) - rrt(349) * density(22) 
  pd(23,22) = pd(23,22) + rrt(349) * density(45) 
  pd(23,45) = pd(23,45) + rrt(349) * density(22) 
  pd(44,22) = pd(44,22) + rrt(349) * density(45) 
  pd(44,45) = pd(44,45) + rrt(349) * density(22) 
  pd(45,22) = pd(45,22) - rrt(349) * density(45) 
  pd(45,45) = pd(45,45) - rrt(349) * density(22) 
  pd(19,20) = pd(19,20) + rrt(350) * density(45) 
  pd(19,45) = pd(19,45) + rrt(350) * density(20) 
  pd(20,20) = pd(20,20) - rrt(350) * density(45) 
  pd(20,45) = pd(20,45) - rrt(350) * density(20) 
  pd(43,20) = pd(43,20) + rrt(350) * density(45) 
  pd(43,45) = pd(43,45) + rrt(350) * density(20) 
  pd(44,20) = pd(44,20) + rrt(350) * density(45) 
  pd(44,45) = pd(44,45) + rrt(350) * density(20) 
  pd(45,20) = pd(45,20) - rrt(350) * density(45) 
  pd(45,45) = pd(45,45) - rrt(350) * density(20) 
  pd(19,21) = pd(19,21) + rrt(351) * density(45) 
  pd(19,45) = pd(19,45) + rrt(351) * density(21) 
  pd(21,21) = pd(21,21) - rrt(351) * density(45) 
  pd(21,45) = pd(21,45) - rrt(351) * density(21) 
  pd(43,21) = pd(43,21) + rrt(351) * density(45) 
  pd(43,45) = pd(43,45) + rrt(351) * density(21) 
  pd(44,21) = pd(44,21) + rrt(351) * density(45) 
  pd(44,45) = pd(44,45) + rrt(351) * density(21) 
  pd(45,21) = pd(45,21) - rrt(351) * density(45) 
  pd(45,45) = pd(45,45) - rrt(351) * density(21) 
  pd(19,22) = pd(19,22) + rrt(352) * density(45) 
  pd(19,45) = pd(19,45) + rrt(352) * density(22) 
  pd(22,22) = pd(22,22) - rrt(352) * density(45) 
  pd(22,45) = pd(22,45) - rrt(352) * density(22) 
  pd(43,22) = pd(43,22) + rrt(352) * density(45) 
  pd(43,45) = pd(43,45) + rrt(352) * density(22) 
  pd(44,22) = pd(44,22) + rrt(352) * density(45) 
  pd(44,45) = pd(44,45) + rrt(352) * density(22) 
  pd(45,22) = pd(45,22) - rrt(352) * density(45) 
  pd(45,45) = pd(45,45) - rrt(352) * density(22) 
  pd(17,20) = pd(17,20) + rrt(353) * density(45) 
  pd(17,45) = pd(17,45) + rrt(353) * density(20) 
  pd(20,20) = pd(20,20) - rrt(353) * density(45) 
  pd(20,45) = pd(20,45) - rrt(353) * density(20) 
  pd(44,20) = pd(44,20) + rrt(353) * density(45) * 2.0d0
  pd(44,45) = pd(44,45) + rrt(353) * density(20) * 2.0d0
  pd(45,20) = pd(45,20) - rrt(353) * density(45) 
  pd(45,45) = pd(45,45) - rrt(353) * density(20) 
  pd(17,21) = pd(17,21) + rrt(354) * density(45) 
  pd(17,45) = pd(17,45) + rrt(354) * density(21) 
  pd(21,21) = pd(21,21) - rrt(354) * density(45) 
  pd(21,45) = pd(21,45) - rrt(354) * density(21) 
  pd(44,21) = pd(44,21) + rrt(354) * density(45) * 2.0d0
  pd(44,45) = pd(44,45) + rrt(354) * density(21) * 2.0d0
  pd(45,21) = pd(45,21) - rrt(354) * density(45) 
  pd(45,45) = pd(45,45) - rrt(354) * density(21) 
  pd(17,22) = pd(17,22) + rrt(355) * density(45) 
  pd(17,45) = pd(17,45) + rrt(355) * density(22) 
  pd(22,22) = pd(22,22) - rrt(355) * density(45) 
  pd(22,45) = pd(22,45) - rrt(355) * density(22) 
  pd(44,22) = pd(44,22) + rrt(355) * density(45) * 2.0d0
  pd(44,45) = pd(44,45) + rrt(355) * density(22) * 2.0d0
  pd(45,22) = pd(45,22) - rrt(355) * density(45) 
  pd(45,45) = pd(45,45) - rrt(355) * density(22) 
  pd(13,13) = pd(13,13) - rrt(356) * density(45) 
  pd(13,45) = pd(13,45) - rrt(356) * density(13) 
  pd(19,13) = pd(19,13) + rrt(356) * density(45) 
  pd(19,45) = pd(19,45) + rrt(356) * density(13) 
  pd(43,13) = pd(43,13) + rrt(356) * density(45) 
  pd(43,45) = pd(43,45) + rrt(356) * density(13) 
  pd(45,13) = pd(45,13) - rrt(356) * density(45) 
  pd(45,45) = pd(45,45) - rrt(356) * density(13) 
  pd(16,16) = pd(16,16) - rrt(357) * density(45) 
  pd(16,45) = pd(16,45) - rrt(357) * density(16) 
  pd(19,16) = pd(19,16) + rrt(357) * density(45) 
  pd(19,45) = pd(19,45) + rrt(357) * density(16) 
  pd(43,16) = pd(43,16) + rrt(357) * density(45) 
  pd(43,45) = pd(43,45) + rrt(357) * density(16) 
  pd(45,16) = pd(45,16) - rrt(357) * density(45) 
  pd(45,45) = pd(45,45) - rrt(357) * density(16) 
  pd(15,15) = pd(15,15) - rrt(358) * density(45) 
  pd(15,45) = pd(15,45) - rrt(358) * density(15) 
  pd(19,15) = pd(19,15) + rrt(358) * density(45) 
  pd(19,45) = pd(19,45) + rrt(358) * density(15) 
  pd(43,15) = pd(43,15) + rrt(358) * density(45) 
  pd(43,45) = pd(43,45) + rrt(358) * density(15) 
  pd(45,15) = pd(45,15) - rrt(358) * density(45) 
  pd(45,45) = pd(45,45) - rrt(358) * density(15) 
  pd(14,14) = pd(14,14) - rrt(359) * density(45) 
  pd(14,45) = pd(14,45) - rrt(359) * density(14) 
  pd(19,14) = pd(19,14) + rrt(359) * density(45) 
  pd(19,45) = pd(19,45) + rrt(359) * density(14) 
  pd(43,14) = pd(43,14) + rrt(359) * density(45) 
  pd(43,45) = pd(43,45) + rrt(359) * density(14) 
  pd(45,14) = pd(45,14) - rrt(359) * density(45) 
  pd(45,45) = pd(45,45) - rrt(359) * density(14) 
  pd(13,13) = pd(13,13) - rrt(360) * density(45) 
  pd(13,45) = pd(13,45) - rrt(360) * density(13) 
  pd(17,13) = pd(17,13) + rrt(360) * density(45) 
  pd(17,45) = pd(17,45) + rrt(360) * density(13) 
  pd(44,13) = pd(44,13) + rrt(360) * density(45) 
  pd(44,45) = pd(44,45) + rrt(360) * density(13) 
  pd(45,13) = pd(45,13) - rrt(360) * density(45) 
  pd(45,45) = pd(45,45) - rrt(360) * density(13) 
  pd(16,16) = pd(16,16) - rrt(361) * density(45) 
  pd(16,45) = pd(16,45) - rrt(361) * density(16) 
  pd(17,16) = pd(17,16) + rrt(361) * density(45) 
  pd(17,45) = pd(17,45) + rrt(361) * density(16) 
  pd(44,16) = pd(44,16) + rrt(361) * density(45) 
  pd(44,45) = pd(44,45) + rrt(361) * density(16) 
  pd(45,16) = pd(45,16) - rrt(361) * density(45) 
  pd(45,45) = pd(45,45) - rrt(361) * density(16) 
  pd(15,15) = pd(15,15) - rrt(362) * density(45) 
  pd(15,45) = pd(15,45) - rrt(362) * density(15) 
  pd(17,15) = pd(17,15) + rrt(362) * density(45) 
  pd(17,45) = pd(17,45) + rrt(362) * density(15) 
  pd(44,15) = pd(44,15) + rrt(362) * density(45) 
  pd(44,45) = pd(44,45) + rrt(362) * density(15) 
  pd(45,15) = pd(45,15) - rrt(362) * density(45) 
  pd(45,45) = pd(45,45) - rrt(362) * density(15) 
  pd(14,14) = pd(14,14) - rrt(363) * density(45) 
  pd(14,45) = pd(14,45) - rrt(363) * density(14) 
  pd(17,14) = pd(17,14) + rrt(363) * density(45) 
  pd(17,45) = pd(17,45) + rrt(363) * density(14) 
  pd(44,14) = pd(44,14) + rrt(363) * density(45) 
  pd(44,45) = pd(44,45) + rrt(363) * density(14) 
  pd(45,14) = pd(45,14) - rrt(363) * density(45) 
  pd(45,45) = pd(45,45) - rrt(363) * density(14) 
  pd(43,43) = pd(43,43) - rrt(364) * density(45) 
  pd(43,45) = pd(43,45) - rrt(364) * density(43) 
  pd(45,43) = pd(45,43) - rrt(364) * density(45) 
  pd(45,45) = pd(45,45) - rrt(364) * density(43) 
  pd(46,43) = pd(46,43) + rrt(364) * density(45) 
  pd(46,45) = pd(46,45) + rrt(364) * density(43) 
  pd(43,43) = pd(43,43) - rrt(365) * density(45) 
  pd(43,45) = pd(43,45) - rrt(365) * density(43) 
  pd(44,43) = pd(44,43) + rrt(365) * density(45) 
  pd(44,45) = pd(44,45) + rrt(365) * density(43) 
  pd(45,43) = pd(45,43) - rrt(365) * density(45) 
  pd(45,45) = pd(45,45) - rrt(365) * density(43) 
  pd(47,43) = pd(47,43) + rrt(365) * density(45) 
  pd(47,45) = pd(47,45) + rrt(365) * density(43) 
  pd(09,09) = pd(09,09) - rrt(366) * density(47) 
  pd(09,47) = pd(09,47) - rrt(366) * density(09) 
  pd(12,09) = pd(12,09) + rrt(366) * density(47) 
  pd(12,47) = pd(12,47) + rrt(366) * density(09) 
  pd(43,09) = pd(43,09) + rrt(366) * density(47) 
  pd(43,47) = pd(43,47) + rrt(366) * density(09) 
  pd(47,09) = pd(47,09) - rrt(366) * density(47) 
  pd(47,47) = pd(47,47) - rrt(366) * density(09) 
  pd(11,11) = pd(11,11) - rrt(367) * density(47) 
  pd(11,47) = pd(11,47) - rrt(367) * density(11) 
  pd(12,11) = pd(12,11) + rrt(367) * density(47) 
  pd(12,47) = pd(12,47) + rrt(367) * density(11) 
  pd(43,11) = pd(43,11) + rrt(367) * density(47) 
  pd(43,47) = pd(43,47) + rrt(367) * density(11) 
  pd(47,11) = pd(47,11) - rrt(367) * density(47) 
  pd(47,47) = pd(47,47) - rrt(367) * density(11) 
  pd(10,10) = pd(10,10) - rrt(368) * density(47) 
  pd(10,47) = pd(10,47) - rrt(368) * density(10) 
  pd(12,10) = pd(12,10) + rrt(368) * density(47) 
  pd(12,47) = pd(12,47) + rrt(368) * density(10) 
  pd(43,10) = pd(43,10) + rrt(368) * density(47) 
  pd(43,47) = pd(43,47) + rrt(368) * density(10) 
  pd(47,10) = pd(47,10) - rrt(368) * density(47) 
  pd(47,47) = pd(47,47) - rrt(368) * density(10) 
  pd(08,09) = pd(08,09) + rrt(369) * density(47) 
  pd(08,47) = pd(08,47) + rrt(369) * density(09) 
  pd(09,09) = pd(09,09) - rrt(369) * density(47) 
  pd(09,47) = pd(09,47) - rrt(369) * density(09) 
  pd(44,09) = pd(44,09) + rrt(369) * density(47) 
  pd(44,47) = pd(44,47) + rrt(369) * density(09) 
  pd(47,09) = pd(47,09) - rrt(369) * density(47) 
  pd(47,47) = pd(47,47) - rrt(369) * density(09) 
  pd(08,11) = pd(08,11) + rrt(370) * density(47) 
  pd(08,47) = pd(08,47) + rrt(370) * density(11) 
  pd(11,11) = pd(11,11) - rrt(370) * density(47) 
  pd(11,47) = pd(11,47) - rrt(370) * density(11) 
  pd(44,11) = pd(44,11) + rrt(370) * density(47) 
  pd(44,47) = pd(44,47) + rrt(370) * density(11) 
  pd(47,11) = pd(47,11) - rrt(370) * density(47) 
  pd(47,47) = pd(47,47) - rrt(370) * density(11) 
  pd(08,10) = pd(08,10) + rrt(371) * density(47) 
  pd(08,47) = pd(08,47) + rrt(371) * density(10) 
  pd(10,10) = pd(10,10) - rrt(371) * density(47) 
  pd(10,47) = pd(10,47) - rrt(371) * density(10) 
  pd(44,10) = pd(44,10) + rrt(371) * density(47) 
  pd(44,47) = pd(44,47) + rrt(371) * density(10) 
  pd(47,10) = pd(47,10) - rrt(371) * density(47) 
  pd(47,47) = pd(47,47) - rrt(371) * density(10) 
  pd(07,07) = pd(07,07) - rrt(372) * density(47) 
  pd(07,47) = pd(07,47) - rrt(372) * density(07) 
  pd(08,07) = pd(08,07) + rrt(372) * density(47) 
  pd(08,47) = pd(08,47) + rrt(372) * density(07) 
  pd(43,07) = pd(43,07) + rrt(372) * density(47) 
  pd(43,47) = pd(43,47) + rrt(372) * density(07) 
  pd(47,07) = pd(47,07) - rrt(372) * density(47) 
  pd(47,47) = pd(47,47) - rrt(372) * density(07) 
  pd(05,05) = pd(05,05) - rrt(373) * density(47) 
  pd(05,47) = pd(05,47) - rrt(373) * density(05) 
  pd(06,05) = pd(06,05) + rrt(373) * density(47) 
  pd(06,47) = pd(06,47) + rrt(373) * density(05) 
  pd(43,05) = pd(43,05) + rrt(373) * density(47) 
  pd(43,47) = pd(43,47) + rrt(373) * density(05) 
  pd(47,05) = pd(47,05) - rrt(373) * density(47) 
  pd(47,47) = pd(47,47) - rrt(373) * density(05) 
  pd(04,05) = pd(04,05) + rrt(374) * density(47) 
  pd(04,47) = pd(04,47) + rrt(374) * density(05) 
  pd(05,05) = pd(05,05) - rrt(374) * density(47) 
  pd(05,47) = pd(05,47) - rrt(374) * density(05) 
  pd(44,05) = pd(44,05) + rrt(374) * density(47) 
  pd(44,47) = pd(44,47) + rrt(374) * density(05) 
  pd(47,05) = pd(47,05) - rrt(374) * density(47) 
  pd(47,47) = pd(47,47) - rrt(374) * density(05) 
  pd(03,03) = pd(03,03) - rrt(375) * density(47) 
  pd(03,47) = pd(03,47) - rrt(375) * density(03) 
  pd(04,03) = pd(04,03) + rrt(375) * density(47) 
  pd(04,47) = pd(04,47) + rrt(375) * density(03) 
  pd(43,03) = pd(43,03) + rrt(375) * density(47) 
  pd(43,47) = pd(43,47) + rrt(375) * density(03) 
  pd(47,03) = pd(47,03) - rrt(375) * density(47) 
  pd(47,47) = pd(47,47) - rrt(375) * density(03) 
  pd(25,26) = pd(25,26) + rrt(376) * density(47) 
  pd(25,47) = pd(25,47) + rrt(376) * density(26) 
  pd(26,26) = pd(26,26) - rrt(376) * density(47) 
  pd(26,47) = pd(26,47) - rrt(376) * density(26) 
  pd(44,26) = pd(44,26) + rrt(376) * density(47) 
  pd(44,47) = pd(44,47) + rrt(376) * density(26) 
  pd(47,26) = pd(47,26) - rrt(376) * density(47) 
  pd(47,47) = pd(47,47) - rrt(376) * density(26) 
  pd(25,28) = pd(25,28) + rrt(377) * density(47) 
  pd(25,47) = pd(25,47) + rrt(377) * density(28) 
  pd(28,28) = pd(28,28) - rrt(377) * density(47) 
  pd(28,47) = pd(28,47) - rrt(377) * density(28) 
  pd(44,28) = pd(44,28) + rrt(377) * density(47) 
  pd(44,47) = pd(44,47) + rrt(377) * density(28) 
  pd(47,28) = pd(47,28) - rrt(377) * density(47) 
  pd(47,47) = pd(47,47) - rrt(377) * density(28) 
  pd(25,27) = pd(25,27) + rrt(378) * density(47) 
  pd(25,47) = pd(25,47) + rrt(378) * density(27) 
  pd(27,27) = pd(27,27) - rrt(378) * density(47) 
  pd(27,47) = pd(27,47) - rrt(378) * density(27) 
  pd(44,27) = pd(44,27) + rrt(378) * density(47) 
  pd(44,47) = pd(44,47) + rrt(378) * density(27) 
  pd(47,27) = pd(47,27) - rrt(378) * density(47) 
  pd(47,47) = pd(47,47) - rrt(378) * density(27) 
  pd(23,26) = pd(23,26) + rrt(379) * density(47) 
  pd(23,47) = pd(23,47) + rrt(379) * density(26) 
  pd(26,26) = pd(26,26) - rrt(379) * density(47) 
  pd(26,47) = pd(26,47) - rrt(379) * density(26) 
  pd(43,26) = pd(43,26) + rrt(379) * density(47) 
  pd(43,47) = pd(43,47) + rrt(379) * density(26) 
  pd(44,26) = pd(44,26) + rrt(379) * density(47) 
  pd(44,47) = pd(44,47) + rrt(379) * density(26) 
  pd(47,26) = pd(47,26) - rrt(379) * density(47) 
  pd(47,47) = pd(47,47) - rrt(379) * density(26) 
  pd(23,28) = pd(23,28) + rrt(380) * density(47) 
  pd(23,47) = pd(23,47) + rrt(380) * density(28) 
  pd(28,28) = pd(28,28) - rrt(380) * density(47) 
  pd(28,47) = pd(28,47) - rrt(380) * density(28) 
  pd(43,28) = pd(43,28) + rrt(380) * density(47) 
  pd(43,47) = pd(43,47) + rrt(380) * density(28) 
  pd(44,28) = pd(44,28) + rrt(380) * density(47) 
  pd(44,47) = pd(44,47) + rrt(380) * density(28) 
  pd(47,28) = pd(47,28) - rrt(380) * density(47) 
  pd(47,47) = pd(47,47) - rrt(380) * density(28) 
  pd(23,27) = pd(23,27) + rrt(381) * density(47) 
  pd(23,47) = pd(23,47) + rrt(381) * density(27) 
  pd(27,27) = pd(27,27) - rrt(381) * density(47) 
  pd(27,47) = pd(27,47) - rrt(381) * density(27) 
  pd(43,27) = pd(43,27) + rrt(381) * density(47) 
  pd(43,47) = pd(43,47) + rrt(381) * density(27) 
  pd(44,27) = pd(44,27) + rrt(381) * density(47) 
  pd(44,47) = pd(44,47) + rrt(381) * density(27) 
  pd(47,27) = pd(47,27) - rrt(381) * density(47) 
  pd(47,47) = pd(47,47) - rrt(381) * density(27) 
  pd(19,26) = pd(19,26) + rrt(382) * density(47) 
  pd(19,47) = pd(19,47) + rrt(382) * density(26) 
  pd(26,26) = pd(26,26) - rrt(382) * density(47) 
  pd(26,47) = pd(26,47) - rrt(382) * density(26) 
  pd(44,26) = pd(44,26) + rrt(382) * density(47) * 2.0d0
  pd(44,47) = pd(44,47) + rrt(382) * density(26) * 2.0d0
  pd(47,26) = pd(47,26) - rrt(382) * density(47) 
  pd(47,47) = pd(47,47) - rrt(382) * density(26) 
  pd(19,28) = pd(19,28) + rrt(383) * density(47) 
  pd(19,47) = pd(19,47) + rrt(383) * density(28) 
  pd(28,28) = pd(28,28) - rrt(383) * density(47) 
  pd(28,47) = pd(28,47) - rrt(383) * density(28) 
  pd(44,28) = pd(44,28) + rrt(383) * density(47) * 2.0d0
  pd(44,47) = pd(44,47) + rrt(383) * density(28) * 2.0d0
  pd(47,28) = pd(47,28) - rrt(383) * density(47) 
  pd(47,47) = pd(47,47) - rrt(383) * density(28) 
  pd(19,27) = pd(19,27) + rrt(384) * density(47) 
  pd(19,47) = pd(19,47) + rrt(384) * density(27) 
  pd(27,27) = pd(27,27) - rrt(384) * density(47) 
  pd(27,47) = pd(27,47) - rrt(384) * density(27) 
  pd(44,27) = pd(44,27) + rrt(384) * density(47) * 2.0d0
  pd(44,47) = pd(44,47) + rrt(384) * density(27) * 2.0d0
  pd(47,27) = pd(47,27) - rrt(384) * density(47) 
  pd(47,47) = pd(47,47) - rrt(384) * density(27) 
  pd(23,24) = pd(23,24) + rrt(385) * density(47) 
  pd(23,47) = pd(23,47) + rrt(385) * density(24) 
  pd(24,24) = pd(24,24) - rrt(385) * density(47) 
  pd(24,47) = pd(24,47) - rrt(385) * density(24) 
  pd(44,24) = pd(44,24) + rrt(385) * density(47) 
  pd(44,47) = pd(44,47) + rrt(385) * density(24) 
  pd(47,24) = pd(47,24) - rrt(385) * density(47) 
  pd(47,47) = pd(47,47) - rrt(385) * density(24) 
  pd(19,24) = pd(19,24) + rrt(386) * density(47) 
  pd(19,47) = pd(19,47) + rrt(386) * density(24) 
  pd(24,24) = pd(24,24) - rrt(386) * density(47) 
  pd(24,47) = pd(24,47) - rrt(386) * density(24) 
  pd(43,24) = pd(43,24) + rrt(386) * density(47) 
  pd(43,47) = pd(43,47) + rrt(386) * density(24) 
  pd(44,24) = pd(44,24) + rrt(386) * density(47) 
  pd(44,47) = pd(44,47) + rrt(386) * density(24) 
  pd(47,24) = pd(47,24) - rrt(386) * density(47) 
  pd(47,47) = pd(47,47) - rrt(386) * density(24) 
  pd(20,20) = pd(20,20) - rrt(387) * density(47) 
  pd(20,47) = pd(20,47) - rrt(387) * density(20) 
  pd(23,20) = pd(23,20) + rrt(387) * density(47) 
  pd(23,47) = pd(23,47) + rrt(387) * density(20) 
  pd(43,20) = pd(43,20) + rrt(387) * density(47) 
  pd(43,47) = pd(43,47) + rrt(387) * density(20) 
  pd(47,20) = pd(47,20) - rrt(387) * density(47) 
  pd(47,47) = pd(47,47) - rrt(387) * density(20) 
  pd(21,21) = pd(21,21) - rrt(388) * density(47) 
  pd(21,47) = pd(21,47) - rrt(388) * density(21) 
  pd(23,21) = pd(23,21) + rrt(388) * density(47) 
  pd(23,47) = pd(23,47) + rrt(388) * density(21) 
  pd(43,21) = pd(43,21) + rrt(388) * density(47) 
  pd(43,47) = pd(43,47) + rrt(388) * density(21) 
  pd(47,21) = pd(47,21) - rrt(388) * density(47) 
  pd(47,47) = pd(47,47) - rrt(388) * density(21) 
  pd(22,22) = pd(22,22) - rrt(389) * density(47) 
  pd(22,47) = pd(22,47) - rrt(389) * density(22) 
  pd(23,22) = pd(23,22) + rrt(389) * density(47) 
  pd(23,47) = pd(23,47) + rrt(389) * density(22) 
  pd(43,22) = pd(43,22) + rrt(389) * density(47) 
  pd(43,47) = pd(43,47) + rrt(389) * density(22) 
  pd(47,22) = pd(47,22) - rrt(389) * density(47) 
  pd(47,47) = pd(47,47) - rrt(389) * density(22) 
  pd(19,20) = pd(19,20) + rrt(390) * density(47) 
  pd(19,47) = pd(19,47) + rrt(390) * density(20) 
  pd(20,20) = pd(20,20) - rrt(390) * density(47) 
  pd(20,47) = pd(20,47) - rrt(390) * density(20) 
  pd(44,20) = pd(44,20) + rrt(390) * density(47) 
  pd(44,47) = pd(44,47) + rrt(390) * density(20) 
  pd(47,20) = pd(47,20) - rrt(390) * density(47) 
  pd(47,47) = pd(47,47) - rrt(390) * density(20) 
  pd(19,21) = pd(19,21) + rrt(391) * density(47) 
  pd(19,47) = pd(19,47) + rrt(391) * density(21) 
  pd(21,21) = pd(21,21) - rrt(391) * density(47) 
  pd(21,47) = pd(21,47) - rrt(391) * density(21) 
  pd(44,21) = pd(44,21) + rrt(391) * density(47) 
  pd(44,47) = pd(44,47) + rrt(391) * density(21) 
  pd(47,21) = pd(47,21) - rrt(391) * density(47) 
  pd(47,47) = pd(47,47) - rrt(391) * density(21) 
  pd(19,22) = pd(19,22) + rrt(392) * density(47) 
  pd(19,47) = pd(19,47) + rrt(392) * density(22) 
  pd(22,22) = pd(22,22) - rrt(392) * density(47) 
  pd(22,47) = pd(22,47) - rrt(392) * density(22) 
  pd(44,22) = pd(44,22) + rrt(392) * density(47) 
  pd(44,47) = pd(44,47) + rrt(392) * density(22) 
  pd(47,22) = pd(47,22) - rrt(392) * density(47) 
  pd(47,47) = pd(47,47) - rrt(392) * density(22) 
  pd(17,20) = pd(17,20) + rrt(393) * density(47) 
  pd(17,47) = pd(17,47) + rrt(393) * density(20) 
  pd(20,20) = pd(20,20) - rrt(393) * density(47) 
  pd(20,47) = pd(20,47) - rrt(393) * density(20) 
  pd(43,20) = pd(43,20) + rrt(393) * density(47) 
  pd(43,47) = pd(43,47) + rrt(393) * density(20) 
  pd(44,20) = pd(44,20) + rrt(393) * density(47) 
  pd(44,47) = pd(44,47) + rrt(393) * density(20) 
  pd(47,20) = pd(47,20) - rrt(393) * density(47) 
  pd(47,47) = pd(47,47) - rrt(393) * density(20) 
  pd(18,18) = pd(18,18) - rrt(394) * density(47) 
  pd(18,47) = pd(18,47) - rrt(394) * density(18) 
  pd(19,18) = pd(19,18) + rrt(394) * density(47) 
  pd(19,47) = pd(19,47) + rrt(394) * density(18) 
  pd(43,18) = pd(43,18) + rrt(394) * density(47) 
  pd(43,47) = pd(43,47) + rrt(394) * density(18) 
  pd(47,18) = pd(47,18) - rrt(394) * density(47) 
  pd(47,47) = pd(47,47) - rrt(394) * density(18) 
  pd(17,18) = pd(17,18) + rrt(395) * density(47) 
  pd(17,47) = pd(17,47) + rrt(395) * density(18) 
  pd(18,18) = pd(18,18) - rrt(395) * density(47) 
  pd(18,47) = pd(18,47) - rrt(395) * density(18) 
  pd(44,18) = pd(44,18) + rrt(395) * density(47) 
  pd(44,47) = pd(44,47) + rrt(395) * density(18) 
  pd(47,18) = pd(47,18) - rrt(395) * density(47) 
  pd(47,47) = pd(47,47) - rrt(395) * density(18) 
  pd(13,13) = pd(13,13) - rrt(396) * density(47) 
  pd(13,47) = pd(13,47) - rrt(396) * density(13) 
  pd(17,13) = pd(17,13) + rrt(396) * density(47) 
  pd(17,47) = pd(17,47) + rrt(396) * density(13) 
  pd(43,13) = pd(43,13) + rrt(396) * density(47) 
  pd(43,47) = pd(43,47) + rrt(396) * density(13) 
  pd(47,13) = pd(47,13) - rrt(396) * density(47) 
  pd(47,47) = pd(47,47) - rrt(396) * density(13) 
  pd(16,16) = pd(16,16) - rrt(397) * density(47) 
  pd(16,47) = pd(16,47) - rrt(397) * density(16) 
  pd(17,16) = pd(17,16) + rrt(397) * density(47) 
  pd(17,47) = pd(17,47) + rrt(397) * density(16) 
  pd(43,16) = pd(43,16) + rrt(397) * density(47) 
  pd(43,47) = pd(43,47) + rrt(397) * density(16) 
  pd(47,16) = pd(47,16) - rrt(397) * density(47) 
  pd(47,47) = pd(47,47) - rrt(397) * density(16) 
  pd(15,15) = pd(15,15) - rrt(398) * density(47) 
  pd(15,47) = pd(15,47) - rrt(398) * density(15) 
  pd(17,15) = pd(17,15) + rrt(398) * density(47) 
  pd(17,47) = pd(17,47) + rrt(398) * density(15) 
  pd(43,15) = pd(43,15) + rrt(398) * density(47) 
  pd(43,47) = pd(43,47) + rrt(398) * density(15) 
  pd(47,15) = pd(47,15) - rrt(398) * density(47) 
  pd(47,47) = pd(47,47) - rrt(398) * density(15) 
  pd(14,14) = pd(14,14) - rrt(399) * density(47) 
  pd(14,47) = pd(14,47) - rrt(399) * density(14) 
  pd(17,14) = pd(17,14) + rrt(399) * density(47) 
  pd(17,47) = pd(17,47) + rrt(399) * density(14) 
  pd(43,14) = pd(43,14) + rrt(399) * density(47) 
  pd(43,47) = pd(43,47) + rrt(399) * density(14) 
  pd(47,14) = pd(47,14) - rrt(399) * density(47) 
  pd(47,47) = pd(47,47) - rrt(399) * density(14) 
  pd(05,05) = pd(05,05) - rrt(400) * density(09) 
  pd(05,09) = pd(05,09) - rrt(400) * density(05) 
  pd(07,05) = pd(07,05) + rrt(400) * density(09) * 2.0d0
  pd(07,09) = pd(07,09) + rrt(400) * density(05) * 2.0d0
  pd(09,05) = pd(09,05) - rrt(400) * density(09) 
  pd(09,09) = pd(09,09) - rrt(400) * density(05) 
  pd(05,05) = pd(05,05) - rrt(401) * density(11) 
  pd(05,11) = pd(05,11) - rrt(401) * density(05) 
  pd(07,05) = pd(07,05) + rrt(401) * density(11) * 2.0d0
  pd(07,11) = pd(07,11) + rrt(401) * density(05) * 2.0d0
  pd(11,05) = pd(11,05) - rrt(401) * density(11) 
  pd(11,11) = pd(11,11) - rrt(401) * density(05) 
  pd(05,05) = pd(05,05) - rrt(402) * density(10) 
  pd(05,10) = pd(05,10) - rrt(402) * density(05) 
  pd(07,05) = pd(07,05) + rrt(402) * density(10) * 2.0d0
  pd(07,10) = pd(07,10) + rrt(402) * density(05) * 2.0d0
  pd(10,05) = pd(10,05) - rrt(402) * density(10) 
  pd(10,10) = pd(10,10) - rrt(402) * density(05) 
  pd(03,03) = pd(03,03) - rrt(403) * density(09) 
  pd(03,09) = pd(03,09) - rrt(403) * density(03) 
  pd(09,03) = pd(09,03) - rrt(403) * density(09) 
  pd(09,09) = pd(09,09) - rrt(403) * density(03) 
  pd(20,03) = pd(20,03) + rrt(403) * density(09) 
  pd(20,09) = pd(20,09) + rrt(403) * density(03) 
  pd(43,03) = pd(43,03) + rrt(403) * density(09) 
  pd(43,09) = pd(43,09) + rrt(403) * density(03) 
  pd(03,03) = pd(03,03) - rrt(404) * density(11) 
  pd(03,11) = pd(03,11) - rrt(404) * density(03) 
  pd(11,03) = pd(11,03) - rrt(404) * density(11) 
  pd(11,11) = pd(11,11) - rrt(404) * density(03) 
  pd(20,03) = pd(20,03) + rrt(404) * density(11) 
  pd(20,11) = pd(20,11) + rrt(404) * density(03) 
  pd(43,03) = pd(43,03) + rrt(404) * density(11) 
  pd(43,11) = pd(43,11) + rrt(404) * density(03) 
  pd(03,03) = pd(03,03) - rrt(405) * density(10) 
  pd(03,10) = pd(03,10) - rrt(405) * density(03) 
  pd(10,03) = pd(10,03) - rrt(405) * density(10) 
  pd(10,10) = pd(10,10) - rrt(405) * density(03) 
  pd(20,03) = pd(20,03) + rrt(405) * density(10) 
  pd(20,10) = pd(20,10) + rrt(405) * density(03) 
  pd(43,03) = pd(43,03) + rrt(405) * density(10) 
  pd(43,10) = pd(43,10) + rrt(405) * density(03) 
  pd(07,09) = pd(07,09) + rrt(406) * density(24) 
  pd(07,24) = pd(07,24) + rrt(406) * density(09) 
  pd(09,09) = pd(09,09) - rrt(406) * density(24) 
  pd(09,24) = pd(09,24) - rrt(406) * density(09) 
  pd(24,09) = pd(24,09) - rrt(406) * density(24) 
  pd(24,24) = pd(24,24) - rrt(406) * density(09) 
  pd(26,09) = pd(26,09) + rrt(406) * density(24) 
  pd(26,24) = pd(26,24) + rrt(406) * density(09) 
  pd(07,11) = pd(07,11) + rrt(407) * density(24) 
  pd(07,24) = pd(07,24) + rrt(407) * density(11) 
  pd(11,11) = pd(11,11) - rrt(407) * density(24) 
  pd(11,24) = pd(11,24) - rrt(407) * density(11) 
  pd(24,11) = pd(24,11) - rrt(407) * density(24) 
  pd(24,24) = pd(24,24) - rrt(407) * density(11) 
  pd(26,11) = pd(26,11) + rrt(407) * density(24) 
  pd(26,24) = pd(26,24) + rrt(407) * density(11) 
  pd(07,10) = pd(07,10) + rrt(408) * density(24) 
  pd(07,24) = pd(07,24) + rrt(408) * density(10) 
  pd(10,10) = pd(10,10) - rrt(408) * density(24) 
  pd(10,24) = pd(10,24) - rrt(408) * density(10) 
  pd(24,10) = pd(24,10) - rrt(408) * density(24) 
  pd(24,24) = pd(24,24) - rrt(408) * density(10) 
  pd(26,10) = pd(26,10) + rrt(408) * density(24) 
  pd(26,24) = pd(26,24) + rrt(408) * density(10) 
  pd(07,09) = pd(07,09) + rrt(409) * density(18) 
  pd(07,18) = pd(07,18) + rrt(409) * density(09) 
  pd(09,09) = pd(09,09) - rrt(409) * density(18) 
  pd(09,18) = pd(09,18) - rrt(409) * density(09) 
  pd(18,09) = pd(18,09) - rrt(409) * density(18) 
  pd(18,18) = pd(18,18) - rrt(409) * density(09) 
  pd(20,09) = pd(20,09) + rrt(409) * density(18) 
  pd(20,18) = pd(20,18) + rrt(409) * density(09) 
  pd(07,11) = pd(07,11) + rrt(410) * density(18) 
  pd(07,18) = pd(07,18) + rrt(410) * density(11) 
  pd(11,11) = pd(11,11) - rrt(410) * density(18) 
  pd(11,18) = pd(11,18) - rrt(410) * density(11) 
  pd(18,11) = pd(18,11) - rrt(410) * density(18) 
  pd(18,18) = pd(18,18) - rrt(410) * density(11) 
  pd(20,11) = pd(20,11) + rrt(410) * density(18) 
  pd(20,18) = pd(20,18) + rrt(410) * density(11) 
  pd(07,10) = pd(07,10) + rrt(411) * density(18) 
  pd(07,18) = pd(07,18) + rrt(411) * density(10) 
  pd(10,10) = pd(10,10) - rrt(411) * density(18) 
  pd(10,18) = pd(10,18) - rrt(411) * density(10) 
  pd(18,10) = pd(18,10) - rrt(411) * density(18) 
  pd(18,18) = pd(18,18) - rrt(411) * density(10) 
  pd(20,10) = pd(20,10) + rrt(411) * density(18) 
  pd(20,18) = pd(20,18) + rrt(411) * density(10) 
  pd(07,09) = pd(07,09) + rrt(412) * density(34) 
  pd(07,34) = pd(07,34) + rrt(412) * density(09) 
  pd(09,09) = pd(09,09) - rrt(412) * density(34) 
  pd(09,34) = pd(09,34) - rrt(412) * density(09) 
  pd(34,09) = pd(34,09) - rrt(412) * density(34) 
  pd(34,34) = pd(34,34) - rrt(412) * density(09) 
  pd(36,09) = pd(36,09) + rrt(412) * density(34) 
  pd(36,34) = pd(36,34) + rrt(412) * density(09) 
  pd(07,11) = pd(07,11) + rrt(413) * density(34) 
  pd(07,34) = pd(07,34) + rrt(413) * density(11) 
  pd(11,11) = pd(11,11) - rrt(413) * density(34) 
  pd(11,34) = pd(11,34) - rrt(413) * density(11) 
  pd(34,11) = pd(34,11) - rrt(413) * density(34) 
  pd(34,34) = pd(34,34) - rrt(413) * density(11) 
  pd(36,11) = pd(36,11) + rrt(413) * density(34) 
  pd(36,34) = pd(36,34) + rrt(413) * density(11) 
  pd(07,10) = pd(07,10) + rrt(414) * density(34) 
  pd(07,34) = pd(07,34) + rrt(414) * density(10) 
  pd(10,10) = pd(10,10) - rrt(414) * density(34) 
  pd(10,34) = pd(10,34) - rrt(414) * density(10) 
  pd(34,10) = pd(34,10) - rrt(414) * density(34) 
  pd(34,34) = pd(34,34) - rrt(414) * density(10) 
  pd(36,10) = pd(36,10) + rrt(414) * density(34) 
  pd(36,34) = pd(36,34) + rrt(414) * density(10) 
  pd(07,09) = pd(07,09) + rrt(415) * density(43) 
  pd(07,43) = pd(07,43) + rrt(415) * density(09) 
  pd(09,09) = pd(09,09) - rrt(415) * density(43) 
  pd(09,43) = pd(09,43) - rrt(415) * density(09) 
  pd(43,09) = pd(43,09) - rrt(415) * density(43) 
  pd(43,43) = pd(43,43) - rrt(415) * density(09) 
  pd(44,09) = pd(44,09) + rrt(415) * density(43) 
  pd(44,43) = pd(44,43) + rrt(415) * density(09) 
  pd(07,11) = pd(07,11) + rrt(416) * density(43) 
  pd(07,43) = pd(07,43) + rrt(416) * density(11) 
  pd(11,11) = pd(11,11) - rrt(416) * density(43) 
  pd(11,43) = pd(11,43) - rrt(416) * density(11) 
  pd(43,11) = pd(43,11) - rrt(416) * density(43) 
  pd(43,43) = pd(43,43) - rrt(416) * density(11) 
  pd(44,11) = pd(44,11) + rrt(416) * density(43) 
  pd(44,43) = pd(44,43) + rrt(416) * density(11) 
  pd(07,10) = pd(07,10) + rrt(417) * density(43) 
  pd(07,43) = pd(07,43) + rrt(417) * density(10) 
  pd(10,10) = pd(10,10) - rrt(417) * density(43) 
  pd(10,43) = pd(10,43) - rrt(417) * density(10) 
  pd(43,10) = pd(43,10) - rrt(417) * density(43) 
  pd(43,43) = pd(43,43) - rrt(417) * density(10) 
  pd(44,10) = pd(44,10) + rrt(417) * density(43) 
  pd(44,43) = pd(44,43) + rrt(417) * density(10) 
  pd(07,07) = pd(07,07) - rrt(418) * density(09) 
  pd(07,09) = pd(07,09) - rrt(418) * density(07) 
  pd(09,07) = pd(09,07) - rrt(418) * density(09) 
  pd(09,09) = pd(09,09) - rrt(418) * density(07) 
  pd(26,07) = pd(26,07) + rrt(418) * density(09) 
  pd(26,09) = pd(26,09) + rrt(418) * density(07) 
  pd(43,07) = pd(43,07) + rrt(418) * density(09) 
  pd(43,09) = pd(43,09) + rrt(418) * density(07) 
  pd(07,07) = pd(07,07) - rrt(419) * density(11) 
  pd(07,11) = pd(07,11) - rrt(419) * density(07) 
  pd(11,07) = pd(11,07) - rrt(419) * density(11) 
  pd(11,11) = pd(11,11) - rrt(419) * density(07) 
  pd(26,07) = pd(26,07) + rrt(419) * density(11) 
  pd(26,11) = pd(26,11) + rrt(419) * density(07) 
  pd(43,07) = pd(43,07) + rrt(419) * density(11) 
  pd(43,11) = pd(43,11) + rrt(419) * density(07) 
  pd(07,07) = pd(07,07) - rrt(420) * density(10) 
  pd(07,10) = pd(07,10) - rrt(420) * density(07) 
  pd(10,07) = pd(10,07) - rrt(420) * density(10) 
  pd(10,10) = pd(10,10) - rrt(420) * density(07) 
  pd(26,07) = pd(26,07) + rrt(420) * density(10) 
  pd(26,10) = pd(26,10) + rrt(420) * density(07) 
  pd(43,07) = pd(43,07) + rrt(420) * density(10) 
  pd(43,10) = pd(43,10) + rrt(420) * density(07) 
  pd(05,05) = pd(05,05) - rrt(421) * density(09) 
  pd(05,09) = pd(05,09) - rrt(421) * density(05) 
  pd(09,05) = pd(09,05) - rrt(421) * density(09) 
  pd(09,09) = pd(09,09) - rrt(421) * density(05) 
  pd(26,05) = pd(26,05) + rrt(421) * density(09) 
  pd(26,09) = pd(26,09) + rrt(421) * density(05) 
  pd(05,05) = pd(05,05) - rrt(422) * density(11) 
  pd(05,11) = pd(05,11) - rrt(422) * density(05) 
  pd(11,05) = pd(11,05) - rrt(422) * density(11) 
  pd(11,11) = pd(11,11) - rrt(422) * density(05) 
  pd(26,05) = pd(26,05) + rrt(422) * density(11) 
  pd(26,11) = pd(26,11) + rrt(422) * density(05) 
  pd(05,05) = pd(05,05) - rrt(423) * density(10) 
  pd(05,10) = pd(05,10) - rrt(423) * density(05) 
  pd(10,05) = pd(10,05) - rrt(423) * density(10) 
  pd(10,10) = pd(10,10) - rrt(423) * density(05) 
  pd(26,05) = pd(26,05) + rrt(423) * density(10) 
  pd(26,10) = pd(26,10) + rrt(423) * density(05) 
  pd(07,09) = pd(07,09) + rrt(424) 
  pd(09,09) = pd(09,09) - rrt(424) 
  pd(43,09) = pd(43,09) + rrt(424) 
  pd(07,11) = pd(07,11) + rrt(425) 
  pd(11,11) = pd(11,11) - rrt(425) 
  pd(43,11) = pd(43,11) + rrt(425) 
  pd(07,10) = pd(07,10) + rrt(426) 
  pd(10,10) = pd(10,10) - rrt(426) 
  pd(43,10) = pd(43,10) + rrt(426) 
  pd(05,07) = pd(05,07) + rrt(427) 
  pd(07,07) = pd(07,07) - rrt(427) 
  pd(43,07) = pd(43,07) + rrt(427) 
  pd(03,07) = pd(03,07) + rrt(428) 
  pd(07,07) = pd(07,07) - rrt(428) 
  pd(44,07) = pd(44,07) + rrt(428) 
  pd(05,07) = pd(05,07) + rrt(429) * density(24) 
  pd(05,24) = pd(05,24) + rrt(429) * density(07) 
  pd(07,07) = pd(07,07) - rrt(429) * density(24) 
  pd(07,24) = pd(07,24) - rrt(429) * density(07) 
  pd(24,07) = pd(24,07) - rrt(429) * density(24) 
  pd(24,24) = pd(24,24) - rrt(429) * density(07) 
  pd(26,07) = pd(26,07) + rrt(429) * density(24) 
  pd(26,24) = pd(26,24) + rrt(429) * density(07) 
  pd(05,05) = pd(05,05) - rrt(430) * density(05) * 4.0d0
  pd(13,05) = pd(13,05) + rrt(430) * density(05) * 2.0d0
  pd(44,05) = pd(44,05) + rrt(430) * density(05) * 2.0d0
  pd(05,05) = pd(05,05) - rrt(431) * density(24) 
  pd(05,24) = pd(05,24) - rrt(431) * density(05) 
  pd(07,05) = pd(07,05) + rrt(431) * density(24) 
  pd(07,24) = pd(07,24) + rrt(431) * density(05) 
  pd(20,05) = pd(20,05) + rrt(431) * density(24) 
  pd(20,24) = pd(20,24) + rrt(431) * density(05) 
  pd(24,05) = pd(24,05) - rrt(431) * density(24) 
  pd(24,24) = pd(24,24) - rrt(431) * density(05) 
  pd(05,05) = pd(05,05) - rrt(432) * density(18) 
  pd(05,18) = pd(05,18) - rrt(432) * density(05) 
  pd(07,05) = pd(07,05) + rrt(432) * density(18) 
  pd(07,18) = pd(07,18) + rrt(432) * density(05) 
  pd(13,05) = pd(13,05) + rrt(432) * density(18) 
  pd(13,18) = pd(13,18) + rrt(432) * density(05) 
  pd(18,05) = pd(18,05) - rrt(432) * density(18) 
  pd(18,18) = pd(18,18) - rrt(432) * density(05) 
  pd(05,05) = pd(05,05) - rrt(433) * density(36) 
  pd(05,36) = pd(05,36) - rrt(433) * density(05) 
  pd(07,05) = pd(07,05) + rrt(433) * density(36) 
  pd(07,36) = pd(07,36) + rrt(433) * density(05) 
  pd(34,05) = pd(34,05) + rrt(433) * density(36) 
  pd(34,36) = pd(34,36) + rrt(433) * density(05) 
  pd(36,05) = pd(36,05) - rrt(433) * density(36) 
  pd(36,36) = pd(36,36) - rrt(433) * density(05) 
  pd(05,05) = pd(05,05) - rrt(434) * density(37) 
  pd(05,37) = pd(05,37) - rrt(434) * density(05) 
  pd(07,05) = pd(07,05) + rrt(434) * density(37) 
  pd(07,37) = pd(07,37) + rrt(434) * density(05) 
  pd(34,05) = pd(34,05) + rrt(434) * density(37) 
  pd(34,37) = pd(34,37) + rrt(434) * density(05) 
  pd(37,05) = pd(37,05) - rrt(434) * density(37) 
  pd(37,37) = pd(37,37) - rrt(434) * density(05) 
  pd(05,05) = pd(05,05) - rrt(435) * density(38) 
  pd(05,38) = pd(05,38) - rrt(435) * density(05) 
  pd(07,05) = pd(07,05) + rrt(435) * density(38) 
  pd(07,38) = pd(07,38) + rrt(435) * density(05) 
  pd(34,05) = pd(34,05) + rrt(435) * density(38) 
  pd(34,38) = pd(34,38) + rrt(435) * density(05) 
  pd(38,05) = pd(38,05) - rrt(435) * density(38) 
  pd(38,38) = pd(38,38) - rrt(435) * density(05) 
  pd(05,05) = pd(05,05) - rrt(436) * density(34) 
  pd(05,34) = pd(05,34) - rrt(436) * density(05) 
  pd(20,05) = pd(20,05) + rrt(436) * density(34) 
  pd(20,34) = pd(20,34) + rrt(436) * density(05) 
  pd(24,05) = pd(24,05) + rrt(436) * density(34) 
  pd(24,34) = pd(24,34) + rrt(436) * density(05) 
  pd(34,05) = pd(34,05) - rrt(436) * density(34) 
  pd(34,34) = pd(34,34) - rrt(436) * density(05) 
  pd(05,05) = pd(05,05) - rrt(437) * density(34) 
  pd(05,34) = pd(05,34) - rrt(437) * density(05) 
  pd(07,05) = pd(07,05) + rrt(437) * density(34) 
  pd(07,34) = pd(07,34) + rrt(437) * density(05) 
  pd(31,05) = pd(31,05) + rrt(437) * density(34) 
  pd(31,34) = pd(31,34) + rrt(437) * density(05) 
  pd(34,05) = pd(34,05) - rrt(437) * density(34) 
  pd(34,34) = pd(34,34) - rrt(437) * density(05) 
  pd(05,05) = pd(05,05) - rrt(438) * density(44) 
  pd(05,44) = pd(05,44) - rrt(438) * density(05) 
  pd(07,05) = pd(07,05) + rrt(438) * density(44) 
  pd(07,44) = pd(07,44) + rrt(438) * density(05) 
  pd(43,05) = pd(43,05) + rrt(438) * density(44) 
  pd(43,44) = pd(43,44) + rrt(438) * density(05) 
  pd(44,05) = pd(44,05) - rrt(438) * density(44) 
  pd(44,44) = pd(44,44) - rrt(438) * density(05) 
  pd(03,05) = pd(03,05) + rrt(439) * density(43) 
  pd(03,43) = pd(03,43) + rrt(439) * density(05) 
  pd(05,05) = pd(05,05) - rrt(439) * density(43) 
  pd(05,43) = pd(05,43) - rrt(439) * density(05) 
  pd(43,05) = pd(43,05) - rrt(439) * density(43) 
  pd(43,43) = pd(43,43) - rrt(439) * density(05) 
  pd(44,05) = pd(44,05) + rrt(439) * density(43) 
  pd(44,43) = pd(44,43) + rrt(439) * density(05) 
  pd(03,05) = pd(03,05) + rrt(440) 
  pd(05,05) = pd(05,05) - rrt(440) 
  pd(43,05) = pd(43,05) + rrt(440) 
  pd(05,05) = pd(05,05) - rrt(441) * density(43) 
  pd(05,43) = pd(05,43) - rrt(441) * density(05) 
  pd(07,05) = pd(07,05) + rrt(441) * density(43) 
  pd(07,43) = pd(07,43) + rrt(441) * density(05) 
  pd(43,05) = pd(43,05) - rrt(441) * density(43) 
  pd(43,43) = pd(43,43) - rrt(441) * density(05) 
  pd(03,03) = pd(03,03) - rrt(442) * density(26) 
  pd(03,26) = pd(03,26) - rrt(442) * density(03) 
  pd(26,03) = pd(26,03) - rrt(442) * density(26) 
  pd(26,26) = pd(26,26) - rrt(442) * density(03) 
  pd(31,03) = pd(31,03) + rrt(442) * density(26) 
  pd(31,26) = pd(31,26) + rrt(442) * density(03) 
  pd(43,03) = pd(43,03) + rrt(442) * density(26) 
  pd(43,26) = pd(43,26) + rrt(442) * density(03) 
  pd(03,03) = pd(03,03) - rrt(443) * density(28) 
  pd(03,28) = pd(03,28) - rrt(443) * density(03) 
  pd(28,03) = pd(28,03) - rrt(443) * density(28) 
  pd(28,28) = pd(28,28) - rrt(443) * density(03) 
  pd(31,03) = pd(31,03) + rrt(443) * density(28) 
  pd(31,28) = pd(31,28) + rrt(443) * density(03) 
  pd(43,03) = pd(43,03) + rrt(443) * density(28) 
  pd(43,28) = pd(43,28) + rrt(443) * density(03) 
  pd(03,03) = pd(03,03) - rrt(444) * density(27) 
  pd(03,27) = pd(03,27) - rrt(444) * density(03) 
  pd(27,03) = pd(27,03) - rrt(444) * density(27) 
  pd(27,27) = pd(27,27) - rrt(444) * density(03) 
  pd(31,03) = pd(31,03) + rrt(444) * density(27) 
  pd(31,27) = pd(31,27) + rrt(444) * density(03) 
  pd(43,03) = pd(43,03) + rrt(444) * density(27) 
  pd(43,27) = pd(43,27) + rrt(444) * density(03) 
  pd(03,03) = pd(03,03) - rrt(445) * density(26) 
  pd(03,26) = pd(03,26) - rrt(445) * density(03) 
  pd(26,03) = pd(26,03) - rrt(445) * density(26) 
  pd(26,26) = pd(26,26) - rrt(445) * density(03) 
  pd(34,03) = pd(34,03) + rrt(445) * density(26) 
  pd(34,26) = pd(34,26) + rrt(445) * density(03) 
  pd(03,03) = pd(03,03) - rrt(446) * density(28) 
  pd(03,28) = pd(03,28) - rrt(446) * density(03) 
  pd(28,03) = pd(28,03) - rrt(446) * density(28) 
  pd(28,28) = pd(28,28) - rrt(446) * density(03) 
  pd(34,03) = pd(34,03) + rrt(446) * density(28) 
  pd(34,28) = pd(34,28) + rrt(446) * density(03) 
  pd(03,03) = pd(03,03) - rrt(447) * density(27) 
  pd(03,27) = pd(03,27) - rrt(447) * density(03) 
  pd(27,03) = pd(27,03) - rrt(447) * density(27) 
  pd(27,27) = pd(27,27) - rrt(447) * density(03) 
  pd(34,03) = pd(34,03) + rrt(447) * density(27) 
  pd(34,27) = pd(34,27) + rrt(447) * density(03) 
  pd(03,03) = pd(03,03) - rrt(448) * density(44) 
  pd(03,44) = pd(03,44) - rrt(448) * density(03) 
  pd(05,03) = pd(05,03) + rrt(448) * density(44) 
  pd(05,44) = pd(05,44) + rrt(448) * density(03) 
  pd(43,03) = pd(43,03) + rrt(448) * density(44) 
  pd(43,44) = pd(43,44) + rrt(448) * density(03) 
  pd(44,03) = pd(44,03) - rrt(448) * density(44) 
  pd(44,44) = pd(44,44) - rrt(448) * density(03) 
  pd(03,03) = pd(03,03) - rrt(449) * density(07) 
  pd(03,07) = pd(03,07) - rrt(449) * density(03) 
  pd(07,03) = pd(07,03) - rrt(449) * density(07) 
  pd(07,07) = pd(07,07) - rrt(449) * density(03) 
  pd(18,03) = pd(18,03) + rrt(449) * density(07) 
  pd(18,07) = pd(18,07) + rrt(449) * density(03) 
  pd(43,03) = pd(43,03) + rrt(449) * density(07) 
  pd(43,07) = pd(43,07) + rrt(449) * density(03) 
  pd(03,03) = pd(03,03) - rrt(450) * density(05) 
  pd(03,05) = pd(03,05) - rrt(450) * density(03) 
  pd(05,03) = pd(05,03) - rrt(450) * density(05) 
  pd(05,05) = pd(05,05) - rrt(450) * density(03) 
  pd(13,03) = pd(13,03) + rrt(450) * density(05) 
  pd(13,05) = pd(13,05) + rrt(450) * density(03) 
  pd(43,03) = pd(43,03) + rrt(450) * density(05) 
  pd(43,05) = pd(43,05) + rrt(450) * density(03) 
  pd(03,03) = pd(03,03) - rrt(451) * density(44) 
  pd(03,44) = pd(03,44) - rrt(451) * density(03) 
  pd(07,03) = pd(07,03) + rrt(451) * density(44) 
  pd(07,44) = pd(07,44) + rrt(451) * density(03) 
  pd(44,03) = pd(44,03) - rrt(451) * density(44) 
  pd(44,44) = pd(44,44) - rrt(451) * density(03) 
  pd(03,03) = pd(03,03) - rrt(452) * density(18) 
  pd(03,18) = pd(03,18) - rrt(452) * density(03) 
  pd(05,03) = pd(05,03) + rrt(452) * density(18) 
  pd(05,18) = pd(05,18) + rrt(452) * density(03) 
  pd(13,03) = pd(13,03) + rrt(452) * density(18) 
  pd(13,18) = pd(13,18) + rrt(452) * density(03) 
  pd(18,03) = pd(18,03) - rrt(452) * density(18) 
  pd(18,18) = pd(18,18) - rrt(452) * density(03) 
  pd(18,18) = pd(18,18) - rrt(453) * density(26) 
  pd(18,26) = pd(18,26) - rrt(453) * density(18) 
  pd(20,18) = pd(20,18) + rrt(453) * density(26) 
  pd(20,26) = pd(20,26) + rrt(453) * density(18) 
  pd(24,18) = pd(24,18) + rrt(453) * density(26) 
  pd(24,26) = pd(24,26) + rrt(453) * density(18) 
  pd(26,18) = pd(26,18) - rrt(453) * density(26) 
  pd(26,26) = pd(26,26) - rrt(453) * density(18) 
  pd(18,18) = pd(18,18) - rrt(454) * density(28) 
  pd(18,28) = pd(18,28) - rrt(454) * density(18) 
  pd(20,18) = pd(20,18) + rrt(454) * density(28) 
  pd(20,28) = pd(20,28) + rrt(454) * density(18) 
  pd(24,18) = pd(24,18) + rrt(454) * density(28) 
  pd(24,28) = pd(24,28) + rrt(454) * density(18) 
  pd(28,18) = pd(28,18) - rrt(454) * density(28) 
  pd(28,28) = pd(28,28) - rrt(454) * density(18) 
  pd(18,18) = pd(18,18) - rrt(455) * density(27) 
  pd(18,27) = pd(18,27) - rrt(455) * density(18) 
  pd(20,18) = pd(20,18) + rrt(455) * density(27) 
  pd(20,27) = pd(20,27) + rrt(455) * density(18) 
  pd(24,18) = pd(24,18) + rrt(455) * density(27) 
  pd(24,27) = pd(24,27) + rrt(455) * density(18) 
  pd(27,18) = pd(27,18) - rrt(455) * density(27) 
  pd(27,27) = pd(27,27) - rrt(455) * density(18) 
  pd(24,26) = pd(24,26) + rrt(456) * density(34) 
  pd(24,34) = pd(24,34) + rrt(456) * density(26) 
  pd(26,26) = pd(26,26) - rrt(456) * density(34) 
  pd(26,34) = pd(26,34) - rrt(456) * density(26) 
  pd(34,26) = pd(34,26) - rrt(456) * density(34) 
  pd(34,34) = pd(34,34) - rrt(456) * density(26) 
  pd(36,26) = pd(36,26) + rrt(456) * density(34) 
  pd(36,34) = pd(36,34) + rrt(456) * density(26) 
  pd(24,28) = pd(24,28) + rrt(457) * density(34) 
  pd(24,34) = pd(24,34) + rrt(457) * density(28) 
  pd(28,28) = pd(28,28) - rrt(457) * density(34) 
  pd(28,34) = pd(28,34) - rrt(457) * density(28) 
  pd(34,28) = pd(34,28) - rrt(457) * density(34) 
  pd(34,34) = pd(34,34) - rrt(457) * density(28) 
  pd(36,28) = pd(36,28) + rrt(457) * density(34) 
  pd(36,34) = pd(36,34) + rrt(457) * density(28) 
  pd(24,27) = pd(24,27) + rrt(458) * density(34) 
  pd(24,34) = pd(24,34) + rrt(458) * density(27) 
  pd(27,27) = pd(27,27) - rrt(458) * density(34) 
  pd(27,34) = pd(27,34) - rrt(458) * density(27) 
  pd(34,27) = pd(34,27) - rrt(458) * density(34) 
  pd(34,34) = pd(34,34) - rrt(458) * density(27) 
  pd(36,27) = pd(36,27) + rrt(458) * density(34) 
  pd(36,34) = pd(36,34) + rrt(458) * density(27) 
  pd(24,26) = pd(24,26) + rrt(459) * density(43) 
  pd(24,43) = pd(24,43) + rrt(459) * density(26) 
  pd(26,26) = pd(26,26) - rrt(459) * density(43) 
  pd(26,43) = pd(26,43) - rrt(459) * density(26) 
  pd(43,26) = pd(43,26) - rrt(459) * density(43) 
  pd(43,43) = pd(43,43) - rrt(459) * density(26) 
  pd(44,26) = pd(44,26) + rrt(459) * density(43) 
  pd(44,43) = pd(44,43) + rrt(459) * density(26) 
  pd(24,28) = pd(24,28) + rrt(460) * density(43) 
  pd(24,43) = pd(24,43) + rrt(460) * density(28) 
  pd(28,28) = pd(28,28) - rrt(460) * density(43) 
  pd(28,43) = pd(28,43) - rrt(460) * density(28) 
  pd(43,28) = pd(43,28) - rrt(460) * density(43) 
  pd(43,43) = pd(43,43) - rrt(460) * density(28) 
  pd(44,28) = pd(44,28) + rrt(460) * density(43) 
  pd(44,43) = pd(44,43) + rrt(460) * density(28) 
  pd(24,27) = pd(24,27) + rrt(461) * density(43) 
  pd(24,43) = pd(24,43) + rrt(461) * density(27) 
  pd(27,27) = pd(27,27) - rrt(461) * density(43) 
  pd(27,43) = pd(27,43) - rrt(461) * density(27) 
  pd(43,27) = pd(43,27) - rrt(461) * density(43) 
  pd(43,43) = pd(43,43) - rrt(461) * density(27) 
  pd(44,27) = pd(44,27) + rrt(461) * density(43) 
  pd(44,43) = pd(44,43) + rrt(461) * density(27) 
  pd(07,26) = pd(07,26) + rrt(462) * density(43) 
  pd(07,43) = pd(07,43) + rrt(462) * density(26) 
  pd(09,26) = pd(09,26) + rrt(462) * density(43) 
  pd(09,43) = pd(09,43) + rrt(462) * density(26) 
  pd(26,26) = pd(26,26) - rrt(462) * density(43) 
  pd(26,43) = pd(26,43) - rrt(462) * density(26) 
  pd(43,26) = pd(43,26) - rrt(462) * density(43) 
  pd(43,43) = pd(43,43) - rrt(462) * density(26) 
  pd(07,28) = pd(07,28) + rrt(463) * density(43) 
  pd(07,43) = pd(07,43) + rrt(463) * density(28) 
  pd(09,28) = pd(09,28) + rrt(463) * density(43) 
  pd(09,43) = pd(09,43) + rrt(463) * density(28) 
  pd(28,28) = pd(28,28) - rrt(463) * density(43) 
  pd(28,43) = pd(28,43) - rrt(463) * density(28) 
  pd(43,28) = pd(43,28) - rrt(463) * density(43) 
  pd(43,43) = pd(43,43) - rrt(463) * density(28) 
  pd(07,27) = pd(07,27) + rrt(464) * density(43) 
  pd(07,43) = pd(07,43) + rrt(464) * density(27) 
  pd(09,27) = pd(09,27) + rrt(464) * density(43) 
  pd(09,43) = pd(09,43) + rrt(464) * density(27) 
  pd(27,27) = pd(27,27) - rrt(464) * density(43) 
  pd(27,43) = pd(27,43) - rrt(464) * density(27) 
  pd(43,27) = pd(43,27) - rrt(464) * density(43) 
  pd(43,43) = pd(43,43) - rrt(464) * density(27) 
  pd(07,26) = pd(07,26) + rrt(465) * 2.0d0
  pd(26,26) = pd(26,26) - rrt(465) 
  pd(07,28) = pd(07,28) + rrt(466) * 2.0d0
  pd(28,28) = pd(28,28) - rrt(466) 
  pd(07,27) = pd(07,27) + rrt(467) * 2.0d0
  pd(27,27) = pd(27,27) - rrt(467) 
  pd(03,03) = pd(03,03) - rrt(468) * density(26) 
  pd(03,26) = pd(03,26) - rrt(468) * density(03) 
  pd(07,03) = pd(07,03) + rrt(468) * density(26) 
  pd(07,26) = pd(07,26) + rrt(468) * density(03) 
  pd(20,03) = pd(20,03) + rrt(468) * density(26) 
  pd(20,26) = pd(20,26) + rrt(468) * density(03) 
  pd(26,03) = pd(26,03) - rrt(468) * density(26) 
  pd(26,26) = pd(26,26) - rrt(468) * density(03) 
  pd(03,03) = pd(03,03) - rrt(469) * density(28) 
  pd(03,28) = pd(03,28) - rrt(469) * density(03) 
  pd(07,03) = pd(07,03) + rrt(469) * density(28) 
  pd(07,28) = pd(07,28) + rrt(469) * density(03) 
  pd(20,03) = pd(20,03) + rrt(469) * density(28) 
  pd(20,28) = pd(20,28) + rrt(469) * density(03) 
  pd(28,03) = pd(28,03) - rrt(469) * density(28) 
  pd(28,28) = pd(28,28) - rrt(469) * density(03) 
  pd(03,03) = pd(03,03) - rrt(470) * density(27) 
  pd(03,27) = pd(03,27) - rrt(470) * density(03) 
  pd(07,03) = pd(07,03) + rrt(470) * density(27) 
  pd(07,27) = pd(07,27) + rrt(470) * density(03) 
  pd(20,03) = pd(20,03) + rrt(470) * density(27) 
  pd(20,27) = pd(20,27) + rrt(470) * density(03) 
  pd(27,03) = pd(27,03) - rrt(470) * density(27) 
  pd(27,27) = pd(27,27) - rrt(470) * density(03) 
  pd(05,05) = pd(05,05) - rrt(471) * density(26) 
  pd(05,26) = pd(05,26) - rrt(471) * density(05) 
  pd(07,05) = pd(07,05) + rrt(471) * density(26) 
  pd(07,26) = pd(07,26) + rrt(471) * density(05) 
  pd(24,05) = pd(24,05) + rrt(471) * density(26) 
  pd(24,26) = pd(24,26) + rrt(471) * density(05) 
  pd(26,05) = pd(26,05) - rrt(471) * density(26) 
  pd(26,26) = pd(26,26) - rrt(471) * density(05) 
  pd(05,05) = pd(05,05) - rrt(472) * density(28) 
  pd(05,28) = pd(05,28) - rrt(472) * density(05) 
  pd(07,05) = pd(07,05) + rrt(472) * density(28) 
  pd(07,28) = pd(07,28) + rrt(472) * density(05) 
  pd(24,05) = pd(24,05) + rrt(472) * density(28) 
  pd(24,28) = pd(24,28) + rrt(472) * density(05) 
  pd(28,05) = pd(28,05) - rrt(472) * density(28) 
  pd(28,28) = pd(28,28) - rrt(472) * density(05) 
  pd(05,05) = pd(05,05) - rrt(473) * density(27) 
  pd(05,27) = pd(05,27) - rrt(473) * density(05) 
  pd(07,05) = pd(07,05) + rrt(473) * density(27) 
  pd(07,27) = pd(07,27) + rrt(473) * density(05) 
  pd(24,05) = pd(24,05) + rrt(473) * density(27) 
  pd(24,27) = pd(24,27) + rrt(473) * density(05) 
  pd(27,05) = pd(27,05) - rrt(473) * density(27) 
  pd(27,27) = pd(27,27) - rrt(473) * density(05) 
  pd(20,24) = pd(20,24) + rrt(474) * density(24) * 2.0d0
  pd(24,24) = pd(24,24) - rrt(474) * density(24) * 4.0d0
  pd(26,24) = pd(26,24) + rrt(474) * density(24) * 2.0d0
  pd(18,20) = pd(18,20) + rrt(475) * density(24) 
  pd(18,24) = pd(18,24) + rrt(475) * density(20) 
  pd(20,20) = pd(20,20) - rrt(475) * density(24) 
  pd(20,24) = pd(20,24) - rrt(475) * density(20) 
  pd(24,20) = pd(24,20) - rrt(475) * density(24) 
  pd(24,24) = pd(24,24) - rrt(475) * density(20) 
  pd(26,20) = pd(26,20) + rrt(475) * density(24) 
  pd(26,24) = pd(26,24) + rrt(475) * density(20) 
  pd(18,21) = pd(18,21) + rrt(476) * density(24) 
  pd(18,24) = pd(18,24) + rrt(476) * density(21) 
  pd(21,21) = pd(21,21) - rrt(476) * density(24) 
  pd(21,24) = pd(21,24) - rrt(476) * density(21) 
  pd(24,21) = pd(24,21) - rrt(476) * density(24) 
  pd(24,24) = pd(24,24) - rrt(476) * density(21) 
  pd(26,21) = pd(26,21) + rrt(476) * density(24) 
  pd(26,24) = pd(26,24) + rrt(476) * density(21) 
  pd(18,22) = pd(18,22) + rrt(477) * density(24) 
  pd(18,24) = pd(18,24) + rrt(477) * density(22) 
  pd(22,22) = pd(22,22) - rrt(477) * density(24) 
  pd(22,24) = pd(22,24) - rrt(477) * density(22) 
  pd(24,22) = pd(24,22) - rrt(477) * density(24) 
  pd(24,24) = pd(24,24) - rrt(477) * density(22) 
  pd(26,22) = pd(26,22) + rrt(477) * density(24) 
  pd(26,24) = pd(26,24) + rrt(477) * density(22) 
  pd(24,24) = pd(24,24) - rrt(478) * density(36) 
  pd(24,36) = pd(24,36) - rrt(478) * density(24) 
  pd(26,24) = pd(26,24) + rrt(478) * density(36) 
  pd(26,36) = pd(26,36) + rrt(478) * density(24) 
  pd(34,24) = pd(34,24) + rrt(478) * density(36) 
  pd(34,36) = pd(34,36) + rrt(478) * density(24) 
  pd(36,24) = pd(36,24) - rrt(478) * density(36) 
  pd(36,36) = pd(36,36) - rrt(478) * density(24) 
  pd(24,24) = pd(24,24) - rrt(479) * density(37) 
  pd(24,37) = pd(24,37) - rrt(479) * density(24) 
  pd(26,24) = pd(26,24) + rrt(479) * density(37) 
  pd(26,37) = pd(26,37) + rrt(479) * density(24) 
  pd(34,24) = pd(34,24) + rrt(479) * density(37) 
  pd(34,37) = pd(34,37) + rrt(479) * density(24) 
  pd(37,24) = pd(37,24) - rrt(479) * density(37) 
  pd(37,37) = pd(37,37) - rrt(479) * density(24) 
  pd(24,24) = pd(24,24) - rrt(480) * density(38) 
  pd(24,38) = pd(24,38) - rrt(480) * density(24) 
  pd(26,24) = pd(26,24) + rrt(480) * density(38) 
  pd(26,38) = pd(26,38) + rrt(480) * density(24) 
  pd(34,24) = pd(34,24) + rrt(480) * density(38) 
  pd(34,38) = pd(34,38) + rrt(480) * density(24) 
  pd(38,24) = pd(38,24) - rrt(480) * density(38) 
  pd(38,38) = pd(38,38) - rrt(480) * density(24) 
  pd(20,24) = pd(20,24) + rrt(481) * density(34) 
  pd(20,34) = pd(20,34) + rrt(481) * density(24) 
  pd(24,24) = pd(24,24) - rrt(481) * density(34) 
  pd(24,34) = pd(24,34) - rrt(481) * density(24) 
  pd(34,24) = pd(34,24) - rrt(481) * density(34) 
  pd(34,34) = pd(34,34) - rrt(481) * density(24) 
  pd(36,24) = pd(36,24) + rrt(481) * density(34) 
  pd(36,34) = pd(36,34) + rrt(481) * density(24) 
  pd(24,24) = pd(24,24) - rrt(482) * density(34) 
  pd(24,34) = pd(24,34) - rrt(482) * density(24) 
  pd(26,24) = pd(26,24) + rrt(482) * density(34) 
  pd(26,34) = pd(26,34) + rrt(482) * density(24) 
  pd(31,24) = pd(31,24) + rrt(482) * density(34) 
  pd(31,34) = pd(31,34) + rrt(482) * density(24) 
  pd(34,24) = pd(34,24) - rrt(482) * density(34) 
  pd(34,34) = pd(34,34) - rrt(482) * density(24) 
  pd(24,24) = pd(24,24) - rrt(483) * density(31) 
  pd(24,31) = pd(24,31) - rrt(483) * density(24) 
  pd(26,24) = pd(26,24) + rrt(483) * density(31) 
  pd(26,31) = pd(26,31) + rrt(483) * density(24) 
  pd(30,24) = pd(30,24) + rrt(483) * density(31) 
  pd(30,31) = pd(30,31) + rrt(483) * density(24) 
  pd(31,24) = pd(31,24) - rrt(483) * density(31) 
  pd(31,31) = pd(31,31) - rrt(483) * density(24) 
  pd(24,24) = pd(24,24) - rrt(484) * density(32) 
  pd(24,32) = pd(24,32) - rrt(484) * density(24) 
  pd(26,24) = pd(26,24) + rrt(484) * density(32) 
  pd(26,32) = pd(26,32) + rrt(484) * density(24) 
  pd(30,24) = pd(30,24) + rrt(484) * density(32) 
  pd(30,32) = pd(30,32) + rrt(484) * density(24) 
  pd(32,24) = pd(32,24) - rrt(484) * density(32) 
  pd(32,32) = pd(32,32) - rrt(484) * density(24) 
  pd(24,24) = pd(24,24) - rrt(485) * density(44) 
  pd(24,44) = pd(24,44) - rrt(485) * density(24) 
  pd(26,24) = pd(26,24) + rrt(485) * density(44) 
  pd(26,44) = pd(26,44) + rrt(485) * density(24) 
  pd(43,24) = pd(43,24) + rrt(485) * density(44) 
  pd(43,44) = pd(43,44) + rrt(485) * density(24) 
  pd(44,24) = pd(44,24) - rrt(485) * density(44) 
  pd(44,44) = pd(44,44) - rrt(485) * density(24) 
  pd(07,24) = pd(07,24) + rrt(486) * density(43) * 2.0d0
  pd(07,43) = pd(07,43) + rrt(486) * density(24) * 2.0d0
  pd(24,24) = pd(24,24) - rrt(486) * density(43) 
  pd(24,43) = pd(24,43) - rrt(486) * density(24) 
  pd(43,24) = pd(43,24) - rrt(486) * density(43) 
  pd(43,43) = pd(43,43) - rrt(486) * density(24) 
  pd(20,24) = pd(20,24) + rrt(487) * density(43) 
  pd(20,43) = pd(20,43) + rrt(487) * density(24) 
  pd(24,24) = pd(24,24) - rrt(487) * density(43) 
  pd(24,43) = pd(24,43) - rrt(487) * density(24) 
  pd(43,24) = pd(43,24) - rrt(487) * density(43) 
  pd(43,43) = pd(43,43) - rrt(487) * density(24) 
  pd(44,24) = pd(44,24) + rrt(487) * density(43) 
  pd(44,43) = pd(44,43) + rrt(487) * density(24) 
  pd(24,24) = pd(24,24) - rrt(488) * density(43) 
  pd(24,43) = pd(24,43) - rrt(488) * density(24) 
  pd(26,24) = pd(26,24) + rrt(488) * density(43) 
  pd(26,43) = pd(26,43) + rrt(488) * density(24) 
  pd(43,24) = pd(43,24) - rrt(488) * density(43) 
  pd(43,43) = pd(43,43) - rrt(488) * density(24) 
  pd(20,24) = pd(20,24) + rrt(489) 
  pd(24,24) = pd(24,24) - rrt(489) 
  pd(43,24) = pd(43,24) + rrt(489) 
  pd(24,24) = pd(24,24) - rrt(490) * density(24) * 4.0d0
  pd(41,24) = pd(41,24) + rrt(490) * density(24) * 2.0d0
  pd(18,18) = pd(18,18) - rrt(491) * density(24) 
  pd(18,24) = pd(18,24) - rrt(491) * density(18) 
  pd(20,18) = pd(20,18) + rrt(491) * density(24) * 2.0d0
  pd(20,24) = pd(20,24) + rrt(491) * density(18) * 2.0d0
  pd(24,18) = pd(24,18) - rrt(491) * density(24) 
  pd(24,24) = pd(24,24) - rrt(491) * density(18) 
  pd(18,20) = pd(18,20) + rrt(492) * density(43) 
  pd(18,43) = pd(18,43) + rrt(492) * density(20) 
  pd(20,20) = pd(20,20) - rrt(492) * density(43) 
  pd(20,43) = pd(20,43) - rrt(492) * density(20) 
  pd(43,20) = pd(43,20) - rrt(492) * density(43) 
  pd(43,43) = pd(43,43) - rrt(492) * density(20) 
  pd(44,20) = pd(44,20) + rrt(492) * density(43) 
  pd(44,43) = pd(44,43) + rrt(492) * density(20) 
  pd(18,21) = pd(18,21) + rrt(493) * density(43) 
  pd(18,43) = pd(18,43) + rrt(493) * density(21) 
  pd(21,21) = pd(21,21) - rrt(493) * density(43) 
  pd(21,43) = pd(21,43) - rrt(493) * density(21) 
  pd(43,21) = pd(43,21) - rrt(493) * density(43) 
  pd(43,43) = pd(43,43) - rrt(493) * density(21) 
  pd(44,21) = pd(44,21) + rrt(493) * density(43) 
  pd(44,43) = pd(44,43) + rrt(493) * density(21) 
  pd(18,22) = pd(18,22) + rrt(494) * density(43) 
  pd(18,43) = pd(18,43) + rrt(494) * density(22) 
  pd(22,22) = pd(22,22) - rrt(494) * density(43) 
  pd(22,43) = pd(22,43) - rrt(494) * density(22) 
  pd(43,22) = pd(43,22) - rrt(494) * density(43) 
  pd(43,43) = pd(43,43) - rrt(494) * density(22) 
  pd(44,22) = pd(44,22) + rrt(494) * density(43) 
  pd(44,43) = pd(44,43) + rrt(494) * density(22) 
  pd(20,20) = pd(20,20) - rrt(495) * density(43) 
  pd(20,43) = pd(20,43) - rrt(495) * density(20) 
  pd(24,20) = pd(24,20) + rrt(495) * density(43) 
  pd(24,43) = pd(24,43) + rrt(495) * density(20) 
  pd(43,20) = pd(43,20) - rrt(495) * density(43) 
  pd(43,43) = pd(43,43) - rrt(495) * density(20) 
  pd(21,21) = pd(21,21) - rrt(496) * density(43) 
  pd(21,43) = pd(21,43) - rrt(496) * density(21) 
  pd(24,21) = pd(24,21) + rrt(496) * density(43) 
  pd(24,43) = pd(24,43) + rrt(496) * density(21) 
  pd(43,21) = pd(43,21) - rrt(496) * density(43) 
  pd(43,43) = pd(43,43) - rrt(496) * density(21) 
  pd(22,22) = pd(22,22) - rrt(497) * density(43) 
  pd(22,43) = pd(22,43) - rrt(497) * density(22) 
  pd(24,22) = pd(24,22) + rrt(497) * density(43) 
  pd(24,43) = pd(24,43) + rrt(497) * density(22) 
  pd(43,22) = pd(43,22) - rrt(497) * density(43) 
  pd(43,43) = pd(43,43) - rrt(497) * density(22) 
  pd(20,20) = pd(20,20) - rrt(498) * density(44) 
  pd(20,44) = pd(20,44) - rrt(498) * density(20) 
  pd(24,20) = pd(24,20) + rrt(498) * density(44) 
  pd(24,44) = pd(24,44) + rrt(498) * density(20) 
  pd(43,20) = pd(43,20) + rrt(498) * density(44) 
  pd(43,44) = pd(43,44) + rrt(498) * density(20) 
  pd(44,20) = pd(44,20) - rrt(498) * density(44) 
  pd(44,44) = pd(44,44) - rrt(498) * density(20) 
  pd(18,20) = pd(18,20) + rrt(499) 
  pd(20,20) = pd(20,20) - rrt(499) 
  pd(43,20) = pd(43,20) + rrt(499) 
  pd(18,21) = pd(18,21) + rrt(500) 
  pd(21,21) = pd(21,21) - rrt(500) 
  pd(43,21) = pd(43,21) + rrt(500) 
  pd(18,22) = pd(18,22) + rrt(501) 
  pd(22,22) = pd(22,22) - rrt(501) 
  pd(43,22) = pd(43,22) + rrt(501) 
  pd(13,13) = pd(13,13) - rrt(502) * density(20) 
  pd(13,20) = pd(13,20) - rrt(502) * density(13) 
  pd(18,13) = pd(18,13) + rrt(502) * density(20) * 2.0d0
  pd(18,20) = pd(18,20) + rrt(502) * density(13) * 2.0d0
  pd(20,13) = pd(20,13) - rrt(502) * density(20) 
  pd(20,20) = pd(20,20) - rrt(502) * density(13) 
  pd(13,13) = pd(13,13) - rrt(503) * density(21) 
  pd(13,21) = pd(13,21) - rrt(503) * density(13) 
  pd(18,13) = pd(18,13) + rrt(503) * density(21) * 2.0d0
  pd(18,21) = pd(18,21) + rrt(503) * density(13) * 2.0d0
  pd(21,13) = pd(21,13) - rrt(503) * density(21) 
  pd(21,21) = pd(21,21) - rrt(503) * density(13) 
  pd(13,13) = pd(13,13) - rrt(504) * density(22) 
  pd(13,22) = pd(13,22) - rrt(504) * density(13) 
  pd(18,13) = pd(18,13) + rrt(504) * density(22) * 2.0d0
  pd(18,22) = pd(18,22) + rrt(504) * density(13) * 2.0d0
  pd(22,13) = pd(22,13) - rrt(504) * density(22) 
  pd(22,22) = pd(22,22) - rrt(504) * density(13) 
  pd(16,16) = pd(16,16) - rrt(505) * density(20) 
  pd(16,20) = pd(16,20) - rrt(505) * density(16) 
  pd(18,16) = pd(18,16) + rrt(505) * density(20) * 2.0d0
  pd(18,20) = pd(18,20) + rrt(505) * density(16) * 2.0d0
  pd(20,16) = pd(20,16) - rrt(505) * density(20) 
  pd(20,20) = pd(20,20) - rrt(505) * density(16) 
  pd(15,15) = pd(15,15) - rrt(506) * density(20) 
  pd(15,20) = pd(15,20) - rrt(506) * density(15) 
  pd(18,15) = pd(18,15) + rrt(506) * density(20) * 2.0d0
  pd(18,20) = pd(18,20) + rrt(506) * density(15) * 2.0d0
  pd(20,15) = pd(20,15) - rrt(506) * density(20) 
  pd(20,20) = pd(20,20) - rrt(506) * density(15) 
  pd(14,14) = pd(14,14) - rrt(507) * density(20) 
  pd(14,20) = pd(14,20) - rrt(507) * density(14) 
  pd(18,14) = pd(18,14) + rrt(507) * density(20) * 2.0d0
  pd(18,20) = pd(18,20) + rrt(507) * density(14) * 2.0d0
  pd(20,14) = pd(20,14) - rrt(507) * density(20) 
  pd(20,20) = pd(20,20) - rrt(507) * density(14) 
  pd(18,20) = pd(18,20) + rrt(508) * density(31) 
  pd(18,31) = pd(18,31) + rrt(508) * density(20) 
  pd(20,20) = pd(20,20) - rrt(508) * density(31) 
  pd(20,31) = pd(20,31) - rrt(508) * density(20) 
  pd(31,20) = pd(31,20) - rrt(508) * density(31) 
  pd(31,31) = pd(31,31) - rrt(508) * density(20) 
  pd(34,20) = pd(34,20) + rrt(508) * density(31) 
  pd(34,31) = pd(34,31) + rrt(508) * density(20) 
  pd(18,21) = pd(18,21) + rrt(509) * density(31) 
  pd(18,31) = pd(18,31) + rrt(509) * density(21) 
  pd(21,21) = pd(21,21) - rrt(509) * density(31) 
  pd(21,31) = pd(21,31) - rrt(509) * density(21) 
  pd(31,21) = pd(31,21) - rrt(509) * density(31) 
  pd(31,31) = pd(31,31) - rrt(509) * density(21) 
  pd(34,21) = pd(34,21) + rrt(509) * density(31) 
  pd(34,31) = pd(34,31) + rrt(509) * density(21) 
  pd(18,22) = pd(18,22) + rrt(510) * density(31) 
  pd(18,31) = pd(18,31) + rrt(510) * density(22) 
  pd(22,22) = pd(22,22) - rrt(510) * density(31) 
  pd(22,31) = pd(22,31) - rrt(510) * density(22) 
  pd(31,22) = pd(31,22) - rrt(510) * density(31) 
  pd(31,31) = pd(31,31) - rrt(510) * density(22) 
  pd(34,22) = pd(34,22) + rrt(510) * density(31) 
  pd(34,31) = pd(34,31) + rrt(510) * density(22) 
  pd(18,20) = pd(18,20) + rrt(511) * density(32) 
  pd(18,32) = pd(18,32) + rrt(511) * density(20) 
  pd(20,20) = pd(20,20) - rrt(511) * density(32) 
  pd(20,32) = pd(20,32) - rrt(511) * density(20) 
  pd(32,20) = pd(32,20) - rrt(511) * density(32) 
  pd(32,32) = pd(32,32) - rrt(511) * density(20) 
  pd(34,20) = pd(34,20) + rrt(511) * density(32) 
  pd(34,32) = pd(34,32) + rrt(511) * density(20) 
  pd(18,20) = pd(18,20) + rrt(512) * density(20) * 2.0d0
  pd(20,20) = pd(20,20) - rrt(512) * density(20) * 4.0d0
  pd(24,20) = pd(24,20) + rrt(512) * density(20) * 2.0d0
  pd(18,21) = pd(18,21) + rrt(513) * density(21) * 2.0d0
  pd(21,21) = pd(21,21) - rrt(513) * density(21) * 4.0d0
  pd(24,21) = pd(24,21) + rrt(513) * density(21) * 2.0d0
  pd(18,22) = pd(18,22) + rrt(514) * density(22) * 2.0d0
  pd(22,22) = pd(22,22) - rrt(514) * density(22) * 4.0d0
  pd(24,22) = pd(24,22) + rrt(514) * density(22) * 2.0d0
  pd(07,07) = pd(07,07) - rrt(515) * density(20) 
  pd(07,20) = pd(07,20) - rrt(515) * density(07) 
  pd(20,07) = pd(20,07) - rrt(515) * density(20) 
  pd(20,20) = pd(20,20) - rrt(515) * density(07) 
  pd(34,07) = pd(34,07) + rrt(515) * density(20) 
  pd(34,20) = pd(34,20) + rrt(515) * density(07) 
  pd(07,07) = pd(07,07) - rrt(516) * density(21) 
  pd(07,21) = pd(07,21) - rrt(516) * density(07) 
  pd(21,07) = pd(21,07) - rrt(516) * density(21) 
  pd(21,21) = pd(21,21) - rrt(516) * density(07) 
  pd(34,07) = pd(34,07) + rrt(516) * density(21) 
  pd(34,21) = pd(34,21) + rrt(516) * density(07) 
  pd(07,07) = pd(07,07) - rrt(517) * density(22) 
  pd(07,22) = pd(07,22) - rrt(517) * density(07) 
  pd(22,07) = pd(22,07) - rrt(517) * density(22) 
  pd(22,22) = pd(22,22) - rrt(517) * density(07) 
  pd(34,07) = pd(34,07) + rrt(517) * density(22) 
  pd(34,22) = pd(34,22) + rrt(517) * density(07) 
  pd(13,20) = pd(13,20) + rrt(518) 
  pd(20,20) = pd(20,20) - rrt(518) 
  pd(44,20) = pd(44,20) + rrt(518) 
  pd(13,21) = pd(13,21) + rrt(519) 
  pd(21,21) = pd(21,21) - rrt(519) 
  pd(44,21) = pd(44,21) + rrt(519) 
  pd(13,22) = pd(13,22) + rrt(520) 
  pd(22,22) = pd(22,22) - rrt(520) 
  pd(44,22) = pd(44,22) + rrt(520) 
  pd(20,20) = pd(20,20) - rrt(521) * density(44) 
  pd(20,44) = pd(20,44) - rrt(521) * density(20) 
  pd(26,20) = pd(26,20) + rrt(521) * density(44) 
  pd(26,44) = pd(26,44) + rrt(521) * density(20) 
  pd(44,20) = pd(44,20) - rrt(521) * density(44) 
  pd(44,44) = pd(44,44) - rrt(521) * density(20) 
  pd(21,21) = pd(21,21) - rrt(522) * density(44) 
  pd(21,44) = pd(21,44) - rrt(522) * density(21) 
  pd(26,21) = pd(26,21) + rrt(522) * density(44) 
  pd(26,44) = pd(26,44) + rrt(522) * density(21) 
  pd(44,21) = pd(44,21) - rrt(522) * density(44) 
  pd(44,44) = pd(44,44) - rrt(522) * density(21) 
  pd(22,22) = pd(22,22) - rrt(523) * density(44) 
  pd(22,44) = pd(22,44) - rrt(523) * density(22) 
  pd(26,22) = pd(26,22) + rrt(523) * density(44) 
  pd(26,44) = pd(26,44) + rrt(523) * density(22) 
  pd(44,22) = pd(44,22) - rrt(523) * density(44) 
  pd(44,44) = pd(44,44) - rrt(523) * density(22) 
  pd(05,05) = pd(05,05) - rrt(524) * density(20) 
  pd(05,20) = pd(05,20) - rrt(524) * density(05) 
  pd(20,05) = pd(20,05) - rrt(524) * density(20) 
  pd(20,20) = pd(20,20) - rrt(524) * density(05) 
  pd(31,05) = pd(31,05) + rrt(524) * density(20) 
  pd(31,20) = pd(31,20) + rrt(524) * density(05) 
  pd(05,05) = pd(05,05) - rrt(525) * density(21) 
  pd(05,21) = pd(05,21) - rrt(525) * density(05) 
  pd(21,05) = pd(21,05) - rrt(525) * density(21) 
  pd(21,21) = pd(21,21) - rrt(525) * density(05) 
  pd(31,05) = pd(31,05) + rrt(525) * density(21) 
  pd(31,21) = pd(31,21) + rrt(525) * density(05) 
  pd(05,05) = pd(05,05) - rrt(526) * density(22) 
  pd(05,22) = pd(05,22) - rrt(526) * density(05) 
  pd(22,05) = pd(22,05) - rrt(526) * density(22) 
  pd(22,22) = pd(22,22) - rrt(526) * density(05) 
  pd(31,05) = pd(31,05) + rrt(526) * density(22) 
  pd(31,22) = pd(31,22) + rrt(526) * density(05) 
  pd(13,18) = pd(13,18) + rrt(527) * density(18) * 2.0d0
  pd(18,18) = pd(18,18) - rrt(527) * density(18) * 4.0d0
  pd(20,18) = pd(20,18) + rrt(527) * density(18) * 2.0d0
  pd(18,18) = pd(18,18) - rrt(528) * density(36) 
  pd(18,36) = pd(18,36) - rrt(528) * density(18) 
  pd(20,18) = pd(20,18) + rrt(528) * density(36) 
  pd(20,36) = pd(20,36) + rrt(528) * density(18) 
  pd(34,18) = pd(34,18) + rrt(528) * density(36) 
  pd(34,36) = pd(34,36) + rrt(528) * density(18) 
  pd(36,18) = pd(36,18) - rrt(528) * density(36) 
  pd(36,36) = pd(36,36) - rrt(528) * density(18) 
  pd(18,18) = pd(18,18) - rrt(529) * density(37) 
  pd(18,37) = pd(18,37) - rrt(529) * density(18) 
  pd(20,18) = pd(20,18) + rrt(529) * density(37) 
  pd(20,37) = pd(20,37) + rrt(529) * density(18) 
  pd(34,18) = pd(34,18) + rrt(529) * density(37) 
  pd(34,37) = pd(34,37) + rrt(529) * density(18) 
  pd(37,18) = pd(37,18) - rrt(529) * density(37) 
  pd(37,37) = pd(37,37) - rrt(529) * density(18) 
  pd(18,18) = pd(18,18) - rrt(530) * density(38) 
  pd(18,38) = pd(18,38) - rrt(530) * density(18) 
  pd(20,18) = pd(20,18) + rrt(530) * density(38) 
  pd(20,38) = pd(20,38) + rrt(530) * density(18) 
  pd(34,18) = pd(34,18) + rrt(530) * density(38) 
  pd(34,38) = pd(34,38) + rrt(530) * density(18) 
  pd(38,18) = pd(38,18) - rrt(530) * density(38) 
  pd(38,38) = pd(38,38) - rrt(530) * density(18) 
  pd(13,18) = pd(13,18) + rrt(531) * density(34) 
  pd(13,34) = pd(13,34) + rrt(531) * density(18) 
  pd(18,18) = pd(18,18) - rrt(531) * density(34) 
  pd(18,34) = pd(18,34) - rrt(531) * density(18) 
  pd(34,18) = pd(34,18) - rrt(531) * density(34) 
  pd(34,34) = pd(34,34) - rrt(531) * density(18) 
  pd(36,18) = pd(36,18) + rrt(531) * density(34) 
  pd(36,34) = pd(36,34) + rrt(531) * density(18) 
  pd(18,18) = pd(18,18) - rrt(532) * density(34) 
  pd(18,34) = pd(18,34) - rrt(532) * density(18) 
  pd(20,18) = pd(20,18) + rrt(532) * density(34) 
  pd(20,34) = pd(20,34) + rrt(532) * density(18) 
  pd(31,18) = pd(31,18) + rrt(532) * density(34) 
  pd(31,34) = pd(31,34) + rrt(532) * density(18) 
  pd(34,18) = pd(34,18) - rrt(532) * density(34) 
  pd(34,34) = pd(34,34) - rrt(532) * density(18) 
  pd(18,18) = pd(18,18) - rrt(533) * density(44) 
  pd(18,44) = pd(18,44) - rrt(533) * density(18) 
  pd(20,18) = pd(20,18) + rrt(533) * density(44) 
  pd(20,44) = pd(20,44) + rrt(533) * density(18) 
  pd(43,18) = pd(43,18) + rrt(533) * density(44) 
  pd(43,44) = pd(43,44) + rrt(533) * density(18) 
  pd(44,18) = pd(44,18) - rrt(533) * density(44) 
  pd(44,44) = pd(44,44) - rrt(533) * density(18) 
  pd(13,18) = pd(13,18) + rrt(534) * density(43) 
  pd(13,43) = pd(13,43) + rrt(534) * density(18) 
  pd(18,18) = pd(18,18) - rrt(534) * density(43) 
  pd(18,43) = pd(18,43) - rrt(534) * density(18) 
  pd(43,18) = pd(43,18) - rrt(534) * density(43) 
  pd(43,43) = pd(43,43) - rrt(534) * density(18) 
  pd(44,18) = pd(44,18) + rrt(534) * density(43) 
  pd(44,43) = pd(44,43) + rrt(534) * density(18) 
  pd(18,18) = pd(18,18) - rrt(535) * density(43) 
  pd(18,43) = pd(18,43) - rrt(535) * density(18) 
  pd(20,18) = pd(20,18) + rrt(535) * density(43) 
  pd(20,43) = pd(20,43) + rrt(535) * density(18) 
  pd(43,18) = pd(43,18) - rrt(535) * density(43) 
  pd(43,43) = pd(43,43) - rrt(535) * density(18) 
  pd(13,18) = pd(13,18) + rrt(536) 
  pd(18,18) = pd(18,18) - rrt(536) 
  pd(43,18) = pd(43,18) + rrt(536) 
  pd(13,13) = pd(13,13) - rrt(537) * density(43) 
  pd(13,43) = pd(13,43) - rrt(537) * density(13) 
  pd(18,13) = pd(18,13) + rrt(537) * density(43) 
  pd(18,43) = pd(18,43) + rrt(537) * density(13) 
  pd(43,13) = pd(43,13) - rrt(537) * density(43) 
  pd(43,43) = pd(43,43) - rrt(537) * density(13) 
  pd(16,16) = pd(16,16) - rrt(538) * density(43) 
  pd(16,43) = pd(16,43) - rrt(538) * density(16) 
  pd(18,16) = pd(18,16) + rrt(538) * density(43) 
  pd(18,43) = pd(18,43) + rrt(538) * density(16) 
  pd(43,16) = pd(43,16) - rrt(538) * density(43) 
  pd(43,43) = pd(43,43) - rrt(538) * density(16) 
  pd(15,15) = pd(15,15) - rrt(539) * density(43) 
  pd(15,43) = pd(15,43) - rrt(539) * density(15) 
  pd(18,15) = pd(18,15) + rrt(539) * density(43) 
  pd(18,43) = pd(18,43) + rrt(539) * density(15) 
  pd(43,15) = pd(43,15) - rrt(539) * density(43) 
  pd(43,43) = pd(43,43) - rrt(539) * density(15) 
  pd(14,14) = pd(14,14) - rrt(540) * density(43) 
  pd(14,43) = pd(14,43) - rrt(540) * density(14) 
  pd(18,14) = pd(18,14) + rrt(540) * density(43) 
  pd(18,43) = pd(18,43) + rrt(540) * density(14) 
  pd(43,14) = pd(43,14) - rrt(540) * density(43) 
  pd(43,43) = pd(43,43) - rrt(540) * density(14) 
  pd(13,13) = pd(13,13) - rrt(541) * density(44) 
  pd(13,44) = pd(13,44) - rrt(541) * density(13) 
  pd(20,13) = pd(20,13) + rrt(541) * density(44) 
  pd(20,44) = pd(20,44) + rrt(541) * density(13) 
  pd(44,13) = pd(44,13) - rrt(541) * density(44) 
  pd(44,44) = pd(44,44) - rrt(541) * density(13) 
  pd(16,16) = pd(16,16) - rrt(542) * density(44) 
  pd(16,44) = pd(16,44) - rrt(542) * density(16) 
  pd(20,16) = pd(20,16) + rrt(542) * density(44) 
  pd(20,44) = pd(20,44) + rrt(542) * density(16) 
  pd(44,16) = pd(44,16) - rrt(542) * density(44) 
  pd(44,44) = pd(44,44) - rrt(542) * density(16) 
  pd(15,15) = pd(15,15) - rrt(543) * density(44) 
  pd(15,44) = pd(15,44) - rrt(543) * density(15) 
  pd(20,15) = pd(20,15) + rrt(543) * density(44) 
  pd(20,44) = pd(20,44) + rrt(543) * density(15) 
  pd(44,15) = pd(44,15) - rrt(543) * density(44) 
  pd(44,44) = pd(44,44) - rrt(543) * density(15) 
  pd(14,14) = pd(14,14) - rrt(544) * density(44) 
  pd(14,44) = pd(14,44) - rrt(544) * density(14) 
  pd(20,14) = pd(20,14) + rrt(544) * density(44) 
  pd(20,44) = pd(20,44) + rrt(544) * density(14) 
  pd(44,14) = pd(44,14) - rrt(544) * density(44) 
  pd(44,44) = pd(44,44) - rrt(544) * density(14) 
  pd(13,13) = pd(13,13) - rrt(545) * density(44) 
  pd(13,44) = pd(13,44) - rrt(545) * density(13) 
  pd(18,13) = pd(18,13) + rrt(545) * density(44) 
  pd(18,44) = pd(18,44) + rrt(545) * density(13) 
  pd(43,13) = pd(43,13) + rrt(545) * density(44) 
  pd(43,44) = pd(43,44) + rrt(545) * density(13) 
  pd(44,13) = pd(44,13) - rrt(545) * density(44) 
  pd(44,44) = pd(44,44) - rrt(545) * density(13) 
  pd(16,16) = pd(16,16) - rrt(546) * density(44) 
  pd(16,44) = pd(16,44) - rrt(546) * density(16) 
  pd(18,16) = pd(18,16) + rrt(546) * density(44) 
  pd(18,44) = pd(18,44) + rrt(546) * density(16) 
  pd(43,16) = pd(43,16) + rrt(546) * density(44) 
  pd(43,44) = pd(43,44) + rrt(546) * density(16) 
  pd(44,16) = pd(44,16) - rrt(546) * density(44) 
  pd(44,44) = pd(44,44) - rrt(546) * density(16) 
  pd(15,15) = pd(15,15) - rrt(547) * density(44) 
  pd(15,44) = pd(15,44) - rrt(547) * density(15) 
  pd(18,15) = pd(18,15) + rrt(547) * density(44) 
  pd(18,44) = pd(18,44) + rrt(547) * density(15) 
  pd(43,15) = pd(43,15) + rrt(547) * density(44) 
  pd(43,44) = pd(43,44) + rrt(547) * density(15) 
  pd(44,15) = pd(44,15) - rrt(547) * density(44) 
  pd(44,44) = pd(44,44) - rrt(547) * density(15) 
  pd(14,14) = pd(14,14) - rrt(548) * density(44) 
  pd(14,44) = pd(14,44) - rrt(548) * density(14) 
  pd(18,14) = pd(18,14) + rrt(548) * density(44) 
  pd(18,44) = pd(18,44) + rrt(548) * density(14) 
  pd(43,14) = pd(43,14) + rrt(548) * density(44) 
  pd(43,44) = pd(43,44) + rrt(548) * density(14) 
  pd(44,14) = pd(44,14) - rrt(548) * density(44) 
  pd(44,44) = pd(44,44) - rrt(548) * density(14) 
  pd(07,07) = pd(07,07) - rrt(549) * density(13) 
  pd(07,13) = pd(07,13) - rrt(549) * density(07) 
  pd(13,07) = pd(13,07) - rrt(549) * density(13) 
  pd(13,13) = pd(13,13) - rrt(549) * density(07) 
  pd(30,07) = pd(30,07) + rrt(549) * density(13) 
  pd(30,13) = pd(30,13) + rrt(549) * density(07) 
  pd(07,07) = pd(07,07) - rrt(550) * density(16) 
  pd(07,16) = pd(07,16) - rrt(550) * density(07) 
  pd(16,07) = pd(16,07) - rrt(550) * density(16) 
  pd(16,16) = pd(16,16) - rrt(550) * density(07) 
  pd(30,07) = pd(30,07) + rrt(550) * density(16) 
  pd(30,16) = pd(30,16) + rrt(550) * density(07) 
  pd(07,07) = pd(07,07) - rrt(551) * density(15) 
  pd(07,15) = pd(07,15) - rrt(551) * density(07) 
  pd(15,07) = pd(15,07) - rrt(551) * density(15) 
  pd(15,15) = pd(15,15) - rrt(551) * density(07) 
  pd(30,07) = pd(30,07) + rrt(551) * density(15) 
  pd(30,15) = pd(30,15) + rrt(551) * density(07) 
  pd(07,07) = pd(07,07) - rrt(552) * density(14) 
  pd(07,14) = pd(07,14) - rrt(552) * density(07) 
  pd(14,07) = pd(14,07) - rrt(552) * density(14) 
  pd(14,14) = pd(14,14) - rrt(552) * density(07) 
  pd(30,07) = pd(30,07) + rrt(552) * density(14) 
  pd(30,14) = pd(30,14) + rrt(552) * density(07) 
  pd(34,36) = pd(34,36) + rrt(553) * density(43) 
  pd(34,43) = pd(34,43) + rrt(553) * density(36) 
  pd(36,36) = pd(36,36) - rrt(553) * density(43) 
  pd(36,43) = pd(36,43) - rrt(553) * density(36) 
  pd(43,36) = pd(43,36) - rrt(553) * density(43) 
  pd(43,43) = pd(43,43) - rrt(553) * density(36) 
  pd(44,36) = pd(44,36) + rrt(553) * density(43) 
  pd(44,43) = pd(44,43) + rrt(553) * density(36) 
  pd(34,37) = pd(34,37) + rrt(554) * density(43) 
  pd(34,43) = pd(34,43) + rrt(554) * density(37) 
  pd(37,37) = pd(37,37) - rrt(554) * density(43) 
  pd(37,43) = pd(37,43) - rrt(554) * density(37) 
  pd(43,37) = pd(43,37) - rrt(554) * density(43) 
  pd(43,43) = pd(43,43) - rrt(554) * density(37) 
  pd(44,37) = pd(44,37) + rrt(554) * density(43) 
  pd(44,43) = pd(44,43) + rrt(554) * density(37) 
  pd(34,38) = pd(34,38) + rrt(555) * density(43) 
  pd(34,43) = pd(34,43) + rrt(555) * density(38) 
  pd(38,38) = pd(38,38) - rrt(555) * density(43) 
  pd(38,43) = pd(38,43) - rrt(555) * density(38) 
  pd(43,38) = pd(43,38) - rrt(555) * density(43) 
  pd(43,43) = pd(43,43) - rrt(555) * density(38) 
  pd(44,38) = pd(44,38) + rrt(555) * density(43) 
  pd(44,43) = pd(44,43) + rrt(555) * density(38) 
  pd(07,36) = pd(07,36) + rrt(556) 
  pd(24,36) = pd(24,36) + rrt(556) 
  pd(36,36) = pd(36,36) - rrt(556) 
  pd(07,37) = pd(07,37) + rrt(557) 
  pd(24,37) = pd(24,37) + rrt(557) 
  pd(37,37) = pd(37,37) - rrt(557) 
  pd(07,38) = pd(07,38) + rrt(558) 
  pd(24,38) = pd(24,38) + rrt(558) 
  pd(38,38) = pd(38,38) - rrt(558) 
  pd(05,05) = pd(05,05) - rrt(559) * density(36) 
  pd(05,36) = pd(05,36) - rrt(559) * density(05) 
  pd(36,05) = pd(36,05) - rrt(559) * density(36) 
  pd(36,36) = pd(36,36) - rrt(559) * density(05) 
  pd(41,05) = pd(41,05) + rrt(559) * density(36) 
  pd(41,36) = pd(41,36) + rrt(559) * density(05) 
  pd(05,05) = pd(05,05) - rrt(560) * density(37) 
  pd(05,37) = pd(05,37) - rrt(560) * density(05) 
  pd(37,05) = pd(37,05) - rrt(560) * density(37) 
  pd(37,37) = pd(37,37) - rrt(560) * density(05) 
  pd(41,05) = pd(41,05) + rrt(560) * density(37) 
  pd(41,37) = pd(41,37) + rrt(560) * density(05) 
  pd(05,05) = pd(05,05) - rrt(561) * density(38) 
  pd(05,38) = pd(05,38) - rrt(561) * density(05) 
  pd(38,05) = pd(38,05) - rrt(561) * density(38) 
  pd(38,38) = pd(38,38) - rrt(561) * density(05) 
  pd(41,05) = pd(41,05) + rrt(561) * density(38) 
  pd(41,38) = pd(41,38) + rrt(561) * density(05) 
  pd(31,34) = pd(31,34) + rrt(562) * density(34) * 2.0d0
  pd(34,34) = pd(34,34) - rrt(562) * density(34) * 4.0d0
  pd(36,34) = pd(36,34) + rrt(562) * density(34) * 2.0d0
  pd(34,34) = pd(34,34) - rrt(563) * density(44) 
  pd(34,44) = pd(34,44) - rrt(563) * density(34) 
  pd(36,34) = pd(36,34) + rrt(563) * density(44) 
  pd(36,44) = pd(36,44) + rrt(563) * density(34) 
  pd(43,34) = pd(43,34) + rrt(563) * density(44) 
  pd(43,44) = pd(43,44) + rrt(563) * density(34) 
  pd(44,34) = pd(44,34) - rrt(563) * density(44) 
  pd(44,44) = pd(44,44) - rrt(563) * density(34) 
  pd(31,34) = pd(31,34) + rrt(564) * density(43) 
  pd(31,43) = pd(31,43) + rrt(564) * density(34) 
  pd(34,34) = pd(34,34) - rrt(564) * density(43) 
  pd(34,43) = pd(34,43) - rrt(564) * density(34) 
  pd(43,34) = pd(43,34) - rrt(564) * density(43) 
  pd(43,43) = pd(43,43) - rrt(564) * density(34) 
  pd(44,34) = pd(44,34) + rrt(564) * density(43) 
  pd(44,43) = pd(44,43) + rrt(564) * density(34) 
  pd(34,34) = pd(34,34) - rrt(565) * density(43) 
  pd(34,43) = pd(34,43) - rrt(565) * density(34) 
  pd(36,34) = pd(36,34) + rrt(565) * density(43) 
  pd(36,43) = pd(36,43) + rrt(565) * density(34) 
  pd(43,34) = pd(43,34) - rrt(565) * density(43) 
  pd(43,43) = pd(43,43) - rrt(565) * density(34) 
  pd(07,34) = pd(07,34) + rrt(566) * density(43) 
  pd(07,43) = pd(07,43) + rrt(566) * density(34) 
  pd(24,34) = pd(24,34) + rrt(566) * density(43) 
  pd(24,43) = pd(24,43) + rrt(566) * density(34) 
  pd(34,34) = pd(34,34) - rrt(566) * density(43) 
  pd(34,43) = pd(34,43) - rrt(566) * density(34) 
  pd(43,34) = pd(43,34) - rrt(566) * density(43) 
  pd(43,43) = pd(43,43) - rrt(566) * density(34) 
  pd(31,34) = pd(31,34) + rrt(567) 
  pd(34,34) = pd(34,34) - rrt(567) 
  pd(43,34) = pd(43,34) + rrt(567) 
  pd(07,34) = pd(07,34) + rrt(568) 
  pd(20,34) = pd(20,34) + rrt(568) 
  pd(34,34) = pd(34,34) - rrt(568) 
  pd(31,31) = pd(31,31) - rrt(569) * density(43) 
  pd(31,43) = pd(31,43) - rrt(569) * density(31) 
  pd(34,31) = pd(34,31) + rrt(569) * density(43) 
  pd(34,43) = pd(34,43) + rrt(569) * density(31) 
  pd(43,31) = pd(43,31) - rrt(569) * density(43) 
  pd(43,43) = pd(43,43) - rrt(569) * density(31) 
  pd(32,32) = pd(32,32) - rrt(570) * density(43) 
  pd(32,43) = pd(32,43) - rrt(570) * density(32) 
  pd(34,32) = pd(34,32) + rrt(570) * density(43) 
  pd(34,43) = pd(34,43) + rrt(570) * density(32) 
  pd(43,32) = pd(43,32) - rrt(570) * density(43) 
  pd(43,43) = pd(43,43) - rrt(570) * density(32) 
  pd(07,31) = pd(07,31) + rrt(571) 
  pd(18,31) = pd(18,31) + rrt(571) 
  pd(31,31) = pd(31,31) - rrt(571) 
  pd(07,07) = pd(07,07) - rrt(572) * density(41) 
  pd(07,41) = pd(07,41) - rrt(572) * density(07) 
  pd(09,07) = pd(09,07) + rrt(572) * density(41) 
  pd(09,41) = pd(09,41) + rrt(572) * density(07) 
  pd(40,07) = pd(40,07) + rrt(572) * density(41) 
  pd(40,41) = pd(40,41) + rrt(572) * density(07) 
  pd(41,07) = pd(41,07) - rrt(572) * density(41) 
  pd(41,41) = pd(41,41) - rrt(572) * density(07) 
  pd(07,41) = pd(07,41) + rrt(573) 
  pd(34,41) = pd(34,41) + rrt(573) 
  pd(41,41) = pd(41,41) - rrt(573) 
  pd(24,41) = pd(24,41) + rrt(574) * 2.0d0
  pd(41,41) = pd(41,41) - rrt(574) 
  pd(05,05) = pd(05,05) - rrt(575) * density(41) 
  pd(05,41) = pd(05,41) - rrt(575) * density(05) 
  pd(41,05) = pd(41,05) - rrt(575) * density(41) 
  pd(41,41) = pd(41,41) - rrt(575) * density(05) 
  pd(42,05) = pd(42,05) + rrt(575) * density(41) 
  pd(42,41) = pd(42,41) + rrt(575) * density(05) 
  pd(43,44) = pd(43,44) + rrt(576) * 2.0d0
  pd(44,44) = pd(44,44) - rrt(576) 
  pd(43,43) = pd(43,43) - rrt(577) * density(43) * 4.0d0
  pd(44,43) = pd(44,43) + rrt(577) * density(43) * 2.0d0
  if( ldensity_constant ) then
    do i = 1, species_max
      if( density_constant(i) ) pd(i,:) = 0.0d0
    enddo
  endif
  if( lgas_heating ) then
    pd(48,1) = eV_to_K * ZDPlasKin_cfg(11)
    pd(48,:) = pd(48,:) * ZDPlasKin_cfg(13)
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
  DOUBLE PRECISION, PARAMETER :: F0 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F1 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F2 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F3 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F4 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F5 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F6 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F7 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F8 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F9 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F10 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F11 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F12 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F13 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F14 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F15 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F16 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F17 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F18 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F19 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F20 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F21 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F22 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F23 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F24 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F25 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F26 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F27 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F28 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F29 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F30 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F31 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F32 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F33 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F34 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F35 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F36 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F37 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F38 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F39 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F40 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F41 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F42 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F43 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F44 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F45 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F46 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F47 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F48 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F49 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F50 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F51 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F52 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F53 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F54 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F55 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F56 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F57 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F58 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F59 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F60 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F61 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F62 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F63 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F64 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F65 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F66 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F67 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F68 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F69 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F70 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F71 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F72 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F73 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F74 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F75 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F76 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F77 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F78 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F79 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F80 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F81 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F82 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F83 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F84 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F85 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F86 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F87 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F88 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F89 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F90 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F91 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F92 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F93 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F94 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F95 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F96 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F97 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F98 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F99 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F100 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F101 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F102 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F103 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F104 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F105 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F106 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F107 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F108 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F109 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F110 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F111 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F112 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F113 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F114 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F115 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F116 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F117 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F118 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F119 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F120 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F121 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F122 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F123 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F124 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F125 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F126 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F127 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F128 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F129 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F130 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F131 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F132 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F133 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F134 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F135 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F136 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F137 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F138 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F139 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F140 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F141 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F142 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F143 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F144 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F145 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F146 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F147 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F148 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F149 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F150 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F151 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F152 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F153 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F154 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F155 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F156 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F157 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F158 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F159 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F160 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F161 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F162 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F163 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F164 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F165 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F166 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F167 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F168 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F169 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F170 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F171 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F172 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F173 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F174 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F175 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F176 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F177 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F178 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F179 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F180 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F181 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F182 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F183 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F184 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F185 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F186 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F187 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F188 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F189 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F190 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F191 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F192 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F193 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F194 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F195 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F196 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F197 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F198 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F199 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F200 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F201 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F202 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F203 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F204 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F205 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F206 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F207 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F208 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F209 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F210 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F211 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F212 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F213 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F214 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F215 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F216 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F217 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F218 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F219 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F220 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F221 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F222 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F223 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F224 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F225 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F226 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F227 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F228 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F229 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F230 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F231 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F232 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F233 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F234 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F235 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F236 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F237 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F238 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F239 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F240 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F241 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F242 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F243 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F244 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F245 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F246 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F247 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F248 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F249 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F250 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F251 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F252 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F253 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F254 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F255 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F256 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F257 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F258 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F259 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F260 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F261 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F262 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F263 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F264 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F265 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F266 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F267 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F268 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F269 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F270 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F271 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F272 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F273 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F274 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F275 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F276 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F277 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F278 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F279 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F280 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F281 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F282 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F283 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F284 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F285 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F286 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F287 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F288 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F289 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F290 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F291 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F292 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F293 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F294 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F295 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F296 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F297 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F298 = 1.000D0
  DOUBLE PRECISION, PARAMETER :: F299 = 1.000D0
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
  rrt(014) = F12*bolsig_rates(bolsig_pointer(14))
  rrt(015) = F12*bolsig_rates(bolsig_pointer(15))
  rrt(016) = F13*bolsig_rates(bolsig_pointer(16))
  rrt(017) = F13*bolsig_rates(bolsig_pointer(17))
  rrt(018) = F13*bolsig_rates(bolsig_pointer(18))
  rrt(019) = F14*bolsig_rates(bolsig_pointer(19))
  rrt(020) = F14*bolsig_rates(bolsig_pointer(20))
  rrt(021) = F14*bolsig_rates(bolsig_pointer(21))
  rrt(022) = F15*bolsig_rates(bolsig_pointer(22))
  rrt(023) = F15*bolsig_rates(bolsig_pointer(23))
  rrt(024) = F15*bolsig_rates(bolsig_pointer(24))
  rrt(025) = F16*bolsig_rates(bolsig_pointer(25))
  rrt(026) = F17*bolsig_rates(bolsig_pointer(26))
  rrt(027) = F18*bolsig_rates(bolsig_pointer(27))
  rrt(028) = F19*bolsig_rates(bolsig_pointer(28))
  rrt(029) = F20*bolsig_rates(bolsig_pointer(29))
  rrt(030) = F21*bolsig_rates(bolsig_pointer(30))
  rrt(031) = F22*bolsig_rates(bolsig_pointer(31))
  rrt(032) = F23*bolsig_rates(bolsig_pointer(32))
  rrt(033) = F23*bolsig_rates(bolsig_pointer(33))
  rrt(034) = F23*bolsig_rates(bolsig_pointer(34))
  rrt(035) = F24*bolsig_rates(bolsig_pointer(35))
  rrt(036) = F24*bolsig_rates(bolsig_pointer(36))
  rrt(037) = F24*bolsig_rates(bolsig_pointer(37))
  rrt(038) = F25*bolsig_rates(bolsig_pointer(38))
  rrt(039) = F25*bolsig_rates(bolsig_pointer(39))
  rrt(040) = F25*bolsig_rates(bolsig_pointer(40))
  rrt(041) = F26*bolsig_rates(bolsig_pointer(41))
  rrt(042) = F26*bolsig_rates(bolsig_pointer(42))
  rrt(043) = F26*bolsig_rates(bolsig_pointer(43))
  rrt(044) = F27*bolsig_rates(bolsig_pointer(44))
  rrt(045) = F28*bolsig_rates(bolsig_pointer(45))
  rrt(046) = F29*bolsig_rates(bolsig_pointer(46))
  rrt(047) = F30*bolsig_rates(bolsig_pointer(47))
  rrt(048) = F31*bolsig_rates(bolsig_pointer(48))
  rrt(049) = F32*bolsig_rates(bolsig_pointer(49))
  rrt(050) = F33*bolsig_rates(bolsig_pointer(50))
  rrt(051) = F33*bolsig_rates(bolsig_pointer(51))
  rrt(052) = F33*bolsig_rates(bolsig_pointer(52))
  rrt(053) = F34*bolsig_rates(bolsig_pointer(53))
  rrt(054) = F34*bolsig_rates(bolsig_pointer(54))
  rrt(055) = F34*bolsig_rates(bolsig_pointer(55))
  rrt(056) = F35*bolsig_rates(bolsig_pointer(56))
  rrt(057) = F35*bolsig_rates(bolsig_pointer(57))
  rrt(058) = F35*bolsig_rates(bolsig_pointer(58))
  rrt(059) = F36*bolsig_rates(bolsig_pointer(59))
  rrt(060) = F36*bolsig_rates(bolsig_pointer(60))
  rrt(061) = F36*bolsig_rates(bolsig_pointer(61))
  rrt(062) = F37*bolsig_rates(bolsig_pointer(62))
  rrt(063) = F37*bolsig_rates(bolsig_pointer(63))
  rrt(064) = F37*bolsig_rates(bolsig_pointer(64))
  rrt(065) = F38*bolsig_rates(bolsig_pointer(65))
  rrt(066) = F38*bolsig_rates(bolsig_pointer(66))
  rrt(067) = F38*bolsig_rates(bolsig_pointer(67))
  rrt(068) = F39*bolsig_rates(bolsig_pointer(68))
  rrt(069) = F40*bolsig_rates(bolsig_pointer(69))
  rrt(070) = F41*bolsig_rates(bolsig_pointer(70))
  rrt(071) = F42*bolsig_rates(bolsig_pointer(71))
  rrt(072) = F43*bolsig_rates(bolsig_pointer(72))
  rrt(073) = F44*bolsig_rates(bolsig_pointer(73))
  rrt(074) = F45*bolsig_rates(bolsig_pointer(74))
  rrt(075) = F45*bolsig_rates(bolsig_pointer(75))
  rrt(076) = F45*bolsig_rates(bolsig_pointer(76))
  rrt(077) = F46*bolsig_rates(bolsig_pointer(77))
  rrt(078) = F46*bolsig_rates(bolsig_pointer(78))
  rrt(079) = F46*bolsig_rates(bolsig_pointer(79))
  rrt(080) = F47*bolsig_rates(bolsig_pointer(80))
  rrt(081) = F47*bolsig_rates(bolsig_pointer(81))
  rrt(082) = F47*bolsig_rates(bolsig_pointer(82))
  rrt(083) = F48*bolsig_rates(bolsig_pointer(83))
  rrt(084) = F48*bolsig_rates(bolsig_pointer(84))
  rrt(085) = F48*bolsig_rates(bolsig_pointer(85))
  rrt(086) = F49*bolsig_rates(bolsig_pointer(86))
  rrt(087) = F49*bolsig_rates(bolsig_pointer(87))
  rrt(088) = F49*bolsig_rates(bolsig_pointer(88))
  rrt(089) = F50*bolsig_rates(bolsig_pointer(89))
  rrt(090) = F51*bolsig_rates(bolsig_pointer(90))
  rrt(091) = F52*bolsig_rates(bolsig_pointer(91))
  rrt(092) = F52*bolsig_rates(bolsig_pointer(92))
  rrt(093) = F52*bolsig_rates(bolsig_pointer(93))
  rrt(094) = F52*bolsig_rates(bolsig_pointer(94))
  rrt(095) = F53*bolsig_rates(bolsig_pointer(95))
  rrt(096) = F53*bolsig_rates(bolsig_pointer(96))
  rrt(097) = F53*bolsig_rates(bolsig_pointer(97))
  rrt(098) = F54*bolsig_rates(bolsig_pointer(98))
  rrt(099) = F54*bolsig_rates(bolsig_pointer(99))
  rrt(100) = F54*bolsig_rates(bolsig_pointer(100))
  rrt(101) = F55*bolsig_rates(bolsig_pointer(101))
  rrt(102) = F55*bolsig_rates(bolsig_pointer(102))
  rrt(103) = F55*bolsig_rates(bolsig_pointer(103))
  rrt(104) = F56*bolsig_rates(bolsig_pointer(104))
  rrt(105) = F56*bolsig_rates(bolsig_pointer(105))
  rrt(106) = F56*bolsig_rates(bolsig_pointer(106))
  rrt(107) = F57*bolsig_rates(bolsig_pointer(107))
  rrt(108) = F57*bolsig_rates(bolsig_pointer(108))
  rrt(109) = F57*bolsig_rates(bolsig_pointer(109))
  rrt(110) = F58*bolsig_rates(bolsig_pointer(110))
  rrt(111) = F58*bolsig_rates(bolsig_pointer(111))
  rrt(112) = F58*bolsig_rates(bolsig_pointer(112))
  rrt(113) = F59*bolsig_rates(bolsig_pointer(113))
  rrt(114) = F59*bolsig_rates(bolsig_pointer(114))
  rrt(115) = F59*bolsig_rates(bolsig_pointer(115))
  rrt(116) = F60*bolsig_rates(bolsig_pointer(116))
  rrt(117) = F61*bolsig_rates(bolsig_pointer(117))
  rrt(118) = F62*bolsig_rates(bolsig_pointer(118))
  rrt(119) = F63*bolsig_rates(bolsig_pointer(119))
  rrt(120) = F64*bolsig_rates(bolsig_pointer(120))
  rrt(121) = F65*bolsig_rates(bolsig_pointer(121))
  rrt(122) = F66*bolsig_rates(bolsig_pointer(122))
  rrt(123) = F67*bolsig_rates(bolsig_pointer(123))
  rrt(124) = F67*bolsig_rates(bolsig_pointer(124))
  rrt(125) = F67*bolsig_rates(bolsig_pointer(125))
  rrt(126) = F68*bolsig_rates(bolsig_pointer(126))
  rrt(127) = F68*bolsig_rates(bolsig_pointer(127))
  rrt(128) = F68*bolsig_rates(bolsig_pointer(128))
  rrt(129) = F69*bolsig_rates(bolsig_pointer(129))
  rrt(130) = F69*bolsig_rates(bolsig_pointer(130))
  rrt(131) = F69*bolsig_rates(bolsig_pointer(131))
  rrt(132) = F70*bolsig_rates(bolsig_pointer(132))
  rrt(133) = F70*bolsig_rates(bolsig_pointer(133))
  rrt(134) = F70*bolsig_rates(bolsig_pointer(134))
  rrt(135) = F71*bolsig_rates(bolsig_pointer(135))
  rrt(136) = F71*bolsig_rates(bolsig_pointer(136))
  rrt(137) = F71*bolsig_rates(bolsig_pointer(137))
  rrt(138) = F72*bolsig_rates(bolsig_pointer(138))
  rrt(139) = F73*bolsig_rates(bolsig_pointer(139))
  rrt(140) = F74*bolsig_rates(bolsig_pointer(140))
  rrt(141) = F75*bolsig_rates(bolsig_pointer(141))
  rrt(142) = F76*bolsig_rates(bolsig_pointer(142))
  rrt(143) = F77*bolsig_rates(bolsig_pointer(143))
  rrt(144) = F77*bolsig_rates(bolsig_pointer(144))
  rrt(145) = F77*bolsig_rates(bolsig_pointer(145))
  rrt(146) = F77*bolsig_rates(bolsig_pointer(146))
  rrt(147) = F78*bolsig_rates(bolsig_pointer(147))
  rrt(148) = F78*bolsig_rates(bolsig_pointer(148))
  rrt(149) = F78*bolsig_rates(bolsig_pointer(149))
  rrt(150) = F78*bolsig_rates(bolsig_pointer(150))
  rrt(151) = F79*bolsig_rates(bolsig_pointer(151))
  rrt(152) = F79*bolsig_rates(bolsig_pointer(152))
  rrt(153) = F79*bolsig_rates(bolsig_pointer(153))
  rrt(154) = F80*bolsig_rates(bolsig_pointer(154))
  rrt(155) = F80*bolsig_rates(bolsig_pointer(155))
  rrt(156) = F80*bolsig_rates(bolsig_pointer(156))
  rrt(157) = F81*bolsig_rates(bolsig_pointer(157))
  rrt(158) = F82*bolsig_rates(bolsig_pointer(158))
  rrt(159) = F82*bolsig_rates(bolsig_pointer(159))
  rrt(160) = F82*bolsig_rates(bolsig_pointer(160))
  rrt(161) = F83*bolsig_rates(bolsig_pointer(161))
  rrt(162) = F83*bolsig_rates(bolsig_pointer(162))
  rrt(163) = F83*bolsig_rates(bolsig_pointer(163))
  rrt(164) = F84*bolsig_rates(bolsig_pointer(164))
  rrt(165) = F84*bolsig_rates(bolsig_pointer(165))
  rrt(166) = F84*bolsig_rates(bolsig_pointer(166))
  rrt(167) = F85*bolsig_rates(bolsig_pointer(167))
  rrt(168) = F86*bolsig_rates(bolsig_pointer(168))
  rrt(169) = F87*bolsig_rates(bolsig_pointer(169))
  rrt(170) = F87*bolsig_rates(bolsig_pointer(170))
  rrt(171) = F88*bolsig_rates(bolsig_pointer(171))
  rrt(172) = F88*bolsig_rates(bolsig_pointer(172))
  rrt(173) = F89*bolsig_rates(bolsig_pointer(173))
  rrt(174) = F89*bolsig_rates(bolsig_pointer(174))
  rrt(175) = F90*1.00D-16
  rrt(176) = F91*1.00D-16
  rrt(177) = F92*1.00D-16
  rrt(178) = F93*1.00D-16
  rrt(179) = F94*1.00D-16
  rrt(180) = F95*1.00D-16
  rrt(181) = F96*1.00D-16
  rrt(182) = F97*1.00D-16
  rrt(183) = F98*1.00D-16
  rrt(184) = F99*1.00D-16
  rrt(185) = F100*1.00D-16
  rrt(186) = F101*1.00D-16
  rrt(187) = F102*1.18D-08*(300./TGAS)**0.5
  rrt(188) = F103*2.42D-08*(300./TGAS)**0.5
  rrt(189) = F104*1.41D-08*(300./TGAS)**0.5
  rrt(190) = F105*2.25D-08*(300./TGAS)**0.5
  rrt(191) = F106*7.88D-09*(300./TGAS)**0.5
  rrt(192) = F107*9.00D-09*(300./TGAS)**0.5
  rrt(193) = F108*1.69D-08*(300./TGAS)**0.5
  rrt(194) = F109*1.00D-08*(300./TGAS)**0.5
  rrt(195) = F110*4.82D-09*(300./TGAS)**0.5
  rrt(196) = F111*2.53D-08*(300./TGAS)**0.5
  rrt(197) = F112*3.23D-08*(300./TGAS)**0.42
  rrt(198) = F113*2.19D-08*(300./TGAS)**0.71
  rrt(199) = F114*3.36D-08*(300./TGAS)**0.71
  rrt(200) = F115*7.70D-09*(300./TGAS)**0.71
  rrt(201) = F116*1,92D-08*(300./TGAS)**0.71
  rrt(202) = F117*1.60D-08*(300./TGAS)**0.71
  rrt(203) = F118*8.98D-09*(300./TGAS)**0.71
  rrt(204) = F119*9.62D-09*(300./TGAS)**0.71
  rrt(205) = F120*8.29D-09*(300./TGAS)**0.71
  rrt(206) = F121*3.43D-08*(300./TGAS)**0.71
  rrt(207) = F122*1.34D-08*(300./TGAS)**0.71
  rrt(208) = F123*4.87D-09*(300./TGAS)**0.71
  rrt(209) = F124*3.17D21/(6.022D23*TE**4.5)
  rrt(210) = F125*3.17D21/(6.022D23*TE**4.5)
  rrt(211) = F126*2.11D-09
  rrt(212) = F127*1.91D-9
  rrt(213) = rrt(212)
  rrt(214) = rrt(212)
  rrt(215) = F128*4.23D-10
  rrt(216) = rrt(215)
  rrt(217) = rrt(215)
  rrt(218) = F129*1.38D-9
  rrt(219) = rrt(218)
  rrt(220) = rrt(218)
  rrt(221) = F130*1.23D-9
  rrt(222) = rrt(221)
  rrt(223) = rrt(221)
  rrt(224) = rrt(221)
  rrt(225) = F131*1.13D-9
  rrt(226) = rrt(225)
  rrt(227) = rrt(225)
  rrt(228) = rrt(225)
  rrt(229) = F132*1.00D-11
  rrt(230) = F133*1.36D-10
  rrt(231) = rrt(230)
  rrt(232) = rrt(230)
  rrt(233) = F134*1.20D-9
  rrt(234) = rrt(233)
  rrt(235) = rrt(233)
  rrt(236) = F135*9.90D-10
  rrt(237) = F136*7.10D-10
  rrt(238) = F137*1.48D-9
  rrt(239) = rrt(238)
  rrt(240) = rrt(238)
  rrt(241) = F138*3.50D-10
  rrt(242) = rrt(241)
  rrt(243) = rrt(241)
  rrt(244) = F139*3.00D-10
  rrt(245) = F140*1.38D-10
  rrt(246) = rrt(245)
  rrt(247) = rrt(245)
  rrt(248) = F141*3.60D-10
  rrt(249) = rrt(248)
  rrt(250) = rrt(248)
  rrt(251) = F142*8.40D-10
  rrt(252) = rrt(251)
  rrt(253) = rrt(251)
  rrt(254) = F143*2.31D-10
  rrt(255) = rrt(254)
  rrt(256) = rrt(254)
  rrt(257) = F144*3.97D-10
  rrt(258) = rrt(257)
  rrt(259) = rrt(257)
  rrt(260) = F145*1.60D-9
  rrt(261) = F146*6.50D-11
  rrt(262) = rrt(261)
  rrt(263) = rrt(261)
  rrt(264) = F147*1.09D-9
  rrt(265) = rrt(264)
  rrt(266) = rrt(264)
  rrt(267) = F148*1.43D-10
  rrt(268) = rrt(267)
  rrt(269) = rrt(267)
  rrt(270) = F149*1.20D-9
  rrt(271) = F150*1.15D-9
  rrt(272) = rrt(271)
  rrt(273) = rrt(271)
  rrt(274) = F151*2.47D-10
  rrt(275) = rrt(274)
  rrt(276) = rrt(274)
  rrt(277) = rrt(274)
  rrt(278) = F152*1.00D-10
  rrt(279) = F153*1.00D-11
  rrt(280) = F154*5.00D-10
  rrt(281) = F155*5.00D-10
  rrt(282) = F156*3.00D-10
  rrt(283) = F157*2.91D-10
  rrt(284) = rrt(283)
  rrt(285) = rrt(283)
  rrt(286) = F158*8.90D-10
  rrt(287) = rrt(286)
  rrt(288) = rrt(286)
  rrt(289) = F159*6.80D-11
  rrt(290) = F160*4.10D-9
  rrt(291) = rrt(290)
  rrt(292) = rrt(290)
  rrt(293) = F161*1.31D-10
  rrt(294) = rrt(293)
  rrt(295) = rrt(293)
  rrt(296) = F162*2.48D-10
  rrt(297) = rrt(296)
  rrt(298) = rrt(296)
  rrt(299) = F163*4.14D-10
  rrt(300) = rrt(299)
  rrt(301) = rrt(299)
  rrt(302) = F164*3.30D-10
  rrt(303) = F165*1.00D-11
  rrt(304) = F166*2.10D-9
  rrt(305) = F167*1.70D-9
  rrt(306) = F168*1.20D-9
  rrt(307) = F169*2.40D-9
  rrt(308) = rrt(307)
  rrt(309) = rrt(307)
  rrt(310) = F170*1.40D-9
  rrt(311) = F171*1.15D-9
  rrt(312) = rrt(311)
  rrt(313) = rrt(311)
  rrt(314) = F172*1.15D-9
  rrt(315) = rrt(314)
  rrt(316) = rrt(314)
  rrt(317) = F173*2.00D-9
  rrt(318) = F174*3.50D-9
  rrt(319) = rrt(318)
  rrt(320) = rrt(318)
  rrt(321) = rrt(318)
  rrt(322) = F175*1.40D-9
  rrt(323) = rrt(322)
  rrt(324) = rrt(322)
  rrt(325) = F176*2.30D-9
  rrt(326) = rrt(325)
  rrt(327) = rrt(325)
  rrt(328) = F177*1.00D-9
  rrt(329) = F178*1.00D-9
  rrt(330) = F179*7.10D-10
  rrt(331) = F180*7.10D-10
  rrt(332) = F181*2.94D-10
  rrt(333) = rrt(332)
  rrt(334) = rrt(332)
  rrt(335) = F182*1.37D-9
  rrt(336) = rrt(335)
  rrt(337) = rrt(335)
  rrt(338) = F183*2.35D-9
  rrt(339) = rrt(338)
  rrt(340) = rrt(338)
  rrt(341) = F184*6.86D-10
  rrt(342) = rrt(341)
  rrt(343) = rrt(341)
  rrt(344) = F185*1.96D-10
  rrt(345) = rrt(344)
  rrt(346) = rrt(344)
  rrt(347) = F186*2.21D-9
  rrt(348) = rrt(347)
  rrt(349) = rrt(347)
  rrt(350) = F187*1.81D-9
  rrt(351) = rrt(350)
  rrt(352) = rrt(350)
  rrt(353) = F188*8.82D-10
  rrt(354) = rrt(353)
  rrt(355) = rrt(353)
  rrt(356) = F189*4.80D-10
  rrt(357) = rrt(356)
  rrt(358) = rrt(356)
  rrt(359) = rrt(356)
  rrt(360) = F190*4.82D-9
  rrt(361) = rrt(360)
  rrt(362) = rrt(360)
  rrt(363) = rrt(360)
  rrt(364) = F191*2.10D-9
  rrt(365) = F192*6.39D-10
  rrt(366) = F193*1.50D-9
  rrt(367) = rrt(366)
  rrt(368) = rrt(366)
  rrt(369) = F194*2.30D-9
  rrt(370) = rrt(369)
  rrt(371) = rrt(369)
  rrt(372) = F195*3.40D-9
  rrt(373) = F196*1.40D-9
  rrt(374) = F197*1.40D-9
  rrt(375) = F198*1.90D-9
  rrt(376) = F199*1.30D-9
  rrt(377) = rrt(376)
  rrt(378) = rrt(376)
  rrt(379) = F200*1.40D-9
  rrt(380) = rrt(379)
  rrt(381) = rrt(379)
  rrt(382) = F201*2.80D-9
  rrt(383) = rrt(382)
  rrt(384) = rrt(382)
  rrt(385) = F202*1.65D-9
  rrt(386) = F203*3.06D-9
  rrt(387) = F204*1.00D-9
  rrt(388) = rrt(387)
  rrt(389) = rrt(387)
  rrt(390) = F205*3.00D-9
  rrt(391) = rrt(390)
  rrt(392) = rrt(390)
  rrt(393) = F206*1.00D-9
  rrt(394) = F207*2.00D-9
  rrt(395) = F208*2.00D-9
  rrt(396) = F209*5.40D-10
  rrt(397) = rrt(396)
  rrt(398) = rrt(396)
  rrt(399) = rrt(396)
  rrt(400) = F210*4.08D-18*TGAS**2*EXP(-4163./TGAS)
  rrt(401) = rrt(400)
  rrt(402) = rrt(400)
  rrt(403) = F211*9.97D-11
  rrt(404) = rrt(403)
  rrt(405) = rrt(403)
  rrt(406) = F212*2.51D-15*(TGAS/298.)**4.14*EXP(-52.55/(R*TGAS))
  rrt(407) = rrt(406)
  rrt(408) = rrt(406)
  rrt(409) = F213*4.26D-15*(TGAS/298.)**4.02*EXP(-22.86/(R*TGAS))
  rrt(410) = rrt(409)
  rrt(411) = rrt(409)
  rrt(412) = F214*3.54D-16*(TGAS/298.)**4.02*EXP(-45.48/(R*TGAS))
  rrt(413) = rrt(412)
  rrt(414) = rrt(412)
  rrt(415) = F215*9.86D-13*(TGAS/298.)**3.00*EXP(-36.67/(R*TGAS))
  rrt(416) = rrt(415)
  rrt(417) = rrt(415)
  rrt(418) = F216*1.33D-10*EXP(-167.00/(R*TGAS))
  rrt(419) = rrt(418)
  rrt(420) = rrt(418)
  rrt(421) = F217*1.90D-12
  rrt(422) = rrt(421)
  rrt(423) = rrt(421)
  rrt(424) = F218*2.40D16*EXP(-52800./TGAS)
  rrt(425) = rrt(424)
  rrt(426) = rrt(424)
  rrt(427) = F219*1.69D-8*EXP(-379./(R*TGAS))
  rrt(428) = F220*6.97D-9*EXP(-345./(R*TGAS))
  rrt(429) = F221*3.00D-44*TGAS**9.10
  rrt(430) = F222*3.32D-10*EXP(-45.98/(R*TGAS))
  rrt(431) = F223*3.01D-11
  rrt(432) = F224*3.01D-11
  rrt(433) = F225*1.61D-15*(TGAS/298.)**3.65*EXP(-29.93/(R*TGAS))
  rrt(434) = rrt(433)
  rrt(435) = rrt(433)
  rrt(436) = F226*3.01D-11
  rrt(437) = F227*3.01D-11
  rrt(438) = F228*8.30D-19*TGAS**2.00*EXP(-3938.65/TGAS)
  rrt(439) = F229*1.00D-11*EXP(7.48/(R*TGAS))
  rrt(440) = F230*6.64D-9*EXP(-348./(R*TGAS))
  rrt(441) = F231*9.96D-10
  rrt(442) = F232*3.00D-11
  rrt(443) = rrt(442)
  rrt(444) = rrt(442)
  rrt(445) = F233*1.14D-29
  rrt(446) = rrt(445)
  rrt(447) = rrt(445)
  rrt(448) = F234*1.79D-10*EXP(-1565.17/TGAS)
  rrt(449) = F235*4.98D-11
  rrt(450) = F236*6.64D-11
  rrt(451) = F237*3.29D-12*TGAS**0.43*EXP(186.21/TGAS)
  rrt(452) = F238*8.30D-11
  rrt(453) = F239*1.46D-13*(TGAS/298.)**3.30*EXP(-43.90/(R*TGAS))
  rrt(454) = rrt(453)
  rrt(455) = rrt(453)
  rrt(456) = F240*1.19D-15*(TGAS/298.)**3.82*EXP(-37.83/(R*TGAS))
  rrt(457) = rrt(456)
  rrt(458) = rrt(456)
  rrt(459) = F241*1.23D-11*(TGAS/298.)**1.50*EXP(-31.01/(R*TGAS))
  rrt(460) = rrt(459)
  rrt(461) = rrt(459)
  rrt(462) = F242*8.97D-20*EXP(-48.64/(R*TGAS))
  rrt(463) = rrt(462)
  rrt(464) = rrt(462)
  rrt(465) = F243*1.80D21*TGAS**-1.24*EXP(-45700./TGAS)
  rrt(466) = rrt(465)
  rrt(467) = rrt(465)
  rrt(468) = F244*1.79D-10*EXP(132.36/TGAS)
  rrt(469) = rrt(468)
  rrt(470) = rrt(468)
  rrt(471) = F245*9.00D-33*TGAS**6.43
  rrt(472) = rrt(471)
  rrt(473) = rrt(471)
  rrt(474) = F246*2.41D-12
  rrt(475) = F247*5.83D-14*(TGAS/298.)**3.13*EXP(-75.33/(R*TGAS))
  rrt(476) = rrt(475)
  rrt(477) = rrt(475)
  rrt(478) = F248*1.61D-15*(TGAS/298.)**3.65*EXP(-38.25/(R*TGAS))
  rrt(479) = rrt(478)
  rrt(480) = rrt(478)
  rrt(481) = F249*1.91D-12
  rrt(482) = F250*2.41D-12
  rrt(483) = F251*1.69D-15*(TGAS/298.)**3.50*EXP(-35.34/(R*TGAS))
  rrt(484) = rrt(483)
  rrt(485) = F252*4.12D-15*(TGAS/298.)**3.60*EXP(-27.77/(R*TGAS))
  rrt(486) = F253*5.99D-11
  rrt(487) = F254*3.32D-12
  rrt(488) = F255*8.65D-7*TGAS**-0.99*EXP(-795.17/TGAS)
  rrt(489) = F256*4.08D12*(TGAS/298.)**1.04*EXP(-154./(R*TGAS))
  rrt(490) = F257*9.55D-12
  rrt(491) = F258*4.42D-11
  rrt(492) = F259*9.00D-10*EXP(-62.36/(R*TGAS))
  rrt(493) = rrt(492)
  rrt(494) = rrt(492)
  rrt(495) = F260*9.68D-12*(TGAS/298.)**1.28*EXP(-5.40/(R*TGAS))
  rrt(496) = rrt(495)
  rrt(497) = rrt(495)
  rrt(498) = F261*9.55D-12
  rrt(499) = F262*4.30D-7*EXP(-404./(R*TGAS))
  rrt(500) = rrt(499)
  rrt(501) = rrt(499)
  rrt(502) = F263*4.00D-11*EXP(-286./(R*TGAS))
  rrt(503) = rrt(502)
  rrt(504) = rrt(502)
  rrt(505) = rrt(502)
  rrt(506) = rrt(502)
  rrt(507) = rrt(502)
  rrt(508) = F264*1.00D-10*EXP(-316./(R*TGAS))
  rrt(509) = rrt(508)
  rrt(510) = rrt(508)
  rrt(511) = rrt(508)
  rrt(512) = F265*8.00D-10*EXP(-299./(R*TGAS))
  rrt(513) = rrt(512)
  rrt(514) = rrt(512)
  rrt(515) = F266*4.23D-18*TGAS**1.60*EXP(-2868.65/TGAS)
  rrt(516) = rrt(515)
  rrt(517) = rrt(515)
  rrt(518) = F267*8.00D12*TGAS**0.44*EXP(-44675.39/TGAS)
  rrt(519) = rrt(518)
  rrt(520) = rrt(518)
  rrt(521) = F268*4.75D-16*EXP(-180./(R*TGAS))
  rrt(522) = rrt(521)
  rrt(523) = rrt(521)
  rrt(524) = F269*5.30D-12*EXP(-2660./TGAS)
  rrt(525) = rrt(524)
  rrt(526) = rrt(524)
  rrt(527) = F270*3.50D-11
  rrt(528) = F271*1.46D-13*(TGAS/298.)**3.30*EXP(-43.90/(R*TGAS))
  rrt(529) = rrt(528)
  rrt(530) = rrt(528)
  rrt(531) = F272*2.01D-12
  rrt(532) = F273*2.01D-12
  rrt(533) = F274*1.61D-13*(TGAS/298.)**2.63*EXP(-35.75/(R*TGAS))
  rrt(534) = F275*1.60D-10
  rrt(535) = F276*1.01D-11*TGAS**0.27*EXP(-140.92/TGAS)
  rrt(536) = F277*2.00D14*EXP(-20000./TGAS)
  rrt(537) = F278*9.30D-12*EXP(-1207.85/TGAS)
  rrt(538) = rrt(537)
  rrt(539) = rrt(537)
  rrt(540) = rrt(537)
  rrt(541) = F279*5.00D-13*EXP(-163./(R*TGAS))
  rrt(542) = rrt(541)
  rrt(543) = rrt(541)
  rrt(544) = rrt(541)
  rrt(545) = F280*4.00D-12*EXP(-272./(R*TGAS))
  rrt(546) = rrt(545)
  rrt(547) = rrt(545)
  rrt(548) = rrt(545)
  rrt(549) = F281*1.58D-5*(TGAS/298.)**-8.58*EXP(-84.81/(R*TGAS))
  rrt(550) = rrt(549)
  rrt(551) = rrt(549)
  rrt(552) = rrt(549)
  rrt(553) = F282*2.19D-180*TGAS**2.54*EXP(-3400.10/TGAS)
  rrt(554) = rrt(553)
  rrt(555) = rrt(553)
  rrt(556) = F283*1.10D17*EXP(-42470./TGAS)
  rrt(557) = rrt(556)
  rrt(558) = rrt(556)
  rrt(559) = F284*4.42D-12
  rrt(560) = rrt(559)
  rrt(561) = rrt(559)
  rrt(562) = F285*2.81D-12
  rrt(563) = F286*3.19D-14*(TGAS/298.)**2.84*EXP(-38.25/(R*TGAS))
  rrt(564) = F287*3.01D-12
  rrt(565) = F288*6.00D-11
  rrt(566) = F289*6.74D-18*TGAS**2.19*EXP(-447.91/TGAS)
  rrt(567) = F290*1.25D17*EXP(-237./(R*TGAS))
  rrt(568) = F291*16.0*(TGAS/298.)**-(10.00)*EXP(-150./(R*TGAS))
  rrt(569) = F292*1.29D-11*(TGAS/298.)**0.51*EXP(-5.15/(R*TGAS))
  rrt(570) = rrt(569)
  rrt(571) = F293*1.28D13*(TGAS/298.)**(-15.70)*EXP(-502./(R*TGAS))
  rrt(572) = F294*8.32D-13*EXP(-56.87/(R*TGAS))
  rrt(573) = F295*8.87D-7*EXP(-180./(R*TGAS))
  rrt(574) = F296*7.84D-6*EXP(-207./(R*TGAS))
  rrt(575) = F297*9.61D-13
  rrt(576) = F298*1.88D-8*(TGAS/298.)**(-1.10)*EXP(-437./(R*TGAS))
  rrt(577) = F299*5.52D-30*TGAS**(-1.00)
  where( lreaction_block(:) ) rrt(:) = 0.0d0
  return
end subroutine ZDPlasKin_reac_rates
!-----------------------------------------------------------------------------------------------------------------------------------
!
! END OF FILE
!
!-----------------------------------------------------------------------------------------------------------------------------------
