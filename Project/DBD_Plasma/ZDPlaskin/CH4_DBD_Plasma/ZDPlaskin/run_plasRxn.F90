program run_plasRxn
!
! declare variables and modules
!
  use ZDPlasKin
  implicit none
  double precision, parameter :: gas_temperature  = 298.15d0, & ! gas temperature, K
                                 reduced_field    = 27.1d0,  & ! reduced electric field, Td
                                 density_ini_ch4   = 2.46d19,  & ! initial CH4 density, cm-3
                                 density_ini_elec = 1.0d0      ! initial electron density, cm-3
  double precision            :: time = 0.0d0, time_end = 1.0d-2, dtime = 1.0d-6 ! times, s
  integer                     :: i
!
! print
!
  write(*,'(/,A)') 'CH4 DBD Plasma Reactor'
!
! initialization of ZDPlasKin
!
  call ZDPlasKin_init()
!
! set the physical conditions of the calculation:
!     the gas temperature and the reduced electric field
!
  call ZDPlasKin_set_conditions(GAS_TEMPERATURE=gas_temperature,REDUCED_FIELD=reduced_field)
!
! set initial densities
!
  call ZDPlasKin_set_density(  'CH4',density_ini_ch4)
  call ZDPlasKin_set_density(    'e',density_ini_elec)
  call ZDPlasKin_set_density('CH4^+',density_ini_elec)
  call ZDPlasKin_set_config(QTPLASKIN_SAVE=.true.)
!
! print column headers and initial values
!
  write(*,'(4(A12))') 'Time_s', ( trim(species_name(i)), i = 1, species_max )
  write(*,'(4(1pe12.4))') time, density(:)
!
! time integration
!
  do while(time .lt. time_end)
    call ZDPlasKin_timestep(time,dtime)
    time = time + dtime
    write(*,'(4(1pe12.4))') time
  enddo
!
! end
!
  write(*,'(/,A,$)') 'PRESS ENTER TO EXIT ...'
  read(*,*)
end program run_plasRxn