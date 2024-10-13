program run_plasRxn
!
! declare variables and modules
!
  use ZDPlasKin
  implicit none
  double precision, parameter :: gas_pressure      = 1.0d0,    & ! pressure, atm
                                 gas_temperature   = 300.0d0,   & ! gas temperature, K
                                 ini_elec_density  = 1.0d8,    & ! initial electron density, cm-3
                                 ini_gas_density   = 101325.0d0 * gas_pressure & 
                                                   / gas_temperature / 1.38d-17 ! gas density, cm-3
  double precision            :: time = 0.0d0, time_end = 16.0d0, dtime = 1.0d-4,  EN = 0.0d0,& ! times, s // Reduced Electric field, Td
                                 t1, t2, tc                       ! calculation time
  integer                     :: i
!
! print
!
  write(*,'(/,A)') 'CH4 DBD Plasma Reactor'
! Initialization of the program
    call CPU_TIME(t1)           ! to get the processing time at the end
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
  call ZDPlasKin_set_density(  'CH4',gas_density)
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
  close(1)
call CPU_TIME(t2)
	tc = (t2 - t1) / 60.
	if(tc .lt. 1.5d0) then
		tc = t2 - t1
		write(*,'(A,1pe12.5,A)') 'Time of calculations:', tc,' (s)'
	else
		tc = floor( (t2-t1) / 60. )
		write(*,'(A,1pe12.5,A)') 'Time of calculations:', tc,' (min)'
	endif

!
  write(*,'(/,A,$)') 'PRESS ENTER TO EXIT ...'
  read(*,*)
end program run_plasRxn