program run_plasRxn
!
! declare variables and modules
!
  use ZDPlasKin
  implicit none
  double precision, parameter :: gas_temperature  = 298.15d0, & ! gas temperature, K
                                 density_ini_ch4   = 2.46d19,  & ! initial CH4 density, cm-3
                                 density_ini_elec = 1.0d0,&      ! initial electron density, cm-3
                                 maxV             = 1.5d04,&      ! maximum electric potential, V
                                 frequency        = 1000,&       ! frequency, Hz
                                 gap_length       = 0.15       ! cm
  double precision            :: time = 0.0d0, time_end = 1.0d0, dtime = 1.0d-5,  EN = 0.0d0,& ! times, s // Reduced Electric field, Td
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
  call ZDPlasKin_set_conditions(GAS_TEMPERATURE=gas_temperature,REDUCED_FIELD=EN)
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
    call ZDPlaskin_set_conditions(REDUCED_FIELD=EN)
    call ZDPlasKin_timestep(time,dtime)
    time = time + dtime
    EN = -maxV * 10 / (10 + 2) / gap_length / density_ini_ch4 * 1.0d17 * sin(3.141592*2*frequency*time)
    EN = abs(EN) + 1.0d-1

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