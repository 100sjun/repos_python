program run_plasRxn
!
! declare variables and modules
!
  use ZDPlasKin
  implicit none
  double precision, parameter :: gas_pressure      = 1.0d0,    & ! pressure, atm
                                 gas_temperature   = 340.15d0,   & ! gas temperature, K
                                 ini_elec_density  = 1.0d6,    & ! initial electron density, cm-3
                                 ini_gas_density   = 101325.0d0 * gas_pressure & 
                                                   / gas_temperature / 1.38d-17,  & ! gas density, cm-3
                                 pd     = 4.317d0 ! Power density, W/cm3
  double precision            :: time, time_end, dtime, EN,& ! times, s // Reduced Electric field, Td
                                 t1, t2, tc, mo, all_neutral                       ! calculation time, mobility
  integer                     :: i
!
! print
!
  write(*,'(/,A)') 'CH4 DBD Plasma Reactor'
! Initialization of the program
    call CPU_TIME(t1)           ! to get the processing time at the end!
    call ZDPlasKin_init()
!
! set the physical conditions of the calculation
!

  call ZDPlasKin_set_conditions(GAS_TEMPERATURE=gas_temperature)
  call ZDPlasKin_set_density(  'CH4',ini_gas_density)
  call ZDPlasKin_set_density(    'e',ini_elec_density)
  call ZDPlasKin_set_density('CH4^+',ini_elec_density)
  call ZDPlasKin_set_config(QTPLASKIN_SAVE=.true.)
!
! initialization of variables
    time = 0.0d0
    time_end = 2.0d2
    dtime = 1.0d-6    ! 초기 dtime 값을 1e-6으로 설정
    EN = sqrt(pd/ini_elec_density/(1.6022d-19)/3.5552d23*ini_gas_density)/ini_gas_density/(1d-17)

! Tells Bolsig the value of the reduced field (to do each time EN changes)
    call ZDPlasKin_set_conditions(REDUCED_FIELD=EN)
!
! write the headers of the output file
    write(*,'(99(A14))') 'Time_s', 'dtime_s', 'EN_Td', 'Elec_cm-3'
    write(1,'(99(A14))', advance='no') 'Time', 'dTime', 'E/N', (trim(species_name(i)), i = 1, species_max)
!
! time integration
!
  WRITE(*, '(A)') "Calculation Start!!"
  do while(time .lt. time_end)
    ! dtime 동적 조정
    if (time .lt. 1.0d-4) then
      dtime = 1.0d-6
    else if (time .lt. 1.0d-3) then
      dtime = 1.0d-5
    else if (time .lt. 1.0d-2) then
      dtime = 1.0d-4
    else if (time .lt. 1.0d-1) then
      dtime = 1.0d-3
    else if (time .lt. 1.0d0) then
      dtime = 1.0d-2
    else
      dtime = 1.0d-1
    end if

    call ZDPlasKin_set_conditions(REDUCED_FIELD=EN)
    call ZDPlasKin_timestep(time,dtime)
    time = time + dtime
    call ZDPlasKin_get_conditions(ELEC_MOBILITY_N=mo)
    call ZDPlasKin_get_density_total(ALL_NEUTRAL = all_neutral)
    EN = sqrt(pd/density(species_electrons)/(1.6022d-19)/mo*all_neutral)/all_neutral/(1d-17)

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