program test_extcirc

    use options
    use ZDPlasKin
    use, intrinsic :: iso_fortran_env
    implicit none

    integer                     ::  i, icount = 0                                   
    double precision            ::  time, dtime, EN, &                              ! time (s), dtime (s), EN (Td)   
                                    density_old(species_max), Vdr, &                ! density (cm-3), Vdr (cm/s)               
                                    t1, t2, tc, V, J,              &                             ! V (V), J(A)
                                    species_source(species_max), all_neutral, &     ! source (cm-3), all_neutral(cm-3)
                                    source_terms(reactions_max), varstep            ! source terms (cm-3/s)                        ! min value for the timestep
    double precision, parameter ::  time_end = 1e-4
    character(len=99)           ::  filename                            
    integer, parameter          ::  icount_save = 1

! Print
    write(*,'(/,A,/)') 'AR DBD PLASMA'
                                 
! Initialization of the program
    call CPU_TIME(t1)           ! to get the processing time at the end
    call set_parameters()       ! read and set the parameters defined in the condition.ini file
    call ZDPlasKin_init()       ! initalize the ZDPlasKin module (species, reaction, etc...))

! open and write the name of the output file
    call file_name(filename)    ! the name is composed of every important parameter of the run
    open(1,file=filename)

! set temperature and initial densities for given pressure
    call ZDPlasKin_set_conditions(GAS_TEMPERATURE=gas_temperature)
    call ZDPlasKin_set_density('Ar',gas_density)
    call ZDPlasKin_set_density('e',elec_density)
    call ZDPlasKin_set_density('Ar^+',elec_density)

! set QTPlasKin configuration
    call ZDPlasKin_set_config(QTPLASKIN_SAVE=.true.)

! initialization of variables
    time    = 0.0d0
    dtime   = 5e-7
    EN      = 0.0d0

! Tells Bolsig the value of the reduced field (to do each time EN changes)
    call ZDPlasKin_set_conditions(REDUCED_FIELD=EN)

! write the headers of the output file
    write(*,'(99(A14))') 'Time_s', 'dtime_s', 'EN_Td', 'Elec_cm-3'
    write(1,'(99(A14))', advance='no') 'Time', 'dTime', 'E/N', (trim(species_name(i)), i = 1, species_max)
    write(1,'(99(i14))', advance='no') (i, i=1, reactions_max)
    write(1,'(99(A14))') 'V', 'J', 'MuN'

! time cycling
    do while(time .lt. time_end)
!   det the new field at each step
        call ZDPlasKin_set_conditions(REDUCED_FIELD=EN)
!	Central routine of the program.
!	Here the master code call ZDPlasKin which will call DVODE knowing
!	the densities, the source terms, the electron temperature, etc...
!	The results returned are given at time = time + dtime
        call ZDPlasKin_timestep(time,dtime)
        time = time + dtime

        if (density(species_electrons) > 1.0d6) then
            call ZDPlasKin_set_density('e', 1.0d6)
            call ZDPlasKin_set_density('Ar^+', 1.0d6)
        endif
!	gives the electron drift velocity
!	(from Bolsig, knowing the field, the composition and 
!	temperature of the gas and the density of electrons)
        call ZDPlasKin_get_conditions(ELEC_DRIFT_VELOCITY=Vdr)

!	Current calculation (considering that the electrons are the only current carriers)
        J  = 1.6d-19 * disk_area * density(species_electrons) * Vdr
!	Knowing the current, deduces the voltage at the plasma
!	considering a simple voltage generator/resistance circuit
        V = maxV * epsr / (epsr + 2)

        EN = -V / gap_length / gas_density * 1.0d17 * sin(3.141592*2*frequency*time)
        EN = abs(EN) ! non-negative    
!------------------------------------

!	get the exact source terms
        call ZDPlasKin_get_rates(SOURCE_TERMS=species_source(:),REACTION_RATES=source_terms(:))
!	output
        if(mod(icount,icount_save) .eq. 0) then
!		  screen print
            write(*,'(99(1pe14.6))') time, dtime, EN, density(species_electrons)
!	write in the output file
            write(1,'(99(1pe14.6))',advance='no') time, dtime, EN,density(:)
            write(1,'(99(1pe14.6))',advance='no') source_terms(:)
            write(1,'(99(1pe14.6))') V, J
        endif

        icount = icount + 1

    enddo 	
!	end of the time loop

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

  write(*,'(A,$)') 'PRESS ENTER TO EXIT ...'
  read(*,*)

end program test_extcirc