program main
!
! decalre variables and modules
!
    use ZDPlasKin
    implicit none
!-----------------------------------------------------------------------------------------------------------------------------------
!
! configuration
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
!
    double precision, parameter     :: gas_temperature  = 400.0d0,  & ! gas temperature, K
                                       density_ini_ar   = 1.83d19,  & ! initial Ar density, cm-3
                                       density_ini_elec = 1.0d0       ! initial elctron density, cm-3
!
! number of pulses
!
    integer, parameter              :: ipulses = 1
!
! filename
!
    character(*), parameter         :: file_in = 'data_EN.dat', file_out = 'out.dat'
!-----------------------------------------------------------------------------------------------------------------------------------
!
! local variables
!
    integer                         :: i, j, idata_max
    double precision                :: time, dtime, EN
    double precision, allocatable   :: data_EN(:,:)
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
! write
!
    write(*,'(/,A)') 'Ar DBD Plasam Case Study'
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
! Load profiles of field and electron from external file
!

    open(1,file=file_in)
    read(1,*)

    idata_max = 0
    do while( .true. )
        read(1,*,end=11) dtime, EN
        idata_max = idata_max + 1
    enddo

11  write(*,'(x,A,i0,A)') 'found ', idata_max, ' points in ' // file_in // 'file'

    if(idata_max .le. 1) then
        write(*,*) 'wrong or missing data in the file'
        close(1)
        goto 999
    else
        allocate(data_EN(idata_max,2))
        rewind(1)
        read(1,*)
        do i = 1, idata_max
            read(1,*) data_EN(i,1:2)
        enddo
        close(1)
    endif

    time = data_EN(1,1)

! initialization of ZDPlasKin
!
    call ZDPlasKin_init()
!
! set the physical conditions of the calculation
!     the gas temperature and the reduced electric field
    call ZDPlasKin_set_conditions(GAS_TEMPERATURE=gas_temperature)
!
! set initial densities
    call ZDPlasKin_set_density(  'Ar',density_ini_ar)
    call ZDPlasKin_set_density(   'e',density_ini_elec)
    call ZDPlasKin_set_density('Ar^+',density_ini_elec)
!
! use QTPlasKin to plot the results
    call ZDPlasKin_set_config(QTPLASKIN_SAVE=.true.)
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
! open file for output and write headers
!
    open(1,file=file_out)
    write(1,'(99A14)') 'TIME', 'E/N', (trim(species_name(i)), i = 1, species_max)
    write(*,'(99A14)') 'TIME'
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
! time integration
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
    do i = 1, ipulses
        do j = 1, idata_max - 1
            dtime = data_EN(j+1,1) - data_EN(j,1)
            EN    = 0.5d0 * (data_EN(j+1,2) + data_EN(j,2))

            call ZDPlasKin_set_conditions(REDUCED_FIELD=EN)
            call ZDPlasKin_timestep(time,dtime)
            time = time + dtime

            write(1,'(99(1pe14.5))') time, EN, density(:)
        enddo

        write(*,'(99(1pe14.5))') time
    enddo
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
! close file and end
!
    close(1)

999 continue
    if (allocated(data_EN) ) deallocate(data_EN)

    write(*,'(/,A,$)') 'PRESS ENTER TO EXIT ...'
    read(*,*)
end program main
                                    