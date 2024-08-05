!
!	In this module all the parameters of the problem are declared
! by reading the condition.ini file and the name of the output file is configured
!

module options
	implicit none
	public
	character(len=12)	::	comment						! comment in the output file name (not more than 12 letters)
	double precision	::	gas_pressure,		&		! pressure (atm)
					    	maxV,        	  	& 		! max voltage (V)
                       		epsr,      		 	&       ! relative permittivity (-)
                       		gap_length,     	&       ! gap length (cm)
							disk_area,			&		! disk area (cm2)
                       		diameter,           &      	! wall diameter (cm)
                       		wall_area,          & 		! wall area(cm2)
                       		gas_temperature,  	&       ! gas temperature (K)
					    	gas_density,		&       ! gas density (cm-3)
                       		elec_density,		&       ! initial electron density (cm-3)
							frequency					! voltage frequency (Hz)

contains

!-------------------------------------------------------------
! read the physical conditions given in the condition.ini file
!-------------------------------------------------------------

subroutine set_parameters()
	implicit none
	integer           :: io, lecteur, compteur = 1
	character(len=12) :: valeur
	open(100,file='conditions.ini')
	open(200,status='scratch')
	do
		read(100,'(1i2,1A12)',iostat=io) lecteur, valeur
		if(io .lt. 0) exit
		if(lecteur .eq. compteur) then
			write(200,'(1A12)',advance='no') trim(adjustl(valeur))
			compteur = compteur + 1
			cycle
		endif
	enddo
	rewind(200)
	read(200,'(1A12,8d12.2)')	comment, gas_pressure, gap_length, diameter, maxV, &
								epsr, gas_temperature, elec_density, frequency

!	screen print before running the code               
	write(*,'((A,1pe9.2))') 'pression, atm			=', gas_pressure,    	&
	                        'gap length, cm         =', gap_length,      	&
	                        'disk diameter, cm      =', diameter,          	&
	                        'Max voltage, V         =', maxV,     		 	&
					        'epsr, -        		=', epsr,      			&
					        'gas temperature, K     =', gas_temperature, 	&
					        'electrin density, cm-3 =', elec_density,		&
							'frequency, Hz			=', frequency  
  write(*,'(A,$)') 'PRESS ENTER TO CONTINUE ...'
	read(*,*)

!	compute the gas density knowing the pressure and the temperature
	disk_area    = 3.1415*diameter**2/4
	gas_density = 1.0d-6 * ( 101325.0d0 * gas_pressure ) / ( 1.38065d-23 * gas_temperature )
	close(200)
	close(100)
end subroutine set_parameters

!-------------------------------------------------------------
! configure the name of output file
!-------------------------------------------------------------

subroutine file_name(x)
	character(len=99), intent(out) :: x
	character(len=14)              :: b, c, d, e, f
	open(10,file='scratch.dat')
	write(10,'(3i14,2f14.1)') floor(maxV), floor(epsr), floor(gas_pressure), gap_length, diameter
	rewind(10)
	read(10,'(5A14)') b, c, d, e, f
	close(10,status='delete')
	comment = adjustl(comment)
	b = adjustl(b)
	c = adjustl(c)
	d = adjustl(d)
	e = adjustl(e)
	f = adjustl(f)
!	here is written the name of the file
	x = trim(comment)//'_'//trim(b)//'V_'//trim(c)//'epsr_'//trim(d)//'P_'//trim(e)//'d_'//trim(f)//'dia.dat'
	return
end subroutine file_name

end module options
