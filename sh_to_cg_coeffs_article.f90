! ----------------------------------------------------------------------------------------
! 
! Commands compile the code :
! > gfortran -Wall -pedantic sh_to_cg_coeffs_article.f90 -o sh_to_cg_coeffs_article
!
! ----------------------------------------------------------------------------------------

! ----------------------------------------------------------------------------------------
!
! Module to define the floating point precision
!
module numbers

    implicit none

    integer, parameter :: float = selected_real_kind( 13, 100 )

end module numbers
!
! ----------------------------------------------------------------------------------------


! ----------------------------------------------------------------------------------------
!
! Module to define the input / output units
!
module iounits

    implicit none

    integer, parameter :: iout = 6

end module iounits
!
! ----------------------------------------------------------------------------------------


! ----------------------------------------------------------------------------------------
!
! Module where the main formulae are implemented
!
! -- Formula to compute the factorial of a non-negative integer number
!    ( function fact )
!
! -- Formula to compute the binomial coefficient of two non-negative integer numbers 
!    ( function binom )
!
! -- Formula for the C_lm_tuv expansion coefficients 
!    ( subroutine clmtuv_formula )
!    that computes the coefficients of the expansion that defines the solid harmonics
!    functions transformation from spherical to Cartesian coordinates
!
module formulae

    use numbers, only : float
    
    implicit none
        
    public :: clmtuv_formula

contains

    ! Function to compute the factorial of a non-negative integer number
    !
    ! Input argument  : 
    !
    ! -- enne  ( non-negative integer )
    !
    ! Output argument :
    !
    ! -- positive real value ( factorial ) that contains the factorial
    !    of the non-negative integer number enne ( n ) --> n!
    ! 
    function fact(enne) result(factorial)
   
        use iounits, only : iout

        implicit none

        integer, intent(in) :: enne
       
        integer :: ind

        real(float) :: factorial
    
        ! check that the input argument is well defined

        if( enne < 0 ) then
            
            write( iout, 10 )
            stop
            
        endif

        ! calculation of the factorial of non-negative integer number n

        factorial = 1.0_float    ! basic cases n = 0 or n = 1 --> 0! = 1! = 1

        if( enne > 1 ) then      ! non-basic cases n > 1 --> n! = n . ( n - 1 ) ... 1

            do ind = 1, enne

                factorial = factorial * real( ind, kind = float )

            enddo

        endif

        return

        10 format(/" Error. Attempt to compute the factorial of a negative number"/)

    end function fact


    ! Function to compute the binomial coefficients of two non-negative integer numbers
    !
    ! Input arguments : 
    !
    ! -- enne  ( upper non-negative integer )
    ! -- kappa ( lower non-negative integer )
    !
    ! Output argument :
    !
    ! -- positive real value ( binomial ) that contains the binomial coefficient
    !    of two non-negative integer numbers enne ( n ) and kappa ( k ) --> nCk
    ! 
    function binom(enne, kappa) result(binomial)
        
        use iounits, only : iout
    
        implicit none

        integer, intent(in) :: enne
        integer, intent(in) :: kappa
       
        logical :: error

        real(float) :: binomial
    
        ! check that the input arguments are well defined

        error = .false.

        if( kappa < 0 .or. enne < 0 ) then

            write( iout, 10)
            error = .true.

        endif
        
        if( kappa > enne ) then
            
            write( iout, 12)
            error = .true.
            
        endif

        if( error ) stop
       
        ! calculation of the binomial coefficient nCk = n! / ( k! (n-k)! )
         
        binomial = fact(enne) / ( fact(kappa) * fact(enne - kappa) )

        return

        10 format(/" Error. Attempt to compute binomial coefficient for k < 0 or n < 0"/)
        12 format(/" Error. Attempt to compute binomial coefficient for k > n"/)

    end function binom


    ! Function to compute the coefficients C_lm_tuv of the expansion
    !
    ! Y_lm(r) = sum_{tuv}^l C_lm_tuv x^t y^u z^v      where t + u + v = l
    !
    ! that permits to convert the solid harmonics Y_lm(r) from spherical to Cartesian coordinates
    !
    ! Input arguments : 
    !
    ! -- elle  ( integer quantum number l )
    ! -- emme  ( integer quantum number m )
    ! -- tcart ( integer number t power of the variable x )
    ! -- ucart ( integer number u power of the variable y )
    ! -- vcart ( integer number v power of the variable z )
    !
    ! Output :
    !
    ! -- expansion coefficient C_lm_tuv
    !    for a given set of integer numbers ( l m t u v ) 
    !    
    function clmtuv_formula(elle, emme, tcart, ucart, vcart) result(clmtuv)
        
        use iounits, only : iout

        implicit none

        integer, intent(in) :: elle
        integer, intent(in) :: emme
        integer, intent(in) :: tcart
        integer, intent(in) :: ucart
        integer, intent(in) :: vcart

        complex(float) :: clmtuv

        logical        :: is_zero
        integer        :: iend, jend
        integer        :: jend_num
        integer        :: lower_bin
        integer        :: iloop, jloop
        real   (float) :: bin_one, bin_two
        complex(float) :: isum, jsum
        complex(float) :: imag_unit, factor

        if( tcart + ucart + vcart .ne. elle ) then

            write( iout, 10 )
            stop

        endif

        imag_unit = ( 0.0_float, 1.0_float )

        iend = ( elle - abs(emme) ) / 2

        jend_num = tcart + ucart - abs(emme)

        if( mod( jend_num, 2 ) .ne. 0 ) then

            clmtuv = ( 0.0_float, 0.0_float )

        else

            jend = jend_num / 2
    
            isum = ( 0.0_float, 0.0_float )
    
            do iloop = 0, iend
    
                jsum = ( 0.0_float, 0.0_float )
    
                do jloop = 0, jend

                    lower_bin = tcart - 2 * jloop

                    is_zero = lower_bin < 0 .or. lower_bin > abs(emme)

                    if( .not. is_zero ) then
                
                        factor = imag_unit ** ( abs(emme) - lower_bin )
   
                        bin_one = binom( jend, jloop )

                        bin_two = binom( abs(emme), lower_bin )

                        jsum = jsum + bin_one * bin_two * factor

                    endif
    
                enddo

                is_zero = jend < 0 .or. jend > iloop

                if( .not. is_zero ) then
                
                    factor = ( - 1.0_float ) ** iloop * fact( 2 * elle - 2 * iloop )
                    
                    factor = factor / fact( elle - abs(emme) - 2 * iloop )

                    bin_one = binom( elle, iloop )

                    bin_two = binom( iloop, jend )
    
                    isum =  isum + ( bin_one * bin_two * factor ) * jsum

                endif
    
            enddo
    
            isum = isum / ( 2.0_float ** elle * fact( elle ) )
            
            clmtuv = isum

        endif

        return

        10 format(/" Error in input arguments : t + u + v differs from elle"/)

    end function clmtuv_formula
    
    
    ! Function to print the coefficients C_lm_tuv computed by function clmtuv_formula
    !
    ! Input arguments : 
    !
    ! -- elle      ( integer quantum number l )
    ! -- emme      ( integer quantum number m )
    ! -- tcart     ( integer number t power of the variable x )
    ! -- ucart     ( integer number u power of the variable y )
    ! -- vcart     ( integer number v power of the variable z )
    ! -- real_harm ( logical true to obtain coefficients for real solid harmonics )
    !
    subroutine print_clmtuv(elle, emme, tcart, ucart, vcart, real_harm)
        
        use iounits, only : iout

        implicit none

        integer, intent(in) :: elle
        integer, intent(in) :: emme
        integer, intent(in) :: tcart
        integer, intent(in) :: ucart
        integer, intent(in) :: vcart
        logical, intent(in) :: real_harm

        complex(float) :: clmtuv

        if( tcart + ucart + vcart .eq. elle ) then

            ! calculation of the coefficient for a given set ( l m t u v )

            clmtuv = clmtuv_formula( elle, emme, tcart, ucart, vcart )
           
            ! print the coefficient for a given set ( l m t u v )

            if( real_harm ) then ! coefficients for real solid harmonics
                
                if( emme .ge. 0 ) then ! m > 0 or m = 0

                    write( iout, 10 ) elle, emme, tcart, ucart, vcart, real(clmtuv)

                else if( emme .lt. 0 ) then ! m < 0

                    write( iout, 10 ) elle, emme, tcart, ucart, vcart, aimag(clmtuv)

                endif

            else ! coefficients for complex solid harmonics

                if( emme .lt. 0 ) then ! m < 0

                    write( iout, 12 )

                else ! m > 0  or  m = 0

                    write( iout, 14 ) elle, emme, tcart, ucart, vcart, clmtuv

                endif

            endif
            
        endif ! t + u + v = l

        return

        10 format(/" D_lm(tuv) coefficient for  ( l m ) = ", 2(i0, 1x), t44, &
                 & " and  ( t u v ) = ", 3(i0, 1x), t69, " :: ", *(es22.14) /)
        12 format(/" Coefficients for complex harmonics are computed for positive emme only")
        14 format(/" C_lm(tuv) coefficient for  ( l m ) = ", 2(i0, 1x), t44, &
                 & " and  ( t u v ) = ", 3(i0, 1x), t69, " :: ", *(es22.14) /)

    end subroutine print_clmtuv


end module formulae
!
! ----------------------------------------------------------------------------------------


! ----------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------
!
program clmtuv_conversion_coeffs

    use iounits , only : iout
    use formulae, only : clmtuv_formula, print_clmtuv

    logical             :: real_harm
    character(len = 32) :: arg
    integer             :: narg
    integer             :: elle_max
    integer             :: elle, emme
    integer             :: tcart, ucart, vcart

    real_harm = .false. ! complex solid harmonics is the default

    ! write the header

    write( iout, 10 ) ! write the header
    
    ! check the correctness of the number of input arguments

    narg = command_argument_count()
   
    if( narg .ne. 1 .and. narg .ne. 2 .and. narg .ne. 5 .and. narg .ne. 6 ) then

        write( iout, 12 ) 
        stop

    endif

    ! case number of arguments = 1 or = 2
    ! --> calculation of all the coefficients in the range l = [ 0 .... elle_max ]
    ! --> if narg == 1 --> coefficients for complex solid harmonics
    ! --> if narg == 2 --> coefficients for real    solid harmonics
    
    if( narg == 1 .or. narg == 2 ) then
   
        ! get the input arguments and check their correctness

        call get_command_argument( 1, arg )
        read( arg, '(i10)' ) elle_max

        if( elle_max < 0 ) then

            write( iout, 14 ) 
            stop

        endif
       
        if( narg == 2 ) then 

            call get_command_argument( 2, arg )

            if( trim( adjustl( arg ) ) == 'real' ) then
                real_harm = .true.
            else
                write( iout, 16 )
                stop
            endif

        endif

        ! calculation and print of the coefficients in the range l = [ 0 .... elle_max ]

        do elle = 0, elle_max
            do emme = -elle, elle
                do tcart = 0, elle
                    do ucart = 0, elle
                        do vcart = 0, elle
                            
                           ! print of the coefficient for a given set ( l m t u v )

                           call print_clmtuv( elle, emme, tcart, ucart, vcart, real_harm )

                        enddo ! vcart
                    enddo ! ucart
                enddo ! tcart
            enddo ! emme
        enddo ! elle

    ! case number of arguments = 5 or = 6
    ! --> calculation of the coefficient for a given set ( l m t u v )
    ! --> if narg == 5 --> coefficients for complex solid harmonics
    ! --> if narg == 6 --> coefficients for real    solid harmonics

    else if( narg == 5 .or. narg == 6 ) then

        ! get the input arguments and check their correctness
        
        call get_command_argument( 1, arg )
        read( arg, '(i10)' ) elle

        if( elle < 0 ) then

            write( iout, 18 )
            stop

        endif
        
        call get_command_argument( 2, arg )
        read( arg, '(i10)' ) emme
        
        if( abs(emme) .gt. elle ) then

            write( iout, 20 )
            stop

        endif
        
        call get_command_argument( 3, arg )
        read( arg, '(i10)' ) tcart
        
        call get_command_argument( 4, arg )
        read( arg, '(i10)' ) ucart
        
        call get_command_argument( 5, arg )
        read( arg, '(i10)' ) vcart

        if( tcart < 0 .or. ucart < 0 .or. vcart < 0 ) then

            write( iout, 22 )
            stop

        endif
        
        if( tcart + ucart + vcart .ne. elle ) then

            write( iout, 24 )
            stop

        endif
        
        if( narg == 6 ) then 

            call get_command_argument( 6, arg )

            if( trim( adjustl( arg ) ) == 'real' ) then
                real_harm = .true.
            else
                write( iout, 26 )
                stop
            endif

        endif
    
        ! print of the coefficient for a given set ( l m t u v )

        call print_clmtuv( elle, emme, tcart, ucart, vcart, real_harm )
    
    endif
    
    write( iout, '( / 120("-") / )' )

10 format(/ 120("-") / 120("-")                                                  / &
   & / 11x, "Calculation of Conversion Coefficients from Spherical to Cartesian ", &
   &  "Coordinates for Solid Harmonics", 11x                                     / &
   & / 120("-") / 120("-"))

12 format(/ " ** Complex Solid Harmonics >> "                                          / &
   & / "    --> coefficient for a given l = elle and m = emme and a triplet ( t u v ) "/ &
   & / "        usage :  ./exe  elle  emme  t  u  v "                                  / &
   & / "    --> coefficients in the range l = [ 0 ... elle_max ] "                     / &
   & / "        usage :  ./exe  elle_max "                                             / &
   & / " ** Real Solid Harmonics >> "                                                  / &
   & / "    --> coefficient for a given l = elle and m = emme and a triplet ( t u v ) "/ &
   & / "        usage :  ./exe  elle  emme  t  u  v  real "                            / &
   & / "    --> coefficients in the range l = [ 0 ... elle_max ] "                     / &
   & / "        usage :  ./exe  elle_max  real "                                       / &
   & / 120("-") /)
 
14 format(/ " Error in input arguments : elle_max must be greater or equal than zero " /)
16 format(/ " Error in input arguments : allowed value of the second argument : real " /)
18 format(/ " Error in input arguments : elle must be greater or equal than zero " /)
20 format(/ " Error in input arguments : |emme| exceeds elle " /)
22 format(/ " Error in input arguments : t u v must be greater or equal than zero " /)
24 format(/ " Error in input arguments : t + u + v differs from elle " /)
26 format(/ " Error in input arguments : allowed value of the second argument : real " /)

end program clmtuv_conversion_coeffs
!
! ----------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------

