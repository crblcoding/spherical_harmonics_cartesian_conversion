# ----------------------------------------------------------------------------------------------------------------------------------
# 
# Commands compile and execute the code :
# > python3 sh_to_cg_coeffs_article.py
#
# ----------------------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------------------
#
# Function to compute the factorial of a non-negative integer number
#
# Input argument  :
#
# -- enne  ( non-negative integer )
#
# Output argument :
#
# -- positive real value ( factorial ) that contains the factorial
#    of the non-negative integer number enne ( n ) --> n!
#
def fact(enne):

    import sys

    if enne < 0:

        print(" Error. Attempt to compute the factorial of a negative number\n")
        sys.exit()

    factorial = 1.0  # basic cases n = 0 or n = 1 --> 0! = 1! = 1

    if enne > 1:     # non-basic cases n > 1 --> n! = n . ( n - 1 ) ... 1

        for index in range(1, enne + 1):

            factorial = factorial * float(index)

    return factorial
#
# ----------------------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------------------
#
# Function to compute the binomial coefficients of two non-negative integer numbers
#
# Input arguments : 
#
# -- enne  ( upper non-negative integer )
# -- kappa ( lower non-negative integer )
#
# Output argument :
#
# -- positive real value ( binomial ) that contains the binomial coefficient
#    of two non-negative integer numbers enne ( n ) and kappa ( k ) --> nCk
# 
#
def binom(enne, kappa):

    import sys

    if kappa < 0 or enne < 0:

        print(" Error. Attempt to compute binomial coefficient for k < 0 or n < 0 \n")
        sys.exit()
    
    if kappa > enne:

        print(" Error. Attempt to compute binomial coefficient for k > n \n")
        sys.exit()
    
    # calculation of the binomial coefficient nCk = n! / ( k! (n-k)! )

    binomial = fact(enne) / ( fact(kappa) * fact(enne - kappa) )

    return binomial
#
# ----------------------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------------------
#
# Function to compute the coefficients C_lm_tuv of the expansion
#
# Y_lm(r) = sum_{tuv}^l C_lm_tuv x^t y^u z^v      where t + u + v = l
#
# that permits to convert the solid harmonics Y_lm(r) from spherical to Cartesian coordinates
#
# Input arguments : 
#
# -- elle  ( integer quantum number l )
# -- emme  ( integer quantum number m )
# -- tcart ( integer number t power of the variable x )
# -- ucart ( integer number u power of the variable y )
# -- vcart ( integer number v power of the variable z )
#
# Output :
#
# -- expansion coefficient C_lm_tuv
#    for a given set of integer numbers ( l m t u v ) 
#    
def clmtuv_formula(elle, emme, tcart, ucart, vcart):

    import sys

    if tcart + ucart + vcart != elle:

        print(" Error in input arguments : t + u + v differs from elle \n")
        sys.exit()
    
    imag_unit = complex( 0.0, 1.0 )

    iend = ( elle - abs(emme) ) // 2

    jend_num = tcart + ucart - abs(emme)
    
    if jend_num % 2 != 0:

        clmtuv = complex( 0.0, 0.0 )

    else:

        jend = jend_num // 2

        isum = complex( 0.0, 0.0 )
        
        for iloop in range(iend + 1):

            jsum = complex( 0.0, 0.0 )

            for jloop in range(jend + 1):

                lower_bin = tcart - 2 * jloop

                is_zero = lower_bin < 0 or lower_bin > abs(emme)

                if not is_zero:

                    factor = imag_unit ** ( abs(emme) - lower_bin )

                    bin_one = binom(jend, jloop)

                    bin_two = binom( abs(emme), lower_bin )

                    jsum = jsum + bin_one * bin_two * factor
            
            is_zero = jend < 0 or jend > iloop

            if not is_zero:

                factor = ( -1.0 ) ** iloop * fact( 2 * elle - 2 * iloop )

                factor = factor / fact( elle - abs(emme) - 2 * iloop )

                bin_one = binom( elle, iloop )

                bin_two = binom( iloop, jend )

                isum = isum + ( bin_one * bin_two * factor ) * jsum
        
        isum = isum / ( 2.0 ** elle * fact(elle) )

        clmtuv = isum
    
    return clmtuv
#
# ----------------------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------------------
#
# Function to print the coefficients C_lm_tuv computed by function clmtuv_formula
#
# Input arguments :
#
# -- elle      ( integer quantum number l )
# -- emme      ( integer quantum number m )
# -- tcart     ( integer number t power of the variable x )
# -- ucart     ( integer number u power of the variable y )
# -- vcart     ( integer number v power of the variable z )
# -- real_harm ( logical true to obtain coefficients for real solid harmonics )
#
def print_clmtuv(elle, emme, tcart, ucart, vcart, real_harm):

    import sys

    if tcart + ucart + vcart == elle:

        # calculation of the coefficient for a given set ( l m t u v )

        clmtuv = clmtuv_formula( elle, emme, tcart, ucart, vcart )
        
        # print the coefficient for a given set ( l m t u v )

        if real_harm:  # coefficients for real solid harmonics

            if emme >= 0:    # m > 0 or m = 0

                print(f" D_lm(tuv) coefficient for  ( l m ) = {elle} {emme}   and  ( t u v ) = {tcart} {ucart} {vcart}   ::   {clmtuv.real:.14e} \n")

            elif emme < 0:   # m < 0
                
                print(f" D_lm(tuv) coefficient for  ( l m ) = {elle} {emme}   and  ( t u v ) = {tcart} {ucart} {vcart}   ::   {clmtuv.imag:.14e} \n")

        else:  # coefficients for complex solid harmonics

            if emme < 0:     # m < 0

                print(" Coefficients for complex harmonics are computed for positive emme only\n")

            else:            # m > 0 or m = 0 

                print(f" C_lm(tuv) coefficient for  ( l m ) = {elle} {emme}   and  ( t u v ) = {tcart} {ucart} {vcart}   ::   {clmtuv.real:.14e}  {clmtuv.imag:.14e}\n")
#
# ----------------------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------
#
# Main program begin
#
def clmtuv_conversion_coeffs():

    import sys

    real_harm = False  # complex solid harmonics is the default

    # write the header

    print("\n" + "-" * 120)
    print("-" * 120)
    print("\n" + " " * 10, "Calculation of Conversion Coefficients from Spherical to Cartesian Coordinates for Solid Harmonics\n")
    print("-" * 120)
    print("-" * 120 + "\n")
    
    # check the correctness of the number of input arguments

    narg = len(sys.argv) - 1

    if narg not in {1, 2, 5, 6}:

        print(" ** Complex Solid Harmonics >> ")
        print("\n    --> coefficient for a given l = elle and m = emme and a triplet ( t u v ) ")
        print("\n        usage :  ./exe  elle  emme  t  u  v ")
        print("\n    --> coefficients in the range l = [ 0 ... elle_max ] ")
        print("\n        usage :  ./exe  elle_max ")
        print("\n ** Real Solid Harmonics >> ")
        print("\n    --> coefficient for a given l = elle and m = emme and a triplet ( t u v ) ")
        print("\n        usage :  ./exe  elle  emme  t  u  v  real ")
        print("\n    --> coefficients in the range l = [ 0 ... elle_max ] ")
        print("\n        usage :  ./exe  elle_max  real \n")
        print("-" * 120 + "\n")

        return

    # case number of arguments = 1 or = 2
    # --> calculation of all the coefficients in the range l = [ 0 .... elle_max ]
    # --> if narg == 1 --> coefficients for complex solid harmonics
    # --> if narg == 2 --> coefficients for real    solid harmonics

    if narg in {1, 2}:

        # get the input arguments and check their correctness

        elle_max = int(sys.argv[1])

        if elle_max < 0:

            print(" Error in input arguments : elle_max must be greater or equal than zero \n")
            return

        if narg == 2:

            if sys.argv[2].strip().lower() == 'real':

                real_harm = True

            else:

                print(" Error in input arguments : allowed value of the second argument : real \n")
                return

        # calculation and print of the coefficients in the range l = [ 0 .... elle_max ]

        for elle in range(elle_max + 1):
            for emme in range(-elle, elle + 1):
                for tcart in range(elle + 1):
                    for ucart in range(elle + 1):
                        for vcart in range(elle + 1):

                            # print of the coefficient for a given set ( l m t u v )

                            print_clmtuv( elle, emme, tcart, ucart, vcart, real_harm )

    # case number of arguments = 5 or = 6
    # --> calculation of the coefficient for a given set ( l m t u v )
    # --> if narg == 5 --> coefficients for complex solid harmonics
    # --> if narg == 6 --> coefficients for real    solid harmonics

    elif narg in {5, 6}:

        # get the input arguments and check their correctness

        elle = int(sys.argv[1])

        if elle < 0:

            print(" Error in input arguments : elle must be greater or equal than zero \n")
            return

        emme = int(sys.argv[2])

        if abs(emme) > elle:

            print(" Error in input arguments : |emme| exceeds elle \n")
            return

        tcart = int(sys.argv[3])
        ucart = int(sys.argv[4])
        vcart = int(sys.argv[5])

        if tcart < 0 or ucart < 0 or vcart < 0:

            print(" Error in input arguments : t u v must be greater or equal than zero \n")
            return

        if tcart + ucart + vcart != elle:

            print(" Error in input arguments : t + u + v differs from elle \n")
            return

        if narg == 6:

            if sys.argv[6].strip().lower() == 'real':

                real_harm = True

            else:

                print(" Error in input arguments : allowed value of the second argument : real \n")
                return

        # print of the coefficient for a given set ( l m t u v )

        print_clmtuv( elle, emme, tcart, ucart, vcart, real_harm )

    print("-" * 120 + "\n")


if __name__ == "__main__":

    clmtuv_conversion_coeffs()

#
# Main program end
#
# ----------------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------
