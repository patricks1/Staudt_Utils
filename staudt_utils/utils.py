import math
import numpy as np
import collections.abc
from astropy import units as u
from IPython.display import display, Latex

def lookup(lookupval,lookuparray,resultarray,threshold):
    '''
    Defining a lookup function that assumes lookuparray is sorted in ascending 
    order
    
    **threshold is the maximum difference between the lookupval and the value
    found in the lookup aray in order for the function to return a result. For
    example, in abundance matching, we would typically set threshold to the
    average step in numden.
    '''
    
    lookuparval,i=float('-inf'),0
    #while the values in the lookup array are less than the lookup value:  
    while lookuparval<lookupval: 
        if i>lookuparray.size-1: 
            #If the loop runs through every element of the lookuparray and 
            #lookuparval is still < lookupval, end the loop, otherwise an error 
            #will occur:
            i-=1
            #Work backwards to find i corresponding to the last time 
            #lookuparval increased:
            while lookuparray[i]==lookuparray[i-1]:
                i-=1
            #Once the nested loop has gotten i back to the point where 
            #lookuparray[i]!=lookuparray[i], break the primary loop:
            break 
        lookuparval=lookuparray[i]
        
        i+=1
        if i>1e4:
            sys.exit("Timeout")
    lowerlookuparval=lookuparray[max(0,i-2)]
    upperlookuparval=lookuparray[max(0,i-1)]
    
    #Check whether the lookuparval below the lookupval or the one above is 
    #closer to the lookupval and use the one that's closer
    upperdiff=np.abs(upperlookuparval-lookupval)
    lowerdiff=np.abs(lookupval-lowerlookuparval)
    if upperdiff>lowerdiff:
        if np.abs(lowerdiff)<threshold:
            result=resultarray[max(0,i-2)]
        else:
            #return #EXCLUDING "NA" FOR NOW
            result=None
    elif np.abs(upperdiff)<threshold:
        result=resultarray[max(0,i-1)]
    else:
        #return #EXCLUDING "NA" FOR NOW
        result=None
    
    return result

def mprint(x, d=3, show=True, order_of_mag=None):
    '''
    Generate a LaTeX formatted version of data in scientific notation

    Parameters
    ----------
    x: float
        The data point to be formatted.
    d: int, default 3
        The number of decimal places to show.
    show: bool, default True
        Whether to display the result or just return it.
    order_of_mag: int, default None
        An order of magnitude nto which the user can force the scientific
        notation. If not specified, the method determines it automatically so
        there is only one leading digit before the decimal.

    Returns
    -------
    string: str
        The LaTeX formatted data point in scientific notation, enclosed in '$' 
        characters.
    '''
    if not isinstance(x,u.quantity.Quantity):
        x = u.quantity.Quantity(x)
    v=x.value

    if order_of_mag is not None:
        exp = order_of_mag
    else:
        exp = int(math.floor(np.log10((abs(v)))))
    if exp != 0:
        order_string = '\\times10^{{{0:d}}}'.format(exp)
    else:
        order_string = ''

    fac = v/10.**exp
    unit_string = '{0:latex}'.format(x.unit)

    #replace '$' so unit_string can be easily incorporated into the larger 
    #string
    unit_string = unit_string.replace('$','') 
    
    #if x is not unitless, add a space before the unit, otherwise just get rid
    #of the string
    if unit_string!='\mathrm{}': 
        unit_string = '\ '+unit_string
    else:
        unit_string = ''

    string = '${0:0.{3}f}{2:s}{1:s}$'.format(fac, unit_string, 
                                             order_string, d)
    if show:
        display(Latex(string))
    return string 

def print_eq_old(lhs, x, d=3, op='='):
    if type(x)==u.quantity.Quantity:
        v=x.value
        unit=x.unit
        fmt='latex'
    else:
        v=x
        unit=''
        fmt='s'
    exp = int(math.floor(np.log10(abs(v))))
    fac = v/10.**exp #factor
    if exp!=0:
        string = '${0:0.{3}f}\\times10^{{{2}}}$ {1:{4}}'.format(fac, unit, exp, 
                                                                d, fmt)
    else:
        string = '${0:0.{3}f}$ {1:{4}}'.format(fac, unit, exp, d, fmt)
    display(Latex('{0} ${2}$ {1}'.format(lhs,string,op)))
    return None 

def print_eq(lhs, x, d=3, op='='):
    '''
    Display a LaTeX equation

    Parameters
    ----------
    lhs: str
        Left hand side of the equation, in LaTeX format but with no $'s
    x: float
        Value to display on right-hand side of the equation
    d: int
        Number of decimals to display
    op: str
        Operator to put between lhs and rhs

    Return
    ------
    None
    '''
    if type(x)==u.quantity.Quantity:
        v=x.value
        unit=x.unit
        fmt='latex'
    else:
        v=x
        unit=''
        fmt='s'
    exp = int(math.floor(np.log10(abs(v))))
    fac = v/10.**exp #factor
    if exp!=0:
        string = '${0:0.{3}f}\\times10^{{{2}}}$ {1:{4}}'.format(fac, unit, exp, 
                                                                d, fmt)
    else:
        string = '${0:0.{3}f}$ {1:{4}}'.format(fac, unit, exp, d, fmt)
    display(Latex('$'+lhs+op+'\:$'+string))
    return None 

def round_up(n, decimals=0):
    multiplier = 10 ** decimals
    return math.ceil(n * multiplier) / multiplier

def determine_er_symmetry(y, maxy, miny):
    '''
    Returns
    -------
    if errors are asymmetric:
        np.array([y, dy_plus, dy_minus]) 
    if errors are symmetric:
        np.array([y, dy]) 
    '''
    dy_plus = maxy-y
    dy_minus = y-miny

    _, dys = sig_figs(y, [dy_plus, dy_minus])     

    if len(dys) == 1:
        # If sig_figs found the errors to be approximately symmetric, use
        # the mean of the +error and the -error for `dy`.
        return np.array([y, np.mean([dy_plus, dy_minus])])
    else:
        return np.array([y, dy_plus, dy_minus])

def log2linear(logy, dlogy, check_symmetry=True):
    '''
    Returns
    -------
    np.array([y, *dy]), shape=(2,) if linear errors are symmetric
                              (3,) if linear errors are asymmetric
        If asymmetric, dy[0] is the plus error; dy[1] is the minus error.
    '''
    y = 10.**logy
    minlogy = logy-dlogy
    maxlogy = logy+dlogy 
    miny = 10.**minlogy
    maxy = 10.**maxlogy
    if check_symmetry:
        return determine_er_symmetry(y, maxy, miny)
    else:
        dy_plus = maxy - y
        dy_minus = y - miny
        return np.array([y, dy_plus, dy_minus])

def linear2log(y, dy, check_symmetry=True):
    '''
    Returns
    -------
    np.array([logy, *dlogy]), shape=(2,) if log-space errors are symmetric
                                    (3,) if log-space errors are asymmetric
        If asymmetric, dy[0] is the plus error; dy[1] is the minus error.
    '''
    logy = np.log10(y)
    maxy = y+dy
    miny = y-dy
    maxlogy = np.log10(maxy)
    minlogy = np.log10(miny)
    if check_symmetry:
        return determine_er_symmetry(logy, maxlogy, minlogy)
    else:
        dlogy_plus = maxlogy - logy
        dlogy_minus = logy - minlogy
        return np.array([logy, dlogy_plus, dlogy_minus])

def sig_figs(y, dys):
    '''
    Determine how many decimal places to show for a value given its 
    uncertainty.

    Parameters
    ----------
    y: float 
	The value of the variable
    dys: float or array-like 
	The uncertainty in the variable. If a float is given, the uncertainty
	is assumed to be symmetric. If an array/list is given, the 0 element
	should be the +uncertainty, the 1 element the -uncertainty.
        The galaxy name string corresponding to an index in df.

    Returns
    -------
    y_string: str
        A string representing the properly rounded target variable
    dy_strings: np.ndarray
        An array of uncertainty in y. If len(dy_strings)==1, the uncertainty
	is symmetrical. If len(dy_strings)==2, the 0 element is the 
	+uncertainty, the 1 element the -uncertainty.
    ''' 

    if not isinstance(y, (int, float)):
        raise ValueError('y should be a float or int but is a {0}. Perhaps you'
                         ' accidentally passed an element of a target vector'
                         ' instead of that element\'s value.'.format(type(y)))

    def formatter(x, decimals):
        x = round(x, decimals)
        fmt_str = '{{0:0.{0:d}f}}'.format(max(0,decimals))
        return fmt_str.format(x)
    
    # The following serves to both work with a copy of dys so we don't
    # modify the original and to make dys into an np.ndarray if it's not
    # already.
    if not isinstance(dys, (np.ndarray, list, tuple)):
        dys = np.array([dys])
    else:
        dys = np.array(dys)
    assert len(dys)<=2
    decimalss = []
    dy_strings = []
    for dy in dys:
        exp = int(math.floor(np.log10(abs(dy))))
        fac = dy/10.**exp
        if str(round(fac, 1))[0] == '1':
            # If the uncertainty would have a '1' as its first digit when we
            # show one decimal... then show one decimal. 
            exp -= 1
        decimals = -exp
        decimalss += [decimals]
        dy = formatter(dy, decimals)
        dy_strings += [dy]
    dy_strings = np.array(dy_strings)
    if len(dy_strings) > 1 and dy_strings[0] == dy_strings[1]:
        # If the uncertainties are symmetric, get rid of one of them
        dy_strings = np.array([dy_strings[0]])
    return formatter(y,max(decimalss)), dy_strings

if __name__=='__main__':
    x = 1.234e8 *u.K/u.m
    mprint(x)
