import math
import numpy as np
import collections.abc
from astropy import units as u
from IPython.display import display, Latex

def lookup(lookupval,lookuparray,resultarray,threshold):
    #defining a lookup function that assumes lookuparray is sorted in ascending 
    #order
    
    #**threshold is the maximum difference between the lookupval and the value
    #found in the lookup aray in order for the function to return a result. For
    #example, in abundance matching, we would typically set threshold to the
    #average step in numden.
    
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

def mprint(x, d=3, show=True):
    if not isinstance(x,u.quantity.Quantity):
        x = u.quantity.Quantity(x)
    v=x.value

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
    dy_plus = maxy-y
    dy_minus = y-miny
    exp_plus = int(math.floor(np.log10(abs(dy_plus))))
    exp_minus = int(math.floor(np.log10(abs(dy_minus))))

    if exp_plus == exp_minus:
        fac_plus = dy_plus/10.**exp_plus
        fac_minus = dy_minus/10.**exp_minus
        if round(fac_plus) == round(fac_minus):
            return np.array([y, np.mean([dy_plus, dy_minus])])
    return np.array([y, dy_plus, dy_minus])

def log2linear(logy, dlogy):
    y = 10.**logy
    minlogy = logy-dlogy
    maxlogy = logy+dlogy 
    miny = 10.**minlogy
    maxy = 10.**maxlogy
    return determine_er_symmetry(y, maxy, miny)

def linear2log(y, dy):
    logy = np.log10(y)
    maxy = y+dy
    miny = y-dy
    maxlogy = np.log10(maxy)
    minlogy = np.log10(miny)
    return determine_er_symmetry(logy, maxlogy, minlogy)

def sig_figs(y, dys):
    def formatter(x, decimals):
        x = round(x, decimals)
        fmt_str = '{{0:0.{0:d}f}}'.format(max(0,decimals))
        return fmt_str.format(x)

    dys = dys.copy()
    if not isinstance(dys, (collections.abc.Sequence, np.ndarray)):
        dys = np.array([dys])
    assert len(dys)<=2
    decimalss = []
    dy_strings = []
    for dy in dys:
        exp = int(math.floor(np.log10(abs(dy))))
        fac = dy/10.**exp
        if round(fac,0) == 1:
            exp -= 1
        decimals = -exp
        decimalss += [decimals]
        dy = formatter(dy, decimals)
        dy_strings += [dy]
    #decimalss = np.array(decimalss)
    #dys = np.array([round(dy,decimals) for dy,decimals in zip(dys,decimalss)])
    dy_strings = np.array(dy_strings)
    return formatter(y,max(decimalss)), dy_strings

if __name__=='__main__':
    x = 1.234e8 *u.K/u.m
    mprint(x)
