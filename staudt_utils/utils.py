import numpy as np
from astropy import units as u
from IPython.display import display, Latex
from math import floor

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

    exp = int(floor(np.log10((abs(v)))))
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
    exp = int(floor(np.log10(abs(v))))
    fac = v/10.**exp #factor
    if exp!=0:
        string = '${0:0.{3}f}\\times10^{{{2}}}$ {1:{4}}'.format(fac, unit, exp, d, fmt)
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
    exp = int(floor(np.log10(abs(v))))
    fac = v/10.**exp #factor
    if exp!=0:
        string = '${0:0.{3}f}\\times10^{{{2}}}$ {1:{4}}'.format(fac, unit, exp, d, fmt)
    else:
        string = '${0:0.{3}f}$ {1:{4}}'.format(fac, unit, exp, d, fmt)
    display(Latex('$'+lhs+op+'\:$'+string))
    return None 

if __name__=='__main__':
    x = 1.234e8 *u.K/u.m
    mprint(x)
