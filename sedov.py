from __future__ import division
from scipy.integrate import quad
import math

#------------------------------------------------------------------------------

gamma = 1.4

E0 = 1.0                    #energy input
rho_0 = 1.0                 #initial constant density

time_end = 4.0              #final time
delta_time = 0.01           #time resolution

rmodel = 2.0                #radius of the model
delta_r = rmodel/1000       #radial resolution

#------------------------------------------------------------------------------
#numbers of the equations are indicated as in the work of James R. Kamm

#equations below correspond to the 'standard' case, constant initial density and
#spherical symmetry

#equations (33)-(37)
a = 1.25*(gamma+1.0)

b = (gamma+1.0)/(gamma-1.0)

c = 2.5*gamma

d = 5.0*(gamma+1.0)/(5.0*(gamma+1.0)-2.0*(2.0+3.0*(gamma-1.0)))

e = 1.0 + 1.5*(gamma-1.0)

#equations (42)-(47)
alpha0 = 0.4

alpha2 = - (gamma-1.0)/(2.0*(gamma-1)+3.0)

alpha1 = 5.0*gamma/(2.0+3.0*(gamma-1.0)) * (6.0*(2.0-gamma)/(gamma*25.0)-alpha2)

alpha3 = 3.0/(2.0*(gamma-1.0)+3.0)

alpha4 = 5.0*alpha1/(2.0-gamma)

alpha5 = -2.0/(2.0-gamma)


V2 = 0.8 / (gamma+1.0) # after equation (17)

V0 = 0.4 / gamma       # equation (23)


def Integrand1(V): # equation (73)
    x1 = a*V
    x2 = b*(c*V-1)
    x3 = d*(1-e*V)
    x4 = b*(1-c*V/gamma)
    return ( -(gamma+1.0)/(gamma-1.0) * V**2 * (alpha0/V + alpha2*c/(c*V-1.0)
        - alpha1*e/(1.0-e*V) )
        * (x1**alpha0 * x2**alpha2 * x3**alpha1)**(-5.0)
        * x2**alpha3 * x3**alpha4 * x4**alpha5 )

def Integrand2(V): # equation (74)
    x1 = a*V
    x2 = b*(c*V-1)
    x3 = d*(1-e*V)
    x4 = b*(1-c*V/gamma)
    return ( -(gamma+1.0)/(2.0*gamma) * V**2 * (c*V-gamma)/(1-c*V)
        * (alpha0/V + alpha2*c/(c*V-1.0) - alpha1*e/(1.0-e*V) )
        * (x1**alpha0 * x2**alpha2 * x3**alpha1)**(-5.0)
        * x2**alpha3 * x3**alpha4 * x4**alpha5 )

result1, err1 = quad(Integrand1,V0+1.0e-16,V2)

result2, err2 = quad(Integrand2,V0+1.0e-16,V2)

J1 = result1    #equation (67)

J2 = result2    #equation (68)

alpha = 2.0*math.pi*result1 + 4.0/(gamma-1.0)*math.pi*result2 #equation (66)

def dm(Vp):
    x1 = a*Vp
    x2 = b*(c*Vp-1)
    x3 = d*(1-e*Vp)
    x4 = b*(1-c*Vp/gamma)
    return ( 4.0*math.pi* r2**3 *rho2 * x1**(-3.0*alpha0)
            * x2**(alpha3-3.0*alpha2)
            * x3**(alpha4-3.0*alpha1) * x4**alpha5
            * (-a*alpha0*x1**(-1) - b*c*alpha2*x2**(-1) + d*e*alpha1*x3**(-1)) )

def m(V): # m(r) = \int_0^r 4 * pi * r'**2 * rho' * dr'
    result_m, err_m = quad(dm,V0+1.0e-10,V)
    return result_m

def diffm(V,mx):
    result_dm, err_dm = quad(dm,V0+1.0e-10,V)
    return result_dm - mx

def r(V): #equation (38)
    x1 = a*V
    x2 = b*(c*V-1)
    x3 = d*(1-e*V)
    return r2 * x1**(-alpha0) * x2**(-alpha2) * x3**(-alpha1)

def diffr(V,rx):
    x1 = a*V
    x2 = b*(c*V-1)
    x3 = d*(1-e*V)
    return r2 * x1**(-alpha0) * x2**(-alpha2) * x3**(-alpha1) - rx

def rho(V): #equation (40)
    x1 = a*V
    x2 = b*(c*V-1)
    x3 = d*(1-e*V)
    x4 = b*(1-c*V/gamma)
    return rho2 * x2**alpha3 * x3**alpha4 * x4**alpha5

#bisection function is based on the one of Praveen Sridhar, taken from
#psbots.blogspot.com
def bisection(function,second_arg,lower_bound,upper_bound,
                                                    error_tolerance=1.0e-16):
    fu=function(upper_bound,second_arg)
    if fu==0.0 :  return upper_bound
    fl=function(lower_bound,second_arg)
    if fl==0.0 :  return lower_bound
    if fu*fl > 0.0:  #if the product of the function values at the two bounds
                     #is positive, then they do not bracket a root
        return None
        print "the given values do not bracket a root"
    else:
        mid=(lower_bound+upper_bound)/2
        error=(upper_bound-lower_bound)/2
        while(error>error_tolerance): #the loop executes till the desired level
                                      #of accuracy is attained
            fm=function(mid,second_arg)
            fu=function(upper_bound,second_arg)
            fl=function(lower_bound,second_arg)

            #if the middle value is indeed the root, it is returned as the root
            if fm==0.0 : return mid

            #if the actual root lies between upper bound and the middle value,
            #the lower bound is set to be the middle value
            elif fu*fm < 0.0 : lower_bound=mid

            #if the actual root lies between lower bound and the middle value,
            #the upper bound is set to be the middle value
            elif fl*fm < 0.0 : upper_bound=mid

            mid=(lower_bound+upper_bound)/2
            error=(upper_bound-lower_bound)/2
        return mid


# here we solve for density as a function of radius in the range between
#                                               0.5*r_shock and r_shock (= r2)
# below 0.5*r_shock density is taken to be zero
# above r_shock density is taken to be initial

outfile = open("sedov_rho_rad.xg","w")

time = 0.0
while time < time_end:
    outfile.write('"Time = ' + str(time) + '\n')
    r2 = (E0/rho_0/alpha)**0.2 * time**0.4   #equation (14)
    rho1 = rho_0   #equation (5)
    rho2 = rho1 * (gamma + 1.0)/(gamma - 1.0)   #equation (13)
    r_coordinate = 0.0
    while r_coordinate < rmodel:
        if r_coordinate < 0.5*r2:
            outfile.write(str(r_coordinate) + ' ' + str(0.0) + '\n')
        elif r_coordinate >= 0.5*r2 and r_coordinate < r2:
            V_x = bisection(diffr,r_coordinate,V0,V2)
            outfile.write(str(r_coordinate) + ' ' + str(rho(V_x)) + '\n')
        elif r_coordinate > r2:
            outfile.write(str(r_coordinate) + ' ' + str(rho1) + '\n')
        r_coordinate = r_coordinate + delta_r
    outfile.write('\n')
    outfile.write('\n')
    time = time + delta_time
outfile.close()
