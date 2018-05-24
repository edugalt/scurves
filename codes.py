import numpy as np
import pylab as pl
from scipy.optimize import  fmin_l_bfgs_b
import string


def xcurve_5(x,x0,x1,x2,a,b):
    '''Solution to the differential equation
       d/dt rho = (x2 - rho)*(a+b*(rho-x1)) if rho in [x1,x2] and 0.0 else;
       with initial condition rho(x0) = 1/2*(y0+y1); s
       this curve goes from x1 to x2;
       note the point x0 is not the intial condition rho(x0) = x1 since this is ill-defined for a=0 (s-curve).
       the special cases are exp-curve (b=0.0) and s-curve (a=0.0)
       
       if a>0, there is a finite time t0: rho(t0)=x1, such that
       rho(t) = x1 for t<t0 and rho(t) = xcurve for t>=t0
       if a==0, rho(t-->-infty) --> x1 <--> t0 =infty
       
       In: x ... array for times at which to compute curve
           x0 ... parameter where rho(x0)=1/2*(y0+y1) 'half the change'
           x1 ... lower asymptotic
           x2 ... upper asymptotic
           a ... external influence
           b ... internal influence
       
    '''
    f_tmp = np.exp((a+b*(x2-x1))*(x-x0))
    if a>0:
        t0 = (x0*(a+b*(x2-x1)) - np.log((2.0*a+b*(x2-x1))/a))/(a+b*(x2-x1))
        y = (x<t0)*x1 + (x>=t0)*(-(a-b*x1)*(x2-x1) + x2*(2.0*a+b*(x2-x1))*f_tmp )/( b*(x2-x1) + (2.0*a+b*(x2-x1))*f_tmp )
    else:
        y = ( -(a-b*x1)*(x2-x1) + x2*(2.0*a+b*(x2-x1))*f_tmp )/( b*(x2-x1) + (2.0*a+b*(x2-x1))*f_tmp )
    return y
    
def fit_mixed(x,y,sigma,theta_0):
    '''fits the mixed model to data (x,y,sigma) with initial parameter theta_0 numerically (local minimum!).
       returns: [estimated parameters, chi_squared of best fit]
    '''
    bnds = [(None,None),(0.0,1.0),(0.0,1.0),(0.0,None),(0.0,None)]
    if sigma.any() == None:
        sigma = np.ones(len(x))
    result_tmp = fmin_l_bfgs_b(ls_fit_mixed,theta_0,args=(x,y,sigma),bounds=bnds, approx_grad=True,disp=0,iprint=-1)
    chi2 = 0
    if result_tmp[2]['warnflag']==0:
        result = result_tmp[0]
        chi2 = ls_fit_mixed(result[:5],x,y,sigma)
    else:
        chi2 = 0
        result = [0,0,0,0,0]
        print('Optimization did not converge - try different initial conditions for the parameters!')
    return [result,chi2]

def ls_fit_mixed(theta,x,y,sigma):
    '''calculates the chi_squared function for the mixed model with parameters theta compared to data (x,y,sigma).
x    '''
    x0 = theta[0]
    x1 = theta[1]
    x2 = theta[2]
    a = theta[3]
    b = theta[4]
    Delta = np.sum(((y - xcurve_5(x,x0,x1,x2,a,b))/sigma)**2)
    return Delta

def read_data(filename):
    ''' import the data'''
    f = open('data/'+filename+'.csv','r')
    header = f.readline()
    x = f.readlines()
    f.close()
    return x


def fit(t,rho,sigma,theta_0):
    ''' Fit the model '''
    # in this example, we formulate the mixed model only for increasing timeseries in rho.
    # if the timeseries is decreasing we fit the transformed timeseries rho' = 1-rho (which is increasing)
    # this is checked by comparing the average of the first 10 versus the last 10 points in the timeseries
    bool_incr = True
    if np.sum(rho[:10])>np.sum(rho[-10:]):
        rho = 1.0 - rho
        bool_incr = False

    # fit the timeseries and extract the parameters
    fit_result = fit_mixed(t,rho,sigma,theta_0) 
    theta = fit_result[0] # estimated parameters [t0,y0,y1,a,b] from best fit
    chi2 = fit_result[1] # chi_squared of the best fit
    
    ## plot the time series and the data
    # timeseries for the mixed model with parameters of the best fit
    t_fit = np.linspace(t[0],t[-1],1000,endpoint=True)
    rho_fit = xcurve_5(t_fit,theta[0],theta[1],theta[2],theta[3],theta[4])
    
    # if the original timeseries rho is decreasing, transform back to rho = 1 - rho'
    if bool_incr == False:
        rho_fit = 1.0-rho_fit
        rho = 1.0-rho
    return [theta,t_fit,rho_fit]
    
