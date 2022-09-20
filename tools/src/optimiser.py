#!/usr/bin/python
# -*- coding: utf-8 -*-

"""INSTALL
Make sure this is in one of the sys.path directories:
> ! ln -vs $PWD/diffTimer.py ~/.ipython

USE

To load:
  import optimiser as o
To use:
  help(o.optimise)
To reload:
  reload(o)
"""

import sys
import subprocess as sp
from time import * #time,sleep
import errno

import numpy as np
m=np.math

#*******************************************************************************
# GLOBAL VARIABLES


#*******************************************************************************
def hyperbolic(x):
  """ Test function : minimum is at x=-(1./3)**.5 """
  return (x**2+1)**.5+.5*x


#*******************************************************************************
def exponential(x):
  """ Test function : minimum is at x=ln(1.5) """
  return np.exp(x)-1.5*x


#*******************************************************************************
def system(x):
  """Test function
  
  See optimise() for an example.
  """
  
  cmd=[system.cmd,str(x)]
  print( '+ '+' '.join(cmd) )
  p=sp.Popen(cmd, stdout=sp.PIPE)
  
  while True:
    l=p.stdout.readline()
    sys.stdout.write( l )
    if l=="":
      break
    lastL=l
  
  r=p.wait()
  
  if r!=0:
    return np.nan
  
  y=float(lastL)
  
  return y


#*******************************************************************************
def dichotomy(Xs,Ys):
  """ Search function  """
  
  i=np.argmin(Ys)
  
  dx0=dx1=-np.inf
  y0=y1=np.inf
  
  if i>0:
    dx0=Xs[i]-Xs[i-1]
    y0=Ys[i-1]
  
  if i<len(Xs)-1:
    dx1=Xs[i+1]-Xs[i]
    y1=Ys[i+1]
  
  if dx0>1.5*dx1:
    i-=1
  elif not dx0*1.5<dx1: # almost equal
    if y0<y1:
      i-=1
  
  x=(Xs[i]+Xs[i+1])*.5
  
  return x


#*******************************************************************************
def n_iteration(i,ix,Xs,Ys):
  """ Stop function  """
  
  if not hasattr(n_iteration,'max'):
    n_iteration.max=10
  
  return i>=n_iteration.max


#*******************************************************************************
def convergence(i,ix,Xs,Ys):
  """ Stop function  """
  
  if not hasattr(convergence,'max'):
    convergence.max=1e-9
  
  x=Xs[ix]
  dx=0
  if ix>=1:
    dxi=x-Xs[ix-1]
    dx=max(dx,dxi)
  if ix<i:
    dxi=Xs[ix+1]-x
    dx=max(dx,dxi)
  
  return dx<=convergence.max*abs(x)


#*******************************************************************************
def optimise(xmin,xmax,get_y,get_next_x,*stops):
  """
  Examples
  ========
  
  Xs,Ys=o.optimise(-1,1,o.hyperbolic,o.dichotomy,o.n_iteration)
  o.n_iteration.max=100
  Xs,Ys=o.optimise(-2,2,o.exponential,o.dichotomy,o.n_iteration)
  
  o.system.cmd='./testCd.sh'
  Xs,Ys=o.optimise(.001,.003,o.system,o.dichotomy,o.n_iteration)
  """
  
  Xs=[xmin,xmax]
  Ys=map(get_y,Xs)
  Is=[0,1]
  
  i=1
  
  stopC=len(stops)
  stop=False
  
  while not stop:
    x=get_next_x(Xs,Ys)
    if x in Xs:
      break
    y=get_y(x)
    
    if np.isnan(y):
      print( 'not coded to handle NAN yet' )
      sys.exit(errno.ENOEXEC)
    
    i+=1
    
    for ix in range(1,len(Xs)):
      if Xs[ix-1]<=x and x<=Xs[ix]:
        break
    Xs.insert(ix,x)
    Ys.insert(ix,y)
    Is.insert(ix,i)
    
    for stopi in range(stopC):
      if stops[stopi](i,ix,Xs,Ys):
        stop=True
        break
  
  Is=np.argsort(Is)
  Xs=np.array(Xs)[Is]
  Ys=np.array(Ys)[Is]
  
  return Xs,Ys


#*******************************************************************************
def main(args):
  """USE
  {0} script xmin xmax

DESCRIPTION
  Optimise the last value printed by script."""
  if len(args)!=4 or args[1]=='-h' or args[1]=='--help':
    print main.func_doc.format(args[0])
    return 0
  
#-------------------------------------------------------------------------------
  # parse arguments
  system.cmd=args[1]
  xmin=float(args[2])
  xmax=float(args[3])
  
#------------------------------------------------------------------------------
  # compute
  Xs,Ys=optimise(xmin,xmax,system,dichotomy,n_iteration)
  
  return 0


#*******************************************************************************
if __name__ == '__main__':
  sys.exit(main(sys.argv))
