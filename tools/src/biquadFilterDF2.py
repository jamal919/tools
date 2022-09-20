#!/usr/bin/python

#http://en.wikipedia.org/wiki/Digital_biquad_filter#Direct_Form_2
#from http://www.musicdsp.org/files/Audio-EQ-Cookbook.txt the transfer function may be LIKE this :
#            b0 + b1*z^-1 + b2*z^-2
#    H(z) = ------------------------
#            a0 + a1*z^-1 + a2*z^-2


""" usage :
"""

import sys #.argv,.exit
from numpy import *
from pylab import *

import scipy.io.netcdf as nc

#*******************************************************************************
def test(args):
  d2s=86400.
  h2s=3600.
  s=1/d2s #sampling frequency
  t0=62.*d2s
  Q=t0/(4.*d2s)
  w0=2*pi/t0 # central frequency in rad/s
  b0=w0/(2*Q*s)
  b1=0;b2=-b0;a1=2*(1-b0)*cos(w0/s);a2=2*b0-1; #Band pass
  print str([t0/h2s,s*h2s,Q,Q/w0*7.*s])
  print str([b0,b1,b2,a1,a2])

  dataDir='/data1/GLORYS/v2'
  path=dataDir+'/sossheig_reg_Glorys2V1_199206-200912.nc'
  path=''
  
  if path:
    print 'reading '+path
    f=nc.netcdf_file(path)
    latName='latitude'
    lonName='longitude'
    lat=f.variables[latName][:]
    lon=f.variables[lonName][:]
    v=f.variables['sossheig']
    t=f.variables['time']
    O=B0=B1=B2=v[0]*0.
    iMax=v.shape[0]
    dPath=dataDir+'/filtering.dat'
  else:
    iMax=6421
    t=range(iMax)
    dPath='test.dat'
    # PLOT WITH: x,y=loadtxt('test.dat',unpack=True)
  
  t0=t[0]
  
  dF=open(dPath,'w')
  O0=B00=B01=B02=0.
  
  for i in range(iMax):
    sys.stderr.write(str(i)+'/'+str(iMax)+'\r')
    
    if path:
      # delay registers
      B2=B1;B1=B0
      B0=B1*a1+B2*a2+v[i]
      O+=(B0*b0+B1*b1+B2*b2)**2
    
    # Dirac impulse response
    B02=B01;B01=B00
    B00=B01*a1+B02*a2+(i==0)
    O0i=B00*b0+B01*b1+B02*b2
    dF.write('{:g} {:g}\n'.format(t[i]-t0,O0i))
    O0+=O0i**2
  
  print str(O0) #=0.00326908704848
  
  if not path:
    return 0
  
  f.close()
  
  path='/data1/FES2012/GLORYS/filtered_reg_Glorys2V1_199206-200912.nc'
  print '\nwriting output to '+path
  f=nc.netcdf_file(path,'w')
  f.createDimension(latName,O.shape[0])
  f.createDimension(lonName,O.shape[1])
  f.createVariable(latName,'f',[latName])[:]=lat
  f.createVariable(lonName,'f',[lonName])[:]=lon
  f.createVariable('rootMeanSquare','f',[latName,lonName])[:]=sqrt(O/iMax)
  f.createVariable('normalisedRootMeanSquare','f',[latName,lonName])[:]=sqrt(O/O0)*2/iMax #=sqrt(O/iMax/(O0*iMax/4))
  f.close()
  
  return 0

#*******************************************************************************
def main(args):
  a=array([.5,.0]);b=array([1.,0.,0.])#Low pass
  x=array([1.,0.,0.])#Dirac impulse
  n=128
  y=empty(n)
  for i in range(n):
    y[i]=sum(x*b)
    x=x[:2]
    x=hstack((sum(x*a),x))
  plot(y)
  show()
  plot(abs(fft(y)))
  show()
  return 0

#*******************************************************************************
if __name__ == '__main__':
  #file:///usr/share/doc/packages/python/html/python-2.7-docs-html/library/sys.html#sys.exit
  sys.exit(test(sys.argv))
