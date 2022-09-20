#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys #.argv,.exit,.std*
import errno
from pylab import *#plot. See http://matplotlib.sourceforge.net/gallery.html

dashLine='-'*80+'\n'
starLine='*'*80

#*******************************************************************************
def reim(m):
  """separates real and imaginary parts
  
  Parameters:
  m : complex matrix
  
  Returns a matrix of column-by-column real and imaginary parts.
  """
  c=size(m,1)
  ri=c_[real(m),imag(m)]
  realCs=arange(c)*2
  return ri[:,r_[realCs,realCs+1]]

#*******************************************************************************
def aliasing(n,s,f,a):
  """test one harmonic analysis
  
  Parameters:
  n : nummber of samples
  s : sampling frequency
  f : signal frequencies
  a : complex amplitudes
  """
  print starLine
  
  print "{} samples at {}".format(n,s)
  t=matrix(arange(0,n)*(1./s)).T
  print 't.T=',t.T
  
  f=array(f)
  print "f=",f
  w=f*(2.*pi)
  a=matrix(a).T
  print "a.T=",a.T
  
  M=exp(1j*t*w)
  print "M=",M
  
  Ma=M*a
  print "imag(M*a).T=",imag(Ma).T
  series=real(Ma)
  print "real(M*a).T=series.T=",series.T
  
#-------------------------------------------------------------------------------
  # real+imaginary solution
  print dashLine+'REAL+IMAGINARY'
  M2=reim(M)
  print "M2=",M2
  
  A=M2.T*M2
  dA=det(A)
  print "A="
  print A,dA
  
  b=M2.T*series
  print "b=",b
  
  pA=pinv(A)
  print 'pinv(A)=',pA
  print 'solution:',(pA*b).T
  return
  
#-------------------------------------------------------------------------------
  # complex solution
  print dashLine+'COMPLEX'
  A=M.H*M
  dA=det(A)
  print "A="
  print A,dA
  
  b=M.H*series
  print "b=",b
  
  pA=pinv(A)
  print 'pinv(A)=',pA
  print 'solution:',(pA*b).T

#*******************************************************************************
def main(args):
  """USE
  {0}

DESCRIPTION
  Tests some harmonic analyses."""
  if len(args)>1:
    print main.func_doc.format(args[0])
    return 0
  
  set_printoptions(precision=4,suppress=True,linewidth=80)
  
  # OK
  aliasing(8,8,[1,2],[1,.5])
  # Still OK
  aliasing(4,4,1,1+1j)
  # self-aliased
  aliasing(2,2,1,1+1j)
  # aliased
  aliasing(9,3,[4,8],[1,1+1j])
  # about OK
  aliasing(9,3,[4,8.25],[1,1+1j])
  
  return 0

#*******************************************************************************
if __name__ == '__main__':
  #file:///usr/share/doc/packages/python/html/python-2.7-docs-html/library/sys.html#sys.exit
  sys.exit(main(sys.argv))
