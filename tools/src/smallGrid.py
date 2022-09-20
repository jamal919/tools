#!/usr/bin/python
# -*- coding: utf-8 -*-

""" usage :
"""

import sys #.argv,.exit
import errno

from numpy import *
import scipy.io.netcdf as nc

class ncvar:
  def __init__(self,n,l):
    """ creates the object """
    self.n = n
    self.l = l
  def __getitem__(self, key):
    """ x.__getitem__(key) <==> x[key] """
    if type(key)==slice:
      r=[]
      if key.step==None:
        s=1
      else:
        s=key.step
      for i in range(key.start,min(key.stop,2),s):
        r+=[self[i]]
      return r
    if key==0:
      return self.n
    if key==1:
      return self.l

def main(args):
  """USE
  {0} FILE

DESCRIPTION
  makes a small grid and saves it into FILE."""
  if len(args)<2 or args[1]=='-h' or args[1]=='--help':
    print main.func_doc.format(args[0])
    return 0
  
  f=nc.netcdf_file(args[1],'w')
  
  doDV=True
  doContent=not doDV
  
  lad=ncvar('lad',5)
  lod=ncvar('lod',5)
  f.createDimension(*lad[:])
  f.createDimension(*lod[:])
  lavd=ncvar('lavd',4)
  loud=ncvar('loud',4)
  f.createDimension(*lavd[:])
  f.createDimension(*loud[:])
  
  if doDV:
    f.createVariable(lad.n,'f',[lad.n]).axis='Y'
    f.createVariable(lod.n,'f',[lod.n]).axis='X'
    f.createVariable(lavd.n,'f',[lavd.n]).axis='Y'
    f.createVariable(loud.n,'f',[loud.n]).axis='X'
  
  lav=f.createVariable('lav','f',['lad','lod']);lav.standard_name='latitude'
  lov=f.createVariable('lov','f',['lad','lod']);lov.standard_name='longitude'
  lauv=f.createVariable('lauv','f',['lad','loud']);lauv.standard_name='latitude_u'
  louv=f.createVariable('louv','f',['lad','loud']);louv.standard_name='longitude_u'
  lavv=f.createVariable('lavv','f',['lavd','lod']);lavv.standard_name='latitude_v'
  lovv=f.createVariable('lovv','f',['lavd','lod']);lovv.standard_name='longitude_v'
  x=f.createVariable('x','f',['lad','lod'])
  u=f.createVariable('u','f',['lad','loud'])
  v=f.createVariable('v','f',['lavd','lod'])
  
  if doContent:
    x.content="YX"
    u.content="YX"
    v.content="YX"
  
  lov[:]=lad.l*[range(lod.l)]
  lav[:]=map(lambda x: [x]*lod.l,range(lad.l))
  x[:]=array(range(lod.l*lad.l)).reshape(x.shape)
  louv[:]=lad.l*[arange(loud.l)+.5]
  lauv[:]=map(lambda x: [x]*loud.l,range(lad.l))
  u[:]=array(range(loud.l*lad.l)).reshape(u.shape)
  lovv[:]=lavd.l*[range(lod.l)]
  lavv[:]=map(lambda x: [x]*lod.l,arange(lavd.l)+.5)
  v[:]=array(range(lod.l*lavd.l)).reshape(v.shape)
  
  mask=float(999.)
  x._FillValue=mask
  x[:][2][2]=mask
  x.scale_factor=.1
  
  f.close()
  
  return errno.ENOEXEC
if __name__ == '__main__':
  #file:///usr/share/doc/packages/python/html/python-2.7-docs-html/library/sys.html#sys.exit
  sys.exit(main(sys.argv))
