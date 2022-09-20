#!/usr/bin/python
# -*- coding: utf-8 -*-

""" usage :
"""

import sys #.argv,.exit
import errno

from numpy import *
import scipy.io.netcdf as nc

#*******************************************************************************
def main(args):
  """USE
  {0} dim1 [ dim2 ... ] input.nc output.nc

DESCRIPTION
  Replace dimensions by time."""
  if len(args)<4 or args[1]=='-h' or args[1]=='--help':
    print main.func_doc.format(args[0])
    return 0
  
  rmds=args[1:-2]
  print repr(rmds)
  
  inF=nc.netcdf_file(args[-2])
  outF=nc.netcdf_file(args[-1],'w')
  
  inDs=inF.dimensions
  outF.createDimension('time',0)
  
  print repr(inDs)
  
  for d in inDs.keys():
    if rmds.count(d):
      continue
    outF.createDimension(d,inDs[d])
    outF.createVariable(d,'i',[d])[:]=range(inDs[d])
  
  for vn in inF.variables.keys():
    if rmds.count(vn):
      continue
    print 'processing '+vn+'...'
    inV=inF.variables[vn]
    inDs=inV.dimensions
    outDs=[]
    outNs=[]
    hasT=1
    for d in inDs:
      if not rmds.count(d):
        outDs.append(d)
        outNs.append(inF.dimensions[d])
      else:
        hasT*=inF.dimensions[d]
    if hasT>1:
      outDs.insert(0,'time')
      outNs.insert(0,hasT)
      print repr(inDs)+'->'+repr(outDs)
      print repr(outNs)
    
    outV=outF.createVariable(vn,inV.typecode(),outDs)
    outV[:]=inV[:].reshape(outNs)
  
  outF.close()
  inF.close()
  
  return errno.ENOEXEC

#*******************************************************************************
if __name__ == '__main__':
  #file:///usr/share/doc/packages/python/html/python-2.7-docs-html/library/sys.html#sys.exit
  sys.exit(main(sys.argv))
