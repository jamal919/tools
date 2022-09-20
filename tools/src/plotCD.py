#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys #.argv,.exit,.std*
import errno
from pylab import *#plot. See http://matplotlib.sourceforge.net/gallery.html

#*******************************************************************************
def verboseSaveFig(fig,path):
  sys.stdout.write('saving to '+path+'...')
  fig.savefig(path,bbox_inches='tight')
  sys.stdout.write(' done.\n')

#*******************************************************************************
def toFloat(s):
  # because fromstring() is just buggy
  try:
    return float(s)
  except ValueError:
    return nan

#*******************************************************************************
def main(args):
  """USE
  {0} OPTIONS

DESCRIPTION
  Plot cumulative distribution.

OPTIONS
  -i  followed by path to input file
  -t  followed by title. Default: input file path with extension removed.
  -o  followed by path to output file. Default: input file path with extension replaced with png
  -x  followed by [Xmin][,[Xmax][,[Xincrement]]]
  -y  followed by Yincrement

EXAMPLE
  altimetry-detidor ... |sed -re 's/[^ ]+ +([^ ]+) correlated with K1/\1/g;t;d' >correlation.log
  plotCD.py -i correlation.log -x .1,1.7,.1 -y 10"""
  n=len(args)
  
  if n<2 or args[1]=='-h' or args[1]=='--help':
    print main.func_doc.format(args[0])
    return 0
  
  titleStr=''
  figP=''
  xO=',,'
  yO=-1
  
  i=1
  while i<n:
    opt=args[i]
    arg=args[i+1]
    i+=2
    if opt=='-i': dataP=arg;continue
    if opt=='-t': titleStr=arg;continue
    if opt=='-o': figP=arg;continue
    if opt=='-x': xO=map(toFloat,arg.split(','));continue
    if opt=='-y': yO=float(arg);continue
    print "*** "+opt+" option not recognised ***"
    print main.func_doc.format(args[0])
    return -1
  
  sys.stdout.write('loading '+dataP+'...')  
  data=loadtxt(dataP)
  sys.stdout.write(' done.\n')
  
  rootP=dataP[:dataP.rfind('.')]
  if not titleStr:
    titleStr=rootP[(rootP.rfind('/')+1):]
    print "title default to: "+titleStr
  if not figP:
    figP=rootP+'.png'
    print "figure path default to: "+figP
  
  n=len(data)
  data=sort(data)
  
  F=figure(1)
  A=F.gca()
  A.plot(data,arange(0,n)*(100./n))
  
  if titleStr:
    A.set_title(titleStr)
  
  xO=xO+[nan]*(3-len(xO))
  print '-x:',xO
  if not isnan(xO[0]):
    A.set_xlim(left=xO[0])
  if not isnan(xO[1]):
    A.set_xlim(right=xO[1])
  if not isnan(xO[2]):
    A.xaxis.set_major_locator(MultipleLocator(xO[2]))
    A.grid(True,axis='x')
  
  A.yaxis.set_major_formatter(FormatStrFormatter('%g%%'))
  if yO>0:
    A.yaxis.set_major_locator(MultipleLocator(yO))
    A.grid(True,axis='y')
  
  verboseSaveFig(F,figP)

#*******************************************************************************
if __name__ == '__main__':
  #file:///usr/share/doc/packages/python/html/python-2.7-docs-html/library/sys.html#sys.exit
  sys.exit(main(sys.argv))
