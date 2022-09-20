#!/usr/bin/python
# -*- coding: utf-8 -*-

""" usage :
plotLogo.py /mnt/srv-ecola/models/global-MR/data/global.nei [name]
"""

import sys #.argv,.exit
from scipy.interpolate import interp1d
from pylab import *#plot. See http://matplotlib.sourceforge.net/gallery.html
import os #os.readable

d2r=pi/180.
#approximately Florent's office
la0,lo0=array([46.75984,1.738281])*d2r
#better look ?
la0,lo0=array([0,-30])*d2r
cla0=cos(la0)
sla0=sin(la0)
mla0=mat([[cla0,-sla0],[sla0,cla0]])

#*******************************************************************************
def float2str(x):
  return '{:g}'.format(x)

#*******************************************************************************
def main(args):
  """USE
  {0} /mnt/srv-ecola/models/global-MR/data/global.nei [name]

DESCRIPTION
  Plot name's logo.
  If name (case insensitive) is not given, plot all logos : Tools,Tugo,xscan,POCViP."""
  if len(args)<2 or args[1]=='-h' or args[1]=='--help':
    print main.func_doc.format(args[0])
    return 0
  
  gF=open(args[1])#grid file
  pC=int(gF.readline())#point count
  nC=int(gF.readline())#neighbour count
  bounds=array(gF.readline().split(),float)
  xyzc=zeros((pC,5))+nan
  figure(figsize=(2.2,2.2))
  hold(True)
  palette=interp1d([0,.4,1],array([[0,1,0],[.5,.5,1],[.125]*3]).transpose())
  for pI in range(pC):
    l=array(gF.readline().split(),float)
    lo,la=l[[1,2]]*d2r-[lo0,0]
    x=cos(la)*cos(lo)
    y=cos(la)*sin(lo)
    z=sin(la)
    x,z=array(mat([x,z])*mla0)[0]
    xyzc[pI]=array([x,y,z,l[4],pI])
    for nI in l[arange(nC)+5]:#for all neigbour indexes
      nI=nI-1 #convert to C-like index
      if 0<=nI and nI<pI and (x>0 or xyzc[nI][0]>0):
        y,z,c=xyzc[[pI,nI]][:,[1,2,3]].transpose()
        plot(y,z,color=palette(min(mean(c)/8e3,1)),linewidth=.3)
    sys.stderr.write(str(pI)+'/'+str(pC)+'\r')
  print
  gF.close()
  x=arange(0,1.0005,.001)*2*pi
  plot(cos(x),sin(x),color='black',linewidth=2)
  A=gca()
  A.set_xlim(-1.1,1.1)
  A.set_ylim(-1.1,1.1)
  A.set_axis_off()
  A.set_position([0,0,1,1])
  root='logo'
  print root
  savefig(root+'.eps')
  savefig(root+'.png',dpi=254)
  for label in (['Tools ',45],['TUGO',40],['xscan',40],['POCViP',28]):
    if len(args)>2 and args[2].lower()!=label[0].lower():
      continue
    root=label[0].replace(' ','')+'-logo'
    print root
    name=text(0,0,label[0],horizontalalignment='center',verticalalignment='center',size=label[1],weight='bold')
    savefig(root+'.eps')
    savefig(root+'.png',dpi=254)
    name.remove()
  return 0

#*******************************************************************************
if __name__ == '__main__':
  #file:///usr/share/doc/packages/python/html/python-2.7-docs-html/library/sys.html#sys.exit
  sys.exit(main(sys.argv))
