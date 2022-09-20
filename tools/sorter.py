#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys #.argv,.exit,.std*
import subprocess as p #Popen,PIPE,call
from datetime import datetime
import struct
import re
import errno

import matplotlib
# TODO: test all interactive backends...
matplotlib.use('TkAgg')
np=matplotlib.numpy
import matplotlib.pyplot as plt #plot. See http://matplotlib.sourceforge.net/gallery.html


#*******************************************************************************
def read1HgResponse(server):
  
  channel, length = struct.unpack('>cI', server.stdout.read(5))
  if channel not in 'oer':
    print "uncoded channel:",channel,' length:',length
    return channel,length
  
  val=server.stdout.read(length)
  
  return channel,val


#*******************************************************************************
def startHgServer():
  
  print 'starting server'
  server = p.Popen(['hg','--config','ui.interactive=True','serve','--cmdserver','pipe'],stdin=p.PIPE,stdout=p.PIPE)
  
  channel,output=read1HgResponse(server)
  print "server said:<",output,">"
  
  return server


#*******************************************************************************
def runHgCommand(server,args,verbose=False):
  
  if verbose:
    print 'sending',repr(args),'to server'
  data='\0'.join(args)
  server.stdin.write('runcommand\n' + struct.pack('>I', len(data)) + data)
  server.stdin.flush()
  
  output=''
  
  while True:
    channel, val = read1HgResponse(server)
    if channel == 'o':
      output+=val
    elif channel == 'e':
      print "error:", repr(val)
      break
    elif channel == 'r':
      if verbose:
        print "exit code:", struct.unpack(">l", val)[0]
      break
  
  return output


#*******************************************************************************
def countOfMains():
  
  server = startHgServer()
  
  output=runHgCommand(server,['grep','-I','*.cpp','-I','src/*.cpp','--all','^ *int +main'],True)
  
  rex=re.compile(':([0-9]+):([+-]):')
  diffdict={}
  for l in output.split('\n'):
    #print l
    if l=="":
      continue
    rev,inc=rex.search(l).groups()
    rev=int(rev)
    if not diffdict.has_key(rev):
      diffdict[rev]=0
    if inc=='+':
      diffdict[rev]+=1
    if inc=='-':
      diffdict[rev]-=1
  
  print repr(diffdict)
  
  output=runHgCommand(server,['log','-r','tip','--template',"{rev}"])
  last=int(output)+1
  diffs=np.zeros(last,int)
  for rev in diffdict.keys():
    diffs[rev]=diffdict[rev]
  
  counts=np.cumsum(diffs)
  print repr(counts)
  
  t=np.zeros(last)
  dateformat='%Y-%m-%d %H:%M:%S'
  for rev in range(last):
    output=runHgCommand(server,['log','-r',str(rev),'--template',"{date(date, '"+dateformat+"')}"])
    t[rev]=matplotlib.dates.date2num(datetime.strptime(output,dateformat))
  
  return t,counts


#*******************************************************************************
def verbose_savefig(F,path,*args,**kwargs):
  print 'saving to '+path
  F.savefig(path,*args,**kwargs)


#*******************************************************************************
def plotCountOfMains():
  t,n=countOfMains()
  
  F=plt.figure(1)
  A=F.gca()
  A.plot_date(t,n,'-')
  A.set_title('number of executables')
  verbose_savefig(F,'ToolsMains.png')


#*******************************************************************************
def main(args):
  """USE
  {0}

DESCRIPTION
  Tools sorting helper.
  
  For the moment, it only produces a figure with the number of main vs. commit date. """
  if len(args)>1:
    print main.func_doc.format(args[0])
    return 0
  
  plotCountOfMains()
  
  return 0

#*******************************************************************************
if __name__ == '__main__':
  #file:///usr/share/doc/packages/python/html/python-2.7-docs-html/library/sys.html#sys.exit
  sys.exit(main(sys.argv))
