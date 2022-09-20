#!/usr/bin/python
# -*- coding: utf-8 -*-

"""INSTALL
Make sure this is in one of the sys.path directories:
> ! ln -s $PWD/diffTimer.py ~/.ipython

USE
To load:
  import diffTimer as dt
To use:
  help(dt)
To reload:
  reload(dt)
"""

import matplotlib
matplotlib.use('TkAgg')
np=matplotlib.numpy
import matplotlib.pyplot as plt

#*******************************************************************************
# DATA
unitPrefix=''

# KEEP THIS CONSTANT
UNIT_TRANSLATION={
  '':1,
  'm':.001,
  'µ':1e-6
  }

#*******************************************************************************
def getTimerStats(path,*args):
  """
  Use
  ---
    dt.getTimerStats(path,substr1[,jid1][,substr2[,jid2][,...]])
  
  Parameters
  ----------
    path
    substr : substring
    jid : proc job id, default: 0
  
  Returns
  -------
    a dict of np.array
  """
  
#-------------------------------------------------------------------------------
  # parse arguments
  
  subs=[]
  i=0
  n=len(args)
  
  while i<n:
    # CODED FROM difftimer_t::destroy()
    sub=args[i]
    
    if i+1<n and type(args[i+1])==int :
      i+=1
      jid=args[i]
    else:
      jid=0
    
    sub+=':difftimer_t(,'+str(jid)+') '
    
    subs+=[sub]
    
    i+=1
  
#-------------------------------------------------------------------------------
  # initialisation
  
  n=len(subs)
  
  outputs={}
  
  for i in range(n):
    outputs[subs[i]]=[]
  
#-------------------------------------------------------------------------------
  # parse input file
  
  f=open(path)
  
  while True:
    
    # read line
    l=f.readline()
    if l=='':
      break
    
    # find match
    for i in range(n):
      sub=subs[i]
      pos=l.find(sub)
      if pos>=0:
        break
    
    if pos<0:
      continue
    
    # extract numbers
    pos=l.find('s = ',pos)
    
    unitPrefix=l[pos-1]
    if unitPrefix==' ':
      unitPrefix=""
    
    row=np.fromstring(l[pos+4:],sep=' ')*UNIT_TRANSLATION[unitPrefix]
    
    outputs[sub]+=[row]
  
  f.close()
  
#-------------------------------------------------------------------------------
  for sub in subs:
    
    output=outputs[sub]
    
    output=np.array(output)
    
    try:
      nj,ni=np.shape(output)
    except:
      print sub+': error with array of shape '+str(np.shape(output))
      continue
    
    print sub+': read {} lines and {} columns'.format(nj,ni)
    outputs[sub]=output
  
  return outputs


#*******************************************************************************
def plotTimerStats(outputs):
  """ plot averages of durations as given by dt.getTimerStats()
  
  Parameter
  ---------
    outputs : from dt.getTimerStats()
  """
  
  for sub,output in outputs.iteritems():
    
    As=np.average(output,axis=0)
    
    N=len(As)
    sumAs=sum(As)
    percents=As*(100./sum(As))
    
    print sub,' '.join(map(lambda x:'{:.4}'.format(x),As))
    
    plt.close(sub)
    F=plt.figure(sub)
    A=F.gca()
    
    if unitPrefix!='':
      As*=1./UNIT_TRANSLATION[unitPrefix]
    
    xs=np.arange(N)
    bars=A.bar(xs,As,align='center')
    
    A.set_title('sum={:.3}'.format(sumAs)+unitPrefix+'s')
    
    A.set_xlim(-.5,N-.5)
    xls=map(lambda x:'{:.3}%'.format(x),percents)
    A.set_xticks(xs)
    A.set_xticklabels(xls)
    
    A.set_ylabel(unitPrefix+'s')
  
  return
