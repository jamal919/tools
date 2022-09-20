#!/usr/bin/python

""" Cleans-up redundant lines in cmdline.* files
usage :
cleanCmdline.py cmdline.*
"""

from glob import * #glob
import sys #.argv,.exit
import os #.environ

#*******************************************************************************
# which
pathDs=os.environ['PATH'].split(':')
def concat(x,y):
  return x+y
def isExecutable(x):
  return os.access(x,os.X_OK)
def which(x):
  return filter(isExecutable,map(concat,pathDs,['/'+x]*len(pathDs)))[0]

#*******************************************************************************
#file:///usr/share/doc/packages/python/html/python-2.7-docs-html/faq/programming.html#how-do-i-find-the-current-module-name
def isNotEmpty(x):
  return x!=''

#*******************************************************************************
def main(args):
  """USE
  {0} [OPTIONS] file1 [file2 ...] [OPTIONS file3 [file4 ...] ]

DESCRIPTION
  Removes redundant entries in the given cmdline files.
  Any non-option argument are considered to be files.
  Directories are replaced by all their cmdline.* files.

OPTIONS
  -h,--help  print this help and exit. When used after any other argument, this will not be considered an option!
  -w,--which  expand programme name as with the which command.
  -b,--basename  reduce programme paths as with the basename command (default).
  -x  execute last command of next file
  """
  if len(args)<2 or args[1]=='-h' or args[1]=='--help':
    print main.func_doc.format(args[0])
    return 0
  
  # replace directories
  iA=1
  while iA<len(args):
    arg=args[iA]
    if not os.path.isdir(arg):
      iA+=1
      continue
    files=glob(arg+'/cmdline.*')
    args=args[:iA]+files+args[1+iA:]
    iA+=len(files)
  
  expand=False
  executeNext=False
  
  for iA in range(1,len(args)):
    try:
      arg=args[iA]
      
      if arg=='-w' or arg=='--which':
        expand=True
        continue
      if arg=='-b' or arg=='--basename':
        expand=False
        continue
      
      if arg=='-x':
        executeNext=True
        continue
      
      # read file
      f=open(arg)
      ls=filter(isNotEmpty,f.read().splitlines())
      f.close()
      
      # pre-process programme paths
      for i1 in range(1,len(ls),2):
        cmd=ls[i1].split(' ',1)
        firstSep=cmd[0].rfind('/')
        if firstSep<0 and expand:
          ls[i1]=which(cmd[0])+' '+cmd[1]
        if firstSep>=0 and not expand:
          ls[i1]=ls[i1][firstSep+1:]
      
      # write file, removing redundancies
      f=open(arg,'w')
      i2=1
      while i2<len(ls):
        i1=i2
        i2=i1+2
        while i2<len(ls) and ls[i1]==ls[i2]:
          i2+=2
        f.write("\n"+",".join(ls[i1-1:i2-2:2])+"\n"+ls[i1]+"\n")
      f.close()
      
      print arg+":\n"+ls[i1]
      
      if executeNext:
        os.system(ls[i1])
        executeNext=False
    except Exception as err:
      print 'error with '+arg
  return 0

#*******************************************************************************
if __name__ == '__main__':
  #file:///usr/share/doc/packages/python/html/python-2.7-docs-html/library/sys.html#sys.exit
  sys.exit(main(sys.argv))
