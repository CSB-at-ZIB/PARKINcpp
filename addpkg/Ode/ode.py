#!/usr/local/bin/python
# KMB 2001 Aug 09,20
# TODO: -r and -a error specification
# Ode friendly interface
# Usage (arguments may be in any order): 
#   ode.py var1.=rhs1 var2.=rhs2 ... var1=iv1 var2=iv2 ... tstart:[number of steps]:tend plotvar1,plotvar2,...
# ("var'=rhs" is also accepted for an ode specification)
#
# e.g. 
#   ode.py x.=-y y.=x x=1 y=0 y 0::10 | p
#   or, instead of:
#     Ode "3*(b-a)" "-a*c+26*a-b" "a*b-c" 0 1 0 0 30 100 | cols -c1,3 | p
#   use:
#     ode.py x.=3*y-3*x y.=-x*z+26*x-y z.=x*y-z x=0 y=1 z=10 x 0::30 | p
#   or:
#     ode.py "x'=3*(y-x)" "y'=-x*z+26*x-y" "z'=x*y-z" x=0 y=1 z=10 x,y,z 0:200:100 
# Backup:
#   mcopy -o ode.py a:.
#   scp ode.py wulf1:.

import re
from sys import argv,stderr,exit
from string import replace,split
from os import system
from commands import getstatusoutput

dbg=0
verbose=0
symbols={}
iv={}
rhs={}
iv={}
ns=0
dn=0
up=1
sp=100
pv=[]
plottime=1
name=     re.compile(r'\b[a-zA-Z]\w*\b(?!\()')
equals=   re.compile(r'=')
dotequals=re.compile(r"[\|.']=")
dollar=   re.compile(r'\$(\d+)')
limits=   re.compile(r'(-?\d+(?:\.\d+)*):(\d*):(-?\d+(?:\.\d+)*)')
plotvars= re.compile(r'\$(\d+)(?:,\$(\d+))*')

for arg in argv[1:]:
  if dbg: print >>stderr,'\narg:',arg
  while 1:
    n=name.search(arg)
    if not n: break
    nm=n.group()
    if dbg: print >>stderr,'found name nm=',nm
    if not nm in symbols.keys():
      ns+=1
      symbols[nm]='$'+str(ns)
    repl=symbols[nm]
    arg=re.compile(r'\b'+nm+r'\b').sub(repl,arg)
    if dbg: print >>stderr,'now arg=',arg
  r=limits.match(arg)
  if r: 
    if dbg: print >>stderr,'got a range specification'
    dn=r.groups()[0]
    sp=r.groups()[1]
    up=r.groups()[2]
    if sp=='': sp='100'
    if dbg: print >>stderr,'dn=',dn,'step=',sp,'up=',up
  else:
    r=dotequals.split(arg)
    if dbg: print >>stderr,'r=',r
    if len(r)>1: 
      if dbg: print >>stderr,'got a d.e.',r
      rhs[r[0]]=r[1]
    else:
      r=equals.split(arg)
      if len(r)>1: 
        if dbg: print >>stderr,'got an i.v.',r
	try:
	  rh=float(eval(r[1]))
        except:
	  print >>stderr,'rhs of initial value equation is invalid'
	  exit(1)
        iv[r[0]]=rh
      else:
        r=split(arg,',')
        if r:  
          if dbg: print >>stderr,'got a plotvars list=',r,
          for p in r: 
            if p!='0': pv.append(p[1:])
            else: plottime=0
          if dbg: print >>stderr,'pv=',pv
        else:
          print >>stderr,'syntax error'
	  exit(1)

if verbose: 
  print >>stderr,'symbols:',symbols
  print >>stderr,'rhs:',rhs
  print >>stderr,'iv:',iv
  print >>stderr,'plotvars:',pv

if len(symbols)<len(rhs):
  print >>stderr,'Too few symbols'
  exit(1)
if len(symbols)<1:
  print >>stderr,'Too few symbols'
  exit(1)
if len(symbols)>20:
  print >>stderr,'Too many symbols'
  exit(1)
if len(symbols)>len(rhs):
  print >>stderr,'Too many symbols:',len(symbols)
  exit(1)
if len(iv)<len(rhs):
  print >>stderr,'Not enough initial values'
  exit(1)
if len(iv)>len(rhs):
  print >>stderr,'Not enough differential equations'
  exit(1)

def dollar2name(x):
  return chr(96+int((x.group()[1:])))

ode='./Ode -r1e-9 -a1e-9 '
kys=rhs.keys()
kys.sort()
for k in kys: # odes
  ode+='"'+dollar.sub(dollar2name,rhs[k])+'" '
for k in kys: # ivs
  ode+=str(iv[k])+' '
ode+=str(dn)+' '+str(up)+' '+str(sp)#+' 2>/dev/null'
if verbose: 
  print >>stderr,'Ode command line: '
  print >>stderr,ode
if dbg: exit()
(status,output)=getstatusoutput(ode)
lines=split(output,'\n')
if not pv: pv=range(1,ns+1)
pv=map(int,pv)
for line in lines:
  if line[0]=='|': # Ode stderr output
    print >>stderr,line
  else:
    cols=split(line)
    if plottime: print cols[0], # 'time'
    for i in pv: 
      if float(cols[i])>=0: print '', # align columns
      print cols[i],
    print
