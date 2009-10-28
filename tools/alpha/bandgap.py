#!/usr/bin/env python

import sys

FILE = open(sys.argv[1],'r')

occupied = []
empty = []
while FILE:
  line = FILE.readline()
  if line=="":
    break
  if line.startswith(" E-fermi"):
    FILE.readline()
    FILE.readline()    
    
    line = FILE.readline()
    while line.startswith(" k-point"):
      occupied.append([])
      empty.append([])
      
      FILE.readline()
      line = FILE.readline()
      data = line.split()
      while len(data) == 3:
        occ = float(data[2])
        if occ > 0.5:
          occupied[-1].append(float(data[1]))
        else:
          empty[-1].append(float(data[1]))
          
        if occ < 0.99 and occ > 0.01:
          print "warning: possibly metallic!" 
                  
        line = FILE.readline()
        data = line.split()
        
      line = FILE.readline()
  
max_occ = [ max(eig) for eig in occupied ]
min_emp = [ min(eig) for eig in empty    ]

print min(min_emp)-max(max_occ)

print len(occupied)
print len(empty)
#print ""
#print occupied
#print ""
#print empty

