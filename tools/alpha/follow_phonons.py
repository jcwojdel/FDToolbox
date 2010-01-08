#!/usr/bin/env python

import sys
import os

sys.path.append(os.path.join(sys.path[0],'../../'))
from fdtoolbox.calculation_set import *
from fdtoolbox.utility import *
from fdtoolbox.linear_expansion import *
from numpy import *

def align(preold_evectors, old_evectors, evalues, evectors):
  if old_evectors is None:
    return

  compatibility_map = abs(dot(old_evectors.T, evectors))
  if preold_evectors is not None:
    compatibility_map += 0.1*abs(dot(preold_evectors.T, evectors))
  new_evalues = []
  new_evectors = []
  for i in range(len(old_evalues)):
    current_row = compatibility_map[i,:]
    index = int( nonzero(current_row == current_row.max())[1] )
    new_evalues.append( evalues[index] )
    new_evectors.append( evectors[:,index] )

    compatibility_map[:,index] = 0.

  evalues[:] = array( new_evalues )
  evectors[:] = array( new_evectors ).T
  
def interpolate(v_list, num, spl = False):
  if num > 1:
    v_list = interpolate(v_list, num-1, spl)

  new_list = []
  for i in range(len(v_list)-1):
    new_list.append(v_list[i])
    if spl is False:
      new_list.append((v_list[i]+v_list[i+1])/2)
    else:
      mid=(v_list[i]+v_list[i+1])/2
      if i == 0:
        mid2=(3*v_list[i]+v_list[i+2])/4
      elif i == len(v_list)-2:
        mid2=(v_list[i-1]+3*v_list[i+1])/4
      else:
        mid2=(v_list[i-1]+v_list[i+2])/2
      new_list.append(mid+(mid-mid2)/6)
  
  new_list.append(v_list[-1])
      
  
  return new_list

preold_evalues = None
preold_evectors = None
old_evalues = None
old_evectors = None
all_evalues = []
all_evectors = []
saxis = [0., 0., 1.]

listofdirs=['thick-0.06', 'thick-0.055', 'thick-0.05', 'thick-0.04', 'thick-0.03', 'thick-0.02', 'thick-0.01', 'thick_limit']
listofdirs.reverse()
#listofdirs=['thick-0.07', 'thick-0.06', 'thick-0.05', 'thick-0.04', 'thick-0.03']
#listofdirs.reverse()

for directory in listofdirs:
  cs=calculation_set.read_directory(directory, saxis, True)

  cs.set_ionic('.*/calc_.._.*[0123]')
  cs.set_lattice('.*/calc_00_...')
  cs.set_groundstate('.*/calc_00_0')
 
  cs.try_fix_displacements()

  lin_exp = linear_expansion( cs )
  lin_exp.calculate_expansion_coefficients()

  fc_matrix, units = lin_exp.force_constant_matrix()
  evalues, evectors = linalg.eig(fc_matrix)
  
  #print mat2str(evectors)
  
  align(preold_evectors, old_evectors, evalues, evectors)

  all_evalues.append( evalues )
  all_evectors.append( evectors )

  preold_evalues = old_evalues
  preold_evectors = old_evectors
  
  old_evalues = evalues
  old_evectors = evectors
  
  
  
  
all_colors = []
for vec in all_evectors:
  vsqr = array(multiply(vec,vec))
  R = sum(vsqr[0:6],0)
  G = sum(vsqr[6:12],0)
  B = sum(vsqr[12:30],0)/3
  RGB=array([sqrt(R),sqrt(G),sqrt(B)])
  #RGB = RGB / sum(RGB)
  all_colors.append(RGB)

all_evalues = interpolate(all_evalues, 4, True)
all_colors = interpolate(all_colors, 4)

for i in range(len(all_colors)):
  c = all_colors[i]
  c = 65536*(255*c[0]).round()+256*(255*c[1]).round()+(255*c[2]).round()
  all_colors[i] = c

#print mat2str( array(all_evalues) )
#print mat2str( array(all_colors), "%d" )

all_evalues = mat(all_evalues)
all_colors = mat(all_colors)

for i in range(all_evalues.shape[1]):
  print mat2str(hstack([ array([range(all_evalues.shape[0])]).T , all_evalues[:,i], all_colors[:,i]]))
                        