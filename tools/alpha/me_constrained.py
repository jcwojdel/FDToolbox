#!/usr/bin/env python

import sys
import os

sys.path.append(os.path.join(sys.path[0],'../../'))
from fdtoolbox.calculation_set import *
from fdtoolbox.utility import *
from fdtoolbox.linear_expansion import *

try:
  directory = sys.argv[1]
except:
  directory = "."
  
try:
  saxis = [float(v) for v in sys.argv[2].split('/')]  
except:
  saxis = [0., 0., 1.]

try:
  usecache = argv[3] != 'nocache'
except:
  usecache = True
  
print 'Reading from directory: %s' % directory
print 'Using saxis: %s' % `saxis`
print 'Using cache: %s' % `usecache`

cs=calculation_set.read_directory(directory, saxis, usecache)

cs.set_ionic('.*/calc_.._.*[0123]')
cs.set_lattice('.*/calc_00_...')
cs.set_groundstate('.*/calc_00_0')

calculation_set.BPH_MAXDIST = 0.25
calculation_set.BPH_CORRECTION = 0.5

cs.try_fix_displacements()  
cs.try_fix_polarizations()

lin_exp = linear_expansion( cs )
lin_exp.calculate_expansion_coefficients()

me, units = lin_exp.magneto_electric_coupling(lattice=False)
print 'Clamped cell ME coupling (%s)'%units
print mat2str( me )
print 'Maximum coupling: %f %s'%(sqrt(max(linalg.eig(dot(me,me.T))[0])), units)

me, units = lin_exp.magneto_electric_coupling()
print 'Full ME coupling (%s)'%units
print mat2str( me )
print 'Maximum coupling: %f %s'%(sqrt(max(linalg.eig(dot(me,me.T))[0])), units)


 
rhombo_cell = cs.groundstate.unit_cell
pcubic_cell = dot( mat([[-1,  1,  1],
                        [ 1, -1,  1],
                        [ 1,  1, -1]]),  rhombo_cell )

planelist = [
             array([[1, 0, 0],
                    [0, 1, 0]]),
             array([[0, 1, 0],
                    [0, 0, 1]]),
             array([[0, 0, 1],
                    [1, 0, 0]]),
             array([[1, 1, 0],
                    [1, 0, 0]]),
             array([[0, 1, 1],
                    [0, 1, 0]]),
             array([[1, 0, 1],
                    [0, 0, 1]]),
             array([[1,-1, 0],
                    [0, 1,-1]])
             ]
namelist = ["[001]","[001]","[001]",
            "[011]","[011]","[011]",
            "[111]"]
for name, plane in zip(namelist,planelist):
  plane = plane*pcubic_cell
  
  for shears in [['xyshear'],['xyshear','oopshear']]:
    print '--------- Plane %s, shears blocked: %s ------------'%(name,' '.join(shears))
    print 'XY plane defined by:'
    print mat2str(plane)
    
    cs.clear_lattice_constraints()
    cs.add_lattice_constraint( plane, shears)
    lin_exp.calculate_expansion_coefficients()
  
    me, units = lin_exp.magneto_electric_coupling()
    print 'Full ME coupling (%s)'%units
    print mat2str( me )
    print 'Maximum coupling: %f %s'%(sqrt(max(linalg.eig(dot(me,me.T))[0])), units)
    
    vx = plane[0,:]
    vy = plane[1,:]
    vz = cross(vx, vy)
    
    a = math.atan2(vz[0,0],vz[0,2])
    r = mat([[ cos(a), 0.,  sin(a)],
             [  0.,    1.,   0.   ],
             [-sin(a), 0.,  cos(a)]])
    
    a = math.atan2((vz*r)[0,1], (vz*r)[0,2])
    r = r * mat([[ 1.,  0.,      0.   ],
                 [ 0., cos(a),  sin(a)],
                 [ 0.,-sin(a),  cos(a)]])
  
    a = math.atan2(((vx+vy)*r)[0,0],((vx+vy)*r)[0,1])
    r = r * mat([[ cos(a), sin(a),   0.],
                 [-sin(a), cos(a),   0.],
                 [  0.,     0.,      1.]])
     
    print 'In pseudocubic cell oriented as follows:'
    print mat2str( dot(pcubic_cell,r) )
    print 'Original [111] axis is now:'
    print mat2str( dot(array([0, 0, 1]),r) )
    print 'Magnetic easy axis is now:'
    print mat2str( dot(array(saxis),r) )
    print 'ME tensor looks like this:'
    print mat2str( dot(linalg.inv(r), dot(me,r) ) )
    print '------------------------------------------'
