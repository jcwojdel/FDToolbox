#!/usr/bin/env python

from fdtoolbox.calculation_set import *
from fdtoolbox.utility import *
from fdtoolbox.linear_expansion import *

import sys
loggable.log_level = LOG_ALLINFO
calculation_set.TRANSLATIONAL_MODE_THRESHOLD = 0.001

saxis = [1., 0., 0.]
usecache = True
directory = "/data/jacek/BiFeO3_100"

  
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

ev=linalg.eig(lin_exp.B_m_n)
print mat2str(ev[0])
print mat2str(ev[1].T)

saxis = [0., 1., 0.]
usecache = True
directory = "/data/jacek/BiFeO3_mono/rotated/thick-0.06"
  
print 'Reading from directory: %s' % directory
print 'Using saxis: %s' % `saxis`
print 'Using cache: %s' % `usecache`

strainedcs=calculation_set.read_directory(directory, saxis, usecache)

strainedcs.set_ionic('.*/calc_.._.*[0123]')
strainedcs.set_groundstate('.*/calc_00_0')
strainedcs.set_lattice('.*/calc_00_...')

strainedcs.try_fix_displacements()  
strainedcs.try_fix_polarizations()

#strainedcs.add_lattice_constraint([[1, 0, 0],
#                                   [0, 1, 0]], ['xyshear', 'oopshear'])
#                                   [0, 1, 0]], ['xyshear'])

strainedlin_exp = linear_expansion( strainedcs )
strainedlin_exp.calculate_expansion_coefficients()

me, units = lin_exp.magneto_electric_coupling(lattice=False)
print 'Clamped cell ME coupling (%s)'%units
print mat2str( me )
print 'Maximum coupling: %f %s'%(sqrt(max(linalg.eig(dot(me,me.T))[0])), units)

me, units = lin_exp.magneto_electric_coupling()
print 'Full ME coupling (%s)'%units
print mat2str( me )
print 'Maximum coupling: %f %s'%(sqrt(max(linalg.eig(dot(me,me.T))[0])), units)


print "Original eigenvalues"
print mat2str( sort( linalg.eig(lin_exp.B_m_n)[0] ) )
print mat2str( sort( linalg.eig(lin_exp.B_j_k)[0] ) )
print mat2str( sort( linalg.eig(lin_exp.Bhat_j_k)[0] ) )


print "Strained eigenvalues"
print mat2str( sort( linalg.eig(strainedlin_exp.B_m_n)[0] ) )
print mat2str( sort( linalg.eig(strainedlin_exp.B_j_k)[0] ) )
print mat2str( sort( linalg.eig(strainedlin_exp.Bhat_j_k)[0] ) )

rotation=array([[  0.707106781186548,   0.460826639679656,   0.536319688396349],
                [ -0.707106781186547,   0.460826639679656,   0.536319688396349],
                [  0.000000000000000,  -0.758470577097829,   0.651707283737789]])

#rotation = dot(rotation, array([[ cos(2*pi/3),-sin(2*pi/3), 0],
#                                [ sin(2*pi/3), cos(2*pi/3), 0],
#                                [ 0,           0,           1]]) )

#rotation = dot(rotation, array([[ 1, 0, 0],
#                                [ 0, -1, 0],
#                                [ 0,           0,           -1]]) )

K = strainedlin_exp.B_m_n.copy()
for i in range(10):
  for j in range(10):
    K[3*i:3*i+3,3*j:3*j+3] = rotatetensor( K[3*i:3*i+3,3*j:3*j+3], linalg.inv(rotation) )
    
C = array(strainedlin_exp.B_j_k.copy()).reshape((3,3,3,3))
C = rotatetensor( C, linalg.inv(rotation) )
C = mat(C.reshape((9,9)))

Bmj = strainedlin_exp.B_m_j.copy()
for i in range(10):
  t = array(Bmj[i*3:i*3+3,:]).reshape((3,3,3))
  t = rotatetensor(t, linalg.inv(rotation))
  Bmj[i*3:i*3+3,:] = t.reshape((3,9))

Bma = strainedlin_exp.B_m_alpha.copy()
for i in range(10):
  Bma[i*3:i*3+3,:] = rotatetensor( Bma[i*3:i*3+3,:], linalg.inv(rotation) )

Baj = array(strainedlin_exp.B_alpha_j.copy()).reshape((3,3,3))
Baj = rotatetensor( Baj, linalg.inv(rotation) )
Baj = mat(Baj.reshape((9,3)))


lin_exp.B_m_n = K
lin_exp.B_j_k = C
lin_exp.B_m_j = Bmj

lin_exp.B_m_alpha = Bma
#lin_exp.B_alpha_j = Baj

lin_exp.calculate_expansion_coefficients(False)

print "Rotated strained eigenvalues"
print mat2str( sort( linalg.eig(lin_exp.B_m_n)[0] ) )
print mat2str( sort( linalg.eig(lin_exp.B_j_k)[0] ) )
print mat2str( sort( linalg.eig(lin_exp.Bhat_j_k)[0] ) )

me, units = lin_exp.magneto_electric_coupling(lattice=False)
print 'Clamped cell ME coupling (%s)'%units
print mat2str( me )
print 'Maximum coupling: %f %s'%(sqrt(max(linalg.eig(dot(me,me.T))[0])), units)

me, units = lin_exp.magneto_electric_coupling()
print 'Full ME coupling (%s)'%units
print mat2str( me )
print 'Maximum coupling: %f %s'%(sqrt(max(linalg.eig(dot(me,me.T))[0])), units)

print "Sanity checks (pairs of the values should be equal):"
print "Force constants for moving atoms along [111]"
for i in range(10):
  dir = mat(zeros((1,30)))
  dir[0,3*i:3*i+3] = strainedcs.groundstate.unit_cell.sum(axis=0)
  dir[0,3*i:3*i+3] /= sqrt(dot(dir[0,3*i:3*i+3],dir[0,3*i:3*i+3].T))
  de = dot( dir, dot( strainedlin_exp.B_m_n, dir.T ) )
  
  dir[0,3*i:3*i+3] = [0., 0., 1.]
  
  print de, "==", dot( dir, dot( K, dir.T ) )
  
  
print "Stress tensor along [111]"
dir = strainedcs.groundstate.unit_cell.sum(axis=0)
dir /= sqrt(dot(dir,dir.T))
dir = dot( dir.T, dir ).reshape(1,9)
de = dot(dir, dot(strainedlin_exp.B_j_k, dir.T) )

dir = mat([[0, 0, 0, 0, 0, 0, 0, 0, 1]])
print de, "==", dot(dir, dot(C, dir.T) )

print "Stress tensor along Z'"
dir = mat([[0, 0, 0, 0, 0, 0, 0, 0, 1]])
de = dot(dir, dot(strainedlin_exp.B_j_k, dir.T) )

dir = mat([[ 0., -3.93100486, 3.37769423]])
dir /= sqrt(dot(dir,dir.T))
dir = dot( dir.T, dir ).reshape(1,9)
print de, "==", dot(dir, dot(C, dir.T) )

print "Ion-relaxed stress tensor along [111]"
dir = strainedcs.groundstate.unit_cell.sum(axis=0)
dir /= sqrt(dot(dir,dir.T))
dir = dot( dir.T, dir ).reshape(1,9)
de = dot(dir, dot(strainedlin_exp.Bhat_j_k, dir.T) )

dir = mat([[0, 0, 0, 0, 0, 0, 0, 0, 1]])
print de, "==", dot(dir, dot(lin_exp.Bhat_j_k, dir.T) )

print "Ion-relaxed stress tensor along Z'"
dir = mat([[0, 0, 0, 0, 0, 0, 0, 0, 1]])
de = dot(dir, dot(strainedlin_exp.Bhat_j_k, dir.T) )

dir = mat([[ 0., -3.93100486, 3.37769423]])
dir /= sqrt(dot(dir,dir.T))
dir = dot( dir.T, dir ).reshape(1,9)
print de, "==", dot(dir, dot(lin_exp.Bhat_j_k, dir.T) )

ev=linalg.eig(lin_exp.B_m_n)
print mat2str(ev[0])
print mat2str(ev[1].T)
