#!/usr/bin/env python

from fdtoolbox.calculation_set import *
from fdtoolbox.linear_expansion import *
from fdtoolbox.utility import *

import sys

try:
  directory = sys.argv[1]
except:
  directory = "."

try:
  usecache = sys.argv[2] != 'nocache'
except:
  usecache = True
  

  
print 'Reading from directory: %s' % directory
print 'Using cache: %s' % `usecache`

saxis = [0., 0., 1.]
cs=calculation_set.read_directory(directory, saxis, usecache)

cs.set_ionic('.*/calc_.._.*[0123]')
cs.set_lattice('.*/calc_00_...')
cs.set_groundstate('.*/calc_00_0')

cs.try_fix_displacements()  

lin_exp = linear_expansion( cs )
lin_exp.calculate_expansion_coefficients()

fc_matrix, units = lin_exp.force_constant_matrix()
print 'Force constant matrix eigenvalues (%s)'%units
print mat2str( sort( linalg.eig(fc_matrix)[0]) )

elastic_c, units = lin_exp.elastic_tensor(ionic=False) 
print 'Clamped-ion elastic constant matrix in Voigt notation (%s)'%units
print mat2str( tensor2voigt(elastic_c, [0,1]), '% 12.4f' )
print 'Eigenvalues: ',mat2str( sort(linalg.eig(elastic_c)[0]) )   

elastic_c, units = lin_exp.elastic_tensor() 
print 'Relaxed-ion elastic constant matrix in Voigt notation (%s)'%units
print mat2str( tensor2voigt(elastic_c, [0,1]), '% 12.4f' )   
print 'Eigenvalues: ', mat2str( sort(linalg.eig(elastic_c)[0]) )   

compl_s, units = lin_exp.compliance_tensor(ionic=False) 
print 'Clamped-ion compliance matrix in Voigt notation (%s)'%units
print mat2str( tensor2voigt(compl_s, [0,1]), '% 12.4f' )
   
compl_s, units = lin_exp.compliance_tensor() 
print 'Relaxed-ion compliance matrix in Voigt notation (%s)'%units
print mat2str( tensor2voigt(compl_s, [0,1]), '% 12.4f' )
