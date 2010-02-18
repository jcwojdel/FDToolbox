#!/usr/bin/env python

import sys
import os

sys.path.append(os.path.join(sys.path[0],'../../'))
from fdtoolbox.calculation_set import *
from fdtoolbox.utility import *
from fdtoolbox.linear_expansion import *

directory1 = sys.argv[1]

directory2 = sys.argv[2]

usecache = True
saxis = [0., 0., 1.]
  
print 'Reading from directory: %s' % directory
print 'Using saxis: %s' % `saxis`
print 'Using cache: %s' % `usecache`

loggable.log_level=LOG_WARNING
cs1=calculation_set.read_directory(directory1, saxis, usecache)

cs1.set_ionic('.*/calc_.._.*[0123]')
cs1.set_lattice('.*/calc_00_...')
cs1.set_groundstate('.*/calc_00_0')

calculation_set.BPH_MAXDIST = 0.25
calculation_set.BPH_CORRECTION = 0.5

cs1.try_fix_displacements()  
cs1.try_fix_polarizations()

cs2=calculation_set.read_directory(directory2, saxis, usecache)

cs2.set_ionic('.*/calc_.._.*[0123]')
cs2.set_lattice('.*/calc_00_...')
cs2.set_groundstate('.*/calc_00_0')

cs2.try_fix_displacements()  
cs2.try_fix_polarizations()


print '==============================================================================='
print '=                             UNPERTURBED STRUCTURE                           ='
print '==============================================================================='

print 'Unit cell (A)'
print mat2str(cs1.groundstate.unit_cell)
print 'a=%f b=%f c=%f, alpha=%f, beta=%f gamma=%f\n'%(cs1.groundstate.cell_a, 
                                                    cs1.groundstate.cell_b,
                                                    cs1.groundstate.cell_c,
                                                    cs1.groundstate.cell_alpha,
                                                    cs1.groundstate.cell_beta,
                                                    cs1.groundstate.cell_gamma)

print 'Atomic positions (A)'
print mat2str(cs1.groundstate.atoms)

pc_calc = copy.copy(cs1.groundstate)
pc_calc.unit_cell = dot(mat("-1 1 1; 1 -1 1; 1 1 -1"), cs1.groundstate.unit_cell )
print 'PC unit cell (A)'
print mat2str( pc_calc.unit_cell )
print 'a=%f b=%f c=%f, alpha=%f, beta=%f gamma=%f\n'%(pc_calc.cell_a, 
                                                    pc_calc.cell_b,
                                                    pc_calc.cell_c,
                                                    pc_calc.cell_alpha,
                                                    pc_calc.cell_beta,
                                                    pc_calc.cell_gamma)




print 'Magnetic moment (mu_B/cell)'
print mat2str(cs1.groundstate.magnetization)


lin_exp1 = linear_expansion( cs1 )
lin_exp1.calculate_expansion_coefficients()

lin_exp2 = linear_expansion( cs2 )
lin_exp2.calculate_expansion_coefficients()

lin_exp = lin_exp1
lin_exp.B_m_alpha = lin_exp2.B_m_alpha
lin_exp.B_alpha_j = lin_exp2.B_alpha_j
lin_exp.calculate_expansion_coefficients(False)

print '==============================================================================='
print '=                         CLAMPED CELL PROPERTIES                             ='
print '==============================================================================='

fc, units = lin_exp.force_constant_matrix()
print 'Force constant matrix (%s)'%units
print mat2str(fc)
ee, _ = linalg.eig(fc)
print 'Matrix eigenvalues'
print mat2str(ee)

bc, units = lin_exp.born_charges()
print 'Born charges (%s)'%units
print mat2str(bc)

ms, units = lin_exp.magnetic_strengths()
print 'Magnetic strengths (1e-3 * %s)'%units
print mat2str(1000.*ms)


me, units = lin_exp.magneto_electric_coupling(lattice=False)
print 'Clamped cell ME coupling (%s)'%units
print mat2str( me )

def explain_me_tensor(me):
  max_v, max_e = linalg.eig(dot(me,me.T))
  
  print 'Coupling strengths: %s'%mat2str(sqrt(max_v))
  for i in range(3):
    print 'E = %s'%mat2str(max_e[:,i].T)
    dm = dot(me.T,max_e[:,i])
    print 'dM = %s'%mat2str(dm.T)
    print '|dM| = %f'%sqrt(dot(dm.T,dm))
    
def decompose_me_tensor(lin_exp, ionic_only = True, ionic_first = True):
  num_phonons = lin_exp.B_m_n.shape[0]
  decomposition = zeros([num_phonons,10])
  
  if ionic_first:
    evals, _ = linalg.eig(lin_exp.B_m_n)
    sort_evals = sort(evals)
      
    for i in range(num_phonons):
      decomposition[i,0] = sort_evals[i]*lin_exp.calcset.groundstate.volume
       
      inv_B_m_n = safe_inv( lin_exp.B_m_n, -3, iterate_components = nonzero(evals<=sort_evals[i])[0] )[0]
      
      decomposition[i,1:] = (1e4*COUPLING_TO_GAUSS*dot(lin_exp.B_m_alpha.T, dot(inv_B_m_n, lin_exp.B_m_mu) )).flatten()
  else:
    pass

  print 'Cumulative ionic eigenmode decomposition'
  print ' eigenval  a_11     a_12     a_13     a_21     a_22     a_23     a_31     a_32     a_33'
  print mat2str(decomposition, "%8.4f")
  
  #decomposition = zeros([9,10])
  decomposition = dot(ones([9,1]), array([decomposition[-1,:]]))
  if ionic_only is False :
    if ionic_first:
      evals, _ = linalg.eig(lin_exp.Bhat_j_k)
      sort_evals=sort(evals)
      
      for i in range(9):
        decomposition[i, 0] = sort_evals[i]*lin_exp.calcset.groundstate.volume
        if decomposition[i, 0]> 1e5:
          decomposition[i, 0]=inf
        inv_B_j_k = safe_inv( lin_exp.Bhat_j_k, lin_exp.calcset.TRANSLATIONAL_MODE_THRESHOLD, iterate_components = nonzero(evals<=sort_evals[i])[0] )[0]
        decomposition[i,1:] = decomposition[i,1:]+(1e4*COUPLING_TO_GAUSS*dot(lin_exp.Bhat_alpha_j.T, dot(inv_B_j_k, lin_exp.Bhat_mu_j) )).flatten()
        
         
    else:
      pass
    
    print 'Cumulative elastic eigenmode decomposition'
    print ' eigenval  a_11     a_12     a_13     a_21     a_22     a_23     a_31     a_32     a_33'
    print mat2str(decomposition, "%8.4f")
 
explain_me_tensor(me)
decompose_me_tensor(lin_exp)

print '==============================================================================='
print '=                         RELAXED CELL PROPERTIES                             ='
print '==============================================================================='


planelist = [
             mat([[1.,0. , 0.],
                  [0., 1., 0.]])
             ]
namelist = ["[001]"]
for name, plane in zip(namelist,planelist):
  
  for shears in [['xyshear'],['xyshear','oopshear']]:
    print '--------- Plane %s, shears blocked: %s ------------'%(name,' '.join(shears))
    print 'XY plane defined by:'
    print mat2str(plane)
    
    cs.clear_lattice_constraints()
    cs.add_lattice_constraint( plane, shears)
    lin_exp.calculate_expansion_coefficients()

    de, units = lin_exp.elastic_tensor(False)
    eval, evec = linalg.eig(de)

    print 'elastic tensor clamped ions (%s)'%units
    print '   xx       xy       xz       yx       yy       yz       zx       zy       zz'
    de[de>1e5]=inf
    de[de<-1e5]=-inf
    print mat2str(de)
    
    print 'tensor\'s eigenvalues (%s)'%units
    eval[eval>1e5]=inf
    eval[eval<-1e5]=-inf
    print mat2str(eval)

    de, units = lin_exp.elastic_tensor()
    eval, evec = linalg.eig(de)

    print 'elastic tensor (%s)'%units
    print '   xx       xy       xz       yx       yy       yz       zx       zy       zz'
    de[de>1e5]=inf
    de[de<-1e5]=-inf
    print mat2str(de)
    
    print 'tensor\'s eigenvalues (%s)'%units
    eval[eval>1e5]=inf
    eval[eval<-1e5]=-inf
    print mat2str(eval)

    de, units = lin_exp.piezoelectric_strain_tensor()
    print 'Piezoelectric strain tensor (%s)'%units
    print '   xx       xy       xz       yx       yy       yz       zx       zy       zz'
    print mat2str(de.T)
    
    de, units = lin_exp.piezomagnetic_strain_tensor()
    print 'Piezomagnetic strain tensor (1e-8 %s)'%units
    print '   xx       xy       xz       yx       yy       yz       zx       zy       zz'
    print mat2str(1e8*de.T)

    me, units = lin_exp.magneto_electric_coupling()
    print 'Full ME coupling (%s)'%units
    print mat2str( me )

    explain_me_tensor(me)
    decompose_me_tensor(lin_exp, ionic_only=False)
    
#    fc, units = lin_exp.force_constant_matrix()
#    inv_B_j_k = invert_with_warning(lin_exp.B_j_k, lin_exp.calcset.TRANSLATIONAL_MODE_THRESHOLD,
#                                    'Inverting elastic constant matrix resulted in %d soft modes', 3)      
#    soft_fc = lin_exp.B_m_n - dot( lin_exp.B_m_j, dot( inv_B_j_k, lin_exp.B_m_j.T ) )
#    soft_fc*= lin_exp.calcset.groundstate.volume
#    
#    e1, _ = linalg.eig( fc )
#    e2, _ = linalg.eig( soft_fc )
#    
#    print 'Eigenvalues of force constant matrix w./w.o. lattice softening (%s)'%units 
#    print mat2str( array([sort(e1),sort(e2)]) )
    
    

    