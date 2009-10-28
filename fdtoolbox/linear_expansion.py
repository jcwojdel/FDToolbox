from numpy import *
from fdtoolbox.calculation_set import *

class linear_expansion(loggable):
  def __init__(self, calcset):
    self.symmetrize = True
    self.calcset = calcset
    
  def calculate_expansion_coefficients(self, update_from_calcset=True):
    if update_from_calcset is True:
      self.volume = self.calcset.groundstate.volume
      self.B_alpha_beta = zeros((3,3))
      self.B_alpha_mu = zeros((3,3))
      self.B_mu_nu = zeros((3,3))
  
      self.B_m_n     = self.calcset.force_constant_matrix('ion', 'ion')[0]
      self.B_j_k     = self.calcset.force_constant_matrix('lat', 'lat')[0]
      self.B_m_j     = 0.5*( self.calcset.force_constant_matrix('lat', 'ion')[0]+
                             self.calcset.force_constant_matrix('ion', 'lat')[0].T )
      if self.symmetrize:
        self.B_m_n     = 0.5*( self.B_m_n + self.B_m_n.T )
        self.B_j_k     = 0.5*( self.B_j_k + self.B_j_k.T )
  
  
      self.B_m_alpha = -self.calcset.electric_polarization_matrix('ion')[0]
      self.B_alpha_j = -self.calcset.electric_polarization_matrix('lat')[0]
      self.B_m_mu    = self.calcset.magnetic_polarization_matrix('ion')[0]
      self.B_mu_j    = self.calcset.magnetic_polarization_matrix('lat')[0]

    inv_B_m_n = invert_with_warning(self.B_m_n, self.calcset.TRANSLATIONAL_MODE_THRESHOLD,
                                    'Inverting force constant matrix resulted in %d soft modes', 3)      

    self.Bhat_alpha_beta = self.B_alpha_beta -    dot(self.B_m_alpha.T, dot(inv_B_m_n, self.B_m_alpha) )
    self.Bhat_mu_nu      = self.B_mu_nu      -    dot(self.B_m_mu.T,    dot(inv_B_m_n, self.B_m_mu) )
    self.Bhat_j_k        = self.B_j_k        -    dot(self.B_m_j.T,     dot(inv_B_m_n, self.B_m_j) )
    self.Bhat_alpha_mu   = self.B_alpha_mu   -    dot(self.B_m_alpha.T, dot(inv_B_m_n, self.B_m_mu) )
    self.Bhat_alpha_j    = self.B_alpha_j    -    dot(self.B_m_j.T,     dot(inv_B_m_n, self.B_m_alpha) )
    self.Bhat_mu_j       = self.B_mu_j       -    dot(self.B_m_j.T,     dot(inv_B_m_n, self.B_m_mu) )

    inv_B_j_k = invert_with_warning(self.Bhat_j_k, self.calcset.TRANSLATIONAL_MODE_THRESHOLD,
                                    'Inverting elastic constant matrix resulted in %d soft modes', 3)      
    
    self.Bgot_alpha_beta = self.Bhat_alpha_beta -    dot(self.Bhat_alpha_j.T, dot(inv_B_j_k, self.Bhat_alpha_j) )
    self.Bgot_mu_nu      = self.Bhat_mu_nu      -    dot(self.Bhat_mu_j.T,    dot(inv_B_j_k, self.Bhat_mu_j) )
    self.Bgot_alpha_mu   = self.Bhat_alpha_mu   -    dot(self.Bhat_alpha_j.T, dot(inv_B_j_k, self.Bhat_mu_j) )

  def born_charges(self):
    return -self.volume*self.B_m_alpha, "|e|"
  
  def piezoelectric_stress_tensor(self, ionic=True):
    if ionic is True:
      rval = self.Bhat_alpha_j
    else:
      rval = self.B_alpha_j
    return -EA2_TO_CM2*rval, "C m**-2"
    
  def piezoelectric_strain_tensor(self, ionic=True):
    if ionic is True:
      sjk = invert_with_warning(self.Bhat_j_k, self.calcset.TRANSLATIONAL_MODE_THRESHOLD,
                              'Inverting elastic constant matrix resulted in %d soft modes', 3)      
      eaj= self.Bhat_alpha_j
    else:
      sjk = invert_with_warning(self.B_j_k, self.calcset.TRANSLATIONAL_MODE_THRESHOLD,
                              'Inverting elastic constant matrix resulted in %d soft modes', 3)      
      eaj= self.B_alpha_j
      
    return -1000*EA2_TO_CM2*dot(sjk, eaj)/EVA3_TO_GPA, "pC N**-1"
  
  def piezomagnetic_stress_tensor(self, ionic=True):
    if ionic is True:
      rval = self.Bhat_mu_j
    else:
      rval = self.B_mu_j
      
    return -BMA3_TO_AM*rval, "A m**-1"
    
  def piezomagnetic_strain_tensor(self, ionic=True):
    if ionic is True:
      sjk = invert_with_warning(self.Bhat_j_k, self.calcset.TRANSLATIONAL_MODE_THRESHOLD,
                                'Inverting elastic constant matrix resulted in %d soft modes', 3)    
      mmj = self.Bhat_mu_j
    else:
      sjk = invert_with_warning(self.B_j_k, self.calcset.TRANSLATIONAL_MODE_THRESHOLD,
                                'Inverting elastic constant matrix resulted in %d soft modes', 3)    
      mmj = self.B_mu_j
       
    return -BMEV_TO_T*dot(sjk, mmj), "T-1"
    
  def elastic_tensor(self, ionic=True):
    if ionic is True:
      rval = self.Bhat_j_k
    else:
      rval = self.B_j_k
    return EVA3_TO_GPA*rval, "GPa" 
    
  def compliance_tensor(self, ionic=True):
    cjk = self.elastic_tensor(ionic)[0]
    sjk = invert_with_warning(cjk, self.calcset.TRANSLATIONAL_MODE_THRESHOLD*EVA3_TO_GPA,
                              'Inverting elastic constant matrix resulted in %d soft modes', 3)
    return 1000*sjk, "TPa**-1"  
  
  def magneto_electric_coupling(self, lattice=True):
    if lattice is True:
      rval = self.Bgot_alpha_mu
    else:
      rval = self.Bhat_alpha_mu

    return -1e4*COUPLING_TO_GAUSS*rval, "1e-4 g.u."
  
  def force_constant_matrix(self):
    return self.calcset.groundstate.volume*self.B_m_n, "eV Andstrom**-2  "
  
  def electric_susceptibility(self, lattice=True):
    if lattice is True:
      rval = self.Bgot_alpha_beta
    else:
      rval = self.Bhat_alpha_beta
    return -EEAEV_TO_EPSILON0*rval, "epsilon0"
      
  def magnetic_susceptibility(self, lattice=True):
    if lattice is True:
      rval = self.Bgot_mu_nu
    else:
      rval = self.Bhat_mu_nu
    return -1e8*BM2A3EV_TO_1MU0*rval, "1e-8 mu0**-1"
      
