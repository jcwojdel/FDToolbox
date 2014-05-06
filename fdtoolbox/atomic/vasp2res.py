#
# This utility reads POSCAR or CONTCAR and generates a .res file
# readable by Materials Studio
# Usage: python vasp2res.py <vasp structure file> <name of .res file> <element list>
# Jorge Iniguez, February 2007
#

import sys
import os
import re
from atomic_set import atomic_set


## Main program

# Read POSCAR
poscar = atomic_set()
poscar.load_from_poscar(sys.argv[1])

print poscar.cell_a, poscar.cell_b, poscar.cell_c, poscar.cell_alpha, poscar.cell_beta, poscar.cell_gamma

if poscar.species is None:
  # Read elements from command line...
  if len(sys.argv)>3:
    elements = sys.argv[3].split()
  # ...or read them from POTCAR in working directory
  elif os.access('POTCAR', os.R_OK):
    elements = []
    with open("POTCAR","r") as pc:
      for line in pc:
        if line.startswith("   VRHFIN"):
          elements.append(re.search('=(.+):',line).group(1))
  else:
    print "Cannot guess atom types - reverting to all 'C'"
    elements = len(poscar.num_per_type)*["C"]

  species = elements
  poscar.species=[]
  for i in poscar.num_per_type:
    poscar.species.extend([species[0]]*i)
    species = species[1:]

print poscar.species

poscar.save_to_res(sys.argv[2]+".res")
