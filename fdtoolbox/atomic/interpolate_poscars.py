import sys
import copy
import optparse
import re 
from atomic_set import atomic_set

parser = optparse.OptionParser(usage='Usage: %prog [options] POSCAR1 POSCAR2')
parser.add_option("-s", "--steps", action='store', type='int', dest='num_steps', help="Interpolate in NUMSTEPS number of steps (default 5)", metavar="NUMSTEPS")
parser.add_option("-b", "--begin", action='store', type='int', dest='begin', help="Start interpolation at step BEG (default 0)", metavar="BEG")
parser.add_option("-e", "--end", action='store', type='int', dest='end', help="Stop interpolation at step END (default NUMSTEPS)", metavar="END")
parser.add_option("-o", "--output", action='store', type='string', dest='output', help="Output file pattern (default POSCAR%02d)")
parser.add_option("-a", "--atoms", action='store_true', dest='atoms', help="Print atomic species line in POSCAR (might break VASP)")
parser.add_option("-d", "--direct", action='store_true', dest='direct', help="Write POSCARs in direct coordinates")
parser.add_option("-r", "--arc", action='store_true', dest='arc', help="Write an ARC file instead of POSCARs")

(options,arguments) = parser.parse_args()

if len(arguments) != 2:
  parser.print_usage()

if options.num_steps is None:
  options.num_steps = 5
  
if options.begin is None:
  options.begin = 0
  
if options.end is None:
  options.end = options.num_steps
    
if options.output is None:
  options.output = 'POSCAR%02d'
  
poscar1=atomic_set()
poscar1.load_from_poscar(arguments[0])

poscar2=atomic_set()
poscar2.load_from_poscar(arguments[1])
  
poscar_int = copy.deepcopy(poscar1)

poscar2.align_to(poscar1)
for i in range(options.begin,options.end+1):
  t = float(i)/options.num_steps
  poscar_int.unit_cell = (1-t)*poscar1.unit_cell + t*poscar2.unit_cell
  poscar_int.atoms = (1-t)*poscar1.atoms + t*poscar2.atoms
  if options.arc:
    poscar_int.save_to_arc(options.output, header=True if i==options.begin else False, comment='Interpolation pos %i'%i)
  else:
    poscar_int.save_to_poscar(options.output%i, direct=options.direct, species_line=options.atoms)
  
