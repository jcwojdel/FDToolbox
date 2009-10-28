SYSTEM   = ( BiMnO3 )*2

ALGO     = N
PREC     = Accurate
ENCUT    = 500.
ADDGRID  = .TRUE.

LREAL    = .FALSE.
LMAXMIX  = 4


ISMEAR   = -5
EDIFF    = 1.e-9
LDAU     = .TRUE.
LDAUTYPE = 2
LDAUL    = -1   2  -1
LDAUU    = 0.0 4.0 0.0
LDAUJ    =  0   0   0

ISPIN    = 2

ISYM     = 0

NELMDL   = -15
NSW      = 0
IBRION   = 1
ISIF     = 3
EDIFFG   = -0.001
LORBIT   = 11

LPLANE  = .TRUE.
LSCALU = .FALSE.
NPAR   = 4
NSIM    = 4

ICHARG = 1
LSORBIT = .TRUE.
SAXIS   = 1 0 0

MAGMOM  = 6*0 0 0 4 0 0 -4 18*0
LWAVE   = .FALSE.
