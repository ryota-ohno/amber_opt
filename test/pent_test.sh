#!/bin/bash
#PJM -L "rscunit=fx"
#PJM -L "rscgrp=fx-small"
#PJM -L "node=1"
#PJM -L "elapse=0:10:00"
#PJM -j
#PJM -S
#PJM "--norestart"

module load amber
parmchk2 -i pent_test.mol2 -f mol2 -o pent_test.frcmod
tleap -f pent_tleap_test.in
sander -O -i FF_calc.in -o pent_test.out -p pent_test.prmtop -c pent_test.inpcrd -r min.rst -ref pent_test.inpcrd 

