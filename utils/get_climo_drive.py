#----------------------------------------------
# This script drive the get_climo computation
# by looping through a number of experiments.
#----------------------------------------------

import numpy as np
import os

#list of all experiments to be processed
exps=[ \
      #'CMIP_RRTMG_UMRad_scat.ne30_ne30.cori-knl', \
      'CMIP_RRTMG_UMRad_scat.ne30_ne30.cori-knl-ens1',\
      #'CMIP_RRTMG_UMRad_scat.ne30_ne30.cori-knl-ens2',\
      #'CMIP_RRTMG_UMRad_scat.ne30_ne30.cori-knl-ens3',\
      #'CMIP_RRTMG_UMRad_abs.ne30_ne30.cori-knl',\
      #'CMIP_RRTMG_UMRad_abs.ne30_ne30.cori-knl-ens1' ] #,\
      'CMIP_RRTMG_UMRad_abs.ne30_ne30.cori-knl-ens2' ] #,\
      #'CMIP_RRTMG_UMRad_abs.ne30_ne30.cori-knl-ens3']

      #'AMIP_RRTMG_UMRad_abs.ne30_ne30.cori-knl',  \
      #'AMIP_RRTMG_UMRad_scat.ne30_ne30.cori-knl']

#command to execute the main python file
cmd = 'python ./get_climo.py'

#write experiment name to file and do calculation one-by-one
for ex in exps:
    fout=open('./exp_id.txt','w')
    fout.write(ex)
    fout.close()
    os.system(cmd)
