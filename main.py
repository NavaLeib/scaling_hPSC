
import math
import numpy as np
import random

from clones_dynamic_functions import *
from analysis import Analysis

for cells in  [ 10 ** 5]:#[10**4, 10 ** 5,10**6]  ## 

    for m in [1]:# [10, 10 ** 2, 10 ** 3]: #mean initial clone size

        for p in [4]:  #dayes to passage

            ratio_arr = []
            WT_Mut_num_cells = dict.fromkeys(['batch' + str(0) + '_num_of_mut','batch' + str(0) + '_num_of_wt'])
            # WT_Mut_num_cells = dict.fromkeys(['rep' + str(0) + '_num_of_mut','rep' + str(0) + '_num_of_wt'])
            i=0
            for y in [10]*1: #paercent to passage

                BIC_Arr = []

                for mut_percent in  [0.00]:

                    print('\n \n paercent to passage:', y)

                    model = Model(number_of_species=int(cells / m), mutant_percent=mut_percent, number_of_mutations=1,
                                  net_growth_rates_array=[0.026, 0,
                                                     0.035, 0],
                                  interaction=0, plating_efficiency=1.0, initial_mean_clone_size=m)

 #

                    X_batches = model.split_batches_emergence_of_variants(number_of_batches=100, number_of_steps=10**8, passaging=p, percent_of_pass=y, var_emerge_rate=2*10**-5)

                    for i in range(100):
                        import matplotlib.pyplot as plt
                        X_X = X_batches[i]

                        ## (num of wt, num of mut)
                        WT_Mut_num_cells['batch' + str(i) + '_num_of_mut']= [int(X_X[i, :int(int(cells/m) * mut_percent)].sum(axis=0)) + int(X_X[i,int(cells):].sum(axis=0)) for i in range(16 + 1)]
                        WT_Mut_num_cells['batch' + str(i) + '_num_of_wt']= [int(X_X[i,int(int(cells/m) * mut_percent):int(cells)].sum(axis=0)) for i in range(16 + 1)]

import pandas as pd
df = pd.DataFrame(WT_Mut_num_cells)

df.to_csv('ScaleOut_'+ 'Batches100_'
           'Cells' + str(cells) +
            '_MeanInitClone' + str(m) +
            '_Days' + str(p) +
            '_PercentY' + str(y) +
            '_mut_percent' + str(mut_percent) + "EveryPassKeepForAnalysisRemainCells_MutationGrow0035_var_emerge_rate2times10toMinus5.csv", index=False)

