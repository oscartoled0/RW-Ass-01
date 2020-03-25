# import necessary libraries
import matplotlib.pyplot as plt
import numpy as np
import BEMmodel as BEM

""" ---- User defined variables ---- """
PR = False
WT = True
N = 100 # annulus spacing
mode = 'cosinus' # 2 possible modes: constant and cosinus distribution of annulus. 
                  # Write 'cosinus' or 'constant' if you want to change the distribution.

[results, CT, CP] = BEM.ExecuteBEM(PR, WT, N, mode)
