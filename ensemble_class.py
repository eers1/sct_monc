#!/usr/bin/env python3.8
# Class for using ppe design 
# R. Sansom 8/8/22

import numpy as np

class Ensemble:
    def __init__(self, em_csv_path, val_csv_path):
        load = lambda path, skiprows : np.loadtxt(path, delimiter=',', skiprows=skiprows)
        
        self.em = load(em_csv_path, 0 if "spin" in em_csv_path else 1)
        self.val = load(val_csv_path, 0 if "spin" in val_csv_path else 1)
        self.design = np.concatenate([self.em, self.val], axis=0)
        self.size = len(self.design)
        self.param_names = ['qv_bl','inv','delt','delq','na','baut']
        self.param_labels = ['$BL~q_{v}$', '$BL~z$', r'$\Delta~\theta$', '$\Delta~q_{v}$', '$BL~N_{a}$', '10^{$b_{aut}$}']
        self.param_mins = [7, 500, 2, -7, 10, 0.0051]
        self.param_maxs = [11, 1300, 21, -1, 500, 0.05]
        self.ax_names = [(7, 11), (500, 1300), (2, 21), (-7, -1), (10, 500), (10**(-2.3), 10**(-1.3))]