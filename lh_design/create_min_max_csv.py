#!/usr/bin/env python
import pandas as pd
import sys


def func(MIN, MAX, LOGTRANSFORM, minmax_extension):
    NAME = ['qv_bl','inv','delt','delq','na','baut']
    
    parameters_data = pd.DataFrame({'NAME':NAME, 'MIN':MIN, 'MAX':MAX, 'LOGTRANSFORM':LOGTRANSFORM})
    parameters_data.to_csv(f"parameters_{minmax_extension}.csv", index=False)

def main():
    minmax_extension = sys.argv[1]

    MIN = [7,500,2,-7,10,-2.3]

    MAX = [11,1300,21,-1,500,-1.3]
    
    LOGTRANSFORM = [0, 0, 0, 0, 0, 1]
    
    func(MIN, MAX, LOGTRANSFORM, minmax_extension)

if __name__=="__main__":
    main()
