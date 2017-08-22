import numpy as np
import pandas as pd
# import json
import matplotlib.pyplot as plt
from io_local.read_func import Time_Intervals_read

"""
Test - unlimited picewise function
"""
#input_file = open('input.dat')
#input_str = input_file.read()
#Param = json.loads(input_str)

#print(list_of_periods)
#print(time_lag)
# time_lag = [0,5,10,15,20]

df = pd.read_csv('Fasten-RR2Pct01.dat', sep=',',header=None)
print(df.values)

time_lag = Time_Intervals_read()
arg_lag = [1,3,2,5,4]
time_exp=np.linspace(0, 48, 100)

def interval_conds(time_exp_local,time_lag_local):
    conds_left = []
    conds_right = []
    for ind in time_lag_local:
        conds_left.append((ind <= time_exp_local))
        conds_right.append((time_exp_local < ind))
        print(ind)
    
    del conds_left[-1]
    del conds_right[0]

    conds = np.logical_and(conds_left, conds_right)
    return conds

def func_zero_order_on_interval(arg_lag_local, time_exp_local):
    from DD_basic_models.zero_order_model import C_zero_order
    funcs = []
    for ind in arg_lag_local:
        funcs.append(lambda time_exp_local: C_zero_order(ind,time_exp_local))
    return funcs

conds2 = interval_conds(time_exp,time_lag)
funcs2 = func_zero_order_on_interval(arg_lag, time_exp)

con_exp2 = np.piecewise(time_exp,conds2, funcs2)

plt.plot(time_exp, con_exp2, '-x')

np.savetxt('data.csv', np.transpose([time_exp, con_exp2]), fmt='%4.2g', delimiter=';  ')

