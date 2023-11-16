#!/usr/bin/env python

# this code is to show energies for adjusting casm fitting
import numpy as np
import os, sys
import pandas as pd
pd.set_option('display.max_rows', 1000)  


def read_file(fname): #read data from file and return as a pandas dataframe (header contains a hashtag)
    with open(fname) as f:
        header_line = f.readline()
    header = header_line.strip().split()[1:]
    df = pd.read_csv(fname, delim_whitespace=True, comment='#',names=header,index_col = 0)
    return df

data_hull = read_file('fit_fit.dat')
data_all = read_file('casm_learn_input')

print('Li3Ta1-xNbxO4   (x = 1 ... 4)')
data = pd.concat([data_hull[['comp(a)','formation_energy','clex(formation_energy)']],data_all[['weight']]],axis =1)
data['Clex-DFT'] = data['clex(formation_energy)'].sub(data['formation_energy'])
data['comp(a)'] = data['comp(a)'].apply(lambda x: x) # convert number of Na to concentration x
data['Clex-DFT'] = data['Clex-DFT'].apply(lambda x: x*1000/4.0) # convert number of unit to meV/f.u.

data.columns = ['x','E_DFT','E_clex','weight', 'E_clex-E_DFT (meV/f.u.)']

#data = data.loc[data['x'] == float(sys.argv[1])]
max_dft = data.loc[data['E_DFT'].idxmax()]
min_dft = data.loc[data['E_DFT'].idxmin()]
max_cx = data.loc[data['E_clex'].idxmax()]
min_cx = data.loc[data['E_clex'].idxmin()]
max_diff = data.loc[data['E_clex-E_DFT (meV/f.u.)'].idxmax()]
min_diff = data.loc[data['E_clex-E_DFT (meV/f.u.)'].idxmin()]

#sort data
data = data.sort_values(by=['E_DFT'],ascending=False)

print(data)

print ("min E_DFT\t", str(min_dft['E_DFT']),'\tat\t',min_dft.name,'\tweight\t',min_dft['weight'])
print ("max E_DFT\t", str(max_dft['E_DFT']),'\tat\t',max_dft.name,'\tweight\t',max_dft['weight'])
print ("min E_CX\t", str(min_cx['E_clex']),'\tat\t',min_cx.name,'\tweight\t',min_cx['weight'])
print ("max E_CX\t", str(max_cx['E_clex']),'\tat\t',max_cx.name,'\tweight\t',max_cx['weight'])
print ("min E_diff\t", str(min_cx['E_clex-E_DFT (meV/f.u.)']),'\tat\t',min_diff.name,'\tweight\t',max_diff['weight'])
print ("max E_diff\t", str(max_cx['E_clex-E_DFT (meV/f.u.)']),'\tat\t',max_diff.name,'\tweight\t',max_diff['weight'])
