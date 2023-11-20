#!/usr/bin/env python
"""
This code contains a set of tools for monte carlo analysis
Author: Zeyu Deng
Email: dengzeyu@gmail.com
Modified by Ziliang & Tara.
"""
import numpy as np
import pandas as pd
import sys,os,re,json
def read_results_file(fname='results.csv'): #deal with results.csv
    with open(fname) as f:
        header_line = f.readline()
    header = header_line.strip().split()[1:]
    df = pd.read_csv(fname, delim_whitespace=True, comment='#',names=header)
    return df
results = read_results_file()
"""
You can change the column to show below in the print function
"""
columns = ['T','param_chem_pot(a)','<formation_energy>','<potential_energy>','<comp(a)>','<clex_hull_dist(ALL,comp)>','susc_n(Nb,Nb)']

with pd.option_context("display.max_rows", 10000):
    print(results[columns])
#print('\n=========== Below are all columns in this table ====================\n')
#print(results.columns.values)
os.system('cp results_table.txt results_table_prev.txt')
results.to_csv('results_table.txt',columns=columns,sep='\t',float_format='%10.4f',index=False)
