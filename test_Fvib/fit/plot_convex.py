#!/usr/bin/env python
"""
This code is used for cluster expansion figure plotting

Author: Zeyu Deng
Email: dengzeyu@gmail.com
"""
import os,json
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib.patches as patches
from matplotlib import rc
from matplotlib.ticker import MultipleLocator

def set_plot_env():
    plt.tick_params(which="major", length=7, width=1.5)
    # set plot LaTeX font
    #plt.rcParams['text.latex.preamble'] = [
    #    r'\usepackage{siunitx}',  # i need upright \micro symbols, but you need...
    #    r'\sisetup{detect-all}',  # ...this to force siunitx to actually use your fonts
    #    r'\usepackage{helvet}',  # set the normal font here
    #    r'\usepackage[eulergreek,EULERGREEK]{sansmath}'  # load up the sansmath so that math -> helvet##
    #    r'\sansmath'  # <- tricky! -- gotta actually tell tex to use!
    #]
    #plt.rcParams['mathtext.fallback_to_cm'] = 'True'
    # plt.rc('font', family='sans-serif')
    # plt.rc('font', family='sans-serif')
    #plt.rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
    #plt.rc('text', usetex=True)
    plt.rcParams['font.family'] = 'Helvetica'
    plt.rcParams['legend.fancybox'] = False
    # plt.rcParams['legend.loc'] = 'upper right'
    # plt.rcParams['legend.numpoints'] = 2
    # plt.rcParams['legend.fontsize'] = 'large'
    plt.rcParams['legend.framealpha'] = None
    # plt.rcParams['legend.scatterpoints'] = 3
    plt.rcParams['legend.edgecolor'] = 'inherit'
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plt.rcParams['xtick.top'] = True
    plt.rcParams['ytick.right'] = True

def read_fit(fname): #read data from 'fit_fit.dat' and return as a pandas dataframe (header contains a hashtag)
    with open(fname) as f:
        header_line = f.readline()
    header = header_line.strip().split()[1:]
    data_all = pd.read_csv(fname, delim_whitespace=True, comment='#',names=header,index_col = 0)
    data_all = data_all.rename({"clex(formation_energy)":"CE fit","formation_energy":"DFT Energies"},axis='columns')
    data_all["E(Clex)-E(DFT)"] = data_all["CE fit"]-data_all["DFT Energies"]
    data_all[['E(Clex)-E(DFT)','hull_dist(ALL,comp)','CE fit','DFT Energies']]=data_all[['E(Clex)-E(DFT)','hull_dist(ALL,comp)','CE fit','DFT Energies']].apply(lambda e:e*1000/4) # convert to meV/fu
    data_all['x']=data_all['comp(a)']#.apply(lambda x:x/4.0)
    dft_hull = (data_all[data_all['hull_dist(ALL,comp)']<1e-5]).sort_values(by='x')
    clex_hull = (data_all[data_all['clex_hull_dist(ALL,comp)']<1e-5]).sort_values(by='x')
    return data_all,dft_hull,clex_hull

def plot_convexhull(plot_mc=False):
    fig,ax = plt.subplots(nrows=1,ncols=2,figsize=(7,3.5))  
    data_all,dft_hull,clex_hull = read_fit('fit_fit.dat')  
    ax[0].axhline(y=0, xmin=0.0, xmax=1, marker='', linestyle='--', linewidth=0.5, color="black", antialiased=True, label="")
    # plot DFT energies and CE fit
    ax[0].scatter(x=data_all['x'], y=data_all['DFT Energies'],s=8,label='DFT',color='#3498db',marker='o',alpha=0.6)
    ax[0].scatter(x=data_all['x'], y=data_all['CE fit'],s=5,label='CE fit',color='orangered',marker='.',alpha=1)
    # plot DFT and CE convex hull
    ax[0].plot(dft_hull['x'],dft_hull['DFT Energies'],'o-',color='#3498db',alpha=1,lw=1,ms=2,label="DFT hull")
    ax[0].plot(clex_hull['x'],clex_hull['CE fit'],'.-',color='orangered',alpha=1,lw=1,ms=1,label="CE hull")
    ax[0].set(xlim=(-0.1,1.1),ylim=(-5,10),xlabel=r'${\rm x}$ in Li$_{\rm 3}$Ta$_{\rm 1-x}$Nb$_{\rm x}$O$_{\rm 4}$',ylabel=r'Formation Energy (meV/f.u.)')
    ax[0].set_xticks(np.arange(0, 1.25, step=0.25))
    ax[0].legend()  
    ax[0].yaxis.set_minor_locator(MultipleLocator(50))
    # plot CE error
    ax[1].scatter(x=data_all['E(Clex)-E(DFT)'], y=data_all['hull_dist(ALL,comp)'],marker="o",color='#95e1d3',s=10,alpha=0.6)
    ax[1].set(xlim=(-2,2),ylim=(0,10),xlabel='Error of CE (meV/f.u.)',ylabel=r'Energy above Convex Hull (meV/f.u.)')  
    # draw 2 boxes for showing error boundaries
    # rect1 = patches.Rectangle((-20,0),40,300,linewidth=0.5,edgecolor="black",facecolor='none',linestyle='--')
    # rect2 = patches.Rectangle((-10,0),20,200,linewidth=0.5,edgecolor="black",facecolor='none',linestyle='--')   
    # ax[1].add_patch(rect1)
    # ax[1].add_patch(rect2)  
    # ax[1].text(0,310,'300 meV',color="black",ha='center')
    # ax[1].text(0,210,'200 meV',color="black",ha='center')
    ax[1].yaxis.set_minor_locator(MultipleLocator(50))
    ax[1].xaxis.set_minor_locator(MultipleLocator(50))

    fig.tight_layout()
    fig.savefig("convex_hull.pdf", format="pdf", bbox_inches='tight')

def draw_stem(ax,group,x,y,color,marker,label):# do a stem plot for "group" with "x" vs "y" using "color" and 'marker' with "label" on "ax"
    markerline, stemlines, baseline=ax.stem(group.get_eff_property_list(x),group.get_eff_property_list(y),markerfmt='.',use_line_collection=True,label=label)
    plt.setp(stemlines,'color',color,'lw',1)
    plt.setp(markerline,'color',color,'lw',1,'marker',marker,'ms',3)
    plt.setp(baseline,'visible',False)   

set_plot_env()
plot_convexhull(plot_mc=False)
