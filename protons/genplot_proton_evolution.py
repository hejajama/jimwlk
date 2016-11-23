#!/usr/bin/python
#-*- coding: UTF-8 -*-
import os
import sys
sys.path.append("/nashome2/hejajama/lib/")
sys.path.append("/home/hejajama/lib/")
sys.path.append("/home/heikki/lib")
sys.path.append("/Users/heikki/lib")
import math
from math import pow
from matplotlibhelper import *
import pylab
import scipy.integrate
import numpy as np
import matplotlib as mp

#import scipy.integrate

from matplotlib import rcParams
#rcParams.update({'figure.autolayout': True})

slides=True
coherent = True

minx=-2.5
miny=-2.5
maxy=2.5
maxx=2.5


minx= -1.5#-5
maxx=1.5 #5
miny=minx
maxy=maxx



rc('text',usetex=True)
rc('text.latex',  preamble='\usepackage{amsmath},\usepackage{amssymb},\usepackage{mathtools}')
textsize=21

rc("xtick", labelsize=textsize)
rc("ytick", labelsize=textsize)

fmgev = 5.068

fig = plt.figure()
ipglasma=True

# colors
cdict = {'red': ((0., 1, 1),
                 (0.05, 1, 1),
                 (0.11, 0, 0),
                 (0.66, 1, 1),
                 (0.89, 1, 1),
                 (1, 0.5, 0.5)),
    'green': ((0., 1, 1),
              (0.05, 1, 1),
              (0.11, 0, 0),
              (0.375, 1, 1),
              (0.64, 1, 1),
              (0.91, 0, 0),
              (1, 0, 0)),
        'blue': ((0., 1, 1),
                 (0.05, 1, 1),
                 (0.11, 1, 1),
                 (0.34, 1, 1),
                 (0.65, 0, 0),
                 (1, 0, 0))}

my_cmap = matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict,256)

#set_cmap('gist_earth_r')
#set_cmap('gist_stern_r')
#grid = AxesGrid(fig, 111,
#               nrows_ncols=(2, 2),
#               axes_pad=0.05,
#               share_all=True,
#               label_mode="L",
#               cbar_location="right",
#               cbar_mode="single",
#               )


#cbar_ax = fig.add_axes([.99, .3, .03, .4])


#alphas=0.14
alphas=0.16

ds=0.0004

cvals=[8,1,4,3] #ipgalsma
cvals=[1,2,3,4,5,6,7]
cvals=[0,1,2,3,4,5]
cvals=[0,1,2,3,4]
#steps=[0,50,100,150]
steps=[0,100,200,250]
steps=[0,100,400,600]


rapidities=[]
for s in steps:
    rapidities.append(s*ds*pi*pi/alphas)

for c in cvals:
    config=-1
    fig, axn = plt.subplots(2, 2, sharex=True, sharey=True)
    for i, ax in enumerate(axn.flat):
        ax.set(adjustable='box-forced', aspect='equal')
        config=config+1
        step=steps[config]
        xvals=[]
        yvals=[]
        xarray=[]
        yarray=[]
        densityarray = []

        # read file
        fname="bp_3.0_bq_0.3_qsfluct_y_1.79_ipglasma_m_0.4_jimwlk_m_0.2/proton_" + str(c) + "_steps_" + str(step)
        print fname
        #fname="ipglasma_nq_5/ipglasma-1.5-0.2-proton-" + str(c)
        f = open(fname, "r")
        lines=f.readlines()
        f.close()
        tmprow = []
        xrow=[]
        yrow=[]
        rowindex=0
        
        ycol=3
        if ipglasma:
            ycol=5
        for i in range(len(lines)):
            s=lines[i].split()
            if (len(s)<2 or lines[i][0]=="#"):
                continue
            x = float(s[1])/fmgev
            y = float(s[0])/fmgev
            d = float(s[ycol])
            
            if len(yvals)==0:
                yvals.append(y)
            #print y, yvals[rowindex],rowindex
            if y == yvals[rowindex]: # continue line
                #print "Add value " + str(d)
                tmprow.append(d)
                xrow.append(x)
                yrow.append(y)
            else:
                #print "Append " + str(y) + " densitylist " + str(tmprow)
                # save rows
                newlist = tmprow[:]
                densityarray.append(newlist)
                newx = xrow[:]
                newy = yrow[:]
                xarray.append(newx)
                yarray.append(newy)
                
                # start new row
                xrow=[x]
                yrow=[y]
                yvals.append(y)
                tmprow = [d]
                rowindex = rowindex + 1

        # save last
        newlist = tmprow[:]
        densityarray.append(newlist)
        newx = xrow[:]
        newy = yrow[:]
        xarray.append(newx)
        yarray.append(newy)

        xvals=list(set(xvals))
        
        #for i in range(len(densityarray)):
        #   for j in range(len(densityarray[0])):
    #       print xarray[i][j],yarray[i][j], densityarray[i][j]
        
    #print densityarray
    #    print yarray
    #   print xarray
        
        #pcolormesh
        im = ax.pcolor(np.array(xarray), np.array(yarray), np.array(densityarray),rasterized=True, cmap=my_cmap)
        ax.axis([minx,maxx,miny,maxy])
        #ax.set_xticks([-2,-1,0,1,2])
        #ax.set_yticks([-2,-1,0,1,2])
        if config == 0 or config == 2:
            ax.set_ylabel(r"$y [\mathrm{fm}]$", fontsize=textsize)
        if config == 2 or config == 3:
            ax.set_xlabel(r"$x [\mathrm{fm}]$", fontsize=textsize)

    cax = fig.add_axes([0.78, 0.1, 0.03, 0.8], )

    bar=fig.colorbar(im, cax=cax)

    if not ipglasma:
        bar.set_ticks([0, 0.05, 0.1, 0.15, 0.2, 0.3])
    else:
        bar.set_ticks([0,0.2,0.4,0.6,0.8,1.0,1.2])
        #fig.suptitle(r"$1.0 - \mathrm{Re} \, \mathrm{Tr} \, V(x,y)$", fontsize=textsize)

    plt.subplots_adjust(left=0., bottom=0.115, wspace=-0.45, hspace=0.1)

    titlestring=r"$y="
    for y in rapidities:
        titlestring= titlestring + "%.2f" % round(y,2)
        if y != rapidities[-1]:
            titlestring += ", "

    titlestring+=r"$"
    fig.suptitle(r"$1.0 - \mathrm{Re} \, \mathrm{Tr} \, V(x,y) $ " + titlestring, fontsize=textsize-2)
                 
    #fig.colorbar(im)

    #fig.tight_layout()

    #plt.show()

    file = "./density_" + str(c) + ".pdf"
    print file
    pp = PdfPages(file)
    savefig(pp, format='pdf', dpi=120, bbox_inches='tight')
    pp.close()
