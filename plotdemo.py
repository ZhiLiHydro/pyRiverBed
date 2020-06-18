# coding: utf-8
"""
Plotting script for pyRiverBed
pyRiverBed: Generate Synthetic Riverbed Topography for Meandering Rivers
MIT License
Author: Zhi Li, Univ. of Illinois Urbana-Champaign
Contact: zhil2[at]illinois[dot]edu

"""

import os
import sys
import numpy as np
import matplotlib
if sys.platform == 'darwin':  ## for tkinter-matplotlib
    matplotlib.use('tkagg')   ## crashing bug in MacOS
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.gridspec as gridspec


def make_figure(x, y, allxyz, cur_ori, cur_flt, cur_lag, s, beck_bed, 
                params, t, oxbowxList, oxbowyList, centerlinexList, centerlineyList):
    """
    Draw and output the figure of computing results.
    
    Parameters
    ----------
    x : ndarray
        x coordinate of river centerline
    y : ndarray
        y coordinate of river centerline
    allxyz : ndarray
        point cloud in 3-column xyz format
    cur_ori : ndarray
        original curvature
    cur_flt : ndarray
        filtered curvature
    cur_lag : ndarray
        filtered curvature with considering phase lag
    s : ndarray 
        streamwise distance
    beck_bed : ndarray
        large matrix storing synthetic bed topography data
    params : tuple
        variables needed in the plotting script
    t : integer
        current # of time step
    oxbowxList : list
        x cooridinates of oxbow lakes
    oxbowyList : list
        y cooridinates of oxbow lakes
    centerlinexList : list
        x cooridinates of centerlines
    centerlineyList : list
        y cooridinates of centerlines
    
    Returns
    -------
    None

    """
    WIDTH, DEPTH, SLOPE, NUM, LAG, FNAME, \
    MIGRATION, DT, TSTEPS, GPRINT, JPG_DIRS = params

    print('+> Plotting...', end='')

##    plt.style.use('ggplot')
    gs = gridspec.GridSpec(8, 8)
    
    plt.figure(figsize=(16, 12))
    
    ## Subplot a: point cloud
    plt.subplot(gs.new_subplotspec((0, 0), colspan=6, rowspan=6))
    if MIGRATION:
        plt.title('(a) Synthetic bed in x-y coordinate\n(exposed sand bar visualization)\nT=' + str(int(t*DT/86400)) + 'd', fontsize=16)
        for oxbowx, oxbowy in zip(oxbowxList, oxbowyList):
            plt.fill(oxbowx, oxbowy, color='tan', label='_nolegend_', zorder=0)
            plt.plot(oxbowx, oxbowy, color='peru', linewidth=3, label='_nolegend_', zorder=0)
    else:
        plt.title('(a) Synthetic bed in x-y coordinate\n(exposed sand bar visualization)', fontsize=16)
    waterlevel = DEPTH/2
    argwet = np.argwhere(allxyz[:, 2]-np.tile((s[-1]-s)*SLOPE, NUM*2+1)<waterlevel)
    argdry = np.argwhere(allxyz[:, 2]-np.tile((s[-1]-s)*SLOPE, NUM*2+1)>=waterlevel)
    wetcolor, drycolor, bgcolor = (70/256,80/256,60/256), (200/256,187/256,160/256), (99/256,120/256,84/256)
    plt.scatter(allxyz[argwet, 0], allxyz[argwet, 1], color=wetcolor, s=2, label='_nolegend_', zorder=1)
    plt.scatter(allxyz[argdry, 0], allxyz[argdry, 1], color=drycolor, s=2, label='_nolegend_', zorder=1)
##    plt.tripcolor(allxyz[:, 0], allxyz[:, 1], allxyz[:, 2]/DEPTH, 
##                cmap='Spectral_r', vmin=-1, vmax=1+SLOPE*np.max(s)/DEPTH, label='_nolegend_')
    plt.plot(x[0], y[0], marker='o', color=wetcolor, markersize=8, linestyle='', label='river channel', zorder=2)
    plt.plot(x[0], y[0], marker='o', color=drycolor, markersize=8, linestyle='', label='exposed sand bar', zorder=2)
    plt.plot(x[0], y[0], marker='o', color='red', markersize=8, linestyle='', label='upstream indicator', zorder=2)
    plt.plot(x, y, color='red', linewidth=.5, label='centerline', zorder=2)
    
    plt.axis('equal')
    plt.xlabel('Longitudinal direction (m)', fontsize=14)
    plt.ylabel('Latitudinal direction (m)', fontsize=14)
    plt.legend(edgecolor='black', facecolor='white', framealpha=1, fontsize=14)
    plt.gca().set_facecolor(bgcolor) # or 'xkcd:dark grass green'

    ## Subplot b: curvature
    plt.subplot(gs.new_subplotspec((6, 0), colspan=6, rowspan=2))
    plt.plot(s/WIDTH, cur_ori*WIDTH, ':', color='darkorange', linewidth=1)
    plt.plot(s/WIDTH, cur_flt*WIDTH, '--', color='dodgerblue', linewidth=1)
    plt.xlim(0, np.max(s)/WIDTH*1.2)
    if MIGRATION:
        plt.ylim(-1, 1)
        plt.gca().yaxis.set_major_formatter(ticker.FormatStrFormatter("%.1f"))
    if LAG == 0:
        plt.legend(['original', 'filtered'], edgecolor='black', facecolor='white', framealpha=1)
    else:
        plt.plot(s/WIDTH, cur_lag*WIDTH, '-', color='orangered', linewidth=1)
        plt.legend(['original', 'filtered', 'filtered & lagged'], edgecolor='black', facecolor='white', framealpha=1)
    plt.xlabel('Dimensionless streamwise distance (normalized by channel width)', fontsize=14)
    plt.ylabel('Dimensionless curvature\n(normalized by channel width)', fontsize=14)
    plt.title('(b) Curvature signal', fontsize=16)

    ## Subplot c
    plt.subplot(gs.new_subplotspec((0, 6), colspan=2, rowspan=8))
    plt.imshow(beck_bed/DEPTH, cmap='gist_earth', aspect=NUM*16/s.size,
               interpolation='bilinear', vmin=-1, vmax=1+SLOPE*np.max(s)/DEPTH)
    plt.xlabel('Transverse direction\n(normalized by channel half-width)', fontsize=14)
    plt.xticks([0, NUM, 2*NUM], ('-1', '0', '1'))
    plt.ylabel('Streamwise direction (normalized by channel length)', fontsize=14)
    plt.yticks(np.linspace(0,s.size,6), ('0', '0.2', '0.4', '0.6', '0.8', '1'))
    cbar = plt.colorbar(extend='both', fraction=0.08)
    if matplotlib.__version__[0] == '3':
        cbar.minorticks_on()
        cbar.ax.set_ylabel('Elevation (normalized by water depth)', fontsize=14)
    plt.title('(c) Synthetic bed in s-n coordinate\n(bilinearly interpolated)', fontsize=16)
    plt.tight_layout()

    if MIGRATION:
        if t == 0:
            extList = ['jpg', 'pdf']
            for ext in extList:
                plt.savefig(FNAME.rsplit('.', 1)[0] + '_pyriverbed.' + ext, dpi=300)
        if not os.path.exists(JPG_DIRS[0]):
            os.makedirs(JPG_DIRS[0])
        plt.savefig(JPG_DIRS[0] + FNAME.rsplit('.', 1)[0] + '_' + str(t) + '_pyriverbed.jpg', quality=80, dpi=150)
    else:
        extList = ['jpg', 'pdf']
        for ext in extList:
            plt.savefig(FNAME.rsplit('.', 1)[0] + '_pyriverbed.' + ext, dpi=300)
    plt.close()

    if MIGRATION:
        plt.figure(figsize=(12,9))
        lw = 1
        for centerlinex, centerliney, i in zip(centerlinexList, centerlineyList, range(len(centerlineyList))):
            if i*GPRINT == 0:
                c = 'red'
            elif i*GPRINT < TSTEPS - 1:
                c = str(1-i*GPRINT/TSTEPS)
            else:
                c = 'dodgerblue'
                lw = 10
            plt.plot(centerlinex, centerliney, color=c, linewidth=lw, label='_nolegend_')
        plt.axis('equal')
        plt.xlabel('Longitudinal direction (m)', fontsize=16)
        plt.ylabel('Latitudinal direction (m)', fontsize=16)
        plt.gca().set_facecolor(bgcolor) # or 'xkcd:dark grass green'
        plt.title('Centerline migration: T=' + str(int(t*DT/86400)) + 'd', fontsize=18)
        plt.tight_layout()
        if not os.path.exists(JPG_DIRS[1]):
            os.makedirs(JPG_DIRS[1])
        plt.savefig(JPG_DIRS[1] + FNAME.rsplit('.', 1)[0] + '_' + str(t) + '_migration.jpg', quality=80, dpi=150) 
        plt.close()
    print(' [done]')

    
