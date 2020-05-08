# coding: utf-8
"""
Main script for pyRiverBed
pyRiverBed: Generate Synthetic Riverbed Topography for Meandering Rivers
MIT License
Author: Zhi Li, Univ. of Illinois Urbana-Champaign
Contact: zhil2[at]illinois[dot]edu

"""

from __future__ import division, print_function
import sys
import os
import shutil
import numpy as np
from scipy.signal import savgol_filter
from scipy.stats import mode
from numba import jit
from tabulate import tabulate
import plotdemo as mf


def print_banner():
    """
    Print the banner of pyRiverBed.
    
    Parameters
    ----------
    None
    
    Returns
    -------
    None

    """
    print('\n\n' + '+' + '-'*60 + '+')
    print('                      - pyRiverBed -')
    print(' Generate Synthetic Riverbed Topography for Meandering Rivers')
    print('                        MIT License')
    print('         Zhi Li, Univ. of Illinois Urbana-Champaign')
    print('                  zhil2[at]illinois[dot]edu')
    print('+' + '-'*60 + '+', end='')
    print("""
                ____  _                ____           _ 
    _ __  _   _|  _ \(_)_   _____ _ __| __ )  ___  __| |
   | '_ \| | | | |_) | \ \ / / _ \ '__|  _ \ / _ \/ _` |
   | |_) | |_| |  _ <| |\ V /  __/ |  | |_) |  __/ (_| |
   | .__/ \__, |_| \_\_| \_/ \___|_|  |____/ \___|\__,_|
   |_|    |___/                                         

   MODE 1: GENERATE KINOSHITA CURVE FROM EQUATION
   MODE 2: READ YOUR OWN RIVER CENTERLINE FROM FILE
    """)


def read_steering():
    """
    Read the steering file to gather user inputs from the GUI of pyRiverBed.

    Parameters are declared as global variables.
    
    Parameters
    ----------
    None
    
    Returns
    -------
    params : tuple
        variables needed in the plotting script

    """
    print('+> Trying to read steering file...', end='')
    
    try:
        d = np.loadtxt('steering.txt', delimiter=',', skiprows=1)
        print(' [done]')
    except IOError:
        print('\nNo steering file found')
        print('Please provide steering file first\n')
        job_done()
        sys.exit(1)
        
    global MODE, NBENDS, LAMBDA, THETA0, JS, JF, WIDTH, DEPTH, SLOPE, DS, \
    NUM, INTERVAL, LAG, LAGSTR, SAVEXYZ, SAVEBOUND, SAVEMESH, FLIPSTRM, \
    FLIPTRANS, MIGRATION, UB0, C0, CF0, FR0, DT, E0, LPRINT, TSTEPS, \
    GPRINT, FPS, ZERO, JPG_DIRS, FNAME, SMOLEV, STCORR
    
    MODE = np.int(d[0])
    NBENDS = np.int(d[1])
    LAMBDA = np.float(d[2])
    THETA0 = np.float(d[3])*np.pi/180
    JS = np.float(d[4])
    JF = np.float(d[5])
    WIDTH = np.float(d[6])
    DEPTH = np.float(d[7])
    SLOPE = np.float(d[8])
    DS = np.float(d[9])
    NUM = np.int(d[10])
    INTERVAL = WIDTH/2/NUM
    LAG = np.int(d[11])
    LAGSTR = d[12]
    SAVEXYZ = np.int(d[13])
    SAVEBOUND = np.int(d[14])
    SAVEMESH = np.int(d[15])
    FLIPSTRM = np.int(d[16])
    FLIPTRANS = np.int(d[17])
    MIGRATION = np.int(d[18])
    UB0 = np.float(d[19])
    C0 = np.float(d[20])
    CF0 = np.float(d[21])
    FR0 = np.float(d[22])
    DT = np.int(d[23])
    E0 = np.float(d[24])
    LPRINT = np.int(d[25])
    TSTEPS = np.int(d[26])
    if MIGRATION == 0:
        TSTEPS = 0
    GPRINT = np.int(d[27])
    FPS = np.int(d[28])
    SMOLEV = np.int(d[29])
    STCORR = d[30]
    ZERO = 1e-8
    JPG_DIRS = ['./jpg1/', './jpg2/']

    with open('steering.txt', 'r') as f:
        lines = f.readlines()
        FNAME = lines[0].rstrip()
    if MODE == 1:
        FNAME = 'kinoshita'
        
    params = WIDTH, DEPTH, SLOPE, NUM, LAG, FNAME, \
             MIGRATION, DT, TSTEPS, GPRINT, JPG_DIRS
    
    return params


def print_para_table(s):
    """
    Print a table displaying parameters read from the steering file. 

    Require 'tabulate' library.
    
    Parameters
    ----------
    s : ndarray
        streamwise distance

    Returns
    -------
    None

    """
    if MODE == 1:
        t = [['Parameter', 'Value', 'Unit'],
             ['Number of bends', NBENDS, '/'], 
             ['Width', WIDTH, 'm'],
             ['Depth', DEPTH, 'm'],
             ['Length', LAMBDA*(NBENDS+1), 'm'],
             ['Arc wavelength', LAMBDA, 'm'],
             ['Slope', SLOPE, '/'],
             ['Streamwise resolution', DS, 'm'],
             ['Transverse resolution', np.around(INTERVAL, decimals=4), 'm'],
             ['Streamwise # of pts', s.size + 2*np.int(LAMBDA/2/DS), '/'],
             ['Transverse # of pts', NUM*2+1, '/']]
    elif MODE == 2:
        if FNAME[0].islower():
            f = FNAME[0].upper() + FNAME[1:]
        else:
            f = FNAME
        t = [['Parameter', 'Value', 'Unit'],
             ['River name', f.rsplit('.', 1)[0], '/'],
             ['Width', WIDTH, 'm'],
             ['Depth', DEPTH, 'm'],
             ['Length', np.round(s[-1], decimals=2), 'm'],
             ['Slope', SLOPE, '/'],
             ['Streamwise resolution', np.round(np.mean(np.diff(s)), decimals=2), 'm'],
             ['Transverse resolution', np.round(INTERVAL, decimals=2), 'm'],
             ['Streamwise # of pts', s.size, '/'],
             ['Transverse # of pts', NUM*2+1, '/']]
    print(tabulate(t, tablefmt='psql', stralign='right', headers='firstrow'))


def print_resamp_table(mean1, median1, mode1, mean2, median2, mode2):
    """
    Print a table displaying mean, median and mode of centerline grid size
    before and after resampling.

    Require 'tabulate' library.
    
    Parameters
    ----------
    mean1, median1, mode1 : ndarray
        before resampling
    mean2, median2, mode2 : ndarray
        after resampling

    Returns
    -------
    None

    """
    t = [['Streamwise\nresolution', 'Before     '
          +'After\nresampling --> resampling', '\nUnit'],
         ['Mean', str(mean1) + ' --> ' + str(mean2), 'm'],
         ['Median', str(median1) + ' --> ' + str(median2), 'm'],
         ['Mode', str(mode1) + ' --> ' + str(mode2), 'm']]
    print(tabulate(t, tablefmt='psql', stralign='center', headers='firstrow'))
    
    
def print_eqn():
    """
    Print Kinoshita Curve equation.

    Only work for Mode 1.
    
    Parameters
    ----------
    None
    
    Returns
    -------
    None

    """
    if sys.stdout.encoding.lower().startswith('utf'):
        if JS != 0 and JF != 0:
            print('Eqn: \u03B8=' + np.str(np.around(THETA0, decimals=6)) +
                  '*sin(2\u03C0s/' + np.str(np.around(LAMBDA, decimals=6)) +
                  ')\n      +' + np.str(np.around(THETA0**3, decimals=6)) +
                  '*[' + np.str(np.around(JS, decimals=6)) + '*cos(6\u03C0s/' +
                  np.str(np.around(LAMBDA, decimals=6)) +
                  ')-' + np.str(np.around(JF, decimals=6)) + '*sin(6\u03C0s/' +
                  np.str(np.around(LAMBDA, decimals=6)) + ')]')
        elif JS == 0 and JF != 0:
            print('Eqn: \u03B8=' + np.str(np.around(THETA0, decimals=6)) +
                  '*sin(2\u03C0s/' + np.str(np.around(LAMBDA, decimals=6)) +
                  ')+' + np.str(np.around(THETA0**3, decimals=6)) + '*[' + 
                  '-' + np.str(np.around(JF, decimals=6)) + '*sin(6\u03C0s/' +
                  np.str(np.around(LAMBDA, decimals=6)) + ')]')
        elif JS != 0 and JF == 0:
            print('Eqn: \u03B8=' + np.str(np.around(THETA0, decimals=6)) +
                  '*sin(2\u03C0s/' + np.str(np.around(LAMBDA, decimals=6)) +
                  ')+' + np.str(np.around(THETA0**3, decimals=6)) +
                  '*[' + np.str(np.around(JS, decimals=6)) + '*cos(6\u03C0s/' +
                  np.str(np.around(LAMBDA, decimals=6)) + ')]')
        elif JS == 0 and JF == 0:
            print('Eqn: \u03B8=' + np.str(np.around(THETA0, decimals=6)) +
                  '*sin(2\u03C0s/' + np.str(np.around(LAMBDA, decimals=6)) + ')')
    else:
        if JS != 0 and JF != 0:
            print('Eqn: THETA=' + np.str(np.around(THETA0, decimals=6)) +
                  '*sin(2PI/' + np.str(np.around(LAMBDA, decimals=6)) +
                  ')\n          +' + np.str(np.around(THETA0**3, decimals=6)) +
                  '*[' + np.str(np.around(JS, decimals=6)) + '*cos(6PI/' +
                  np.str(np.around(LAMBDA, decimals=6)) +
                  ')-' + np.str(np.around(JF, decimals=6)) + '*sin(6PI/' +
                  np.str(np.around(LAMBDA, decimals=6)) + ')]')
        elif JS == 0 and JF != 0:
            print('Eqn: THETA=' + np.str(np.around(THETA0, decimals=6)) +
                  '*sin(2PI/' + np.str(np.around(LAMBDA, decimals=6)) +
                  ')+' + np.str(np.around(THETA0**3, decimals=6)) + '*[' + 
                  '-' + np.str(np.around(JF, decimals=6)) + '*sin(6PI/' +
                  np.str(np.around(LAMBDA, decimals=6)) + ')]')
        elif JS != 0 and JF == 0:
            print('Eqn: THETA=' + np.str(np.around(THETA0, decimals=6)) +
                  '*sin(2PI/' + np.str(np.around(LAMBDA, decimals=6)) +
                  ')+' + np.str(np.around(THETA0**3, decimals=6)) +
                  '*[' + np.str(np.around(JS, decimals=6)) + '*cos(6PI/' +
                  np.str(np.around(LAMBDA, decimals=6)) + ')]')
        elif JS == 0 and JF == 0:
            print('Eqn: THETA=' + np.str(np.around(THETA0, decimals=6)) +
                  '*sin(2PI/' + np.str(np.around(LAMBDA, decimals=6)) + ')')


def build_kinoshita():
    """
    Build Kinoshita Curve (non-computational part).

    Only work for Mode 1.
    
    Parameters
    ----------
    None
    
    Returns
    -------
    s : ndarray
        streamwise distance
    x : ndarray
        x coordinate of Kinoshita Curve
    y : ndarray
        y coordinate of Kinoshita Curve
    cur : ndarray
        curvature
    theta : ndarray
        angular ampitude

    """
    if MODE != 1:
        return [], [], [], [], []
    print('MODE 1: GENERATE KINOSHITA CURVE FROM EQUATION is selected')
    print('Kinoshita Curve parameters are read from steering file:')
    print_eqn()
    s = np.linspace(0, NBENDS*LAMBDA, np.int(NBENDS*LAMBDA/DS) + 1)
    print_para_table(s)
    print('+> Calculating Kinoshita Curve...', end='')
    s, x, y, cur, theta = compute_kinoshita(s)
    print(' [done]')
    return s, x, y, cur, theta


@jit(nopython=True)
def compute_kinoshita(s):
    """
    Build Kinoshita Curve (computational part).

    Numba nopyton mode is on.

    Only work for Mode 1.
    
    Parameters
    ----------
    s : ndarray 
        streamwise distance
    
    Returns
    -------
    s : ndarray 
        streamwise distance
    x : ndarray
        x coordinate of Kinoshita Curve
    y : ndarray
        y coordinate of Kinoshita Curve
    cur : ndarray
        curvature
    theta : ndarray
        angular ampitude

    """
    length = np.int(NBENDS*LAMBDA/DS) + 1
    x = np.zeros(length)
    y = np.zeros(length)
    cur = np.zeros(length+1)
    theta = THETA0*np.sin(2*np.pi*s/LAMBDA) \
            + THETA0**3*(JS*np.cos(6*np.pi*s/LAMBDA) \
            - JF*np.sin(6*np.pi*s/LAMBDA))
    theta[np.abs(theta)<ZERO] = 0
    for i in range(length):
        cossum, sinsum = 0, 0
        for j in range(i):
            cossum += DS*np.cos(theta[j])
            sinsum += DS*np.sin(theta[j])
        x[i] = 0 if np.abs(cossum) < ZERO else cossum
        y[i] = 0 if np.abs(sinsum) < ZERO else sinsum
    x = np.concatenate((x, np.array([x[-1]+x[1]-x[0]])))
    y = np.concatenate((y, np.array([y[-1]+y[1]-y[0]])))
    s = np.concatenate((s, np.array([s[-1]+DS])))
    theta = np.concatenate((theta, np.array([theta[-1]])))
    if FLIPSTRM:
        x = x[::-1]
        y = y[::-1]
        theta = np.concatenate((theta[::-1][1:], np.array([theta[0]])))
    for i in range(1, length):
        cur[i] = (theta[i]-theta[i-1])/DS
        cur[i] = 0 if np.abs(cur[i]) < ZERO else cur[i]
    cur[0], cur[-1] = cur[-2], cur[1]
    return s, x, y, cur, theta


def read_centerline(s, x, y, cur, theta):
    """
    Read river centerline coordinates from user-prepared centerline file.

    Centerline is then resampled to prevent ununiform spacing.
    
    Only work for Mode 2.
    
    Parameters
    ----------
    s : ndarray 
        streamwise distance (empty)
    x : ndarray
        x coordinate of river centerline read from file (empty)
    y : ndarray
        y coordinate of river centerline read from file (empty)
    cur : ndarray
        curvature (empty)
    theta : ndarray
        angular ampitude (empty)
    
    Returns
    -------
    s : ndarray 
        streamwise distance
    x : ndarray
        x coordinate of river centerline read from file
    y : ndarray
        y coordinate of river centerline read from file
    cur : ndarray
        curvature
    theta : ndarray
        angular ampitude

    """
    if MODE == 2:
        print('MODE 2: READ YOUR OWN RIVER CENTERLINE FROM FILE is selected')
        try:
            centerlinexy = np.loadtxt(FNAME)
        except IOError:
            print('\'' + FNAME + '\' not found')
            print('Please place \'' + FNAME + '\' in working directory\n')
            job_done()
            sys.exit(1)
    else:
        return s, x, y, cur, theta
    x = centerlinexy[:, 0]
    y = centerlinexy[:, 1]
    if FLIPSTRM:
        x = x[::-1]
        y = y[::-1]
#    if np.mean(np.abs(x)) > 1e6 or np.mean(np.abs(y)) > 1e6:
#        print('!!! centerline X/Y too large, forced to shift toward (0, 0) !!!')
#        print('shifting vector: ('+str(-np.mean(x))+', '+str(-np.mean(y))+')')
#        x -= np.mean(x)
#        y -= np.mean(y)
    length = x.size
    s = np.zeros(length)
    for j in range(1, x.size):
        s[j] = s[j-1] + np.sqrt((x[j]-x[j-1])**2 + (y[j]-y[j-1])**2)
    mean1 = np.around(np.mean(np.diff(s)), decimals=2)
    median1 = np.around(np.median(np.diff(s)), decimals=2)
    mode1 = np.around(mode(np.diff(s))[0][0], decimals=2)
    print('+> Resampling centerline & Calculating curvature...', end='')
    s, x, y, cur, theta = resample_centerline(s, x, y)
    print(' [done]')
    mean2 = np.around(np.mean(np.diff(s)), decimals=2)
    median2 = np.around(np.median(np.diff(s)), decimals=2)
    mode2 = np.around(mode(np.diff(s))[0][0], decimals=2)
    print_resamp_table(mean1, median1, mode1, mean2, median2, mode2)
    print_para_table(s)
    return s, x, y, cur, theta


def extend_centerline(s, x, y, cur, theta):
    """
    Extend centerline to have straight channels at both ends.

    Parameters
    ----------
    s : ndarray 
        streamwise distance
    x : ndarray
        x coordinate of centerline
    y : ndarray
        y coordinate of centerline
    cur : ndarray
        curvature
    theta : ndarray
        angular ampitude
    
    Returns
    -------
    s : ndarray 
        extended streamwise distance
    x : ndarray
        extended x coordinate of centerline
    y : ndarray
        extended y coordinate of centerline
    cur : ndarray
        extended curvature
    theta : ndarray
        extended angular ampitude

    """
    print('+> Extending centerline to have straight channels at both ends...', end='')
    if MODE == 1:
        extlength = LAMBDA/10
        d = DS     
    elif MODE == 2:
        extlength = WIDTH
        d = INTERVAL
    num = np.int(extlength/d)
        
    coshead = (x[1] - x[0])/d
    sinhead = (y[1] - y[0])/d
    headx = np.linspace(x[0] - extlength*coshead, x[0] - d*coshead, num)
    heady = np.linspace(y[0] - extlength*sinhead, y[0] - d*sinhead, num)

    costail = (x[-1] - x[-2])/d
    sintail = (y[-1] - y[-2])/d
    tailx = np.linspace(x[-1] + d*costail, x[-1] + extlength*costail, num)
    taily = np.linspace(y[-1] + d*sintail, y[-1] + extlength*sintail, num)

    x = np.concatenate((headx, x, tailx), axis=0)
    y = np.concatenate((heady, y, taily), axis=0)
    s, x, y = smooth_centerline(x, y)
    s, x, y, cur, theta = resample_centerline(s, x, y)
    print(' [done]')
    return s, x, y, cur, theta


@jit(nopython=True)
def resample_centerline(s, x, y):
    """
    Resample centerline coordinates.
    
    Parameters
    ----------
    s : ndarray 
        streamwise distance
    x : ndarray
        x coordinate of river centerline read from file
    y : ndarray
        y coordinate of river centerline read from file
    
    Returns
    -------
    sn : ndarray 
        streamwise distance after resampling
    xn : ndarray
        x coordinate of river centerline after resampling
    yn : ndarray
        y coordinate of river centerline after resampling

    """
    N = int(s[-1]//INTERVAL) + 1
    ds = s[-1]/N
    length = s.size
    xn = np.zeros(N+1)
    yn = np.zeros(N+1)
    xn[0], yn[0], xn[-1], yn[-1] = x[0], y[0], x[-1], y[-1]
    for i in range(1, N):
        curr_loc = i*ds
        for j in range(length-1):
            if s[j] < curr_loc and s[j+1] > curr_loc:
                r = (curr_loc-s[j])/(s[j+1]-s[j])
                break
        xn[i] = x[j] + (x[j+1]-x[j])*r
        yn[i] = y[j] + (y[j+1]-y[j])*r
    sn = np.zeros(xn.size)
    for i in range(1,xn.size):
        sn[i] = sn[i-1] + np.sqrt((xn[i]-xn[i-1])**2 + (yn[i]-yn[i-1])**2)
    curn, thetan = tan2curv(sn, xn, yn)
    return sn, xn, yn, curn, thetan


def smooth_centerline(x, y):
    """
    Smooth centerline using Savitzky–Golay filter 
    (5-point quadratic polynomial method).
    
    Number of passes is equal to:
        SMOOTHINGLEVEL (SMOOTHINGLEVEL < 39)
        int(1.1^SMOOTHINGLEVEL) (SMOOTHINGLEVEL >= 39)
    
    Parameters
    ----------
    x : ndarray
        x coordinate of river centerline before smoothing
    y : ndarray
        y coordinate of river centerline before smoothing
    
    Returns
    -------
    s : ndarray 
        streamwise distance after smoothing
    x : ndarray
        x coordinate of river centerline after smoothing
    y : ndarray
        y coordinate of river centerline after smoothing
    
    """
    n = SMOLEV if SMOLEV < 39 else np.int(np.around(1.1**SMOLEV,decimals=0))
    xa, xb, ya, yb = x[0], x[-1], y[0], y[-1]
    for i in range(n):
        x = savgol_filter(x, 5, 2, mode='nearest')
        y = savgol_filter(y, 5, 2, mode='nearest')
        x[0], x[-1], y[0], y[-1] = xa, xb, ya, yb
    s = np.zeros(x.size)
    for i in range(1, x.size):
        s[i] = s[i-1] + np.sqrt((x[i]-x[i-1])**2 + (y[i]-y[i-1])**2)
    return s, x, y


def filter_curvature(cur, t):
    """
    Filter the curvature signal by the Savitzky–Golay filter 
    (5-point quadratic polynomial method) and 5-point moving
    average.
    
    Parameters
    ----------
    cur : ndarray
        curvature
    t : integer
        current # of time step
    
    Returns
    -------
    cur : ndarray
        filtered curvature

    """
    if np.mod(t, LPRINT) == 0:
        print('+> Filtering curvature (5-pt Savitzky–Golay + 5-pt Moving Average)...', end='')
    cur = savgol_filter(savgol_filter(cur, 5, 2, mode='nearest'), 5, 1, mode='nearest')
    if np.mod(t, LPRINT) == 0:
        print(' [done]')
    return cur


def lag(s, cur, t):
    """
    Impose a phase lag to the curvature signal by replacing the local 
    curvature with the upstream-wise moving averaged curvature.
    
    Parameters
    ----------
    s : ndarray 
        streamwise distance
    cur : ndarray
        curvature
    t : integer
        current # of time step
    
    Returns
    -------
    s : ndarray 
        streamwise distance with considering phase lag
    cur : ndarray
        curvature with considering phase lag

    """
    if LAG == 0:
        return cur
    else:
        if MODE == 1:
            num = np.int(WIDTH*LAGSTR/DS)
        elif MODE == 2:
            num = np.int(WIDTH*LAGSTR/np.mean(np.diff(s)))
    if np.mod(t, LPRINT) == 0:
        print('+> Adding phase lag to curvature...', end='')
    cur = compute_lag(cur, num)
    if np.mod(t, LPRINT) == 0:
        print(' [done]')
    return cur


@jit(nopython=True)
def compute_lag(cur, num):
    """
    Compute phase lag.
    
    Numba nopyton mode is on.
    
    Parameters
    ----------
    cur : ndarray
        curvature
    num : integer
        moving average windows size
    
    Returns
    -------
    cur : ndarray
        curvature with considering phase lag
    
    """
    length = cur.size
    cur0 = np.copy(cur)
    for i in range(2, length):
        M = i if i < num else num
        c = 0
        for j in range(M):
            c += (2/M-j*2/M/(M-1))*cur0[i-j]
        cur[i] = c
    return cur


@jit(nopython=True)
def tan2curv(s, x, y):
    """
    Compute curvature using 'arctan2' method.
    
    Parameters
    ----------
    s : ndarray 
        streamwise distance
    x : ndarray
        x coordinate of river centerline
    y : ndarray
        y coordinate of river centerline
    
    Returns
    -------
    cur : ndarray
        curvature
    forw : ndarray
        angular ampitude

    """
    length = x.size
    cur = np.zeros(length)
    forw = np.zeros(length)
    back = np.zeros(length)
    for i in range(1, length-1):
        forw[i] = np.arctan2(y[i+1]-y[i], x[i+1]-x[i])
        back[i] = np.arctan2(y[i]-y[i-1], x[i]-x[i-1])
        angle_atan2 = forw[i] - back[i]
        cur[i] = angle_atan2/(s[i+1]-s[i-1])*2
        if np.abs(cur[i]) < ZERO:
            cur[i] = 0
    for i in range(1, length-1):
        ave = (cur[i-1]+cur[i+1])/2
        if np.abs(cur[i]-ave) > 5*np.abs(cur[i-1]-cur[i+1]):
            cur[i] = ave
    forw[0], forw[-1] = back[1], forw[-2]
    return cur, forw


@jit(nopython=True)
def coscurv(s, x, y):
    """
    Compute curvature using 'law of cosine' method.
    
    Parameters
    ----------
    s : ndarray 
        streamwise distance
    x : ndarray
        x coordinate of river centerline
    y : ndarray
        y coordinate of river centerline
    
    Returns
    -------
    cur : ndarray
        curvature

    """
    length = x.size
    cur = np.zeros(length)
    for i in range(1, length-1):
        a = np.array([x[i+1]-x[i], y[i+1]-y[i]])
        b = np.array([x[i]-x[i-1], y[i]-y[i-1]])
        c = np.array([1, 0])
        flag = 1
        if flag == 1 and a[1] < 0:
            flag = -1
        elif flag == -1 and a[1] <= 0:
            flag = 1
        angle_cos = flag \
            *(np.arccos(np.vdot(a, c)/np.linalg.norm(a)/np.linalg.norm(c)) \
            - np.arccos(np.vdot(b, c)/np.linalg.norm(b)/np.linalg.norm(c)))
        cur[i] = angle_cos/(s[i+1]-s[i-1])*2
        if np.abs(cur[i]) < ZERO:
            cur[i] = 0
    for i in range(1, length-1):
        ave = (cur[i-1]+cur[i+1])/2
        if np.abs(cur[i]-ave) > 5*np.abs(cur[i-1]-cur[i+1]):
            cur[i] = ave
    return cur


@jit(nopython=True)
def threeptscurv(x, y):
    """
    Compute curvature using 'triangle's circumscribed circle ' method.
    
    Parameters
    ----------
    s : ndarray 
        streamwise distance
    x : ndarray
        x coordinate of river centerline
    y : ndarray
        y coordinate of river centerline
    
    Returns
    -------
    cur : ndarray
        curvature

    """
    length = x.size
    R = np.zeros(length)
    cur = np.zeros(length)
    for i in range(1, length-1):
        a = np.sqrt((x[i+1]-x[i])**2 + (y[i+1]-y[i])**2)
        b = np.sqrt((x[i+1]-x[i-1])**2 + (y[i+1]-y[i-1])**2)
        c = np.sqrt((x[i]-x[i-1])**2 + (y[i]-y[i-1])**2)
        p = (a+b+c)/2
        R[i] = a*b*c/4/np.sqrt(p*(p-a)*(p-b)*(p-c))
        cur[i] = 1/R[i]
        if R[i] > 1/ZERO or np.isnan(R[i]):
            cur[i] = 0
    return cur


def build_beck(cur, s, t):
    """
    Build synthetic bed topography (non-computational part).
    
    Parameters
    ----------
    cur : ndarray
        curvature
    s : ndarray 
        streamwise distance
    t : integer
        current # of time step

    Returns
    -------
    beck_bed : ndarray
        synthetic bed topography data  

    """
    if np.mod(t, LPRINT) == 0:
        print('+> Calculating synthetic riverbed topography...', end='')
    beck_bed = compute_beck(cur, s)
    beck_bed[np.abs(beck_bed)<ZERO] = 0
    if np.mod(t, LPRINT) == 0:
        print(' [done]')
    return beck_bed


@jit(nopython=True)
def compute_beck(cur, s):
    """
    Build synthetic bed topography using Beck-1988 formula (computational part).

    Numba nopyton mode is on.
    
    Parameters
    ----------
    cur : ndarray
        curvature
    s : ndarray 
        streamwise distance
    
    Returns
    -------
    beck_bed : ndarray
        synthetic bed topography data  

    """
    halfwidth = WIDTH/2
    A = 3.8*(1+halfwidth/6.96/DEPTH*np.exp(-6.96*DEPTH/halfwidth))
    st = -A*DEPTH*cur*STCORR
    length = cur.size
    hc = np.ones(length)
    for i in range(length):
        if np.abs(st[i]) < ZERO:
            st[i] = ZERO
        hc[i] = (4*halfwidth*DEPTH*np.abs(st[i])-st[i]**2*halfwidth**2) \
                /(2*halfwidth*np.abs(st[i])+2*DEPTH-2*DEPTH \
                *np.exp(-np.abs(st[i])*halfwidth/DEPTH))
    beck_bed = np.zeros((length, 2*NUM+1))
    slope = (np.max(s)-s)*SLOPE
    for j in range(2*NUM+1):
        if j == NUM:
            beck_bed[:, NUM] = hc - slope
            continue
        n = -halfwidth + j*INTERVAL
        for i in range(length):
            beck_bed[i, j] = (1 - hc[i]/st[i]/n)*np.maximum(-st[i]*n, 0) \
                + hc[i]/st[i]/n*np.exp(-st[i]*n/DEPTH)*np.maximum(st[i]*n, 0) \
                - slope[i]
    beck_bed = DEPTH - beck_bed
    if FLIPTRANS:
        beck_bed = beck_bed.T[::-1].T
    return beck_bed


@jit(nopython=True)
def offset(x, y, L):
    """
    Compute left and right offset polylines of centerline 
    with an offset distance of L.
    
    Parameters
    ----------
    x : ndarray
        x coordinate of river centerline
    y : ndarray
        y coordinate of river centerline
    L : float
        offset distance
    
    Returns
    -------
    offsetx : ndarray
        x coordinates of left and right offset polylines
    offsety : ndarray
        y coordinates of left and right offset polylines  

    """
    length = x.size
    Cx = np.zeros(length)
    Cy = np.zeros(length)
    avek = np.zeros(length)
    theta = np.zeros(length)
    R = np.zeros(length)
    offsetx = np.zeros((length, 2))
    offsety = np.zeros((length, 2))
    for i in range(1, length-1):
        a = np.sqrt((x[i+1]-x[i])**2 + (y[i+1]-y[i])**2)
        b = np.sqrt((x[i+1]-x[i-1])**2 + (y[i+1]-y[i-1])**2)
        c = np.sqrt((x[i]-x[i-1])**2 + (y[i]-y[i-1])**2)
        p = a + b + c
        area = 0.5*np.abs((x[i-1]*y[i] - x[i]*y[i-1] + x[i]*y[i+1]
            - x[i+1]*y[i] + x[i+1]*y[i-1] - x[i-1]*y[i+1]))
        Cx[i] = (x[i-1]*a + x[i]*b + x[i+1]*c)/p
        Cy[i] = (y[i-1]*a + y[i]*b + y[i+1]*c)/p
        if area/WIDTH**2 < ZERO and np.abs(y[i] - y[i-1]) < ZERO:
            avek[i] = 1/ZERO
        elif area/WIDTH**2 < ZERO and np.abs(x[i] - x[i-1]) < ZERO:
            avek[i] = 0
        elif area/WIDTH**2 < ZERO:
            avek[i] = -(x[i] - x[i-1])/(y[i] - y[i-1])
        elif np.abs(x[i] - Cx[i]) < ZERO:
            avek[i] = 1/ZERO
        elif np.abs(y[i] - Cy[i]) < ZERO:
            avek[i] = 0
        else:
            avek[i] = (y[i] - Cy[i])/(x[i] - Cx[i])
        m = np.array([x[i+1] - x[i], y[i+1] - y[i]])
        n = np.array([x[i-1] - x[i], y[i-1] - y[i]])
        if np.abs(np.vdot(m, n)/np.linalg.norm(m)/np.linalg.norm(n)+1) < ZERO:
            theta[i] = np.pi
        elif np.abs(np.vdot(m, n)/np.linalg.norm(m)/np.linalg.norm(n)-1) < ZERO:
            theta[i] = 0
        elif np.abs(np.vdot(m, n)/np.linalg.norm(m)/np.linalg.norm(n)) < ZERO:
            theta[i] = np.pi/2
        else:
            theta[i] = np.arccos(np.vdot(m, n)/np.linalg.norm(m)/np.linalg.norm(n))
        R[i] = L/np.sin(theta[i]/2)
        if np.isnan(R[i]) or np.isinf(R[i]):
            R[i] = L
        if np.abs(avek[i]) < ZERO:
            offsetx[i, 0] = x[i] + R[i]
            offsetx[i, 1] = x[i] - R[i]
            offsety[i, 0] = y[i]
            offsety[i, 1] = y[i]
        elif np.isinf(np.abs(avek[i])) or np.abs(avek[i]) > 1/ZERO:
            offsetx[i, 0] = x[i]
            offsetx[i, 1] = x[i]
            offsety[i, 0] = y[i] + R[i]
            offsety[i, 1] = y[i] - R[i]
        else:
            offsetx[i, 0] = R[i]/np.sqrt(avek[i]**2 + 1) + x[i]
            offsetx[i, 1] = -R[i]/np.sqrt(avek[i]**2 + 1) + x[i]
            offsety[i, 0] = avek[i]*(offsetx[i, 0] - x[i]) + y[i]
            offsety[i, 1] = avek[i]*(offsetx[i, 1] - x[i]) + y[i]
        if np.sqrt((offsetx[i, 0]-offsetx[i-1, 0])**2+(offsety[i, 0]-offsety[i-1, 0])**2) > R[i]*2:
            offsetx[i, 0], offsetx[i, 1] = offsetx[i, 1], offsetx[i, 0]
            offsety[i, 0], offsety[i, 1] = offsety[i, 1], offsety[i, 0]
    coshead = (x[1] - x[0])/np.sqrt((x[1]-x[0])**2+(y[1]-y[0])**2)
    sinhead = (y[1] - y[0])/np.sqrt((x[1]-x[0])**2+(y[1]-y[0])**2)
    costail = (x[-1] - x[-2])/np.sqrt((x[-1]-x[-2])**2+(y[-1]-y[-2])**2)
    sintail = (y[-1] - y[-2])/np.sqrt((x[-1]-x[-2])**2+(y[-1]-y[-2])**2)
    offsetx[0, 0] = offsetx[1, 0] - np.sqrt((x[1]-x[0])**2+(y[1]-y[0])**2)*coshead
    offsetx[0, 1] = offsetx[1, 1] - np.sqrt((x[1]-x[0])**2+(y[1]-y[0])**2)*coshead
    offsety[0, 0] = offsety[1, 0] - np.sqrt((x[1]-x[0])**2+(y[1]-y[0])**2)*sinhead
    offsety[0, 1] = offsety[1, 1] - np.sqrt((x[1]-x[0])**2+(y[1]-y[0])**2)*sinhead
    offsetx[-1, 0] = offsetx[-2, 0] + np.sqrt((x[-1]-x[-2])**2+(y[-1]-y[-2])**2)*costail
    offsetx[-1, 1] = offsetx[-2, 1] + np.sqrt((x[-1]-x[-2])**2+(y[-1]-y[-2])**2)*costail
    offsety[-1, 0] = offsety[-2, 0] + np.sqrt((x[-1]-x[-2])**2+(y[-1]-y[-2])**2)*sintail
    offsety[-1, 1] = offsety[-2, 1] + np.sqrt((x[-1]-x[-2])**2+(y[-1]-y[-2])**2)*sintail
    return offsetx, offsety


def offset_all(x, y, beck_bed, t):
    """
    Compute the offset polylines of centerline.

    Merge coordinates data (x & y information) with bed topography 
    data (z information) to form a point cloud dataset in 3-column xyz format.
    
    Parameters
    ----------
    x : ndarray
        x coordinate of river centerline
    y : ndarray
        y coordinate of river centerline
    beck_bed : ndarray
        synthetic bed topography data
    t : integer
        current # of time step
    
    Returns
    -------
    allxyz : ndarray
        point cloud in 3-column xyz format

    """
    length = x.size
    xyz1 = np.zeros((length, 3))
    xyz2 = np.zeros((length, 3))
    xyz1[:, 0] = np.copy(x)
    xyz1[:, 1] = np.copy(y)
    xyz1[:, 2] = np.copy(beck_bed[:, NUM])
    allxyz = np.copy(xyz1)
    offsetx = np.zeros((length, 2))
    offsety = np.zeros((length, 2))
    for i in range(NUM-1, -1, -1):
        """Offset distance L is looping from INTERVAL to B."""
        if np.mod(t, LPRINT) == 0:
            if i == NUM - 1:
                extr = '...(innermost)'
            elif i == 0:
                extr = '...(outermost)'
            else:
                extr = '...'
            print('+> Offsetting Polyline #'
                  + np.str(i+1) + ' & #' + np.str(2*NUM+1-i) + extr, end='')
        offsetx, offsety = offset(x, y, WIDTH/2-i*INTERVAL)
        if i == 0 and SAVEBOUND and t == 0:
            t1 = np.copy(offsetx)
            t2 = np.copy(offsetx)
            t1[:,0] = np.copy(offsetx[:, 0])
            t1[:,1] = np.copy(offsety[:, 0])
            t2[:,0] = np.copy(offsetx[:, 1])
            t2[:,1] = np.copy(offsety[:, 1])
            t3 = np.concatenate((t1, t2[::-1], np.array([t1[0, :]])), axis=0)
            np.savetxt(FNAME.rsplit('.', 1)[0] + '_boundary.i2s', t3, fmt='%.6e')
        xyz1[:, 0] = offsetx[:, 0]
        xyz1[:, 1] = offsety[:, 0]
        xyz1[:, 2] = beck_bed[:, -1-i]
        xyz2[:, 0] = offsetx[:, 1]
        xyz2[:, 1] = offsety[:, 1]
        xyz2[:, 2] = beck_bed[:, i]
        allxyz = np.concatenate((allxyz, xyz1, xyz2), axis=0)
        if np.mod(t, LPRINT) == 0:
            print(' [done]')
        if i == 0 and np.mod(t, LPRINT) == 0:
            print('   * Note: Polyline #' + np.str(NUM + 1) + ' is centerline')
    return allxyz


def write_xyz_file(allxyz):
    """
    Write the point cloud of riverbed topography data.
    
    Parameters
    ----------
    allxyz : ndarray
        point cloud in 3-column xyz format

    Returns
    -------
    None
    
    Notes
    -----
    xyz : can be loaded in Blue Kenue.

    """
    if SAVEXYZ:
        print('+> Saving riverbed topography file...', end='')
        if MODE == 1:
            np.savetxt('kinoshita_topo.xyz', allxyz, fmt='%.6e')
        elif MODE == 2:
            np.savetxt(FNAME.rsplit('.', 1)[0] + '_topo.xyz', allxyz, fmt='%.6e')
        print(' [done]')


def write_mesh_file(allxyz, beck_bed):
    """
    Build and write the finite element mesh (non-compuational).
    
    Parameters
    ----------
    allxyz : ndarray
        point cloud in 3-column xyz format
    beck_bed : ndarray
        synthetic bed topography data
    
    Returns
    -------
    None
    
    Notes
    -----
    Two mesh formats:
    t3s : can be loaded in Blue Kenue and used to built a Selafin 
      file for TELEMAC modeling.
    dat : can be loaded in Tecplot for visualization.
    
    Two boundary condition formats:
    cli : can be used directly for TELEMAC modeling.
    bc2 : can be loaded in BlueKenue to check/modify BC info and metadata.

    """
    if SAVEMESH:
        print('+> Saving finite element mesh files...', end='')
        fname = FNAME.rsplit('.', 1)[0]
        ncol = beck_bed[0,:].size
        nrow = beck_bed[:,0].size
        nele = (nrow-1)*(ncol-1)*2
        d = compute_mesh(nrow, ncol, nele)
        h = ':NodeCount ' + str(allxyz[:,0].size) + '\n:ElementCount ' \
            + str(nele) + '\n#\n:EndHeader\n'
        with open(fname + '_mesh.t3s', 'w') as f: 
            f.write(h)
        with open(fname + '_mesh.t3s', 'a') as f:
            np.savetxt(f, allxyz, fmt='%.6e')
            np.savetxt(f, d, fmt='%d')
            f.write('\n\n')
        h = 'TITLE = \"' + fname \
            + '_mesh\"\nVARIABLES = \"X\", \"Y\", \"' + fname \
            + '_mesh\"\nZONE NODES=' + str(allxyz[:,0].size) + ', ELEMENTS=' \
            + str(nele) + ', DATAPACKING=POINT, ZONETYPE=FETRIANGLE\n'
        with open(fname + '_mesh.dat', 'w') as f: 
            f.write(h)
        with open(fname + '_mesh.dat', 'a') as f:
            np.savetxt(f, allxyz, fmt='%.6e')
            np.savetxt(f, d, fmt='%d')
            f.write('\n\n')
        inlet = np.zeros((ncol,), dtype=int)
        outlet = np.zeros((ncol,), dtype=int)
        for i in range(ncol):
            inlet[i] = 1 + i*nrow
            outlet[i] = (1 + i)*nrow
        left = np.zeros((nrow-2,), dtype=int)
        right = np.zeros((nrow-2,), dtype=int)
        for i in range(1, nrow-1):
            left[i-1] = (ncol-2)*nrow + i + 1
            right[i-1] = (ncol-1)*nrow + i + 1
        cli = np.zeros((2*(nrow+ncol-2), 13))
        cli[:,:2] = 2
        cli[:,7] = 2
        cli[:,11] = np.concatenate((inlet, outlet, left, right))
        cli[:,12] = np.arange(2*(nrow+ncol-2)) + 1
        cli[:ncol,0] = 4
        cli[:ncol,1] = 5
        cli[:ncol,2] = 5
        cli[:ncol,7] = 4
        cli[ncol:2*ncol,0] = 5
        cli[ncol:2*ncol,1] = 4
        cli[ncol:2*ncol,2] = 4
        cli[ncol:2*ncol,7] = 4
        np.savetxt(fname + '_BC_tmp.cli', cli, fmt='%d')
        with open(fname + '_BC.cli', 'w') as out_f:
            with open(fname + '_BC_tmp.cli', 'r') as in_f:
                for i, line in enumerate(in_f):
                    if i < ncol:
                        s = ' #Inlet'
                    elif i >= ncol and i < 2*ncol:
                        s = ' #Outlet'
                    else:
                        s = ' #'
                    out_f.write(line.rstrip('\n') + s + '\n')
            out_f.write('\n')
        os.remove(fname + '_BC_tmp.cli')
        h = ':FileType bc2  ASCII  EnSim 1.0' \
            + '\n:NodeCount ' + str(allxyz[:,0].size) \
            + '\n:ElementCount ' + str(nele) \
            + '\n:ElementType  T3' \
            + '\n:BoundarySegmentCount 2' \
            + '\n# id  code  sectionCount startNode1 endNode1 startNode2 endNode2 tracerCode name' \
            + '\n:BoundarySegment 1  455  1 1 ' + str(ncol) + ' 1 1  4  \"Inlet\"' \
            + '\n:BoundarySegment 2  544  1 ' + str(ncol+1) + ' ' + str(2*ncol) + ' 1 1  4  \"Outlet\"' \
            + '\n:ShorelineCount 1' \
            + '\n:ShorelineNodeCount ' + str(2*(nrow+ncol-2)) \
            + '\n:EndHeader' \
            + '\n:BeginNodes ' + str(allxyz[:,0].size) + '\n'
        with open(fname + '_BC.bc2', 'w') as f: 
            f.write(h)
        with open(fname + '_BC.bc2', 'a') as f:
            xyz = np.copy(allxyz)
            xyz[:,2] = 0
            np.savetxt(f, xyz, fmt='%.6e')
            f.write(':EndNodes\n:BeginElements ' + str(nele) + '\n')
            np.savetxt(f, d, fmt='%d')
            f.write(':EndElements\n:BeginTable ' + str(2*(nrow+ncol-2)) + ' 15\n')
            with open(fname + '_BC.cli', 'r') as g:
                lines = g.read()
            f.write(lines[:-1])
            f.write(':EndTable\n\n')
        print(' [done]')


@jit(nopython=True)
def compute_mesh(nrow, ncol, nele):
    """
    Build the finite element mesh (compuational).

    Numba nopyton mode is on.
    
    Parameters
    ----------
    nrow : int
        number of rows of beck_bed
    ncol : int
        number of columns of beck_bed
    nele : int
        number of triangular elements in mesh
    
    Returns
    -------
    tri_index : ndarray
        index information of finite element triangles

    """
    tri_index = np.zeros((nele, 3))
    for i in range(nrow-1):
        for j in range(NUM):
            if j == 0:
                tri_index[i*4*NUM+j*4, 0] = (i+1)+(2*j+1)*nrow
                tri_index[i*4*NUM+j*4, 1] = (i+1)
                tri_index[i*4*NUM+j*4, 2] = (i+2)

                tri_index[i*4*NUM+j*4+1, 0] = (i+1)+(2*j+1)*nrow
                tri_index[i*4*NUM+j*4+1, 1] = (i+2)
                tri_index[i*4*NUM+j*4+1, 2] = (i+2)+(2*j+1)*nrow
            else:
                tri_index[i*4*NUM+j*4, 0] = (i+1)+(2*j+1)*nrow
                tri_index[i*4*NUM+j*4, 1] = (i+1)+(2*j-1)*nrow
                tri_index[i*4*NUM+j*4, 2] = (i+2)+(2*j-1)*nrow

                tri_index[i*4*NUM+j*4+1, 0] = (i+1)+(2*j+1)*nrow
                tri_index[i*4*NUM+j*4+1, 1] = (i+2)+(2*j-1)*nrow
                tri_index[i*4*NUM+j*4+1, 2] = (i+2)+(2*j+1)*nrow
            
            tri_index[i*4*NUM+j*4+2, 0] = (i+1)+2*j*nrow
            tri_index[i*4*NUM+j*4+2, 1] = (i+1)+2*(j+1)*nrow
            tri_index[i*4*NUM+j*4+2, 2] = (i+2)+2*(j+1)*nrow

            tri_index[i*4*NUM+j*4+3, 0] = (i+1)+2*j*nrow
            tri_index[i*4*NUM+j*4+3, 1] = (i+2)+2*(j+1)*nrow
            tri_index[i*4*NUM+j*4+3, 2] = (i+2)+2*j*nrow
    return tri_index


def migration(s, x, y, cur_flt, cur_lag, theta, t):
    """
    Channel migration
    
    Parameters
    ----------
    s : ndarray 
        streamwise distance
    x : ndarray
        x coordinate of river centerline
    y : ndarray
        y coordinate of river centerline
    cur_flt : ndarray
        filtered curvature
    cur_lag : ndarray
        filtered curvature with considering phase lag
    theta : ndarray
        angular ampitude
    t : integer
        current # of time step
    
    Returns
    -------
    s : ndarray 
        streamwise distance after migration
    x : ndarray
        x coordinate of river centerline after migration
    y : ndarray
        y coordinate of river centerline after migration

    """
    if MIGRATION:
        if np.mod(t, LPRINT) == 0 and np.mod(t, GPRINT) != 0:
            print('\n[Time Step '+str(t)+'/'+str(TSTEPS)
                  +']\n+> Running channel migration...', end='')
        elif np.mod(t, LPRINT) == 0 and np.mod(t, GPRINT) == 0:
            print('\n[Time Step with graphic printout '+str(t)+'/'+str(TSTEPS)
                  +']\n+> Running channel migration...', end='')

        beta = WIDTH/2/DEPTH
        A = 3.8*(1+beta/6.96*np.exp(-6.96/beta))

        s = s/WIDTH
        x = x/WIDTH
        y = y/WIDTH
        cur_flt = cur_flt*WIDTH
        cur_lag = cur_lag*WIDTH

        chi = (np.sqrt((x[0]-x[-1])**2 + (y[0]-y[-1])**2)/s[-1])**(1/3)
        a1 = UB0*(2*np.random.random(1)[0] - 1) + chi*C0
        a2 = 2*CF0*beta*chi
        a3 = -chi
        a4 = CF0*beta*(chi**5*FR0**2 + (A+1)*chi**2 + 5*np.sqrt(CF0)*(A + chi**2*FR0**2))
        ub = a1*np.exp(-a2*s) + a3*cur_flt + a4*cur_lag 

        x += E0*ub*DT*np.sin(theta)
        y -= E0*ub*DT*np.cos(theta)

        x = x*WIDTH
        y = y*WIDTH
        
        for j in range(1, x.size):
            s[j] = np.sqrt((x[j]-x[j-1])**2 + (y[j]-y[j-1])**2) + s[j-1]
        if np.mod(t, LPRINT) == 0:
            print(' [done]')
    return s, x, y


@jit(nopython=True)
def cutoff(s, x, y):
    """
    Find neck cutoff. If found, remake centerline.

    Criterion: 
    distance between 2 nodes > WIDTH 
    and 
    index difference between the 2 nodes > 4*NUM

    Parameters
    ----------
    s : ndarray 
        streamwise distance
    x : ndarray
        x coordinate of river centerline
    y : ndarray
        y coordinate of river centerline
    
    Returns
    -------
    s : ndarray 
        streamwise distance after cutoff (if found cutoff)
    x : ndarray
        x coordinate of river centerline after cutoff (if found cutoff)
    y : ndarray
        y coordinate of river centerline after cutoff (if found cutoff)
    oxbowx : ndarray
        x coordinate of oxbow lake after cutoff (if found cutoff)
    oxbowy : ndarray
        y coordinate of oxbow lake after cutoff (if found cutoff)
    found_cutoff : boolean

    """
    oxbowx, oxbowy = np.zeros(0), np.zeros(0)
    found_cutoff = False
    if MIGRATION:
        for i in range(1, s.size): 
            for j in range(1, s.size):
                if j-i > 4*NUM and np.sqrt((x[i]-x[j])**2+(y[i]-y[j])**2) < WIDTH:
                    oxbowx, oxbowy = np.copy(x[i+1:j]), np.copy(y[i+1:j])
                    x = np.concatenate((x[:i+1], x[j:]), axis=0)
                    y = np.concatenate((y[:i+1], y[j:]), axis=0)
                    found_cutoff = True
                    s = np.zeros(x.size)
                    for j in range(1, x.size):
                        s[j] = s[j-1] + np.sqrt((x[j]-x[j-1])**2 + (y[j]-y[j-1])**2)
                    return s, x, y, oxbowx, oxbowy, found_cutoff
    return s, x, y, oxbowx, oxbowy, found_cutoff

    
def make_gif():
    """
    Make channel migration movie in gif format.
    Inspired from:
    https://gist.github.com/engineersportal/246155e07556c7ae2dce6a37593b1f31

    """
    if MIGRATION:
        import imageio
        for n, JPG_DIR in enumerate(JPG_DIRS):
            images, image_file_names = [], []
            for file_name in os.listdir(JPG_DIR):
                if file_name.endswith('.jpg'):
                    image_file_names.append(file_name)       
            sorted_files = sorted(image_file_names, key=lambda y: int(y.split('_')[1]))
            for i in range(len(sorted_files)):       
                file_path = os.path.join(JPG_DIR, sorted_files[i])
                images.append(imageio.imread(file_path))
            imageio.mimsave(FNAME.rsplit('.', 1)[0] + '_migration' + str(n) + '.gif', images, 'GIF', loop=1, fps=FPS)


def job_done():
    """
    Clean cache if there is any. Print job done.
    
    Parameters
    ----------
    None
    
    Returns
    -------
    None

    """
    try:
        shutil.rmtree('__pycache__')
    except OSError:
        pass
    print('+> My job is done\n')
    input('Press <Enter> to quit\n')


def main():
    """
    Execute the workflow of pyRiverBed.
    
    Parameters
    ----------
    None
    
    Returns
    -------
    None

    """
    print_banner()
    params = read_steering()
    s, x, y, cur, theta = build_kinoshita()
    s, x, y, cur, theta = read_centerline(s, x, y, cur, theta)
    s, x, y, cur, theta = extend_centerline(s, x, y, cur, theta)
    for t in range(TSTEPS+1):
        cur, theta = tan2curv(s, x, y)
        cur_ori = np.copy(cur)
        cur = filter_curvature(cur, t)
        cur_flt = np.copy(cur)
        cur = lag(s, cur, t)
        cur_lag = np.copy(cur)
        beck_bed = build_beck(cur, s, t)
        allxyz = offset_all(x, y, beck_bed, t)
        if t == 0:
            write_xyz_file(allxyz)
            write_mesh_file(allxyz, beck_bed)
            oxbowxList, oxbowyList = [], []
            centerlinexList, centerlineyList = [], []
        if np.mod(t, GPRINT) == 0:
            centerlinexList.append(x)
            centerlineyList.append(y)
            mf.make_figure(x, y, allxyz, cur_ori, cur_flt, cur_lag, s, beck_bed,
                           params, t, oxbowxList, oxbowyList, centerlinexList, centerlineyList)
        if t == TSTEPS:
            break
        s, x, y = migration(s, x, y, cur_flt, cur_lag, theta, t)
        s, x, y, oxbowx, oxbowy, found_cutoff = cutoff(s, x, y)
        s, x, y = smooth_centerline(x, y)
        s, x, y, cur, theta = resample_centerline(s, x, y)
        if found_cutoff:
            oxbowxList.append(oxbowx)
            oxbowyList.append(oxbowy)
    make_gif()
    job_done()


if __name__ == '__main__':
    main()
    
