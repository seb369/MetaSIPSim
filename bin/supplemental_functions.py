# supplemental_functions.py
## Writen by Samuel Barnett with most code from or adapted from SIPSim created by Nicholas Youngblut
### Youngblut ND, Barnett SE, Buckley DH. SIPSim: A modeling toolkit to predict accuracy and aid design of DNA-SIP experiments. ### Frontiers in Microbiology. 2018;9. doi:10.3389/fmicb.2018.00570.
## This file contains functions used within SIPSim_metagenome. Most relating to calculating BD of DNA fragments.

import numpy as np
import pandas as pd
import scipy
import math
from scipy.integrate import quad
from scipy.stats import norm as normdist

###############################################################################################################################
################################################ Functions for main function ##################################################
###############################################################################################################################


# functions for trig with angle in degrees
cos_d = lambda d : np.cos(np.deg2rad(d))
sin_d = lambda d : np.sin(np.deg2rad(d))
asin_d = lambda x : np.arcsin(x) * 180/np.pi  
acos_d = lambda x : np.arccos(x) * 180/np.pi

def divide(x, divide=1):
    return x / divide

def BD2distFromAxis(BD, D, B, w, I):
    # Convert buoyant density (BD) to distance from axis of rotation (units = cm). If distance > r_max, nan is returned.
    x = np.sqrt(((BD - D) * 2.0*B/w) + I**2)
    return x

def _SphVol(t, r, p2, R12):
    # helper function for axisDist2angledTubeVol
    v1 = t*((2*r)-t)/2
    v2 = 2*np.pi*((p2-t)/R12)
    v3 = np.sin(2*np.pi*((p2-t)/R12))
    return v1 * (v2 - v3)

def _CylWedVol(t, r, b, h):
    # helper function for axisDist2angledTubeVol
    return 2*(h*(t-r+b)/ b) * np.sqrt(r**2-t**2)

def axisDist2angledTubeVol(x, A, r, r_max):
    # Convert distance from axis of rotation to volume of gradient where the BD is >= to the provided BD.
    # return nan if nan provided
    
    a = np.deg2rad(A)
    p1 = r-(r*np.cos(a))
    p2 = r+(r*np.cos(a))
    R12 = p2-p1
    d = r_max-x
    D1 = r_max-p1
    D2 = r_max-p2 
    if x < D2:
        if a == 0:
            z = 1
        else:
            z = np.sin(a)
        h1 = (D2-x)/z
        h2 = (D1-x)/z
        volume1 = (2/3.0)*np.pi*r**3
        volume2 = (0.5)*np.pi*r**2*(h1+h2)
        volume = volume1 + volume2
    elif D1 >= x >= D2:
        volume1 = (1/3.0)*np.pi*p1**2*(3*r-p1)
        volume2 = scipy.integrate.quad(_SphVol, p1, d, args=(r, p2, R12))
        b = (d-p1)/np.cos(a)
        if a == 0:
            h = b
        else:
            h = b/np.tan(a)
        volume3 = scipy.integrate.quad(_CylWedVol, r-b, r, args=(r, b, h))
        volume = volume1+volume2[0]+volume3[0]
    elif r_max >= x > D1:
        volume = (1/3.0)*np.pi*d**2*(3*r-d)
    elif x > r_max:  ####### FIX THIS!!!!!!!!!!!!!!!!!
        print(' '.join(['Error: estimated distance of BD from axis', str(x), 
                        'is greater than maximum distance', str(r_max), 
                        '\nIncrease centrifugation speed or average gradient density']))
        volume = 0
    return volume

def _cylVol2height(v, r): 
    # v = volume (ml)
    # r = tube radius (cm)
    h1 = v / (np.pi * r**2)
    return h1


def _sphereCapVol2height(v, r):
    # v = volume (ml)
    # r = tube radius (cm)
    # height = h**3 - 3*r*h**2 + (3v / pi) = 0
    f = lambda x : x**3 - 3*r*x**2 + 3*v/np.pi
    h2 = scipy.optimize.brentq(f, 0, r*2, maxiter=1000)
    return h2

def tubeVol2vertTubeHeight(v, r):
    # Convert angled tube volume (see axisDist2angledTubeVol) to height in the vertical tube.
    sphere_half_vol = (4.0/3.0 * np.pi * r**3)/2.0

    if v <= sphere_half_vol:
        # height does not extend to cylinder
        h = _sphereCapVol2height(v, r)
    else:
        # height = sphere_cap_height (r) + cylinder_height
        #sphere_cap_height = r #_sphereCapVol2height(sphere_half_vol, r)
        h =  r + _cylVol2height(v - sphere_half_vol, r)
    return h

def vertTubePos_BD_fit(BDs, vert_tube_pos):
    # Making a continuous function: BD ~ vert_tube_pos. This function will be used to convert angle_tube_pos to the BD in the vertical tube gradient. Using polynomial curve fit to define function.
    deg = 3
    
    # trimming NA values
    vert_tube_pos = vert_tube_pos[~np.isnan(vert_tube_pos)]
    BDs = BDs[~np.isnan(vert_tube_pos)]

    # fitting variables (polynomial fitting)
    fit = np.polyfit(vert_tube_pos, BDs, deg=deg)
    return  np.poly1d(fit)

# Calculate diffusion standard deviation

def diffusion_SD(BD, length, diffVar):
    sig = math.sqrt((BD/length)*diffVar)
    return sig

# The following functions allow for calculations of the diffusive boundery layer for a given fragment.

def axisDist2angledTubePos(x, r_max, A, r):
    # Convert distance from axis of rotation to angled tube position
    if np.isnan(x):
        return (x, x)

    # low position
    if(x >= r_max - (r * cos_d(A)) - r):
        # band in rounded bottom of cfg tube
        d = x - (r_max - r)
        a = A - asin_d(d / r)
        LowH = r - r * cos_d(a)
    else:
        # band in cylinder of cfg tube
        d = r_max - (r * cos_d(A)) - r - x
        h_c = d/sin_d(A)
        LowH = r + h_c
    # high position
    if(x > r_max - (r - r * cos_d(A))):
        # Equation for finding the upper band
        d = x - (r_max - r)
        a = A - (180 - asin_d(d/r))
        HighH = r - r * cos_d(a)
    else:
        # This band will be in the cylinder part        
        d = r_max - (r - r * cos_d(A)) - x
        if A == 0:
            h_c = d
        else:
            h_c = d/sin_d(A)
        HighH = r + h_c
        
    return(LowH, HighH)

def angle_tube_pos_to_BD(angle_tube_pos, VTP2BD_func):
    # Converting angle_tube_pos (min/max) to vertical tube BD min/max (span of DBL), converting span to uniform distribution function (min/max of DBL), then making index: (BD : uniform_dist_func)
    low_pos_BD = VTP2BD_func(angle_tube_pos[0])
    high_pos_BD = VTP2BD_func(angle_tube_pos[1])

    return (low_pos_BD, high_pos_BD)

def DBL_prop(DBL_range,  maxBd,  minBd):
    # returns the proportion of the DBL abundance of the fragment that can be found in the given window.
    DBL_sampled = (maxBd-minBd)/(DBL_range[0] - DBL_range[1]) # proportion of DBL found in this window
    return DBL_sampled



# The following functions assist in the final function
def calc_isoconc_point(r_min, r_max):
    # Convert min/max distance from axis of rotation to the isoconcentration point of the gradient (units = cm).
    I = np.sqrt((r_min**2.0 + r_min * r_max + r_max**2.0)/3.0)
    return I

def absolute_to_relative_abundance(row, total_abundance):
    # assigns relative abundance to each read in the library
    relative_abundance = row['abundance']/total_abundance
    return relative_abundance

###############################################################################################################################
############################################## Functions for calculating mamory ###############################################
###############################################################################################################################

def mem_usage(pandas_obj):
    if isinstance(pandas_obj,pd.DataFrame):
        usage_b = pandas_obj.memory_usage(deep=True).sum()
    else: # we assume if not a df it's a series
        usage_b = pandas_obj.memory_usage(deep=True)
    usage_mb = usage_b / 1024 ** 2 # convert bytes to megabytes
    return "{:03.2f} MB".format(usage_mb)
