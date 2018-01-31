# -*- coding: utf-8 -*-
"""
Last updated on Jan 31 2018

@author: Mengyao XUE
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.path import Path
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import Longitude
from astropy.table import Table

def get_gal(c):
    gl_all = Longitude(c.galactic.l)
    gl_all.wrap_angle = 180 * u.deg

    gl_start=np.array(gl_all.rad[np.argmax(gl_all)+1:])
    gl_end=np.array(gl_all.rad[:np.argmax(gl_all)])
    gb_start=np.array(c.galactic.b.rad[np.argmax(gl_all)+1:])
    gb_end=np.array(c.galactic.b.rad[:np.argmax(gl_all)])
    gl=np.concatenate((gl_start,gl_end))*(-1)
    gb=np.concatenate((gb_start,gb_end))

    gl0=np.full(90,gl[0])
    gl1=np.full(90,gl[len(gl)-1])
    gb0=np.linspace(-1.57,gb[0],90)
    gb1=np.linspace(gb[len(gb)-1],-1.57,90)
    for i in range(len(gl)-1):
        if (gl[i]*gl[i+1]<0 and gl[i]>3):
            print gl[i],gl[i+1]
            gl2=np.full(90,gl[i])
            gb2=np.linspace(gb[i],1.57,90)
            gl3=np.full(90,gl[i+1])
            gb3=np.linspace(1.57,gb[i+1],90)
            gl=np.concatenate((gl0,gl[:i],gl2,gl3,gl[i+1:],gl1))
            gb=np.concatenate((gb0,gb[:i],gb2,gb3,gb[i+1:],gb1))
            break

    print i
    return gl,gb


def get_path(glc,gbc):

    verts=zip(glc,gbc)
    codes = np.ones(len(glc), int) * Path.LINETO
    codes[0] = Path.MOVETO
    codes[len(glc)-1] = Path.CLOSEPOLY
    sPath=Path(verts, codes)
    return sPath


if __name__ == "__main__":
    data_psr=open('psrcatlog1.56.csv').readlines()
    ra_S=[]
    dec_S=[]
    ra_N=[]
    dec_N=[]
    for line_psr in data_psr:
        if line_psr.startswith ('#'):
            continue
        words=line_psr.split(',')
        if float(words[3]) < 30:
            ra_S.append(float(words[2]))
            dec_S.append(float(words[3]))

        else:
            ra_N.append(float(words[2]))
            dec_N.append(float(words[3]))

    data_MWA=open('MWA_PSR_20170665.csv').readlines()
    ra_MWA=[]
    dec_MWA=[]

    for line_MWA in data_MWA:
        if line_MWA.startswith ('#'):
            continue
        words_MWA=line_MWA.split(',')
        ra_MWA.append(float(words_MWA[2]))
        dec_MWA.append(float(words_MWA[3]))


    # To plot the celestial equator in galactic coordinates
    alpha = np.linspace(-180.,180.,3600.)
    delta0 = np.zeros(len(alpha))
    cline_0 = SkyCoord(ra=alpha*u.degree, dec=(delta0-10)*u.degree)
    gl_0, gb_0 = get_gal(cline_0)
    cline_1 = SkyCoord(ra=alpha*u.degree, dec=(delta0+30)*u.degree)
    gl_1, gb_1 = get_gal(cline_1)
    sPath_1 = get_path(gl_1,gb_1)
    cline_2 = SkyCoord(ra=alpha*u.degree, dec=(delta0-50)*u.degree)
    gl_2, gb_2 = get_gal(cline_2)
    sPath_2 = get_path(gl_2,gb_2)

    cpsr_S = SkyCoord(ra=ra_S*u.deg, dec=dec_S*u.deg)
    gl_S = Longitude(cpsr_S.galactic.l)
    gl_S.wrap_angle = 180 * u.deg
    cpsr_N = SkyCoord(ra=ra_N*u.deg, dec=dec_N*u.deg)
    gl_N = Longitude(cpsr_N.galactic.l)
    gl_N.wrap_angle = 180 * u.deg
    cpsr_MWA = SkyCoord(ra=ra_MWA*u.deg, dec=dec_MWA*u.deg)
    gl_MWA = Longitude(cpsr_MWA.galactic.l)
    gl_MWA.wrap_angle = 180 * u.deg



    plt.figure(figsize=(10,6),dpi=300)
    bx = plt.subplot(111, projection = 'mollweide')
    #bx = plt.subplot(111, projection = 'aitoff')

    p_S=bx.scatter(-gl_S.rad, cpsr_S.galactic.b.rad, 1.5, lw=0, marker='o', color ='mediumblue', zorder=2, label='All cataloged pulsars (PSRCAT v1.56)')
    p_N=bx.scatter(-gl_N.rad, cpsr_N.galactic.b.rad, 1.5, lw=0, marker='o', color ='gray', zorder=2)
    p_MWA=bx.scatter(-gl_MWA.rad, cpsr_MWA.galactic.b.rad, 2.1, marker='o', color ='red', zorder=3, label='MWA-VCS incoherent-sum detected pulsars')

    p0 = bx.plot(gl_0, gb_0, linewidth=0.8, dashes=[5, 5], color='dimgray', zorder=4, label='Declination limit of the observable sky from LOFAR')
    p1 = bx.plot(gl_1, gb_1, linewidth=0.8,color='gray', zorder=6)
    p2 = bx.plot(gl_2, gb_2, linewidth=0.8,color='gray', zorder=6)

    spch_1 = patches.PathPatch(sPath_1, facecolor='k', edgecolor='none', alpha=0.1, zorder=1, label='Observable sky from the MWA')
    spch_2 = patches.PathPatch(sPath_2, facecolor='gray', edgecolor='none', alpha=0.5, zorder=0, label='Exclusive sky of the MWA at 80-300 MHz')
    bx.add_patch(spch_1)
    bx.add_patch(spch_2)

    handles, labels = bx.get_legend_handles_labels()
    handles = handles[::-1]
    labels = labels[::-1]
    bx.legend(handles, labels, bbox_to_anchor=(0.65, 1.02,0.34,0.2), loc=3, scatterpoints=1 #, numpoints=1
             ,ncol=1, mode="expand", borderaxespad=0., fontsize=6,handlelength=3)

    #xtick_labels = ['10h', '8h', '6h', '4h', '2h', '0h', '22h', '20h', '18h', '16h', '14h']
    d = u'\N{DEGREE SIGN}'
    xtick_labels = ['150'+d,'120'+d,'90'+d,'60'+d,'30'+d,'0'+d,'330'+d,'300'+d,'270'+d,'240'+d,'210'+d]
    #bx.set_xticklabels(xtick_labels)
    bx.set_xticklabels([],fontsize='xx-small')
    bx.annotate('0'+d, xy=(0.005, 1.05), xycoords='data',fontsize=8)
    bx.annotate('300'+d, xy=(0.93, 1.05), xycoords='data',fontsize=8)
    bx.annotate('60'+d, xy=(-1.051, 1.05), xycoords='data',fontsize=8)
    bx.annotate('240'+d, xy=(2.0, 1.05), xycoords='data',fontsize=8)
    bx.annotate('120'+d, xy=(-2.095, 1.05), xycoords='data',fontsize=8)
    #bx.set_yticklabels(bx.get_yticklabels(),fontsize='xx-small')
    bx.set_xlabel('l',fontstyle='italic')
    bx.set_ylabel('b',fontstyle='italic')
    plt.grid(True, color='gray', lw=0.5)
    plt.rcParams['font.size']= 8
    plt.rcParams['legend.fontsize']= 4
    print plt.rcParams.keys()
    print plt.rcParams['font.size']

    plt.savefig('gal_skycover_201801.png',figsize=(20,12),dpi=600)
    #plt.savefig('gal_skycover.ps')
    #plt.show()
