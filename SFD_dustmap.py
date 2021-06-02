
import sfdmap
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable

"""
Query the Schlegel, Finkbeiner & Davis (1998) dust map FITS files.

By default, a scaling of 0.86 is applied to the map values to reflect
the recalibration by Schlafly & Finkbeiner (2011)

By default the coordinates are assumed to be in degrees in the ICRS
coordinate system (e.g., "J2000")

SFD dust map package: https://github.com/kbarbary/sfdmap
"""

# Name for the region, (ra, de) center, and radius (arcmin)
cl_list = (
    ('Haffner 14',   116.2125,  -28.3667,   5.0),
    ('Ruprecht 41',  118.4500,  -26.9611,   5.0),
    ('Ruprecht 42',  119.4000,  -25.9167,   5.0),
    ('Ruprecht 44',  119.7125,  -28.5833,   5.0),
    ('Ruprecht 152', 118.6167,  -38.2372,   5.0),
    )
#     ('NGC_188', 11.798  , 85.244),
#     ('King_2',  12.741  , 58.188),
#     ('IC_166',  28.094  , 61.857),
#     ('NGC_1193',    46.486  , 44.383),
#     ('Berkeley_12', 71.1    , 42.691),
#     ('Berkeley_14', 74.935  , 43.488),
#     ('NGC_1798',    77.914  , 47.691),
#     ('Berkeley_17', 80.13   , 30.574),
#     ('Berkeley_18', 80.531  , 45.442),
#     ('Berkeley_19', 81.014  , 29.575),
#     ('Berkeley_70', 81.457  , 41.951),
#     ('Berkeley_21', 87.93   , 21.812),
#     ('NGC_2141',    90.734  , 10.451),
#     ('NGC_2158',    91.862  , 24.099),
#     ('NGC_2243',    97.395  , -31.282),
#     ('Berkeley_23', 98.318  , 20.535),
#     ('Trumpler_5',  99.126  , 9.465),
#     ('Collinder_110',   99.677  , 2.069),
#     ('Berkeley_32', 104.53  , 6.433),
#     ('Tombaugh_2',  105.773 , -20.82),
#     ('Melotte_66',  111.573 , -47.685),
#     ('Berkeley_39', 116.702 , -4.665),
#     ('NGC_2506',    120.01  , -10.773),
#     ('NGC_2682',    132.846 , 11.814),
#     ('Collinder_261',   189.519 , -68.377),
#     ('Trumpler_20', 189.882 , -60.637),
#     ('NGC_6253',    254.778 , -52.712),
#     ('NGC_6791',    290.221 , 37.778),
#     ('NGC_6819',    295.327 , 40.19),
#     ('NGC_7142',    326.29  , 65.782),
#     ('King_11', 356.912 , 68.636),
#     ('NGC_7789',    359.334 , 56.726)
# )

# Length of the box (degrees)
length = 1.
# Interpolation step
step = 0.01

# Load fits files
dmap = sfdmap.SFDMap('sfddata-master')


def querySFD(cx, cy, length):
    # RA coords
    x = np.arange(cx - length, cx + length, step)
    # DE coords
    y = np.arange(cy - length, cy + length, step)
    # Query dust maps
    e_bv = []
    for xi in x:
        e_bv.append(dmap.ebv(xi, y))
    e_bv = np.array(e_bv).T

    return e_bv


for cluster in cl_list:

    # Name for the region, and (ra, de) center
    name, cx, cy, rad = cluster
    # name, cx, cy = cluster
    # rad = 5

    rad_arcmin = rad / 60.
    e_bv_rad = querySFD(cx, cy, rad_arcmin)
    e_bv_median, e_BV_min, e_BV_max = np.median(e_bv_rad), e_bv_rad.min(),\
        e_bv_rad.max()
    print("{}: median={:.3f}, min={:.3f}, max={:.3f}".format(
        name, e_bv_median, e_BV_min, e_BV_max))

    # For plotting
    e_bv = querySFD(cx, cy, length)

    gs_unit = 5
    gs_x, gs_y = 2, 1
    fig = plt.figure(figsize=(gs_unit * gs_x, gs_unit * gs_y))
    gs = gridspec.GridSpec(gs_y, gs_x)

    ax = plt.subplot(gs[0:1, 0:1])
    plt.title("median={:.3f}, min={:.3f}, max={:.3f}".format(
        e_bv_median, e_BV_min, e_BV_max))
    ax.minorticks_on()
    xmin, xmax = -length, length
    ymin, ymax = -length, length
    # xmin, xmax = np.amin(x), np.amax(x)
    # ymin, ymax = np.amin(y), np.amax(y)
    im = plt.imshow(e_bv, origin='lower',  # cmap="RdYlBu_r",
                    extent=(xmin, xmax, ymin, ymax))

    plt.scatter(0., 0., marker='x', s=60, c='r')
    circ = plt.Circle((0, 0), rad_arcmin, color='r', fill=False)
    ax.add_patch(circ)
    plt.grid(c='k', ls=":", which='both')

    ax.set_xlabel(r"$\alpha^{*}$")
    ax.set_ylabel(r"$\delta^{*}$")
    ax.set_ylim(ymin, ymax)
    ax.set_xlim(xmax, xmin)
    ax.set_xticks((xmin, 0, xmax))
    ax.set_yticks((ymin, 0, ymax))

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.05)
    cbar = plt.colorbar(im, cax=cax)
    # cbar.set_ticks([
    #     np.min(kde), (np.max(kde) + np.min(kde)) * .5, np.max(kde)])
    cbar.ax.minorticks_off()
    cbar.set_label(r'E$_{(B-V)}$')

    fig.tight_layout()
    out_fig = 'out/' + name + "_SFD.png"
    plt.savefig(out_fig, dpi=300, bbox_inches='tight')
    plt.clf()
    plt.close("all")
