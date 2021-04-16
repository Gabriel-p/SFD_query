
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

# Name for the region, and (ra, de) center
name, cx, cy = "NGC2516", 119.5166667, -60.7533332

# Length of the box
length = 2.
# Interpolation step
step = 0.01

# Load fits files
dmap = sfdmap.SFDMap('sfddata-master')

# RA coords
x = np.arange(cx - length, cx + length, step)
# DE coords
y = np.arange(cy - length, cy + length, step)

# Query dust maps
e_bv = []
for xi in x:
    e_bv.append(dmap.ebv(xi, y))
e_bv = np.array(e_bv).T

# idx_max = np.unravel_index(e_bv.argmax(), e_bv.shape)
# print(idx_max, e_bv[idx_max[0]][idx_max[1]])

gs_unit = 5
gs_x, gs_y = 2, 1
fig = plt.figure(figsize=(gs_unit * gs_x, gs_unit * gs_y))
gs = gridspec.GridSpec(gs_y, gs_x)

ax = plt.subplot(gs[0:1, 0:1])
ax.minorticks_on()
xmin, xmax = -length, length
ymin, ymax = -length, length
# xmin, xmax = np.amin(x), np.amax(x)
# ymin, ymax = np.amin(y), np.amax(y)
im = plt.imshow(e_bv, origin='lower',  # cmap="RdYlBu_r",
                extent=(xmin, xmax, ymin, ymax))

plt.scatter(0., 0., marker='x', s=60, c='r')
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
out_fig = name + "_SFD.png"
plt.savefig(out_fig, dpi=300, bbox_inches='tight')
