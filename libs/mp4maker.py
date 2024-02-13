import os
import re
import sys
import pandas as pd
import xarray as xr
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.axes_grid1 import Divider, Size

filename = sys.argv[1]

# load ensemble csv
csvpath = '/Users/jonniebarnsley/Code/phd/Postprocessing/Pliocene_control/AIS_PPE_control_ensemble_SLC.csv'
data = pd.read_csv(csvpath, index_col=0, dtype=float)
data = data[data.index < 10_000]
data = data[data.index > 0]

# get slc timeseries
run_num = re.search('run(\d+)_control_thickness_0lev.nc', filename).group(1)
run = data[run_num]
slc = run.values

# check if already processed
mp4path = '/Users/jonniebarnsley/data/phd/Control/mp4s_timeseries/run{}_control_thickness.mp4'.format(run_num)
if os.path.exists(mp4path):
    sys.exit()

# load netcdf
dir = '/Users/jonniebarnsley/data/phd/Control/thickness'
filepath = os.path.join(dir, filename)
ds = xr.open_dataset(filepath)

x=ds.x
y=ds.y
t=ds.t

# plot options
font = {'weight' : 'normal',
        'size'   : 8}
mpl.rc('font', **font)

fig = plt.figure(figsize=(8, 5))

# spatial axis
h = [Size.Fixed(0.0), Size.Fixed(4.5)]
v = [Size.Fixed(0.0), Size.Fixed(4.5)]
divider = Divider(fig, (0, 0, 1, 1), h, v, aspect=False)
ax1 = fig.add_axes(divider.get_position(),
                  axes_locator=divider.new_locator(nx=1, ny=1))

# timeseries axis
h = [Size.Fixed(5), Size.Fixed(2.5)]
v = [Size.Fixed(1.25), Size.Fixed(2.)]
divider = Divider(fig, (0, 0, 1, 1), h, v, aspect=False)
ax2 = fig.add_axes(divider.get_position(),
                  axes_locator=divider.new_locator(nx=1, ny=1))

ax1.set_aspect('equal')

ax2.set_xlim([0, 10000])
ax2.spines['bottom'].set_position('zero')
ax2.spines['left'].set_position('zero')
ax2.spines['top'].set_color('none')
ax2.spines['right'].set_color('none')
ax2.set_xlabel('time ($years$)')
ax2.set_ylabel('sea level contribution ($m$)')

# plot timeseries and initial marker location
ax2.plot(t, slc)
point, = ax2.plot(30, 0, 'o', c='r')

def update(frame):

    '''
    Updates the figure at each frame with the appropriate thickness data
    and point location.
    '''
    subsetted_ds = ds.sel(t=frame)
    thk = subsetted_ds.thickness

    # thickness data
    ax1.clear()
    plot = ax1.pcolormesh(x, y, thk, shading='auto', cmap='Blues')
    ax1.set_title(f'year = {frame}')
    ax1.set_axis_off()

    # timeseries
    sl = run.loc[run.index == frame].values
    point.set_data([frame], [sl])

    return plot, point

ani = animation.FuncAnimation(
    fig, 
    update, 
    frames = np.arange(30, 10020, 30), 
    blit=True
    )

writervideo = animation.FFMpegWriter(fps=30)
ani.save(f'/Users/jonniebarnsley/data/phd/Control/mp4s_timeseries/run{run_num}_thickness.mp4', writer=writervideo)
plt.close(fig)