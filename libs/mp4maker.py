import argparse
import xarray as xr
import matplotlib.pyplot as plt
from pathlib import Path
from xarray import Dataset
from matplotlib.animation import FuncAnimation, FFMpegWriter
from utils import forceNamingConvention

def make_animation(ds: Dataset, variable: str, cmap: str='viridis') -> FuncAnimation:

    x=ds.x
    y=ds.y
    time=ds.time

    fig = plt.figure(figsize = (8, 8), dpi=300)
    ax = plt.axes()
    ax.set_aspect('equal')
    

    def update(frame):

        timeslice = ds.sel(time=frame)
        data = timeslice[variable]

        # plot data
        ax.clear()
        plot = ax.pcolormesh(x, y, data, shading='auto', cmap=cmap)
        ax.set_title(f'year = {frame}')
        ax.set_axis_off()

        return plot

    animation = FuncAnimation(
        fig, 
        update, 
        frames = time.values, 
        blit=False
        )
    
    return animation

def main(args) -> None:

    netcdf = Path(args.netcdf)
    variable = args.variable
    outfile = Path(args.outfile)
    cmap = args.cmap if args.cmap else 'Blues'

    if outfile.is_file() and not args.overwrite:
        print(f'{outfile.name} already exists')
        return

    ds = xr.open_dataset(netcdf)
    ds = forceNamingConvention(ds)
    animation = make_animation(ds, variable, cmap=cmap)
    writervideo = FFMpegWriter(fps=30, bitrate=5000)
    animation.save(outfile, writer=writervideo, dpi=300)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Process inputs and options"
        )
    
    # add arguments
    parser.add_argument("netcdf", type=str, help="netcdf to turn into an mp4") 
    parser.add_argument("variable", type=str, help="variable to extract from netcdf")
    parser.add_argument("outfile", type=str, help="save path for output mp4")

    # add optional arguments
    parser.add_argument("--cmap", type=int, help="colormap for pcolormesh")
    parser.add_argument("--overwrite", action="store_true", 
                        help="Will overwrite outfile if it already exists")

    args = parser.parse_args()
    main(args)