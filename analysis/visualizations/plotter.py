import numpy as np
import os
import time
import glob
import re
from PIL import Image

import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use("Agg")
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

import imageio.v2 as imageio
from moviepy.video.io.ImageSequenceClip import ImageSequenceClip


class Plotter:
    def __init__(self, figsize=(6, 6), xsize=1, ysize=1, fontsize=15, colorbar=False, bgcolor=None, ftcolor=None, right=0.88):
        
        self.fig, self.ax = plt.subplots(xsize, ysize, figsize=figsize, facecolor=bgcolor, edgecolor=ftcolor)
        
        if colorbar:
            self.fig.subplots_adjust(right=right)
        
        if isinstance(self.ax, np.ndarray):
            self.axes = self.ax.ravel()
        else:
            self.axes = [self.ax]
        self.figsize = figsize
        
        reference_area = 6 * 6
        figure_area = figsize[0] * figsize[1]

        scale = np.sqrt(figure_area / reference_area)

        #Prevent font from becoming too small or large
        scale = np.clip(scale, 0.7, 1.6)

        self.fontsize = fontsize * scale

        self.fontsizes = {
            "label": self.fontsize,
            "tick": self.fontsize * 0.80,
            "legend": self.fontsize * 0.70,
            "title": self.fontsize * 1.10,
            "colorbar": self.fontsize,
            "colorbar_tick": self.fontsize * 0.80,
        }

        self.colorbar = colorbar
        self.config = {"bgcolor": bgcolor, "ftcolor": ftcolor}

        if bgcolor and ftcolor:
            self.fig.patch.set_facecolor(bgcolor)
            for ax in self.axes:

                ax.set_facecolor(bgcolor)

                ax.tick_params(colors=ftcolor)

                ax.xaxis.label.set_color(ftcolor)
                ax.yaxis.label.set_color(ftcolor)

                ax.title.set_color(ftcolor)

                for spine in ax.spines.values():
                    spine.set_color(ftcolor)

        self.legends = {}

    def get_ax(self, i):
        return self.axes[i]
    
    def add_details(self, ax=None, xlabel='', ylabel='', xscale='linear', yscale='linear',
                xlim=None, ylim=None, colorbar=False, norm=None, cmap=None, h=None, cbarlabel='',
                minorticks=True, legend=False, legendloc='upper right', grid=False, 
                invert_xaxis=False, invert_yaxis=False):
        
        axes = self.axes if ax is None else [ax]

        for a in axes:
            a.set(xlabel=xlabel, ylabel=ylabel, xscale=xscale, yscale=yscale)
            if xlim: a.set_xlim(xlim)
            if ylim: a.set_ylim(ylim)
            a.tick_params(axis='both', which='both', direction='in', labelsize=self.fontsizes["tick"])
            if minorticks:
                a.minorticks_on()
            else:
                a.minorticks_off()

            a.xaxis.label.set_size(self.fontsizes["label"])
            a.yaxis.label.set_size(self.fontsizes["label"])

            if colorbar:
                divider = make_axes_locatable(a)
                cax = divider.append_axes("right", size="5%", pad=0.05)
                if norm is not None and cmap is not None:
                    sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
                    sm.set_array([])
                    cbar = self.fig.colorbar(sm, cax=cax, label=cbarlabel, fontsize=self.fontsizes["colorbar"])
                    cbar.ax.tick_params(labelsize=self.fontsizes["colorbar_tick"])
                else:
                    cbar = self.fig.colorbar(h, cax=cax, label=cbarlabel, fontsize=self.fontsizes["colorbar"])
                    cbar.ax.tick_params(labelsize=self.fontsizes["colorbar_tick"])
                cbar.ax.set_aspect('auto')

            if legend:
                bgcolor = self.config.get("bgcolor", None)
                leg = a.legend(loc=legendloc, frameon=False, facecolor=bgcolor, fontsize=self.fontsizes["legend"])
                ftcolor = self.config.get("ftcolor", None)

                if ftcolor is not None:
                    plt.setp(leg.get_texts(), color=ftcolor)

            if grid:
                a.grid(True, zorder=0, alpha=0.5, which="both")
            if invert_xaxis:
                a.invert_xaxis()
            if invert_yaxis:
                a.invert_yaxis()
    
    def add_legends(self, ax=None, type="line", color="", label="", linestyle="-", marker="o", markersize=5, lw=3, alpha=1.0, edgecolor=None):
        
        if ax is None:
            ax = self.axes[0]
        legend = None
        if type == "line":
            legend = Line2D([0], [0], color=color, lw=lw, markerfacecolor=color, label=label, linestyle=linestyle, markersize=markersize)
        elif type == "scatter":
            legend = Line2D([0], [0], marker=marker, color=color, label=label, markerfacecolor=color, markersize=5, linestyle="None")
        elif type == "errorbar":
            legend = Line2D([0], [0], marker=marker, color=color, label=label, markerfacecolor=color, markeredgecolor="k", markersize=10, linestyle="None")
        elif type == "patch":
            legend = Patch(facecolor=color, edgecolor=edgecolor, label=label, alpha=alpha, linestyle=linestyle, lw=lw)
        if legend is not None:
            if ax not in self.legends:
                self.legends[ax] = []
            self.legends[ax].append(legend)

    def save(self, filename, pad=0.4, h_pad=None, w_pad=None, close=True, legend_loc="upper left", legend_divide=True, fontsize=12):
        
        if not self.colorbar:
            self.fig.tight_layout(pad=pad, h_pad=h_pad, w_pad=w_pad)

        for ax, handles in self.legends.items():
            if len(handles) > 6 and legend_divide == True:
                print(f"Number of legends {len(handles)} is too large, dividing the legends in multiple columns.")
                leg = ax.legend(handles=handles, frameon=False, loc=legend_loc, fontsize=self.fontsizes["legend"], ncol=2)
                
                ftcolor = self.config.get("ftcolor", None)

                if ftcolor is not None:
                    plt.setp(leg.get_texts(), color=ftcolor)
            else:
                ax.legend(handles=handles, frameon=False, loc=legend_loc, fontsize=self.fontsizes["legend"], ncol=1)

                ftcolor = self.config.get("ftcolor", None)

                if ftcolor is not None:
                    plt.setp(leg.get_texts(), color=ftcolor)
        
        format = os.path.splitext(filename)[1][1:]
        self.fig.savefig(os.path.join("Images/", filename), dpi=300, format=format, bbox_inches="tight")
        if close:
            plt.close(self.fig)
        print(f"Plot Created at loc: Images/{filename}")
      
    def create_colorbar(self, cmap, vmin, vmax, label="", pad=0.1, left=0.91, bottom=0.1, width=0.025, height=0.80):
        norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        
        # [left, bottom, width, height]
        cax = self.fig.add_axes([left, bottom, width, height])
        
        cbar = self.fig.colorbar(sm, cax=cax, pad=pad)
        cbar.set_label(label, fontsize=self.fontsizes["label"])
        cbar.ax.tick_params(labelsize=self.fontsizes["colorbar_tick"])
        return cbar


class GridPlot(Plotter):

    def __init__(self, figsize=(6, 6), nx=1, ny=1, colorbar=False, bgcolor=None, ftcolor=None, fontsize=15):

        self.fig = plt.figure(figsize=figsize)

        self.nx = nx
        self.ny = ny
        self.colorbar = colorbar
        reference_area = 6 * 6
        figure_area = figsize[0] * figsize[1]

        scale = np.sqrt(figure_area / reference_area)

        #Prevent font from becoming too small or large
        scale = np.clip(scale, 0.7, 1.6)

        self.fontsize = fontsize * scale

        self.fontsizes = {
            "label": self.fontsize,
            "tick": self.fontsize * 0.80,
            "legend": self.fontsize * 0.80,
            "title": self.fontsize * 1.10,
            "colorbar": self.fontsize,
            "colorbar_tick": self.fontsize * 0.80,
        }

        self.config = {"bgcolor": bgcolor, "ftcolor": ftcolor}

        # allocate axes array
        self.ax = np.empty((ny, nx), dtype=object)

        x_length = 1.0
        if colorbar:
            x_length = 0.88  # leave space for colorbar

        y_length = 1.0

        for j in range(ny):        # rows (y)
            for i in range(nx):    # columns (x)

                ax = self.fig.add_axes([
                    x_length * i / nx,
                    y_length * (ny - j - 1) / ny,   # <-- flip y (important)
                    x_length / nx,
                    y_length / ny
                ])

                # store axis
                self.ax[j, i] = ax

                # remove inner labels
                if i != 0:
                    ax.set_yticklabels([])
                if j != ny - 1:
                    ax.set_xticklabels([])

        self.legends = {}

class Video:
    def __init__(self, input_dir, output_dir=None):
        """
        Parameters
        ----------
        input_dir : str
            Directory containing the rendered frames.

        output_dir : str, optional
            Directory where the video will be written.
            If None, uses input_dir.
        """
        self.input_dir = input_dir
        self.output_dir = output_dir or input_dir

        os.makedirs(self.output_dir, exist_ok=True)

    ### Frame discovery

    def get_frames(self, prefix="frame", pattern="*.png"):
        """
        Find and sort images frames.
        """

        frames = glob.glob(os.path.join(self.input_dir, f"{prefix}_{pattern}"))

        if not frames:
            raise FileNotFoundError(f"No frames found in {self.input_dir} matching '{prefix}_{pattern}'")

        frames.sort(key=lambda x: int(
            re.search(
            rf"{re.escape(prefix)}_(\d+)", os.path.basename(x)
            ).group(1))
        )

        return frames


    ### Create video

    def make_video(self, prefix="frame", pattern="*.png", output_name="movie.mp4", fps=20, typ="mp4", frame_size=800):

        """
        Assemble rendered frames into a video.

        Parameters
        ----------
        pattern : str
            Filename pattern for input frames.

        output_name : str
            Name of the output video.

        fps : int or float
            Frames per second.

        typ : str
            'mp4' or 'gif'
        """
        frame_paths = self.get_frames(prefix, pattern)

        images = [imageio.imread(path) for path in frame_paths]

        valid_frames = []

        print(f"Found {len(valid_frames)} frames.")

        for i, image in enumerate(images):
            if np.shape(image) == (frame_size, frame_size, 4):
                valid_frames.append(image)

        output_path = os.path.join(self.output_dir, output_name)

        if typ == "gif":
            imageio.mimsave(output_path, valid_frames, duration=1 / fps, palettesize=256, subrectangles=True)

        elif typ == 'mp4':
            clip = ImageSequenceClip(valid_frames, fps=fps)
            clip.write_videofile(output_path, codec='libx264', fps=fps)

        else:
            raise ValueError(
                f"Unknown video type '{typ}'. "
                "Use 'mp4' or 'gif'."
            )

        print(f"Video created at: {output_path}")

        return output_path
'''
For a dark background slide
Background	#07111F
Text	#E8EDF7
Accent	#F6B73C
Blue	#56B4E9
Orange	#E69F00
Green	#009E73
Purple	#CC79A7
Yellow	#F0E442
'''