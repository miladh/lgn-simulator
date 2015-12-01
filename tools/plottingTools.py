import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import axes3d
from mpl_toolkits.axes_grid1 import make_axes_locatable
import time
import colormaps as cmaps
import numpy as np

def simpleaxis(ax):
    """
    Removes axis lines
    """
    # plt.axis('off')
    ax.set_axis_off()
    # ax.spines['top'].set_visible(False)
    # ax.get_xaxis().tick_bottom()


def raster(spike_times, figsize = (12,8), ax = None, ylabel = "Cell Id"):
    """
    Raster plot
    """

    if not isinstance(spike_times, list):
        spike_times  = [spike_times]

    num_cells = len(spike_times)
    if ax ==None:
        f = plt.figure(figsize=figsize)
        ax = f.add_subplot(111)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    ax.tick_params(axis='x', labelsize=12)

    yticks = []
    yticks.append("")
    for i in range(num_cells):
        for t in spike_times[i][0]:
            ax.vlines(t, i + .9, i + 1.1)
        yticks.append(spike_times[i][1])

        plt.ylim(0.,num_cells+0.5)
        plt.xlabel('t [s]', fontsize = 16)
        plt.ylabel(ylabel, fontsize = 16)

    plt.yticks(range(num_cells+1), yticks, fontsize = 12)
    plt.tight_layout()
    return ax

def animate3dPlots(data, figsize = (8,6), cmap = cmaps.viridis, resolution = 0,
                        save_animation = False, animation_name = "unnamed" ):

    num_subplots = len(data)
    imshowPlots = []
    nStates = data[0].shape[0]
    nSpatialPoints = data[0].shape[1]
    resolution = 2**resolution

    num_cols = 2 if num_subplots >= 2 else num_subplots
    num_rows = int(np.ceil(num_subplots/2.))

    plt.ion()
    fig = plt.figure(figsize=figsize)
    xs = np.linspace(-1, 1, nSpatialPoints/resolution)
    ys = np.linspace(-1, 1, nSpatialPoints/resolution)
    X, Y = np.meshgrid(xs, ys)


    iplot = [num_rows, num_cols, 0]
    axes = []
    frames = [None]* num_subplots
    framesOld = [None]* num_subplots
    for i in range(num_subplots):
        iplot[2] += 1
        axes.append(fig.add_subplot(iplot[0], iplot[1], iplot[2], projection='3d'))


    tstart = time.time()
    for i in range(nStates):
        for j in range(num_subplots):
            data3d = data[j][i,::resolution,::resolution]
            framesOld[j] = frames[j]
            frame = axes[j].plot_surface(X, Y, data3d,
                rstride=2, cstride=2, alpha = 0.9,
                linewidth=0.1, antialiased=True, cmap=cmap)

            if framesOld[j] is not None:
                axes[j].collections.remove(framesOld[j])
        fig.canvas.draw()
        plt.pause(0.0001)

    plt.show()



def animateImshowPlots(data, dt = None, figsize = (8,15), cmap = cmaps.viridis,
                        save_animation = False, colorbar = False, animation_name = "unnamed" ):



    num_subplots = len(data)
    imshowPlots = []
    nStates = np.array(data[0][0]).shape[0]
    Nx = np.array(data[0][0]).shape[1]
    Ny = np.array(data[0][0]).shape[2]
    dt = 1 if dt==None else dt

    num_cols = 2 if num_subplots >= 2 else num_subplots
    num_rows = int(np.ceil((num_subplots-1)/2.))+1

    fig = plt.figure(figsize=figsize)

    # extent = self.screenSizeX * np.arctan(1./self.screenDist) * 180./np.pi


    def init():
        for i in range(num_subplots):
            imshowPlots[i].set_data(data[i][0][1,:,:])
        ttl.set_text("")
        return imshowPlots, ttl

    def animate(j):
        for i in range(num_subplots):
            imshowPlots[i].set_data(data[i][0][j,:,:])
        t = j*dt
        ttl.set_text("Time = " + str('%.2f' % (t,)) + "s")
        return imshowPlots, ttl



    ttl = plt.suptitle("",fontsize = 16)
    # ttl = plt.suptitle("",fontsize = 16, x = 0.45, horizontalalignment = "left")

    #Plot stimulus
    colspanStim = 2
    ax = plt.subplot2grid((num_rows, num_cols),(0, 0), colspan = colspanStim)
    simpleaxis(ax)
    ax.set_title(data[0][1])
    imshowPlots.append(ax.imshow(data[0][0][0,:,:], cmap="gray", interpolation="None"))
    if(colorbar):
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(imshowPlots[-1], ax=ax, orientation='vertical', cax = cax)

    k = 1
    for i in range(1, num_rows):
        for j in range(num_cols):
            if(k>= num_subplots):
                break
            ax = plt.subplot2grid((num_rows, num_cols),(i, j))
            simpleaxis(ax)
            ax.set_title(data[k][1])
            imshowPlots.append(ax.imshow(data[k][0][0,:,:], cmap=cmap,
            interpolation="None"))
            if(colorbar):
                divider = make_axes_locatable(ax)
                cax = divider.append_axes("right", size="5%", pad=0.05)
                plt.colorbar(imshowPlots[-1], ax=ax, orientation='vertical',
                cax = cax)


            k+=1

    plt.tight_layout()
    plt.subplots_adjust(top=0.95)

    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                   frames=nStates, interval=20, blit=False)

    if(save_animation):
        anim.save(animation_name + ".mp4",fps=30,  writer="avconv", codec="libx264")

    plt.show()




if __name__ == "__main__":
    import h5py
    from glob import glob
    import Simulation as sim


    #outputFilePath = "/home/milad/Dropbox/projects/edog/extendedDOG/eDOG/DATA/*.h5"
    outputFilePath = "/home/milad/kurs/*.h5"
    outputFile = glob(outputFilePath)[0]
    f = h5py.File(outputFile, "r")
    exp = sim.Simulation(f)

    # spikeTrain = exp.spikeTrain("ganglion", 64, 64, num_trails = 1)
    # spikeTrain2 = exp.spikeTrain("relay", 64, 64, num_trails = 1)

    data = [
     [exp.stimulus["spatioTemporal"], "Stimulus"]
    ,[exp.ganglion["response"]["spatioTemporal"], "Ganglion cell response"]
    ,[exp.ganglion["impulseResponse"]["spatioTemporal"], "Ganglion cell impulse response"]
    ,[exp.interneuron["response"]["spatioTemporal"], "Interneuron"]
    ,[exp.interneuron["impulseResponse"]["spatioTemporal"], "Interneuron"]
    ,[exp.relay["response"]["spatioTemporal"], "Relay cell response"]
    ,[exp.relay["impulseResponse"]["spatioTemporal"], "Relay cell impulse response"]
    ,[exp.cortical["response"]["spatioTemporal"], "Cortical cell response"]
    ,[exp.cortical["impulseResponse"]["spatioTemporal"], "Cortical impulse response"]
    ]

    # raster([[spikeTrain, "Ganglion"], [spikeTrain2, "Relay"]] )
    # print exp.ganglion["impulseResponse"]
    print exp.dt*2**8
    animateImshowPlots(data,exp.dt, colorbar = True, save_animation = True, animation_name = "rat1")
    # animate3dPlots(data, resolution = 3)
    plt.show()
