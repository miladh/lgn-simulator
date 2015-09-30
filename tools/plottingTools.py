import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import axes3d
import time
import colormaps as cmaps
import numpy as np

def plotResponse(data, figsize = (15,8)):
    fig = plt.figure(figsize=figsize)
    plt.plot(data)
    plt.show()


def animate3dPlots(data, figsize = (8,6), cmap = cmaps.viridis,
                        save_animation = False, animation_name = "unnamed" ):

    num_subplots = len(data)
    imshowPlots = []
    nStates = data[0].shape[0]
    nSpatialPoints = data[0].shape[1]
    resolution = 1

    num_cols = 3 if num_subplots >= 3 else num_subplots
    num_rows = int(np.ceil(num_subplots/3.))

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



def animateImshowPlots(data, figsize = (15,8), cmap = cmaps.viridis,
                        save_animation = False, animation_name = "unnamed" ):
    num_subplots = len(data)
    imshowPlots = []
    nStates = data[0].shape[0]

    num_cols = 3 if num_subplots >= 3 else num_subplots
    num_rows = int(np.ceil(num_subplots/3.))

    fig = plt.figure(figsize=figsize)


    def init():
        for i in range(num_subplots):
            imshowPlots[i].set_data(data[i][0,:,:])
        ttl.set_text("")
        return imshowPlots, ttl

    def animate(j):
        for i in range(num_subplots):
            imshowPlots[i].set_data(data[i][j,:,:])
        ttl.set_text("timestep = " + str(j) + " s")
        return imshowPlots, ttl


    ttl = plt.suptitle("",fontsize = 16)
    iplot = [num_rows, num_cols, 0]
    for i in range(num_subplots):
        iplot[2] += 1
        ax = fig.add_subplot(iplot[0], iplot[1], iplot[2])
        imshowPlots.append(ax.imshow(data[i][0,:,:], origin='lower', cmap=cmap,
        interpolation="None"))
        plt.colorbar(imshowPlots[-1], ax=ax, orientation='vertical')
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


    outputFilePath = "/home/milad/Dropbox/projects/edog/extendedDOG/eDOG/DATA/*.h5"
    outputFile = glob(outputFilePath)[0]
    f = h5py.File(outputFile, "r")
    exp = sim.Simulation(f)

    data = [
    exp.stimuli["spatioTemporal"], exp.stimuli["spatioTemporal"]
    ]
    animateImshowPlots(data)
    # animate3dPlots(data)
