import matplotlib.pyplot as plt
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


def raster(spike_times,
           ax = None,
           figsize = (12,8),
           ylabel = "Cell Id",
           title = None):
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
        for t in spike_times[i][:]:
            ax.vlines(t, i + .9, i + 1.1)

    plt.ylim(0.8, num_cells+0.2)
    plt.xlabel('t [s]', fontsize = 16)
    plt.ylabel(ylabel, fontsize = 16)
    plt.yticks(range(1,num_cells+1), fontsize = 12)

    if not title==None:
        ax.set_title(title)
    plt.tight_layout()
    return ax


def animateImshowPlots(data,
                       dt = None,
                       figsize = (8,15),
                       cmap = cmaps.inferno,
                       colorbar = False,
                       save_animation = False,
                       animation_name = "unnamed" ):

    """
    Animate imshow plots
    """
    import matplotlib.animation as animation
    from mpl_toolkits.mplot3d import axes3d
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    import time

    num_subplots = len(data)
    imshowPlots = []
    nStates = np.array(data[0][0]).shape[0]
    Nx = np.array(data[0][0]).shape[1]
    Ny = np.array(data[0][0]).shape[2]
    dt = 1 if dt==None else dt

    num_cols = 2 if num_subplots >= 2 else num_subplots
    num_rows = int(np.ceil((num_subplots-1)/2.))+1

    fig = plt.figure(figsize=figsize)

    def init():
        for i in range(num_subplots):
            imshowPlots[i].set_data(data[i][0][0,:,:])
        ttl.set_text("")
        return imshowPlots, ttl

    def animate(j):
        for i in range(num_subplots):
            imshowPlots[i].set_data(data[i][0][j,:,:])
            # imshowPlots[i].autoscale()
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
            interpolation="None",)) #vmin=0.0, #vmax=-1

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



def imshowPlotsOfImpulseResponses(data,
                                  x_imshow=True,
                                  y_imshow=True,
                                  idx=0,
                                  idy=0,
                                  figsize=(14,8),
                                  cmap=cmaps.inferno,
                                  colorbar=True,
                                  save_figure=False,
                                  figure_name="unnamed"):

    """
    Imshow plots of impulse response functions
    """
    num_cols = len(data)
    num_rows = int(x_imshow) + int(y_imshow)

    fig, axarr=plt.subplots(num_rows, num_cols, figsize=figsize, sharex=True, sharey=True)

    for j in range(num_cols):
        axarr[0,j].set_title(data[j][1])
        i=0
        if(x_imshow):
            axarr[i,j].set_adjustable('box-forced')
            im = axarr[i,j].imshow(data[j][0][:,idy,:], cmap=cmap, origin="lower")
            axarr[i,0].set_xlabel(r"$x(\theta)$")
            axarr[-1,j].set_ylabel(r"$\tau(ms)$")
            if(colorbar):
                fig.colorbar(im, ax=axarr[i,j],orientation='horizontal')
            i+=1
        if(y_imshow):
            axarr[i,j].set_adjustable('box-forced')
            axarr[i,j].imshow(data[j][0][:,:,idx], cmap=cmap, origin="lower")
            axarr[i,0].set_xlabel(r"$y(\theta)$")
            axarr[-1,j].set_ylabel(r"$\tau(ms)$")
            if(colorbar):
                fig.colorbar(im, ax=axarr[i,j], orientation='horizontal')
            i+=1
    if(save_figure):
        fig.save(animation_name + ".svg")
    plt.show()






def plot3dOfImpulseResponses(data,
                             x_3d=True,
                             y_3d=True,
                             idx=0,
                             idy=0,
                             figsize=(15,10),
                             cmap=cmaps.inferno,
                             colorbar=False,
                             save_figure=False,
                             figure_name="unnamed"):

    """
    3D plots of impulse response functions
    """
    from mpl_toolkits.mplot3d import Axes3D
    num_cols = len(data)
    num_rows = int(x_3d) + int(y_3d)

    Nt = np.array(data[0][0]).shape[0]
    Nx = np.array(data[0][0]).shape[1]
    Ny = np.array(data[0][0]).shape[2]

    X = range(0,Nx)
    T = range(0,Nt)
    T,X = np.meshgrid(X, T)

    fig = plt.figure(figsize=figsize)
    p=2
    for j in range(num_cols):
        i=0
        if(x_3d):
            ax = plt.subplot2grid((num_rows, num_cols),(i,j), projection='3d')
            surf = ax.plot_surface(X[::p,::p],T[::p,::p], data[j][0][::p,Ny/2,::p],
                                cmap=cmap, edgecolors="k", alpha=0.9,  shade=False,
                                rstride=1, cstride=1, linewidth=0.0, antialiased=False)
            ax.set_title(data[j][1])
            ax.set_ylabel(r"$x(\theta)$")
            ax.set_xlabel(r"$\tau(ms)$")
            ax.set_zlabel(r"$W$")
            ax.view_init(elev=46., azim=130)
            if(colorbar):
                cbar = plt.colorbar(surf, ax=ax, orientation='horizontal')
            i+=1
        if(y_3d):
            ax = plt.subplot2grid((num_rows, num_cols),(i,j), projection='3d')
            surf = ax.plot_surface(X[::p,::p],T[::p,::p], data[j][0][::p,::p,Nx/2],
            cmap=cmap, edgecolors="k", alpha=0.9,  shade=False,
            rstride=1, cstride=1, linewidth=0.0, antialiased=False)
            ax.set_ylabel(r"$y(\theta)$")
            ax.set_xlabel(r"$\tau(ms)$")
            ax.set_zlabel(r"$W$")
            ax.view_init(elev=46., azim=130)
            if(colorbar):
                fig.colorbar(surf, ax=ax, orientation='horizontal')
            i+=1


    fig.tight_layout()
    if(save_figure):
        fig.save(animation_name + ".svg")
    plt.show()




def line3dPlotsOfImpulseResponses(data,
                                  x_line3d=True,
                                  y_line3d=False,
                                  idx=0,
                                  idy=0,
                                  figsize = (14,8),
                                  cmap = cmaps.inferno,
                                  colorbar = False,
                                  save_figure = False,
                                  figure_name = "unnamed"):

    """
    Imshow plots of impulse response functions
    """
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib.collections import PolyCollection
    num_cols = len(data)
    num_rows = int(x_line3d) + int(y_line3d)

    Nt = np.array(data[0][0]).shape[0]
    Nx = np.array(data[0][0]).shape[1]
    Ny = np.array(data[0][0]).shape[2]

    # fig, axarr = plt.subplots(num_rows, num_cols, figsize=figsize)
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    x_vec = np.linspace(-2.56/2, 256/2-1, 128)
    t_vec = np.linspace(0, 51.2, 256)

    zs = [0]
    for j in range(num_cols):
        # axarr[0,j].set_title(data[j][1])
        i=0
        if(x_line3d):
            verts = []
            verts.append(list(zip(t_vec[::100] , data[j][0][::100,0,idy]  )))
            # verts.append(list(zip(t_vec[::100] , data[j][0][::100,100,idy]  )))
            poly = PolyCollection(verts)
            ax.add_collection3d(poly , zdir='x')
            ax.set_ylabel(r"$x(\theta)$")
            ax.set_xlabel(r"$\tau(ms)$")



    plt.tight_layout()
    if(save_figure):
        fig.save(animation_name + ".svg")
    plt.show()



if __name__ == "__main__":
    import h5py
    from glob import glob
    import Simulation as sim

    outputFilePath =  "/home/milad/Dropbox/projects/lgn/code/lgn-simulator/apps/firingSynchrony/firingSynchrony.h5"
    outputFile = glob(outputFilePath)[0]
    f = h5py.File(outputFile, "r")
    exp = sim.Simulation(f)

    data = [
     [exp.stimulus.spatioTemporal, "Stimulus"]
    # ,[exp.ganglion.response["spatioTemporal"], "Ganglion cell response"]
    # ,[exp.ganglion.impulseResponse["spatioTemporal"], "Ganglion cell impulse response"]
    # ,[exp.interneuron.response["spatioTemporal"], "Interneuron"]
    # ,[exp.interneuron.impulseResponse["spatioTemporal"], "Interneuron"]
    ,[exp.relay.response["spatioTemporal"], "Relay cell response"]
    ,[exp.relay.impulseResponse["spatioTemporal"], "Relay cell impulse response"]
    ,[exp.cortical.response["spatioTemporal"], "Cortical cell response"]
    ,[exp.cortical.impulseResponse["spatioTemporal"], "Cortical impulse response"]
    ]

    idx = exp.integrator.nPointsSpatial/2
    idy = idx

    # animateImshowPlots(data, exp.integrator.temporalResolution, colorbar = True, save_animation = False, animation_name = "rat")
    # animate3dPlots(data, resolution = 3)

    data = [
    [exp.ganglion.response["spatioTemporal"], "Ganglion"]
    ,[exp.relay.response["spatioTemporal"], "Relay"]
    ,[exp.cortical.response["spatioTemporal"], "Cortical"]
    ]


    imshowPlotsOfImpulseResponses(data)
    # line3dPlotsOfImpulseResponses(data)
    # plot3dOfImpulseResponses(data, colorbar=True)

    plt.show()
