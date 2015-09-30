import matplotlib.pyplot as plt
import matplotlib.animation as animation
import colormaps as cmaps
import numpy as np

def animateImshowPlots(data, figsize = (15,8), cmap = cmaps.viridis):
    num_subplots = len(data)
    imshowPlots = []
    nStates = data[0].shape[0]

    if(num_subplots ==1):
        num_cols = 1
        num_rows = 1
    elif(num_subplots ==2):
        num_cols = 2
        num_rows = 1
    else:
        num_cols = 3
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
        cellType = "No name"
        ttl.set_text(cellType + "\n" + "t = " + str(j) + " s")
        return imshowPlots, ttl


    ttl = plt.suptitle("")
    iplot = [num_rows, num_cols, 0]
    for i in range(num_subplots):
        iplot[2] += 1
        ax = fig.add_subplot(iplot[0], iplot[1], iplot[2])
        imshowPlots.append(ax.imshow(data[i][0,:,:], origin='lower', cmap=cmap,
        interpolation="None"))
        #plt.colorbar(data[i+j][0,:,:], ax=axarr[i,j], orientation='horizontal')


    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                   frames=nStates, interval=20, blit=False)
    plt.tight_layout()
    plt.show()
    # # anim.save('basic_animation.mp4',fps=30,  writer="avconv", codec="libx264")



if __name__ == "__main__":
    import h5py
    from glob import glob
    import Simulation as sim


    outputFilePath = "/home/milad/Dropbox/projects/edog/extendedDOG/eDOG/DATA/*.h5"
    outputFile = glob(outputFilePath)[0]
    f = h5py.File(outputFile, "r")
    exp = sim.Simulation(f)

    data = [exp.stimuli["spatioTemporal"], exp.stimuli["spatioTemporal"], exp.stimuli["spatioTemporal"], exp.stimuli["spatioTemporal"]
    ]
    animateImshowPlots(data)
