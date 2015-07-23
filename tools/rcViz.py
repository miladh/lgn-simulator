import h5py
from glob import glob
from pylab import*


import matplotlib.pyplot as plt
import matplotlib.animation as animation

# User commands:
cellType = "cortical"
outputFilePath = "/home/milad/Dropbox/projects/edog/extendedDOG/eDOG/DATA/*.h5"
# outputFilePath = "/home/milad/Dropbox/projects/edog/extendedDOG/eDOG/tools/fb/*.h5"
outputFile = glob(outputFilePath)[0]

f = h5py.File(outputFile, "r")
states = f.get("/")
dt = states.attrs["dt"]

R = []
G = []
S = []

for stateId, state in enumerate(states):
    dataset = states.get(state)
    S.append(array(dataset.get("stimuli/real")))
    R.append(array(dataset.get(cellType+"/response/real")).clip(min =0))
    G.append(array(dataset.get(cellType+"/impulseResponse/real")))

    nStates = len(states)
    print "Number of states: ", nStates

2


#####################################################################
def init():
    S_im.set_data(S[0])
    G_im.set_data(G[0])
    R_im.set_data(R[0])
    ttl.set_text("")
    return [S_im, G_im, R_im], ttl


def animate(i):
    S_im.set_array(S[i])
    G_im.set_array(G[i])
    R_im.set_array(R[i])
    ttl.set_text(cellType + "\n" +"t = " + str(i*dt) + " s")
#    print Rc[i].max(), " - ", Rc[i].min()
    return [S_im, G_im, R_im], ttl


fig, axarr = plt.subplots(1,3,figsize=(15,8))
tight_layout()
S_im=axarr[0].imshow(S[0], origin='lower', cmap='gray')
colorbar(S_im, ax = axarr[0], orientation='horizontal')

G_im=axarr[1].imshow(G[0], origin='lower', cmap='jet')
colorbar(G_im, ax = axarr[1], orientation='horizontal')

R_im=axarr[2].imshow(R[0], origin='lower', cmap='jet')
colorbar(R_im, ax = axarr[2], orientation='horizontal')

ttl = plt.suptitle("")
axarr[0].set_title("Stimuli")
axarr[1].set_title("Impulse Response")
axarr[2].set_title("Response")


anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=nStates, interval=20, blit=False)
show()
#anim.save('basic_animation.mp4',fps=30,  writer="avconv", codec="libx264")
