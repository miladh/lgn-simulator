import h5py
from glob import glob
from pylab import*


import matplotlib.pyplot as plt
import matplotlib.animation as animation

outputFilePath = "/home/milad/Dropbox/projects/edog/extendedDOG/eDOG/DATA/*.h5"
outputFile = glob(outputFilePath)[0]

f = h5py.File(outputFile, "r")             
states = f.get("/")

dt = states.attrs["dt"]


Rc = [] 
for stateId, state in enumerate(states):
    dataset = states.get(state)
    Rc.append(array(dataset.get("complex")))

nStates = len(states)
print "Number of states: ", nStates

f.close()


#####################################################################

def init():
    im.set_data(Rc[0])
    return [im]


def animate(i):
    im.set_array(Rc[i])
    plt.title("t= " + str(i*dt) + " s")
#    print Rc[i].max(), " - ", Rc[i].min() 
    return [im]


fig = plt.figure()
im=plt.imshow(Rc[0], origin='lower', cmap='jet')

colorbar()

anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=nStates, interval=20, blit=True)

anim.save('basic_animation.mp4',fps=30,  writer="avconv", codec="libx264")