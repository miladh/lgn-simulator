import matplotlib.pyplot as plt
import matplotlib.animation as animation




def init():
    S_im.set_data(S[0,:,:])
    G_im.set_data(G[0,:,:])
    R_im.set_data(R[0,:,:])
    ttl.set_text("")
    return [S_im, R_im, G_im], ttl


def animate(i):
    S_im.set_array(S[i,:,:])
    G_im.set_array(G[i,:,:])
    R_im.set_array(R[i,:,:])
    ttl.set_text(cellType + "\n" + "t = " + str(i) + " s")
#    print Rc[i].max(), " - ", Rc[i].min()
    return [S_im, R_im, G_im], ttl




fig, axarr = plt.subplots(1, 3, figsize=(15, 8))
tight_layout()
S_im = axarr[0].imshow(S[0,:,:], origin='lower', cmap="gray", interpolation="None")
colorbar(S_im, ax=axarr[0], orientation='horizontal')

#cmaps.plasma
R_im = axarr[1].imshow(R[0,:,:], origin='lower', cmap=cmaps.viridis, interpolation="None")
colorbar(R_im, ax=axarr[1], orientation='horizontal', norm=mpl.colors.Normalize(vmin=-10, vmax=10))

G_im = axarr[2].imshow(G[1,:,:], origin='lower', cmap=cmaps.viridis)  # NOTE G[1]!!!!!NOTE
colorbar(G_im, ax=axarr[2], orientation='horizontal')

ttl = plt.suptitle("")
axarr[0].set_title("Stimuli")
axarr[1].set_title("Response")
axarr[2].set_title("Impulse Response")


anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=nStates, interval=20, blit=False)
show()
# # anim.save('basic_animation.mp4',fps=30,  writer="avconv", codec="libx264")
