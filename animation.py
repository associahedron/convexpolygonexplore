import numpy as np
import matplotlib.pyplot as plt
from codewords import *
from sys import exit, argv


if __name__ == '__main__':
    if len(argv) < 3:
        print("Usage: python animation.py <N> <frames_per_step>")
        exit(0)
    N = int(argv[1])
    frames_per_step = int(argv[2])

    plt.figure()
    ax1 = plt.subplot2grid((1, 10), (0, 0))
    ax2 = plt.subplot2grid((1, 10), (0, 1))
    ax3 = plt.subplot2grid((1, 10), (0, 2), colspan=8)
    ax3.set_xticks([])
    ax3.set_yticks([])


    a = Associahedron(N, {"draw_tree":True}, ax=ax3)
    C = np.array([c["c"].w for c in a.codewords])
    S = np.array([c["s"] for c in a.codewords])[:, 1:-1]
    ax1.imshow(C, aspect='auto', interpolation='none')
    S = S%2
    SDisp = np.array(S)
    SDisp[SDisp == 1] = 2
    ax2.imshow(SDisp, aspect='auto', interpolation='none', cmap='gray', vmin=-2, vmax=2)
    ax1.set_xticks([])
    ax1.set_yticks([])
    ax2.set_xticks([])
    ax2.set_yticks([])

    ly = [0, 0]
    line1, = ax1.plot([0, C.shape[1]-1], ly, c='C1', zorder=200)
    line2, = ax2.plot([0, S.shape[1]-1], ly, c='C1', zorder=200)

    tcirc = np.linspace(0, 2*np.pi, 100)
    circ = np.zeros((tcirc.size, 2))
    circ[:, 0] = np.cos(tcirc)
    circ[:, 1] = np.sin(tcirc)
    scirc = ax3.scatter(circ[:, 0], circ[:, 1], zorder=400)

    pidx = 0
    x1 = 0
    y1 = 0
    i = 0

    NORMAL_TRAVERSAL = 0
    MOVE_TO_STACK = 1
    state = NORMAL_TRAVERSAL

    ax1.set_title("{}".format(a.codewords[0]["c"].w))

    while i < len(a.codewords):
        c = a.codewords[i]
        x2 = c["x"]
        y2 = c["y"]
        if i > 0 and not np.allclose(S[i-1], S[i]):
            d = np.where(S[i-1] != S[i])[0][-1]+1
            if state == NORMAL_TRAVERSAL:
                ax3.set_title("Moving to new stack of dim {}".format(d))
                state = MOVE_TO_STACK
                x2 += 1.5*d
            elif state == MOVE_TO_STACK:
                ax3.set_title("")
                state = NORMAL_TRAVERSAL
        line1.set_ydata([i-1, i-1])
        line2.set_ydata([i-1, i-1])
        for t in np.linspace(0, 1, frames_per_step):
            x = (1-t)*x1 + t*x2
            y = (1-t)*y1 + t*y2
            circ = np.zeros((tcirc.size, 2))
            circ[:, 0] = x + np.cos(tcirc)
            circ[:, 1] = y + np.sin(tcirc)
            scirc.set_offsets(circ)
            ax3.set_xlim(x-3, x+3)
            ax3.set_ylim(y-3, y+3)
            plt.savefig("Codeword{}.png".format(pidx))
            pidx += 1
        x1, y1 = x2, y2
        ax1.set_title("{}".format(c["c"].w))
        if state == NORMAL_TRAVERSAL:
            i += 1

