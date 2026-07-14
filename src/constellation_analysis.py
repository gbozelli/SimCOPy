import numpy as np
import matplotlib.pyplot as plt

def constellation_analysis(s):
    plt.figure()
    plt.plot(np.real(s),np.imag(s),'.')

    H, xedges, yedges = np.histogram2d(np.real(s),np.imag(s),100)
    #plt.figure()
    #plt.contourf(counts)

    fig = plt.figure()
    #ax = fig.add_subplot(131, title='imshow: square bins')
    plt.imshow(np.log10(H), interpolation='nearest', origin='lower', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],cmap="hot")

    fig = plt.figure()
    #ax = fig.add_subplot(131, title='imshow: square bins')
    plt.imshow(H,cmap="hot")

    fig = plt.figure()
    #ax = fig.add_subplot(131, title='imshow: square bins')
    plt.imshow(np.log10(H), interpolation='nearest', origin='lower', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],cmap="Blues")

    fig = plt.figure()
    #ax = fig.add_subplot(131, title='imshow: square bins')
    plt.imshow(H,cmap="Blues")

