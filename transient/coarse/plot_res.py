from func import *
import matplotlib.pyplot as plt


def flip_dict(d, x_var=''):
    k = list(d[0].keys())
    if x_var == '':
        x_var = k[0]

    out = {}
    for ki in k:
        out[ki] = []
    
    for di in d:
        for ki in k:
            out[ki].append( di[ki] )
    
    return out


def plot_all():
    d = get_perf(verb=False, time='all')
    d = flip_dict(d)

    fig, ax = plt.subplots(2)

    ax[0].plot(d['time'], d['thrust_ratio'])
    ax[0].set_xlabel('time (s)')
    ax[0].set_ylabel('thrust_ratio')

    ax[1].plot(d['time'], d['bypass_ratio'])
    ax[1].set_xlabel('time (s)')
    ax[1].set_ylabel('bypass_ratio')

    plt.show()
    


if __name__=='__main__':
    plot_all()