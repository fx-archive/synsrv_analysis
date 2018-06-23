
import numpy as np

from brian2.units import mV, ms, second


def conductance_mult_trace(ax, tr, crun='run_00000000', tmin=0*second, tmax='None'):

    df = tr.crun.GExc_stat

    if tmax=='None':
        tmax=tr.T

    js = np.array((range(len(df.t))))
    j_tmin, j_tmax = js[df.t>=(tmin/second)][0], js[df.t<=(tmax/second)][-1]
    
    ge_max, gi_min = [],[]
    
    for i in df.record:
        gi_min.append(1.1*np.min(df.gi[i,j_tmin:j_tmax]))
        ax.plot(df.t[j_tmin:j_tmax],
                np.sum(ge_max)-np.sum(gi_min)+df.ge[i,j_tmin:j_tmax],
                color='blue')
        ax.plot(df.t[j_tmin:j_tmax],
                np.sum(ge_max)-np.sum(gi_min)+df.gi[i,j_tmin:j_tmax],
                color='red')

        ge_max.append(1.1*np.max(df.ge[i,j_tmin:j_tmax]))

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_title('Conductance Traces')
    ax.set_xlabel('time [s]')
    ax.set_ylabel('conductance')
