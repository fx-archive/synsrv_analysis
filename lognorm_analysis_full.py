
import sys, os, itertools, matplotlib, math
matplotlib.use('Agg')
import matplotlib.pyplot as pl
from matplotlib import rc

rc('text', usetex=True)
pl.rcParams['text.latex.preamble'] = [
    r'\usepackage{tgheros}',    # helvetica font
    r'\usepackage{sansmath}',   # math-font matching helvetica
    r'\sansmath'                # actually tell tex to use it!
    r'\usepackage{siunitx}',    # micro symbols
    r'\sisetup{detect-all}',    # force siunitx to use the fonts
]  

import numpy as np

from brian2.units import mV, ms, second

from pypet.trajectory import Trajectory
from pypet.brian2.parameter import Brian2Parameter, Brian2MonitorResult

from multiprocessing import Pool

from plot_methods import raster_plot, voltage_traces, isi_distribution, \
                         conductance_mult_trace, firing_rate_distribution_exc, \
                         firing_rate_distribution_inh, synapse_active_trace, \
                         synapse_weight_distribution_t, print_membrane_params, \
                         synapse_lifetime_distribution, Apre_traces, \
                         Apost_traces, membrane_threshold_distribution_t, \
                         membrane_threshold_traces, synapse_weight_traces, \
                         synapse_lifetime_distribution_loglog, \
                         synapse_lifetime_distribution_loglin, \
                         synapse_lifetime_distribution_loglog_linbin, \
                         synapse_weight_distribution_log_t, \
                         synapse_weight_distribution_log_t_faulty, \
                         synapse_weight_change_on_spike, \
                         synapse_weight_change_on_spike_log, \
                         synapse_weight_change_on_spike_loglog, \
                         print_netw_params, print_stdp_params, \
                         synapse_deathtime_distribution, ax_off, \
                         synapse_weight_distribution_loglog_t




def load_trajectory(fname):
    
    tr = Trajectory(name='tr1',
                    add_time=False, 
                    filename=fname,
                    dynamic_imports=[Brian2MonitorResult, Brian2Parameter])

    # pypet.pypetconstants.LOAD_NOTHING  --> 0
    # pypet.pypetconstants.LOAD_SKELETON --> 1
    # pypet.pypetconstants.LOAD_DATA     --> 2
    tr.f_load(load_parameters=2, load_derived_parameters=2, load_results=1)
    tr.v_auto_load = True

    return tr



def lognorm_analysis_figure(idx):

    global fname, ftitle
    
    tr = load_trajectory(fname)
    tr.v_idx = idx
    crun = tr.v_crun
    
    pl.close()
    fig = pl.figure()

    ax_lines, ax_cols = 6,8
    s = 150
    fig.set_size_inches(1920/s*ax_cols/4,1080/s*ax_lines/3)

    axs = {}
    for x,y in itertools.product(range(ax_lines),range(ax_cols)):
        axs['%d,%d'%(x+1,y+1)] = pl.subplot2grid((ax_lines, ax_cols), (x, y))

    midT = int(tr.T/tr.dt/2)

    # if tr.T/second > 2:
    raster_plot(axs['1,1'], tr, crun,
                tmin=0*second, tmax=2*second)
    raster_plot(axs['1,2'], tr, crun,
                tmin=tr.T/2-1*second, tmax=tr.T/2+1*second)
    raster_plot(axs['1,3'], tr, crun,
                tmin=tr.T-2*second, tmax=tr.T)

    voltage_traces(axs['2,1'], tr, crun, tmin=0*second, tmax=0.5*second)
    voltage_traces(axs['2,2'], tr, crun, tmin=tr.T/2-0.25*second,
                tmax=tr.T/2+0.25*second)
    voltage_traces(axs['2,3'], tr, crun, tmin=tr.T-0.25*second,
                   tmax=tr.T)

    print_netw_params(axs['1,4'], tr, crun)
    print_stdp_params(axs['2,4'], tr, crun)
    
    firing_rate_distribution_exc(axs['3,1'], tr, crun, steps=10)
    firing_rate_distribution_inh(axs['3,2'], tr, crun, steps=15)

    isi_distribution(axs['3,3'], tr, crun, bins=50)
    ax_off(axs['3,4'])

    synapse_weight_traces(axs['4,1'], tr, crun,
                          tmin=0*second, tmax=2*second)
    synapse_weight_traces(axs['4,2'], tr, crun,
                          tmin=tr.T/2-1*second, tmax=tr.T/2+1*second)
    synapse_weight_traces(axs['4,3'], tr, crun,
                          tmin=tr.T-2*second, tmax=tr.T)
    synapse_weight_traces(axs['4,4'], tr, crun)


    # conductance_mult_trace(axs['3,1'], tr, crun, tmin=0,tmax=2000)

    # conductance_mult_trace(axs['3,2'], tr, crun, tmin=midT, tmax=midT+2000)
    # conductance_mult_trace(axs['3,3'], tr, crun, tmin=-2001,tmax=-1)


    
    synapse_weight_distribution_t(axs['5,1'], tr, crun, tstep=0)
    synapse_weight_distribution_t(axs['5,2'], tr, crun, tstep=1)
    synapse_weight_distribution_t(axs['5,3'], tr, crun, tstep=2)
    synapse_weight_distribution_t(axs['5,4'], tr, crun, tstep=-1)

    synapse_weight_distribution_log_t(axs['6,1'], tr, crun, tstep=1)
    synapse_weight_distribution_log_t(axs['6,2'], tr, crun, tstep=2)
    synapse_weight_distribution_log_t(axs['6,3'], tr, crun, tstep=3)
    synapse_weight_distribution_log_t(axs['6,4'], tr, crun, tstep=-1)

    synapse_weight_distribution_log_t(axs['1,5'], tr, crun, tstep=1, bins=50)
    synapse_weight_distribution_log_t(axs['1,6'], tr, crun, tstep=2, bins=50)
    synapse_weight_distribution_log_t(axs['1,7'], tr, crun, tstep=3, bins=50)
    synapse_weight_distribution_log_t(axs['1,8'], tr, crun, tstep=-1, bins=50)

    synapse_weight_distribution_log_t(axs['2,5'], tr, crun, tstep=1, bins=100)
    synapse_weight_distribution_log_t(axs['2,6'], tr, crun, tstep=2, bins=100)
    synapse_weight_distribution_log_t(axs['2,7'], tr, crun, tstep=3, bins=100)
    synapse_weight_distribution_log_t(axs['2,8'], tr, crun, tstep=-1, bins=100)

    synapse_weight_distribution_log_t(axs['3,5'], tr, crun, tstep=1, bins=250)
    synapse_weight_distribution_log_t(axs['3,6'], tr, crun, tstep=2, bins=250)
    synapse_weight_distribution_log_t(axs['3,7'], tr, crun, tstep=3, bins=250)
    synapse_weight_distribution_log_t(axs['3,8'], tr, crun, tstep=-1, bins=250)

    synapse_weight_distribution_loglog_t(axs['4,5'], tr, crun, tstep=1, bins=50)
    synapse_weight_distribution_loglog_t(axs['4,6'], tr, crun, tstep=2, bins=50)
    synapse_weight_distribution_loglog_t(axs['4,7'], tr, crun, tstep=3, bins=50)
    synapse_weight_distribution_loglog_t(axs['4,8'], tr, crun, tstep=-1, bins=50)

    synapse_weight_distribution_loglog_t(axs['5,5'], tr, crun, tstep=1, bins=100)
    synapse_weight_distribution_loglog_t(axs['5,6'], tr, crun, tstep=2, bins=100)
    synapse_weight_distribution_loglog_t(axs['5,7'], tr, crun, tstep=3, bins=100)
    synapse_weight_distribution_loglog_t(axs['5,8'], tr, crun, tstep=-1, bins=100)

    synapse_weight_distribution_loglog_t(axs['6,5'], tr, crun, tstep=1, bins=250)
    synapse_weight_distribution_loglog_t(axs['6,6'], tr, crun, tstep=2, bins=250)
    synapse_weight_distribution_loglog_t(axs['6,7'], tr, crun, tstep=3, bins=250)
    synapse_weight_distribution_loglog_t(axs['6,8'], tr, crun, tstep=-1, bins=250)
    
    
    
    # synapse_weight_distribution_t(axs['4,2'], tr, crun, tstep=-1)
    # synapse_weight_distribution_log_t(axs['4,3'], tr, crun, bins=30,
    #                                  tstep=-1)
    # synapse_weight_distribution_log_t(axs['4,4'], tr, crun, bins=30,
    #                                   a_min=tr.a_insert, tstep=-1)    
    # synapse_weight_distribution_log_t_faulty(axs['4,5'], tr, crun, bins=30)
        
    
    # synapse_active_trace(axs['5,1'], tr, crun)
    # synapse_deathtime_distribution(axs['5,2'], tr, crun)
    # synapse_lifetime_distribution(axs['5,3'], tr, crun)
    # synapse_lifetime_distribution_loglog_linbin(axs['5,4'], tr, crun)
    # synapse_lifetime_distribution_loglog(axs['5,5'], tr, crun)


    # synapse_weight_traces(axs['6,1'], tr, crun)
    # synapse_weight_traces(axs['7,1'], tr, crun, tmin=tr.T-1*second, tmax=tr.T)
    # synapse_weight_change_on_spike(axs['6,2'], tr, crun)
    # synapse_weight_change_on_spike_log(axs['6,3'], tr, crun)
    # synapse_weight_change_on_spike_loglog(axs['6,4'], tr, crun)

    # synapse_lifetime_distribution_loglin(axs['6,5'], tr, crun)
    
    # # Apre_traces(axs['6,4'], tr, crun, tmin=-2001, tmax=-1)
    # # Apost_traces(axs['6,5'], tr, crun, tmin=-2001, tmax=-1)

    # membrane_threshold_distribution_t(axs['7,2'], tr, crun, tstep=0)
    # membrane_threshold_distribution_t(axs['7,3'], tr, crun, tstep=-1)
    # membrane_threshold_traces(axs['7,4'], tr, crun)
    # print_membrane_params(axs['7,5'], tr, crun)
            
    pl.tight_layout()

    subdir = "analysis/lognorm_full".format(ftitle)
    if not os.path.exists(subdir):
        os.makedirs(subdir)
        
    pl.savefig(subdir+"/log_norm_full_{:s}_{:s}.png".format(ftitle, crun),
               dpi=100, bbox_inches='tight')



    
if __name__ == "__main__":

    fname = sys.argv[1]
    ftitle = os.path.basename(os.path.splitext(fname)[0])
    
    n_pool, run_pool =  int(sys.argv[2]), bool(int(sys.argv[3]))

    tr = load_trajectory(fname)
    num_runs = len(list(tr.f_iter_runs()))

    if run_pool:
        print("Using pool")
        
        pool = Pool(n_pool)
        pool.map(lognorm_analysis_figure, range(num_runs))
        pool.close()
        pool.join()

    else:
        print("Not using multiprocessing.")
        for idx in range(num_runs):
        # for idx in range(1):
            lognorm_analysis_figure(idx)
    
        
