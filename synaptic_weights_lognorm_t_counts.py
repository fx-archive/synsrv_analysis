
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

    synapse_weight_distribution_log_t(axs['1,1'], tr, crun, tstep=1, bins=50, density=False)
    synapse_weight_distribution_log_t(axs['1,2'], tr, crun, tstep=2, bins=50, density=False)
    synapse_weight_distribution_log_t(axs['1,3'], tr, crun, tstep=3, bins=50, density=False)
    synapse_weight_distribution_log_t(axs['1,4'], tr, crun, tstep=4, bins=50, density=False)
    synapse_weight_distribution_log_t(axs['1,5'], tr, crun, tstep=5, bins=50, density=False)
    synapse_weight_distribution_log_t(axs['1,6'], tr, crun, tstep=6, bins=50, density=False)
    synapse_weight_distribution_log_t(axs['1,7'], tr, crun, tstep=-1, bins=50, density=False)

    fit=True
    synapse_weight_distribution_log_t(axs['2,1'], tr, crun, tstep=1, bins=50,
                                      low_bound=0.00000000000000001*tr.ATotalMax, fit=fit, density=False)
    synapse_weight_distribution_log_t(axs['2,2'], tr, crun, tstep=2, bins=50,
                                      low_bound=0.00000000000000001*tr.ATotalMax, fit=fit, density=False)
    synapse_weight_distribution_log_t(axs['2,3'], tr, crun, tstep=3, bins=50,
                                      low_bound=0.00000000000000001*tr.ATotalMax, fit=fit, density=False)
    synapse_weight_distribution_log_t(axs['2,4'], tr, crun, tstep=4, bins=50,
                                      low_bound=0.00000000000000001*tr.ATotalMax, fit=fit, density=False)
    synapse_weight_distribution_log_t(axs['2,5'], tr, crun, tstep=5, bins=50,
                                      low_bound=0.00000000000000001*tr.ATotalMax, fit=fit, density=False)
    synapse_weight_distribution_log_t(axs['2,6'], tr, crun, tstep=6, bins=50,
                                      low_bound=0.00000000000000001*tr.ATotalMax, fit=fit, density=False)
    synapse_weight_distribution_log_t(axs['2,7'], tr, crun, tstep=-1, bins=50,
                                      low_bound=0.00000000000000001*tr.ATotalMax, fit=fit, density=False)

    synapse_weight_distribution_log_t(axs['3,1'], tr, crun, tstep=1, bins=50,
                                      low_bound=0.00000000001*tr.ATotalMax, fit=fit, density=False)
    synapse_weight_distribution_log_t(axs['3,2'], tr, crun, tstep=2, bins=50,
                                      low_bound=0.00000000001*tr.ATotalMax, fit=fit, density=False)
    synapse_weight_distribution_log_t(axs['3,3'], tr, crun, tstep=3, bins=50,
                                      low_bound=0.00000000001*tr.ATotalMax, fit=fit, density=False)
    synapse_weight_distribution_log_t(axs['3,4'], tr, crun, tstep=4, bins=50,
                                      low_bound=0.00000000001*tr.ATotalMax, fit=fit, density=False)
    synapse_weight_distribution_log_t(axs['3,5'], tr, crun, tstep=5, bins=50,
                                      low_bound=0.00000000001*tr.ATotalMax, fit=fit, density=False)
    synapse_weight_distribution_log_t(axs['3,6'], tr, crun, tstep=6, bins=50,
                                      low_bound=0.00000000001*tr.ATotalMax, fit=fit, density=False)
    synapse_weight_distribution_log_t(axs['3,7'], tr, crun, tstep=-1, bins=50,
                                      low_bound=0.00000000001*tr.ATotalMax, fit=fit, density=False)
    

    synapse_weight_distribution_log_t(axs['4,1'], tr, crun, tstep=1, bins=50,
                                      low_bound=0.00000001*tr.ATotalMax, fit=fit, density=False)
    synapse_weight_distribution_log_t(axs['4,2'], tr, crun, tstep=2, bins=50,
                                      low_bound=0.00000001*tr.ATotalMax, fit=fit, density=False)
    synapse_weight_distribution_log_t(axs['4,3'], tr, crun, tstep=3, bins=50,
                                      low_bound=0.00000001*tr.ATotalMax, fit=fit, density=False)
    synapse_weight_distribution_log_t(axs['4,4'], tr, crun, tstep=4, bins=50,
                                      low_bound=0.00000001*tr.ATotalMax, fit=fit, density=False)
    synapse_weight_distribution_log_t(axs['4,5'], tr, crun, tstep=5, bins=50,
                                      low_bound=0.00000001*tr.ATotalMax, fit=fit, density=False)
    synapse_weight_distribution_log_t(axs['4,6'], tr, crun, tstep=6, bins=50,
                                      low_bound=0.00000001*tr.ATotalMax, fit=fit, density=False)
    synapse_weight_distribution_log_t(axs['4,7'], tr, crun, tstep=-1, bins=50,
                                      low_bound=0.00000001*tr.ATotalMax, fit=fit, density=False)

    synapse_weight_distribution_log_t(axs['5,1'], tr, crun, tstep=1, bins=50,
                                      low_bound=0.00001*tr.ATotalMax, fit=fit, density=False)
    synapse_weight_distribution_log_t(axs['5,2'], tr, crun, tstep=2, bins=50,
                                      low_bound=0.00001*tr.ATotalMax, fit=fit, density=False)
    synapse_weight_distribution_log_t(axs['5,3'], tr, crun, tstep=3, bins=50,
                                      low_bound=0.00001*tr.ATotalMax, fit=fit, density=False)
    synapse_weight_distribution_log_t(axs['5,4'], tr, crun, tstep=4, bins=50,
                                      low_bound=0.00001*tr.ATotalMax, fit=fit, density=False)
    synapse_weight_distribution_log_t(axs['5,5'], tr, crun, tstep=5, bins=50,
                                      low_bound=0.00001*tr.ATotalMax, fit=fit, density=False)
    synapse_weight_distribution_log_t(axs['5,6'], tr, crun, tstep=6, bins=50,
                                      low_bound=0.00001*tr.ATotalMax, fit=fit, density=False)
    synapse_weight_distribution_log_t(axs['5,7'], tr, crun, tstep=-1, bins=50,
                                      low_bound=0.00001*tr.ATotalMax, fit=fit, density=False)

    synapse_weight_distribution_log_t(axs['6,1'], tr, crun, tstep=1, bins=50,
                                      low_bound=0.001*tr.ATotalMax, fit=fit, density=False)
    synapse_weight_distribution_log_t(axs['6,2'], tr, crun, tstep=2, bins=50,
                                      low_bound=0.001*tr.ATotalMax, fit=fit, density=False)
    synapse_weight_distribution_log_t(axs['6,3'], tr, crun, tstep=3, bins=50,
                                      low_bound=0.001*tr.ATotalMax, fit=fit, density=False)
    synapse_weight_distribution_log_t(axs['6,4'], tr, crun, tstep=4, bins=50,
                                      low_bound=0.001*tr.ATotalMax, fit=fit, density=False)
    synapse_weight_distribution_log_t(axs['6,5'], tr, crun, tstep=5, bins=50,
                                      low_bound=0.001*tr.ATotalMax, fit=fit, density=False)
    synapse_weight_distribution_log_t(axs['6,6'], tr, crun, tstep=6, bins=50,
                                      low_bound=0.001*tr.ATotalMax, fit=fit, density=False)
    synapse_weight_distribution_log_t(axs['6,7'], tr, crun, tstep=-1, bins=50,
                                      low_bound=0.001*tr.ATotalMax, fit=fit, density=False)

    
    print_netw_params(axs['1,8'], tr, crun)
    print_stdp_params(axs['2,8'], tr, crun)
    ax_off(axs['3,8'])
    ax_off(axs['4,8'])
    ax_off(axs['5,8'])
    ax_off(axs['6,8'])

            
    pl.tight_layout()

    subdir = "analysis/synaptic_weight_lognorm_t_counts"
    if not os.path.exists(subdir):
        os.makedirs(subdir)
        
    pl.savefig(subdir+"/syn_lnt_{:s}_{:s}.png".format(ftitle, crun),
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
    
        
