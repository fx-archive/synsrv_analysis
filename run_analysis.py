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

from post_processing import extract_lifetimes


def life_time_distribution_compare(tr, bin_w, ftitle):

    print(''' Not discarding any ts. Need to fix or reconsider,
              as early (t<5s) synapse dynamics may distort 
              lifetime distribution ''')
    
    pl.close()
    fig = pl.figure()
    fig.set_size_inches(12,9)
    
    ax = fig.add_subplot(111)

    for crun in tr.f_iter_runs():

        _lt, _dt = extract_lifetimes(tr.crun.turnover, tr.N_e)
        life_t, death_t = _lt*second, _dt*second

        counts, edges = np.histogram(life_t/ms,
                                     bins=np.arange(tr.dt/ms,tr.T/ms,bin_w))

        centers = (edges[:-1] + edges[1:])/2.
        ax.plot(centers, counts, '.', markersize=3., label=tr.taupre/ms)

    ax.legend(loc='lower left')

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_title('Lifetime distribution')
    ax.set_xlabel('time [ms]')

    directory = "default_output/{:s}".format(ftitle)
    if not os.path.exists(directory):
        os.makedirs(directory)
        
    pl.savefig(directory+"/{:s}__ltcompare.png".format(ftitle),
               dpi=100, bbox_inches='tight')




if __name__ == "__main__":

    fname = sys.argv[1]
    ftitle = os.path.basename(os.path.splitext(fname)[0])

    tr = Trajectory(name='tr1',
                    add_time=False, 
                    filename=fname,
                    dynamic_imports=[Brian2MonitorResult, Brian2Parameter])

    # pypet.pypetconstants.LOAD_NOTHING  --> 0
    # pypet.pypetconstants.LOAD_SKELETON --> 1
    # pypet.pypetconstants.LOAD_DATA     --> 2
    tr.f_load(load_parameters=2, load_derived_parameters=2, load_results=1)
    tr.v_auto_load = True

    life_time_distribution_compare(tr, 10, ftitle)
    
