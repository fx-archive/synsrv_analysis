
import numpy as np
from scipy.stats import norm, lognorm
#from decimal import Decimal

from brian2.units import mV, ms, second




def n_active_synapses(ax, tr, crun='run_00000000'):

    df = tr.crun.SynEE_a

    active_at_t = np.sum(df['syn_active'], axis=0)
    all_at_t = np.shape(df['syn_active'])[0]

    assert all_at_t == tr.N_e*(tr.N_e-1)

    ax.plot(df.t, active_at_t/all_at_t, lw=2)

    # ax.set_title('E'+r'$\to$'+'E weights at t=\SI{'+ \
    #              str(df.t[tstep])+'}{s}')

    ax.set_xlabel('time [s]')
    ax.set_ylabel('fraction of synapses active')
    
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')



def weight_distribution_t(ax, tr, crun='run_00000000', bins=50,
                          tstep=0, low_bound=-1):

    df = tr.crun.SynEE_a

    weight_at_t = df['a'][:,tstep]
    states_at_t = df['syn_active'][:,tstep]

    active_weight_at_t = weight_at_t[states_at_t==1]

    ax.hist(active_weight_at_t, bins=bins)
    
    ax.set_title('E'+r'$\to$'+'E weights at t=\SI{'+ \
                 str(df.t[tstep])+'}{s}')

    ax.set_xlabel('synaptic weight')
    
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')



def weight_distribution_log_t(ax, tr, crun='run_00000000', bins=50,
                              tstep=0, low_bound=-1, fit=True):

    df = tr.crun.SynEE_a

    weight_at_t = df['a'][:,tstep]
    states_at_t = df['syn_active'][:,tstep]

    active_weight_at_t = weight_at_t[states_at_t==1]

    log_weights = np.log(active_weight_at_t)
    
    ax.hist(log_weights, bins=bins, density=True)

    if fit:
        floc, fscale = norm.fit(log_weights)
        f_rv = norm(loc=floc, scale=fscale)
        xs = np.linspace(start=np.min(log_weights),
                         stop=np.max(log_weights),
                         num = 1000)
        ax.plot(xs, f_rv.pdf(xs), lw=2, color='red',
                linestyle='-')
    
    ax.set_title('E'+r'$\to$'+'E weights at t=\SI{'+ \
                 str(df.t[tstep])+'}{s}')

    ax.set_xlabel('log of synaptic weight')

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')


    

def synapse_weight_distribution_log_t(ax, tr, crun='run_00000000',
                                      bins=30, a_min=-1, a_max=-1,
                                      tstep=-1, low_bound = -1,
                                      fit=False, fit_adapt_ylim=True,
                                      ax_text=False, density=True):

    df = tr.crun.SynEE_a
    data = df['a'][:,tstep]
    a_vals = data[df['syn_active'][:,tstep]==1]

    a_cut = a_vals[a_vals>0]

    if low_bound!=-1:
        a_cut = a_cut[a_cut>low_bound]

    # if a_min==-1:
    #     a_min = tr.prn_thrshld
    # if a_max==-1:
    #     try:
    #         a_max = np.max(a_cut)
    #     except ValueError:
    #         a_max = 1

    a_min = np.min(a_cut)
    a_max = np.max(a_cut)

    val, bins, _ = ax.hist(a_cut, 10**np.linspace(np.log10(a_min),
                                              np.log10(a_max),bins),
                           density=density)

    if fit:
        fs, floc, fscale = lognorm.fit(a_cut, floc=0)
        f_rv = lognorm(fs, loc=0, scale=fscale)
        xs = np.logspace(start=np.log10(a_min),
                         stop=np.log10(a_max),
                         base=10., num=5000)
        ax.plot(xs, f_rv.pdf(xs))

        if fit_adapt_ylim:
            ax.set_ylim(0, 1.2*np.max(val))
        

    ax.set_xscale('log')

    # ax.axvline(tr.a_insert, color='red')
    # ax.axvline(tr.prn_thrshld, color='red')

    ax.set_title('E-E weights at t=\SI{'+ \
                 str(df.t[tstep])+'}{s}')

    if ax_text:
        
        ax.text(1.0, 0.875,
                r'view thrs '+'%.2E' % Decimal(low_bound/tr.ATotalMax) + \
                r'$\, a_{\mathrm{TotalMax}}$',
                horizontalalignment='right',
                verticalalignment='top',
                bbox={'boxstyle': 'square, pad=0.3', 'facecolor':'white',
                      'alpha':1, 'edgecolor':'none'},
                transform = ax.transAxes)

    
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')

def view_thrshld_ax_text(ax, tr, low_bound, tstep=-1):

    df = tr.crun.SynEE_a
    data = df['a'][:,tstep]
    a_vals = data[df['syn_active'][:,tstep]==1]

    a_zeroplus = a_vals[a_vals>0]

    if low_bound!=-1:
        a_cut = a_zeroplus[a_zeroplus>low_bound]

    percent_shown = len(a_cut)/len(a_vals)
    
    ax.text(0.0, 0.75,
            r'view thrs '+'%.2E' % Decimal(low_bound/tr.ATotalMax) + \
            r'$\, a_{\mathrm{TotalMax}}$' + \
            '\n'+'%f of weights shown' % percent_shown + \
            '\n'+'%f above zero' % (len(a_zeroplus)/len(a_vals)),
            horizontalalignment='left',
            verticalalignment='top',
            bbox={'boxstyle': 'square, pad=0.3', 'facecolor':'white',
                  'alpha':1, 'edgecolor':'none'},
            transform = ax.transAxes, fontsize=16)

    ax.axis('off')      
    

def synapse_weight_distribution_loglog_t(ax, tr, crun='run_00000000',
                                         bins=30, a_min=-1, a_max=-1,
                                         tstep=-1, low_bound=-1):

    df = tr.crun.SynEE_a
    data = df['a'][:,tstep]
    a_vals = data[df['syn_active'][:,tstep]==1]

    a_cut = a_vals[a_vals>0]

    if low_bound!=-1:
        a_cut = a_cut[a_cut>low_bound]


    a_min = np.min(a_cut)
    a_max = np.max(a_cut)

    ax.hist(a_cut, 10**np.linspace(np.log10(a_min),
                                   np.log10(a_max),bins),
            density=True)

    ax.set_xscale('log')
    ax.set_yscale('log')

    # ax.axvline(tr.a_insert, color='red')
    # ax.axvline(tr.prn_thrshld, color='red')

    ax.set_title('SynEE Weights at t=\SI{'+ \
                 str(df.t[tstep])+'}{s}')
    
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')



def synapse_weight_distribution_log_t_faulty(ax, tr, crun='run_00000000',
                                             bins=30, a_min=-1, a_max=-1,
                                             tstep=-1):

    df = tr.crun.SynEE_a
    a_vals = df['a'][:,tstep]
    a_cut = a_vals[a_vals>tr.prn_thrshld]

    if a_min==-1:
        a_min = tr.prn_thrshld
    if a_max==-1:
        try:
            a_max = np.max(a_cut)
        except ValueError:
            a_max = 1

    ax.hist(np.log10(a_cut), bins=25, normed=True, color='gray')
  
    ax.axvline(np.log10(tr.a_insert), color='red')
    ax.axvline(np.log10(tr.prn_thrshld), color='red')

    ax.set_title('SynEE Weights at t=\SI{'+ \
                 str(df.t[tstep])+'}{s} (invalid)')


def synapse_active_trace(ax, tr,  crun='run_00000000'):

    if not len(tr.crun.turnover)==0:
        _asc_t, asc_n = extract_active_synapse_count(tr.crun.turnover)
        asc_t = _asc_t*second

        ax.plot(asc_t/ms,asc_n/(tr.N_e*tr.N_e-1))
        
    ax.set_xlabel('time [ms]')
    ax.set_ylabel('frac. of syn. active')
    ax.set_title('active synapses')
    
    
def synapse_weight_traces(ax, tr, crun='run_00000000', tmin=0*ms, tmax=-1*ms,
                          plot_thresholds=False):

    if tmax==-1*ms:
        tmax = tr.T

    df = tr.crun.SynEE_stat

    if df == []:
        pass
    
    else:    
    
        for i in range(len(df.a)):
            indx = np.logical_and(df.t > tmin/second,
                                  df.t < tmax/second)
            ax.plot(df.t[indx],df.a[i,:][indx], color='grey')
            ax.text(0.47, 0.95,
            'ascale='+'%.2E' % Decimal(tr.ascale),
            horizontalalignment='left',
            verticalalignment='top',
            bbox={'boxstyle': 'square, pad=0.3', 'facecolor':'white',
                  'alpha':1, 'edgecolor':'none'},
            transform = ax.transAxes)

        if plot_thresholds:
            ax.axhline(tr.a_insert, linestyle='dashed', color='grey')
            ax.axhline(tr.prn_thrshld, linestyle='dashed', color='grey')
            ax.axhline(tr.amax, color='red')

        ax.set_title('Synaptic Weight Traces')
        ax.set_xlabel('time [s]')

        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')


        
def synapse_lifetime_distribution(ax, tr, crun='run_00000000', discard_t=0.):
    ''' 
    discard all values until discard_t
    '''
    if discard_t!=0.:
        raise NotImplementedError

    if not len(tr.crun.turnover) == 0:
        _lt, _dt = extract_lifetimes(tr.crun.turnover, tr.N_e)
        life_t, death_t = _lt*second, _dt*second

        ax.hist(life_t/ms)

    ax.set_title('Lifetime distribution')
    ax.set_xlabel('time [ms]')

def synapse_deathtime_distribution(ax, tr, crun, discard_t=0.):
    ''' 
    discard all values until discard_t
    '''
    if discard_t!=0.:
        raise NotImplementedError

    if not len(tr.crun.turnover) == 0:
        _lt, _dt = extract_lifetimes(tr.crun.turnover, tr.N_e)
        life_t, death_t = _lt*second, _dt*second
        ax.hist(death_t/ms)
        
    ax.set_title('Deathtime distribution')
    ax.set_xlabel('time [ms]')    
    

def synapse_lifetime_distribution_loglog(ax, tr, crun='run_00000000',bins=50,
                                         discard_t=0.):
    ''' 
    discard all values until discard_t
    '''
    if discard_t!=0.:
        raise NotImplementedError
    else:
        print("not discarding any ts")

    if not len(tr.crun.turnover) == 0:
        _lt, _dt = extract_lifetimes(tr.crun.turnover, tr.N_e)
        life_t, death_t = _lt*second, _dt*second

        if len(life_t) == 0:
            print('No recorded lifetimes, not plotting distribution')
            ax.set_title('No recorded lifetimes')
        else:
            b_min, b_max = tr.dt/ms, np.max(life_t/ms)
            bins = np.linspace(np.log10(b_min), np.log10(b_max), bins)
            ax.hist(life_t/ms, 10**bins, log=True, normed=True)
            ax.set_title('Lifetime distribution')
            ax.set_xlabel('time [ms]')
            ax.set_xscale('log')
            ax.set_xlim(b_min,b_max)


def synapse_lifetime_distribution_loglog_linbin(ax, tr, crun='run_00000000',
                                                bin_w=10., discard_t=0.):
    ''' 
    discard all values until discard_t
    '''
    if discard_t!=0.:
        raise NotImplementedError
    else:
        print("not discarding any ts")

    if not len(tr.crun.turnover) == 0:
    
        _lt, _dt = extract_lifetimes(tr.crun.turnover, tr.N_e)
        life_t, death_t = _lt*second, _dt*second

        counts, edges = np.histogram(life_t/ms,
                                     bins=np.arange(tr.dt/ms,tr.T/ms,bin_w))
        centers = (edges[:-1] + edges[1:])/2.
        ax.plot(centers, counts, '.', markersize=0.5)

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_title('Lifetime distribution')
    ax.set_xlabel('time [ms]')


def synapse_lifetime_distribution_loglin(ax, tr, crun='run_00000000',bins=50,                                            discard_t=0.):
    ''' 
    discard all values until discard_t
    '''
    if discard_t!=0.:
        raise NotImplementedError

    if not len(tr.crun.turnover) == 0:
        _lt, _dt = extract_lifetimes(tr.crun.turnover, tr.N_e)
        life_t, death_t = _lt*second, _dt*second

        if len(life_t) == 0:
            print("No recorded lifetimes. Not plotting.")
            ax.set_title("No recorded lifetimes")
        else:
            b_min, b_max = tr.dt/ms, np.max(life_t/ms)
            bins = np.linspace(np.log10(b_min), np.log10(b_max), bins)
            ax.hist(life_t/ms, 10**bins, normed=False, histtype='step')
            ax.set_title('Lifetime distribution')
            ax.set_xlabel('time [ms]')
            ax.set_xscale('log')
            ax.set_xlim(b_min,b_max)

def get_dA_prepost_minmax(dA_on_pre, dA_on_post):
    '''
    returns 
    '''
    try:
        dA_on_pre_min = np.min(dA_on_pre)
    except ValueError:
        dA_on_pre_min = 0.

    try:
        dA_on_pre_max = np.max(dA_on_pre)
    except ValueError:
        dA_on_pre_max = 0.

    try:
        dA_on_post_min = np.min(dA_on_post)
    except ValueError:
        dA_on_post_min = 0.

    try:
        dA_on_post_max = np.min(dA_on_post)
    except ValueError:
        dA_on_post_max = 0.

    dA_max = np.max((-1*dA_on_pre_min, dA_on_post_max))
    dA_min = np.min((-1*dA_on_pre_max, dA_on_post_min))
        
    return dA_max, dA_min

def synapse_weight_change_on_spike(ax, tr, crun='run_00000000', bins=100):
    '''
    '''
    spk_rgst = tr.crun.spk_register
    if len(spk_rgst)==0:
        print("no recorded spikes, not plotting weight changes")
        ax.set_title("No recorded weight changes on spike")
    else:
        dA_on_pre, dA_on_post = extract_delta_a_on_spike(spk_rgst)

        dA_max, dA_min = get_dA_prepost_minmax(dA_on_pre, dA_on_post)

        dA_post_bins = np.linspace(dA_min, dA_max, bins)
        dA_pre_bins = -1*np.flip(dA_post_bins,0)

        ax.hist(dA_on_post, bins=dA_post_bins)
        ax.hist(dA_on_pre, bins=dA_pre_bins)
        ax.set_title('Weight Changes through STDP')
        ax.set_xlabel(r'$\Delta a$')

def synapse_weight_change_on_spike_log(ax, tr, crun='run_00000000', bins=100):
    '''
    '''
    spk_rgst = tr.crun.spk_register
    if len(spk_rgst)==0:
        print("no recorded spikes, not plotting weight changes")
        ax.set_title("No recorded weight changes on spike")

    else:
        dA_on_pre, dA_on_post = extract_delta_a_on_spike(spk_rgst)

        dA_max, dA_min = get_dA_prepost_minmax(dA_on_pre, dA_on_post)
        
        dA_post_bins = np.linspace(dA_min, dA_max, bins)
        dA_pre_bins = -1*np.flip(dA_post_bins,0)

        ax.hist(dA_on_post, bins=dA_post_bins, log=True)
        ax.hist(dA_on_pre, bins=dA_pre_bins, log=True)

        ax.set_title('Weight Changes through STDP')
        ax.set_xlabel(r'$\Delta a$')

def synapse_weight_change_on_spike_loglog(ax, tr, crun='run_00000000',
                                          bins=100):
    '''
    '''
    spk_rgst = tr.crun.spk_register
    if len(spk_rgst)==0:
        print("no recorded spikes, not plotting weight changes")
        ax.set_title("No recorded weight changes on spike")

    else:
        dA_on_pre, dA_on_post = extract_delta_a_on_spike(spk_rgst)

        dA_max, dA_min = get_dA_prepost_minmax(dA_on_pre, dA_on_post)
        
        dA_post_bins = np.linspace(dA_min, dA_max, bins)
        dA_pre_bins = -1*np.flip(dA_post_bins,0)


        counts, edges = np.histogram(dA_on_post, bins=dA_post_bins)
        centers = (edges[:-1] + edges[1:])/2.
        ax.plot(centers, counts, '.', markersize=1.)

        counts, edges = np.histogram(-1.*dA_on_pre, bins=dA_post_bins)
        centers = (edges[:-1] + edges[1:])/2.
        ax.plot(centers, counts, '.', markersize=2.)

        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_title('Weight Changes through STDP')
        ax.set_xlabel(r'$\Delta a$')



def Apre_traces(ax, tr, crun='run_00000000', tmin=0, tmax=-1):
    df = tr.crun.SynEE_stat
    for i in range(len(df.a)):
        ax.plot(df.t[tmin:tmax]/ms,df.Apre[i,tmin:tmax])
    ax.set_title('Apre Traces')
    ax.set_xlabel('time [ms]')

def Apost_traces(ax, tr, crun='run_00000000', tmin=0, tmax=-1):
    df = tr.crun.SynEE_stat
    for i in range(len(df.a)):
        ax.plot(df.t[tmin:tmax]/ms,df.Apost[i,tmin:tmax])
    ax.set_title('Apost Traces')
    ax.set_xlabel('time [ms]')

def membrane_threshold_distribution_t(ax, tr, crun='run_00000000', tstep=-1):

    df = tr.crun.GExc_vts

    if df == []:
        pass

    else:
    
        for i,t in enumerate(df.t):
            print('''Assuming unit of t as second here, 
                     but need to confirm and/ or fix!''')
            #if t*second == tr.sim.T:
            if t == df.t[tstep]:
                ax.hist(df.Vt[:,i]/mV, label='Exc')
                
        ax.set_title('Mem. Thresholds at t=\SI{'+ \
                     str(df.t[tstep])+'}{s}')
        ax.set_xlabel('threshold [mV]')
        ax.set_ylabel('frequency')
        ax.legend(loc='upper right', frameon=False)

def membrane_threshold_traces(ax, tr, crun='run_00000000', tmax=-1):

    df = tr.crun.GExc_stat

    if df == []:
        pass

    try:
        for i in df.record:
            ax.plot(df.t[:tmax]/ms, df.Vt[i,:tmax]/mV)
        ax.ticklabel_format(useOffset=False)

        ax.set_title('Mem. Threshold Traces')
        ax.set_xlabel('time [ms]')
        ax.set_ylabel('voltage [mV]')

    except AttributeError:
        ax.axis('off')      
    
   
    
        
