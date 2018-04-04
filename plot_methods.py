import matplotlib, math
matplotlib.use('Agg')
import matplotlib.pyplot as pl

import numpy as np
from decimal import Decimal

from brian2.units import mV, ms, second

from post_processing import extract_lifetimes, extract_active_synapse_count, \
                            extract_delta_a_on_spike


def raster_plot(ax, tr, crun='run_00000000', tmin=0.*second, tmax=-1.*second):
    if tmax == -1.*second:
        tmax=tr.T
        
    df_e, df_i = tr.crun.GExc_spks, tr.crun.GInh_spks
    
    if (df_e == []) or (df_i == []):
        pass
    
    else:
        
        try:
            indx = np.logical_and(df_e.t/ms>tmin/ms, df_e.t/ms<tmax/ms)
            ax.plot(df_e.t[indx]/second, df_e.i[indx], marker='.',
                    color='blue', markersize=.5, linestyle='None')
        except AttributeError:
            print(crun, " no exc. spikes")
            
        try:
            indx = np.logical_and(df_i.t/ms>tmin/ms, df_i.t/ms<tmax/ms)
            ax.plot(df_i.t[indx]/second, df_i.i[indx]+tr.N_e, marker='.',
                    color='red', markersize=.5, linestyle='None')
        except AttributeError:
            print(crun, " no inh. spikes")

        ax.set_ylim(0, tr.N_e + tr.N_i)
        ax.set_title('T='+str(tr.T/second)+' s, ' +\
                     '$\sigma$=$\sqrt{\mathrm{' +\
                     '%.2E' % Decimal((tr.sigma/mV)**2)+'}}$ mV')
        ax.set_xlabel('time [s]')

        
def firing_rate_distribution_exc(ax, tr, crun='run_00000000', steps=25):

    df_e = tr.crun.GExc_spks

    if df_e == []:
        pass
    
    else:
    
        try:
            spk_counts = np.bincount(df_e.i)
            delta_spk = math.ceil(0.05*(np.max(spk_counts)-np.min(spk_counts)))
            bins = np.linspace(np.min(spk_counts)-delta_spk,
                             np.max(spk_counts)+delta_spk,
                             steps)/tr.T/second
            ax.hist(spk_counts/tr.T/second, bins=bins, color='blue')
        except AttributeError:
            print(crun, " no exc. spikes")
            
        text, textcol = 'inactive', 'red'

        if bool(tr.it_active):
            text, textcol = 'active','green'

        ax.text(0.02, 0.95,
                'it ' + text, color=textcol,
                horizontalalignment='left',
                verticalalignment='top',
                bbox={'boxstyle': 'square, pad=0.3', 'facecolor':'white',
                      'alpha':1, 'edgecolor':'none'},
                transform = ax.transAxes)        

        ax.set_title('Firing Rate Distribution (Exc.)')
        ax.set_xlabel('firing rate [Hz]')
        ax.set_ylabel('counts')

        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')


def firing_rate_distribution_inh(ax, tr, crun='run_00000000', steps=25):
    df_i = tr.crun.GInh_spks
    try:
        spk_counts = np.bincount(df_i.i)
        delta_spk = math.ceil(0.05*(np.max(spk_counts)-np.min(spk_counts)))
        bins = np.linspace(np.min(spk_counts)-delta_spk,
                         np.max(spk_counts)+delta_spk,
                           steps)/tr.T/second
        ax.hist(spk_counts/(tr.T/second), bins=bins, color='red')
    except AttributeError:
        print(crun, " no inh. spikes")
    ax.set_title('Firing Rate Distribution (Inh.)')
    ax.set_xlabel('firing rate [Hz]')
    ax.set_ylabel('counts')

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')



def voltage_traces(ax, tr, crun='run_00000000', tmin=0.*second, tmax=-1.*second):
    
    if tmax == -1.*second:
        tmax=tr.T

    df = tr.crun.GExc_stat

    if df == []:
        pass
    
    else:
    
        print("Assuming t in second. It's not clear why df.t doesn't",
              "\ncarry units and this should be fixed eventually!")

        indx = np.logical_and(df.t>tmin/second, df.t<tmax/second)

        for i in df.record:
            ax.plot(df.t[indx]/second, df.V[i,:][indx]/mV)

        ax.set_ylim(tr.Vr_e/mV-5, tr.Vt_e/mV+5)
        ax.set_title('Membrane Voltage Traces')
        ax.set_xlabel('time [s]')
        ax.set_ylabel('voltage [mV]')

        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')
        

        
def isi_distribution(ax, tr, crun='run_00000000', bins=25):
    df_e, df_i = tr.crun.GExc_spks, tr.crun.GInh_spks
    exc_isi = []
    try:
        for k in range(tr.N_e):
            exc_isi += list(np.diff(df_e.t[df_e.i==k])/ms)
        ax.hist(exc_isi, bins=bins)
    except AttributeError:
        print(crun, " no exc. spikes")
    ax.set_title('ISI Distribution (Exc.)')
    ax.set_xlabel('interspike interval [ms]')
    ax.set_ylabel('counts')

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')


    # try:
    #     spk_counts = np.bincount(df_e.i)
    #     delta_spk = math.ceil(0.1*(np.max(spk_counts)-np.min(spk_counts)))
    #     bins = np.arange(np.min(spk_counts)-delta_spk,
    #                      np.max(spk_counts)+step+delta_spk,
    #                      step)/tr.T/second
    #     ax.hist(spk_counts/tr.T/second, bins=bins, color='blue')
    # except AttributeError:
    #     print(crun, " no exc. spikes")
    # ax.set_title('Firing Rate Distribution (Exc.)')
    # ax.set_xlabel('firing rate [Hz]')
    # ax.set_ylabel('counts')
    #return df_e

def conductance_traces(ax, tr, crun='run_00000000', tmin=0, tmax=-1):
    df = tr.crun.GExc_stat 
    for i in df.record:
        if i==2:
            ax.plot(df.t[tmin:tmax]/ms, df.ge[i,tmin:tmax])
            ax.plot(df.t[tmin:tmax]/ms, df.gi[i,tmin:tmax])
    #ax.set_ylim(tr.Vr_e/mV-5, tr.Vt_e/mV+5)
    ax.set_title('Conductance Traces')
    ax.set_xlabel('time [ms]')
    ax.set_ylabel('conductance')

def conductance_mult_trace(ax, tr, crun='run_00000000', tmin=0, tmax=-1):

    df = tr.crun.GExc_stat

    if df == []:
        pass
    
    else:
    
        ge_max, gi_min = [],[]

        for i in df.record:
            gi_min.append(1.1*np.min(df.gi[i,tmin:tmax]))
            ax.plot(df.t[tmin:tmax]/ms,
                    np.sum(ge_max)-np.sum(gi_min)+df.ge[i,tmin:tmax],
                    color='blue')
            ax.plot(df.t[tmin:tmax]/ms,
                    np.sum(ge_max)-np.sum(gi_min)+df.gi[i,tmin:tmax],
                    color='red')

            ge_max.append(1.1*np.max(df.ge[i,tmin:tmax]))

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_title('Conductance Traces')
        ax.set_xlabel('time [ms]')
        ax.set_ylabel('conductance')
    

def synapse_weight_distribution(ax, tr, crun='run_00000000', bins=50):
    df = tr.crun.SynEE_a
    bins = np.linspace(0, tr.amax, num=bins)
    for i,t in enumerate(df.t):
        ax.hist(np.nan_to_num(df.a[:,i].flatten()), bins=bins)
    ax.set_title('Synaptic Weight Distribution')

def synapse_weight_distribution_t(ax, tr, crun='run_00000000', bins=25, tstep=0, low_bound=-1):
    df = tr.crun.SynEE_a

    data = df['a'][:,tstep]
    data = data[df['syn_active'][:,tstep]==1]

    if low_bound!=-1:
        data = data[data>low_bound]
    
    print('max: ', np.max(data), '\t min: ', np.min(data))
    # data = data[data>tr.prn_thrshld]

    # try:
    #     amax = np.max(data)
    # except ValueError:
    #     amax = 4*tr.a_insert
    # ax.hist(data, bins=np.linspace(tr.prn_thrshld, 4*tr.a_insert, bins))

    ax.hist(data, bins=bins)
    
    ax.text(0.50, 0.95,
            r'$a_{\text{max}}$='+'%.2E' % Decimal(tr.amax) + '\n' + \
            'prn=' + '%.2E' % Decimal(tr.prn_thrshld) +'\n' +\
            'ins=' + '%.2E' % Decimal(tr.a_insert),
            horizontalalignment='left',
            verticalalignment='top',
            bbox={'boxstyle': 'square, pad=0.3', 'facecolor':'white',
                  'alpha':1, 'edgecolor':'none'},
            transform = ax.transAxes)
            
    # ax.axvline(tr.a_insert, color='red')
    # ax.axvline(tr.prn_thrshld, color='red')

    text, textcol = 'inactive', 'red'
    if bool(tr.scl_active):
        text, textcol = 'active', 'green'
    if tstep==0:
        ax.text(0.02, 0.95,
                'scl ' + text, color=textcol, 
                horizontalalignment='left',
                verticalalignment='top',
                bbox={'boxstyle': 'square, pad=0.3', 'facecolor':'white',
                      'alpha':1, 'edgecolor':'none'},
                transform = ax.transAxes)        

    ax.set_title('SynEE Weights at t=\SI{'+ \
                 str(df.t[tstep])+'}{s}')

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')


def synapse_weight_distribution_log_t(ax, tr, crun='run_00000000',
                                      bins=30, a_min=-1, a_max=-1,
                                      tstep=-1, low_bound = -1):

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

    ax.hist(a_cut, 10**np.linspace(np.log10(a_min),
                                   np.log10(a_max),bins),
            density=True)

    ax.set_xscale('log')

    # ax.axvline(tr.a_insert, color='red')
    # ax.axvline(tr.prn_thrshld, color='red')

    ax.set_title('SynEE Weights at t=\SI{'+ \
                 str(df.t[tstep])+'}{s}')
    
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')

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
            ax.plot(df.t[indx],df.a[i,:][indx])
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

    else:
    
        for i in df.record:
            ax.plot(df.t[:tmax]/ms, df.Vt[i,:tmax]/mV)
        ax.ticklabel_format(useOffset=False)
        ax.set_title('Mem. Threshold Traces')
        ax.set_xlabel('time [ms]')
        ax.set_ylabel('voltage [mV]')

        
def print_membrane_params(ax, tr, crun='run_00000000'):
    ax.axis('off')
    ax.text(0.0, 1.0,
            '$V_t^e=$\SI{'+str(tr.Vt_e/mV)+'}{mV}' + \
            '\n$V_t^i=$\SI{'+str(tr.Vt_i/mV)+'}{mV}\n' +\
            r'$\eta_{\text{IP}}=$\SI{' +\
            '%.2E' % Decimal(tr.eta_ip/(mV*ms))+'}{mV ms}'+\
            '\n'+r'$\Delta_t^{\mathrm{it}}$='+str(tr.it_dt/ms)+' ms',
            horizontalalignment='left',
            verticalalignment='top',
            bbox={'boxstyle': 'square, pad=0.3', 'facecolor':'white',
                  'alpha':1, 'edgecolor':'none'},
            transform = ax.transAxes)

def print_netw_params(ax, tr, crun='run_00000000'):
    ax.axis('off')
    ax.text(0.0, 1.0,
            '$N_e$=' +str(int(tr.N_e)) + ', $N_i$=%d' %(int(tr.N_i)) +\
            ', T=%.2f s' %(tr.T/second) +
            '\n$\sigma$=$\sqrt{\mathrm{' +\
            '%.2E' % Decimal((tr.sigma/mV)**2)+'}}$ mV' +\
            '\n $'+r'\tau'+'$=\SI{%.2f}{ms}' %(tr.tau/ms) +\
            '\n$E_l$= \SI{%.0f}{mV}' % (tr.El/mV) +\
            '$, E_e$= \SI{%.0f}{mV}' % (tr.Ee/mV) +\
            '$, E_i$= \SI{%.0f}{mV}' % (tr.Ei/mV) +\
            '\n $'+r'\tau_e'+'$=\SI{%.2f}{ms}' % (tr.tau_e/ms) +\
            ', $'+r'\tau_i'+'$=\SI{%.2f}{ms}' % (tr.tau_i/ms) + \
            '\n $V^r_{e}$=\SI{%.0f}{mV}' % (tr.Vr_e/mV) +\
            ', $V^r_{i}$=\SI{%.0f}{mV}' % (tr.Vr_i/mV) +\
            '\n $V^T_{e}$=\SI{%.0f}{mV}' % (tr.Vt_e/mV) +\
            ', $V^T_{i}$=\SI{%.0f}{mV}' % (tr.Vt_i/mV) +\
            '\n $p_{ee}$ = %.2f' % (tr.p_ee) +\
            ', $p_{ie}$ = %.2f' % (tr.p_ie) +\
            ', $p_{ei}$ = %.2f' % (tr.p_ei) +\
            ', $p_{ii}$ = %.2f' % (tr.p_ii) +\
            '\n $a_{ie}$ = %.2f $a_{\mathrm{scale}}$' %(tr.a_ie/tr.ascale) +\
            ', $a_{ei}$ = %.2f $a_{\mathrm{scale}}$' % (tr.a_ei/tr.ascale) +\
            ', $a_{ii}$ = %.2f $a_{\mathrm{scale}}$' % (tr.a_ii/tr.ascale),
            horizontalalignment='left',
            verticalalignment='top',
            linespacing = 1.95,
            bbox={'boxstyle': 'square, pad=0.3', 'facecolor':'white',
                  'alpha':1, 'edgecolor':'none'},
            transform = ax.transAxes)

def print_stdp_params(ax, tr, crun='run_00000000'):
    ax.axis('off')
    ax.text(0.0, 1.0,
            r'$\tau'+'_{\mathrm{pre}}$= \SI{%.2f}{ms}' %(tr.taupre/ms) +\
            r', $\tau'+'_{\mathrm{post}}$= \SI{%.2f}{ms}' %(tr.taupost/ms) +\
            '\n $A_{\mathrm{plus}}$= $%f\, a_{\mathrm{scale}}$' %(tr.Aplus/tr.ascale) + \
            '\n $A_{\mathrm{minus}}$= $%f\, a_{\mathrm{scale}}$' %(tr.Aminus/tr.ascale) + \
            '\n $a_{\mathrm{scale}}$ = %.4E' %(Decimal(tr.ascale)) +\
            '\n------------------------------------------' +\
            '\n $a_{\mathrm{TotalMax}}$ = $%.4f \, a_{\mathrm{scale}}$' %(tr.ATotalMax/tr.ascale) +\
            '\n $\Delta_{\mathrm{scaling}}$ = \SI{%.2f}{ms}' %(tr.dt_synEE_scaling/ms) +\
            '\n $\eta_{\mathrm{scaling}}$ = %.4f' %(tr.eta_scaling) +\
            '\n------------------------------------------' +\
            '\n $a_{\mathrm{insert}}$ = $%.2f\, a_\mathrm{scale}$' %(tr.a_insert/tr.ascale) +\
            ', $a_{\mathrm{thrshld}}$ = $%.2f\, a_\mathrm{scale}$' %(tr.prn_thrshld/tr.ascale) +\
            '\n $\Delta_{\mathrm{strct}}$ = \SI{%.2f}{ms}' %(tr.strct_dt/ms) +\
            ', $p_{\mathrm{insert}}$ = %.4f' %(tr.insert_P),
            horizontalalignment='left',
            verticalalignment='top',
            linespacing = 1.75,
            bbox={'boxstyle': 'square, pad=0.3', 'facecolor':'white',
                  'alpha':1, 'edgecolor':'none'},
            transform = ax.transAxes)


def ax_off(ax):
    ax.axis('off')
