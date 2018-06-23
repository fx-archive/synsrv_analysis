
from decimal import Decimal
from brian2.units import mV, ms, second

def netw_params_display(ax, tr, crun='run_00000000'):

    ax.axis('off')

    text = r'\textbf{Network configuration}'

    ax.text(0., 1.0, text,
            horizontalalignment='left',
            verticalalignment='top',
            linespacing = 1.95,
            fontsize=13,
            bbox={'boxstyle': 'square, pad=0.3', 'facecolor':'white',
                  'alpha':1, 'edgecolor':'none'},
            transform = ax.transAxes)


    text = '$N_e$=' +str(int(tr.N_e)) + ', $\, N_i$=%d' %(int(tr.N_i)) +\
           '\n $p_{ee}$ = %.2f' % (tr.p_ee) +\
           ', $\, p_{ie}$ = %.2f' % (tr.p_ie) +\
           '\n $p_{ei}$ = %.2f' % (tr.p_ei) +\
           ', $\, p_{ii}$ = %.2f' % (tr.p_ii)

    ax.text(0., 0.85, text,
            horizontalalignment='left',
            verticalalignment='top',
            linespacing = 1.95,
            fontsize=12,
            bbox={'boxstyle': 'square, pad=0.3', 'facecolor':'white',
                  'alpha':1, 'edgecolor':'none'},
            transform = ax.transAxes)

    
    text = r'\textbf{Simulation}'

    ax.text(0., 0.365, text,
            horizontalalignment='left',
            verticalalignment='top',
            linespacing = 1.95,
            fontsize=13,
            bbox={'boxstyle': 'square, pad=0.3', 'facecolor':'white',
                  'alpha':1, 'edgecolor':'none'},
            transform = ax.transAxes)

    
    text = 'T=%.2f s' %(tr.T/second)


    ax.text(0., 0.215, text,
            horizontalalignment='left',
            verticalalignment='top',
            linespacing = 1.95,
            fontsize=12,
            bbox={'boxstyle': 'square, pad=0.3', 'facecolor':'white',
                  'alpha':1, 'edgecolor':'none'},
            transform = ax.transAxes)

    


def neuron_params_display(ax, tr, crun='run_00000000'):

    ax.axis('off')

    text = r'\textbf{Neuron parameters}'

    ax.text(0., 1.0, text,
            horizontalalignment='left',
            verticalalignment='top',
            linespacing = 1.95,
            fontsize=13,
            bbox={'boxstyle': 'square, pad=0.3', 'facecolor':'white',
                  'alpha':1, 'edgecolor':'none'},
            transform = ax.transAxes)


    
    text = '$'+r'\tau'+'$=\SI{%.0f}{ms}' %(tr.tau/ms) +\
           ', $\,E_l$= \SI{%.0f}{mV}' % (tr.El/mV) +\
           '\n $\, E_e$= \SI{%.0f}{mV}' % (tr.Ee/mV) +\
           '$,\, E_i$= \SI{%.0f}{mV}' % (tr.Ei/mV) +\
           '\n $'+r'\,\,\,\,\,\, \tau_e'+'$=\SI{%.0f}{ms}' % (tr.tau_e/ms) +\
           ', $'+r'\tau_i'+'$=\SI{%.0f}{ms}' % (tr.tau_i/ms)

    ax.text(0., 0.85, text,
            horizontalalignment='left',
            verticalalignment='top',
            linespacing = 1.95,
            fontsize=12,
            bbox={'boxstyle': 'square, pad=0.3', 'facecolor':'white',
                  'alpha':1, 'edgecolor':'none'},
            transform = ax.transAxes)
        

       
  
    text = r'\textbf{excitatory, inhibitory}'

    ax.text(0., 0.365, text,
            horizontalalignment='left',
            verticalalignment='top',
            linespacing = 1.95,
            fontsize=13,
            bbox={'boxstyle': 'square, pad=0.3', 'facecolor':'white',
                  'alpha':1, 'edgecolor':'none'},
            transform = ax.transAxes)

    
    text = '$\sigma_e$=$\sqrt{\mathrm{' +\
           '%.2E' % Decimal((tr.sigma_e/mV)**2)+'}}$ mV' +\
           ', $\, \sigma_i$=$\sqrt{\mathrm{' +\
           '%.2E' % Decimal((tr.sigma_i/mV)**2)+'}}$ mV' +\
           '\n $V^T_{e}$=\SI{%.2f}{mV}' % (tr.Vt_e/mV) +\
           ', $V^T_{i}$=\SI{%.2f}{mV}' % (tr.Vt_i/mV) +\
           '\n $V^r_{e}$=\SI{%.0f}{mV}' % (tr.Vr_e/mV) +\
           ', $V^r_{i}$=\SI{%.0f}{mV}' % (tr.Vr_i/mV) 


    ax.text(0., 0.2, text,
            horizontalalignment='left',
            verticalalignment='top',
            linespacing = 1.95,
            fontsize=12,
            bbox={'boxstyle': 'square, pad=0.3', 'facecolor':'white',
                  'alpha':1, 'edgecolor':'none'},
            transform = ax.transAxes)



    


def synapse_params_display(ax, tr, crun='run_00000000'):

    ax.axis('off')

    text = r'\textbf{Synapse parameters}'

    ax.text(0., 0.9, text,
            horizontalalignment='left',
            verticalalignment='top',
            linespacing = 1.95,
            fontsize=13,
            bbox={'boxstyle': 'square, pad=0.3', 'facecolor':'white',
                  'alpha':1, 'edgecolor':'none'},
            transform = ax.transAxes)


    
    text = '$a_{ee}=$ %.2E ' %(Decimal(tr.a_ee/tr.ascale)) +\
           '$a_{\mathrm{scale}}$' +\
           '\n$a_{ie}=$ %.2E ' %(Decimal(tr.a_ie/tr.ascale)) +\
           '$a_{\mathrm{scale}}$' +\
           '\n $a_{ei}$ = %.2E ' % (Decimal(tr.a_ei/tr.ascale)) +\
           '$a_{\mathrm{scale}}$' +\
           '\n $a_{ii}$ = %.2E '  % Decimal((tr.a_ii/tr.ascale))  +\
           '$a_{\mathrm{scale}}$' +\
           '\n$a_{\mathrm{scale}} =$ %.3E' % Decimal(tr.ascale)

    ax.text(0., 0.75, text,
            horizontalalignment='left',
            verticalalignment='top',
            linespacing = 1.95,
            fontsize=12,
            bbox={'boxstyle': 'square, pad=0.3', 'facecolor':'white',
                  'alpha':1, 'edgecolor':'none'},
            transform = ax.transAxes)
        

       
  

def stdp_params_display(ax, tr, crun='run_00000000'):

    ax.axis('off')

    text = r'\textbf{STDP parameters}'

    ax.text(0., 1.0, text,
            horizontalalignment='left',
            verticalalignment='top',
            linespacing = 1.95,
            fontsize=13,
            bbox={'boxstyle': 'square, pad=0.3', 'facecolor':'white',
                  'alpha':1, 'edgecolor':'none'},
            transform = ax.transAxes)
    
    
    text = r'$\tau'+'_{\mathrm{pre}} =\,$ \SI{%.2f}{ms}' %(tr.taupre/ms) +\
            r',$\,\tau'+'_{\mathrm{post}}=\,$ \SI{%.2f}{ms}' %(tr.taupost/ms) +\
            '\n $A_{\mathrm{plus}}\,\,\,$= $\,\,\,\,%f$' %(tr.Aplus) + \
            '\n $A_{\mathrm{minus}}$= $%f$' %(tr.Aminus) +\
            '\n $a_{\mathrm{max}} = '+'%f' %(tr.amax) +'$'

    ax.text(0., 0.85, text,
            horizontalalignment='left',
            verticalalignment='top',
            linespacing = 1.95,
            fontsize=12,
            bbox={'boxstyle': 'square, pad=0.3', 'facecolor':'white',
                  'alpha':1, 'edgecolor':'none'},
            transform = ax.transAxes)
        



def sn_params_display(ax, tr, crun='run_00000000'):

    ax.axis('off')

    text = r'\textbf{Synaptic scaling parameters}'

    ax.text(0., 1.0, text,
            horizontalalignment='left',
            verticalalignment='top',
            linespacing = 1.95,
            fontsize=13,
            bbox={'boxstyle': 'square, pad=0.3', 'facecolor':'white',
                  'alpha':1, 'edgecolor':'none'},
            transform = ax.transAxes)


    text = r'$a_{\mathrm{TotalMax}} = %f $' %(tr.ATotalMax) +\
           '\n $\Delta_{\mathrm{scaling}} =\,$' +\
           '\SI{%.2f}{ms}' %(tr.dt_synEE_scaling/ms) + \
           '\n $\eta_{\mathrm{scaling}} = %f$' %(tr.eta_scaling) 


    ax.text(0., 0.85, text,
            horizontalalignment='left',
            verticalalignment='top',
            linespacing = 1.95,
            fontsize=12,
            bbox={'boxstyle': 'square, pad=0.3', 'facecolor':'white',
                  'alpha':1, 'edgecolor':'none'},
            transform = ax.transAxes)
        

    

def sn_params_display(ax, tr, crun='run_00000000'):

    ax.axis('off')

    text = r'\textbf{Synaptic scaling parameters}'

    ax.text(0., 1.0, text,
            horizontalalignment='left',
            verticalalignment='top',
            linespacing = 1.95,
            fontsize=13,
            bbox={'boxstyle': 'square, pad=0.3', 'facecolor':'white',
                  'alpha':1, 'edgecolor':'none'},
            transform = ax.transAxes)


    text = r'$a_{\mathrm{TotalMax}} = %f $' %(tr.ATotalMax) +\
           '\n $\Delta_{\mathrm{scaling}} =\,$' +\
           '\SI{%.2f}{ms}' %(tr.dt_synEE_scaling/ms) + \
           '\n $\eta_{\mathrm{scaling}} = %f$' %(tr.eta_scaling) 


    ax.text(0., 0.85, text,
            horizontalalignment='left',
            verticalalignment='top',
            linespacing = 1.95,
            fontsize=12,
            bbox={'boxstyle': 'square, pad=0.3', 'facecolor':'white',
                  'alpha':1, 'edgecolor':'none'},
            transform = ax.transAxes)
        

    
def sn_params_display(ax, tr, crun='run_00000000'):

    ax.axis('off')

    text = r'\textbf{Synaptic scaling parameters}'

    ax.text(0., 1.0, text,
            horizontalalignment='left',
            verticalalignment='top',
            linespacing = 1.95,
            fontsize=13,
            bbox={'boxstyle': 'square, pad=0.3', 'facecolor':'white',
                  'alpha':1, 'edgecolor':'none'},
            transform = ax.transAxes)


    text = r'$a_{\mathrm{TotalMax}} = %f $' %(tr.ATotalMax) +\
           '\n $\Delta_{\mathrm{scaling}} =\,$' +\
           '\SI{%.2f}{ms}' %(tr.dt_synEE_scaling/ms) + \
           '\n $\eta_{\mathrm{scaling}} = %f$' %(tr.eta_scaling) 


    ax.text(0., 0.85, text,
            horizontalalignment='left',
            verticalalignment='top',
            linespacing = 1.95,
            fontsize=12,
            bbox={'boxstyle': 'square, pad=0.3', 'facecolor':'white',
                  'alpha':1, 'edgecolor':'none'},
            transform = ax.transAxes)
        

    
def sn_params_display(ax, tr, crun='run_00000000'):

    ax.axis('off')

    text = r'\textbf{Synaptic scaling parameters}'

    ax.text(0., 1.0, text,
            horizontalalignment='left',
            verticalalignment='top',
            linespacing = 1.95,
            fontsize=13,
            bbox={'boxstyle': 'square, pad=0.3', 'facecolor':'white',
                  'alpha':1, 'edgecolor':'none'},
            transform = ax.transAxes)


    text = r'$a_{\mathrm{TotalMax}} = %f $' %(tr.ATotalMax) +\
           '\n $\Delta_{\mathrm{scaling}} =\,$' +\
           '\SI{%.2f}{ms}' %(tr.dt_synEE_scaling/ms) + \
           '\n $\eta_{\mathrm{scaling}} = %f$' %(tr.eta_scaling) 


    ax.text(0., 0.85, text,
            horizontalalignment='left',
            verticalalignment='top',
            linespacing = 1.95,
            fontsize=12,
            bbox={'boxstyle': 'square, pad=0.3', 'facecolor':'white',
                  'alpha':1, 'edgecolor':'none'},
            transform = ax.transAxes)
        

    
def sn_params_display(ax, tr, crun='run_00000000'):

    ax.axis('off')

    text = r'\textbf{Synaptic scaling parameters}'

    ax.text(0., 1.0, text,
            horizontalalignment='left',
            verticalalignment='top',
            linespacing = 1.95,
            fontsize=13,
            bbox={'boxstyle': 'square, pad=0.3', 'facecolor':'white',
                  'alpha':1, 'edgecolor':'none'},
            transform = ax.transAxes)


    text = r'$a_{\mathrm{TotalMax}} = %f $' %(tr.ATotalMax) +\
           '\n $\Delta_{\mathrm{scaling}} =\,$' +\
           '\SI{%.2f}{ms}' %(tr.dt_synEE_scaling/ms) + \
           '\n $\eta_{\mathrm{scaling}} = %f$' %(tr.eta_scaling) 


    ax.text(0., 0.85, text,
            horizontalalignment='left',
            verticalalignment='top',
            linespacing = 1.95,
            fontsize=12,
            bbox={'boxstyle': 'square, pad=0.3', 'facecolor':'white',
                  'alpha':1, 'edgecolor':'none'},
            transform = ax.transAxes)
        

    
def sn_params_display(ax, tr, crun='run_00000000'):

    ax.axis('off')

    text = r'\textbf{Synaptic scaling parameters}'

    ax.text(0., 1.0, text,
            horizontalalignment='left',
            verticalalignment='top',
            linespacing = 1.95,
            fontsize=13,
            bbox={'boxstyle': 'square, pad=0.3', 'facecolor':'white',
                  'alpha':1, 'edgecolor':'none'},
            transform = ax.transAxes)


    text = r'$a_{\mathrm{TotalMax}} = %f $' %(tr.ATotalMax) +\
           '\n $\Delta_{\mathrm{scaling}} =\,$' +\
           '\SI{%.2f}{ms}' %(tr.dt_synEE_scaling/ms) + \
           '\n $\eta_{\mathrm{scaling}} = %f$' %(tr.eta_scaling) 


    ax.text(0., 0.85, text,
            horizontalalignment='left',
            verticalalignment='top',
            linespacing = 1.95,
            fontsize=12,
            bbox={'boxstyle': 'square, pad=0.3', 'facecolor':'white',
                  'alpha':1, 'edgecolor':'none'},
            transform = ax.transAxes)
        

    

def strct_params_display(ax, tr, crun='run_00000000'):

    ax.axis('off')

    text = r'\textbf{Structural plasticity parameters}'

    ax.text(0., 1.0, text,
            horizontalalignment='left',
            verticalalignment='top',
            linespacing = 1.95,
            fontsize=13,
            bbox={'boxstyle': 'square, pad=0.3', 'facecolor':'white',
                  'alpha':1, 'edgecolor':'none'},
            transform = ax.transAxes)


    text = '$p_{\mathrm{insert}} = %f$' %(tr.insert_P) +\
           '\n$a_{\mathrm{insert}} = %f$' %(tr.a_insert) +\
           '\n $a_{\mathrm{thrshld}} = %f$' %(tr.prn_thrshld) +\
           '\n $\Delta_{\mathrm{strct}} =\,$ \SI{%.2f}{ms}' %(tr.strct_dt/ms)


    ax.text(0., 0.85, text,
            horizontalalignment='left',
            verticalalignment='top',
            linespacing = 1.95,
            fontsize=12,
            bbox={'boxstyle': 'square, pad=0.3', 'facecolor':'white',
                  'alpha':1, 'edgecolor':'none'},
            transform = ax.transAxes)
        

    
