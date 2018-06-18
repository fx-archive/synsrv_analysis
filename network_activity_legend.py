
from decimal import Decimal
from brian2.units import mV, ms, second

def netw_params_table(ax, tr, crun='run_00000000'):

    ax.axis('off')

    text =  '$N_e$=' +str(int(tr.N_e)) + ', $N_i$=%d' %(int(tr.N_i)) +\
            ', T=%.2f s' %(tr.T/second) +\
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
            ', $a_{ii}$ = %.2f $a_{\mathrm{scale}}$' % (tr.a_ii/tr.ascale)

    ax.text(0.0, 1.0, text,
            horizontalalignment='left',
            verticalalignment='top',
            linespacing = 1.95,
            bbox={'boxstyle': 'square, pad=0.3', 'facecolor':'white',
                  'alpha':1, 'edgecolor':'none'},
            transform = ax.transAxes)
