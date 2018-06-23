
import numpy as np

from elephant.conversion import BinnedSpikeTrain
from elephant.spike_train_correlation import corrcoef, cch
from brian2.units import second, ms
import quantities as pq
import neo

from mpl_toolkits.axes_grid1 import make_axes_locatable



def correlation_matrix(fig, ax, tr, crun='run_00000000', N=50, nbin=None):

    if nbin==None:
        nbin=int(tr.T/(50*ms))

    df = tr.crun.GExc_spks
    xt, xi = df.t, df.i

    sts = [neo.SpikeTrain(xt[xi==i]/second*pq.s,
                          t_stop=tr.T/second*pq.s) for i in range(N)]

    x = corrcoef(BinnedSpikeTrain(sts, num_bins = nbin))


    x[np.diag_indices(N)]=0

    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.05)

    im = ax.imshow(x)
    fig.colorbar(im, cax=cax, orientation='vertical')

# fig.savefig('this.png')


# with open("sEE_adj.p", "rb") as pfile:
#     sEE_src, sEE_tar = pickle.load(pfile)

# adj = np.zeros((400,400))
# for src,tar in zip(sEE_src, sEE_tar):
#     adj[src,tar]=1 

# adj_cut = adj[:N,:N]

# print(adj_cut[x>0.3])
# #for srcs, tars in zip(*np.where(x>0.3)):


    
# fig = pl.figure()
# ax = fig.add_subplot(111)

# x[np.diag_indices(N)]=0
# cax = ax.imshow(x)
# cbar = fig.colorbar(cax)

# fig.savefig('this.png')
