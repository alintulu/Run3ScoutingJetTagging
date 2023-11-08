import numpy as np
import awkward as ak
import mplhep

def mean(h):
    total = 0
    for i, (center, value) in enumerate(zip(h.axes[0].centers, h.values())):
        total += center * value
    return total / sum(h.values())

def median(h, returnValue=False):
    centers = h.axes[0].centers
    
    total = 0
    median_index = (sum(h.values()) + 1) / 2
    for i, value in enumerate(h.values()):
        total += value
        if total > median_index:
            if returnValue:
                return centers[i]
            return i, value

def get_median(xvals, yvals, bin_edges, Neff):
    ''' Calculate median and median error (assuming Gaussian distribution).
    This is the binned median, not the real data median
    Extrapolation withing bins is performed.
    '''
    yvals_cumsum = np.cumsum(yvals)
    N = np.sum(yvals)

    # once adding weights, Neff appears to be ~1/4 - 1/3 of N when not using weights,
    # so changing limits to match the both cases
    if np.abs(np.sum(yvals)-Neff)/Neff<1e-5:
        N_min_limit=200
    else:
        N_min_limit=50

    if Neff>N_min_limit:
        med_bin = np.nonzero(yvals_cumsum>N/2)[0][0]
        median = bin_edges[med_bin] + (N/2 - yvals_cumsum[med_bin-1])/yvals[med_bin]*(bin_edges[med_bin+1]
                                                                                 - bin_edges[med_bin])
    else:
        median = 0

    hist_mean = np.sum(xvals*yvals)/sum(yvals) 
    hist_rms = np.sqrt(np.sum(yvals*((hist_mean-xvals)**2))/sum(yvals))
    medianstd = 1.253 * hist_rms/np.sqrt(Neff)
    
    return median, medianstd
        
def confidence(h, confLevel = 0.68):
    
    def get_width(i, edges):
        return edges[i + 1] - edges[i]

    values = h.values()
    centers = h.axes[0].centers
    edges = h.axes[0].edges
    ix, nsum = median(h)
    ixlow = ix
    ixhigh = ix
    nb = len(centers)
    ntot = sum(values)
    width = get_width(ix, edges)
    if ntot == 0:
        return 0
    
    while (nsum < confLevel * ntot):
        nlow = values[ixlow - 1] if ixlow > 0 else 0
        nhigh = values[ixhigh + 1] if ixhigh < nb else 0
        
        if (nsum + max(nlow, nhigh) < confLevel * ntot):
            if (nlow >= nhigh and ixlow > 0):
                nsum += nlow
                ixlow -= 1
                width += get_width(ixlow, edges)
            elif ixhigh < nb:
                nsum += nhigh
                ixhigh +=1
                width += get_width(ixhigh, edges)
            else: raise ValueError('BOOM')
        else:
            if (nlow > nhigh):
                width += get_width(ixlow-1, edges) * (confLevel * ntot - nsum) / nlow
            else:
                width += get_width(ixhigh+1, edges) * (confLevel * ntot - nsum) / nhigh
            nsum = ntot
            
    return width / (2*0.9945)

def err(meanhist):
    return np.sqrt(meanhist.variances() * meanhist.counts())

def err_err(meanhist):
    return np.sqrt(err(meanhist)**2 / (2 * meanhist.counts()))

def confidence_unc(sigma, counts):
    return np.sqrt(sigma**2 / (2 * counts))

def err_prop(val_1, val_2, err_1, err_2):
    
    err = []
    
    for i in range(len(val_1)):
        tmp = (val_1[i] / val_2[i]) * np.sqrt((err_1[i] / val_1[i])**2 + (err_2[i] / val_2[i])**2)
        err.append(tmp)
        
    return err

from iminuit import cost, Minuit
from numba_stats import norm , crystalball
import numpy as np
from scipy.stats import multivariate_normal as mvnorm
from scipy import interpolate

def crystalball_cdf(x, alpha, n, mean, sigma):
    return crystalball.cdf(x, alpha, n, mean, sigma)

def gauss_cdf(x, mean, sigma):
    return norm.cdf(x, mean, sigma)

def get_centers_edges(c):
    n = c._masked[..., 0] if c._bztrafo else c._masked
    ne = (c._masked[..., 1] if c._bztrafo else c._masked) ** 0.5
    xe = c.xe
    cx = 0.5 * (xe[1:] + xe[:-1])
    if c.mask is not None:
        cx = cx[c.mask]
        
    return cx, xe

def get_prediction(c, args):
    return c.prediction(args)

def get_xerr(centers):
    return [
            [np.abs(centers[i] - centers[i-1]) / 2 if i > 0 else np.abs(centers[i+1] - centers[i]) / 2 for i, _ in enumerate(centers)],
            [np.abs(centers[i+1] - centers[i]) / 2 if i < len(centers) - 1 else np.abs(centers[i] - centers[i-1]) / 2 for i, _ in enumerate(centers)]
        ]

def set_title(is_data, ax):
    if is_data:
        label= mplhep.cms.label(ax=ax, data=is_data, year="2022", com=13.6, label="Preliminary", lumi=34)
    else:
        label = mplhep.cms.label(ax=ax, data=is_data, year="2022", com=13.6, label="Preliminary")
    return label