import numpy as np
import awkward as ak
import random
import os
import uproot

def run_deltar_matching(obj1, obj2, radius=0.4): # NxM , NxG arrays
    _, obj2 = ak.unzip(ak.cartesian([obj1, obj2], nested=True)) # Obj2 is now NxMxG
    obj2['dR'] = obj1.delta_r(obj2)  # Calculating delta R
    t_index = ak.argmin(obj2.dR, axis=-2) # Finding the smallest dR (NxG array)
    s_index = ak.local_index(obj1.eta, axis=-1) #  NxM array
    _, t_index = ak.unzip(ak.cartesian([s_index, t_index], nested=True)) 
    obj2 = obj2[s_index == t_index] # Pairwise comparison to keep smallest delta R

    # Cutting on delta R
    obj2 = obj2[obj2.dR < radius] # Additional cut on delta R, now a NxMxG' array 
    return obj2

def require_back2back(obj1, obj2, phi=2.7):
    return ak.to_numpy((abs(obj1.delta_phi(obj2)) > phi))

def require_3rd_jet(jets, pt_type="pt"):
    jet = jets[:, 2]
    pt_ave = (jets[:, 0][pt_type] + jets[:, 1][pt_type]) / 2

    return ak.to_numpy(((jet[pt_type] / pt_ave) < 0.05))

def require_eta(jets):
    return (
        (abs(jets[:, 0].eta) < 1.3) 
        | (abs(jets[:, 1].eta) < 1.3)
    )

def require_n(jet1, jet2=None, two=True):
    if two:
        jet_cut = np.array((ak.num(jet1) == 2))
        if jet2 is not None:
            jet_cut *= np.array((ak.num(jet2) == 2))
            return jet1[jet_cut], jet2[jet_cut]
    else:
        jet_cut = np.array((ak.num(jet1) > 2))
        if jet2 is not None:
            jet_cut *= np.array((ak.num(jet2) > 2))
            return jet1[jet_cut], jet2[jet_cut]

    return jet1[jet_cut]

def criteria_one(jet1, jet2=None, phi=2.7):
    b2b = require_back2back(jet1[:, 0], jet1[:, 1], phi)

    if jet2 is not None:
        b2b2 = require_back2back(jet2[:, 0], jet2[:, 1], phi)
        b2b *= b2b2
        return jet1[b2b], jet2[b2b]

    return jet1[b2b]

def criteria_n(jet1, jet2=None, phi=2.7):
    third_jet = require_3rd_jet(jet1, pt_type)
    b2b = require_back2back(jet1[:, 0], jet1[:, 1], phi)

    if jet2 is not None:
        third_jet2 = require_3rd_jet(jet2, pt_type)
        third_jet *= third_jet2
        
        b2b2 = require_back2back(jet2[:, 0], jet2[:, 1], phi)
        b2b *= b2b2
        return jet1[(third_jet & b2b)], jet2[(third_jet & b2b)]

    return jet1[(third_jet & b2b)]

def cut_ratio(ratio, h4=False):
    if h4:
        return ratio, np.ones(len(ratio), dtype=bool)
    cut = (np.abs(ratio - 1) < 2) 

    return ratio[cut], cut

def select_eta(eta, jet1, jet2=None):
    if eta not in {"barrel", "endcap", 0, 1}:
        raise ValueError("The requested eta region has to be either 'barrel', 'endcap', '0' or '1'")
        
    if eta == "barrel" or eta == 0:
        eta_cut = ak.to_numpy((
            (np.abs(jet1[:,0].eta) < 1.3)
            & (np.abs(jet1[:,1].eta) < 1.3)
        ))
        if jet2 is not None:
            eta_cut *= ak.to_numpy((
                (np.abs(jet2[:,0].eta) < 1.3)
                & (np.abs(jet2[:,1].eta) < 1.3)
                ))
            
    elif eta == "endcap" or eta == 1:
        eta_cut = ak.to_numpy((
            ((np.abs(jet1[:,0].eta) > 1.3) & (np.abs(jet1[:,0].eta) < 2.5))
            & ((np.abs(jet1[:,1].eta) > 1.3) & (np.abs(jet1[:,1].eta) < 2.5))
        ))
        if jet2 is not None:
            eta_cut *= ak.to_numpy((
            ((np.abs(jet2[:,0].eta) > 1.3) & (np.abs(jet2[:,0].eta) < 2.5))
            & ((np.abs(jet2[:,1].eta) > 1.3) & (np.abs(jet2[:,1].eta) < 2.5))
        ))
    return eta_cut

def compute_asymmetry(jets, pt_type="pt"):
    shuffle = random.choices([-1, 1], k=len(jets[:,0]))
    shuffle_opp = [s * -1 for s in shuffle]

    asymmetry = (shuffle * jets[:,0][pt_type] + shuffle_opp * jets[:,1][pt_type]) / (jets[:,0][pt_type] + jets[:,1][pt_type])
    asymmetry_cut = (np.abs(asymmetry) < 0.5)

    return asymmetry[asymmetry_cut], asymmetry_cut


def add_pileup_weight(events: ak.Array, pileup_profile: str = None) -> ak.Array:
    path_pileup_profile = os.path.join(
        os.path.dirname(__file__),
        "../data/pileup_histograms/PileupHistogram_Data2022.root",
    )
    pileup_profile = uproot.open(path_pileup_profile)["pileup"]
    pileup_profile = pileup_profile.to_numpy()[0]
    pileup_profile /= pileup_profile.sum()

    pileup_MC = np.histogram(ak.to_numpy(events.Pileup.nPU), bins=100, range=(0, 100))[0].astype("float64")
    # avoid division by zero later
    pileup_MC[pileup_MC == 0.] = 1
    pileup_MC /= pileup_MC.sum()

    pileup_correction = pileup_profile / pileup_MC
    # remove large MC reweighting factors to prevent artifacts
    pileup_correction[pileup_correction > 10] = 1

    weight_pileup = pileup_correction[ak.to_numpy(events.Pileup.nPU)]
    events["weight_pileup"] = weight_pileup
    
    return events
