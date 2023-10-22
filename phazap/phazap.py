import numpy as np
import h5py

import gwphase
import gw_utils as gwutils
import tension_utils as tension

def phazap_all_metrics_one_ordering(event_name_1,event_name_2,fbest=40.0,fhigh=100.0,flow=20.0,dir_phase = 'phazap_phases_o4/'):

    parameters_1, parameters_2, det_phases_1, det_phases_2, tau_phases_1, tau_phases_2, Dphi_f_1, Dphi_f_2, above_below = phases_events(event_name_1,event_name_2,fbest,fhigh,flow,dir_phase)

    nphases = 3

    #Record number of effective phases
    prior_range = np.pi/2
    neff = tension.n_eff_pair(det_phases_1,det_phases_2,nphases,prior_range)

    #Compute volumes
    vol_1 = tension.volume_phases(parameters_1,nphases)
    vol_2 = tension.volume_phases(parameters_2,nphases)        
    vol_phases_1, vol_phases_2 = tension.volume_only_phases(det_phases_1,det_phases_2)

    #Compute distances
    #We divide the samples in two groups, above and below the plane
    #only if the fraction of samples in each group is between 0.05 and 0.95
    frac_samples_below = len(above_below[above_below<0]) / len(above_below)
    frac_limit = 0.05
    if (frac_samples_below>frac_limit) & (frac_samples_below < 1 - frac_limit):
        #All phases
        dist_M = tension.distance_phases_with_shift(parameters_1,parameters_2[above_below<0],nphases)
        dist_P = tension.distance_phases_with_shift(parameters_1,parameters_2[above_below>0],nphases)
        dist_all = np.minimum(dist_M,dist_P)
        #Time delay phases
        dist_Tphases_M = tension.distance(tau_phases_1,tau_phases_2[above_below<0])
        dist_Tphases_P = tension.distance(tau_phases_1,tau_phases_2[above_below>0])
        dist_Tphases = np.minimum(dist_Tphases_M,dist_Tphases_P)
        #Detector phases
        dist_phases_M = tension.distance_only_phases_with_shift(det_phases_1,det_phases_2[above_below<0])
        dist_phases_P = tension.distance_only_phases_with_shift(det_phases_1,det_phases_2[above_below>0])
        dist_phases = np.minimum(dist_phases_M,dist_phases_P)
    else:
        dist_all = tension.distance_phases_with_shift(parameters_1,parameters_2,nphases)
        dist_phases = tension.distance_only_phases_with_shift(det_phases_1,det_phases_2)
        dist_Tphases = tension.distance(tau_phases_1,tau_phases_2)

    dist_Dphi = tension.distance1D(Dphi_f_1,Dphi_f_2)

    return neff, vol_1, vol_2, vol_phases_1, vol_phases_2, dist_all, dist_phases, dist_Tphases, dist_Dphi

def phazap_one_ordering(event_name_1,event_name_2,fbest=40.0,fhigh=100.0,flow=20.0,dir_phase = 'phazap_phases_o4/'):

    parameters_1, parameters_2, det_phases_1, det_phases_2, tau_phases_1, tau_phases_2, Dphi_f_1, Dphi_f_2, above_below = phases_events(event_name_1,event_name_2,fbest,fhigh,flow,dir_phase)

    nphases = 3
    #Compute volume phases     
    vol_phases_1, vol_phases_2 = tension.volume_only_phases(det_phases_1,det_phases_2)

    #Compute distances
    #We divide the samples in two groups, above and below the plane
    #only if the fraction of samples in each group is between 0.05 and 0.95
    frac_samples_below = len(above_below[above_below<0]) / len(above_below)
    frac_limit = 0.05
    if (frac_samples_below>frac_limit) & (frac_samples_below < 1 - frac_limit):
        #All phases
        dist_M = tension.distance_phases_with_shift(parameters_1,parameters_2[above_below<0],nphases)
        dist_P = tension.distance_phases_with_shift(parameters_1,parameters_2[above_below>0],nphases)
        dist_all = np.minimum(dist_M,dist_P)
    else:
        dist_all = tension.distance_phases_with_shift(parameters_1,parameters_2,nphases)

    #ADD NEFF

    return vol_phases_1, dist_all

def phazap(event_name_1,event_name_2,fbest=40.0,fhigh=100.0,flow=20.0,dir_phase = 'phazap_phases_o4/'):
    #Compute volumes and distances for both orderings
    vol_phases_12, dist_12 = phazap_one_ordering(event_name_1,event_name_2,fbest,fhigh,flow,dir_phase)
    vol_phases_21, dist_21 = phazap_one_ordering(event_name_2,event_name_1,fbest,fhigh,flow,dir_phase)

    dist_21_swap = np.array([dist_21[0],dist_21[3],dist_21[4],dist_21[1],dist_21[2]])
    D_J_n = np.maximum(dist_12,dist_21_swap)
    D_J = np.min(D_J_n)
    
    vol_J = vol_phases_12 + vol_phases_21

    #Phase shift
    phase_shifts = np.array([0,1,2,-1,-2])*np.pi/2
    phase_shift = phase_shifts[np.argmin(D_J_n)]

    #ADD P-VALUE

    return D_J, vol_J, phase_shift, D_J_n

def phazap_summary(event_name_1,event_name_2,fbest=40.0,fhigh=100.0,flow=20.0,dir_phase = 'phazap_phases_o4/'):
    D_J, vol_J, phase_shift, D_J_n = phazap(event_name_1,event_name_2,fbest,fhigh,flow,dir_phase)

    return D_J, vol_J, phase_shift