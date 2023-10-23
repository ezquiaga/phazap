import numpy as np
import os

from . import gwphase
from . import gw_utils as gwutils
from . import tension_utils as tension
from .postprocess_phase import postprocess_phase, PostprocessedPhase, _variables_event_1, _variables_event_2

def phases_events(event1_postprocessed_phase, event2_postprocessed_phase):
    assert event1_postprocessed_phase.fbest == event2_postprocessed_phase.fbest, "The two sets of phases have different fbest"
    fbest = event1_postprocessed_phase.fbest # Either one should be fine
    assert event1_postprocessed_phase.flow == event2_postprocessed_phase.flow, "The two sets of phases have different flow"
    flow = event1_postprocessed_phase.flow
    assert event1_postprocessed_phase.fhigh == event2_postprocessed_phase.fhigh, "The two sets of phases have different fhigh"
    fhigh = event1_postprocessed_phase.fhigh

    # Define variables *explicitly* so that linter will be satisfied
    # Read from memory
    geocent_time_1 = event1_postprocessed_phase.dataset['geocent_time']
    phase_H_1 = event1_postprocessed_phase.dataset['phase_H']
    phase_L_1 = event1_postprocessed_phase.dataset['phase_L']
    phase_V_1 = event1_postprocessed_phase.dataset['phase_V']
    tau_HL_1 = event1_postprocessed_phase.dataset['tau_HL']
    tau_HV_1 = event1_postprocessed_phase.dataset['tau_HV']
    Dphi_f_1 = event1_postprocessed_phase.dataset['Dphi_f']

    geocent_time_2 = event2_postprocessed_phase.dataset['geocent_time']
    ra_2 = event2_postprocessed_phase.dataset['ra']
    dec_2 = event2_postprocessed_phase.dataset['dec']
    psi_2 = event2_postprocessed_phase.dataset['psi']
    phase_fbest_2 = event2_postprocessed_phase.dataset['phase_fbest']
    a22_fbest_2 = event2_postprocessed_phase.dataset['a22_fbest']
    zeta_fbest_2 = event2_postprocessed_phase.dataset['zeta_fbest']
    Dphi_f_2 = event2_postprocessed_phase.dataset['Dphi_f']


    Nsamples_2 = len(phase_fbest_2)
    tc_shift_21 = np.mean(geocent_time_2) - np.mean(geocent_time_1)
    
    above_below = np.zeros(Nsamples_2)

    tau_H_2_wrt_1 = np.zeros(Nsamples_2)
    tau_L_2_wrt_1 = np.zeros(Nsamples_2)
    tau_V_2_wrt_1 = np.zeros(Nsamples_2)
    for l in range(Nsamples_2):
        tau_H_2_wrt_1[l] = 2.*np.pi*fbest*gwutils.time_delay_det("H1", ra_2[l], dec_2[l], geocent_time_2[l]-tc_shift_21)
        tau_L_2_wrt_1[l] = 2.*np.pi*fbest*gwutils.time_delay_det("L1", ra_2[l], dec_2[l], geocent_time_2[l]-tc_shift_21)
        tau_V_2_wrt_1[l] = 2.*np.pi*fbest*gwutils.time_delay_det("V1", ra_2[l], dec_2[l], geocent_time_2[l]-tc_shift_21)
         #Define quadrant of the source
        above_below[l] = gwutils.N_dot_cross_d123(ra_2[l],dec_2[l],geocent_time_2[l],gwutils.H1_vertex_meters,gwutils.L1_vertex_meters,gwutils.V1_vertex_meters)
    
    #Arrival time phase difference between detectors
    tau_HL_2 = tau_H_2_wrt_1-tau_L_2_wrt_1
    tau_HV_2 = tau_H_2_wrt_1-tau_V_2_wrt_1
    tau_LV_2 = tau_L_2_wrt_1-tau_V_2_wrt_1
    
    #Detector phases at reference frame of event 1
    Fp_H_2, Fx_H_2 = gwutils.FpFx("H1",ra_2, dec_2, psi_2, geocent_time_2 - tc_shift_21)
    Fp_L_2, Fx_L_2 = gwutils.FpFx("L1",ra_2, dec_2, psi_2, geocent_time_2 - tc_shift_21)
    Fp_V_2, Fx_V_2 = gwutils.FpFx("V1",ra_2, dec_2, psi_2, geocent_time_2 - tc_shift_21)
    
    phase_H_2 = gwphase.phase_d(phase_fbest_2,a22_fbest_2,zeta_fbest_2,Fp_H_2,Fx_H_2)
    phase_L_2 = gwphase.phase_d(phase_fbest_2,a22_fbest_2,zeta_fbest_2,Fp_L_2,Fx_L_2)
    phase_V_2 = gwphase.phase_d(phase_fbest_2,a22_fbest_2,zeta_fbest_2,Fp_V_2,Fx_V_2)
    
    #Sets of parameters
    #All phases
    parameters_1 = np.stack(
        (
            phase_H_1,
            phase_L_1,
            phase_V_1,
            tau_HL_1,
            tau_HV_1,
            Dphi_f_1
        ),
        axis=1
    )
    
    parameters_2 = np.stack(
        (
            phase_H_2,
            phase_L_2,
            phase_V_2,
            tau_HL_2,
            tau_HV_2,
            Dphi_f_2
        ),
        axis=1
    )
    
    #Detector phases
    det_phases_1 = np.stack(
        (
            phase_H_1,
            phase_L_1,
            phase_V_1
        ),
        axis=1
    )
    
    det_phases_2 = np.stack(
        (
            phase_H_2,
            phase_L_2,
            phase_V_2
        ),
        axis=1
    )
    
    #Time delay phases
    tau_phases_1 = np.stack(
        (
            tau_HL_1,
            tau_HV_1
        ),
        axis=1
    )
    
    tau_phases_2 = np.stack(
        (
            tau_HL_2,
            tau_HV_2
        ),
        axis=1
    )
    
    return parameters_1, parameters_2, det_phases_1, det_phases_2, tau_phases_1, tau_phases_2, Dphi_f_1, Dphi_f_2, above_below

def phazap_all_metrics_one_ordering(event1_postprocessed_phase, event2_postprocessed_phase):
    parameters_1, parameters_2, det_phases_1, det_phases_2, tau_phases_1, tau_phases_2, Dphi_f_1, Dphi_f_2, above_below = phases_events(event1_postprocessed_phase, event2_postprocessed_phase)

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

def phazap_one_ordering(event1_postprocessed_phase, event2_postprocessed_phase):

    parameters_1, parameters_2, det_phases_1, det_phases_2, tau_phases_1, tau_phases_2, Dphi_f_1, Dphi_f_2, above_below = phases_events(event1_postprocessed_phase, event2_postprocessed_phase)

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


def _phazap(event1_postprocessed_phase, event2_postprocessed_phase):
    #Compute volumes and distances for both orderings
    vol_phases_12, dist_12 = phazap_one_ordering(event1_postprocessed_phase, event2_postprocessed_phase)
    vol_phases_21, dist_21 = phazap_one_ordering(event1_postprocessed_phase, event2_postprocessed_phase)

    dist_21_swap = np.array([dist_21[0],dist_21[3],dist_21[4],dist_21[1],dist_21[2]])
    D_J_n = np.maximum(dist_12,dist_21_swap)
    D_J = np.min(D_J_n)
    
    vol_J = vol_phases_12 + vol_phases_21

    #Phase shift
    phase_shifts = np.array([0,1,2,-1,-2])*np.pi/2
    phase_shift = phase_shifts[np.argmin(D_J_n)]

    #ADD P-VALUE

    return D_J, vol_J, phase_shift, D_J_n

def phazap(event_1, event_2):
    def check_if_postprocessed(x):
        if type(x) is PostprocessedPhase:
            # x is indeed a PostprocessedPhase object, return it directly
            return x
        if type(x) is str:
            if os.path.exists(x):
                # x is a file path
                try:
                    # Maybe it is already postprocessed?
                    return PostprocessedPhase.from_file(x)
                except:
                    try:
                        # Maybe it is a summary file?
                        return postprocess_phase(x)
                    except:
                        raise ValueError(f"Does not understand {x}")
            else:
                raise NotImplemented("Currently only support a file path as the value")
    
    postprocessed_phase_1 = check_if_postprocessed(event_1)
    postprocessed_phase_2 = check_if_postprocessed(event_2)

    return _phazap(postprocessed_phase_1, postprocessed_phase_2)

def phazap_summary(event_1, event_2):
    D_J, vol_J, phase_shift, D_J_n = phazap(event_1, event_2)

    return D_J, vol_J, phase_shift
