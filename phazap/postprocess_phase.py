import os
import numpy as np
import h5py
import bilby
from tqdm import tqdm

from . import logger, _default_postprocessed_phase_dir, _default_postprocessed_phase_filename_str
from . import gwphase
from . import gw_utils as gwutils
from .pe_input import ParameterEstimationInput

_internal_counter = 1

def postprocess_phase(
        pe_result,
        flow=20.,
        fhigh=100.,
        fbest=40.,
        label=None,
        output_dir=_default_postprocessed_phase_dir,
        output_filename=None,
    ):
    """
    #sname is the superevent name
    #fbest is the pivotal frequency
    #fhigh is the high frequency for the phase evolution
    #flow is the low frequency for the phase evolution
    """

    if type(pe_result) is ParameterEstimationInput:
        ip = pe_result # No need to do anything
    elif type(pe_result) is str:
        try:
            # Maybe it is a bilby result file?
            ip = ParameterEstimationInput.from_bilby_result_file(pe_result)
        except:
            try:
                # Maybe it is a PESummary file?
                ip = ParameterEstimationInput.from_PESummary_file(pe_result)
            except:
                raise ValueError(f"Does not understand {pe_result}")

    posterior_samples = ip.posterior_samples
    fref = ip.reference_frequency
    duration = ip.duration 
    sampling_frequency = ip.sampling_frequency
    ifos = ip.ifo_list
    logger.info(f"Detectors online {ifos}")
    approx = ip.waveform_approximant
    logger.info(f"Waveform approximant {approx}")

    #15D
    phiRef = posterior_samples['phase']
    theta_jn = posterior_samples['theta_jn']
    ra = posterior_samples['ra']
    dec = posterior_samples['dec']
    psi = posterior_samples['psi']
    geocent_time = posterior_samples['geocent_time']
    a1 = posterior_samples['a_1'] 
    a2 = posterior_samples['a_2'] 
    tilt_1 = posterior_samples['tilt_1'] 
    tilt_2 = posterior_samples['tilt_2'] 
    phi_jl = posterior_samples['phi_jl'] 
    phi_12 = posterior_samples['phi_12'] 
    mass_1 = posterior_samples['mass_1'] 
    mass_2 = posterior_samples['mass_2']
    luminosity_distance = posterior_samples['luminosity_distance'] 

    N_posteriors = len(mass_1)

    modes = [[2,2],[2,-2]]

    waveform_arguments = dict(
            waveform_approximant=approx,
            reference_frequency=fref,
            minimum_frequency=min(fref,flow) - 2./duration,
            maximum_frequency=fhigh + 2./duration,
            mode_array = modes
        )

    waveform_generator = bilby.gw.WaveformGenerator(
            duration=duration,
            sampling_frequency=sampling_frequency,
            frequency_domain_source_model=bilby.gw.source.lal_binary_black_hole,
            parameter_conversion=bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters,
            waveform_arguments=waveform_arguments,
        )

    phase_fbest = np.zeros(N_posteriors)
    a22_fbest = np.zeros(N_posteriors)
    zeta_fbest = np.zeros(N_posteriors)
    r_fbest = np.zeros(N_posteriors)
    phase_fhigh = np.zeros(N_posteriors)
    phase_flow = np.zeros(N_posteriors)

    fs = bilby.core.utils.series.create_frequency_series(sampling_frequency,duration)
    rang_fs = (fs >= flow) & (fs <= fhigh)
    i_best = np.argmin(abs(fs[rang_fs]-fbest))

    for i in tqdm(range(N_posteriors)):
        phase_fbest[i], a22_fbest[i], zeta_fbest[i], r_fbest[i], phase_fhigh[i], phase_flow[i] = gwphase.phases_fs(waveform_generator,
                                        detec = "H1",
                                        rang_fs = rang_fs,
                                        i_best = i_best,
                                        mass_1=mass_1[i],
                                        mass_2=mass_2[i],
                                        a_1=a1[i],
                                        a_2=a2[i],
                                        tilt_1=tilt_1[i],
                                        tilt_2=tilt_2[i],
                                        phi_12=phi_12[i],
                                        phi_jl=phi_jl[i],
                                        luminosity_distance=luminosity_distance[i],
                                        theta_jn=theta_jn[i],
                                        psi = psi[i],
                                        phase_ref= phiRef[i],
                                        geocent_time= geocent_time[i],
                                        ra= ra[i],
                                        dec= dec[i],
                                        duration = duration,
                                        sampling_frequency = sampling_frequency)

    time_delay_H = np.zeros(len(geocent_time))
    time_delay_L = np.zeros(len(geocent_time))
    time_delay_V = np.zeros(len(geocent_time))
    for i in range(len(geocent_time)):
        time_delay_H[i] = gwutils.time_delay_det("H1", ra[i], dec[i], geocent_time[i])
        time_delay_L[i] = gwutils.time_delay_det("L1", ra[i], dec[i], geocent_time[i])
        time_delay_V[i] = gwutils.time_delay_det("V1", ra[i], dec[i], geocent_time[i])

    #Arrival time phase at each detector
    tau_H = 2.*np.pi*fbest*(geocent_time + time_delay_H)
    tau_L = 2.*np.pi*fbest*(geocent_time + time_delay_L)
    tau_V = 2.*np.pi*fbest*(geocent_time + time_delay_V)

    #Arrival time phase difference between detectors
    tau_HL = tau_H-tau_L
    tau_HV = tau_H-tau_V
    tau_LV = tau_L-tau_V

    #Detector phases
    Fp_H, Fx_H = gwutils.FpFx("H1",ra, dec, psi, geocent_time)
    Fp_L, Fx_L = gwutils.FpFx("L1",ra, dec, psi, geocent_time)
    Fp_V, Fx_V = gwutils.FpFx("V1",ra, dec, psi, geocent_time)

    phase_H = gwphase.phase_d(phase_fbest,a22_fbest,zeta_fbest,Fp_H,Fx_H)
    phase_L = gwphase.phase_d(phase_fbest,a22_fbest,zeta_fbest,Fp_L,Fx_L)
    phase_V = gwphase.phase_d(phase_fbest,a22_fbest,zeta_fbest,Fp_V,Fx_V)

    #Detector phase evolution in frequency
    Dphi_f = phase_fhigh - phase_flow

    #Post-processed phases
    variables_event_1 = ['phase_H', 'phase_L', 'phase_V', 'Dphi_f', 'tau_HL', 'tau_HV', 'tau_LV']

    #For event 2 we need to compute the phases in the frame of event 1
    variables_event_2 = ['phase_fbest','a22_fbest','zeta_fbest','r_fbest','phase_fhigh','phase_flow','tau_H','tau_L','tau_V',
                         'ra', 'dec', 'psi', 'geocent_time']

    variables = [*variables_event_1, *variables_event_2]

    # Saving the data
    if label is None:
        global _internal_counter # Not a very good practice
        label = f"event{_internal_counter}"
        _internal_counter += 1

    if output_filename is None:
        output_filename = _default_postprocessed_phase_filename_str.format(
            label,
            fbest,
            fhigh,
            flow
        )

    output_fullpath = os.path.join(output_dir, output_filename)
    with h5py.File(output_fullpath, "w") as f:
        for var in variables:
            dset = f.create_dataset(var, data=eval(var))
        
        f.attrs["fbest"] = fbest
        f.attrs["fhigh"] = fhigh
        f.attrs["flow"] = flow
        f.attrs["label"] = label
    
    logger.info(f"Postprocessing completed and saved to {output_fullpath}")