import numpy as np
import lal
import bilby

DETECTORS = {'H1': lal.CachedDetectors[lal.LHO_4K_DETECTOR],
             'L1': lal.CachedDetectors[lal.LLO_4K_DETECTOR],
             'V1': lal.CachedDetectors[lal.VIRGO_DETECTOR]}

def time_delay_det(detec,ra, dec, geocent_time):
    #det: detector
    #ra: right ascension
    #dec: declination
    #geocent_time: geocentric time
    return lal.TimeDelayFromEarthCenter(DETECTORS[detec].location, ra, dec, geocent_time)

def FpFx(det,ra, dec, psi, geocent_time):
    #det: detector
    #ra: right ascension
    #dec: declination
    #psi: polarization
    #geocent_time: geocentric time
    gmst = [lal.GreenwichMeanSiderealTime(t) for t in geocent_time]      
    return np.transpose([lal.ComputeDetAMResponse(DETECTORS[det].response, r, d, p, g)
        for r, d, p, g in np.broadcast(ra, dec, psi, gmst)])

#https://docs.ligo.org/lscsoft/lalsuite/lal/group___detector_constants.html
#Position of the detectors vertex in meters with respect to geocenter
H1_vertex_meters = np.array([lal.LHO_4K_VERTEX_LOCATION_X_SI,
                            lal.LHO_4K_VERTEX_LOCATION_Y_SI,
                            lal.LHO_4K_VERTEX_LOCATION_Z_SI])

L1_vertex_meters = np.array([lal.LLO_4K_VERTEX_LOCATION_X_SI,
                            lal.LLO_4K_VERTEX_LOCATION_Y_SI,
                            lal.LLO_4K_VERTEX_LOCATION_Z_SI])

V1_vertex_meters = np.array([lal.VIRGO_VERTEX_LOCATION_X_SI,
                            lal.VIRGO_VERTEX_LOCATION_Y_SI,
                            lal.VIRGO_VERTEX_LOCATION_Z_SI])

def Nvector(source_right_ascension_radians,source_declination_radians,gpstime):
    #compute the unit vector pointing from the the source to geocenter
    #in LAL they compute -N
    
    gmst = lal.GreenwichMeanSiderealTime(gpstime)
    greenwich_hour_angle = gmst - source_right_ascension_radians
    
    mN_x = np.cos(source_declination_radians) * np.cos(greenwich_hour_angle)
    mN_y = np.cos(source_declination_radians) * (-np.sin(greenwich_hour_angle))
    mN_z = np.sin(source_declination_radians)
    return np.array([-mN_x,-mN_y,-mN_z])

def N_dot_cross_d123(ra,dec,gpstime,d1,d2,d3):
    N = Nvector(ra,dec,gpstime)
    d12 = d1 - d2
    d13 = d1 - d3
    return np.dot(N, np.cross(d12,d13))

def N_dot_d12(ra,dec,gpstime,d1,d2):
    N = Nvector(ra,dec,gpstime)
    d12 = d1 - d2
    return np.dot(N, d12) / lal.C_SI

"""Generate waveform"""
def hp_hx(binary_parameters,waveform_generator):
    generate_h = waveform_generator.frequency_domain_strain(binary_parameters)

    hp = generate_h['plus']
    hx = generate_h['cross']
    
    return hp, hx
