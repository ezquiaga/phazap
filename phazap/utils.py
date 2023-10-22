import numpy as np

""" Compute mode """

#1D Maximum likelihood value from posterior samples
def mode(posterior,bins=30):
    bins, param = np.histogram(posterior,bins=bins)
    return param[np.argmax(bins)]

def modes(posteriors,bins=30):
    length = np.shape(posteriors)[-1]
    
    modes_c = np.zeros(length)
    for i in range(length):
        modes_c[i] = mode(posteriors[:,i],bins)
    return modes_c

""" Wrap phases around their modes """

def wrap_phases(parameters_1,parameters_2,nphases):
    #Wrap phases around their modes 
    #to avoid discontinuities from periodic boundaries
    phases_1 = parameters_1[:,:nphases]
    phases_2 = parameters_2[:,:nphases]
    
    mode_1 = modes(phases_1)
    mode_2 = modes(phases_2)
   
    phases_1_wrap = np.mod(phases_1 - mode_1 + np.pi,2*np.pi)
    phases_2_wrap = np.mod(phases_2 - mode_2 + np.pi,2*np.pi) + mode_2 - mode_1
    
    #Put with rest of parameters
    parameters_1_wrap = np.concatenate((phases_1_wrap,parameters_1[:,nphases:]),axis=1)
    parameters_2_wrap = np.concatenate((phases_2_wrap,parameters_2[:,nphases:]),axis=1)
    
    mean_1_wrap = np.mean(parameters_1_wrap,axis=0)
    mean_2_wrap = np.mean(parameters_2_wrap,axis=0)
    
    #Make sure mean is within 0, 2pi
    parameters_1_mod = np.concatenate((np.mod(phases_1_wrap,2*np.pi),parameters_1[:,nphases:]),axis=1)
    parameters_2_mod = np.concatenate((np.mod(phases_2_wrap,2*np.pi),parameters_2[:,nphases:]),axis=1)
    
    mean_1_mod = np.mean(parameters_1_mod,axis=0)
    mean_2_mod = np.mean(parameters_2_mod,axis=0)
    
    parameters_1_wrap_mod = parameters_1_wrap - mean_1_wrap + mean_1_mod
    parameters_2_wrap_mod = parameters_2_wrap - mean_2_wrap + mean_2_mod
    
    return parameters_1_wrap_mod, parameters_2_wrap_mod

def wrap_only_phases(phases_1,phases_2):
    #Wrap phases around their modes 
    #to avoid discontinuities from periodic boundaries
    
    mode_1 = modes(phases_1)
    mode_2 = modes(phases_2)
   
    phases_1_wrap = np.mod(phases_1 - mode_1 + np.pi,2*np.pi)
    phases_2_wrap = np.mod(phases_2 - mode_2 + np.pi,2*np.pi) + mode_2 - mode_1   
    
    mean_2_wrap = np.mean(phases_2_wrap,axis=0)
    
    #Make sure mean is within 0, 2pi
    phases_2_mod = np.mod(phases_2_wrap,2*np.pi)
    
    mean_2_mod = np.mean(phases_2_mod,axis=0)
    
    phases_2_wrap_mod = phases_2_wrap - mean_2_wrap + mean_2_mod
    
    return phases_1_wrap, phases_2_wrap_mod

def wrap_phase(phases_1):
    #Wrap phases around their modes 
    #to avoid discontinuities from periodic boundaries
    
    mode_1 = mode(phases_1)
   
    phases_1_wrap = np.mod(phases_1 - mode_1 + np.pi,2*np.pi)
    
    return phases_1_wrap