import numpy as np

from . import utils
from . import tension_utils as tension

#PLOTS
import getdist
from getdist import plots
import matplotlib.pyplot as plt
#%matplotlib inline
#%config InlineBackend.figure_format='retina'
#mpl.rcParams['text.usetex'] = True
#mpl.rcParams['font.family'] = 'serif'

def phazap_plot(event_name_1,event_name_2,fbest=40.0,fhigh=100.0,flow=20.0,dir_phase = 'phazap_phases_o4/',dir_out = 'phazap_plots/'):

    parameters_1, parameters_2, det_phases_1, det_phases_2, tau_phases_1, tau_phases_2, Dphi_f_1, Dphi_f_2, above_below = phazap.phases_events(event_name_1,event_name_2,fbest,fhigh,flow,dir_phase)

    nphases = 3
    length = np.shape(parameters_1)[1]
    
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

    dist = dist_all

    phase_shifts = np.array([0,1,2,-1,-2])*np.pi/2
    phase_shift = phase_shifts[np.argmin(dist)]*np.hstack((np.ones(nphases),np.zeros(length-nphases)))

    #Event 2 is shifted by the phase shift that minimizes the distance
    parameters_2_shifted = parameters_2 + phase_shift
    parameters_2_shifted[:,:nphases] = np.mod(parameters_2_shifted[:,:nphases],2*np.pi)

    parameters_1_wrap, parameters_2_wrap = utils.wrap_phases(parameters_1,parameters_2_shifted,nphases)

    #labels
    labels = [r'\phi_\mathrm{H}',
              r'\phi_\mathrm{L}',
              r'\phi_\mathrm{V}',
              r'\Delta\phi_f',
             r'\tau_\mathrm{HL}',
                 r'\tau_\mathrm{HV}']

    label_1 = event_name_1
    label_2 = event_name_2+' (shifted by $%s\pi$)' % np.round(phase_shifts[np.argmin(dist)]/np.pi,1)

    color_1 = 'blue'
    color_2 = 'green'

    # Get the getdist MCSamples objects for the samples, specifying same parameter
    # names and labels; if not specified weights are assumed to all be unity
    names = ["x%s"%i for i in range(length)]
    #labels =  ["x_%s"%i for i in range(ndim)]
    samples_1 = getdist.MCSamples(samples=parameters_1_wrap,names=names, labels = labels,label=label_1)
    samples_2 = getdist.MCSamples(samples=parameters_2_wrap,names=names, labels = labels,label=label_2)


    # Start plotting
    g = plots.get_subplot_plotter(width_inch=6, scaling=False, rc_sizes=True)

    # Settings
    g.settings.figure_legend_frame = False

    # Triangle plot
    g.triangle_plot([samples_1, samples_2], 
                    filled=True, 
                    colors = [color_1, color_2],
                   line_args = [{'color':color_1},{'color':color_2}])

    plt.suptitle(r'Minimum distance =%s' % np.round(np.min(dist),2)+r', volume = %s' % np.round(np.min(vol_phases_1),2), va='bottom')

    plt.savefig(dir_out+'phazap_'+event_name_1+'_'+event_name_2+'_fbest_%s_fhigh_%s_flow_%s.pdf' % (int(fbest),int(fhigh),int(flow)), bbox_inches='tight', transparent=True)


    return print(event_name_1,event_name_2,' plot done!')