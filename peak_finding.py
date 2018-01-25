import functools
import operator
import numpy as np
import scipy.stats as sci
from scipy import signal
import seaborn as sns
import pandas as pd
import matplotlib.gridspec as gridspec
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from datetime import datetime
import peakutils.peak as pk


def measure_len_runs_of_nans_array(bits):
    # make sure all runs of ones are well-bounded
    bounded = np.hstack(([0], bits, [0]))
    # get 1 at run starts and -1 at run ends
    difs = np.diff(bounded)
    run_starts, = np.where(difs > 0)
    run_ends, = np.where(difs < 0)
    return zip(run_starts,run_ends)

def return_clean_malariatherapy_trajectories(data,lapse):
    df = data
    asexual_trajectories = {}
    gametocyte_trajectories = {}

    for id,measurements in df.groupby('patient'):
        try:
            asex = np.array([np.nan if x == '.' else int(x) for x in np.array(measurements.Asexual)])
            gam = np.array([np.nan if x == '.' else int(x) for x in np.array(measurements.Gametocytes)])

            mask_as = np.isnan(asex)
            mask_gam = np.isnan(gam)

            asex[mask_as] = np.interp(np.flatnonzero(mask_as),np.flatnonzero(~mask_as),asex[~mask_as])
            gam[mask_gam] = np.interp(np.flatnonzero(mask_gam), np.flatnonzero(~mask_gam), gam[~mask_gam])

            masker_gam = np.append((np.diff(np.flatnonzero(mask_gam), n=2) == 0), 0)
            masker_as = np.append((np.diff(np.flatnonzero(mask_as), n=2) == 0), 0)

            nan_edges_gam = measure_len_runs_of_nans_array(masker_gam)
            nan_edges_as = measure_len_runs_of_nans_array(masker_as)

            bad_list_gam = []
            bad_list_as = []

            for i in range(len(nan_edges_gam)):
                if nan_edges_gam[i][1] - nan_edges_gam[i][0] > lapse:
                    bad_list_gam += range(nan_edges_gam[i][0], nan_edges_gam[i][1])
                else:
                    pass

            for i in range(len(nan_edges_as)):
                if nan_edges_as[i][1] - nan_edges_as[i][0] > lapse:
                    bad_list_as += range(nan_edges_as[i][0], nan_edges_as[i][1])
                else:
                    pass
            # These are the indices to hide of the array
            hide_these_indices_gam = bad_list_gam
            hide_these_indices_as = bad_list_as

            if hide_these_indices_gam:
                gam[hide_these_indices_gam] = np.nan
            if hide_these_indices_as:
                asex[hide_these_indices_as] = np.nan

            # mask_check = measurements[pd.isnull(measurements.Asexual)]
            asexual_trajectories[id] = asex
            gametocyte_trajectories[id] = gam
        except:
            print(id,sum(gam))
    return [asexual_trajectories,gametocyte_trajectories]


def smooth(x, window_len=11, window='hanning'):
    """smooth the data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.

    input:
        x: the input signal
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal

    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)

    see also:

    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter

    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."

    if window_len < 3:
        return x

    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"

    s = np.r_[x[window_len - 1:0:-1], x, x[-2:-window_len - 1:-1]]
    # print(len(s))
    if window == 'flat':  # moving average
        w = np.ones(window_len, 'd')
    else:
        w = eval('np.' + window + '(window_len)')

    y = np.convolve(w / w.sum(), s, mode='valid')
    return y


if __name__ == '__main__':
    mt_data = pd.read_table(r'C:\Users\jorussell\Dropbox (IDM)\Malaria Team Folder\data\Malariatherapy\Malariatherapy.txt')
    [asexuals,gametocytes] = return_clean_malariatherapy_trajectories(mt_data,15)

    for individual in asexuals.keys():
        try:
            asexual_densities = smooth(asexuals[individual],window_len=4)
            gametocyte_densities = smooth(gametocytes[individual],window_len=4)
# find the peak indices using a wavelet convolution
            asexual_peakind = signal.find_peaks_cwt(asexual_densities,widths = np.arange(1,30),max_distances=np.arange(1,30), min_length = 7)
            gametocyte_peakind = signal.find_peaks_cwt(gametocyte_densities,widths =  np.arange(1, 30),max_distances=np.arange(1,30), min_length=7)
            plt.figure(1)
            plt.plot(asexual_densities,color = 'red',alpha = 0.5)
            plt.plot(asexual_peakind, asexual_densities[asexual_peakind], linestyle='None', marker='o', color='red', alpha=1)
            plt.savefig(r'C:\Users\jorussell\Dropbox (IDM)\Malaria Team Folder\projects\updated_infection_and_immunity\malariatherapy\peak_identification\output_images\asexual_parasite_trace_for_patient_'+str(individual)+'.png')
            plt.clf()

            plt.figure(1)
            plt.plot(gametocyte_densities, color='blue', alpha =0.5)
            plt.plot(gametocyte_peakind,gametocyte_densities[gametocyte_peakind], linestyle='None', marker='o', color='blue', alpha=1)
            plt.savefig(
                r'C:\Users\jorussell\Dropbox (IDM)\Malaria Team Folder\projects\updated_infection_and_immunity\malariatherapy\peak_identification\output_images\gametocyte_trace_for_patient_' + str(
                    individual) + '.png')

            plt.clf()
        except:
            print(individual)
            plt.clf()
print('test')