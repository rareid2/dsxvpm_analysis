import matplotlib.pyplot as plt 
import numpy as np
import datetime as dt 
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from parula_colormap import parula
import scipy.signal

def plot_TNT(subplot, tnttimes, durations, startf, stopf): 
    ax = subplot
    for di, dur in enumerate(durations):
        tim = tnttimes[di]
        start = startf[di]
        stop = stopf[di]
        if dur == 150:
            ax.plot([tim, tim + dt.timedelta(microseconds=75e3)],[start/1e3,(start+100)/1e3],'b')
            ax.plot([tim + dt.timedelta(microseconds=75e3),tim + dt.timedelta(microseconds=150e3)],[(start+100)/1e3,(start-100)/1e3],'b')
        else:
            ax.plot([tim, tim + dt.timedelta(microseconds=dur*1e3)],[start/1e3,stop/1e3],'b')

    ax.set_ylabel('Frequency [kHz]')


def plot_burst(fig, subplot, burst):
    cfg = burst['config']
    E_FD = subplot

    system_delay_samps_FD = 200;
    fs = 80000;
    cm = parula();  # This is a mockup of the current Matlab colormap (which is proprietary)

    e_clims = np.array([-40, 10]) # for newly calibrated data
    E_coef = burst['CAL'] # calibrate into uV/m units 

    # Construct the appropriate time and frequency axes
    # Get the equivalent sample rate, if decimated
    if cfg['DECIMATE_ON']==1:
        fs_equiv = 80000./cfg['DECIMATION_FACTOR']
    else:
        fs_equiv = 80000.

    if cfg['SAMPLES_OFF'] == 0:
        max_ind = max(len(burst['E']), len(burst['B']))
        t_axis = np.arange(max_ind)/fs_equiv
    else:
        # Seconds from the start of the burst
        t_axis = np.array([(np.arange(cfg['SAMPLES_ON']))/fs_equiv +\
                        (k*(cfg['SAMPLES_ON'] + cfg['SAMPLES_OFF']))/fs_equiv for k in range(cfg['burst_pulses'])]).ravel()

    # GPS timestamps are taken at the end of each contiguous recording.
    # (I think "samples on" is still undecimated, regardless if decimation is being used...)
    start_timestamp = dt.datetime.utcfromtimestamp(burst['G'][0]['timestamp']) - dt.timedelta(seconds=float(cfg['SAMPLES_ON']/fs))

    # the "samples on" and "samples off" values are counting at the full rate, not the decimated rate.
    sec_on  = cfg['SAMPLES_ON']/fs
    sec_off = cfg['SAMPLES_OFF']/fs
    
    nfft=1024;
    overlap = 0.5
    window = 'hanning'

    if cfg['SAMPLES_OFF'] == 0:
        E_td_spaced = E_coef*burst['E']
    else:
        # Insert nans into vector to account for "off" time sections
        E_td_spaced = []

        for k in np.arange(cfg['burst_pulses']):
            E_td_spaced.append(E_coef*burst['E'][k*cfg['SAMPLES_ON']:(k+1)*cfg['SAMPLES_ON']])
            E_td_spaced.append(np.ones(cfg['SAMPLES_OFF'])*np.nan)

        E_td_spaced = np.concatenate(E_td_spaced).ravel()

    # E spectrogram -- "spectrum" scaling -> V^2; "density" scaling -> V^2/Hz
    ff,tt, FE = scipy.signal.spectrogram(E_td_spaced, fs=fs_equiv, window=window,
                nperseg=nfft, noverlap=nfft*overlap,mode='psd',scaling='density') # changed to density
    E_S_mag = 20*np.log10(np.sqrt(FE))
    E_S_mag[np.isinf(E_S_mag)] = -100

    # change time axis
    dt_axis = [start_timestamp + dt.timedelta(seconds=t) for t in tt]

    pe = E_FD.pcolormesh(dt_axis, ff/1000, E_S_mag, cmap = cm,  vmin=e_clims[0], vmax=e_clims[1])

    left, bottom, width, height = E_FD.get_position().bounds
    cax = fig.add_axes([left*7.28, bottom, width * 0.02, height])
    ce = fig.colorbar(pe, orientation='vertical', cax=cax)
    ce.set_label(f'dB[(uV/m)^2/Hz]')
    E_FD.set_ylabel('Frequency [kHz]')
    E_FD.set_xlabel('Time')

    return start_timestamp