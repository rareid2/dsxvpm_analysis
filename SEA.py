
# break up the burst into chuncks
# first let's look @ TNT
import numpy as np 
import matplotlib.pyplot as plt
import datetime as dt 
import os
from os import listdir
from os.path import isfile, join
cwd = os.getcwd()
from read_files import find_TNT, TNT_log, read_burst_XML
from plots import plot_TNT, plot_burst, plot_shift, plot_dist, plot_burst_spaced
import time
from shifts import get_lshell_sats
os.chdir(cwd)
import matplotlib.pyplot as plt 
import numpy as np
import datetime as dt 
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
import matplotlib as mpl 
from parula_colormap import parula
import scipy.signal

burst_dir = '817'
onlyfiles = [f for f in listdir(burst_dir) if isfile(join(burst_dir, f))] 

fname_burst = onlyfiles[0]
# get burst data
burst = read_burst_XML(burst_dir + '/' + fname_burst)

# get TNT log data
fname_tnt = find_TNT(fname_burst,'TNT_logs/')
tnttimes, durations, startf, stopf = TNT_log(fname_tnt)

start_timestamp = dt.datetime(2020,8,17,21,20,30, tzinfo=dt.timezone.utc)

for ti, tt in enumerate(tnttimes):
    tt = tt.replace(tzinfo=dt.timezone.utc) # cant be timezone naive
    start_timestamp = start_timestamp.replace(tzinfo=dt.timezone.utc)

    tdiff = start_timestamp - tt
    if tdiff < dt.timedelta(seconds=1):
        tnt_start = tt
        tnt_start_ind = ti
        break

#stop_timestamp = start_timestamp + dt.timedelta(seconds=burst_dur)
stop_timestamp = start_timestamp + dt.timedelta(seconds=100)

for ti, tt in enumerate(tnttimes):
    tt = tt.replace(tzinfo=dt.timezone.utc)
    start_timestamp = start_timestamp.replace(tzinfo=dt.timezone.utc)

    tdiff = tt - stop_timestamp
    if tdiff > dt.timedelta(seconds=0.5):
        tnt_stop = tt
        tnt_stop_ind = ti
        break

# for testing
#tnt_stop_ind = tnt_start_ind + 15 + 2
#tnt_start_ind = tnt_start_ind + 15

# narrow down lists
tnt_times_shift = tnttimes[tnt_start_ind:tnt_stop_ind]
dur_shift = durations[tnt_start_ind:tnt_stop_ind]
startf_shift = startf[tnt_start_ind:tnt_stop_ind]
stopf_shift = stopf[tnt_start_ind:tnt_stop_ind]


# period is 8 seconds!
# chunk the data into groups of 8 seconds!

burst = burst[0]
cfg = burst['config']

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
if len(burst['G']) > 1:
    start_timestamp = dt.datetime.utcfromtimestamp(burst['G'][0]['timestamp']) - dt.timedelta(seconds=float(cfg['SAMPLES_ON']/fs))
else:
    start_timestamp = dt.datetime.utcfromtimestamp(burst['header_timestamp']) - dt.timedelta(seconds=float(cfg['SAMPLES_ON']/fs))

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

# break up and average


intchecks = [1]
for intcheck in intchecks:
    fig, ax2 = plt.subplots(nrows=1, figsize=(10,6), sharex=True, sharey=True)
    E_FD = ax2
    allE = np.zeros(8*fs)
    for bi, br in enumerate(range(0,8*(intcheck+1)*fs,8*fs)):
        E_chunk = E_td_spaced[br:br+(8*fs)]

        if bi > 11:
            print('hey')
            an_array = np.empty(4*fs)
            an_array[:] = np.NaN
            E_chunk = list(E_chunk)
            [E_chunk.append(an) for an in an_array]
            E_chunk = np.array(E_chunk)

        allE = np.nansum(np.dstack((allE,E_chunk)),2)

    final_ae = [ae/(bi+1) for ae in allE[0][:fs*8]]
    if intcheck ==12:
        final_ae = [ae/(bi+1) for ae in allE[0][:fs*7]]

        final_ae2 = [ae/(bi) for ae in allE[0][fs*7:fs*8]]

        [final_ae.append(ae) for ae in final_ae2]
    E_td_spacedb = np.array(final_ae)
    # E spectrogram -- "spectrum" scaling -> V^2; "density" scaling -> V^2/Hz
    ff,tt, FE = scipy.signal.spectrogram(E_td_spacedb, fs=fs_equiv, window=window,
                nperseg=nfft, noverlap=nfft*overlap,mode='psd',scaling='density') # changed to density
    E_S_mag = 20*np.log10(np.sqrt(FE))
    E_S_mag[np.isinf(E_S_mag)] = -100

    # change time axis
    dt_axis = [start_timestamp + dt.timedelta(seconds=t) for t in tt]

    pe = E_FD.pcolormesh(dt_axis, ff/1000, E_S_mag, cmap = cm,  vmin=e_clims[0], vmax=e_clims[1])

    left, bottom, width, height = E_FD.get_position().bounds
    cax = fig.add_axes([left*7.28, bottom, width * 0.02, height])
    ce = fig.colorbar(pe, orientation='vertical', cax=cax)
    ce.set_label('dB[(uV/m)^2/Hz]')
    E_FD.set_ylabel('Frequency [kHz]')
    E_FD.set_xlabel('Time')

    ax = ax2

    for di, dur in enumerate(dur_shift[:4]):
        tim = tnt_times_shift[di]
        start = startf_shift[di]
        stop = stopf_shift[di]
        tim = tim  + dt.timedelta(seconds=6.8)

        ax.plot([tim, tim + dt.timedelta(microseconds=dur*1e3)],[start/1e3,stop/1e3],'white')
    ax2.set_ylim([7,10])

    labels = [item.get_text() for item in ax.get_xticklabels()]
    labels[:7] = [1,2,3,4,5,6,7]

    ax.set_xticklabels(labels)
    ax2.set_xlabel('seconds')
    ax2.set_title(str(intcheck)+' peroids')
    plt.savefig(str(intcheck)+'.png')
    plt.close()