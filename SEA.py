
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
import time
from inspect import getmembers, isclass
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

os.chdir(cwd)
import matplotlib.pyplot as plt 
import numpy as np
import datetime as dt 
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
import matplotlib as mpl 
from parula_colormap import parula
import scipy.signal

# script to run supoerposed epoch analysis on the burst data

burst_dir = '817'
onlyfiles = [f for f in listdir(burst_dir) if isfile(join(burst_dir, f))] 

fname_burst = onlyfiles[0]
# get burst data
burst = read_burst_XML(burst_dir + '/' + fname_burst)

# get TNT log data
fname_tnt = find_TNT(fname_burst,'TNT_logs/')
tnttimes, durations, startf, stopf = TNT_log(fname_tnt)

start_timestamp = dt.datetime(2020,8,17,21,20,30, tzinfo=dt.timezone.utc)

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

from matplotlib import gridspec


# plot it
import seaborn as sns
                           
fig = plt.figure(figsize=(12,6))

gs = gridspec.GridSpec(2, 2, wspace=0.18,hspace=0.28)

ax1 = plt.subplot(gs[0:2,1])
ax2 = plt.subplot(gs[0,0])
ax3 = plt.subplot(gs[1,0])

plt.tight_layout()

intchecks = [1]
for intcheck in intchecks:
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
    E_FD.set_ylabel('Frequency [kHz]',fontsize=12)
    #E_FD.set_xlabel('Time')

    ax = ax2

    ax2.set_ylim([7,10])
    ax2.set_xlim([dt_axis[0]-dt.timedelta(seconds=0.01),dt_axis[-1]+dt.timedelta(seconds=0.01)])

    labels = [0,1,2,3,4,5,6,7,8]
    ax.set_xticklabels(labels)
    #ax2.set_xlabel('seconds')
    ax2.set_title(str(intcheck)+' period',fontsize=12)
    #plt.show()
    #plt.savefig(str(intcheck)+'.png')
    #plt.close()

intchecks = [12]
for intcheck in intchecks:
    #fig, ax2 = plt.subplots(nrows=1, figsize=(10,6), sharex=True, sharey=True)
    E_FD = ax3
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

    E_FD.set_ylabel('Frequency [kHz]',fontsize=12)
    #E_FD.set_xlabel('Time')

    ax = ax3

    ax.set_ylim([7,10])
    ax.set_xlim([dt_axis[0]-dt.timedelta(seconds=0.01),dt_axis[-1]+dt.timedelta(seconds=0.01)])

    labels = [0,1,2,3,4,5,6,7,8]
    ax.set_xticklabels(labels)
    ax.set_xlabel('Seconds',fontsize=12)
    ax.set_title(str(intcheck)+' periods',fontsize=12)



tnt_times, durations, start_f, stop_f = TNT_log('/home/rileyannereid/workspace/dsxvpm_analysis/TNT_logs/TNT_HPT_Catalog_2020-08-17_210503-2020-08-17_213457.txt')
for ti,tt in enumerate(tnt_times):
    if tt - dt.datetime(2020,8,17,21,20,35) < dt.timedelta(seconds=0.1) and tt > dt.datetime(2020,8,17,21,20,34):
        tt_start = ti
    if dt.datetime(2020,8,17,21,22,5) - tt < dt.timedelta(seconds=0.1) and tt > dt.datetime(2020,8,17,21,22,0):
        tt_end = ti-1
        break

t_axis_tnt = tnt_times[tt_start:tt_end]

for tam, dur, stf, spf in zip(t_axis_tnt, durations[tt_start:tt_end], start_f[tt_start:tt_end], stop_f[tt_start:tt_end]):
    tam = tam - tnt_times[tt_start]
    tam = tam.total_seconds()
    ax1.plot([tam, tam+(dur/1e3)],[stf/1e3,spf/1e3],'r')

E_FD = ax1
# E spectrogram -- "spectrum" scaling -> V^2; "density" scaling -> V^2/Hz
ff,tt, FE = scipy.signal.spectrogram(E_td_spaced, fs=fs_equiv, window=window,
            nperseg=nfft, noverlap=nfft*overlap,mode='psd',scaling='density') # changed to density
E_S_mag = 20*np.log10(np.sqrt(FE))
E_S_mag[np.isinf(E_S_mag)] = -100

pe = ax1.pcolormesh(tt, ff/1000, E_S_mag, cmap = cm,  vmin=e_clims[0], vmax=e_clims[1])

left, bottom, width, height = E_FD.get_position().bounds
cax = fig.add_axes([left*7.28, bottom, width * 0.02, height])
ce = fig.colorbar(pe, orientation='vertical', cax=cax)
ce.set_label(r"$dB[(uV/m)^2/Hz]$",fontsize=12)
ax1.set_ylabel('Frequency [kHz]', fontsize=12)
ax1.set_xlim([0,90])
ax1.set_xlabel('Seconds',fontsize=12)

fig.suptitle('VPM Burst-Mode Data 08-17-2020 21:20:35 UTC',fontsize=18)
rasterize_list = [ax1,ax2, ax3]
rasterize_and_save('SEA_forpub.svg', rasterize_list, dpi=500)
plt.close()