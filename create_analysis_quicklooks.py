import numpy as np 
import matplotlib.pyplot as plt
import datetime as dt 
import os
cwd = os.getcwd()
from read_files import find_TNT, TNT_log, read_burst_XML
from plots import plot_TNT, plot_burst, plot_shift

# make this a loop for bursts!
fname_burst = 'dsx_burst/VPM_burst_TD_2020-04-06_2204.xml'

# get burst data
burst = read_burst_XML(fname_burst)

# get TNT log data
fname_tnt = find_TNT(fname_burst,'TNT_logs/')
tnttimes, durations, startf, stopf = TNT_log(fname_tnt)
burst_samples = float(burst[0]['config']['SAMPLES_ON']) / float(burst[0]['config']['SAMPLES_OFF'])
burst_dur = float(burst[0]['config']['burst_pulses']) * 10 + (float(burst[0]['config']['burst_pulses']) * 10 / burst_samples) 

# ------- ------------------ plotting ------------------- -------
fig, ax2 = plt.subplots(nrows=1, sharex=True, sharey=True)

start_timestamp = plot_burst(fig, ax2, burst[0])
plot_TNT(ax2,tnttimes,durations,startf,stopf)

# some set up for plotting the shifted pulses
# find 5 seconds before start_timestamp and 5 seconds after -- minimize how much we have to raytrace
for ti, tt in enumerate(tnttimes):
    tt = tt.replace(tzinfo=dt.timezone.utc) # cant be timezone naive
    start_timestamp = start_timestamp.replace(tzinfo=dt.timezone.utc)

    tdiff = start_timestamp - tt
    if tdiff < dt.timedelta(seconds=5):
        tnt_start = tt
        tnt_start_ind = ti
        break

stop_timestamp = start_timestamp + dt.timedelta(seconds=burst_dur)
for ti, tt in enumerate(tnttimes):
    tt = tt.replace(tzinfo=dt.timezone.utc)
    start_timestamp = start_timestamp.replace(tzinfo=dt.timezone.utc)

    tdiff = tt - stop_timestamp
    if tdiff > dt.timedelta(seconds=5):
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

nrays = int(1e3)
rayfile_directory = '/home/rileyannereid/workspace/SR_output'

# plot the shifted pulses - runs the raytracer and calcualtes time
plot_shift(ax2, nrays, rayfile_directory, tnt_times_shift, dur_shift, startf_shift, stopf_shift)

# also plot 3D view of the vectors in doppler calculation !

# ------- --------------------- formatting ------------------------------ -----
os.chdir(cwd)
time_str = dt.datetime.strftime(start_timestamp, '%Y-%m-%d %H:%M:%S')
fig.suptitle('VPM Burst ' + time_str)
plt.savefig('VPM Burst' + time_str)
plt.show()

# add MLT to plots - or just get the plots from the burst dir