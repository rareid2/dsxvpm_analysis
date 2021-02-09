import numpy as np 
import matplotlib.pyplot as plt
import datetime as dt 

from read_files import TNT_log, read_burst_XML
from plots import plot_TNT, plot_burst

fname_tnt = 'TNT_logs/TNT_HPT_Catalog_2020-05-20_181503-2020-05-20_184457.txt'
tnttimes, durations, startf, stopf = TNT_log(fname_tnt)

fname_burst = 'dsx_burst/VPM_burst_TD_2020-05-20_1827.xml'
burst = read_burst_XML(fname_burst)

# ------- ------ plotting ------ -------
fig, (ax1, ax2) = plt.subplots(nrows=2, sharex=True, sharey=True)

plot_TNT(ax1, tnttimes, durations, startf, stopf)
start_timestamp = plot_burst(fig, ax2, burst[0])

# ------- ----- formatting ------- -----
time_str = dt.datetime.strftime(start_timestamp, '%Y-%m-%d %H:%M:%S')
fig.suptitle('VPM Burst ' + time_str)
plt.savefig('VPM Burst' + time_str)
plt.show()


# fix coordinates first!!
# then create a ray tracing plot (separate plot) -- maybe 3 snapshots? 
# make sure to add MLT
# add time delay! and doppler!