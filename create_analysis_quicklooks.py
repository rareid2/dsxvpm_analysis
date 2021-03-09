import numpy as np 
import matplotlib.pyplot as plt
import datetime as dt 
import os
from os import listdir
from os.path import isfile, join
cwd = os.getcwd()
from read_files import find_TNT, TNT_log, read_burst_XML
from plots import plot_TNT, plot_burst, plot_shift, plot_dist, plot_burst_spaced

burst_dir = 'focus'
onlyfiles = [f for f in listdir(burst_dir) if isfile(join(burst_dir, f))] 
for fi, fname_burst in enumerate(onlyfiles):

    # get burst data
    burst = read_burst_XML(burst_dir + '/' + fname_burst)

    # get TNT log data
    fname_tnt = find_TNT(fname_burst,'TNT_logs/')
    #if fname_tnt == 0:
    #    continue
    #tnttimes, durations, startf, stopf = TNT_log(fname_tnt)
    burst_samples = float(burst[0]['config']['SAMPLES_ON']) / float(burst[0]['config']['SAMPLES_OFF'])
    burst_dur = float(burst[0]['config']['burst_pulses']) * 10 + (float(burst[0]['config']['burst_pulses']) * 10 / burst_samples) 

    # ------- ------------------ plotting ------------------- -------
    output = 'quicklooks/'+fname_burst[13:28]
    try:  
        os.mkdir(output)  
    except OSError as error:  
        print(error)
    for section in range(0,5):
        for chun in range(0,10,2):
            fig, ax2 = plt.subplots(nrows=1, figsize=(6,6), sharex=True, sharey=True)

            start_timestamp = plot_burst_spaced(fig, ax2, burst[0], section, chun)
            
            # some set up for plotting the shifted pulses
            # find 5 seconds before start_timestamp and 5 seconds after -- minimize how much we have to raytrace
            """
            for ti, tt in enumerate(tnttimes):
                tt = tt.replace(tzinfo=dt.timezone.utc) # cant be timezone naive
                start_timestamp = start_timestamp.replace(tzinfo=dt.timezone.utc)

                tdiff = start_timestamp - tt
                if tdiff < dt.timedelta(seconds=1):
                    tnt_start = tt
                    tnt_start_ind = ti
                    break

            #stop_timestamp = start_timestamp + dt.timedelta(seconds=burst_dur)
            stop_timestamp = start_timestamp + dt.timedelta(seconds=10)
            
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

            nrays = int(10)
            rayfile_directory = '/home/rileyannereid/workspace/SR_output'
            
            savef = plot_TNT(ax2,tnt_times_shift,dur_shift,startf_shift,stopf_shift)

            # plot the shifted pulses - runs the raytracer and calcualtes time
            #lhr = plot_shift(ax2, nrays, rayfile_directory, tnt_times_shift, dur_shift, startf_shift, stopf_shift, 
            #    plot_distribution=False, plot_HR=True)

            #ax2.plot([start_timestamp + dt.timedelta(seconds=i*0.5) for i in range(30)], np.ones(30)*lhr/1e3, 'w--', )
            if savef !=0:
                ax2.set_ylim([savef - 2, savef + 2])
            ax2.set_ylim([2,5])
            """
            # ------- --------------------- formatting ------------------------------ -----
            os.chdir(cwd)
            time_str = dt.datetime.strftime(start_timestamp, '%Y-%m-%d %H:%M:%S')
            
            fig.suptitle('VPM Burst ' + time_str)
            plt.savefig(output+'/VPMBurst' + time_str+'_section'+str(section))
            #plt.show()
            plt.close()