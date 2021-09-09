import matplotlib.pyplot as plt 
import numpy as np
import datetime as dt 
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
import matplotlib as mpl 
from parula_colormap import parula
import scipy.signal
import matplotlib.gridspec as GS

# some helper functions for plotting VPM burst data etc. 

def plot_TNT(subplot, tnttimes, durations, startf, stopf): 
    ax = subplot
    for di, dur in enumerate(durations):
        tim = tnttimes[di]
        start = startf[di]
        stop = stopf[di]
        if dur == 150:
            ax.plot([tim, tim + dt.timedelta(microseconds=75e3)],[start/1e3,(start+100)/1e3],'white')
            ax.plot([tim + dt.timedelta(microseconds=75e3),tim + dt.timedelta(microseconds=150e3)],[(start+100)/1e3,(start-100)/1e3],'red')
            savef = start/1e3

    ax.tick_params(labelsize=7)
    if len(durations) < 1:
        savef = 0
    return savef 


def calc_LUHR(ray): # WE CAN OUTPUT STOP CONDITIONS
    EPS0 = 8.854187817e-12        # C^2/Nm^2

    # need info at last point (VPM alt)
    ind = ray['time'].index[-1]
    tt = ray['time']
    t = len(tt) - 1
    
    B   =  [ray['B0'].x[ind], ray['B0'].y[ind], ray['B0'].z[ind]]
    Bmag = np.linalg.norm(B)

    Q    = np.abs(np.array(ray['qs'].iloc[t,:]))
    M    = np.array(ray['ms'].iloc[t,:]) # elec, proton, He, O
    Ns   = np.array(ray['Ns'].iloc[t,:])

    Wcs   = Q*Bmag/M
    Wps2  = Ns*pow(Q,2)/EPS0/M
    Wps   = np.sqrt(Wps2)

    Mefflist = [ns_ion/m_ion for ns_ion, m_ion in zip(Ns[1:],M[1:])]
    Meffsum = sum(Mefflist)
    Meff_inv = Meffsum*M[0]/Ns[0]

    w_LHR2 = 1 / (1/(Wcs[0]*Wcs[-1]) + 1/(np.sqrt(Wps)) ) 
    w_LHR = np.sqrt(w_LHR2)
    f_LHR = w_LHR / (2*np.pi)
    print(f_LHR)
    
    # calc upper
    w_UHR2 = Wps2[0] + Wcs[0]**2
    w_UHR = np.sqrt(w_UHR2)
    f_UHR = w_UHR / (2*np.pi)

    # ask Bob -- LHR with several species
    return f_LHR, f_UHR


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
    ce.set_label('dB[(uV/m)^2/Hz]')
    E_FD.set_ylabel('Frequency [kHz]')
    E_FD.set_xlabel('Time')

    return start_timestamp

def plot_shift(subplot, nrays, rayfile_directory, tnt_times_shift, dur_shift, startf_shift, stopf_shift, plot_distribution=False, plot_HR=False):
    import os
    cwd = os.getcwd()
    from shifts import dopp_delay
    os.chdir(cwd)

    tnt_dop, tnt_t, alldop, allsec, allthetas, ray = dopp_delay(nrays, rayfile_directory, tnt_times_shift, dur_shift, startf_shift, stopf_shift)
    #print(allthetas)
    if plot_distribution==True:
        plot_dist(alldop, allsec, allthetas)
    
    ax = subplot

    if plot_HR == True:
        lhr, uhr = calc_LUHR(ray)
    else:
        lhr = 0
        
    for dopps, tdelays in zip(tnt_dop, tnt_t):
        for p_d, p_t in zip(dopps, tdelays):
            ax.scatter(p_t, p_d, color='red',s=2)

    ax.set_ylim([0,40])

    return lhr 

def is_outlier(points, thresh=3.5):

    if len(points.shape) == 1:
        points = points[:,None]
    median = np.median(points, axis=0)
    diff = np.sum((points - median)**2, axis=-1)
    diff = np.sqrt(diff)
    med_abs_deviation = np.median(diff)

    modified_z_score = 0.6745 * diff / med_abs_deviation

    return modified_z_score > thresh

def find_colors(bin_edges, hist_array, filtered_th, N, pltype):
    
    allbinthetas = []

    # find color avg for each bin
    for bi,bb in enumerate(bin_edges):
        bin_thetas = []
        current_bin = [bb, bin_edges[bi+1]]
        for di, dd in enumerate(hist_array):
            if current_bin[0] < dd and dd < current_bin[1]:
                bin_thetas.append(np.abs(filtered_th[di]))
        if len(bin_thetas) > 0:
            if pltype == 'dist':
                allbinthetas.append((sum(bin_thetas)/len(bin_thetas))-60)
            else:
                allbinthetas.append(np.std(np.array(bin_thetas)))
        else:
            allbinthetas.append(0)
        if bi == len(bin_edges) - 2:
            break

    my_cmap = plt.get_cmap("viridis",N)
    finalbinthetas = []
    for alls in allbinthetas:
        if alls > N-1:
            finalbinthetas.append(N-1)
        else:
            finalbinthetas.append(alls)

    th_colors =[my_cmap.colors[int(th)] for th in finalbinthetas]

    return th_colors, my_cmap

def find_bin_centers(bin_edges):
    bars_locs = []
    widths = []

    for bi,bb in enumerate(bin_edges):
        bars_locs.append((bb + bin_edges[bi+1]) /2)
        widths.append(0.95*(bb - bin_edges[bi+1]))
        if bi == len(bin_edges) - 2:
            break

    return bars_locs, widths

# -------------------------------------------------------------------
def plot_dist(tnt_dop, tnt_t, allthetas):
    pltypes = ['dist', 'std']
    for pltype in pltypes:

        fig, axs = plt.subplots(1, 3, gridspec_kw={'width_ratios':[1,1,0.1]}, tight_layout=True)
        n_bins=50

        # make arrays
        tnt_dop = np.array(tnt_dop)
        tnt_t = np.array(tnt_t)
        allthetas = np.array(allthetas)

        # filter 95%
        filtered_dop = tnt_dop[~is_outlier(tnt_dop)]
        filtered_t = tnt_t[~is_outlier(tnt_t)]

        filtered_thd = allthetas[~is_outlier(tnt_dop)]
        filtered_tht = allthetas[~is_outlier(tnt_t)]

        counts_d, dbin_edges = np.histogram(filtered_dop, bins=n_bins)
        counts_t, tbin_edges = np.histogram(filtered_t, bins=n_bins)

        # find corresponding color (by theta)
        if pltype=='dist':
            cn = 30
        else:
            cn = 15

        thd_colors, my_cmap = find_colors(dbin_edges, filtered_dop, filtered_thd, cn, pltype)
        tht_colors, my_cmap = find_colors(tbin_edges, filtered_t, filtered_tht, cn, pltype)

        # finally, get bin locations
        blocs_d, width_d = find_bin_centers(dbin_edges)
        blocs_t, width_t = find_bin_centers(tbin_edges)

        # add to subplot
        axs[0].bar(blocs_d, counts_d, align='center', width=width_d, color=thd_colors)
        axs[0].set_xlabel('doppler shift [Hz]')

        axs[1].bar(blocs_t, counts_t, align='center', width=width_t, color=tht_colors)
        axs[1].set_xlabel('seconds')

        # get color bar
        norm = mpl.colors.Normalize(vmin=0, vmax=cn)
        cb1 = mpl.colorbar.ColorbarBase(axs[2], cmap=my_cmap,
                                        norm=norm,
                                        orientation='vertical')
        if pltype=='dist':
            cb1.set_label('avg magnitude of initial polar wavenormal')
            cb1.ax.set_yticklabels([r'$60^{\circ}$',r'$65^{\circ}$',r'$70^{\circ}$',r'$75^{\circ}$',r'$80^{\circ}$',r'$85^{\circ}$',r'$90^{\circ}$'])
        else:
            cb1.set_label('std')
        
        axs[0].title.set_text('Doppler Shift')
        axs[1].title.set_text('Time Delay')
        
        fig.suptitle('8.834kHz Pulse at 2020-08-17 21:20:42')
        plt.savefig('dop_time_dist_20200817212042_'+pltype+'.png')
        plt.show()
        plt.close()

# -------------------------------------------------------------------
def plot_burst_spaced(fig, subplot, burst, section, chun, chunked=True):
    cfg = burst['config']
    E_FD = subplot

    system_delay_samps_FD = 200
    fs = 80000
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
    
    sections = []

    if cfg['SAMPLES_OFF'] == 0:
        E_td_spaced = E_coef*burst['E']
    else:
        # Insert nans into vector to account for "off" time sections
        E_td_spaced = []

        for k in np.arange(cfg['burst_pulses']):
            E_td_spaced.append(E_coef*burst['E'][k*cfg['SAMPLES_ON']:(k+1)*cfg['SAMPLES_ON']])
            sections.append(E_coef*burst['E'][k*cfg['SAMPLES_ON']:(k+1)*cfg['SAMPLES_ON']])
            E_td_spaced.append(np.ones(cfg['SAMPLES_OFF'])*np.nan)

        E_td_spaced = np.concatenate(E_td_spaced).ravel()
        
    if chunked==True:
        ax = E_FD
        nfft = int(2**12)
        overlap = 0.75
        window = 'hanning'
        try:
            # E spectrogram -- "spectrum" scaling -> V^2; "density" scaling -> V^2/Hz
            ff,tt, FE = scipy.signal.spectrogram(sections[section][int(chun*fs_equiv):int((chun+2)*fs_equiv)], fs=fs_equiv, window=window,
                    nperseg=nfft, noverlap=nfft*overlap,mode='psd', nfft=nfft, scaling='density') # changed to density
        except:
            print('short burst')
        E_S_mag = 20*np.log10(np.sqrt(FE))
        E_S_mag[np.isinf(E_S_mag)] = -100

        # change time axis
        section_timestamp = start_timestamp + dt.timedelta(seconds=section*20) + dt.timedelta(seconds=chun)
        dt_axis = [section_timestamp + dt.timedelta(seconds=t) for t in tt]
        
        try:
            pe = ax.pcolormesh(dt_axis, ff/1000, E_S_mag, cmap = cm,  vmin=e_clims[0], vmax=e_clims[1])
            left, bottom, width, height = ax.get_position().bounds
            cax = fig.add_axes([left*7.28, bottom, width * 0.02, height])
            ce = fig.colorbar(pe, orientation='vertical', cax=cax)
            ce.set_label('dB[(uV/m)^2/Hz]')
            ax.set_ylabel('Frequency [kHz]')
            ax.set_xlabel('Time')
        except:
            print('oops')


    else:
        nfft=2048
        overlap = 0.5
        window = 'hanning'
        # E spectrogram -- "spectrum" scaling -> V^2; "density" scaling -> V^2/Hz
        ff,tt, FE = scipy.signal.spectrogram(sections[section], fs=fs_equiv, window=window,
                    nperseg=nfft, noverlap=nfft*overlap,mode='psd',scaling='density') # changed to density
        E_S_mag = 20*np.log10(np.sqrt(FE))
        E_S_mag[np.isinf(E_S_mag)] = -100

        # change time axis
        section_timestamp = start_timestamp + dt.timedelta(seconds=section*20)
        dt_axis = [section_timestamp + dt.timedelta(seconds=t) for t in tt]

        pe = E_FD.pcolormesh(dt_axis, ff/1000, E_S_mag, cmap = cm,  vmin=e_clims[0], vmax=e_clims[1])

        left, bottom, width, height = E_FD.get_position().bounds
        cax = fig.add_axes([left*7.28, bottom, width * 0.02, height])
        ce = fig.colorbar(pe, orientation='vertical', cax=cax)
        ce.set_label('dB[(uV/m)^2/Hz]')
        E_FD.set_ylabel('Frequency [kHz]')
        E_FD.set_xlabel('Time')

    return section_timestamp

def plot_survey(S_data):
    # colormap -- parula is a clone of the Matlab colormap; also try plt.cm.jet or plt.cm.viridis
    cm = parula(); #plt.cm.viridis;

    # Sort by header timestamps
    S_data = sorted(S_data, key = lambda f: f['header_timestamp'])

    # Subset of data with GPS stamps included.
    # We need these for the line plots, regardless if we're using payload or bus timestamps.
    # Also confirm that we have at least one field from BESTPOS and BESTVEL messages,
    # since on rare occasions we miss one or the other.
    S_with_GPS = list(filter(lambda x: (('GPS' in x) and 
                                        ('timestamp' in x['GPS'][0]) and
                                        ('lat' in x['GPS'][0]) and
                                        ('horiz_speed' in x['GPS'][0])), S_data))
    S_with_GPS = sorted(S_with_GPS, key = lambda f: f['GPS'][0]['timestamp'])

    T_gps = np.array([x['GPS'][0]['timestamp'] for x in S_with_GPS])
    dts_gps = np.array([dt.datetime.fromtimestamp(x, tz=dt.timezone.utc) for x in T_gps])

    # Build arrays
    E = []
    B = []
    T = []
    F = np.arange(512)*40/512;
    
    # # Only plot survey data if we have GPS data to match
    bus_timestamps=True
    if bus_timestamps:
        # Sort using bus timestamp (finer resolution, but 
        # includes transmission error from payload to bus)
        for S in S_data:
            T.append(S['header_timestamp'])
            E.append(S['E_data'])
            B.append(S['B_data'])


            print(S)
            #print(gain)
    else:
        # Sort using payload GPS timestamp (rounded to nearest second.
        # Ugh, why didn't we just save a local microsecond counter... do that on CANVAS please)
        for S in S_with_GPS:
            T.append(S['GPS'][0]['timestamp'])
            E.append(S['E_data'])
            B.append(S['B_data'])
    T = np.array(T)

    dates = np.array([dt.datetime.utcfromtimestamp(t) for t in T])
    # -----------------------------------
    # Spectrograms
    # -----------------------------------
    E = np.array(E); B = np.array(B); T = np.array(T);
    fig, _ = plt.subplots()

    gs_data = GS.GridSpec(2, 2, width_ratios=[20, 1], wspace = 0.05, hspace = 0.05, figure = fig)
    ax1 = fig.add_subplot(gs_data[0,0])
    e_cbax = fig.add_subplot(gs_data[0,1])
    b_cbax = fig.add_subplot(gs_data[1,1])

    e_clims = [50,255] #[0,255] #[-80,-40]
    b_clims = [150,255] #[0,255] #[-80,-40]

    date_edges = np.insert(dates, 0, dates[0] - dt.timedelta(seconds=26))

    # Insert columns of NaNs wherever we have gaps in data (dt > 27 sec)
    per_sec = 26 # Might want to look this up for the shorter survey modes
    gaps = np.where(np.diff(date_edges) > dt.timedelta(seconds=(per_sec+2)))[0]

    d_gapped = np.insert(dates, gaps, dates[gaps] - dt.timedelta(seconds=per_sec + 3))
    E_gapped = np.insert(E.astype('float'), gaps - 1, np.nan*np.ones([1,512]), axis=0)

    # Plot E data
    p1 = ax1.pcolormesh(d_gapped,F,E_gapped.T, vmin=e_clims[0], vmax=e_clims[1], shading='flat', cmap = cm);
    cb1 = fig.colorbar(p1, cax = e_cbax)
    cb1.set_label(f'Raw value [{e_clims[0]}-{e_clims[1]}]')

    # # vertical lines at each edge (kinda nice, but messy for big plots)
    # g1 = ax1.vlines(dates, 0, 40, linewidth=0.2, alpha=0.5, color='w')
    # g2 = ax2.vlines(dates, 0, 40, linewidth=0.2, alpha=0.5, color='w')

    ax1.set_xticklabels([])
    ax1.set_ylim([0,40])

    fig.autofmt_xdate()
    ax1.set_xlabel("Time (H:M:S) on \n%s"%dt.datetime.utcfromtimestamp(T[0]).strftime("%Y-%m-%d"))
    # ax2.set_xlabel("Time (H:M:S)")

    ax1.set_ylabel('E channel\nFrequency [kHz]')