import matplotlib.pyplot as plt 
import numpy as np
import datetime as dt 
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
import matplotlib as mpl 
from parula_colormap import parula
import scipy.signal

def plot_TNT(subplot, tnttimes, durations, startf, stopf): 
    ax = subplot
    for di, dur in enumerate(durations):
        tim = tnttimes[di]
        start = startf[di]
        stop = stopf[di]
        if dur == 150:
            ax.plot([tim, tim + dt.timedelta(microseconds=75e3)],[start/1e3,(start+100)/1e3],'red')
            ax.plot([tim + dt.timedelta(microseconds=75e3),tim + dt.timedelta(microseconds=150e3)],[(start+100)/1e3,(start-100)/1e3],'red')
            savef = start/1e3
        elif dur == 250:
            savef = start/1e3
            ax.plot([tim, tim + dt.timedelta(microseconds=dur*1e3)],[start/1e3,stop/1e3],'red')
        else:
            savef = 0
            ax.plot([tim, tim + dt.timedelta(microseconds=dur*1e3)],[start/1e3,stop/1e3],'red')

    ax.set_ylabel('Frequency [kHz]')
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

    Mefflist = [ns_ion/m_ion for ns_ion, m_ion in zip(Ns[1:],M[1:])]
    Meffsum = sum(Mefflist)
    Meff_inv = Meffsum*M[0]/Ns[0]

    w_LHR2 = 1 / (1/(Wcs[0]*Wcs[-1]) + 1/(np.sqrt(Wps2[0])*np.sqrt(Wps2[-1])) ) 
    w_LHR = np.sqrt(w_LHR2)
    f_LHR = w_LHR / (2*np.pi)
    
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
    
    if plot_distribution==True:
        plot_dist(alldop, allsec, allthetas)
    
    ax = subplot

    if plot_HR == True:
        lhr, uhr = calc_LUHR(ray)
    else:
        lhr = 0
        
    for dopps, tdelays in zip(tnt_dop, tnt_t):
        for p_d, p_t in zip(dopps, tdelays):
            ax.scatter(p_t, p_d, color='lightsalmon',s=2)

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

def find_colors(bin_edges, hist_array, filtered_th, N):
    allbinthetas = []
    allstdthetas = []

    # find color avg for each bin
    for bi,bb in enumerate(bin_edges):
        bin_thetas = []
        current_bin = [bb, bin_edges[bi+1]]
        for di, dd in enumerate(hist_array):
            if current_bin[0] < dd and dd < current_bin[1]:
                bin_thetas.append(np.abs(filtered_th[di]))
        if len(bin_thetas) > 0:
            allbinthetas.append(sum(bin_thetas)/len(bin_thetas))
            allstdthetas.append(np.std(np.array(bin_thetas)))
        else:
            allstdthetas.append(0)
            allbinthetas.append(0)
        if bi == len(bin_edges) - 2:
            break

    my_cmap = plt.get_cmap("viridis",N)
    finalstdthetas = []
    for alls in allstdthetas:
        if alls > N-1:
            finalstdthetas.append(N-1)
        else:
            finalstdthetas.append(alls)

    th_colors =[my_cmap.colors[int(th)] for th in finalstdthetas]

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
    cn=15
    thd_colors, my_cmap = find_colors(dbin_edges, filtered_dop, filtered_thd, cn)
    tht_colors, my_cmap = find_colors(tbin_edges, filtered_t, filtered_tht, cn)

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
    cb1.set_label('avg magnitude of initial polar wavenormal')
    cb1.set_label('std')
    
    #cb1.ax.set_yticklabels([r'$0^{\circ}$',r'$10^{\circ}$',r'$20^{\circ}$',r'$30^{\circ}$',r'$40^{\circ}$',r'$50^{\circ}$',r'$60^{\circ}$', r'$70^{\circ}$', r'$80^{\circ}$',r'$90^{\circ}$'])
    
    axs[0].title.set_text('Doppler Shift')
    axs[1].title.set_text('Time Delay')
    
    fig.suptitle('8.834kHz Pulse at 2020-08-17 21:20:42')
    #plt.savefig('dop_time_dist_20200817212042.png')
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
        # E spectrogram -- "spectrum" scaling -> V^2; "density" scaling -> V^2/Hz
        ff,tt, FE = scipy.signal.spectrogram(sections[section][int(chun*fs_equiv):int((chun+2)*fs_equiv)], fs=fs_equiv, window=window,
                    nperseg=nfft, noverlap=nfft*overlap,mode='psd', nfft=nfft, scaling='density') # changed to density
        E_S_mag = 20*np.log10(np.sqrt(FE))
        E_S_mag[np.isinf(E_S_mag)] = -100

        # change time axis
        section_timestamp = start_timestamp + dt.timedelta(seconds=section*20) + dt.timedelta(seconds=chun)
        dt_axis = [section_timestamp + dt.timedelta(seconds=t) for t in tt]
        
        pe = ax.pcolormesh(dt_axis, ff/1000, E_S_mag, cmap = cm,  vmin=e_clims[0], vmax=e_clims[1])

        left, bottom, width, height = ax.get_position().bounds
        cax = fig.add_axes([left*7.28, bottom, width * 0.02, height])
        ce = fig.colorbar(pe, orientation='vertical', cax=cax)
        ce.set_label('dB[(uV/m)^2/Hz]')
        ax.set_ylabel('Frequency [kHz]')
        ax.set_xlabel('Time')
        #ax.set_ylim([7,10])

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
