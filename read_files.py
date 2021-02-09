import numpy as np 
import datetime as dt
import xml.etree.ElementTree as ET

def TNT_log(fname):
	f = open(fname, "r")
	datalines = [line for line in f]
	datalines = [line.strip() for line in datalines[12:]]
	datalines = [line.split() for line in datalines]
	
	tnt_times = []
	durations = []
	start_f = []
	stop_f = []

	for line in datalines:
		ymd = dt.datetime.strptime(line[0],'%Y-%m-%d')
		hms = dt.datetime.strptime(line[1],'%H:%M:%S.%f')
		tnt_time = ymd + dt.timedelta(hours=hms.hour,minutes=hms.minute,seconds=hms.second,microseconds=hms.microsecond)
		
		tnt_times.append(tnt_time)
		durations.append(float(line[2]))
		start_f.append(float(line[3]))
		stop_f.append(float(line[4]))

	return tnt_times, durations, start_f, stop_f


def read_burst_XML(filename):   
    ''' Reads burst elements from an xml file. '''

    # Open it
    with open(filename,'r') as f:
        tree = ET.parse(f)
    outs = []
    
    # Process all "burst" elements
    for S in tree.findall('burst'):
        d = dict()

        # Load the iso-formatted UTC string, and cast it to a Unix timestamp
        # (This is consistent with how we're storing times internally)
        header_timestamp_isoformat = S.attrib['header_timestamp']
        d['header_timestamp'] = dt.datetime.fromisoformat(header_timestamp_isoformat).replace(tzinfo=dt.timezone.utc).timestamp()
        
        if 'footer_timestamp' in S.attrib:
            footer_timestamp_isoformat = S.attrib['footer_timestamp']
            d['footer_timestamp'] = dt.datetime.fromisoformat(footer_timestamp_isoformat).replace(tzinfo=dt.timezone.utc).timestamp()    

        if 'experiment_number' in S.attrib:
            d['experiment_number'] = int(S.attrib['experiment_number'])

        if 'GAIN' in S.attrib:
            d['GAIN'] = str(S.attrib['GAIN'])     
        if 'CAL' in S.attrib:
            d['CAL'] = float(S.attrib['CAL']) 
        if 'FILT' in S.attrib:
            d['FILT'] = str(S.attrib['FILT'])      

        # Load burst configuration
        d['config'] = dict()
        for el in S.find('burst_config'):
            # print(cfg)
            # print(el.tag, el.text)
            # # for el in cfg:
                # print(el.name, el.text)   
            if el.tag in ['str', 'BINS']:
                d['config'][el.tag] = el.text
            else:
                d['config'][el.tag] = int(el.text)

        if S.find('bbr_config') is not None:
            # Load bbr configuration
            d['bbr_config'] = dict()
            for el in S.find('bbr_config'):
                d['bbr_config'][el.tag] = int(el.text)

        TD_FD_SELECT = d['config']['TD_FD_SELECT']

        # Load data fields -- as dtype "float" to preserve NaNs
        if TD_FD_SELECT == 1:
            # Time domain
            d['E'] = np.fromstring(S.find('E_data').text, dtype='float', sep=',')
            d['B'] = np.fromstring(S.find('B_data').text, dtype='float', sep=',')

        elif TD_FD_SELECT == 0:
            # Frequency domain
            ER = np.fromstring(S.find('E_data').find('real').text, dtype='float', sep=',')
            EI = np.fromstring(S.find('E_data').find('imag').text, dtype='float', sep=',')
            d['E'] = ER + 1j*EI
            
            BR = np.fromstring(S.find('B_data').find('real').text, dtype='float', sep=',')
            BI = np.fromstring(S.find('B_data').find('imag').text, dtype='float', sep=',')
            d['B'] = BR + 1j*BI

        # Load GPS data
        d['G'] = []
        for g in S.find('GPS'):
            tmp_dict = dict()
            for el in g:
                try:
                    tmp_dict[el.tag] = int(el.text)
                except:
                    tmp_dict[el.tag] = float(el.text)
            d['G'].append(tmp_dict)
        outs.append(d)

    # Return a list of dicts
    return outs