import numpy as np 
import matplotlib.pyplot as plt
import datetime as dt 
import os
from os import listdir
from os.path import isfile, join
cwd = os.getcwd()
from read_files import find_TNT, TNT_log, read_burst_XML, read_survey_XML
from plots import plot_TNT, plot_burst, plot_shift, plot_dist, plot_burst_spaced, plot_survey
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

sdata = read_survey_XML('817/VPM_survey_data_2020-08-17.xml')


plot_survey(sdata)
plt.show()
plt.close()