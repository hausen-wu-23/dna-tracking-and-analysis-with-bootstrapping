from cProfile import label
from re import X
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit

import os

# Import all the folders that include the csv files
path ="Y:\Hausen W\Data\DNA Data\MSD\C+KN"
path1 ="Y:\Hausen W\Data\DNA Data\MSD\C-KN"
path2 ="Y:\Hausen W\Data\DNA Data\MSD\L+KN"
path3 ="Y:\Hausen W\Data\DNA Data\MSD\L-KN"
path4 ="Y:\Hausen W\Data\DNA Data\MSD\C-KN_SPT"

# create list to save the alpha, effective coefficient and the labels
a_l = []
d_l = []
label_l = []

# for fitting
def func(x, a, b):
    return a*x**b


def make_graph(file_p, label, color):
	# Store the path of all the files found
	filelist = []

	for root, dirs, files in os.walk(file_p):
		for file in files:
			#append the file name to the list
			filelist.append(os.path.join(root,file))

	# read the first csv in the file list and create a dataframe in the same format
	df = pd.read_csv(filelist[0])
	df = pd.DataFrame(columns = df.columns)

	# concatinate all of the csv into one single dataframe
	for i in filelist:
		df = pd.concat([df, pd.read_csv(i)])
	
	# calculate msd/time
	df['msd'] = df['msd']/df['lagt']

	# turn lagt into integer so pandas can group them easily*
	# *floats cannot be directly compared
	df['lagt'] = df['lagt'] * 10

	# sort lag time (unnecessary but it's nice)
	df = df.sort_values('lagt')

	# explicitly declare types to prevent rounding errors
	df['msd'] = df['msd'].astype(float)
	df['lagt'] = df['lagt'].astype(int)

	# group all the same lag times together to average/std all the msd on the same lag time
	grouped = df.groupby('lagt', as_index=False).msd.mean()
	error = df.groupby('lagt', as_index=False).msd.std()['msd']

	# divide by 10 to get orginal lagtime
	grouped['lagt'] = grouped['lagt']/10

	# calculate standard error
	error = error / np.sqrt(len(filelist))
	
	x = grouped['lagt']
	y = grouped['msd']
	
	# uncomment following to graph the datapoints with errorbar
	# plt.errorbar(x, y, yerr = error, fmt='%so' % color, alpha=.2,ms=5, label=label)
	
	# fitting lagt to msd
	popt, pcov = curve_fit(func, x[:10], y[:10], maxfev=12000)
	x_lin = np.linspace(x[1], x[9], 500)

	# uncomment following to enable fit line of first half
	# plt.plot(x_lin, func(x_lin, *popt), '%s-' % color, label='%s: $y = %.5f * x^{%.5f}$' % (label, popt[0], popt[1]))
	a_l.append(popt[1])
	d_l.append(popt[0])
	label_l.append(label)
	
	popt2, pcov2 = curve_fit(func, x[10:], y[10:], maxfev=12000)
	x_lin2 = np.linspace(x[10], x[99], 500)

	# uncomment following to enable fit line of second half
	# plt.plot(x_lin2, func(x_lin2, *popt2), '%s-' % color, label='%s: $y = %.5f * x^{%.5f}$' % (label, popt2[0], popt2[1]))
	a_l.append(popt2[1])
	d_l.append(popt2[0])
	label_l.append('%s-2' % label)
	
	

ax = plt.axes()

# uncomment following two lines to use loglog scale
#ax.set_xscale("log")
#ax.set_yscale("log")
make_graph(path, 'C+KN', 'y')
make_graph(path1, 'C-KN DDM', 'b')
make_graph(path4, 'C-KN SPT', 'm')
make_graph(path2, 'L+KN', 'k')		
make_graph(path3, 'L-KN', 'r')	

# comment out following line to disable bar graph
plt.bar(label_l, d_l)

plt.legend()
plt.show()
