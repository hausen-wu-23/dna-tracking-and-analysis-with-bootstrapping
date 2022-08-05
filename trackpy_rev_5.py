from __future__ import division, unicode_literals, print_function

import matplotlib as mpl
import matplotlib.pyplot as plt

# tweaking style
mpl.rc('figure',  figsize=(10, 5))
mpl.rc('image', cmap='gray')

# import necessary libraries
import numpy as np
import pandas as pd

import pims
import trackpy as tp

# options for fitting
def func(x, a, b, c):
    return a*x**b + c

def func2(x, a, b):
    return a*x**b 

def func3(x, a, c):
    return a*x**alpha + c

# create a file list for all the tif files needed to be processed 
file_l = []

# 0618
for i in range(1,17):
    file_l.append('Y://Hausen W//Data//DNA Data//20220618//EXP2//FOS25C_488MT+KN__%i//FOS25C_488MT+KN__%i_MMStack_Default.ome.tif' % (i, i))

for i in range(2,12):
    file_l.append('Y://Hausen W//Data//DNA Data//20220618//FOS25C+KN_SPT//FOS25C_488MT+KN_SPT__%i_MMStack_Default.ome.tif' % (i))

for i in range(1,11):
    file_l.append('Y://Hausen W//Data//DNA Data//20220618//FOS25L+KN_SPT//FOS25L_488MT+KN_SPT__%i_MMStack_Default.ome.tif' % (i))
    
# 0619
for i in range(1,11):
    file_l.append('Y://Hausen W//Data//DNA Data//20220619//FOS25C-KN_SPT//%i//FOS25C-KN_SPT__%i_MMStack_Default.ome.tif' % (i, i))

# 0622
for i in range(1,11):
    file_l.append('Y://Hausen W//Data//DNA Data//20220622//FOS-25C_DDM_-KN//%i//FOS25C-KN_DDM__%i_MMStack_Default.ome.tif' % (i, i))
    
for i in range(1,11):
    file_l.append('Y://Hausen W//Data//DNA Data//20220622//FOS-25L_SPT_-KN//%i//FOS25L-KN_SPT__%i_MMStack_Default.ome.tif' % (i, i))

print(file_l)

# frame rate
fps = 10

# calibration of microscope
pixel_per_micron = 1/.256

# whether to invert the data for tracking, false for tracking bright spots, true for tracking dark spots
inv = False

# particle diameter and minmass, determined beforehand
par_d = 15
min_mass = 1200

# iterate through each file and run tracking
for file_n in file_l:
	# opening the file
	file_dir = file_n
	work_dir = '/'.join(file_dir.split('/')[:-1]) + '/'
	frames = pims.open(file_dir)

	# start tracking
	tp.quiet(suppress=False)
	f = tp.batch(frames[:], par_d, minmass=min_mass, invert=inv, processes=1)

	# link the paths
	t = tp.link_df(f, 7, memory=10)

	# filter out very short trajectories (less than 50 frames)
	t1 = tp.filter_stubs(t, 50)

	# Compare the number of particles in the unfiltered and filtered data.
	print('Before:', t['particle'].nunique())
	print('After:', t1['particle'].nunique())

	# filter by eccentricity
	t2 = t1[t1['ecc'] < .78]
    
	# ensemble msd
	em1 = tp.emsd(t2, 1/pixel_per_micron, fps)

	# create msd graph and save to the same folder as the tif file
	plt.loglog(em1.index, em1.values, 'o', markersize=8, alpha=0.2)
	plt.ylabel(r'$\langle \Delta r^2 \rangle$ [$\mu$m$^2$]',fontsize=14)
	plt.xlabel('lag time $(s)$', fontsize=14);
	plt.xticks(fontsize=14)
	plt.yticks(fontsize=14)
	plt.legend()
	plt.savefig(work_dir + "output%s.png" % (file_dir.split('/')[-1]))
	em1.to_csv(work_dir + 'emsd%s.csv'%(file_dir.split('/')[-1]), mode='w')
    
###################################################################################
#                                  bootstrapping                                  #
###################################################################################

	from scipy.optimize import curve_fit
	# progressbar is unnecessary, if you wish to use it to visualise progress, use `pip install progressbar2`
	import progressbar

	# list of effective diffusion coefficient and list of alpha
	lod = []
	loa = []

	# create progressbar widget, comment out if you don't want it
	widgets = [' [',
			 progressbar.Timer(format= 'elapsed time: %(elapsed)s'),
			 '] ',
			   progressbar.Bar('*'),' (',
			   progressbar.ETA(), ') ',
			  ]
	bar = progressbar.ProgressBar(max_value=39, 
								  widgets=widgets).start()

	# Run bootstrapping
	for i in range(39):
		# update progressbar
		bar.update(i)

		# select a random set of 80 trajectories, same trajectory can be picked multiple times, will be handled later
		rand = np.random.choice(t2['particle'].unique(), 80)

		# create dataframe according to the columns of original trajectories dataframe
		t3 = pd.DataFrame(columns=t2.columns)

		# save a list of already seen particles, for handling duplicate trajectories
		# trackpy will not process the same particle - we need to generate random seed to spoof trackpy into thinking it's a different particle
		alr_seen = []
		for i in rand:
			# find the trajectory of particle
			df_t = t2[t2['particle'] == i]

			# if not already seen, concatinate the trajectory
			if i not in alr_seen:
				t3 = pd.concat([t3, df_t])

				# add to alrady seen
				alr_seen.append(i)
			else:
				# generate random seed
				k = np.random.randint(1000, 3939)

				# run until the seed is unique
				while k in alr_seen:
					k = np.random.randint(1000, 3939)

				# assign new particle id
				df_t = df_t.assign(particle = k)

				# add to datagrame
				t3 = pd.concat([t3, df_t])

				# add to already seen
				alr_seen.append(k)
		
		# sort the values
		t3 = t3.sort_values(by='frame', ascending=True)

		# get the new msd
		em = tp.emsd(t3, 1/pixel_per_micron, fps)

		# fitting
		lim1 = 5
		lim2 = 65
		x = em.index[lim1:lim2]
		y = em.values[lim1:lim2]
		popt, pcov = curve_fit(func, x, y)
		lod.append(popt[0]/4)
		loa.append(popt[1])

	# save the fitting result to a file
	file1 = open('Y://Hausen W//Data//DNA Data//analysis.dat', "a")  # append mode
	file1.write("%s\n%.6f\n%.6f\n%.6f\n%.6f\n******\n" % (file_n, np.average(lod), np.average(loa),np.std(lod),np.std(loa)))
	file1.close()
	
	# print results, can be disabled
	print(file_n)
	print(lod)
	print(loa)

	print(np.average(lod))
	print(np.average(loa))

	print(np.std(lod))
	print(np.std(loa))
	print('********************************************')