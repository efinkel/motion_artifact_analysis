
import numpy as np
import pandas as pd
import matplotlib.pylab as mpl
import matplotlib.patches as patches
from matplotlib import gridspec
import scipy

def plot_zscore_psth(data, xlims, bin_size, colors, labels = None):
	mpl.close('all')
	
	num_subplots = len(list(data.keys()))
	fig = mpl.figure(figsize = (16,4))

	gs1 = gridspec.GridSpec(1,num_subplots)
	gs2 = gridspec.GridSpec(1,num_subplots)

	gs1.update(bottom=0.25, top=0.8, left = 0.15, right = 0.9, wspace=0.05)
	gs2.update(bottom = 0.81, top=0.87, left = 0.15, right = 0.9, wspace=0.05)
	
	axs = [mpl.subplot(gs1[0, 0])]
	for num in range(1,num_subplots):
		axs.append(mpl.subplot(gs1[0, num],sharey=axs[0]))

	xvals = np.arange(xlims[0],xlims[1],bin_size)
	## -1 is the first timepoint for all trials.
	## need to subtract 1 bin from both ends because np.histogram counts spikes occuring from and including
	## the bin label up to but not including the label for bin+1 - this would make it look like activity 
	## starts before the stimulus onset.
	xrange = np.array([abs(-1 - xvals[0]), abs(1+xvals[-1])])
	xrange = (xrange-bin_size)/bin_size
	inds = np.arange(int(round(xrange[0])),int(round(xrange[1])),1)
	
	stim_color = ['C0', 'C1']
	for i, col in enumerate(list(data.keys())):
		trial_types = data[col]
		ax = mpl.subplot(gs2[0, i],sharex=axs[i])
		ax.add_patch(patches.Rectangle((0,0), 0.15, 2, facecolor = stim_color[i], alpha = 0.5))
		ax.set_ylim(0,2)
		ax.axis('off')

		for tt_num, tt in enumerate(trial_types):
			z_map = tt[inds]
			m_psth = z_map.mean(axis = 0)
			
			sem_psth = scipy.stats.sem(z_map, axis = 0)
			axs[i].plot(xvals[:-1],m_psth, color = colors[tt_num])
			axs[i].fill_between(xvals[0:-1], m_psth-sem_psth, m_psth+sem_psth,alpha=0.5,color = colors[tt_num])
			axs[i].set_xticks(np.arange(xlims[0],xlims[1],abs(xlims[0])))
			axs[i].set_xlim(xvals[0],xvals[-1])
			axs[i].spines['right'].set_visible(False)
			axs[i].spines['top'].set_visible(False)
			axs[i].set_xlabel('Time from\nstim onset (s)')
			if i>0:
				axs[i].axes.get_yaxis().set_visible(False)
				axs[i].spines['left'].set_visible(False)
			
	axs[0].set_ylabel('Mean Z-score')
	
	for j,label in enumerate(labels):
		axs[-1].text(.35, .65+.1*j, label , transform=axs[i].transAxes, color = colors[0])
	ylim_hit_psth = axs[0].get_ylim()

	return fig