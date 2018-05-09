import matplotlib.pylab as mpl
import matplotlib.patches as patches
from matplotlib import gridspec
import numpy as np
import pandas as pd
from scipy import stats

def plot_rasters(T_rasters, V_rasters, window, bin_size, modality = 'touch', ylim_r = None, ylim_p = None, size_mult = 1):
    
    fig = mpl.figure(figsize=(8*size_mult, 7*size_mult))
    first_raster = T_rasters[0]
    gs1 = gridspec.GridSpec(1,1)
    gs2 = gridspec.GridSpec(1,1)
    gs3 = gridspec.GridSpec(1,1)
    gs4 = gridspec.GridSpec(1,1)
    gs1.update(bottom = 0.88, top=0.95, left = 0.2, right = 0.83)
    gs3.update(bottom=0.15, top=0.41, left = 0.2, right = 0.83)
    gs2.update(bottom=0.45, top=0.88, left = 0.2, right = 0.83)
    gs4.update(bottom=0.45, top=0.88, left = 0.83, right = 0.9)
    
    ax1 = mpl.subplot(gs1[0, 0])
    ax2 = mpl.subplot(gs2[0, 0])
    ax3 = mpl.subplot(gs3[0, 0])
    patch_ax = mpl.subplot(gs4[0, 0], sharey = ax2)

    trial_type = 0
    trial_total = 1
    hists = []
    colors= ['C4', 'C7','C2']
    
    if modality == 'visual':
        rasters = V_rasters
        blocks =['Visual\nblock', 'Touch\nblock']
        ax1.add_patch(patches.Rectangle((0,0), 0.15, 1, facecolor = 'C1', alpha = 0.5))
        block_colors = ['C1', 'C0']
    else:
        rasters = T_rasters
        blocks =['Touch\nblock','Visual\nblock']
        ax1.add_patch(patches.Rectangle((0,0), 0.15, 1, facecolor = 'C0', alpha = 0.5))
        block_colors = ['C0', 'C1']
    block_lims = []
    for i in range(len(rasters)):
        ras = rasters[trial_type]
        spike_counts = []
        for trial, spike in enumerate(ras['spike_times(stim_aligned)']):
            spike = spike[(spike>window[0]) & (spike<=window[1])]
            ax2.vlines(spike, trial + trial_total - .5, trial + trial_total +.5)
            ax2.vlines(ras.iloc[trial]['first_lick'], trial + trial_total - .5, trial +
                       trial_total + .5, color = colors[i], linewidth = 5)
            
            spike = spike[(spike>window[0]) & (spike<=window[1])]
            edges = np.arange(window[0], window[1]+bin_size*2, bin_size)
            count, _ = np.histogram(spike,edges)
            spike_counts.append(count)

        patch_ax.add_patch(patches.Rectangle((0,trial_total -.5), 1,
                                       trial+1, facecolor = colors[trial_type], alpha = 0.5))
        if trial_type in [0,1]:
            patch_ax.plot([0,2],[trial_total-.5, trial_total-.5], '--k')
            block_lims.append(trial_total)
        elif trial_type == 3:
            patch_ax.plot([0,2],[trial_total+ trial, trial_total+trial], '--k')
            block_lims.append(trial_total+ trial+1)

        trial_total = trial_total + trial +1
        trial_type += 1
        average_hist = np.mean(spike_counts, axis=0)/bin_size
        SE_hist = stats.sem(spike_counts)/bin_size
        
        ax3.plot(edges[1:], average_hist, color = colors[i])
        ax3.fill_between(edges[1:], average_hist-SE_hist, average_hist+SE_hist, alpha = 0.5, color = colors[i])
    
    if modality == 'Touch':
        text_loc = 0.1
    else:
        text_loc = 0.1
    ax3.text(text_loc, .57, "CR", transform=ax3.transAxes, color = colors[0])
#     ax3.text(text_loc, .67, "False alarms", transform=ax3.transAxes, color = colors[1])
    ax3.text(text_loc, .70, "Misses", transform=ax3.transAxes, color = colors[1])
    ax3.text(text_loc, .84, "Hits", transform=ax3.transAxes, color = colors[2])

    for ax in [ax1,ax2,ax3, patch_ax]:
        ax.set_xlim(window[0],window[1]-bin_size)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        
    ax1.axis('off')
    ax1.set_ylim(0,2)

    ax2.spines['bottom'].set_visible(False)
    ax2.set_ylabel('Trials')    
    ax2.axes.get_xaxis().set_ticks([])
    ax2.set_ylim(-1, trial_total+.5)

    ax3.set_xlabel('Time(s)')
    ax3.set_ylabel('Firing\nrate (Hz)') 
    
    og_ylim = ax2.get_ylim()
    if ylim_r != None:
        ax2.set_ylim(ylim_r)
        ax2.spines['left'].set_bounds(0, trial_total)
        ax2.set_yticks(np.arange(0, og_ylim[1], 20))
    if ylim_p != None:
        ax3.set_ylim(ylim_p)

    patch_ax.axis('off')
    patch_ax.set_xlim(0,2)
    patch_ax.text(1.2, (block_lims[1] - block_lims[0])/2.5, blocks[1], color = block_colors[1])
    patch_ax.text(1.2,  (block_lims[1]/3) + block_lims[1], blocks[0],  color = block_colors[0])
    
    return fig
	
def plot_unit(log_df, uni_id, x_min, x_max, bin_size, ylim_p=None, modality = 'touch', size_mult = 1):

    ## need to subtract 1 bin from both ends because np.histogram counts spikes occuring from and including
    ## the bin label up to but not including the label for bin+1 - this would make it look like activity 
    ## starts before the stimulus onset.
    # x_min, xmax = x_min - bin_size, x_max-bin_size

    x_min, x_max = x_min-bin_size, x_max-bin_size

    current_cell = log_df[log_df['uni_id'] == uni_id]
    tac_stim = 'Stim_Som'
    vis_stim = 'Stim_Vis'


    cell_hit_v = current_cell[(current_cell['response'] != 0) &
                           (current_cell['trial_type'] == vis_stim) & (current_cell['correct'] == 1)]
    cell_miss_v = current_cell[(current_cell['trial_type'] == vis_stim) &
                          (current_cell['response'] == 0) & (current_cell['block_type'] == 'Visual')]
    cell_cr_v = current_cell[(current_cell['trial_type'] == vis_stim)&
                          (current_cell['response'] == 0) & (current_cell['block_type'] == 'Whisker')]


    cell_hit_t = current_cell[(current_cell['response'] != 0) &
                           (current_cell['trial_type'] == tac_stim) & (current_cell['correct'] == 1)]
    cell_miss_t = current_cell[(current_cell['trial_type'] == tac_stim)&
                          (current_cell['response'] == 0) & (current_cell['block_type'] == 'Whisker')]
    cell_cr_t = current_cell[(current_cell['trial_type'] == tac_stim)&
                          (current_cell['response'] == 0) & (current_cell['block_type'] == 'Visual')]
#     from IPython.core.debugger import Tracer; Tracer()() 
    t_rasters = [cell_cr_t, cell_miss_t, cell_hit_t]
    v_rasters = [cell_cr_v, cell_miss_v, cell_hit_v]
    max_ylim_r = current_cell['trial_type'].value_counts().max()

    fig = plot_rasters(t_rasters, v_rasters, [x_min, x_max], bin_size, ylim_p = ylim_p, ylim_r = (0,max_ylim_r), modality = modality, size_mult = size_mult)

    return fig

mpl.close('all')

mpl.close('all')

def plot_auc(uni_id, cp_df, t_long_cp, v_long_cp, t_short_cp, v_short_cp, window, bin_size, size_mult=1):

    ## need to subtract 1 bin from both ends because np.histogram counts spikes occuring from and including
    ## the bin label up to but not including the label for bin+1 - this would make it look like activity 
    ## starts before the stimulus onset.
    offset = np.abs(-1/bin_size)
    start = int(offset + (window[0]/bin_size))
    stop = int(offset + (window[1]/bin_size))


    edges = np.arange(window[0],window[1],bin_size)
    fig = mpl.figure(figsize=(18*size_mult, 4*size_mult))
    fig.subplots_adjust(bottom=0.2, left = 0.2, right = 0.83, wspace = 1.05)


    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)


    unit_ind = t_long_cp['uni_id'] == uni_id
    unit_row = cp_df[cp_df['uni_id'] == uni_id]

    t_long_cp_aucs = t_long_cp.loc[unit_ind, t_long_cp.columns.str.contains('auc')].as_matrix()[0]
    t_long_cp_upper = t_long_cp.loc[unit_ind, t_long_cp.columns.str.contains('up')].as_matrix()[0]
    t_long_cp_lower = t_long_cp.loc[unit_ind, t_long_cp.columns.str.contains('low')].as_matrix()[0]
    v_long_cp_aucs = v_long_cp.loc[unit_ind, v_long_cp.columns.str.contains('auc')].as_matrix()[0]
    v_long_cp_upper = v_long_cp.loc[unit_ind, v_long_cp.columns.str.contains('up')].as_matrix()[0]
    v_long_cp_lower = v_long_cp.loc[unit_ind, v_long_cp.columns.str.contains('low')].as_matrix()[0]

    t_short_cp_aucs = t_short_cp.loc[unit_ind, t_short_cp.columns.str.contains('auc')].as_matrix()[0]
    t_short_cp_upper = t_short_cp.loc[unit_ind, t_short_cp.columns.str.contains('up')].as_matrix()[0]
    t_short_cp_lower = t_short_cp.loc[unit_ind, t_short_cp.columns.str.contains('low')].as_matrix()[0]
    v_short_cp_aucs = v_short_cp.loc[unit_ind, v_short_cp.columns.str.contains('auc')].as_matrix()[0]
    v_short_cp_upper = v_short_cp.loc[unit_ind, v_short_cp.columns.str.contains('up')].as_matrix()[0]
    v_short_cp_lower = v_short_cp.loc[unit_ind, v_short_cp.columns.str.contains('low')].as_matrix()[0]

    t_late_onset = list(unit_row['touchL_cp_onset'])
    t_late_onset_1cyc = list(unit_row['touchS_cp_onset'])
    v_late_onset = list(unit_row['visL_cp_onset'])
    v_late_onset_1cyc = list(unit_row['visS_cp_onset'])

    ax1.add_patch(patches.Rectangle((0,0), 0.15, 1, facecolor = 'C7', alpha = 0.4))

    # from IPython.core.debugger import Tracer; Tracer()() 

    ax1.plot(edges, t_long_cp_aucs[start:stop], 'o-',color = 'C0', markersize = 5*size_mult)
    ax1.plot([t_late_onset]*2, [0,1], '--', color = 'C0')
    ax1.fill_between(edges, t_long_cp_lower[start:stop].astype(float),
                     t_long_cp_upper[start:stop].astype(float), alpha = 0.5,
                     color = 'C0')

    ax2.add_patch(patches.Rectangle((0,0), 0.15, 1, facecolor = 'C7', alpha = 0.4))
    ax2.plot(edges, v_long_cp_aucs[start:stop], 'o-', color = 'C1',  markersize = 5*size_mult)
    ax2.plot([v_late_onset]*2, [0,1],color = 'C1')
    ax2.fill_between(edges, v_long_cp_lower[start:stop].astype(float),
                     v_long_cp_upper[start:stop].astype(float), alpha = 0.5,
                     color = 'C1')

    ax1.plot(edges, t_short_cp_aucs[start:stop], 'o-',color = '#1d91c0', markersize = 5*size_mult)
    ax1.plot([t_late_onset_1cyc]*2, [0,1],color = 'c')
    ax1.fill_between(edges, t_short_cp_lower[start:stop].astype(float),
                     t_short_cp_upper[start:stop].astype(float), alpha = 0.5,
                     color = 'c')

    ax2.plot(edges, v_short_cp_aucs[start:stop], 'o-', color = '#fe9929', markersize = 5*size_mult)
    ax2.plot([v_late_onset_1cyc]*2, [0,1],color = '#feb24c')
    ax2.fill_between(edges, v_short_cp_lower[start:stop].astype(float),
                     v_short_cp_upper[start:stop].astype(float), alpha = 0.5,
                     color = '#feb24c')

    for ax in [ax1, ax2]:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        ax.set_xlim(window[0]-bin_size, window[1]-bin_size)
        ax.set_ylim([0,1])
        ax.set_xlabel('Time from stim onset(s)')
        ax.set_ylabel('AUC score')
        ax.plot(ax.get_xlim(), [0.5, 0.5], linestyle = '--', color = 'k' )
        ax.set_xticks([0,0.5])
    ax1.text(.3, .26, "Touch-lick/no lick",transform=ax1.transAxes, color = 'C0')
    ax1.text(.3, .15, "Short touch-lick/no lick",transform=ax1.transAxes, color = 'c')

    ax2.text(.4, .21, "Visual-lick/no lick",transform=ax2.transAxes, color = 'C1')
    ax2.text(.4, .1, "Short visual-lick/no lick",transform=ax2.transAxes, color = '#feb24c')
    return fig, ax1, ax2

