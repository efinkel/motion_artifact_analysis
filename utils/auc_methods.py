import scipy as sp
import scipy.io
import numpy as np
import pandas as pd
from tqdm import tnrange, tqdm_notebook
import matplotlib.pylab as mpl
from scipy import stats

def find_sig_bins(aucs):

	auc_col_names = aucs.columns.str.contains('auc')
	low_conf_col_names = aucs.columns.str.contains('low')
	up_conf_col_names = aucs.columns.str.contains('up')
	
	raw_sig_AUC = pd.DataFrame(~((aucs.loc[:,low_conf_col_names].fillna(0.5) >= 0.5).as_matrix() & 
						(aucs.loc[:,up_conf_col_names].fillna(0.5) <= 0.5).as_matrix())*1, index = aucs['uni_id'])

	auc_dir = (aucs.loc[:,auc_col_names] >= 0.5)*1 + (aucs.loc[:,auc_col_names] < 0.5)*-1
	auc_dir.index = aucs['uni_id']; auc_dir.columns = raw_sig_AUC.columns
	raw_sig_AUC = raw_sig_AUC*auc_dir
	smooth_sig_AUC = raw_sig_AUC.apply(lambda y: np.convolve(y, [1,1,1], 'same'), axis = 1)
	
	return smooth_sig_AUC


def get_cp_groups(auc_array, unit_key_df, bin_size):

	smooth_tac_sig_AUC = find_sig_bins(auc_array[0])
	smooth_vis_sig_AUC = find_sig_bins(auc_array[1])

	### to get the first bin of 3 that are significantly different from 0 need to subtract 1 from all_sig_tac/vis since convolve 
	### was used with 'same' setting (otherwise would need to subtract 2 if setting was 'full')
	all_sig_tac = (smooth_tac_sig_AUC.isin([3,-3])*1).apply(lambda y: np.where(y), axis = 1)
	all_sig_tac = all_sig_tac.apply(lambda y: (y[0]-1)*bin_size-1 if len(y[0])>0 else [])

	all_sig_vis = (smooth_vis_sig_AUC.isin([3,-3])*1).apply(lambda y: np.where(y), axis = 1)
	all_sig_vis = all_sig_vis.apply(lambda y: (y[0]-1)*bin_size-1 if len(y[0])>0 else [])

	all_first_sig_tac = all_sig_tac.apply(lambda y:  y[(y>=0) & (y<1)] if len(y)!=0 else [])
	all_first_sig_tac = all_first_sig_tac.apply(lambda y:  np.nan if len(y)==0 else np.min(y))

	all_first_sig_vis = all_sig_vis.apply(lambda y:  y[(y>=0) & (y<1)] if len(y)!=0 else [])
	all_first_sig_vis = all_first_sig_vis.apply(lambda y:  np.nan if len(y)==0 else np.min(y))

	sig_right_lick_only = ~np.isnan(all_first_sig_tac) & np.isnan(all_first_sig_vis)
	sig_left_lick_only = np.isnan(all_first_sig_tac) & ~np.isnan(all_first_sig_vis)
	sig_bidirec_lick = ~np.isnan(all_first_sig_tac) & ~np.isnan(all_first_sig_vis)

	sig_right_lick_only = list(sig_right_lick_only[sig_right_lick_only].index)
	sig_left_lick_only = list(sig_left_lick_only[sig_left_lick_only].index)
	sig_bidirec_lick = list(sig_bidirec_lick[sig_bidirec_lick].index)

	lick_group_inds = [sig_right_lick_only,sig_left_lick_only,sig_bidirec_lick]
	lick_group = []
	for group in lick_group_inds:
		late_onsets = pd.concat([all_first_sig_tac[group], all_first_sig_vis[group]], axis = 1).reset_index()
		late_onsets.columns = [['uni_id', 'Tac_late_onset', 'Vis_late_onset']]
		
		unit_inds = unit_key_df['uni_id'].isin(group)
		first_sig = pd.merge(unit_key_df.loc[unit_inds, ['uni_id','mouse_name', 'date', 'cluster_name',
														  'RT_median_TLR','RT_median_VLL','RT_mean_TLR','RT_mean_VLL',
														  'RT_std_TLR','RT_std_VLL', 'RT_num_TLR', 'RT_num_VLL', 
														 'activation_resp_touch', 'activation_resp_vis']], late_onsets, on = ['uni_id'])

		
		lick_group.append(first_sig.reset_index(drop=True))

	touch_lick_units, vis_lick_units, bimodal_lick_units = lick_group[0], lick_group[1], lick_group[2]
	
	return touch_lick_units, vis_lick_units, bimodal_lick_units
	