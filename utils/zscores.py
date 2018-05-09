import pandas as pd
import numpy as np
import scipy as sp
from tqdm import tnrange, tqdm_notebook


def label_trials(log_df, conds):
	whisker = log_df['block_type'] == 'Whisker'
	visual = log_df['block_type'] == 'Visual'
	correct = log_df['correct'] == 1
	lick = log_df['response'] != 0
	lick_right = log_df['response'] == 1
	lick_left = log_df['response'] == 2
	t_stim = log_df['trial_type'].str.contains('Som')
	v_stim = log_df['trial_type'].str.contains('Vis')
	short = log_df['trial_type'].str.contains('1Cyc')
	
	log_df['trial_label'] = None
	
	if np.any(['Hit' in cond for cond in conds]):
		log_df.loc[(t_stim  & whisker & correct), 'trial_label'] = 'Touch Stim Hit'
		log_df.loc[(t_stim & visual & correct), 'trial_label'] = 'Touch Stim CR'
		log_df.loc[(t_stim & lick & ~correct), 'trial_label'] = 'Touch Stim FA'
		log_df.loc[(t_stim & whisker & ~lick), 'trial_label'] = 'Touch Stim Miss'

		log_df.loc[(v_stim & visual & correct), 'trial_label'] = 'Visual Stim Hit'
		log_df.loc[(v_stim & lick & ~correct), 'trial_label'] = 'Visual Stim FA'
		log_df.loc[(v_stim & whisker & correct), 'trial_label'] = 'Visual Stim CR'
		log_df.loc[(v_stim & visual & ~lick), 'trial_label'] = 'Visual Stim Miss'
	elif np.any(['Left' in cond for cond in conds]):
		log_df.loc[(t_stim  & lick_right), 'trial_label'] = 'Touch Stim Lick Right'
		log_df.loc[(t_stim & lick_left), 'trial_label'] = 'Touch Stim Lick Left'
		log_df.loc[(t_stim & ~lick), 'trial_label'] = 'Touch Stim No Lick'
		# from IPython.core.debugger import Tracer; Tracer()()
		log_df.loc[(v_stim  & lick_right), 'trial_label'] = 'Visual Stim Lick Right'
		log_df.loc[(v_stim & lick_left), 'trial_label'] = 'Visual Stim Lick Left'
		log_df.loc[(v_stim & ~lick), 'trial_label'] = 'Visual Stim No Lick'
	elif np.any(['Lick' in cond for cond in conds]):
		log_df.loc[(t_stim  & ~short & lick), 'trial_label'] = 'Touch Stim Lick'
		log_df.loc[(t_stim & ~short & ~lick), 'trial_label'] = 'Touch Stim No Lick'
		log_df.loc[(v_stim & ~short & lick), 'trial_label'] = 'Visual Stim Lick'
		log_df.loc[(v_stim & ~short & ~lick), 'trial_label'] = 'Visual Stim No Lick'

		log_df.loc[(t_stim & short & lick), 'trial_label'] = 'Short Touch Stim Lick'
		log_df.loc[(t_stim & short & ~lick), 'trial_label'] = 'Short Touch Stim No Lick'
		log_df.loc[(v_stim & short & lick), 'trial_label'] = 'Short Visual Stim Lick'
		log_df.loc[(v_stim & short & ~lick), 'trial_label'] = 'Short Visual Stim No Lick'
	else:
		raise Exception('invalid conds')
	
	return log_df


def calculate_mean_FRs(log_df, unit_key_df, conds):

    def calc_unit_mean_FRs(uni_id):
        unit_rows = log_df[log_df['uni_id'] == uni_id]
		
        def calc_cond_FRs(cond):
            cond_rows = unit_rows.loc[unit_rows['trial_label'] == cond, 'spike_counts(stim_aligned)']
            if cond_rows.shape[0] < 1:
               bin_means = np.array([np.nan]*159)
            else:
               bin_means = np.mean(np.stack(cond_rows.as_matrix(), axis = 0), axis=0)       
               return bin_means
		
        mean_cond_FRs = pd.Series({c : calc_cond_FRs(c) for c in conds}) 
		
        return mean_cond_FRs
		
    mean_FRs = pd.DataFrame({uni_id : calc_unit_mean_FRs(uni_id) for uni_id in tqdm_notebook(unit_key_df['uni_id'])}).T
    mean_FRs.index.rename('uni_id', inplace = True)
    mean_FRs.reset_index(inplace=True)
    unit_key_df = unit_key_df.merge(mean_FRs, on = ['uni_id'])
	
    return unit_key_df

def calc_activation_resp(unit_key_df, bin_size, window, z_conds):
	act_resp_window = np.round(([0,.500] + np.abs(window[0]))/bin_size) -1

	##need to do this in two parts otherwise python crashes
	if np.any(['Hit' in cond for cond in z_conds]):
		act_resp_touch = unit_key_df.loc[:,'Touch Stim Hit(z_score)'].apply(lambda y: y[int(act_resp_window[0]):int(act_resp_window[1])])
		act_resp_vis = unit_key_df.loc[:,'Visual Stim Hit(z_score)'].apply(lambda y: y[int(act_resp_window[0]):int(act_resp_window[1])])
		
	elif np.any(['Left' in cond for cond in z_conds]):
		act_resp_touch = unit_key_df.loc[:,'Touch Stim Lick Right(z_score)'].apply(lambda y: y[int(act_resp_window[0]):int(act_resp_window[1])])		
		act_resp_vis = unit_key_df.loc[:,'Visual Stim Lick Right(z_score)'].apply(lambda y: y[int(act_resp_window[0]):int(act_resp_window[1])])

	elif np.any(['Lick' in cond for cond in z_conds]):
		act_resp_touch = unit_key_df.loc[:,'Touch Stim Lick(z_score)'].apply(lambda y: y[int(act_resp_window[0]):int(act_resp_window[1])])
		act_resp_vis = unit_key_df.loc[:,'Visual Stim Lick(z_score)'].apply(lambda y: y[int(act_resp_window[0]):int(act_resp_window[1])])
		
	act_resp_touch  = act_resp_touch.apply(lambda y: np.nanmean(y))
	act_resp_vis  = act_resp_vis.apply(lambda y: np.nanmean(y))

	unit_key_df['activation_resp_touch'] = act_resp_touch
	unit_key_df['activation_resp_vis'] = act_resp_vis
	
	return unit_key_df
	
def calc_z_values(unit_key_df, conds):
	z_scores_df = unit_key_df[conds].apply(lambda y: (y - unit_key_df['FR_mean'])/unit_key_df['FR_std'])
	
	for cond in conds:
		nans = z_scores_df[cond].isnull()
		if nans.any() :
			z_scores_df.loc[nans, cond] = [[np.nan]*159]
		
	return z_scores_df

	
def calc_z_scores(log_df, unit_key_df, bin_size, window, conds=None):
	
    if conds == None:
        conds = ['Touch Stim Hit', 'Touch Stim CR', 'Touch Stim Miss','Touch Stim FA',
				 'Visual Stim Hit','Visual Stim CR', 'Visual Stim Miss', 'Visual Stim FA']
		
        z_conds = ['Touch Stim Hit(z_score)', 'Touch Stim CR(z_score)', 'Touch Stim Miss(z_score)', 'Touch Stim FA(z_score)',
				   'Visual Stim Hit(z_score)','Visual Stim CR(z_score)','Visual Stim Miss(z_score)', 'Visual Stim FA(z_score)']
    else:
        z_conds = [cond+'(z_score)' for cond in conds]
	
    log_df = label_trials(log_df, conds)
    unit_key_df = calculate_mean_FRs(log_df, unit_key_df, conds)
    z_score_df = calc_z_values(unit_key_df, conds)
    z_score_df.columns = z_conds
    unit_key_df = pd.concat([unit_key_df,z_score_df], axis = 1)
    unit_key_df = calc_activation_resp(unit_key_df, bin_size, window, z_conds)
	
    return unit_key_df