function [fig, rms_summary] = plotPreStimMove(mouse_name, filtSignal, filtPerch, time, a, range)

badIndexes = abs(filtSignal) > 0.2;
filtSignal(badIndexes) = NaN;

t_hits = a.intanTrials_SsSb(cellfun(@(x) x.response == 1, a.BcontTrials_SsSb));
t_CR = a.intanTrials_SsVb(cellfun(@(x) x.response == 0, a.BcontTrials_SsVb));
t_misses = a.intanTrials_SsSb(cellfun(@(x) x.response == 0, a.BcontTrials_SsSb));

ts_hits = a.intanTrials_SscSb(cellfun(@(x) x.response == 1, a.BcontTrials_SscSb));
ts_CR = a.intanTrials_SscVb(cellfun(@(x) x.response == 0, a.BcontTrials_SscVb));
ts_misses = a.intanTrials_SscSb(cellfun(@(x) x.response == 0, a.BcontTrials_SscSb));

v_hits = a.intanTrials_VsVb(cellfun(@(x) x.response == 2, a.BcontTrials_VsVb));
v_CR = a.intanTrials_VsSb(cellfun(@(x) x.response == 0, a.BcontTrials_VsSb));
v_misses = a.intanTrials_VsVb(cellfun(@(x) x.response == 0, a.BcontTrials_VsVb));

vs_hits = a.intanTrials_VscVb(cellfun(@(x) x.response == 2, a.BcontTrials_VscVb));
vs_CR = a.intanTrials_VscSb(cellfun(@(x) x.response == 0, a.BcontTrials_VscSb));
vs_misses = a.intanTrials_VscVb(cellfun(@(x) x.response == 0, a.BcontTrials_VscVb));



fig = figure; hold on
set(fig, 'Visible', 'off');

touch_trials = {t_hits,ts_hits, t_CR, ts_CR, t_misses, ts_misses};
visual_trials = {v_hits, vs_hits, v_CR,vs_CR, v_misses,vs_misses};
tt = {touch_trials, visual_trials};
hit_color = [76, 167, 51]/256;
cr_color = [145, 104, 191]/256;
miss_color = [129, 129, 129]/256;
line_colors = {hit_color, hit_color, cr_color, cr_color, miss_color, miss_color};
stim_colors = {[51, 119, 182]/256, [246, 130, 0]/256};

pre_stim_rms_all = [];
rms_labels_all = [];
mod_label = [];
subplots = [];
line_heights = [];
for m = 1:2
    sp = subaxis(1,3,m, 'Margintop', 0.06, 'Marginbottom', 0.06); hold on;
    subplots = [subplots, sp];
    trial_types = tt{m};
    stim_col = stim_colors{m};
    count = 0;
    pre_stim_rms_mod = [];
    rms_labels_mod = [];
    for i=1:numel(trial_types)
        trial_starts = cellfun(@(x) x.rawTime(1), trial_types{i});
        stim_onsets = cellfun(@(x) x.stimOnsetTime, trial_types{i});
        aligned_so = stim_onsets+trial_starts;

        trial_start_inds = find(ismember(time,trial_starts));
        stim_onset_inds = stim_onsets*30000 + trial_start_inds;

        window_starts = stim_onset_inds + range(1)*30000;
        window_ends = stim_onset_inds + range(2)*30000;
        
        windows = arrayfun(@(x,y) {filtSignal(x:y), filtPerch(x:y)}, window_starts, window_ends, 'uni', 0);
        pre_stim = arrayfun(@(x,y) {filtSignal(x:y), filtPerch(x:y)}, window_starts, stim_onset_inds, 'uni', 0);
        post_stim = arrayfun(@(x,y) {filtSignal(x:y), filtPerch(x:y)}, stim_onset_inds, window_ends, 'uni', 0);
        
        
        pre_stim_rms = cellfun(@(x) x{1}(~isnan(x{1})), pre_stim, 'uni', 0);
        pre_stim_rms = cellfun(@(x) rms(x), pre_stim_rms);
        [pre_stim_rms, sorted_inds] = sort(pre_stim_rms);

        windows = windows(sorted_inds);
        pre_stim = pre_stim(sorted_inds);
        post_stim = post_stim(sorted_inds);
        
        pre_stim_rms_mod = [pre_stim_rms_mod, pre_stim_rms];
        mod_label = [mod_label, ones([1,numel(pre_stim_rms)])*m];
        rms_labels_mod = [rms_labels_mod, ones([1,numel(pre_stim_rms)])*i];


        offset = 0.15;
        for j = 1:numel(trial_starts)
            if mod(i,2)==0
                stim_len = 0.15;
            else
                stim_len = 0.075;
            end

            patch([0,stim_len,stim_len,0], [offset*count - offset/2, offset*count - offset/2,...
                offset*count + offset/2, offset*count + offset/2],...
                stim_col,'edgecolor', stim_col);
            line_y = offset*count;
            line_heights = [line_heights, line_y];
            plot([range(1):1/30000:range(2)], windows{j}{1}+line_y, 'k')
            plot([stim_len-1/30000:1/30000:range(2)], post_stim{j}{1}(stim_len*30000:end)+offset*count,...
                'color', line_colors{i})
            
            plot([range(1):1/30000:range(2)], windows{j}{2}*2000+line_y, 'r')
            plot([stim_len-1/30000:1/30000:range(2)], post_stim{j}{2}(stim_len*30000:end)*2000+offset*count,...
                'color', 'r')
            count = count + 1;
        end
        count = count + 2;
    end
    ylim([-offset*2, offset*(count+3)])
    pre_stim_rms_all = [pre_stim_rms_all, pre_stim_rms_mod];
    rms_labels_all = [rms_labels_all, rms_labels_mod];
    xlabel('Time (s)', 'FontSize', 16 , 'FontWeight', 'Bold');
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    box off
end
trial_labels = cell([1, numel(rms_labels_all)]);
trial_labels(ismember(rms_labels_all, [1,2])) = {'hit'};
trial_labels(ismember(rms_labels_all, [3,4])) = {'CR'};
trial_labels(ismember(rms_labels_all, [5,6])) = {'miss'};

trial_labels(mod_label==1) = cellfun(@(x) ['touch_',x], trial_labels(mod_label==1),'uni', 0);
trial_labels(mod_label==2) = cellfun(@(x) ['visual_',x], trial_labels(mod_label==2),'uni', 0);


touch_90p = prctile(pre_stim_rms_all(mod_label == 1),90);
visual_90p = prctile(pre_stim_rms_all(mod_label == 2),90);
overall_90p = prctile(pre_stim_rms_all,90);

rms_summary = {mouse_name, a.session_name, pre_stim_rms_all, trial_labels, overall_90p};

over_90_touch = line_heights(mod_label == 1 & pre_stim_rms_all > overall_90p);
over_90_visual = line_heights(mod_label == 2 & pre_stim_rms_all > overall_90p);

scatter(subplots(1), -ones([1,numel(over_90_touch)]), over_90_touch, 'ro', 'filled')
scatter(subplots(2), -ones([1,numel(over_90_visual)]), over_90_visual, 'ro', 'filled')


subaxis(1,3,3, 'Margintop', 0.06, 'Marginbottom', 0.06)
handles = plotSpread(pre_stim_rms_all,'distributionIdx', mod_label, 'categoryIdx', rms_labels_all,...
   'binWidth',0.001, 'xyori', 'flipped', 'spreadwidth', 1.25);
plot([touch_90p, touch_90p], [0.9, 1.1], 'k')
plot([visual_90p, visual_90p], [1.9, 2.1], 'k')
plot([overall_90p, overall_90p], [0.8, 2.2], 'r')


ylim([0.4, 2.7])
curr_xlim = get(gca, 'XLim');
xlim([0, 0.03])

text(0.05, 0.2, 'Touch Stim', 'color', stim_colors{1}, 'Units','normalized', 'rotation', 90, 'fontsize', 16)
text(0.05, 0.66, 'Visual Stim', 'color', stim_colors{2},'Units', 'normalized', 'rotation', 90, 'fontsize', 16)

text(0.7, 0.82, 'Hits', 'color', line_colors{1}, 'Units', 'normalized','fontsize', 16)
text(0.7, 0.85, 'CR', 'color', line_colors{3},'Units', 'normalized', 'fontsize', 16)
text(0.7, 0.88, 'Miss', 'color', line_colors{5}, 'Units', 'normalized','fontsize', 16)

text(0.7, 0.5, '90th percentile', 'color', 'r', 'Units', 'normalized','fontsize', 16)
text(0.7, 0.47, '90th percentile(within modality)', 'color', 'k', 'Units', 'normalized','fontsize', 16)


set(gcf, 'Position', [0 0 1920 1060])
set(gcf,'color','w');
box off
xlabel('rms', 'FontSize', 16 , 'FontWeight', 'Bold');
set(gca,'ytick',[])
set(gca,'yticklabel',[])

markerColors = [line_colors, line_colors];
for h=1:numel(handles{1})
    set(handles{1}(h), 'MarkerSize', 20)
    set(handles{1}(h), 'MarkerFaceColor', 'k')
    set(handles{1}(h), 'MarkerEdgeColor', markerColors{h})
end
