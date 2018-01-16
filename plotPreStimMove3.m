function fig = plotPreStimMove3(sorted_signal, range, sampleRate)

fig = figure; hold on
set(fig, 'Visible', 'off');

numTT = numel(sorted_signal);
touch_trials = sorted_signal(1:(numTT/2));
visual_trials = sorted_signal((numTT/2)+1:end);
tt = {touch_trials, visual_trials};
hit_color = [76, 167, 51]/256;
cr_color = [145, 104, 191]/256;
miss_color = [129, 129, 129]/256;
line_colors = {hit_color, cr_color, miss_color, hit_color, cr_color, miss_color};
stim_colors = {[51, 119, 182]/256, [246, 130, 0]/256};
subplots=[];
line_heights = [];

for modality=1:numel(tt)
    
sp = subaxis(1,3,modality, 'Margintop', 0.06, 'Marginbottom', 0.06); hold on;
subplots = [subplots, sp];
modality_signal = tt{modality};
stim_col = stim_colors{modality};
count = 0;

line_heights_mod = [];
    for i = 1:numel(modality_signal)
        offset = 1;
        current_trial_type = modality_signal{i};
        [~,rms_sorted_inds] = sort(cell2mat(current_trial_type(:,3)));
        trial_type_signal = current_trial_type(rms_sorted_inds, 2);
        for j = 1:size(trial_type_signal,1)
            if ismember(i,1:3)
                stim_len = 0.15;
            else
                stim_len = 0.075;
            end

            patch([0,stim_len,stim_len,0], [offset*count - offset/2, offset*count - offset/2,...
                offset*count + offset/2, offset*count + offset/2],...
                stim_col,'edgecolor', stim_col);
            line_y = offset*count;
            line_heights_mod = [line_heights_mod, line_y];
            xvals = range(1):1/sampleRate:range(2);
            plot(xvals, trial_type_signal{j,:}*20+line_y, 'k')
            
            postStim_xvals = stim_len- 1/sampleRate:1/sampleRate:range(2);
            postStim = trial_type_signal{j}((abs(range(1)) + stim_len)*30000:end);
            plot(postStim_xvals, postStim*20 + line_y, 'color', line_colors{i})
            count = count + 1;
        end
        count = count + 2;
    end
    line_heights = [line_heights, {line_heights_mod}];
    ylim([-offset*2, offset*(count+3)])
    xlim(range)

%     pre_stim_rms_all = [pre_stim_rms_all, pre_stim_rms_mod];
%     rms_labels_all = [rms_labels_all, rms_labels_mod];
    xlabel('Time (s)', 'FontSize', 16 , 'FontWeight', 'Bold');
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    box off
end

% 
% trial_labels = cell([1, numel(rms_labels_all)]);
% trial_labels(ismember(rms_labels_all, [1,2])) = {'hit'};
% trial_labels(ismember(rms_labels_all, [3,4])) = {'CR'};
% trial_labels(ismember(rms_labels_all, [5,6])) = {'miss'};
% 
% trial_labels(mod_label==1) = cellfun(@(x) ['touch_',x], trial_labels(mod_label==1),'uni', 0);
% trial_labels(mod_label==2) = cellfun(@(x) ['visual_',x], trial_labels(mod_label==2),'uni', 0);

touch_rms = cellfun(@(x) sort(cell2mat(x(:,3))), touch_trials, 'uni',0);
touch_rms = vertcat(touch_rms{:});

vis_rms = cellfun(@(x) sort(cell2mat(x(:,3))), visual_trials, 'uni',0);
vis_rms = vertcat(vis_rms{:});

touch_90p = prctile(touch_rms,90);
visual_90p = prctile(vis_rms,90);
overall_90p = prctile([touch_rms; vis_rms],90);

% rms_summary = {mouse_name, a.session_name, pre_stim_rms_all, trial_labels, overall_90p};

over_90_touch = line_heights{1}(touch_rms > overall_90p);
over_90_visual = line_heights{2}(vis_rms > overall_90p);

scatter(subplots(1), -ones([1,numel(over_90_touch)]), over_90_touch, 'ro', 'filled')
scatter(subplots(2), -ones([1,numel(over_90_visual)]), over_90_visual, 'ro', 'filled')

mod_label = [ones(size(touch_rms)); ones(size(vis_rms))*2];
rms_labels = num2cell([1,2,3,1,2,3,1,2,3,1,2,3])';
rms_labels = rms_labels(1:numTT);
rms_labels_all = cell2mat(cellfun(@(x,y) ones(size(x(:,1)))*y, sorted_signal, rms_labels, 'uni', 0));
subaxis(1,3,3, 'Margintop', 0.06, 'Marginbottom', 0.06)
handles = plotSpread([touch_rms; vis_rms],'distributionIdx', mod_label, 'categoryIdx', rms_labels_all,...
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