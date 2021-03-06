function fig = plotBinnedPerch(binnedPerch, range, bin, sampleRate, clims, rts)

fig = figure; hold on
% set(fig, 'Visible', 'off');

numTT = numel(binnedPerch);
touch_trials = binnedPerch(1:(numTT/2));
visual_trials = binnedPerch((numTT/2)+1:end);
tt = {touch_trials, visual_trials};
hit_color = [76, 167, 51]/256;
cr_color = [145, 104, 191]/256;
miss_color = [129, 129, 129]/256;
line_colors = {hit_color, cr_color, miss_color, hit_color, cr_color, miss_color};
stim_colors = {[51, 119, 182]/256, [246, 130, 0]/256};
subplots=[];
line_heights = [];

rts = [{horzcat(rts{1:6})}, {horzcat(rts{6:12})}];
for modality=1:numel(tt)

    sp = subaxis(1,3,modality, 'Margintop', 0.06, 'Marginbottom', 0.06); hold on;
    subplots = [subplots, sp];
    modality_perch = tt{modality};
    stim_col = stim_colors{modality};
    count = 0;
    
    modality_perch = vertcat(modality_perch{:});
    modality_perch = cell2mat(modality_perch);
    
    xvals = [range(1):bin:range(2)];
    imagesc(xvals, 1:size(modality_perch,1), flipud(modality_perch))
    cmocean('balance', 'negative')
    set(gca,'YDir','reverse');
    set(gca, 'TickDir', 'out')
    caxis(clims)
    box off
    hold on
    
    xlim([xvals(1), xvals(end)])
    ylim([0,size(modality_perch,1)])
    xlabel('Time (s)', 'FontSize', 16 , 'FontWeight', 'Bold');
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    box off
    
    plot([0,0], [0,size(modality_perch,1)], ':w')
    plot([0.15,0.15], [0,size(modality_perch,1)], ':w')
    
    rts_mod = rts{modality};
    num_trials = size(modality_perch,1);
    for tn=1:numel(rts{modality})
        scatter(rts_mod(tn), num_trials+1-tn, 'wo', 'markerfacecolor', 'w')
    end
end
cb = colorbar([.64 0.06 .02 0.87]);
ylabel(cb, 'a.u.')
% 
% touch_Perch = cellfun(@(x) sort(cell2mat(x(:,3))), touch_trials, 'uni',0);
% touch_Perch = vertcat(touch_Perch{:});
% vis_Perch = cellfun(@(x) sort(cell2mat(x(:,3))), visual_trials, 'uni',0);
% vis_Perch = vertcat(vis_Perch{:});
% 
% touch_90p = prctile(touch_Perch,90);
% visual_90p = prctile(vis_Perch,90);
% overall_90p = prctile([touch_Perch; vis_Perch],90);
% over_90_touch = line_heights{1}(touch_Perch > overall_90p);
% over_90_visual = line_heights{2}(vis_Perch > overall_90p);
% 
% scatter(subplots(1), -ones([1,numel(over_90_touch)]), over_90_touch, 'ro', 'filled')
% scatter(subplots(2), -ones([1,numel(over_90_visual)]), over_90_visual, 'ro', 'filled')
% 
% mod_label = [ones(size(touch_Perch)); ones(size(vis_Perch))*2];
% Perch_labels = num2cell([1,2,3,1,2,3,1,2,3,1,2,3])';
% Perch_labels = Perch_labels(1:numTT);
% Perch_labels_all = cell2mat(cellfun(@(x,y) ones(size(x(:,1)))*y, sorted_perch, Perch_labels, 'uni', 0));
% 
% subaxis(1,3,3, 'Margintop', 0.06, 'Marginbottom', 0.06)
% handles = plotSpread([touch_Perch; vis_Perch],'distributionIdx', mod_label, 'categoryIdx', Perch_labels_all,...
%    'binWidth',0.001, 'xyori', 'flipped', 'spreadwidth', 1.25);
% plot([touch_90p, touch_90p], [0.9, 1.1], 'k')
% plot([visual_90p, visual_90p], [1.9, 2.1], 'k')
% plot([overall_90p, overall_90p], [0.8, 2.2], 'r')
% 
% ylim([0.4, 2.7])
% curr_xlim = get(gca, 'XLim');
% xlim([0, 0.03])
% 
% text(0.05, 0.2, 'Touch Stim', 'color', stim_colors{1}, 'Units','normalized', 'rotation', 90, 'fontsize', 16)
% text(0.05, 0.66, 'Visual Stim', 'color', stim_colors{2},'Units', 'normalized', 'rotation', 90, 'fontsize', 16)
% text(0.7, 0.82, 'Hits', 'color', line_colors{1}, 'Units', 'normalized','fontsize', 16)
% text(0.7, 0.85, 'CR', 'color', line_colors{3},'Units', 'normalized', 'fontsize', 16)
% text(0.7, 0.88, 'Miss', 'color', line_colors{5}, 'Units', 'normalized','fontsize', 16)
% text(0.7, 0.5, '90th percentile', 'color', 'r', 'Units', 'normalized','fontsize', 16)
% text(0.7, 0.47, '90th percentile(within modality)', 'color', 'k', 'Units', 'normalized','fontsize', 16)
% 
set(gcf, 'Position', [0 0 1920 1060])
set(gcf,'color','w');
box off
% xlabel('Perch', 'FontSize', 16 , 'FontWeight', 'Bold');
% set(gca,'ytick',[])
% set(gca,'yticklabel',[])
% 
% markerColors = [line_colors, line_colors];
% for h=1:numel(handles{1})
%     set(handles{1}(h), 'MarkerSize', 20)
%     set(handles{1}(h), 'MarkerFaceColor', 'k')
%     set(handles{1}(h), 'MarkerEdgeColor', markerColors{h})
% end
