
touch_hit = cellfun(@(x,y,z) sum(x(ismember(y, 'touch_hit'))>z)/numel(x(ismember(y, 'touch_hit'))), rms_all(:,3), rms_all(:,4), rms_all(:,5));
touch_CR = cellfun(@(x,y,z) sum(x(ismember(y, 'touch_CR'))>z)/numel(x(ismember(y, 'touch_CR'))), rms_all(:,3), rms_all(:,4), rms_all(:,5));
touch_miss = cellfun(@(x,y,z) sum(x(ismember(y, 'touch_miss'))>z)/numel(x(ismember(y, 'touch_miss'))), rms_all(:,3), rms_all(:,4), rms_all(:,5));
visual_hit = cellfun(@(x,y,z) sum(x(ismember(y, 'visual_hit'))>z)/numel(x(ismember(y, 'visual_hit'))), rms_all(:,3), rms_all(:,4), rms_all(:,5));
visual_CR = cellfun(@(x,y,z) sum(x(ismember(y, 'visual_CR'))>z)/numel(x(ismember(y, 'visual_CR'))), rms_all(:,3), rms_all(:,4), rms_all(:,5));
visual_miss = cellfun(@(x,y,z) sum(x(ismember(y, 'visual_miss'))>z)/numel(x(ismember(y, 'visual_miss'))), rms_all(:,3), rms_all(:,4), rms_all(:,5));

frac_above_90 = [touch_hit, touch_CR, touch_miss, visual_hit, visual_CR, visual_miss];
meanFrac = mean(frac_above_90, 1);
stim_colors = {[51, 119, 182]/256, [246, 130, 0]/256};

fig = figure;
subaxis(2,1,1, 'Margintop', 0.08, 'Marginbottom', 0.08)
hold on

for d = 1:size(frac_above_90,1)
	scatter([1,2,3], frac_above_90(d,1:3), 65, stim_colors{1}, 'filled');
    scatter([5,6,7], frac_above_90(d,4:6), 65, stim_colors{2}, 'filled')
    
    plot([1,2,3,5,6,7], frac_above_90(d,:), 'color', [0.8,0.8, 0.8])
end

tt_pos = [1,2,3,5,6,7];
for tt=1:numel(meanFrac)
    plot([tt_pos(tt)-0.33, tt_pos(tt)+0.33] , [meanFrac(tt), meanFrac(tt)], 'k')
end
set(gca,'xtick',[1,2,3,5,6,7])
set(gca,'xticklabel',{'Hit', 'CR', 'Miss', 'Hit', 'CR', 'Miss'}, 'fontsize', 12)
ylabel('Fraction above overall 90th percentile', 'fontsize', 12)
box off
set(gca, 'TickDir', 'out')

%%
%%fraction of trials above 90th percentile for rms are of a certain trial
%%type

frac_touch_hit = cellfun(@(x,y,z) sum(x(ismember(y, 'touch_hit'))>z)/numel(x(x>z)), rms_all(:,3), rms_all(:,4), rms_all(:,5));
frac_touch_CR = cellfun(@(x,y,z) sum(x(ismember(y, 'touch_CR'))>z)/numel(x(x>z)), rms_all(:,3), rms_all(:,4), rms_all(:,5));
frac_touch_miss = cellfun(@(x,y,z) sum(x(ismember(y, 'touch_miss'))>z)/numel(x(x>z)), rms_all(:,3), rms_all(:,4), rms_all(:,5));
frac_visual_hit = cellfun(@(x,y,z) sum(x(ismember(y, 'visual_hit'))>z)/numel(x(x>z)), rms_all(:,3), rms_all(:,4), rms_all(:,5));
frac_visual_CR = cellfun(@(x,y,z) sum(x(ismember(y, 'visual_CR'))>z)/numel(x(x>z)), rms_all(:,3), rms_all(:,4), rms_all(:,5));
frac_visual_miss = cellfun(@(x,y,z) sum(x(ismember(y, 'visual_miss'))>z)/numel(x(x>z)), rms_all(:,3), rms_all(:,4), rms_all(:,5));

comp_above_90 = [frac_touch_hit, frac_touch_CR, frac_touch_miss, ones([numel(frac_touch_hit), 1])*0.15,...
    frac_visual_hit, frac_visual_CR, frac_visual_miss];

hit_color = [76, 167, 51]/256;
cr_color = [145, 104, 191]/256;
miss_color = [129, 129, 129]/256;
block_colors = {hit_color, cr_color, miss_color, 'w', miss_color, cr_color, hit_color};

subaxis(2,1,2, 'Margintop', 0.08, 'Marginbottom', 0.08, 'PaddingTop', 0.08)
title('Composition of trials above 90th percentile','fontsize', 12)
hold on
for d = 1:size(comp_above_90,1)
    x_spread = 0;
    for tt = 1:numel(comp_above_90(1,:))
        xvals = [x_spread, x_spread+ comp_above_90(d,tt), x_spread+ comp_above_90(d,tt), x_spread];
        yvals = [d, d, d+0.8, d+0.8];
        patch(xvals, yvals, block_colors{tt}, 'edgecolor',block_colors{tt})
        x_spread = x_spread + comp_above_90(d,tt);
    end
end
plot([0.7,0.8], [0.5,0.5], 'k')
text(0.74, 0.15, '0.1', 'color', 'k', 'fontsize', 12)

text(0.15, 0.95, 'Touch Stim', 'color', stim_colors{1}, 'Units','normalized', 'fontsize', 12)
text(0.65, 0.95, 'Visual Stim', 'color', stim_colors{2},'Units', 'normalized', 'fontsize', 12)
xlim([-0.1, 1.35])
% ylim([-2, size(comp_above_90,1)+2])

text(0.9, 0.7, 'Hits', 'color', block_colors{1}, 'Units', 'normalized','fontsize', 12)
text(0.9, 0.75, 'CR', 'color', block_colors{2},'Units', 'normalized', 'fontsize', 12)
text(0.9, 0.8, 'Miss', 'color', block_colors{3}, 'Units', 'normalized','fontsize', 12)

set(gcf, 'Position', [50 50 1000 1000])
set(gcf,'color','w');
box off
set(gca,'ytick',[0:size(comp_above_90,1)], 'fontsize', 12)
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'xcolor','w')
ylabel('Day', 'fontsize', 12)
set(gca, 'TickDir', 'out')
