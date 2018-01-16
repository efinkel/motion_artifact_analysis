function fig = plotPreStimMove(mouse_name, filtSignal, perch, filtPerch, time, a, clust, range, plotType)

opengl software
%%% identify lick artifacts and replace them with NaN values
badIndexes = abs(filtSignal) > 0.2;
filtSignal(badIndexes) = NaN;

%%% group trials by type
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

%%
touch_trials = {t_hits,ts_hits, t_CR, ts_CR, t_misses, ts_misses}; 
touch_trials = {t_hits, t_CR, t_misses}; 

visual_trials = {v_hits, vs_hits, v_CR,vs_CR, v_misses,vs_misses};
visual_trials = {v_hits, v_CR, v_misses};


tt = {touch_trials, visual_trials};

hit_color = [76, 167, 51]/256;
cr_color = [145, 104, 191]/256;
miss_color = [129, 129, 129]/256;
line_colors = {hit_color, hit_color, cr_color, cr_color, miss_color, miss_color};
line_colors = {hit_color, cr_color, miss_color};

stim_colors = {[51, 119, 182]/256, [246, 130, 0]/256};

perch_all = [];
maps = cell(1,2);
all_RTs = {};
axs= [];
sp = [1,3];
hists_all = [];
for m = 1:2
    axs = [axs, subaxis(2,4,sp(m), 'Margintop', 0.09, 'Marginbottom', 0.09, 'Marginleft', 0.06, 'MarginRight', 0.06, 'padding', 0.01)];
    hold on
    trial_types = tt{m};
    stim_col = stim_colors{m};
    count = 0;
    perch_mod = [];
    hists_mod = [];
   
    %%% cut perch, motion artifact, and spikeing activity into trials
    for i=1:numel(trial_types)
        trial_starts = cellfun(@(x) x.rawTime(1), trial_types{i});
        stim_onsets = cellfun(@(x) x.stimOnsetTime, trial_types{i});
        licks = cellfun(@(x) [x.rPreciseTrialLickTimes; x.lPreciseTrialLickTimes], trial_types{i}, 'uni', 0);
        alicks = cellfun(@(x,y) x-y, licks, num2cell(stim_onsets), 'uni', 0);
        aligned_first_licks = cellfun(@(x) min(x(x>0.1 & x<2)), alicks, 'uni', 0);

        trial_start_inds = find(ismember(time,trial_starts));
        stim_onset_inds = stim_onsets*30000 + trial_start_inds;

        window_starts = stim_onset_inds + range(1)*30000;
        window_ends = stim_onset_inds + range(2)*30000;
        
        windows = arrayfun(@(x,y) filtSignal(x:y), window_starts, window_ends, 'uni', 0);
        pre_stim = arrayfun(@(x,y) filtSignal(x:y), window_starts, stim_onset_inds, 'uni', 0);
        post_stim = arrayfun(@(x,y) filtSignal(x:y), stim_onset_inds, window_ends, 'uni', 0);
        windowsPerch = arrayfun(@(x,y) filtPerch(x:y), window_starts, window_ends, 'uni', 0)';
        perch_window = arrayfun(@(x,y) perch(x:y), window_starts, window_ends, 'uni', 0)';
        
        %%% sort all signals by trial RT
        if i == 1
            [aligned_first_licks, sorted_inds] = sort(horzcat(aligned_first_licks{:}));
            all_RTs = [all_RTs, {aligned_first_licks}];
        else
            sorted_inds = 1:numel(aligned_first_licks);
        end

        windows = windows(sorted_inds);
        pre_stim = pre_stim(sorted_inds);
        post_stim = post_stim(sorted_inds);
        perch_window = perch_window(sorted_inds);
        
        windowsPerch_sort = cell2mat(windowsPerch(sorted_inds));
        maps{m} = [maps{m}; windowsPerch_sort];
        perch_mod = [perch_mod; {perch_window}];

        offset = 1;
        
        %%%plot spikes
        hist = [];
        for j = 1:numel(trial_starts)
           
            line_y = offset*count;
            trial = trial_types{i}{j};
            
            if strcmp(plotType,'rasters')
                fields = fieldnames(trial);
                tetrodes = fields(end-7:end);
                clusters = cellfun(@(x) getfield(trial,x), tetrodes, 'uni', 0);
                clusters = horzcat(clusters{:});

                aligned_spikes = clusters{clust}-stim_onsets(j);
                window_spikes = aligned_spikes(aligned_spikes > range(1) & aligned_spikes < range(2));
                trial_spRate = histc(window_spikes, range(1):0.025:range(2));
                
                hist = [hist; reshape(trial_spRate, [1,numel(trial_spRate)])];

                arrayfun(@(x) plot([x,x] , [count-0.5, count+0.5], 'k'), window_spikes)
            elseif strcmp(plotType,'artifact')
                 plot([range(1):1/30000:range(2)], windows{j}*5+line_y, 'k')
            end
            
            scatter(range(2), line_y,50,'o', 'filled', 'markerfacecolor', line_colors{i});
            if i==1
                scatter(aligned_first_licks(j), offset*count+.5, 25, 'r','v', 'filled')
            end
            count = count + 1;
        end
        hists_mod = [hists_mod; {hist}];

    end
    hists_all = [hists_all, {hists_mod}];
    perch_all = [perch_all, {perch_mod}];
    ylim([0, offset*(count)])
    xlim(range)
    xlabel('Time (s)', 'FontSize', 16 , 'FontWeight', 'Bold');
    box off
end


maps2 = [{cell2mat(vertcat(perch_all{1}{:}))}, {cell2mat(vertcat(perch_all{2}{:}))}];
sp = [2,4];
C = [1,2];

for h=1:numel(maps2)
    axs = [axs, subaxis(2,4,sp(h), 'Margintop', 0.09, 'Marginbottom', 0.09, 'Marginleft', 0.06, 'MarginRight', 0.06, 'padding', 0.01)];
    xvals = [range(1):1/30000:range(2)];
    imagesc(xvals, 1:size(maps2{h},1), flipud(maps2{h}*C(h)), [-0.3,0.3])
    cmocean('balance', 'negative')
    set(gca,'YDir','reverse');
    set(gca, 'TickDir', 'out')
    box off
    hold on
    numTrials = size(maps2{h},1);

%     for trial=1:numel(perch_all{h})
%         mod_trials = (perch_all{h});
%         tr = mod_trials{trial};
%         baseline_mean = mean(tr(1:abs(range(1))*30000));
%         plot(xvals, (tr- baseline_mean)*C(h)+(trial) , 'color', [0.4,0.4,0.4], 'linewidth', 0.1)
%     end
    xlim([xvals(1), xvals(end)])
    ylim([0,size(maps2{h},1)])
end

for modality=1:numel(all_RTs)
    axes(axs(2+modality))
    hold on
    xvals = (all_RTs{modality});
    yvals = size(maps{modality},1):-1:size(maps{modality},1)-numel(all_RTs{modality}) + .5; %%
%     yvals = 1:numel(all_RTs{modality});
    scatter(xvals, yvals, 25, 'r','v', 'filled')
end

sp2 = [5,7];
ymax = [];
for m = 1:2
    
    axs = [axs, subaxis(2,4,sp2(m), 'Margintop', 0.09, 'Marginbottom', 0.09, 'Marginleft', 0.06, 'MarginRight', 0.06, 'padding', 0.01)]; hold on
   
    for tt = 1:numel(hists_all{m})
        meanFR = mean(hists_all{m}{tt}/0.025,1);
        semFR = (std(hists_all{m}{tt}/0.025,1)/sqrt(size(hists_all{m}{tt},1)));
        xvals = range(1):0.025:range(2);
        if strcmp(plotType,'rasters')
            plot(xvals, meanFR, 'color', line_colors{tt})
            fill([xvals,fliplr(xvals)],[meanFR-semFR,fliplr(meanFR+semFR)], 'k', 'faceColor',...
                line_colors{tt}, 'faceAlpha', 0.4, 'EdgeColor', line_colors{tt}, 'EdgeAlpha', 0.4);
        end
    end
    
    xlim(range)
    y_limits = get(gca, 'ylim');
    ymax = [ymax, y_limits(2)];
    
    
    axs = [axs, subaxis(2,4,sp2(m)+1, 'Margintop', 0.09, 'Marginbottom', 0.09, 'Marginleft', 0.06, 'MarginRight', 0.06, 'padding', 0.01)]; hold on
    
    for tt = 1:numel(hists_all{m})
        trials = cell2mat(perch_all{m}{tt});
        meanFR = mean(trials,1);
        semFR = (std(trials,1)/sqrt(size(trials,1)));
        xvals = range(1):01/30000:range(2);
        if strcmp(plotType,'rasters')
            plot(xvals, meanFR, 'color', line_colors{tt})
            fill([xvals';flipud(xvals')],[meanFR'-semFR';flipud(meanFR'+semFR')], 'k', 'faceColor',...
                line_colors{tt}, 'faceAlpha', 0.4, 'EdgeColor', line_colors{tt}, 'EdgeAlpha', 0.4);
        end
    end
    
    xlim(range)

end
linkaxes(axs([end-3,end-1]), 'y')
set(axs(end-1), 'ylim', [0,max(ymax)])

for ax = 1:numel(axs)
    stim_kind = 1;
    if mod(ax,2) == 0
        stim_kind = 2;
    end
    stim_len = 0.15;
    axes(axs(ax))
    plot([0,0], ylim, 'color', stim_colors{stim_kind},  'linewidth', 2)
    plot([0.15,0.15], ylim, 'color', stim_colors{stim_kind}, 'linewidth', 2)
end
set(gcf, 'Position', [20 300 1900 600])
set(gcf,'color','w');
