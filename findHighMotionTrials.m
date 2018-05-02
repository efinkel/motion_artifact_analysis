%%% Script that batch processes motion artifact collected from intan
%%% recordings. Will output a list of trials that have high pre stimulus
%%% motion artifact so that they can be excluded in later analyses. Also
%%% plots figures that visualize the recorded motion artifact for every
%%% trial in a session.

clear
directories = {'D:\EF0074', 'D:\EF0076', 'M:\EF0088', 'M:\EF0089', 'M:\EF0094', 'M:\EF0098', 'M:\EF0099',...
    'E:\EF0077-EF0090\EF0077','E:\EF0077-EF0090\EF0079','E:\EF0077-EF0090\EF0081','E:\EF0077-EF0090\EF0083',...
    'E:\EF0077-EF0090\EF0084', 'E:\EF0077-EF0090\EF0085'};

BcontDir = ['C:\Users\PC\Documents\Current Data\BCont_data', '\'];
trialsToExclude = [];
for d = 12:numel(directories)
    cd(directories{d})
    [~, mouse_name, ~] = fileparts(pwd);
    mouse_sessions = arrayfun(@(x) x.name(1:(end)), dir(['*_FT']),'UniformOutput',false);
    
    for s = 1:numel(mouse_sessions)
        
        % identify session date, direction
        session_dir = [directories{d}, '/' mouse_sessions{s}];
        cd(session_dir)
        [~, seshName, ~] = fileparts(pwd);
        seshDate = seshName(1:end-3);
        
        % get motion artifact signal and filter it
        [rawSignal, sampleRate] = preProcessMotion.get_motion_artifact(session_dir, 'perch', 0);
        filteredSignal = preProcessMotion.filterSignal(rawSignal, sampleRate);
        
        

        % load behavior and intan data that will be used to identify where 
        % to chop up and how to sort filtered signal
        b1 = preProcessMotion.processBeh(BcontDir,mouse_name, seshDate);
        try
            cd('D:\NeuralData')
            load([mouse_name, '_', seshDate, '.mat'])
        catch
            cd('C:\Users\PC\Documents\Current Data\working folder\NeuralData')
            load([mouse_name, '_', seshDate, '.mat'])
        end
        a.combineData(b1)

        
        
        %chop up and sort filtered signal
        xrange = [-1,1];
        parsedSignal = preProcessMotion.parseSignalTrials(a, filteredSignal, sampleRate, xrange);
        sortedTrials  = preProcessMotion.sortByTrialType(a);
        
        sortedTrialNums = sortedTrials(:,2);
        signalTrialNums = cell2mat(parsedSignal(:,1));
        
        sortedSignal_inds = cellfun(@(x) find(ismember(signalTrialNums, x)), sortedTrialNums, 'uni', 0);
        sortedSignal = cellfun(@(x)  parsedSignal(x,:), sortedSignal_inds, 'uni', 0);
        
% %         range should be identical to range used in parsedSignalTrials
% %         function
%         fig = plotPreStimMove3(sortedSignal, xrange, sampleRate);
%         cd('C:\Users\PC\Documents\Current Data\motion_artifact_analysis\Figures\trial_exclusion_criteria')
%         savefig(fig, [mouse_name, '_', seshDate, '.fig'])
%         set(gcf,'PaperPositionMode','auto')
%         print(fig, [mouse_name, '_', seshDate],'-dtiffn','-r0')
%         close(fig)
%         
        %identify which trials to exclude
        rms_trialNums = cell2mat(parsedSignal(:,[1,3]));
        cutoff =  prctile(rms_trialNums(:,2),90);
        cut_trials = rms_trialNums(rms_trialNums(:,2) > cutoff, 1);
        
        % append to cell table
        mouse_name_rep = cell(size(cut_trials)); mouse_name_rep(:) = {mouse_name};
        date_rep = cell(size(cut_trials)); date_rep(:) = {seshDate};
        section = [mouse_name_rep,date_rep, num2cell(cut_trials)];
        trialsToExclude = [trialsToExclude; section];
    end
end

    
    
