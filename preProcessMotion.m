classdef preProcessMotion
    
    methods(Static)
        
        function [rawSignal, sampleRate] = get_motion_artifact(directory, varargin)
            %export and concatenate median of raw motion artifact
            p = inputParser;
            addRequired(p,'directory',@isstr)
            addParamValue(p, 'perch', 0, @isnumeric)
            parse(p, directory, varargin{:})
            
            directory = p.Results.directory;
            perch = p.Results.perch;
            
            cd(directory)
            rawfiles = arrayfun(@(x) x.name(1:(end)), dir('*.rhd'),'UniformOutput',false);
            
            rawSignal = [];
            session_perch = [];
            for k=1:numel(rawfiles)
                part = Intan.Sweep(char(rawfiles(k)));
                rawSignal = [rawSignal,median(part.rawSignal)];
                if perch == 1
                    session_perch = [session_perch, part.adcInRawSignal([6,7],:)];
                end
            end
            sampleRate = part.sampleRate;
        end
        
        function filtered = filterSignal(rawSignal, sampleRate)
            
            Wn = [300, 800] / (sampleRate/2);
            [b,a]=butter(2,Wn,'bandpass');
            filtered = filtfilt(b, a, rawSignal);
            
            %licking causes big transients that will bias rms calculation.
            %NaN licks to prevent bias.
            badIndexes = abs(filtered) > 0.1;
            filtered(badIndexes) = NaN;
        end
        
        function sortedTrials  = sortByTrialType(a)
            
            t_hits = {a.intanTrials_SsSb(cellfun(@(x) x.response ~= 0 & x.trialCorrect == 1, a.BcontTrials_SsSb))};
            t_CR = {a.intanTrials_SsVb(cellfun(@(x) x.response == 0, a.BcontTrials_SsVb))};
            t_misses = {a.intanTrials_SsSb(cellfun(@(x) x.response == 0, a.BcontTrials_SsSb))};
            
            v_hits = {a.intanTrials_VsVb(cellfun(@(x) x.response ~= 0 & x.trialCorrect == 1, a.BcontTrials_VsVb))};
            v_CR = {a.intanTrials_VsSb(cellfun(@(x) x.response == 0, a.BcontTrials_VsSb))};
            v_misses = {a.intanTrials_VsVb(cellfun(@(x) x.response == 0, a.BcontTrials_VsVb))};
            
            % some mice didnt get the single cycle stimulus
            try
                ts_hits = {a.intanTrials_SscSb(cellfun(@(x) x.response == 1, a.BcontTrials_SscSb))};
                ts_CR = {a.intanTrials_SscVb(cellfun(@(x) x.response == 0, a.BcontTrials_SscVb))};
                ts_misses = {a.intanTrials_SscSb(cellfun(@(x) x.response == 0, a.BcontTrials_SscSb))};

                vs_hits = {a.intanTrials_VscVb(cellfun(@(x) x.response == 2, a.BcontTrials_VscVb))};
                vs_CR = {a.intanTrials_VscSb(cellfun(@(x) x.response == 0, a.BcontTrials_VscSb))};
                vs_misses = {a.intanTrials_VscVb(cellfun(@(x) x.response == 0, a.BcontTrials_VscVb))};
                
                ttLabels = [{'t_hits'}; {'t_CR'}; {'t_misses'}; {'ts_hits'}; {'ts_CR'}; {'ts_misses'}; {'v_hits'};...
                {'v_CR'}; {'v_misses'}; {'vs_hits'}; {'vs_CR'}; {'vs_misses'}];
            catch
                ts_hits = {};
                ts_CR = {};
                ts_misses = {};

                vs_hits = {};
                vs_CR = {};
                vs_misses = {};
                ttLabels = [{'t_hits'}; {'t_CR'}; {'t_misses'};{'v_hits'};...
                {'v_CR'}; {'v_misses'};];
            end

            sortedtrialTypes = [t_hits, t_CR, t_misses, ts_hits, ts_CR, ts_misses, v_hits, v_CR, v_misses, vs_hits, vs_CR, vs_misses]';
            
            sortedTrialNums = [];
            for tt = 1:numel(sortedtrialTypes)
                tn = {cellfun(@(x) x.trialNum, sortedtrialTypes{tt})};
                sortedTrialNums = [sortedTrialNums; tn];
            end
            
            sortedTrials = [ttLabels, sortedTrialNums, sortedtrialTypes];   
        end
        
        function parsedSignal =  parseSignalTrials(a, filtSignal, sampleRate, range)
            
            %get timestamp of the start of each trial and each trial's stim
            %onset relative to the start of the session.
            trialNums = cellfun(@(x) x.trialNum(1), a.correctedIntanTrials);
            trial_starts = cellfun(@(x) x.rawTime(1), a.correctedIntanTrials);
            stim_onsets = cellfun(@(x) x.stimOnsetTime, a.correctedIntanTrials);
            aligned_so = stim_onsets+trial_starts;
            
            %mutliply the stim onset time by the sampling rate will give
            %the index of that timestamp
            stim_onset_inds = aligned_so*sampleRate;
            stim_onset_inds = uint64(stim_onset_inds);
            
            %find the indeces of every timepoint around a fixed interval
            %around each stim.
            window_starts = stim_onset_inds + range(1)*30000;
            window_ends = stim_onset_inds + range(2)*30000;
            pSignal = arrayfun(@(x,y) filtSignal(x:y), window_starts, window_ends, 'uni', 0);
            pre_stim = arrayfun(@(x,y) filtSignal(x:y), window_starts, stim_onset_inds, 'uni', 0);
            
            %calculate pre-stim rms
            pre_stim_rms = cellfun(@(x) x(~isnan(x)), pre_stim, 'uni', 0);
            pre_stim_rms = cellfun(@(x) rms(x), pre_stim_rms);
            
            parsedSignal = [num2cell(trialNums)', pSignal', num2cell(pre_stim_rms)'];
        end
        
        function b1 = processBeh(BcontDir,mouse_name, seshDate)
            cd([BcontDir, mouse_name])
            rawSeshDate = seshDate([7,8,1,2,4,5]);
            bcont_type = cell2mat(arrayfun(@(x) x.name(1:(end)), dir(['*',num2str(rawSeshDate),'*']),'UniformOutput',false));

            if strfind(bcont_type, 'switch')
                b1 = Solo.EFcross3_switchArray([BcontDir, mouse_name, '\data_@EFcross3_switchobj_', mouse_name,'_', rawSeshDate, 'a'],...
                    [mouse_name, '_', rawSeshDate, 'a']);
            elseif strfind(bcont_type, 'EFcross3rev')
                b1 = Solo.EFcross3revArray([BcontDir, mouse_name, '\data_@EFcross3revobj_', mouse_name,'_', rawSeshDate, 'a'],...
                    [mouse_name, '_', rawSeshDate, 'a']);
            elseif strfind(bcont_type, 'EFcross3')
                b1 = Solo.efcross3Array([BcontDir, mouse_name, '\data_@EFcross3obj_', mouse_name,'_', rawSeshDate, 'a'],...
                    [mouse_name, '_', rawSeshDate, 'a']);
            elseif strfind(bcont_type, 'EFcross2')
                b1 = Solo.EFcross2Array([BcontDir, mouse_name, '\data_@EFcross2obj_', mouse_name,'_', rawSeshDate, 'a'],...
                    [mouse_name, '_', rawSeshDate, 'a']);
            else
                error('Do not have appropriate script to process behavior data type')
            end
        end
        function binnedRMS = motionRMS(sortedSignal, xrange, bin, samplingRate)
            %input sorted motion artifact recordings, bin, calculate and
            %return rms of binned motion artifact aligned to stim onset
            
            edgeInds = [1:bin*samplingRate:(xrange(2)+abs(xrange(1))*samplingRate)];
            binnedRMS = [];
            for tt = 1:numel(sortedSignal)
                
                trialType = sortedSignal{tt};
                binnedRMS_trialType = [];
                for t = 1:size(trialType,1)
                    trial = trialType{t,2};
                    signalBins = arrayfun(@(x,y) trial(x:y), edgeInds(1:end-1), edgeInds(2:end), 'uni', 0);
                    binnedRMS_trial = cellfun(@(x) rms(x), signalBins, 'uni',0);
                    binnedRMS_trialType = [binnedRMS_trialType; binnedRMS_trial];

                end
                binnedRMS = [binnedRMS; {binnedRMS_trialType}];
            end
        end           
    end
end
