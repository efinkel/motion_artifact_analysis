%%% To do:
%%% -Load example session's raw data
%%% -export and process median motion artifact signal
%%% -cut up signal into trials using previously processed intan/bcont data
%%% -plot motion artifact
%%% -see if calculating moving binned rms is a good idea (bin it at the
%%%     same specs as neural data
%%% -see if CP activity precedes rms

example_dir = 'E:\EF0077-EF0090\EF0079\03-09-16_FT\';

% identify session date, direction
cd(example_dir)
[~, seshName, ~] = fileparts(pwd);
seshDate = seshName(1:end-3);
cd ../
[~, mouse_name, ~] = fileparts(pwd);
cd(example_dir)


% get motion artifact signal and filter it
[rawSignal, sampleRate] = Intan.preProcessMotion.get_motion_artifact(session_dir, 'perch', 0);
filteredSignal = Intan.preProcessMotion.filterSignal(rawSignal, sampleRate);

% load behavior and intan data that will be used to identify where 
% to chop up and how to sort filtered signal
b1 = Intan.preProcessMotion.processBeh(BcontDir,mouse_name, seshDate);

cd('D:\NeuralData')
load([mouse_name, '_', seshDate, '.mat'])
a.combineData(b1)

%chop up and sort filtered signal
xrange = [-1,1];
parsedSignal = Intan.preProcessMotion.parseSignalTrials(a, filteredSignal, sampleRate, xrange);
sortedTrials  = Intan.preProcessMotion.sortByTrialType(a);

sortedTrialNums = sortedTrials(:,2);
signalTrialNums = cell2mat(parsedSignal(:,1));

sortedSignal_inds = cellfun(@(x) find(ismember(signalTrialNums, x)), sortedTrialNums, 'uni', 0);
sortedSignal = cellfun(@(x)  parsedSignal(x,:), sortedSignal_inds, 'uni', 0);

% range should be identical to range used in parsedSignalTrials
% function
fig = plotPreStimMove3(sortedSignal, xrange, sampleRate);
cd('C:\Users\PC\Documents\Current Data\motion_artifact_analysis\Figures\trial_exclusion_criteria')
savefig(fig, [mouse_name, '_', seshDate, '.fig'])
set(gcf,'PaperPositionMode','auto')
print(fig, [mouse_name, '_', seshDate],'-dtiffn','-r0')
close(fig)
        