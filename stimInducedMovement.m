%%% To do:
%%% -Load example session's raw data
%%% -export and process median motion artifact signal
%%% -cut up signal into trials using previously processed intan/bcont data
%%% -plot motion artifact
%%% -see if calculating moving binned rms is a good idea (bin it at the
%%%     same specs as neural data
%%% -see if CP activity precedes rms

example_dir = 'C:\Users\PC\Documents\Current Data\YT015_FT\09-26-17_FT\';
BcontDir = ['C:\Users\PC\Documents\Current Data\BCont_data', '\'];

% identify session date, direction
cd(example_dir)
[~, seshName, ~] = fileparts(pwd);
seshDate = seshName(1:end-3);
cd ../
[~, mouse_name, ~] = fileparts(pwd);
cd(example_dir)

% get motion artifact signal and perch. Filter them
[rawSignal, sampleRate, rawPerch] = preProcessMotion.get_motion_artifact(example_dir, 'perch', 1);
filteredSignal = preProcessMotion.filterSignal(rawSignal, sampleRate);
filteredPerch = Intan.PawFuncPerch.FiltMovement(rawPerch(2,:), sampleRate);
firstDerivPerch = gradient(filteredPerch);
normedPerch = (firstDerivPerch-min(firstDerivPerch))/(max(firstDerivPerch) - min (firstDerivPerch));

% load behavior and intan data that will be used to identify where 
% to chop up and how to sort filtered signal
b1 = preProcessMotion.processBeh(BcontDir,mouse_name, seshDate);

cd('D:\NeuralData')
load([mouse_name, '_', seshDate, '.mat'])
a.combineData(b1)

%chop up and sort filtered signal
xrange = [-1,1];
parsedSignal = preProcessMotion.parseSignalTrials(a, filteredSignal, sampleRate, xrange);
parsedPerch = preProcessMotion.parseSignalTrials(a, normedPerch, sampleRate, xrange);
sortedTrials  = preProcessMotion.sortByTrialType(a);

sortedTrialNums = sortedTrials(:,2);
signalTrialNums = cell2mat(parsedSignal(:,1));

sortedSignal_inds = cellfun(@(x) find(ismember(signalTrialNums, x)), sortedTrialNums, 'uni', 0);
sortedSignal = cellfun(@(x)  parsedSignal(x,:), sortedSignal_inds, 'uni', 0);
sortedPerch = cellfun(@(x)  parsedPerch(x,:), sortedSignal_inds, 'uni', 0);


%calculate rms of binned signal
bin_size = 0.01;
binnedRMS = preProcessMotion.motionRMS(sortedSignal, xrange, bin_size, sampleRate);

% range should be identical to range used in parsedSignalTrials
% function
set(0,'defaultAxesFontSize',18)
[~, preStimRMS_inds] = cellfun(@(x) sort(cell2mat(x(:,3))), sortedSignal, 'uni',0);

binnedRMS = cellfun(@(x,y) x(y,:), binnedRMS, preStimRMS_inds, 'uni',0);
sortedPerch = cellfun(@(x,y) x(y,:), sortedPerch, preStimRMS_inds, 'uni',0);

rts = preProcessMotion.get_reactionTimes(sortedTrials);
[rts, rt_inds] = cellfun(@(x) sort(x), rts, 'uni',0);

binnedRMS([1,7]) = cellfun(@(x,y) x(y,:), binnedRMS([1,7]), rt_inds([1,7]), 'uni',0);
sortedPerch([1,7]) = cellfun(@(x,y) x(y,:), sortedPerch([1,7]), rt_inds([1,7]), 'uni',0);


fig = plotBinnedRMS(binnedRMS, xrange, bin_size, sampleRate, [0,0.05], rts);
fig = plotBinnedPerch(sortedPerch, xrange, 1/30000, sampleRate, [0.35,0.65], rts);
fig = plotPreStimMove3(sortedSignal, xrange, sampleRate);

cd('C:\Users\PC\Documents\Current Data\motion_artifact_analysis\Figures\trial_exclusion_criteria')
savefig(fig, [mouse_name, '_', seshDate, '.fig'])
set(gcf,'PaperPositionMode','auto')
print(fig, [mouse_name, '_', seshDate],'-dtiffn','-r0')
close(fig)
        