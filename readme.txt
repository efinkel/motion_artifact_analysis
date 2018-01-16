to do:

- starting input: mouse main directories
- list all sessions for each mouse

- extract and concatenate median electrode signal, sampling rate, and maybe perch signal
- filter median signal and maybe perch signal
- identify intant data object and behavior data object. Open both and combine.
- cut up signals according to trial structures
- calculate rms of median signal for each pre stimulus period of each trial

####
- plot summary figure for each day using available script and save
- find trial numbers of trials that fall within the top 10% of pre trial_rms
- append a row for each trial within a session to a cell array containing following info: mouse_name, session_date, trial_num
