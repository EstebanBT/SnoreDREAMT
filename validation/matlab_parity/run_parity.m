addpath('/home/esteban/toolbox/fieldtrip-20230716'); ft_defaults;
addpath('/home/esteban/BreathREM');
load(fullfile(fileparts(mfilename('fullpath')), 'seg.mat'));
data = []; data.label = {'F4-M1'; 'C4-M1'}; data.fsample = 100; data.trial = {X};
data.time = {(0:size(X,2)-1)/100}; data.sampleinfo = [1 size(X,2)]; data.scoring = {scoring};
cfg = []; cfg.bpfreq = [0.3 1.5]; cfg.slStages = [2 3]; cfg.thresCh = data.label; cfg.detectCh = data.label;
cfg.thresType = 'channelwise'; cfg.filtType = 'but'; cfg.criterionLen = [0.8 2]; cfg.criterionParam = 'both';
cfg.criterionCenter = 'median'; cfg.criterionVar = 'scaleCenter'; cfg.criterionScale = [2.5 2.5];
cfg.criterionPadding = [0 0]; cfg.findEvtFree = 0; cfg.findEvtFreeWin = 300; cfg.raw_trough = false;
so = slythm.detectSOs_splitStages(cfg, data);
cfg = []; cfg.bpfreq = [12 16]; cfg.slStages = [2 3]; cfg.thresCh = data.label; cfg.detectCh = data.label;
cfg.thresType = 'channelwise'; cfg.filtType = 'but'; cfg.envelopeType = 'rms'; cfg.envelopeWin = 0.2;
cfg.criterionLen = [0.5 3]; cfg.criterionCenter = 'median'; cfg.criterionVar = 'scaleCenter'; cfg.criterionScale = 1.5;
cfg.criterionPadding = [0 0]; cfg.searchIndivEvtPeak = 0; cfg.searchMargin = 1; cfg.mergeEvts = 0; cfg.mergeEvtsGap = 0.5;
cfg.doubleThres = 0; cfg.doubleThresScale = 10; cfg.findEvtFree = 0; cfg.findEvtFreeWin = 300;
sp = slythm.detectSpindles(cfg, data);
out = struct();
for c = 1:2
  out.(sprintf('so_tr_%d', c)) = so.EvtInfo(c).minTime;
  out.(sprintf('so_sta_%d', c)) = so.EvtInfo(c).staTime;
  out.(sprintf('sp_sta_%d', c)) = sp.EvtInfo(c).staTime;
  out.(sprintf('sp_end_%d', c)) = sp.EvtInfo(c).endTime;
end
out.so_thr = so.EvtThres; out.so_thr2 = so.EvtThres2; out.sp_thr = sp.EvtThres;
save(fullfile(fileparts(mfilename('fullpath')), 'matlab_out.mat'), '-struct', 'out');
