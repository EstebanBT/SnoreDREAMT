% CircStat reference values for the Python port (snoredreamt/circstat.py). Run from this folder:
%   matlab -batch "run('validation/circstat_parity/run_parity.m')"
addpath('/home/esteban/toolbox/CircStat2012a');
here = fileparts(mfilename('fullpath'));
load(fullfile(here, 'inputs.mat'));   % a1, a2 (column vectors, radians)
[p_ww, tab] = circ_wwtest(a1, a2); F_ww = tab{2, 5};
[p_cm, med_cm, P_cm] = circ_cmtest(a1, a2);
med1 = circ_median(a1);
[h_mt, mu_mt, ul_mt, ll_mt] = circ_mtest(a1, 0);
p_med = circ_medtest(a1, 0);
save(fullfile(here, 'matlab_out.mat'), 'p_ww', 'F_ww', 'p_cm', 'P_cm', 'med1', 'h_mt', 'mu_mt', 'ul_mt', 'll_mt', 'p_med');
