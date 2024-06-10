%%%%%%%%%%%%%%%%%%%
%Either run from directory containing all raw TGS data files for this exposure, or change directory here (or both :) )
pname = 'example_data\tungsten_depth_study';
%%%%%%%%%%%%%%%%%%%
close all;

grat=3.7; % Initial grating spacing guess. Can be (nominally, MIT) 1.6, 1.9, 2.2, 2.5, 5.8, 6.4, 7.0, 7.6, 8.2, 8.8, 9.4
                                                                 % 2.8, 3.1, 3.4, 3.7, 4.0, 4.2, 4.4, 4.6, 4.9, 5.2, 5.5 (layout of phase mask B)
posstr='Tungsten_Calibration-2022-05-19-06.40um-spot1-POS-1.txt';
negstr='Tungsten_Calibration-2022-05-19-06.40um-spot1-NEG-1.txt';
%Baseline handling if no baseline arguments are provided, the code still works - it just doesn't do any baseline correction
baselineBool = 0; % 1 for "yes do baseline subtraction", 0 for "no don't do that" 
POSbaselineStr = pname + string('\Tungsten_Calibration-2022-05-19-03.70um-baseline') + string('-POS-1.txt'); %ADJUST
NEGbaselineStr = pname + string('\Tungsten_Calibration-2022-05-19-03.70um-baseline') + string('-NEG-1.txt'); %ADJUST

[freq_final,freq_error,speed,diffusivity,diffusivity_err,tau, tauErr, paramrA, AErr, ParamBeta, BetaErr, paramB, BErr, paramTheta, thetaErr, paramC, CErr, file_date_time] = TGSPhaseAnalysis(strcat(pname,posstr),strcat(pname,negstr),grat,2,0,baselineBool, POSbaselineStr, NEGbaselineStr);
disp('Frequency = ');
disp(freq_final);
disp('Frequency error = ');
disp(freq_error);
disp('=> grating spacing on tungsten sample(m): ');
disp(2665.9/freq_final);
disp('with and error of (m): ');
disp(2665.9*freq_error/(freq_final*freq_final));
disp('Thermal diffusivity = ');
disp(diffusivity);
disp('Thermal diffusivity error = ');
disp(diffusivity_err);

close all % comment this out if you want to keep the plot windows at the end
