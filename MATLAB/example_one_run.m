%close all current plots
close all;
% directory location
pname = 'C:\Users\apcwy\OneDrive - Massachusetts Institute of Technology\MIT_postdoc\TGS\Data\PI3\YBCO_uncoated\YBCO_Tape_3'; %ADJUST
%First part of file name
str_base='YBCO_uncoated-YBCO_Tape_3-05.80um-spotForIrradiation_2';  %ADJUST
% Calibrated grating spacings
grat = 5.8; %in um %ADJUST

% baseline handling
blbool = 1;
baselinePOS = pname + string('\') + 'YBCO_uncoated-YBCO_Tape_3-05.20um-baselinePreIrradiation-POS-1.txt';
baselineNEG = pname + string('\') + 'YBCO_uncoated-YBCO_Tape_3-05.20um-baselinePreIrradiation-NEG-1.txt';

printFile = pname + string('\') + str_base + string('_postprocessing.txt');
fid1 = fopen( printFile, 'wt' );

posstr = strcat(pname,string('\'),str_base,'-POS-1.txt');
negstr = strcat(pname,string('\'),str_base,'-NEG-1.txt');

[freq_final,freq_error,speed,diffusivity,diffusivity_err,tau, tauErr, ParamA, AErr, ParamBeta, BetaErr, paramB, BErr, paramTheta, thetaErr, paramC, CErr] = TGSPhaseAnalysis(posstr,negstr,grat,2,0,blbool,baselinePOS, baselineNEG);
disp('Frequency = ');
disp(freq_final);
disp('Frequency error = ');
disp(freq_error);
disp('Thermal diffusivity = ');
disp(diffusivity);
disp('Thermal diffusivity error = ');
disp(diffusivity_err);

fid1 = fopen( printFile, 'a' );
    naming_str=strcat(str_base,num2str(i));
    fprintf(fid1, string('\n%s %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g'), naming_str, grat, freq_final, freq_error, ParamA, AErr, diffusivity, diffusivity_err, ParamBeta, BetaErr, paramB, BErr, paramTheta, thetaErr, tau(3), tauErr, paramC, CErr);
    fclose(fid1);
%close all % comment this out if you want to keep the plot windows at the end