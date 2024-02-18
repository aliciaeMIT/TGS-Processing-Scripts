% directory location
pname = 'example_data\tungsten_depth_study'; %ADJUST
%First part of file name and number of spots
str_base='Tungsten_Calibration-2022-05-19-0';  %ADJUST
%Change filename modifying number here
spots=5:20; %ADJUST
%Baseline handling if no baseline arguments are provided, the code still works - it just doesn't do any baseline correction
baselineBool = 1; % 1 for "yes do baseline subtraction", 0 for "no don't do that" 

% Calibrated grating spacings - this script ran through a bunch of nominal
% grating spacings for tungsten to check depth effects. gratArr was changed
% each time.
gratArr = [5.2250,5.2968, 5.0289]; %, 5.8382, 5.8371, 5.8309, 6.4742, 6.4832, 6.4587, 7.0421, 7.0282, 7.0426, 7.6541, 7.6670, 7.6535, 9.5139, 9.5135, 9.5195]; %in um %ADJUST

% Initialize all outputs
freq=0;
freq_err=0;
speed=0;
diff=0;
diff_err=0;
tau=zeros(3,length(spots));
tau_err = 0;
A = 0;
AErr = 0;
beta = 0;
betaErr = 0;
B = 0;
BErr = 0;
theta = 0;
thetaErr = 0;
C = 0;
CErr = 0;
nominalSpacing = 0;
formatSpec = '%.1f';
%Set up the output files with their headers, clears the file at beginning.
printFile = pname + string('\') + str_base + string('grating_depth_study_3.7um.txt');
fid1 = fopen( printFile, 'wt' );
fprintf(fid1,'%s', 'grating_value[um] SAW_freq[Hz] SAW_freq_error[Hz] A[Wm^-2] A_err[Wm^-2] alpha[m^2s^-1] alpha_err[m2s-1] beta[s^0.5] beta_err[s^0.5] B[Wm^-2] B_err[Wm^-2] theta theta_err tau[s] tau_err[s] C[Wm^-2] C_err[Wm^-2]');
fclose(fid1);

%Run through each of the files, sweeping across gratings and peak-skips
for i=1:length(gratArr)
    
    if i<4
        nominalSpacing = 3.7;
    elseif i<7
        nominalSpacing = 5.8;
    elseif i<10
        nominalSpacing = 6.4;
    elseif i<13
        nominalSpacing = 7.0;
    elseif i<16
        nominalSpacing = 7.6;
    else
        nominalSpacing = 9.4;
    end
    
    close all
    pos_str=strcat(pname,string('\'),str_base,num2str(nominalSpacing,formatSpec),'0um-spot',num2str(mod(i-1,3)),'-POS-1.txt');
    neg_str=strcat(pname,string('\'),str_base,num2str(nominalSpacing,formatSpec),'0um-spot',num2str(mod(i-1,3)),'-NEG-1.txt');
    POSbaselineStr = strcat(pname,'\',str_base,num2str(nominalSpacing,formatSpec),'0um-baseline-POS-1.txt'); %ADJUST
    NEGbaselineStr = strcat(pname,'\',str_base,num2str(nominalSpacing,formatSpec),'0um-baseline-NEG-1.txt'); %ADJUST
    
    disp(['current run is: ', pos_str, '   with grating: ', num2str(gratArr(i)),' um']);
    disp(['current baseline is: ', POSbaselineStr]);
    [freq,freq_err,speed,diff,diff_err,tau(:,i),tauErr, A, AErr, beta, betaErr, B, BErr, theta, thetaErr, C, CErr]=TGSPhaseAnalysis(pos_str,neg_str,gratArr(i),2,0,1, POSbaselineStr, NEGbaselineStr);
    
    fid1 = fopen( printFile, 'a' );
    fprintf(fid1, string('\n%.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g'), gratArr(i), freq, freq_err, A, AErr, diff, diff_err, beta, betaErr, B, BErr, theta, thetaErr, tau(3,i), tauErr, C, CErr);
    fclose(fid1);
    
end

disp('F L A W L E S S');
