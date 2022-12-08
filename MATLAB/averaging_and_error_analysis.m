%Inputs go here
calibration_filename_base = 'C:\Users\apcwy\Dropbox (The University of Manchester)\Physics PhD\MIT_2022\TGS_Data\old_data_2019\Cr_single_crystal\Tungsten-Calibration-2019-08-22-2019-08-22-06.40um-2ndLastIrrSpot0spot01';
pname = 'C:\Users\apcwy\Dropbox (The University of Manchester)\Physics PhD\MIT_2022\TGS_Data\old_data_2019\Cr_single_crystal'; %ADJUST
outputFileName = 'Cr_damaged_TGS_Outputs.txt';
%First part of file name and number of spots
str_base='Cr-single-crystal-2019-08-22-06.40um-5.487909e-5C40degSpot0'; 
spots=1:10;
dir=cd(pname);
cd(dir)
%Baseline handling - you always have to supply some form of baseline anyway for this to work, even if it is not used!
baselineBool = 0; % 1 for "yes do baseline subtraction", 0 for "no don't do that" 
POSbaselineStr = pname + string('\Tungsten_Irradiation_5-2022-01-13-06.40um-baseline') + string('-POS-1.txt'); %ADJUST
NEGbaselineStr = pname + string('\Tungsten_Irradiation_5-2022-01-13-06.40um-baseline') + string('-NEG-1.txt'); %ADJUST
initial_grating_guess = 6.4; %in um, use the nominal value here

%Figure out true grating size from calibration run
%initialise calibration outputs
freq_cal = 0;
freq_err_cal = 0;
speed_cal = 0;
diff_cal = 0;
diff_err_cal = 0;
tau_cal=zeros(3,1);
tau_err_cal = 0;
A_cal = 0;
AErr_cal = 0;
beta_cal = 0;
betaErr_cal = 0;
B_cal = 0;
BErr_cal = 0;
theta_cal = 0;
thetaErr_cal = 0;
C_cal = 0;
CErr_cal = 0;
% Calibrated grating spacing
pos_str_cal=strcat(calibration_filename_base,'-POS-1.txt');
neg_str_cal=strcat(calibration_filename_base,'-NEG-1.txt');
[freq_cal,freq_err_cal,speed_cal,diff_cal,diff_err_cal,tau_cal(:,1),tauErr_cal,A_cal,AErr_cal,beta_cal,betaErr_cal,B_cal,BErr_cal,theta_cal,thetaErr_cal,C_cal,CErr_cal]=TGSPhaseAnalysis(pos_str_cal,neg_str_cal,initial_grating_guess,2,0,'overlay1','overlay2', 0, POSbaselineStr, NEGbaselineStr);
grat=1E6*2665.9/freq_cal; %in um, speed of sound in tungsten is 2665.9 m/s
disp(['calculated grating size: ',num2str(grat,8),' um']);

% Initialize all outputs
freq=zeros(1,length(spots));
freq_err=zeros(1,length(spots));
speed=zeros(1,length(spots));
diff=zeros(1,length(spots));
diff_err=zeros(1,length(spots));
tau=zeros(3,length(spots));
tau_err = zeros(1,length(spots));
A = zeros(1,length(spots));
AErr = zeros(1,length(spots));
beta = zeros(1,length(spots));
betaErr = zeros(1,length(spots));
B = zeros(1,length(spots));
BErr = zeros(1,length(spots));
theta = zeros(1,length(spots));
thetaErr = zeros(1,length(spots));
C = zeros(1,length(spots));
CErr = zeros(1,length(spots));

%Calculate fits for every data set and save output in arrays
adder = 0;
for jj=1:length(spots)
%     if jj>7
%         adder = 1;
%     else
%         adder = 0;
%     end
    if jj==10
        str_base=str_base(1:end-1);
    end
    close all
    pos_str=strcat(pname,'\',str_base,num2str(spots(jj)+adder),'-POS-1.txt');
    neg_str=strcat(pname,'\',str_base,num2str(spots(jj)+adder),'-NEG-1.txt');
    disp(['current run is: ', pos_str]);
    [freq(jj),freq_err(jj),speed(jj),diff(jj),diff_err(jj),tau(:,jj),tauErr(jj), A(jj), AErr(jj), beta(jj), betaErr(jj), B(jj), BErr(jj), theta(jj), thetaErr(jj), C(jj), CErr(jj)]=TGSPhaseAnalysis(pos_str,neg_str,grat,2,0,'overlay1','overlay2', baselineBool, POSbaselineStr, NEGbaselineStr);
end

printFile1 = pname + string('\') + outputFileName;
%outputting to txt
fid1 = fopen( printFile1, 'a' );
%Calculate mean and combined experimental+statistical error for these data.
[freq_mean, freq_error_final] = meanAndError(freq,freq_err);
[A_mean, A_error_final] =  meanAndError(A, AErr);
[diff_mean, diff_error_final] = meanAndError(diff, diff_err);
[beta_mean, beta_error_final] = meanAndError(beta, betaErr);
[B_mean, B_error_final] = meanAndError(B, BErr);
[theta_mean, theta_error_final] = meanAndError(theta, thetaErr);
[tau_mean, tau_error_final] = meanAndError(tau(3), tauErr);
[C_mean, C_error_final] = meanAndError(C, CErr);
fprintf(fid1, string('\n%s %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g'), str_base, grat, freq_mean, freq_error_final, A_mean, A_error_final, diff_mean, diff_error_final, beta_mean, beta_error_final, B_mean, B_error_final, theta_mean, theta_error_final, tau_mean, tau_error_final, C_mean, C_error_final);
fclose(fid1);

close all
disp('finished');
