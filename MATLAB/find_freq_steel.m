function [freq_final,freq_err,flat] = find_freq_steel(pos_file,neg_file,grat)

plotty=0;
plottx=1;

hdr_len=16;

if nargin<5
    rangefrac=0.1;
end

%read in data files for this procedure
pos=dlmread(pos_file,'',hdr_len,0);
neg=dlmread(neg_file,'',hdr_len,0);
% neg(:,1)=neg(:,1)-neg(1,1); %%%%% Comment when cables are same length
% posbas=dlmread('UO2_100-2021_10_01-04.20um-baseline-POS-1.txt','',hdr_len,0);
% negbas=dlmread('UO2_100-2021_10_01-04.20um-baseline-POS-1.txt','',hdr_len+36,0);

% %%%%%%%%%%% UNCOMMENT FOR F'd UP BASELINE
% pos(:,2)=pos(:,2)-posbas(:,2);
% neg(:,2)=neg(:,2)-negbas(:,2);

%normalize each set of data to the zero level before the pump impulse
pos(:,2)=pos(:,2)-mean(pos(1:50,2));
neg(:,2)=neg(:,2)-mean(neg(1:50,2));


[~,time_index]=max(neg(1:1000,2));
time_naught=neg(time_index,1);


%If recorded traces differ in length, fix them to shorter of the two
if length(pos(:,1))>length(neg(:,1))
    pos=pos(1:length(neg(:,1)),:);
%     posbas=posbas(1:length(neg(:,1)),:);
elseif length(neg(:,1))>length(pos(:,1))
    neg=neg(1:length(pos(:,1)),:);
%     negbas=negbas(1:length(pos(:,1)),:);
end

fixed_short=[pos(:,1) pos(:,2)-neg(:,2)];

[~,fix_index]=max(fixed_short(:,2));
fixed_short=fixed_short(fix_index:end,:);

if plotty
    figure()
    plot(fixed_short(:,1),fixed_short(:,2),'r')
    title('this is fixed short');
end

%normalize to initial level before pump pulse
fixed_short(:,2)=fixed_short(:,2)-mean(fixed_short(end-50:end,2));

fixed_short(:,1)=fixed_short(:,1)-fixed_short(1,1); %slide peak back to 0 to fit more easily

%Take a fit to a purely erfc() decay to subtract a decay profile before
%taking fft
q=2*pi/(grat*10^(-6));
LB=[0 0 0];
UB=[10 5*10^-4 0.1];
ST=[0.05 5*10^-5 0];
OPS=fitoptions('Method','NonLinearLeastSquares','Lower',LB,'Upper',UB,'Start',ST);
TYPE=fittype('A.*erfc(q*sqrt(k*(x)))+c;','options',OPS,'problem','q','coefficients',{'A','k','c'});

[f0,~]=fit(fixed_short(:,1),fixed_short(:,2),TYPE,'problem',q);

if plotty
    %plot fit to see how well it matches
    figure()
    plot(fixed_short(:,1),f0(fixed_short(:,1)),'b',fixed_short(:,1),fixed_short(:,2),'r');
    title('this is an erfc(..) fit of fixed short');
end

%this is the thermal decay filtering of the signal recorded vs time. This will clean up the DC end of the power spectrum
flat1=[fixed_short(:,1) fixed_short(:,2)-f0(fixed_short(:,1))];

% Shorten 'flat1' to 'flat' so the FFT samples better data
newlength=ceil(numel(fixed_short(:,1))*rangefrac);
flat=zeros(newlength,2);
flat(:,1)=flat1(1:newlength,1);
flat(:,1)=flat(:,1)*1e8;
flat(:,2)=flat1(1:newlength,2)/max(flat1(1:newlength,2));


if plottx
    figure()
    plot(flat(:,1),flat(:,2),'b-')
    title('this is fixed short minus the erfc(..) fit and shortened');
end

%%%%%%%%%%% Fit 4th order polynomial for noise
STnoise=[0 0 0 0 0];
OPS=fitoptions('Method','NonLinearLeastSquares','Start',STnoise);
TYPEnoise=fittype('A*(x)^4+B*(x)^3+C*(x)^2+D*(x)-E;','options',OPS,'coefficients',{'A','B','C','D','E'});
[fnoise,~]=fit(flat(:,1),flat(:,2),TYPEnoise);
flat(:,2) = flat(:,2) - fnoise(flat(:,1));


%% Fit damped sine wave
STsteel=[.5 0.25 3.5 1.5 0];
% LB=[0 0 3 1 -0.001];
% UB=[1 0.5 4 2 0.001];
% OPS=fitoptions('Method','NonLinearLeastSquares','Lower',LB,'Upper',UB,'Start',STsteel);
OPS=fitoptions('Method','NonLinearLeastSquares','Start',STsteel);
TYPEsteel=fittype('A.*exp(-b*x)*(sin(6.283*c*(x)+d))+e;','options',OPS,'coefficients',{'A','b','c','d','e'});
[fsteel,~]=fit(flat(:,1),flat(:,2),TYPEsteel);

display(fsteel.c)
error=confint(fsteel,0.95);
freq_final = fsteel.c*1e8;
freq_err=(freq_final-error(1,2)*1e8)/2;

% xx=0:1e-3:6e-1;
% yy=spline(flat(:,1),flat(:,2),xx);
% 
% clear flat
% flat=zeros(length(xx),2);
% flat(:,1) = xx;
% flat(:,2) = yy;
% 
% %% Find local maxima
% [psor,lsor] = findpeaks(flat(:,2),flat(:,1),'SortStr','descend');
% % findpeaks(flat(:,2),flat(:,1))
% % text(lsor+.02,psor,num2str((1:numel(psor))'))
% f1 = lsor(2)-0;
% f2 = lsor(3)-lsor(2);
% 
% %% Find local minima
% [psob,lsob] = findpeaks(0-flat(:,2),flat(:,1),'SortStr','descend');
% % findpeaks(0-flat(:,2),flat(:,1))
% % text(lsob+.02,psob,num2str((1:numel(psob))'))
% f3 = lsob(2)-lsob(1);
% f4 = lsob(3)-lsob(2);

% freq_final = 1e8/(0.25*(f1+f2+f3+f4));


if plottx
    figure()
    plot(flat(:,1),flat(:,2),'b-')
    hold on
    plot(flat(:,1),fsteel(flat(:,1)))
    title('Fit of damped sine wave');
    display(fsteel)
end


end
