function [freq_final,freq_error,speed,diffusivity,diffusivity_err,tau] = TGSPhaseAnalysis(pos_file,neg_file,grat,start_phase,two_mode,end_time,time_index,overlay1,overlay2)

%     Function to determine thermal diffusivity from phase grating TGS data
%   Data is saved in two files, positive (with one heterodyne phsae) and
%       negative (with another), must provide both files
%   pos_file: positive phase TGS data file
%   neg_file: negative phae TGS data file
%   grat: calibrated grating spacing in um
%   start_phase: provide integer between 1 and 4 to pick the null-point start
%               from which fit will begin
%   two_mode: a boolean value, default value is 0 if only one acoustic mode
%               is present in the measurement. If two are present, provide 1 for both
%               frequencies and speeds to be output.
%   end_time: a shortened fit end time if you do not want ot fit the whole profile.
%               If argument not given, will be set to default for 200ns data collection window

%%%%%%%%%%%%%%%%%%%%%%%%%
%   Extra parameter: time_index: Define where the actual fit starts.
%   Make it negative if you want the script to find it. Specify it if not.
%   Experimental, working on a physical way to define this.
%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%Write this to include a sine variation in the fit by default, but to
%%%%start the fits from a fixed null point, not time, relative to
%%%% the initial SAW maximum.

%Settings for various plotting and output options to be set by boolean arguments
find_max=0;
plotty=0;
plot_trace=0;
plot_psd=1;
plot_final=1;
print_final_fit=0;
two_detectors=1;
q=2*pi/(grat*10^(-6));
tstep=5e-11; %Set by scope used for data collection
no_pre_calc=0;
amp_grat=0;

%if on, record the amplitude ratio of A/B in final fit at the end of tau
%vector (so it is four elements instead of 3)
record_amp_ratio=0;

spot_size=140e-6; %measured pump spot size of 140um - MIT
% spot_size=170e-6; %measured pump spot size of 170um - Sandia

derivative=0;

%%%%%%%%%%%%%%%%%%%
%Tau selection block. Allows you to select the final model using a variety
%of time constant schemes
%%%%%%%%%%%%%%%%%%%
%if fitting tau need to select one of three options for start point
tau_final_fit=1;
start_constant=1;
start_lorentz=0;
start_walkoff=0;
%if not fitting tau, need to select one of two constant schemes
fixed_tau=0;
fixed_lorentz=0;
fixed_walkoff=0;
%check schemes to make sure that the code will always run. default to
%fitting with constant values
if tau_final_fit==0 && fixed_tau==0
    tau_final_fit=1;
elseif tau_final_fit==1 && fixed_tau==1
    tau_final_fit=0;
end
if tau_final_fit==1 && start_constant==0 && start_lorentz==0 && start_walkoff==0
    start_constant=1;
end
if fixed_tau==1 && fixed_lorentz==0 && fixed_walkoff==0
    fixed_walkoff=1;
end
%%%%%%%%%%%%%%%%%%%

%How far from guessed values for diffusivity and beta do you vary in the
%end, d for diffusivity b for beta. Diffusivity is most important
percent_range_d=0.45;
percent_range_b=2.25;
if tau_final_fit && ~start_constant
    if start_walkoff
        percent_range_t=2.5;
    elseif start_lorentz
        percent_range_t=0.9;
    end
end

%Output tau will have multiple values, initialize variable.
%tau(1) is the lorentzian decay time
%tau(2) is the walk-off time
%tau(3) is the final value either used in or optimized by fit
tau=zeros(1,3);

if nargin<5
    two_mode=0;
end

%Difference in file write format based on newer or older acquisition. hdr_len should be 16 for the Ge dataset
if two_detectors
    hdr_len=16;
else
    hdr_len=15;
end

%Generate filtered power spectrum using TGS_phase_fft() and find the peak frequency from that profile using lorentzian_peak()
[fft]=TGS_phase_fft(pos_file,neg_file,grat,1);
[freq_final,freq_error,peak_fwhm,tau_lorentz]=lorentzian_peak_fit(fft,two_mode,plot_psd);
speed=freq_final*grat*10^(-6);
tau(1)=tau_lorentz(1);

walk_off=spot_size/(speed(1)*2); %this is the walk off time in sec
tau(2)=walk_off;

peak_freq=freq_final(1);

%read in data files for this procedure
pos=dlmread(pos_file,'',hdr_len,0);
neg=dlmread(neg_file,'',hdr_len,0);

%sometimes written data is off by one time step at the end, chop that off if they do not match
if length(pos(:,1))>length(neg(:,1))
    pos=pos(1:length(neg(:,1)),:);
elseif length(neg(:,1))>length(neg(:,1))
    neg=neg(1:length(pos(:,1)),:);
end

%normalize each set of data to the zero level before the pump impulse
pos(:,2)=pos(:,2)-mean(pos(1:50,2));
neg(:,2)=neg(:,2)-mean(neg(1:50,2));

%%%%%Time indexing block, important to keep track of%%%%%%%%
%%%%%%%%%%%%%%%%%% now it's not!!! %%%%%%%%%%%%%%%%%%%%%%%%%
 if time_index<0
    [time_index,end_time]=findTimeIndex(pos_file,neg_file);
    display(['Time Index is: ',num2str(time_index)])
 end

% time_index=186;
% time_index=235;

time_naught=neg(time_index,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end_index=floor(end_time/tstep);

%re-normalize data to end signal decayed state if grating decays entirely during collection window, if not, do not re-normalize
if grat<8
    base_index=floor(length(pos(:,2))/5);
    long_base=mean((pos(1:time_index,2)-neg(1:time_index,2)));
%     long_base=mean((pos(end-base_index:end,2)-neg(end-base_index:end,2)));
else
    long_base=0;
end

if plot_trace
    figure()
    plot(neg(:,1)*10^9,(pos(:,2)-neg(:,2)-long_base)*10^3,'-','Color',[0 0 0.75],'LineWidth',1.25)
    hold on
%     plot([neg(1,1) neg(end,1)]*10^9,[0 0],'k--','LineWidth',1.5)
    xlim([0 (end_time/2)*10^9])
    set(gca,...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',16,...
        'FontName','Helvetica',...
        'LineWidth',1.25)
    ylabel({'Amplitude [mV]'},...
        'FontUnits','points',...
        'FontSize',20,...
        'FontName','Helvetica')
    xlabel({'Time [ns]'},...
        'FontUnits','points',...
        'FontSize',20,...
        'FontName','Helvetica')
    
    txt3 = {['Time Index = ',num2str(time_index)],['End Time = ',num2str(end_time)]};
    text(neg(end/2-1000,1)*10^9,(max(pos(:,2)-neg(:,2)-long_base))*.75*10^3,txt3,'FontSize',24)
    saveas(gcf,"TGS_Trace.png")
end

if start_phase==0
    diffusivity=0;
    diffusivity_err=0;
else
    
    if amp_grat==0

        fixed_short=[pos(time_index:end_index,1)-time_naught pos(time_index:end_index,2)-neg(time_index:end_index,2)-long_base];
        if derivative
            der_len=length(fixed_short(:,1))-1;
            fixed_derivative=zeros(1,der_len);
            for jj=1:der_len
                fixed_derivative(jj)=(fixed_short(jj+1,2)-fixed_short(jj,2))/tstep;
            end
            figure()
            plot(fixed_short(1:der_len,1),fixed_derivative,'k-')
            title('This is the derivative of fixed short')
        end

        %if you don't want to automatically find the peak t_0 of the profile, setting find_max to true above will allow
        %you to select a region on the plot within with to search. Useful if there are initial transients.
        if find_max || plotty
            figure()
            plot(fixed_short(:,1),fixed_short(:,2),'k-')
            xlim([0 1.5*10^-7])
            title('this is fixed short');
            if find_max
                hold on
                [x_cord,~]=ginput(2);
                neg_x_cord=x_cord(1);
                pos_x_cord=x_cord(2);
                pos_x_ind=floor(pos_x_cord/tstep);
                neg_x_ind=floor(neg_x_cord/tstep);
                [~,max_index]=max(fixed_short(neg_x_ind:pos_x_ind,2));
                time_max=fixed_short(max_index+neg_x_ind,1);
                close(gcf)
            end
        end

        %Otherwise, find t_0 from the profile directly
        if ~find_max
            time_offset_index=20; %was 20 before
            [~,max_index]=max(fixed_short(time_offset_index:end,2));
            max_index=max_index+time_offset_index-1;
            time_max=fixed_short(max_index,1);

            start_time_phase=find_start_phase(fixed_short(max_index:end,1),fixed_short(max_index:end,2),start_phase,grat);
%             start_time_phase=find_start_phase(fixed_short(:,1),fixed_short(:,2),start_phase,grat);
            start_index_master=round(start_time_phase/tstep)+1;
            start_time_master=fixed_short(start_index_master,1);

            %Fitting parameters for initial naive fit
            LB=[0 0];
            UB=[1 5*10^-4];   % Increased to account for silver-alloys, which are REALLY HIGH thermal diffusivity
            ST=[.05 5*10^-6];

            OPS=fitoptions('Method','NonLinearLeastSquares','Lower',LB,'Upper',UB,'Start',ST);
            TYPE=fittype('A.*erfc(q*sqrt(k*(x+time_max)))','options',OPS,'problem',{'q','time_max'},'coefficients',{'A','k'});
            [f0,gof]=fit(fixed_short(:,1),fixed_short(:,2),TYPE,'problem',{q,time_max});

            diffusivity=f0.k;
            error=confint(f0,0.95);
            %factor of 2 makes the 1 sigma confidence interval come out
            diffusivity_err=[diffusivity-error(1,2) error(2,2)-diffusivity]/2;

            if plotty
                figure()
                plot(fixed_short(:,1),fixed_short(:,2),fixed_short(:,1),f0(fixed_short(:,1)))
                hold on
                title('First naive fit')
            end

            %We'll call the parameter beta the ratio of the amplitudes of the
            %displacement versus temperature grating. beta should be a small number.

            for jj=1:10

                beta=q*sqrt(diffusivity/pi)*(q^2*diffusivity+1/(2*time_max))^(-1);

                start_time=start_time_master;
                start_index=start_index_master;

                %Conduct initial parameter estimation without using an sin(x) contribution to the fit

                LB1=[0 0];
                UB1=[1 5*10^-4];   % Increased to account for silver-alloys, which are REALLY HIGH thermal diffusivity
                ST1=[.05 10^-5];

                OPS1=fitoptions('Method','NonLinearLeastSquares','Lower',LB1,'Upper',UB1,'Start',ST1);
                TYPE1=fittype('A.*(erfc(q*sqrt(k*(x+start_time)))-beta*exp(-q^2*k*(x+start_time))./sqrt((x+start_time)))','options',OPS1,'problem',{'q','beta','start_time'},'coefficients',{'A','k'});
                [f1,gof]=fit(fixed_short(start_index:end,1),fixed_short(start_index:end,2),TYPE1,'problem',{q,beta,start_time});

                diffusivity=f1.k;
                error=confint(f1,0.95);
                %factor of 2 makes the 1 sigma confidence interval come out
                diffusivity_err=[diffusivity-error(1,2) error(2,2)-diffusivity]/2;

                if plotty
                    figure()
                    plot(fixed_short(:,1),fixed_short(:,2))
                    hold on
                    title(strcat('Fit number ',num2str(jj+1),' - fixed beta'))
                end

            end

            %If you've elected not to pre-compute, provide hard initial guesses for diffusivity and beta based on
            %the bulk value of Ge diffusivity. Set ranges for final fit.
            if no_pre_calc
                diffusivity=0.3636*10^-4;
                beta=2e-5;
                low_bound=[1e-5 0];
                up_bound=[1e-3 1e-4];
            else
                low_bound=[diffusivity*(1-percent_range_d) beta*(1-percent_range_b)];
                up_bound=[diffusivity*(1+percent_range_d) beta*(1+percent_range_b)];
                if percent_range_d>1
                    low_bound(1)=0;
                end
                if percent_range_b>1
                    low_bound(2)=0;
                end
            end

            start_time2=start_time_master;
            start_index2=start_index_master;

            if tau_final_fit
                if start_constant
                    low_t=5e-9;
                    up_t=7e-7;
                    start_tau=1e-8;
                elseif start_lorentz
                    low_t=tau(1)*(1-percent_range_t);
                    up_t=tau(1)*(1+percent_range_t);
                    if percent_range_t>1
                        low_t(1)=0;
                    end
                    start_tau=tau(1);
                    LB2=[0 low_bound(1) low_bound(2) 0 -2*pi low_t -5e-3];
                    UB2=[1 up_bound(1) up_bound(2) 10 2*pi up_t 5e-3];
                    ST2=[.05 diffusivity beta 0.05 0 tau(1) 0];
                elseif start_walkoff
                    low_t=tau(2)*(1-percent_range_t);
                    up_t=tau(2)*(1+percent_range_t);
                    if percent_range_t>1
                        low_t=0;
                    end
                    start_tau=tau(2);
                end

    %            LB2=[0 low_bound(1) low_bound(2) 0 -2*pi low_t -5e-3];
    %            UB2=[1 up_bound(1) up_bound(2) 10 2*pi up_t 5e-3];
                LB2=[0 low_bound(1) low_bound(2) 0 -2*pi low_t -5e-3];
                UB2=[1 up_bound(1) up_bound(2) 10 2*pi up_t 5e-3];
                ST2=[.05 diffusivity beta 0.05 0 start_tau 0];

                OPS2=fitoptions('Method','NonLinearLeastSquares','Lower',LB2,'Upper',UB2,'Start',ST2);
                TYPE2=fittype('A.*(erfc(q*sqrt(k*(x+start_time)))-beta*exp(-q^2*k*(x+start_time))./sqrt((x+start_time)))+B.*sin(2*pi*(peak_freq)*(x+start_time)+p)*exp(-(x+start_time)/t)+D','options',OPS2,'problem',{'q','start_time','peak_freq'},'coefficients',{'A','k','beta','B','p','t','D'});
                [f2,gof]=fit(fixed_short(start_index2:end,1),fixed_short(start_index2:end,2),TYPE2,'problem',{q,start_time2,peak_freq});

            elseif fixed_tau
                if fixed_lorentz
                    pp_tau=tau(1);
                elseif fixed_walkoff
                    pp_tau=tau(2);
                end
                LB2=[0 low_bound(1) low_bound(2) 0 -2*pi -5e-3];
                UB2=[1 up_bound(1) up_bound(2) 10 2*pi 5e-3];
                ST2=[.05 diffusivity beta 0.05 0 0];

                OPS2=fitoptions('Method','NonLinearLeastSquares','Lower',LB2,'Upper',UB2,'Start',ST2);
                TYPE2=fittype('A.*(erfc(q*sqrt(k*(x+start_time)))-beta*exp(-q^2*k*(x+start_time))./sqrt((x+start_time)))+B.*sin(2*pi*(peak_freq)*(x+start_time)+p)*exp(-(x+start_time)/t)+D','options',OPS2,'problem',{'q','start_time','peak_freq','t'},'coefficients',{'A','k','beta','B','p','D'});
                [f2,gof]=fit(fixed_short(start_index2:end,1),fixed_short(start_index2:end,2),TYPE2,'problem',{q,start_time2,peak_freq,pp_tau});
%                 [f2,gof]=fit(fixed_short((start_index2+time_index):end,1),fixed_short((start_index2+time_index):end,2),TYPE2,'problem',{q,start_time2,peak_freq,pp_tau});
            end
%             display(['Start time: ', num2str(start_time2)])
%             display(['Start index: ', num2str(start_index2)])
%             display(['Start index + Time_index: ', num2str(start_index2+time_index)])

            if print_final_fit
                display(f2)
            end

            diffusivity=f2.k;

            error=confint(f2,0.95);
            %factor of 2 makes the 1 sigma confidence interval come out
            diffusivity_err=(diffusivity-error(1,2))/2;

            %final fit (on constant provided) version of tau
            tau(3)=f2.t;

            if record_amp_ratio
                gamma=f2.A/f2.B;
                tau=[tau gamma];
            end

            %first checks that diffusivity has not pegged to fit bounds, second checks that beta
            %has not pegged. If either do, it is a bad fit
            bad_alpha=isnan(diffusivity_err(1));
            bad_beta=isnan(error(1,3));
            if bad_alpha && bad_beta
                display(strcat('Bad fit for: ',pos_file,'~re: tau (likely)'))
            elseif bad_alpha && ~bad_beta
                display(strcat('Bad fit for: ',pos_file,'~re: alpha (thermal diffusivity)'))
            elseif ~bad_alpha && bad_beta
                display(strcat('Bad fit for: ',pos_file,'~re: beta (ratio of reflectivity to displacement)'))
            end

            if plot_final
                %Plotting factor for generation of traces in Figure 6
                amp_factor=1;
                %%%%%%%%%%%%%%%
                %Block to reconstruct the best-fit model without the sinusoidal
                %contribution, for comparison
                if tau_final_fit
                    f_remove_sine=cfit(TYPE2,f2.A,f2.k,f2.beta,f2.B,f2.p,0,f2.D,q,start_time2,peak_freq);
                elseif fixed_tau
                    f_remove_sine=cfit(TYPE2,f2.A,f2.k,f2.beta,f2.B,f2.p,f2.D,q,start_time2,peak_freq,0);
                end
                %%%%%%%%%%%%%%%
                figure()
                plot((neg(:,1)-time_naught)*10^9,(pos(:,2)-neg(:,2)-long_base)*10^3/amp_factor,'k-','LineWidth',5,'DisplayName','Raw TGS Trace')
                hold on
                %plot vertical line at start time
                plot([fixed_short(start_index2,1) fixed_short(start_index2,1)]*10^9,ylim,'g--','LineWidth',5,'DisplayName','Start Index')
                hold on
                plot([neg(time_index,1)-time_naught neg(time_index,1)-time_naught]*10^9,ylim,'c--','LineWidth',5,'DisplayName','Time Index')
                hold on
                plot(fixed_short(start_index2:end,1)*10^9,(f2(fixed_short(start_index2:end,1)))*10^3/amp_factor,'r--','LineWidth',5,'DisplayName','Full Functional Fit')
                hold on
                plot(fixed_short(start_index2:end,1)*10^9,(f_remove_sine(fixed_short(start_index2:end,1)))*10^3/amp_factor,'-','Color',[0 0 0.75],'LineWidth',5,'DisplayName','Thermal Fit')
                hold on
                xlim([-5 end_time*10^9])
%                 xlim([-5 50])
%                 ylim([10 30])
                set(gcf,'Position',[0 0 1920 1080])
		if nargin>7	%% Overlays come last, only place if specified. This stops throwing errors when no grapn overlays are given.
	             annotation('textbox',[0.02 0.01 0.5 0.03],'String',overlay2,'FontSize',25,'FontName','Arial','FontWeight','bold','LineStyle','none')
	             annotation('textbox',[0.6 0.01 0.35 0.03],'String',overlay1,'FontSize',25,'FontName','Arial','FontWeight','bold','LineStyle','none')
		end
                hold on
                set(gca,...
                    'FontUnits','points',...
                    'FontWeight','normal',...
                    'FontSize',30,...
                    'FontName','Helvetica',...
                    'LineWidth',5)
                ylabel({'Signal Amplitude [mV]'},...
                    'FontUnits','points',...
                    'FontSize',40,...
                    'FontName','Helvetica')
                xlabel({'Time [ns]'},...
                    'FontUnits','points',...
                    'FontSize',40,...
                    'FontName','Helvetica')
            legend('Location','best')

    % Display the already-made FFT as an inset image

            axes('pos',[.48 .48 .5 .5])
            imshow('TGS_FFT.png')
            saveas(gcf,strcat(pos_file,"_TGS_Final_Fit.png"))
            end

        end
%         display(['Start Index is: ',num2str(start_index2+time_index)])

    else
        fixed_short=[pos(time_index:end_index,1)-time_naught pos(time_index:end_index,2)-neg(time_index:end_index,2)-long_base];
    %     start_time_phase=find_start_phase(fixed_short(:,1),fixed_short(:,2),start_phase,grat);
    %     start_index=round(start_time_phase/tstep)+1;
    %     start_time=fixed_short(start_index,1);

        LBamp=[0 0 0];
        UBamp=[1 5*10^-3 1];
        STamp=[.1 5*10^-5 .1];
        start_time=time_index*tstep;

        OPSamp=fitoptions('Method','NonLinearLeastSquares','Lower',LBamp,'Upper',UBamp,'Start',STamp);
        TYPEamp=fittype('(A./sqrt((x+start_time)))*exp(-q^2*k*(x+start_time))+D','options',OPSamp,'problem',{'q','start_time'},'coefficients',{'A','k','D'});
        [famp,gof]=fit(fixed_short(time_index:end,1),fixed_short(time_index:end,2),TYPEamp,'problem',{q,start_time});

        diffusivity=famp.k;
        amp_factor=1;

        if plot_final
            
            figure()
            plot((neg(:,1)-time_naught)*10^9,(pos(:,2)-neg(:,2)-long_base)*10^3/amp_factor,'k-','LineWidth',5,'DisplayName','Raw TGS Trace')
            hold on
            plot(fixed_short(time_index:end,1)*10^9,(famp(fixed_short(time_index:end,1)))*10^3/amp_factor,'r--','LineWidth',5,'DisplayName','Full Functional Fit')
            hold on
            xline(fixed_short(time_index,1)*10^9,'-b','LineWidth',5,'DisplayName','time index');
            xlim([-5 end_time*10^9])
            set(gcf,'Position',[0 0 1920 1080])
	    if nargin>7     %% Overlays come last, only place if specified. This stops throwing errors when no grapn overlays are given.
                annotation('textbox',[0.02 0.01 0.5 0.03],'String',overlay2,'FontSize',25,'FontName','Arial','FontWeight','bold','LineStyle','none')
                annotation('textbox',[0.6 0.01 0.35 0.03],'String',overlay1,'FontSize',25,'FontName','Arial','FontWeight','bold','LineStyle','none')
	    end
            hold on
            set(gca,...
                'FontUnits','points',...
                'FontWeight','normal',...
                'FontSize',30,...
                'FontName','Helvetica',...
                'LineWidth',5)
            ylabel({'Signal Amplitude [mV]'},...
                'FontUnits','points',...
                'FontSize',40,...
                'FontName','Helvetica')
            xlabel({'Time [ns]'},...
                'FontUnits','points',...
                'FontSize',40,...
                'FontName','Helvetica')
            legend('Location','northeast')
        end
        
        diffusivity=0;
        diffusivity_err=0;
        
        if print_final_fit
                display(famp)
        end
    end
end

 fileID = fopen('Analysis/Compiled-Analysis.csv','a');
 fprintf(fileID, '%s',pos_file);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% % Uncomment this to include the time_index for easier looping/testing purposes
% %

 fprintf(fileID, ',%E', time_index);

% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%
% 
 fprintf(fileID, ',%E', freq_final);
 fprintf(fileID, ',%E', freq_error);
 fprintf(fileID, ',%E', diffusivity);
 fprintf(fileID, ',%E', diffusivity_err);
 fprintf(fileID, ',%E,%E,%E\n', tau);
 fclose(fileID);
