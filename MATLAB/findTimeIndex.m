function [time_index,end_time] = findTimeIndex(pos_file,neg_file)
%%%%
%%%% Finds the time_index based on the index at which the 2nd derivative of
%%%% fixed_short is a maximum
%%%% 
%%%%
%%%%
%   pos_file: positive phase data pos_file
%   neg_file: negative phase data neg_file

plot_ti=0;   % Show the peaks of raw trace and time index point
plot_derv=0; % Show second derivative and peak fits
hdr_len=16;
off_20ns=23;
off_50ns=19;


%%%%%%%%%%%%% Read in data files and create fixed_short %%%%%%%%%%%%%%%%
pos=dlmread(pos_file,'',hdr_len,0);
neg=dlmread(neg_file,'',hdr_len,0);

pos(:,2)=pos(:,2)-mean(pos(1:50,2));
neg(:,2)=neg(:,2)-mean(neg(1:50,2));

%If recorded traces differ in length, fix them to shorter of the two
if length(pos(:,1))>length(neg(:,1))
    pos=pos(1:length(neg(:,1)),:);
elseif length(pos(:,1))<length(neg(:,1))
    neg=neg(1:length(pos(:,1)),:);
end

%%%%%%%%% Make the Fixed Short Variable
fixed_short=[pos(:,1) pos(:,2)-neg(:,2)];




%%%%%%%%%%%%%%%%%%%%%%%%%%% Methods for finding time_index %%%%%%%%%%%%%%%%

% %%%%%%%%% Find Time Index from above mean(floor) + 10*std(floor)
% flat_index=0;
% hzcond=0;
% toohigh=mean(fixed_short(1:100,2))+10*std(fixed_short(1:100,2));
% while hzcond < toohigh
%     flat_index=flat_index+1;
%     hzcond=fixed_short(flat_index,2);
% end

%%%%%%%%% Find Time Index from Trace MAX
[max_trace,max_time]=max(fixed_short(1:600,2));


%%%%%%%%% Find Time Index from Second Derivative
first = gradient(fixed_short(:,2)) ;   % derivative of fixed_short
second_derv = gradient(first);         % second derivative of fixed_short
[~,derv_index] = max(second_derv(1:max_time));
minprom=5*max(second_derv(1:50));
[pos_val,pos_loc]=findpeaks(second_derv(1:max_time),'MinPeakProminence',minprom);

%%%%%%%%% Make sure the peak it is using is the first after the rise of the
%%%%%%%%% raw trace
if length(pos_loc)>0
    if pos_loc(1)<derv_index
        derv_index=pos_loc(1);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%% Assign time_index and end time values %%%%%%%%%%%%%%%%

%%%%%%%%% Determine indices for different methods and assign end times
time_len=length(pos(:,1))/1000;
if time_len < 5
    end_time=2e-7; %for 20ns base on scope
%     flat_index=flat_index-off_20ns;
%     max_index = max_time-47;
    derv_index = derv_index-off_20ns;
else
    end_time=5e-7; %for 50ns base on scope
%     flat_index=flat_index-off_50ns;
%     max_index = max_time-42;
    derv_index = derv_index-off_50ns;
end

% derv_index = derv_index-off_50ns; %%%% comment out for 20ns timebase
time_index = derv_index;




%%%%%%%%% Plot
if plot_ti
    figure()
    plot(fixed_short(:,2),'k','LineWidth',3)
    xlim([0 max_time+50])
    ylim([-0.005 max_trace+0.005])
    hold on
%     plot([max_index max_index],[-0.005 max_trace+0.005],'g','LineWidth',3)
%     hold on
    plot([derv_index derv_index],[-0.005 max_trace+0.005],'r--','LineWidth',3)
    hold on
%     plot([max_time max_time],[-0.005 max_trace+0.005],'g','LineWidth',3)
    title('Raw Trace with Calculated Time Index');
    ylabel({'Signal Amplitude [mV]'})
    xlabel({'Time Index'})
    
%     txt = ['Time Index_{max} = ',num2str(max_index),'\rightarrow '];
%     text(max_index,(max_trace*0.33),txt,'HorizontalAlignment','right','FontSize',14)
%     txt2 = ['Time Index_{2derv} = ',num2str(derv_index), '\rightarrow '];
%     text(derv_index,(max_trace*0.66),txt2,'HorizontalAlignment','right','FontSize',16)
    txt3 = {['Time Index = ',num2str(time_index)],['End Time = ',num2str(end_time)]};
    text(10,0.9*max_trace,txt3,'FontSize',20)
    
    set(gca,...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',24,...
        'FontName','Helvetica',...
        'LineWidth',3)
    ylabel({'Raw Trace'},...
        'FontUnits','points',...
        'FontSize',24,...
        'FontName','Helvetica')
    xlabel({'Time Index'},...
        'FontUnits','points',...
        'FontSize',24,...
        'FontName','Helvetica')
    
    if plot_derv
        yyaxis right
        hold on
        plot(second_derv*10^4,'b','LineWidth',3)
        hold on
        plot(pos_loc,pos_val*10^4,'o','LineWidth',3)
        ylabel({'2^{nd} Derivative'},...
            'FontUnits','points',...
            'FontSize',24,...
            'FontName','Helvetica')
    end
end

end