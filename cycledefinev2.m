
% A more straightforward algorithm for identifying the starting points of 
% each cycle in an oscillating "sawtooth-like" signal in a robust manner. 
% Expecially good for integrated bursting signals, such as those recovered
% from respiratory EMG recordings. Best performance when the cycles are 
% fairly consistent, or when a large signal-noise ratio.
% Assumes sampling rate of 1000Hz

%Input: 
% signal = oscillating signal
% cutoff = filter cutoff, typically set at 1. Can be raised to suit the 
    % rhythm of signal and the amount of noise present 

%Written by Mitchell Prostebby, 2023


function [ind,pkind,pkamp] = cycledefinev2(signal,cutoff)

sigd1 = diff(signal);
[b,a] = butter(1,cutoff./500);  
sigd1 = filtfilt(b,a,sigd1);
sigd2 = diff(sigd1);

%find peaks
[pkamp,pkind,~] = findpeaks(signal,'MinPeakDistance',5);
%thresh 1
pkind(pkamp<0) = []; pkamp(pkamp<0) = [];
pkamp(isnan(pkind)) = []; pkind(isnan(pkind)) = [];

nBinDET = 2;
pkthresh1 = figure;
disp('Threshold the peak amp')
histogram(log(pkamp),round(length(pkamp)./nBinDET,0)); 
xlabel('Signal Peak prominences (log) {Discard below selection]'); [x,~] = ginput(1);
pkind(pkamp<exp(x)) = []; pkamp(pkamp<exp(x)) = []; 
close(pkthresh1)
% pkamp = signal(pkind);
%remove first cycle
pkind(1) = []; pkamp(1) = [];

%finding cycle start
sigd2(signal<0) = 0;
sigd2(sigd1<0) = 0;
[d_pkamp,ind,~,d_pkprom] = findpeaks(sigd2,'MinPeakDistance',5);
ind(d_pkamp<0) = []; d_pkprom(d_pkamp<0) = []; d_pkamp(d_pkamp<0) = []; 

t = 1;
while t==1
    disp('Threshold the amp of the derivative');
    c = figure; histogram(log(d_pkamp),length(d_pkamp).*2); xlabel('amp of signal derivative [Discard below selection]');
    [x,~] = ginput(1);
    % if ~isempty(d_pkamp<exp(x))
    %     t = true;
    % else
    %     t = false;
    % end
    ind(d_pkamp<exp(x)) = []; d_pkprom(d_pkamp<exp(x)) = []; d_pkamp(d_pkamp<exp(x)) = []; 
    ind = ind+1; %resetting due to diff function. 
    close(c)
    disp('Threshold the prominence of the derivative');
    c = figure; histogram(log(d_pkprom),length(d_pkprom).*2); xlabel('Prominence of signal derivative [Discard below selection]');
    [x,~] = ginput(1);
    ind(d_pkprom<exp(x)) = []; d_pkamp(d_pkprom<exp(x)) = []; d_pkprom(d_pkprom<exp(x)) = [];
    ind = ind+1; %resetting due to diff function.
    close(c)
    f = figure; plot(sigd1); hold on; scatter(ind,zeros(length(ind),1),'r');
    t = input('Rethreshold?(1=y,2=no)=');
    close(f)
end


ind(1) = []; d_pkamp(1) = []; d_pkprom(1) = [];
% correction = figure;
% plot(signal(ind(1)-1000:ind(1)+1000)); title('Choose cycle start correction point');
% [c,~] = ginput(1); c = (c-1000); close(correction);
% ind = round(ind+c,0);

%% checkin each peak
if pkind(1)<ind(1)
    pkamp(1) = []; pkind(1) = [];
end
y_zoom = [min(signal) max(signal)];
r_ind = [];
for i = (2:length(pkind))-1
    firstpeakstart = find(ind<pkind(i),1,'last');
    secondpeakstart = find(ind<pkind(i+1),1,'last');
    if firstpeakstart==secondpeakstart

        chek = figure;
        plot(signal(pkind(i-1):pkind(i+2)),'k'); hold on
        scatter(pkind(i-1:i+2)-pkind(i-1),ones(4,1).*mean(signal),'g');
        scatter(ind(firstpeakstart)-pkind(i-1),ones(1,1).*mean(signal),'r')
        edit_value = input('add (1) or remove (2)?');
        if edit_value == 1
            disp('Choose between next 2 peaks a peak cycle starting point');
            [d,~] = ginput(1); d = d+pkind(i-1); %absolute index
            ind(firstpeakstart+2:end+1) = ind(firstpeakstart+1:end);
            ind(firstpeakstart+1) = round(d,0); %Inserts a new index after 

        elseif edit_value == 2
            r = input('Which one to remove?(1/2)');
            % r_ind(end+1) = pkind(i+r-1);
            r_ind(end+1) = i+r-1;
        end
        close(chek);
    end 
end
pkamp(r_ind) = []; pkind(r_ind) = []; 


%Checking the cycle starts: each cycle start that does not have a
%corresponding peak, one of them should be deleted
for i = 1:(length(ind)-1)
   firstcyclepeak = find(pkind>ind(i),1,'first');
   secondcyclepeak = find(pkind>ind(i+1),1,'first');
   if firstcyclepeak==secondcyclepeak
      disp('Removing false cycle starts');
      % figure(main); xlim([(strt_ind(i+1)./1000)-1 (strt_ind(i+1)/1000)+1]);
      chek = figure;
      plot(signal(ind(i-1):ind(i+2)),'k'); hold on
      scatter(ind(i-1:i+2)-ind(i-1),ones(4,1).*mean(signal),'r');
      scatter(pkind(firstcyclepeak)-ind(i-1),1.*mean(signal),'g');
      r = input('Which one to remove (1/2) [0 to add]');
      if r==0
          [x,~] = ginput(1); x = x + ind(i-1); %absolute_index
          pkind(firstcyclepeak+1:end+1) = pkind(firstcyclepeak:end);
          pkamp(firstcyclepeak+1:end+1) = pkamp(firstcyclepeak:end);
          pkind(firstcyclepeak) = round(x);
          pkamp(firstcyclepeak) = signal(pkind(firstcyclepeak));
      end
      close(chek)
      ind(i+r-1) = NaN;
   end
end
ind(isnan(ind)) = [];

disp(strcat('Pks-starts =',num2str(length(pkind)-length(ind))));

%Correction factor: may improve the quality of the cycle start
for i = 1:length(ind)
    ind(i) = ind(i) + find(signal(ind(i):end)>0,1,'first');
end
%remove last one
ind(end) = []; d_pkamp(end) = []; d_pkprom(end) = []; pkind(end) = []; pkamp(end) = [];
end