% Function used to remove segments of data from a set of timeseries for the
% purposes of artifact removal and data segementation.

% Input up to 5 time-locked time-series, all of which should be the same
% length (and sample rate). Function prompts user for regions of the
% signals to be censored and puts those segments to NaN.

% 'removeopt' (optional) = 'delete', or []
    %'delete' deletes censored regions instead of setting to NaN. Alters
    %the length of inputted signals.

% Written by Mitchell Prostebby, 2023    

function [sigN,sig2N,sig3N,sig4N,sig5N] = datacensor(removeopt,sig,sig2,sig3,sig4,sig5)
sigN = [];
if nargin == 2
    sig2 = []; 
    sig3 = []; 
    sig4 = []; 
    sig5 = []; 
elseif nargin == 3
    sig3 = []; 
    sig4 = []; 
    sig5 = []; sig2N = [];
elseif nargin == 4
    sig4 = []; sig2N = [];
    sig5 = []; sig3N = [];
elseif nargin == 5
    sig5 = []; sig2N = [];
               sig3N = [];
               sig4N = [];
elseif nargin == 6
               sig2N = [];
               sig3N = [];
               sig4N = [];
               sig5N = [];
end

f = figure;
ax1 = subplot(2,1,1); plot(sig);
ax2 = subplot(2,1,2); plot(sig3);
linkaxes([ax1 ax2],'x');

segflag = 'y';
while strcmp(segflag,'y')
    start = input('(Zoom to artifact) Beginning of Artifact (Sample)=');
    % if isempty(start)
    %     disp('Choose Zoom range 1')
    %     [x,~] = ginput(2);
    %     xlim([x(1) x(2)]);
    %     ylim([min(sig(x(1):x(2)))*1.1 max(sig(x(1):x(2)))*1.1]); 
        [start,~] = ginput(1); start = round(start,0);
        [fin,~] = ginput(1); fin = round(fin,0);
        % xlim([1 length(sig)]);
        hold on; scatter(start,mean(sig),'r');
        hold on; scatter(fin,mean(sig),'r');
    % end
    % fin = input('() End of Artifact (Sample) =');
    % if strcmp(fin,'end')
    %     fin = length(sig);
    % elseif isempty(fin)
    %     disp('Choose Zoom range 2')
    %     [x,~] = ginput(2);
    %     xlim([x(1) x(2)]);
    %     ylim([min(sig(x(1):x(2)))*1.1 max(sig(x(1):x(2)))*1.1]); 
    %     [fin,~] = ginput(1); fin = round(fin,0);
    %     xlim([1 length(sig)]);
    % end
    
    sig(start:fin) = NaN; sigN = sig;
    
    if ~isempty(sig2)
        sig2(start:fin) = NaN; sig2N = sig2;
        if ~isempty(sig3)
            sig3(start:fin) = NaN; sig3N = sig3;
            if ~isempty(sig4)
                sig4(start:fin) = NaN; sig4N = sig4;
                if ~isempty(sig5)
                    sig5(start:fin) = NaN; sig5N = sig5;
                end
            end
        end
    end
    segflag = input('Remove more artifacts (y/n)?:','s');
end
if strcmp(removeopt,'delete')
    sigN(isnan(sigN)) = [];
    sig2N(isnan(sig2N)) = [];
    sig3N(isnan(sig3N)) = [];
    sig4N(isnan(sig4N)) = [];
    sig5N(isnan(sig5N)) = [];
end
close(f)
end
