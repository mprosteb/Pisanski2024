function [DATA,INFO] = LCdataload()
%Loads data into a structure format from a lab chart .mat export 
% Written by Mitchell Prostebby, 2023

disp('Please choose a labchart file')
[file,path] = uigetfile;
load([path file])
if sum(firstsampleoffset,'all','omitnan')>0
    disp('Error: First sample of lab chart file has offset');
    return
end

%Setting up the data structure
%Channels first
INFO.Channels = titles; %#ofchannels
INFO.SampleRates = samplerate;
INFO.Name = convertCharsToStrings(file);
if exist('comtext')~=0
   INFO.COMMENT = comtext;
   INFO.COMMENTTIMES = [];
end
    
nCh = size(samplerate,1);
nBlock = size(dataend,2);
blockdata= [];
blocklen = [0];
commenttimes = [];
for i = 1:nCh
    for j = 1:nBlock
        blocklen(j+1,1) = dataend(i,j)-datastart(i,j)+1;%may cause issues if first sample offset isnt 0 (hench line 9)
        blockdata = vertcat(blockdata,data(datastart(i,j)-firstsampleoffset(i,j):dataend(i,j)-firstsampleoffset(i,j))');
        if exist('com')~=0
            if sum(com(:,2)==j)>0 & i==1
                commenttimes = sum(blocklen(1:j,1))+com([com(:,2)==j],3);
                INFO.COMMENTTIMES = vertcat(INFO.COMMENTTIMES,commenttimes);
                commenttimes = [];
            end
        end
        
    end
    DATA(:,i) = blockdata;
    blockdata = [];
end
end