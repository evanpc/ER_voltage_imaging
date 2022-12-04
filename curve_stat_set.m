
%% Inputs %% 

%Trace is a normalized time series with an upward deflection response

%Stim start is the colum number when the singular stimulation occurs

%Peak_range is a reasonable number of columns **after the stiulous onset 
%   (and not within the global number of columns)** that the maximal 
%   response occurs within



%% Outputs %%

%curve.max(1) is the maximal resonse recorded

%curve.max(2) is the index of maximal response

%curve.hwfm


%To note: peak_range is arbitrary in this function, and is set at user's
%   discretion. 
%To note: start and peak_range have to be in # of cells not time units





function curve=curve_stat_set(trace,stim_start,peak_range)


% subtracting off start value from fluorescence
%trace=trace-trace(stim_start); 
trace=trace+1;
% calculating max (in units of trace) and max cell %
[curve.max(1),curve.max(2)]=max(trace(stim_start:(stim_start+peak_range)));
curve.max(2)=curve.max(2)+stim_start-1;%setting to full x not just peak range
curve.max(1)=mean(trace(curve.max(2):(curve.max(2)+3)));


end



%maybe for fwhm the bounds shoud be determined by findin the avarage point
%   of half maxb before and after peak? (instead of first and last).  This
%   may be helpful for particularly noisy or small responses

