%% Details


%This function compiles multiple replicate stimulations in one recording to
%generate the average, stdev, and stderr of the responses in one trace



%To note: stim times need to be listed in a 1xn row

%To note: stim times, pre, and post, need to be numbers of cells

%To note: THIS FUNCTION NORMALIZES THE COMPILED TRACES TO SET PERIOD OF






function compiled=bksub_cmple(trace,stim_times,pre,post)

for i=1:size(stim_times,2)
    compiled.traces(i,:)=(trace(stim_times(i)-pre:stim_times(i)+post))-...
        (mean(trace((stim_times(i)-15):(stim_times(i)-5))));
end

    compiled.stats.ave=mean(compiled.traces,1);
    compiled.stats.stdev=std(compiled.traces,1);
    compiled.stats.stderr=compiled.stats.stdev/(sqrt(3));
    compiled.stats.t_stderr=compiled.stats.ave+compiled.stats.stderr;
    compiled.stats.b_stderr=compiled.stats.ave-compiled.stats.stderr;
end