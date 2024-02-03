function [per,dur]=taskprogress(i,itot)

%% MONITORING TASK PROGRESS

%% INPUT:
%  i = index task progression
%  itot = index of task completion

%% OUTPUT: 
%  per = percentage of already executed task
%  dur = expected time to completion of task

persistent t0  

frac = i/itot;
per  = 100*frac;

t = clock();
if isempty(t0)
    t0 = t;
else
    if (frac < 1)
        for j=1:43
            fprintf('\b')
        end
        fprintf('%05.2f%%',per)  
        dur = etime(clock(),t0) *(itot-i)/(itot*frac);
        fprintf(' | Time to completion: %s',time2str(dur));
        
    else
        for j=1:43
            fprintf('\b')
        end
        fprintf('%06.2f%%',100)  
        fprintf(' | Time to completion: %s\n\n',0);
    end
end

end % End of function taskprogress