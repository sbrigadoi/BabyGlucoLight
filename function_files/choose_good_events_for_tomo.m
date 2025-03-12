function [event_good] = choose_good_events_for_tomo(eventType,subjectN,time_window);


load("event_stats_"+eventType+"_PD"+subjectN+"_"+time_window+".mat")

%1 = good, 0 = bad
%row 1: check N of goof channels
%row 2: % of T pts, MA train
event_good_bad = zeros(2,size(events_stats.eventN_allE,2));
event_good = zeros(1,size(events_stats.eventN_allE,2));

for i=1:size(events_stats.eventN_allE,2)
    %check N of good channels
    if nnz(events_stats.goodch_idx_allE(:,i))>23
        event_good_bad(1,i) = 1; %good due to N good chs
    end
    %check % of M.A T Pts (out of all time)
    if events_stats.pc_Tpts_MA_train_allE(i,1) < 0.5
        event_good_bad(2,i) = 1; %good due to N good chs
    end
end

for i=1:size(events_stats.eventN_allE,2)
    if sum(event_good_bad(:,i)) == 2
        event_good(1,i) = events_stats.eventN_allE(1,i);
    end
end
event_good = nonzeros(event_good)';



end