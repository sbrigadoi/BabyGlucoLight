function PD_data = block_avg_data(PD_data,baseline)
    PD_data.dod3 = PD_data.dod2(:,:) - mean(PD_data.dod2(baseline(1):baseline(2),:));
end