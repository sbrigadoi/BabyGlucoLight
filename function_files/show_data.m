function show_data(PD_data, dodConv,ch_idx,data_type)


if data_type == 1;
    figure()
    plot(PD_data.t,    dodConv(:,ch_idx));
    xlabel('t / s');
    ylabel('DOD / A.U');
    title('DODConv');
    ylim([-2 2])
end

if data_type == 2;
    figure()
    plot(PD_data.t_ds, PD_data.dod(:,ch_idx));
    xlabel('t / min');
    ylabel('DOD / A.U');
    title('DODConv');
    ylim([-2 2])
end


if data_type == 3;
    figure()
    plot(PD_data.t,    dodConv(:,ch_idx));
    xlabel('t / s');
    ylabel('d HbO /M');
    title('d HbO');
    %ylim([-2 2])
end


if data_type == 4;
    figure()
    plot(PD_data.t,    dodConv(:,ch_idx));
    xlabel('t / s');
    ylabel('d Hb /M');
    title('d Hb');
    %ylim([-2 2])
end

if data_type == 5;
    figure()
    plot(PD_data.t,    dodConv(:,ch_idx));
    xlabel('t / s');
    ylabel('d HbT /M');
    title('d HbT');
    %ylim([-2 2])
end



end