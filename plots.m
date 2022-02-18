function plots(time_steps, dtau, t_gc, t_ac, t_abc, t_fc, t_ic, t_amb,t_in, t_out, Q_dot,eff, selected_node, n_nodes)
    time_seconds=(1:time_steps)*dtau;
    subplot(2,2,1)
    plot(time_seconds,t_gc, ...
        time_seconds,t_ac, ...
        time_seconds,t_abc, ...
        time_seconds,t_fc, ...
        time_seconds,t_ic, ...
        time_seconds,t_amb(1:time_steps))
    legend({'t\_gc','t\_ac', 't\_abc', 't\_fc', 't\_ic','t\_amb'},'Location','eastoutside')
    ylabel('Temperature (K)')
    xlabel('time (s)')
    subplot(2,2,2)
    plot(time_seconds,t_in(1:time_steps), ...
        time_seconds,t_out, ...
        time_seconds,t_fc)
    legend({'t\_in','t\_out', 't\_fc'},'Location','eastoutside')
    ylabel('Temperature (K)')
    xlabel('time (s)')
    subplot(2,2,3)
    plot(time_seconds,Q_dot)
    ylabel('Q\_dot')
    xlabel('time (s)')
    subplot(2,2,4)
    %plot(T,G_r(1:time_steps))
    %ylabel('G\_r')
    plot(time_seconds,eff)
    ylabel("efficiency")
    xlabel('time (s)')
    sgtitle(['Selected node: ' num2str(selected_node) ' of ' num2str(n_nodes)]) 
end

