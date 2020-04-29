for offset = 1:2
    
    if offset == 1
        F = load('final_test');
        G = load('reweighting/results');
    else
        F = load('final_test_half');
        G = load('reweighting/results_half');
    end
    algos = {'syn. $$\ell_1$$ (\texttt{none})',...
             'ana. $$\ell_1$$ (\texttt{none})',...
             'syn. $$\ell_1$$ (\texttt{abs})',...
             'ana. $$\ell_1$$ (\texttt{energy})',...
             'ana. $$\ell_1$$ (\texttt{energy}) + tdc',...
             'syn. $$\ell_1$$ (\texttt{iterative})',...   
             'ana. $$\ell_1$$ (\texttt{iterative})',...   
             'S-SPAIN H',...
             'A-SPAIN',...
             'Janssen'};

    gaps = 5:5:50;     
    figure;

    %            DR        CP        wDR       wCP       wCPtdc    reDR      reCP      S-SPAIN   A-SPAIN,  Janssen
    colors  = [ .5 .5 .5; .5 .5 .5;  0  0  0;  0  0  0;  0  0  0;  0  1  1;  0  1  1;  0  0  1;  0  0  1;  1  0  0 ];
    styles  = { '--',      '-',      '--',     '-',      '-.',      '--',     '-',      '--',     '-',      '-'     };
    markers = { 'none',    'none',   'none',   'none',   'none',   'none',   'none',   'none',   'none',   'none'  };
    widths  = [ 1,         1,        1,        1,        1,        1,        1,        1,        1,        1       ];

    K = 10;

    %% SNRs
    subplot(2,2,1);
    for algo = 1:10
        if algo < 6
            avg_per_algo = mean(F.SNRs(:,algo,:,:),1,'omitnan');
        else
            if algo < 8
                avg_per_algo = mean(G.SNRs(:,algo-3,:,:),1,'omitnan');
            else
                avg_per_algo = mean(F.SNRs(:,algo-2,:,:),1,'omitnan');
            end            
        end
        avg_per_algo_and_gap = mean(avg_per_algo,4);   
        plot(gaps,avg_per_algo_and_gap(:),...
            'color',colors(algo,:),...
            'linestyle',styles{algo},...
            'marker',markers{algo},...
            'linewidth',widths(algo));
        hold on;
    end
    legend(algos(1:10),'location','southwest','numcolumns',2);
    xlabel('gap length [ms]');
    ylabel('SNR [dB]');
    ylim([0 18]);
    xlim([5 50]);
    grid on;
    ax = gca;
    ax.XTick = 5:5:50;
    ax.YTick = 0:2:18;

    %% ODGs
    subplot(2,2,2);
    for algo = 1:K
        if algo < 6
            avg_per_algo = mean(F.ODGs(:,algo,:,:),1,'omitnan');
        else
            if algo < 8
                avg_per_algo = mean(G.ODGs(:,algo-3,:,:),1,'omitnan');
            else
                avg_per_algo = mean(F.ODGs(:,algo-2,:,:),1,'omitnan');
            end            
        end
        avg_per_algo_and_gap = mean(avg_per_algo,4);   
        plot(gaps,avg_per_algo_and_gap(:),...
            'color',colors(algo,:),...
            'linestyle',styles{algo},...
            'marker',markers{algo},...
            'linewidth',widths(algo));
        hold on;
    end
    legend(algos(1:K),'location','southwest','numcolumns',2);
    xlabel('gap length [ms]');
    ylabel('ODG');
    xlim([5 50]);
    grid on;
    ax = gca;
    ax.XTick = 5:5:50;

    %% PSMs
    subplot(2,2,3);
    for algo = 1:K
        if algo < 6
            avg_per_algo = mean(F.PSMs(:,algo,:,:),1,'omitnan');
        else
            if algo < 8
                avg_per_algo = mean(G.PSMs(:,algo-3,:,:),1,'omitnan');
            else
                avg_per_algo = mean(F.PSMs(:,algo-2,:,:),1,'omitnan');
            end            
        end
        avg_per_algo_and_gap = mean(avg_per_algo,4);   
        plot(gaps,avg_per_algo_and_gap(:),...
            'color',colors(algo,:),...
            'linestyle',styles{algo},...
            'marker',markers{algo},...
            'linewidth',widths(algo));
        hold on;
    end
    legend(algos,'location','southwest');
    xlabel('gap length [ms]');
    ylabel('PSM');
    xlim([5 50]);
    grid on;
    ax = gca;
    ax.XTick = 5:5:50;

    %% PSMts
    subplot(2,2,4);
    for algo = 1:K
        if algo < 6
            avg_per_algo = mean(F.PSMts(:,algo,:,:),1,'omitnan');
        else
            if algo < 8
                avg_per_algo = mean(G.PSMts(:,algo-3,:,:),1,'omitnan');
            else
                avg_per_algo = mean(F.PSMts(:,algo-2,:,:),1,'omitnan');
            end            
        end
        avg_per_algo_and_gap = mean(avg_per_algo,4);   
        plot(gaps,avg_per_algo_and_gap(:),...
            'color',colors(algo,:),...
            'linestyle',styles{algo},...
            'marker',markers{algo},...
            'linewidth',widths(algo));
        hold on;
    end
    legend(algos,'location','southwest');
    xlabel('gap length [ms]');
    ylabel('PSM$$_t$$');
    xlim([5 50]);
    grid on;
    ax = gca;
    ax.XTick = 5:5:50;

    if offset == 1
        sgtitle('Overall comparison of the methods with full offset');
    else
        sgtitle('Overall comparison of the methods with half offset');
    end
    
end