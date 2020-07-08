%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             plotting offset                             %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Date: 08/07/2020
% By Ondrej Mokry
% Brno University of Technology
% Contact: ondrej.mokry@mensa.cz

colors = [  1  0  0;...
            0  0  0;...
            0  1  1;...
            0  0  1;...
            0  0 .5];
        
% alpha = 0.5;
alpha = 1;

weightingtypes = {'\texttt{none}','\texttt{supp}','\texttt{abs}','\texttt{norm}','\texttt{energy}'};
leg = cell(5,1);
figure
h = gobjects(5,1);
for model = 1:2
    XX = zeros(800,5);
    YY = zeros(800,5);
    
    if model == 1
        str = 'analysis';
    else
        str = 'synthesis';
    end
    subplot(1,2,model)
    hold on    
    for weighting = 1:5
        
        % reading the data
        load(['global_test_',str,'_full'])
        x = SNRs.weighted(:,:,:,weighting);
        load(['global_test_',str,'_half'])
        y = SNRs.weighted(:,:,:,weighting);
        
        % scatter
        h(weighting) = scatter(x(:),y(:),[],colors(weighting,:),'MarkerEdgeAlpha',alpha);
        
        % saving the data for the future computation of the percentage
        XX(:,weighting) = x(:);
        YY(:,weighting) = y(:);
        
        % computing and saving the percentage for given weighting
        perc = round(100*sum(y(:)>x(:))/length(x(:)));
        leg{weighting} = [weightingtypes{weighting}, ' (', num2str(perc), ' \%)'];
        
    end
    xlabel('SNR [dB] for full offset')
    ylabel('SNR [dB] for half offset')
    grid on
    box on
    
    % computing the overall percentage
    perc = round(100*sum(YY(:)>XX(:))/length(XX(:)));
    title({[str,' model'],[num2str(perc),' \% above the diagonal line']})
    
    % plotting the diagonal line
    X = get(gca,'XLim');
    Y = get(gca,'YLim');
    maximum = max([X(:); Y(:)]);
    line([-5, maximum+5],[-5, maximum+5],'Color',[.5 .5 .5],'LineWidth',2)
    xlim([0, maximum])
    ylim([0, maximum])
    legend(h,leg,'location','southeast')
end
sgtitle('Comparison of the two offset approaches')