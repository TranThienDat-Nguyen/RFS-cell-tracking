function plot_cardinality(truth, est, identifier)
    %% Plot cardinality
    figure('Name', ['Cardinality_', identifier]); 
    box on; hold on;
    stairs(1:truth.K,truth.N,'k'); 
    plot(1:truth.K,est.N,'k.');

    grid on;
    legend(gca,'True','Estimated');
    set(gca, 'XLim',[1 truth.K]); set(gca, 'YLim',[0 max([max(truth.N), max(est.N)])+1]);
    xlabel('Time'); ylabel('Cardinality');
end