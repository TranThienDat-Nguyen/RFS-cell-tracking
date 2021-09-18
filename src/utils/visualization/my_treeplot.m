function my_treeplot(p, id)
figure('Name', ['Cell_lineages_', id]);
hold on
LW = 1.5 ; 
x_offset_all = zeros(1,length(p)) ; y_min = 10000 ; 
cnt = 0 ; 
for fidx = 1 : length(p)
    if length(p{fidx})>1
        [x,y,~]=treelayout(p{fidx});
        f = find(p{fidx}~=0);
        pp = p{fidx}(f);
        X{fidx} = [x(f); x(pp); NaN(size(f))];
        Y{fidx} = [y(f); y(pp); NaN(size(f))];
        X{fidx} = X{fidx}(:);
        Y{fidx} = Y{fidx}(:);
        if min(Y{fidx})<y_min
            y_min = min(Y{fidx}) ; 
        end
        if ~isempty(X{fidx})
            x_offset_all(fidx) = max(max(X{fidx})) ; 
        end
        cnt = cnt + 1 ;
    end
end
color_new = distinguishable_colors(cnt) ; 
idx = 1 ; 
for fidx = 1 : length(p)
    if length(p{fidx})>1
        y_offset = min(Y{fidx}) - y_min ; 
        x_offset = sum(x_offset_all(1:fidx))+0.5 ; 
        for n = 1 : 3 : length(X{fidx})
            if ~isnan(X{fidx}(n))
                plot([x_offset+X{fidx}(n), x_offset+X{fidx}(n+1)],  [Y{fidx}(n+1)-y_offset, Y{fidx}(n+1)-y_offset], 'LineWidth', LW, 'Color', color_new(idx,:));
                plot([x_offset+X{fidx}(n), x_offset+X{fidx}(n)],  [Y{fidx}(n)-y_offset, Y{fidx}(n+1)-y_offset], 'LineWidth', LW, 'Color', color_new(idx,:));
            end
        end
        idx = idx + 1 ;
    end
end
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
title('Cell lineages')
