function movie_out = gen_cell_movie(img_path, img_ext, model, est)
    disp('Generating tracking movie ...')
    % check images info
    img_info = dir(fullfile(img_path, ['*.', img_ext])) ; 
    K = length(img_info) ; 
    labels = est.L ; 
    state = est.X ;

    % plotting properties
    MS = 20 ; % marker size

    labelcount = countlabels (labels) ;
    colorarray = makecolorarray(labelcount);
    for k = 1 : K
      f = figure('visible' , 'off') ;
      img = imread(fullfile(img_path, img_info(k).name)) ; 
      imshow(img) ; 
      hold on ;
      if ~isempty(labels{k}) 
          for t = 1 : length(labels{k})
              Pt = state{k}(model.xy_pos,t) ; 
              color_opt = colorarray.rgb(assigncolor(labels{k}(:,t)),:) ;
              plot( Pt(1),Pt(2),'ks','MarkerSize' , MS ,'Color' , color_opt  ); 
          end
      end
      movie_out(k) = getframe(f) ; 
      close(f) ; 
    end


    function ca= makecolorarray(nlabels)
        lower= 0.1;
        upper= 0.9;
        rrr= rand(1,nlabels)*(upper-lower)+lower;
        ggg= rand(1,nlabels)*(upper-lower)+lower;
        bbb= rand(1,nlabels)*(upper-lower)+lower;
        ca.rgb= [rrr; ggg; bbb]';
        ca.lab= cell(nlabels,1);
        ca.cnt= 0;   
    end

    function idx= assigncolor(label)
        str= sprintf('%.0f-' ,label{1});
        str= str(1:end-1);
        tmp= strcmp(str,colorarray.lab);
        if any(tmp)
            idx= find(tmp);
        else
            colorarray.cnt= colorarray.cnt + 1;
            colorarray.lab{colorarray.cnt}= str;
            idx= colorarray.cnt;
        end
    end

    function count= countlabels(labels)
        strlabel= [];
        for k = 1 : K
            tmp_str_label = cell(1 , length(labels{k}) ) ;
            for lidx = 1 : length(labels{k})
                tmp_str_label{lidx}= sprintf('%.0f,' ,labels{k}{lidx}');
            end
            strlabel = [strlabel tmp_str_label] ;
        end
        [c,~,~]= unique(strlabel ,'rows');
        count=size(c,2);
    end

end

