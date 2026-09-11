function plotOptiLandscape(gridNames, gridVecs, gridVar1, gridVar2, metric, h5folder,exportpath,gridResList)
        f= figure;
        t = tiledlayout(1,1);
        f.Position = [100 100 700 500];
        ax1 = axes(t);
        scaling1 = 1;
        scaling2 = 1;
        if gridVar1>=4
            scaling1 = 1/degree;
        end
        if gridVar2>=4
            scaling2 = 1/degree;f
        end  

            imagesc(gridVecs(:,gridVar1)*scaling1,gridVecs(:,gridVar2)*scaling2,metric,'Parent',ax1); 
        if gridVar1>=4
            xlabel(["\Delta" + sprintf(" %s, °", gridNames(gridVar1))]);
        else
            xlabel(["\Delta" + sprintf(" %s", gridNames(gridVar1))]);
        end
        if gridVar2>=4
            ylabel("\Delta" + sprintf(" %s, °", gridNames(gridVar2)));
        else
            ylabel("\Delta" + sprintf(" %s", gridNames(gridVar2)));
        end
        ax1.YAxis.Exponent = 0;
        ax1.XAxis.Exponent = 0;
        % xtickformat('%.3f')
        % ytickformat('%.3f')
        % axis equal xy
        % yyaxis right;
        % xxaxis top;
        ax2 = axes(t);
        % plot(ax2,x2,y2,'-k')
        imagesc(gridVecs(:,gridVar1)/gridResList(gridVar1),gridVecs(:,gridVar2)/gridResList(gridVar1),metric,'Parent',ax2);
        % axis equal xy
        ax2.XAxisLocation = 'top';
        ax2.YAxisLocation = 'right';
        xlabel("\Delta" + sprintf(" %s, px", gridNames(gridVar1)),'Parent',ax2);
        ylabel("\Delta" + sprintf(" %s, px", gridNames(gridVar2)),'Parent',ax2);
       

        % end
        % else
            % xlabel("\Delta" + sprintf(" %s, px", gridNames(gridVar2)),'Parent',ax2);
        % end
        % twiny()
        % xsecondarylabel("\Delta" + sprintf("%s, px", gridNames(gridVar2)))
        colorbar
        exportfolder = fullfile(h5folder,exportpath);
        mkdir(exportfolder);
        exportgraphics(f,fullfile(exportfolder,sprintf("%s_%s.png", gridNames(gridVar1), gridNames(gridVar2))));
        print(f,fullfile(exportfolder,sprintf("%s_%s.svg", gridNames(gridVar1), gridNames(gridVar2))),'-dsvg');
        %%
end