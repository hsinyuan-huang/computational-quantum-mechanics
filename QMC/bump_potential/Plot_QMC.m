function Plot_QMC()
N = 1;
d = 2;
thermal_steps = 50000;
stable_steps = 20000;

fileID = fopen('QMC_system.txt','r');
qmc_sys = fscanf(fileID, '%f');
fclose(fileID);

fileID = fopen('walkerN.txt','r');
walkerN = fscanf(fileID, '%f');
fclose(fileID);

accu_Nw = 0;
for i = 1:thermal_steps
    Nw = walkerN(i);
    posi = qmc_sys(N * d * accu_Nw + 1 : N * d * (accu_Nw + Nw));
    posi = reshape(posi, d, Nw * N);
    
    accu_Nw = Nw + accu_Nw;
end

frame_cnt = 1;

final_posi = [];
for i = 1:stable_steps
    Nw = walkerN(i + thermal_steps);
    posi = qmc_sys(N * d * accu_Nw + 1 : N * d * (accu_Nw + Nw));
    posi = reshape(posi, d, Nw * N);
    final_posi = [final_posi, posi];
    
    totalNum = size(final_posi(1, :), 2);
    
    if (i < 750) || (mod(i, 50) == 0)
        map = [];
        poss_list = linspace(9999, max(10000-totalNum, 0), 1000);
     
        for x = poss_list
            map = [map; [x / 10000, 1, x / 10000]];
        end
        colormap(map);
        
        xedge = linspace(-7, 7, 50);
        yedge = linspace(-7, 7, 50);
        if (i < 1000)
            h = histogram2(final_posi(1, :), final_posi(2, :), xedge, yedge, 'Normalization', 'count', 'EdgeAlpha', 1, 'FaceAlpha', 0.8);
        else
            h = histogram2(final_posi(1, :), final_posi(2, :), xedge, yedge, 'Normalization', 'count');
        end
        h.FaceColor = 'flat';
        hold on
    
        if (i < 750)
            h2 = scatter3(posi(1, :), posi(2, :), 0 * ones(size(posi(2, :))), 100, 'blue', 'o', 'filled', 'MarkerFaceAlpha', 1);
            uistack(h2, 'top');
            if i < 600
                axis([-7 7 -7 7 0 100]);
            else
                axis([-7 7 -7 7 0 (100 + 700 * (i - 600) / 150)]);
            end
        else
            axis([-7 7 -7 7 0 800]);
        end
    
        set(gcf, 'Position', [100, 100, 1049, 895]);
        
        theta = 60;
        view(-37.5, theta);
        drawnow;
        M(frame_cnt) = getframe(gcf);
        frame_cnt = frame_cnt + 1;
    
        if i ~= stable_steps
            delete(h);
            if i < 750
                delete(h2);
            end
        end
    end
    
    accu_Nw = Nw + accu_Nw;
end

v = VideoWriter('illustration');
open(v);
writeVideo(v, M);
close(v);

figure
xedge = linspace(-7, 7, 50);
yedge = linspace(-7, 7, 50);
h = histogram2(final_posi(1, :), final_posi(2, :), xedge, yedge, 'Normalization', 'count');
h.FaceColor = 'flat';

end