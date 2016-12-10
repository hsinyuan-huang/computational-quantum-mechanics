function Plot_QMC()
N = 4;
d = 3;
stable_steps = 400;

fileID = fopen('QMC_system.txt','r');
qmc_sys = fscanf(fileID, '%f');
fclose(fileID);

fileID = fopen('walkerN.txt','r');
walkerN = fscanf(fileID, '%f');
fclose(fileID);

accu_Nw = 0;

final_posi = [];
for i = 1:stable_steps
    Nw = walkerN(i);
    posi = qmc_sys(N * d * accu_Nw + 1 : N * d * (accu_Nw + Nw));
    posi = reshape(posi, d, N * Nw);
    final_posi = [final_posi, posi(1:d, :)];
    
    accu_Nw = Nw + accu_Nw;
end

figure
if d == 1
    histogram(final_posi(1, :), 'FaceColor', 'blue');
    axis([-8 8 0 inf]);
elseif d == 2
    histogram2(final_posi(1, :), final_posi(2, :), 'FaceColor', 'flat');
else
    [threeDH, mid, ~, ~] = histcn(final_posi', 50, 50, 50);

    iso = [60, 300];
    color = {[0.5 0.5 0.5], 'blue'};
    for i = 1:2
        p = patch(isosurface(mid{1}, mid{2}, mid{3}, threeDH, iso(i)));
        isonormals(mid{1}, mid{2}, mid{3}, threeDH, p);
        set(p, 'FaceColor', 'none', 'EdgeColor', color{i}, 'EdgeAlpha', 0.2);
        daspect([1,1,1]);
    end
    box on;
    grid on;
    axis tight;
    %figure
    %scatter3(final_posi(1, :), final_posi(2, :), final_posi(3, :));
end

end