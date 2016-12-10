function Plot_eigenstate()
topK = 500; % number of excited states from C++ code
printK = 150; % number of excited states to be drawn
n = 40; % base grid (-n ~ 0 ~ n)
N = 2 * n + 1; % Number of grid points
dx = 0.05;

dont_plot = 0; % draw the excited states or not?

% Plot the potential
x = ((-n:n) * dx);
y = ((-n:n) * dx);
[X, Y] = meshgrid(x, y);
R = (X.*X + Y.*Y).^0.5;

%Z = 0.5 * (X .* X + Y .* Y); %SHO
Z = -600 .* 0.3 ./ (abs(X) + 0.3)  - 600 .* 0.3 ./ (abs(Y) + 0.3) - 600 .* 0.3 ./ (abs(R - 0.8) + 0.3); %Weird Potential

surf(X, Y, Z, 'EdgeColor', 'None', 'facecolor', 'interp');
view(2); axis equal;
axis off;

map = [];
for i = 0:999
    map = [map; [i / 1000, i / 1000, i / 1000]];
end
colormap(map);

saveas(gcf, 'Potential.png')
RemoveWhiteSpace([], 'file', 'Potential.png');

% Plot the eigenfunctions
if dont_plot
    return;
end

fileID = fopen('top_eigenstates.txt','r');
estate = fscanf(fileID, '%f');
fclose(fileID);

estate = reshape(estate, N, N * topK);

for K = 0:(printK-1) 
    [Xq, Yq] = meshgrid((-n:0.2:n) * dx);
    Zq = interp2(X, Y, estate(1:N, N*K+1:N*K+N), Xq, Yq, 'cubic');
    
    surf(Xq, Yq, Zq, 'EdgeColor', 'None', 'facecolor', 'interp');
    shading interp
    
    minPSI = min(min(estate(1:N, N*K+1:N*K+N)));
    maxPSI = max(max(estate(1:N, N*K+1:N*K+N)));
    
    maxABS = max(abs(minPSI), abs(maxPSI)) * 0.5;
    
    map = [];
    for i = 0:999
        curPSI = minPSI + (maxPSI - minPSI) / 1000.0 * i;
        white = min((abs(curPSI) - maxABS) / maxABS, 1);
        
        if curPSI < 0
            if abs(curPSI) > maxABS
                map = [map; [white, 1, white]]; % Green
            else
                map = [map; [0, abs(curPSI) / maxABS, 0]]; % Green
            end
        else
            if abs(curPSI) > maxABS
                map = [map; [1, 1, white]]; % 
            else
                map = [map; [abs(curPSI) / maxABS, abs(curPSI) / maxABS, 0]];
            end
        end
    end
        
    colormap(map);
    view(2);
    axis equal;
    axis off;
    
    saveas(gcf, ['Eigenstate_', num2str(K), '.png']);
    RemoveWhiteSpace([], 'file', ['Eigenstate_', num2str(K), '.png']);
end

end