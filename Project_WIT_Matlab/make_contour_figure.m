function make_contour_figure(nodeMatrix, Cu,Cv)
    
    minn = min(Cv);
    maxx = max(Cv);
    
    % Plot the color plot
    figure;
    scatter(nodeMatrix(:,1), -nodeMatrix(:,2), 200, Cv, 'filled');
    xlabel('X Coordinate');
    ylabel('Y Coordinate');
    caxis([minn maxx]);
    colorbar;
    zlabel('Value');
    title('Cv gros')
    axis equal

    minn = min(Cu);
    maxx = max(Cu);
    
    % Plot the color plot
    figure;
    scatter(nodeMatrix(:,1), -nodeMatrix(:,2), 200, Cu, 'filled');
    xlabel('X Coordinate');
    ylabel('Y Coordinate');
    caxis([minn maxx]);
    colorbar;
    zlabel('Value');
    title('Cu')
    axis equal

    writematrix(Cu,'C:\Users\tdieg\team_01\semi_circle\cu.txt')
    writematrix(Cv,'C:\Users\tdieg\team_01\semi_circle\cv.txt')
    writematrix(nodeMatrix,'C:\Users\tdieg\team_01\semi_circle\finalnodes.txt')
end