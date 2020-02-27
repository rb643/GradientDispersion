%% TOPOLOGY
matrix = thresholdMatrix;
weighted = 1;

for topology_figures = 1
       
    % degree and strength distribution
    [x1, y1] = ksdensity(sum(matrix>0));
    [h, pValue, KSstatistic, criticalValue] = kstest(sum(matrix>0))

    [x2, y2] = ksdensity(mean(matrix));
    [h, pValue, KSstatistic, criticalValue] = kstest(mean(matrix))

    hubness = mean(mean(matrix)) + (2*(std(mean(matrix))));
    hub = mean(matrix);
    hub(hub<hubness)=0;

    % create random matrices
    rand_matrix = randomize_matrix(matrix);

    
    % basic topology
    if weighted == 1
            Cp = clustering_coef_matrix(matrix, 'O');
            rand_Cp = clustering_coef_matrix(rand_matrix, 'O');
            
            Gmatrix = sparse(matrix);
            rand_Gmatrix = sparse(rand_matrix);

    else
            Cp = clustering_coef_bu(matrix);
            rand_Cp = clustering_coef_bu(rand_matrix);
            
            Gmatrix = sparse(matrix>0);
            rand_Gmatrix = sparse(rand_matrix>0);

    end
    
    dist_matrix = graphallshortestpaths(Gmatrix);
    Lp = mean(dist_matrix);  
    rand_dist_matrix = graphallshortestpaths(rand_Gmatrix);
    rand_Lp = mean(rand_dist_matrix);
    
    % normalised topology
    Cnorm = Cp ./ rand_Cp;
    Lnorm = Lp ./ rand_Lp;

    f=figure;

    % degree distribution
    a(1) = axes('position',[0.125 .75 0.15 .2]); 
    plot(y1,x1,'Color', 'r', 'LineWidth', 3);

    % strength distribution
    a(2) = axes('position',[0.325 .75 0.15 .2]); 
    plot(y2,x2,'Color', 'r', 'LineWidth', 3);

    % nodal strength on surface  
    uparcel = unique(nodes);
    for this_parcel = 1:length(uparcel)
       strength_whole(nodes==uparcel(this_parcel)) = mean(matrix(this_parcel,:),2);
    end

    positions = [   0.7 0.8 0.2 0.2; 0.5 0.8 0.2 0.2; ...
                    0.7 0.65 0.2 0.2; 0.5 0.65 0.2 0.2]; 
    extent=[min(strength_whole) max(strength_whole)];
    ax = 1:4;            
    BoSurfStat_calibrate4Views(strength_whole, Slo, positions, ax, extent, (colorbrewer.seq.YlOrRd{1,9}./256))

    a(3) = axes('position',[0.9 .775 0.03 .1]); 
    imagesc(linspace(min(strength_whole), max(strength_whole), 100)',[min(strength_whole) max(strength_whole)]), axis off
    colormap(a(3),flipud((colorbrewer.seq.YlOrRd{1,9}./256)))

    uparcel = unique(nodes);
    for this_parcel = 1:length(uparcel)
       Cp_whole(nodes==uparcel(this_parcel)) = Cnorm(this_parcel);
       Lp_whole(nodes==uparcel(this_parcel)) = Lnorm(this_parcel);
    end

    % nodal clustering on surface     
    positions = [   0.3 0.5 0.2 0.2; 0.1 0.5 0.2 0.2; ...
                    0.3 0.35 0.2 0.2; 0.1 0.35 0.2 0.2]; 
    extent=[min(Cp_whole) max(Cp_whole)];
    ax = 1:4;            
    BoSurfStat_calibrate4Views(Cp_whole, Slo, positions, ax, extent, (colorbrewer.seq.Purples{1,9}./256))

    a(3) = axes('position',[0.5 .45 0.03 .1]); 
    imagesc(linspace(min(strength_whole), max(strength_whole), 100)',[min(strength_whole) max(strength_whole)]), axis off
    colormap(a(3),flipud((colorbrewer.seq.Purples{1,9}./256)))
    
    % nodal path length on surface     
    positions = [   0.7 0.5 0.2 0.2; 0.5 0.5 0.2 0.2; ...
                    0.7 0.35 0.2 0.2; 0.5 0.35 0.2 0.2]; 
    extent=[min(Lp_whole) max(Lp_whole)];
    ax = 1:4;            
    BoSurfStat_calibrate4Views(Lp_whole, Slo, positions, ax, extent, (colorbrewer.seq.Purples{1,9}./256))

    %a(4) = axes('position',[0.95 .425 0.03 .1]); 
    %imagesc(linspace(min(Cp_whole), max(Cp_whole), 100)',[min(Cp_whole) max(Cp_whole)]), axis off
    %colormap(a(4),flipud(inferno))

    % distance
    distance_matrix = matrix .* double(geoD);
    distance_vector = mean(distance_matrix);

    %distance_matrix_bin = (Z>0) .* double(geoD);
    %distance_vector_bin = mean(distance_matrix_bin);

    [x3, y3] = ksdensity(distance_vector);
    a(4) = axes('position',[0.575 .1 0.15 .2]); 
    plot(y3,x3,'Color', 'b', 'LineWidth', 3);


    uparcel = unique(nodes);
    for this_parcel = 1:length(uparcel)
       distance_whole(nodes==uparcel(this_parcel)) = distance_vector(this_parcel);
    end


    positions = [   0.3 0.2 0.2 0.2; 0.1 0.2 0.2 0.2; ...
                    0.3 0.05 0.2 0.2; 0.1 0.05 0.2 0.2]; 
    extent=[min(distance_whole) max(distance_whole)];
    ax = 1:4;            
    BoSurfStat_calibrate4Views(distance_whole, Slo, positions, ax, extent, (colorbrewer.seq.Blues{1,9}./256))

    a(5) = axes('position',[0.5 .1 0.03 .1]); 
    imagesc(linspace(min(distance_whole), max(distance_whole), 100)',[min(distance_whole) max(distance_whole)]), axis off
    colormap(a(5),flipud(colorbrewer.seq.Blues{1,9}./256))

    exportfigbo(f,[UPATH 'Figure3_topology.png'],'png',10);

end 