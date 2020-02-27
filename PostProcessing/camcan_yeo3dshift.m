% main function to quantify gradient shift
function [] = camcan_yeo3dshift(subjectlist, gradientlist, yeomap)
  % function: Short description
  %% camcan_yeo3dshift
  %% 1. Take in first three gradients and compute centre coordinates of gravity
  %% for each yeo network
  %% 2. Generate plot in 3D and triangular 2D and compare controlling

  for s=1:length(subjectlist)
    % grab only first three gradients
    gradient = gradient(:,1:3);
    % get unique yeo points
    for i=1:length(unique(yeomap))
      % grab the points from that yeo network
      yeomask = yeomap == i;
      tempG = gradient(yeomask,:);
      % get the median of each gradient to obtain centroid coordinates
      G3D(;,) = median(tempG);
    end

  end
end

%% maybe a plotting function that animates the  gradientshift by a specific ordering variable
function [ out ] = camcan_yeo3dshift_plot(ordervariable,subjectlist,gradientlist,yeomap)
  % function: Short description
  % IN: an ordering variable (e.g. age)
  %     a bin size for the animation (age-bins in this case)
  %     a list of gradients
  %     a list of yeo classes of the same dimension as the gradient rows


end  % function
