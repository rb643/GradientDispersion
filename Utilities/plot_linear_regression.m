function plot_linear_regression(x,y)

colour = [0.9412, 0.2314, 0.1255];

[p, S] = polyfit(x, y, 1);
[y_fit,delta] = polyval(p,x,S);
[xsort, idx] = sort(x);
ysort=y_fit(idx);
plot(xsort, ysort+2*delta, '--', 'Color', colour); hold on
plot(xsort, ysort-2*delta, '--', 'Color', colour)
plot(xsort,ysort, 'LineWidth', 2, 'Color', colour);