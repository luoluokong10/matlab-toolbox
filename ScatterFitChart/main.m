f = uifigure;
g = uigridlayout(f, [1 1]);
x = randn(1000, 1);
y = 2 * x + 1 + randn(size(x));
s = ScatterFitChart('parent',g, 'XData', x,'YData', y);