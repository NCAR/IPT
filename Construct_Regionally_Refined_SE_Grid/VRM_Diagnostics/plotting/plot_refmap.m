function plot_refmap(file)

% Load the whole refmap
X = load(file);

contourf(X);