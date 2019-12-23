function [pdf, centers, edges] = computePDF(centers, values, varargin)
% Computes the empirical probability density function underlying the values at
% the points in edges. Edges are the edges of each bin, ie, the first it the
% left side of the first bin and the last the right side of the last bin. The
% vector corresponding to bins centers is returned along the pdf.
%
% INPUTS :
% ------
% edges : bins edges at which the pdf will be computed.
% values : data
% 'method', value : method to estimade pdf : 'hist' (default, simple
% normalized histogram) or 'kde' (kernel smoothing density).
% 'param', cell : additional arguments for histcount or ksdensity
% functions.
%
% OUTPUTS :
% -------
% pdf : empirical probability density function.
% centers : bin centers, pdf should be plotted against this.

% --- Check input
p = inputParser;
p.addRequired('edges', @isnumeric);
p.addRequired('values', @isnumeric);
p.addParameter('method', 'hist', @(x) ischar(x)||isstring(x));
p.addParameter('param', {}, @iscell);
p.parse(centers, values, varargin{:});

values = p.Results.values;
method = p.Results.method;
param = p.Results.param;

if size(centers, 1) ~= 1
    flipcenters = true;
    centers = centers';
else
    flipcenters = false;
end
if size(values, 1) ~= 1
    flipvalues = true;
    values = values';
else
    flipvalues = false;
end

% --- Processing
switch method
    case 'hist'
        d = diff(centers)/2;
        edges = [centers(1)-d(1), centers(1:end-1)+d, centers(end)+d(end)];
        edges(2:end) = edges(2:end)+eps(edges(2:end));
        pdf = histcounts(values, edges, 'Normalization', 'pdf', param{:});
    case 'kde'
        pdf = ksdensity(values, centers, param{:});
end

if flipcenters
    centers = centers';
end
if flipvalues
    pdf = pdf';
end
