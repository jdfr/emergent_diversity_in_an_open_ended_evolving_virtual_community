function [varargout] = drawtreeP(P, idx, varargin)

[varargout{1:nargout}]= drawtree(P.raster{idx},P.dimensions(idx,:),P.offset(idx,:),varargin{:});
