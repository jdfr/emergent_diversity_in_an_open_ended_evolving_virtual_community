function [varargout] = cellfunc(varargin)
[varargout{1:nargout}] = cellfun(varargin{:}, 'uniformoutput', false);