function [varargout] = arrayfunc(varargin)
[varargout{1:nargout}] = arrayfun(varargin{:}, 'uniformoutput', false);