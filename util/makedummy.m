function newfun = makedummy(fun)
newfun = @(varargin)dummy(fun, varargin);

function d = dummy(fun, args)
d = [];
if iscell(fun)
  for k=1:numel(fun)
    fun{k}(args{:});
  end
else
  fun(args{:});
end
