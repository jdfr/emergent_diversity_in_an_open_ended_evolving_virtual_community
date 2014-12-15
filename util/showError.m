%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str = showError(error)
if isempty(error)
    str='<EMPTY ERROR>';
elseif (~isstruct(error)) && (~isa(error, 'MException')) %&& (~all(ismember({'message', 'identifier', 'stack'}, lower(fieldnames(error)))))
    str=['<NO STRUCT ERROR: ' any2str(error) '>'];
else
    lines = cell(size(error.stack));
    for m=1:numel(error.stack)
        lines{m} = sprintf('    File: %s\n  Name: %s\n  Line: %g\n', error.stack(m).file, error.stack(m).name, error.stack(m).line);
    end
    str = horzcat(...
        sprintf('ERROR:\n  Message:    <<%s>>\n  Identifier: <<%s>>\n  Stack:\n', error.message, error.identifier), ...
        lines{:});
end
str = strrep(str, '\', '\\');
