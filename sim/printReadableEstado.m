%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function printReadableEstado(filename, ACIParams, name)
f = fopen(filename, 'w');
printStruct(f, ACIParams, name);
fclose(f);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function printStruct(f, estructura, name)
fields = fieldnames(estructura);
for k=1:numel(fields)
  if isstruct(estructura.(fields{k}))
    printStruct(f, estructura.(fields{k}), [name '.' fields{k}]);
  else
    numChars = 60-numel(name)-numel(fields{k})-1;
    if numChars<0; numChars = 0; end
    fprintf(f, [name '.' fields{k} repmat(' ', 1, numChars) ' = %s;\n'], any2str(estructura.(fields{k})));
  end
end
