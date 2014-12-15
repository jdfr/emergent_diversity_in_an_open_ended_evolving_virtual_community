function l = getLastLine(file)

if not(exist(file, 'file'))
  error('This file does not exist: <%s>', file);
end
f = fopen(file, 'r');
l='';
while not(feof(f))
  l1 = fgetl(f);
  if not(isempty(l1))
    l = l1;
  end
end
fclose(f);