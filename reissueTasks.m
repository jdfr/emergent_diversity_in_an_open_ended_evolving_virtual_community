function job = reissueTasks(jm, path, tasks, independent)

if not(exist('independent', 'var'))
  independent = false;
end

if not(iscell(tasks))
  tasks = array2cell(tasks);
end


if independent
  
  job = cellfunc(@(x)reissueTasks(x, false), tasks);
  job = [job{:}];
  
else
  
  job = createJob(jm, 'Tag', 'reissueTask', 'PathDependencies', genpath(path));
  for k=1:numel(tasks)
    t = tasks{k};
    createTask(job, t.Function, t.NumberOfOutputArguments, t.InputArguments);
  end
  submit(job);
  
end
