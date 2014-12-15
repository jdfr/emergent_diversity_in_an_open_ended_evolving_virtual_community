%RECOGERTRABAJOS 
%   Crea un jobManager, y una lista con todos los trabajos de ese tag.
%   Devuelve el estado y los valores de salida. Imprimirá tantos mensajes
%   de error como tareas hayan finalizado con error. El vector de valores
%   de salida devuelve la cadena 'No terminado' en aquellos trabajos no
%   terminados, y un cell vacío en los terminados con error.
%   Guarda la salida obtenida en un fichero llamado salida_temp.

if ~(exist('tag', 'var'))
    error('Por favor, crea la variable tag con el nombre del Tag');
end
if ~(ischar(tag))
    error('La variable Tag debe ser de tipo cadena de caracteres');
end

jm = findResource('scheduler', 'configuration', 'terclus')

myJobs = findJob(jm, 'Tag', tag)

estado = get(myJobs, 'state')
salida = cell(0);
for p = 1 : length(myJobs)
    if strcmp(estado{p},'finished')
        salida(p) = {getAllOutputArguments(myJobs(p))};
        if isempty(salida{p})
            disp(get(myJobs(p).Tasks, {'ErrorMessage'}))
        end
    else
        salida(p) = {'No terminado'};
    end
end
salida

if length(salida) > 0
    save salida_temp salida
    disp('  Salida guardada');
end