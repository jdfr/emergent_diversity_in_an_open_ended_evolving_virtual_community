function [varargout] = ls2(regla,T,doSanitize,returnSanitized)
% rules are delimited by spaces and sides by '=>', like this 'F=>FF G=>F[+G][-G]F[+G][-G]FG'

% L_SYSTEM_2D simulates the L-system or string-rewriting-system
% in 2 dimensions.
% The recursive principle is the following:
%
%   axiom: initial string
%       |
%       |
%       -----> input string -----> production rules ----
%                   ^                                   |
%                   |                                   |
%                   --------output string <-------------
%
% The construction rules can be specified/varied in the file, as
% well as the initial string and the number of iterations.
%
% The results are plotted with a chosen angle between the axioms
% and a chosen length of each axiom.
% Brackets '[' save the current position and angle, ']' returns
% to the saved position.
%
% The default setting draws a 'branch', a tree-like structure
% plotted with the colours brown and green with 5 iterations.

% This file was generated by students as a partial fulfillment
% for the requirements of the course "Fractals", Winter term
% 2004/2005, Stuttgart University.   
   
% Author : Marko Grob, Wiebke Heidelberg
% Date   : Feb 2005
% Version: 1.0

% 0. PART - INPUT VARIABLES
%============================

%   input rules for the construction of the picture
%   rule(1).pre = 'F';
%   rule(1).pos = 'FF';
%   rule(2).pre = 'G';
%   rule(2).pos = 'F[+G][-G]F[+G][-G]FG';


% Chequear que los corchetes est�n balanceados
if(~corchetesBalanceados(regla))
    error('ls2: regla con corchetes no balanceados: %s', regla);
end


% CONSIDERAR APARTE LOS SIMBOLOS RESERVADOS, NO COMO REGLAS


%ProductionSystem = ['[=>[ ]=>] +=>+ -=>- G=>',regla];
ProductionSystem = ['G=>',regla];

% %BEWARE: THIS CODE OUTPUTS EXACTLY THE SAME RESULTS AS THE OLD CODE, BUT IS
% %NOT PREPARED FOR L-SYSTEMS WITH SEVERAL RULES, SINCE THE REPLACEMENT IS
% %NOT CONCURRENT BUT RECURSIVE!!!!!!!!!!!!!
% % rules are separated by space character, and pre and post sides are
% % separtated by '=>': 'F=>FF G=>F[+G][-G]F[+G][-G]FG'
% rules = regexp(ProductionSystem, '([^ ])=>([^ ]*)', 'tokens'); %1xnrules cell array of 1x2 cell arrays
% n_Rules = numel(rules);
% for i=1:n_Repeats
%   for j=1:n_Rules
%     %for rule j, replace the occurences of the precedent in the string by
%     %the consequent
%     axiom=strrep(axiom, rules{j}{:});
%   end
% end

%BEWARE: THIS IMPLEMENTATION AIMS TO BE EFFICIENT AT THE COST OF REQUIRING
%INPUT STRINGS TO BE UPPER CASE, AND THE PRECEDENTS TO BE LETTERS. IF THESE
%REQUERIMENTS ARE NOT MET, THIS FUNCTION WILL YIELD BADLY WRONG RESULTS

if any(upper(ProductionSystem)~=ProductionSystem)
  error('This code only can handle production systems whose letters are uppercase');
end

% rules are separated by space character, and pre and post sides are
% separatated by '=>': 'F=>FF G=>F[+G][-G]F[+G][-G]FG'
rules = regexp(lower(ProductionSystem), '([^ ])=>([^ ]*)', 'tokens'); %1xnrules cell array of 1x2 cell arrays
n_Rules = numel(rules);

rules = vertcat(rules{:});
pres  = upper(horzcat(rules{:,1})); %precedents  are uppercase
randomExpression = isstruct(T) && (isfield(T, 'randomize')) && logical(T.randomize);
if (nargin<3) || doSanitize
  posts = cellfun(@(x)sanitizeString(x,randomExpression), rules(:,2), 'uniformoutput', false);                 %consequents are lowercase
else
  posts = rules(:,2);                 %consequents are lowercase
end
if ~all(isletter(pres)) || any(pres=='�') || any(pres=='�') % '�' and '�' are deemed letters by matlab, but they have not uppercase counterparts
  error('This code only can handle rules whose precedent are letters!');
end

% starting string (axiom)
axiom = 'G';

% iterations (choose only from 1 to 7, >= 8 critical,
% depends on the string and on the computer !!
n_Repeats = 3;

% 1. PART - CALCULATE THE STRING
%=================================
for i=1:n_Repeats
  for j=1:n_Rules
    %for rule j, replace the occurences of the precedent in the string by
    %the consequent. As consequents are lowercase while axiom and
    %precendents are uppercase, there will be no collisions
    axiom=strrep(axiom, pres(j), posts{j});
  end
  %turn all the string uppercase for the next expansion
  axiom=upper(axiom);
end


varargout = cell(1,nargout);
%[varargout{:}] = tree(axiom, T);
if (nargin>3)&&returnSanitized
  [varargout{1:end-1}] = treeRaster(axiom, T);
  varargout{end} = upper(posts{1});
else
  [varargout{:}] = treeRaster(axiom, T);
end

end

