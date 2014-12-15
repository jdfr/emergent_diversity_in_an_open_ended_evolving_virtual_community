function entornoRasterLatex(bd, nms, scale, range)
%bd = 'lsys\sc0'; names = dir([bd filesep '*.png']); names = {names.name}'; entornoRasterLatex(bd, names, '0.1');  

%example: entornoRasterLatex('D:/G/pepedavid/geb/thesis/t/lsys/sc/', '0.17', 1:20);

%cf = @(varargin)cellfun(varargin{:}, 'uniformoutput', false);
%af = @(varargin)arrayfun(varargin{:}, 'uniformoutput', false);

%bd = 'lsys\0.75\';adasd010.png
% d = dir([bd '*.png']);
% dn = {d.name}';
% dnn = sort(dn(not([d.isdir])));
% if exist('range', 'var')
%   dnn = dnn(range);
% end
%imfs = cf(@(x)imfinfo([bd x], 'png'), dnn);
%[ws hs] = cellfun(@(x)deal(x.Width, x.Height), imfs);
%nms = cf(@(x)x(end-6:end-4), dnn);

f = fopen([bd filesep 'include.tex'], 'w');

fprintf(f, '\\noindent\n\\begin{center}');
bds = strrep(bd, '\', '/');
indexes = 1:numel(nms);
% idxs = [
% 001
% 002
% 003
% 005
% 004
% 006
% 007
% 008
% 009
% 010
% 011
% 012
% 013
% 015
% 021
% 014
% 017
% 020
% 016
% 018
% 025
% 019
% 022
% 023
% 024
% 026
% 027
% 029
% 030
% 028
% 031
% 032
% 033
% 034
% 035
% 036
% 037
% 038
% 039
% 040
% 041
% 042
% 043
% 044
% 045
% 046
% 047
% 048
% 049
% 055
% 056
% 050
% 052
% 051
% 053
% 054
% 057
% 058
% 059
% 060
% 061
% 062
% 063
% 064
% 065
% 066
% 067
% 068
% 072
% 069
% 070
% 071
% 073
% 074
% 075
% 076
% 077
% 078
% 079
% 080
% 081
% 082
% 083
% 084
% 085
% 086
% 087
% 088
% 089
% 090
% 091
% 092
% 093
% 094
% 095
% 096
% 097
% 098
% 099
% 100
% 101
% 102
% 103
% 104
% 105
% 106
% 107
% 108
% 109
% 110
% 111
% 112
% 113
% 114
% 115
% 116
% 117
% 118
% 119
% 120
% 121
% 122
% 123
% 124
% 125
% 126
% 127
% 128
% ];
% indexes = idxs';
for k=indexes
  p = find(nms{k}=='.');
  if numel(p)>1
    newnm = nms{k};
    newnm(p(1:end-1)) = '_';
    movefile([bd filesep nms{k}], [bd filesep newnm]);
    nms{k} = newnm;
  end
    n = lower(nms{k});
    sub = find(n=='_');
    alpha = n(3:6);
    alpha(alpha=='_') = '.';
    alpha = mat2str(str2num(alpha));
    ng = n(34:end);
    ng = ng(1:find(ng=='_', 1)-1);
    p = strfind(n, '_sc500');
    pn = p+6;
    if isempty(p)
      p = strfind(n, '_sc');
      pn = p+3;
    end
    pixs = n(pn:end);
    pix1 = pixs(1:find(pixs=='x', 1)-1);
    pix2 = pixs(find(pixs=='x', 1)+1:find(pixs=='.')-1);
    ratio = (str2double(pix2)/str2double(pix1));
    fprintf(f, '\\begin{center}\\fbox{\\includegraphics[%s]{%s/%s}} \\nopagebreak \\\\ \\nopagebreak idx %03d, $\\alpha=%s$, last generation: %s, image width: %s, image height: %s \\end{center} \n\n', scale, bds, nms{k}, k, alpha, ng, pix1, pix2);
    %fprintf(f, '\\vbox{\\hbox{\\fbox{\\includegraphics[%s]{%s/%s}}}} \\parbox{\\stretch}{ idx %03d, $\\alpha=%s$, last generation: %s, image width: %s, image height: %s } \n\n', scale, bds, nms{k}, k, alpha, ng, pix1, pix2);


%    %fprintf(f, '\\begin{center}\\includegraphics[%s]{%s/%s}\\end{center} \n\n', scale, bds, nms{k});
%    fprintf(f, '\\begin{center}\\includegraphics[%s,bb={60bp 5bp 1054bp 527bp},clip]{%s/%s}\\end{center} \n\n', scale, bds, nms{k});
end
fprintf(f, '\\end{center}');

fclose(f);
