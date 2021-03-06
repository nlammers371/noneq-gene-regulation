function find_duplicate(fname)

P=path;
P=strsplit(P, pathsep());
% mydir='/home/myusername/matlabdir';
% P=P(strncmpi(mydir,P,length(mydir)));
P=cellfun(@(x) what(x),P,'UniformOutput',false);
P=vertcat(P{:});
Q=arrayfun(@(x) x.mat,P,'UniformOutput',false); % Q is a cell of cells of strings
Q=vertcat(Q{:});
R=arrayfun(@(x) repmat({x.path},size(x.m)),P,'UniformOutput',false); % R is a cell of cell of strings
R=vertcat(R{:});
[C,ia,ic]=unique(Q);
f_index = find(ismember(C,{fname}));
for c = f_index
    ind=strcmpi(C{c},Q);
%    if sum(ind)>1
       fprintf('duplicate %s at paths\n\t',C{c});
       fprintf('%s\n\t',R{ind});
       fprintf('\n');
%    end
end

end