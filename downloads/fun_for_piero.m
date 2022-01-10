function cl = fun_for_piero(ID)
[a,~,c] = unique(ID);
ii = reshape(c,size(ID));
k = 1;
i1 = 1;
cl{k} = [];
c = 1:numel(a);
while ~isempty(ii)
  t = any(ismember(ii,i1),2);
  if ~any(t)
    k = k + 1;
    cl{k} = [];
    i1 = c(1);
  else
    p = unique(ii(t,:));
    cl{k} = [cl{k};sort(a(p(:)))];
    ii = ii(~t,:);
    c = setdiff(c,p);
    i1 = p(:);
  end
end