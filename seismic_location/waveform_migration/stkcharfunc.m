function stk=stkcharfunc(data)
% This function is used to stack the characteristic function to image the
% earthquake.
%
% INPUT--------------------------------------------------------------------
% data: the input charactersitic function, 2D array, shape: nt*nre;
%
% OUTPUT:------------------------------------------------------------------
% stk: scalar, the final stacked value for the input characteristic function;


stk=sum(data,'all');



end