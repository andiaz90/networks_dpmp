function [imax] = imax(MatVar)
% Finds maximum element of a vector or an array
% If row or column vector, returns a scalar

index = size(MatVar);
Index = index*0;
M = max(MatVar(:));
A = find(MatVar==max(MatVar(:)),1);
for i = 1:length(index)
   Index(i) = mod(ceil(A),index(i));
   A = A/index(i);
end
Index(Index==0)=index(Index==0);

imax=Index;


% If row or column vector, returns a scalar
if size(MatVar,1)*size(MatVar,2)==numel(MatVar)
    if size(MatVar,1)==1 && size(MatVar,2)>1
        imax=imax(2);
    end
    if size(MatVar,1)>1 && size(MatVar,2)==1
        imax=imax(1);
    end
end
        

end

