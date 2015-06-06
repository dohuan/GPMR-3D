function A = circulant_block_3D(X)
%% X: 3-D base cubic for the block-circulant matrix
% Required: BlockCirculant.m

n1 = size(X,1);
n2 = size(X,2);
n3 = size(X,3);

%A = zeros(n1*n2*n3,n1*n2*n3);
U = [];
for j=1:n3
    M = []; % fundamental block: ABCD
    for i=1:n2
        C_temp = circulant(X(:,i,j)');
        M = [M C_temp];
    end
    % --- Apply BlockCirculant.m function to produce:
    %   ABCD
    %   DABC
    %   CDAB
    %   BCDA
    A_temp_1 = BlockCirculant(M,n2);
    U = [U double(A_temp_1)];
end
% --- treat U as n3 block of circulant matrices
A_temp_2 = BlockCirculant(U,n3);
A = double(A_temp_2);

% % --- treat U as a base matrix (n1xn2)x(n3)
% clear M C_temp
% M = [];
% for i=1:n3
%     C_temp = circulant(U(:,i)');
%     M = [M C_temp];
% end
% A_temp_2 = BlockCirculant(M,n3);
% A = double(A_temp_2);

end

function C = circulant(vec,direction)
% Original file refer to circulant.m
% error checks
if (nargin<1) || (nargin > 2)
    error('circulant takes only one or two input arguments')
end

if (nargin < 2) || isempty(direction)
    direction = 1;
elseif ~ismember(direction,[1,-1])
    error('direction must be either +1 or -1 if it is supplied')
end

if ~isvector(vec)
    error('vec must be a vector')
elseif length(vec) == 1
    C = vec;
    return
end
n = length(vec);
n1 = n-1;

if direction == -1
    C = repmat(0:n1,n,1);
    C = vec(mod(C+C',n)+1);
else
    if size(vec,1) == 1
        rind = 1:n;
        cind =  n + 2 - rind' ;
        cind(cind == (n+1)) = 1;
    else
        cind = (1:n)';
        rind =  n + 2 - cind';
        rind(rind == (n+1)) = 1;
    end
    C = vec(toeplitz(cind,rind));
    
end

end
