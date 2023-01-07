function [ThreeDout]=ThreeDkron(x,y)
% this function gets 3-dimensional arrays and outputs the kronecker product page-wise. meaning-- kron(x(:,:,i),y(:,:,i))   
[m,n,r] = size(x); %// Get size
N = m*n; %// number of elements in one 3D slice

%// ------------- PART 1: Get indices for each 3D slice

%// Get the first mxm block of kron-corresponding indices and then add to
%// each such block for the indices corresponding to the kron multiplications
%// of each iteration
a1 = bsxfun(@plus,reshape((0:N-1)*N+1,m,m),permute((0:N-1),[1 3 2]));

%// Now, a1 is a 3D array, we need to make 2D array out of it.
%// So, concatenate along rows to make it a "slimish" 2D array
a2 = reshape(permute(a1,[1 3 2]),size(a1,1)*size(a1,3),[]);

%// Cut after every N rows to make it a square 2D array.
%// These are the indices for each frontal tensor of kron muliplications
slice_idx = reshape(permute(reshape(a2,N,size(a2,1)/N,[]),[1 3 2]),N,N);

%// -------------  PART 2: Get kron equivalent output

%// Perform x:(Nx1) x y:(1xN) multiplications
vals = bsxfun(@times,reshape(x,m*n,1,r),reshape(y,1,m*n,r)); %//multiplications

%// Get indices for all 3D slices and then index into those multiplications
%// with these for the final kron equivalent output
all_idx=bsxfun(@plus,slice_idx,permute((0:r-1)*m*m*n*n,[1 3 2])); %//all indices
ThreeDout = vals(all_idx); %// final output of kron equivalent multiplications
end

