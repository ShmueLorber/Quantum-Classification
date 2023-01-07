function [Hamiltonian] = H_matrix(H,n,k,q) 
%formation of the Hamiltonain with proper constraints:
l=n+1;
H_mat=reshape(H,n,n);
Di=diag(H_mat);     %diagonal of H_mat
H_mat(1:k,1:k)=0;                %input-input interactions
H_mat(1:k,l-k:end)=0;            %input-output without layers
%H_mat(1:k+layer_size,l-k:end)=0;             %input-intermediate-output interactions for layers
%H_mat(1:k,l-k-layer_size:end)=0;             %input-intermediate interactions for layers
H_mat(l-q:end,l-q:end)=0;    %outpu-output interactions

H_mat=H_mat.*~eye(size(H_mat));
H_mat=H_mat+diag(Di);        %comment out for zero diagonal
H_mat=triu(H_mat)+triu(H_mat,1)';
H_mat=[zeros(n,1),H_mat];                    %empty state
H_mat=[zeros(1,l);H_mat];

Hamiltonian = H_mat;
end