function [ParCostFunction]=parallelCF(H,k,n,q,gammain,gammadep,gammaout,psi)
     %This part is intended to create a cost function which finds the most discriminative current, means J_vec=[J_tot;0],[0;J_tot] 
     %notice the actual cost at the end, commonly needs to get modified
     %also notice sometimes modification of the hamiltonian is neccesary
    l=n+1;%total number of sites%
    I=eye(l);
    J_vec=zeros(q,size(psi,1));
    for inx=1:size(psi,1)     %runs on each csi
        Vout=zeros(l,l,q);
        Vin=zeros(l,l,k);%empty array of k lxl for Vin matrices%
        %Vin formation:
        Vhelp=eye(l);
        Vhelp=reshape(Vhelp,[1 l l]);
        Vin(:,1,:)=Vhelp(:,:,2:k+1);
        psihelp=reshape(psi(inx,:),[1 1 k]);
        Vin=Vin.*psihelp;
        Vin=sum(Vin,3);
        Lin=sqrt(gammain).*diss3(Vin,l);
        %Vdep formation:
        if gammadep~=0 %Vdep is pretty big (it's 3D is as big as l-(k+q+2)) so trying to reduce computation time
            Vdep=bsxfun(@times,pagetranspose(Vhelp(:,:,k+2:l-q)),Vhelp(:,:,k+2:l-q));
            Ldep=sqrt(gammadep).*sum(diss3(Vdep,l),3);
        else
            Ldep=zeros(l^2);
        end
        %Vout formation:
        Vhelp=flip(eye(l));
        Vhelp=reshape(Vhelp,[1 l l]);
        Vout(1,:,:)=Vhelp(:,:,1:q);
        gammaouthelp=reshape(gammaout,[1 1 q]);
        Lout=gammaouthelp.*diss3(Vout,l);       
        %final Lindblad equation 
        H_mat=H_matrix(H,n,k,q); %formation of the Hamiltonain with proper constraints:
        VonNeuman=-1i*(kron(I,H_mat)-kron(H_mat.',I)); %commutation of H and rho%
        source=Lin; %dissipator for lindblad in operator%
        sink=sum(Lout,3); %summing of vectorizations for Lout%
        L=VonNeuman+source+Ldep+sink; %full lindbladian%
        %finding rho of the steady state
        rho_vector=null(L);%finding the vector form of rho%
        %assuring trace=1
        diag_index=1:l+1:l^2;
        trace_of_rho=sum(rho_vector(diag_index,:));
        rho_vector_normalized=bsxfun(@rdivide,rho_vector,trace_of_rho);
        rho=reshape(rho_vector_normalized,l,l,size(rho_vector_normalized,2));
        %assuring real solutions%
        helper=rho_vector_normalized(diag_index,:);
        helper=round(helper,9); %sensitivity level for imaginary numbers
        i=helper==real(helper);
        i=all(i);
        rho_physical=rho(:,:,i);
        if size(rho_physical,3)>1
            rho_physical=rho_physical(:,:,1);
        end
        rho_physical=squeeze(rho_physical);
        if isempty(rho_physical)  % in case no real rho was found (numerical issue)
            rho_physical=eye(l);
        end
%         % current calculation: 
%         helper=diag(rho_physical);
%         J_vec(:,inx)=helper(l-q+1:l);
%         J_tot=sum(J_vec(:,inx));
%         J_vec(:,inx)=J_vec(:,inx)./J_tot; %normalization of J_vec
        Vout=flip(Vout,3);
        for indx=1:q              %extraction of the current for every output
            out_lindbladian=diss3(Vout(:,:,indx),l)*reshape(rho_physical,l^2,1);
            out_lindbladian=reshape(out_lindbladian,l,l);
            n_matrix=zeros(l);
            n_matrix(l-q+indx,l-q+indx)=1;
            J=trace(n_matrix*out_lindbladian);
            J_vec(indx,inx)=-J;
        end
    end
    %cost function%
    J_TS=eye(q);    %test currents
    ParCostFunction=sum(vecnorm(J_vec(:,1)-J_TS(:,1),2,1).^2)+...
        sum(vecnorm(J_vec(:,2)-J_TS(:,2),2,1).^2);
    % FIX = sum(vecnorm(J_vec(:,1:23)-J_TS(:,1),2,1).^2)+...
    %     sum(vecnorm(J_vec(:,24:44)-J_TS(:,2),2,1).^2)+...
    %     sum(vecnorm(J_vec(:,45:end)-J_TS(:,3),2,1).^2);
end
