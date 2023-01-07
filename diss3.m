function[u3]=diss3(A,l)
% 3D version of the lindblad superoperator
I=repmat(eye(l),[1 1 size(A,3)]);
u3=pagemtimes(ThreeDkron(A,I),ThreeDkron(I,conj(A)))-0.5.*(pagemtimes(ThreeDkron(I,pagetranspose(A)),ThreeDkron(I,conj(A)))+pagemtimes(ThreeDkron(conj(pagetranspose(A)),I),ThreeDkron(A,I)));
% u3=ThreeDkron((A').',A)-0.5.*(pagemtimes(ThreeDkron(I,pagetranspose(A)),ThreeDkron(I,conj(A)))+pagemtimes(ThreeDkron(conj(pagetranspose(A)),I),ThreeDkron(A,I)));

end