%% compares the method that generates the reduced operators
%% via the full operator and the reduction matrix to the 
%% improved method that defines the operators directly


num_l = 8;	% set to different system sizes
sum_errors = 0;
t1=0;
t2=0;

for l_mv = 1:num_l
	
	% first method
	tic;
	N=2^num_l-1;	
	V=generate_V(num_l,l_mv);
	Vp=pinv(V);
	marked_f=zeros(N,1);
	marked_f(2^(l_mv-1),1)=1;
	initial_f=1/sqrt(N).*ones(N,1);
	marked_r=V*marked_f;
	initial_r=V*initial_f;
	L=generate_L(num_l);
	H=L-diag(marked_f);
	Hb=V*H*Vp;
	t1=t1+toc;
	
	% improved method
	tic;	
	[Lb,imv,w,s]=generate_reduced_L_directly(num_l,l_mv);
    t2=t2+toc;

	% compare
	sum_errors = sum_errors + norm(Hb-Lb+diag(w))+norm(marked_r-w)+norm(initial_r-s);

end

fprintf('sum of all errors: %e\n',sum_errors)
fprintf('elapsed time for first method: %f\n',t1)
fprintf('elapsed time for improved method: %f\n',t2)
fprintf('(the main improvement actually is about memory...)\n')
