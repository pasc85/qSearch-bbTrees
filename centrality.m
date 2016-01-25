%%% this script computes several centrality measures of the tree

%% parameters to set
n = 12;
beta = 0.5;
% system size
N = 2^n - 1;

%% construct adjacency matrix
A = zeros( N );
for k=2:2:N-1
    A(k/2,k)=1;
    A(k/2,k+1)=1;
    A(k,k/2)=1;
    A(k+1,k/2)=1;
end

%% 0. degree centrality
Df=sum(A)';
Dr = zeros(n,1);
for k = 1:n
    Dr(k,1) = Df(2^(k-1),1);
end
D = zeros(5,1);
D(1,1)=Dr(1,1);
D(2,1)=Dr(n/4,1);
D(3,1)=Dr(n/2,1);
D(4,1)=Dr(3*n/4,1);
D(5,1)=Dr(n,1);

%% 1. communicability
% compute matrix exp and communicability
unit = ones(N,1);
CM = expm(beta.*A);
Cf = CM*unit;
% reduce communicability vector to levels
Cr = zeros(n,1);
for k = 1:n
    Cr(k,1) = Cf(2^(k-1),1);
end
% reduce communicability vector to levels 1,n/4,n/2,3n/4,n
C = zeros(5,1);
C(1,1)=Cr(1,1);
C(2,1)=Cr(n/4,1);
C(3,1)=Cr(n/2,1);
C(4,1)=Cr(3*n/4,1);
C(5,1)=Cr(n,1);

%% 2. subgraph centrality
Sf = diag(CM);
Sr = zeros(n,1);
for k = 1:n
    Sr(k,1) = Sf(2^(k-1),1);
end
S = zeros(5,1);
S(1,1)=Sr(1,1);
S(2,1)=Sr(n/4,1);
S(3,1)=Sr(n/2,1);
S(4,1)=Sr(3*n/4,1);
S(5,1)=Sr(n,1);

%% 3. eigenvector centrality
[Ef,lambda]=eigs(A,1,'la');
Ef=abs(Ef);
Er = zeros(n,1);
for k = 1:n
    Er(k,1) = Ef(2^(k-1),1);
end
E = zeros(5,1);
E(1,1)=Er(1,1);
E(2,1)=Er(n/4,1);
E(3,1)=Er(n/2,1);
E(4,1)=Er(3*n/4,1);
E(5,1)=Er(n,1);

%% 4. Katz centrality
alpha=1/(2*lambda(1,1));
KM=inv(eye(N)-alpha.*A);
Kf=KM*unit;
Kr = zeros(n,1);
for k = 1:n
    Kr(k,1) = Kf(2^(k-1),1);
end
K = zeros(5,1);
K(1,1)=Kr(1,1);
K(2,1)=Kr(n/4,1);
K(3,1)=Kr(n/2,1);
K(4,1)=Kr(3*n/4,1);
K(5,1)=Kr(n,1);


%% output
fprintf('\ndegree centrality, communicability and subgraph centrality (both with beta = %.1f), eigenvector centrality, and Katz centrality (with alpha=1/(2*lambda_1)) at levels 1,n/4,n/2,3n/4,n for n = %i:\n',beta,n);
D
C
S
E
K
