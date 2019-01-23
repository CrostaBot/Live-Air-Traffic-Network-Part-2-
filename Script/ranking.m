close all
clear all
clc

%% Import data

G = importdata('dataset_3.txt', '\t', 4);

%% Adjacency matrix

G.data = G.data + 1;
N = max(max(G.data));
A = sparse(G.data(:,2),G.data(:,1),ones(size(G.data,1),1),N,N);
Au = 1*(A+A'>0);
clear G;

%% PageRank

M = A*sparse(diag(1./sum(A)));
c = 0.85;
q = ones(N,1)/N;
disp('Computing PageRank - linear system solution')
tic
r = sparse((eye(N)-c*M)/(1-c))\q;
r = r/sum(r);
toc

disp('Computing PageRank - power iteration')
tic
p0 = ones(N,1)/N;
s = [];
for k = 1:35
    p00 = p0;
    p0 = c*M*p0+(1-c)*q;
    p0 = p0/sum(p0);
    s(k) = norm(p0-r)/sqrt(N);
end
toc 
    
disp('Extracting PageRank eigenvalues')
tic
eig = eigs(M,size(M,1));
toc

figure(1)
set(0,'defaultTextInterpreter','latex')
ref = (c*abs(eig(end-1))).^(1:35);
semilogy([s;ref/ref(end)*s(end)]')
grid
legend('power iteration','second eigenvalue')
xlabel('k [iteration \#]')
ylabel('$\|r_k - r_\infty\|$')
title('PageRank convergence')

figure(2)
plot(eig,'x')
hold on
plot(exp(2i*pi*(0:0.001:1)))
hold off
grid
title('PageRank eigenvalues')

%% HITS - authorities

M = A*A';
disp('Computing HITS - eigenvalue extraction')
tic
[pp,ee] = eigs(M,2);
toc
p = -pp(:,1)/norm(pp(:,1));

disp('Computing HITS - power iteration')
tic
N = size(M,1);
p0 = ones(N,1)/sqrt(N);
s = [];
for k = 1:35
    p00 = p0;
    p0 = A*(A'*p0);
    p0 = p0/norm(p0);
    s(k) = norm(p0-p00)/sqrt(N);
end
toc 

figure(3)
ref = (ee(2,2)/ee(1,1)).^(1:35);
semilogy([s;ref*s(end)/ref(end)]')
grid
legend('power iteration','second eigenvalue')
xlabel('k [iteration \#]')
ylabel('$\|r_k - r_\infty\|$')
title('HITS convergence - authorities')

%% HITS - hubs

A = A';
M = A*A';
disp('Computing HITS - eigenvalue extraction')
tic
[pp,ee] = eigs(M,2);
toc
h = -pp(:,1)/norm(pp(:,1));

disp('Computing HITS - power iteration')
tic
N = size(M,1);
p0 = ones(N,1)/sqrt(N);
s = [];
for k = 1:35
    p00 = p0;
    p0 = A*(A'*p0);
    p0 = p0/norm(p0);
    s(k) = norm(p0-p00)/sqrt(N);
end
toc 

figure(4)
ref = (ee(2,2)/ee(1,1)).^(1:35);
semilogy([s;ref*s(end)/ref(end)]')
grid
legend('power iteration','second eigenvalue')
xlabel('k [iteration \#]')
ylabel('$\|r_k - r_\infty\|$')
title('HITS convergence - hubs')


%% Comparison 

figure(5)
plot([h/sum(h),-r/sum(r)])
grid
legend('HITS','PageRank')
title('PageRank vs HITS ')

figure(6)
plot(h/sum(h),r/sum(r),'x')
grid
xlabel('HITS score')
ylabel('Pagerank score')
title('PageRank vs HITS ')

%% Node score

disp(' ')
disp('PageRank')
for i= 1:10
    [~,node] = max(r);
    disp([num2str(node) ' & ' num2str(max(r)) ' \\']);
    r(node)= [];
end

disp(' ')
disp('HITS - authorities')
for i= 1:10
    [~,node] = max(p);
    disp([num2str(node) ' & ' num2str(max(p)) ' \\']);
    p(node)= [];
end

disp(' ')
disp('HITS - hubs')
for i= 1:10
    [~,node] = max(h);
    disp([num2str(node) ' & ' num2str(max(h)) ' & ']);
    h(node)= [];
end
