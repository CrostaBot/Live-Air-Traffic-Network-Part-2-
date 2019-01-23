close all
clear all
clc
tic 
%% Import data

G = importdata('dataset_3.txt', '\t', 4) 

%% Adjacency matrix

G.data = G.data + 1;
N = max(max(G.data));
A = sparse(G.data(:,2),G.data(:,1),ones(size(G.data,1),1),N,N);
Au = 1*(A+A'>0); 
Au = Au - diag(diag(Au)); 

%% Extract eigensystem info

d = full(sum(Au));
Di = spdiags(1./sqrt(d'),0,N,N);
L = spdiags(ones(N,1),0,N,N) - Di*Au*Di;

if N<2e3 
    la = eigs(L,N);  
    figure(1)
    plot(la,'x')
    grid
    title('eigenvalues (of the normalized Laplacian)')
end

[V,DD] = eigs(L,3,'SA');
Vv = Di*V;
v1 = Vv(:,2)/norm(Vv(:,2));
v2 = Vv(:,3)/norm(Vv(:,3)); 

%% Sweep wrt the ordering identified by v1

[v1s,pos] = sort(v1);
Au1 = Au(pos,pos);
a = sum(triu(Au1));
b = sum(tril(Au1));
d = a+b;
D = sum(d);
assoc = cumsum(d);
assoc = min(assoc,D-assoc);
cut = cumsum(b-a);
conduct = cut./assoc;
conduct = conduct(1:end-1);

figure(2)
plot(conduct,'x-')
grid
title('conductance')

[~,mpos] = min(conduct);
threshold = mean(v1s(mpos:mpos+1));
disp(['Minimum conductance: ' num2str(conduct(mpos))])
disp(['   Cheeger''s upper bound: ' num2str(sqrt(2*DD(2,2)))])
disp(['   # of links: ' num2str(D/2)])
disp(['   Cut value: ' num2str(cut(mpos))])
disp(['   Assoc value: ' num2str(assoc(mpos))])
disp(['   Community size #1: ' num2str(mpos)])
disp(['   Community size #2: ' num2str(N-mpos)])

figure(3)
plot(v1,v2,'.')
grid
hold on
plot(threshold*[1,1],ylim,'r-')
hold off
title('communities')

if D<1e4 
    figure(4)
    plot(v1,v2,'.')
    grid
    [I,J,~] = find(Au);
    hold on
    plot([v1(I),v1(J)]',[v2(I),v2(J)]')
    plot(threshold*[1,1],ylim,'k-')
    hold off
    title('communities (with links)')
end
 
if (mpos>=N-mpos)
    A = A(pos(1:mpos),pos(1:mpos));
else
    A = A(pos(mpos+1:end),pos(mpos+1:end));
end
save('previous_community','A')

toc
