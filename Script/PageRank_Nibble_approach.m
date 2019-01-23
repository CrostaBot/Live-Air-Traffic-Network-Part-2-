% close all
clear all
clc

%% Import data

G = importdata('dataset_3.txt', '\t', 4) 

%% pre-processing

G.data = G.data + 1;
N = max(max(G.data));
A = sparse(G.data(:,2),G.data(:,1),ones(size(G.data,1),1),N,N);
Au = 1*(A+A'>0); 
Au = Au - diag(diag(Au)); 

d = full(sum(Au)); 
D = sum(d); 
I = spdiags(ones(N,1),0,N,N); 
Di = spdiags(1./sqrt(d'),0,N,N); 
L = I - Di*Au*Di;
M = Au*Di*Di; 

%% Spectral approach

[V,DD] = eigs(L,2,'SA');
Vv = Di*V;
v1 = Vv(:,2)/norm(Vv(:,2)); 
[v1s,pos] = sort(v1,'descend');
Au1 = Au(pos,pos);
a = sum(triu(Au1));
b = sum(tril(Au1));
assoc = cumsum(a+b);
assoc = min(assoc,D-assoc);
cut = cumsum(b-a);
conduct = cut./assoc;
conduct = conduct(1:end-1);
[~,mpos] = min(conduct);
threshold = mean(v1s(mpos:mpos+1));
disp('spectral approach')
disp(['   Minimum conductance: ' num2str(conduct(mpos))])
disp(['   Cheeger''s upper bound: ' num2str(sqrt(2*DD(2,2)))])
disp(['   # of links: ' num2str(D/2)])
disp(['   Cut value: ' num2str(cut(mpos))])
disp(['   Assoc value: ' num2str(assoc(mpos))])
disp(['   Community size #1: ' num2str(mpos)])
disp(['   Community size #2: ' num2str(N-mpos)])
disp([' '])

%% PageRank-nibble approach

if mpos<N-mpos 
    i = pos(1); 
else
    i = pos(end);
end
q = zeros(N,1);
q(i) = 1;
c = 0.85;
r = (I-c*M)\((1-c)*q);
ep = 1e-3; 

u = zeros(N,1); 
v = q; 
th = full(ep*d/D)'; 
count = 0; 
complexity = 0;
ii = i;
while (count<N)
    if v(ii)>th(ii)
        tmp = v(ii);
        u(ii) = u(ii)+(1-c)*tmp;
        v(ii) = 0;
        v = v + c*M(:,ii)*tmp;    
        complexity = complexity + d(ii); 
        count = 0;
    else 
        count = count + 1; 
        ii = mod(ii,N)+1; 
    end
end

[u1s,pos2] = sort(u,'descend');
Nmax = find(u1s>0,1,'last');
Au1 = Au(pos2,pos2(1:Nmax));
a = sum(triu(Au1));
b = sum(tril(Au1));
assoc = cumsum(a+b);
assoc = min(assoc,D-assoc);
cut = cumsum(b-a);
conduct = cut./assoc;
conduct = conduct(1:Nmax-1);
[~,mpos2] = min(conduct);
threshold2 = mean(u1s(mpos2:mpos2+1));
disp('PageRank-nibble approach')
disp(['   complexity/D: ' num2str((complexity/D))])
disp(['   epsilon: ' num2str(ep)])
disp(['   prec: ' num2str(norm(r-u,1))])
disp(['   Minimum conductance: ' num2str(conduct(mpos2))])
disp(['   # of links: ' num2str(D/2)])
disp(['   Cut value: ' num2str(cut(mpos2))])
disp(['   Assoc value: ' num2str(assoc(mpos2))])
disp(['   Community size #1: ' num2str(mpos2)])
disp(['   Community size #2: ' num2str(N-mpos2)])

figure(2)
plot(conduct)
grid
ylabel('conductance')
title('sweep choice')

figure(1)
plot(u,v1,'k.')
hold on
plot(threshold2*[1,1],ylim,'g-')
plot(xlim,threshold*[1,1],'r-')
hold off
grid
ylabel('Fiedler''s eigenvector value')
xlabel('PageRank value')
title('communities')
