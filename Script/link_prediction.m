close all
clear all
clc

%% Import data

G = importdata('p2p-Gnutella25.txt', '\t', 4);

%% Adjacency matrix

G.data = G.data + 1;
N = max(max(G.data));
A = sparse(G.data(:,2),G.data(:,1),ones(size(G.data,1),1),N,N);
clear G;

%% Setup        

G = digraph(A);
total_links = numedges(G);
train_size = round(total_links*0.3);

for i = 1:train_size
    removal_edge = randi(total_links);
    G = rmedge(G, removal_edge);
    total_links = numedges(G);
end 

A_train = adjacency(G);
A_test = A - A_train;

AUC_vector = [];

%% Neighbour technique

% common neighbours
similarity = A_train*A_train;
AUC_vector = [AUC(A_train, A_test, similarity)];

% Adamic Adar
A_train_1 = A_train ./ repmat(log(sum(A_train,2)),[1,size(A_train,1)]); 
A_train_1(isnan(A_train_1)) = 0; 
A_train_1(isinf(A_train_1)) = 0;  
similarity = A_train * A_train_1;   
clear train1;
AUC_vector = [AUC_vector, AUC(A_train, A_test, similarity)];

% resource allocation
A_train_1 = A_train ./ repmat(sum(A_train,2),[1,size(A_train,1)]); 
A_train_1(isnan(A_train_1)) = 0; 
A_train_1(isinf(A_train_1)) = 0;
similarity = A_train * A_train_1;  
clear train1;
AUC_vector = [AUC_vector, AUC(A_train, A_test, similarity)];

%% Path technique

% local path
beta = 0.0001
similarity = A_train^2 + beta*A_train^3;
AUC_vector = [AUC_vector, AUC(A_train, A_test, similarity)];
  
% Kats 
beta = 0.001
similarity = inv(sparse(eye(size(A_train,1))) - beta*A_train);
similarity = similarity - sparse(eye(size(A_train,1)));
AUC_vector = [AUC_vector, AUC(A_train, A_test, similarity)];

%% Random walk techinique
beta = 0.85
steps = 3

% SWR 
deg = repmat(sum(A_train,2),[1,size(A_train,2)]);
A_train = max(A_train ./ deg,0); 
clear deg;
I = sparse(eye(size(A_train,1)));                                
tempsimilarity = I;                            
stepi = 0; similarity = sparse(size(A_train,1),size(A_train,2));  

while(stepi < steps)
	tempsimilarity = (1-beta)*I + beta * A_train' * tempsimilarity;
	stepi = stepi + 1;
	similarity = similarity + tempsimilarity;
end
    
similarity = similarity+similarity';                        
A_train = spones(A_train); 
AUC_vector = [AUC_vector, AUC(A_train, A_test, similarity)];
    
% LWR  
deg = repmat(sum(A_train,2),[1,size(A_train,2)]);
A_train = max(A_train ./ deg,0); clear deg;                                
I = sparse(eye(size(A_train,1)));                                 
similarity = I;
stepi = 0;

while(stepi < steps)                                     
    similarity = (1-beta)*I + beta * A_train' * similarity;
    stepi = stepi + 1;
end

similarity = similarity+similarity';
A_train = spones(A_train); 
AUC_vector = [AUC_vector, AUC(A_train, A_test, similarity)];

% RWR      
deg = repmat(sum(A_train,2),[1,size(A_train,2)]);
A_train = max(A_train ./ deg,0); 	
clear deg;
I = sparse(eye(size(A_train,1)));                                
similarity = (1 - beta) * inv(I- beta * A_train') * I;
similarity = similarity+similarity';                           
A_train = spones(A_train);  
AUC_vector = [AUC_vector, AUC(A_train, A_test, similarity)];

%% fun CALCAUC

function [auc] = AUC(A_train, A_test, similarity)
    similarity = triu(similarity - similarity.*A_train) - diag(diag(similarity));
    non = 1 - A_train - A_test - eye(max(size(A_train,1),size(A_train,2)));
    A_test = triu(A_test);
    non = triu(non);
    A_test_num = nnz(A_test);
    non_num = nnz(non);
    A_test_pre = similarity .* A_test;
    non_pre = similarity .* non;
    A_test_data =  A_test_pre( A_test ~= 0 )';  
    non_data =  non_pre( non ~= 0 )';   
    labels = [ones(1,size(A_test_data,2)), zeros(1,size(non_data,2))];
    scores = [A_test_data, non_data];
    [X,Y,T,auc] = perfcurve(labels, scores, 1);
end
