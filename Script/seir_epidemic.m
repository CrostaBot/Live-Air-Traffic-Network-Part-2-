close all
clear all
clc

%% Import data

G = importdata('dataset_3.txt', '\t', 4) 

%% Adjacency matrix

G.data = G.data + 1 
N = max(max(G.data)) 
A = sparse(G.data(:,2),G.data(:,1),ones(size(G.data,1),1),N,N) 
clear G 

%% Setup

disease = [0,0,0.025,0.075,0.175,0.225,0.250,0.250]; 

initial_states = zeros(1,N);
for i = 1:15
    initial_states(randi(N)) = 1;
end

results = epidemic(initial_states,A,disease);
figure(1)
plot(results);
title('Epidemic Simulation')
legend('Subsceptible','Infected','Recovered');
axis([0 max(size(results)) 0 N]);

%% Epidemic step        

function new_states = epidemic_step(old_states, graph, disease);

infection = zeros(length(old_states),1);     
for individual = 1:length(old_states)
    if (old_states(individual) > 0)
        infection(individual) = disease(old_states(individual));
    end;
end;

probability = (graph*infection);

for individual = 1:length(old_states)
  if (old_states(individual) > 0)
      if (old_states(individual) == length(disease))
          new_states(individual) = -1;
      else
          new_states(individual) = old_states(individual) + 1;
          end;
  else
      if (old_states(individual) == 0)
          if (rand<probability(individual))
              new_states(individual) = 1;
          else
              new_states(individual) = 0;
              end;
      else
          new_states(individual) = -1;
          end;
      end;
end;
end

%% History

function history = epidemic(initial_states,graph,disease);

states = initial_states;
count = sir(states);
history = count;

while (count(2)>0)
    states = epidemic_step(states,graph,disease);
    count = sir(states);
    history(length(history(:,1))+1,:) = count;
end;
end

%% SEIR

function seir = sir(states);

seir = zeros(1,3);

for individual = 1:length(states)
    switch states(individual)
        case -1
            seir(3) = seir(3)+1;
        case 0
            seir(1) = seir(1)+1;
        otherwise
            seir(2) = seir(2)+1;
    end;
end;
end
