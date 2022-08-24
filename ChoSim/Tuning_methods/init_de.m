function [M,GC,CR,F,nd,popsize,pop,bests,bests_ind,count,pop_c,flag_s,cost,new_cost,trial,new_pop] = init_de(mina, maxa, u_length, num_u, batch_opt)
%Function to select the constant parameter of the diagonal suppresion matrix
%using diferential evolution
%Copyright Carlos Alberto Duran Villalobos, February 2022-University of Manchester


popsize = 4;         %size of the population (NP)
GC = [maxa;mina];          %Gene constraints
CR = 0.8;           %Cross-over Rate
F = 0.5;             %Mutation factor
nd = size(GC,2);      %Dimensionality of the problem 1 for maxa could be 2 for maxa+mina
%DE Initialization

%DE Initialization
pop=ones(popsize,nd);
for i=1:nd
    pop(:,i)=unifrnd(GC(2*i),GC(2*i-1),popsize,1); %create population
end
trial=zeros(popsize,nd);            % Define trial vector
new_pop=ones(popsize,nd);           %Define new selected population
maxIter = round(batch_opt/popsize); %Maximum number of iterations could be less for more dimensionality /popsize*n
bests=zeros(maxIter,1);          % Create the output matrix "bests"
bests_ind=1;                       % index of bests
count=maxIter;
pop_c=1;                         %Set the counter for the population element evaluated in the simulation
flag_s=1;                        %Set a flag for the DE iteration: 1=Initial population 2=trial population 3=best for each iteration
cost=ones(popsize,nd);           %Create fitness values matrix
new_cost=ones(popsize,nd);       %Create new fitness values matrix

M=pop(pop_c)*ones(u_length*num_u,1);
end

