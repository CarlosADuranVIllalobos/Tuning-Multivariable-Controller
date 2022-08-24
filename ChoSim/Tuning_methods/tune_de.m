function [M, flag_s, pop_c, pop, trial, new_pop, cost, bests, bests_ind] = tune_de(dy, flag_s, pop_c, popsize, CR, F, pop, trial, new_pop, cost, bests, bests_ind, nd)
%Function to select the constant parameter of the diagonal suppresion matrix
%using diferential evolution
%Copyright Carlos Alberto Duran Villalobos, February 2022-University of Manchester
if (flag_s==1) %first population
    cost(pop_c)=dy;
    if (pop_c < popsize)
        pop_c = pop_c + 1;
        M=pop(pop_c);
    else %max popsize
        [val,~]=max(cost); %min for minimum
        bests(bests_ind)=val;
        bests_ind=bests_ind+1;
        pop_c=1;
        flag_s=2;
    end
end

if (flag_s==2) %trials
    if pop_c==1  %First new mutation
        %%%Mutation
        rands=rand(popsize,3); %Randomly choose 3 individuals for each member of the population
        rnd_in=round(rands*(popsize-1)+1);
        donor=pop(rnd_in(:,1),:)+F*(pop(rnd_in(:,2),:)-pop(rnd_in(:,3),:));%Donor vectors for each individual
        %%%Recombination
        random=rands(:,1);  %create random vector and I_rand:
        Irand=round(1+(nd-1)*rands(:,1:nd));
        
        for ii=1:popsize
            for jj=1:nd
                if random(ii)<=CR || jj==Irand(ii,jj)
                    trial(ii,:)=donor(ii,:);
                elseif random(ii)>CR && jj~=Irand(ii,jj)
                    trial(ii,:)=pop(ii,:);
                end
            end
        end
        M=trial(pop_c);
        pop_c = pop_c + 1;
    else %next new mutations
        new_cost(pop_c-1)=dy;
        if (pop_c-1 < popsize)
            M=trial(pop_c);
            pop_c = pop_c + 1;
        else %max popsize
            %%%Selection  
            for kk=1:popsize
                % For maximization (if you want to minimize, change inequality as "<")
                if cost(kk)>=new_cost(kk)
                    new_pop(kk,:)=pop(kk,:);
                else
                    new_pop(kk,:)=trial(kk,:);
                end
            end
            pop_c=1;
            %flag_s=3;
            flag_s=1;
            M=new_pop(pop_c);
            pop=new_pop; 
        end
    end
end

% if (flag_s==3) %finding best
%     if pop_c==1
%         M=new_pop(pop_c);
%         pop_c = pop_c + 1;
%     else
%         cost(pop_c-1)=dy;
%         if (pop_c < popsize)
%             M=new_pop(pop_c);
%             pop_c = pop_c + 1;
%             
%         else %max popsize
%             [val,~]=max(cost); %min for minimum
%             bests(bests_ind)=val;
%             pop=new_pop;
%             pop_c=1;
%             flag_s=1;
%             M=pop(pop_c);
%             bests_ind=bests_ind+1;
%         end
%         
%         
%     end
%     
%     
% end
% 
% 
% end

