% This is a script for the Gillespie Model  %
% that tries to understand relation between %
% robustness and evolvability               %
clc;
clear all;

w = 20;
% Parameters for the model
mu = 0.9;  % Probability of mutation
q  = 0.9;  % Probability that the mutation is neutral
K  = 70; %[5,15,40,70,100];      % 30,100 % Number of accessible phenotype(phenotypic neighborhood) for a genotype
P  = 100;    % Number of phenotypes in the fitness landscape
N  = 5000;  % Number of individuals

qvec = linspace(0.01,0.99,20); %Vector of robustness values

itime  = 10000; % Iteration time for data collection
nrep = 100; % Number of replications
timelimit = 40000; % length of running time

tim = zeros(itime + 1, nrep); %Time vector
ss = 100;
%Vector of species number over time and replications
NA = zeros(itime + 1, nrep);
NB = zeros(itime + 1, nrep);
NC = zeros(itime + 1, nrep);

NAA  = zeros(ss+1, nrep);
NBB  = zeros(ss+1, nrep); 

%Mean vectors over replications
MeanA = zeros(itime + 1, length(qvec));
MeanB = zeros(itime + 1, length(qvec));
MeanC = zeros(itime + 1, length(qvec));

%Vector for getting adaptation time over each replication
%and over all replications
adtime = zeros(length(qvec), nrep);
ngeneration = zeros(length(qvec), nrep);
cumadtime = zeros(length(qvec), 1);
firstgeneration = zeros(length(qvec), nrep);
mean_firstgeneration = zeros(length(qvec), 1);


% Loop over all robustness values
for qi = 1:length(qvec)
    % Rate constants for the individual reactions
    r0 = 0;
    r1 = mu*(1-(K/P))*qvec(qi);
    r2 = mu*(K/P)*qvec(qi);
    r3 = mu*(K/P)*qvec(qi);
    r4 = mu*(1-(K/P))*qvec(qi);
    r5 = mu*(1-qvec(qi))*(1/K);
    sumadtime = 0; % Variable to hold the sum of adaptation times
    sumgeneration = 0;
    sumfirstgeneration = 0;
    %#############################################


    %Loop over ensembles
    for rep = 1:nrep
        NAA(1, rep) = 0.5*N; 
        NBB(1, rep) = 0.5*N; 
        
        % taking the A and B population to a steady state
        
        for u = 1:ss
        
           sum =  r0*(NAA(u, rep) + NBB(u, rep)) + ...
                    r1*NAA(u,rep) + r2*NAA(u,rep) + ...
                    r3*NBB(u,rep) + r4*NBB(u,rep);
            
            % Reaction propensities
            s1 = r0*NAA(u, rep)/sum;
            s2 = s1 + r0*NBB(u, rep)/sum;
            s3 =  s2 + r1*NAA(u, rep)/ sum;
            s4 =  s3 + r2*NAA(u, rep) / sum;
            s5 = s4 + r3*NBB(u, rep) / sum;
            s6 =  s5 + r4*NBB(u, rep) / sum;
            
            tau = -(1/sum) *(log(rand(1)));
            Pr = rand(1);

            %Reaction event A -> A
            if Pr <= s3 && Pr > s2
                NAA(u+1, rep) = NAA(u, rep) ;
                NBB(u+1, rep) = NBB(u, rep);
            end
            
            
           %Reaction event B -> B
            if Pr <= s5 && Pr > s4
                NAA(u+1, rep) = NAA(u, rep);
                NBB(u+1, rep) = NBB(u, rep);
            end 
            
            
            
           %Reaction event A -> B
            if Pr <= s4 && Pr > s3
                NAA(u+1, rep) = NAA(u, rep) - 1;
                NBB(u+1, rep) = NBB(u, rep) + 1;
            end   
	    %Reaction event B -> A
            if Pr <= s6 && Pr > s5 && NBB(u, rep) > 0
                NAA(u+1, rep) = NAA(u, rep) + 1;
                NBB(u+1, rep) = NBB(u, rep) - 1;
            end
                     
        end
         
        NA(1, rep) = NAA(u,rep); 
        NB(1, rep) = NBB(u,rep); 
  
        
	adapt = [];
	generation = [];
	mod_rep = nrep; %Modified counter if no C class is created
	oldit = 0;
	%Loop over time
        for it = 1:itime
            sum = r0*(NA(it, rep) + NB(it, rep)) + ...
                    r1*NA(it,rep) + r2*NA(it,rep) + ...
                    r3*NB(it,rep) + r4*NB(it,rep) + r5*NB(it,rep); 
            % Reaction propensities
            s1 = r0*NA(it, rep)/sum;
            s2 = s1 + r0*NB(it, rep)/sum;
            s3 =  s2 + r1*NA(it, rep)/ sum;
            s4 =  s3 + r2*NA(it, rep) / sum;
            s5 = s4 + r3*NB(it, rep) / sum;
            s6 =  s5 + r4*NB(it, rep) / sum;
            s7 =  s6 + r5*NB(it, rep) / sum;

            tau = -(1/sum) *(log(rand(1)));
            Pr = rand(1);

	    %Reaction event A -> B
            if Pr <= s4 && Pr > s3
                NA(it+1, rep) = NA(it, rep) - 1;
                NB(it+1, rep) = NB(it, rep) + 1;
                NC(it+1, rep) = NC(it, rep);
	    %Reaction event B -> A
            elseif Pr <= s6 && Pr > s5 && NB(it, rep) > 0
                NA(it+1, rep) = NA(it, rep) + 1;
                NB(it+1, rep) = NB(it, rep) - 1;
                NC(it+1, rep) = NC(it, rep);
	    %Reaction event B -> C
            elseif Pr <= s7 && Pr >= s6 && NB(it, rep) > 0
                NA(it+1, rep) = NA(it, rep);
                NB(it+1, rep) = NB(it, rep) - 1;
                NC(it+1, rep) = NC(it, rep) + 1;
		tm = tim(it, rep) + tau;
		adapt = [adapt tm];
		generation = [generation it-oldit];
		oldit = it;
	    %Everything else
            else
                NA(it+1, rep) = NA(it, rep);
                NB(it+1, rep) = NB(it, rep);
                NC(it+1, rep) = NC(it, rep);
            end

            tim(it+1, rep) = tim(it, rep) + tau;
        end %End of iterattion over time
	%Find adapation time from last two (steady state)
	if length(adapt) > 1
		adtime(qi, rep) = adapt(end) - adapt(end-1);
		ngeneration(qi, rep) = generation(end);
		firstgeneration(qi, rep) = generation(1);
	elseif length(adapt) == 1
		adtime(qi, rep) = adapt;
		ngeneration(qi, rep) = generation;
		firstgeneration(qi, rep) = generation;
	else
		adtime(qi, rep) = 0;
		ngeneration(qi, rep) = 0;
		mod_rep = mod_rep - 1;
	end
	sumadtime = sumadtime + adtime(qi, rep);
	sumgeneration = sumgeneration + ngeneration(qi, rep)*NB(itime, rep)/N;
	sumfirstgeneration = sumfirstgeneration + firstgeneration(qi, rep);
    end % End of loop over ensembles
    %Calculate cumulative adaptation time averaged over ensembles
    cumadtime(qi) = sumadtime/mod_rep;
    cumgeneration(qi) = sumgeneration/mod_rep;
    mean_firstgeneration(qi) = sumfirstgeneration/mod_rep;
    %Average numbers of A and B over ensembles
    MeanA(:, qi) = mean(NA, 2);
    MeanB(:, qi) = mean(NB, 2);
    MeanC(:, qi) = mean(NC, 2);
end %End of loop over robustness values


%Plot cumulative adaptation time
f1 = figure;
plot(qvec, cumadtime);
xlabel('Robustness');
ylabel('Mean adaptation time');
saveas(f1, 'Adaptation_Time.jpeg');

f2 = figure;
plot(qvec, mean_firstgeneration);
xlabel('Robustness');
ylabel('Mean generations');
saveas(f2, 'Mean_Generation.jpeg');
%Plot cumulative adaptation generation
f3 = figure;
plot(qvec, cumgeneration);
xlabel('Robustness');
ylabel('Mean adaptation generation');
saveas(f3, 'Adaptation_Generation.jpeg');
%set(gca, 'YScale','log');


%Plot individual numbers
% %%
% mean_t = mean(tim, 2);
% f4 = figure;
% plot(mean_t, MeanA(:, 10),'r','linewidth',2);
% hold on 
% txt1 = {'A'};
% text(10,2000,txt1,'Color','red','FontSize',15)
% 
% plot(mean_t, MeanB(:, 10),'g','linewidth',2);
% hold on
% txt2 = {'B'};
% text(10,3300,txt2,'Color','green','FontSize',15)
% 
% plot( mean_t, MeanC(:, 10),'b','linewidth',2);
% hold on
% txt3 = {'C'};
% text(10,300,txt3,'Color','blue','FontSize',15)
% 
% plot(mean_t,MeanA(:, 10) + MeanB(:, 10) + MeanC(:, 10),'k','linewidth',2);
% hold on
% txt4 = {'Total'};
% text(10,5200,txt4,'Color','black','FontSize',15)
% 
% title('For K/P = 70/100','FontSize',15)
% xlabel('Time');
% ylabel('Population');
% set(gca,'FontSize',15)
