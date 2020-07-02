%LP solver for multiple given set control
%in this case, distributed policy can be used to get the same result as
%centralized with reduced computation
%10-node case for sandy
clear;close all;
%% parameter
N=10;
%only given set of nodes have controller
aflag=ones(N,1);
aflag(9:10)=0;%aflag(2)=1;
%aflag(3)=1;%aflag(4)=1;%we can only control node 1 and 2
%aflag(6)=0;aflag(5)=0;aflag(4)=0;aflag(3)=0;
%the probablity and load parameter
beta=1;
%reward are given based on the real flow
%R=ones(N,1).*[1,2,3,4]';%reward is 5
% R=ones(N,1).*[1,2,3,4,5,6]';
% R=R';
% R=100*R/sum(R);

cost=ones(N+1,1); %cost of taking action
%cost=-cost/100;
cost=-2*cost;
cost(N+1)=0;


%% REAL 
%source: http://web.mta.info/nyct/facts/ridership/ridership_sub_annual.htm
%R(1)=5222096;
R(1)=11492780;
R(2)=5340581;
R(3)=1855650;
R(4)=2209212;
%R(6)=2344141;
%R(7)=4021679;
R(5)=11792114;
R(6)=10650508;
R(7)=15364366;
R(8)=1280063;
R(9)=7034940;
R(10)=7990789;
%normalize
R=100*R/sum(R);
%%


%other parameter
gama=0.9;
P_prior=0.01; %prior failure probability

%initilize
Prob=zeros(N,1);
position=zeros(N,1);
ind=0;
inde=0;
kkke=0;
index2=0;
kkk=0;

%% G is the topology
G=eye(N,N);
Gtotal=zeros(N,N);
%%
%ring
% G(1,2)=1;
% G(2,3)=1;
% G(3,4)=1;
% G(4,5)=1;
% G(5,6)=1;
% G(6,1)=1;
% 
% G(2,1)=1;
% G(3,2)=1;
% G(4,3)=1;
% G(1,4)=1;

% %star order1
% G(2,1)=1;
% G(3,1)=1;
% G(4,1)=1;

% G(1,2)=1;
% G(1,3)=1;
% G(1,4)=1;

%star order2
% G(1,4)=1;
% G(2,4)=1;
% G(3,4)=1;

% G(4,1)=1;
% G(4,2)=1;
% G(4,3)=1;
%star order3
% G(1,2)=1;
% G(3,2)=1;
% G(4,2)=1;

% G(3,1)=1;
% G(3,2)=1;
% G(3,4)=1;


%line oneway or twoway
% G(1,2)=1;
% G(2,3)=1;
% %G(3,4)=1;
% 
% G(2,1)=1;
% G(3,2)=1;
% %G(4,3)=1;

%random
% G(2,3)=1;
% G(3,1)=1;
% G(2,4)=1;
% G(3,4)=1;

%% REAl
G(1,2)=1;G(1,4)=1;
G(2,3)=1;G(2,10)=1;
G(3,9)=1;
G(4,5)=1;
G(5,6)=1;
G(6,7)=1;
G(7,8)=1;G(7,10)=1;
G(8,9)=1;
G(9,10)=1;

%when G is bi-direction

    G=G | G';
%%
cflag=ones(1,N);%flag of usage of C0{i}

degree=zeros(1,N);
degreeout=zeros(1,N);
Gnew=zeros(N,N);

for i=1:N
    degree(i)=sum(G(:,i));
    degreeout(i)=sum(G(i,:))-1;
end

%construct GG showing the scope without ei
GG=zeros(N,N);
for i=1:N
    for j=1:N
if G(i,j)==1
    GG(:,i)=GG(:,i) | G(:,j);
end
    end
end


%% physical setting
P_piror=0.01;%this should a rare event. without generalization, this P is same for all node.
%C=[70 70 70 70];
%C=[60 70 75]; %capability of five nodes, could define a C_average and make it random for each
%rng('default');
%P_old=C.*rand(1,5); %the net injection of each node [0 C(i)]  a+(b-a)*rand(m,n)
%P_old=ones(N,1);
P_old=R;
C=1.2*P_old;

rng(0);
%weight of basis function
w = sdpvar(N+1,1); %the last one is w0

%% set up table of failure probability P(i'| i,\Omega_i) 
%basis function's demension should be N*2^(1+|neibor|)
C0 = cell(N,1);
for i=1:N
    C0{i}=rand(1,2^degree(i));
end

%for each action we generate new constraints
%a=N+1 is the default action
Constraints=[];

%every nodes have control, use brute of a=1:2^N or factor action space also
%we can have given set(#m) of node have action 2^m. Make a=vector
%have a given number of node to control, how to choose?
for a=1:2^sum(aflag)
    ahln=de2bi(a-1,sum(aflag));%1->100
        ahln=fliplr(ahln); %1->001
        %put ahln the length as aflag so that we can compare later in bit
        ahlnflag=zeros(N,1);
        yyy=find(aflag==1);
        ahlnflag(yyy)=ahln;
        
tic
%% put the corresonding probablity inside C{i} considering the load-probablity relation.
for i=1:N
    for j=1:2^degree(i)
        
        hln=de2bi(j-1,degree(i));%1->100
        hln=fliplr(hln); %1->001
        %find G(i,i)'s rank in all 1 in the column
        geshu=find(G(1:i,i)==1);
        geshu3=find(G(:,i)==1);
        position(i)=length(geshu);
        geshu2=find(hln==0);
        
        if aflag(i)==1 && 1==ahlnflag(i) %first controllable  %then, we do control/repair i
        C0{i}(j)=1;
        elseif position(i)~=0 && hln(position(i))==0 %i fail and not repair i. xi is 0 here.
        C0{i}(j)=0;
        else %i good but other fails
            %should for extra
            P=P_old;
            for l=1:length(find(hln==0)); %for this else, xi can only be 1
                xxx=geshu3(geshu2(l));
            P(i)=P(i)+P_old(xxx)/degreeout(xxx);
            end
            %compute Probablity
            if P(i)<=C(i)
                Prob(i)=1-P_prior; %Prob is the probability of true
            else
                Prob(i)=1+(1/P_prior-1)*exp(-beta*(P(i)-C(i)));
                Prob(i)=1/Prob(i);
                
                Prob(i)=1-Prob(i); %Prob is the probability of true
            end
            %put h=probability into C0
            C0{i}(j)=Prob(i);
           
        end
        
        %c=gama*g-h
            C0{i}(j)=gama*C0{i}(j);
            if hln(position(i))==1
            C0{i}(j)=C0{i}(j)-1;
            end
            
    end
end

for i=1:N
    C0{i}=C0{i}*w(i);
end

%% variable elimination in LP constraints
C1=cell(N,1);
E1=cell(N,1);
Gtotal(:,1)=GG(:,1);
for i=1:N %first define all the variable. Remanber E1 are all varibles.
   %first elimination do not need to consider e, only consider GG. Gnew=GG(,1) 
   C1{i}(2,:)=sdpvar(1,2^(sum(Gtotal(:,i))));
   C1{i}(1,:)=zeros(2^(sum(Gtotal(:,i))),1);
   
       %generate new LP variable each step
       Gtotal(i,i)=0;
   E1{i} = sdpvar(2^(sum(Gtotal(:,i))),1);
   
   %fill in ei %the elimination variable is always at the first
       C1{i}(2,1:2^(sum(Gtotal(:,i))))=E1{i};
       C1{i}(2,2^(sum(Gtotal(:,i)))+1:2^(sum(Gtotal(:,i))+1))=E1{i};
   
   %Gnew show the scope of e1
   %each time eliminate a variable xj, search for weather ei involves xj
  if i+1<=N 
   index_e=find(Gtotal(i+1,:)==1);
   for j=index_e
   Gtotal(:,i+1)=zeros(N,1) | Gtotal(:,j);
   %Gnew(:,j) has been used, so should be zero now
   Gtotal(:,j)=0;
   end
   Gtotal(:,i+1)=Gtotal(:,i+1) | GG(:,i+1);
   Gtotal(i,i+1)=0;
  end
     
   GG(i,:)=0;

end
%recover the oringial GG
GG=zeros(N,N);
for i=1:N
    for j=1:N
if G(i,j)==1
    GG(:,i)=GG(:,i) | G(:,j);
end
    end
end
%


Gtotal(:,1)=GG(:,1);
%elimination order from 1 to N
for i=1:N
    
Gnew(:,i)=Gtotal(:,i);
Gnew(i,i)=0;
    
   %fill in the form
   for j=1:2^(sum(Gtotal(:,i)))
       %fill in Reward 
        hln=de2bi(j-1,sum(Gtotal(:,i)));%1->100
        hln=fliplr(hln); %1->001
       if hln(1)==1 %??the variable to eliminte is always the first??
       C1{i}(1,j)=C1{i}(1,j)+R(i);
       end
       %fill in cost
       %C1{i}(1,j)=C1{i}(1,j)+cost(1:N)'*ahlnflag;
       %fill in C0{i}s of the i
       if cflag(i)~=0
           GGG=Gtotal(:,i)+G(:,i); %should be Gnew here, which show the scope of the whole equation
           GGG(GGG==0)=[];
           pos=find(GGG==2);
       hln1=hln(pos); 
       hln1=fliplr(hln1);
       hln1=bi2de(hln1)+1; %the index for the C0{i}
       C1{i}(1,j)=C1{i}(1,j)+C0{i}(hln1);
       end
       %fill in others
            for k=i+1:N %start from i will be enough?
                if G(i,k)==1 && cflag(k)~=0
                %abstract the index  
                 index=0;
                for kk=1:N
                    if Gtotal(kk,i)==1 && G(kk,k)==1 %GG->Gnew
                        index=index+1;
                        flag(index)=index; %flag will be the index of abstraction
                    elseif Gtotal(kk,i)+G(kk,k)==1
                     index=index+1;
                    end
                    
                end
       flag(flag==0)=[];
       hln2=hln(flag); 
       hln2=fliplr(hln2);
       hln2=bi2de(hln2)+1; %the index for the C0{i}
       C1{i}(1,j)=C1{i}(1,j)+C0{k}(hln2);
              flag=[];  
             ind=ind+1;
              kkk(ind)=k;
                end
            end
         %fill in all the e1 term containing the variable i

   index_e=find(Gnew(i,:)==1);
   for jj=index_e
   Gtotal(:,i+1)=zeros(N,1) | Gnew(:,jj);
   %Gnew(:,j) has been used, so should be zero now
   
         for jjj=i:N
                 if Gtotal(jjj,i)==1 && Gnew(jjj,jj)==1
                      index2=index2+1;
                        flag2(index2)=index2; %flag will be the index of abstraction                 
                 elseif Gtotal(jjj,i)+Gnew(jjj,jj)==1
                     index2=index2+1;
                 end    
         end
                 flag2(flag2==0)=[];
       hln3=hln(flag2); 
       hln3=fliplr(hln3);
       hln3=bi2de(hln3)+1; %the index for the C0{i}
       C1{i}(1,j)=C1{i}(1,j)+E1{jj}(hln3); %%jj is the index of new LP variable e
       flag2=[];
       index2=0;
       
       inde=inde+1;
       kkke(inde)=jj;
   end
  if i+1<=N
   Gtotal(:,i+1)=GG(:,i+1) | Gtotal(:,i);
   Gtotal(i,i+1)=0; %xi is eliminated
  end
   end  
     if kkk~=0
      cflag(kkk)=0;%mark the C0{i} after usage 
     end
        cflag(i)=0;%mark the C0{i} after usage
        
    %put Gnew=0 after usage 
    if kkke~=0
    Gnew(:,kkke)=0;
    end
    
   GG(i,:)=0;
end


%% solve LP by yalmip


%final one
Constraints=Constraints + (E1{N}-0.1*w(N+1)+cost(1:N)'*ahlnflag<=0);

%constrains generated at each step
for m=1:N
    for n=1:size(C1{m},2)
Constraints = Constraints + (C1{m}(2,n)>=C1{m}(1,n));
    end
end

%Constraints = Constraints + [norm(w,2)<=5];
%GG return to the previous one because we will have another a
%recover the oringial GG
GG=zeros(N,N);
for i=1:N
    for j=1:N
if G(i,j)==1
    GG(:,i)=GG(:,i) | G(:,j);
end
    end
end

%recover all flags
cflag=ones(1,N);%flag of usage of C0{i}
kkk=0;
ind=0;
index2=0;
inde=0;
kkke=0;
Gtotal=zeros(N,N);
toc
end
%Constraints=Constraints+[w>=0];
Objective = sum(w)+w(N+1); %default is min, our objective is to min
%.*[1,2,3,4]'


%options = sdpsettings('verbose',1,'solver','quadprog','quadprog.maxiter',100);
options = sdpsettings('solver','gurobi');
% Solve the problem
sol = optimize(Constraints,Objective,options);
% Analyze error flags
if sol.problem == 0
 % Extract and display value
 solution = value(w)
else
 display('Hmm, something went wrong!');
 sol.info
 yalmiperror(sol.problem)
end

%% distributed way to compute policy
policy=cell(N,1); %state to action_1,...action_N

 for actionindex=1:N %compute a_actionindex one by one
        gindex=find(G(:,actionindex)==1);
        policy{actionindex}=zeros(2^length(gindex),length(gindex)+1);
     for ii=1:2^length(gindex) %for each state 
         

        hln=de2bi(ii-1,length(gindex));%1->100
        hln=fliplr(hln); %1->001   
        
            policy{actionindex}(ii,1:length(gindex))=hln;
        
  
        
        for l=1:length(hln)
            if hln(l)==0
       P(actionindex)=P(actionindex)+P_old(gindex(l))/degreeout(gindex(l));
            end
        end
        
          %compute Probablity
            if P(actionindex)<=C(actionindex)
                Prob(actionindex)=1-P_prior; %Prob is the probability of true
            else
                Prob(actionindex)=1+(1/P_prior-1)*exp(-beta*(P(actionindex)-C(actionindex)));
                Prob(actionindex)=1/Prob(actionindex);
                
                Prob(actionindex)=1-Prob(actionindex); %Prob is the probability of true
            end
        
             %if a_i=1
        %compare=gama*value(w(actionindex))+cost(actionindex);
        %if a_i=0
 compare=(1-Prob(actionindex))*gama*value(w(actionindex))+cost(actionindex);
 if compare>0  %a_i=1
  policy{actionindex}(ii,1+length(gindex))=1;
 end

    end
 end

%% distributed way to compute policy
tic
policy=cell(N,1); %state to action_1,...action_N

 for actionindex=1:N %compute a_actionindex one by one
        gindex=find(G(:,actionindex)==1);
        policy{actionindex}=zeros(2^length(gindex),length(gindex)+1);
     for ii=1:2^length(gindex) %for each state 
         

        hln=de2bi(ii-1,length(gindex));%1->100
        hln=fliplr(hln); %1->001   
        
            policy{actionindex}(ii,1:length(gindex))=hln;
        
  %if node itself is down, then g=0
        PPindex=find(gindex==actionindex);
        if hln(PPindex)==0
            Prob(actionindex)=0;
        else
            
        P=P_old;
        
        
        for l=1:length(hln)
            if hln(l)==0 && l~=PPindex
       P(actionindex)=P(actionindex)+P_old(gindex(l))/degreeout(gindex(l));
            end
        end
        
          %compute Probablity
            if P(actionindex)<=C(actionindex)
                Prob(actionindex)=1-P_prior; %Prob is the probability of true
            else
                Prob(actionindex)=1+(1/P_prior-1)*exp(-beta*(P(actionindex)-C(actionindex)));
                Prob(actionindex)=1/Prob(actionindex);
                
                Prob(actionindex)=1-Prob(actionindex); %Prob is the probability of true
            end
        
             %if a_i=1
        %compare=gama*value(w(actionindex))+cost(actionindex);
        %if a_i=0
         end
 compare=(1-Prob(actionindex))*gama*value(w(actionindex))+cost(actionindex);
 if compare>0  %a_i=1
  policy{actionindex}(ii,1+length(gindex))=1;
 end

    end
end
toc
 
% %% optimal policy
% tic
% policy=zeros(2^N,N+sum(aflag));
% Reward=zeros(2^sum(aflag),1);
% Reward=Reward+value(w(N+1));%reward with the constant w0
% 
% for ii=2^N:-1:1 %for each state
%     for action=1:2^sum(aflag) %find the best action
%         
%     ahln=de2bi(action-1,sum(aflag));%1->100
%         ahln=fliplr(ahln); %1->001
%         %put ahln the length as aflag so that we can compare later in bit
%         ahlnflag=zeros(N,1);
%         yyy=find(aflag==1);
%         ahlnflag(yyy)=ahln;
%         
%         %% recompute C0
%         
% for i=1:N
%     for j=1:2^degree(i)
%         
%         hln=de2bi(j-1,degree(i));%1->100
%         hln=fliplr(hln); %1->001
%         %find G(i,i)'s rank in all 1 in the column
%         geshu=find(G(1:i,i)==1);
%         geshu3=find(G(:,i)==1);
%         position(i)=length(geshu);
%         geshu2=find(hln==0);
%         
%         if aflag(i)==1 && 1==ahlnflag(i) %first controllable  %then, we do control/repair i
%         C0{i}(j)=1;
%         elseif position(i)~=0 && hln(position(i))==0 %i fail and not repair i. xi is 0 here.
%         C0{i}(j)=0;
%         else %i good but other fails
%             %should for extra
%             P=P_old;
%             for l=1:length(find(hln==0)); %for this else, xi can only be 1
%                 xxx=geshu3(geshu2(l));
%             P(i)=P(i)+P_old(xxx)/degreeout(xxx);
%             end
%           %compute Probablity
%             if P(i)<=C(i)
%                 Prob(i)=1-P_prior; %Prob is the probability of true
%             else
%                 Prob(i)=1+(1/P_prior-1)*exp(-beta*(P(i)-C(i)));
%                 Prob(i)=1/Prob(i);
%                 
%                 Prob(i)=1-Prob(i); %Prob is the probability of true
%             end
%             %put h=probability into C0
%             C0{i}(j)=Prob(i);
%             
%         end
%         
%         %c=gama*g-h
%             C0{i}(j)=gama*C0{i}(j);
%             if hln(position(i))==1
%             C0{i}(j)=C0{i}(j)-1;
%             end
%             
%     end
% end
% 
% for i=1:N
%     C0{i}=C0{i}*value(w(i));
% end
%         
%         
%         
%       %%  
%       
%         hln=de2bi(ii-1,N);%1->100
%         hln=fliplr(hln); %1->001
%         %add reward
%         Reward(action)=Reward(action)+hln*R';
%         Reward(action)=Reward(action)+cost(1:N)'*ahlnflag;
%         
%             for kk=1:N %addup all C0{k}
%                 indexa=hln(G(:,kk)==1);
%                 indexa=fliplr(indexa);
%                 indexa=bi2de(indexa)+1;%1->100
%                Reward(action)=Reward(action)+C0{kk}(indexa); 
%                %fill in cost
%             end
%         
%     end
%     %Reward
%            [maxvalue,policyindex]=max(Reward);
%            
%         hlnpolicy=de2bi(policyindex-1,sum(aflag));%1->100
%         hlnpolicy=fliplr(hlnpolicy); %1->001
%            
%            policy(ii,N+1:N+sum(aflag))=hlnpolicy;
%            policy(ii,1:N)=hln;
% %restart for a new one           
% Reward=zeros(2^sum(aflag),1);
% Reward=Reward+value(w(N+1));%reward with the constant w0
% end
% 
% 
% %% result procession
% Objective
% 
% toc
sum(value(w))
