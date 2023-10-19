% By Abdollah Ghaffari Sheshjavani
% the code is for the paper in title "Content Caching in Shared Medium Networks with Non-Uniform and User-Dependent Demands" IEEE ICC, 2022

clear all
Policy=1; % 1=duplicate select files is consider for not cached files(uncoded requests)   else=  duplicate select files is not considered for not cached files
ShowDetailResults=false;
Time=2000; % rounds of runs
clustering=0  % 0=without clustering    1=with clustering
clusterCachePartition=2 % 1= all cache allocated for local       2= cache partition into two part a)global b)local 
Method=1  %1=our hybrid   2=purecoded  3=pureUncoded
HuristicType=1% 0= generalPopularity(independent)        1= Huristic1(dependent)
% input parameters 
K=10; % number of SBS
N=1000;%1000; % number of video files
M=50;%100; % SBS cache size
GLBeta=1/M;

% these bellow parameters used for evaluating different distribution of users inside SBSs (here is 10 SBSs)
Z1=[1 1 1 1 1 5 15 20 25 30 ];
Z2=[0 2 2 3 7 11 14 16 20 25 ];
Z3=[0 2 4 6 9 11 14 16 18 20 ];
Z4=[1 3 5 7 9 11 13 15 17 19 ];
Z5=[2 4 6 8 9 11 12 14 16 18 ];
Z6=[3 5 7 9 9 11 11 13 15 17 ];
Z7=[4 6 9 9 9 10 11 12 14 16 ];
Z8=[5 7 9 9 9 10 11 12 13 15 ];
Z9=[6 8 9 9 9 10 11 12 12 14 ];
Z10=[8 9 9 9 9 10 11 11 12 12 ];
Z11=[10 10 10 10 10 10 10 10 10 10];
Z12=[1 1 1 1 1 1 1 1 1 1 ];

Z=Z11; % selecting distribution of users (the Number of users in the coverage of each SBS) for simulation and analysis

%////////////////////////////////////////////popularity set//////////////
% popularityArray=[ [0.3  0.2  0.5  0  ]   % this array is used only for evaluation optimal placement in tiny scale
%                  [0.2  0.3  0.5  0  ]
%                  [0.3  0.2  0  0.5  ]
%                  [0.2  0.3  0  0.5  ]]
group_count = 5;   % number of different popularity groups of contents
ZipfParameter = [1, 1, 1, 1, 1];   % popularity parameter for each group of contents 
group_interest = [0.5, 0.5, 0.0];  % percent of request for local popular group of contents, percent of request for global popular contents, percent of request for other group
%GlobalPopularContentNumber=200;
GlobalPopularContentNumber=0.2*N;
%//////////////////////////////////////////END popularity set//////////////

resultsize=10;
resultIndex=1;

Analyticalresult=zeros(1,resultsize);
Simulationresult=zeros(1,resultsize);
StandardDeviationR=zeros(1,resultsize);
TimeComplexity=zeros(1,resultsize);
progress = waitbar(0,'Please wait...');
evaluation_param0=[0.02 0.04 0.06 0.08 0.1 0.2 0.4 0.5  1]; % its only for GLBeta 

for evaluation_param= 0.5:0.1:1.5  %1:1:10
   % ######## from below parameter uncommented only one for evaluating corresponding feature (and change the value of evaluation_param in above for to valid values) other parameters are set above ##############
    %GLBeta = evaluation_param0(evaluation_param)
    %N=evaluation_param
    %M=0.05*N
    %GlobalPopularContentNumber=.2*N
    %Z=Z12*evaluation_param
    %M=evaluation_param
    ZipfParameter = [evaluation_param, evaluation_param, evaluation_param, evaluation_param, evaluation_param]
    %group_interest = [evaluation_param, 1-evaluation_param, 0.0]
   %############# end of selecting evaluating feature ###################### 

    popularityArray = groupingPreferenceMaker(N, K, ZipfParameter, group_count, group_interest);
    tic;
    
    placement=clusterHuristicFunction(Z,K,N,M,popularityArray,GlobalPopularContentNumber,Method,clustering,HuristicType,clusterCachePartition);% 1 = hybrid coded-uncoded, 2=pure coded, 3=pure uncoded
   
    TimeComplexity(1,resultIndex)=toc;
    Analyticalresult(1,resultIndex)=placement(K+1,1);
    %placement
    
%///////////////////////////////////////////////////////END Analytical/////////////////////////////////////////////
%///////////////////////////////////////////////////////Start Simuation/////////////////////////////////////////////
    Zmax=max(Z);
    Allusers_count=sum(Z);
    Zvar=0;
    for i=1:K
       Zvar=Zvar+(Z(i)-(Allusers_count/K))^2;
    end
    Zvar=Zvar/(K-1);
    ZstandardDev=sqrt(Zvar);
    Zvar
    allratecoded=0;
    allratecodedFormula=0;
    allrateuncoded=0;
    allrate=0;
    allrateFormula=0;
    sumprobe=zeros(K,N);
    for i=1:K
        for idp=1:N
            sumprobe(i,idp)=sum(popularityArray(i,1:idp));
        end
    end
    %sumprobe
     %/////////////////////////////////////////// finding cluster Informations ////////////
    if clustering==0
        allclusters=0;
    else
        allclusters=placement(K+1,3);
    end
    MG=placement(K+1,2);
    ML=M-MG;
    M1=zeros(1,allclusters+1);
    NL=zeros(1,allclusters+1); 
    Kcluster=zeros(1,allclusters+1);
    cluster_cells=zeros(allclusters+1,K); 
    
    for clusterindex=1:allclusters+1 
        if clusterindex<=allclusters    
            Ktemp=0; 
            for content_index=1:N
                selected=false;
                for cell_index=1:K
                    if placement(cell_index,content_index)==-clusterindex
                        if Kcluster(1,clusterindex)==0
                            Ktemp=Ktemp+1;
                            cluster_cells(clusterindex,Ktemp)=cell_index;
                        end
                        if selected==false
                             NL(1,clusterindex)= NL(1,clusterindex)+1; % contents that cached with coded scheme are same for all cells
                             selected=true;
                        end
                    elseif placement(cell_index,content_index)==clusterindex+2                  
                        M1(1,clusterindex)= M1(1,clusterindex)+1; % contents that cached with uncoded scheme may be not same for all cell           
                    end
                end
                Kcluster(1,clusterindex)=Ktemp;
            end
            M1(1,clusterindex)=M1(1,clusterindex)/Kcluster(1,clusterindex);% we count these contents in all cell caches 
            M2=ML-M1(1,clusterindex);
        elseif clusterindex==allclusters+1 
            Kcluster(1,clusterindex)=K;
            for cell_index=1:K
               cluster_cells(clusterindex,cell_index)=cell_index;
            end
            for content_index=1:N
                if placement(1,content_index)==1
                    M1(1,clusterindex)=M1(1,clusterindex)+1; % M1 Global
                elseif(placement(1,content_index)==2)
                    NL(1,clusterindex)=NL(1,clusterindex)+1;
                end
            end 
         
        end
    end
    %cluster_cells
     %///////////////////////////////////////// END finding cluste Informations////////////
     RateResults=zeros(Time);
     for t=1 : Time
        progress = waitbar(((resultIndex-1)*Time+ t)/(Time*resultsize),progress,...
        ['All Simulation Progress =',num2str(((resultIndex-1)*Time+ t)/(Time*resultsize)*100,'%4.1f'),'%' ]);
        SBSRequests=zeros(K,Zmax);
        SBSGlobalcodedQueu=zeros(K,Zmax);
        GlobalrequestQueusize=zeros(1,K);
        SBSLocalcodedQueu=zeros(K,Zmax);
        LocalrequestQueusize=zeros(1,K);
        uncodedRequest=zeros(1,K);
        cacheMiss=zeros(1,Zmax);
        rateuncoded=0;
    % \\\\\\\\\\\\\\\\\\\\\this is request generator based-on popularity distribution
        for i=1 : K  
                for idx2=1: Z(i)
                   temp=rand;
                   for idx=1: N
                       if temp<=sumprobe(i,idx);
                            SBSRequests(i,idx2)=idx;
                            break;
                        end
                    end
                end
    % \\\\\\\\\\\\\\\\\\\\end of request generator based-on popularity distribution
           for j=1 : Z(i)
              SBSduplicateRequest=false;
              MBSUncodedduplicateRequest=false;
              if j>1 
                  for n=1 : j-1
                       if SBSRequests(i,j)==SBSRequests(i,n)
                           SBSduplicateRequest=true;
                       end
                  end
              end
              if i>1
                  for n=1 : i-1
                      for l=1: Z(i)
                           if SBSRequests(i,j)==SBSRequests(n,l)&& placement(n,SBSRequests(n,l))==0
                               MBSUncodedduplicateRequest=true;
                           end
                      end
                  end
              end
              if placement(i,SBSRequests(i,j))==2 &&  SBSduplicateRequest==false %SBSRequests(i,j)>M1
                  GlobalrequestQueusize(i)= GlobalrequestQueusize(i)+1;
                  SBSGlobalcodedQueu(i,GlobalrequestQueusize(i))=SBSRequests(i,j);
              elseif placement(i,SBSRequests(i,j))<0 &&  SBSduplicateRequest==false %SBSRequests(i,j)>M1
                  LocalrequestQueusize(i)= LocalrequestQueusize(i)+1;
                  SBSLocalcodedQueu(i,LocalrequestQueusize(i))=SBSRequests(i,j);
              elseif placement(i,SBSRequests(i,j))==0 %SBSRequests(i,j)>N1 
                  if Policy==1
                      if SBSduplicateRequest==false && MBSUncodedduplicateRequest==false;
                         uncodedRequest(i)=uncodedRequest(i)+1; 
                      end
                  else
                      uncodedRequest(i)=uncodedRequest(i)+1;
                  end
              end
                   %end
           end  
        end 

    %     placement
    %     SBSRequests
    %     uncodedRequest
        for p=1: K
            rateuncoded=rateuncoded+uncodedRequest(p);  
        end
    %/////////////////////////////////End computing r2\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        
        ratecoded=zeros(1,allclusters+1);
        ratecoded2=zeros(1,allclusters+1);
       
        for clusterindex=1:allclusters+1 
            M2=ML-M1(clusterindex);
            if clusterindex==allclusters+1
                M2=MG-M1(clusterindex);
            end 
            T=(Kcluster(1,clusterindex)*M2)/NL(1,clusterindex);
            a=zeros(0);
            for f=1 :Kcluster(1,clusterindex)
                a=cat(1,a,f);
            end
            if M2>0
                ComKT=combination(Kcluster(1,clusterindex),T);
                CombKT1=0;
                if Kcluster(1,clusterindex)>T
                    CombKT1=combination(Kcluster(1,clusterindex),T+1);
                end
                slices=nchoosek(a,T);
                sizeTchooseK=size(slices(:,1));
                allXOR=nchoosek(a,T+1);
      %//////////////////Coded cache placment phase//////////////////                
                 if ShowDetailResults ==true
                    for l=1:Kcluster(1,clusterindex)
                        index=1;
                        for j=1 : sizeTchooseK
                            for i=1 :T
                                if (slices(j,i)==l)
                                   for h=1:T
                                      SBSCodedcache(index,h,l)=slices(j,h);
                                   end
                                   index=index+1;
                                end
                            end
                        end
                    end
                    size(slices)
                    slices
                    SBSCodedcache
                end
       %/////////////////// end of coded cache placement Phase//////////////     
               
                for q=1: Zmax
                   ki=0;
                   index=0;
                   for p=1 : Kcluster(1,clusterindex)
                       if  clusterindex <=allclusters 
                           if LocalrequestQueusize(cluster_cells(clusterindex,p))>0
                                ki=ki+1;
                                LocalrequestQueusize(cluster_cells(clusterindex,p))=LocalrequestQueusize(cluster_cells(clusterindex,p))-1;
                           end
                       elseif clusterindex ==allclusters+1 
                           if GlobalrequestQueusize(cluster_cells(clusterindex,p))>0
                                ki=ki+1;
                                GlobalrequestQueusize(cluster_cells(clusterindex,p))=GlobalrequestQueusize(cluster_cells(clusterindex,p))-1;
                           end
                       end
                   end
                   if ki>0
                    % calculating rate
                        for i=1: size(allXOR(:,1))
                           flag=false;
                           for j=1: T+1
                               if  clusterindex <=allclusters 
                                   if SBSLocalcodedQueu(cluster_cells(clusterindex,allXOR(i,j)),q)> 0
                                       flag=true;
                                   end
                               elseif clusterindex ==allclusters+1
                                   if SBSGlobalcodedQueu(cluster_cells(clusterindex,allXOR(i,j)),q)> 0
                                       flag=true;
                                   end
                               end
                           end 
                           if flag==true
                              index=index+1;
                              if ShowDetailResults ==true
                                sendSlice(index,:)=allXOR(i,:);
                              end
                              ratecoded(1,clusterindex)=ratecoded(1,clusterindex)+(1/sizeTchooseK(1));
                           end
                        end
                        if ShowDetailResults ==true
                           showSend=zeros(index,T+1);
                           for i=1: index 
                             showSend(i,:)=sendSlice(i,:);
                           end
                           showSend
                        end 
                       % The end of calculating rate 

                   % calculating rate with use formula
                    notsending=0;
                    if Kcluster(1,clusterindex)-ki>=T+1
                        notsending=combination(Kcluster(1,clusterindex)-ki,T+1);
                    end
                    %EKi(q)=EKi(q)+ki;
                    ratecoded2(1,clusterindex)=ratecoded2(1,clusterindex)+min ((CombKT1-notsending)/ComKT,(ki-((ki*M2)/(NL(1,clusterindex)-M1(clusterindex))))*((NL(1,clusterindex)-M1(clusterindex))/ki));% CombKT1 is combination(K,T+1)
                  % The end of calculating rate with use formula
                   end
                end
            end
            if ShowDetailResults ==true
                SBSRequests
                SBScodedQueu
            end
            allratecoded=allratecoded+ratecoded(1,clusterindex);
            allratecodedFormula=allratecodedFormula+ratecoded2(1,clusterindex);
        end
        
        RateResults(t)=sum(ratecoded)+rateuncoded;
        allrateuncoded=allrateuncoded+rateuncoded;
        allrate=allratecoded+allrateuncoded;
        allrateFormula=allratecodedFormula+allrateuncoded;
        
    
    end
% calculating Average of results
    allratecoded=allratecoded/Time;
    allratecodedFormula=allratecodedFormula/Time;
    allrateuncoded=allrateuncoded/Time;
    allrate=allrate/Time;
    allrateFormula=allrateFormula/Time;
    allrateuncoded
    allratecoded
    allrate
    for i=1:Time
        StandardDeviationR(1,resultIndex)=StandardDeviationR(1,resultIndex)+(RateResults(i)-allrate)^2; 
    end
    StandardDeviationR(1,resultIndex)=sqrt(StandardDeviationR(1,resultIndex)/(Time-1));

%     allratecodedFormula
%     allrateFormula
%   
%///////////////////////////////////////////////////////END Simuation////////////////////////////////////////////
    Simulationresult(1,resultIndex)=allrate;
    resultIndex=resultIndex+1;
% lll = '======================================='
 
end
close(progress);
for i=1:resultIndex
    confidenceInt95L(i)=Simulationresult(1,i)-1.96*(StandardDeviationR(1,i)/sqrt(Time));
    confidenceInt95H(i)=Simulationresult(1,i)+1.96*(StandardDeviationR(1,i)/sqrt(Time));
end
%placement
Analyticalresult
Simulationresult
StandardDeviationR
confidenceInt95L
confidenceInt95H
TimeComplexity
lll = '======================================='