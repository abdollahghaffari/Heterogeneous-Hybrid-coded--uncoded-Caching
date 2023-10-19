% Abdollah Ghaffari sheshjavani 4/5/2020 , 19/01/1399
function out=clusterHuristicFunction(Z,K,N,M,popularityArray,GlobalPopularContentNumber,Method,clustering,HuristicType,clusterCachePartition)

%Method=1; % 1=our hybrid   2=purecoded  3=pureUncoded
% Cntr+R  Cntrl+T
out=zeros(K+1,N);
Zmax=max(Z);
%///////////////////////////////////////////////////////// Pre computation Section //////////////////////////////////////
myfactorial(1)=1;
for l=1 : K
    myfactorial(l+1)=l*myfactorial(l); % notice ==> factorial(k)=myfactorial(k+1)
end
% /////////this section compute combination(Z,j) or factorial(n)/(factorial(k)*factorial(n-k))with dynomic programing
for i=1 : K   
  temp1(i)=factorial(Z(i));
 % temp1(i)
end
combZ=zeros(K,Zmax);
for i=1 : K
    temp2=1;
    temp3=temp1(i);
    combZ(i,1)=1;
    for j=1 : Z(i)
        temp2=temp2*j;
        temp3=temp3/(Z(i)-j+1);
        combZ(i,(j+1))=temp1(i)/(temp2*temp3);%combination(Z,j);
    end
end
% ////////the end of computing combination(Z,j)

% ///////pre computation for huristic and .... placement
most_popularityindex = zeros(K, N);
for k_itrator=1:K
    [dummy,most_popularityindex(k_itrator,:)] =sort(popularityArray(k_itrator,:),'descend');
end
%most_popularityindex
sum_popularityArray = zeros(1,N);
for j=1 :K
    for i=1 : N
        sum_popularityArray(1,i) = sum_popularityArray(1,i)+ popularityArray(j,i);
    end
end
%sum_popularityArray= sum_popularityArray/K;
sorted_sum_popularityArray=zeros(1,N);
[dummy,sorted_sum_popularityArray] = sort(sum_popularityArray,'descend');


MGstart=0;
MGEnd=M;
if clustering==0     
   GlobalPopularContentNumber=N;
   MGstart=M;
elseif clusterCachePartition==1
    GlobalPopularContentNumber=0;
    MGEnd=0;   
end
globalPopularityscore=zeros(1,N);
for contentindex=1 : N
    for cellindex=1 : K
        if contentindex<= 2*GlobalPopularContentNumber %NG
            globalPopularityscore(1,most_popularityindex(cellindex,contentindex))=globalPopularityscore(1,most_popularityindex(cellindex,contentindex))+1;
        end
    end 
end
sortedglobalcontents=zeros(1,N);
for N_itrator=1:N
     [dummy,sortedglobalcontents(1,:)] =sort(globalPopularityscore(1,:),'descend');
end

GlobalPopularArray = zeros(K, N);
LocalPopularArray = zeros(K, N);

% GlobalPopularArray
% LocalPopularArray 
        

%///////////////////////////////////////////////////////// END of Pre computation Section //////////////////////////////////////

%popularityArray
bestRate=-1;
bestMG=-1;
bestClusterNumber=-1;
Allrate=0;

progress = waitbar(0,'Please wait...');
for MG=MGstart : MGEnd 
    progress = waitbar(MG/M,progress,...
    ['Analytical Progress =',num2str(MG*100/M,'%4.1f'),'%' ]);
    if MG==0 
        for cellindex=1 : K
            for contentindex=1 : N
                LocalPopularArray(cellindex,sortedglobalcontents(1,contentindex))=popularityArray(cellindex,sortedglobalcontents(1,contentindex)); 
            end     
        end 
    else
        LocalPopularArray = zeros(K, N);
        for cellindex=1 : K
            for contentindex=1 : N
              if contentindex<= GlobalPopularContentNumber%NG
                GlobalPopularArray(cellindex,sortedglobalcontents(1,contentindex))=popularityArray(cellindex,sortedglobalcontents(1,contentindex));
              else
                LocalPopularArray(cellindex,sortedglobalcontents(1,contentindex))=popularityArray(cellindex,sortedglobalcontents(1,contentindex));
              end
            end     
        end 
    end
        placementTemp=zeros(K,N);
        placementGlobal=subclusterHuristicFunction(Z,K,N,MG,GlobalPopularArray,Method,HuristicType);
        globalRate=placementGlobal(K+1,1);
        %globalRate
        bestMGClusterNumber=-1;
%         for x=1: K
%             for y=1: N
%                  placementTemp(x,y)=placementGlobal(x,y);
%             end
%         end 
        bestClusterRate=-1;
        for i=1:K
            placementClusterTemp=zeros(K,N);
            localplacement=zeros(K,N);
            clusterR1Rate=0;
            d=0;
            if M-MG>0 
                try 
                      [clusters,c,d] = kmeans(LocalPopularArray,i,'MaxIter',300,'Replicates',10); 
                      %clusters
                      %d
                catch
                  %errorclusters=i
                  continue;
                end
                for j=1:i
                   clustersize=0;
                   clusterPopularityArray=zeros(0,N);
                   clusterZ=zeros(0);
                   for k=1:K
                     if clusters(k)==j
                         clusterZ=cat(1,clusterZ,Z(k));
                         clustersize=clustersize+1;
                         clusterPopularityArray=cat(1,clusterPopularityArray,LocalPopularArray(k,:));
                     end
                   end
                   %clusterZ
                  % p=M-MG
                   %clusterPopularityArray
                  % clusterPopularityArray 
                  if clustersize>0
                       placementcluster=subclusterHuristicFunction(clusterZ,clustersize,N,M-MG,clusterPopularityArray,Method,HuristicType);
                       clusterR1Rate=clusterR1Rate+placementcluster(clustersize+2,1);
                       %clusterR1Rate
                       %placementcluster
                       SBSindex=0;
                       for x=1: K
                         if clusters(x)==j
                             SBSindex=SBSindex+1;
                              for y=1: N
                                  if placementcluster(SBSindex,y)==2 %placementcluster(SBSindex,y)==0 ||placementcluster(SBSindex,y)==2
                                      localplacement(x,y)=-j;
                                  elseif placementcluster(SBSindex,y)==1
                                      localplacement(x,y)=j+2; % for distinguish M1 of clusters
                                  end
                                  if(placementcluster(SBSindex,y)>0) % for computing r2
                                       placementClusterTemp(x,y)=1;
                                  end
                              end
                         end
                       end
                  end
                end
                %localplacement
            end
            %placementClusterTemp
         % ///////////////////////////computing r2//////////////////////////////////////////////
            localR2=0;
            for n=1 : N
                temp_r2=1;
                for c=1 : K 
                   temp_r2=temp_r2*(1-(LocalPopularArray(c,n)*(1-placementClusterTemp(c,n))))^Z(c);
                end
                temp_r2=1-temp_r2;
                %r2=r2+ (1-X(1,n))*temp_r2;%(1-(1-(((1/n)^ziph_parameter)/allziph))^sigmaZ); %(Z*K));
                localR2=localR2+temp_r2;  
            end

      % ///////////////////////////END computing r2//////////////////////////////////////////////     clusterRate=clusterR1Rate+clusterR2Rate;
            clusterRate= clusterR1Rate+ localR2; 
            %MG
            %NumberoFClusters=i
            %RateToTal=clusterRate+globalRate
            if bestClusterRate==-1 || clusterRate<bestClusterRate
                bestClusterRate=clusterRate;
                bestMGClusterNumber=i;
                %clusterR1Rate
                %localR2
                placementTemp=zeros(K,N);
                for x=1: K
                     for y=1: N   
                         placementTemp(x,y)=placementGlobal(x,y)+localplacement(x,y); 
                      end
                end
            end
            if d==0
               break;
            end
        end
        %placementTemp
        Allrate=globalRate+bestClusterRate;
        if bestRate==-1 || Allrate<bestRate
            bestRate=Allrate;
            bestMG=MG;
            bestClusterNumber=bestMGClusterNumber;
            for x=1: K
                 for y=1: N
                     %if(placementTemp(x,y)==0)
                        out(x,y)=placementTemp(x,y);
                     %end
                  end
            end
        end   
   % end
end
close(progress);
%out 
bestMG
bestClusterNumber
bestRate

out(K+1,1)= bestRate; 
out(K+1,2)= bestMG;
out(K+1,3)= bestClusterNumber;
%out

end
