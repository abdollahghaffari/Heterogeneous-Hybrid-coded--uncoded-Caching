% Abdollah Ghaffari sheshjavani 4/5/2020 , 19/01/1399
function out=subclusterHuristicFunction(Z,K,N,M,popularityArray,Method,popularityType)

%Method=1; % 1=our hybrid   2=purecoded  3=pureUncoded
%popularityType; =0= generalPopularity        =1= Huristic1
% Cntr+R  Cntrl+T
out=zeros(K+2,N);
cachceckstart=0;
cachchecksize=0;
Zmax=max(Z);
if Method==1
    cachchecksize=M;
end
if Method==3
    cachceckstart=M;
    cachchecksize=M;
end
% this is variables
BestBandwidth=-1;
bestN1=-1; % it is the best N*
bestM1=-1;
base=0;
bestbase=-1;
best_r2=-1;
MBSBackhaulOverhead=0;
%///////////////////////////////////////////////////////// Pre computation Section //////////////////////////////////////
% myfactorial(1)=1;
% for l=1 : K
%     myfactorial(l+1)=l*myfactorial(l); % notice ==> factorial(k)=myfactorial(k+1)
% end
% /////////this section compute combination(Z,j) or factorial(n)/(factorial(k)*factorial(n-k))with dynomic programing
%com_max=max(K,Zmax);
computed_comb=zeros(K,K);
for i=1 : K
    computed_comb(i,1)=i;
   for j=2 : K
       if i<j
            computed_comb(i,j)=0;
       else
           computed_comb(i,j)=computed_comb(i,j-1)*(i-j+1)/j;
       end
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
%sorted_sum_popularityArray
%///////////////////////////////////////////////////////// END of Pre computation Section //////////////////////////////////////
%popularityArray
for N1=M : N 
    zarib1=0;
    for M1=cachceckstart : cachchecksize
        Flag=false;
        if M1==M 
            if N1==M
                Flag=true;
            end
            if N1>M
                break;
            end
        end
        if mod(K*(M-M1),(N1-M1))==0 || Flag==true  
        
            X=zeros(1,N);
            Y=zeros(K,N);      
  
     %//////////////////placement (huristic and ...) methods to determined X and Y ...
            
            if popularityType==1 % Huristic1
                for M1_itrator=1 : M1
                    for k_itrator=1:K
                        Y(k_itrator,most_popularityindex(k_itrator,M1_itrator))= 1;
                    end
                end

                N1_counter=0;
                if N1>M1
                    for content_itrator=1 : N
                        coded_content=true;
                        for k_itrator=1:K
                            if Y(k_itrator,sorted_sum_popularityArray(1,content_itrator))==1
                                coded_content=false;
                            end
                        end
                        if coded_content == true
                            X(1,sorted_sum_popularityArray(1,content_itrator))=1;
                            N1_counter=N1_counter+1;
                        end
                        if N1_counter>= (N1-M1)
                            break;
                        end
                    end
                    if N1_counter<(N1-M1)
                        errrrrrrrorrrrrrrrrrrrrrrrrrrr=N1    
                        break;
                    end
                end
            elseif popularityType==0
                for M1_itrator=1 : M1
                    for k_itrator=1:K
                        Y(k_itrator,sorted_sum_popularityArray(1,M1_itrator))= 1;
                    end
                end
                if N1>M1
                    for content_itrator=M1+1 : N1
                       X(1,sorted_sum_popularityArray(1,content_itrator))=1; 
                    end 
                end     
            end 
   %///////////////END placement (huristic and general.) methods to determined X and Y ...          
  % ///////////////////////////computing r2//////////////////////////////////////////////
           r2=0;
            for n=1 : N
                temp_r2=1;
                for c=1 : K 
                   temp_r2=temp_r2*(1-(popularityArray(c,n)*(1-Y(c,n))))^Z(c);
                end
                temp_r2=1-temp_r2;
                r2=r2+ (1-X(1,n))*temp_r2;%(1-(1-(((1/n)^ziph_parameter)/allziph))^sigmaZ); %(Z*K));
                %order1=order1+1;
            end
 % ///////////////////////////END computing r2////////////////////////////////////////////// 
          if Flag==true
              MBSBackhaulOverhead=r2;
              T=0;
          else
               T=K*(M-M1)/(N1-M1);
               base=Zmax*(K*(N1-M)/((N1-M1)+K*(M-M1))); % K-T/T+1
               MBSBackhaulOverhead=base+r2;
          end
           expectation=0;
           if T>0
               %zarib1=1/(combination(K,T));
               %zarib1=1/ (myfactorial(K+1)/(myfactorial(T+1)*myfactorial(K-T+1)));
               zarib1=1/ computed_comb(K,T);
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                PZ= zeros(K,Z(c)+1, Z(c)+1);
                for c=1: K
                   SumofP_Codedfiles=0;
                   for n=1 : N
                       SumofP_Codedfiles=SumofP_Codedfiles+ X(1,n)*popularityArray(c,n);
                   end  
                   qcn(1)= SumofP_Codedfiles;
                   
                   EP_codedfiles=SumofP_Codedfiles/(N1-M1);                 
                   for j=2: Z(c)+1
                           qcn(j)= SumofP_Codedfiles-((j-1)*EP_codedfiles); 
                           if (j>N1-M1) 
                                qcn(j)=0;  
                           end
                        %end
                   end
%                    
                  
%                    for j=2: Z(c)+1
%                        qcn(j)= qcn(j-1)*((N1-M1-j+1)/(N1-M1)); 
%                        if (j>N1-M1) 
%                             qcn(j)=0;  
%                        end   
%                    end
%                    
                   
%                    qcn(1)= SumofP_Codedfiles;
%                    codedindex=1;
%                    for j=2: Z(c)+1
%                        if (j>N1-M1) 
%                             qcn(j)=0;  
%                        else
%                            while X(1,most_popularityindex(c,codedindex))==0
%                                codedindex=codedindex+1;
%                            end
%                            qcn(j)=qcn(j-1)-popularityArray(c,most_popularityindex(c,codedindex));
%                            codedindex=codedindex+1;
%                        end
%                    end
%                    qcn
     %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\            
                   expectation=0;
                   PZ(c,1,1)= 1.0;
                   for s=2: Z(c)+1
                        for qipz=1: s
                            if qipz == 1 
                                PZ(c,s,qipz)=PZ(c,s-1,qipz)*(1 - qcn(1));
                            else 
                                PZ(c,s,qipz) = (PZ(c,s-1,qipz) * (1 -qcn(qipz))+ PZ(c,s-1,qipz-1) * qcn(qipz-1)); 
                            end
                        end
                   end
               end
% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
               for i=1: Zmax
                   for ki=0: K
                       PQ = zeros(K, K+1);
                       PQ(1,1)=1.0;
                       for c=1: K  
                            Pzi(c)=0;
                            if i<=Z(c)
                               for j=i :Z(c)   
                                    Pzi(c)=Pzi(c)+PZ(c,Z(c)+1,j+1);
                                end
                            end
                %  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                            for qi=0: c
                                if qi == 0 
                                    PQ(c+1,qi+1)=PQ(c,qi+1)*(1 - Pzi(c));
                                else 
                                    PQ(c+1,qi+1) = (PQ(c,qi+1) * (1 - Pzi(c))+ PQ(c,qi) * Pzi(c)); 
                                end
                            end
                        end
                        if K-ki >=T+1  
                             %expectation=expectation+PQ(K+1,ki+1)*(myfactorial(K-ki+1)/(myfactorial(T+1+1)*myfactorial(K-ki-T)));
                             expectation=expectation+PQ(K+1,ki+1)*computed_comb(K-ki,T+1);
                        end
                   end
               end
           end
            
           %N1
          % M1
          % r2
          % expectation
         %  zarib1
          % MBSBackhaulOverhead
           MBSBackhaulOverhead=MBSBackhaulOverhead-(zarib1*expectation);
           if MBSBackhaulOverhead < 0
                errorr=1
           end
          % MBSBackhaulOverhead
          if BestBandwidth==-1
              BestBandwidth=MBSBackhaulOverhead;
              bestN1=N1;
              bestM1=M1;
              for x=1: K
                 for y=1: N
                     out(x,y)=0;
                     if(X(1,y)==1)
                         out(x,y)=2;
                     end
                     if(Y(x,y)==1)
                         out(x,y)=1;
                     end
                 end
              end
               best_r2=r2;
               bestbase=base;
              %bestN1
              %bestM1
              %BestBandwidth
          else if MBSBackhaulOverhead<=BestBandwidth
               BestBandwidth=MBSBackhaulOverhead;
               bestN1=N1;
               bestM1=M1;
               for x=1: K
                 for y=1: N
                     out(x,y)=0;
                     if(X(1,y)==1)
                         out(x,y)=2;
                     end
                     if(Y(x,y)==1)
                         out(x,y)=1;
                     end
                 end
              end
               best_r2=r2;
               bestbase=base;
               %bestN1
               %bestM1
               %BestBandwidth
              end
          end
        end
end
end
% Show Results
%Zmax
%ziph_parameter
%bestN1
%bestM1
%best_r2
%BestBandwidth
%AnalyticalRate=BestBandwidth
out(K+1,1)=BestBandwidth;
out(K+2,1)=BestBandwidth-best_r2;
%out
%best_r2
%bestbase
%BestBandwidth
%x = '=======================================';
%x
end
%order1
%order2
%order3