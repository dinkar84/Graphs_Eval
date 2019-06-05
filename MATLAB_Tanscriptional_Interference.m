clc; clear all; % v=50bp/s; s=70; b=40; b+2f=40+2*(20)=80;

Ta = 0;
Tb = 40;
Tc = 0;
s = 70;


TCLMean = 0;
TCRMean = 0;

TCLTotal = 0;
TCRTotal = 0;

Rep = 1;

TCLAddRep = zeros(1,Rep);
TCRAddRep = zeros(1,Rep);


k = input('Enter the rate constants for the left and right promoters: ');

r = k(1) + k(2);


T = 1000;


for i1=1:1:Rep
    
    i1
    
    TLO = []; % array containing time of transcription initiation from the left promoter    
    TRO = []; % array containing time of transcription initiation from the right promoter
    
    TL = []; % array containing time of transcription initiation of transcribing RNAPs from the left promoter   
    TR = []; % array containing time of transcription initiation of transcribing RNAPs from the right promoter
    
    TCL = []; % array containing time of transcription completion of RNAPs originating from the left promoter 
    TCR = []; % array containing time of transcription completion of RNAPs originating from the right promoter
                                                                                                              
    TNCL_time = []; % array containing time of collision of RNAPs initiated from the left promoter
    TNCR_time = []; % array containing time of collision of RNAPs initiated from the right promoter
    
    TNCL_len = []; % array containing the length of the truncated transcript of the gene on the left
    TNCR_len = []; % array containing the length of the truncated transcript of the gene on the right
    
    TNCL_init_time = []; % array containing time of transcription initiation from the left promoter of those RNAPs that underwent collision
    TNCR_init_time = []; % array containing time of transcription initiation from the right promoter of those RNAPs that underwent collision
        
      
     
    Time = zeros (1,1);
 
    m=0; mO=0; n=0; nO=0; u=1; v=1; wL=1; wR=1; t=1;
    
    tao = 0; 
    Ts1 = 0;
    Ts2 = 0;
               
     
% disp('While loop starts');

 while 1
 %    disp ('Time is: ');
 %    Time(1,t)
    if (Time(1,t) > T)
     % disp('Program ends');
        break;
                
    else 
          B = rand();
          C = log(B);
          tao = (-1)*C/r;
                    
          Time (1,t+1) = Time (1,t) + tao;
          t=t+1;
          
          D = rand();
          
          if (D*r <= k(1))  
              
              if (tao<Ts1)
                  Ts1 = Ts1 - tao;
                  Ts2 = Ts2 - tao;
                  continue;
                  
              else
              Ts1 = 1/50;  
              Ts2 = Ts2 - tao;              
              m=length(TL)+1;
              TL (1,m) = Time(1,t); 
              mO=mO+1;              
              TLO (1,mO) = Time (1,t);              
              end
                        
%              disp ('TL = ');
%              disp('Left promoter');
%              TL
             
          else  
              
               if (tao<Ts2)                   
                   Ts1 = Ts1 - tao;
                   Ts2 = Ts2 - tao;
                   continue;
                  
               else   
               Ts1 = Ts1 - tao;
               Ts2 = 1/50;               
               n=length(TR)+1;
               TR (1,n) = Time (1,t);
               nO=nO+1;              
               TRO (1,nO) = Time (1,t);
               end
    end
%              disp ('TR = ');
%              disp('Right promoter');
%              TR
%              
          end
    
    
               if (isempty(TL) || isempty(TR)) 
                    
                   if (isempty(TL))
                       if (length(TR) == 1)                               
                          %  disp('TR has only one element, so continued');
                          continue;
                                                                    
                       else
                                 c=0; % disp('TL empty TR more than one RNAP');  
                             for j= 1:1:length(TR)-1
                                 if ((TR(length(TR)) - TR(j)) > Tc + Tb)                                  
                                 TCR(1,v) = TR(j) + Tc + Tb;
                             %   TCR(1,v) = TR(length(TR)); 
                                 v=v+1;
                                 c=c+1;
                                 end
                             
                             end
                            % removing RNAPs which have completed transcription  
                             for i= 1:1:c
                             TR(1) = []; 
                             end
                             n=n-c;                            
                             continue;
                             
                       end
                            
                   else          %(isempty(TR))
                                     
                                     
                    if (length(TL) == 1)                       
                       %  disp('TL has only one element, so continued');
                       continue;
                       
                    else
                           c=0; % disp('TR empty TL more than one RNAP');  
                       for j= 1:1:length(TL)-1
                            if  ((TL(length(TL)) - TL(j)) > Ta + Tb)  
                                TCL(1,u) = TL(j) + Ta + Tb;
                            %   TCL(1,u) = TL(length(TL));
                                u=u+1;
                                c=c+1; 
                            end
                            
                       end
                       % removing RNAPs which have completed transcription  
                           for i= 1:1:c
                           TL(1) = []; 
                           end 
                           m=m-c;                          
                           %%
                           continue;
                            
                    end
                   end   
                   
                   
               else             
%                            disp ('TL1');
%                            disp(TL);
%                            disp(length(TL));
%                            
%                            disp ('TR1');
%                            disp(TR);                           
%                            disp(length(TR));
                                                    
                             
                             iL=1;                               
                             while (1)   
                                                                  
                      
                                    iR=1;                                    
                                    while (1)                                         
                                         z=0;
                                                                                                                                               
                                       if (    (   (  ( ((TL(iL) + Ta) - (TR(iR) + Tc))  >= 0 ) && ( ((TL(iL) + Ta) - (TR(iR) + Tc))  < Tb )  ) || (  ( ((TR(iR) + Tc) - (TL(iL) + Ta))  >= 0 ) &&  ( ((TR(iR) + Tc) - (TL(iL) + Ta)) < Tb )  )   ) &&  ( (Time(1,t) - TL(iL)) > Ta ) && ( (Time(1,t) - TR(iR)) > Tc ) && (length(TL)~=0) && (length(TR)~=0)    )                               
                                           
                                           z=1;
                                           TNCL_len(wL)=(TL(iL)+TR(iR)+Ta+Tb+Tc)/2 - TL(iL); % calculates the length of the truncated transcript from the left promoter
                                           TNCL_time(wL)=(TL(iL)+TR(iR)+Ta+Tb+Tc)/2;
                                         % TNCL_time(wL)=Time(1,t);
                                         
                                           TNCR_len(wR)=(TL(iL)+TR(iR)+Ta+Tb+Tc)/2 - TR(iR); % calculates the length of the truncated transcript from the right promoter
                                           TNCR_time(wR)=(TL(iL)+TR(iR)+Ta+Tb+Tc)/2;
                                         % TNCR_time(wR)=Time(1,t);
                                                                                               
                                           wL=wL+1;
                                           TL(iL)=[];
                                           iL=iL-1;
                                                                        
                                           wR=wR+1;
                                           TR(iR)=[];
                                           iR=iR-1;
                                                                                    
                                           
                                           if (length(TL)==0 || length(TR)==0)                                               
                                               break;                                                                                   
                                           end           
                                                       
                                       end              
                                                                                                                                                                   
                                       if (iL==length(TL) || iR==length(TR))                                                    
                                           break;                                           
                                       end                                                                                 
                                             
                                          if (z==1)
                                              iL=iL+1;                                              
                                          end
                                          
                                          iR=iR+1;
                                                                                                                                                                                 
                                    end           
                                                   
                                                              
                                 if (length(TL)==0 || length(TR)==0)                                    
                                     break;                                                                                   
                                 end         
                                             
                                 if (iL==length(TL))                                 
                                     break;         
                                 end                       
                                                    
                                    iL=iL+1;
                                                                                                                         
                             end          
                             
                             
%                            disp ('TL2');
%                            disp(TL);
%                                                     
%                            disp ('TR2');
%                            disp(TR);
                                                                             
                                            
                             
                            
                           if (length(TL)~=0)
                           
                             iL=1;                          
                             while (1)                                
                                                                     
                                % if (  ( (TL(iL) + Ta + Tb) < (TR(iR) + Tc) ) &&  ( (Time(1,t) - TL(iL)) > (Ta + Tb) )  ) 
                                  if ( (Time(1,t) - TL(iL)) > (Ta + Tb) )  
                                      TCL(1,u) = TL(iL) + Ta + Tb;    
                                    % TCL(1,u) = Time(1,t);
                                      u=u+1;
                                      TL(iL)=[];
                                      iL=iL-1;
                                     
                                    
                                    if (length(TL)==0)
%                                         disp('(TL)==0');
                                        break;                                                                            
                                    end
                                    
                                  end
                                                        
                                                                                                                       
                                  if (iL==length(TL))
                                      break;
                                  end
                                  
                                  iL=iL+1;
                                  
                             end 
                             
                           end 
                           
%                            disp('TL ended');
%                              
%                            
%                            disp ('TL3');
%                            disp(TL);
%                            disp(length(TL));
%                            
%                            disp ('TR3');
%                            disp(TR);
%                            disp(length(TR));
                           
                            
                           if (length(TR) ~= 0)
                             % RNAPs completing transcription from the right promoter 
                             % using while loop instead of for loop due to non-constant size of the array TR     
                                                              
                               iR=1;                          
                               while (1)                                 
                                        
                                                                                                        
                                % if (  ( (TL(iL) + Ta) > (TR(iR) + Tc + Tb) ) &&  ( (Time(1,t) - TR(iR)) > (Tc + Tb) )  )  
                                  if ( (Time(1,t) - TR(iR)) > (Tc + Tb) )         
                                      TCR(1,v) = TR(iR) + Tc + Tb;    
                                    % TCR(1,v) = Time(1,t);
                                      v=v+1;
                                      TR(iR)=[];
                                      iR=iR-1;
                                      
                                      
                                    if (length(TR)==0)
%                                         disp('(TR)==0');
                                        break;
                                                                               
                                    end
                                    
                                  end
                                  
                                  
                                 
                                  % last RNAP cannot complete transcription as it has just bound the promoter!
                            
                                  if (iR==length(TR))
                                      break;
                                  end
                                  
                                  iR=iR+1;
                                  
                               end 
                               
                           end
                           
%                            disp('TR ended');
                             
               end
              
                             
              
                             

%disp ('Time of RNAP completing the transcription from the left promoter');

%disp(TCL)

%disp ('Time of RNAP completing the transcription from the right promoter');
%disp(TCR);



%disp ('Number of RNAPs from the left promoter: ');
%disp (length(TCL));
%disp ('Number of RNAPs from the right promoter: ');
%disp (length(TCR));

 end

TCLAddRep(1,i1) = length(TCL);
TCRAddRep(1,i1) = length(TCR);

end

TCLTotal = 0;
TCRTotal = 0;

for i = 1:1:Rep
    
TCLTotal = TCLTotal + TCLAddRep(1,i1);
TCRTotal = TCRTotal + TCRAddRep(1,i1);

end   


TCLMean = TCLTotal/Rep;
TCRMean = TCRTotal/Rep;

TCLMean
TCRMean


TCLVarT = 0;
TCRVarT = 0;

for i = 1:1:Rep
    
TCLVarT = TCLVarT + power((TCLAddRep(1,i)-TCLMean),2);
TCRVarT = TCRVarT + power((TCRAddRep(1,i)-TCRMean),2);

end   
      
TCLVar = TCLVarT/Rep;
TCRVar = TCRVarT/Rep;

disp ('Left SD: ');
sqrt(TCLVar)

disp ('Right SD: ');
sqrt(TCRVar) 

% disp('length(TL)');
% disp(length(TL));
% 
% disp('length(TR)');
% disp(length(TR));

        
      
% TLO
% TL
% TCL
% TNCL_time
% TNCL_len
% TNCL_init_time = TNCL_time - TNCL_len 
% 
% 
% TRO
% TR
% TCR
% TNCR_time
% TNCR_len
% TNCR_init_time = TNCR_time - TNCR_len