function [LTSINR,positions,r0,r,Losses]=long_term_SINR(Pt,type,centers,coordinates) 
    %pairnei eisodo ta kentra twn cells kai tis syntetagmenes olwn twn simeiwn

    %% Calculate distance matrix r0
    totaleNBs=length(centers);      
    totalUEs=length(coordinates);   
    UEspercell=length(coordinates)/length(centers);
    Pnoise=-174; %dBm
    interferingchannels=totaleNBs;
    
    positions=zeros(totalUEs,4);
    
    for i=1:totalUEs
        positions(i,1)=coordinates(i,1);
        positions(i,2)=coordinates(i,2);
    end
    
    for j=0:totaleNBs-1
        for i=1:UEspercell
            positions(j*UEspercell+i,3)=centers(j+1,1);
            positions(j*UEspercell+i,4)=centers(j+1,2);
        end
    end
    
    r0=zeros(totalUEs,1);
    for i=1:totalUEs
        r0(i,1)=sqrt((positions(i,1)-positions(i,3))^2+(positions(i,2)-positions(i,4))^2);
        
    end
    
    %% Calculate distance matrix r
    
    r=zeros(totalUEs,interferingchannels);
    
 
    %An exw interference apo ola ta alla eNBs
    for i=1:totalUEs
        for j=1:interferingchannels
            if centers(j,1)~=positions(i,3) || centers(j,2)~=positions(i,4)
                r(i,j)=sqrt((positions(i,1)-centers(j,1))^2+(positions(i,2)-centers(j,2))^2);
            end
        end
    end
        
   %% Computing Losses of UE
   
   Losses0 = zeros(totalUEs,1);
   for i=1:totalUEs
       Losses0(i,1)=PathLoss(Pt,type,r0(i,1),1);    
       Losses0(i,1) = 10^(-Losses0(i,1)/10);    %linear times
       %Losses0(i,2)=shadowing;
   end

   Losses = zeros(totalUEs,interferingchannels); 
   
   
   %An exw interference apo ola ta alla eNBs
   for i=1:totalUEs
       for j=1:interferingchannels
           if r(i,j)~=0
               Losses(i,j)=PathLoss(Pt,type,r(i,j),1);
               Losses(i,j) = 10^(-Losses(i,j)/10);
           end
       end
   end
    
   
   Pnoise_mW=10^(Pnoise/10);
   Pt_mW=10^(Pt/10);
   
   SINR=zeros(totalUEs,1);
   for i=1:totalUEs
       numerator = (Losses0(i,1))*Pt_mW;        %((10^((Losses0(i,1)+Losses0(i,2))/10))^(-1)*1000)*Pt_mW;
       sum=0;
       for j=1:interferingchannels
           if Losses(i,j)~=0
                sum=sum+Losses(i,j);
           end
       end
       denominator = Pnoise_mW + Pt_mW*sum;
       SINR(i,1) =  numerator / denominator; %((L0^(-1))*Pt_mW)/(Pnoise_mW + (L^(-1))*Pt_mW );
   end
    
      
    LTSINR=zeros(UEspercell,totaleNBs,'single');
    j=1; i=1;
    for q=1:totalUEs
        %sinr1(j,i)=10*log10(tempsinr1(1,q));
        %sinr2(j,i)=10*log10(tempsinr2(1,q));
        LTSINR(j,i)=SINR(q,1);
        j=j+1;
        if mod(q,UEspercell)==0
            i=i+1;
            j=1;
        end
    end
end
