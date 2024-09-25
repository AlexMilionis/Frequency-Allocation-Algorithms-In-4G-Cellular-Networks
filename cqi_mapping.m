function [cqi,modulation_order,coding_rate,cqisum]=cqi_mapping(sinr1,sinr2,UEspercell,totaleNBs,tti,cqisum)
            
    %% arxikopoihseis
    cqi=zeros(UEspercell,totaleNBs,2,'single');
    modulation_order=zeros(UEspercell,totaleNBs,2,'single');
    coding_rate=zeros(UEspercell,totaleNBs,2,'single');
    
    
    %% sinr1
    for j=1:UEspercell
        for i=1:totaleNBs
            if sinr1(j,i)<-5.147
                cqi(j,i,1)=1;
                modulation_order(j,i,1)=2; 
                coding_rate(j,i,1)=78;
                %sumcqi(1,1)=sumcqi(1,1)+1;
            elseif sinr1(j,i)>=-5.147 && sinr1(j,i)<-3.18 
                cqi(j,i,1)=2;     
                modulation_order(j,i,1)=2;
                coding_rate(j,i,1)=120;
                %sumcqi(2,2)=sumcqi(2,2)+1;
            elseif sinr1(j,i)>=-3.18 && sinr1(j,i)<-1.254  
                cqi(j,i,1)=3;
                modulation_order(j,i,1)=2;
                coding_rate(j,i,1)=193;
                %sumcqi(3,2)=sumcqi(3,2)+1;
            elseif sinr1(j,i)>=-1.254 && sinr1(j,i)<0.761    
                cqi(j,i,1)=4;  
                modulation_order(j,i,1)=2;
                coding_rate(j,i,1)=308;
                %sumcqi(4,2)=sumcqi(4,2)+1;
            elseif sinr1(j,i)>=0.761 && sinr1(j,i)<2.70   
                cqi(j,i,1)=5;
                modulation_order(j,i,1)=2;
                coding_rate(j,i,1)=449;
                %sumcqi(5,2)=sumcqi(5,2)+1;
            elseif sinr1(j,i)>=2.70 && sinr1(j,i)<4.697   
                cqi(j,i,1)=6;
                modulation_order(j,i,1)=2;
                coding_rate(j,i,1)=602;
                %sumcqi(6,2)=sumcqi(6,2)+1;
            elseif sinr1(j,i)>=4.697 && sinr1(j,i)<6.528  
                cqi(j,i,1)=7;
                modulation_order(j,i,1)=4;
                coding_rate(j,i,1)=378;
                %sumcqi(7,2)=sumcqi(7,2)+1;
            elseif sinr1(j,i)>=6.528 && sinr1(j,i)<8.576    
                cqi(j,i,1)=8;
                modulation_order(j,i,1)=4;
                coding_rate(j,i,1)=490;
                %sumcqi(8,2)=sumcqi(8,2)+1;
            elseif sinr1(j,i)>=8.576 && sinr1(j,i)<10.37 
                cqi(j,i,1)=9;
                modulation_order(j,i,1)=4;
                coding_rate(j,i,1)=616;
                %sumcqi(9,2)=sumcqi(9,2)+1;
            elseif sinr1(j,i)>=10.37 && sinr1(j,i)<12.3 
                cqi(j,i,1)=10;
                modulation_order(j,i,1)=6;
                coding_rate(j,i,1)=466;
                %sumcqi(10,2)=sumcqi(10,2)+1;
            elseif sinr1(j,i)>=12.3 && sinr1(j,i)<14.18 
                cqi(j,i,1)=11;
                modulation_order(j,i,1)=6;
                coding_rate(j,i,1)=567;
                %sumcqi(11,2)=sumcqi(11,2)+1;
            elseif sinr1(j,i)>=14.18 && sinr1(j,i)<15.89 
                cqi(j,i,1)=12;
                modulation_order(j,i,1)=6;
                coding_rate(j,i,1)=666;
                %sumcqi(12,2)=sumcqi(12,2)+1;
            elseif sinr1(j,i)>=15.89 && sinr1(j,i)<17.82 
                cqi(j,i,1)=13;
                modulation_order(j,i,1)=6;
                coding_rate(j,i,1)=772;
                %sumcqi(13,2)=sumcqi(13,2)+1;
            elseif sinr1(j,i)>=17.82 && sinr1(j,i)<19.83  
                cqi(j,i,1)=14;
                modulation_order(j,i,1)=6;
                coding_rate(j,i,1)=873;
                %sumcqi(14,2)=sumcqi(14,2)+1;
            elseif sinr1(j,i)>=19.83
                cqi(j,i,1)=15;
                modulation_order(j,i,1)=6;
                coding_rate(j,i,1)=948;
                %sumcqi(15,2)=sumcqi(16,2)+1;
            end
            cqisum(tti,cqi(j,i,1),1)=cqisum(tti,cqi(j,i,1),1)+1;
            coding_rate(j,i,1)=coding_rate(j,i,1)/1024;
            modulation_order(j,i,1)=modulation_order(j,i,1);
        end
    end
    
    %% sinr2
    for j=1:UEspercell
        for i=1:totaleNBs
            if sinr2(j,i)<-5.147
                cqi(j,i,2)=1;
                modulation_order(j,i,2)=2;
                coding_rate(j,i,2)=78;
                %sumcqi(1,1)=sumcqi(1,1)+1;
            elseif sinr2(j,i)>=-5.147 && sinr2(j,i)<-3.18 
                cqi(j,i,2)=2;
                modulation_order(j,i,2)=2;
                coding_rate(j,i,2)=120;
                %sumcqi(2,2)=sumcqi(2,2)+1;
            elseif sinr2(j,i)>=-3.18 && sinr2(j,i)<-1.254  
                cqi(j,i,2)=3; 
                modulation_order(j,i,2)=2;
                coding_rate(j,i,2)=193;
                %sumcqi(3,2)=sumcqi(3,2)+1;
            elseif sinr2(j,i)>=-1.254 && sinr2(j,i)<0.761    
                cqi(j,i,2)=4;
                modulation_order(j,i,2)=2;
                coding_rate(j,i,2)=308;
                %sumcqi(4,2)=sumcqi(4,2)+1;
            elseif sinr2(j,i)>=0.761 && sinr2(j,i)<2.70   
                cqi(j,i,2)=5;
                modulation_order(j,i,2)=2;
                coding_rate(j,i,2)=449;
                %sumcqi(5,2)=sumcqi(5,2)+1;
            elseif sinr2(j,i)>=2.70 && sinr2(j,i)<4.697   
                cqi(j,i,2)=6;
                modulation_order(j,i,2)=2;
                coding_rate(j,i,2)=602;
                %sumcqi(6,2)=sumcqi(6,2)+1;
            elseif sinr2(j,i)>=4.697 && sinr2(j,i)<6.528  
                cqi(j,i,2)=7;
                modulation_order(j,i,2)=4;
                coding_rate(j,i,2)=378;
                %sumcqi(7,2)=sumcqi(7,2)+1;
            elseif sinr2(j,i)>=6.528 && sinr2(j,i)<8.576    
                cqi(j,i,2)=8;
                modulation_order(j,i,2)=4;
                coding_rate(j,i,2)=490;
                %sumcqi(8,2)=sumcqi(8,2)+1;
            elseif sinr2(j,i)>=8.576 && sinr2(j,i)<10.37 
                cqi(j,i,2)=9;
                modulation_order(j,i,2)=4;
                coding_rate(j,i,2)=616;
                %sumcqi(9,2)=sumcqi(9,2)+1;
            elseif sinr2(j,i)>=10.37 && sinr2(j,i)<12.3 
                cqi(j,i,2)=10;
                modulation_order(j,i,2)=6;
                coding_rate(j,i,2)=466;
                %sumcqi(10,2)=sumcqi(10,2)+1;
            elseif sinr2(j,i)>=12.3 && sinr2(j,i)<14.18 
                cqi(j,i,2)=11;
                modulation_order(j,i,2)=6;
                coding_rate(j,i,2)=567;
                %sumcqi(11,2)=sumcqi(11,2)+1;
            elseif sinr2(j,i)>=14.18 && sinr2(j,i)<15.89 
                cqi(j,i,2)=12;
                modulation_order(j,i,2)=6;
                coding_rate(j,i,2)=666;
                %sumcqi(12,2)=sumcqi(12,2)+1;
            elseif sinr2(j,i)>=15.89 && sinr2(j,i)<17.82 
                cqi(j,i,2)=13;
                modulation_order(j,i,2)=6;
                coding_rate(j,i,2)=772;
                %sumcqi(13,2)=sumcqi(13,2)+1;
            elseif sinr2(j,i)>=17.82 && sinr2(j,i)<19.83  
                cqi(j,i,2)=14;
                modulation_order(j,i,2)=6;
                coding_rate(j,i,2)=873;
                %sumcqi(14,2)=sumcqi(14,2)+1;
            elseif sinr2(j,i)>=19.83
                cqi(j,i,2)=15;
                modulation_order(j,i,2)=6;
                coding_rate(j,i,2)=948;
                %sumcqi(15,2)=sumcqi(16,2)+1;
            end
            cqisum(tti,cqi(j,i,2),2)=cqisum(tti,cqi(j,i,2),2)+1;
            coding_rate(j,i,2)=coding_rate(j,i,2)/1024;
            modulation_order(j,i,2)=modulation_order(j,i,2);
        end
    end
    
        
end