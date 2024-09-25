%Round robin scheduler
function [instant_throughput,usagepertti,usage,allocation,rrcounter]=round_robin_scheduler(sectoring,tti,totaleNBs,UEspercell,instant_throughput,usagepertti,...
usage,allocation,rrcounter,modulation_order,coding_rate)

%% Arxikopoihseis
totalRBs=100;

%% round robin scheduling
if sectoring==0
    for i=1:totaleNBs
        for m=1:totalRBs
            usagepertti(rrcounter(1,i),i)=usagepertti(rrcounter(1,i),i)+1;
            allocation{tti,i}(rrcounter(1,i),m)=1;
            if rrcounter(1,i)<UEspercell
                rrcounter(1,i)=rrcounter(1,i)+1;
            else
                rrcounter(1,i)=1;
            end
        end
    end
    for i=1:totaleNBs
        for j=1:UEspercell
            instant_throughput(j,i)=throughputcalc(tti,usagepertti(j,i),modulation_order(j,i,1),modulation_order(j,i,2),coding_rate(j,i,1),coding_rate(j,i,2));
            usage(j,i)=usage(j,i)+usagepertti(j,i);
        end
    end
    
elseif sectoring==1
    for i=1:totaleNBs
        for m=1:totalRBs
            usagepertti(rrcounter(1,i),i)=usagepertti(rrcounter(1,i),i)+1;
            allocation{tti,i}(rrcounter(1,i),m)=1;
            if rrcounter(1,i)<UEspercell(1,i)
                rrcounter(1,i)=rrcounter(1,i)+1;
            else
                rrcounter(1,i)=1;
            end
        end
    end
    for i=1:totaleNBs %totalsectors
        for j=1:UEspercell(1,i)  %UEspersector
            tf=isnan(instant_throughput(j,i));  %tf=1 an einai NaN, tf=0 an oxi
            if tf==0
                instant_throughput(j,i)=throughputcalc(tti,usagepertti(j,i),modulation_order(j,i,1),modulation_order(j,i,2),coding_rate(j,i,1),coding_rate(j,i,2));
                usage(j,i)=usage(j,i)+usagepertti(j,i);
            end
        end
    end
end
