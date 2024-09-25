%Best CQI scheduler
function [usagepertti,instant_throughput,usage,allocation]=best_cqi_scheduler(sectoring,usagepertti,tti,UEspercell,...
totaleNBs,totalRBs,instant_throughput,usage,modulation_order,coding_rate,allocation)

if sectoring==0
    metric=zeros(UEspercell,totaleNBs,'single');
    for i=1:totaleNBs
        for j=1:UEspercell
            if usage(j,i)==0
                metric(j,i)=throughputcalc(tti,5,modulation_order(j,i,1),modulation_order(j,i,2),coding_rate(j,i,1),coding_rate(j,i,2));
            else
                metric(j,i)=throughputcalc(tti,1,modulation_order(j,i,1),modulation_order(j,i,2),coding_rate(j,i,1),coding_rate(j,i,2));
            end
        end
        for m=1:totalRBs
            %[maxthr,maxind]=max(metric(:,i));
            maxthr=max(metric(:,i));
            maxindexes = find(metric(:,i)==maxthr);
            if length(maxindexes)>1
                maxind=maxindexes(randi(numel(maxindexes)));
            else
                maxind=maxindexes;
            end
            usagepertti(maxind,i)=usagepertti(maxind,i)+1;
            usage(maxind,i)=usage(maxind,i)+1;
            allocation{tti,i}(maxind,m)=1;
            %scheduled_once(maxind,i)=1;
            metric(maxind,i)=throughputcalc(tti,usagepertti(maxind,i),modulation_order(maxind,i,1),modulation_order(maxind,i,2),...
            coding_rate(maxind,i,1),coding_rate(maxind,i,2));
        end
        for j=1:UEspercell
            instant_throughput(j,i)=throughputcalc(tti,usagepertti(j,i),modulation_order(j,i,1),modulation_order(j,i,2),...
            coding_rate(j,i,1),coding_rate(j,i,2));
            %usage(j,i)=usage(j,i)+usagepertti{tti,1}(j,i);
        end
    end

elseif sectoring==1
    metric=zeros(max(UEspercell),totaleNBs,'single');
    for i=1:totaleNBs
        for j=UEspercell(1,i)+1:max(UEspercell)
            metric(j,i)=NaN;
        end
    end
    for i=1:totaleNBs
        for j=1:UEspercell(1,i)
            tf=isnan(instant_throughput(j,i));  %tf=1 an einai NaN, tf=0 an oxi
            if usage(j,i)==0 && tf==0
                metric(j,i)=throughputcalc(tti,5,modulation_order(j,i,1),modulation_order(j,i,2),coding_rate(j,i,1),coding_rate(j,i,2));
            elseif usage(j,i)~=0 && tf==0
                metric(j,i)=throughputcalc(tti,1,modulation_order(j,i,1),modulation_order(j,i,2),coding_rate(j,i,1),coding_rate(j,i,2));
            end
        end
        for m=1:totalRBs
            %[maxthr,maxind]=max(metric(:,i));
            maxthr=max(metric(:,i));
            maxindexes = find(metric(:,i)==maxthr);
            if length(maxindexes)>1
                maxind=maxindexes(randi(numel(maxindexes)));
            else
                maxind=maxindexes;
            end
            usagepertti(maxind,i)=usagepertti(maxind,i)+1;
            usage(maxind,i)=usage(maxind,i)+1;
            allocation{tti,i}(maxind,m)=1;
            %scheduled_once(maxind,i)=1;
            metric(maxind,i)=throughputcalc(tti,usagepertti(maxind,i),modulation_order(maxind,i,1),modulation_order(maxind,i,2),...
            coding_rate(maxind,i,1),coding_rate(maxind,i,2));
        end
        for j=1:UEspercell(1,i)
            tf=isnan(instant_throughput(j,i));  %tf=1 an einai NaN, tf=0 an oxi
            if tf==0
                instant_throughput(j,i)=throughputcalc(tti,usagepertti(j,i),modulation_order(j,i,1),modulation_order(j,i,2),...
                coding_rate(j,i,1),coding_rate(j,i,2));
                %usage(j,i)=usage(j,i)+usagepertti{tti,1}(j,i);
            end
        end
    end
end
    