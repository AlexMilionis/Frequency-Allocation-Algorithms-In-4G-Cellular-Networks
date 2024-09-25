%Proportionally Fair scheduler
function [usagepertti,instant_throughput,usage,allocation,metric,num,den]=proportional_fair_scheduler(sectoring,usagepertti,tti,UEspercell,totaleNBs,...
totalRBs,instant_throughput,usage,modulation_order,coding_rate,allocation,metric,num,den)

%% arxikopoihseis
k=10;    %previous tti memory
a=1-1/k;

if sectoring==0
    for i=1:totaleNBs
        if tti==1
            for j=1:UEspercell
                %num(j,i)=throughputcalc(tti,5,modulation_order(j,i,1),modulation_order(j,i,2),coding_rate(j,i,1),coding_rate(j,i,2));
                den(j,i)=1;
                %metric(j,i,tti)=num(j,i)/den(j,i);
            end
        end
        for m=1:totalRBs
            for j=1:UEspercell
                num(j,i)=throughputcalc(tti,1,modulation_order(j,i,1),modulation_order(j,i,2),coding_rate(j,i,1),coding_rate(j,i,2));
                metric(j,i,tti)=num(j,i)/den(j,i);
            end
            maxthr=max(metric(:,i,tti));
            maxindexes = find(metric(:,i,tti)==maxthr);
            if length(maxindexes)>1
                maxind=maxindexes(randi(numel(maxindexes)));
            else
                maxind=maxindexes;
            end
            usagepertti(maxind,i)=usagepertti(maxind,i)+1;
            allocation{tti,i}(maxind,m)=1;
            %scheduled_once(maxind,i)=1;
            den(maxind,i)=a*den(maxind,i)+(1-a)*num(maxind,i);
            %metric(maxind,i,tti)=num(maxind,i)/den(maxind,i);
            for j=1:UEspercell
                %if j~=maxind && usagepertti{tti,1}(j,i)==0
                 if j~=maxind   
                    %num(j,i)=throughputcalc(tti,1,modulation_order(j,i,1),modulation_order(j,i,2),coding_rate(j,i,1),coding_rate(j,i,2));
                    den(j,i)=a*den(j,i);
                    %metric(j,i,tti)=num(j,i)/den(j,i);
                end
            end
        end
        for j=1:UEspercell
            if usagepertti(j,i)==0
                instant_throughput(j,i)=0;
            else
                instant_throughput(j,i)=throughputcalc(tti,usagepertti(j,i),modulation_order(j,i,1),modulation_order(j,i,2),...
                coding_rate(j,i,1),coding_rate(maxind,i,2));
            end
            usage(j,i)=usage(j,i)+usagepertti(j,i);
        end
    end
    
elseif sectoring==1
    for i=1:totaleNBs
        if tti==1 
            for j=1:UEspercell(1,i)
                tf=isnan(instant_throughput(j,i));  %tf=1 an einai NaN, tf=0 an oxi
                if tf==0
                    den(j,i)=1;
                end
            end
        end
        for m=1:totalRBs
            for j=1:UEspercell(1,i)
                tf=isnan(metric(j,i));  %tf=1 an einai NaN, tf=0 an oxi
                if tf==0
                    num(j,i)=throughputcalc(tti,1,modulation_order(j,i,1),modulation_order(j,i,2),coding_rate(j,i,1),coding_rate(j,i,2));
                    metric(j,i,tti)=num(j,i)/den(j,i);
                end
            end
            maxthr=max(metric(:,i,tti));
            maxindexes = find(metric(:,i,tti)==maxthr);
            if length(maxindexes)>1
                maxind=maxindexes(randi(numel(maxindexes)));
            else
                maxind=maxindexes;
            end
            usagepertti(maxind,i)=usagepertti(maxind,i)+1;
            allocation{tti,i}(maxind,m)=1;
            den(maxind,i)=a*den(maxind,i)+(1-a)*num(maxind,i);
            for j=1:UEspercell(1,i)
                tf=isnan(metric(j,i));  %tf=1 an einai NaN, tf=0 an oxi
                if j~=maxind && tf==0
                    den(j,i)=a*den(j,i);
                end
            end
        end
        for j=1:UEspercell(1,i)
            tf=isnan(instant_throughput(j,i));
            if tf==0
                if usagepertti(j,i)==0
                    instant_throughput(j,i)=0;
                else
                    instant_throughput(j,i)=throughputcalc(tti,usagepertti(j,i),modulation_order(j,i,1),modulation_order(j,i,2),...
                    coding_rate(j,i,1),coding_rate(maxind,i,2));
                end
                usage(j,i)=usage(j,i)+usagepertti(j,i);
            end
        end
    end
end

end
    