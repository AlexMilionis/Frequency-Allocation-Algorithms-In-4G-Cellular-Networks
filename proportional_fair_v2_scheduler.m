%Proportionally Fair v2 scheduler
function [usagepertti,instant_throughput,usage,allocation,metric,num,den,starvcounter]=proportional_fair_v2_scheduler(type,sectoring,usagepertti,tti,UEspercell,...
totaleNBs,totalRBs,instant_throughput,usage,modulation_order,coding_rate,allocation,metric,num,den,starvcounter)

%% arxikopoihseis
k=10;    %previous tti memory
a=1-1/k;
if type=="urban"
    death=14;
elseif type=="suburban"
    death=4;
end

%%
if sectoring==0
    for i=1:totaleNBs
        if tti==1
            for j=1:UEspercell
                den(j,i)=1;
            end
        end
        for m=1:totalRBs
            x=find(starvcounter(:,i)>=death,1); %briskw posa einai starved se auto to tti
            tf=isempty(x);  %tf=1 kanena starved, tf=0 exw starved
            if tf==1  %den exw starved UEs 
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
                starvcounter(maxind,i)=0;
                usagepertti(maxind,i)=usagepertti(maxind,i)+1;
                allocation{tti,i}(maxind,m)=1;
                den(maxind,i)=a*den(maxind,i)+(1-a)*num(maxind,i);
                for j=1:UEspercell
                     if j~=maxind   
                        den(j,i)=a*den(j,i);
                        %starvcounter(j,i)=starvcounter(j,i)+1;
                    end
                end
            elseif tf==0 %exw starved UEs
                [maxstarved,maxind]=max(starvcounter(:,i));
                usagepertti(maxind,i)=usagepertti(maxind,i)+1;
                allocation{tti,i}(maxind,m)=1;
                den(maxind,i)=a*den(maxind,i)+(1-a)*num(maxind,i);
                starvcounter(maxind,i)=0;
                for j=1:UEspercell
                    if j~=maxind
                        den(j,i)=a*den(j,i);
                        %starvcounter(j,i)=starvcounter(j,i)+1;
                    end
                end   
            end
        end
        for j=1:UEspercell
            if usagepertti(j,i)==0
                instant_throughput(j,i)=0;
                starvcounter(j,i)=starvcounter(j,i)+1;
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
                den(j,i)=1;
            end
        end
        for m=1:totalRBs
            x=find(starvcounter(:,i)>=death,1); %briskw posa einai starved se auto to tti
            tf=isempty(x);  %tf=1 kanena starved, tf=0 exw starved
            if tf==1  %den exw starved UEs 
                for j=1:UEspercell(1,i)
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
                starvcounter(maxind,i)=0;
                usagepertti(maxind,i)=usagepertti(maxind,i)+1;
                allocation{tti,i}(maxind,m)=1;
                den(maxind,i)=a*den(maxind,i)+(1-a)*num(maxind,i);
                for j=1:UEspercell(1,i)
                     if j~=maxind   
                        den(j,i)=a*den(j,i);
                        %starvcounter(j,i)=starvcounter(j,i)+1;
                    end
                end
            elseif tf==0 %exw starved UEs
                [maxstarved,maxind]=max(starvcounter(:,i));
                usagepertti(maxind,i)=usagepertti(maxind,i)+1;
                allocation{tti,i}(maxind,m)=1;
                den(maxind,i)=a*den(maxind,i)+(1-a)*num(maxind,i);
                starvcounter(maxind,i)=0;
                for j=1:UEspercell(1,i)
                    if j~=maxind
                        den(j,i)=a*den(j,i);
                        %starvcounter(j,i)=starvcounter(j,i)+1;
                    end
                end   
            end
        end
        for j=1:UEspercell(1,i)
            tf=isnan(instant_throughput(j,i));
            if tf==0
                if usagepertti(j,i)==0
                    instant_throughput(j,i)=0;
                    starvcounter(j,i)=starvcounter(j,i)+1;
                else
                    instant_throughput(j,i)=throughputcalc(tti,usagepertti(j,i),modulation_order(j,i,1),modulation_order(j,i,2),...
                    coding_rate(j,i,1),coding_rate(maxind,i,2));
                end
                usage(j,i)=usage(j,i)+usagepertti(j,i);
            end
        end
    end
    
end