%% epilogi propagation environment & algorithmou scheduling

answer = questdlg('Select propagation environment','Thesis','Urban','Suburban','Urban');
switch answer
    case 'Urban'
        type="urban";
        answer3=questdlg('Sectoring','Sectoring','Yes','No','Yes');
        switch answer3
            case 'Yes'
                sectoring=true(1);
            case 'No'
                sectoring=false(1);
        end
    
    case 'Suburban'
        type="suburban";
        sectoring=false(1);  
end
disp(type);

tic;

%% Arxikopoihsh parametrwn
Pt=30; %dBm
Pt_W=10^(Pt/10)/1000;
noise_var_dBm=-174; %dBm
noise_var=10^(noise_var_dBm/10); %linear
radius = path_loss_model(Pt,type,0,0); %(Pt,type,d,value)   bazo d=value=0; epeidi thelo na brw radius, type=1 urban,type=2 suburban 

%% diagramma isxyos-apostasis
f = 2140;   %MHz 
hr = 1.5;
hb = 12.5;
if type=="urban"
    x=linspace(1,radius,200);
    y=zeros(1,200);
    for i=1:200
        if x(1,i)<20
            y(1,i)=-35.4 + 26*log10(x(1,i)) + 20*log10(f) + 4;
        elseif x(1,i)>=20
            y(1,i)=-55.9 + 38*log10(x(1,i)) +(24.5+1.5*(f/925))*log10(f) + 8.5;
        end
    end
elseif type=="suburban"
    d0=100;
    lamda=physconst('LightSpeed')/(f*10^6);
    A=20*log10(4*pi*d0/lamda);
    Xf=6*log10(f/2000);
    Xh=-10.8*log10(hr/2000);
    a=4; b=0.0065; c=17.1;
    n=a-b*hb+c/hb;
    x=linspace(1,radius,200);
    y=A + 10*n*log10(x/d0) + Xf + Xh + 4;
end
figure;
plot(x,y);
xlabel('Distance (m)');
ylabel('Path Loss (dB)');

%% Frame
t = linspace(0,2*pi,7);
frame=1000;
x1=-frame/2; x2=frame/2;    %akra tou frame
y1=-frame/2; y2=frame/2;
x = [x1, x2, x2, x1, x1];
y = [y1, y1, y2, y2, y1];
figure;
plot(x,y,'r', 'LineWidth',1);
xlim([-frame/2-radius, frame/2+radius]);
ylim([-frame/2-radius, frame/2+radius]);
hold on;

%% Geografikes parametroi

if type=="urban"
    % pairno dimo athinon opou zoun peripou 17044 anthropoi/km^2. Theoro 1 UE
    % ana user, ara exw 17044 UEs
    %hexagon_surface=((3*sqrt(3))/2)*(radius^2);
    totaleNBs=17;
    %totaleNBs = ceil(frame^2 / hexagon_surface);    %number of eNBs fitting in the frame
    UEspercell=ceil(17044/totaleNBs);
    totalUEs=UEspercell*totaleNBs;
    %UEspercell = (totalUEs/totaleNBs); %number of UEs per cell
elseif type=="suburban"
    % pairno dimo naksou opou zoun peripou 576 anthropoi/km^2. Theoro 1 UE
    % ana user, ara exw 576 UEs
    %hexagon_surface=((3*sqrt(3))/2)*(radius^2);
    totaleNBs=11; %mporei 13 na kanw to sxima
    %totaleNBs = ceil(frame^2 / hexagon_surface);    %number of eNBs fitting in the frame
    UEspercell=ceil(576/totaleNBs); %53
    totalUEs=UEspercell*totaleNBs;  %583
end 
centers=zeros(totaleNBs,2);

%% Cell Map
%Center cell
x1 = 0 + radius*cos(t);
y1 = 0 + radius*sin(t);
c_x1 = 0; %BS Location Center Cell x-axis
c_y1 = 0; %BS Location Center Cell y-axis
plot(x1,y1,'blue');
plot(c_x1,c_y1,'r+');
hold on;
coord=subscribers_coordinates(c_x1,c_y1,radius,UEspercell);
centers(1,1)=c_x1;  centers(1,2)=c_y1;
coordinates=zeros(totalUEs,2);
counter=0;
for j=1:UEspercell
    coordinates((counter*UEspercell)+j,1)=coord(j,1);
    coordinates((counter*UEspercell)+j,2)=coord(j,2);
end
counter=counter+1;

% Loop
sumoutside=0;
n=1;
dx=1.5*radius;
dy=(sqrt(3)/2)*radius;
i=2;
while sumoutside ~= 6*n
    sumoutside=0;
    cx=c_x1 + 0;
    cy=c_y1 + n*2*dy;
    
    if (cx>=-frame/2) && (cx<=frame/2) && (cy>=-frame/2) && (cy<=frame/2)
        x = cx + radius*cos(t);
        y = cy + radius*sin(t);
        plot(x,y,'blue');
        plot(cx,cy,'r+');
        hold on;
        coord=subscribers_coordinates(cx,cy,radius,UEspercell);
        centers(i,1)=cx;
        centers(i,2)=cy;
        i=i+1;
        for j=1:UEspercell
            coordinates((counter*UEspercell)+j,1)=coord(j,1);
            coordinates((counter*UEspercell)+j,2)=coord(j,2);
        end
        counter=counter+1;
    else
        sumoutside=sumoutside+1;
    end
    
    for k=1:n
        cx = cx + dx;
        cy = cy - dy;
        if (cx>=-frame/2) && (cx<=frame/2) && (cy>=-frame/2) && (cy<=frame/2)
            x = cx + radius*cos(t);
            y = cy + radius*sin(t);
            plot(x,y,'blue');
            plot(cx,cy,'r+');
            hold on;
            coord=subscribers_coordinates(cx,cy,radius,UEspercell);
            centers(i,1)=cx;
            centers(i,2)=cy;
            i=i+1;
            for j=1:UEspercell
                coordinates((counter*UEspercell)+j,1)=coord(j,1);
                coordinates((counter*UEspercell)+j,2)=coord(j,2);
            end
            counter=counter+1;
        else
            sumoutside=sumoutside+1;
        end
    end
    
    for k=1:n
        cy=cy-2*dy;
        if (cx>=-frame/2) && (cx<=frame/2) && (cy>=-frame/2) && (cy<=frame/2)
            x = cx + radius*cos(t);
            y = cy + radius*sin(t);
            plot(x,y,'blue');
            plot(cx,cy,'r+');
            hold on;
            coord=subscribers_coordinates(cx,cy,radius,UEspercell);
            centers(i,1)=cx;
            centers(i,2)=cy;
            i=i+1;
            for j=1:UEspercell
                coordinates((counter*UEspercell)+j,1)=coord(j,1);
                coordinates((counter*UEspercell)+j,2)=coord(j,2);
            end
            counter=counter+1;
        else
            sumoutside=sumoutside+1;
        end
    end
    
    for k=1:n
        cx=cx-dx;
        cy=cy-dy;
        if (cx>=-frame/2) && (cx<=frame/2) && (cy>=-frame/2) && (cy<=frame/2)
            x = cx + radius*cos(t);
            y = cy + radius*sin(t);
            plot(x,y,'blue');
            plot(cx,cy,'r+');
            hold on;
            coord=subscribers_coordinates(cx,cy,radius,UEspercell);
            centers(i,1)=cx;
            centers(i,2)=cy;
            i=i+1;
            for j=1:UEspercell
                coordinates((counter*UEspercell)+j,1)=coord(j,1);
                coordinates((counter*UEspercell)+j,2)=coord(j,2);
            end
            counter=counter+1;
        else
            sumoutside=sumoutside+1;
        end
    end
    
    for k=1:n
        cx=cx-dx;
        cy=cy+dy;
        if (cx>=-frame/2) && (cx<=frame/2) && (cy>=-frame/2) && (cy<=frame/2)
            x = cx + radius*cos(t);
            y = cy + radius*sin(t);
            plot(x,y,'blue');
            plot(cx,cy,'r+');
            hold on;
            coord=subscribers_coordinates(cx,cy,radius,UEspercell);
            centers(i,1)=cx;
            centers(i,2)=cy;
            i=i+1;
            for j=1:UEspercell
                coordinates((counter*UEspercell)+j,1)=coord(j,1);
                coordinates((counter*UEspercell)+j,2)=coord(j,2);
            end
            counter=counter+1; 
        else
            sumoutside=sumoutside+1;
        end
    end
    
    for k=1:n
        cy = cy + 2*dy;
        if (cx>=-frame/2) && (cx<=frame/2) && (cy>=-frame/2) && (cy<=frame/2)
            x = cx + radius*cos(t);
            y = cy + radius*sin(t);
            plot(x,y,'blue');
            plot(cx,cy,'r+');
            hold on;
            coord=subscribers_coordinates(cx,cy,radius,UEspercell);
            centers(i,1)=cx;
            centers(i,2)=cy;
            i=i+1;
            for j=1:UEspercell
                coordinates((counter*UEspercell)+j,1)=coord(j,1);
                coordinates((counter*UEspercell)+j,2)=coord(j,2);
            end
            counter=counter+1;
        else
            sumoutside=sumoutside+1;
        end
    end
    
    for k=1:n
        cx = cx + dx;
        cy = cy + dy;
        if (cx>=-frame/2) && (cx<=frame/2) && (cy>=-frame/2) && (cy<=frame/2)
            x = cx + radius*cos(t);
            y = cy + radius*sin(t);
            plot(x,y,'blue');
            plot(cx,cy,'r+');
            hold on;
            if k~=n
                coord=subscribers_coordinates(cx,cy,radius,UEspercell);
                centers(i,1)=cx;
                centers(i,2)=cy;
                i=i+1;
                for j=1:UEspercell
                    coordinates((counter*UEspercell)+j,1)=coord(j,1);
                    coordinates((counter*UEspercell)+j,2)=coord(j,2);
                end
                counter=counter+1;
            end
        else
            sumoutside=sumoutside+1;
        end
    end
    sumoutside=sumoutside-1;
    if sumoutside==6*n
        break;
    end
    n=n+1;
    sumoutside=0;
end

%% sectoring
if type=="urban" && sectoring==1
    coordinatespercell=zeros(UEspercell,totaleNBs,2);
    totalsectors=3*totaleNBs;
    UEspersector=zeros(1,totalsectors);

    j=1; i=1;
    for q=1:totalUEs  %bazw syntetagmenenes x,y se pinaka (j,i)
        coordinatespercell(j,i,1)=coordinates(q,1);
        coordinatespercell(j,i,2)=coordinates(q,2);
        j=j+1;
        if mod(q,UEspercell)==0
            i=i+1;
            j=1;
        end
    end

    s=0;
    for i=1:totaleNBs  %sximatizw tis grammes twn sectors
        cx=centers(i,1);    cy=centers(i,2);
        x1=cx-radius;   y1=cy; %A
        x2=cx-radius/2; y2=cy+sqrt(3)*radius/2; %B
        x3=cx+radius/2; y3=cy+sqrt(3)*radius/2; %C
        x4=cx+radius; y4=cy; %D
        x5=cx+radius/2; y5=cy-sqrt(3)*radius/2; %E
        x6=cx-radius/2; y6=cy-sqrt(3)*radius/2; %F

        t=x1:cx;
        y = cy + (cy-y1)/(cx-x1) * (t-cx);
        plot(t,y,'--k');
        hold on

        t=cx:x3;
        y = cy + (cy-y3)/(cx-x3) * (t-cx);
        plot(t,y,'--k');
        hold on

        t=cx:x3;
        y = cy + (cy-y5)/(cx-x5) * (t-cx);
        plot(t,y,'--k');
        hold on

        for j=1:UEspercell %briskw posa UE exei o kathe sector
            if inpolygon(coordinatespercell(j,i,1),coordinatespercell(j,i,2),[x1,x2,x3,cx,x1],[y1,y2,y3,cy,y1])==1
                UEspersector(1,s+1)=UEspersector(1,s+1)+1;
            elseif inpolygon(coordinatespercell(j,i,1),coordinatespercell(j,i,2),[cx,x3,x4,x5,cx],[cy,y3,y4,y5,cy])==1
                UEspersector(1,s+2)=UEspersector(1,s+2)+1;
            elseif inpolygon(coordinatespercell(j,i,1),coordinatespercell(j,i,2),[x1,cx,x5,x6,x1],[y1,cy,y5,y6,y1])==1
                UEspersector(1,s+3)=UEspersector(1,s+3)+1;
            end
        end
        s=s+totalsectors/totaleNBs;
    end
end

%% Long Term SINR
[Long_Term_SINR,positions,r0,r,Losses]=long_term_SINR(Pt,type,centers,coordinates); %epistrefei linear sinr

%% clear variables
clearvars c_x1 c_y1 counter cx cy dx dy i j k n sumoutside t x x1 x2 y y1 y2

%% Channel Settings
% Settings
Nt=4; %transmit antennas  2x4 MIMO
Nr=2; %receive antennas
e=0.85; %time correlation
t=0.5;  %spatial correlation
tti=1000;  %tti=1ms
%----------------------------------------------------------------------------------------------------------------------
Hbar=cell(totalUEs,totaleNBs);      %iid circularly symmetric complex gaussian random variables that are uncorrelated
N=cell(totalUEs,totaleNBs);         %noise matrix
R=cell(totalUEs,totaleNBs);         %transmit correlation matrix
H=cell(totalUEs,totaleNBs);         %channel matrix
Rn=cell(totalUEs,2);                %covariance matrix

% Htest=cell(tti,1);
% vartest=zeros(tti,1);
% Ntest=cell(tti,1);

P=precoding(Pt);   %16 precoding matrices, size 4x2
norms=zeros(16,1,'single');
Popt=cell(1,totalUEs);  %matrix pou krata to optimal precoder gia kathe link enos UE me to eNB tou

tempsinr1=zeros(1,totalUEs,'single');
tempsinr2=zeros(1,totalUEs,'single');
sinr1=zeros(UEspercell,totaleNBs,'single');        %SINR1 gia kathe link
sinr2=zeros(UEspercell,totaleNBs,'single');        %SINR2 gia kathe link
sinr1_dB=zeros(UEspercell,totaleNBs,'single');
sinr2_dB=zeros(UEspercell,totaleNBs,'single');
cqisum=zeros(tti,15,2);
avg_cqi=zeros(UEspercell,totaleNBs);

%----------------------------------------------------------------------------------------------------------------------
%xwris sectoring
totalRBs=100;
instant_throughput=cell(tti,4);
usagepertti=cell(tti,4);
usage=cell(1,4);
for alg=1:4
    usage{1,alg}=zeros(UEspercell,totaleNBs); 
    for k=1:tti
        instant_throughput{k,alg}=zeros(UEspercell,totaleNBs);
        usagepertti{k,alg}=zeros(UEspercell,totaleNBs);  
    end
end

allocation=cell(tti,totaleNBs);
for k=1:tti
    for i=1:totaleNBs
        allocation{k,i}=zeros(UEspercell,100);
    end
end
scheduled_once=zeros(UEspercell,totaleNBs,'single');    %1 an exei ginei schedule mia fora, allios 0
%-----------------------------------------------------------------------------------------------------
if sectoring==1  
    instant_throughputsc=cell(tti,4);
    usageperttisc=cell(tti,4);
    usagesc=cell(1,alg);
    for alg=1:4
        usagesc{1,alg}=zeros(max(UEspersector),totalsectors);
        for k=1:tti
            instant_throughputsc{k,alg}=zeros(max(UEspersector),totalsectors);
            usageperttisc{k,alg}=zeros(max(UEspersector),totalsectors);  
        end
    end
    for alg=1:4
        for i=1:totalsectors
            for j=(UEspersector(1,i)+1):max(UEspersector)
                usagesc{1,alg}(j,i)=NaN;
            end
        end
        for k=1:tti
            for i=1:totalsectors
                for j=(UEspersector(1,i)+1):max(UEspersector)
                    instant_throughputsc{k,alg}(j,i)=NaN;
                    usageperttisc{k,alg}(j,i)=NaN;
                end
            end
        end
    end
%     instant_throughputsc=zeros(max(UEspersector),totalsectors,tti,'single');
%     usagesc=zeros(max(UEspersector),totalsectors,'single');
%     usageperttisc=cell(tti,1);
    allocationsc=cell(tti,totalsectors);
    for k=1:tti
        usageperttisc{k,1}=zeros(max(UEspersector),totalsectors);
        for i=1:totalsectors
            allocationsc{k,i}=zeros(max(UEspersector),100);
        end
    end
    %bazw NaN sta stoixeia pou perisevoun
    for i=1:totalsectors
        for j=(UEspersector(1,i)+1):max(UEspersector)
%             usagesc(j,i)=NaN;
%             instant_throughputsc(j,i,:)=NaN;
            for k=1:tti
%                 usageperttisc{k,1}(j,i)=NaN;
                for m=1:totalRBs
                    allocationsc{k,i}(j,m)=NaN;
                end
            end
        end
    end    
end
%-----------------------------------------------------------------------------------------------------
if sectoring==0
    rrcounter=ones(1,totaleNBs);
else
    rrcounter=ones(1,totalsectors);
end
if sectoring==0
    metricPF1=zeros(UEspercell,totaleNBs,tti,'single');
    numPF1=zeros(UEspercell,totaleNBs);
    denPF1=zeros(UEspercell,totaleNBs);
    metricPF2=zeros(UEspercell,totaleNBs,tti,'single');
    numPF2=zeros(UEspercell,totaleNBs);
    denPF2=zeros(UEspercell,totaleNBs);
    starvcounter=zeros(UEspercell,totaleNBs);
elseif sectoring==1
    metricscPF1=zeros(max(UEspersector),totalsectors,tti,'single');
    numscPF1=zeros(max(UEspersector),totalsectors);
    denscPF1=zeros(max(UEspersector),totalsectors);
    metricscPF2=zeros(max(UEspersector),totalsectors,tti,'single');
    numscPF2=zeros(max(UEspersector),totalsectors);
    denscPF2=zeros(max(UEspersector),totalsectors);
    starvcountersc=zeros(max(UEspersector),totalsectors);
 
    for i=1:totalsectors
        for j=(UEspersector(1,i)+1):max(UEspersector)
            metricscPF1(j,i,:)=NaN;
            numscPF1(j,i)=NaN;
            denscPF1(j,i)=NaN;
            metricscPF2(j,i,:)=NaN;
            numscPF2(j,i)=NaN;
            denscPF2(j,i)=NaN;
            starvcountersc(j,i)=NaN;
        end
    end
end


%% R (correlation matrix), den allazei me to xrono (4x4)
for q=1:totalUEs
    for i=1:totaleNBs
        if centers(i,1)~=positions(q,3) || centers(i,2)~=positions(q,4)
            tqi = 0;    %den iparxei correlation me tous interferer eNBs
        else
            fi=(2*pi-0).*rand(1)+0;
            tqi = t*exp(1i*fi); % gia ton eNB poy anikei to UE
        end
        R{q,i}=[1 tqi tqi^2 tqi^3 ; conj(tqi) 1 tqi tqi^2 ;(conj(tqi))^2 conj(tqi) 1 tqi ; (conj(tqi))^3 (conj(tqi))^2 conj(tqi) 1];
    end
end

%% TTI Loop
for k=1:tti
    for q=1:totalUEs
        for i=1:totaleNBs
            %% Noise matrix (2x4)
            %N{q,i} = sqrt(1/2)*(randn(Nr,Nt)+1i*randn(Nr,Nt));
            N{q,i} = sqrt(1/2)*(normrnd(0,0.5,[Nr,Nt])+1i*normrnd(0,1,[Nr,Nt]));
            
            if centers(i,1)==positions(q,3) && centers(i,2)==positions(q,4)
                var = 10^(-path_loss_model(Pt,type,r0(q,1),1)/10);    %PL tou q-UE me ton eNB tou cell tou
            elseif centers(i,1)~=positions(q,3) || centers(i,2)~=positions(q,4)
                var = 10^(-path_loss_model(Pt,type,r(q,i),1)/10);  %PL tou q-UE me eNB diaforetikou cell
            end
            Hbar{q,i} = sqrt(1/2)*(normrnd(0,0.5,[Nr,Nt])+1i*normrnd(0,1,[Nr,Nt]));

            
            %% H (channel matrix) (2x4)
            H{q,i} = sqrt(var)*Hbar{q,i}*(sqrtm(R{q,i}));
              
            %% Popt (optimal precoding matrix) (4x2)
            %tora tha brw gia kathe UE poio P(i} tairiazei kalitera apo ta 16 me ton eNB tou cell tou,
            %diladi poio dinei max ||H*P||
            if r(q,i)==0
                maxindex=0;
                maximum=0;
                for n=1:16  %16 sindiasmoi precoding gia kathe link enos UE me ton enb tou cell tou
                    HP=H{q,i}*P{n,1};  %to H.P einai 2x2 matrix!
                    %norms(n,1)=norm(HP,'fro')^2;
                    norms(n,1)=norm(HP);
                    if norms(n,1)>maximum
                        maximum=norms(n,1);
                        maxindex=n;
                    end
                end
                Popt{1,q}=P{maxindex,1};
            end
%             Hbar{1,q,i}=Hbar{2,q,i};
%             Hbar{2,q,i}=zeros(Nr,Nt);
        end
        %Popt{1,q}=P{maxindex,1};
    end
    
    %% Rn (covariance matrix) (2x2)
    for q=1:totalUEs
        suma1=zeros(2,2,'single');    %sum gia to proto stream tou protou orou tou Rn
        suma2=zeros(2,2,'single');    %sum gia to deutero stream tou protou orou tou Rn
        sumb=0;              %sum deuterou orou tou 1Rn
        for i=1:totaleNBs
            if r(q,i)==0
                suma1=(H{q,i}*Popt{1,q}(:,1)) * (H{q,i}*Popt{1,q}(:,1))';
                suma2=(H{q,i}*Popt{1,q}(:,2)) * (H{q,i}*Popt{1,q}(:,2))';
            end
            if r(q,i)~=0   %q user den anikei sto i eNB
                n=randi(16);
                sumb = sumb + (H{q,i}*P{n,1})*(H{q,i}*P{n,1})';
            end
        end
        Rn{q,1} = suma1 + sumb + noise_var*eye(2);    %covariance matrix tou 1st stream
        Rn{q,2} = suma2 + sumb + noise_var*eye(2);    %covariance matrix tou 2nd stream
    end
    
    %% Instantaneous SINR 
    for q=1:totalUEs
        for i=1:totaleNBs
            if r(q,i)==0
                tempsinr1(1,q)=abs((H{q,i}*Popt{1,q}(:,1))' * (Rn{q,2}^(-1)) * (H{q,i}*Popt{1,q}(:,1)));
                tempsinr2(1,q)=abs((H{q,i}*Popt{1,q}(:,2))' * (Rn{q,1}^(-1)) * (H{q,i}*Popt{1,q}(:,2)));

            end
        end
    end
    j=1; i=1;
    for q=1:totalUEs
        sinr1(j,i)=tempsinr1(1,q);
        sinr2(j,i)=tempsinr2(1,q);
        sinr1_dB(j,i)=10*log10(sinr1(j,i));
        sinr2_dB(j,i)=10*log10(sinr2(j,i));
        j=j+1;
        if mod(q,UEspercell)==0
            i=i+1;
            j=1;
        end
    end
    
    
    %% CQI
    [cqi,modulation_order,coding_rate,cqisum]=cqi_mapping(sinr1_dB,sinr2_dB,UEspercell,totaleNBs,k,cqisum);
    for i=1:totaleNBs
        for j=1:UEspercell
            avg_cqi(j,i)=avg_cqi(j,i)+(cqi(j,i,1)+cqi(j,i,2))/2;
        end
    end
    
    if sectoring==1
        cqisc=zeros(max(UEspersector),totalsectors,2,'single');
        modulation_ordersc=zeros(max(UEspersector),totalsectors,2,'single');
        coding_ratesc=zeros(max(UEspersector),totalsectors,2,'single');
        for i=1:totalsectors
            for j=(UEspersector(1,i)+1):max(UEspersector)
                cqisc(j,i,:)=NaN;
                modulation_order(j,i,:)=NaN;
                coding_rate(j,i,:)=NaN;
            end
        end
        s=0;
        for i=1:totaleNBs
            for a=1:2
                x=UEspersector(1,s+1);
                y=UEspersector(1,s+2);
                z=UEspersector(1,s+3);
                cqisc(1:x,s+1,a)=cqi(1:x,i,a);
                cqisc(1:y,s+2,a)=cqi(x+1:x+y,i,a);
                cqisc(1:z,s+3,a)=cqi(x+y+1:x+y+z,i,a);
                modulation_ordersc(1:x,s+1,a)=modulation_order(1:x,i,a);
                modulation_ordersc(1:y,s+2,a)=modulation_order(x+1:x+y,i,a);
                modulation_ordersc(1:z,s+3,a)=modulation_order(x+y+1:x+y+z,i,a);
                coding_ratesc(1:x,s+1,a)=coding_rate(1:x,i,a);
                coding_ratesc(1:y,s+2,a)=coding_rate(x+1:x+y,i,a);
                coding_ratesc(1:z,s+3,a)=coding_rate(x+y+1:x+y+z,i,a);
            end
            s=s+3;
        end
    end
    
    %% Scheduling
    X=['tti=',num2str(k)];
    disp(X);
    if sectoring==0
        temp1=instant_throughput{k,1}; temp2=usagepertti{k,1}; temp3=usage{1,1};
        [temp1,temp2,temp3,allocation,rrcounter]=round_robin_scheduler(sectoring,k,totaleNBs,UEspercell,temp1,temp2,temp3,allocation,rrcounter,modulation_order,coding_rate);
        instant_throughput{k,1}=temp1; usagepertti{k,1}=temp2; usage{1,1}=temp3;
        
        temp1=instant_throughput{k,2}; temp2=usagepertti{k,2}; temp3=usage{1,2};
        [temp2,temp1,temp3,allocation,metricPF1,numPF1,denPF1]=proportional_fair_scheduler(sectoring,temp2,k,UEspercell,totaleNBs,totalRBs,temp1,temp3,modulation_order,coding_rate,allocation,metricPF1,numPF1,denPF1);
        instant_throughput{k,2}=temp1; usagepertti{k,2}=temp2; usage{1,2}=temp3;
        
        temp1=instant_throughput{k,3}; temp2=usagepertti{k,3}; temp3=usage{1,3};
        [temp2,temp1,temp3,allocation,metricPF2,numPF2,denPF2,starvcounter]=proportional_fair_v2_scheduler(type,sectoring,temp2,k,UEspercell,totaleNBs,totalRBs,temp1,temp3,modulation_order,coding_rate,allocation,metricPF2,numPF2,denPF2,starvcounter);
        instant_throughput{k,3}=temp1; usagepertti{k,3}=temp2; usage{1,3}=temp3;
        
        temp1=instant_throughput{k,4}; temp2=usagepertti{k,4}; temp3=usage{1,4};
        [temp2,temp1,temp3,allocation] =best_cqi_scheduler(sectoring,temp2,k,UEspercell,totaleNBs,totalRBs,temp1,temp3,modulation_order,coding_rate,allocation);
        instant_throughput{k,4}=temp1; usagepertti{k,4}=temp2; usage{1,4}=temp3;
        
    elseif sectoring==1
        temp1=instant_throughputsc{k,1}; temp2=usageperttisc{k,1}; temp3=usagesc{1,1};
        [temp1,temp2,temp3,allocationsc,rrcounter]=round_robin_scheduler(sectoring,k,totalsectors,UEspersector,temp1,temp2,temp3,allocationsc,rrcounter,modulation_ordersc,coding_ratesc);
        instant_throughputsc{k,1}=temp1; usageperttisc{k,1}=temp2; usagesc{1,1}=temp3;
        
        temp1=instant_throughputsc{k,2}; temp2=usageperttisc{k,2}; temp3=usagesc{1,2};
        [temp2,temp1,temp3,allocationsc,metricscPF1,numscPF1,denscPF1]=proportional_fair_scheduler(sectoring,temp2,k,UEspersector,totalsectors,totalRBs,temp1,temp3,modulation_ordersc,coding_ratesc,allocationsc,metricscPF1,numscPF1,denscPF1);
        instant_throughputsc{k,2}=temp1;  usageperttisc{k,2}=temp2; usagesc{1,2}=temp3;
        
        temp1=instant_throughputsc{k,3}; temp2=usageperttisc{k,3}; temp3=usagesc{1,3};
        [temp2,temp1,temp3,allocationsc,metricscPF2,numscPF2,denscPF2,starvcountersc]=proportional_fair_v2_scheduler(type,sectoring,temp2,k,UEspersector,totalsectors,totalRBs,temp1,temp3,modulation_ordersc,coding_ratesc,allocationsc,metricscPF2,numscPF2,denscPF2,starvcountersc);
        instant_throughputsc{k,3}=temp1; usageperttisc{k,3}=temp2; usagesc{1,3}=temp3;
        
        temp1=instant_throughputsc{k,4}; temp2=usageperttisc{k,4}; temp3=usagesc{1,4};
        [temp2,temp1,temp3,allocationsc] =bestcqischeduler(sectoring,temp2,k,UEspersector,totalsectors,totalRBs,temp1,temp3,modulation_ordersc,coding_ratesc,allocationsc);
        instant_throughputsc{k,4}=temp1; usageperttisc{k,4}=temp2; usagesc{1,4}=temp3;
    end
  
end

%% post-scheduling results
if sectoring==1
    s=0;
    for alg=1:4
        for i=1:totaleNBs
            x=UEspersector(1,s+1);
            y=UEspersector(1,s+2);
            z=UEspersector(1,s+3);
            usage{1,alg}(1:x,i)=usagesc{1,alg}(1:x,s+1);
            usage{1,alg}(x+1:x+y,i)=usagesc{1,alg}(1:y,s+2);
            usage{1,alg}(x+y+1:x+y+z,i)=usagesc{1,alg}(1:z,s+3);
            s=s+3;
        end
        s=0;
    end
    
    s=0;
    for alg=1:4
        for k=1:tti
            for i=1:totaleNBs
                x=UEspersector(1,s+1);
                y=UEspersector(1,s+2);
                z=UEspersector(1,s+3);
                instant_throughput{k,alg}(1:x,i)=instant_throughputsc{k,alg}(1:x,s+1);
                instant_throughput{k,alg}(x+1:x+y,i)=instant_throughputsc{k,alg}(1:y,s+2);
                instant_throughput{k,alg}(x+y+1:x+y+z,i)=instant_throughputsc{k,alg}(1:UEspersector(1,s+3),s+3);
                usagepertti{k,alg}(1:x,i)=usageperttisc{k,alg}(1:x,s+1);
                usagepertti{k,alg}(x+1:x+y,i)=usageperttisc{k,alg}(1:y,s+2);
                usagepertti{k,alg}(x+y+1:x+y+z,i)=usageperttisc{k,alg}(1:z,s+3);
                for m=1:totalRBs
                    allocation{k,i}(1:x,m)=allocationsc{k,s+1}(1:x,m);
                    allocation{k,i}(x+1:x+y,m)=allocationsc{k,s+2}(1:y,m);
                    allocation{k,i}(x+y+1:x+y+z,m)=allocationsc{k,s+3}(1:z,m); 
                end
                s=s+3;
            end
            s=0;
        end
        s=0;
    end 
end

%------------------------------------------------------------------------------------------------------------------------
avg_throughput=cell(1,4);
for alg=1:4
    avg_throughput{1,alg}=zeros(UEspercell,totaleNBs);
end
avgthr_eNB=cell(1,4);
for alg=1:4
    avgthr_eNB{1,alg}=zeros(1,totaleNBs);
end
total_cellthr=cell(1,4);
for alg=1:4
    total_cellthr{1,alg}=zeros(1,totaleNBs);
end
avgthr_final=zeros(1,4);
cellthr_final=zeros(1,4);
for alg=1:4
    for i=1:totaleNBs
        for j=1:UEspercell
            for k=1:tti
                avg_throughput{1,alg}(j,i)=avg_throughput{1,alg}(j,i)+instant_throughput{k,alg}(j,i);
            end
            avg_throughput{1,alg}(j,i)=avg_throughput{1,alg}(j,i)/tti;        % bits ana tti
            avg_throughput{1,alg}(j,i)=avg_throughput{1,alg}(j,i)*1000;       % bits ana sec
            avg_throughput{1,alg}(j,i)=avg_throughput{1,alg}(j,i)*(10^(-6));  % Mbits ana sec
        end
        avgthr_eNB{1,alg}(1,i)=sum(avg_throughput{1,alg}(:,i))/UEspercell;
        total_cellthr{1,alg}(1,i)=sum(avg_throughput{1,alg}(:,i));
    end
    avgthr_final(1,alg)=sum(avgthr_eNB{1,alg}(1,:))/totaleNBs;
    cellthr_final(1,alg)=sum(total_cellthr{1,alg}(1,:))/totaleNBs;
end

%% graphs

%cdf RR,PF1,PF2,BCQI
figure(3);
h1=cdfplot(avg_throughput{1,1}(:));
h1.Color='blue'; 
h1.LineWidth=2;
hold on;
h2=cdfplot(avg_throughput{1,2}(:));
h2.Color='magenta';
h2.LineWidth=2;
%h2.LineStyle='--';
hold on;
h3=cdfplot(avg_throughput{1,3}(:));
h3.Color='green';
h3.LineWidth=2;
%h3.LineStyle='--';
hold on;
h4=cdfplot(avg_throughput{1,4}(:));
h4.Color='red';
h4.LineWidth=2;
hold on;
legend('Round Robin','Proportional Fair','Proportional Fair v2','Best CQI');
hold off;

%cdf RR,PF1,PF2
figure(4);
h1=cdfplot(avg_throughput{1,1}(:));
h1.Color='blue'; 
h1.LineWidth=2;
hold on;
h2=cdfplot(avg_throughput{1,2}(:));
h2.Color='magenta';
h2.LineWidth=2;
%h2.LineStyle='--';
hold on;
h3=cdfplot(avg_throughput{1,3}(:));
h3.Color='green';
h3.LineWidth=2;
%h3.LineStyle='--';
hold on;
legend('Round Robin','Proportional Fair','Proportional Fair v2','Best CQI');
hold off;

%throughput barchart
figure(5);
bar(avgthr_final);
x={'Round Robin';'Proportional Fair';'Proportional Fair v2';'Best CQI'};
set(gca,'xticklabel',x);
xlabel('Schedulers');
ylabel('Average throughput (Mbit/s)');

%cell throughput
figure(6);
bar(cellthr_final);
x={'Round Robin';'Proportional Fair';'Proportional Fair v2';'Best CQI'};
set(gca,'xticklabel',x);
xlabel('Schedulers');
ylabel('Cell throughput (Mbit/s)');

%fairness
fairness=zeros(1,4);
for alg=1:4
    sum1=0; sum2=0;
    for i=1:totaleNBs
        for j=1:UEspercell
            sum1=sum1+avg_throughput{1,alg}(j,i);
            sum2=sum2+(avg_throughput{1,alg}(j,i))^2;
        end
    end
    fairness(1,alg)=(sum1^2)/(UEspercell*totaleNBs*sum2);
end
figure(7);
bar(fairness);
x={'Round Robin';'Proportional Fair';'Proportional Fair v2';'Best CQI'};
set(gca,'xticklabel',x);
xlabel('Schedulers');
ylabel('Fairness');

for i=1:totaleNBs
    for j=1:UEspercell
        avg_cqi(j,i)=avg_cqi(j,i)/tti;
    end
end


