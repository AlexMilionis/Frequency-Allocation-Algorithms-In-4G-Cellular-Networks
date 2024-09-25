function x=path_loss_model(Pt,type,d,value)
    %d einai h apostasi UE-eNB kai value= 0 or 1 = cell radius or losses
    
    %% System specs
    %Pt = 30; %dBm
    f = 2140;   %MHz 
    hb = 12.5; %m (10-80) transmitter antenna
    hr = 1.5; %m receiver antenna
    Gt = 10; %dBi transmission antenna gain
    Gr = 0; %dBi receiver antenna gain
    Lcoax = (0.1/hb)*30; %dB apwleies kalodiou
    Pthr = -90; %dBm
    
    %% Urban environment: COST 231 Walfisch Ikegami model
    
    if type=="urban"
        if value==0 && d==0  %Calculating cell radius 

            % PL = -55.9+38logd +(24.5+1.5*(f/925))*logf + S;
            PL = Pt + Gt + Gr - Lcoax - Pthr;  %ypologismos pathloss L
            S = 8.5; %dB
            x=10^((PL + 55.9 - (24.5 + 1.5*f/925)*log10(f)-S)/38);
            return;

        elseif value==1    % Calculating Losses
            d0=20;   %m
            if d>=d0
                S = 8.5; %dB
                PL_dB = -55.9 + 38*log10(d) +(24.5+1.5*(f/925))*log10(f) + S; %gia d>20m theoro NLOS
            else
                S=4;
                PL_dB = -35.4 + 26*log10(d) + 20*log10(f) + S;  %gia d<20m theoro proseggistika LOS
            end
            shadowing_dB = S*randn(1); %dB
            Losses_dB = PL_dB + shadowing_dB - Gt;
            %Losses = 10^(-Losses_dB/10);
            x = Losses_dB;            
        end
        
    end
    %% Suburban environment: %%SUI
       
    if type=="suburban"
        d0=100;
        lamda=physconst('LightSpeed')/(f*10^6);
        A=20*log10(4*pi*d0/lamda);
        Xf=6*log10(f/2000);
        Xh=-10.8*log10(hr/2000);
        a=4; b=0.0065; c=17.1;
        n=a-b*hb+c/hb;
        S=4;
        
        if value==0 && d==0  %Calculating cell radius
            
            % PL = A + 10*n*log10(d/d0) + Xf + Xh + S;
            PL = Pt + Gt + Gr - Lcoax - Pthr;  %ypologismos pathloss PL
            exponential=(PL-A-Xf-Xh-S)/(10*n);
            x=d0*(10^exponential);
            return;
    
        elseif value==1    % Calculating Losses
            PL_dB = A + 10*n*log10(d/d0) + Xf + Xh + S;
            shadowing_dB = S*randn(1); %dB
            Losses_dB = PL_dB + shadowing_dB - Gt;
            x=Losses_dB;
        end
    end
    

       
end