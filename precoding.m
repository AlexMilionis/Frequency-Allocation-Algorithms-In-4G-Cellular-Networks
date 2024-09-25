function P=precoding(Pt)

%% Precoding parametres
% theoro streams=layers=sinoliko data stream ston receiver: 
%U=2;
% thero antenna ports 
%p=4;

Nt=4;
Nr=2;

%% %% Beamforming matrix W (p x u, u=2 layers, p=4 antenna ports)

%W = I - 2*u*(u^H)/((u^H)*u)

u=cell(16,1);
u{1,1}=[1 -1 -1 -1].';                           %0
u{2,1}=[1 -1i 1 1i].';                           %1
u{3,1}=[1 1 -1 1].';                             %2
u{4,1}=[1 1i 1 -1i].';                           %3
u{5,1}=[1 (-1-1i)/sqrt(2) -1i (1-1i)/sqrt(2)].'; %4
u{6,1}=[1 (1-1i)/sqrt(2) 1i (-1-1i)/sqrt(2)].';  %5
u{7,1}=[1 (1+1i)/sqrt(2) -1i (-1+1i)/sqrt(2)].'; %6
u{8,1}=[1 (-1+1i)/sqrt(2) 1i (1+1i)/sqrt(2)].';  %7  
u{9,1}=[1 -1 1 1].';                             %8
u{10,1}=[1 -1i -1 -1i].';                        %9
u{11,1}=[1 1 1 -1].';                            %10
u{12,1}=[1 1i -1 1i].';                          %11
u{13,1}=[1 -1 -1 1].';                           %12
u{14,1}=[1 -1 1 -1].';                           %13
u{15,1}=[1 1 -1 -1].';                           %14
u{16,1}=[1 1 1 1].';                             %15

Wn=cell(16,1);
for i=1:16
    Wn{i,1}=eye(4)- (2*u{i,1}*(u{i,1}'))/((u{i,1}')*u{i,1});
end


W=cell(16,1);
tempW=zeros(4,2);
for i=1:16
    if i==1 || i==5 || i==6 || i==10
        tempWn=Wn{i,1};
        tempW(:,1)=tempWn(:,1);
        tempW(:,2)=tempWn(:,4);
        W{i,1}=tempW/sqrt(2);
    end
    if i==2 || i==3 || i==4 || i==9 || i==13 || i==16
        tempWn=Wn{i,1};
        tempW(:,1)=tempWn(:,1);
        tempW(:,2)=tempWn(:,2);
        W{i,1}=tempW/sqrt(2);
    end
    if i==7 || i==8 || i==11 || i==12 || i==14 || i==15
        tempWn=Wn{i,1};
        tempW(:,1)=tempWn(:,1);
        tempW(:,2)=tempWn(:,3);
        W{i,1}=tempW/sqrt(2);
    end
end

%% Power allocation matrix S (2x2)

S = [sqrt(Pt/Nt/2) 0 ; 0 sqrt(Pt/Nt/2)];

%% Precoding matrix P  (4x2)

P=cell(16,1);
for i=1:16
    P{i,1} = W{i,1}*sqrtm(S);
end

% to thema omws einai pio apo ta 16 W(i) tha epilekso
% h 3gpp den proteinei tropo gia non CDD precoding, opote paw gia H*P->1
    


end