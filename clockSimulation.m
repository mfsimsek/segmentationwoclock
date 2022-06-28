%%2-D simulation using a hexagonal structure%%
clear;clc;
clockamp=20; % clockamp=20 for clock, 0 for w/o clock
drug=1.5; % fraction of inhibitor effect on ppp % drug=1.5 for treatment, 1 for control
pulse=6; % pulse freq in each cycle

%%space-time settings%%
L=40; % PSM Size in Cell #
prec=30;
perc=1.333;
delt=prec/71/perc;
tcell=6;
x=1:L;
R = 4; % PSM Width in Cell #
tailbud = 7; % # of Tailbud Cells
T = 2218; % Simulation Duration

randDirection = zeros(10^8, 1); % pre-allocates random direction matrix

%%Clock Matrix Settings%%
CRMatrix = zeros(L,T); % empty clock matrix
phi=zeros(L,T); %clock phase matrix
f=1./(1+exp(5*(x-35)/L));
f=2*pi/prec*f;
for i=1:T
    
    CRMatrix(:,i)=clockamp*(sin(phi(:,i))+1)/2; % clock signal
    phi(:,i+1)=phi(:,i)+delt*f(:);
    
    if mod(i,round(tcell/delt,0))==0
        phi=circshift(phi,1,1);
        phi(1,:)=phi(2,:);
    end
end

%%Simulation%%
tic

pb = .01; % Ligand Receptor Binding Rate
Diff = 2; % Ligand Diffusion Speed


ppd  = 1.0; % Dephosphorylation Rate due to Repressor
pcd = 0.5; % Dephosphorylation Rate due to Clock
ppp = 1.0; % Phosphorylation Rate

pppm=ppp*ones(1,T);
for t=1:T
    if mod(t,71*perc)<floor(71*perc/pulse)
        pppm(t)=pppm(t)/drug;
    end
end

prp = 1.2; % Repressor Protein Translation Rate
prr = -0.5; % Repressor Protein Decay / Translation Rate


actdel = 14; % Activation Delay from Complex to ERK
trlatedel = 14; % Trl+Secr Delay of FGF
reprdel = 35; % Repression (Trx+Trl) Delay

randDirection = randi(3, [10^8 1]); % direction for neighborFinder, only enough for 1 parameter set
dirCounter = 0; % Chooses which random num

Vg = 0.07; % Tailbud Growth Speed

pbd = -1; % Comp Decay Rate
rd = -0.025; % Ligand RNA Decay Rate
pfp = 1; % Free Ligand Protein Translation Rate
pfd = -0.25; % Free Ligand Protein Decay Rate

pd = -0.1; % Inactive Signaling Protein Decay Rate
pp = 12.5; % Inactive Signaling Protein Translation Rate
rp = 10; % Ligand RNA Transcription Rate

r0 = 0; % Initial Ligand RNA
pf0 = 0; % Initial Free Ligand Protein
pb0 = 0; % Initial Bound Ligand Protein
pr0 = 0; % Initial Repressor Protein
pp0 = 0; % Initial Active Signaling Protein
p0 = 0; % Initial Inactive Signaling Protein

RMatrix = r0 * ones(L, T, R); % Ligand RNA Matrix (mLIG)
PFMatrix = r0 * ones(L, T, R); % Ligand Protein Matrix (LIG)
PBMatrix = r0 * ones(L, T, R); % Bound Ligand Receptor Protein Matrix (COMP)
PRMatrix = r0 * ones(L, T, R); % Repressor Protein Matrix (INH)
PMatrix = r0 * ones(L, T, R); % Inactive Signaling Protein Matrix (SIG-)
PPMatrix = r0 * ones(L, T, R); % Active Signaling Protein Matrix (SIG+)
StatMatrix = PPMatrix; % Stationary Matrix

ry0matrix = r0 * ones(L, 1, R);
pfy0matrix = pf0 * ones(L, 1, R);
pby0matrix = pb0 * ones(L, 1, R);
pry0matrix = pr0 * ones(L, 1, R);
ppy0matrix = pp0 * ones(L, 1, R);
y0matrix = p0 * ones(L, 1, R);

for ti = 0:T-1
    %Euler's Method for dy/dt=f(y,t)
    N = 100;
    for j = 1:R
        for k = 1:tailbud
            
            tf = ti + 1; % initial and final ti
            h = (tf - ti) / N; % Time step
            ry = zeros(N + 1, 1);
            pfy = zeros(N + 1, 1);
            pby = zeros(N + 1, 1);
            pry = zeros(N + 1, 1);
            py = zeros(N + 1, 1);
            ppy = zeros(N + 1, 1);
            
            if ti == 0
                
                ry(1) = ry0matrix(k, 1, j);
                pfy(1) = pfy0matrix(k, 1, j);
                pby(1) = pby0matrix(k, 1, j);
                pry(1)=pry0matrix(k,1,j);
                py(1) = y0matrix(k, 1, j);
                ppy(1) = ppy0matrix(k, 1, j);
                
            else
                
                ry(1) = RMatrix(k, ti, j);
                pfy(1) = PFMatrix(k, ti, j);
                pby(1) = PBMatrix(k, ti, j);
                pry(1) = PRMatrix(k, ti, j);
                py(1) = PMatrix(k, ti, j);
                ppy(1) = PPMatrix(k, ti, j);
                
            end
            
            for n = 1:N
                
                ry(n + 1) = ry(n) + h * rp + h * rd * ry(n);
                
                if ti > actdel && ti > trlatedel
                    
                    pfy(n + 1) = pfy(n) + h * (pfp * RMatrix(k, ti - trlatedel, j) + pfd * pfy(n) - pb * PFMatrix(k, ti - actdel, j));
                    pby(n + 1) = pby(n) + h * (pb * PFMatrix(k, ti - actdel, j) + pbd * pby(n));
                    
                else
                    
                    pfy(n + 1) = pfy(n);
                    pby(n + 1) = pby(n);
                    
                end
                
                if ti > reprdel
                    
                    pry(n + 1) = pry(n) + h * (prp * PPMatrix(k, ti - reprdel, j) + prr * prp * pry(n));
                    
                    if py(n) < ppy(n)
                        
                        ppy(n + 1) = ppy(n) + h * (pppm(ti) * pby(n) * py(n) + pd * ppy(n) - ppd * pry(n) * ppy(n) - pcd * CRMatrix(k, ti) * ppy(n));
                        py(n + 1) = py(n) + ppy(n) + h * (pp + pd * (py(n) + ppy(n))) - ppy(n + 1);
                        
                    else
                        
                        py(n + 1) = py(n) + h * (pp + pcd * CRMatrix(k, ti) * ppy(n) + ppd * pry(n) * ppy(n) + pd * py(n) - pppm(ti) * pby(n) * py(n));
                        ppy(n + 1) = ppy(n) + py(n) + h * (pp + pd * (ppy(n) + py(n))) - py(n + 1);
                        
                    end
                    
                    
                    if py(n + 1) < 0 || py(n + 1) > 50000 || ppy(n + 1) < 0 || ppy(n + 1) > 50000
                        
                        py(n + 1) = 0;
                        ppy(n + 1) = 0;
                        
                    end
                    
                else
                    
                    pry(n + 1) = pry(n);
                    ppy(n + 1) = ppy(n);
                    py(n + 1) = py(n);
                    
                end
                
            end
            
            ry0 = ry(N + 1);
            pfy0 = pfy(N + 1);
            pby0 = pby(N + 1);
            pry0 = pry(N + 1);
            py0 = py(N + 1);
            ppy0 = ppy(N + 1);
            
            ry0matrix(k, 1, j) = ry0;
            pfy0matrix(k, 1, j) = pfy0;
            pby0matrix(k, 1, j) = pby0;
            pry0matrix(k, 1, j) = pry0;
            y0matrix(k, 1, j) = py0;
            ppy0matrix(k, 1, j) = ppy0;
            
            RMatrix(k, tf, j) = ry0;
            PFMatrix(k, tf, j) = pfy0;
            PBMatrix(k, tf, j) = pby0;
            PRMatrix(k, tf, j) = pry0;
            PMatrix(k, tf, j) = py0;
            PPMatrix(k, tf, j)=ppy0;
        end
    end
end


r0 = ry0matrix(tailbud, 1, j);
pf0 = pfy0matrix(tailbud, 1, j);
pb0 = pby0matrix(tailbud, 1, j);
pr0 = pry0matrix(tailbud, 1, j);
pp0 = ppy0matrix(tailbud, 1, j);
p0 = y0matrix(tailbud, 1, j);


for j = 1:R 
    for k = tailbud + 1:L
        
        N = 100; % Number of ti steps
        ry0 = r0;
        pfy0 = pf0;
        pby0 = pb0;
        pry0 = pr0;
        py0 = p0;
        ppy0 = pp0;
        
        past = floor((k - tailbud) / Vg);
        for ti = 1:past
            
            %Euler's Method for dy/dt=f(y,t)
            tf = ti + 1; % initial and final ti
            h = (tf - ti) / N; % Time step
            ry = zeros(N + 1, 1);
            pfy = zeros(N + 1, 1);
            pby = zeros(N + 1, 1);
            pry = zeros(N + 1, 1);
            py = zeros(N + 1, 1);
            ppy = zeros(N + 1, 1);
            
            ry(1) = ry0;
            pfy(1) = pfy0;
            pby(1) = pby0;
            pry(1) = pry0;
            py(1) = py0;
            ppy(1) = ppy0;
            
            for n = 1:N
                
                ry(n + 1) = ry(n) + h * rd * ry(n);
                
                if ti > actdel && ti > trlatedel
                    
                    pfy(n + 1) = pfy(n) + h * (pfp * RMatrix(k, ti - trlatedel, j) + pfd * pfy(n) - pb * PFMatrix(k, ti - actdel, j));
                    pby(n + 1) = pby(n) + h * (pb * PFMatrix(k, ti - actdel, j) + pbd * pby(n));
                    
                else
                    
                    pfy(n + 1) = pfy(n);
                    pby(n + 1) = pby(n);
                    
                end
                
                if ti > reprdel
                    
                    pry(n + 1) = pry(n) + h * (prp * PPMatrix(k, ti - reprdel, j) + prr * prp * pry(n));
                    
                    if py(n) < ppy(n)
                        
                        ppy(n + 1) = ppy(n) + h * (pppm(ti) * pby(n) * py(n) + pd * ppy(n) - ppd * pry(n) * ppy(n) - pcd * CRMatrix(k, ti) * ppy(n));
                        py(n + 1) = py(n) + ppy(n) + h * (pp + pd * (py(n) + ppy(n))) - ppy(n + 1);
                        
                    else
                        
                        py(n + 1) = py(n) + h * (pp + pcd * CRMatrix(k, ti) * ppy(n) + ppd * pry(n) * ppy(n) + pd * py(n) - pppm(ti) * pby(n) * py(n));
                        ppy(n + 1) = ppy(n) + py(n) + h * (pp + pd * (ppy(n) + py(n))) - py(n + 1);
                        
                    end
                else
                    
                    pry(n + 1) = pry(n);
                    ppy(n + 1) = ppy(n);
                    py(n + 1) = py(n);
                    
                end
                
            end
            
            ry0 = ry(N + 1);
            pfy0 = pfy(N + 1);
            pby0 = pby(N + 1);
            pry0 = pry(N + 1);
            py0 = py(N + 1);
            ppy0 = ppy(N + 1);
            
            ry0matrix(k, 1, j) = ry0;
            pfy0matrix(k, 1, j) = pfy0;
            pby0matrix(k, 1, j) = pby0;
            pry0matrix(k, 1, j) = pry0;
            y0matrix(k, 1, j) = py0;
            ppy0matrix(k, 1, j) = ppy0;
            
            RMatrix(k, tf, j) = ry0;
            PFMatrix(k, tf, j) = pfy0;
            PBMatrix(k, tf, j) = pby0;
            PRMatrix(k, tf, j) = pry0;
            PMatrix(k, tf, j) = py0;
            PPMatrix(k, tf, j) = ppy0;
            
        end
    end
end

for ti = 0:T - 1
    %Euler's Method for dy/dt=f(y,t)
    
    for j = 1:R
        for k = 1:L - 1
            
            tf = ti + 1; % initial and final ti
            h = (tf - ti) / N; % Time step
            ry = zeros(N + 1, 1);
            pfy = zeros(N + 1, 1);
            pby = zeros(N + 1, 1);
            pry = zeros(N + 1, 1);
            py = zeros(N + 1, 1);
            ppy = zeros(N + 1, 1);
            
            if ti == 0
                
                ry(1) = ry0matrix(k, 1, j);
                pfy(1) = pfy0matrix(k, 1, j);
                pby(1) = pby0matrix(k, 1, j);
                pry(1) = pry0matrix(k, 1, j);
                py(1) = y0matrix(k, 1, j);
                ppy(1) = ppy0matrix(k, 1, j);
                
            else
                
                ry(1) = RMatrix(k, ti, j);
                pfy(1) = PFMatrix(k, ti, j);
                pby(1) = PBMatrix(k, ti, j);
                pry(1) = PRMatrix(k,ti,j);
                py(1) = PMatrix(k, ti, j);
                ppy(1) = PPMatrix(k, ti, j);
            end
            
            for n = 1:N
                if k <= tailbud
                    ry(n + 1) = ry(n) + h * (rp + rd * ry(n));
                else
                    ry(n + 1) = ry(n) + h * rd * ry(n);
                end
                
                if ti > trlatedel && ti > actdel
                    value = pfy(n);
                    
                    dirCounter = dirCounter+1;
                    direction = randDirection(dirCounter);
                    c = ceil((k - 1) / L);
                    
                    switch direction % determines adjacent cells for hexagonal matrix
                        case 1 % 12 o'clock/6 o'clock
                            
                            if k == 1
                                neighbor1 = 0;
                            else
                                neighbor1 = PFMatrix(k - 1, ti, j); % 12 o'clock
                            end
                            
                            if k >= L
                                neighbor2 = 0;
                            else
                                neighbor2 = PFMatrix(k + 1, ti, j); % 6 o'clock
                            end
                            
                        case 2 % 2 o'clock/8 o'clock
                            
                            if (mod(j, 2) == 0 && k == 1) || j == R
                                neighbor1 = 0; % right (2 o'clock) not possible
                            elseif mod(j,2) == 0
                                neighbor1 = PFMatrix(k - 1, ti, j + 1); % 2 o'clock (even j)
                            else
                                neighbor1 = PFMatrix(k, ti, j + 1); % 2 o'clock (odd j)
                            end
                            
                            if j == 1 || (mod(j, 2) ~= 0 && k >= L)
                                neighbor2 = 0; % left (8 o'clock) not possible
                            elseif mod(j, 2) == 0
                                neighbor2 = PFMatrix(k, ti, j - 1); % 8 o'clock (even j)
                            else
                                neighbor2 = PFMatrix(k + 1, ti, j-1); % 8 o'clock (odd j)
                            end
                            
                        case 3 % 4 o'clock/10 o'clock
                            
                            if (mod(j, 2) ~= 0 && k >= L) || j == R
                                neighbor1 = 0; % right (4 o'clock) not possible
                            elseif mod(j, 2) == 0
                                neighbor1 = PFMatrix(k, ti, j + 1); % 4 o'clock (even j)
                            else
                                neighbor1 = PFMatrix(k + 1, ti, j + 1); % 4 o'clock (odd j)
                            end
                            
                            if j == 1 || (mod(j, 2) == 0 && k == 1)
                                neighbor2 = 0; % left (10 o'clock) not possible
                            elseif mod(j, 2) == 0
                                neighbor2 = PFMatrix(k - 1, ti, j - 1); % 10 o'clock (even j)
                            else
                                neighbor2 = PFMatrix(k, ti, j - 1); % 10 o'clock (odd j)
                            end
                    end
                    
                    switch c % k==1
                        case 0
                            if j == 1
                                switch direction
                                    case 1
                                        x = ((2 * neighbor2) - 2 * value);
                                    otherwise
                                        x = ((2 * neighbor1) - 2 * value);
                                end
                            elseif j == R
                                if mod(j, 2) == 0 && direction == 3
                                    x = 0;
                                else
                                    x = ((2 * neighbor2) - 2 * value);
                                end
                            else
                                if mod(j, 2) == 0
                                    switch direction
                                        case 3
                                            x = ((2 * neighbor1) - 2 *value);
                                        otherwise
                                            x = ((2 * neighbor2) - 2 * value);
                                    end
                                else
                                    switch direction
                                        case 1
                                            x = ((2 * neighbor2) - 2 * value);
                                        otherwise
                                            x = (neighbor1 + neighbor2 - 2 * value);
                                    end
                                end
                            end
                            
                        case 2 % k > L (k>M)
                            if j == 1
                                switch direction
                                    case 3
                                        x = (-2 * value);
                                    otherwise
                                        x = (neighbor1 - 2 * value);
                                end
                            elseif j == R
                                if mod(j, 2) == 0
                                    switch direction
                                        case 1
                                            x = (neighbor1 - 2 * value);
                                        otherwise
                                            x = (neighbor2 - 2 * value);
                                    end
                                else
                                    switch direction
                                        case 1
                                            x = (neighbor1 - 2 * value);
                                        case 2
                                            x = (-2 * value);
                                        case 3
                                            x = (neighbor2 - 2 * value);
                                    end
                                end
                            else
                                if mod(j, 2) == 0
                                    switch direction
                                        case 1
                                            x = (neighbor1 - 2 * value);
                                        otherwise
                                            x = (neighbor1 + neighbor2 - 2 * value);
                                    end
                                else
                                    switch direction
                                        case 3
                                            x = (neighbor2 - 2 * value);
                                        otherwise
                                            x = (neighbor1 - 2 * value);
                                    end
                                end
                            end
                            
                        otherwise % k is between 1 and L-1
                            if j == 1
                                switch direction
                                    case 1
                                        x = (neighbor1 + neighbor2 - 2 * value);
                                    otherwise
                                        x = (2 * neighbor1 - 2 * value);
                                end
                            elseif j == R
                                switch direction
                                    case 1
                                        x = (neighbor1 + neighbor2 - 2 * value);
                                    otherwise
                                        x = (2 * neighbor2 - 2 * value);
                                end
                            else
                                x = (neighbor1 + neighbor2 - 2 * value);
                            end
                            
                    end
                    pfy(n + 1) = pfy(n) + h * (pfp * RMatrix(k, ti - trlatedel, j) + pfd * pfy(n) - pb * PFMatrix(k, ti - actdel, j) + Diff * x);
                    
                else
                    pfy(n + 1) = pfy(n);
                end
                
                if ti > actdel
                    pby(n + 1) = pby(n) + h * (pb * PFMatrix(k, ti - actdel, j) + pbd * pby(n));
                else
                    pby(n + 1) = pby(n);
                end
                
                if ti > reprdel
                    pry(n + 1) = pry(n) + h * (prp * PPMatrix(k, ti - reprdel, j) + prr*prp * pry(n));
                    
                    if py(n) < ppy(n)
                        ppy(n+1)=ppy(n)+h*(pppm(ti)*pby(n)*py(n)+pd*ppy(n)-ppd*pry(n)*ppy(n) - pcd * CRMatrix(k,ti) * ppy(n));
                        py(n+1)=py(n) + ppy(n)+h*(pp+pd*(py(n)+ppy(n))) - ppy(n+1);
                        
                    else
                        py(n+1)=py(n)+h*(pp + pcd * CRMatrix(k,ti) * ppy(n) + ppd*pry(n)*ppy(n)+pd*py(n)-pppm(ti)*pby(n)*py(n));
                        ppy(n+1)=ppy(n)+py(n)+h*(pp + pd*(ppy(n)+py(n))) - py(n+1);
                    end
                    
                else
                    pry(n + 1) = pry(n);
                    ppy(n + 1) = ppy(n);
                    py(n + 1) = py(n);
                end
                
            end
            
            
            
            ry0 = ry(N + 1);
            pfy0 = pfy(N + 1);
            pby0 = pby(N + 1);
            pry0 = pry(N + 1);
            py0 = py(N + 1);
            ppy0 = ppy(N + 1);
            
            ry0matrix(k, 1, j) = ry0;
            pfy0matrix(k, 1, j) = pfy0;
            pby0matrix(k, 1, j) = pby0;
            pry0matrix(k, 1, j) = pry0;
            y0matrix(k, 1, j) = py0;
            ppy0matrix(k, 1, j) = ppy0;
            
            RMatrix(k, tf, j) = ry0;
            PFMatrix(k, tf, j) = pfy0;
            PBMatrix(k, tf, j) = pby0;
            PRMatrix(k, tf, j) = pry0;
            PMatrix(k, tf, j) = py0;
            PPMatrix(k, tf, j) = ppy0;
            
        end
        k = L;
        if ti > 0
            tf = ti + 1; % initial and final ti
            h = (tf - ti) / N; % Time step
            ry = zeros(N + 1, 1);
            pfy = zeros(N + 1, 1);
            pby = zeros(N + 1, 1);
            pry = zeros(N + 1, 1);
            py = zeros(N + 1, 1);
            ppy = zeros(N + 1, 1);
            
            
            ry(1) = RMatrix(k, ti, j);
            pfy(1) = PFMatrix(k, ti, j);
            pby(1) = PBMatrix(k, ti, j);
            pry(1) = PRMatrix(k, ti, j);
            py(1) = PMatrix(k, ti, j);
            ppy(1) = PPMatrix(k, ti, j);
            
            
            for n = 1:N
                
                ry(n + 1) = ry(n) + h * rd * ry(n);
                
                if ti > trlatedel && ti > actdel
                    value = pfy(n);
                    dirCounter = dirCounter + 1;
                    direction = randDirection(dirCounter);
                    c = ceil((k - 1) / L);
                    
                    switch direction % determines adjacent cells for hexagonal matrix
                        case 1 % 12 o'clock/6 o'clock
                            
                            if k == 1
                                neighbor1 = 0;
                            else
                                neighbor1 = PFMatrix(k - 1, ti, j); % 12 o'clock
                            end
                            
                            if k >= L
                                neighbor2 = 0;
                            else
                                neighbor2 = PFMatrix(k + 1, ti, j); % 6 o'clock
                            end
                            
                        case 2 % 2 o'clock/8 o'clock
                            
                            if (mod(j, 2) == 0 && k == 1) || j == R
                                neighbor1 = 0; % right (2 o'clock) not possible
                            elseif mod(j, 2) == 0
                                neighbor1 = PFMatrix(k - 1, ti, j + 1); % 2 o'clock (even j)
                            else
                                neighbor1 = PFMatrix(k, ti, j + 1); % 2 o'clock (odd j)
                            end
                            
                            if j == 1 || (mod(j, 2) ~= 0 && k >= L)
                                neighbor2 = 0; % left (8 o'clock) not possible
                            elseif mod(j, 2) == 0
                                neighbor2 = PFMatrix(k, ti, j - 1); % 8 o'clock (even j)
                            else
                                neighbor2 = PFMatrix(k + 1, ti, j-1); % 8 o'clock (odd j)
                            end
                            
                        case 3 % 4 o'clock/10 o'clock
                            
                            if (mod(j, 2) ~= 0 && k >= L) || j == R
                                neighbor1 = 0; % right (4 o'clock) not possible
                            elseif mod(j, 2) == 0
                                neighbor1 = PFMatrix(k, ti, j + 1); % 4 o'clock (even j)
                            else
                                neighbor1 = PFMatrix(k + 1, ti, j + 1); % 4 o'clock (odd j)
                            end
                            
                            if j == 1 || (mod(j, 2) == 0 && k == 1)
                                neighbor2 = 0; % left (10 o'clock) not possible
                            elseif mod(j, 2) == 0
                                neighbor2 = PFMatrix(k - 1, ti, j - 1); % 10 o'clock (even j)
                            else
                                neighbor2 = PFMatrix(k, ti, j - 1); % 10 o'clock (odd j)
                            end
                    end
                    
                    switch c % k==1
                        case 0
                            if j == 1
                                switch direction
                                    case 1
                                        x = ((2 * neighbor2) - 2 * value);
                                    otherwise
                                        x = ((2 * neighbor1) - 2 * value);
                                end
                            elseif j == R
                                if mod(j, 2) == 0 && direction == 3
                                    x = 0;
                                else
                                    x = ((2 * neighbor2) - 2 * value);
                                end
                            else
                                if mod(j, 2) == 0
                                    switch direction
                                        case 3
                                            x = ((2 * neighbor1) - 2 * value);
                                        otherwise
                                            x = ((2 * neighbor2) - 2 * value);
                                    end
                                else
                                    switch direction
                                        case 1
                                            x = ((2 * neighbor2) - 2 * value);
                                        otherwise
                                            x = (neighbor1 + neighbor2 - 2 * value);
                                    end
                                end
                            end
                            
                        case 2 % k > L (k>M)
                            if j == 1
                                switch direction
                                    case 3
                                        x = (-2 * value);
                                    otherwise
                                        x = (neighbor1 - 2 * value);
                                end
                            elseif j == R
                                if mod(j, 2) == 0
                                    switch direction
                                        case 1
                                            x = (neighbor1 - 2 * value);
                                        otherwise
                                            x = (neighbor2 - 2 * value);
                                    end
                                else
                                    switch direction
                                        case 1
                                            x = (neighbor1 - 2 * value);
                                        case 2
                                            x = (-2 * value);
                                        case 3
                                            x = (neighbor2 - 2 * value);
                                    end
                                end
                            else
                                if mod(j, 2) == 0
                                    switch direction
                                        case 1
                                            x = (neighbor1 - 2 * value);
                                        otherwise
                                            x = (neighbor1 + neighbor2 - 2 * value);
                                    end
                                else
                                    switch direction
                                        case 3
                                            x = (neighbor2 - 2 * value);
                                        otherwise
                                            x = (neighbor1 - 2 * value);
                                    end
                                end
                            end
                            
                        otherwise % k is between 1 and L-1
                            if j == 1
                                switch direction
                                    case 1
                                        x = (neighbor1 + neighbor2 - 2 * value);
                                    otherwise
                                        x = (2 * neighbor1 - 2 * value);
                                end
                            elseif j == R
                                switch direction
                                    case 1
                                        x = (neighbor1 + neighbor2 - 2 * value);
                                    otherwise
                                        x = (2 * neighbor2 - 2 * value);
                                end
                            else
                                x = (neighbor1 + neighbor2 - 2 * value);
                            end
                            
                    end
                    pfy(n + 1) = pfy(n) + h * (pfp * RMatrix(k, ti - trlatedel, j) + pfd * pfy(n) - pb * PFMatrix(k, ti - actdel, j) + Diff * x);
                    
                else
                    pfy(n + 1) = pfy(n);
                end
                
                if ti > actdel
                    pby(n +1 ) = pby(n) + h * (pb * PFMatrix(k, ti - actdel, j) + pbd * pby(n));
                else
                    pby(n + 1) = pby(n);
                end
                
                if ti > reprdel
                    pry(n + 1) = pry(n) + h * (prp * PPMatrix(k, ti - reprdel, j) + prr*prp * pry(n));
                    
                    if py(n) < ppy(n)
                        ppy(n+1)=ppy(n)+h*(pppm(ti)*pby(n)*py(n)+pd*ppy(n)-ppd*pry(n)*ppy(n) - pcd * CRMatrix(k,ti) * ppy(n));
                        py(n+1)=py(n) + ppy(n)+h*(pp+pd*(py(n)+ppy(n))) - ppy(n+1);
                        
                    else
                        py(n+1)=py(n)+h*(pp + pcd * CRMatrix(k,ti) * ppy(n) + ppd*pry(n)*ppy(n)+pd*py(n)-pppm(ti)*pby(n)*py(n));
                        ppy(n+1)=ppy(n)+py(n)+h*(pp + pd*(ppy(n)+py(n))) - py(n+1);
                    end
                    
                else
                    pry(n + 1) = pry(n);
                    ppy(n + 1) = ppy(n);
                    py(n + 1) = py(n);
                end
                
            end
            
            
            
            ry0 = ry(N + 1);
            pfy0 = pfy(N + 1);
            pby0 = pby(N + 1);
            pry0 = pry(N + 1);
            py0 = py(N + 1);
            ppy0 = ppy(N + 1);
            
            ry0matrix(k, 1, j) = ry0;
            pfy0matrix(k, 1, j) = pfy0;
            pby0matrix(k, 1, j) = pby0;
            pry0matrix(k, 1, j) = pry0;
            y0matrix(k, 1, j) = py0;
            ppy0matrix(k, 1, j) = ppy0;
            
            RMatrix(k, tf, j) = ry0;
            PFMatrix(k, tf, j) = pfy0;
            PBMatrix(k, tf, j) = pby0;
            PRMatrix(k, tf, j) = pry0;
            PMatrix(k, tf, j) = py0;
            PPMatrix(k, tf, j) = ppy0;
            
        end
    end
    if mod(tf, floor(1 / Vg)) == 0
        StatMatrix(:, tf + 1 - floor(1 / Vg):tf, :) = PPMatrix(:, tf + 1 - floor(1 / Vg):tf, :);
        RMatrix = circshift(RMatrix, [1 0 0]); % shifts matrix 1 position along 1st dimension
        PFMatrix = circshift(PFMatrix, [1 0 0]);
        PBMatrix = circshift(PBMatrix, [1 0 0]);
        PRMatrix = circshift(PRMatrix, [1 0 0]);
        PMatrix = circshift(PMatrix, [1 0 0]);
        PPMatrix = circshift(PPMatrix, [1 0 0]);
        
        PFMatrix(1, :, :) = PFMatrix(2,:,:); % accounts for matrix shift
        RMatrix(1, :, :)=RMatrix(2,:,:);
        PRMatrix(1, :, :)=PRMatrix(2,:,:);
        PBMatrix(1, :, :)=PBMatrix(2,:,:);
        PMatrix(1, :, :)=PMatrix(2,:,:);
        PPMatrix(1, :, :)=PPMatrix(2,:,:);
    end
end

A = PPMatrix(1, T, 1) / 2;

if PPMatrix(tailbud+1, T, 1) > A && PPMatrix(37, T, 1) < A
    zmatrix = zeros(L + 1 - tailbud, T);
    for i = 1:L + 1 - tailbud
        for j = 1:T
            zmatrix(i, j + 1) = PPMatrix(tailbud + i - 1, j)- 0.05 * PPMatrix(tailbud, j);
        end
    end
    xnew = linspace(1, L + 1 - tailbud, 6 * (L + 1 - tailbud) + 1);
    znew = interp1(zmatrix, xnew, 'linear');
    
    graph=znew(:, T, 1)';
    
end

toc
figure
mesh(StatMatrix(:,T-380:T,1))

figure
PPMatrix=zmatrix;
for i=1:size(PPMatrix,2)
    PPMatrix(:,i)=PPMatrix(:,i)/max(PPMatrix(1,:));
end

for j=T+1:-1:T-380
    PPMatrix(1:floor((T+1-j)*Vg),j)=NaN;
end

mesh(PPMatrix(:,T-380:T))

figure
NormMatrix=zmatrix;
for i=1:size(NormMatrix,2)
    NormMatrix(:,i)=NormMatrix(:,i)/NormMatrix(1,i);
end

for j=T+1:-1:T-380
    NormMatrix(1:floor((T+1-j)*Vg),j)=NaN;
end

mesh(NormMatrix(:,T-380:T))

figure
SFC=zmatrix;
SFC=SFC*0;
for i=1:size(zmatrix,2)
    j=1;
    while SFC(j,i)<=1 && j<size(zmatrix,1)
        j=j+1;
        SFC(j,i)=zmatrix(j-1,i)/zmatrix(j,i)-1;
    end
    SFC(1,i)=SFC(2,i);
    siz=size(SFC(SFC(:,i)<=0,i));
    x=NaN(siz);
    SFC(SFC(:,i)<=0,i)=x;
    
    siz1=size(SFC(SFC(:,i)>=1,i));
    y=ones(siz1);
    SFC(SFC(:,i)>=1,i)=y;
    
end

for j=T+1:-1:T-380
    SFC(1:floor((T+1-j)*Vg),j)=NaN;
end

mesh(SFC(:,T-380:T))
