clear;
Time = xlsread('data.xls','A2:A1000');
pPSM_ppERK = xlsread('data.xls','B2:B1000');
pPSM_clock = xlsread('data.xls','C2:C1000');
opol = 10;
[p,s,mu] = polyfit(Time,pPSM_ppERK,opol);
f_y = polyval(p,Time,[],mu);
dt_erk = pPSM_ppERK - f_y;
plot(Time(1:645),dt_erk(1:645)), grid
title 'Detrended ERK Activity', ylabel 'Signal (a.u.)'
xlim([0 5])

figure
[upp,lop]=envelope(dt_erk,99,'rms');
normdt_erk=(dt_erk-lop)./(upp-lop);
[upp2,lop2]=envelope(normdt_erk,44,'peak');
normdt_erk2=(normdt_erk-lop2)./(upp2-lop2);
plot(Time(1:645),normdt_erk2(1:645)), grid
title 'Normalized Detrended ERK Activity', ylabel 'Signal (a.u.)'
ylim([-0.05 1.05])
xlim([0 5])

figure
[c,s,mu] = polyfit(Time,pPSM_clock,opol);
g_y = polyval(c,Time,[],mu);
dt_clock = pPSM_clock - g_y;
plot(Time(1:645),dt_clock(1:645)), grid
title 'Detrended Her7-Venus', ylabel 'Signal (a.u.)'
xlim([0 5])

figure
[upc,loc]=envelope(dt_clock,99,'rms');
normdt_clock=(dt_clock-loc)./(upc-loc);
[upc2,loc2]=envelope(normdt_clock,44,'peak');
normdt_clock2=(normdt_clock-loc2)./(upc2-loc2);

plot(Time(1:645),normdt_clock2(1:645)), grid
title 'Normalized Detrended Her7-Venus', ylabel 'Signal (a.u.)'
ylim([-0.05 1.05])
xlim([0 5])

figure
plot(Time(1:645),normdt_erk2(1:645),'g--*',Time(1:645),normdt_clock2(1:645),'m--x'), grid
ylim([-0.1 1.1])
xlim([0 5])

erktrend=f_y;
clocktrend=(g_y-50)/0.062271;

waveerk=2*normdt_erk2-1;
waveerk=min(1,waveerk);
waveerk=max(-1,waveerk);
phaseerk=acos(waveerk);

switchphase1=fix(max(asin(waveerk.^2)-pi()/2==phaseerk,asin(waveerk.^2)+pi()/2==phaseerk));
switchphase1=(circshift(switchphase1,-1)>switchphase1);

phaseerk2=phaseerk;
flip=0;
for k=2:length(Time)
    if switchphase1(k-1)
        flip=~flip;
    end
    if flip
        phaseerk2(k)=asin(waveerk(k))-pi()/2;
    end
end


waveclock=2*normdt_clock2-1;
waveclock=min(1,waveclock);
waveclock=max(-1,waveclock);
phaseclock=acos(waveclock);

switchphase1=fix(max(asin(waveclock.^2)-pi()/2==phaseclock,asin(waveclock.^2)+pi()/2==phaseclock));
switchphase1=(circshift(switchphase1,-1)>switchphase1);

phaseclock2=phaseclock;
flip=0;
for k=2:length(Time)
    if switchphase1(k-1)
        flip=~flip;
    end
    if flip
        phaseclock2(k)=asin(waveclock(k))-pi()/2;
    end
end
figure
plot(Time,phaseerk2+pi(),'g',Time,phaseclock2+pi(),'m')

phasediff=(phaseerk2-phaseclock2);
phasediff=phasediff-(sign(phasediff)-1)*pi()
figure
[deltaph,frame]=resample(phasediff,Time,10);
plot(frame,deltaph,'ro',Time,phasediff,'g-')
ylim([0,2*pi()])
