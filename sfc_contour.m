abcds=SFC(:,T-380:T-14)*100;
[X Y] = meshgrid(T-380:T-14,1:34);
[xx yy] = meshgrid(T-380:71/30:T-14,1:1/8:34);
abcs=interp2(X,Y,abcds,xx,yy);
figure
colormap(jet(64));
[Ms,cs]=contourf(abcs,[22 22]);

ylabel('PSM (um)');

% Create xlabel
xlabel('Time (min)');
% Create title
%title(strcat('SFC with clock T=',num2str(30*perc),'min'),'FontName','Arial', 'FontSize', 28);
title(strcat('SFC inh=',num2str(round(su-1,1)),' pulse=',num2str(round(30*perc/pulse,1)),"'-",num2str(round(30*perc*(1-1/pulse),1)),"'"),'FontName','Arial', 'FontSize', 22);


thres=0.1;
abcdp=PPMatrix(:,T-380:T-14);
abcp=interp2(X,Y,abcdp,xx,yy);
figure
colormap(jet(64));
[Mp,cp]=contourf(abcp,[thres thres]);

ylabel('PSM (um)');

% Create xlabel
xlabel('Time (min)');
% Create title
%title(strcat('ppERK with clock T=',num2str(30*perc),'min'),'FontName','Arial', 'FontSize', 28);
title(strcat('ppERK EC',num2str(thres*100),' inh=',num2str(round(su-1,1)),' pulse=',num2str(round(30*perc/pulse,1)),"'-",num2str(round(30*perc*(1-1/pulse),1)),"'"),'FontName','Arial', 'FontSize', 22);


thres=0.5;
abcdp=PPMatrix(:,T-380:T-14);
abcp=interp2(X,Y,abcdp,xx,yy);
figure
colormap(jet(64));
[Mp,cp]=contourf(abcp,[thres thres]);

ylabel('PSM (um)');

% Create xlabel
xlabel('Time (min)');
% Create title
%title(strcat('ppERK with clock T=',num2str(30*perc),'min'),'FontName','Arial', 'FontSize', 28);
title(strcat('ppERK EC',num2str(thres*100),' inh=',num2str(round(su-1,1)),' pulse=',num2str(round(30*perc/pulse,1)),"'-",num2str(round(30*perc*(1-1/pulse),1)),"'"),'FontName','Arial', 'FontSize', 22);




%{

MovingClock = nan * ones(K, T); 
for tf = 1:T
    if mod(tf, floor(1 / Vg)) == 0
            MovingClock(:, tf + 1 - floor(1 / Vg):tf)=CRMatrix(:, tf + 1 - floor(1 / Vg):tf);
            MovingClock = circshift(MovingClock, [1 0 0]);
            MovingClock(1, :, :)=MovingClock(2,:,:);
    end
end

cmatrix = zeros(K + 1 - tailbud, T);
        for i = 1:K + 1 - tailbud
            for j = 1:T
                cmatrix(i, j + 1) = MovingClock(tailbud + i - 1, j);
            end
        end

    MovingClock=cmatrix;
for i=1:size(MovingClock,2)
MovingClock(:,i)=MovingClock(:,i)/max(MovingClock(1,:));
end

for j=T+1:-1:T-380
MovingClock(1:floor((T+1-j)*Vg),j)=NaN;
end

abcdl=MovingClock(:,T-380:T-14);
abcl=interp2(X,Y,abcdl,xx,yy);

figure
colormap(jet(64));
[Mc,cc]=contourf(abcl,[0.95 0.95]);


ylabel('PSM (um)');

% Create xlabel
xlabel('Time (min)');
% Create title
title(strcat('Clock T=',num2str(30*perc),'min'),'FontName','Arial', 'FontSize', 28);

%}
%{
for t=T-380:T-1
    if(max(abcd(abs(MovingClock(:,t)-1)<0.05,t-T+381)>22))
        ;
end
%}