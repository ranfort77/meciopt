clear all; clc; close all; 

fname = 'test01_h2s/h2s_avdz_hlf_oldg.out'; select_k = 4;

% fname = 'test02_ethylene/c2h4_etzm_hlf_oldg.out'; select_k = 6;
% -------------------------------------------------------------------
ap = 3.5:0.5:14.0;
bt = 0.001:0.001:0.020;
[AP, BT] = meshgrid(ap, bt);
cut_bt = 0.02;
cut_ap = 3.5;
% -------------------------------------------------------------------
f = fopen(fname, 'r');
k = 0;
while ~feof(f)
    line = fgetl(f);

    if strfind(line, 'calculate gradients using central')
        k = k + 1;
        if k == select_k
            for n = 1:6; fgetl(f); end
            info = textscan(f, '%f%f%f%f');
            sg = info{3};
            gm = info{4};
            fgetl(f); 
            info = textscan(f, '%s%f%f%f%f');
            fgetl(f);
            info = textscan(f, '%s%f%f%f');
            sgp = info{4};
            fgetl(f);
            info = textscan(f, '%s%f%f%f');
            gmp = info{4};
        end
    end
end
% -------------------------------------------------------------------
% scalar functions
% -------------------------------------------------------------------
funP = @(gm, bt) gm^2 /(gm + bt);
funPp = @(gm,gmp,bt) (gmp*(gm^2 + 2*bt*gm)) / ((gm + bt)^2);
funL = @(sg, gm, ap, bt) sg + ap*funP(gm, bt);
funLp = @(gm, sgp, gmp, ap, bt) sgp + ap*funPp(gm,gmp,bt);
funU = @(gm,gmp,bt) funPp(gm,gmp,bt) / norm(funPp(gm,gmp,bt));
funpara = @(Lp, U, ap) dot(Lp, U) / ap;
funperp = @(Lp, U) norm(Lp - dot(Lp, U)*U);


% -------------------------------------------------------------------
% calculation data
% -------------------------------------------------------------------
L = nan*ones(size(AP));
P = nan*ones(size(AP));
PARA = nan*ones(size(AP));
PERP = nan*ones(size(AP));
for n = 1:numel(AP)
    P(n) = funP(gm, BT(n));
    L(n) = funL(sg, gm, AP(n), BT(n));
    elp = funLp(gm, sgp, gmp, AP(n), BT(n));
    u = funU(gm, gmp, BT(n));
    PARA(n) = funpara(elp, u, AP(n));
    PERP(n) = funperp(elp, u);
end
% -------------------------------------------------------------------
PARA = abs(PARA);
THRESH = 0.005*ones(size(PARA));
figure('color','w','pos',[100   300   400*3   400]); 
subplot(131)
mesh(AP, BT, PARA,'edgecolor','b'); hold on; box on;
mesh(AP, BT, PERP,'edgecolor','r'); hold on; box on;
mesh(AP, BT, THRESH, 'edgecolor','k','linestyle','-'); hold on; box on;
xlim([ap(1), ap(end)])
ylim([bt(1), bt(end)])
zlim([0, max([PARA(:); PERP(:)])])
xlabel('\alpha')
ylabel('\beta')
set(gca,'fontsize',14)
view([-51,20])
% index calc: NOTICE!!!!!!! meshgrid index
ii = find(abs(ap - cut_ap) < 1e-12);
jj = find(abs(bt - cut_bt) < 1e-12);
% frozen beta, alpha VS. para, perp
c1ap = AP(jj,:);
c1bt = BT(jj,:);
c1para = PARA(jj,:);
c1perp = PERP(jj,:);
c1thresh = THRESH(jj,:);
% frozen alpha, beta VS. para, perp
c2ap = AP(:,ii);
c2bt = BT(:,ii);
c2para = PARA(:,ii);
c2perp = PERP(:,ii);
c2thresh = THRESH(:,ii);
% plot3
h(1) = plot3(c1ap,c1bt,c1para,'bo-','MarkerFaceColor','b','MarkerSize',4);
h(2) = plot3(c1ap,c1bt,c1perp,'ro-','MarkerFaceColor','r','MarkerSize',4);
h(1) = plot3(c2ap,c2bt,c2para,'b^-','MarkerFaceColor','b','MarkerSize',4);
h(2) = plot3(c2ap,c2bt,c2perp,'r^-','MarkerFaceColor','r','MarkerSize',4);

subplot(132)
h(1) = plot(c1ap,c1para,'bo-','MarkerFaceColor','b','MarkerSize',4); hold on; grid on;
h(2) = plot(c1ap,c1perp,'ro-','MarkerFaceColor','r','MarkerSize',4);
h(2) = plot(c1ap,c1thresh,'k-','MarkerFaceColor','r','MarkerSize',4);
xlabel('\alpha')
set(gca,'xlim', [ap(1), ap(end)]);
set(gca,'fontsize',14)
title(sprintf('Fronzen \\beta = %8.4f', c1bt(1)))
legend('para','perp')

subplot(133)
h(1) = plot(c2bt,c2para,'bo-','MarkerFaceColor','b','MarkerSize',4); hold on; grid on;
h(2) = plot(c2bt,c2perp,'ro-','MarkerFaceColor','r','MarkerSize',4);
h(2) = plot(c2bt,c2thresh,'k-','MarkerFaceColor','r','MarkerSize',4);
xlabel('\beta')
set(gca,'xlim', [bt(1), bt(end)]);
set(gca,'fontsize',14)
title(sprintf('Fronzen \\alpha = %8.4f', c2ap(1)))
legend('para','perp')






