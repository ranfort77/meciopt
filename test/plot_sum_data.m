clear all; clc; close all; 

% prefix = 'test01_h2s/h2s_avdz_hlf_oldg'; ylim3 = [0, 0.006]; ylim1 = [-398.6613, -398.6595]; ylim2 = [0, 0.03];

prefix = 'test02_ethylene/c2h4_etzm_hlf_oldg'; ylim3 = [0, 0.025];


au2ev = 27.2107;

% -------------------------------------------------------------------
% parsing thresholds
% -------------------------------------------------------------------
fname = sprintf('%s.out', prefix);
fid = fopen(fname);
while ~feof(fid)
    sline = fgetl(fid);
    if strfind(sline, 'THRESH_STEP');
        s = regexp(sline, '=', 'split');
        thresh_step = str2double(s{2});
        sline = fgetl(fid);
        s = regexp(sline, '=', 'split');
        thresh_grad = str2double(s{2});
        sline = fgetl(fid);
        s = regexp(sline, '=', 'split');
        thresh_egap = str2double(s{2});        
        break
    end
end
fclose(fid);
% -------------------------------------------------------------------
% parsing summary data
% -------------------------------------------------------------------
fname = sprintf('%s_sum.out', prefix);
fid = fopen(fname);
while ~feof(fid)   
    sline = fgets(fid);
    if strfind(sline, 'Summary Table 1')
        header = fgets(fid);
        tab1 = textscan(fid, '%d%f%f%f%f%f%f');
    end
    
    if strfind(sline, 'Summary Table 2')
        header = fgets(fid);
        tab2 = textscan(fid, '%d%f%f%f%d%f%f%f%f%s%s%s%s');
    end    
end
fclose(fid);
% -------------------------------------------------------------------
% parsing thresholds
% -------------------------------------------------------------------
k   = tab1{1};
sg  = tab1{2};
gm  = tab1{3};
L  = tab1{4};
P  = tab1{5};
ap  = tab1{6};
bt  = tab1{7};
pk = tab2{2};
dr = tab2{3};
sl = tab2{4};
nl = tab2{5};
para = tab2{6};
perp = tab2{7};
step = tab2{8};
egap = tab2{9};
lw = sg - gm/2;
up = sg + gm/2;
kmax = max(k);

markersize = 5; fontsize = 12;
xlim1 = [0, max(k) + 1];

figure('color','w', 'pos', [100, 100, 420, 420*2]); 
subplot(511)
clear h
h(1) = plot(k, lw, 'ko-'); hold on; grid on;
h(2) = plot(k, up, 'ko-');
h(3) = plot(k, sg, 'bo-');
h(4) = plot(k, L, 'r.-');
set(h, 'MarkerFaceColor', 'w', 'MarkerSize', markersize);
set(gca, 'FontSize', fontsize);
set(gca, 'xlim', xlim1)
if exist('ylim1', 'var'); 
    set(gca, 'ylim', ylim1); 
else
    set(gca, 'ylim', [min(lw), max(up)])
end
title('E_{low}, E_{upp}, \color{blue}\Sigma, and \color{red}L')

subplot(512)
clear h
h(1) = plot(k, gm*au2ev, 'bo-'); hold on; grid on;
plot([0,kmax+1], thresh_egap*ones(size([0,kmax+1]))*au2ev, 'k--');
set(h, 'MarkerFaceColor', 'w', 'MarkerSize', markersize);
set(gca, 'FontSize', fontsize);
set(gca, 'xlim', xlim1)
if exist('ylim2', 'var'); set(gca, 'ylim', ylim2); end
ylabel('[eV]')
title('\Gamma')

subplot(513)
clear h
h(1) = plot(k, abs(para), 'bo-'); hold on; grid on;
h(2) = plot(k, perp, 'ro-'); hold on; grid on;
plot([0,kmax+1], thresh_grad*ones(size([0,kmax+1])), 'k--');
set(h, 'MarkerFaceColor', 'w', 'MarkerSize', markersize);
set(gca, 'FontSize', fontsize);
set(gca, 'xlim', xlim1)
if exist('ylim3', 'var'); set(gca, 'ylim', ylim3); end
title('|grad,para| and grad,perp')
legend('|para|','perp')

subplot(514)
clear h
h(1) = plot(k(2:end), abs(step(2:end)*au2ev), 'bo-'); hold on; grid on;
plot([0,kmax+1], thresh_step*ones(size([0,kmax+1]))*au2ev, 'k--');
set(h, 'MarkerFaceColor', 'w', 'MarkerSize', markersize);
set(gca, 'FontSize', fontsize);
set(gca, 'xlim', xlim1)
set(gca, 'ylim', [0, 20*thresh_step*au2ev])
ylabel('[eV]')
title('step = | L(k) - L(k-1) |')

subplot(515)
clear h
h(1) = plot(k(1:end-1), nl(1:end-1), 'bo-'); hold on; grid on;
set(h, 'MarkerFaceColor', 'w', 'MarkerSize', markersize);
set(gca, 'FontSize', fontsize);
set(gca, 'xlim', xlim1)
title('#LineSearch')






