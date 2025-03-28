% Automate the plotting process for comparison of results.
% Currently only two different results are compared.

% X. B. 2025-02-22 Initial work

%% Load results
file = dir(fullfile('NucleocytoplasmicTransport*.mat'));
S1 = load(file(1).name);
t1 = deal(S1.t); z1 = deal(S1.z); 
NC1 = deal(S1.NCcell); Vcy1 = deal(S1.Vcy); Vn1 = deal(S1.Vn);
phintf1 = deal(S1.phintf); phiran1 = deal(S1.phiran); p1 = deal(S1.p);
S2 = load(file(2).name);
t2 = deal(S2.t); z2 = deal(S2.z); 
NC2 = deal(S2.NCcell); Vcy2 = deal(S2.Vcy); Vn2 = deal(S2.Vn);
phintf2 = deal(S2.phintf); phiran2 = deal(S2.phiran); p2 = deal(S2.p);
clear S1 S2
% Total number of Ran
zran1 = z1(:,2)+z1(:,3)+z1(:,6)+z1(:,9)+z1(:,10)+z1(:,13)+z1(:,17);
zran2 = z2(:,2)+z2(:,3)+z2(:,6)+z2(:,9)+z2(:,10)+z2(:,13)+z2(:,17);
% Number of Ran in the nucleus
zrann1 = z1(:,9)+z1(:,10)+z1(:,13)+z1(:,17);
zrann2 = z2(:,9)+z2(:,10)+z2(:,13)+z2(:,17);
% Fraction of Ran in the nucleus
frann1 = zrann1./zran1;
frann2 = zrann2./zran2;
% Total number of cargo proteins
zn1 = z1(:,4)+z1(:,5)+z1(:,11)+z1(:,12);
zn2 = z2(:,4)+z2(:,5)+z2(:,11)+z2(:,12);
% Total number of cargo proteins in the nucleus
znn1 = z1(:,11)+z1(:,12);
znn2 = z2(:,11)+z2(:,12);
% Fraction of cargo proteins in the nucleus
fnn1 = znn1./zn1;
fnn2 = znn2./zn2;
% Total number of NTF2
zntf1 = z1(:,1)+z1(:,3)+z1(:,8)+z1(:,10);
zntf2 = z2(:,1)+z2(:,3)+z2(:,8)+z2(:,10);
% Total number of NTF2 in the nucleus
zntfn1 = z1(:,8)+z1(:,10);
zntfn2 = z2(:,8)+z2(:,10);
% Fraction of NTF2 in the nucleus
fntfn1 = zntfn1./zntf1;
fntfn2 = zntfn2./zntf2;
% Total number of Importin
zt1 = z1(:,4)+z1(:,6)+z1(:,7)+z1(:,11)+z1(:,13)+z1(:,14);
zt2 = z2(:,4)+z2(:,6)+z2(:,7)+z2(:,11)+z2(:,13)+z2(:,14);
% Total number of Importin in the nucleus
ztn1 = z1(:,11)+z1(:,13)+z1(:,14);
ztn2 = z2(:,11)+z2(:,13)+z2(:,14);
% Fraction of Importin in the nucleus
ftn1 = ztn1./zt1;
ftn2 = ztn2./zt2;
% Total number of Nucleocytoplasmic Machinery proteins
zm1 = sum(z1(:,1:14),2)+z1(:,17);
zm2 = sum(z2(:,1:14),2)+z2(:,17);
% Total number of Nucleocytoplasmic Machinery proteins in the nucleus
zmn1 = sum(z1(:,8:14),2)+z1(:,17);
zmn2 = sum(z2(:,8:14),2)+z2(:,17);
% Fraction of Nucleocytoplasmic Machinery proteins in the nucleus
fmn1 = zmn1./zm1;
fmn2 = zmn2./zm2;

%% Plots
% %{
figure(1)
plot(t1,p1.Vref*z1(:,1)./Vcy1,'Linewidth',2)
hold on
plot(t2,p2.Vref*z2(:,1)./Vcy2,'Linewidth',2)
hold off
legend(['\phi_{NTF}=' num2str(phintf1,'%.2f') ', \phi_{Ran}=' num2str(phiran1,'%.3f') ', C_{Ran}=' num2str(p1.Cran,'%.2f')],...
['\phi_{NTF}=' num2str(phintf2,'%.2f') ', \phi_{Ran}=' num2str(phiran2,'%.3f') ', C_{Ran}=' num2str(p2.Cran,'%.2f')])
xlabel('Time')
ylabel('Cytoplasmic NTF2 (nM)')
figure(2)
plot(t1,p1.Vref*z1(:,2)./Vcy1,'Linewidth',2)
hold on
plot(t2,p2.Vref*z2(:,2)./Vcy2,'Linewidth',2)
hold off
legend(['\phi_{NTF}=' num2str(phintf1,'%.2f') ', \phi_{Ran}=' num2str(phiran1,'%.3f') ', C_{Ran}=' num2str(p1.Cran,'%.2f')],...
['\phi_{NTF}=' num2str(phintf2,'%.2f') ', \phi_{Ran}=' num2str(phiran2,'%.3f') ', C_{Ran}=' num2str(p2.Cran,'%.2f')])
xlabel('Time')
ylabel('Cytoplasmic RanGDP (nM)')
figure(3)
plot(t1,p1.Vref*z1(:,3)./Vcy1,'Linewidth',2)
hold on
plot(t2,p2.Vref*z2(:,3)./Vcy2,'Linewidth',2)
hold off
legend(['\phi_{NTF}=' num2str(phintf1,'%.2f') ', \phi_{Ran}=' num2str(phiran1,'%.3f') ', C_{Ran}=' num2str(p1.Cran,'%.2f')],...
['\phi_{NTF}=' num2str(phintf2,'%.2f') ', \phi_{Ran}=' num2str(phiran2,'%.3f') ', C_{Ran}=' num2str(p2.Cran,'%.2f')])
xlabel('Time')
ylabel('Cytoplasmic NTF2-RanGDP complex (nM)')
figure(4)
plot(t1,p1.Vref*z1(:,4)./Vcy1,'Linewidth',2)
hold on
plot(t2,p2.Vref*z2(:,4)./Vcy2,'Linewidth',2)
hold off
legend(['\phi_{NTF}=' num2str(phintf1,'%.2f') ', \phi_{Ran}=' num2str(phiran1,'%.3f') ', C_{Ran}=' num2str(p1.Cran,'%.2f')],...
['\phi_{NTF}=' num2str(phintf2,'%.2f') ', \phi_{Ran}=' num2str(phiran2,'%.3f') ', C_{Ran}=' num2str(p2.Cran,'%.2f')])
xlabel('Time')
ylabel('Cytoplasmic importin-cargo complex (nM)')
figure(5)
plot(t1,p1.Vref*z1(:,5)./Vcy1,'Linewidth',2)
hold on
plot(t2,p2.Vref*z2(:,5)./Vcy2,'Linewidth',2)
hold off
legend(['\phi_{NTF}=' num2str(phintf1,'%.2f') ', \phi_{Ran}=' num2str(phiran1,'%.3f') ', C_{Ran}=' num2str(p1.Cran,'%.2f')],...
['\phi_{NTF}=' num2str(phintf2,'%.2f') ', \phi_{Ran}=' num2str(phiran2,'%.3f') ', C_{Ran}=' num2str(p2.Cran,'%.2f')])
xlabel('Time')
ylabel('Cytoplasmic cargo (nM)')
figure(6)
plot(t1,p1.Vref*z1(:,6)./Vcy1,'Linewidth',2)
hold on
plot(t2,p2.Vref*z2(:,6)./Vcy2,'Linewidth',2)
hold off
legend(['\phi_{NTF}=' num2str(phintf1,'%.2f') ', \phi_{Ran}=' num2str(phiran1,'%.3f') ', C_{Ran}=' num2str(p1.Cran,'%.2f')],...
['\phi_{NTF}=' num2str(phintf2,'%.2f') ', \phi_{Ran}=' num2str(phiran2,'%.3f') ', C_{Ran}=' num2str(p2.Cran,'%.2f')])
xlabel('Time')
ylabel('Cytoplasmic importin-RanGTP complex (nM)')
figure(7)
plot(t1,p1.Vref*z1(:,7)./Vcy1,'Linewidth',2)
hold on
plot(t2,p2.Vref*z2(:,7)./Vcy2,'Linewidth',2)
hold off
legend(['\phi_{NTF}=' num2str(phintf1,'%.2f') ', \phi_{Ran}=' num2str(phiran1,'%.3f') ', C_{Ran}=' num2str(p1.Cran,'%.2f')],...
['\phi_{NTF}=' num2str(phintf2,'%.2f') ', \phi_{Ran}=' num2str(phiran2,'%.3f') ', C_{Ran}=' num2str(p2.Cran,'%.2f')])
xlabel('Time')
ylabel('Cytoplasmic importin (nM)')
figure(8)
plot(t1,p1.Vref*z1(:,8)./Vn1,'Linewidth',2)
hold on
plot(t2,p2.Vref*z2(:,8)./Vn2,'Linewidth',2)
hold off
legend(['\phi_{NTF}=' num2str(phintf1,'%.2f') ', \phi_{Ran}=' num2str(phiran1,'%.3f') ', C_{Ran}=' num2str(p1.Cran,'%.2f')],...
['\phi_{NTF}=' num2str(phintf2,'%.2f') ', \phi_{Ran}=' num2str(phiran2,'%.3f') ', C_{Ran}=' num2str(p2.Cran,'%.2f')])
xlabel('Time')
ylabel('Nuclear NTF2 (nM)')
figure(9)
plot(t1,p1.Vref*z1(:,9)./Vn1,'Linewidth',2)
hold on
plot(t2,p2.Vref*z2(:,9)./Vn2,'Linewidth',2)
hold off
legend(['\phi_{NTF}=' num2str(phintf1,'%.2f') ', \phi_{Ran}=' num2str(phiran1,'%.3f') ', C_{Ran}=' num2str(p1.Cran,'%.2f')],...
['\phi_{NTF}=' num2str(phintf2,'%.2f') ', \phi_{Ran}=' num2str(phiran2,'%.3f') ', C_{Ran}=' num2str(p2.Cran,'%.2f')])
xlabel('Time')
ylabel('Nuclear RanGTP (nM)')
figure(10)
plot(t1,p1.Vref*z1(:,10)./Vn1,'Linewidth',2)
hold on
plot(t2,p2.Vref*z2(:,10)./Vn2,'Linewidth',2)
hold off
legend(['\phi_{NTF}=' num2str(phintf1,'%.2f') ', \phi_{Ran}=' num2str(phiran1,'%.3f') ', C_{Ran}=' num2str(p1.Cran,'%.2f')],...
['\phi_{NTF}=' num2str(phintf2,'%.2f') ', \phi_{Ran}=' num2str(phiran2,'%.3f') ', C_{Ran}=' num2str(p2.Cran,'%.2f')])
xlabel('Time')
ylabel('Nuclear NTF2-RanGDP complex (nM)')
figure(11)
plot(t1,p1.Vref*z1(:,11)./Vn1,'Linewidth',2)
hold on
plot(t2,p2.Vref*z2(:,11)./Vn2,'Linewidth',2)
hold off
legend(['\phi_{NTF}=' num2str(phintf1,'%.2f') ', \phi_{Ran}=' num2str(phiran1,'%.3f') ', C_{Ran}=' num2str(p1.Cran,'%.2f')],...
['\phi_{NTF}=' num2str(phintf2,'%.2f') ', \phi_{Ran}=' num2str(phiran2,'%.3f') ', C_{Ran}=' num2str(p2.Cran,'%.2f')])
xlabel('Time')
ylabel('Nuclear importin-cargo complex (nM)')
figure(12)
plot(t1,p1.Vref*z1(:,12)./Vn1,'Linewidth',2)
hold on
plot(t2,p2.Vref*z2(:,12)./Vn2,'Linewidth',2)
hold off
legend(['\phi_{NTF}=' num2str(phintf1,'%.2f') ', \phi_{Ran}=' num2str(phiran1,'%.3f') ', C_{Ran}=' num2str(p1.Cran,'%.2f')],...
['\phi_{NTF}=' num2str(phintf2,'%.2f') ', \phi_{Ran}=' num2str(phiran2,'%.3f') ', C_{Ran}=' num2str(p2.Cran,'%.2f')])
xlabel('Time')
ylabel('Nuclear cargo (nM)')
figure(13)
plot(t1,p1.Vref*z1(:,13)./Vn1,'Linewidth',2)
hold on
plot(t2,p2.Vref*z2(:,13)./Vn2,'Linewidth',2)
hold off
legend(['\phi_{NTF}=' num2str(phintf1,'%.2f') ', \phi_{Ran}=' num2str(phiran1,'%.3f') ', C_{Ran}=' num2str(p1.Cran,'%.2f')],...
['\phi_{NTF}=' num2str(phintf2,'%.2f') ', \phi_{Ran}=' num2str(phiran2,'%.3f') ', C_{Ran}=' num2str(p2.Cran,'%.2f')])
xlabel('Time')
ylabel('Nuclear importin-RanGTP complex (nM)')
figure(14)
plot(t1,p1.Vref*z1(:,14)./Vn1,'Linewidth',2)
hold on
plot(t2,p2.Vref*z2(:,14)./Vn2,'Linewidth',2)
hold off
legend(['\phi_{NTF}=' num2str(phintf1,'%.2f') ', \phi_{Ran}=' num2str(phiran1,'%.3f') ', C_{Ran}=' num2str(p1.Cran,'%.2f')],...
['\phi_{NTF}=' num2str(phintf2,'%.2f') ', \phi_{Ran}=' num2str(phiran2,'%.3f') ', C_{Ran}=' num2str(p2.Cran,'%.2f')])
xlabel('Time')
ylabel('Nuclear importin (nM)')
figure(15)
plot(t1,p1.Vref*z1(:,17)./Vn1,'Linewidth',2)
hold on
plot(t2,p2.Vref*z2(:,17)./Vn2,'Linewidth',2)
hold off
legend(['\phi_{NTF}=' num2str(phintf1,'%.2f') ', \phi_{Ran}=' num2str(phiran1,'%.3f') ', C_{Ran}=' num2str(p1.Cran,'%.2f')],...
['\phi_{NTF}=' num2str(phintf2,'%.2f') ', \phi_{Ran}=' num2str(phiran2,'%.3f') ', C_{Ran}=' num2str(p2.Cran,'%.2f')])
xlabel('Time')
ylabel('Nuclear RanGDP (nM)')
% %}
figure(16)
plot(t1,NC1,'Linewidth',2)
hold on
plot(t2,NC2,'Linewidth',2)
hold off
legend(['\phi_{NTF}=' num2str(phintf1,'%.2f') ', \phi_{Ran}=' num2str(phiran1,'%.3f') ', C_{Ran}=' num2str(p1.Cran,'%.2f')],...
['\phi_{NTF}=' num2str(phintf2,'%.2f') ', \phi_{Ran}=' num2str(phiran2,'%.3f') ', C_{Ran}=' num2str(p2.Cran,'%.2f')])
xlabel('Time')
ylabel('Nuclear-to-cell ratio')
% %{
figure(17)
plot(t1,frann1,'Linewidth',2)
hold on
plot(t2,frann2,'Linewidth',2)
hold off
legend(['\phi_{NTF}=' num2str(phintf1,'%.2f') ', \phi_{Ran}=' num2str(phiran1,'%.3f') ', C_{Ran}=' num2str(p1.Cran,'%.2f')],...
['\phi_{NTF}=' num2str(phintf2,'%.2f') ', \phi_{Ran}=' num2str(phiran2,'%.3f') ', C_{Ran}=' num2str(p2.Cran,'%.2f')])
xlabel('Time')
ylabel('Fraction of Ran in the nucleus')
figure(18)
plot(t1,fnn1,'Linewidth',2)
hold on
plot(t2,fnn2,'Linewidth',2)
hold off
legend(['\phi_{NTF}=' num2str(phintf1,'%.2f') ', \phi_{Ran}=' num2str(phiran1,'%.3f') ', C_{Ran}=' num2str(p1.Cran,'%.2f')],...
['\phi_{NTF}=' num2str(phintf2,'%.2f') ', \phi_{Ran}=' num2str(phiran2,'%.3f') ', C_{Ran}=' num2str(p2.Cran,'%.2f')])
xlabel('Time')
ylabel('Fraction of Cargo in the nucleus')
figure(19)
plot(t1,fntfn1,'Linewidth',2)
hold on
plot(t2,fntfn2,'Linewidth',2)
hold off
legend(['\phi_{NTF}=' num2str(phintf1,'%.2f') ', \phi_{Ran}=' num2str(phiran1,'%.3f') ', C_{Ran}=' num2str(p1.Cran,'%.2f')],...
['\phi_{NTF}=' num2str(phintf2,'%.2f') ', \phi_{Ran}=' num2str(phiran2,'%.3f') ', C_{Ran}=' num2str(p2.Cran,'%.2f')])
xlabel('Time')
ylabel('Fraction of NTF2 in the nucleus')
figure(20)
plot(t1,ftn1,'Linewidth',2)
hold on
plot(t2,ftn2,'Linewidth',2)
hold off
legend(['\phi_{NTF}=' num2str(phintf1,'%.2f') ', \phi_{Ran}=' num2str(phiran1,'%.3f') ', C_{Ran}=' num2str(p1.Cran,'%.2f')],...
['\phi_{NTF}=' num2str(phintf2,'%.2f') ', \phi_{Ran}=' num2str(phiran2,'%.3f') ', C_{Ran}=' num2str(p2.Cran,'%.2f')])
xlabel('Time')
ylabel('Fraction of Importin in the nucleus')
figure(21)
plot(t1,fmn1,'Linewidth',2)
hold on
plot(t2,fmn2,'Linewidth',2)
hold off
legend(['\phi_{NTF}=' num2str(phintf1,'%.2f') ', \phi_{Ran}=' num2str(phiran1,'%.3f') ', C_{Ran}=' num2str(p1.Cran,'%.2f')],...
['\phi_{NTF}=' num2str(phintf2,'%.2f') ', \phi_{Ran}=' num2str(phiran2,'%.3f') ', C_{Ran}=' num2str(p2.Cran,'%.2f')])
xlabel('Time')
ylabel('Fraction of NCT Machinery in the nucleus')
figure(22)
plot(t1,p1.Vref*zrann1./Vn1,'Linewidth',2)
hold on
plot(t2,p2.Vref*zrann2./Vn2,'Linewidth',2)
hold off
legend(['\phi_{NTF}=' num2str(phintf1,'%.2f') ', \phi_{Ran}=' num2str(phiran1,'%.3f') ', C_{Ran}=' num2str(p1.Cran,'%.2f')],...
['\phi_{NTF}=' num2str(phintf2,'%.2f') ', \phi_{Ran}=' num2str(phiran2,'%.3f') ', C_{Ran}=' num2str(p2.Cran,'%.2f')])
xlabel('Time')
ylabel('Total nuclear Ran (nM)')
% %}
