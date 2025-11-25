%% photodiode

data = readtable("photodiode_test.csv");

figure
hold on
plot(data.time_s,data.PD1)
plot(data.time_s,data.PD2)
plot(data.time_s,data.PD3)
legend('PD1','PD2','PD3')

%%
LED_dat = readtable("LED_RLI.csv"); %wavelength vs RLI (%)
PD_dat = readtable("OPT101_PR.csv"); %wavelendth vs Photodiode Responsivity (A/W)
CIE_dat = readtable("CIE_sle_photopic.csv"); %https://doi.org/10.25039/CIE.DS.dktna2s3

%interpolate onto same range
dl = 1;
lambd = (200:dl:1100)';
LED_dat = interp1(LED_dat.Var1,LED_dat.Var2./sum(LED_dat.Var2),lambd,'linear',1e-12);
PD_dat = interp1(PD_dat.Var1,PD_dat.Var2,lambd,'linear',1e-12);
CIE_dat = interp1(CIE_dat.Var1,CIE_dat.Var2,lambd,'linear',1e-12);

LED_intensity = 38000; %mcd
viewing_angle = 15; %deg

LED_lm = LED_intensity * (2*pi*(1-cosd(viewing_angle/2)));

srf_wm = (LED_lm * LED_dat)./(683 * CIE_dat);

R_eff = sum(srf_wm .* PD_dat) * dl / (sum(srf_wm) * dl);

figure
plot(data.time_s,(data.PD1/(1e6*R_eff)))
hold on
plot(data.time_s,(data.PD2/(1e6*R_eff)))
plot(data.time_s,(data.PD3/(1e6*R_eff)))
legend('PD1','PD2','PD3')

%% electrometer
data = readtable("electrometer_test.csv",'NumHeaderLines',1);
firstRow = readcell('electrometer_test.csv', 'Range', '1:1');
Fs = sscanf(firstRow{1,1}, 'Sample rate: %f S/s');

R_fb = 1.01e9;
eps_0 = 8.8541878188e-12;
r_e = 6e-3;


BE1 = movmean(data.BE1,floor(Fs/10));
BE2 = movmean(data.BE2,floor(Fs/10));
BE3 = movmean(data.BE3,floor(Fs/10));
BE4 = movmean(data.BE4,floor(Fs/10));

BE1_ER = - BE1 / (4 * pi * r_e^2 * R_fb * eps_0);
BE2_ER = - BE2 / (4 * pi * r_e^2 * R_fb * eps_0);
BE3_ER = - BE3 / (4 * pi * r_e^2 * R_fb * eps_0);
BE4_ER = - BE4 / (4 * pi * r_e^2 * R_fb * eps_0);

figure
tiledlayout(4,1)
nexttile
yyaxis left
plot(data.time_s,data.BE4,'Color',[0.30,0.75,0.93])
hold on
plot(data.time_s,BE4,'-')
yyaxis right
plot(data.time_s,BE4_ER)
%legend('BE4 Voltage','0.1 s Moving Average','$\frac{dE}{dt}$','Location','east')
xticklabels('')

nexttile
yyaxis left
plot(data.time_s,data.BE3,'Color',[0.30,0.75,0.93])
hold on
plot(data.time_s,BE3,'-')
yyaxis right
plot(data.time_s,BE3_ER)
xticklabels('')

nexttile
yyaxis left
plot(data.time_s,data.BE2,'Color',[0.30,0.75,0.93])
hold on
plot(data.time_s,BE2,'-')
yyaxis right
plot(data.time_s,BE2_ER)
xticklabels('')

nexttile
yyaxis left
plot(data.time_s,data.BE1,'Color',[0.30,0.75,0.93])
hold on
plot(data.time_s,BE1,'-')
ylabel('Voltage [V]')
yyaxis right
plot(data.time_s,BE1_ER)
xlabel('Time [s]')
ylabel('$\frac{dE}{dt}$ [V/m/s]')

a4 = decimate(BE4_ER,floor(Fs/10));
a3 = decimate(BE3_ER,floor(Fs/10));
a2 = decimate(BE2_ER,floor(Fs/10));
a1 = decimate(BE1_ER,floor(Fs/10));

t = decimate(data.time_s,floor(Fs/10));

E4 = 0.5 * (a4(1:end-1)+a4(2:end)) * (t(2)-t(1));
E4 = [NaN;E4];
E3 = 0.5 * (a3(1:end-1)+a3(2:end)) * (t(2)-t(1));
E3 = [NaN;E3];
E2 = 0.5 * (a2(1:end-1)+a2(2:end)) * (t(2)-t(1));
E2 = [NaN;E2];
E1 = 0.5 * (a1(1:end-1)+a1(2:end)) * (t(2)-t(1));
E1 = [NaN;E1];
figure
tiledlayout(4,1)
nexttile(1)
plot(t,E4)
nexttile(2)
plot(t,E3)
nexttile(3)
plot(t,E2)
nexttile(4)
plot(t,E1)
%% windspeed
data = readtable("anemometer_zeros.csv",'NumHeaderLines',1);

TIN_zw = mean(data.TInW);
Vert_zw =  mean(data.VertW);

data = readtable("anemometer_zeros2.csv",'NumHeaderLines',1);

TIN_zw(2) = mean(data.TInW);
Vert_zw(2) =  mean(data.VertW);

TIN_zw = mean(TIN_zw);
Vert_zw = mean(Vert_zw);

data = readtable("anemometer_test.csv",'NumHeaderLines',1);
firstRow = readcell('anemometer_test.csv', 'Range', '1:1');
Fs = sscanf(firstRow{1,1}, 'Sample rate: %f S/s');

TInW = movmean(data.TInW,floor(Fs/2));
TInT = movmean(data.TInT,floor(Fs/2));
VertW = movmean(data.VertW,floor(Fs/2));
VertT = movmean(data.VertT,floor(Fs/2));

temp = convert_to_temperature(TInT);
winds = convert_to_winsdpeed(TInW,TIN_zw,temp,'kph');
figure(1)
plot(data.time_s, winds)
hold on

temp = convert_to_temperature(VertT);
winds = convert_to_winsdpeed(VertW,Vert_zw,temp,'kph');
figure(2)
plot(data.time_s, winds)
hold on

data = readtable("anemometer_test2.csv",'NumHeaderLines',1);
firstRow = readcell('anemometer_test2.csv', 'Range', '1:1');
Fs = sscanf(firstRow{1,1}, 'Sample rate: %f S/s');

TInW = movmean(data.TInW,floor(Fs/2));
TInT = movmean(data.TInT,floor(Fs/2));
VertW = movmean(data.VertW,floor(Fs/2));
VertT = movmean(data.VertT,floor(Fs/2));

temp = convert_to_temperature(TInT);
winds = convert_to_winsdpeed(TInW,TIN_zw,temp,'kph');
figure(1)
plot(data.time_s, winds)

temp = convert_to_temperature(VertT);
winds = convert_to_winsdpeed(VertW,Vert_zw,temp,'kph');
figure(2)
plot(data.time_s, winds)


%%
data = readtable("psd_from_image.csv");
histogram(data.size_um,200,'Normalization','pdf')
xlabel("Particle Size [$\mu$m]")
ylabel("Probability Density [-]")
set(gca,"XScale","log")
pd = fitdist(data{:,1}, 'Lognormal');

x_values = linspace(0, 5000, 1000);
y_values = pdf(pd, x_values);

hold on
plot(x_values, y_values, 'r-', 'LineWidth', 2);


%%
function out = convert_to_temperature(input)
 out = (input - 0.400)./(0.0195);
end

function out = convert_to_winsdpeed(input_w, cal_f, temp,unit)
 out = (((input_w - cal_f)./(3.038517 .* (temp.^(0.115157))))./(0.087288)).^(3.009364);
 if strcmp(unit,'mph')
    out(out<0) = 0;
    out = out;
 elseif strcmp(unit,'ms')
    out(out<0) = 0;
    out = out * 0.44704; %mph to m/s
 elseif strcmp(unit,'kph')
    out(out<0) = 0;
    out = out * 1.609344; %mph to kph
 else
     out = out;
 end
end