data = readtable("photodiode_test.csv");

figure
hold on
plot(data.time_s,data.PD1)
plot(data.time_s,data.PD2)
plot(data.time_s,data.PD3)
legend('PD1','PD2','PD3')

%%
data = readtable("electrometer_test.csv",'NumHeaderLines',1);
firstRow = readcell('electrometer_test.csv', 'Range', '1:1');
Fs = sscanf(firstRow{1,1}, 'Sample rate: %f S/s');

BE1 = movmean(data.BE1,floor(Fs/10));
BE2 = movmean(data.BE2,floor(Fs/10));
BE3 = movmean(data.BE3,floor(Fs/10));
BE4 = movmean(data.BE4,floor(Fs/10));

figure
tiledlayout(4,1)
nexttile
plot(data.time_s,data.BE4)
hold on
%plot(data.time_s,smooth(data.BE4,floor(Fs)/4,'lowess'))
plot(data.time_s,BE4)
legend('BE4','10 s Moving Average','Location','east')
xticklabels('')

nexttile
plot(data.time_s,data.BE3)
hold on
%plot(data.time_s,smooth(data.BE4,floor(Fs)/4,'lowess'))
plot(data.time_s,BE3)
legend('BE3','10 s Moving Average','Location','east')
xticklabels('')

nexttile
plot(data.time_s,data.BE2)
hold on
%plot(data.time_s,smooth(data.BE4,floor(Fs)/4,'lowess'))
plot(data.time_s,BE2)
legend('BE2','10 s Moving Average','Location','west')
xticklabels('')

nexttile
plot(data.time_s,data.BE1)
hold on
%plot(data.time_s,smooth(data.BE4,floor(Fs)/4,'lowess'))
plot(data.time_s,BE1)
legend('BE1','10 s Moving Average','Location','west')
xlabel('Time [s]')

%%
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
histogram(data.size_um,500)

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