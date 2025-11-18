data = readtable("sand_data.csv");
nick_data = readtable("nick_iceland_data.csv",NumHeaderLines=1);
data2 = readtable("airport.csv");
data3 = readtable("Ajman.csv");
data4 = readtable("CentralAlAin.csv");
data5 = readtable("MGS1.csv");
data6 = readtable("MGS1c.csv");

bins = [data.SizeClass;data.x__m_(end)];
bins2 = [data2.SizeClass;data2.x__m_(end)];
bins3 = [data3.SizeClass;data3.x_m_(end)];
bins4 = [data4.SizeClass;data4.x__m_(end)];
bins5 = [data5.SizeClass;data5.x__m_(end)];
bins6 = [data6.SizeClass;data6.x__m_(end)];

cum_dist = cumsum(data.p3___);
cum_dist2 = cumsum(data2.p3___);
cum_dist3 = cumsum(data3.p3___);
cum_dist4 = cumsum(data4.p3___);
cum_dist5 = cumsum(data5.p3___);
cum_dist6 = cumsum(data6.p3___);

bin_mid = movmean(bins,2);
bin_mid = bin_mid(2:end);
bin_mid2 = movmean(bins2,2);
bin_mid2 = bin_mid2(2:end);
bin_mid3 = movmean(bins3,2);
bin_mid3 = bin_mid3(2:end);
bin_mid4 = movmean(bins4,2);
bin_mid4 = bin_mid4(2:end);
bin_mid5 = movmean(bins5,2);
bin_mid5 = bin_mid5(2:end);
bin_mid6 = movmean(bins6,2);
bin_mid6 = bin_mid6(2:end);

colours = lines(6);
figure
semilogx(bin_mid,cum_dist./100,'k')
hold on
semilogx(nick_data.X,nick_data.Y','--','Color',colours(1,:))
semilogx(nick_data.X_1,nick_data.Y_1,'Color',colours(1,:))

semilogx(nick_data.X_2,nick_data.Y_2,'--','Color',colours(4,:))
semilogx(nick_data.X_3,nick_data.Y_3,'Color',colours(4,:))

semilogx(bin_mid2,cum_dist2./100,'Color',colours(3,:))
semilogx(bin_mid3,cum_dist3./100,'--','Color',colours(3,:))
semilogx(bin_mid4,cum_dist4./100,':','Color',colours(3,:))

semilogx(bin_mid5,cum_dist5./100,'Color',colours(2,:))
semilogx(bin_mid6,cum_dist6./100,'--','Color',colours(2,:))

legend('Sand','Holasandur (10cm depth)','Holasandur (surface)', ...
    'Hverfjall (10cm depth)','Hverfjall (surface)', 'Airport', ...
    'Ajman','Central Al Ain','MGS-1','MGS-1C')

xlabel('Particle Size [$\mu$m]')
ylabel('Cumulative Volume Fraction [-]')