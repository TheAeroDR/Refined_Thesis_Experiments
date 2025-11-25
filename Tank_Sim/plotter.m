filestring = sprintf('./field_t*.h5');

files = dir(filestring);

% Sort files by name to ensure correct temporal order
[~, idx] = sort({files.name});
files = files(idx);
 
for j = 1:120
    x_all(:,j) = h5read(fullfile(files(j).folder,files(j).name),'/x');
    y_all(:,j) = h5read(fullfile(files(j).folder,files(j).name),'/y');
    z_all(:,j) = h5read(fullfile(files(j).folder,files(j).name),'/z');
    Ex_all(:,j) = h5read(fullfile(files(j).folder,files(j).name),'/Ex');
    Ey_all(:,j) = h5read(fullfile(files(j).folder,files(j).name),'/Ey');
    Ez_all(:,j) = h5read(fullfile(files(j).folder,files(j).name),'/Ez');
    Bx_all(:,j) = h5read(fullfile(files(j).folder,files(j).name),'/Bx');
    By_all(:,j) = h5read(fullfile(files(j).folder,files(j).name),'/By');
    Bz_all(:,j) = h5read(fullfile(files(j).folder,files(j).name),'/Bz');
end

grid_size_x = length(unique(x_all(:,1)));
grid_size_y = length(unique(y_all(:,1)));

%%
x_reshaped = reshape(x_all(:,1), grid_size_x, grid_size_y)';
x_unique = x_reshaped(1,:);
y_reshaped = reshape(y_all(:,1), grid_size_x, grid_size_y)'; 
y_unique = y_reshaped(:,1);


[X,Y] = meshgrid(x_unique, y_unique);
for i = 1:size(Bz_all,2)
    Bx(:,:,i) = reshape(Bx_all(:,i), grid_size_x, grid_size_y)';
    By(:,:,i) = reshape(By_all(:,i), grid_size_x, grid_size_y)';
    Bz(:,:,i) = reshape(Bz_all(:,i), grid_size_x, grid_size_y)';
    Ex(:,:,i) = reshape(Ex_all(:,i), grid_size_x, grid_size_y)';
    Ey(:,:,i) = reshape(Ey_all(:,i), grid_size_x, grid_size_y)';
    Ez(:,:,i) = reshape(Ez_all(:,i), grid_size_x, grid_size_y)';
end

%%
figure(1)

bz_min = -2e-12;
bz_max = 2e-12;
ez_min = -6e6;
ez_max = 6e6;

s = surf(X, Y, Bx(:,:,1), 'EdgeColor', 'none');
colormap jet;
c = colorbar;
xlabel('X [m]'); ylabel('Y [m]'); zlabel('Bz');
axis tight;
axis equal;
xlim([-0.6,0.6])
ylim([-0.6,0.6])
clim([bz_min, bz_max]);
c.Label.Interpreter = "latex";
c.Label.String = "$B_X$ [T]";
view(2);
hold on
radius = 0.505;
theta = 0:0.01:2*pi;
z_height = 3e-12;
x = radius * cos(theta);
y = radius * sin(theta);
z = ones(size(theta)) * z_height;

plot3(x, y, z, 'b', 'LineWidth', 2); 

%%
for i = 1:size(Bz,3)
    set(s,'ZData',Bx(:,:,i))
    drawnow;  % Make sure the figure is updated
    pause(0.05);
end

%%

figure(2)
tiledlayout(1,3)
nexttile(1)
surf(X, Y, Bx(:,:,100), 'EdgeColor', 'none');
colormap jet;
c = colorbar;
xlabel('X [m]'); ylabel('Y [m]'); zlabel('Bz');
axis tight;
axis equal;
xlim([-0.6,0.6])
ylim([-0.6,0.6])
clim([bz_min, bz_max]);
c.Label.Interpreter = "latex";
c.Label.String = "$B_X$ [T]";
view(2);
hold on
plot3(x, y, z, 'b', 'LineWidth', 2); 

nexttile(2)
surf(X, Y, By(:,:,100), 'EdgeColor', 'none');
colormap jet;
c = colorbar;
xlabel('X [m]'); ylabel('Y [m]'); zlabel('Bz');
axis tight;
axis equal;
xlim([-0.6,0.6])
ylim([-0.6,0.6])
clim([bz_min, bz_max]);
c.Label.Interpreter = "latex";
c.Label.String = "$B_Y$ [T]";
view(2);
hold on
plot3(x, y, z, 'b', 'LineWidth', 2); 

nexttile(3)
s = surf(X, Y, Bz(:,:,100), 'EdgeColor', 'none');
colormap jet;
c = colorbar;
xlabel('X [m]'); ylabel('Y [m]'); zlabel('Bz');
axis tight;
axis equal;
xlim([-0.6,0.6])
ylim([-0.6,0.6])
clim([bz_min, bz_max]);
c.Label.Interpreter = "latex";
c.Label.String = "$B_Z$ [T]";
view(2);
hold on
plot3(x, y, z, 'b', 'LineWidth', 2); 