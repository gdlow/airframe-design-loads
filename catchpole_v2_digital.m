clear all
close all
clc

% R_t = 0.75
fid = fopen('R_t_0.75.txt');
vals_1 = fscanf(fid, '%f %f ',[2 Inf]);
catchpole_1 = transpose(vals_1);
fclose(fid);

% R_t = 1.00
fid = fopen('R_t_1.00.txt');
vals_2 = fscanf(fid, '%f %f ',[2 Inf]);
catchpole_2 = transpose(vals_2);
fclose(fid);

% R_t = 1.25
fid = fopen('R_t_1.25.txt');
vals_3 = fscanf(fid, '%f %f ',[2 Inf]);
catchpole_3 = transpose(vals_3);
fclose(fid);

% R_t = 1.50
fid = fopen('R_t_1.50.txt');
vals_4 = fscanf(fid, '%f %f ',[2 Inf]);
catchpole_4 = transpose(vals_4);
fclose(fid);

% R_t = 1.75
fid = fopen('R_t_1.75.txt');
vals_5 = fscanf(fid, '%f %f ',[2 Inf]);
catchpole_5 = transpose(vals_5);
fclose(fid);

% R_t = 1.875
fid = fopen('R_t_1.875.txt');
vals_6 = fscanf(fid, '%f %f ',[2 Inf]);
catchpole_6 = transpose(vals_6);
fclose(fid);

% R_t = 2.00
fid = fopen('R_t_2.00.txt');
vals_7 = fscanf(fid, '%f %f ',[2 Inf]);
catchpole_7 = transpose(vals_7);
fclose(fid);

% R_t = 2.25
fid = fopen('R_t_2.25.txt');
vals_8 = fscanf(fid, '%f %f ',[2 Inf]);
catchpole_8 = transpose(vals_8);
fclose(fid);

% R_t = 2.50
fid = fopen('R_t_2.50.txt');
vals_9 = fscanf(fid, '%f %f ',[2 Inf]);
catchpole_9 = transpose(vals_9);
fclose(fid);

% R_t = 3.00
fid = fopen('R_t_3.00.txt');
vals_10 = fscanf(fid, '%f %f ',[2 Inf]);
catchpole_10 = transpose(vals_10);
fclose(fid);

figure
hold on
plot(catchpole_1(:,1),catchpole_1(:,2),'r-')
plot(catchpole_2(:,1),catchpole_2(:,2),'b-')
plot(catchpole_3(:,1),catchpole_3(:,2),'k-')
plot(catchpole_4(:,1),catchpole_4(:,2),'r-.')
plot(catchpole_5(:,1),catchpole_5(:,2),'b-.')
plot(catchpole_6(:,1),catchpole_6(:,2),'k-.')
plot(catchpole_7(:,1),catchpole_7(:,2),'r--')
plot(catchpole_8(:,1),catchpole_8(:,2),'b--')
plot(catchpole_9(:,1),catchpole_9(:,2),'k--')
plot(catchpole_10(:,1),catchpole_10(:,2),'r:')
