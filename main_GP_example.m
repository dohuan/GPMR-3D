close all;
clear all;
clc;
%% Generate the field
addpath(genpath('./gpml'))

% --- Config field 1
Psi1 = [4 4 10]; % psi_2x, psi_2y, psi_2f
sig_w1 = 1e-5; % small number to avoid numerical problem
covfunc1 = {'covSum', {'covSEard','covNoise'}};
loghyper1 = [log(Psi1(1));log(Psi1(2));log(Psi1(3));log(sig_w1)];

% --- Config field 2
Psi2 = [7 7 5]; % psi_2x, psi_2y, psi_2f
sig_w2 = 1e-5; % small number to avoid numerical problem
covfunc2 = {'covSum', {'covSEard','covNoise'}};
loghyper2 = [log(Psi2(1));log(Psi2(2));log(Psi2(3));log(sig_w2)];

num_grid = 51;
XMIN = -15; XMAX = 15;
YMIN = -15; YMAX = 15;
x_mesh = linspace(XMIN,XMAX,num_grid)';
y_mesh = linspace(YMIN,YMAX,num_grid)';
[S1,S2] = meshgrid(x_mesh,y_mesh);
s_grid = [S1(:),S2(:)];

Z_grid1 = chol(feval(covfunc1{:},loghyper1,s_grid))'*randn(num_grid^2,1);
Z_field1 = reshape(Z_grid1,num_grid,num_grid);

Z_grid2 = chol(feval(covfunc2{:},loghyper2,s_grid))'*randn(num_grid^2,1);
Z_field2 = reshape(Z_grid2,num_grid,num_grid);


figure('name','field');
subplot(2,2,1)
surf(Z_field1,'edgecolor','none');
title('true field 1')

subplot(2,2,3)
surf(Z_field2,'edgecolor','none');
title('true field 2')


%% Create sample locations of the fields
nt = 40;

index_temp = round(1+(size(s_grid,1)-1).*rand(nt,1));
sample_loc1 = s_grid(index_temp,:);
sample_field1 = Z_grid1(index_temp,:);

index_temp = round(1+(size(s_grid,1)-1).*rand(nt,1));
sample_loc2 = s_grid(index_temp,:);
sample_field2 = Z_grid2(index_temp,:);

%% Run GMPL to estimate the fields

covFunc = @covSEard;
likFunc = @likGauss;

% --- Estimate field 1
hyp1.cov = log([2 2 6]);
hyp1.lik = log(5e-5);

hyp1 = minimize(hyp1,@gp,-50,@infExact,[],covFunc,likFunc,...
                        sample_loc1,sample_field1);

[est1,var1] = gp(hyp1,@infExact,[],covFunc,likFunc,sample_loc1,sample_field1,s_grid);
est_field1 = reshape(est1,num_grid,num_grid);

fprintf('True hyper-parameter field 1: %d %d %d\n',Psi1(1),Psi1(2),Psi1(3));
fprintf('Estimated hyper-parameter field 1: %.1f %.1f %.1f\n',...
    exp(hyp1.cov(1)),exp(hyp1.cov(2)),exp(hyp1.cov(3)));

% --- Estimate field 2
hyp2.cov = log([5 5 7]);
hyp2.lik = log(5e-5);

hyp2 = minimize(hyp2,@gp,-50,@infExact,[],covFunc,likFunc,...
                        sample_loc2,sample_field2);

[est2,var2] = gp(hyp2,@infExact,[],covFunc,likFunc,sample_loc2,sample_field2,s_grid);
est_field2 = reshape(est2,num_grid,num_grid);

fprintf('True hyper-parameter field 2: %d %d %d\n',Psi2(1),Psi2(2),Psi2(3));
fprintf('Estimated hyper-parameter field 2: %.1f %.1f %.1f\n',...
    exp(hyp2.cov(1)),exp(hyp2.cov(2)),exp(hyp2.cov(3)));



subplot(2,2,2)
surf(est_field1,'edgecolor','none');
title('estimated field 1')

subplot(2,2,4)
surf(est_field2,'edgecolor','none');
title('estimated field 2')

figure('name','squared error')
subplot(1,2,1)
title('field 1')
serr_field1 = (est_field1-Z_field1).^2;
imagesc(serr_field1);

subplot(1,2,2)
title('field 2')
serr_field2 = (est_field2-Z_field2).^2;
imagesc(serr_field2);
colorbar

figure('name','sampled positions')
subplot(1,2,1)
title('field 1')
imagesc(x_mesh,y_mesh,est_field1);
hold on
plot(sample_loc1(:,1),sample_loc1(:,2),'kx');

subplot(1,2,2)
title('field 2')
imagesc(x_mesh,y_mesh,est_field2);
hold on
plot(sample_loc2(:,1),sample_loc2(:,2),'kx');



