n=128;
mid=n/2+1;
h=1/n;
gd=h*(0:n-1);

load('test1.mat');
ugd=reshape(real(u),[n,n,n]);
range_u=min(abs(max(ugd(:))),abs(min(ugd(:))));
vel=vel(:,:,mid);
ugd=ugd(:,:,mid);

figure; imagesc(gd,gd,vel); colorbar; axis equal; axis tight; colormap jet;
set(gca, 'FontSize', 16);
bb=get(gca);
set(bb.XLabel, 'FontSize', 16);
set(bb.YLabel, 'FontSize', 16);
set(bb.ZLabel, 'FontSize', 16);
set(bb.Title, 'FontSize', 16);
print(gcf, '-depsc', 'H31c');

figure; imagesc(gd,gd,ugd); colorbar; axis equal; axis tight; colormap jet;
caxis([-range_u,range_u]);
set(gca, 'FontSize', 16);
bb=get(gca);
set(bb.XLabel, 'FontSize', 16);
set(bb.YLabel, 'FontSize', 16);
set(bb.ZLabel, 'FontSize', 16);
set(bb.Title, 'FontSize', 16);
print(gcf, '-depsc', 'H31u');

load('test2.mat');
ugd=reshape(real(u),[n,n,n]);
range_u=min(abs(max(ugd(:))),abs(min(ugd(:))));
vel=vel(:,:,mid);
ugd=ugd(:,:,mid);

figure; imagesc(gd,gd,vel); colorbar; axis equal; axis tight; colormap jet;
set(gca, 'FontSize', 16);
bb=get(gca);
set(bb.XLabel, 'FontSize', 16);
set(bb.YLabel, 'FontSize', 16);
set(bb.ZLabel, 'FontSize', 16);
set(bb.Title, 'FontSize', 16);
print(gcf, '-depsc', 'H32c');

figure; imagesc(gd,gd,ugd); colorbar; axis equal; axis tight; colormap jet;
caxis([-range_u,range_u]);
set(gca, 'FontSize', 16);
bb=get(gca);
set(bb.XLabel, 'FontSize', 16);
set(bb.YLabel, 'FontSize', 16);
set(bb.ZLabel, 'FontSize', 16);
set(bb.Title, 'FontSize', 16);
print(gcf, '-depsc', 'H32u');

load('test3.mat');
ugd=reshape(real(u),[n,n,n]);
range_u=min(abs(max(ugd(:))),abs(min(ugd(:))));
V=V(:,:,mid);
ugd=ugd(:,:,mid);

figure; imagesc(gd,gd,V); colorbar; axis equal; axis tight; colormap jet;
set(gca, 'FontSize', 16);
bb=get(gca);
set(bb.XLabel, 'FontSize', 16);
set(bb.YLabel, 'FontSize', 16);
set(bb.ZLabel, 'FontSize', 16);
set(bb.Title, 'FontSize', 16);
print(gcf, '-depsc', 'S31V');

figure; imagesc(gd,gd,ugd); colorbar; axis equal; axis tight; colormap jet;
caxis([-range_u,range_u]);
set(gca, 'FontSize', 16);
bb=get(gca);
set(bb.XLabel, 'FontSize', 16);
set(bb.YLabel, 'FontSize', 16);
set(bb.ZLabel, 'FontSize', 16);
set(bb.Title, 'FontSize', 16);
print(gcf, '-depsc', 'S31u');

load('test4.mat');
ugd=reshape(real(u),[n,n,n]);
range_u=min(abs(max(ugd(:))),abs(min(ugd(:))));
V=V(:,:,mid);
ugd=ugd(:,:,mid);

figure; imagesc(gd,gd,V); colorbar; axis equal; axis tight; colormap jet;
set(gca, 'FontSize', 16);
bb=get(gca);
set(bb.XLabel, 'FontSize', 16);
set(bb.YLabel, 'FontSize', 16);
set(bb.ZLabel, 'FontSize', 16);
set(bb.Title, 'FontSize', 16);
print(gcf, '-depsc', 'S32V');

figure; imagesc(gd,gd,ugd); colorbar; axis equal; axis tight; colormap jet;
caxis([-range_u,range_u]);
set(gca, 'FontSize', 16);
bb=get(gca);
set(bb.XLabel, 'FontSize', 16);
set(bb.YLabel, 'FontSize', 16);
set(bb.ZLabel, 'FontSize', 16);
set(bb.Title, 'FontSize', 16);
print(gcf, '-depsc', 'S32u');