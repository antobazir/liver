%test_init liver
clear all
load('liver_snap','-ascii');
x = find(liver_snap(:,1)==1);
snap1 = liver_snap(x(1):x(2),:);
snap2 = liver_snap(x(2)+1:x(3),:);
snap3 = liver_snap(x(3)+1:x(4),:);
snap4 = liver_snap(x(4)+1:x(5),:);
snap5 = liver_snap(x(5)+1:x(6),:);
snap6 = liver_snap(x(6)+1:x(7),:);
snap7 = liver_snap(x(7)+1:x(8),:);
snap8 = liver_snap(x(8)+1:x(9),:);
snap9 = liver_snap(x(9)+1:x(10),:);


%figure
%scatter3(liver_snap(:,13),liver_snap(:,14),liver_snap(:,15),18,liver_snap(:,17),'filled')

%figure
%plot3(liver_snap(:,13),liver_snap(:,14),liver_snap(:,15),18,liver_snap(:,17))

figure
scatter3(snap1(:,13),snap1(:,14),snap1(:,15),18,snap1(:,17),'filled')

figure
scatter3(snap2(:,13),snap2(:,14),snap2(:,15),18,snap2(:,17),'filled')

figure
scatter3(snap3(:,13),snap3(:,14),snap3(:,15),18,snap3(:,17),'filled')


figure
scatter3(snap3(:,13),snap3(:,14),snap3(:,15),18,snap3(:,18),'filled')
colorbar

figure
scatter3(snap3(:,13),snap3(:,14),snap3(:,15),18,snap3(:,19),'filled')
colorbar
