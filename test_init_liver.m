%test_init liver

load('liver_snap','-ascii');

figure
scatter3(liver_snap(:,13),liver_snap(:,14),liver_snap(:,15),18,liver_snap(:,17),'filled')

figure
plot3(liver_snap(:,13),liver_snap(:,14),liver_snap(:,15),18,liver_snap(:,17))

