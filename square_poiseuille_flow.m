close all

v_full = zeros(35,35);

sz = 9;

%velocity field
v = zeros(sz,sz);

h = 18;
l = 18;

%y and z axes
y = linspace(0,h,sz);
z = linspace(0,l,sz);



%pressure gradient as given by the infinite plate poiseuille flow
%150 is linear blood velocity
% blood viscosity : 3 cP -> 0.03 P -> 0.1 Pa s = 1 P -> 0.003 Pa s
% vessel width = 10 µm
G = 2*150*8*0.003/18/18;

for i=1:sz
  for j=1:sz
    %first term
    v(i,j) = G/(2*0.003)*(y(i)*(h-y(i)));

    for n=1:100
        beta = (2*n-1)*pi/h;
       v(i,j) = v(i,j) - (4*G*h^2)/(0.003*pi^3)*1/((2*n-1)^3)*((sinh(beta*z(j))+ sinh(beta*(l-z(j))))/(sinh(beta*l)))*sin(beta*y(i));
    endfor
  endfor
endfor

figure
imagesc(y,z,v);

figure
plot(y,v(:,round(sz/2)))

figure
plot(z,v(round(sz/2),:))


v_full(4:12,4:12) = v;
v_full(v_full<0)=0;
v_full = vec(v_full);

filename = "velocity_field";
fid = fopen (filename, "w");
fwrite(fid,v_full,"double")
fclose(fid);
%save("velocity_field","v_full","-ascii");

