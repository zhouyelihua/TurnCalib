T=[-0.22891259728559904	15.847950494563172	734.5971220956136]'
r=[0.018710803963259883	0.96961672900583951	-0.24391249383142491]
theta =20/180*pi
R =rodrigues(r .* theta)
for i=1:2
    theta_i=theta*i;
    R =rodrigues(r .* theta_i);
    eval(['fid=fopen(''output' num2str(i) '.asc'',''w'');']);
    eval(['[data_X,data_Y,data_Z,B,G,yyyy]=textread(''' num2str(i) '.asc'',''%f %f %f %f %f %f'');']);
    sizeX=size(data_X);
    n_points=sizeX(1);
    for i=1:n_points
       coord=[data_X(i);data_Y(i);data_Z(i)];
       temp=coord-T;
       temp=R*temp;
       temp=temp+T;
       fprintf(fid,'%f\t%f\t%f\t\r\n',temp(1),temp(2),temp(3));
    end
    fclose(fid);
end
