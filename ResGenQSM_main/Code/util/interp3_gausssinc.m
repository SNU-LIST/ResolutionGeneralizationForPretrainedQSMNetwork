function result = interp3_gausssinc(img,matrix_size_in,matrix_size_out,voxel_size_in,voxel_size_out,stdev)
if(~exist('stdev','var'))
    stdev = inf;
end

% original points
x = [-matrix_size_in(1)/2:matrix_size_in(1)/2-1]*voxel_size_in(1);
y = [-matrix_size_in(2)/2:matrix_size_in(2)/2-1]*voxel_size_in(2);
z = [-matrix_size_in(3)/2:matrix_size_in(3)/2-1]*voxel_size_in(3);
% interpolation query points
xq = [-matrix_size_out(1)/2:matrix_size_out(1)/2-1]*voxel_size_out(1);
yq = [-matrix_size_out(2)/2:matrix_size_out(2)/2-1]*voxel_size_out(2);
zq = [-matrix_size_out(3)/2:matrix_size_out(3)/2-1]*voxel_size_out(3);

[Xq,X] = ndgrid(xq,x);
[Yq,Y] = ndgrid(yq,y);
[Zq,Z] = ndgrid(zq,z);


a = pagemtimes(sinc((Xq-X)/voxel_size_in(1)).*exp(-1/2*((Xq-X)/voxel_size_in(1)/stdev).^2),img);
b = permute(pagemtimes(sinc((Yq-Y)/voxel_size_in(2)).*exp(-1/2*((Yq-Y)/voxel_size_in(2)/stdev).^2),permute(a,[2,1,3])),[2,1,3]);
result = permute(pagemtimes(sinc((Zq-Z)/voxel_size_in(3)).*exp(-1/2*((Zq-Z)/voxel_size_in(3)/stdev).^2),permute(b,[3,2,1])),[3,2,1]);
end