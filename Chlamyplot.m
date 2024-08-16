function [S, E]= Chlamyplot(xc,yc,zc,length,thickness,bodyrot,bodyangle,eyerot,eyeangle,transparent,ratio)

[x, y, z] = ellipsoid(0,0,0,length/2,thickness/2,thickness/2,30);
t = hgtransform;
S = surfl(x, y, z);
set(S,'parent',t);
% rotate(S,bodyrot,bodyangle)
R=makehgtform('axisrotate',bodyrot,bodyangle);
Tx = makehgtform('translate',[xc yc zc]);
t.Matrix = Tx*R;

set(S,'FaceColor',[0.4 1 0],'FaceAlpha',transparent,'EdgeColor','none'); 

