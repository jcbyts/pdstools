
%% Replace quiver arrowheads using text annotation

[x,y]=meshgrid(-15:2:15);
mask=mvnpdf([x(:) y(:)], [0 0], diag([40 40]));
u=10;
v=5;

figure(1); clf
qh=quiver(x(:), y(:), u.*mask, v.*mask);
hold on
headwidth=150;
headlength=150;
clr=(mask(:)/max(mask))*[.8 0 .8]; % clr can be 1x3 or nx3
pdsa.quiverArrowFix(qh, headwidth, headlength, clr, 'headstyle', 'deltoid')