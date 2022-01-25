import three;
import graph3;
import palette;
import x11colors;
currentprojection=orthographic(3,8,1.5);
settings.outformat = "pdf";
settings.render = 16;
size(200);

void draw_grid(triple origin, pen p=currentpen)
{
    triple pa = origin;
    int n = 9;
    triple pb = origin + (1,0,0);
    triple pc = pb + (0,1,0);
    triple pd = pc + (-1,0,0);
    draw(pa--pb--pc--pd--cycle);
    for (int i = 0; i < n; ++i)
    {
        draw((origin+(i*1.0/(n-1),0,0))--(origin+(i*1.0/(n-1),1,0)), dotted);
        draw((origin+(0,i*1.0/(n-1),0))--(origin+(1,i*1.0/(n-1),0)), dotted);
    }
}

real f(pair z) { return 0.3*((z.x - .5)^2 + (z.y - .5)^2); }
surface s = surface(f, (0,-0.5), (1.5, 1.0),nx=12,ny=12);
pen[] p = {MidnightBlue,LightSeaGreen,Yellow,Tomato};
for (var pp : p) {
    p = p + opacity(1);
}
s.colors(palette(s.map(zpart), Gradient(100 ... p)));
draw(s, meshpen=black+thick(), light=nolight);

//draw((0.5, 0.5, 0)--(.5, .5 ,-.8), dotted);

draw_grid((0,0,0));
draw_grid((0.5,-0.5,0));

//pair[] samples = {(1,0), (1, 1)};
//for (var s : samples) {
//    real fx;
//    fx = f(s);
//    draw((s.x, s.y, 0)--(s.x, s.y, fx));
//    dot((s.x, s.y, fx));
//}

// test connection
draw((1,0,0)..(1,0,-.5)..(1,0.5,0), EndArrow3);


