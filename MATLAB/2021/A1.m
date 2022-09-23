deltat=0.125;
t=0;
Tf=10;
m0=1000;
vi=0;
ve=-2000;
cd=0.1;
g=9.81;
dVdt=0;
vpoints=[vi];
mpoints=[];
tpoints=[];

 while t <= Tf
    if t <= 0.4*Tf
        me = 100;
        mpoint = m0-(me*(t));
        dVdt = (((me/(mpoint)) * (vpoints(end)-(ve)))) - (g) - ((cd/(mpoint)) * (abs(vpoints(end)) * vpoints(end)));
        vpoint = (deltat*(dVdt)) + vpoints(end);  
    elseif t > .4*Tf
        dVdt = - (g) - ((cd/(mpoint)) * (abs(vpoints(end)) * vpoints(end)));
        vpoint = (deltat*(dVdt)) + vpoints(end); 
    end
    
    tpoints(end+1) = t;
    vpoints(end+1) = vpoint;
    mpoints(end+1) = mpoint;
    t = t + deltat;
    
 end
 
figure(1)
plot(tpoints,mpoints,'DisplayName','f()x','LineWidth',2);
grid on;

vpoints(end)=[];
figure(2)
plot(tpoints,vpoints,'DisplayName','f()x','LineWidth',2);
grid on;
