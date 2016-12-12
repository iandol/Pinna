function Pinna_grating_match(condition,match_time)
global sM w wrect ifi waitFrames  ang1 ang2  black ...
	spa spb spc ok esc cKey onFrames approach r1Origin abandon ...
  xc xxc yc ovalRect r1Match  eleTexMatch1 txtColorMat ...
	shiftAng1 shiftAng2 numChoice ppd breakLoop

abandon = 1;
flag = 0;
num_rings = 3;
% expansion = 1;
eachConditionSecs = match_time;  %>10s abandon
colorMat = txtColorMat(:,:,1); % all white
%%%%%move_speed = 1
move_speed = 1 * ppd;
r_dis = r1Match-r1Origin;
onFrames = fix((r_dis + sqrt(r_dis*r_dis + 4*r_dis*move_speed*ifi))/(2*move_speed*ifi));

%%%%%shiftAng = 20deg
shiftAng(1:3) = 0;
if condition == 1  %exp
    approach = 1;
    shiftAng(1) =  -pi*(20.*ifi)/180;%-shiftAng1;
    shiftAng(2) = 0;
    shiftAng(3) = pi*(20.*ifi)/180;%shiftAng1;
end
if condition == 2  %con
    approach = 0;
    shiftAng(1) = -pi*(20.*ifi)/180;%-shiftAng1;
    shiftAng(2) = 0;
    shiftAng(3) = pi*(20.*ifi)/180;%shiftAng1;
end
if condition == 3
    approach = 2;
    %     shiftAng(1:3) =  shiftAng1; % CW
    shiftAng(1) = pi*(20.*ifi)/180;%shiftAng1;
    shiftAng(2) = pi*(20.*ifi)/180;%shiftAng1;
    shiftAng(3) = pi*(20.*ifi)/180;%shiftAng1;
end
if condition == 4
    approach = 2;
    %     shiftAng(1:3) =  shiftAng2; % CCW
    shiftAng(1) = -pi*(20.*ifi)/180;%shiftAng2;
    shiftAng(2) = -pi*(20.*ifi)/180;%shiftAng2;
    shiftAng(3) = -pi*(20.*ifi)/180;%shiftAng2;
end
flag = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%expansion
i = 1;  %when only has one ring
ii = 1;   %%%%%%%%%%%%%%%%%%%
jj = 1;
kk = 1;
ll = 1;
mm = 1;
nn = 1;
oo = 1;
pp = 1;
qq = 1;
rr = 1;
j = 1;
k = 1;
l = 1;
m = 1;
n = 1;
o = 1;
p = 1;
q = 1;
r = 1;
if num_rings ==2
    j = round(onFrames/2); %2、
    %     jj = j;
end
if num_rings == 3
    j = round(onFrames/3); %li:入 确保i j k为整数 保证轨迹相同，且能实现i == onFrames
    %     jj = j;
    k = round(2*onFrames/3);
    %     kk = k;
end
if num_rings == 4
    j = round(onFrames/4);
    %     jj = j;
    k = round(2*onFrames/4);
    %     kk = k;
    l = round(3*onFrames/4);
    %     ll = l;
end
if num_rings == 5
    j = round(onFrames/5);
    %     jj = j;
    k = round(2*onFrames/5);
    %     kk = k;
    l = round(3*onFrames/5);
    %     ll = l;
    m = round(4*onFrames/5);
    %     mm = m;
end
if num_rings == 6
    j = round(onFrames/6);
    k = round(2*onFrames/6);
    l = round(3*onFrames/6);
    m = round(4*onFrames/6);
    n = round(5*onFrames/6);
    %     jj = j;
    %     kk = k;
    %     ll = l;
    %     mm = m;
    %     nn = n;
end
if num_rings == 7
    j = round(onFrames/7);
    k = round(2*onFrames/7);
    l = round(3*onFrames/7);
    m = round(4*onFrames/7);
    n = round(5*onFrames/7);
    o = round(6*onFrames/7);
    %     jj = j;
    %     kk = k;
    %     ll = l;
    %     mm = m;
    %     nn = n;
    %     oo = o;
end
if num_rings == 8
    j = round(onFrames/8);
    k = round(2*onFrames/8);
    l = round(3*onFrames/8);
    m = round(4*onFrames/8);
    n = round(5*onFrames/8);
    o = round(6*onFrames/8);
    p = round(7*onFrames/8);
    %     jj = j;
    %     kk = k;
    %     ll = l;
    %     mm = m;
    %     nn = n;
    %     oo = o;
    %     pp = p;
end
if num_rings == 9
    j = round(onFrames/9);
    k = round(2*onFrames/9);
    l = round(3*onFrames/9);
    m = round(4*onFrames/9);
    n = round(5*onFrames/9);
    o = round(6*onFrames/9);
    p = round(7*onFrames/9);
    q = round(8*onFrames/9);
    %     jj = j;
    %     kk = k;
    %     ll = l;
    %     mm = m;
    %     nn = n;
    %     oo = o;
    %     pp = p;
    %     qq = q;
end
if num_rings == 10
    j = round(onFrames/10);
    k = round(2*onFrames/10);
    l = round(3*onFrames/10);
    m = round(4*onFrames/10);
    n = round(5*onFrames/10);
    o = round(6*onFrames/10);
    p = round(7*onFrames/10);
    q = round(8*onFrames/10);
    r = round(9*onFrames/10);
    %     jj = j;
    %     kk = k;
    %     ll = l;
    %     mm = m;
    %     nn = n;
    %     oo = o;
    %     pp = p;
    %     qq = q;
    %     rr = r;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if approach == 2   %%每个trial ,ii ,jj ……初始化
    %     ii = i;
    %     jj = j;
    %     kk = k;
    %     ll = l;
    %     mm = m;
    %     nn = n;
    %     oo = o;
    %     pp = p;
    %     qq = q;
    %     rr = r;
    
end

vbl = Screen('Flip',w);
vblendtime = vbl + eachConditionSecs;
max_r = r1Match;
min_size = 5;
max_size = 35;
r_a = 2*(max_r-r1Origin)/(onFrames)^2;
size_a = 2*(max_size-min_size)/(onFrames)^2;
if condition==1 || condition ==2
    while vbl < vblendtime
        
        %%%%%%%%%%1st ring
        r1 = round(r1Origin + 0.5*r_a*i^2);
        side1P(1) = round(min_size + 0.5*size_a*i^2);
        side1P(2) = round(min_size + 0.5*size_a*i^2);
        %area1
        dstRect1a  = [xxc(1)+r1*sin(mod(ang1+(ii-1)*shiftAng(1),2*pi))-side1P(2)/2;yc(1)-r1*cos(mod(ang1+(ii-1)*shiftAng(1),2*pi))-side1P(1)/2;...
            xxc(1)+r1*sin(mod(ang1+(ii-1)*shiftAng(1),2*pi))+side1P(2)/2;yc(1)-r1*cos(mod(ang1+(ii-1)*shiftAng(1),2*pi))+side1P(1)/2];
        %area2
        dstRect1b  = [xxc(2)+r1*sin(mod(ang1+(ii-1)*shiftAng(2),2*pi))-side1P(2)/2;yc(1)-r1*cos(mod(ang1+(ii-1)*shiftAng(2),2*pi))-side1P(1)/2;...
            xxc(2)+r1*sin(mod(ang1+(ii-1)*shiftAng(2),2*pi))+side1P(2)/2;yc(1)-r1*cos(mod(ang1+(ii-1)*shiftAng(2),2*pi))+side1P(1)/2];
        %area3
        dstRect1c  = [xxc(3)+r1*sin(mod(ang1+(ii-1)*shiftAng(3),2*pi))-side1P(2)/2;yc(1)-r1*cos(mod(ang1+(ii-1)*shiftAng(3),2*pi))-side1P(1)/2;...
            xxc(3)+r1*sin(mod(ang1+(ii-1)*shiftAng(3),2*pi))+side1P(2)/2;yc(1)-r1*cos(mod(ang1+(ii-1)*shiftAng(3),2*pi))+side1P(1)/2];
        
        %%%%%%%%%2 ring
        if num_rings >=2
            r2 = round(r1Origin + 0.5*r_a*j^2);
            side2P(1) = round(min_size + 0.5*size_a*j^2);
            side2P(2) = round(min_size + 0.5*size_a*j^2);
            %area1
            dstRect2a = [xxc(1)+r2*sin(mod(ang2+(jj-1)*shiftAng(1),2*pi))-side2P(2)/2;yc(1)-r2*cos(mod(ang2+(jj-1)*shiftAng(1),2*pi))-side2P(1)/2;...
                xxc(1)+r2*sin(mod(ang2+(jj-1)*shiftAng(1),2*pi))+side2P(2)/2;yc(1)-r2*cos(mod(ang2+(jj-1)*shiftAng(1),2*pi))+side2P(1)/2];
            %area2
            dstRect2b  = [xxc(2)+r2*sin(mod(ang2+(jj-1)*shiftAng(2),2*pi))-side2P(2)/2;yc(1)-r2*cos(mod(ang2+(jj-1)*shiftAng(2),2*pi))-side2P(1)/2;...
                xxc(2)+r2*sin(mod(ang2+(jj-1)*shiftAng(2),2*pi))+side2P(2)/2;yc(1)-r2*cos(mod(ang2+(jj-1)*shiftAng(2),2*pi))+side2P(1)/2];
            %area3
            dstRect2c  = [xxc(3)+r2*sin(mod(ang2+(jj-1)*shiftAng(3),2*pi))-side2P(2)/2;yc(1)-r2*cos(mod(ang2+(jj-1)*shiftAng(3),2*pi))-side2P(1)/2;...
                xxc(3)+r2*sin(mod(ang2+(jj-1)*shiftAng(3),2*pi))+side2P(2)/2;yc(1)-r2*cos(mod(ang2+(jj-1)*shiftAng(3),2*pi))+side2P(1)/2];
            
        end
        %%%%%%%%%3 ring
        if num_rings >=3
            r3 = round(r1Origin + 0.5*r_a*k^2);
            side3P(1) = round(min_size + 0.5*size_a*k^2);
            side3P(2) = round(min_size + 0.5*size_a*k^2);
            %area1
            dstRect3a = [xxc(1)+r3*sin(mod(ang1+(kk-1)*shiftAng(1),2*pi))-side3P(2)/2;yc(1)-r3*cos(mod(ang1+(kk-1)*shiftAng(1),2*pi))-side3P(1)/2;...
                xxc(1)+r3*sin(mod(ang1+(kk-1)*shiftAng(1),2*pi))+side3P(2)/2;yc(1)-r3*cos(mod(ang1+(kk-1)*shiftAng(1),2*pi))+side3P(1)/2];
            %area2
            dstRect3b  = [xxc(2)+r3*sin(mod(ang1+(kk-1)*shiftAng(2),2*pi))-side3P(2)/2;yc(1)-r3*cos(mod(ang1+(kk-1)*shiftAng(2),2*pi))-side3P(1)/2;...
                xxc(2)+r3*sin(mod(ang1+(kk-1)*shiftAng(2),2*pi))+side3P(2)/2;yc(1)-r3*cos(mod(ang1+(kk-1)*shiftAng(2),2*pi))+side3P(1)/2];
            %area3
            dstRect3c  = [xxc(3)+r3*sin(mod(ang1+(kk-1)*shiftAng(3),2*pi))-side3P(2)/2;yc(1)-r3*cos(mod(ang1+(kk-1)*shiftAng(3),2*pi))-side3P(1)/2;...
                xxc(3)+r3*sin(mod(ang1+(kk-1)*shiftAng(3),2*pi))+side3P(2)/2;yc(1)-r3*cos(mod(ang1+(kk-1)*shiftAng(3),2*pi))+side3P(1)/2];
            %     %area4
            
        end
        %%%%%%%%%4 ring
        if num_rings >= 4
            r4 = round(r1Origin + 0.5*r_a*l^2);
            side4P(1) = round(min_size + 0.5*size_a*l^2);
            side4P(2) = round(min_size + 0.5*size_a*l^2);
            dstRect4a = [xxc(1)+r4*sin(mod(ang2+(ll-1)*shiftAng(1),2*pi))-side4P(2)/2;yc(1)-r4*cos(mod(ang2+(ll-1)*shiftAng(1),2*pi))-side4P(1)/2;...
                xxc(1)+r4*sin(mod(ang2+(ll-1)*shiftAng(1),2*pi))+side4P(2)/2;yc(1)-r4*cos(mod(ang2+(ll-1)*shiftAng(1),2*pi))+side4P(1)/2];
            %area2
            dstRect4b  = [xxc(2)+r4*sin(mod(ang2+(ll-1)*shiftAng(2),2*pi))-side4P(2)/2;yc(1)-r4*cos(mod(ang2+(ll-1)*shiftAng(2),2*pi))-side4P(1)/2;...
                xxc(2)+r4*sin(mod(ang2+(ll-1)*shiftAng(2),2*pi))+side4P(2)/2;yc(1)-r4*cos(mod(ang2+(ll-1)*shiftAng(2),2*pi))+side4P(1)/2];
            %area3
            dstRect4c  = [xxc(3)+r4*sin(mod(ang2+(ll-1)*shiftAng(3),2*pi))-side4P(2)/2;yc(1)-r4*cos(mod(ang2+(ll-1)*shiftAng(3),2*pi))-side4P(1)/2;...
                xxc(3)+r4*sin(mod(ang2+(ll-1)*shiftAng(3),2*pi))+side4P(2)/2;yc(1)-r4*cos(mod(ang2+(ll-1)*shiftAng(3),2*pi))+side4P(1)/2];
            
        end
        %%%%%%%%%5 ring
        if num_rings >=5
            r5 = round(r1Origin + 0.5*r_a*m^2);
            side5P(1) = round(min_size + 0.5*size_a*m^2);
            side5P(2) = round(min_size + 0.5*size_a*m^2);
            dstRect5a = [xxc(1)+r5*sin(mod(ang1+(mm-1)*shiftAng(1),2*pi))-side5P(2)/2;yc(1)-r5*cos(mod(ang1+(mm-1)*shiftAng(1),2*pi))-side5P(1)/2;...
                xxc(1)+r5*sin(mod(ang1+(mm-1)*shiftAng(1),2*pi))+side5P(2)/2;yc(1)-r5*cos(mod(ang1+(mm-1)*shiftAng(1),2*pi))+side5P(1)/2];
            %area2
            dstRect5b  = [xxc(2)+r5*sin(mod(ang1+(mm-1)*shiftAng(2),2*pi))-side5P(2)/2;yc(1)-r5*cos(mod(ang1+(mm-1)*shiftAng(2),2*pi))-side5P(1)/2;...
                xxc(2)+r5*sin(mod(ang1+(mm-1)*shiftAng(2),2*pi))+side5P(2)/2;yc(1)-r5*cos(mod(ang1+(mm-1)*shiftAng(2),2*pi))+side5P(1)/2];
            %area3
            dstRect5c  = [xxc(3)+r5*sin(mod(ang1+(mm-1)*shiftAng(3),2*pi))-side5P(2)/2;yc(1)-r5*cos(mod(ang1+(mm-1)*shiftAng(3),2*pi))-side5P(1)/2;...
                xxc(3)+r5*sin(mod(ang1+(mm-1)*shiftAng(3),2*pi))+side5P(2)/2;yc(1)-r5*cos(mod(ang1+(mm-1)*shiftAng(3),2*pi))+side5P(1)/2];
            
        end
        %%%%%%%%%6 ring
        if num_rings >=6
            r6 = round(r1Origin + 0.5*r_a*n^2);
            side6P(1) = round(min_size + 0.5*size_a*n^2);
            side6P(2) = round(min_size + 0.5*size_a*n^2);
            dstRect6a = [xxc(1)+r6*sin(mod(ang2+(nn-1)*shiftAng(1),2*pi))-side6P(2)/2;yc(1)-r6*cos(mod(ang2+(nn-1)*shiftAng(1),2*pi))-side6P(1)/2;...
                xxc(1)+r6*sin(mod(ang2+(nn-1)*shiftAng(1),2*pi))+side6P(2)/2;yc(1)-r6*cos(mod(ang2+(nn-1)*shiftAng(1),2*pi))+side6P(1)/2];
            %area2
            dstRect6b  = [xxc(2)+r6*sin(mod(ang2+(nn-1)*shiftAng(2),2*pi))-side6P(2)/2;yc(1)-r6*cos(mod(ang2+(nn-1)*shiftAng(2),2*pi))-side6P(1)/2;...
                xxc(2)+r6*sin(mod(ang2+(nn-1)*shiftAng(2),2*pi))+side6P(2)/2;yc(1)-r6*cos(mod(ang2+(nn-1)*shiftAng(2),2*pi))+side6P(1)/2];
            %area3
            dstRect6c  = [xxc(3)+r6*sin(mod(ang2+(nn-1)*shiftAng(3),2*pi))-side6P(2)/2;yc(1)-r6*cos(mod(ang2+(nn-1)*shiftAng(3),2*pi))-side6P(1)/2;...
                xxc(3)+r6*sin(mod(ang2+(nn-1)*shiftAng(3),2*pi))+side6P(2)/2;yc(1)-r6*cos(mod(ang2+(nn-1)*shiftAng(3),2*pi))+side6P(1)/2];
            
        end
        
        
        %%%%%%%%%7 ring
        if num_rings >=7
            r7 = round(r1Origin + 0.5*r_a*o^2);
            side7P(1) = round(min_size + 0.5*size_a*o^2);
            side7P(2) = round(min_size + 0.5*size_a*o^2);
            dstRect7a = [xxc(1)+r7*sin(mod(ang1+(oo-1)*shiftAng(1),2*pi))-side7P(2)/2;yc(1)-r7*cos(mod(ang1+(oo-1)*shiftAng(1),2*pi))-side7P(1)/2;...
                xxc(1)+r7*sin(mod(ang1+(oo-1)*shiftAng(1),2*pi))+side7P(2)/2;yc(1)-r7*cos(mod(ang1+(oo-1)*shiftAng(1),2*pi))+side7P(1)/2];
            %area2
            dstRect7b  = [xxc(2)+r7*sin(mod(ang1+(oo-1)*shiftAng(2),2*pi))-side7P(2)/2;yc(1)-r7*cos(mod(ang1+(oo-1)*shiftAng(2),2*pi))-side7P(1)/2;...
                xxc(2)+r7*sin(mod(ang1+(oo-1)*shiftAng(2),2*pi))+side7P(2)/2;yc(1)-r7*cos(mod(ang1+(oo-1)*shiftAng(2),2*pi))+side7P(1)/2];
            %area3
            dstRect7c  = [xxc(3)+r7*sin(mod(ang1+(oo-1)*shiftAng(3),2*pi))-side7P(2)/2;yc(1)-r7*cos(mod(ang1+(oo-1)*shiftAng(3),2*pi))-side7P(1)/2;...
                xxc(3)+r7*sin(mod(ang1+(oo-1)*shiftAng(3),2*pi))+side7P(2)/2;yc(1)-r7*cos(mod(ang1+(oo-1)*shiftAng(3),2*pi))+side7P(1)/2];
            
        end
        %%%%%%%%%8 ring
        if num_rings >=8
            r8 = round(r1Origin + 0.5*r_a*p^2);
            side8P(1) = round(min_size + 0.5*size_a*p^2);
            side8P(2) = round(min_size + 0.5*size_a*p^2);
            dstRect8a = [xxc(1)+r8*sin(mod(ang2+(pp-1)*shiftAng(1),2*pi))-side8P(2)/2;yc(1)-r8*cos(mod(ang2+(pp-1)*shiftAng(1),2*pi))-side8P(1)/2;...
                xxc(1)+r8*sin(mod(ang2+(pp-1)*shiftAng(1),2*pi))+side8P(2)/2;yc(1)-r8*cos(mod(ang2+(pp-1)*shiftAng(1),2*pi))+side8P(1)/2];
            %area2
            dstRect8b  = [xxc(2)+r8*sin(mod(ang2+(pp-1)*shiftAng(2),2*pi))-side8P(2)/2;yc(1)-r8*cos(mod(ang2+(pp-1)*shiftAng(2),2*pi))-side8P(1)/2;...
                xxc(2)+r8*sin(mod(ang2+(pp-1)*shiftAng(2),2*pi))+side8P(2)/2;yc(1)-r8*cos(mod(ang2+(pp-1)*shiftAng(2),2*pi))+side8P(1)/2];
            %area3
            dstRect8c  = [xxc(3)+r8*sin(mod(ang2+(pp-1)*shiftAng(3),2*pi))-side8P(2)/2;yc(1)-r8*cos(mod(ang2+(pp-1)*shiftAng(3),2*pi))-side8P(1)/2;...
                xxc(3)+r8*sin(mod(ang2+(pp-1)*shiftAng(3),2*pi))+side8P(2)/2;yc(1)-r8*cos(mod(ang2+(pp-1)*shiftAng(3),2*pi))+side8P(1)/2];
            
            
            %%%%%%%%%9 ring
            if num_rings >=9
                r9 = round(r1Origin + 0.5*r_a*q^2);
                side9P(1) = round(min_size + 0.5*size_a*q^2);
                side9P(2) = round(min_size + 0.5*size_a*q^2);
                dstRect9a = [xxc(1)+r9*sin(mod(ang1+(qq-1)*shiftAng(1),2*pi))-side9P(2)/2;yc(1)-r9*cos(mod(ang1+(qq-1)*shiftAng(1),2*pi))-side9P(1)/2;...
                    xxc(1)+r9*sin(mod(ang1+(qq-1)*shiftAng(1),2*pi))+side9P(2)/2;yc(1)-r9*cos(mod(ang1+(qq-1)*shiftAng(1),2*pi))+side9P(1)/2];
                %area2
                dstRect9b  = [xxc(2)+r9*sin(mod(ang1+(qq-1)*shiftAng(2),2*pi))-side9P(2)/2;yc(1)-r9*cos(mod(ang1+(qq-1)*shiftAng(2),2*pi))-side9P(1)/2;...
                    xxc(2)+r9*sin(mod(ang1+(qq-1)*shiftAng(2),2*pi))+side9P(2)/2;yc(1)-r9*cos(mod(ang1+(qq-1)*shiftAng(2),2*pi))+side9P(1)/2];
                %area3
                dstRect9c  = [xxc(3)+r9*sin(mod(ang1+(qq-1)*shiftAng(3),2*pi))-side9P(2)/2;yc(1)-r9*cos(mod(ang1+(qq-1)*shiftAng(3),2*pi))-side9P(1)/2;...
                    xxc(3)+r9*sin(mod(ang1+(qq-1)*shiftAng(3),2*pi))+side9P(2)/2;yc(1)-r9*cos(mod(ang1+(qq-1)*shiftAng(3),2*pi))+side9P(1)/2];
                
                
                %%%%%%%%%10 ring
                if num_rings>= 10
                    r10 = round(r1Origin + 0.5*r_a*r^2);
                    side10P(1) = round(min_size + 0.5*size_a*r^2);
                    side10P(2) = round(min_size + 0.5*size_a*r^2);
                    dstRect10a = [xxc(1)+r10*sin(mod(ang2+(rr-1)*shiftAng(1),2*pi))-side10P(2)/2;yc(1)-r10*cos(mod(ang2+(rr-1)*shiftAng(1),2*pi))-side10P(1)/2;...
                        xxc(1)+r10*sin(mod(ang2+(rr-1)*shiftAng(1),2*pi))+side10P(2)/2;yc(1)-r10*cos(mod(ang2+(rr-1)*shiftAng(1),2*pi))+side10P(1)/2];
                    %area2
                    dstRect10b  = [xxc(2)+r10*sin(mod(ang2+(rr-1)*shiftAng(2),2*pi))-side10P(2)/2;yc(1)-r10*cos(mod(ang2+(rr-1)*shiftAng(2),2*pi))-side10P(1)/2;...
                        xxc(2)+r10*sin(mod(ang2+(rr-1)*shiftAng(2),2*pi))+side10P(2)/2;yc(1)-r10*cos(mod(ang2+(rr-1)*shiftAng(2),2*pi))+side10P(1)/2];
                    %area3
                    dstRect10c  = [xxc(3)+r10*sin(mod(ang2+(rr-1)*shiftAng(3),2*pi))-side10P(2)/2;yc(1)-r10*cos(mod(ang2+(rr-1)*shiftAng(3),2*pi))-side10P(1)/2;...
                        xxc(3)+r10*sin(mod(ang2+(rr-1)*shiftAng(3),2*pi))+side10P(2)/2;yc(1)-r10*cos(mod(ang2+(rr-1)*shiftAng(3),2*pi))+side10P(1)/2];
                    
                end %10
            end  %9
        end   %8
        
        dstRect1 = [dstRect1a dstRect1b dstRect1c];
        Screen('DrawTextures',w,eleTexMatch1,[],dstRect1,[(180/pi)*(mod(ang1+(ii-1)*shiftAng(1),2*pi)) (180/pi)*(mod(ang1+(ii-1)*shiftAng(2),2*pi)) (180/pi)*(mod(ang1+(ii-1)*shiftAng(3),2*pi))],...
            1,[],[],[],[],[]); % 0 for nearest neighboring filtering, 1 for bilinear filtering
        
        
        if num_rings>=2
            dstRect2 = [dstRect2a dstRect2b dstRect2c];
            Screen('DrawTextures',w,eleTexMatch1,[],dstRect2,[(180/pi)*(mod(ang2+(jj-1)*shiftAng(1),2*pi)) (180/pi)*(mod(ang2+(jj-1)*shiftAng(2),2*pi)) (180/pi)*(mod(ang2+(jj-1)*shiftAng(3),2*pi))],...
                1,[],[],[],[],[]); % 0 for nearest neighboring filtering, 1 for bilinear filtering
        end
        if num_rings>=3
            dstRect3 = [dstRect3a dstRect3b dstRect3c];
            Screen('DrawTextures',w,eleTexMatch1,[],dstRect3,[(180/pi)*(mod(ang1+(kk-1)*shiftAng(1),2*pi)) (180/pi)*(mod(ang1+(kk-1)*shiftAng(2),2*pi)) (180/pi)*(mod(ang1+(kk-1)*shiftAng(3),2*pi))],...
                1,[],[],[],[],[]); % 0 for nearest neighboring filtering, 1 for bilinear filtering
        end
        if num_rings>=4
            dstRect4 = [dstRect4a dstRect4b dstRect4c];
            Screen('DrawTextures',w,eleTexMatch1,[],dstRect4,[(180/pi)*(mod(ang2+(ll-1)*shiftAng(1),2*pi)) (180/pi)*(mod(ang2+(ll-1)*shiftAng(2),2*pi)) (180/pi)*(mod(ang2+(ll-1)*shiftAng(3),2*pi))],...
                1,[],[],[],[],[]); % 0 for nearest neighboring filtering, 1 for bilinear filtering
        end
        if num_rings>=5
            dstRect5 = [dstRect5a dstRect5b dstRect5c];
            Screen('DrawTextures',w,eleTexMatch1,[],dstRect5,[(180/pi)*(mod(ang1+(mm-1)*shiftAng(1),2*pi)) (180/pi)*(mod(ang1+(mm-1)*shiftAng(2),2*pi)) (180/pi)*(mod(ang1+(mm-1)*shiftAng(3),2*pi))],...
                1,[],[],[],[],[]); % 0 for nearest neighboring filtering, 1 for bilinear filtering
        end
        if num_rings>=6
            dstRect6 = [dstRect6a dstRect6b dstRect6c];
            Screen('DrawTextures',w,eleTexMatch1,[],dstRect6,[(180/pi)*(mod(ang2+(nn-1)*shiftAng(1),2*pi)) (180/pi)*(mod(ang2+(nn-1)*shiftAng(2),2*pi)) (180/pi)*(mod(ang2+(nn-1)*shiftAng(3),2*pi))],...
                1,[],[],[],[],[]); % 0 for nearest neighboring filtering, 1 for bilinear filtering
        end
        if num_rings >=7
            dstRect7 = [dstRect7a dstRect7b dstRect7c];
            Screen('DrawTextures',w,eleTexMatch1,[],dstRect7,[(180/pi)*(mod(ang1+(oo-1)*shiftAng(1),2*pi)) (180/pi)*(mod(ang1+(oo-1)*shiftAng(2),2*pi)) (180/pi)*(mod(ang1+(oo-1)*shiftAng(3),2*pi))],...
                1,[],[],[],[],[]); % 0 for nearest neighboring filtering, 1 for bilinear filtering
        end
        if num_rings >=8
            dstRect8 = [dstRect8a dstRect8b dstRect8c];
            Screen('DrawTextures',w,eleTexMatch1,[],dstRect8,[(180/pi)*(mod(ang2+(pp-1)*shiftAng(1),2*pi)) (180/pi)*(mod(ang2+(pp-1)*shiftAng(2),2*pi)) (180/pi)*(mod(ang2+(pp-1)*shiftAng(3),2*pi))],...
                1,[],[],[],[],[]); % 0 for nearest neighboring filtering, 1 for bilinear filtering
            if num_rings >=9
                dstRect9 = [dstRect9a dstRect9b dstRect9c];
                Screen('DrawTextures',w,eleTexMatch1,[],dstRect9,[(180/pi)*(mod(ang1+(qq-1)*shiftAng(1),2*pi)) (180/pi)*(mod(ang1+(qq-1)*shiftAng(2),2*pi)) (180/pi)*(mod(ang1+(qq-1)*shiftAng(3),2*pi))],...
                    1,[],[],[],[],[]); % 0 for nearest neighboring filtering, 1 for bilinear filtering
                if num_rings >=10
                    dstRect10 = [dstRect10a dstRect10b dstRect10c];
                    Screen('DrawTextures',w,eleTexMatch1,[],dstRect10,[(180/pi)*(mod(ang2+(rr-1)*shiftAng(1),2*pi)) (180/pi)*(mod(ang2+(rr-1)*shiftAng(2),2*pi)) (180/pi)*(mod(ang2+(rr-1)*shiftAng(3),2*pi))],...
                        1,[],[],[],[],[]); % 0 for nearest neighboring filtering, 1 for bilinear filtering
                    
                end  %10
            end %9
        end %8
        %              Screen('BlendFunction',w,GL_ONE_MINUS_SRC_ALPHA,GL_SRC_ALPHA,[1 1 1 0]);
        %              Screen('DrawTexture',w,masktex0);
        %              Screen('BlendFunction',w,GL_ONE,GL_ZERO,[1 1 1 1]);
        
        if approach ==1
            i = i+1;
            %             ii = i;
            if i == onFrames
                i = 1;
                %                 ii = i;
            end
            if num_rings>=2
                j = j+1;
                %                 jj = j;
                if j == onFrames
                    j = 1;
                    %                     jj = j;
                end
            end
            if num_rings >=3
                k = k+1;
                %                 kk = k;
                if k == onFrames
                    k = 1;
                    %                     kk = k;
                end
            end
            if num_rings >=4
                l = l+1;
                %                 ll = l;
                if l == onFrames
                    l = 1;
                    %                     ll = l;
                end
            end
            if num_rings >=5
                m = m+1;
                %                 mm = m;
                if m == onFrames
                    m = 1;
                    %                     mm = m;
                end
            end
            if num_rings>=6
                n = n+1;
                %                 nn = n;
                if n == onFrames
                    n = 1;
                    %                     nn = n;
                end
            end
            if num_rings >=7
                o = o+1;
                %                 oo = o;
                if o == onFrames
                    o = 1;
                    %                     oo = o;
                end
                if num_rings >=8
                    p = p +1;
                    %                     pp = p;
                    if p == onFrames
                        p = 1;
                        %                         pp = p;
                    end
                    if num_rings >=9
                        q = q +1;
                        %                         qq = q;
                        if q == onFrames
                            q = 1;
                            %                             qq = q;
                        end
                        if num_rings >= 10
                            r = r+1;
                            %                             rr = r;
                            if r == onFrames
                                r = 1;
                                %                                 rr = r;
                            end
                            
                        end  %10
                    end   %9
                end  %8
            end%7
            
        end
        if approach == 0
            i = i-1;
            if i == 0
                i = onFrames;
            end
            if num_rings>=2
                j = j-1;
                if j == 0
                    j = onFrames;
                end
                if num_rings >=3
                    k = k-1;
                    if k == 0
                        k = onFrames;
                    end
                    if num_rings >=4
                        l = l-1;
                        if l == 0
                            l = onFrames;
                        end
                        if num_rings >=5
                            m = m-1;
                            if m == 0
                                m = onFrames;
                            end
                            if num_rings >=6
                                n = n-1;
                                if n == 0
                                    n = onFrames;
                                end
                                if num_rings >= 7
                                    o = o-1;
                                    if o == 0
                                        o = onFrames;
                                    end
                                    if num_rings >= 8
                                        p = p -1;
                                        if p==0
                                            p = onFrames;
                                        end
                                        
                                        if num_rings >=9
                                            q = q-1;
                                            if q == 0
                                                q = onFrames;
                                            end
                                            
                                            if num_rings >=10
                                                r = r-1;
                                                if r == 0
                                                    r =  onFrames;
                                                end
                                                
                                            end %10
                                        end  %9
                                    end  %8
                                end  %7
                            end   %6
                        end
                    end
                end
            end
            %             ii = i;
            %             jj = j;
            %             kk = k;
            %             ll = l;
            %             mm = m;
            %             nn = n;
            %             oo = o;
            %             pp = p;
            %             qq = q;
            %             rr = r;
            
        end
        
        %         if approach == 2   %  只为了静止时，真实旋转可以实现
        ii = ii+1;
        jj = jj+1;
        kk = kk+1;
        ll = ll+1;
        mm = mm+1;
        nn = nn+1;
        oo = oo +1;
        pp = pp+1;
        qq = qq +1;
        rr = rr +1;
        
        %         end
        vbl = Screen('Flip',w,vbl+(waitFrames-0.5)*ifi);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        sM.drawCross(0.3,[0 0 0 1]);
        
        Screen('DrawText',w,'a',wrect(3)/4,wrect(4)/2+200,colorMat(1,:));
        
        Screen('DrawText',w,'b',wrect(3)*2/4,wrect(4)/2+200,colorMat(2,:));
        
        Screen('DrawText',w,'c',wrect(3)*3/4,wrect(4)/2+200,colorMat(3,:));
        
        %     Screen('DrawText',w,'4',wrect(3)*7/8,0+20,colorMat(4,:));
        %
        %     Screen('DrawText',w,'5',wrect(3)/8,wrect(4)/2+20,colorMat(5,:));
        %
        %     Screen('DrawText',w,'6',wrect(3)*3/8,wrect(4)/2+20,colorMat(6,:));
        %
        %     Screen('DrawText',w,'7',wrect(3)*5/8,wrect(4)/2+20,colorMat(7,:));
        
        [keyisdown,secs,keycode] = KbCheck(-1);
				
				if keycode(cKey)
					abandon = 1;
					flag = 0;
					stopRecording(eL);
					setOffline(eL);
					trackerSetup(eL);
					WaitSecs('YieldSecs',1);
					break
				end
        
        if keycode(spa)
            numChoice = 1;  %real angle speed ,no CCW and CW
            flag = 1;
						abandon = 0;
        end
        if keycode(spb)
            numChoice = 2;  %real angle speed ,no CCW and CW
            flag = 1;
						abandon = 0;
        end
        if keycode(spc)
            numChoice = 3;  %real angle speed ,no CCW and CW
            flag = 1;
						abandon = 0;
        end
        %     if keycode(spd)
        %         number_area_speed = angSpeed1(4);  %real angle speed ,no CCW and CW
        %         flag = 1;
        %     end
        %     if keycode(spe)
        %         number_area_speed = angSpeed1(5);  %real angle speed ,no CCW and CW
        %         flag = 1;
        %     end
        %     if keycode(spf)
        %         number_area_speed = angSpeed1(6);  %real angle speed ,no CCW and CW
        %         flag = 1;
        %     end
        %     if keycode(spg)
        %         number_area_speed = angSpeed1(7);  %real angle speed ,no CCW and CW
        %         flag = 1;
        %     end
        
        if keycode(ok) && flag == 1
            break
        end
        
        if keycode(esc)
            breakLoop = true;
            break;
        end
        
		end
elseif condition == 3 || condition == 4
    i1 = 1; i2 = 1; i3 = 1;
    ii1 = 1; ii2 = 1; ii3 = 1;
    jj1 = 1; jj2 = 1; jj3 = 1;
    kk1 = 1; kk2 = 1; kk3 = 1;
    ll1 = 1; ll2 = 1; ll3 = 1;
    mm1 = 1; mm2 = 1; mm3 = 1;
    nn1 = 1; nn2 = 1; nn3 = 1;
    oo1 = 1; oo2 = 1; oo3 = 1;
    pp1 = 1; pp2 = 1; pp3 = 1;
    qq1 = 1; qq2 = 1; qq3 = 1;
    rr1 = 1; rr2 = 1; rr3 = 1;
    j1 = 1; j2 = 1; j3 = 1;
    k1 = 1; k2 = 1; k3 = 1;
    l1 = 1; l2 = 1; l3 = 1;
    m1 = 1; m2 = 1; m3 = 1;
    n1 = 1; n2 = 1; n3 = 1;
    o1 = 1; o2 = 1; o3 = 1;
    p1 = 1; p2 = 1; p3 = 1;
    q1 = 1; q2 = 1; q3 = 1;
    r1 = 1; r2 = 1; r3 = 1;
    if num_rings == 2
        j1 = round(onFrames/2);j2 = round(onFrames/2);j3 = round(onFrames/2);
        %         jj1 = j1; jj2 = j2; jj3 = j3;
    end
    if num_rings == 3
        j1 = round(onFrames/3);i2 = round(2*onFrames/5);j2 = round(3*onFrames/5);j3 = round(onFrames/3);
        k1 = round(2*onFrames/3); k2 = round(4*onFrames/5); k3 = round(2*onFrames/3);
        
        %         jj1 = j1; jj2 = j2; jj3 = j3;
        %         kk1 = k1;  kk2 = k2;  kk3 = k3;
    end
    if num_rings == 4
        j1 = round(onFrames/4);j2 = round(onFrames/4);j3 = round(onFrames/4);
        k1 = round(2*onFrames/4); k2 = round(2*onFrames/4); k3 = round(2*onFrames/4);
        l1 = round(3*onFrames/4); l2 = round(3*onFrames/4); l3 = round(3*onFrames/4);
        
        %         jj1 = j1; jj2 = j2; jj3 = j3;
        %         kk1 = k1;  kk2 = k2;  kk3 = k3;
        %         ll1 = l1; ll2 = l2; ll3 = l3;
    end
    if num_rings == 5
        j1 = round(onFrames/5);j2 = round(onFrames/5);j3 = round(onFrames/5);
        k1 = round(2*onFrames/5); k2 = round(2*onFrames/5); k3 = round(2*onFrames/5);
        l1 = round(3*onFrames/5); l2 = round(3*onFrames/5); l3 = round(3*onFrames/5);
        m1 = round(4*onFrames/5);m2 = round(4*onFrames/5);m3 = round(4*onFrames/5);
        
        %         jj1 = j1; jj2 = j2; jj3 = j3;
        %         kk1 = k1;  kk2 = k2;  kk3 = k3;
        %         ll1 = l1; ll2 = l2; ll3 = l3;
        %         mm1 = m1; mm2 = m2; mm3 = m3;
    end
    if num_rings == 6
        j1 = round(onFrames/6);j2 = round(onFrames/6);j3 = round(onFrames/6);
        k1 = round(2*onFrames/6); k2 = round(2*onFrames/6); k3 = round(2*onFrames/6);
        l1 = round(3*onFrames/6); l2 = round(3*onFrames/6); l3 = round(3*onFrames/6);
        m1 = round(4*onFrames/6);m2 = round(4*onFrames/6);m3 = round(4*onFrames/6);
        n1 = round(5*onFrames/6);n2 = round(5*onFrames/6);n3 = round(5*onFrames/6);
        
        %         jj1 = j1; jj2 = j2; jj3 = j3;
        %         kk1 = k1;  kk2 = k2;  kk3 = k3;
        %         ll1 = l1; ll2 = l2; ll3 = l3;
        %         mm1 = m1; mm2 = m2; mm3 = m3;
        %         nn1 = n1;nn2 = n2;nn3 = n3;
    end
    if num_rings == 6
        j1 = round(onFrames/6);j2 = round(onFrames/6);j3 = round(onFrames/6);
        k1 = round(2*onFrames/6); k2 = round(2*onFrames/6); k3 = round(2*onFrames/6);
        l1 = round(3*onFrames/6); l2 = round(3*onFrames/6); l3 = round(3*onFrames/6);
        m1 = round(4*onFrames/6);m2 = round(4*onFrames/6);m3 = round(4*onFrames/6);
        n1 = round(5*onFrames/6);n2 = round(5*onFrames/6);n3 = round(5*onFrames/6);
        
        %         jj1 = j1; jj2 = j2; jj3 = j3;
        %         kk1 = k1;  kk2 = k2;  kk3 = k3;
        %         ll1 = l1; ll2 = l2; ll3 = l3;
        %         mm1 = m1; mm2 = m2; mm3 = m3;
        %         nn1 = n1;nn2 = n2;nn3 = n3;
    end
    if num_rings == 7
        j1 = round(onFrames/7);j2 = round(onFrames/7);j3 = round(onFrames/7);
        k1 = round(2*onFrames/7); k2 = round(2*onFrames/7); k3 = round(2*onFrames/7);
        l1 = round(3*onFrames/7); l2 = round(3*onFrames/7); l3 = round(3*onFrames/7);
        m1 = round(4*onFrames/7);m2 = round(4*onFrames/7);m3 = round(4*onFrames/7);
        n1 = round(5*onFrames/7);n2 = round(5*onFrames/7);n3 = round(5*onFrames/7);
        o1 = round(6*onFrames/7); o2 = round(6*onFrames/7); o3= round(6*onFrames/7);
        %         jj1 = j1; jj2 = j2; jj3 = j3;
        %         kk1 = k1;  kk2 = k2;  kk3 = k3;
        %         ll1 = l1; ll2 = l2; ll3 = l3;
        %         mm1 = m1; mm2 = m2; mm3 = m3;
        %         nn1 = n1;nn2 = n2;nn3 = n3;
        %         oo1 = o1; oo2 = o2; oo3 = o3;
    end
    if num_rings == 9
        j1 = round(onFrames/9);j2 = round(onFrames/9);j3 = round(onFrames/9);
        k1 = round(2*onFrames/9);k2 = round(2*onFrames/9);k3 = round(2*onFrames/9);
        l1 = round(3*onFrames/9);l2 = round(3*onFrames/9);l3 = round(3*onFrames/9);
        m1 = round(4*onFrames/9); m2 = round(4*onFrames/9); m3 = round(4*onFrames/9);
        n1 = round(5*onFrames/9);n2 = round(5*onFrames/9);n3 = round(5*onFrames/9);
        o1 = round(6*onFrames/9);o2 = round(6*onFrames/9);o3 = round(6*onFrames/9);
        p1 = round(7*onFrames/9);p2 = round(7*onFrames/9);p3 = round(7*onFrames/9);
        q1 = round(8*onFrames/9);q2 = round(8*onFrames/9);q3 = round(8*onFrames/9);
        %         jj1 = j1;jj2 = j2;jj3 = j3;
        %         kk1 = k1; kk2 = k2; kk3 = k3;
        %         ll1 = l1; ll2 = l2; ll3 = l3;
        %         mm1 = m1;mm2 = m2;mm3 = m3;
        %         nn1 = n1; nn2 = n2; nn3 = n3;
        %         oo1 = o1;oo2 = o2;oo3 = o3;
        %         pp1 = p1;pp2 = p2;pp3 = p3;
        %         qq1 = q1;qq2 = q2;qq3 = q3;
    end
    if num_rings == 10
        j1 = round(onFrames/10);j2 = round(onFrames/10);j3 = round(onFrames/10);
        k1 = round(2*onFrames/10);k2 = round(2*onFrames/10);k3 = round(2*onFrames/10);
        l1 = round(3*onFrames/10); l2 = round(3*onFrames/10); l3 = round(3*onFrames/10);
        m1 = round(4*onFrames/10);m2 = round(4*onFrames/10);m3 = round(4*onFrames/10);
        n1 = round(5*onFrames/10); n2 = round(5*onFrames/10); n3 = round(5*onFrames/10);
        o1 = round(6*onFrames/10); o2 = round(6*onFrames/10); o3 = round(6*onFrames/10);
        p1 = round(7*onFrames/10);p2 = round(7*onFrames/10);p3 = round(7*onFrames/10);
        q1 = round(8*onFrames/10); q2 = round(8*onFrames/10); q3 = round(8*onFrames/10);
        r1 = round(9*onFrames/10); r2 = round(9*onFrames/10); r3 = round(9*onFrames/10);
        %         jj1 = j1; jj2 = j2; jj3 = j3;
        %         kk1 = k1; kk2 = k2; kk3 = k3;
        %         ll1 = l1;  ll2 = l2;  ll3 = l3;
        %         mm1 = m1; mm2 = m2; mm3 = m3;
        %         nn1 = n1;nn2 = n2;nn3 = n3;
        %         oo1 = o1; oo2 = o2; oo3 = o3;
        %         pp1 = p1;  pp2 = p2;  pp3 = p3;
        %         qq1 = q1;qq2 = q2;qq3 = q3;
        %         rr1 = r1;rr2 = r2;rr3 = r3;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     if approach ==2   %%每个trial ,ii ,jj ……初始化
    %         ii1 = i1;ii2 = i2;ii3 = i3;
    %         jj1 = j1;jj2 = j2;jj3 = j3;
    %         kk1 = k1; kk2 = k2; kk3 = k3;
    %         ll1 = l1; ll2 = l2; ll3 = l3;
    %         mm1 = m1; mm2 = m2; mm3 = m3;
    %         nn1 = n1;nn2 = n2;nn3 = n3;
    %         oo1 = o1;oo2 = o2;oo3 = o3;
    %         pp1 = p1; pp2 = p2; pp3 = p3;
    %         qq1 = q1;qq2 = q2;qq3 = q3;
    %         rr1 = r1; rr2 = r2; rr3 = r3;
    %
    %     end
    while vbl < vblendtime
        
        %%%%%%%%%%1st ring
        
        %area1
        radius1 = round(r1Origin + 0.5*r_a*i1^2);
        side1P(1) = round(min_size + 0.5*size_a*i1^2);
        side1P(2) = round(min_size + 0.5*size_a*i1^2);
        dstRect1a  = [xxc(1)+radius1*sin(mod(ang1+(ii1-1)*shiftAng(1),2*pi))-side1P(2)/2;yc(1)-radius1*cos(mod(ang1+(ii1-1)*shiftAng(1),2*pi))-side1P(1)/2;...
            xxc(1)+radius1*sin(mod(ang1+(ii1-1)*shiftAng(1),2*pi))+side1P(2)/2;yc(1)-radius1*cos(mod(ang1+(ii1-1)*shiftAng(1),2*pi))+side1P(1)/2];
        %area2
        radius1 = round(r1Origin + 0.5*r_a*i2^2);
        side1P(1) = round(min_size + 0.5*size_a*i2^2);
        side1P(2) = round(min_size + 0.5*size_a*i2^2);
        dstRect1b  = [xxc(2)+radius1*sin(mod(ang1+(ii2-1)*shiftAng(2),2*pi))-side1P(2)/2;yc(1)-radius1*cos(mod(ang1+(ii2-1)*shiftAng(2),2*pi))-side1P(1)/2;...
            xxc(2)+radius1*sin(mod(ang1+(ii2-1)*shiftAng(2),2*pi))+side1P(2)/2;yc(1)-radius1*cos(mod(ang1+(ii2-1)*shiftAng(2),2*pi))+side1P(1)/2];
        %area3
        radius1 = round(r1Origin + 0.5*r_a*i3^2);
        side1P(1) = round(min_size + 0.5*size_a*i3^2);
        side1P(2) = round(min_size + 0.5*size_a*i3^2);
        dstRect1c  = [xxc(3)+radius1*sin(mod(ang1+(ii3-1)*shiftAng(3),2*pi))-side1P(2)/2;yc(1)-radius1*cos(mod(ang1+(ii3-1)*shiftAng(3),2*pi))-side1P(1)/2;...
            xxc(3)+radius1*sin(mod(ang1+(ii3-1)*shiftAng(3),2*pi))+side1P(2)/2;yc(1)-radius1*cos(mod(ang1+(ii3-1)*shiftAng(3),2*pi))+side1P(1)/2];
        
        %%%%%%%%%2 ring
        if num_rings>=2
            %area1
            radius2 = round(r1Origin + 0.5*r_a*j1^2);
            side2P(1) = round(min_size + 0.5*size_a*j1^2);
            side2P(2) = round(min_size + 0.5*size_a*j1^2);
            dstRect2a = [xxc(1)+radius2*sin(mod(ang2+(jj1-1)*shiftAng(1),2*pi))-side2P(2)/2;yc(1)-radius2*cos(mod(ang2+(jj1-1)*shiftAng(1),2*pi))-side2P(1)/2;...
                xxc(1)+radius2*sin(mod(ang2+(jj1-1)*shiftAng(1),2*pi))+side2P(2)/2;yc(1)-radius2*cos(mod(ang2+(jj1-1)*shiftAng(1),2*pi))+side2P(1)/2];
            %area2
            radius2 = round(r1Origin + 0.5*r_a*j2^2);
            side2P(1) = round(min_size + 0.5*size_a*j2^2);
            side2P(2) = round(min_size + 0.5*size_a*j2^2);
            dstRect2b  = [xxc(2)+radius2*sin(mod(ang2+(jj2-1)*shiftAng(2),2*pi))-side2P(2)/2;yc(1)-radius2*cos(mod(ang2+(jj2-1)*shiftAng(2),2*pi))-side2P(1)/2;...
                xxc(2)+radius2*sin(mod(ang2+(jj2-1)*shiftAng(2),2*pi))+side2P(2)/2;yc(1)-radius2*cos(mod(ang2+(jj2-1)*shiftAng(2),2*pi))+side2P(1)/2];
            %area3
            radius2 = round(r1Origin + 0.5*r_a*j3^2);
            side2P(1) = round(min_size + 0.5*size_a*j3^2);
            side2P(2) = round(min_size + 0.5*size_a*j3^2);
            dstRect2c  = [xxc(3)+radius2*sin(mod(ang2+(jj3-1)*shiftAng(3),2*pi))-side2P(2)/2;yc(1)-radius2*cos(mod(ang2+(jj3-1)*shiftAng(3),2*pi))-side2P(1)/2;...
                xxc(3)+radius2*sin(mod(ang2+(jj3-1)*shiftAng(3),2*pi))+side2P(2)/2;yc(1)-radius2*cos(mod(ang2+(jj3-1)*shiftAng(3),2*pi))+side2P(1)/2];
        end
        %%%%%%%%%3 ring
        if num_rings >=3
            %area1
            radius3 = round(r1Origin + 0.5*r_a*k1^2);
            side3P(1) = round(min_size + 0.5*size_a*k1^2);
            side3P(2) = round(min_size + 0.5*size_a*k1^2);
            dstRect3a = [xxc(1)+radius3*sin(mod(ang1+(kk1-1)*shiftAng(1),2*pi))-side3P(2)/2;yc(1)-radius3*cos(mod(ang1+(kk1-1)*shiftAng(1),2*pi))-side3P(1)/2;...
                xxc(1)+radius3*sin(mod(ang1+(kk1-1)*shiftAng(1),2*pi))+side3P(2)/2;yc(1)-radius3*cos(mod(ang1+(kk1-1)*shiftAng(1),2*pi))+side3P(1)/2];
            %area2
            radius3 = round(r1Origin + 0.5*r_a*k2^2);
            side3P(1) = round(min_size + 0.5*size_a*k2^2);
            side3P(2) = round(min_size + 0.5*size_a*k2^2);
            dstRect3b  = [xxc(2)+radius3*sin(mod(ang1+(kk2-1)*shiftAng(2),2*pi))-side3P(2)/2;yc(1)-radius3*cos(mod(ang1+(kk2-1)*shiftAng(2),2*pi))-side3P(1)/2;...
                xxc(2)+radius3*sin(mod(ang1+(kk2-1)*shiftAng(2),2*pi))+side3P(2)/2;yc(1)-radius3*cos(mod(ang1+(kk2-1)*shiftAng(2),2*pi))+side3P(1)/2];
            %area3
            radius3 = round(r1Origin + 0.5*r_a*k3^2);
            side3P(1) = round(min_size + 0.5*size_a*k3^2);
            side3P(2) = round(min_size + 0.5*size_a*k3^2);
            dstRect3c  = [xxc(3)+radius3*sin(mod(ang1+(kk3-1)*shiftAng(3),2*pi))-side3P(2)/2;yc(1)-radius3*cos(mod(ang1+(kk3-1)*shiftAng(3),2*pi))-side3P(1)/2;...
                xxc(3)+radius3*sin(mod(ang1+(kk3-1)*shiftAng(3),2*pi))+side3P(2)/2;yc(1)-radius3*cos(mod(ang1+(kk3-1)*shiftAng(3),2*pi))+side3P(1)/2];
        end
        %%%%%%%%%4 ring
        if num_rings>=4
            r4 = round(r1Origin + 0.5*r_a*l1^2);
            side4P(1) = round(min_size + 0.5*size_a*l1^2);
            side4P(2) = round(min_size + 0.5*size_a*l1^2);
            dstRect4a = [xxc(1)+r4*sin(mod(ang2+(ll1-1)*shiftAng(1),2*pi))-side4P(2)/2;yc(1)-r4*cos(mod(ang2+(ll1-1)*shiftAng(1),2*pi))-side4P(1)/2;...
                xxc(1)+r4*sin(mod(ang2+(ll1-1)*shiftAng(1),2*pi))+side4P(2)/2;yc(1)-r4*cos(mod(ang2+(ll1-1)*shiftAng(1),2*pi))+side4P(1)/2];
            %area2
            r4 = round(r1Origin + 0.5*r_a*l2^2);
            side4P(1) = round(min_size + 0.5*size_a*l2^2);
            side4P(2) = round(min_size + 0.5*size_a*l2^2);
            dstRect4b  = [xxc(2)+r4*sin(mod(ang2+(ll2-1)*shiftAng(2),2*pi))-side4P(2)/2;yc(1)-r4*cos(mod(ang2+(ll2-1)*shiftAng(2),2*pi))-side4P(1)/2;...
                xxc(2)+r4*sin(mod(ang2+(ll2-1)*shiftAng(2),2*pi))+side4P(2)/2;yc(1)-r4*cos(mod(ang2+(ll2-1)*shiftAng(2),2*pi))+side4P(1)/2];
            %area3
            r4 = round(r1Origin + 0.5*r_a*l3^2);
            side4P(1) = round(min_size + 0.5*size_a*l3^2);
            side4P(2) = round(min_size + 0.5*size_a*l3^2);
            dstRect4c  = [xxc(3)+r4*sin(mod(ang2+(ll3-1)*shiftAng(3),2*pi))-side4P(2)/2;yc(1)-r4*cos(mod(ang2+(ll3-1)*shiftAng(3),2*pi))-side4P(1)/2;...
                xxc(3)+r4*sin(mod(ang2+(ll3-1)*shiftAng(3),2*pi))+side4P(2)/2;yc(1)-r4*cos(mod(ang2+(ll3-1)*shiftAng(3),2*pi))+side4P(1)/2];
        end
        %%%%%%%%%5 ring
        if num_rings >=5
            r5 = round(r1Origin + 0.5*r_a*m1^2);
            side5P(1) = round(min_size + 0.5*size_a*m1^2);
            side5P(2) = round(min_size + 0.5*size_a*m1^2);
            dstRect5a = [xxc(1)+r5*sin(mod(ang1+(mm1-1)*shiftAng(1),2*pi))-side5P(2)/2;yc(1)-r5*cos(mod(ang1+(mm1-1)*shiftAng(1),2*pi))-side5P(1)/2;...
                xxc(1)+r5*sin(mod(ang1+(mm1-1)*shiftAng(1),2*pi))+side5P(2)/2;yc(1)-r5*cos(mod(ang1+(mm1-1)*shiftAng(1),2*pi))+side5P(1)/2];
            %area2
            r5 = round(r1Origin + 0.5*r_a*m2^2);
            side5P(1) = round(min_size + 0.5*size_a*m2^2);
            side5P(2) = round(min_size + 0.5*size_a*m2^2);
            dstRect5b  = [xxc(2)+r5*sin(mod(ang1+(mm2-1)*shiftAng(2),2*pi))-side5P(2)/2;yc(1)-r5*cos(mod(ang1+(mm2-1)*shiftAng(2),2*pi))-side5P(1)/2;...
                xxc(2)+r5*sin(mod(ang1+(mm2-1)*shiftAng(2),2*pi))+side5P(2)/2;yc(1)-r5*cos(mod(ang1+(mm2-1)*shiftAng(2),2*pi))+side5P(1)/2];
            %area3
            r5 = round(r1Origin + 0.5*r_a*m3^2);
            side5P(1) = round(min_size + 0.5*size_a*m3^2);
            side5P(2) = round(min_size + 0.5*size_a*m3^2);
            dstRect5c  = [xxc(3)+r5*sin(mod(ang1+(mm3-1)*shiftAng(3),2*pi))-side5P(2)/2;yc(1)-r5*cos(mod(ang1+(mm3-1)*shiftAng(3),2*pi))-side5P(1)/2;...
                xxc(3)+r5*sin(mod(ang1+(mm3-1)*shiftAng(3),2*pi))+side5P(2)/2;yc(1)-r5*cos(mod(ang1+(mm3-1)*shiftAng(3),2*pi))+side5P(1)/2];
        end
        %%%%%%%%%6 ring
        if num_rings>=6
            r6 = round(r1Origin + 0.5*r_a*n1^2);
            side6P(1) = round(min_size + 0.5*size_a*n1^2);
            side6P(2) = round(min_size + 0.5*size_a*n1^2);
            dstRect6a = [xxc(1)+r6*sin(mod(ang2+(nn1-1)*shiftAng(1),2*pi))-side6P(2)/2;yc(1)-r6*cos(mod(ang2+(nn1-1)*shiftAng(1),2*pi))-side6P(1)/2;...
                xxc(1)+r6*sin(mod(ang2+(nn1-1)*shiftAng(1),2*pi))+side6P(2)/2;yc(1)-r6*cos(mod(ang2+(nn1-1)*shiftAng(1),2*pi))+side6P(1)/2];
            %area2
            r6 = round(r1Origin + 0.5*r_a*n2^2);
            side6P(1) = round(min_size + 0.5*size_a*n2^2);
            side6P(2) = round(min_size + 0.5*size_a*n2^2);
            dstRect6b  = [xxc(2)+r6*sin(mod(ang2+(nn2-1)*shiftAng(2),2*pi))-side6P(2)/2;yc(1)-r6*cos(mod(ang2+(nn2-1)*shiftAng(2),2*pi))-side6P(1)/2;...
                xxc(2)+r6*sin(mod(ang2+(nn2-1)*shiftAng(2),2*pi))+side6P(2)/2;yc(1)-r6*cos(mod(ang2+(nn2-1)*shiftAng(2),2*pi))+side6P(1)/2];
            %area3
            r6 = round(r1Origin + 0.5*r_a*n3^2);
            side6P(1) = round(min_size + 0.5*size_a*n3^2);
            side6P(2) = round(min_size + 0.5*size_a*n3^2);
            dstRect6c  = [xxc(3)+r6*sin(mod(ang2+(nn3-1)*shiftAng(3),2*pi))-side6P(2)/2;yc(1)-r6*cos(mod(ang2+(nn3-1)*shiftAng(3),2*pi))-side6P(1)/2;...
                xxc(3)+r6*sin(mod(ang2+(nn3-1)*shiftAng(3),2*pi))+side6P(2)/2;yc(1)-r6*cos(mod(ang2+(nn3-1)*shiftAng(3),2*pi))+side6P(1)/2];
        end
        
        %%%%%%%%%7 ring
        if num_rings >=7
            r7 = round(r1Origin + 0.5*r_a*o1^2);
            side7P(1) = round(min_size + 0.5*size_a*o1^2);
            side7P(2) = round(min_size + 0.5*size_a*o1^2);
            dstRect7a = [xxc(1)+r7*sin(mod(ang1+(oo1-1)*shiftAng(1),2*pi))-side7P(2)/2;yc(1)-r7*cos(mod(ang1+(oo1-1)*shiftAng(1),2*pi))-side7P(1)/2;...
                xxc(1)+r7*sin(mod(ang1+(oo1-1)*shiftAng(1),2*pi))+side7P(2)/2;yc(1)-r7*cos(mod(ang1+(oo1-1)*shiftAng(1),2*pi))+side7P(1)/2];
            %area2
            r7 = round(r1Origin + 0.5*r_a*o2^2);
            side7P(1) = round(min_size + 0.5*size_a*o2^2);
            side7P(2) = round(min_size + 0.5*size_a*o2^2);
            dstRect7b  = [xxc(2)+r7*sin(mod(ang1+(oo2-1)*shiftAng(2),2*pi))-side7P(2)/2;yc(1)-r7*cos(mod(ang1+(oo2-1)*shiftAng(2),2*pi))-side7P(1)/2;...
                xxc(2)+r7*sin(mod(ang1+(oo2-1)*shiftAng(2),2*pi))+side7P(2)/2;yc(1)-r7*cos(mod(ang1+(oo2-1)*shiftAng(2),2*pi))+side7P(1)/2];
            %area3
            r7 = round(r1Origin + 0.5*r_a*o3^2);
            side7P(1) = round(min_size + 0.5*size_a*o3^2);
            side7P(2) = round(min_size + 0.5*size_a*o3^2);
            dstRect7c  = [xxc(3)+r7*sin(mod(ang1+(oo3-1)*shiftAng(3),2*pi))-side7P(2)/2;yc(1)-r7*cos(mod(ang1+(oo3-1)*shiftAng(3),2*pi))-side7P(1)/2;...
                xxc(3)+r7*sin(mod(ang1+(oo3-1)*shiftAng(3),2*pi))+side7P(2)/2;yc(1)-r7*cos(mod(ang1+(oo3-1)*shiftAng(3),2*pi))+side7P(1)/2];
        end
        %%%%%%%%%8 ring
        if num_rings >=8
            r8 = round(r1Origin + 0.5*r_a*p1^2);
            side8P(1) = round(min_size + 0.5*size_a*p1^2);
            side8P(2) = round(min_size + 0.5*size_a*p1^2);
            dstRect8a = [xxc(1)+r8*sin(mod(ang2+(pp1-1)*shiftAng(1),2*pi))-side8P(2)/2;yc(1)-r8*cos(mod(ang2+(pp1-1)*shiftAng(1),2*pi))-side8P(1)/2;...
                xxc(1)+r8*sin(mod(ang2+(pp1-1)*shiftAng(1),2*pi))+side8P(2)/2;yc(1)-r8*cos(mod(ang2+(pp1-1)*shiftAng(1),2*pi))+side8P(1)/2];
            %area2
            r8 = round(r1Origin + 0.5*r_a*p2^2);
            side8P(1) = round(min_size + 0.5*size_a*p2^2);
            side8P(2) = round(min_size + 0.5*size_a*p2^2);
            dstRect8b  = [xxc(2)+r8*sin(mod(ang2+(pp2-1)*shiftAng(2),2*pi))-side8P(2)/2;yc(1)-r8*cos(mod(ang2+(pp2-1)*shiftAng(2),2*pi))-side8P(1)/2;...
                xxc(2)+r8*sin(mod(ang2+(pp2-1)*shiftAng(2),2*pi))+side8P(2)/2;yc(1)-r8*cos(mod(ang2+(pp2-1)*shiftAng(2),2*pi))+side8P(1)/2];
            %area3
            r8 = round(r1Origin + 0.5*r_a*p3^2);
            side8P(1) = round(min_size + 0.5*size_a*p3^2);
            side8P(2) = round(min_size + 0.5*size_a*p3^2);
            dstRect8c  = [xxc(3)+r8*sin(mod(ang2+(pp3-1)*shiftAng(3),2*pi))-side8P(2)/2;yc(1)-r8*cos(mod(ang2+(pp3-1)*shiftAng(3),2*pi))-side8P(1)/2;...
                xxc(3)+r8*sin(mod(ang2+(pp3-1)*shiftAng(3),2*pi))+side8P(2)/2;yc(1)-r8*cos(mod(ang2+(pp3-1)*shiftAng(3),2*pi))+side8P(1)/2];
            
            
            %%%%%%%%%9 ring
            if num_rings >=9
                r9 = round(r1Origin + 0.5*r_a*q1^2);
                side9P(1) = round(min_size + 0.5*size_a*q1^2);
                side9P(2) = round(min_size + 0.5*size_a*q1^2);
                dstRect9a = [xxc(1)+r9*sin(mod(ang1+(qq1-1)*shiftAng(1),2*pi))-side9P(2)/2;yc(1)-r9*cos(mod(ang1+(qq1-1)*shiftAng(1),2*pi))-side9P(1)/2;...
                    xxc(1)+r9*sin(mod(ang1+(qq1-1)*shiftAng(1),2*pi))+side9P(2)/2;yc(1)-r9*cos(mod(ang1+(qq1-1)*shiftAng(1),2*pi))+side9P(1)/2];
                %area2
                r9 = round(r1Origin + 0.5*r_a*q2^2);
                side9P(1) = round(min_size + 0.5*size_a*q2^2);
                side9P(2) = round(min_size + 0.5*size_a*q2^2);
                dstRect9b  = [xxc(2)+r9*sin(mod(ang1+(qq2-1)*shiftAng(2),2*pi))-side9P(2)/2;yc(1)-r9*cos(mod(ang1+(qq2-1)*shiftAng(2),2*pi))-side9P(1)/2;...
                    xxc(2)+r9*sin(mod(ang1+(qq2-1)*shiftAng(2),2*pi))+side9P(2)/2;yc(1)-r9*cos(mod(ang1+(qq2-1)*shiftAng(2),2*pi))+side9P(1)/2];
                %area3
                r9 = round(r1Origin + 0.5*r_a*q3^2);
                side9P(1) = round(min_size + 0.5*size_a*q3^2);
                side9P(2) = round(min_size + 0.5*size_a*q3^2);
                dstRect9c  = [xxc(3)+r9*sin(mod(ang1+(qq3-1)*shiftAng(3),2*pi))-side9P(2)/2;yc(1)-r9*cos(mod(ang1+(qq3-1)*shiftAng(3),2*pi))-side9P(1)/2;...
                    xxc(3)+r9*sin(mod(ang1+(qq3-1)*shiftAng(3),2*pi))+side9P(2)/2;yc(1)-r9*cos(mod(ang1+(qq3-1)*shiftAng(3),2*pi))+side9P(1)/2];
                
                
                %%%%%%%%%10 ring
                if num_rings>= 10
                    r10 = round(r1Origin + 0.5*r_a*r1^2);
                    side10P(1) = round(min_size + 0.5*size_a*r1^2);
                    side10P(2) = round(min_size + 0.5*size_a*r1^2);
                    dstRect10a = [xxc(1)+r10*sin(mod(ang2+(rr1-1)*shiftAng(1),2*pi))-side10P(2)/2;yc(1)-r10*cos(mod(ang2+(rr1-1)*shiftAng(1),2*pi))-side10P(1)/2;...
                        xxc(1)+r10*sin(mod(ang2+(rr1-1)*shiftAng(1),2*pi))+side10P(2)/2;yc(1)-r10*cos(mod(ang2+(rr1-1)*shiftAng(1),2*pi))+side10P(1)/2];
                    %area2
                    r10 = round(r1Origin + 0.5*r_a*r2^2);
                    side10P(1) = round(min_size + 0.5*size_a*r2^2);
                    side10P(2) = round(min_size + 0.5*size_a*r2^2);
                    dstRect10b  = [xxc(2)+r10*sin(mod(ang2+(rr2-1)*shiftAng(2),2*pi))-side10P(2)/2;yc(1)-r10*cos(mod(ang2+(rr2-1)*shiftAng(2),2*pi))-side10P(1)/2;...
                        xxc(2)+r10*sin(mod(ang2+(rr2-1)*shiftAng(2),2*pi))+side10P(2)/2;yc(1)-r10*cos(mod(ang2+(rr2-1)*shiftAng(2),2*pi))+side10P(1)/2];
                    %area3
                    r10 = round(r1Origin + 0.5*r_a*r3^2);
                    side10P(1) = round(min_size + 0.5*size_a*r3^2);
                    side10P(2) = round(min_size + 0.5*size_a*r3^2);
                    dstRect10c  = [xxc(3)+r10*sin(mod(ang2+(rr3-1)*shiftAng(3),2*pi))-side10P(2)/2;yc(1)-r10*cos(mod(ang2+(rr3-1)*shiftAng(3),2*pi))-side10P(1)/2;...
                        xxc(3)+r10*sin(mod(ang2+(rr3-1)*shiftAng(3),2*pi))+side10P(2)/2;yc(1)-r10*cos(mod(ang2+(rr3-1)*shiftAng(3),2*pi))+side10P(1)/2];
                    
                    
                end %10
            end  %9
        end   %8
        
        
        dstRect1 = [dstRect1a dstRect1b dstRect1c];
        Screen('DrawTextures',w,eleTexMatch1,[],dstRect1,[(180/pi)*(mod(ang1+(ii1-1)*shiftAng(1),2*pi)) (180/pi)*(mod(ang1+(ii2-1)*shiftAng(2),2*pi)) (180/pi)*(mod(ang1+(ii3-1)*shiftAng(3),2*pi))],...
            1,[],[],[],[],[]); % 0 for nearest neighboring filtering, 1 for bilinear filtering
        if num_rings>=2
            dstRect2 = [dstRect2a dstRect2b dstRect2c];
            Screen('DrawTextures',w,eleTexMatch1,[],dstRect2,[(180/pi)*(mod(ang2+(jj1-1)*shiftAng(1),2*pi)) (180/pi)*(mod(ang2+(jj2-1)*shiftAng(2),2*pi)) (180/pi)*(mod(ang2+(jj3-1)*shiftAng(3),2*pi))],...
                1,[],[],[],[],[]); % 0 for nearest neighboring filtering, 1 for bilinear filtering
        end
        if num_rings>=3
            dstRect3 = [dstRect3a dstRect3b dstRect3c];
            Screen('DrawTextures',w,eleTexMatch1,[],dstRect3,[(180/pi)*(mod(ang1+(kk1-1)*shiftAng(1),2*pi)) (180/pi)*(mod(ang1+(kk2-1)*shiftAng(2),2*pi)) (180/pi)*(mod(ang1+(kk3-1)*shiftAng(3),2*pi))],...
                1,[],[],[],[],[]); % 0 for nearest neighboring filtering, 1 for bilinear filtering
        end
        if num_rings>=4
            dstRect4 = [dstRect4a dstRect4b dstRect4c];
            Screen('DrawTextures',w,eleTexMatch1,[],dstRect4,[(180/pi)*(mod(ang2+(ll1-1)*shiftAng(1),2*pi)) (180/pi)*(mod(ang2+(ll2-1)*shiftAng(2),2*pi)) (180/pi)*(mod(ang2+(ll3-1)*shiftAng(3),2*pi))],...
                1,[],[],[],[],[]); % 0 for nearest neighboring filtering, 1 for bilinear filtering
        end
        if num_rings>=5
            dstRect5 = [dstRect5a dstRect5b dstRect5c];
            Screen('DrawTextures',w,eleTexMatch1,[],dstRect5,[(180/pi)*(mod(ang1+(mm1-1)*shiftAng(1),2*pi)) (180/pi)*(mod(ang1+(mm2-1)*shiftAng(2),2*pi)) (180/pi)*(mod(ang1+(mm3-1)*shiftAng(3),2*pi))],...
                1,[],[],[],[],[]); % 0 for nearest neighboring filtering, 1 for bilinear filtering
        end
        if num_rings>=6
            dstRect6 = [dstRect6a dstRect6b dstRect6c];
            Screen('DrawTextures',w,eleTexMatch1,[],dstRect6,[(180/pi)*(mod(ang2+(nn1-1)*shiftAng(1),2*pi)) (180/pi)*(mod(ang2+(nn2-1)*shiftAng(2),2*pi)) (180/pi)*(mod(ang2+(nn3-1)*shiftAng(3),2*pi))],...
                1,[],[],[],[],[]); % 0 for nearest neighboring filtering, 1 for bilinear filtering
        end
        if num_rings >=7
            dstRect7 = [dstRect7a dstRect7b dstRect7c];
            Screen('DrawTextures',w,eleTexMatch1,[],dstRect7,[(180/pi)*(mod(ang1+(oo1-1)*shiftAng(1),2*pi)) (180/pi)*(mod(ang1+(oo2-1)*shiftAng(2),2*pi)) (180/pi)*(mod(ang1+(oo3-1)*shiftAng(3),2*pi))],...
                1,[],[],[],[],[]); % 0 for nearest neighboring filtering, 1 for bilinear filtering
        end
        if num_rings >=8
            dstRect8 = [dstRect8a dstRect8b dstRect8c];
            Screen('DrawTextures',w,eleTexMatch1,[],dstRect8,[(180/pi)*(mod(ang2+(pp1-1)*shiftAng(1),2*pi)) (180/pi)*(mod(ang2+(pp2-1)*shiftAng(2),2*pi)) (180/pi)*(mod(ang2+(pp3-1)*shiftAng(3),2*pi))],...
                1,[],[],[],[],[]); % 0 for nearest neighboring filtering, 1 for bilinear filtering
            if num_rings >=9
                dstRect9 = [dstRect9a dstRect9b dstRect9c];
                Screen('DrawTextures',w,eleTexMatch1,[],dstRect9,[(180/pi)*(mod(ang1+(qq1-1)*shiftAng(1),2*pi)) (180/pi)*(mod(ang1+(qq2-1)*shiftAng(2),2*pi)) (180/pi)*(mod(ang1+(qq3-1)*shiftAng(3),2*pi))],...
                    1,[],[],[],[],[]); % 0 for nearest neighboring filtering, 1 for bilinear filtering
                if num_rings >=10
                    dstRect10 = [dstRect10a dstRect10b dstRect10c];
                    Screen('DrawTextures',w,eleTexMatch1,[],dstRect10,[(180/pi)*(mod(ang2+(rr1-1)*shiftAng(1),2*pi)) (180/pi)*(mod(ang2+(rr2-1)*shiftAng(2),2*pi)) (180/pi)*(mod(ang2+(rr3-1)*shiftAng(3),2*pi))],...
                        1,[],[],[],[],[]); % 0 for nearest neighboring filtering, 1 for bilinear filtering
                    
                end  %10
            end %9
        end %8
        %              Screen('BlendFunction',w,GL_ONE_MINUS_SRC_ALPHA,GL_SRC_ALPHA,[1 1 1 0]);
        %              Screen('DrawTexture',w,masktex0);
        %              Screen('BlendFunction',w,GL_ONE,GL_ZERO,[1 1 1 1]);
        
        %%%%%%Exp  area3
        i3 = i3+1;
        %         ii3 = i3;
        if i3 == onFrames
            i3= 1;
            %             ii3 = i3;
        end
        if num_rings >=2
            j3 = j3+1;
            %             jj3 = j3;
            if j3 == onFrames
                j3 = 1;
                %                 jj3 = j3;
            end
        end
        if num_rings>=3
            k3 = k3+1;
            %             kk3 = k3;
            if k3 == onFrames
                k3 = 1;
                %                 kk3 = k3;
            end
        end
        if num_rings>=4
            l3 = l3+1;
            %             ll3 = l3;
            if l3 == onFrames
                l3 = 1;
                %                 ll3 = l3;
            end
        end
        if num_rings>=5
            m3 = m3+1;
            %             mm3 = m3;
            if m3 == onFrames
                m3 = 1;
                %                 mm3 = m3;
            end
        end
        if num_rings>=6
            n3 = n3+1;
            %             nn3 = n3;
            if n3 == onFrames
                n3 = 1;
                %                 nn3 = n3;
            end
        end
        if num_rings >=7
            o3 = o3+1;
            %             oo3 = o3;
            if o3 == onFrames
                o3 = 1;
                %                 oo3 = o3;
            end
            if num_rings >=8
                p3 = p3 +1;
                %                 pp3 = p3;
                if p3 == onFrames
                    p3 = 1;
                    %                     pp3 = p3;
                end
                if num_rings >=9
                    q3 = q3 +1;
                    %                     qq3 = q3;
                    if q3 == onFrames
                        q3 = 1;
                        %                         qq3 = q3;
                    end
                    if num_rings >= 10
                        r3 = r3+1;
                        %                         rr3 = r3;
                        if r3 == onFrames
                            r3 = 1;
                            %                             rr3 = r3;
                        end
                        
                    end  %10
                end   %9
            end  %8
        end %7
        
        %%%%Con  area1
        i1 = i1-1;
        if i1 == 0
            i1 = onFrames;
        end
        if num_rings>=2
            j1 = j1-1;
            if j1 == 0
                j1 = onFrames;
            end
            if num_rings >=3
                k1 = k1-1;
                if k1 == 0
                    k1 = onFrames;
                end
                if num_rings >=4
                    l1 = l1-1;
                    if l1 == 0
                        l1 = onFrames;
                    end
                    if num_rings >=5
                        m1 = m1-1;
                        if m1 == 0
                            m1 = onFrames;
                        end
                        if num_rings >=6
                            n1 = n1-1;
                            if n1 == 0
                                n1 = onFrames;
                            end
                            if num_rings >= 7
                                o1 = o1-1;
                                if o1 == 0
                                    o1 = onFrames;
                                end
                                if num_rings >= 8
                                    p1 = p1 -1;
                                    if p1==0
                                        p1 = onFrames;
                                    end
                                    
                                    if num_rings >=9
                                        q1 = q1-1;
                                        if q1 == 0
                                            q1 = onFrames;
                                        end
                                        
                                        if num_rings >=10
                                            r1 = r1-1;
                                            if r1 == 0
                                                r1 =  onFrames;
                                            end
                                            
                                        end %10
                                    end  %9
                                end  %8
                            end  %7
                        end   %6
                    end
                end
            end
        end
        %         ii1 = i1;
        %         jj1 = j1;
        %         kk1 = k1;
        %         ll1 = l1;
        %         mm1 = m1;
        %         nn1 = n1;
        %         oo1 = o1;
        %         pp1 = p1;
        %         qq1 = q1;
        %         rr1 = r1;
        
        ii1 = ii1+1;
        jj1 = jj1+1;
        kk1 = kk1+1;
        ll1 = ll1+1;
        mm1 = mm1+1;
        nn1 = nn1+1;
        oo1 = oo1+1;
        pp1 = pp1+1;
        qq1 = qq1+1;
        rr1 = rr1+1;
        
        ii3 = ii3+1;
        jj3 = jj3+1;
        kk3 = kk3+1;
        ll3 = ll3+1;
        mm3 = mm3+1;
        nn3 = nn3+1;
        oo3 = oo3+1;
        pp3 = pp3+1;
        qq3 = qq3+1;
        rr3 = rr3+1;
        
        %  只为了静止时，真实旋转可以实现
        ii2 = ii2+1;
        jj2 = jj2+1;
        kk2 = kk2+1;
        ll2 = ll2+1;
        mm2 = mm2+1;
        nn2 = nn2+1;
        oo2 = oo2+1;
        pp2 = pp2+1;
        qq2 = qq2+1;
        rr2 = rr2+1;
        
        
        vbl = Screen('Flip',w,vbl+(waitFrames-0.5)*ifi);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        sM.drawCross(0.3,[0 0 0 1]);

        Screen('DrawText',w,'a',wrect(3)/4,wrect(4)/2+200,colorMat(1,:));
        
        Screen('DrawText',w,'b',wrect(3)*2/4,wrect(4)/2+200,colorMat(2,:));
        
        Screen('DrawText',w,'c',wrect(3)*3/4,wrect(4)/2+200,colorMat(3,:));
        
        
        [keyisdown,secs,keycode] = KbCheck(-1);
        
        if keycode(spa)
            numChoice = 1;  %real angle speed ,no CCW and CW
            flag = 1;
						abandon = 0;
        end
        if keycode(spb)
            numChoice = 2;  %real angle speed ,no CCW and CW
            flag = 1;
						abandon = 0;
        end
        if keycode(spc)
            numChoice = 3;  %real angle speed ,no CCW and CW
            flag = 1;
						abandon = 0;
				end
        
        if keycode(ok) && flag == 1
            break;
        end
        
        if keycode(esc)
            breakLoop = true;
            break;
        end
        
    end
     
end

if flag==0
    abandon = 1;   %>10s no select,restart
end

%     i = i+1;
%     if keycode(ok) && flag == 1
%         break;
%     end
%  end
return;