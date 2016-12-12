function Pinna_grating_match(condition,match_time)
global sM w wrect ifi waitFrames ang1 ang2 ...
	spa spb ok esc cKey onFrames approach r1Origin abandon ...
  xc yc r1Match  eleTexMatch1 txtColorMat ...
	numChoice ppd breakLoop

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
shiftAng(1:2) = 0;

approach = 2;

shiftAng(1) = -pi*(20.*ifi)/180;
shiftAng(2) = pi*(20.*ifi)/180;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%expansion
i = 1;  %when only has one ring
ii = 1; %同角速度 只用一个变量即可
j = 1;
k = 1;
l = 1;
m = 1;
n = 1;
if num_rings ==2
    j = round(onFrames/2); %2、
end
if num_rings == 3
    j = round(onFrames/3); %li:入 确保i j k为整数 保证轨迹相同，且能实现i == onFrames
    k = round(2*onFrames/3);
end
if num_rings == 4
    j = round(onFrames/4);
    k = round(2*onFrames/4);
    l = round(3*onFrames/4);
end
if num_rings == 5
    j = round(onFrames/5);
    k = round(2*onFrames/5);
    l = round(3*onFrames/5);
    m = round(4*onFrames/5);
end
if num_rings == 6
    j = round(onFrames/6);
    k = round(2*onFrames/6);
    l = round(3*onFrames/6);
    m = round(4*onFrames/6);
    n = round(5*onFrames/6);
end

vbl = Screen('Flip',w);
vblendtime = vbl + eachConditionSecs;
max_r = r1Match;
min_size = 5;
max_size = 35;
r_a = 2*(max_r-r1Origin)/(onFrames)^2;
size_a = 2*(max_size-min_size)/(onFrames)^2;

while vbl < vblendtime
    
    %%%%%%%%%%1st ring
    %area1
    radius1 = round(r1Origin + 0.5*r_a*i^2);
    side1P(1) = round(min_size + 0.5*size_a*i^2);
    side1P(2) = round(min_size + 0.5*size_a*i^2);
    dstRect1a  = [xc(1)+radius1*sin(mod(ang1+(ii-1)*shiftAng(1),2*pi))-side1P(2)/2;yc(1)-radius1*cos(mod(ang1+(ii-1)*shiftAng(1),2*pi))-side1P(1)/2;...
        xc(1)+radius1*sin(mod(ang1+(ii-1)*shiftAng(1),2*pi))+side1P(2)/2;yc(1)-radius1*cos(mod(ang1+(ii-1)*shiftAng(1),2*pi))+side1P(1)/2];
    %area2
    dstRect1b  = [xc(2)+radius1*sin(mod(ang1+(ii-1)*shiftAng(2),2*pi))-side1P(2)/2;yc(1)-radius1*cos(mod(ang1+(ii-1)*shiftAng(2),2*pi))-side1P(1)/2;...
        xc(2)+radius1*sin(mod(ang1+(ii-1)*shiftAng(2),2*pi))+side1P(2)/2;yc(1)-radius1*cos(mod(ang1+(ii-1)*shiftAng(2),2*pi))+side1P(1)/2];
    
    %%%%%%%%%2 ring
    if num_rings>=2
        %area1
        radius2 = round(r1Origin + 0.5*r_a*j^2);
        side2P(1) = round(min_size + 0.5*size_a*j^2);
        side2P(2) = round(min_size + 0.5*size_a*j^2);
        dstRect2a = [xc(1)+radius2*sin(mod(ang2+(ii-1)*shiftAng(1),2*pi))-side2P(2)/2;yc(1)-radius2*cos(mod(ang2+(ii-1)*shiftAng(1),2*pi))-side2P(1)/2;...
            xc(1)+radius2*sin(mod(ang2+(ii-1)*shiftAng(1),2*pi))+side2P(2)/2;yc(1)-radius2*cos(mod(ang2+(ii-1)*shiftAng(1),2*pi))+side2P(1)/2];
        %area2
        dstRect2b  = [xc(2)+radius2*sin(mod(ang2+(ii-1)*shiftAng(2),2*pi))-side2P(2)/2;yc(1)-radius2*cos(mod(ang2+(ii-1)*shiftAng(2),2*pi))-side2P(1)/2;...
            xc(2)+radius2*sin(mod(ang2+(ii-1)*shiftAng(2),2*pi))+side2P(2)/2;yc(1)-radius2*cos(mod(ang2+(ii-1)*shiftAng(2),2*pi))+side2P(1)/2];
        
    end
    %%%%%%%%%3 ring
    if num_rings >=3
        %area1
        radius3 = round(r1Origin + 0.5*r_a*k^2);
        side3P(1) = round(min_size + 0.5*size_a*k^2);
        side3P(2) = round(min_size + 0.5*size_a*k^2);
        dstRect3a = [xc(1)+radius3*sin(mod(ang1+(ii-1)*shiftAng(1),2*pi))-side3P(2)/2;yc(1)-radius3*cos(mod(ang1+(ii-1)*shiftAng(1),2*pi))-side3P(1)/2;...
            xc(1)+radius3*sin(mod(ang1+(ii-1)*shiftAng(1),2*pi))+side3P(2)/2;yc(1)-radius3*cos(mod(ang1+(ii-1)*shiftAng(1),2*pi))+side3P(1)/2];
        %area2
        dstRect3b  = [xc(2)+radius3*sin(mod(ang1+(ii-1)*shiftAng(2),2*pi))-side3P(2)/2;yc(1)-radius3*cos(mod(ang1+(ii-1)*shiftAng(2),2*pi))-side3P(1)/2;...
            xc(2)+radius3*sin(mod(ang1+(ii-1)*shiftAng(2),2*pi))+side3P(2)/2;yc(1)-radius3*cos(mod(ang1+(ii-1)*shiftAng(2),2*pi))+side3P(1)/2];
    end
    %%%%%%%%%4 ring
    if num_rings>=4
        r4 = round(r1Origin + 0.5*r_a*l^2);
        side4P(1) = round(min_size + 0.5*size_a*l^2);
        side4P(2) = round(min_size + 0.5*size_a*l^2);
        dstRect4a = [xc(1)+r4*sin(mod(ang2+(ii-1)*shiftAng(1),2*pi))-side4P(2)/2;yc(1)-r4*cos(mod(ang2+(ii-1)*shiftAng(1),2*pi))-side4P(1)/2;...
            xc(1)+r4*sin(mod(ang2+(ii-1)*shiftAng(1),2*pi))+side4P(2)/2;yc(1)-r4*cos(mod(ang2+(ii-1)*shiftAng(1),2*pi))+side4P(1)/2];
        %area2
        dstRect4b  = [xc(2)+r4*sin(mod(ang2+(ii-1)*shiftAng(2),2*pi))-side4P(2)/2;yc(1)-r4*cos(mod(ang2+(ii-1)*shiftAng(2),2*pi))-side4P(1)/2;...
            xc(2)+r4*sin(mod(ang2+(ii-1)*shiftAng(2),2*pi))+side4P(2)/2;yc(1)-r4*cos(mod(ang2+(ii-1)*shiftAng(2),2*pi))+side4P(1)/2];
    end
    %%%%%%%%%5 ring
    if num_rings >=5
        r5 = round(r1Origin + 0.5*r_a*m^2);
        side5P(1) = round(min_size + 0.5*size_a*m^2);
        side5P(2) = round(min_size + 0.5*size_a*m^2);
        dstRect5a = [xc(1)+r5*sin(mod(ang1+(ii-1)*shiftAng(1),2*pi))-side5P(2)/2;yc(1)-r5*cos(mod(ang1+(ii-1)*shiftAng(1),2*pi))-side5P(1)/2;...
            xc(1)+r5*sin(mod(ang1+(ii-1)*shiftAng(1),2*pi))+side5P(2)/2;yc(1)-r5*cos(mod(ang1+(ii-1)*shiftAng(1),2*pi))+side5P(1)/2];
        %area2
        dstRect5b  = [xc(2)+r5*sin(mod(ang1+(ii-1)*shiftAng(2),2*pi))-side5P(2)/2;yc(1)-r5*cos(mod(ang1+(ii-1)*shiftAng(2),2*pi))-side5P(1)/2;...
            xc(2)+r5*sin(mod(ang1+(ii-1)*shiftAng(2),2*pi))+side5P(2)/2;yc(1)-r5*cos(mod(ang1+(ii-1)*shiftAng(2),2*pi))+side5P(1)/2];
    end
    %%%%%%%%%6 ring
    if num_rings>=6
        r6 = round(r1Origin + 0.5*r_a*n^2);
        side6P(1) = round(min_size + 0.5*size_a*n^2);
        side6P(2) = round(min_size + 0.5*size_a*n^2);
        dstRect6a = [xc(1)+r6*sin(mod(ang2+(ii-1)*shiftAng(1),2*pi))-side6P(2)/2;yc(1)-r6*cos(mod(ang2+(ii-1)*shiftAng(1),2*pi))-side6P(1)/2;...
            xc(1)+r6*sin(mod(ang2+(ii-1)*shiftAng(1),2*pi))+side6P(2)/2;yc(1)-r6*cos(mod(ang2+(ii-1)*shiftAng(1),2*pi))+side6P(1)/2];
        %area2
        dstRect6b  = [xc(2)+r6*sin(mod(ang2+(ii-1)*shiftAng(2),2*pi))-side6P(2)/2;yc(1)-r6*cos(mod(ang2+(ii-1)*shiftAng(2),2*pi))-side6P(1)/2;...
            xc(2)+r6*sin(mod(ang2+(ii-1)*shiftAng(2),2*pi))+side6P(2)/2;yc(1)-r6*cos(mod(ang2+(ii-1)*shiftAng(2),2*pi))+side6P(1)/2];
	 end
    
	 %Screen('BlendFunction',w,GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA,[1 1 1 1]);
    dstRect1 = [dstRect1a dstRect1b];
    Screen('DrawTextures',w,eleTexMatch1,[],dstRect1,[(180/pi)*(mod(ang1+(ii-1)*shiftAng(1),2*pi)) (180/pi)*(mod(ang1+(ii-1)*shiftAng(2),2*pi))],...
        1,[],[],[],[],[]); % 0 for nearest neighboring filtering, 1 for bilinear filtering
    if num_rings>=2
        dstRect2 = [dstRect2a dstRect2b];
        Screen('DrawTextures',w,eleTexMatch1,[],dstRect2,[(180/pi)*(mod(ang2+(ii-1)*shiftAng(1),2*pi)) (180/pi)*(mod(ang2+(ii-1)*shiftAng(2),2*pi))],...
            1,[],[],[],[],[]); % 0 for nearest neighboring filtering, 1 for bilinear filtering
    end
    if num_rings>=3
        dstRect3 = [dstRect3a dstRect3b];
        Screen('DrawTextures',w,eleTexMatch1,[],dstRect3,[(180/pi)*(mod(ang1+(ii-1)*shiftAng(1),2*pi)) (180/pi)*(mod(ang1+(ii-1)*shiftAng(2),2*pi))],...
            1,[],[],[],[],[]); % 0 for nearest neighboring filtering, 1 for bilinear filtering
    end
    if num_rings>=4
        dstRect4 = [dstRect4a dstRect4b];
        Screen('DrawTextures',w,eleTexMatch1,[],dstRect4,[(180/pi)*(mod(ang2+(ii-1)*shiftAng(1),2*pi)) (180/pi)*(mod(ang2+(ii-1)*shiftAng(2),2*pi))],...
            1,[],[],[],[],[]); % 0 for nearest neighboring filtering, 1 for bilinear filtering
    end
    if num_rings>=5
        dstRect5 = [dstRect5a dstRect5b];
        Screen('DrawTextures',w,eleTexMatch1,[],dstRect5,[(180/pi)*(mod(ang1+(ii-1)*shiftAng(1),2*pi)) (180/pi)*(mod(ang1+(ii-1)*shiftAng(2),2*pi))],...
            1,[],[],[],[],[]); % 0 for nearest neighboring filtering, 1 for bilinear filtering
    end
    if num_rings>=6
        dstRect6 = [dstRect6a dstRect6b];
        Screen('DrawTextures',w,eleTexMatch1,[],dstRect6,[(180/pi)*(mod(ang2+(ii-1)*shiftAng(1),2*pi)) (180/pi)*(mod(ang2+(ii-1)*shiftAng(2),2*pi))],...
            1,[],[],[],[],[]); % 0 for nearest neighboring filtering, 1 for bilinear filtering
	 end
	 %Screen('BlendFunction',w,GL_ONE,GL_ZERO,[1 1 1 1]);
    
    ii = ii+1;
    
		sM.drawCross(0.3,[0 0 0 1]);
    Screen('DrawText',w,'a',wrect(3)/3,wrect(4)/2+200,colorMat(1,:));
    Screen('DrawText',w,'b',wrect(3)*2/3,wrect(4)/2+200,colorMat(2,:));

    vbl = Screen('Flip',w,vbl+(waitFrames-0.5)*ifi);

    [~,~,keycode] = KbCheck(-1);
		
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
    
    if keycode(ok) & flag == 1
       break
    end
    
    if keycode(esc)
        breakLoop = true;
        break
		end
end

if flag==0
    abandon = 1;   %>10s no select,restart
end