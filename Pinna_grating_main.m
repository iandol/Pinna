function Pinna_grating_main(angle_pattern,move_speed_i,angle_speed_i,ResultDir,one_trials,duration,match_time,is_binary_mask,mask_diameter,mask_xpos,mask_ypos)

global sM ana w wrect ifi waitFrames ang1 ang2 white black gray  approach spa spb r1Origin abandon move_speed...
	ok  xc yc ovalRect r1Match  eleTexMatch1  txtColorMat shiftAng1  onFrames numChoice angSpeed1 ppd

PsychDefaultSetup(2);
Screen('Preference', 'SkipSyncTests', 2)
KbName('UnifyKeyNames');
esc = KbName('escape');
if ispc; ok = KbName('return'); else; ok = KbName('enter'); end
spa = KbName('uparrow');
spb = KbName('downarrow');

% viewing parameters ------------------------------------------------------
screenNumber = max(Screen('Screens'));%-1;
pixelsPerCm = 35;
distance = 36.5;
windowed = [];
backgroundColor = [0.5 0.5 0.5];
approach = 1;%[0 1]; % simulate approaching (1) or leaving (0)
% directions = [0 1]; % whether the inner ring has CW (0) or CCW (1) rotational direction when approaching, equivalently CCW (0) or CW (1) when leaving
% speeds = 200;  % translation speed, in mm/sec
% angSpeed = 0; % rotation speeds, in degrees/sec, +: CW, -: CCW
allAngle1 =  angle_pattern;%[45 135]; %angle of inclination  ,absolute value  no + -
num_rings = 10;
is_mask_outer = 1;
is_mask_inner = 1;
f = 0.05; %grating frequency
grating_size = 15; %micropattern size (pixel)
maskinner_radius = 2; %deg
maskouter_radius = 18; %deg
eachConditionSecs = duration; %sec
% number of elements
num1 = 20;

%for 1st ring
offset1 = 0;
ang1 = mod(offset1+linspace(0,360-360/num1,num1),360);  %
ang1 = pi*ang1/180; %degtorad(ang1); %%
offset2 = 180/num1;  % when different rings ,second ring's position is different
ang2 = mod(offset2+linspace(0,360-360/num1,num1),360);  %
ang2 = pi*ang2/180; %degtorad(ang2); %%

% radius of ring of elements
r1 = 0.4; % in degree

% fixation point rectangle
fixSide = 8;% radius of fixation point is 8 pixel
fixColor = 0;  %black


% %%%%%%%%%%%%%%%arrange angle ,10*5 = 50
% AngleArray = [];
numAngle1 = size(allAngle1,2);
num_move_speed = length(move_speed_i);
num_angle_speed = length(angle_speed_i);
if num_move_speed==num_angle_speed
	trials = numAngle1 * one_trials * num_move_speed;  %be sure one by one
else
	msgbox('The number of two kinds of speed is different', 'warning');
end

%%%%randperm
column_rank = randperm(trials);  %rand(50)

%%%%save marker number 1 2 3.. angles and move_speed and angle_speed
pattern_angle_index = ones(1,trials); %value is 1~num_angle
% condition_index = ones(1,trials); %value is 1~num_conditions
rotation_speed_index = ones(1,trials); %value is 1~num_speed
move_speed_index = ones(1,trials);
for i = 1:trials
	%     pattern_angle_index(i) = fix((column_rank(i)-1)/(one_trials * num_move_speed * num_angle_speed))+1; %value is 1~num_angle
	%     rotation_speed_index(i) = fix(rem(column_rank(i)-1,(one_trials * num_move_speed * num_angle_speed))/(one_trials * num_move_speed)) + 1;   %value is 1~num_rotation_speed
	%     move_speed_index(i) = rem(column_rank(i)-1,num_move_speed)+1;  %value is 1~num_move_speed
	pattern_angle_index(i) = fix((column_rank(i)-1)/(one_trials * num_move_speed))+1; %value is 1~num_angle
	rotation_speed_index(i) = fix(rem(column_rank(i)-1,(one_trials * num_move_speed))/(one_trials)) + 1;   %value is 1~num_rotation_speed
	move_speed_index = rotation_speed_index;%value is 1~num_move_speed And num_rotation_speed = num_move_speed; because one-to-one
end
%------------
%global struct ana to save para
ana.angle_index = pattern_angle_index;
ana.move_speed_index = move_speed_index;
ana.angle_speed_index = rotation_speed_index;
ana.randArray = column_rank;
ana.result = zeros(1,trials);
ana.allAngle = allAngle1;
ana.one_trials = one_trials;
ana.moving_speed = move_speed_i;
ana.angle_speed = angle_speed_i;
ana.sti_duration = duration;
ana.sti_match_time = match_time;
ana.is_binary_mask = is_binary_mask;
ana.mask_diameter = mask_diameter;
ana.mask_xpos = mask_xpos;
ana.mask_ypos = mask_ypos;
%--------------------------------------------------------------------------
try
	%-----------------------open the PTB screens------------------------
	%Screen('BlendFunction', w, GL_ONE, GL_ONE);
	sM = screenManager('verbose',false,'blend',true,'screen',screenNumber,...
		'pixelsPerCm',pixelsPerCm,...
		'distance',distance,'bitDepth','FloatingPoint32BitIfPossible',...
		'debug',false,'antiAlias',0,'nativeBeamPosition',0, ...
		'srcMode','GL_SRC_ALPHA','dstMode','GL_ONE_MINUS_SRC_ALPHA',...
		'windowed',windowed,'backgroundColour',[backgroundColor 0],...
		'gammaTable', calibrationFile); %use a temporary screenManager object
	screenVals = open(sM); %open PTB screen
	w = sM.win;
	wrect = sM.winRect;
	ifi = sM.screenVals.ifi;
	ppd = sM.ppd;
	
	% screen
	white = sM.screenVals.white;
	black = sM.screenVals.black;
	gray = (white+black)/2; % index for white, black and gray
	if gray == white
		gray = white/2;
	end
	inc = white-gray;
	
	% screen areas
	xc = (1:2).*wrect(3)/3; %match -only have A and B
	yc = wrect(4)/2;
	r1Match = round(45.0/sizePixel);%30mm/0.25=120pixels ;max radius at match
	
	%procedural gabor
	sizeGabor = 301;
	phase = 0;
	sc = 10.0;
	freq = .05;
	contrast = 10;
	aspectratio = 1.0;
	%mypars = repmat([phase+180, freq, sc, contrast, aspectratio, 0, 0, 0]', 1, ngabors);
	gabortex = CreateProceduralGabor(w, sizeGabor, sizeGabor, 1);
	
	%%%%%%% 3 fixPoint
	fixSideMatch = round(1.7/sizePixel);%6;
	ovalRect = zeros(4,2);
	ovalRect(1,1:2) = xc-0.5*fixSideMatch;
	ovalRect(2,1:2) = yc-0.5*fixSideMatch;
	ovalRect(3,1:2) = xc+0.5*fixSideMatch;
	ovalRect(4,1:2) = yc+0.5*fixSideMatch;
	
	r1Origin = r1 * ppd;% convert degrees to pixels
	
	if is_binary_mask
		mask_diameter = mask_diameter * ppd;
		radius = mask_diameter/2; %diameter to radius
		masksize = [wrect(4),wrect(3)];
		mask = ones(masksize(1),masksize(2),2)*gray;
		[x1,y1] = meshgrid(0:masksize(2)-1,0:masksize(1)-1);
		mask(:,:,2) = white*((x1-mask_xpos).^2 + (y1-mask_ypos).^2 <=(radius).^2);
		masktex11 = Screen('MakeTexture', w, mask);
	end
	
	max_r = 800;
	min_size = 12;
	max_size = 120;
	r_dis = max_r-r1Origin;
	
	width_line = 5; %width of line
	m = ones(61,61)*128; %line length = 50
	m(:,28:32) = 0;
	eleTexMatch1 = Screen('MakeTexture',w,m);
	
	% color matrix for letters
	txtColorMat = 255.*ones(10,3,11);
	for i = 1:10
		txtColorMat(i,2:3,i+1) = 0;
	end
	%----------------------------------------------------------------------
	AllElements = [];
	%make grating texture
	f=f*2*pi; %used for angle1 and angle2
	for jj = 1:numAngle1   %number of the MakeTexture is 20 ,at least
		[x,y]=meshgrid(-30:30,-30:30);  %used for angle1 and angle2
		angle1=(-allAngle1(jj)-90)*pi/180;
		a1=cos(angle1)*f;
		b1=sin(angle1)*f;
		m1=exp(-((x/grating_size).^2)-((y/grating_size).^2)).*sin(a1*x+b1*y-5);
		AllElements(1,jj).eleTex1P = Screen('MakeTexture',w,gray+inc*m1);    %CCW
		
	end
	%%%%%%%%%%%%%%%%%%%%%%outside mask & inside mask
	
	maskinner_radius = maskinner_radius * ppd;% convert degrees to pixels
	maskouter_radius = maskouter_radius * ppd;% convert degrees to pixels
	
	if is_mask_outer == 1
		masksize = [wrect(4)/4,wrect(3)/4];
		mask = ones(masksize(1),masksize(2),2)*gray;
		[x,y] = meshgrid(0:masksize(2)-1,0:masksize(1)-1);
		mask(:, :, 2)=white *(exp(-((x-xCen).^2+(y-yCen).^2)/(maskouter_radius)^2)); % gaussian mask
		mask(mask>white)=white;
		masktex = Screen('MakeTexture', w, mask);
	end
	if is_mask_inner == 1
		masksize = [wrect(4)/4,wrect(3)/4];
		mask2 = ones(masksize(1),masksize(2),2)*gray;
		[x,y] = meshgrid(0:masksize(2)-1,0:masksize(1)-1);
		mask2(:, :, 2)=white * (1-exp(-((x-xCen).^2+(y-yCen).^2)/((maskinner_radius)*2)^2)); % gaussian mask
		mask2(mask2>white)=white;
		mask2(mask2<30) = black;  %
		masktex2 = Screen('MakeTexture', w, mask2);
	end
	%----------------------
	
	HideCursor;
	Screen('FillRect',w,gray,[]);
	normBoundsRect = Screen('TextBounds', w, 'Please fix at the central point.');
	Screen('DrawText',w,'Please fix at the central point.',xCen-normBoundsRect(3)/2,yCen-normBoundsRect(4)/2);
	Screen('Flip',w);
	WaitSecs(2);
	
	Screen('FillOval',w,fixColor,[xCen-fixSide/2,yCen-fixSide/2,xCen+fixSide/2,yCen+fixSide/2]);
	Screen('Flip',w);
	WaitSecs(1);
	iii = 1;
	breakLoop = false;
	%     for iii = 1:trials
	while ~breakLoop
		if iii<=trials
			waitFrames = 1;
			angSpeed1 = angle_speed_i(rotation_speed_index(iii));
			shiftAng1 = pi*(angSpeed1.*ifi)/180;  %degtorad(angSpeed1.*ifi);
			move_speed = move_speed_i(move_speed_index(iii)) * ppd;
			% radius of ring of elements
			if move_speed>0 %Exp And Exp+CW And Exp+CCW
				speeds = move_speed;
				onFrames = fix((r_dis + sqrt(r_dis*r_dis + 4*r_dis*speeds*ifi))/(2*speeds*ifi));
				approach = 1;
				shiftAng = shiftAng1;%shiftAng1>0isCW ;shiftAng1<0 is CCW;shiftAng1=0 is no rotation
				condition = 1;
			end
			if move_speed<0 %Con And Con+CW And Con+CCW
				speeds = -move_speed;
				onFrames = fix((r_dis + sqrt(r_dis*r_dis + 4*r_dis*speeds*ifi))/(2*speeds*ifi));
				approach = 0;
				shiftAng = shiftAng1;
				condition = 2;
			end
			if angSpeed1>0&&move_speed==0 %CW
				speeds = 5*ppd;
				onFrames = fix((r_dis + sqrt(r_dis*r_dis + 4*r_dis*speeds*ifi))/(2*speeds*ifi));
				approach = 2;  %when speed = 0,  approach = 2；
				shiftAng = shiftAng1; %degtorad(angSpeed*ifi); % degree of rotation per frame  %%真实旋转
				condition = 3;
			end
			if angSpeed1<0&&move_speed==0  %CCW
				speeds = 5*ppd;
				onFrames = fix((r_dis + sqrt(r_dis*r_dis + 4*r_dis*speeds*ifi))/(2*speeds*ifi));
				approach = 2;  %when speed = 0,  approach = 2；
				shiftAng = shiftAng1; %degtorad(angSpeed*ifi); % degree of rotation per frame  %%真实旋转
				condition = 4;
			end
			
			%%%%%Part3 is rotation_speed!=0 && move_speed_index !=0 的情况
			
			row = 1;
			column = pattern_angle_index(iii);
			
			%------------------------
			
			i = 1;  %when only has one ring
			ii = 1;
			jj = 1;
			kk = 1;
			ll = 1;
			mm = 1;
			nn = 1;
			oo = 1;
			pp = 1;
			qq = 1;
			rr = 1;
			ss = 1;
			tt = 1;
			j = 1;
			k = 1;
			l = 1;
			m = 1;
			n = 1;
			o = 1;
			p = 1;
			q = 1;
			r = 1;
			s = 1;
			t = 1;
			
			if num_rings == 8
				j = round(onFrames/8);
				k = round(2*onFrames/8);
				l = round(3*onFrames/8);
				m = round(4*onFrames/8);
				n = round(5*onFrames/8);
				o = round(6*onFrames/8);
				p = round(7*onFrames/8);
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
				
			end
			if num_rings == 11
				j = round(onFrames/11);
				k = round(2*onFrames/11);
				l = round(3*onFrames/11);
				m = round(4*onFrames/11);
				n = round(5*onFrames/11);
				o = round(6*onFrames/11);
				p = round(7*onFrames/11);
				q = round(8*onFrames/11);
				r = round(9*onFrames/11);
				s = round(10*onFrames/11);
				
			end
			
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			if approach ==2   %%每个trial ,ii ,jj ……初始化
				ii = 1;
				jj = 1;
				kk = 1;
				ll = 1;
				mm = 1;
				nn = 1;
				oo = 1;
				pp = 1;
				qq = 1;
				rr = 1;
				ss = 1;
				tt = 1;
				
			end
			
			vbl = Screen('Flip',w);
			vblendtime = vbl + eachConditionSecs;
			r_a = 2*(max_r-r1Origin)/(onFrames)^2;
			size_a = 2*(max_size-min_size)/(onFrames)^2;
			
			while vbl < vblendtime
				
				r1 = round(r1Origin + 0.5*r_a*i^2);
				side1P(1) = round(min_size + 0.5*size_a*i^2);
				side1P(2) = round(min_size + 0.5*size_a*i^2);
				dstRect1 = [xCen+r1*sin(mod(ang1+(ii-1)*shiftAng,2*pi))-side1P(2)/2;yCen-r1*cos(mod(ang1+(ii-1)*shiftAng,2*pi))-side1P(1)/2;...
					xCen+r1*sin(mod(ang1+(ii-1)*shiftAng,2*pi))+side1P(2)/2;yCen-r1*cos(mod(ang1+(ii-1)*shiftAng,2*pi))+side1P(1)/2];
				
				r2 = round(r1Origin + 0.5*r_a*j^2);
				side2P(1) = round(min_size + 0.5*size_a*j^2);
				side2P(2) = round(min_size + 0.5*size_a*j^2);
				dstRect2 = [xCen+r2*sin(mod(ang2+(jj-1)*shiftAng,2*pi))-side2P(2)/2;yCen-r2*cos(mod(ang2+(jj-1)*shiftAng,2*pi))-side2P(1)/2;...
					xCen+r2*sin(mod(ang2+(jj-1)*shiftAng,2*pi))+side2P(2)/2;yCen-r2*cos(mod(ang2+(jj-1)*shiftAng,2*pi))+side2P(1)/2];
				
				
				r3 = round(r1Origin + 0.5*r_a*k^2);
				side3P(1) = round(min_size + 0.5*size_a*k^2);
				side3P(2) = round(min_size + 0.5*size_a*k^2);
				dstRect3 = [xCen+r3*sin(mod(ang1+(kk-1)*shiftAng,2*pi))-side3P(2)/2;yCen-r3*cos(mod(ang1+(kk-1)*shiftAng,2*pi))-side3P(1)/2;...
					xCen+r3*sin(mod(ang1+(kk-1)*shiftAng,2*pi))+side3P(2)/2;yCen-r3*cos(mod(ang1+(kk-1)*shiftAng,2*pi))+side3P(1)/2];
				
				r4 = round(r1Origin + 0.5*r_a*l^2);
				side4P(1) = round(min_size + 0.5*size_a*l^2);
				side4P(2) = round(min_size + 0.5*size_a*l^2);
				dstRect4 = [xCen+r4*sin(mod(ang2+(ll-1)*shiftAng,2*pi))-side4P(2)/2;yCen-r4*cos(mod(ang2+(ll-1)*shiftAng,2*pi))-side4P(1)/2;...
					xCen+r4*sin(mod(ang2+(ll-1)*shiftAng,2*pi))+side4P(2)/2;yCen-r4*cos(mod(ang2+(ll-1)*shiftAng,2*pi))+side4P(1)/2];
				
				r5 = round(r1Origin + 0.5*r_a*m^2);
				side5P(1) = round(min_size + 0.5*size_a*m^2);
				side5P(2) = round(min_size + 0.5*size_a*m^2);
				dstRect5 = [xCen+r5*sin(mod(ang1+(mm-1)*shiftAng,2*pi))-side5P(2)/2;yCen-r5*cos(mod(ang1+(mm-1)*shiftAng,2*pi))-side5P(1)/2;...
					xCen+r5*sin(mod(ang1+(mm-1)*shiftAng,2*pi))+side5P(2)/2;yCen-r5*cos(mod(ang1+(mm-1)*shiftAng,2*pi))+side5P(1)/2];
				
				r6 = round(r1Origin + 0.5*r_a*n^2);
				side6P(1) = round(min_size + 0.5*size_a*n^2);
				side6P(2) = round(min_size + 0.5*size_a*n^2);
				dstRect6 = [xCen+r6*sin(mod(ang2+(nn-1)*shiftAng,2*pi))-side6P(2)/2;yCen-r6*cos(mod(ang2+(nn-1)*shiftAng,2*pi))-side6P(1)/2;...
					xCen+r6*sin(mod(ang2+(nn-1)*shiftAng,2*pi))+side6P(2)/2;yCen-r6*cos(mod(ang2+(nn-1)*shiftAng,2*pi))+side6P(1)/2];
				
				r7 = round(r1Origin + 0.5*r_a*o^2);
				side7P(1) = round(min_size + 0.5*size_a*o^2);
				side7P(2) = round(min_size + 0.5*size_a*o^2);
				dstRect7 = [xCen+r7*sin(mod(ang1+(oo-1)*shiftAng,2*pi))-side7P(2)/2;yCen-r7*cos(mod(ang1+(oo-1)*shiftAng,2*pi))-side7P(1)/2;...
					xCen+r7*sin(mod(ang1+(oo-1)*shiftAng,2*pi))+side7P(2)/2;yCen-r7*cos(mod(ang1+(oo-1)*shiftAng,2*pi))+side7P(1)/2];
				if num_rings >=8
					r8 = round(r1Origin + 0.5*r_a*p^2);
					side8P(1) = round(min_size + 0.5*size_a*p^2);
					side8P(2) = round(min_size + 0.5*size_a*p^2);
					dstRect8 = [xCen+r8*sin(mod(ang2+(pp-1)*shiftAng,2*pi))-side8P(2)/2;yCen-r8*cos(mod(ang2+(pp-1)*shiftAng,2*pi))-side8P(1)/2;...
						xCen+r8*sin(mod(ang2+(pp-1)*shiftAng,2*pi))+side8P(2)/2;yCen-r8*cos(mod(ang2+(pp-1)*shiftAng,2*pi))+side8P(1)/2];
					if num_rings >=9
						r9 = round(r1Origin + 0.5*r_a*q^2);
						side9P(1) = round(min_size + 0.5*size_a*q^2);
						side9P(2) = round(min_size + 0.5*size_a*q^2);
						dstRect9 = [xCen+r9*sin(mod(ang1+(qq-1)*shiftAng,2*pi))-side9P(2)/2;yCen-r9*cos(mod(ang1+(qq-1)*shiftAng,2*pi))-side9P(1)/2;...
							xCen+r9*sin(mod(ang1+(qq-1)*shiftAng,2*pi))+side9P(2)/2;yCen-r9*cos(mod(ang1+(qq-1)*shiftAng,2*pi))+side9P(1)/2];
						
						if num_rings>= 10
							r10 = round(r1Origin + 0.5*r_a*r^2);
							side10P(1) = round(min_size + 0.5*size_a*r^2);
							side10P(2) = round(min_size + 0.5*size_a*r^2);
							dstRect10 = [xCen+r10*sin(mod(ang2+(rr-1)*shiftAng,2*pi))-side10P(2)/2;yCen-r10*cos(mod(ang2+(rr-1)*shiftAng,2*pi))-side10P(1)/2;...
								xCen+r10*sin(mod(ang2+(rr-1)*shiftAng,2*pi))+side10P(2)/2;yCen-r10*cos(mod(ang2+(rr-1)*shiftAng,2*pi))+side10P(1)/2];
							if num_rings >=11
								r11 = round(r1Origin + 0.5*r_a*s^2);
								side11P(1) = round(min_size + 0.5*size_a*s^2);
								side11P(2) = round(min_size + 0.5*size_a*s^2);
								dstRect11 = [xCen+r11*sin(mod(ang1+(ss-1)*shiftAng,2*pi))-side11P(2)/2;yCen-r11*cos(mod(ang1+(ss-1)*shiftAng,2*pi))-side11P(1)/2;...
									xCen+r11*sin(mod(ang1+(ss-1)*shiftAng,2*pi))+side11P(2)/2;yCen-r11*cos(mod(ang1+(ss-1)*shiftAng,2*pi))+side11P(1)/2];
								
							end  %11
						end %10
					end  %9
				end   %8
				
				mypars = repmat([phase+180, freq, sc, contrast, aspectratio, 0, 0, 0]', 1, length(ang1));
				
				Screen('DrawTextures', w, gabortex, [], dstRect1, 180*(mod(ang1+(ii-1)*shiftAng,2*pi))/pi, [], [], [], [], kPsychDontDoRotation, mypars);
				%Screen('DrawTextures',w,AllElements(row,column).eleTex1P,[],dstRect1,180*(mod(ang1+(ii-1)*shiftAng,2*pi))/pi,1,[],[],[],[],[]); % 0 for nearest neighboring filtering, 1 for bilinear filtering
% 				Screen('DrawTextures',w,AllElements(row,column).eleTex1P,[],dstRect2,180*(mod(ang2+(jj-1)*shiftAng,2*pi))/pi,1,[],[],[],[],[]); % 0 for nearest neighboring filtering, 1 for bilinear filtering
% 				Screen('DrawTextures',w,AllElements(row,column).eleTex1P,[],dstRect3,180*(mod(ang1+(kk-1)*shiftAng,2*pi))/pi,1,[],[],[],[],[]); % 0 for nearest neighboring filtering, 1 for bilinear filtering
% 				Screen('DrawTextures',w,AllElements(row,column).eleTex1P,[],dstRect4,180*(mod(ang2+(ll-1)*shiftAng,2*pi))/pi,1,[],[],[],[],[]); % 0 for nearest neighboring filtering, 1 for bilinear filtering
% 				Screen('DrawTextures',w,AllElements(row,column).eleTex1P,[],dstRect5,180*(mod(ang1+(mm-1)*shiftAng,2*pi))/pi,1,[],[],[],[],[]); % 0 for nearest neighboring filtering, 1 for bilinear filtering
% 				Screen('DrawTextures',w,AllElements(row,column).eleTex1P,[],dstRect6,180*(mod(ang2+(nn-1)*shiftAng,2*pi))/pi,1,[],[],[],[],[]); % 0 for nearest neighboring filtering, 1 for bilinear filtering
% 				Screen('DrawTextures',w,AllElements(row,column).eleTex1P,[],dstRect7,180*(mod(ang1+(oo-1)*shiftAng,2*pi))/pi,1,[],[],[],[],[]); % 0 for nearest neighboring filtering, 1 for bilinear filtering
% 				if num_rings >=8
% 					Screen('DrawTextures',w,AllElements(row,column).eleTex1P,[],dstRect8,180*(mod(ang2+(pp-1)*shiftAng,2*pi))/pi,1,[],[],[],[],[]); % 0 for nearest neighboring filtering, 1 for bilinear filtering
% 					if num_rings >=9
% 						Screen('DrawTextures',w,AllElements(row,column).eleTex1P,[],dstRect9,180*(mod(ang1+(qq-1)*shiftAng,2*pi))/pi,1,[],[],[],[],[]); % 0 for nearest neighboring filtering, 1 for bilinear filtering
% 						if num_rings >=10
% 							Screen('DrawTextures',w,AllElements(row,column).eleTex1P,[],dstRect10,180*(mod(ang2+(rr-1)*shiftAng,2*pi))/pi,1,[],[],[],[],[]); % 0 for nearest neighboring filtering, 1 for bilinear filtering
% 							if num_rings >=11
% 								Screen('DrawTextures',w,AllElements(row,column).eleTex1P,[],dstRect11,180*(mod(ang1+(ss-1)*shiftAng,2*pi))/pi,1,[],[],[],[],[]); % 0 for nearest neighboring filtering, 1 for bilinear filtering
% 							end %11
% 						end  %10
% 					end %9
% 				end %8
				
% 				if is_mask_outer == 1
% 					Screen('BlendFunction',w,GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA,[1 1 1 0]);
% 					Screen('DrawTexture',w,masktex);
% 					Screen('BlendFunction',w,GL_ONE,GL_ZERO,[1 1 1 1]);
% 				end
% 				if is_mask_inner == 1
% 					Screen('BlendFunction',w,GL_ONE_MINUS_SRC_ALPHA,GL_SRC_ALPHA,[1 1 1 0]);
% 					Screen('DrawTexture',w,masktex2);
% 					Screen('BlendFunction',w,GL_ONE,GL_ZERO,[1 1 1 1]);
% 				end
% 				if is_binary_mask
% 					Screen('BlendFunction',w,GL_ONE_MINUS_SRC_ALPHA,GL_SRC_ALPHA,[1 1 1 0]);
% 					Screen('DrawTexture',w,masktex11);
% 					Screen('BlendFunction',w,GL_ONE,GL_ZERO,[1 1 1 1]);
% 				end
				Screen('FillOval',w,fixColor,[xCen-0.5*fixSide,yCen-0.5*fixSide,xCen+0.5*fixSide,yCen+0.5*fixSide]);
				Screen('DrawingFinished', w);
				
				if approach ==1
					i = i+1;
					if i == onFrames
						i = 1;
					end
					
					j = j+1;
					if j == onFrames
						j = 1;
					end
					
					k = k+1;
					if k == onFrames
						k = 1;
					end
					
					l = l+1;
					if l == onFrames
						l = 1;
					end
					
					m = m+1;
					if m == onFrames
						m = 1;
					end
					n = n+1;
					if n == onFrames
						n = 1;
					end
					o = o+1;
					if o == onFrames
						o = 1;
					end
					if num_rings >=8
						p = p +1;
						if p == onFrames
							p = 1;
						end
						if num_rings >=9
							q = q +1;
							if q == onFrames
								q = 1;
							end
							if num_rings >= 10
								r = r+1;
								if r == onFrames
									r = 1;
								end
								
								if num_rings >=11
									s = s+1;
									if s == onFrames
										s = 1;
									end
									
									if num_rings >=12
										t = t+1;
										if t == onFrames
											t = 1;
										end
									end %12
								end %11
							end  %10
						end   %9
					end  %8
					
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
														
														if num_rings >=11
															s = s-1;
															if s == 0
																s = onFrames;
															end
															
															if num_rings >=12
																t = t-1;
																if t == 0
																	t = onFrames;
																end
															end %12
														end %11
													end %10
												end  %9
											end  %8
										end  %7
									end   %6
								end
							end
						end
					end
					
				end
				
				%                 if approach == 2   %  只为了静止时，真实旋转可以实现
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
				ss = ss +1;
				tt = tt +1;
				%                 end
				vbl = Screen('Flip',w,vbl+(waitFrames-0.5)*ifi);
			end
			
			[~,~,keycode] = KbCheck(-1);
			if keycode(esc)
				breakLoop = true;
				break;  %break while 1
			end
			
			Pinna_grating_match(condition,match_time);
			if abandon==0   %else restart,and don't save data
				ana.result(1,iii) = numChoice;
				iii = iii+1;
			end
			
			Screen('FillRect',w,gray,[]);
			Screen('FillOval',w,fixColor,[xCen-fixSide/2,yCen-fixSide/2,xCen+fixSide/2,yCen+fixSide/2]);
			Screen('Flip',w);
			WaitSecs(1); % wait for 1 seconds
			
		else
			break;
		end
	end
	
	save([ResultDir '\' 'data'],'ana');
	Screen('CloseAll');
	ShowCursor;
catch ME
	Screen('CloseAll');sca
	Priority(0);
	psychrethrow(psychlasterror);
	ple(ME)
end

return;