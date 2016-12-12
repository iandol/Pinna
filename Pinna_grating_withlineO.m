function Pinna_grating_withlineO(is_with_line, benchmark)
if ~exist('is_with_line','var'); is_with_line=false; end
if ~exist('benchmark','var'); benchmark=false; end

PsychDefaultSetup(0);
KbName('UnifyKeyNames');
esc = KbName('escape');
Screen('Preference', 'SkipSyncTests', 0)
% screen
screenID = max(Screen('Screens'));

% viewing parameters ------------------------------------------------------
stdDis = 600; %mm
% one_trials = 5;
% approach = ones(trials); %[0 1]; % simulate approaching (1) or leaving (0)
approach = 1;
% directions = [0 1]; % whether the inner ring has CW (0) or CCW (1) rotational direction when approaching, equivalently CCW (0) or CW (1) when leaving
speeds = 200;  % translation speed, in mm/sec
angSpeed = 0; % rotation speeds, in degrees/sec, +: CW, -: CCW
allAngle1 =  45; %angle of inclination  ,absolute value  no + -
num_rings = 10;
is_mask_outer = 1;
is_mask_inner = 0;
f = 0.05;
grating_size = 15;
maskinner_radius = 2;
maskouter_radius = 10;
eachConditionSecs = 3; %1s
% number of elements
num1 = 20;

%for 1st ring
offset1 = 0;
ang1 = mod(offset1+linspace(0,360-360/num1,num1),360);  %
ang1 = degtorad(ang1); %%
offset2 = 180/num1;  % when different rings ,second ring's position is different
ang2 = mod(offset2+linspace(0,360-360/num1,num1),360);  %
ang2 = degtorad(ang2); %%

% radius of ring of elements
r1 = 0.2; % in degree  %% 初始ring的半径

% fixation point rectangle
fixSide = 8;% radius of fixation point is 8 pixel
fixColor = 0;  %black

sizePixel = 0.25;   %in mm, calculated from different Screen  ,has different value

%--------------------------------------------------------------------------
try
	white = WhiteIndex(screenID);
	black = BlackIndex(screenID);
	gray = (white+black)/2; % index for white, black and gray
	if gray == white
		gray = white/2;
	end
	inc = white-gray;
	[w,wrect] = Screen('OpenWindow',screenID,gray,[]);
	ifi = Screen('GetFlipInterval',w);
	halfifi = ifi/2;
	xCen = wrect(3)/2;
	yCen = wrect(4)/2;
	ScreenWidth = round(wrect(3)*sizePixel);
	ppd = wrect(3)/2/atand(ScreenWidth/2/stdDis);
	
	%make grating texture
	f=f*2*pi; %used for angle1 and angle2
	for jj = 1:1   %number of the MakeTexture is 20 ,at least
		[x,y]=meshgrid(-30:30,-30:30);  %used for angle1 and angle2
		angle1=(-allAngle1-90)*pi/180;
		a1=cos(angle1)*f;
		b1=sin(angle1)*f;
		m1=exp(-((x/grating_size).^2)-((y/grating_size).^2)).*sin(a1*x+b1*y-5);
		eleTex1P = Screen('MakeTexture',w,gray+inc*m1);    %CCW
		
		%         angle2=(allAngle1(jj)-90)*pi/180;   %
		%         a2=cos(angle2)*f;
		%         b2=sin(angle2)*f;
		%         m2=exp(-((x/grating_size).^2)-((y/grating_size).^2)).*sin(a2*x+b2*y-5);
		%         AllElements(2,jj).eleTex1P = Screen('MakeTexture',w,gray+inc*m2);    %CW
		
	end
	%%%%%%%%%%%%%%%%%%%%%%outside mask & inside mask
	
	maskinner_radius = maskinner_radius * ppd;% convert degrees to pixels
	maskouter_radius = maskouter_radius * ppd;% convert degrees to pixels
	
	if is_mask_outer == 1
		masksize = [wrect(4),wrect(3)];
		mask = ones(masksize(1),masksize(2),2)*gray;
		[x,y] = meshgrid(0:masksize(2)-1,0:masksize(1)-1);
		mask(:, :, 2)=white *(exp(-((x-xCen).^2+(y-yCen).^2)/(maskouter_radius)^2)); % gaussian mask
		mask(mask>white)=white;
		masktex = Screen('MakeTexture', w, mask);
	end
	if is_mask_inner == 1
		masksize = [wrect(4),wrect(3)];
		mask2 = ones(masksize(1),masksize(2),2)*gray;
		[x,y] = meshgrid(0:masksize(2)-1,0:masksize(1)-1);
		mask2(:, :, 2)=white * (1-exp(-((x-xCen).^2+(y-yCen).^2)/((maskinner_radius)*2)^2)); % gaussian mask
		mask2(mask2>white)=white;
		mask2(mask2<30) = black;  %
		masktex2 = Screen('MakeTexture', w, mask2);
	end
	%----------------------
	
	r1Origin = r1 * ppd;% convert degrees to pixels
	
	
	% radius of ring of elements
	onFrames = fix(stdDis/(speeds*ifi) +1);  %li:为了必须超出screen范围，取最大的帧数，同时不能超过半径公式的极大值， controllde through d0s and speeds
	waitFrames = 1;
	shiftAng = degtorad(angSpeed*ifi); % degree of rotation per frame  %%真实旋转
	
	if speeds == 0
		onFrames = 200;  %
		approach = 2;  %when speed = 0,  approach = 2；
	end
	
	HideCursor;
	Screen('FillRect',w,gray,[]);
	%     normBoundsRect = Screen('TextBounds', w, 'Please fix at the central point.');
	%     Screen('DrawText',w,'Please fix at the central point.',xCen-normBoundsRect(3)/2,yCen-normBoundsRect(4)/2);
	%     Screen('Flip',w);
	%     WaitSecs(2);
	
	Screen('FillOval',w,fixColor,[xCen-fixSide/2,yCen-fixSide/2,xCen+fixSide/2,yCen+fixSide/2]);
	Screen('Flip',w);
	WaitSecs(1);
	
	breakLoop = false;
	
	while ~breakLoop
		
		%      %%%%%%%judge row and column ,get element
		%      if randAngleArray(iii)>= 0   %CCW
		%         row = 1;
		%      else
		%         row = 2;
		%      end
		%      column = find(allAngle1 == abs(randAngleArray(iii)));
		
		%------------------------
		
		i = 1;  %when only has one ring
		ii = 1;   %%%%%%%%%%%%%%%%%%%ii,jj,kk,ll,mm,nn只针对真实旋转的而递增的变量，除了了speed=0 时 与i,j,k,l,m,n一致，且speed = 0时 一直递增！！！
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
			jj = j;
			kk = k;
			ll = l;
			mm = m;
			nn = n;
			oo = o;
			pp = p;
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
			jj = j;
			kk = k;
			ll = l;
			mm = m;
			nn = n;
			oo = o;
			pp = p;
			qq = q;
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
			jj = j;
			kk = k;
			ll = l;
			mm = m;
			nn = n;
			oo = o;
			pp = p;
			qq = q;
			rr = r;
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
			jj = j;
			kk = k;
			ll = l;
			mm = m;
			nn = n;
			oo = o;
			pp = p;
			qq = q;
			rr = r;
			ss = s;
		end
		
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		if approach ==2   %%每个trial ,ii ,jj ……初始化
			ii = i;
			jj = j;
			kk = k;
			ll = l;
			mm = m;
			nn = n;
			oo = o;
			pp = p;
			qq = q;
			rr = r;
			ss = s;
			tt = t;
		end
		
		vbl = Screen('Flip',w);
		vblendtime = vbl + eachConditionSecs;
		ts = vbl;
		count = 0;
		
		max_r = 800;
		min_size = 12;
		max_size = 120;
		r_a = 2*(max_r-r1Origin)/(onFrames)^2;
		size_a = 2*(max_size-min_size)/(onFrames)^2;
		
		%%%%%%
		xy1 = zeros(2,2*num1);
		xy2 = zeros(2,2*num1);
		xy1(1,2:2:end) = 800 * sin(ang1);
		xy1(2,2:2:end) = -800 * cos(ang1);
		xy2(1,2:2:end) = 800 * sin(ang2);
		xy2(2,2:2:end) = -800 * cos(ang2);
		
		%-----------
		while GetSecs < vblendtime
			
			count = count + 1;
			
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
			
			%%%%%%withline
			if is_with_line
				Screen('Drawlines',w,xy1,1,200,[xCen yCen],[]);
				Screen('Drawlines',w,xy2,1,200,[xCen yCen],[]);
			end
			
			Screen('DrawTextures',w,eleTex1P,[],dstRect1,radtodeg(mod(ang1+(ii-1)*shiftAng,2*pi)),1,[],[],[],[],[]); % 0 for nearest neighboring filtering, 1 for bilinear filtering
			Screen('DrawTextures',w,eleTex1P,[],dstRect2,radtodeg(mod(ang2+(jj-1)*shiftAng,2*pi)),1,[],[],[],[],[]); % 0 for nearest neighboring filtering, 1 for bilinear filtering
			Screen('DrawTextures',w,eleTex1P,[],dstRect3,radtodeg(mod(ang1+(kk-1)*shiftAng,2*pi)),1,[],[],[],[],[]); % 0 for nearest neighboring filtering, 1 for bilinear filtering
			Screen('DrawTextures',w,eleTex1P,[],dstRect4,radtodeg(mod(ang2+(ll-1)*shiftAng,2*pi)),1,[],[],[],[],[]); % 0 for nearest neighboring filtering, 1 for bilinear filtering
			Screen('DrawTextures',w,eleTex1P,[],dstRect5,radtodeg(mod(ang1+(mm-1)*shiftAng,2*pi)),1,[],[],[],[],[]); % 0 for nearest neighboring filtering, 1 for bilinear filtering
			Screen('DrawTextures',w,eleTex1P,[],dstRect6,radtodeg(mod(ang2+(nn-1)*shiftAng,2*pi)),1,[],[],[],[],[]); % 0 for nearest neighboring filtering, 1 for bilinear filtering
			Screen('DrawTextures',w,eleTex1P,[],dstRect7,radtodeg(mod(ang1+(oo-1)*shiftAng,2*pi)),1,[],[],[],[],[]); % 0 for nearest neighboring filtering, 1 for bilinear filtering
			if num_rings >=8
				Screen('DrawTextures',w,eleTex1P,[],dstRect8,radtodeg(mod(ang2+(pp-1)*shiftAng,2*pi)),1,[],[],[],[],[]); % 0 for nearest neighboring filtering, 1 for bilinear filtering
				if num_rings >=9
					Screen('DrawTextures',w,eleTex1P,[],dstRect9,radtodeg(mod(ang1+(qq-1)*shiftAng,2*pi)),1,[],[],[],[],[]); % 0 for nearest neighboring filtering, 1 for bilinear filtering
					if num_rings >=10
						Screen('DrawTextures',w,eleTex1P,[],dstRect10,radtodeg(mod(ang2+(rr-1)*shiftAng,2*pi)),1,[],[],[],[],[]); % 0 for nearest neighboring filtering, 1 for bilinear filtering
						if num_rings >=11
							Screen('DrawTextures',w,eleTex1P,[],dstRect11,radtodeg(mod(ang1+(ss-1)*shiftAng,2*pi)),1,[],[],[],[],[]); % 0 for nearest neighboring filtering, 1 for bilinear filtering
						end %11
					end  %10
				end %9
			end %8
			
			
			%%%%%%%%************************** add mask
			if is_mask_outer == 1
				Screen('BlendFunction',w,GL_ONE_MINUS_SRC_ALPHA,GL_SRC_ALPHA,[1 1 1 0]);
				Screen('DrawTexture',w,masktex);
				Screen('BlendFunction',w,GL_ONE,GL_ZERO,[1 1 1 1]);
			end
			if is_mask_inner == 1
				Screen('BlendFunction',w,GL_ONE_MINUS_SRC_ALPHA,GL_SRC_ALPHA,[1 1 1 0]);
				Screen('DrawTexture',w,masktex2);
				Screen('BlendFunction',w,GL_ONE,GL_ZERO,[1 1 1 1]);
			end
			%---------------------------------
			Screen('FillOval',w,fixColor,[xCen-0.5*fixSide,yCen-0.5*fixSide,xCen+0.5*fixSide,yCen+0.5*fixSide]);
			
			
			
			if approach ==1
				i = i+1;
				ii = i;
				if i == onFrames
					i = 1;
					ii = i;
				end
				
				j = j+1;
				jj = j;
				if j == onFrames
					j = 1;
					jj = j;
				end
				
				k = k+1;
				kk = k;
				if k == onFrames
					k = 1;
					kk = k;
				end
				
				l = l+1;
				ll = l;
				if l == onFrames
					l = 1;
					ll = l;
				end
				
				m = m+1;
				mm = m;
				if m == onFrames
					m = 1;
					mm = m;
				end
				n = n+1;
				nn = n;
				if n == onFrames
					n = 1;
					nn = n;
				end
				o = o+1;
				oo = o;
				if o == onFrames
					o = 1;
					oo = o;
				end
				if num_rings >=8
					p = p +1;
					pp = p;
					if p == onFrames
						p = 1;
						pp = p;
					end
					if num_rings >=9
						q = q +1;
						qq = q;
						if q == onFrames
							q = 1;
							qq = q;
						end
						if num_rings >= 10
							r = r+1;
							rr = r;
							if r == onFrames
								r = 1;
								rr = r;
							end
							
							if num_rings >=11
								s = s+1;
								ss = s;
								if s == onFrames
									s = 1;
									ss = s;
								end
								
								if num_rings >=12
									t = t+1;
									tt = t;
									if t == onFrames
										t = 1;
										tt = t;
									end
								end %12
							end %11
						end  %10
					end   %9
				end  %8
				
			end
			
			
			if approach == 2   %  只为了静止时，真实旋转可以实现
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
			end
			if benchmark == 0
				[vbl,~,~,missed] = Screen('Flip',w,vbl+halfifi);
				if missed>0;fprintf('---!!! Missed frame !!!---\n');end
			else
				Screen('Flip',w,0,2,2);
			end
			[~,~,keycode] = KbCheck(-1);
			if keycode(esc)
				breakLoop = true;
				break;
			end
		end
		%       Pinna_grating_match(row);
		%       ana.result(1,iii) = number_area_speed;
		%result = ;实时赋值
		%     Screen('FillRect',w,gray,[]);
		%     Screen('FillOval',w,fixColor,[xCen-fixSide/2,yCen-fixSide/2,xCen+fixSide/2,yCen+fixSide/2]);
		%     Screen('Flip',w);
		%     WaitSecs(1); % wait for 2 seconds
		if benchmark > 0
			avgfps = count / (GetSecs - ts);
			fprintf('---> The average FPS was: %f fps.\n',avgfps);
		end
	end
	%save result .mat
	%  save(file.sta_fileName,'ana');
	%  save('ana');
	Screen('FillRect',w,gray,[]);
	Screen('FillOval',w,fixColor,[xCen-fixSide/2,yCen-fixSide/2,xCen+fixSide/2,yCen+fixSide/2]);
	Screen('Flip',w);
	WaitSecs(1); % wait for 2 seconds
	Screen('CloseAll');
	ShowCursor;
catch
	Screen('CloseAll');
	Priority(0);
	psychrethrow(psychlasterror);
	
end

return;