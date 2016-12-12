function Pinna_grating_withlineP(is_with_line, benchmark)
if ~exist('is_with_line','var'); is_with_line=false; end
if ~exist('benchmark','var'); benchmark=false; end

PsychDefaultSetup(2);
KbName('UnifyKeyNames');
esc = KbName('escape');
Screen('Preference', 'SkipSyncTests', 1)

% viewing parameters ------------------------------------------------------
screenNumber = max(Screen('Screens'));%-1;
pixelsPerCm = 35;
sizePixel = 10/pixelsPerCm;   %in mm, calculated from different Screen  ,has different value
distance = 56.5;
stdDis = distance*10; %mm
windowed = [];
backgroundColour = [0.5 0.5 0.5 1];
% one_trials = 5;
% approach = ones(trials); %[0 1]; % simulate approaching (1) or leaving (0)
approach = 1;
% directions = [0 1]; % whether the inner ring has CW (0) or CCW (1) rotational direction when approaching, equivalently CCW (0) or CW (1) when leaving
speeds = 200;  % translation speed, in mm/sec
angSpeed = 10; % rotation speeds, in degrees/sec, +: CW, -: CCW
allAngle1 =  45; %angle of inclination  ,absolute value  no + -
num_rings = 10;
is_mask_outer = 0;
is_mask_inner = 0;
f = 0.05;
maskinner_radius = 2;
maskouter_radius = 10;
eachConditionSecs = 3; %1s
% number of elements
num1 = 20;

%for 1st ring
offset1 = 0;
ang1 = mod(offset1+linspace(0,360-360/num1,num1),360);  %
ang1 = deg2rad(ang1); %%
offset2 = 180/num1;  % when different rings ,second ring's position is different
ang2 = mod(offset2+linspace(0,360-360/num1,num1),360);  %
ang2 = deg2rad(ang2); %%

% radius of ring of elements
r1 = 0.4; % in degree  

%--------------------------------------------------------------------------
try
	%-----------------------open the PTB screens------------------------
	sM = screenManager('verbose',false,'screen',screenNumber,...
		'pixelsPerCm',pixelsPerCm,...
		'distance',distance,'bitDepth','8bit',...
		'debug',false,'antiAlias',0, ...
		'blend',true,'srcMode','GL_ONE','dstMode','GL_ZERO',...
		'windowed',windowed,'backgroundColour',[backgroundColour],...
		'gammaTable', []); %use a temporary screenManager object
	screenVals = open(sM); %open PTB screen
	w = sM.win;
	wrect = sM.winRect;
	xCen = sM.xCenter;
	yCen = sM.yCenter;
	ifi = sM.screenVals.ifi;
	halfifi = ifi/2;
	ppd = sM.ppd;
	[os,od,oc]=Screen('BlendFunction',w);
	
	% screen
	white = sM.screenVals.white;
	black = sM.screenVals.black;
	gray = sM.screenVals.gray;
	
	%==============procedural gabor===========================
	degsize = 3;
	sizeGabor = [degsize*ppd degsize*ppd];
	phase = 90;
	sc = 8.0;
	freq = 2/ppd; %cycles per pixel
	contrast = 0.99;
	aspectratio = 1.0;
	nonsymmetric = false;
	backgroundColorOffset = [0.5 0.5 0.5 1];
	disableNorm = true;
	contrastPreMultiplicator = 0.5;
	mypars = [phase, freq, sc, contrast, aspectratio, 0, 0, 0]';
	allPars = repmat(mypars,1,num_rings*num1);
	gabortex = CreateProceduralGabor(w, sizeGabor(1), sizeGabor(2), nonsymmetric, backgroundColorOffset, disableNorm, contrastPreMultiplicator);
	
	srcRect = Screen('Rect',gabortex);
	dstRect1 = ScaleRect(srcRect, 0.25, 0.25);
	dstRect2 = srcRect;
	dstRect3 = ScaleRect(srcRect, 4, 4);
	dstRect1 = CenterRectOnPointd(dstRect1, xCen-300, yCen);
	dstRect2 = CenterRectOnPointd(dstRect2, xCen, yCen);
	dstRect3 = CenterRectOnPointd(dstRect3, xCen+300, yCen);
	dstRects = [dstRect1', dstRect2', dstRect3'];
	
	auxP = [phase, freq, sc, contrast, aspectratio, 0, 0, 0]';
	auxP = repmat(auxP,1,3);
	
	fprintf('===>>> Draw example procedural gabors\n');
	Screen('DrawTextures', w, gabortex, [], dstRects, [45, 45, 45], ...
		[], [], [], [], kPsychDontDoRotation, auxP);
	Screen('DrawText', w, 'These are the source procedural gabors, scaled x0.25 x1 and x2',0,0);
	sM.flip;
	WaitSecs('Yieldsecs',2);
	
% 	% We create a Luminance+Alpha matrix for use as transparency mask:
% 	% Layer 1 (Luminance) is filled with luminance value 'gray' of the
% 	% background.
% 	ms = 100;
% 	transLayer = 2;
% 	[x,y] = meshgrid(-ms:ms, -ms:ms);
% 	maskblob = ones(2*ms+1, 2*ms+1, transLayer) * gray;
% 	size(maskblob);
% 	
% 	% Layer 2 (Transparency aka Alpha) is filled with gaussian transparency
% 	% mask.
% 	xsd = ms/2.0;
% 	ysd = ms/2.0;
% 	maskblob(:,:,transLayer)= white - exp(-((x/xsd).^2)-((y/ysd).^2))*white;
% 	
% 	% Build a single transparency mask texture
% 	masktex=Screen('MakeTexture', w, maskblob);
	
	maskinner_radius = maskinner_radius * ppd;% convert degrees to pixels
	maskouter_radius = maskouter_radius * ppd;% convert degrees to pixels
	
	aLayer = 4;
	if is_mask_outer == 1
		masksize = [wrect(4)/2,wrect(3)/2];
		mask = ones(masksize(1),masksize(2),aLayer);
		mask(:,:,1) = 0.4;
		[x,y] = meshgrid(0:masksize(2)-1,0:masksize(1)-1);
		mask(:, :, aLayer) = white *(exp(-((x-xCen).^2+(y-yCen).^2)/(maskouter_radius)^2)); % gaussian mask
		masktex = Screen('MakeTexture', w, mask);
	end
	if is_mask_inner == 1
		masksize = [wrect(4),wrect(3)];
		mask2 = ones(masksize(1),masksize(2),aLayer)*gray;
		[x,y] = meshgrid(0:masksize(2)-1,0:masksize(1)-1);
		mask2(:, :, aLayer)=white * (1-exp(-((x-xCen).^2+(y-yCen).^2)/((maskinner_radius)*2)^2)); % gaussian mask
		masktex2 = Screen('MakeTexture', w, mask2);
	end
	
	r1Origin = r1 * ppd;% convert degrees to pixels
	
	% radius of ring of elements
	onFrames = fix(stdDis/(speeds*ifi) +1);  %li:为了必须超出screen范围，取最大的帧数，同时不能超过半径公式的极大值， controllde through d0s and speeds
	shiftAng = deg2rad(angSpeed*ifi); % degree of rotation per frame  %%真实旋转
	
	if speeds == 0
		onFrames = 200;  %
		approach = 2;  %when speed = 0,  approach = 2；
	end
	
	breakLoop = false;
	Priority(MaxPriority(w));
	Screen('FillRect',w,gray,[]);
	sM.drawCross(0.4,[0 0 0 1]);
	Screen('Flip',w);
	fprintf('\n===>>> Start Pinna display...\n');
	WaitSecs(1);
	
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
		ii = 1;		jj = 1;         kk = 1;		ll = 1;        mm = 1;
		nn = 1;        oo = 1;		pp = 1;        qq = 1;
		rr = 1;        ss = 1;		tt = 1;        j = 1;
		k = 1;        l = 1;		m = 1;        n = 1;
		o = 1;        p = 1;		q = 1;        r = 1;
		s = 1;        t = 1;
		
		switch num_rings
			case 8
				j = round(onFrames/8);
				k = round(2*onFrames/8);
				l = round(3*onFrames/8);
				m = round(4*onFrames/8);
				n = round(5*onFrames/8);
				o = round(6*onFrames/8);
				p = round(7*onFrames/8);
				jj = j;            kk = k;
				ll = l;            mm = m;
				nn = n;            oo = o;
				pp = p;
			case 9
				j = round(onFrames/9);
				k = round(2*onFrames/9);
				l = round(3*onFrames/9);
				m = round(4*onFrames/9);
				n = round(5*onFrames/9);
				o = round(6*onFrames/9);
				p = round(7*onFrames/9);
				q = round(8*onFrames/9);
				jj = j;            kk = k;
				ll = l;            mm = m;
				nn = n;            oo = o;
				pp = p;            qq = q;
			case 10
				j = round(onFrames/10);
				k = round(2*onFrames/10);
				l = round(3*onFrames/10);
				m = round(4*onFrames/10);
				n = round(5*onFrames/10);
				o = round(6*onFrames/10);
				p = round(7*onFrames/10);
				q = round(8*onFrames/10);
				r = round(9*onFrames/10);
				jj = j;            kk = k;
				ll = l;            mm = m;
				nn = n;            oo = o;
				pp = p;            qq = q;
				rr = r;
			case 11
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
				jj = j;            kk = k;
				ll = l;            mm = m;
				nn = n;            oo = o;
				pp = p;            qq = q;
				rr = r;            ss = s;
		end
		
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		if approach ==2   %%每个trial ,ii ,jj ……初始化
			ii = i;            jj = j;
			kk = k;            ll = l;
			mm = m;            nn = n;
			oo = o;            pp = p;
			qq = q;            rr = r;
			ss = s;            tt = t;
		end
		
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
		

		vbl = Screen('Flip',w);
		ts = vbl;
		count = 0;
		vblendtime = vbl + eachConditionSecs;
		%-----------
		while GetSecs < vblendtime
			
			count = count + 1;
			
			a1 = rad2deg(mod(ang1+(ii-1)*shiftAng,2*pi));
			r1 = r1Origin + 0.5*r_a*i^2;
			side1P(1) = min_size + 0.5*size_a*i^2;
			side1P(2) = side1P(1);
			dstRect1 = [xCen+r1*sin(mod(ang1+(ii-1)*shiftAng,2*pi))-side1P(2)/2;yCen-r1*cos(mod(ang1+(ii-1)*shiftAng,2*pi))-side1P(1)/2;...
				xCen+r1*sin(mod(ang1+(ii-1)*shiftAng,2*pi))+side1P(2)/2;yCen-r1*cos(mod(ang1+(ii-1)*shiftAng,2*pi))+side1P(1)/2];
			
			a2 = rad2deg(mod(ang2+(jj-1)*shiftAng,2*pi));
			r2 = r1Origin + 0.5*r_a*j^2;
			side2P(1) = min_size + 0.5*size_a*j^2;
			side2P(2) = side2P(1);
			dstRect2 = [xCen+r2*sin(mod(ang2+(jj-1)*shiftAng,2*pi))-side2P(2)/2;yCen-r2*cos(mod(ang2+(jj-1)*shiftAng,2*pi))-side2P(1)/2;...
				xCen+r2*sin(mod(ang2+(jj-1)*shiftAng,2*pi))+side2P(2)/2;yCen-r2*cos(mod(ang2+(jj-1)*shiftAng,2*pi))+side2P(1)/2];
			
			a3 = rad2deg(mod(ang1+(kk-1)*shiftAng,2*pi));
			r3 = r1Origin + 0.5*r_a*k^2;
			side3P(1) = min_size + 0.5*size_a*k^2;
			side3P(2) = side3P(1);
			dstRect3 = [xCen+r3*sin(mod(ang1+(kk-1)*shiftAng,2*pi))-side3P(2)/2;yCen-r3*cos(mod(ang1+(kk-1)*shiftAng,2*pi))-side3P(1)/2;...
				xCen+r3*sin(mod(ang1+(kk-1)*shiftAng,2*pi))+side3P(2)/2;yCen-r3*cos(mod(ang1+(kk-1)*shiftAng,2*pi))+side3P(1)/2];
			
			a4 = rad2deg(mod(ang2+(ll-1)*shiftAng,2*pi));
			r4 = r1Origin + 0.5*r_a*l^2;
			side4P(1) = min_size + 0.5*size_a*l^2;
			side4P(2) = side4P(1);
			dstRect4 = [xCen+r4*sin(mod(ang2+(ll-1)*shiftAng,2*pi))-side4P(2)/2;yCen-r4*cos(mod(ang2+(ll-1)*shiftAng,2*pi))-side4P(1)/2;...
				xCen+r4*sin(mod(ang2+(ll-1)*shiftAng,2*pi))+side4P(2)/2;yCen-r4*cos(mod(ang2+(ll-1)*shiftAng,2*pi))+side4P(1)/2];
			
			a5 = rad2deg(mod(ang1+(mm-1)*shiftAng,2*pi));
			r5 = round(r1Origin + 0.5*r_a*m^2);
			side5P(1) = round(min_size + 0.5*size_a*m^2);
			side5P(2) = side5P(1);
			dstRect5 = [xCen+r5*sin(mod(ang1+(mm-1)*shiftAng,2*pi))-side5P(2)/2;yCen-r5*cos(mod(ang1+(mm-1)*shiftAng,2*pi))-side5P(1)/2;...
				xCen+r5*sin(mod(ang1+(mm-1)*shiftAng,2*pi))+side5P(2)/2;yCen-r5*cos(mod(ang1+(mm-1)*shiftAng,2*pi))+side5P(1)/2];
			
			a6 = rad2deg(mod(ang2+(nn-1)*shiftAng,2*pi));
			r6 = r1Origin + 0.5*r_a*n^2;
			side6P(1) = min_size + 0.5*size_a*n^2;
			side6P(2) = side6P(1);
			dstRect6 = [xCen+r6*sin(mod(ang2+(nn-1)*shiftAng,2*pi))-side6P(2)/2;yCen-r6*cos(mod(ang2+(nn-1)*shiftAng,2*pi))-side6P(1)/2;...
				xCen+r6*sin(mod(ang2+(nn-1)*shiftAng,2*pi))+side6P(2)/2;yCen-r6*cos(mod(ang2+(nn-1)*shiftAng,2*pi))+side6P(1)/2];
			
			a7 = rad2deg(mod(ang1+(oo-1)*shiftAng,2*pi));
			r7 = r1Origin + 0.5*r_a*o^2;
			side7P(1) = min_size + 0.5*size_a*o^2;
			side7P(2) = side7P(1);
			dstRect7 = [xCen+r7*sin(mod(ang1+(oo-1)*shiftAng,2*pi))-side7P(2)/2;yCen-r7*cos(mod(ang1+(oo-1)*shiftAng,2*pi))-side7P(1)/2;...
				xCen+r7*sin(mod(ang1+(oo-1)*shiftAng,2*pi))+side7P(2)/2;yCen-r7*cos(mod(ang1+(oo-1)*shiftAng,2*pi))+side7P(1)/2];
						
			dstRectA = [dstRect1, dstRect2, dstRect3, dstRect4, dstRect5, dstRect6, dstRect7];
			aA = [a1 a2 a3 a4 a5 a6 a7];
			
			if num_rings >=8
				a8 = rad2deg(mod(ang2+(pp-1)*shiftAng,2*pi));
				r8 = r1Origin + 0.5*r_a*p^2;
				side8P(1) = min_size + 0.5*size_a*p^2;
				side8P(2) = side8P(1);
				dstRect8 = [xCen+r8*sin(mod(ang2+(pp-1)*shiftAng,2*pi))-side8P(2)/2;yCen-r8*cos(mod(ang2+(pp-1)*shiftAng,2*pi))-side8P(1)/2;...
					xCen+r8*sin(mod(ang2+(pp-1)*shiftAng,2*pi))+side8P(2)/2;yCen-r8*cos(mod(ang2+(pp-1)*shiftAng,2*pi))+side8P(1)/2];
				dstRectA = [dstRectA, dstRect8];
				aA = [aA, a8];
				if num_rings >=9
					a9 = rad2deg(mod(ang1+(qq-1)*shiftAng,2*pi));
					r9 = r1Origin + 0.5*r_a*q^2;
					side9P(1) = min_size + 0.5*size_a*q^2;
					side9P(2) = side9P(1);
					dstRect9 = [xCen+r9*sin(mod(ang1+(qq-1)*shiftAng,2*pi))-side9P(2)/2;yCen-r9*cos(mod(ang1+(qq-1)*shiftAng,2*pi))-side9P(1)/2;...
						xCen+r9*sin(mod(ang1+(qq-1)*shiftAng,2*pi))+side9P(2)/2;yCen-r9*cos(mod(ang1+(qq-1)*shiftAng,2*pi))+side9P(1)/2];
					dstRectA = [dstRectA, dstRect9];
					aA = [aA, a9];
					if num_rings>= 10
						a10 = rad2deg(mod(ang2+(rr-1)*shiftAng,2*pi));
						r10 = r1Origin + 0.5*r_a*r^2;
						side10P(1) = min_size + 0.5*size_a*r^2;
						side10P(2) = side10P(1);
						dstRect10 = [xCen+r10*sin(mod(ang2+(rr-1)*shiftAng,2*pi))-side10P(2)/2;yCen-r10*cos(mod(ang2+(rr-1)*shiftAng,2*pi))-side10P(1)/2;...
							xCen+r10*sin(mod(ang2+(rr-1)*shiftAng,2*pi))+side10P(2)/2;yCen-r10*cos(mod(ang2+(rr-1)*shiftAng,2*pi))+side10P(1)/2];
						dstRectA = [dstRectA, dstRect10];
						aA = [aA, a10];
						if num_rings >=11
							a11 = rad2deg(mod(ang1+(ss-1)*shiftAng,2*pi));
							r11 = r1Origin + 0.5*r_a*s^2;
							side11P(1) = round(min_size + 0.5*size_a*s^2);
							side11P(2) = side11P(1);
							dstRect11 = [xCen+r11*sin(mod(ang1+(ss-1)*shiftAng,2*pi))-side11P(2)/2;yCen-r11*cos(mod(ang1+(ss-1)*shiftAng,2*pi))-side11P(1)/2;...
								xCen+r11*sin(mod(ang1+(ss-1)*shiftAng,2*pi))+side11P(2)/2;yCen-r11*cos(mod(ang1+(ss-1)*shiftAng,2*pi))+side11P(1)/2];
							dstRectA = [dstRectA, dstRect11];
							aA = [aA, a11];
						end  %11
					end %10
				end  %9
			end   %8
			
			%%%%%%withline
			if is_with_line
				Screen('Drawlines',w,xy1,1,200,[xCen yCen],[]);
				Screen('Drawlines',w,xy2,1,200,[xCen yCen],[]);
			end
			
			%Screen('BlendFunction',w,GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
			Screen('DrawTextures', w, gabortex, [], dstRectA, aA, [], [], [], [], kPsychDontDoRotation, allPars);
			if is_mask_outer == 1
				Screen('BlendFunction',w,GL_ONE,GL_ONE);
				Screen('DrawTexture',w,masktex);
				Screen('BlendFunction',w,GL_ONE,GL_ZERO,[1 1 1 1]);
			end
			if is_mask_inner == 1
				Screen('BlendFunction',w,GL_ONE_MINUS_SRC_ALPHA,GL_SRC_ALPHA);
				Screen('DrawTexture',w,masktex2);
				Screen('BlendFunction',w,GL_ONE,GL_ZERO);
			end
			sM.drawCross(0.4,[0 0 0 1]);
			Screen('DrawingFinished', w);
			
			if approach ==1
				i = i+1; 				ii = i;
				if i == onFrames;	i = 1; ii = i; end
				
				j = j+1;				jj = j;
				if j == onFrames;	j = 1; jj = j; end
				
				k = k+1; kk = k;
				if k == onFrames; k = 1; kk = k; end
				
				l = l+1; ll = l;
				if l == onFrames;	l = 1; ll = l; end
				
				m = m+1;				mm = m;
				if m == onFrames;	m = 1; mm = m; end
				
				n = n+1;	nn = n;
				if n == onFrames;	n = 1; nn = n; end
				
				o = o+1;				oo = o;
				if o == onFrames; o = 1; oo = o; end
				
				if num_rings >=8
					p = p +1; pp = p;
					if p == onFrames; p = 1; pp = p; end
					
					if num_rings >=9
						q = q +1; qq = q;
						if q == onFrames;	q = 1; qq = q; end
						
						if num_rings >= 10
							r = r+1; rr = r;
							if r == onFrames;	r = 1; rr = r; end
							
							if num_rings >=11
								s = s+1;		ss = s;
								if s == onFrames;	s = 1; ss = s; end
								
								if num_rings >=12
									t = t+1;	tt = t;
									if t == onFrames; t = 1; tt = t;end
								end %12
							end %11
						end  %10
					end   %9
				end  %8
				
			elseif approach == 2   %  只为了静止时，真实旋转可以实现
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
		if benchmark > 0
			avgfps = count / (GetSecs - ts);
			fprintf('---> The average FPS was: %f fps.\n',avgfps);
		end
	end
	%save result .mat
	%  save(file.sta_fileName,'ana');
	%  save('ana');
	sM.drawBackground();
	sM.drawCross(0.4,[0 0 0 1]);
	Priority(0);
	Screen('Flip',w);
	WaitSecs('YieldSecs',1); 
	Screen('Close', gabortex); Screen('Close', masktex);
	sM.close();
	ShowCursor;
catch ME
	Screen('CloseAll');
	Priority(0);
	getReport(ME);
end

return;