function Pinna_grating_main(angle_pattern,radial_speed_i,...
	rotation_speed_i,ResultDir,one_trials,duration,...
	match_time,is_binary_mask,mask_diameter,mask_xpos,...
	mask_ypos, match_value, use_eyelink, use_staircase, ...
	fix_radius,staircase_use_radial,calib_file)

global eL sM ana w wrect ifi waitFrames ang1 ang2 ...
	white black gray approach spa spb spc esc ok ckey r1Origin ...
	abandon move_speed xc xxc yc ovalRect r1Match  eleTexMatch1 ...
	txtColorMat shiftAng1 onFrames numChoice rot_speed breakLoop ...
	prior45 prior90 prior135;

if nargin == 0; error('No parameters passed');end

ins = inputdlg({'Subject Name','Comments (room, lights etc.)'});
if isempty(ins); return; end
ana=[]; %this holds the experiment data and parameters
ana.subject = ins{1};
ana.comments = ins{2};
ana.date = datestr(datetime);
ana.version = Screen('Version');
ana.computer = Screen('Computer');
ana.calibFile = calib_file;
useEyeLink = use_eyelink;
isDummy = false;
PsychDefaultSetup(0);
Screen('Preference', 'SkipSyncTests', 0);
KbName('UnifyKeyNames');
esc = KbName('q');
ok = KbName('uparrow');
spa = KbName('leftarrow');
spb = KbName('rightarrow');
spc = KbName('downarrow');
ckey = KbName('c');

%--------------fix parameters
fixX = 0;
fixY = 0;
firstFixInit = 1;
firstFixTime = 1;
firstFixRadius = fix_radius;
strictFixation = true;

%----------------Make a name for this run-----------------------
pf='PINNA_';
nameExp = [pf ana.subject];
c = sprintf(' %i',fix(clock()));
c = regexprep(c,' ','_');
nameExp = [nameExp c];

%---------------------- viewing parameters -------------------------------
screenID = max(Screen('Screens'));%-1;
pixelsPerCm = 35;
sizePixel = 10 / pixelsPerCm;   %in mm, calculated from different Screen  ,has different value
distance = 56.5;
windowed = [];
backgroundColor = [127.5 127.5 127.5];
stdDis = distance*10; %mm
approach = 1;%[0 1]; % simulate approaching (1) or leaving (0)
% directions = [0 1]; % whether the inner ring has CW (0) or CCW (1) rotational direction when approaching, equivalently CCW (0) or CW (1) when leaving
% speeds = 200;  % translation speed, in mm/sec
% angSpeed = 0; % rotation speeds, in degrees/sec, +: CW, -: CCW
gaborAngles =  angle_pattern;%[45 135]; %angle of inclination  ,absolute value  no + -
num_rings = 10;
is_mask_outer = 1;
is_mask_inner = 1;
f = 0.05; %grating frequency
grating_size = 15; %micropattern size (pixel)
maskinner_radius = 2.5; %deg
maskouter_radius = 9; %deg
eachConditionSecs = duration; %sec
% number of elements
numElements = 20;

%for 1st ring
offset1 = 0;
ang1 = mod(offset1+linspace(0,360-360/numElements,numElements),360);  %
ang1 = deg2rad(ang1);
offset2 = 180/numElements;  % when different rings ,second ring's position is different
ang2 = mod(offset2+linspace(0,360-360/numElements,numElements),360);  %
ang2 = deg2rad(ang2);

% radius of ring of elements
r1 = 0.4; % in degree

%============================SET UP VARIABLES=====================================
numgaborAngles = size(gaborAngles,2);
num_move_speed = length(radial_speed_i);
num_angle_speed = length(rotation_speed_i);
% if num_move_speed==num_angle_speed
% 	trials = numgaborAngles * one_trials * num_move_speed;  %be sure one by one
% else
% 	msgbox('The number of two kinds of speed is different', 'warning');
% end

if use_staircase
	setupTrial();
	pattern_angle_index = ones(1,stopRule*3);
	column_rank = randperm(stopRule*3);
	for i = 1:stopRule*3
		pattern_angle_index(i) = fix((column_rank(i)-1)/stopRule)+1;
	end
else
	trials = numgaborAngles * one_trials * num_move_speed;
	column_rank = randperm(trials);  %rand(50)
	pattern_angle_index = ones(1,trials); %value is 1~num_angle
	% condition_index = ones(1,trials); %value is 1~num_conditions
	rotation_speed_index = ones(1,trials); %value is 1~num_speed
	radial_speed_index = ones(1,trials);
	for i = 1:trials
		pattern_angle_index(i) = fix((column_rank(i)-1)/(one_trials * num_move_speed))+1; %value is 1~num_angle
		rotation_speed_index(i) = fix(rem(column_rank(i)-1,(one_trials * num_move_speed))/(one_trials)) + 1;   %value is 1~num_rotation_speed
		radial_speed_index = rotation_speed_index;%value is 1~num_move_speed And num_rotation_speed = num_move_speed; because one-to-one
	end
end
%============================================================================

try
	white = WhiteIndex(screenID);
	black = BlackIndex(screenID);
	gray = (white+black)/2; % index for white, black and gray
	if gray == white
		gray = white/2;
	end
	inc = white-gray;
	PsychImaging('PrepareConfiguration');
	PsychImaging('AddTask', 'General', 'UseFastOffscreenWindows');
	%PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible');
	[w,wrect] = PsychImaging('OpenWindow',screenID,gray,windowed);
	ana.ifi = Screen('GetFlipInterval',w);
	ana.halfifi = ana.ifi/2;
	xCen = wrect(3)/2;
	yCen = wrect(4)/2;
	ScreenWidth = round(wrect(3)*sizePixel);
	ana.ppd = wrect(3)/2/atand(ScreenWidth/2/stdDis);
	
	sM = screenManager;
	sM.screen = screenID;
	sM.pixelsPerCm = pixelsPerCm;
	sM.distance = distance;
	sM.backgroundColour = backgroundColor;
	sM.forceWin(w);
	
	%==============================setup eyelink==========================
	if useEyeLink == true
		eL = eyelinkManager('IP',[]);
		%eL.verbose = true;
		eL.isDummy = isDummy; %use dummy or real eyelink?
		eL.name = nameExp;
		eL.saveFile = [nameExp '.edf'];
		eL.recordData = true; %save EDF file
		eL.sampleRate = 500;
		eL.remoteCalibration = false; % manual calibration?
		eL.calibrationStyle = 'HV5'; % calibration style
		eL.modify.calibrationtargetcolour = [0 0 0];
		eL.modify.calibrationtargetsize = 0.5;
		eL.modify.calibrationtargetwidth = 0.05;
		eL.modify.waitformodereadytime = 500;
		eL.modify.devicenumber = -1; % -1 = use any keyboard
		% X, Y, FixInitTime, FixTime, Radius, StrictFix
		updateFixationValues(eL, fixX, fixY, firstFixInit, firstFixTime, firstFixRadius, strictFixation);
		initialise(eL, sM); %use sM to pass screen values to eyelink
		setup(eL); % do setup and calibration
		WaitSecs('YieldSecs',0.25);
		getSample(eL); %make sure everything is in memory etc.
	else
		eL = [];
	end
	
	%===================================================
	%global struct ana to save parameters
	saveMetaData();
	%======================================================
	
	% screen areas
	xc = (1:2).*wrect(3)/3; %match -only have A and B
	xxc = (1:3).*wrect(3)/4; %match -only have A and B
	yc = wrect(4)/2;
	r1Match = round(45.0/sizePixel); %30mm/0.25=120pixels ;max radius at match
	
	% 3 fixPoint
	%fixRect = [xCen-fixSide/2,yCen-fixSide/2,xCen+fixSide/2,yCen+fixSide/2];
	fixSideMatch = round(1.7/sizePixel);%6;
	ovalRect = zeros(4,2);
	ovalRect(1,1:2) = xc-0.5*fixSideMatch;
	ovalRect(2,1:2) = yc-0.5*fixSideMatch;
	ovalRect(3,1:2) = xc+0.5*fixSideMatch;
	ovalRect(4,1:2) = yc+0.5*fixSideMatch;
	
	r1Origin = r1 * ana.ppd;% convert degrees to pixels
	
	if is_binary_mask
		mask_diameter = mask_diameter * ana.ppd;
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
	
	%=====================MAKE GABORS===============================
	AllElements = [];
	%make grating texture
	f=f*2*pi; %used for angle1 and angle2
	for jj = 1:numgaborAngles   %number of the MakeTexture is 20 ,at least
		[x,y]=meshgrid(-30:30,-30:30);  %used for angle1 and angle2
		angle1=(-gaborAngles(jj)-90)*pi/180;
		a1=cos(angle1)*f;
		b1=sin(angle1)*f;
		m1=exp(-((x/grating_size).^2)-((y/grating_size).^2)).*sin(a1*x+b1*y-5);
		AllElements(1,jj).eleTex1P = Screen('MakeTexture',w,gray+inc*m1);    %CCW
	end
	
	%-----------------outside mask & inside mask-------------------------
	maskinner_radius = maskinner_radius * ana.ppd;% convert degrees to pixels
	maskouter_radius = maskouter_radius * ana.ppd;% convert degrees to pixels
	
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
	%------------------------------------------------------------
	
	figH = figure('Position',[100 100 700 600],'NumberTitle','off','Name',...
		['Subject: ' ana.subject  ' started ' ana.date ' | ' ana.comments]);
	box on; grid on; grid minor; ylim([0 3]);
	xlabel('Trials ')
	ylabel('Response')
	title('Pinna Psychophysics')
	drawnow; WaitSecs('YieldSecs',0.25);
	
	if ~useEyeLink
		Screen('FillRect',w,gray,[]);
		normBoundsRect = Screen('TextBounds', w, 'Please fix at the central point.');
		Screen('DrawText',w,'Please fix at the central point.',xCen-normBoundsRect(3)/2,yCen+normBoundsRect(4));
		drawCross(sM,0.3,[0 0 0 1],fixX,fixY);
		Screen('Flip',w);
		WaitSecs('YieldSecs',2);
	end
	
	% initialise our trial variables
	iii = 1;
	breakLoop = false;
	fixated = 'no';
	response = NaN;
	numChoice = -1;
	
	%%%%%%%%%%%%%%%%%%%%%%--WE LOOP BABY--%%%%%%%%%%%%%%%%%%%%%%%%%
	while ~breakLoop
		fixated = 'no';
		response = NaN;
		numChoice = -1;
		sM.drawBackground();
		Screen('Flip',w);
		WaitSecs('YieldSecs',1);
		Priority(MaxPriority(w));
		ifi = ana.ifi;
		
		%======================SET UP THIS TRIAL PARAMATERS========================
		if iii<=trials %initial fixation held
			if use_staircase
				if staircase_use_radial
					waitFrames = 1;
					rot_speed = rotation_speed_i(1);
					shiftAng1 = deg2rad( rot_speed .* ifi);
					if pattern_angle_index(iii) == 1
						move_speed = staircase45.xCurrent * ana.ppd;
					elseif pattern_angle_index(iii) == 2
						move_speed = staircase90.xCurrent * ana.ppd;
					elseif pattern_angle_index(iii) == 3
						move_speed = staircase135.xCurrent * ana.ppd;
					end
				else
					waitFrames = 1;
					if pattern_angle_index(iii) == 1
						rot_speed = staircase45.xCurrent;
						shiftAng1 = deg2rad( rot_speed .* ifi);
					elseif pattern_angle_index(iii) == 2
						rot_speed = staircase90.xCurrent;
						shiftAng1 = deg2rad( rot_speed .* ifi);
					elseif pattern_angle_index(iii) == 3
						rot_speed = staircase135.xCurrent;
						shiftAng1 = deg2rad( rot_speed .* ifi);
					end
					move_speed = radial_speed_i(1) * ana.ppd;
				end
			else
				waitFrames = 1;
				rot_speed = rotation_speed_i(rotation_speed_index(iii));
				shiftAng1 = deg2rad( rot_speed .* ifi);
				move_speed = radial_speed_i(radial_speed_index(iii)) * ana.ppd;
			end
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
			if rot_speed>0&&move_speed==0 %CW
				speeds = 5*ana.ppd;
				onFrames = fix((r_dis + sqrt(r_dis*r_dis + 4*r_dis*speeds*ifi))/(2*speeds*ifi));
				approach = 2;  %when speed = 0,  approach = 2?�?�?
				shiftAng = shiftAng1; %degtorad(angSpeed*ifi); % degree of rotation per frame  %%?�?��?�?��?
				condition = 3;
			end
			if rot_speed<0&&move_speed==0  %CCW
				speeds = 5*ana.ppd;
				onFrames = fix((r_dis + sqrt(r_dis*r_dis + 4*r_dis*speeds*ifi))/(2*speeds*ifi));
				approach = 2;  %when speed = 0,  approach = 2?�?�?
				shiftAng = shiftAng1; %degtorad(angSpeed*ifi); % degree of rotation per frame  %%?�?��?�?��?
				condition = 4;
			end
			fprintf('==>> Trial = %i rot_speed = %g | move_speed = %g | speed = %g | approach = %g | Condition = %g\n', iii,rot_speed,move_speed,speeds,approach,condition);
			
			row = 1;
			column = pattern_angle_index(iii);
			i = 1;  %when only has one ring
			ii = 1;			jj = 1;			kk = 1;			ll = 1;
			mm = 1;			nn = 1;			oo = 1;			pp = 1;
			qq = 1;			rr = 1;			ss = 1;			tt = 1;
			j = 1;			k = 1;			l = 1;			m = 1;
			n = 1;			o = 1;			p = 1;			q = 1;
			r = 1;			s = 1;			t = 1;
			
			switch num_rings
				case 8
					j = round(onFrames/8);
					k = round(2*onFrames/8);
					l = round(3*onFrames/8);
					m = round(4*onFrames/8);
					n = round(5*onFrames/8);
					o = round(6*onFrames/8);
					p = round(7*onFrames/8);
				case 9
					j = round(onFrames/9);
					k = round(2*onFrames/9);
					l = round(3*onFrames/9);
					m = round(4*onFrames/9);
					n = round(5*onFrames/9);
					o = round(6*onFrames/9);
					p = round(7*onFrames/9);
					q = round(8*onFrames/9);
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
					
			end
			
			if approach ==2   %%ÿ?�?�trial ,ii ,jj ?�?�?�?�?�?��?�?�?
				ii = 1;			jj = 1;			kk = 1;			ll = 1;
				mm = 1;			nn = 1;			oo = 1;			pp = 1;
				qq = 1;			rr = 1;			ss = 1;			tt = 1;
			end
			
			r_a = 2*(max_r-r1Origin)/(onFrames)^2;
			size_a = 2*(max_size-min_size)/(onFrames)^2;
			
			%===============================SUBJECT FIXATE===========================
			if useEyeLink
				resetFixation(eL);
				updateFixationValues(eL, fixX, fixY, firstFixInit, firstFixTime, firstFixRadius, strictFixation);
				trackerClearScreen(eL);
				trackerDrawFixation(eL); %draw fixation window on eyelink computer
				edfMessage(eL,'V_RT MESSAGE END_FIX END_RT');  %this 3 lines set the trial info for the eyelink
				edfMessage(eL,['TRIALID ' num2str(iii)]);  %obj.getTaskIndex gives us which trial we're at
				startRecording(eL);
				statusMessage(eL,'INITIATE FIXATION...');
				fixated = '';
				ListenChar(2);
				syncTime(eL);
				%fprintf('===>>> INITIATE FIXATION Trial = %i\n', iii);
				while ~strcmpi(fixated,'fix') && ~strcmpi(fixated,'breakfix')
					drawCross(sM,0.3,[0 0 0 1],fixX,fixY);
					Screen('Flip',w); %flip the buffer
					getSample(eL);
					fixated=testSearchHoldFixation(eL,'fix','breakfix');
					[keyIsDown, ~, keyCode] = KbCheck(-1);
					if keyIsDown == 1
						rchar = KbName(keyCode); if iscell(rchar);rchar=rchar{1};end
						switch lower(rchar)
							case {'c'}
								fixated = 'breakfix';
								stopRecording(eL);
								setOffline(eL);
								trackerSetup(eL);
								WaitSecs('YieldSecs',2);
							case {'d'}
								fixated = 'breakfix';
								stopRecording(eL);
								driftCorrection(eL);
								WaitSecs('YieldSecs',2);
							case {'escape'}
								fixated = 'breakfix';
								response = -100;
								breakLoop = true;
						end
						%ListenChar(0);
					end
				end
				if strcmpi(fixated,'breakfix')
					response = -1;
					fprintf('===>>> BROKE INITIATE FIXATION Trial = %i\n', iii);
					statusMessage(eL,'Subject Broke Initial Fixation!');
					edfMessage(eL,'MSG:BreakFix');
					resetFixation(eL);
					stopRecording(eL);
					edfMessage(eL,['TRIAL_RESULT ' num2str(response)]);
					setOffline(eL);
					continue
				end
			else %no eyetracker, simple show cross
				fixated = 'fix';
				drawCross(sM,0.3,[0 0 0 1],fixX,fixY);
				Screen('Flip',w);
				WaitSecs('YieldSecs',0.75);
			end
			%if we lost fixation then
			if ~strcmpi(fixated,'fix'); continue; end
			
			%=========================Our actual stimulus drawing loop==========================
			if useEyeLink; edfMessage(eL,'END_FIX'); statusMessage(eL,'Show Stimulus...'); end
			vbl = Screen('Flip',w);
			vblendtime = vbl + eachConditionSecs;
			while vbl < vblendtime
				
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
				
				Screen('DrawTextures',w,AllElements(row,column).eleTex1P,[],dstRectA,aA,1,[],[],[],[],[]); % 0 for nearest neighboring filtering, 1 for bilinear filtering
				
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
				if is_binary_mask
					Screen('BlendFunction',w,GL_ONE_MINUS_SRC_ALPHA,GL_SRC_ALPHA,[1 1 1 0]);
					Screen('DrawTexture',w,masktex11);
					Screen('BlendFunction',w,GL_ONE,GL_ZERO,[1 1 1 1]);
				end
				
				drawCross(sM,0.3,[0 0 0 1],fixX,fixY); %Screen('FillOval',w,fixColor,fixRect);
				Screen('DrawingFinished', w);
				
				if approach == 1
					i = i+1; if i == onFrames;	i = 1; end
					j = j+1;	if j == onFrames;	j = 1; end
					k = k+1;	if k == onFrames;	k = 1; end
					l = l+1;	if l == onFrames;	l = 1; end
					m = m+1; if m == onFrames;	m = 1; end
					n = n+1; if n == onFrames;	n = 1; end
					o = o+1; if o == onFrames;	o = 1; end
					if num_rings >=8
						p = p +1;if p == onFrames;	p = 1;end
						if num_rings >=9
							q = q +1;if q == onFrames;q = 1;end
							if num_rings >= 10
								r = r+1;if r == onFrames;r = 1;end
								if num_rings >=11
									s = s+1;if s == onFrames;s = 1;end
									if num_rings >=12
										t = t+1;if t == onFrames;t = 1;end
									end %12
								end %11
							end  %10
						end   %9
					end  %8
				elseif approach == 0
					i = i-1;if i == 0;i = onFrames;end
					if num_rings>=2
						j = j-1;if j == 0;j = onFrames;end
						if num_rings >=3
							k = k-1;if k == 0;k = onFrames;end
							if num_rings >=4
								l = l-1;if l == 0;l = onFrames;end
								if num_rings >=5
									m = m-1;if m == 0;m = onFrames;end
									if num_rings >=6
										n = n-1;if n == 0;n = onFrames;end
										if num_rings >= 7
											o = o-1;if o == 0;o = onFrames;end
											if num_rings >= 8
												p = p -1;if p==0;p = onFrames;end
												if num_rings >=9
													q = q-1;if q == 0;q = onFrames;end
													if num_rings >=10
														r = r-1;if r == 0;r = onFrames;end
														if num_rings >=11
															s = s-1;if s == 0;s = onFrames;end
															if num_rings >=12
																t = t-1;if t == 0;t = onFrames;end
															end %12
														end %11
													end %10
												end %9
											end %8
										end %7
									end %6
								end %5
							end %4
						end %3
					end %2
				end %approach
				
				%                 if approach == 2   %  ֻΪ?��?�ֹ�?�?�?�?��?�?��?�?�?�?��?�?�?
				ii = ii+1;				jj = jj+1;
				kk = kk+1;				ll = ll+1;
				mm = mm+1;				nn = nn+1;
				oo = oo +1;				pp = pp+1;
				qq = qq +1;				rr = rr +1;
				ss = ss +1;				tt = tt +1;
				%                 end
				
				if useEyeLink
					getSample(eL);
					isfix = isFixated(eL);
					if ~isfix
						fixated = 'breakfix';
						break %break the while loop
					end
				end
				[vbl,~,~,missed] = Screen('Flip', w, vbl + ana.halfifi);
				if missed>0 && isempty(windowed); fprintf('---!!! Missed frame !!!---\n'); end
				
				%------escape li
				[~, ~, keycode] = KbCheck(-1);
				if keycode(esc)
					breakLoop = true;
					break
				end
				%----
			end % end main vblendtime draw loop
			
			sM.drawBackground();
			Screen('Flip',w);
			thisTrial = iii;
			% check if we lost fixation
			if ~strcmpi(fixated,'fix')
				fprintf('===>>> BROKE FIXATION Trial = %i\n', iii);
				response = -1;
				statusMessage(eL,'Subject Broke Fixation!');
				edfMessage(eL,'MSG:BreakFix');
				resetFixation(eL);
				stopRecording(eL);
				edfMessage(eL,['TRIAL_RESULT ' num2str(response)]);
				setOffline(eL);
				continue
			end
			
			% get subject reponse
			if match_value == 1
				Pinna_grating_match(condition,match_time);
			else
				Pinna_grating_match3(condition,match_time);
			end
			if abandon == 0   %else restart,and don't save data
				ana.result(1,iii) = numChoice;
				response = numChoice;
				%---------------------update staircase data
				if use_staircase
					if pattern_angle_index(iii) == 1
						staircase45 = PAL_AMPM_updatePM(staircase45, response);
					elseif pattern_angle_index(iii) == 2
						staircase90 = PAL_AMPM_updatePM(staircase90, response);
					elseif pattern_angle_index(iii) == 3
						staircase135 = PAL_AMPM_updatePM(staircase135, response);
					end
				end
				doPlot();
				iii = iii+1;
			end
			fprintf('==>> Trial = %i RESPONSE = %.2g\n', thisTrial,response);
			%log response to eyelink
			if useEyeLink
				resetFixation(eL);
				stopRecording(eL);
				edfMessage(eL,['TRIAL_RESULT ' num2str(response)]);
				setOffline(eL);
			end
			%------------escape li
			if breakLoop
				breakLoop = true;
				break
			end
		end % END iii <= trials
		if iii>trials; breakLoop = true; end
	end % breakLoop
	close(sM);
	ListenChar(0);ShowCursor;Priority(0);Screen('CloseAll');
	if exist(ResultDir,'dir')
		cd(ResultDir);
	end
	disp(['==>> SAVE, saved current data to: ' pwd]);
	if use_staircase
		save([nameExp '.mat'],'ana','staircase90','staircase45','staircase135','eL');
	else
		save([nameExp '.mat'],'ana','eL');
	end
	savefig(figH, [nameExp '.fig']);
	if useEyeLink == true; close(eL); end
catch ME
	close(sM);
	ListenChar(0);ShowCursor;Priority(0);Screen('CloseAll');
	getReport(ME)
end

	function saveMetaData()
		ana.match_value = match_value;
		ana.angle_pattern = angle_pattern;
		ana.angle_index = pattern_angle_index;
		if ~use_staircase
			ana.move_speed_index = radial_speed_index;
			ana.angle_speed_index = rotation_speed_index;
		end
		ana.randArray = column_rank;
% 		ana.result = zeros(1,trials);
		ana.allAngle = gaborAngles;
		ana.one_trials = one_trials;
		ana.moving_speed = radial_speed_i;
		ana.angle_speed = rotation_speed_i;
		ana.sti_duration = duration;
		ana.sti_match_time = match_time;
		ana.is_binary_mask = is_binary_mask;
		ana.mask_diameter = mask_diameter;
		ana.mask_xpos = mask_xpos;
		ana.mask_ypos = mask_ypos;
		ana.ResultDir =  ResultDir;
		ana.backgroundColor = backgroundColor;
		ana.stdDis = stdDis;
		ana.xCen = xCen;
		ana.yCen = yCen;
		ana.num_rings = num_rings;
		ana.numElements = numElements;
		
		ana.useEyeLink = useEyeLink;
		ana.isDummy = isDummy;
		ana.pixelsPerCm = pixelsPerCm; %26 for Dorris lab,32=Lab CRT -- 44=27"monitor or Macbook Pro
		ana.distance = distance; %64.5 in Dorris lab;
		ana.windowed = windowed;
		ana.eL = eL;
		
		ana.fixX = fixX;
		ana.fixY = fixY;
		ana.firstFixInit = firstFixInit;
		ana.firstFixTime = firstFixTime;
		ana.firstFixRadius = firstFixRadius;
		ana.strictFixation = strictFixation;
	end

	function setupTrial()
		
		stopCriterion = 'trials';
		trials = 90;
		stopRule = 30;
		usePriors = false;
		grain = 100;
		
		% 		task = stimulusSequence();
		% 		task.name = nameExp;
		% 		task.nVar(1).name = 'angle';
		% 		task.nVar(1).stimulus = [1];
		% 		task.nVar(1).values = [45 90 135];
		% 		task.nBlocks = length(task.nVar(1).values) * stopRule;
		% 		initialise(task);
		
		if staircase_use_radial
			stims = linspace(min(radial_speed_i),max(radial_speed_i),grain);
			priorAlpha90 = linspace(min(radial_speed_i), max(radial_speed_i),grain);
			priorAlpha45 = linspace(min(radial_speed_i), max(radial_speed_i),grain);
			priorAlpha135 = linspace(min(radial_speed_i), max(radial_speed_i),grain);
			priorBeta = linspace(0.1, 25,grain); %our slope
			priorGammaRange = 0.02;  %fixed value (using vector here would make it a free parameter)
			priorLambdaRange = 0.02; %ditto
		else
			stims = linspace(min(rotation_speed_i),max(rotation_speed_i),grain);
			priorAlpha90 = linspace(min(rotation_speed_i),max(rotation_speed_i),grain);
			priorAlpha45 = linspace(min(rotation_speed_i),max(rotation_speed_i),grain);
			priorAlpha135 = linspace(min(rotation_speed_i),max(rotation_speed_i),grain);
			priorBeta = linspace(0.1, 25,grain); %our slope
			priorGammaRange = 0.02;  %fixed value (using vector here would make it a free parameter)
			priorLambdaRange = 0.02; %ditto
		end
		
		staircase90 = PAL_AMPM_setupPM('stimRange',stims,'PF',@PAL_Logistic,...
			'priorAlphaRange', priorAlpha90, 'priorBetaRange', priorBeta,...
			'priorGammaRange',priorGammaRange, 'priorLambdaRange',priorLambdaRange,...
			'numTrials', stopRule,'marginalize','lapse');
		
		staircase45 = PAL_AMPM_setupPM('stimRange',stims,'PF',@PAL_Logistic,...
			'priorAlphaRange', priorAlpha45, 'priorBetaRange', priorBeta,...
			'priorGammaRange',priorGammaRange, 'priorLambdaRange',priorLambdaRange,...
			'numTrials', stopRule,'marginalize','lapse');
		
		staircase135 = PAL_AMPM_setupPM('stimRange',stims,'PF',@PAL_Logistic,...
			'priorAlphaRange', priorAlpha135, 'priorBetaRange', priorBeta,...
			'priorGammaRange',priorGammaRange, 'priorLambdaRange',priorLambdaRange,...
			'numTrials', stopRule,'marginalize','lapse');
		
		if usePriors
			prior90 = PAL_pdfNormal(staircase90.priorAlphas,0,20).*PAL_pdfNormal(staircase90.priorBetas,2,1);
			prior45 = PAL_pdfNormal(staircase45.priorAlphas,0,20).*PAL_pdfNormal(staircase45.priorBetas,2,1);
			prior135 = PAL_pdfNormal(staircase135.priorAlphas,0,20).*PAL_pdfNormal(staircase135.priorBetas,2,1);
% 			figure;
% 			subplot(1,3,1);imagesc(staircase90.priorBetaRange,staircase90.priorAlphaRange,prior90);axis square
% 			ylabel('Threshold');xlabel('Slope');title('Initial Bayesian Priors 90deg')
% 			subplot(1,3,2);imagesc(staircase45.priorBetaRange,staircase45.priorAlphaRange,prior45); axis square
% 			ylabel('Threshold');xlabel('Slope');title('Initial Bayesian Priors 45deg')
% 			subplot(1,3,3);imagesc(staircase135.priorBetaRange,staircase135.priorAlphaRange,prior135); axis square
% 			ylabel('Threshold');xlabel('Slope');title('Initial Bayesian Priors 135deg')
			staircase90 = PAL_AMPM_setupPM(staircase90,'prior',prior90);
			staircase45 = PAL_AMPM_setupPM(staircase45,'prior',prior45);
			staircase135 = PAL_AMPM_setupPM(staircase135,'prior',prior135);
		end
	end

	function doPlot()
		if ~isfield(ana,'result') || isempty(ana.result)
			return
		end
		if use_staircase
			figure(figH);
			subplot(3,1,1)
			hold on
			if pattern_angle_index(iii) == 1
				plot(iii,ana.result(end)+4,'ro','MarkerFaceColor','r','MarkerSize',8);
			elseif pattern_angle_index(iii) == 2
				plot(iii,ana.result(end)+2,'go','MarkerFaceColor','g','MarkerSize',8);
			elseif pattern_angle_index(iii) == 3
				plot(iii,ana.result(end),'bo','MarkerFaceColor','b','MarkerSize',8);
			end
			tit = sprintf('NEXT TRIAL:%i', iii);
			title(tit);
			box on; grid on; grid minor; ylim([0 6])
			xlabel('Total Trials');
			ylabel('Subject Response');
			hold off;
			subplot(3,1,2)
			if staircase_use_radial
				if length(staircase45.x)>1 && length(staircase90.x)>1 && length(staircase135.x)>1
					[SL45, NP45, OON45] = PAL_PFML_GroupTrialsbyX(staircase45.x(1:length(staircase45.x)-1),staircase45.response,ones(size(staircase45.response)));
					for SR45 = 1:length(SL45(OON45~=0))
						plot(SL45(SR45),NP45(SR45)/OON45(SR45),'ro','markerfacecolor','r','markersize',20*sqrt(OON45(SR45)./sum(OON45)))
						hold on
					end
					[SL90, NP90, OON90] = PAL_PFML_GroupTrialsbyX(staircase90.x(1:length(staircase90.x)-1),staircase90.response,ones(size(staircase90.response)));
					for SR90 = 1:length(SL90(OON90~=0))
						plot(SL90(SR90),NP90(SR90)/OON90(SR90),'go','markerfacecolor','g','markersize',20*sqrt(OON90(SR90)./sum(OON90)))
						hold on
					end
					[SL135, NP135, OON135] = PAL_PFML_GroupTrialsbyX(staircase135.x(1:length(staircase135.x)-1),staircase135.response,ones(size(staircase135.response)));
					for SR135 = 1:length(SL135(OON135~=0))
						plot(SL135(SR135),NP135(SR135)/OON135(SR135),'bo','markerfacecolor','b','markersize',20*sqrt(OON135(SR135)./sum(OON135)))
						hold on
					end
					plot(min(radial_speed_i):0.01:max(radial_speed_i),...
						staircase45.PF([staircase45.threshold(length(staircase45.threshold)), staircase45.slope(length(staircase45.threshold)), 0.02, 0.02],...
						min(radial_speed_i):0.01:max(radial_speed_i)),'r-','linewidth',2)
					hold on
					plot(min(radial_speed_i):0.01:max(radial_speed_i),...
						staircase90.PF([staircase90.threshold(length(staircase90.threshold)), staircase90.slope(length(staircase90.threshold)), 0.02, 0.02],...
						min(radial_speed_i):0.01:max(radial_speed_i)),'g-','linewidth',2)
					hold on
					plot(min(radial_speed_i):0.01:max(radial_speed_i),...
						staircase135.PF([staircase135.threshold(length(staircase135.threshold)), staircase135.slope(length(staircase135.threshold)), 0.02, 0.02],...
						min(radial_speed_i):0.01:max(radial_speed_i)),'b-','linewidth',2)
					hold off
				end
			else
				if length(staircase45.x)>1 && length(staircase90.x)>1 && length(staircase135.x)>1
					[SL45, NP45, OON45] = PAL_PFML_GroupTrialsbyX(staircase45.x(1:length(staircase45.x)-1),staircase45.response,ones(size(staircase45.response)));
					for SR45 = 1:length(SL45(OON45~=0))
						plot(SL45(SR45),NP45(SR45)/OON45(SR45),'ro','markerfacecolor','r','markersize',20*sqrt(OON45(SR45)./sum(OON45)))
						hold on
					end
					[SL90, NP90, OON90] = PAL_PFML_GroupTrialsbyX(staircase90.x(1:length(staircase90.x)-1),staircase90.response,ones(size(staircase90.response)));
					for SR90 = 1:length(SL90(OON90~=0))
						plot(SL90(SR90),NP90(SR90)/OON90(SR90),'go','markerfacecolor','g','markersize',20*sqrt(OON90(SR90)./sum(OON90)))
						hold on
					end
					[SL135, NP135, OON135] = PAL_PFML_GroupTrialsbyX(staircase135.x(1:length(staircase135.x)-1),staircase135.response,ones(size(staircase135.response)));
					for SR135 = 1:length(SL135(OON135~=0))
						plot(SL135(SR135),NP135(SR135)/OON135(SR135),'bo','markerfacecolor','b','markersize',20*sqrt(OON135(SR135)./sum(OON135)))
						hold on
					end
					plot(min(rotation_speed_i):0.01:max(rotation_speed_i),...
						staircase45.PF([staircase45.threshold(length(staircase45.threshold)), staircase45.slope(length(staircase45.threshold)), 0.02, 0.02],...
						min(rotation_speed_i):0.01:max(rotation_speed_i)),'r-','linewidth',2)
					hold on
					plot(min(rotation_speed_i):0.01:max(rotation_speed_i),...
						staircase90.PF([staircase90.threshold(length(staircase90.threshold)), staircase90.slope(length(staircase90.threshold)), 0.02, 0.02],...
						min(rotation_speed_i):0.01:max(rotation_speed_i)),'g-','linewidth',2)
					hold on
					plot(min(rotation_speed_i):0.01:max(rotation_speed_i),...
						staircase135.PF([staircase135.threshold(length(staircase135.threshold)), staircase135.slope(length(staircase135.threshold)), 0.02, 0.02],...
						min(rotation_speed_i):0.01:max(rotation_speed_i)),'b-','linewidth',2)
					hold off
				end
			end
			subplot(3,1,3)
			if length(staircase90.threshold) > 1
				analysisCore.areabar(1:length(staircase90.threshold),staircase90.threshold,staircase90.seThreshold,[0.5 1 0.5],[],'g-');
				plot(1:length(staircase90.x),staircase90.x,'g');
			end
			if length(staircase45.threshold) >1
				hold on
				analysisCore.areabar(1:length(staircase45.threshold),staircase45.threshold,staircase45.seThreshold,[1 0.5 0.5],[],'r-');
				plot(1:length(staircase45.x),staircase45.x,'r');
			end
			if length(staircase135.threshold) > 1
				analysisCore.areabar(1:length(staircase135.threshold),staircase135.threshold,staircase135.seThreshold,[0.5 0.5 1],[],'b-');
				plot(1:length(staircase135.x),staircase135.x,'b');
			end
			hold off
			drawnow
		else
			figure(figH);
			xpl = 1:length(ana.result);
			plot(xpl, ana.result,'ro','MarkerFaceColor','r','MarkerSize',8);
			tit = sprintf('NEXT TRIAL:%i', iii);
			title(tit);
			box on; grid on; grid minor; ylim([0 3])
			xlabel('Total Trials');
			ylabel('Subject Response');
			hold off;
			drawnow;
		end
	end
end
