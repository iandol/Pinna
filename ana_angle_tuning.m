function ana_angle_tuning(ResultDir,CommonResultDir)
global ana prefer_angle_data name
%%%%
index = strfind(ResultDir,'\');
index_end = index(1,end);
student_name = ResultDir(index_end+1:end);  %can use chinese
%-----------
cd(ResultDir);
cwd=pwd;
load('data.mat')
allAngle1 =  ana.allAngle;
% angSpeed1 = 0:5:30;
one_trials = ana.one_trials;
% trials = one_trials * size(allAngle1,1);
ana.randAngleArray = abs(ana.randAngleArray);
num_angle = size(allAngle1,2);
all_speed = zeros(one_trials,num_angle);
for i = 1:num_angle
    index = find(ana.randAngleArray == allAngle1(i));
    angle_trials = length(index);
    speed = ana.result(index);
    all_speed(:,i) = speed';
    ave_speed(i) = sum(speed)/angle_trials;
end
speed_err = std(all_speed,0,1)/sqrt(angle_trials);
figure;
h1 = errorbar(allAngle1, ave_speed, speed_err, 'o-','MarkerSize',5);
set(h1,'color','b','linewidth',2)
% plot(allAngle1,ave_speed,'linewidth',2);

xlabel('angle of inclination ','fontsize',10);
ylabel('real angle speed ','fontsize',10);
ylim([0 35]);
legend(student_name,1);

max_speed = max(ave_speed);
prefer_angle_index = find(ave_speed == max_speed);
prefer_angle = allAngle1(prefer_angle_index');
number = length(prefer_angle);
max_speed1 = [];
for i=1:number
   max_speed1 = [max_speed1,max_speed]; 
end
A= [prefer_angle;max_speed1];
hold on;
for i =1:number
 plot([prefer_angle(i) prefer_angle(i)],[0 max_speed],'r--');
 plot([0 prefer_angle(i)],[max_speed max_speed],'r--');
end
saveas(gcf,strcat('angle_tuning_curve'), 'fig')


%%%%%%%%%%%population
% fileCommonPath = 'E:\PinnaGratingPreferAngleTrends';
cd(CommonResultDir);
% load('prefer_angle_data.mat');
%%%%------------------------
if isempty(prefer_angle_data)  %first save data 
    for i = 1:number
       name{i} =  student_name;
    end
%    figure;
else
    for i = 1:number
      name = [name,student_name];
    end
%     open prefer_angle_trend.fig;
end
%--------------------------

prefer_angle_data = [prefer_angle_data,A];  %row 1 is angle data;row 2 is speed data;
save([CommonResultDir '\' 'prefer_angle_data'],'prefer_angle_data','name');
figure;
rgb= [0 0 0;1 0 0;0 1 0;0 0 1;1 1 0;1 0 1;0 1 1;1 1 1;0.5 0 0;0 0.5 0;0 0 0.5;0.5 0.5 0;0.5 0 0.5;0 0.5 0.5;
      0.5 0.5 0.5;1 0.5 0;1 0 0.5;1 0.5 0.5;0.5 1 0;0 1 0.5;0.5 1 0.5;0 0.5 1;0.5 0 1;0.5 0.5 1;1 1 0.5;1 0.5 1;
      0.5 1 1;1 0 0;1 0 0;1 0 0;1 0 0;1 0 0;1 0 0;1 0 0;1 0 0;1 0 0;1 0 0;1 0 0;1 0 0;1 0 0;1 0 0];   %27+14
for j = 1:length(name)
   hold on;
   plot(prefer_angle_data(1,j),prefer_angle_data(2,j),'p','Color',rgb(j,:),'MarkerSize',10);
end
legend(name,-1);
xlabel('prefer angle','fontsize',10);
ylabel('max real speed ','fontsize',10);
ylim([0 35]);
xlim([0 90]);
saveas(gcf,strcat('prefer_angle_trend'), 'fig')
cd(cwd);
