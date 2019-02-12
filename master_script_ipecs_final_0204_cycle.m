%This script works with both Comparison and Prototype data FOR SUBJECTS AA02 AND AA04 ONLY, but only MATLAB
%data from ONE of the two protheses may be in the same path as the script at a
%time.

close
clc
clearvars -except scalef offset %clears all variables from previous iteration, except for calibration that carries over
cal=input('Do you need to run 2-point calibration? Press "y" or "n" ','s');
side=input('Is the prothesis on the right or left? Press "r" for right, "l" for left ','s'); %some axes definitions will need to be flipped depending on the input here
weight=input('What is the subjects weight in kg? ','s'); %used to normalize moment data by weight of subject

%% Determination of Offsets and Scaling Factors
if strcmp(cal,'y')==1
    load iPECS_Cal.mat %loading calibration file

    force_cal=iPECS_Cal.Analog.Data(13:20,1:20000); %pull data from necessary 8 channels

    for n=1:8 %sort channel data by the low or high voltage input
        force_cal_low=[];
        force_cal_high=[];
        
        for j=1:length(force_cal)
            if force_cal(n,j)<-2 %sorting by picking threshold (+/- 2V) 
                force_cal_low=[force_cal_low; force_cal(n,j)];
            elseif force_cal(n,j)>2
                force_cal_high=[force_cal_high; force_cal(n,j)];
            end
        end

        low_avg(n)=mean(force_cal_low); %average the low voltage readings
        high_avg(n)=mean(force_cal_high); %average high readings
        plot(force_cal(n,:),'k') %ploting calibration data, followed by the low and high voltage values determined by the script for a visual check
        hold on
        plot([1 length(force_cal)],[low_avg(n) low_avg(n)],'r')
        plot([1 length(force_cal)],[high_avg(n) high_avg(n)],'b')
        title(num2str(n))
        g=axis;
        axis([g(1) g(2) -3 3])
        pause %allows user to let the script proceed only if the calibration values are determined to be accurate
        clf
        clearvars force_cal_low force_cal_high %clear high/low voltages to allow for subsequent channel to overwrite them

        scalef(n)=2.45/(low_avg(n)-high_avg(n));%calibration formula given in iPECS manual
        offset(n)=1.225-(low_avg(n)*scalef(n)); %offset formula given in manual
    end
end

%% Calculation of Forces and Moments
D=dir('*.mat'); %create directory for each condition's matlab file
% change this to user-prompted

cal_matrix=[35.693526 -19.833401 39.455188 -278.22993 37.023151 211.599345 35.277711 7.139588; %establish calibration matrix from iPECS
    -4.128547 11.220743 -0.169142 24.752593 233.393184 17.232101 -274.45351 6.098171;
    69.204972 -5.297664 79.234431 8.532263 -25.149022 22.450788 -9.112311 688.268961;
    1108.95636 1134.2613 -44.205883 45.72884 -84.968904 46.387532 -17.059699 274.634798;
    115.530562 -1138.4244 -977.89314 -101.16652 72.40158 -34.431581 74.364339 -252.79711;
    -37.491058 9.286816 -27.48374 177.715892 129.73555 177.496754 126.394302 -12.441065];

    scalef_mat=[scalef(1) 0 0 0 0 0 0 0; %create scaling matrix with each channel's scaling factor
        0 scalef(2) 0 0 0 0 0 0;
        0 0 scalef(3) 0 0 0 0 0;
        0 0 0 scalef(4) 0 0 0 0;
        0 0 0 0 scalef(5) 0 0 0;
        0 0 0 0 0 scalef(6) 0 0;
        0 0 0 0 0 0 scalef(7) 0;
        0 0 0 0 0 0 0 scalef(8)];
    
   offset_mat=[offset(1) 0 0 0 0 0 0 0; %create offset matrix with each channel's offset
        0 offset(2) 0 0 0 0 0 0;
        0 0 offset(3) 0 0 0 0 0;
        0 0 0 offset(4) 0 0 0 0;
        0 0 0 0 offset(5) 0 0 0;
        0 0 0 0 0 offset(6) 0 0;
        0 0 0 0 0 0 offset(7) 0;
        0 0 0 0 0 0 0 offset(8)]; 

%[d,c] = butter(4,(50*0.89568)/(1200/2),'low'); % Low pass filter with
%effective order=8, cutoff at 50 Hz. This has been turned off for now, but
%left in the code in case it is determined to be useful by the user.

for filei=1:length(D) %run for loop for each condition
    
    filename=D(filei).name; %load condition, store data, clear data matrix after pulling data from it
    eval(['load ' filename]);
    eval(['data = ' filename(1:(length(filename)-4)) ';']);
    eval(['clear ' filename(1:(length(filename)-4)) ';']);
    
    voltage=data.Analog.Data(13:20,:); %pull voltage data from 8 channels of respective condition
    
    voltage=voltage.'; %transpose for easier calculation
    
    offset_new=ones(size(voltage))*offset_mat; %augment offset matrix
    IPECS_measure=cal_matrix*(scalef_mat*(voltage)'+offset_new'); %find measured forces/moments in lbs and in*lbs for each channel
    
    Fipecs_filt=IPECS_measure(1:3,:)'; %separating forces into separate array    
    Mipecs_filt=IPECS_measure(4:6,:)'; %separating moments into separate array
    
    %The next few lines filter all force and moment data, but again have
    %been turned off at this time.
    
%     Fipecs_filt(:,1)=filtfilt(d,c,Fipecs_filt(:,1));
%     Fipecs_filt(:,2)=filtfilt(d,c,Fipecs_filt(:,2));
%     Fipecs_filt(:,3)=filtfilt(d,c,Fipecs_filt(:,3)); 
%     Mipecs_filt(:,1)=filtfilt(d,c,Mipecs_filt(:,1)); 
%     Mipecs_filt(:,2)=filtfilt(d,c,Mipecs_filt(:,2));
%     Mipecs_filt(:,3)=filtfilt(d,c,Mipecs_filt(:,3));
    
    Values{filei,1}=filename; %creating cell to hold forces/moments for each condition
    Values{filei,2}=Fipecs_filt*4.44822; %converting from pounds (lb) to newtons (N)
    Values{filei,3}=Mipecs_filt*0.112984829; %converting from in*lb to N*m
    
    clear voltage Fipecs_filt Mipecs_filt IPECS_measure offset_new %clearing all unnessary variables that will be written over when analyzing the next channel
        
end

%These next few lines remove files that will not be analyzed, the names of
%which differ depending on if the prothesis is Comparison or Prototype.
%AnkleStatic data is searched for in both cases, but if it is not present
%in the set of files being analzyed, it will not affect the code.
    
r=find(strncmp('Ramp',Values(:,1),4));
Values(r,:)=[];
D(r)=[];
s=find(strncmp('Static',Values(:,1),4));
Values(s,:)=[];
D(s)=[];
ic=find(strncmp('iPECS',Values(:,1),4));
Values(ic,:)=[];
D(ic)=[];
as=find(strncmp('Ankle',Values(:,1),4));
Values(as,:)=[];
D(as)=[];

clear r s ic as

%% Determination of Ankle/Knee Markers for Comparison Foot
%similar to above, this section will only run for the Comparison foot as
%the method of determining ankle/knee markers is different

load Static.mat %loading static trial

%removing unnessary fourth value from each vector
Static.Trajectories.Labeled.Data(:,4,:)=[]; 

%for the Comparison foot, one coordinate system can be used to find all
%necessary markers. This coordinate system will be based on the shank,
%and thus will be heavily dependant on the side of the prothesis.
if strcmp(side,'r')==1
    lk=find(strcmp('R_Knee_Lat',Static.Trajectories.Labeled.Labels)); %finding indices of important markers that are side-dependant
    mk=find(strcmp('R_Knee_Med',Static.Trajectories.Labeled.Labels));
    la=find(strcmp('R_Ank_Lat',Static.Trajectories.Labeled.Labels));
    ma=find(strcmp('R_Ank_Med',Static.Trajectories.Labeled.Labels));
    s=find(strcmp('R_Shank',Static.Trajectories.Labeled.Labels));
elseif strcmp(side,'l')==1
    lk=find(strcmp('L_Knee_Lat',Static.Trajectories.Labeled.Labels)); 
    mk=find(strcmp('L_Knee_Med',Static.Trajectories.Labeled.Labels));
    la=find(strcmp('L_Ank_Lat',Static.Trajectories.Labeled.Labels));
    ma=find(strcmp('L_Ank_Med',Static.Trajectories.Labeled.Labels));
    s=find(strcmp('L_Shank',Static.Trajectories.Labeled.Labels)); 
end

im=find(strcmp('iPECS_Med',Static.Trajectories.Labeled.Labels)); %finding indices of side-independant markers
il=find(strcmp('iPECS_Lat',Static.Trajectories.Labeled.Labels));
ip=find(strcmp('iPECS_Post',Static.Trajectories.Labeled.Labels));

for n=1:length(Static.Trajectories.Labeled.Data)  
    %creating intermediate shank coordinate system defined per the
    %specifications in the OrthoTrak manual. As none of the axes
    %depend on the relationship between medial and lateral markers,
    %axes do not have to be flipped in definiition depending on the
    %side of the prosthesis.
    lat_knee_lat_shank_vec(n,:)=normalize(Static.Trajectories.Labeled.Data(s,:,n)-Static.Trajectories.Labeled.Data(lk,:,n));
    shank_inter_x(n,:)=normalize(Static.Trajectories.Labeled.Data(la,:,n)-Static.Trajectories.Labeled.Data(lk,:,n)); ... %X axis of coord. system
    shank_inter_y(n,:)=cross(lat_knee_lat_shank_vec(n,:),shank_inter_x(n,:)); %Y axis found by same method as previous coord. system(s)
    shank_inter_z(n,:)=cross(shank_inter_x(n,:),shank_inter_y(n,:)); %Z axis found by same method as previous coord. system(s)

    med_ankle(n,:)=Static.Trajectories.Labeled.Data(ma,:,n); %defining important marker positions and placing them in their own arrays
    lat_ankle(n,:)=Static.Trajectories.Labeled.Data(la,:,n);
    lat_knee(n,:)=Static.Trajectories.Labeled.Data(lk,:,n);
    med_knee(n,:)=Static.Trajectories.Labeled.Data(mk,:,n);

    shank_inter_x(n,:)=normalize(shank_inter_x(n,:)); %normalizing vectors (again) for use in matrices
    shank_inter_y(n,:)=normalize(shank_inter_y(n,:));
    shank_inter_z(n,:)=normalize(shank_inter_z(n,:));

    med_ankle_t(:,n)=[med_ankle(n,1);med_ankle(n,2);med_ankle(n,3); 1]; %augmenting and transposing important markers for use in transformation matrix
    lat_ankle_t(:,n)=[lat_ankle(n,1);lat_ankle(n,2);lat_ankle(n,3); 1];
    ankle_cent_t(:,n)=(lat_ankle_t(:,n)+med_ankle_t(:,n))./2;
    lat_knee_t(:,n)=[lat_knee(n,1);lat_knee(n,2);lat_knee(n,3); 1];
    med_knee_t(:,n)=[med_knee(n,1);med_knee(n,2);med_knee(n,3); 1];
    knee_cent_t(:,n)=(lat_knee_t(:,n)+med_knee_t(:,n))./2;

    T_inter_shank(:,:,n)=[shank_inter_x(n,1) shank_inter_y(n,1) shank_inter_z(n,1) lat_knee_t(1,n); %forming intermediate shank coord. system transformation matrix, where the lat knee marker is the origin
                    shank_inter_x(n,2) shank_inter_y(n,2) shank_inter_z(n,2) lat_knee_t(2,n); 
                    shank_inter_x(n,3) shank_inter_y(n,3) shank_inter_z(n,3) lat_knee_t(3,n);
                    0 0 0 1];

    med_ankle_inter(:,n)=T_inter_shank(:,:,n)\med_ankle_t(:,n); %finding important marker positions in intermediate shank coord. system
    lat_ankle_inter(:,n)=T_inter_shank(:,:,n)\lat_ankle_t(:,n);
    med_knee_inter(:,n)=T_inter_shank(:,:,n)\med_knee_t(:,n);
    lat_knee_inter(:,n)=T_inter_shank(:,:,n)\lat_knee_t(:,n);

 %finding ankle and knee centers in intermediate coord. system; to
 %debug, these positions can be compared against ankle/knee centers
 %found by calculating the midpoint of the med-lat markers in lab.
 %system
    ankle_center_inter(:,n)=T_inter_shank(:,:,n)\ankle_cent_t(:,n);
    knee_center_inter(:,n)=T_inter_shank(:,:,n)\knee_cent_t(:,n);
end

med_ankle_inter=mean(med_ankle_inter,2); %averaging inter. shank coord. system positions of important markers, as they should not be moving regardless
lat_ankle_inter=mean(lat_ankle_inter,2);
med_knee_inter=mean(med_knee_inter,2);
lat_knee_inter=mean(lat_knee_inter,2);

ankle_center_inter=mean(ankle_center_inter,2);
knee_center_inter=mean(knee_center_inter,2);

%% Translation to Dynamic Trials
%This section is ran by both protheses, but various sections will either be
%ran or skipped over depending on the prothesis. The coordinate systems
%defined in this section are identical to the systems previously defined
%(unless noted), so their definitions will not be explained.
 
for filei=1:length(D) %run for loop for each condition
    
    filename=D(filei).name; %load condition, store data, clear data matrix after pulling data from it
    eval(['load ' filename]);
    eval(['data = ' filename(1:(length(filename)-4)) ';']);
    eval(['clear ' filename(1:(length(filename)-4)) ';']);

    data.Trajectories.Labeled.Data(:,4,:)=[]; %removing unnessessary fourth data value

    %clearing marker positions remaining from previous conditions, 
    %as they will need to be overwritten by their positions in the dynamic 
    %trials--only intermediate marker positions should remain
    clearvars -except scalef offset D Values weight cal side filei data ankle_center_inter knee_center_inter ... 
        lat_ankle_inter lat_knee_inter med_ankle_inter med_knee_inter
        
    if strcmp(side,'r')==1
        al=find(strcmp('R_Ank_Lat',data.Trajectories.Labeled.Labels)); %finding indices of important markers, the side of which will depend on the side of the prothesis 
        kl=find(strcmp('R_Knee_Lat',data.Trajectories.Labeled.Labels));
        s=find(strcmp('R_Shank',data.Trajectories.Labeled.Labels));
        h=find(strcmp('R_Heel',data.Trajectories.Labeled.Labels));
        t=find(strcmp('R_Toe',data.Trajectories.Labeled.Labels));
        im=find(strcmp('iPECS_Med',data.Trajectories.Labeled.Labels));
        il=find(strcmp('iPECS_Lat',data.Trajectories.Labeled.Labels));
        ip=find(strcmp('iPECS_Post',data.Trajectories.Labeled.Labels));
        ala=find(strcmp('R_Ank_Lat_Ant',data.Trajectories.Labeled.Labels));
        alp=find(strcmp('R_Ank_Lat_Post',data.Trajectories.Labeled.Labels));
        ama=find(strcmp('R_Ank_Med_Ant',data.Trajectories.Labeled.Labels));
        amp=find(strcmp('R_Ank_Med_Post',data.Trajectories.Labeled.Labels));
    elseif strcmp(side,'l')==1
        al=find(strcmp('L_Ank_Lat',data.Trajectories.Labeled.Labels)); 
        kl=find(strcmp('L_Knee_Lat',data.Trajectories.Labeled.Labels));
        s=find(strcmp('L_Shank',data.Trajectories.Labeled.Labels));
        h=find(strcmp('L_Heel',data.Trajectories.Labeled.Labels));
        t=find(strcmp('L_Toe',data.Trajectories.Labeled.Labels));
        im=find(strcmp('iPECS_Med',data.Trajectories.Labeled.Labels));
        il=find(strcmp('iPECS_Lat',data.Trajectories.Labeled.Labels));
        ip=find(strcmp('iPECS_Post',data.Trajectories.Labeled.Labels));
        ala=find(strcmp('L_Ank_Lat_Ant',data.Trajectories.Labeled.Labels));
        alp=find(strcmp('L_Ank_Lat_Post',data.Trajectories.Labeled.Labels));
        ama=find(strcmp('L_Ank_Med_Ant',data.Trajectories.Labeled.Labels));
        amp=find(strcmp('L_Ank_Med_Post',data.Trajectories.Labeled.Labels));
    end
    
    %This newly added section will hopefully correct for marker dropout by
    %finding where the marker drops and removing that section from ALL
    %marker positions.
    markers=[al kl s h t im il ip ala alp ama amp];
    j=0;
    for k=1:length(markers)
        bad=isnan(data.Trajectories.Labeled.Data(markers(k),1,:));
        badf=find(bad);
        badfs=size(badf);
        if badfs(1)==1
            badf=[];
        end
        if isempty(badf)~=1
            j=j+1;
            first=badf(1);
            last=badf(end);
            if j==1
                left_cut=first;
                right_cut=last;
            else
                if first<left_cut
                    left_cut=first;
                end
                if last>right_cut
                    right_cut=last;
                end
            end
        end
        clear bad badf badfs first last
    end
    clear j
    
    if exist('left_cut')==1
        data.Trajectories.Labeled.Data(:,:,(left_cut:right_cut))=[];
    end
    
    %as denoted above, this section contains the coord. systems used to
    %define ankle/knee positions for Comparison. The same coord. systems
    %will be used to find these marker positions in the lab coord. system
    %by the reverse process
    for n=1:length(data.Trajectories.Labeled.Data)
        lat_knee_lat_shank_vec(n,:)=normalize(data.Trajectories.Labeled.Data(s,:,n)-data.Trajectories.Labeled.Data(kl,:,n));
        shank_inter_x(n,:)=normalize(data.Trajectories.Labeled.Data(al,:,n)-data.Trajectories.Labeled.Data(kl,:,n)); 
        shank_inter_y(n,:)=cross(lat_knee_lat_shank_vec(n,:),shank_inter_x(n,:));
        shank_inter_z(n,:)=cross(shank_inter_x(n,:),shank_inter_y(n,:));

        lat_knee(n,:)=data.Trajectories.Labeled.Data(kl,:,n); %pulling lateral knee and ankle markers for comparison to calculated locations (for debugging purposes)
        lat_ankle(n,:)=data.Trajectories.Labeled.Data(al,:,n);

        shank_inter_x(n,:)=normalize(shank_inter_x(n,:)); 
        shank_inter_y(n,:)=normalize(shank_inter_y(n,:));
        shank_inter_z(n,:)=normalize(shank_inter_z(n,:));

        T_inter_shank(:,:,n)=[shank_inter_x(n,1) shank_inter_y(n,1) shank_inter_z(n,1) lat_knee(n,1); 
                    shank_inter_x(n,2) shank_inter_y(n,2) shank_inter_z(n,2) lat_knee(n,2);
                    shank_inter_x(n,3) shank_inter_y(n,3) shank_inter_z(n,3) lat_knee(n,3);
                    0 0 0 1];

        med_ankle_lab(n,:)=(T_inter_shank(:,:,n)*med_ankle_inter)'; %finding positions of each marker in lab frame by multiplying known intermediate position
        lat_ankle_lab(n,:)=(T_inter_shank(:,:,n)*lat_ankle_inter)'; %by transformation matrix that changes upon each new position
        med_knee_lab(n,:)=(T_inter_shank(:,:,n)*med_knee_inter)';
        lat_knee_lab(n,:)=(T_inter_shank(:,:,n)*lat_knee_inter)';

        ankle_center_lab(n,:)=(lat_ankle_lab(n,:)+med_ankle_lab(n,:))/2; %finding ankle and knee centers by finding midpoint of lateral and medial markers
        knee_center_lab(n,:)=(lat_knee_lab(n,:)+med_knee_lab(n,:))/2;

        ankle_center_lab_check(:,n)=T_inter_shank(:,:,n)*ankle_center_inter; %checking positions by transforming known joint center poition in intermediate frame
        knee_center_lab_check(:,n)=T_inter_shank(:,:,n)*knee_center_inter;

        %the iPECS coord. system is defined here, even though it is
        %not used to find ankle or knee markers, as the position of
        %the ankle center MUST be found in its coordinate system to
        %have accurate moment arm lengths (from iPECS center to ankle center) 
        %for force/moment calculations that need to be done later.

        iPECS_med(n,:)=data.Trajectories.Labeled.Data(im,:,n); 
        iPECS_lat(n,:)=data.Trajectories.Labeled.Data(il,:,n);
        iPECS_post(n,:)=data.Trajectories.Labeled.Data(ip,:,n);
        iPECS_mid(n,:)=(iPECS_med(n,:)+iPECS_lat(n,:))/2;

        if strcmp(side,'r')==1
            iPECS_inter_x(n,:)=normalize(iPECS_lat(n,:)-iPECS_med(n,:));
        else
            iPECS_inter_x(n,:)=normalize(iPECS_med(n,:)-iPECS_lat(n,:));
        end
        iPECS_y(n,:)=normalize(iPECS_mid(n,:)-iPECS_post(n,:));
        iPECS_z(n,:)=cross(iPECS_inter_x(n,:),iPECS_y(n,:));
        iPECS_x(n,:)=cross(iPECS_y(n,:),iPECS_z(n,:));

        iPECS_x(n,:)=normalize(iPECS_x(n,:));
        iPECS_y(n,:)=normalize(iPECS_y(n,:));
        iPECS_z(n,:)=normalize(iPECS_z(n,:));

        R_iPECS(:,:,n)=[iPECS_x(n,1) iPECS_y(n,1) iPECS_z(n,1); 
                    iPECS_x(n,2) iPECS_y(n,2) iPECS_z(n,2);
                    iPECS_x(n,3) iPECS_y(n,3) iPECS_z(n,3)];

        T_iPECS(:,:,n)=[iPECS_x(n,1) iPECS_y(n,1) iPECS_z(n,1) iPECS_mid(n,1); 
                    iPECS_x(n,2) iPECS_y(n,2) iPECS_z(n,2) iPECS_mid(n,2);
                    iPECS_x(n,3) iPECS_y(n,3) iPECS_z(n,3) iPECS_mid(n,3);
                    0 0 0 1];

    end

    med_ankle_lab(:,4)=[]; %removing unnessesary fourth value that is added when marker positions are found in lab. coord system
    lat_ankle_lab(:,4)=[];
    med_knee_lab(:,4)=[];
    lat_knee_lab(:,4)=[];
    ankle_center_lab(:,4)=[];
    knee_center_lab(:,4)=[];
    
    %this next section of code is run regardless of the prothesis being
    %analyzed, as the coordinate systems defined are used in later
    %calculations and can be defined on both protheses.
    for n=1:length(data.Trajectories.Labeled.Data) %this loop defines a new shank coordinate system that differs slightly from the conventions in OrthoTrak.
        if strcmp(side,'r')==1 %temporary Y axis is defined running med-lat or lat-med depending on side of prosthesis
            shank_prin_inter_y(n,:)=normalize(med_ankle_lab(n,:)-lat_ankle_lab(n,:));
        else
            shank_prin_inter_y(n,:)=normalize(lat_ankle_lab(n,:)-med_ankle_lab(n,:));
        end
        shank_principle_z(n,:)=normalize(knee_center_lab(n,:)-ankle_center_lab(n,:)); %Z axis defined vertically, running from ankle center to knee center
        shank_principle_x(n,:)=cross(shank_prin_inter_y(n,:),shank_principle_z(n,:)); %X axis defined running post-ant through cross product
        shank_principle_y(n,:)=cross(shank_principle_z(n,:),shank_principle_x(n,:)); %new Y axis created by cross product, matching coord. system conventions

        shank_principle_x(n,:)=normalize(shank_principle_x(n,:)); 
        shank_principle_y(n,:)=normalize(shank_principle_y(n,:));
        shank_principle_z(n,:)=normalize(shank_principle_z(n,:));

        %creating rotation matrix to convert forces and moments later in 
        %code by creating a matrix to move from iPECS coord. system to the shank coord. system
        R_shank(:,:,n)=[shank_principle_x(n,1) shank_principle_y(n,1) shank_principle_z(n,1);  
                    shank_principle_x(n,2) shank_principle_y(n,2) shank_principle_z(n,2);
                    shank_principle_x(n,3) shank_principle_y(n,3) shank_principle_z(n,3)];
                
        R_iPECS_shank(:,:,n)=R_iPECS(:,:,n)*R_shank(:,:,n);

        %This next piece of code creates an anatomical foot coordinate
        %system defined using OrthoTrak conventions. It is not used
        %anywhere else in the code, except in the "rpy_angles" function that
        %does not seem to produce agreeable data. For this reason, the
        %coordinate system has been turned off but left in the script in
        %case this method of finding ankle angles is developed further.
        
        toe(n,:)=data.Trajectories.Labeled.Data(t,:,n);
        heel(n,:)=data.Trajectories.Labeled.Data(h,:,n);
        
        anklecen_kneecen_vec(n,:)=normalize(knee_center_lab(n,:)-ankle_center_lab(n,:));
        foot_principle_y(n,:)=normalize(data.Trajectories.Labeled.Data(t,:,n)-data.Trajectories.Labeled.Data(h,:,n));
%         foot_principle_z(n,:)=cross(foot_principle_y(n,:),anklecen_kneecen_vec(n,:)); 
%         foot_principle_x(n,:)=cross(foot_principle_y(n,:),foot_principle_z(n,:)); 
% 
%         foot_principle_x(n,:)=normalize(foot_principle_x(n,:));
%         foot_principle_y(n,:)=normalize(foot_principle_y(n,:));
%         foot_principle_z(n,:)=normalize(foot_principle_z(n,:));
% 
%         R_foot(:,:,n)=[foot_principle_x(n,1) foot_principle_y(n,1) foot_principle_z(n,1); 
%                     foot_principle_x(n,2) foot_principle_y(n,2) foot_principle_z(n,2);
%                     foot_principle_x(n,3) foot_principle_y(n,3) foot_principle_z(n,3)];

%         R_shank_foot(:,:,n)=R_shank(:,:,n)*R_foot(:,:,n);
%         rpy_angles(n,:)=get_rpy(R_shank_foot(:,:,n)); %Sara's function to find roll, pitch, yaw angles; does not work for me for reasons I do not understand

        %these lines find the ankle center in the coord. system of the
        %iPECS (which changes during any dynamic trial), which will be used
        %to calculate moment arms. The sign is reversed because the coordinate
        %position of the iPECS center is actually what is of interest,
        %which is simply the same coordinate with signs flipped.
        ankle_center_lab_t(:,n)=[ankle_center_lab(n,1); ankle_center_lab(n,2); ankle_center_lab(n,3); 1];
        ankle_center_ipecs_inter(:,n)=T_iPECS(:,:,n)\ankle_center_lab_t(:,n);
        ankle_center_ipecs_inter(:,n)=-ankle_center_ipecs_inter(:,n);
                
    end
    
%Gait Cycle Analysis%

    %The first piece of this section pulls raw iPECS force data, decimates
    %it to match frequency of marker data, and filters it to find the
    %times of heel strike/toe off gait events
    cycle_info_raw(1,:)=Values{filei,2}(:,3);
    cycle_info_d=cycle_info_raw(1:10:length(cycle_info_raw));
    [d,c] = butter(4,(50*0.89568)/(1200/2),'low');
    cycle_info=filtfilt(d,c,cycle_info_d);
    
    %The following line removes the corresponding indices of force/moment
    %data so the lengths of both said data and the marker data are synched
    if exist('left_cut')
        cycle_info(left_cut:right_cut)=[];
    end

    %defines force threshold corresponding to existance of gait event and
    %marks the time at which the event occurs
    for n=2:length(cycle_info)   
        thres=40;
        if cycle_info(n)>thres && cycle_info(n-1)<thres %viewing indices directly before and after gait event
            event(n)=n;
        elseif cycle_info(n)<thres && cycle_info(n-1)>thres
            event(n)=n;
        end
    end

    event=event(event~=0); %remove zeros
    
    %This loop is run differently depending on whether data recording is
    %begun while subject is in stance or swing phase, in order to insure
    %that the gait cycle taken from the data is stance phase
    for n=1:length(event)-1 %creates vector of the length of each stance phase recorded by amputee
        if cycle_info(1)<thres
            if mod(n,2)~=0
                period(n)=event(n+1)-event(n);
            end
        else
            if mod(n,2)==0
                period(n)=event(n+1)-event(n);
            end
        end
    end
    
    %this "period" vector stores the length of each stance phase recorded,
    %and will have exactly the same amount of elements as stance phases
    %seen during recording if the code is working correctly. The vector
    %can be viewed for debugging purposes
    period=period(period~=0); %remove zeros

    for n=1:length(cycle_info) %my method of finding ankle angle over time, may not be entirely correct, rpy_angles would likely be a better substitute
        anklecen_kneecen_vec(n,:)=normalize(anklecen_kneecen_vec(n,:));

        cos_angle(n)=dot(foot_principle_y(n,:),anklecen_kneecen_vec(n,:)); %formula to find angle between two vectors
        angle(n)=acosd(cos_angle(n));

        angle_norm(n)=angle(n); %changed from normalized angle to absolute angle, but the vector name has remained
    end

%The following lines create a splined plot of ankle angle as a function of % gait cycle where each vector has 101 points, can be turned on or off for inspection:    
    angle_all=[];
    for n=1:length(period)
        if cycle_info(1)<thres
             high_bound=event(2*n)-1;
             low_bound=event(2*n-1);
        else
             high_bound=event(2*n+1)-1;
             low_bound=event(2*n);
        end
        angle_norm_cycle=angle_norm(low_bound:high_bound); 
        xx=0:0.01:1;
        percent_cycle=linspace(0,1,length(angle_norm_cycle));
        angle_splined=spline(percent_cycle,angle_norm_cycle(1,:),xx);
        angle_all=[angle_all; angle_splined];
        %plot(0:100,90-(angle_splined-90),'k','linewidth',2); %comment here to
        %pause(1)                                         %turn plot off
        %drawnow
        %hold on
        clear angle_norm_cycle percent_cycle xx angle_splined high_bound low_bound
    end
    
    %Force/Moment Translation/Transformation%
    
    %The following code takes the moment arms in each direction (X,Y,Z) from ankle 
    %center to iPECS center in the iPECS coord. system and multiplies them
    %by the forces in each direction recorded by the iPECS according to formulas
    %in literature to find the cooresponding forces and moments STILL IN THE
    %iPECS COORDNIATE SYSTEM, but at the ankle center. Any variable with two 
    %modifiers in the following code follows this notation: 1st is location, 2nd is
    %current coordinate system

    %finding moment arms and converting to proper units (mm to m)
    dx=(ankle_center_ipecs_inter(1,:)*(0.001));
    dy=(ankle_center_ipecs_inter(2,:)*(0.001));
    dz=(ankle_center_ipecs_inter(3,:)*(0.001));

    for n=1:3   %separating forces/moments into separate arrays and decimating
        f(n,:)=Values{filei,2}(:,n);
        f_ipecs_ipecs(n,:)=f(n,1:10:length(f(n,:)));
        m(n,:)=Values{filei,3}(:,n);
        m_ipecs_ipecs(n,:)=m(n,1:10:length(m(n,:)));
    end
    
    %Removing indices cooresponding to removed marker indices
    if exist('left_cut')==1
        f_ipecs_ipecs(:,left_cut:right_cut)=[];
        m_ipecs_ipecs(:,left_cut:right_cut)=[];
    end

    %finding forces at ankle
    f_ank_ipecs(1,:)=-f_ipecs_ipecs(1,:); 
    f_ank_ipecs(2,:)=-f_ipecs_ipecs(2,:);
    f_ank_ipecs(3,:)=-f_ipecs_ipecs(3,:);
    
    %finding moments at ankle
    m_ank_ipecs(1,:)=-m_ipecs_ipecs(1,:)-f_ipecs_ipecs(3,:).*dy+f_ipecs_ipecs(2,:).*dz; 
    m_ank_ipecs(2,:)=-m_ipecs_ipecs(2,:)+f_ipecs_ipecs(3,:).*dx-f_ipecs_ipecs(1,:).*dz;
    m_ank_ipecs(3,:)=-m_ipecs_ipecs(3,:)-f_ipecs_ipecs(2,:).*dx+f_ipecs_ipecs(1,:).*dy; 

    %this loop then transforms forces/moments at ankle center from iPECS
    %coord. system to final shank coord. system
    for n=1:length(f_ipecs_ipecs)
        f_ank_shank(n,:)=(inv(R_iPECS_shank(:,:,n))*f_ank_ipecs(:,n))'; 
        m_ank_shank(n,:)=(inv(R_iPECS_shank(:,:,n))*m_ank_ipecs(:,n))';
    end
    
    Values{filei,4}=f_ank_shank; %storing forces and moments in Values cell, normalizing moments by weight of subject
    Values{filei,5}=m_ank_shank./(str2double(weight));
    
    %creates plot of ankle angle vs. moment at ankle w.r.t lateral/medial
    %axis (Y axis in this case) that will update upon each gait cycle
    m_all=[];
    angle_all=[];
    for n=1:length(period) 
       if cycle_info(1)<thres
           high_bound=event(2*n)-1;
           low_bound=event(2*n-1);
       else
           high_bound=event(2*n+1)-1;
           low_bound=event(2*n);
       end
       angle_norm_cycle=angle_norm(low_bound:high_bound); 
       m_cycle=Values{filei,5}(low_bound:high_bound,2);
       xx=0:0.01:1;
       percent_cycle=linspace(0,1,length(angle_norm_cycle));
       angle_splined=spline(percent_cycle,angle_norm_cycle(1,:),xx);
       m_splined=spline(percent_cycle,m_cycle(:,1),xx);
       angle_all=[angle_all; angle_splined];
       m_all=[m_all; m_splined];
       %plot(90-(angle_splined-90),m_splined,'k','linewidth',2); %turn here to "hold on" off to stop plotting
       %title('Ankle Angle vs. Ankle Moment')
       %xlabel('Ankle Angle (°)')
       %ylabel('Ankle Moment (ft-lb/lbs)')
       %pause(1)                               
       %drawnow
       %hold on
       clear angle_norm_cycle m_cycle percent_cycle xx angle_splined m_splined high_bound low_bound
    end
    
    %clear ALL unnessary variables in preparation for next condition (i.e
    %Dn5 --> Level)
    clearvars -except angle_all m_all scalef offset cal D data filei side weight Values ...
        ankle_center_inter knee_center_inter lat_ankle_inter lat_knee_inter med_ankle_inter med_knee_inter period
    
    %NEW PLOT GENERATION%
    close all
    size(angle_all);
    for n=1:ans(1)
        plot(90-(angle_all(n,:)-90),-m_all(n,:))
        hold on
    end
    saveas(gcf,Values{filei,1}(1:length(Values{filei,1})-4),'png')
    %END NEW PLOT GENGERATION%
    
    
    close
    angle_mean=mean(angle_all); %averaging angle and ankle moment data across gait cycles into one vector representative of all cycles
    m_mean=mean(m_all);
    plot(90-(angle_mean-90),-m_mean,'linewidth',2); %avg. ankle angle vs. avg. ankle moment plot that contains all conditions
    a=axis;
    hold on
    
    figure %separate plot of avg. ankle angle vs. avg. ankle moment for the condition being tested, along with a plot of both data sets against gait cycle
    subplot 121
    plot(90-(angle_mean-90),-m_mean,'linewidth',2);
    subplot 122
    plot(90-(angle_mean-90))
    hold on
    yyaxis right
    plot(-m_mean)
    g=axis;
    axis([0 100 g(3) g(4)])
    pause
    close
    
    %This newly-added section will find the equilibrium point of each trial
    %(i.e at what angle the ankle moment is zero) and average each trial to
    %find averages of each respective condition.
%     for n=1:length(period)
%         equil=intersections(90-(angle_all(n,:)-90),-m_all(n,:),[min(90-(angle_all(n,:)-90)) max(90-(angle_all(n,:)-90))],[0 0]);
%         equil_a(n)=equil(1);
%         clear equil
%     end
%     equil_m=zeros(1,length(equil_a));
%     
%     figure
%     scatter(equil_a,equil_m,20,[0.64 0.83 0.93],'filled')
%     hold on
%     a_avg=mean(equil_a);
%     m_avg=mean(equil_m);
%     scatter(a_avg,m_avg,50,[1 0 0],'filled')
%     axis([a(1) a(2) a(3) a(4)])
%     pause
%     close
    
    figure
    
end