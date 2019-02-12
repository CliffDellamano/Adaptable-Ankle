%This script works with both Comparison and Prototype data, but only MATLAB
%data from ONE of the two protheses may be in the same path as the script at a
%time.

close
clc
clearvars -except scalef offset %clears all variables from previous iteration, except for calibration that carries over
fname=input('Which files to process? Press "c" for Comparison, "p" for Prototype ','s'); %these inputs tell the script which sections to run
side=input('Is the prothesis on the right or left? Press "r" for right, "l" for left ','s'); %some axes definitions will need to be flipped depending on the input here
weight=input('What is the subjects weight in kg? ','s'); %used to normalize moment data by weight of subject

%% Instantiating Forces and Moments
D=dir('*.mat'); %create directory for each condition's matlab file
bad=[];

for filei=1:length(D) %run for loop for each condition
    
    filename=D(filei).name; %load condition, store data, clear data matrix after pulling data from it
    eval(['load ' filename]);
    eval(['data = ' filename(1:(length(filename)-4)) ';']);
    eval(['clear ' filename(1:(length(filename)-4)) ';']);
    
    if D(filei).name(1) == 'S' || D(filei).name(1) == 'i' || D(filei).name(1) == 'A'
        bad=[bad; filei];
    elseif D(filei).name(1) ~= 'S' && D(filei).name(1) ~= 'i' && D(filei).name(1) ~= 'A'
        if D(filei).name(1) == 'D' % Walking downhill is opposite direction
            [Lcopx Lcopy Lcopz FL ML_top] = get_bertec_cop_w_markers(data,1,500);
            [Rcopx Rcopy Rcopz FR MR_top] = get_bertec_cop_w_markers(data,2,500);
        else
            [Lcopx Lcopy Lcopz FL ML_top] = get_bertec_cop_w_markers(data,2,500);
            [Rcopx Rcopy Rcopz FR MR_top] = get_bertec_cop_w_markers(data,1,500);
        end

        % Change to forces/moments on the foot instead of on the plate (Newton's 3rd
        % law - equal and opposite reaction forces) 

        Values{filei,1}=filename; %creating cell to hold forces/moments for each condition
        Values{filei,2}=FL; %storing force data, used to apply conversion from lb to N in iPECS code but unnessesary for Bertec
        Values{filei,3}=ML_top; %storing moment data
        Values{filei,4}=Lcopx;
        Values{filei,5}=Lcopy;
        Values{filei,6}=Lcopz;
        Values{filei,7}=FR;
        Values{filei,8}=MR_top;
        Values{filei,9}=Rcopx;
        Values{filei,10}=Rcopy;
        Values{filei,11}=Rcopz;

        clear FL ML_top Lcopx Lcopy Lcopz FR MR_top Rcopx Rcopy Rcopz %clearing all unnessary variables that will be written over when analyzing the next channel
    end
        
end

%These next few lines remove files that will not be analyzed, the names of
%which differ depending on if the prothesis is Comparison or Prototype.
%AnkleStatic data is searched for in both cases, but if it is not present
%in the set of files being analzyed, it will not affect the code.

D(bad)=[];
Values(cellfun('isempty',Values(:,1)),:)=[];

%% Determination of Ankle/Knee Markers for Prototype Foot
%This next section of code will only run if the set of files being analyzed
%belongs to a Prototype prothesis. The process for determining ankle/knee marker
%positions is different for this prothesis than the Comparison.

if strcmp(fname,'p')==1
    load AnkleStatic.mat %load static trial
    
    %removing unnessary fourth value from each vector
    AnkleStatic.Trajectories.Labeled.Data(:,4,:)=[];
    
    al=find(strcmp('Ank_Lat',AnkleStatic.Trajectories.Labeled.Labels)); %finding indices of important markers in dataset
    am=find(strcmp('Ank_Med',AnkleStatic.Trajectories.Labeled.Labels));
    ala=find(strcmp('Ank_Lat_Ant',AnkleStatic.Trajectories.Labeled.Labels));
    alp=find(strcmp('Ank_Lat_Post',AnkleStatic.Trajectories.Labeled.Labels));
    ama=find(strcmp('Ank_Med_Ant',AnkleStatic.Trajectories.Labeled.Labels));
    amp=find(strcmp('Ank_Med_Post',AnkleStatic.Trajectories.Labeled.Labels));
    
    for n=1:length(AnkleStatic.Trajectories.Labeled.Data) %axes of an ankle coordinate system are determined here, with the X axis running from posterior to anterior
        ankle_x(n,:)=normalize(((AnkleStatic.Trajectories.Labeled.Data(ala,:,n)+AnkleStatic.Trajectories.Labeled.Data(ama,:,n))./2)-...
            ((AnkleStatic.Trajectories.Labeled.Data(alp,:,n)+AnkleStatic.Trajectories.Labeled.Data(amp,:,n))./2));
        if strcmp(side,'r')==1 %temporary Y axis of ankle coordinate system, running medial to lateral on the right side/lateral to medial on left
            ankle_inter_y(n,:)=normalize(((AnkleStatic.Trajectories.Labeled.Data(ala,:,n)+AnkleStatic.Trajectories.Labeled.Data(alp,:,n))./2)-...
                ((AnkleStatic.Trajectories.Labeled.Data(ama,:,n)+AnkleStatic.Trajectories.Labeled.Data(amp,:,n))./2));
        else
            ankle_inter_y(n,:)=normalize(((AnkleStatic.Trajectories.Labeled.Data(ama,:,n)+AnkleStatic.Trajectories.Labeled.Data(amp,:,n))./2)-...
                ((AnkleStatic.Trajectories.Labeled.Data(ala,:,n)+AnkleStatic.Trajectories.Labeled.Data(alp,:,n))./2));
        end
        ankle_z(n,:)=cross(ankle_inter_y(n,:),ankle_x(n,:)); %defining Z axis, running vertical, by the cross product of the previous two axes
        ankle_y(n,:)=cross(ankle_z(n,:),ankle_x(n,:)); %new Y axis, running in the opposite direction as the intermediate, redefined to match conventions for axes (cross of x and y equals z)
        
        ankle_med(n,:)=AnkleStatic.Trajectories.Labeled.Data(am,:,n); %finding all important markers from ankle static trial,
        ankle_lat(n,:)=AnkleStatic.Trajectories.Labeled.Data(al,:,n); %in lab coordinate system
        ankle_lat_post(n,:)=AnkleStatic.Trajectories.Labeled.Data(alp,:,n);
    
        ankle_x(n,:)=normalize(ankle_x(n,:)); %normalizing vectors for use in matrices
        ankle_y(n,:)=normalize(ankle_y(n,:));
        ankle_z(n,:)=normalize(ankle_z(n,:));

        ankle_med_t(:,n)=[ankle_med(n,1);ankle_med(n,2);ankle_med(n,3); 1]; %augmenting important marker positions and transposing them for use in the following transformation matrix
        ankle_lat_t(:,n)=[ankle_lat(n,1);ankle_lat(n,2);ankle_lat(n,3); 1];
        ankle_cent_t(:,n)=(ankle_lat_t(:,n)+ankle_med_t(:,n))./2;
        ankle_lat_post_t(:,n)=[ankle_lat_post(n,1);ankle_lat_post(n,2);ankle_lat_post(n,3); 1];

        T_inter_ankle(:,:,n)=[ankle_x(n,1) ankle_y(n,1) ankle_z(n,1) ankle_lat_post_t(1,n); %forming ankle coord. system transformation matrix, with the lateral posterior ankle marker as the center
                    ankle_x(n,2) ankle_y(n,2) ankle_z(n,2) ankle_lat_post_t(2,n); 
                    ankle_x(n,3) ankle_y(n,3) ankle_z(n,3) ankle_lat_post_t(3,n);
                    0 0 0 1];
            
        med_ankle_inter(:,n)=T_inter_ankle(:,:,n)\ankle_med_t(:,n); %finding important marker positions in ankle coord. system
        lat_ankle_inter(:,n)=T_inter_ankle(:,:,n)\ankle_lat_t(:,n);
        ankle_center_inter(:,n)=T_inter_ankle(:,:,n)\ankle_cent_t(:,n);
    end
    
    med_ankle_inter=nanmean(med_ankle_inter,2); %averaging the ankle coord. system positions of each marker, as they should be completely static regardless
    lat_ankle_inter=nanmean(lat_ankle_inter,2);
    ankle_center_inter=nanmean(ankle_center_inter,2);
    
    load Static.mat
    %removing unnessary fourth value from each vector
    Static.Trajectories.Labeled.Data(:,4,:)=[];

    im=find(strcmp('iPECS_Med',Static.Trajectories.Labeled.Labels)); %finding indices of important markers in dataset
    il=find(strcmp('iPECS_Lat',Static.Trajectories.Labeled.Labels));
    ip=find(strcmp('iPECS_Post',Static.Trajectories.Labeled.Labels));
    if strcmp(side,'r')==1
        mk=find(strcmp('R_Knee_Med',Static.Trajectories.Labeled.Labels)); %the choice of marker here is determined by side of prothesis because the axis needs to point in the same direction regardless
        lk=find(strcmp('R_Knee_Lat',Static.Trajectories.Labeled.Labels));
    else
        mk=find(strcmp('L_Knee_Med',Static.Trajectories.Labeled.Labels));
        lk=find(strcmp('L_Knee_Lat',Static.Trajectories.Labeled.Labels));
    end  
    
    for n=1:length(Static.Trajectories.Labeled.Data)
        
        iPECS_med(n,:)=Static.Trajectories.Labeled.Data(im,:,n); %defining important marker positions and placing them in their own arrays
        iPECS_lat(n,:)=Static.Trajectories.Labeled.Data(il,:,n);
        iPECS_post(n,:)=Static.Trajectories.Labeled.Data(ip,:,n);
        
        iPECS_mid(n,:)=(iPECS_med(n,:)+iPECS_lat(n,:))/2;
        iPECS_mid_t(:,n)=[iPECS_mid(n,1);iPECS_mid(n,2);iPECS_mid(n,3); 1];
    
        med_knee(n,:)=Static.Trajectories.Labeled.Data(mk,:,n); %defining the medial knee marker position and placing it in its own array
        med_knee_t(:,n)=[med_knee(n,1);med_knee(n,2);med_knee(n,3); 1]; %augmenting and transposing in same manner as previous coord. system
        
        lat_knee(n,:)=Static.Trajectories.Labeled.Data(lk,:,n);
        lat_knee_t(:,n)=[lat_knee(n,1);lat_knee(n,2);lat_knee(n,3); 1];
        
        if strcmp(side,'r')==1 %an iPECS coordinate system is created here, with the temporary X axis running med-lat or lat-med
            iPECS_inter_x(n,:)=normalize(iPECS_lat(n,:)-iPECS_med(n,:)); %the order of markers here is decided in the same way as above for the ankle coord. system
        else
            iPECS_inter_x(n,:)=normalize(iPECS_med(n,:)-iPECS_lat(n,:));
        end
        iPECS_y(n,:)=normalize(iPECS_mid(n,:)-iPECS_post(n,:)); %Y axis defined as running posterior-anterior, using the midpoint of the med-lat markers as the anterior 
        iPECS_z(n,:)=cross(iPECS_inter_x(n,:),iPECS_y(n,:)); %Z axis created by cross product of temporary X and Y axes
        iPECS_x(n,:)=cross(iPECS_y(n,:),iPECS_z(n,:)); %new X axis created by cross product of Y and Z axes, matching coordinate system conventions
    
        iPECS_x(n,:)=normalize(iPECS_x(n,:)); %normalizing vectors once more
        iPECS_y(n,:)=normalize(iPECS_y(n,:));
        iPECS_z(n,:)=normalize(iPECS_z(n,:));

        T_iPECS(:,:,n)=[iPECS_x(n,1) iPECS_y(n,1) iPECS_z(n,1) iPECS_mid_t(1,n); %forming iPECS coord. system transformation matrix, with the midpoint of the med-lat markers as the origin
                        iPECS_x(n,2) iPECS_y(n,2) iPECS_z(n,2) iPECS_mid_t(2,n); 
                        iPECS_x(n,3) iPECS_y(n,3) iPECS_z(n,3) iPECS_mid_t(3,n);
                        0 0 0 1];

        med_knee_inter(:,n)=T_iPECS(:,:,n)\med_knee_t(:,n); %finding medial knee marker position in ankle coord. system
        lat_knee_inter(:,n)=T_iPECS(:,:,n)\lat_knee_t(:,n);
    end
    
    med_knee_inter=nanmean(med_knee_inter,2); %averaging ankle coord. system position(s) of important marker(s)
    lat_knee_inter=nanmean(lat_knee_inter,2);
    
end 
%% Determination of Ankle/Knee Markers for Comparison Foot
%similar to above, this section will only run for the Comparison foot as
%the method of determining ankle/knee markers is different

if strcmp(fname,'c')==1
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
    
    med_ankle_inter=nanmean(med_ankle_inter,2); %averaging inter. shank coord. system positions of important markers, as they should not be moving regardless
    lat_ankle_inter=nanmean(lat_ankle_inter,2);
    med_knee_inter=nanmean(med_knee_inter,2);
    lat_knee_inter=nanmean(lat_knee_inter,2);
    
    ankle_center_inter=nanmean(ankle_center_inter,2);
    knee_center_inter=nanmean(knee_center_inter,2);

end
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
    clearvars -except scalef offset D Values weight cal fname side filei data ankle_center_inter knee_center_inter ... 
        lat_ankle_inter lat_knee_inter med_ankle_inter med_knee_inter
        
    if strcmp(side,'r')==1
        al=find(strcmp('R_Ank_Lat',data.Trajectories.Labeled.Labels)); %finding indices of important markers, the side of which will depend on the side of the prothesis 
        kl=find(strcmp('R_Knee_Lat',data.Trajectories.Labeled.Labels));
        s=find(strcmp('R_Shank',data.Trajectories.Labeled.Labels));
        h=find(strcmp('R_Heel',data.Trajectories.Labeled.Labels));
        t=find(strcmp('R_Toe',data.Trajectories.Labeled.Labels));
        ala=find(strcmp('R_Ank_Lat_Ant',data.Trajectories.Labeled.Labels));
        alp=find(strcmp('R_Ank_Lat_Post',data.Trajectories.Labeled.Labels));
        ama=find(strcmp('R_Ank_Med_Ant',data.Trajectories.Labeled.Labels));
        amp=find(strcmp('R_Ank_Med_Post',data.Trajectories.Labeled.Labels));
        im=find(strcmp('iPECS_Med',data.Trajectories.Labeled.Labels));
        il=find(strcmp('iPECS_Lat',data.Trajectories.Labeled.Labels));
        ip=find(strcmp('iPECS_Post',data.Trajectories.Labeled.Labels));
    elseif strcmp(side,'l')==1
        al=find(strcmp('L_Ank_Lat',data.Trajectories.Labeled.Labels)); 
        kl=find(strcmp('L_Knee_Lat',data.Trajectories.Labeled.Labels));
        s=find(strcmp('L_Shank',data.Trajectories.Labeled.Labels));
        h=find(strcmp('L_Heel',data.Trajectories.Labeled.Labels));
        t=find(strcmp('L_Toe',data.Trajectories.Labeled.Labels));
        ala=find(strcmp('L_Ank_Lat_Ant',data.Trajectories.Labeled.Labels));
        alp=find(strcmp('L_Ank_Lat_Post',data.Trajectories.Labeled.Labels));
        ama=find(strcmp('L_Ank_Med_Ant',data.Trajectories.Labeled.Labels));
        amp=find(strcmp('L_Ank_Med_Post',data.Trajectories.Labeled.Labels));
        im=find(strcmp('iPECS_Med',data.Trajectories.Labeled.Labels));
        il=find(strcmp('iPECS_Lat',data.Trajectories.Labeled.Labels));
        ip=find(strcmp('iPECS_Post',data.Trajectories.Labeled.Labels));
    end
    
    %as denoted above, this section contains the coord. systems used to
    %define ankle/knee positions for Comparison. The same coord. systems
    %will be used to find these marker positions in the lab coord. system
    %by the reverse process
    if strcmp(fname,'c')==1
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
            
        end

        med_ankle_lab(:,4)=[]; %removing unnessesary fourth value that is added when marker positions are found in lab. coord system
        lat_ankle_lab(:,4)=[];
        med_knee_lab(:,4)=[];
        lat_knee_lab(:,4)=[];
        ankle_center_lab(:,4)=[];
        knee_center_lab(:,4)=[];
        
    end
    
    %as done above, the coord. systems used in the static trials for the
    %Prototype foot are repeated here to find each important marker
    %position in the lab frame during dynamic trials.
    if strcmp(fname,'p')==1
        for n=1:length(data.Trajectories.Labeled.Data)
            ankle_x(n,:)=normalize(((data.Trajectories.Labeled.Data(ala,:,n)+data.Trajectories.Labeled.Data(ama,:,n))./2)-...
            ((data.Trajectories.Labeled.Data(alp,:,n)+data.Trajectories.Labeled.Data(amp,:,n))./2));
            if strcmp(side,'r')==1
                ankle_inter_y(n,:)=normalize(((data.Trajectories.Labeled.Data(ala,:,n)+data.Trajectories.Labeled.Data(alp,:,n))./2)-...
                    ((data.Trajectories.Labeled.Data(ama,:,n)+data.Trajectories.Labeled.Data(amp,:,n))./2));
            else
                ankle_inter_y(n,:)=normalize(((data.Trajectories.Labeled.Data(ama,:,n)+data.Trajectories.Labeled.Data(amp,:,n))./2)-...
                    ((data.Trajectories.Labeled.Data(ala,:,n)+data.Trajectories.Labeled.Data(alp,:,n))./2));
            end
            ankle_z(n,:)=cross(ankle_inter_y(n,:),ankle_x(n,:));
            ankle_y(n,:)=cross(ankle_z(n,:),ankle_x(n,:));
         
            ankle_lat_post(n,:)=data.Trajectories.Labeled.Data(alp,:,n);
            
            ankle_x(n,:)=normalize(ankle_x(n,:)); 
            ankle_y(n,:)=normalize(ankle_y(n,:));
            ankle_z(n,:)=normalize(ankle_z(n,:));
            
            T_inter_ankle(:,:,n)=[ankle_x(n,1) ankle_y(n,1) ankle_z(n,1) ankle_lat_post(n,1); 
                ankle_x(n,2) ankle_y(n,2) ankle_z(n,2) ankle_lat_post(n,2); 
                ankle_x(n,3) ankle_y(n,3) ankle_z(n,3) ankle_lat_post(n,3);
                0 0 0 1];
            
            med_ankle_lab(n,:)=(T_inter_ankle(:,:,n)*med_ankle_inter)'; %finding positions of each marker in lab frame by multiplying known intermediate position
            lat_ankle_lab(n,:)=(T_inter_ankle(:,:,n)*lat_ankle_inter)';
            ankle_center_lab(n,:)=(lat_ankle_lab(n,:)+med_ankle_lab(n,:))/2;
            
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
            
            T_iPECS(:,:,n)=[iPECS_x(n,1) iPECS_y(n,1) iPECS_z(n,1) iPECS_mid(n,1); 
                        iPECS_x(n,2) iPECS_y(n,2) iPECS_z(n,2) iPECS_mid(n,2);
                        iPECS_x(n,3) iPECS_y(n,3) iPECS_z(n,3) iPECS_mid(n,3);
                        0 0 0 1];            
            
            med_knee_lab(n,:)=(T_iPECS(:,:,n)*med_knee_inter)';
            lat_knee_lab(n,:)=(T_iPECS(:,:,n)*lat_knee_inter)';
            knee_center_lab(n,:)=(lat_knee_lab(n,:)+med_knee_lab(n,:))/2;

        end
        
        med_ankle_lab(:,4)=[]; 
        lat_ankle_lab(:,4)=[];
        med_knee_lab(:,4)=[];
        ankle_center_lab(:,4)=[];
        knee_center_lab(:,4)=[];
        
    end
    
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
                
    end
    
%Gait Cycle Analysis%

    %The first piece of this section pulls raw Bertec force data, decimates
    %it to match frequency of marker data, and filters it to find the
    %times of heel strike/toe off gait events
    if strcmp(side,'r')==1
        cycle_info_d(1,:)=-Values{filei,7}(3,:);
        cycle_info_opposite_d(1,:)=-Values{filei,2}(3,:);
    elseif strcmp(side,'l')==1
        cycle_info_d(1,:)=-Values{filei,2}(3,:);
        cycle_info_opposite_d(1,:)=-Values{filei,7}(3,:);
    end

    [d,c] = butter(4,(50*0.89568)/(1200/2),'low');
    cycle_info=filtfilt(d,c,cycle_info_d);
    cycle_info_opposite=filtfilt(d,c,cycle_info_opposite_d);

    %defines force threshold corresponding to existance of gait event and
    %marks the time at which the event occurs
    figure;
    plot(cycle_info);
    hold on
    plot(cycle_info_opposite);
    title('Select force threshold for gait events.');
    [~,y]=ginput(1);
    thres=round(y);
    close;
    
    for n=2:length(cycle_info)   
        if cycle_info(n)>thres && cycle_info(n-1)<thres %viewing indices directly before and after gait event
            HS(n)=n;
        elseif cycle_info(n)<thres && cycle_info(n-1)>thres
            TO(n)=n;
        elseif cycle_info_opposite(n)>thres && cycle_info_opposite(n-1)<thres
            OHS(n)=n;
        elseif cycle_info_opposite(n)<thres && cycle_info_opposite(n-1)>thres
            OTO(n)=n;
        end
    end

    HS=HS(HS~=0); %remove zeros
    TO=TO(TO~=0);
    OHS=OHS(OHS~=0);
    OTO=OTO(OTO~=0);
    
    if length(OHS)>length(OTO)
        OHS(1)=[];
    elseif length(OHS)<length(OTO)
        OTO(end)=[];
    elseif OHS(1)<OTO(1)
        OHS(1)=[];
        OTO(end)=[];
    end

    for n=1:length(cycle_info) %my method of finding ankle angle over time, may not be entirely correct, rpy_angles would likely be a better substitute
        anklecen_kneecen_vec(n,:)=normalize(anklecen_kneecen_vec(n,:));

        cos_angle(n)=dot(foot_principle_y(n,:),anklecen_kneecen_vec(n,:)); %formula to find angle between two vectors
        angle(n)=acosd(cos_angle(n));

        angle_norm(n)=angle(n); %changed from normalized angle to absolute angle, but the vector name has remained
    end

%The following lines create a splined plot of ankle angle as a function of % gait cycle where each vector has 101 points, can be turned on or off for inspection:    
    angle_all=[];
    for n=1:length(OHS)
        low_bound=OTO(n);
        high_bound=OHS(n);
        angle_norm_cycle=angle_norm(low_bound:high_bound); 
        xx=0:0.01:1;
        percent_cycle=linspace(0,1,length(angle_norm_cycle));
        angle_splined=spline(percent_cycle,angle_norm_cycle(1,:),xx);
        angle_all=[angle_all; angle_splined];
%         plot(0:100,90-(angle_splined-90),'k','linewidth',2); %comment here to
%         pause(1)                                         %turn plot off
%         drawnow
%         hold on
        clear angle_norm_cycle percent_cycle xx angle_splined high_bound low_bound
    end
    
    %Force/Moment Translation/Transformation%
    
    %The following code takes the moment arms in each direction (X,Y,Z) from the location 
    %of the center of pressure (COP) to ankle center in the shank coord. system and 
    %multiplies them by the forces in each direction recorded by the force plates according to 
    %formulas both in literature and the 'find_rocker_shapes'.m file to find the 
    %cooresponding forces and moments at the ankle center. Any variable with two 
    %modifiers in the following code follows this notation: 1st is location, 2nd is
    %current coordinate system
   
    
    % Change to forces/moments on the foot instead of on the plate (Newton's 3rd
    % law - equal and opposite reaction forces)
    FLx=-Values{filei,2}(1,:);
    FLy=-Values{filei,2}(2,:);
    FLz=-Values{filei,2}(3,:);
    FRx=-Values{filei,7}(1,:);
    FRy=-Values{filei,7}(2,:);
    FRz=-Values{filei,7}(3,:);
    
    MLx=-Values{filei,3}(1,:);
    MLy=-Values{filei,3}(2,:);
    MLz=-Values{filei,3}(3,:);
    MRx=-Values{filei,8}(1,:);
    MRy=-Values{filei,8}(2,:);
    MRz=-Values{filei,8}(3,:);
    
    Lcopx=Values{filei,4};
    Lcopy=Values{filei,5};
    Lcopz=Values{filei,6};
    Rcopx=Values{filei,9};
    Rcopy=Values{filei,10};
    Rcopz=Values{filei,11};
   
    if strcmp(side,'r')==1
        m_global(1,:)=MRx+FRz.*(0.001*(Rcopy-ankle_center_lab(:,2)'))+FRy.*(0.001*(ankle_center_lab(:,3)'-Rcopz));
        m_global(2,:)=MRy-FRz.*(0.001*(Rcopx-ankle_center_lab(:,1)'))-FRx.*(0.001*(ankle_center_lab(:,3)'-Rcopz));
        m_global(3,:)=MRz+FRy.*(0.001*(Rcopx-ankle_center_lab(:,1)'))-FRx.*(0.001*(Rcopy-ankle_center_lab(:,2)'));
    else
        m_global(1,:)=MLx+FLz.*(0.001*(Lcopy-ankle_center_lab(:,2)'))+FLy.*(0.001*(ankle_center_lab(:,3)'-Lcopz));
        m_global(2,:)=MLy-FLz.*(0.001*(Lcopx-ankle_center_lab(:,1)'))-FLx.*(0.001*(ankle_center_lab(:,3)'-Lcopz));
        m_global(3,:)=MLz+FLy.*(0.001*(Lcopx-ankle_center_lab(:,1)'))-FLx.*(0.001*(Lcopy-ankle_center_lab(:,2)'));
    end
    
    for n=1:length(R_shank)
        m_ankle_sagital(:,n)=inv(R_shank(:,:,n))*m_global(:,n);
    end
    
    Values{filei,12}=m_ankle_sagital./(str2double(weight));
    
    %creates plot of ankle angle vs. moment at ankle w.r.t lateral/medial
    %axis (Y axis in this case) that will update upon each gait cycle
    m_all=[];
    angle_all=[];
    for n=1:length(OHS) 
       low_bound=OTO(n);
       high_bound=OHS(n);
       angle_norm_cycle=angle_norm(low_bound:high_bound);
       m_cycle=Values{filei,12}(2,low_bound:high_bound);
       if sum(isnan(m_cycle))==0       
           xx=0:0.01:1;
           percent_cycle=linspace(0,1,length(angle_norm_cycle));
           angle_splined=spline(percent_cycle,angle_norm_cycle(1,:),xx);
           m_splined=spline(percent_cycle,m_cycle(1,:),xx);
           angle_all=[angle_all; angle_splined];
           m_all=[m_all; m_splined];
    %        plot(0:100,m_splined,'k','linewidth',2);       %comment here to
    %        pause(1)                                         %turn plots off
    %        drawnow
    %        hold on       
    %        plot(90-(angle_splined-90),m_splined,'k','linewidth',2);
    %        title('Ankle Angle vs. Ankle Moment')
    %        xlabel('Ankle Angle (°)')
    %        ylabel('Ankle Moment (ft-lb/lbs)')
    %        pause(1)                               
    %        drawnow
    %        hold on
           clear angle_norm_cycle m_cycle percent_cycle xx angle_splined m_splined high_bound low_bound
       end
    end
    
    %clear ALL unnessary variables in preparation for next condition (i.e
    %Dn5 --> Level)
    clearvars -except angle_all m_all scalef offset cal D data filei fname side weight Values ...
        ankle_center_inter knee_center_inter lat_ankle_inter lat_knee_inter med_ankle_inter med_knee_inter period
    
    close
    angle_mean=nanmean(angle_all); %averaging angle and ankle moment data across gait cycles into one vector representative of all cycles
    m_mean=nanmean(m_all);
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
%     a_avg=nanmean(equil_a);
%     m_avg=nanmean(equil_m);
%     scatter(a_avg,m_avg,50,[1 0 0],'filled')
%     axis([a(1) a(2) a(3) a(4)])
%     pause
%     close
    
    figure
    
end