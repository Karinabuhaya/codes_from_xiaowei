function [Patient_reorder, PatientNameDate, varargout] = ListFolderContents_SPK_LFP_RAW_Coherent(varargin) 
% D:\Matlab\Test_data\patients

dbstop if error;

PatientNamePath = 'D:\Matlab\Test_data\patients';

file_name = dir(PatientNamePath);
num_file = length(file_name);

 for jj = 3:num_file
%for jj = 3
    FileNamePath = fullfile(PatientNamePath, file_name(jj).name);
    FolderContent = dir(FileNamePath);
    if length(FolderContent)<5 
        if strcmp(FolderContent(end).name(end), 'R')
            cd(FileNamePath);
            RightSide = dir('*_R');
            RightSidePath = fullfile(FileNamePath, RightSide.name);
            PatientMatPath(jj-2).RightSidePath = RightSidePath;
            PatientMatPath(jj-2).LeftSidePath = NaN;
            PatientNameDate(jj-2).RightSide = [file_name(jj).name,'_',RightSide.name];
            PatientNameDate(jj-2).LeftSide = NaN;
        elseif strcmp(FolderContent(end).name(end), 'L')
            cd(FileNamePath);
            LeftSide = dir('*_L');
            LeftSidePath = fullfile(FileNamePath, LeftSide.name);
            PatientMatPath(jj-2).RightSidePath = NaN;
            PatientMatPath(jj-2).LeftSidePath = LeftSidePath;
            PatientNameDate(jj-2).RightSide = NaN;
            PatientNameDate(jj-2).LeftSide = [file_name(jj).name,'_',LeftSide.name];
        end
    else
        cd(FileNamePath);
        RightSide = dir('*_R');
        LeftSide = dir('*_L');
        RightSidePath = fullfile(FileNamePath, RightSide.name);
        LeftSidePath = fullfile(FileNamePath, LeftSide.name);
        PatientMatPath(jj-2).RightSidePath = RightSidePath;
        PatientMatPath(jj-2).LeftSidePath = LeftSidePath;
        PatientNameDate(jj-2).RightSide = [file_name(jj).name,'_',RightSide.name];
        PatientNameDate(jj-2).LeftSide = [file_name(jj).name,'_',LeftSide.name];
    end
end
if nargin>0
    Patients_Indx = varargin{1};
else
    Patients_Indx = 1:length(PatientMatPath);
end

%%
for jj = Patients_Indx % patients index
    if ~isnan(PatientMatPath(jj).RightSidePath)
        % right side analysis of every patient
        cd(PatientMatPath(jj).RightSidePath);
        RightSideDepth = dir('*.mat');
        % find the multiple recording
        RightSideDepth = FindAndRemove_Multiple_Recording(RightSideDepth);
        % now, RightSideDepth.name include []
        
        if nargin>0
            if strcmp(varargin{2}(1),'R')
                Depth_Indx = 1:length(RightSideDepth);
            else
                Depth_Indx = [];
            end
        else
            Depth_Indx = 1:length(RightSideDepth);
        end
        % extract the depth number which is included in the RightSideDepth.name
        n = 0;
        for ii = Depth_Indx % depth index
            
            RightSideDepthName = RightSideDepth(ii).name;
            
            if isempty(RightSideDepthName)
                continue;
            else
                n = n+1;
            end
            
            Fi = strfind(RightSideDepthName, 'F');
            RightSideDepthPath = fullfile(RightSideDepth(ii).folder, RightSideDepth(ii).name);
            
            if nargin>0
                if str2double(RightSideDepthName(5:(Fi(1)-1))) == varargin{3}
                    [LFP_SPK_PSD_RMS, LFP_SPK_Coherence, FH] = LoadPath_SPK_LFP_RAW_Coherent(RightSideDepthPath,varargin{2});
                else
                    LFP_SPK_PSD_RMS = struct('LFP_PSD_1',[],'LFP_fre_1',[],'LFP_RMS_1',[],'LFP_PSD_2',[],'LFP_fre_2',[],'LFP_RMS_2',[],...
                        'SPK_PSD_1',[],'SPK_fre_1',[],'SPK_RMS_1',[],'SPK_PSD_2',[],'SPK_fre_2',[],'SPK_RMS_2',[]);
                    LFP_SPK_Coherence = struct('LFP_SPK_Coherence_1',[],'LFP_SPK_Coherence_2',[],'Fre',[]);
                end
            else
                [LFP_SPK_PSD_RMS, LFP_SPK_Coherence] = LoadPath_SPK_LFP_RAW_Coherent(RightSideDepthPath);
            end
            %{
            Patients_Recording_Time(jj).RightSide(n).Depth = str2double(RightSideDepthName(5:(Fi(1)-1)));
            Patients_Recording_Time(jj).RightSide(n).t_1 = LFP_SPK_PSD_RMS;
            Patients_Recording_Time(jj).RightSide(n).t_2 = LFP_SPK_Coherence;
            %}
            %
            Patient(jj).RightSideDepth(n).LocationNum = str2double(RightSideDepthName(5:(Fi(1)-1)));
            
            Patient(jj).RightSideDepth(n).LFP_PSD_1 = LFP_SPK_PSD_RMS.LFP_PSD_1;
            Patient(jj).RightSideDepth(n).LFP_Fre_1 = LFP_SPK_PSD_RMS.LFP_fre_1;
            Patient(jj).RightSideDepth(n).LFP_RMS_1 = LFP_SPK_PSD_RMS.LFP_RMS_1;
            Patient(jj).RightSideDepth(n).LFP_PSD_2 = LFP_SPK_PSD_RMS.LFP_PSD_2;
            Patient(jj).RightSideDepth(n).LFP_Fre_2 = LFP_SPK_PSD_RMS.LFP_fre_2;
            Patient(jj).RightSideDepth(n).LFP_RMS_2 = LFP_SPK_PSD_RMS.LFP_RMS_2;
            
            Patient(jj).RightSideDepth(n).SPK_PSD_1 = LFP_SPK_PSD_RMS.SPK_PSD_1;
            Patient(jj).RightSideDepth(n).SPK_Fre_1 = LFP_SPK_PSD_RMS.SPK_fre_1;
            Patient(jj).RightSideDepth(n).SPK_RMS_1 = LFP_SPK_PSD_RMS.SPK_RMS_1;
            Patient(jj).RightSideDepth(n).SPK_PSD_2 = LFP_SPK_PSD_RMS.SPK_PSD_2;
            Patient(jj).RightSideDepth(n).SPK_Fre_2 = LFP_SPK_PSD_RMS.SPK_fre_2;
            Patient(jj).RightSideDepth(n).SPK_RMS_2 = LFP_SPK_PSD_RMS.SPK_RMS_2;
            
            Patient(jj).RightSideDepth(n).LFP_SPK_Coherence_1 = LFP_SPK_Coherence.LFP_SPK_Coherence_1;
            Patient(jj).RightSideDepth(n).LFP_SPK_Coherence_2 = LFP_SPK_Coherence.LFP_SPK_Coherence_2;
            Patient(jj).RightSideDepth(n).Fre = LFP_SPK_Coherence.Fre;
            %}
        end
    else
        Patient(jj).RightSideDepth = NaN;
        %Patients_Recording_Time(jj).RightSide = NaN;
    end
  
    
    
    %%%%%%%%%%%%% left side analysis of every patient %%%%%%%%%%%%%%%%
    
    if ~isnan(PatientMatPath(jj).LeftSidePath)
        cd(PatientMatPath(jj).LeftSidePath);
        LeftSideDepth = dir('*.mat');
        % find the multiple recording
        LeftSideDepth = FindAndRemove_Multiple_Recording(LeftSideDepth);
        % now, LeftSideDepth.name include []
        
        if nargin>0
            if strcmp(varargin{2}(1),'L')
                Depth_Indx = 1:length(LeftSideDepth);
            else
                Depth_Indx = [];
            end
        else
            Depth_Indx = 1:length(LeftSideDepth);
        end

        %extract the depth number which is included in the LeftSideDepth.name
        m = 0;
        for ii = Depth_Indx
            
            LeftSideDepthName = LeftSideDepth(ii).name;
            
            if isempty(LeftSideDepthName)
                continue;
            else
                m = m+1;
            end
            
            Fi = strfind(LeftSideDepthName, 'F');
            LeftSideDepthPath = fullfile(LeftSideDepth(ii).folder, LeftSideDepth(ii).name);

            if nargin>0
                if str2double(LeftSideDepthName(5:(Fi(1)-1))) == varargin{3}
                    [LFP_SPK_PSD_RMS, LFP_SPK_Coherence, FH] = LoadPath_SPK_LFP_RAW_Coherent(LeftSideDepthPath,varargin{2});
                else
                    LFP_SPK_PSD_RMS = struct('LFP_PSD_1',[],'LFP_fre_1',[],'LFP_RMS_1',[],'LFP_PSD_2',[],'LFP_fre_2',[],'LFP_RMS_2',[],...
                        'SPK_PSD_1',[],'SPK_fre_1',[],'SPK_RMS_1',[],'SPK_PSD_2',[],'SPK_fre_2',[],'SPK_RMS_2',[]);
                    LFP_SPK_Coherence = struct('LFP_SPK_Coherence_1',[],'LFP_SPK_Coherence_2',[],'Fre',[]);
                end
            else
                [LFP_SPK_PSD_RMS, LFP_SPK_Coherence] = LoadPath_SPK_LFP_RAW_Coherent(LeftSideDepthPath);
            end
            %}
  
            %{
            Patients_Recording_Time(jj).LeftSide(m).Depth = str2double(LeftSideDepthName(5:(Fi(1)-1)));
            Patients_Recording_Time(jj).LeftSide(m).t_1 = LFP_SPK_PSD_RMS;
            Patients_Recording_Time(jj).LeftSide(m).t_2 = LFP_SPK_Coherence;
            %}
            %
            Patient(jj).LeftSideDepth(m).LocationNum = str2double(LeftSideDepthName(5:(Fi(1)-1)));
            
            Patient(jj).LeftSideDepth(m).LFP_PSD_1 = LFP_SPK_PSD_RMS.LFP_PSD_1;
            Patient(jj).LeftSideDepth(m).LFP_Fre_1 = LFP_SPK_PSD_RMS.LFP_fre_1;
            Patient(jj).LeftSideDepth(m).LFP_RMS_1 = LFP_SPK_PSD_RMS.LFP_RMS_1;
            Patient(jj).LeftSideDepth(m).LFP_PSD_2 = LFP_SPK_PSD_RMS.LFP_PSD_2;
            Patient(jj).LeftSideDepth(m).LFP_Fre_2 = LFP_SPK_PSD_RMS.LFP_fre_2;
            Patient(jj).LeftSideDepth(m).LFP_RMS_2 = LFP_SPK_PSD_RMS.LFP_RMS_2;
            
            Patient(jj).LeftSideDepth(m).SPK_PSD_1 = LFP_SPK_PSD_RMS.SPK_PSD_1;
            Patient(jj).LeftSideDepth(m).SPK_Fre_1 = LFP_SPK_PSD_RMS.SPK_fre_1;
            Patient(jj).LeftSideDepth(m).SPK_RMS_1 = LFP_SPK_PSD_RMS.SPK_RMS_1;
            Patient(jj).LeftSideDepth(m).SPK_PSD_2 = LFP_SPK_PSD_RMS.SPK_PSD_2;
            Patient(jj).LeftSideDepth(m).SPK_Fre_2 = LFP_SPK_PSD_RMS.SPK_fre_2;
            Patient(jj).LeftSideDepth(m).SPK_RMS_2 = LFP_SPK_PSD_RMS.SPK_RMS_2;

            Patient(jj).LeftSideDepth(m).LFP_SPK_Coherence_1 = LFP_SPK_Coherence.LFP_SPK_Coherence_1;
            Patient(jj).LeftSideDepth(m).LFP_SPK_Coherence_2 = LFP_SPK_Coherence.LFP_SPK_Coherence_2;
            Patient(jj).LeftSideDepth(m).Fre = LFP_SPK_Coherence.Fre;
            %}
        end
    else
        Patient(jj).LeftSideDepth = NaN;
        %Patients_Recording_Time(jj).LeftSide = NaN;
    end
end

%save('D:\Matlab\Test_data\test code\The Data of Avariables in Workspace\Patients_Recording_Time', 'Patients_Recording_Time');

%%
% reorder the result of 'patient' based on the field 'LocationNum'
if nargout==2
    for jj = 1: length(Patient)
        if isstruct(Patient(jj).RightSideDepth)
            [~, I_R] = sortrows([Patient(jj).RightSideDepth(:).LocationNum].', 1, 'descend');
            Patient_reorder(jj).RightSideDepth = Patient(jj).RightSideDepth(I_R);
            clear I_R;
        else
            Patient_reorder(jj).RightSideDepth = NaN;
        end
        if isstruct(Patient(jj).LeftSideDepth)
            [~, I_L] = sortrows([Patient(jj).LeftSideDepth(:).LocationNum].', 1, 'descend');
            Patient_reorder(jj).LeftSideDepth = Patient(jj).LeftSideDepth(I_L);
            clear I_L;
        else
            Patient_reorder(jj).LeftSideDepth = NaN;
        end
    end
end

if nargout==3
    Patient_reorder = [];
    varargout{1} = FH;
end


end

%% subfunction: find and remove replicated recordings

function Depth_Unique = FindAndRemove_Multiple_Recording(Depth)
Depth_Unique = Depth;
Depth_Cell = struct2cell(Depth(:));
Name_Cell = Depth_Cell(1, :);
%Fields = {'name', 'folder', 'date', 'bytes', 'isdir', 'datenum'};
for ii = 1:length(Depth)-1
    Fi = strfind(Depth(ii).name, 'F');
    if isempty(Fi)
        error('there is no F in Depth(ii).name');
    end
    if strcmp(Depth(ii).name(5:(Fi(1)-1)), Depth(ii+1).name(5:(Fi(1)-1)))
        IF = contains(Name_Cell, Depth(ii).name(5:(Fi(1)-1)));
        Indx = find(IF);
        [~, Max_Indx] = max([Depth(Indx).bytes]);
        for jj = 1:length(Indx)
            if jj == Max_Indx
                Depth_Unique(Indx(Max_Indx)).name = Depth(Indx(Max_Indx)).name;
            else
                Depth_Unique(Indx(jj)).name = [];
            end
        end
    end
end
end

