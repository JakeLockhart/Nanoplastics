classdef MassSpectrometryDataAnalysis
    % <Documentation>
        % MassSpectrometryDataAnalysis()
        %   Load and process nanoplastics mass spectrometry data, perform quality
        %   control, generate summary statistics and plots. Data collected from the 
        %   Proteomics and Mass Spectrometry Core Facility - PSU Huck.
        %   Created by: jsl5865
        %
        % Syntax:
        %   Analyzer = MassSpectrometryDataAnalysis
        %
        % Description:
        %   This class (via constructor) initializes an object `Analyzer` that handles
        %   mass spectrometry datasets related to nanoplastic analysis. It performs:
        %     • Data loading from .csv or .xlsx files
        %     • Preprocessing
        %       - Adding normalized spectral abundancy factor to data
        %       - Updating original data tables
        %     • User selection of file(s) (Gels) and parameter(s) of interest
        %     • Concatentated data display
        %       - Generating new tables with only parameters of interest
        %       - Combining files to compare parameters across gels
        %     • Visualization
        %       - Generating protein fold enrichment plots 
        %   Note that all files must be contained within a parent folder. Constructor
        %   looks for all .csv files within a single folder selected upon startup.
        %
        % Input:
        %   Class properties are created during constructor, no input variables to class.
        %   Most instance methods require UI interface as input.
        %
        % Output:
        %   Class does not have output.
        %   Instance methods have outputs.
        %
        % Methods (public):
        %   TotalProteins()
        %   GelParameters(SelectionMode)
        %   CreateExcelFile(GelData, FileName)
        %   PlotFoldEnrchiment()
    % <End Documentation>

    properties
        Variables = [
                     "Checked"
                     "Master"
                     "Accession"
                     "Description"
                     "Coverage [%]"
                     "# Peptides"
                     "# PSMs"
                     "# Unique Peptides"
                     "# AAs"
                     "MW [kDa]"
                     "calc. pI"
                     "Score Sequest HT: Sequest HT"
                     "# Peptides (by Search Engine): Sequest HT"
                     "Contaminant"
                     "Biological Process"
                     "Cellular Component"
                     "Molecular Function"
                     "Pfam IDs"
                     "Entrez Gene ID"
                     "Gene Symbol"
                     "Gene ID"
                     "Ensembl Gene ID"
                     "WikiPathways"
                     "Reactome Pathways"
                     "# Protein Pathway Groups"
                     "Found in Sample: [S13] F13: Sample"
                     "# Protein Groups"
                    ];
        Variables_Numeric = [ 
                            "Coverage [%]"
                            "# Peptides"
                            "# PSMs"
                            "# Unique Peptides"
                            "# AAs"
                            "MW [kDa]"
                            "calc. pI"
                            "Score Sequest HT: Sequest HT"
                            "# Peptides (by Search Engine): Sequest HT"
                            "Entrez Gene ID"
                            "# Protein Pathway Groups"
                            "# Protein Groups"
                           ];
        Variables_Text = [
                         "Checked"
                         "Master"
                         "Accession"
                         "Description"
                         "Contaminant"
                         "Biological Process"
                         "Cellular Component"
                         "Molecular Function"
                         "Pfam IDs"
                         "Gene Symbol"
                         "Gene ID"
                         "Ensembl Gene ID"
                         "WikiPathways"
                         "Reactome Pathways"
                         "Found in Sample: [S13] F13: Sample"
                        ];
        Keys = [
                "Accession"
                "Description"
               ]
    end
    
    properties
        DirectoryInfo struct = struct()
        MatLabVariables struct = struct()
        Data struct = struct()
    end

    %% Constructor
    methods
        function obj = MassSpectrometryDataAnalysis()
            obj.DirectoryInfo = MassSpectrometryDataAnalysis.FindFiles("xlsx", "SingleFolder");
            obj.MatLabVariables = MassSpectrometryDataAnalysis.MakeValidVariables(obj);
            obj = ReadData(obj);
            obj = Add_SAF_NSAF(obj);
        end
    end

    %% Callable functions
    methods (Access = public)
        function TotalProteins(obj)
            Fields = string(fieldnames(obj.Data));
            Alignment = max(strlength(Fields));
            for i = 1:length(Fields)
                ProteinCount = height(obj.Data.(Fields(i)).Description);
                fprintf("%*s   -   %i\n",Alignment, Fields(i), ProteinCount)
            end
        end

        function [GelData, Display] = GelParameters(obj, SelectionMode)
            arguments
                obj
                SelectionMode string {mustBeMember(SelectionMode, ["single", "multiple"])} = "single";
            end
            Display = obj.UI_FilesandVariables(SelectionMode);
            GelData = [];
            switch SelectionMode
                case "single"
                    GelData = SingleGelProperties(obj, Display);
                case "multiple"
                    GelData = MultipleGelProperties(obj, Display);
            end
        end

        function CreateExcelFile(~, GelData, FileName)
            FileName = FileName + ".xlsx";
            writetable(GelData, FileName);
        end

        function [GelData, Proteins, MissingProteins, PointsOfInterest] = PlotFoldEnrchiment(obj)
            GelSet = UI_DefineInputUnboundElution(obj);
            Parameter = UI_GetVariables(obj, GelSet);
            Fields = fieldnames(GelSet);
            for i = 1:length(fieldnames(GelSet))     
                Display.Gel(i) = string(GelSet.(Fields{i}));
            end
            Display.Vars = Parameter;
            Display.Vars = union(obj.Keys, Display.Vars, "stable");

            GelData = MultipleGelProperties(obj, Display);
            Elutions = obj.ExtractValuesFromMergedTable(GelData, Fields);
            Proteins = obj.FoldEnrichment(GelData, Elutions, Fields);
            MissingProteins = obj.MissingFoldEnrichment(Proteins);
            PointsOfInterest = obj.ScatterProteins(Proteins);
            obj.ScatterMissingProteins(MissingProteins);
            obj.UniProtIdentification(PointsOfInterest);
        end

        function [ProteinInfo] = ProteinIdentification(obj, ProteinIDs)
            [ProteinInfo] = UniProtIdentification(obj, ProteinIDs);
        end
    end

    %% Initializtion functions
    methods (Static, Access = private)
        function DirectoryInfo = FindFiles(FileType, SearchMode, ConstantAddress)
            arguments
                FileType char {mustBeMember(FileType, {'csv', 'xlsx', 'txt', 'tiff', 'mdf'})} = 'csv'
                SearchMode char {mustBeMember(SearchMode, {'SingleFolder', 'AllSubFolders', 'TroubleShoot'})} = 'SingleFolder'
                ConstantAddress char = ''
            end

            DirectoryInfo = FileLookup(FileType, SearchMode, ConstantAddress);

            function Lookup = FileLookup(FileType, SearchMode, ConstantAddress)
                %% Define the file structure
                    Lookup.FileType = strcat('*.', FileType);   % Choose file type
                
                %% Select the folder
                    switch SearchMode                                                               % SearchMode = Constant filepath for testing
                        case 'TroubleShoot'                                                             % User input: Troubleshoot
                            if ~isfolder(ConstantAddress)                                                   % Validate that user defined path exists and is valid
                                error("For 'TroubleShoot' mode, you must provide a" + ...                   % Throw error if invalid
                                    " valid folder path as the third argument.");      
                            end
                                Lookup.FolderAddress = ConstantAddress;                                     % Set folder address to user defined path
                                searchPattern = fullfile(Lookup.FolderAddress, Lookup.FileType);            % Create a search pattern based on folder address
                        case 'SingleFile'                                                               % User input: SingleFile
                            [FileName, FolderPath] = uigetfile(Lookup.FileType, 'Select a file');           % Prompt user to select a single file
                            if isequal(FileName, 0)                                                         % Validate file selection 
                                error("No file selected");                                                  % Throw error if invalid
                            end
                            Lookup.FolderAddress = FolderPath;                                              % Set folder address to user defined path
                            searchPattern = fullfile(Lookup.FolderAddress, FileName);                       % Create a search pattern based on folder address
                        case 'SingleFolder'                                                             % User input: SingleFolder
                            Lookup.FolderAddress = uigetdir('*.*', 'Select a folder');                      % Prompt user to select a folder
                            if isequal(Lookup.FolderAddress, 0)                                             % Validate folder selection
                                error('No folder selected');                                                % Throw error if invalid
                            end
                            searchPattern = fullfile(Lookup.FolderAddress, Lookup.FileType);                % Create a search pattern based on an individual folder
                        case 'AllSubFolders'                                                            % User input: Troubleshoot
                            Lookup.FolderAddress = uigetdir('*.*', 'Select a folder');                      % Prompt user to select a folder
                            if isequal(Lookup.FolderAddress, 0)                                             % Validate folder selection
                                error('No folder selected');                                                % Throw error if invalid
                            end
                            searchPattern = fullfile(Lookup.FolderAddress, '**', Lookup.FileType);          % Create a search pattern based on all sub folders
                    end

                %% Find All FileType within defined folder
                    Lookup.AllFiles = searchPattern;                                                                    % Create general file path
                    Lookup.FolderInfo = dir(searchPattern);                                                             % Identify the folder directory
                    Lookup.FileCount = length(Lookup.FolderInfo);                                                       % Determine the number of files in this folder
                    Lookup.FolderCount = length(unique({Lookup.FolderInfo.folder}));                                    % Determine the number of folders 
                    [~, Lookup.CurrentFolder] = fileparts(Lookup.FolderAddress);                                        % Collect folder information
                    Lookup.Path = arrayfun(@(x) fullfile(x.folder, x.name), Lookup.FolderInfo, 'UniformOutput', false); % Identify the file path
            end

        end

        function MatLabVariables = MakeValidVariables(obj)
            VariableList = ["Variables" "Variables_Numeric" "Variables_Text"];
            
            for i = 1:length(VariableList)
                Temp_Variables = obj.(VariableList(i));
                MatLabVariables.(VariableList(i)) = matlab.lang.makeValidName(Temp_Variables);
            end
        end
    end

    %% Create default data set
    methods (Access = private)        
        function obj = ReadData(obj)
            Dir = obj.DirectoryInfo;
            Var = obj.MatLabVariables.Variables;
            NVar = obj.MatLabVariables.Variables_Numeric;

            for i = 1:Dir.FileCount
                GelBand = erase(Dir.FolderInfo(i).name, ".xlsx");
                Temp_Dir = fullfile(Dir.FolderInfo(i).folder, GelBand);
                Temp_File = readtable(Temp_Dir, "VariableNamingRule", "preserve");
                for j = 1:length(Var)
                    GelBand = matlab.lang.makeValidName(GelBand);
                    if ismember(Var(j), NVar)
                        obj.Data.(GelBand).(Var(j)) = table2array(Temp_File(:,j));
                    else
                        obj.Data.(GelBand).(Var(j)) = Temp_File(:,j);
                    end
                end
            end
        end

        function obj = Add_SAF_NSAF(obj)
            GelBands = fieldnames(obj.Data);
            for i = 1:length(fieldnames(obj.Data))
                obj.Data.(GelBands{i}).SAF = Calculate_SAF(obj, GelBands{i});
                obj.Data.(GelBands{i}).NSAF = Calculate_NSAF(obj, GelBands{i});
            end
            obj = UpdateVariables(obj, "SAF");
            obj = UpdateVariables(obj, "NSAF");
        end
    end

    %% Helper functions
    methods (Access = private)
        function obj = UpdateVariables(obj, NewVariable)
            VariableList = ["Variables" "Variables_Numeric" "Variables_Text"];
            NewVariable = matlab.lang.makeValidName(NewVariable);

            for i = 1:length(VariableList)
                obj.(VariableList(i)) = union(obj.(VariableList(i)), NewVariable);
                obj.MatLabVariables.(VariableList(i)) = union(obj.(VariableList(i)), NewVariable);
            end
        end

        function SAF = Calculate_SAF(obj, GelBand)
            SAF = obj.Data.(GelBand).x_PSMs ./ obj.Data.(GelBand).x_AAs;
        end

        function NSAF = Calculate_NSAF(obj, GelBand)
            SAF = obj.Data.(GelBand).SAF;
            NSAF = SAF ./ sum(SAF);
        end

        function Display = UI_FilesandVariables(obj, SelectionMode)
            arguments
                obj
                SelectionMode string {mustBeMember(SelectionMode, ["single", "multiple"])} = "single";
            end

            Display.Gel = listdlg("PromptString", "Choose a Gel Band", "SelectionMode", SelectionMode, "ListString", fieldnames(obj.Data));
            Temp_Gel = fieldnames(obj.Data);
            Display.Gel = string(Temp_Gel(Display.Gel));

            Display.Vars = listdlg("PromptString", "Select properties to view", "SelectionMode", "multiple", "ListString", fieldnames(obj.Data.(Display.Gel(1))));
            Temp_Vars = fieldnames(obj.Data.(Display.Gel(1)));
            Display.Vars = Temp_Vars(Display.Vars);
            Display.Vars = union(obj.Keys, Display.Vars, "stable");
            

            fprintf("Gel band(s) selected:\n")
            for i = 1:length(Display.Gel)
                fprintf("\t%s\n", Display.Gel(i))
            end
            fprintf("Property(s) selected: \n")
            for i = 1:length(Display.Vars)
                fprintf("\t%s\n", Display.Vars(i))
            end
        end

        function GelTypes = UI_DefineInputUnboundElution(obj)
            Figure = figure("Name", "Define gel samples {Input, Unbound, Bound}", ...
                            'NumberTitle', 'off', 'MenuBar', 'none', 'ToolBar', 'none', ...
                            'Resize', 'off', 'Position', [500, 500, 300, 200]);

            Labels = {"Choose input gel", "Choose unbound gel", "Choose bound gel"};
            OptionList = fieldnames(obj.Data);

            Field.LabelWidth = 150;
            Field.DropDownWidth = 150;
            Field.RowHeight = 30;
            Field.TopMargin = 50;
            Field.Spacing = 10;

            DropdownHandles = gobjects(1, length(Labels));

            for i = 1:length(Labels)
                FieldPosition = Field.TopMargin + (length(Labels)-i) * (Field.RowHeight + Field.Spacing);

                uicontrol("Style", "text", "Parent", Figure, ...
                          "String", Labels{i}, 'HorizontalAlignment', 'left', ...
                          'Position', [20, FieldPosition, Field.LabelWidth, Field.RowHeight]);

                DropdownHandles(i) = uicontrol("Style", "popupmenu", "Parent", Figure, ...
                                               "String", OptionList, ...
                                               'Position', [130, FieldPosition, Field.DropDownWidth, Field.RowHeight]);
            end

            uicontrol('Style', 'pushbutton', 'String', 'OK', ...
                      'Position', [50, 10, 80, 30], ...
                      'Callback', @(src, event) onOK());

            uicontrol('Style', 'pushbutton', 'String', 'Cancel', ...
                      'Position', [170, 10, 80, 30], ...
                      'Callback', @(src, event) close(Figure));

            uiwait(Figure);

            function onOK()
                selections = cell(1, length(DropdownHandles));
                for j = 1:length(DropdownHandles)
                    idx = DropdownHandles(j).Value;
                    selections{j} = OptionList{idx};
                end

                if numel(unique(selections)) < numel(selections)
                    errordlg('Each selection must be unique. Please choose different gels.', ...
                             'Duplicate Selection');
                    return; 
                end

                GelTypes = struct( ...
                                  'InputGel', selections{1}, ...
                                  'UnboundGel', selections{2}, ...
                                  'BoundGel', selections{3});
                close(Figure);
            end
        end

        function Parameter = UI_GetVariables(obj, GelTypes)
            Parameter = listdlg("PromptString", "Select properties to view", "SelectionMode", "single", "ListString", fieldnames(obj.Data.(GelTypes.InputGel)));
            Temp_Parameter = fieldnames(obj.Data.(GelTypes.InputGel));
            Parameter = string(Temp_Parameter(Parameter));
        end

        function GelData = SingleGelProperties(obj, Display)
            GelData = obj.Data.(Display.Gel).Description;
            for i = 1:length(Display.Vars)
                switch Display.Vars(i)
                    case "Description"
                        continue
                    otherwise
                        ParameterData = obj.Data.(Display.Gel).(Display.Vars(i));
                        if istable(ParameterData)
                            ParameterData.Properties.VariableNames = {char(Display.Vars(i))};
                        else
                            ParameterData = table(ParameterData, 'VariableNames', {char(Display.Vars(i))});
                        end
                        GelData = [GelData, ParameterData];
                end
            end
        end

        function [GelData] = MultipleGelProperties(obj, Display)
            DataSet = cell(1, length(Display.Gel));
            for i = 1:length(Display.Gel)
                DataSet{i} = obj.Data.(Display.Gel(i));
            end

            MergedTable = obj.MergeTables(DataSet);
            MergedParameters = obj.MergeTableParameters(MergedTable, Display);
            ReducedParameters = obj.ShowDesiredParameters(MergedParameters, Display);
            GelData = obj.ExpandCells(ReducedParameters, Display);
        end

        function CellData = NormalizeTableVars(~, GelData)
            CellData = GelData;
            Fields = fieldnames(CellData);
            for i = 1:length(Fields)
                value = GelData.(Fields{i});
                switch class(CellData.(Fields{i}))
                    case {'cell', 'double'}
                        CellData.(Fields{i}) = value;
                    case 'table'
                        TempVar = value.Properties.VariableNames{1};
                        CellData.(Fields{i}) = value.(TempVar);
                end
            end
            CellData.Description = string(CellData.Description);
            CellData = struct2table(CellData);
        end

        function CombinedTable = MergeTables(obj, DataSet)
            Tables = cellfun(@(d) obj.NormalizeTableVars(d), DataSet, 'UniformOutput', false);
        
            for i = 1:numel(Tables)
                varsToRename = setdiff(Tables{i}.Properties.VariableNames, obj.Keys);
                newNames = varsToRename + "_T" + string(i);
                Tables{i} = renamevars(Tables{i}, varsToRename, newNames);
            end
        
            CombinedTable = Tables{1};
            for i = 2:numel(Tables)
                CombinedTable = outerjoin(CombinedTable, Tables{i}, ...
                    'Keys', obj.Keys, 'MergeKeys', true, 'Type', 'full');
            end
        end

        function DataSet = MergeTableParameters(~, DataSet, Display)
            AllParameters = DataSet.Properties.VariableNames;
            Suffix = '_T\d+$';
            Parameters = unique(regexprep(AllParameters, Suffix, ''));
            for i = 1:length(Parameters)
                BaseParameter = Parameters{i};
                pattern = "^" + BaseParameter + "_T\d+$";
                FauxParameters = AllParameters(~strcmp(AllParameters, BaseParameter) & ~cellfun('isempty', regexp(AllParameters, pattern)));
                if isempty(FauxParameters)
                    continue
                end
                CombinedValue = cell(height(DataSet), 1);
                for row = 1:height(DataSet)
                    CombinedRow = cell(1, length(Display.Gel));
                    for column = 1:length(FauxParameters)
                        ColumnName = FauxParameters{column};
                        Value = DataSet.(ColumnName)(row);
                        if iscell(Value)
                            if isempty(Value) || (numel(Value) == 1 && isempty(Value{1}))
                                CombinedRow{column} = NaN;
                            else
                                CombinedRow{column} = Value;
                            end
                        elseif isempty(Value)
                            CombinedRow{column} = NaN;
                        else
                            CombinedRow{column} = Value;
                        end
                    end
                    if numel(FauxParameters) < length(Display.Gel)
                        CombinedRow(numel(FauxParameters)+1:length(Display.Gel)) = {NaN};
                    end
                    CombinedValue{row} = CombinedRow;
                end
                DataSet.(BaseParameter) = CombinedValue;
                DataSet(:, FauxParameters) = []; 
                AllParameters = DataSet.Properties.VariableNames;
            end
        end

        function SelectedDataSet = ShowDesiredParameters(~, DataSet, Display)
            SelectedDataSet = DataSet(:, Display.Vars);
        end

        function DisplayDataSet = ExpandCells(obj, DataSet, Display)
            for i = 1:length(Display.Vars)
                Parameter = Display.Vars{i};
                if ~ismember(Parameter, obj.Keys)
                    ParameterData = obj.NaN2Inf(DataSet.(Parameter));
                    DataSet.(Parameter) = cellfun(@(x) strjoin(string(x), sprintf(' ')), ParameterData, "UniformOutput", false);
                end
            end
            DisplayDataSet = DataSet(:,Display.Vars);
        end

        function DataSet = NaN2Inf(~, DataSet)
            for i = 1:length(DataSet)
                for j = 1:length(DataSet{i})
                    if isnan(DataSet{i}{j}(:))
                        DataSet{i}{j} = "NaN";
                    end
                end
            end
        end

        function Elutions = ExtractValuesFromMergedTable(~, GelData, Fields)
            TotalGels = length(split(GelData{1,:}{3}));
            Values = zeros(height(GelData), TotalGels);
            for i = 1:height(GelData)
                Values(i,:,:,:) = double(string(split(GelData{i,:}{3})));
            end

            for i = 1:length(Fields)
                Elutions.(Fields{i}) = Values(:,i);
            end
            disp(Elutions')
        end

        function [Proteins] = FoldEnrichment(~, GelData, Elutions, Fields)
            Proteins.Reference = Elutions.(Fields{1});
            Proteins.Unbound = Elutions.(Fields{2});
            Proteins.Bound = Elutions.(Fields{3});

            Proteins.BindingRatio = Proteins.Bound ./ Proteins.Unbound;
            Proteins.Accession = GelData.Accession;
        end

        function [PointsOfInterest] = ScatterProteins(~, Proteins)
            fig = figure;

            ValidIndex = ~isnan(Proteins.Reference) | ~isnan(Proteins.Bound) | ~isnan(Proteins.Unbound);
            Proteins.Reference = Proteins.Reference(ValidIndex);
            Proteins.BindingRatio = Proteins.BindingRatio(ValidIndex);
            Proteins.Accession = Proteins.Accession(ValidIndex);

            hold on;
            PlotInfo = scatter(Proteins.Reference, Proteins.BindingRatio, 144, 'black.');
            yline(2, '--r', '2 Fold or Fold ↑', 'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'top')
            set(gca, 'XScale', 'log')
            xlabel('Gel Input (NSAF)')
            ylabel('Bound / UnBound (NSAF)')
            title("Protein Fold Enrichment")

            PlotInfo.CData = repmat([0 0 0], length(Proteins.Reference), 1);
            PlotInfo.UserData.X = Proteins.Reference;
            PlotInfo.UserData.Y = Proteins.BindingRatio;
            PlotInfo.UserData.Accession = Proteins.Accession;
            PlotInfo.UserData.TextHandles = gobjects(length(Proteins.Reference),1);

            PlotInfo.PickableParts = 'all';
            PlotInfo.ButtonDownFcn = @(src, event) Toggle(src, event);

            function Toggle(src, event)
                ClickPoint = event.IntersectionPoint(1:2);

                X = src.UserData.X;
                Y = src.UserData.Y;
                A = src.UserData.Accession;
                C = src.CData;
                T = src.UserData.TextHandles;

                [~, idx] = min((log10(X) - log10(ClickPoint(1))).^2 + (Y - ClickPoint(2)).^2);

                if all(C(idx,:) == [0 0 0])
                    C(idx,:) = [1 0 0];
                    T(idx) = text(X(idx), Y(idx), A{idx}, ...
                                'FontSize', 8, 'HorizontalAlignment', 'left', ...
                                'VerticalAlignment', 'bottom', 'PickableParts','none');
                else
                    C(idx,:) = [0 0 0];
                    if isvalid(T(idx))
                        delete(T(idx));
                    end
                    T(idx) = gobjects(1);
                end

                src.CData = C;
                src.UserData.TextHandles = T;
            end

            PointsOfInterest = {};
            fig.KeyPressFcn = @(src, event) onKey(src, event, PlotInfo);

            function onKey(~, event, scatterObj)
                if strcmp(event.Key, 'return')
                    redIdx = all(scatterObj.CData == [1 0 0], 2);
                    selected = scatterObj.UserData.Accession(redIdx);

                    PointsOfInterest = selected;

                    if ~isempty(selected)
                        disp('Selected proteins:');
                        disp(selected(:));
                    else
                        disp('No proteins selected.');
                    end

                    uiresume(fig);
                end
            end

            uiwait(fig);
        end

        function [MissingProteins] = MissingFoldEnrichment(~, Proteins)
            MissingProteins.Reference = Proteins.Reference;
            LowerBound = min(Proteins.Reference);
            for i = 1:length(Proteins.Reference)
                if isnan(Proteins.Reference(i))
                    MissingProteins.Accession{i} = Proteins.Accession{i};
                    MissingProteins.Reference(i) = LowerBound;
                    MissingProteins.BindingRatio(i) = -1;
                elseif isnan(Proteins.Bound(i))
                    MissingProteins.Accession{i} = Proteins.Accession{i};
                    MissingProteins.BindingRatio(i) = -2;
                elseif isnan(Proteins.Unbound(i))
                    MissingProteins.Accession{i} = Proteins.Accession{i};
                    MissingProteins.BindingRatio(i) = -3;
                else 
                    continue
                end
            end
        end

        function [PointsOfInterest] = ScatterMissingProteins(~, MissingProteins)
            figure
            Jitter = 0.2;
            MissingProteins.BindingRatio = MissingProteins.BindingRatio + (rand(size(MissingProteins.BindingRatio)) -0.5) * 2 * Jitter;
            scatter(MissingProteins.Reference, MissingProteins.BindingRatio, 144, 'black.')

            ax = gca;
            set(gca, 'XScale', 'log')
            axis tight; 

            OldLabels = [-1, -2, -3];
            NewLabels = ["NaN Input", "NaN Bound", "NaN UnBound"];
            
            DisplayLabels = OldLabels(ismember(OldLabels, ax.YTick));
            DisplayLabels = sort(DisplayLabels);
            ax.YTick = DisplayLabels;
            Labels = strings(size(DisplayLabels));
            for i = 1:length(DisplayLabels)
                Index = OldLabels == DisplayLabels(i);
                Labels(i) = NewLabels(Index);
            end
            ax.YTickLabel = Labels;
            xlabel("Gel Input (NSAF)")
            ylim([min(DisplayLabels) - 0.5, max(DisplayLabels) + 0.5])
            
            text(1, 1, 'Jitter applied vertically', 'Units', 'normalized', ...
                 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top')

            for i = 1:length(MissingProteins.Reference)
                text(MissingProteins.Reference(i), MissingProteins.BindingRatio(i), string(MissingProteins.Accession{i}), ...
                     'FontSize', 8, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom')
            end

            % Clickable function
            % PointsOfInterest printout function
            PointsOfInterest = 1;
        end

        function [Info] = UniProtIdentification(~, ProteinIDs)
            Info.Accession    = ProteinIDs(:);
            Info.ProteinName  = repmat({''}, length(ProteinIDs), 1);
            Info.Function     = repmat({''}, length(ProteinIDs), 1);

            Header = ["Accept", "application/json" 
                      "Accept-Encoding", "identity"
                     ];
            Options = weboptions("ContentType", 'json', ...
                                 "Timeout", 15, ...
                                 "HeaderFields", Header ...
                                );

            for i = 1:length(ProteinIDs)
                UniProtID = ProteinIDs{i};
                url = ['https://rest.uniprot.org/uniprotkb/' UniProtID '.json'];

                try 
                    WebData = webread(url, Options);
                    Description = WebData.proteinDescription;
                    if isfield(Description, 'recommendedName')
                        Info.ProteinName{i} = Description.recommendedName.fullName.value;
                    elseif isfield(Description, 'submissionNames')
                        Info.ProteinName{i} = Description.submissionNames.fullName.value;
                    else
                        Info.ProteinName{i} = 'N/A';
                    end

                    Info.Function{i} = 'N/A';
                    if isfield(WebData, 'comments') && ~isempty(WebData.comments)
                        comments = WebData.comments;

                        if isstruct(comments)
                            if isscalar(comments)
                                comments = {comments};
                            else
                                comments = num2cell(comments);
                            end
                        end

                        for j = 1:length(comments)
                            comment = comments{j};
                            if strcmp(comment.commentType, 'FUNCTION')
                                if isfield(comment.texts, 'value')
                                    Info.Function{i} = comment.texts.value;
                                    break;
                                end
                            end
                        end
                    end

                catch ME
                    Info.ProteinName{i} = 'Error';
                    Info.Function{i} = ['Error: ' ME.message];
                end
            end
            disp([Info.Accession, Info.ProteinName, Info.Function])
        end

    end
end