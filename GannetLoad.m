function MRS_struct = GannetLoad(varargin)
% Gannet 3.1 GannetLoad
% Started by RAEE Nov. 5, 2012
% Updates by MGS, MM, GO 2016-2020

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Workflow summary
%   1. Pre-initialise
%   2. Determine data parameters from headers
%   3. Load data from files
%   4. Reconstruction of coil-sensitivity maps (PRIAM only)
%   5. Apply appropriate pre-processing
%   6. Build GannetLoad output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MRS_struct.version.Gannet = '3.2.0-rc.1';
MRS_struct.version.load   = '201124';
VersionCheck(0, MRS_struct.version.Gannet);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   -1. Determine if the GannetGUI was used and adjust inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

runWithGUI = 0;

if isnumeric(varargin{1})
    runWithGUI = 1;
    configFilePath = varargin{2};
    varargin = varargin(3:end);
    nargin = length(varargin);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   0. Parse the input arguments and check for typos
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

metabfile = varargin{1};
MRS_struct.metabfile = metabfile;
missing = 0;
for filecheck = 1:length(metabfile)
    if ~exist(metabfile{filecheck}, 'file')
        fprintf('\nThe file ''%s'' (%d) is missing. Typo?\n', metabfile{filecheck}, filecheck);
        missing = 1;
    end
end

if nargin > 1 && ~isempty(varargin{2})
    waterfile = varargin{2};
    MRS_struct.waterfile = waterfile;
    for filecheck = 1:length(waterfile)
        if ~exist(waterfile{filecheck}, 'file')
            fprintf('\nThe water reference file ''%s'' (%d) is missing. Typo?\n', waterfile{filecheck}, filecheck);
            missing = 1;
        end
    end
end

if missing
    fprintf('\n');
    error('Not all the input files can be found. Exiting...');
end

if nargin == 3
    mode = varargin{3};
    if ~any(strcmp(mode, {'batch','join'}))
        error('The third input argument must be either ''batch'' or ''join''.');
    end
end

if ~exist('mode','var')
    mode = 'batch';
end

if nargin == 4
    if isnumeric(varargin{4})
        trimAvgs = varargin{4};
    else
        [~,~,ext] = fileparts(varargin{4});
        assert(strcmpi(ext,'.csv'), 'The fourth argument must be a .csv file.');
        trimAvgs = readtable(varargin{4});
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   1. Pre-initialise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~runWithGUI
    MRS_struct = GannetPreInitialise(MRS_struct);
else
    try
        MRS_struct = GannetPreInitialiseGUIVersion(MRS_struct, configFilePath);
    catch
        warning(strcat("Error reading configuration file: ", configFilePath, ". Make sure that it was created using the GUI or has all fields properly named. Stopping GannetLoad."))
        return
    end
end

CheckTargets(MRS_struct);

if MRS_struct.p.PRIAM
    vox = MRS_struct.p.vox;
else
    vox = MRS_struct.p.vox(1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   2. Determine data parameters from header
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Discern input data format
MRS_struct = DiscernDataType(metabfile{1}, MRS_struct);

% Determine number of provided water-suppressed files in the batch
switch mode
    case 'batch'
        numscans = numel(metabfile);
        % For Siemens RDA, each acquisition has two files so correct the number
        if strcmp(MRS_struct.p.vendor,'Siemens_rda')
            numscans = numscans/2;
        end
        numfilesperscan = 1;
    case 'join'
        numscans = 1;
        numfilesperscan = numel(metabfile);
        % For Siemens RDA, each acquisition has two files so correct the number
        if strcmp(MRS_struct.p.vendor,'Siemens_rda')
            numfilesperscan = numfilesperscan/2;
        end
        fprintf('\nRunning GannetLoad in ''join'' mode. Joining FIDs from %i separate files...', numfilesperscan);
end

% Determine number of provided water-unsuppressed files in the batch
if exist('waterfile','var')
    MRS_struct.p.Reference_compound = 'H2O';
    numwaterscans = numel(waterfile);
    switch mode
        case 'batch'
            if numwaterscans ~= numscans
                error ('Number of water-unsuppressed files does not match number of water-suppressed files.');
            end
    end
else
    MRS_struct.p.Reference_compound = 'Cr';
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   3. Load data from files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MRS_struct.p.numscans = numscans;
errorReport = cell(1);
tryCatchInd = 1;

warning('off','stats:nlinfit:ModelConstantWRTParam');
warning('off','stats:nlinfit:IllConditionedJacobian');
warning('off','MATLAB:rankDeficientMatrix');

for ii = 1:numscans % Loop over all files in the batch (from metabfile)
    
    [~,b,c] = fileparts(metabfile{ii});
    if strcmp(mode,'join')
        fprintf('\nLoading %s and %i other files...', [b c], numfilesperscan - 1);
    else
        fprintf('\nLoading %s...', [b c]);
    end
    
    f = dir(metabfile{ii});
    if f.bytes/1e6 > 150
        fprintf('\nLarge file detected (%.2f MB). Please wait...', f.bytes/1e6);
    end
    
    try % pass to next dataset if errors occur
        
        MRS_struct.ii = ii;
        
        switch MRS_struct.p.vendor
            
            case 'GE'
                
                MRS_struct = GERead(MRS_struct, metabfile{ii});
                MRS_struct.p.Reference_compound = 'H2O';
                %MRS_struct.fids.data = MRS_struct.fids.data * MRS_struct.p.nrows(ii)/MRS_struct.p.Navg(ii);                
                
            case 'Siemens_twix'
                
                if exist('waterfile','var')
                    if strcmp(mode,'batch')
                        MRS_struct = SiemensTwixRead(MRS_struct, metabfile{ii}, waterfile{ii});
                    else
                        % Load each input file and append the FIDs
                        MRS_struct = SiemensTwixRead(MRS_struct, metabfile{ii}, waterfile{ii});
                        for kk = 2:numfilesperscan
                            sub_MRS_struct = SiemensTwixRead(MRS_struct, metabfile{kk}, waterfile{ii});
                            MRS_struct.fids.data = [MRS_struct.fids.data sub_MRS_struct.fids.data];
                        end
                    end
                    % Correct the total number of averages
                    MRS_struct.p.nrows = MRS_struct.p.nrows * numfilesperscan;
                    MRS_struct.p.Navg = MRS_struct.p.Navg * numfilesperscan;
                else
                    if strcmp(mode,'batch')
                        MRS_struct = SiemensTwixRead(MRS_struct, metabfile{ii});
                    else
                        % Load each input file and append the FIDs
                        MRS_struct = SiemensTwixRead(MRS_struct, metabfile{1});
                        for kk = 2:numfilesperscan
                            sub_MRS_struct = SiemensTwixRead(MRS_struct, metabfile{kk});
                            MRS_struct.fids.data = [MRS_struct.fids.data sub_MRS_struct.fids.data];
                        end
                    end
                    % Correct the total number of averages
                    MRS_struct.p.nrows = MRS_struct.p.nrows * numfilesperscan;
                    MRS_struct.p.Navg = MRS_struct.p.Navg * numfilesperscan;
                end
                
            case 'Siemens_dicom'
                
                if exist('waterfile','var')
                    MRS_struct = SiemensDICOMRead(MRS_struct,metabfile{ii}, waterfile{ii});
                else
                    MRS_struct = SiemensDICOMRead(MRS_struct, metabfile{ii});
                end
                
            case 'dicom'
                
                if exist('waterfile','var')
                    MRS_struct = DICOMRead(MRS_struct,metabfile{ii}, waterfile{ii});
                else
                    MRS_struct = DICOMRead(MRS_struct, metabfile{ii});
                end
                
                switch MRS_struct.p.ONOFForder
                    % Not sure whether this is always the case, but the CMRR
                    % sequence appears to go OFF-OFF-ON-ON in the DICOM
                    % sorting?! Fixing this hard for now.
                    case 'onfirst'
                        %if strcmp(MRS_struct.p.seq,'""%CustomerSeq%\eja_svs_mpress""')
                        %    MRS_struct.fids.ON_OFF = repmat([1 1 0 0],[1 MRS_struct.p.Navg(ii)/4]);
                        %    MRS_struct.fids.ON_OFF = MRS_struct.fids.ON_OFF(:).';
                        %else
                            MRS_struct.fids.ON_OFF = repmat([1 0],[1 size(MRS_struct.fids.data,2)/2]);
                            MRS_struct.fids.ON_OFF = MRS_struct.fids.ON_OFF(:).';
                        %end
                    case 'offfirst'
                        %if strcmp(MRS_struct.p.seq,'""%CustomerSeq%\eja_svs_mpress""')
                        %    MRS_struct.fids.ON_OFF = repmat([0 0 1 1],[1 MRS_struct.p.Navg(ii)/4]);
                        %    MRS_struct.fids.ON_OFF = MRS_struct.fids.ON_OFF(:).';
                        %else
                            MRS_struct.fids.ON_OFF = repmat([0 1],[1 size(MRS_struct.fids.data,2)/2]);
                            MRS_struct.fids.ON_OFF = MRS_struct.fids.ON_OFF(:).';
                        %end
                end
                
            case 'Siemens_rda'
                
                if exist('waterfile','var')
                    if strcmp(mode,'batch')
                        MRS_struct = SiemensRead(MRS_struct, metabfile{ii*2}, metabfile{ii*2-1}, waterfile{ii});
                    else
                        % Load each input file and append the FIDs
                        MRS_struct = SiemensRead(MRS_struct, metabfile{2}, metabfile{1}, waterfile{ii});
                        for kk = 2:numfilesperscan
                            sub_MRS_struct = SiemensRead(MRS_struct, metabfile{kk*2}, metabfile{kk*2-1}, waterfile{ii});
                            MRS_struct.fids.data = [MRS_struct.fids.data sub_MRS_struct.fids.data];
                        end
                    end
                    % Correct the total number of averages
                    MRS_struct.p.Navg = MRS_struct.p.Navg * numfilesperscan;
                    MRS_struct.p.Nwateravg = 1;
                else
                    if strcmp(mode,'batch')
                        MRS_struct = SiemensRead(MRS_struct, metabfile{ii*2}, metabfile{ii*2-1});
                    else
                        % Load each input file and append the FIDs
                        MRS_struct = SiemensRead(MRS_struct, metabfile{2}, metabfile{1});
                        for kk = 2:numfilesperscan
                            sub_MRS_struct = SiemensRead(MRS_struct, metabfile{kk*2}, metabfile{kk*2-1});
                            MRS_struct.fids.data = [MRS_struct.fids.data sub_MRS_struct.fids.data];
                        end
                    end
                    % Correct the total number of averages
                    MRS_struct.p.Navg = MRS_struct.p.Navg * numfilesperscan;
                end
                
            case 'Philips'
                
                if exist('waterfile','var')
                    MRS_struct = PhilipsRead(MRS_struct, metabfile{ii}, waterfile{ii});
                else
                    MRS_struct = PhilipsRead(MRS_struct, metabfile{ii});
                end
                
            case 'Philips_data'
                
                % If a water reference scan is acquired, it is saved as a mix
                % in the DATA/LIST files. Later: add option to provide an additional
                % water reference file (i.e. short-TE).
                if nargin > 1
                    MRS_struct = PhilipsRead_data(MRS_struct, metabfile{ii}, waterfile{ii});
                else
                    MRS_struct = PhilipsRead_data(MRS_struct, metabfile{ii});
                end
                if isfield(MRS_struct.fids, 'data_water')
                    MRS_struct.p.Reference_compound = 'H2O';
                else
                    MRS_struct.p.Reference_compound = 'Cr';
                end
                
            case 'Philips_raw'
                
                MRS_struct = PhilipsRawLoad(MRS_struct,metabfile{ii},3,0);
                MRS_struct.fids.data = conj(squeeze(MRS_struct.multivoxel.allsignals(:,:,1,:)));
                if exist('waterfile','var')
                    MRS_struct.p.Reference_compound = 'H2O';
                end
                
        end % end of vendor switch loop for data load
        
        if nargin == 4
            if isnumeric(varargin{4})
                if size(trimAvgs,1) > 1
                    t_start = trimAvgs(ii,1);
                    t_end   = trimAvgs(ii,2);
                else
                    t_start = trimAvgs(1,1);
                    t_end   = trimAvgs(1,2);
                end
            else
                fileInd = find(strcmpi(metabfile(ii), trimAvgs{:,1}));
                t_start = trimAvgs{fileInd,2};
                t_end   = trimAvgs{fileInd,3};
            end
            % Make sure t_start is an odd number and t_end is an even number
            if ~rem(t_start,2)
                t_start = t_start - 1;
            end
            if rem(t_end,2)
                t_end = t_end + 1;
            end
            MRS_struct.fids.data  = MRS_struct.fids.data(:,t_start:t_end);
            MRS_struct.p.Navg(ii) = size(MRS_struct.fids.data,2);
        end
        
        if ~strcmp(MRS_struct.p.vendor, 'dicom')
            % Determine order of ON and OFF acquisitions
            MRS_struct = SpecifyOnOffOrder(MRS_struct);
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %   4. Reconstruction of coil-sensitivity maps
        %      (PRIAM only)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % if a PRIAM dataset is processed, load the coil reference scan and
        % calculate the SENSE reconstruction matrix here
        if MRS_struct.p.PRIAM
            MRS_struct = senseRecon(MRS_struct);
            PRIAMData = zeros(length(MRS_struct.p.vox),MRS_struct.p.Navg,MRS_struct.p.npoints);
            PRIAMWaterData = zeros(length(MRS_struct.p.vox),MRS_struct.p.Nwateravg,MRS_struct.p.npoints);
            for kk = 1:MRS_struct.p.Navg
                PRIAMData(:,kk,:) = MRS_struct.p.SENSE.U * squeeze(MRS_struct.fids.data(:,kk,:));
                % Phase by multiplying with normalized complex conjugate of first point
                conj_norm = conj(PRIAMData(:,kk,1)) ./ abs(conj(PRIAMData(:,kk,1)));
                PRIAMData(:,kk,:) = PRIAMData(:,kk,:) .* repmat(conj_norm, [1 1 MRS_struct.p.npoints]);
            end
            for kk = 1:MRS_struct.p.Nwateravg
                PRIAMWaterData(:,kk,:) = MRS_struct.p.SENSE.U * squeeze(MRS_struct.fids.data_water(:,kk,:));
                % Phase by multiplying with normalized complex conjugate of first point
                conj_norm = conj(PRIAMWaterData(:,kk,1)) ./ abs(conj(PRIAMWaterData(:,kk,1)));
                PRIAMWaterData(:,kk,:) = PRIAMWaterData(:,kk,:) .* repmat(conj_norm, [1 1 MRS_struct.p.npoints]);
            end
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %   5. Apply appropriate pre-processing
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        for kk = 1:length(vox) % loop over number of voxels
            
            % Select data from first voxel
            if MRS_struct.p.PRIAM
                MRS_struct.fids.data = squeeze(-PRIAMData(kk,:,:))';
                MRS_struct.fids.data_water = squeeze(PRIAMWaterData(kk,:,:))';
            end
            
            % Zero-fill to obtain nominal spectral resolution of 0.061 Hz/point
            MRS_struct.p.ZeroFillTo(ii) = round(32768 / 2000 * MRS_struct.p.sw(ii));
            MRS_struct.p.zf(ii) = MRS_struct.p.ZeroFillTo(ii) / MRS_struct.p.npoints(ii);
            time = (1:size(MRS_struct.fids.data,1)) / MRS_struct.p.sw(ii);
            
            % Finish processing water data
            if strcmp(MRS_struct.p.Reference_compound,'H2O')
                
                % Some data formats have not averaged the water signal
                % (change this later so that it happens during data loading instead)
                if strcmp(MRS_struct.p.vendor,'dicom')
                    MRS_struct.fids.data_water = mean(MRS_struct.fids.data_water,2);
                elseif strcmp(MRS_struct.p.vendor,'Philips_raw')
                    MRS_struct.fids.data_water = mean(MRS_struct.fids.data_water(kk,:,:),2);
                end
                
                % Performing eddy current corrrection on the water-suppressed data
                if MRS_struct.p.metab_ECC
                    MRS_struct.fids.data = EddyCurrentCorrection(MRS_struct.fids.data, MRS_struct.fids.data_water);
                end
                
                % Performing eddy current corrrection on the unsuppressed water data
                if MRS_struct.p.water_ECC
                    MRS_struct.fids.data_water = EddyCurrentCorrection(MRS_struct.fids.data_water, MRS_struct.fids.data_water);
                end
                
                % Line-broadening, zero-filling and FFT
                % Water data may have different bandwidth
                if isfield(MRS_struct.p, 'sw_water')
                    time_water = (1:size(MRS_struct.fids.data_water,1)) / MRS_struct.p.sw_water(ii);
                else
                    time_water = (1:size(MRS_struct.fids.data_water,1)) / MRS_struct.p.sw(ii);
                end
                MRS_struct.fids.data_water = MRS_struct.fids.data_water .* exp(-time_water' * MRS_struct.p.LB * pi);
                MRS_struct.spec.(vox{kk}).water(ii,:) = fftshift(fft(MRS_struct.fids.data_water, MRS_struct.p.ZeroFillTo(ii), 1),1);
                
            end % end of H2O reference loop
            
            % Global zero-order phase correction
            if ~MRS_struct.p.phantom
                MRS_struct.fids.data = PhaseCorrection(MRS_struct.fids.data, MRS_struct);
            end
            
            % Line-broadening, zero-filling and FFT
            AllFramesFT = MRS_struct.fids.data .* repmat(exp(-time' * MRS_struct.p.LB * pi), [1 size(MRS_struct.fids.data,2)]);
            AllFramesFT = fftshift(fft(AllFramesFT, MRS_struct.p.ZeroFillTo(ii), 1),1);
            
            % Work out ppm axis
            freqRange = MRS_struct.p.sw(ii) / MRS_struct.p.LarmorFreq(ii);
            if MRS_struct.p.phantom
                F0 = 4.8;
            else
                F0 = 4.68;
            end
            MRS_struct.spec.freq = (MRS_struct.p.ZeroFillTo(ii) + 1 - (1:1:MRS_struct.p.ZeroFillTo(ii))) / MRS_struct.p.ZeroFillTo(ii) * freqRange + F0 - freqRange/2;
            
            MRS_struct.p.df(ii)             = abs(MRS_struct.spec.freq(1) - MRS_struct.spec.freq(2));
            MRS_struct.p.SpecRes(ii)        = MRS_struct.p.sw(ii) / MRS_struct.p.npoints(ii);
            MRS_struct.p.SpecResNominal(ii) = MRS_struct.p.sw(ii) / MRS_struct.p.ZeroFillTo(ii);
            MRS_struct.p.Tacq(ii)           = 1/MRS_struct.p.SpecRes(ii);
            
            % Frame-by-frame determination of frequency of residual water or Cr (if HERMES/HERCULES or GSH editing)
            if MRS_struct.p.HERMES || any(strcmp(MRS_struct.p.target,'GSH'))
                F0freqRange = MRS_struct.spec.freq - 3.02 >= -0.15 & MRS_struct.spec.freq - 3.02 <= 0.15;
            else
                F0freqRange = MRS_struct.spec.freq - F0 >= -0.2 & MRS_struct.spec.freq - F0 <= 0.2;
            end
            [~,FrameMaxPos] = max(abs(real(AllFramesFT(F0freqRange,:))),[],1);
            F0freqRange = MRS_struct.spec.freq(F0freqRange);
            MRS_struct.spec.F0freq{ii} = F0freqRange(FrameMaxPos);
            
            % Estimate average amount of F0 offset
            if MRS_struct.p.HERMES || any(strcmp(MRS_struct.p.target,'GSH'))
                MRS_struct.out.AvgDeltaF0(ii) = mean(F0freqRange(FrameMaxPos) - 3.02);
            elseif any(strcmp(MRS_struct.p.vendor,{'Siemens_rda','Siemens_twix','Siemens_dicom'}))
                MRS_struct.out.AvgDeltaF0(ii) = mean(F0freqRange(FrameMaxPos) - 4.7); % Siemens assumes 4.7 ppm as F0
            else
                MRS_struct.out.AvgDeltaF0(ii) = mean(F0freqRange(FrameMaxPos) - F0);
            end
            
            % Use frame-by-frame frequency of Cr for RobustSpectralRegistration
            F0freqRange = MRS_struct.spec.freq - 3.02 >= -0.15 & MRS_struct.spec.freq - 3.02 <= 0.15;
            [~,FrameMaxPos] = max(abs(real(AllFramesFT(F0freqRange,:))),[],1);
            F0freqRange = MRS_struct.spec.freq(F0freqRange);
            MRS_struct.spec.F0freq2{ii} = F0freqRange(FrameMaxPos);
            
            % Frame-by-frame alignment
            switch MRS_struct.p.alignment
                case {'Cr','Cho','NAA'}
                    [AllFramesFTrealign, MRS_struct] = AlignUsingPeak(AllFramesFT, MRS_struct);
                case 'H2O'
                    [AllFramesFTrealign, MRS_struct] = AlignUsingH2O(AllFramesFT, MRS_struct);
                case 'SpecReg'
                    [AllFramesFTrealign, MRS_struct] = SpectralRegistration(MRS_struct,0);
                case 'SpecRegDual'
                    %Dual-channel Spectral Registration is applied separately to ON and OFF and they are coregistered after...
                    [AllFramesFTrealign, MRS_struct] = SpectralRegistration(MRS_struct,0,1);
                case 'SpecRegHERMES'
                    [AllFramesFTrealign, MRS_struct] = SpectralRegistrationHERMES(MRS_struct);
                case 'RobustSpecReg'
                    [AllFramesFTrealign, MRS_struct] = RobustSpectralRegistration(MRS_struct);
                case 'none'
                    % do nothing
                    AllFramesFTrealign = AllFramesFT;
                    MRS_struct.out.reject{ii} = zeros(1,size(AllFramesFT,2));
                otherwise
                    filepath = fullfile(fileparts(which(mfilename('fullpath'))), 'GannetPreInitialise.m');
                    msg = 'FPC parameter in GannetPreInitialise.m not recognized. Check spelling.';
                    msg = hyperlink(['matlab: opentoline(''' filepath ''', 22, 0)'], 'FPC parameter in GannetPreInitialise.m not recognized', msg);
                    error(msg);
            end
            
            MRS_struct.spec.AllFramesFT = AllFramesFT;
            MRS_struct.spec.AllFramesFTrealign = AllFramesFTrealign;
            
            % Average subspectra and generate DIFF spectra
            MRS_struct = SignalAveraging(MRS_struct, AllFramesFT, AllFramesFTrealign, ii, kk, vox);
            
            % Remove residual water from diff and diff_noalign spectra using HSVD
            if MRS_struct.p.water_removal
                
                for jj = 1:length(MRS_struct.p.target)
                    if jj == 1
                        fprintf('\nRemoving the residual water signal using HSVD...\n');
                    end
                    
                    % Convert DIFF spectra to time domain, apply water filter, convert back to frequency domain
                    fids.diff = WaterRemovalHSVD(ifft(ifftshift(MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{jj}).diff(ii,:).')), ...
                        MRS_struct.p.sw(ii)/1e3, 8, -0.08, 0.08, 0, 2048);
                    MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{jj}).diff(ii,:) = fftshift(fft(fids.diff));
                    
                    fids.diff_noalign = WaterRemovalHSVD(ifft(ifftshift(MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{jj}).diff_noalign(ii,:).')), ...
                        MRS_struct.p.sw(ii)/1e3, 8, -0.08, 0.08, 0, 2048);
                    MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{jj}).diff_noalign(ii,:) = fftshift(fft(fids.diff_noalign));
                    
                    % Need to perform baseline correction on filtered data
                    freqbounds = MRS_struct.spec.freq <= 8 & MRS_struct.spec.freq >= 7;
                    baseMean_diff = mean(real(MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{jj}).diff(ii,freqbounds)));
                    baseMean_diffnoalign = mean(real(MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{jj}).diff_noalign(ii,freqbounds)));
                    
                    MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{jj}).diff(ii,:) = MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{jj}).diff(ii,:) - baseMean_diff;
                    MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{jj}).diff_noalign(ii,:) = MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{jj}).diff_noalign(ii,:) - baseMean_diffnoalign;
                end
                
            end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %   6. Build GannetLoad output
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if ishandle(101)
                clf(101);
            end
            h = figure(101);
            % Open figure in center of screen
            scr_sz = get(0,'ScreenSize');
            fig_w = 1000;
            fig_h = 707;
            set(h,'Position',[(scr_sz(3)-fig_w)/2, (scr_sz(4)-fig_h)/2, fig_w, fig_h]);
            set(h,'Color',[1 1 1]);
            figTitle = 'GannetLoad Output';
            set(h,'Name',figTitle,'Tag',figTitle,'NumberTitle','off');
            
            % Top left
            if length(MRS_struct.p.target) == 3
                subplot(5,2,1:2:5);
            else
                subplot(2,2,1);
            end
            PlotPrePostAlign(MRS_struct, vox, ii, kk);
            
            % Top right
            if MRS_struct.p.phantom
                if MRS_struct.p.HERMES
                    F0 = 3.02;
                else
                    F0 = 4.8;
                end
            elseif MRS_struct.p.HERMES || any(strcmp(MRS_struct.p.target,'GSH'))
                F0 = 3.02;
            else
                F0 = 4.68;
            end
            
            subplot(2,2,2);
            if ~MRS_struct.p.weighted_averaging && size(MRS_struct.fids.data,2) >= 4
                rejectframesplot = (1./MRS_struct.out.reject{ii}) .* MRS_struct.spec.F0freq{ii};
            end
            hold on;
            plot([1 size(MRS_struct.fids.data,2)], [F0 F0], '-k')
            plot([1 size(MRS_struct.fids.data,2)], [F0-0.04 F0-0.04], '--k')
            plot([1 size(MRS_struct.fids.data,2)], [F0+0.04 F0+0.04], '--k');
            plot(1:size(MRS_struct.fids.data,2), MRS_struct.spec.F0freq{ii}', 'Color', 'b');
            if ~MRS_struct.p.weighted_averaging && size(MRS_struct.fids.data,2) >= 4
                plot(1:size(MRS_struct.fids.data,2), rejectframesplot, 'ro');
            end
            hold off;
            if MRS_struct.p.HERMES || any(strcmp(MRS_struct.p.target,'GSH'))
                text(size(MRS_struct.fids.data,2) + 0.025*size(MRS_struct.fids.data,2), F0, {'Nominal','Cr freq.'}, 'FontSize', 8);
            else
                text(size(MRS_struct.fids.data,2) + 0.025*size(MRS_struct.fids.data,2), F0, {'Nominal','water freq.'}, 'FontSize', 8);
            end
            set(gca,'TickDir','out','box','off','XLim',[1 size(MRS_struct.fids.data,2)], ...
                'YLim',[min([F0-0.06 MRS_struct.spec.F0freq{ii}-0.005]) max([F0+0.06 MRS_struct.spec.F0freq{ii}+0.005])]);
            if size(MRS_struct.fids.data,2) == 2
                set(gca,'XTick',[1 2]);
            end
            xlabel('average');
            ylabel('ppm');
            if MRS_struct.p.HERMES || any(strcmp(MRS_struct.p.target,'GSH'))
                title('Cr Frequency');
            else
                title('Water Frequency');
            end
            
            % Bottom left
            if length(MRS_struct.p.target) == 3
                subplot(3,2,5);
            else
                subplot(2,2,3);
            end
            if ~strcmp(MRS_struct.p.alignment, 'none')
                CrFitLimLow = 2.72;
                CrFitLimHigh = 3.12;
                plotrange = MRS_struct.spec.freq <= CrFitLimHigh & MRS_struct.spec.freq >= CrFitLimLow;
                CrFitRange = sum(plotrange);
                plotrealign = [real(AllFramesFT(plotrange,:)); real(AllFramesFTrealign(plotrange,:))];
                % Don't display rejects
                if ~MRS_struct.p.weighted_averaging && size(MRS_struct.fids.data,2) >= 4
                    plotrealign(CrFitRange+1:end,(MRS_struct.out.reject{ii} == 1)) = min(plotrealign(:));
                end
                imagesc(plotrealign);
                colormap('parula');
                title({'Cr Frequency','(pre- and post-alignment)'});
                xlabel('average');
                ylabel('ppm');
                set(gca, 'YTick', [CrFitRange * (CrFitLimHigh - 3.02) / (CrFitLimHigh - CrFitLimLow) ...
                                   CrFitRange ...
                                   CrFitRange + CrFitRange * (CrFitLimHigh - 3.02) / (CrFitLimHigh - CrFitLimLow) ...
                                   CrFitRange * 2], ...
                    'YTickLabel', [3.02 CrFitLimLow 3.02 CrFitLimLow], ...
                    'XLim', [1 size(MRS_struct.fids.data,2)], ...
                    'YLim', [1 CrFitRange * 2], ...
                    'TickDir','out','box','off');
                if size(MRS_struct.fids.data,2) == 2
                    set(gca, 'XTick', [1 2]);
                end
                % Add in labels for pre/post
                text(size(plotrealign,2)/18*17, 0.4*size(plotrealign,1), 'PRE', 'Color', [1 1 1], 'HorizontalAlignment', 'right');
                text(size(plotrealign,2)/18*17, 0.9*size(plotrealign,1), 'POST', 'Color', [1 1 1], 'HorizontalAlignment', 'right');
            else
                tmp = 'No alignment';
                text(0.5, 0.5, tmp, 'FontName', 'Arial', 'FontSize', 20, 'HorizontalAlignment', 'center');
                set(get(gca,'XAxis'),'Visible','off');
                set(get(gca,'YAxis'),'Visible','off');
            end
            
            % Bottom right
            subplot(2,2,4);
            axis off;
            
            if strcmp(MRS_struct.p.vendor, 'Siemens_rda')
                [~,tmp,tmp2] = fileparts(MRS_struct.metabfile{ii*2-1});
            else
                [~,tmp,tmp2] = fileparts(MRS_struct.metabfile{ii});
            end
            fname = [tmp tmp2];
            if length(fname) > 30
                fname = [fname(1:12) '...' fname(end-11:end)];
            end
            text(0.25, 0.9, 'Filename: ', 'FontName', 'Arial', 'FontSize', 13, 'HorizontalAlignment', 'right');
            text(0.275, 0.9, fname, 'FontName', 'Arial', 'FontSize', 13, 'Interpreter', 'none');
            
            tmp = [num2str(MRS_struct.p.Navg(ii)) ' averages'];
            text(0.25, 0.8, 'N_{avg}: ', 'FontName', 'Arial', 'FontSize', 13, 'HorizontalAlignment', 'right');
            text(0.275, 0.8, tmp, 'FontName', 'Arial', 'FontSize', 13);
            
            tmp = [num2str(prod(MRS_struct.p.voxdim(ii,:))/1e3) ' mL'];
            text(0.25, 0.7, 'Volume: ', 'FontName', 'Arial', 'FontSize', 13, 'HorizontalAlignment', 'right');
            text(0.275, 0.7, tmp, 'FontName', 'Arial', 'FontSize', 13);
            
            text(0.25, 0.6, 'Alignment: ', 'FontName', 'Arial', 'FontSize', 13, 'HorizontalAlignment', 'right');
            if strcmp(MRS_struct.p.alignment,'RobustSpecReg') && MRS_struct.p.use_prealign_ref
                text(0.275, 0.6, [MRS_struct.p.alignment ' (PreAlignRef)'], 'FontName', 'Arial', 'FontSize', 13);
            else
                text(0.275, 0.6, MRS_struct.p.alignment, 'FontName', 'Arial', 'FontSize', 13);
            end
            
            tmp = [num2str(MRS_struct.p.LB) ' Hz'];
            text(0.25, 0.5, 'LB: ', 'FontName', 'Arial', 'FontSize', 13, 'HorizontalAlignment', 'right');
            text(0.275, 0.5, tmp, 'FontName', 'Arial', 'FontSize', 13);
            
            text(0.25, 0.4, 'Rejects: ', 'FontName', 'Arial', 'FontSize', 13, 'HorizontalAlignment', 'right');
            if MRS_struct.p.weighted_averaging
                text(0.275, 0.4, 'n/a - weighted averaging used', 'FontName', 'Arial', 'FontSize', 13);
            else
                text(0.275, 0.4, num2str(sum(MRS_struct.out.reject{ii})), 'FontName', 'Arial', 'FontSize', 13);
            end
            
            text(0.25, 0.3, 'LoadVer: ', 'FontName', 'Arial', 'FontSize', 13, 'HorizontalAlignment', 'right');
            text(0.275, 0.3, MRS_struct.version.load, 'FontName', 'Arial', 'FontSize', 13);
            
            % Gannet logo
            Gannet_logo = fullfile(fileparts(which('GannetLoad')), 'Gannet3_logo.png');
            I = imread(Gannet_logo,'png','BackgroundColor',[1 1 1]);
            axes('Position', [0.825, 0.05, 0.125, 0.125]);
            imshow(I);
            text(0.9, 0, MRS_struct.version.Gannet, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'left');
            axis off;
            axis square;
            
            % Gannet documentation
            axes('Position', [(1-0.9)/2, 0.025, 0.9, 0.15]);
            str = 'For complete documentation, please visit: https://markmikkelsen.github.io/Gannet-docs';
            text(0.5, 0, str, 'FontName', 'Arial', 'FontSize', 11, 'HorizontalAlignment', 'center');
            axis off;
            axis square;
            
            % For Philips .data
            if strcmp(MRS_struct.p.vendor,'Philips_data')
                fullpath = MRS_struct.metabfile{ii};
                fullpath = regexprep(fullpath, '.data', '_data');
                fullpath = regexprep(fullpath, '\', '_');
                fullpath = regexprep(fullpath, '/', '_');
            end
            
            if strcmp(MRS_struct.p.vendor,'Siemens_rda')
                [~,metabfile_nopath] = fileparts(MRS_struct.metabfile{ii*2-1});
            else
                [~,metabfile_nopath] = fileparts(MRS_struct.metabfile{ii});
            end
            
            if any(strcmp(listfonts,'Arial'))
                set(findall(h,'-property','FontName'),'FontName','Arial');
            end
            set(findall(h,'-property','XColor'),'XColor',[0 0 0]);
            set(findall(h,'-property','YColor'),'YColor',[0 0 0]);
            
            % Create output folder
            if ~exist(fullfile(pwd, 'GannetLoad_output'),'dir')
                mkdir(fullfile(pwd, 'GannetLoad_output'));
            end
            
            % Save PDF
            set(h,'PaperUnits','inches');
            set(h,'PaperSize',[11 8.5]);
            set(h,'PaperPosition',[0 0 11 8.5]);
            
            if strcmp(MRS_struct.p.vendor,'Philips_data')
                pdfname = fullfile(pwd, 'GannetLoad_output', [fullpath '_' vox{kk} '_load.pdf']);
            else
                pdfname = fullfile(pwd, 'GannetLoad_output', [metabfile_nopath '_' vox{kk} '_load.pdf']);
            end
            saveas(h, pdfname);
            
            % Export the processed data into an SDAT file
            if MRS_struct.p.sdat
                if strcmpi(MRS_struct.p.vendor,'Philips')
                    % Set up filenames
                    sdat_G_name = ['GannetLoad_output/' metabfile_nopath  '_' vox{kk} '_G.sdat'];
                    spar_G_name = ['GannetLoad_output/' metabfile_nopath  '_' vox{kk} '_G.spar'];
                    % Make file copies for SDAT/SPAR files
                    copyfile(metabfile{ii}, sdat_G_name);
                    sparname = [metabfile{ii}(1:end-4) MRS_struct.p.spar_string];
                    copyfile(sparname, spar_G_name);
                    % Write DIFF data into the SDAT file
                    sdat_diff_out = conj(ifft(fftshift(MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{1}).diff(ii,:),2),[],2));
                    sdat_diff_out = sdat_diff_out(1:MRS_struct.p.npoints(ii));
                    % Also write out OFF data
                    sdat_off_out = conj(ifft(fftshift(MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{1}).off(ii,:),2),[],2));
                    sdat_off_out = sdat_off_out(1:MRS_struct.p.npoints(ii));
                    fileid  = fopen(sdat_G_name,'w','ieee-le');
                    ff(:,1:2:2*MRS_struct.p.npoints(ii)) = real(sdat_diff_out);
                    ff(:,2:2:2*MRS_struct.p.npoints(ii)) = imag(sdat_diff_out);
                    gg(:,1:2:2*MRS_struct.p.npoints(ii)) = real(sdat_off_out);
                    gg(:,2:2:2*MRS_struct.p.npoints(ii)) = imag(sdat_off_out);
                    fwriteVAXD(fileid,[ff.' gg.'],'float');
                    fclose(fileid);
                else
                    warning('Only Philips SDAT files can be exported! No data exported.');
                end
            end
            
            % Reorder structure
            if isfield(MRS_struct, 'waterfile')
                structorder = {'version', 'ii', 'metabfile', ...
                    'waterfile', 'p', 'fids', 'spec', 'out'};
            else
                structorder = {'version', 'ii', 'metabfile', ...
                    'p', 'fids', 'spec', 'out'};
            end
            MRS_struct = orderfields(MRS_struct, structorder);
            
            if ii == numscans && MRS_struct.p.mat % save MRS_struct as mat file
                mat_name = fullfile(pwd, ['MRS_struct_' vox{kk} '.mat']);
                if exist(mat_name, 'file')
                    fprintf('\nUpdating results in %s\n', ['MRS_struct_' vox{kk} '.mat...']);
                else
                    fprintf('\nSaving results to %s\n', ['MRS_struct_' vox{kk} '.mat...']);
                end
                save(mat_name, 'MRS_struct');
            end
            
        end % end of output loop over voxels
        
    catch ME
        
        fprintf('\n');
        if ispc
            fname = regexprep(MRS_struct.metabfile{ii}, '\', '\\\');
        else
            fname = MRS_struct.metabfile{ii};
        end
        warning('********** An error occured while loading dataset: ''%s''. Check data. Skipping to next dataset in batch **********', fname);
        errorReport{tryCatchInd} = sprintf(['Filename: ' fname '\n\n' getReport(ME,'extended','hyperlinks','on')]);
        tryCatchInd = tryCatchInd + 1;
        
    end % end of load-and-processing loop over datasets
    
    if ~isempty(errorReport{1}) && ii == numscans
        for ll = 1:size(errorReport,2)
            if ll == 1
                fprintf('\n\n********** ERROR REPORT **********');
            end
            fprintf(['\n\n' errorReport{ll} '\n']);
        end
        fprintf('\n');
    end
    
end

warning('on','stats:nlinfit:ModelConstantWRTParam');
warning('on','stats:nlinfit:IllConditionedJacobian');
warning('on','MATLAB:rankDeficientMatrix');



