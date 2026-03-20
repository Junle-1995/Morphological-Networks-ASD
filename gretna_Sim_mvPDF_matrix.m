function gretna_Sim_mvPDF_matrix(Data_path, Output_path, PDFType, Npoints, SIMType, varargin)

%==========================================================================
% This function is used to calculate similarity among different probability
% distribution functions based on the Jensen-Shannon (JS) divergence or
% Kolmogorov-Smirnov (KS) distance. NB, the probability distribution functions
% are estimated from high-dimensional (>=2) variables.
%
% Syntax: function gretna_Sim_mvPDF_matrix(Data_path, Output_path, PDFType, Npoints, SIMType, varargin)
%
% Inputs:
%         Data_path:
%                   A 1-by-n struct with n denoting the number of variables.
%                   The struct contains two fields with one named 'Path' and
%                   the other named 'File_filter'. The 'Path' field specifies
%                   the directory for each variable where the data files are
%                   sorted, and the 'File_filter' field specifies the prefix
%                   of corresponding data files to be processed.
%       Output_path:
%                   The directory where the resultant files are sorted.
%           PDFType:
%                   'akde';
%                   'mvksdensity'.
%           Npoints:
%                   The number of sampling points (default 2^6).
%           SIMType:
%                   'JS': Jensen¨CShannon distance-based similarity
%                         (default);
%                   'KL': Kullback¨CLeibler divergence-based similarity.
%          varargin:
%                   See gretna_mvPDF_mvksdensity.m or gretna_mvPDF_akde.m
%                   for details.
%
% Example:
%        Suppose we have extracted two regional features of cortical thickness
%        and gyrification, based on which we want to estimate interregional
%        similarity. The 'Data_path' argument can be set as follows:
%
%        Data_path(1).Path        = 'D:\ROI_signal\thickness';
%        Data_path(2).Path        = 'D:\ROI_signal\gyrification';
%        Data_path(1).File_filter = 'CT_signal';
%        Data_path(2).File_filter = 'GI_signal';
%
%        gretna_Sim_mvPDF_matrix(Data_path, 'D:\SIM', 'mvksdensity', 64, 'JS')
%
% Ningkai Wang, HZNU, Hangzhou, 2018/12/07, ningkai.wang.1993@gmail.com
% Jinhui Wang,  HZNU, Hangzhou, 2018/12/07, jinhui.Wang.1982@gmail.com.
%==========================================================================

if nargin < 2
    error('At least two arguments are needed!'); end

if nargin == 2
    PDFType = 'mvksdensity'; Npoints = 2^6; SIMType = 'JS'; end

if nargin == 3
    Npoints = 2^6; SIMType = 'JS'; end

if nargin == 4
    SIMType = 'JS'; end

if nargin > 6
    error('At most five arguments are permitted!'); end

Num_var  = length(Data_path);
Sig_file = cell(Num_var,1);
Npoints  = repmat(Npoints, 1, Num_var./length(Npoints));
Len_Data = prod(Npoints);

for iDim = 1:Num_var
    Sig_file{iDim} = dir([Data_path(iDim).Path filesep Data_path(iDim).File_filter '*.mat']);
end

for isub = 1:size(Sig_file{1},1)
    
    fprintf('Calculating similarity matrix for Subject %s\n', num2str(isub));
    
    MV_Data = cell(Num_var,1);
    
    for iDim = 1:Num_var
        ind_data        = load([Data_path(iDim).Path filesep deblank(Sig_file{iDim}(isub).name)]);
        MV_Data{iDim,1} = ind_data.Sig_roi;
    end
    
    [Num_imgs, Num_regs] = size(ind_data.Sig_roi);
    
    PDF = zeros(Len_Data, Num_regs);
    SIM = zeros(Num_regs, Num_regs, Num_imgs);
    
    for iimg = 1:Num_imgs
        
        for ireg = 1:Num_regs
            
            Data = zeros(length(ind_data.Sig_roi{iimg,ireg}),Num_var);
            
            for iDim = 1:Num_var
                Data(:,iDim) = MV_Data{iDim}{iimg,ireg};
            end
            
            switch lower(PDFType)
                case 'mvksdensity'
                    if isempty(varargin)
                        [PDF_ireg, ~] = gretna_mvPDF_mvksdensity(Data, Npoints);
                    else
                        [PDF_ireg, ~] = gretna_mvPDF_mvksdensity(Data, Npoints, varargin);
                    end
                case 'akde'
                    if isempty(varargin)
                        [PDF_ireg, ~] = gretna_mvPDF_akde(Data, Npoints);
                    else
                        [PDF_ireg, ~] = gretna_mvPDF_akde(Data, Npoints, varargin);
                    end
                otherwise
                    error('Undefined PDFType ''%s''.',PDFType);
            end
            
            PDF(:,ireg)   = PDF_ireg(:);
        end
        
        switch upper(SIMType)
            case 'JS'
                Sim = gretna_JSDs(PDF);
                
            case 'KS'
                Sim = gretna_KSDs(PDF);
                
            otherwise
                error('Undefined SIMType ''%s''.',SIMType)
        end
        
        SIM(:,:,iimg) = Sim;
    end
    
    save([Output_path filesep upper(SIMType) '_' strrep(num2str(Npoints), ' ', '') '_' Sig_file{1}(isub).name], 'SIM');
    
    fprintf('Calculating similarity matrix for Subject %s ...... is done\n', num2str(isub));
end

return