function [PDF, Grid] = gretna_mvPDF_mvksdensity(Data, Npoints, varargin)

%==========================================================================
% This function is used to calculate the probability distribution function
% of high-dimensional (>=2) data based on an optimized kernel density estimate.
%
%
% Syntax: function [PDF, Grid] = gretna_mvPDF_mvksdensity(Data, Npoints, varargin)
%
% Inputs:
%        Data:
%             An n-by-d matrix with n denoting the number of data points
%             and d denoting the number of dimensions.
%     Npoints:
%             A scalar or 1-by-d vector for the size of resultant grid. E.g.,
%             for a 3-D data, the output pdf will be a 64-64-64 matrix if
%             Npoints = 64, and a 60-80-100 matrix if Npoints = [60 80 100].
%             For huge data, 64 is recommended.
%    varargin:
%             By calling MATLAB build-in function 'mvksdensity', this
%             argument specifies parameter name/value pairs to control the
%             density estimation. Valid name/value pairs are as follows.
%             'Kernel' -- Type of kernel smoother:
%                      'normal' (default);
%                      'box';
%                      'triangle';
%                      'epanechnikov'.
%             'Support' -- Support for the density:
%                      'unbounded' (default) allows the density to extend
%                          over the whole real line;
%                      'positive' restricts the density to positive values;
%                      2-by-d matrix specifies the finite lower (first row)
%                          and upper (second row) bounds for the support of
%                          the density. Each column contains the limits for
%                          one dimension of Data. NB. For data such as volume
%                          and thickness, 'positive' is recommended.
%             'Weights':
%                      A n-by-1 vector assigning weights to each data point
%                      in Data (default equal weight).
%             'BandWidth':
%                      A scalar value or 1-by-d vector assigning the
%                      bandwidth of the kernel smoothing window for each
%                      dimension in Data. The default bandwidth is
%                      calculated based on Scott's rule.
%             'Function' -- Function to estimate:
%                      'pdf' (default) for probability density function;
%                      'cdf' for cumulative distribution function;
%                      'survivor' for survivor function.
% Outputs:
%         PDF:
%             The resulant probability distribution function.
%        Grid:
%             Optimized points at which the PDF are estimated.
%
% Examples:
%        For a 3-D Data,
%        1. [PDF, Grid] = gretna_mvPDF(Data, 64, 'Kernel', 'normal', 'Function', 'pdf')
%        returns a 64-64-64 PDF matrix with a normal kernel smoother.
%        2. [PDF, Grid] = gretna_mvPDF(Data, [60, 80, 100]);
%        returns a 60-80-100 PDF matrix estimated with a customized grid
%        size.
%        3. [PDF, Grid] = gretna_mvPDF(Data, 64, 'BandWidth', [.3 .4 .6])
%        returns a 64-64-64 PDF matrix with a customized bandWidth for each
%        dimension.
%
% Reference:
%   Scott, D. W., & Sain, S. R. (2005). Multidimensional density estimation.
%   Handbook of statistics, 24, 229-261.
%
% Ningkai WANG,HZNU, Hangzhou, 2018/10/20, ningkai.wang.1993@gmail.com
% Jinhui WANG, SCNU, Guangzhou,2018/10/20, jinhui.wang.1982@gmail.com
%==========================================================================

matlabVersionDetection = isempty(which('mvksdensity'));
if matlabVersionDetection
    error([' This function depends on MATLAB built-in function %s,'...
        ' which is released in the Statistics and Machine Learning Toolbox'...
        ' of MATLAB R2016a and later.'],'''mvksdensity''');
end

[n, d]  = size(Data);
Npoints = repmat(Npoints, 1, d./length(Npoints));

if sum(strcmpi(varargin, 'Support'))
    supportType = varargin{find(strcmpi(varargin, 'Support'))+1};
else
    supportType = [];
end

if d == 1
    error('Please use gretna_PDF.m for 1-D data!');

elseif d >= 2 && rem(d, 1) == 0
    Sig  = median(abs(Data - repmat(median(Data), n, 1)))/0.6745;
    % The estimation of the standard deviation

    Grid = makeGrid(Data, Npoints, d, n, Sig, supportType);

    if d == 2
        [f, Grid] = ksdensity(Data, Grid, varargin{:});

    elseif d >= 3
        whereBW = strcmpi(varargin, 'BandWidth');
        isBW    = sum(whereBW);

        if ~isBW
            BW = Sig.*(n.^(-1./(d+4)));
            % The multivariate rule-of-thumb for the bandwidth of the
            % kernel smoothing window (normal product kernel).
            
            f  = mvksdensity(Data, Grid, 'BandWidth', BW, varargin{:});
        else
            f  = mvksdensity(Data, Grid, varargin{:});
        end
        
    end

    f(f <= 0) = eps;
    f         = reshape(f, ones(1, d).* Npoints);
    PDF       = f./sum(f(:));
end

end

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
function Grid = makeGrid(Data, Npoints, d, n, Sig, supportType)
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

if min(Sig) <= 0
    Sig(Sig<=0) = max(Data, [], 1) - min(Data, [], 1);
end

if min(Sig) > 0
    Data_u = Sig .* (4/(3*n))^0.2;
else
    Data_u = ones(1, d);
end

if isempty(supportType) || strcmpi(supportType, 'unbounded')
    MIN = min(Data, [], 1) - 3*Data_u;
    MAX = max(Data, [], 1) + 3*Data_u;

elseif ~isempty(supportType)
    if ~ischar(supportType)
        MIN = supportType(2,:) - 3*Data_u;
        MAX = max(Data, [], 1) + 3*Data_u;

    elseif ischar(supportType)
        if strcmpi(supportType, 'positive')
            MIN = zeros(1, d);
            MAX = max(Data, [], 1) + 3*Data_u;

        else
            error('Unidentified ''Support'' type.');

        end
    end
end

scaling = MAX - MIN;
Grid    = computeNdgridMatrix(Npoints, MIN, scaling, d);

end


%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
function Grid = computeNdgridMatrix(Npoints, MIN, scaling, d)
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

inputs  = cell(1, d);
outputs = cell(1, d);
scaling = scaling./(Npoints-1);

for iDim = 1:d
    add_Scaling  = bsxfun(@times, scaling(:, iDim), (0:Npoints(iDim)-1)');
    A            = repmat(MIN(iDim), [Npoints(iDim), 1]) + add_Scaling;
    inputs(iDim) = mat2cell(A, Npoints(iDim), ones(1,1));
end

[outputs{:}] = ndgrid(inputs{:});

for iDim = 1:d
    outputs{iDim} = outputs{iDim}(:);
end

Grid = cell2mat(outputs);

end