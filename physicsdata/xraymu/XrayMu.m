function [mus,xray_energies,ztable] = XrayMu(chem_spec,egys,varargin)
% function [mus,xray_energies,ztable] = xraymu(chem_spec,egys,varargin)
%  X-ray attenuation coefficients for a specified compound from tabulated data
%  returns NaN's for energies beyond tabulated data range
%  inputs:
%    chem_spec: (case-sensitive string) chemical formula for the compound e.g. NaSO4 or H2O
%               or an n by 2 matrix, 1 row per element with columns [(atomic number of element), fraction by weight]
%    egys: a vector of x-ray energy (keV) 
%  optional arguments
%    'energy_absorption' if specified the function returns energy absorption coefficents (cm^2/gm)
%    'weight_fraction' if specified expects the formula to specify weight fraction
%         e.g. for water H.111O.888
%  outputs:
%    mus: a column vector of the attenuation coefficients at the specified energies (cm^2/gm)
%    xray_energies: the x-ray energies for the mus values (keV)
%    ztable: a 2 column matrix with atomic numbers in 1st column and fraction by weight in 2nd
%
% for more information see http://aprendtech.com/wordpress/?p=45
%
%  REA 5/28/07-2/10/09 from PhotonAttenuation.m by Jaroslaw Tuszynski
%  REA 7/31/10 add ztable output

nreqargs = 2;
assert(nargin>=nreqargs);

mutype = 'att';
do_weight_fraction = false;

if(nargin>nreqargs)
  i=1;
  while(i<=size(varargin,2))
     switch lower(varargin{i})
     case 'energy_absorption';   mutype = 'abs';
     case 'weight_fraction';     do_weight_fraction = true; 
    otherwise
        error('Unknown argument %s given',varargin{i});
     end
     i=i+1;
  end
end

assert(~isempty(egys) );
assert(all(egys>0.0));

if ischar(chem_spec)
  [elems,ns] = ParseChemicalFormula(chem_spec);  
  [atomic_numbers,atomic_weights] = ChemicalSymbols2AtomicNumbers(elems);
elseif isnumeric(chem_spec) && (size(chem_spec,2)==2) && (all(chem_spec(:,1)<92)) && (all(chem_spec(:,2)<=1)) 
   do_weight_fraction = true;
   atomic_numbers = chem_spec(:,1);
   ns =  chem_spec(:,2);  % this var used for wgt fractions below if do_weight_fraction==true
else
   error('chem_spec format not recognized');
end

   % get the attenuation coefficients for all the elements used 
musall = MusElementsLoc(atomic_numbers,egys,mutype);

   % prepare the weighting
if do_weight_fraction % formula is weight fraction e.g. H.111O.888
  atomic_weightsall = ns(:);
else % formula is atomic fraction e.g. H2O
  atomic_weightsall = ns(:).*atomic_weights(:); 
  atomic_weightsall = (1/sum(atomic_weightsall))*atomic_weightsall;
end
ztable = [atomic_numbers(:) (1/sum(atomic_weightsall))*atomic_weightsall];

mus = sum( musall.* repmat(atomic_weightsall',numel(egys),1), 2); % sum along rows
xray_energies = egys(:);


function musall = MusElementsLoc(atomic_numbers,egys,type)
%  function musall = MusElementsLoc(atomic_numbers,egys)
%  Returns a 2D array of the attenuation coefficients 
%  of the elements with atomic numbers atomic_numbers at energies egys
%  inputs:
%    atomic_numbers: atomic numbers of the elements 1:93
%    egys: the x-ray energies (keV)
%    type: (string) 'att' (mass attenuation coeffs), 'abs' energy absorption
%  outputs:
%    musall: a matrix with each column the mass attentuation coefficients for the 
%            corresponding element in atomic_numbers at the specified energies
%  REA 5/28/07

    % get the tables of data
[weights,symbols,xdatas] = XrayData;

  % check inputs
assert(all(floor(atomic_numbers(:)) == atomic_numbers(:)));    % make sure they are integers      
mn = min(atomic_numbers(:)); mx = max(atomic_numbers(:));     % and in proper range
assert( (mn>0) && (mn<=length(xdatas) ) );
assert( (mx>0) && (mx<=length(xdatas)) );
assert(all(egys>0));
egys = 0.001*egys; % convert to MeV

    % prep log-log interpolation from NIST data
musall = zeros(numel(egys),numel(atomic_numbers));
for kz = 1:numel(atomic_numbers)
  z = atomic_numbers(kz);
  xdata = xdatas{z};
  etable = log(xdata.PhotonEnergy);
  switch type
    case 'att'      
      mus = log( xdata.MAC' );
    case 'abs'
      mus = log( xdata.MEAC' );
     otherwise
      error('MusElementsLoc: absorption type %s unknown',type)
  end
         % in NIST table, energies are duplicated at absorption edges
         % but this gives interp1 problems, so nudge them just below and above edge 
  d = 5*eps;            % offset energies this much at absorption edges (MeV)
  k = find(diff(etable)==0);% find edges
  etable(k) = etable(k)-d;   % give them some width
  etable(k+1) = etable(k+1)+d;
  a = interp1(etable, mus, log(egys), 'PCHIP'); % cubic Hermite polynomial
  a(a>9.267) = 9.267;  % extrapolated values can get high especially close to edges. Clip them.
  musall(:,kz) = exp(a(:));
end

