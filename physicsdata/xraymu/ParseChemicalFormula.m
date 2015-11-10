function [symbols,ns] = ParseChemicalFormula(formula)
%  function [symbols,ns] = ParseChemicalFormula(formula)
%  Returns the elements and number of atoms for a chemical formula
%  inputs:
%    formula: a case sensitive string e.g. Na2SO4
%  outputs:
%    symbols: a cell array containing the element symbols e.g. {'Na'},{'S'},{'O'}
%    ns: an array containing the number of atoms or atomic fraction (can be a float)
%      e.g. [2,1,4]
% for more information see http://aprendtech.com/wordpress/?p=45
%  REA 5/8/07

assert(ischar(formula));

  % pre-allocate empty cell array with max number of elements possible
symbols = cell(numel(formula),1);

  % init other variables
nsymbols = 0;
ns = [];
   % regular expression for a float number
numexpr = '(\+|-)?([0-9]+\.?[0-9]*|\.[0-9]+)([eE](\+|-)?[0-9]+)?';

while ~isempty(formula) && (length(formula)>1)
    % look for 2 char element symbol
  s = regexp(formula(1:2),'[A-Z][a-z]','match');
  if ~isempty(s) % found a 2 char symbol so process it
    nsymbols = nsymbols + 1;
    symbols{nsymbols} = s;
    formula = formula(3:end);
  else % did not find 2 element symbol, so next letter in string must be single char symbol
    nsymbols = nsymbols + 1;
    symbols{nsymbols} = formula(1);
    formula = formula(2:end);
  end
  
    % look for number of atoms
  [startpos endpos] = regexp(formula,numexpr);
  if (~isempty(startpos)) && (startpos(1) ==1) % found an atom number
    ns(end+1) = str2num(formula(startpos:endpos));
    formula = formula((endpos+1):end);
  else % assume it's a 1
      ns(end+1) = 1;
  end
end

if ~isempty(formula) % handle case where last symbol is a single symbol element e.g. H2O)
  assert(length(formula)==1);
  nsymbols = nsymbols+1;
  symbols{nsymbols} = formula(1);
  ns(end+1) = 1;
end

symbols = symbols(1:nsymbols); % clip off empty cells
ns = ns(:); % column vector