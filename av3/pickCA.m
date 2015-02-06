function crd = pickCA(pdbFileName)

% Read CA protein data from pdb text file.

fid = fopen(pdbFileName);   % Open coordinate file
Eline = fgetl(fid);         % Skip till you get to ATOM lines
word = words(Eline);        % word-break the line to words separated by space

% ~ is a logical not. Read next line if ~ TER and ~ ATOM
while (and(~strcmp(word{1},'TER'), ~strcmp(word{1},'ATOM'))) 
   Fline = fgetl(fid);   % The first atom is N
   word = words(Fline);
end

% if TER end of chain
if (strcmp(word{1},'TER')) break; end;

Fline = fgetl(fid);     % The second atom is CA
word = words(Fline);

%special protein data bank cases are handled below
if or(strcmp(word{5},'R'), strcmp(word{5},'A'))
   x=7; y=8; z=9;
else
   x=6; y=7; z=8;
end

i=0;
while (~strcmp(word{1},'TER'))
   i = i+1;
   coords(i,:) =  [str2double(word{x}),str2double(word{y}),str2double(word{z})];
   while (~strcmp(word{1},'TER'))
      Fline = fgetl(fid);    % Find the CA of the next amino acid
      word = words(Fline);
      if strcmp(word{1},'TER') break; end;
      if ~strcmp(word{1},'ATOM') break; end;
      if strcmp(word{3},'CA') break; end;
   end
   if or(strcmp(word{1},'TER'),~strcmp(word{1},'ATOM')) break; end;
end
fclose(fid);

crd = coords(1:i,:);


function word = words(string)

% Break a string into whitespace delimited words.
% Produces a cell-array of strings.  Individual strings are 
% accessed via s{1}, s{2}, and s{3}.  Note the curly braces.

s = size(string,2);                    % Length of input string
blanks = find(isspace(string));        % Find indices of all whitespace
zblanks = [0 blanks];                  % 0 counts as whitespace, too
n = size(zblanks,2);                   % Number of words (includes length 0)
lengths = [blanks s+1] - zblanks - 1;  % Lengths of all words
result = cell(1,n);                    % Reserve space
for i = 1:n                            % Copy all words (includes length 0)
   result{1,i} = string(zblanks(i)+1:zblanks(i)+lengths(i));
end
word = result(lengths > 0);     % Remove length 0 words