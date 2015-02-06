function tom_pdbwrite(pdb_name,pdb_st)
%PDBWRITE writes pdb file
% 
%   PDBWRITE(pdb_name, pdb_st) writes the contents of pdb_st to
%   
%   Example:
%
%    pdb=tom_pdbread2('test.pdb');   
%    tom_pdbwrite()   
%
%   See also tom_pdbread2
%
%   Copyright 2006-2007 The MathWorks, Inc.
%   $Revision: 1.1.6.11 $  $Date: 2010/12/22 16:17:38 $

%   Reference:
%   The PDB format contents guide version 2.3.
%   http://www.wwpdb.org/documentation/format23/v2.3.html


% Validate input data
if nargin < 2
    error('no enough input vals');
end

% Check pdb_name opens
fid = fopen(pdb_name, 'wt');

if fid == (-1)
    [theDir, theFile, theExtension] = fileparts(pdb_name);
    if ~isempty(theDir)
        error(message('bioinfo:pdbwrite:CouldNotOpenFileinDir', [ theFile, theExtension ], theDir));
    else
        error(message('bioinfo:pdbwrite:CouldNotOpenFileinCurrentDir', pdb_name));
    end
end

% Check that the record has all the unconditional mandatory fields,
% otherwise throw a warning.
% MandatoryFields = {'Header',...
%                'Compound',...
%                'Source',...
%                'Authors',...
%                'RevisionDate',...
%                'Remark2',...
%                'Remark3',...
%                'Cryst1',...
%                'OriginX',...
%                'Scale',...
%                'Master'};
MandatoryFields = {''};
           
m = isfield(pdb_st, MandatoryFields); 

% if ~all(m)
%     missingFields = MandatoryFields(~m);
%     if numel(missingFields)== numel(MandatoryFields)
%         warning(message('bioinfo:pdbwrite:MissingAllMandatoryFields')); 
%     else
%         mf = ''; % list with mandatory fields
%         for i=1:numel(missingFields)
%             mf = [mf missingFields{i} ', ']; %#ok<AGROW>
%         end
%         mf = mf(1:end-2); % remove extra comma and space from the end
%         warning(message('bioinfo:pdbwrite:MissingMandatoryField', mf));
%         
%     end
% end

% Writes pdbrecords
try
    pdbrecord = getPDBRecord(pdb_st);
    
    numRecords = size(pdbrecord, 1);
    for i=1:numRecords
        fprintf(fid,'%s\n', pdbrecord(i,:));
    end
    
    fclose(fid);
    
catch le
    % close and delete the file if there was an error
    fclose(fid);
    delete(pdb_name);
    %rethrow the exception
    rethrow(le);
end

%-------------------------------------------------------------------------
function pdbout = getPDBRecord(pdb_struct)
% Returns a character array in which each row corresponds to a line in a
% PDB record.

pdbout = [];
%write HEADER
if(isfield(pdb_struct,'Header'))
    pdbout = [pdbout;writeHeader(pdb_struct.Header)];
end
%write OBSLTE
if(isfield(pdb_struct,'Obsolete'))
    pdbout = [pdbout;writeOBS_SPR(pdb_struct.Obsolete,'OBSLTE')];
end
%write TITLE
if(isfield(pdb_struct,'Title'))
    pdbout = [pdbout;writeSingleRecMultipleLine(pdb_struct.Title,'TITLE')];
end
%write CAVEAT
if(isfield(pdb_struct,'Caveat'))
    pdbout = [pdbout;writeCaveat(pdb_struct.Caveat)];
end
%write COMPND
if(isfield(pdb_struct,'Compound'))
    pdbout = [pdbout;writeSingleRecMultipleLine(pdb_struct.Compound,'COMPND')];
end
%write SOURCE
if(isfield(pdb_struct,'Source'))
    pdbout = [pdbout;writeSingleRecMultipleLine(pdb_struct.Source,'SOURCE')];
end
%write KEYWDS
if(isfield(pdb_struct,'Keywords'))
    pdbout = [pdbout;writeSingleRecMultipleLine(pdb_struct.Keywords,'KEYWDS')];
end
%write EXPDTA
if(isfield(pdb_struct,'ExperimentData'))
    pdbout = [pdbout;writeSingleRecMultipleLine(pdb_struct.ExperimentData,'EXPDTA')];
end
%write AUTHOR
if(isfield(pdb_struct,'Authors'))
    pdbout = [pdbout;writeSingleRecMultipleLine(pdb_struct.Authors,'AUTHOR')];
end
%write REVDAT
if(isfield(pdb_struct,'RevisionDate'))
    pdbout = [pdbout;writeREVDAT(pdb_struct.RevisionDate)];
end
%write SPRSDE
if(isfield(pdb_struct,'Superseded'))
    pdbout = [pdbout;writeOBS_SPR(pdb_struct.Superseded,'SPRSDE')];
end
%write JRNL
if(isfield(pdb_struct,'Journal'))
    pdbout = [pdbout;writeJournal(pdb_struct.Journal)];
end
%write REMARK
remarks = regexp(fieldnames(pdb_struct),'Remark[0-9]*','match');
remarks =[remarks{:}];
for i = 1:length(remarks)
    remarkNo = sscanf(remarks{i},'Remark%d');
    switch remarkNo
        case 1
            if(isfield(pdb_struct.Remark1,'JournalEntry'))
                pdbout = [pdbout;writeRemark1(pdb_struct.Remark1.JournalEntry)]; %#ok
            end
        case 2
            pdbout = [pdbout;writeRemark2(pdb_struct.Remark2)];%#ok
        case 3
            pdbout = [pdbout;writeRemark3(pdb_struct.Remark3)];%#ok
        otherwise
            if(~isempty(remarkNo))
                pdbout = [pdbout;writeRemarkN(pdb_struct.(remarks{i}),remarkNo)];%#ok
            end;
    end
end

%write DBREF
if(isfield(pdb_struct,'DBReferences'))
    pdbout = [pdbout;writeDBREF(pdb_struct.DBReferences)];
end
%write SEQADV
if(isfield(pdb_struct,'SequenceConflicts'))
    pdbout = [pdbout;writeSeqAdv(pdb_struct.SequenceConflicts)];
end
%write SEQRES
if(isfield(pdb_struct,'Sequence'))
    pdbout = [pdbout;writeSEQRES(pdb_struct.Sequence)];
end

%write FTNOTE
footnotes = regexp(fieldnames(pdb_struct),'Footnote[0-9]*','match');
footnotes =[footnotes{:}];
for i = 1:length(footnotes)
    footnoteNo = sscanf(footnotes{i},'Footnote%d');
    pdbout = [pdbout;writeFootnoteN(pdb_struct.(footnotes{i}),footnoteNo)];%#ok
end

%write MODRES
if(isfield(pdb_struct,'ModifiedResidues'))
    pdbout = [pdbout;writeMODRES(pdb_struct.ModifiedResidues)];
end
%write HET
if(isfield(pdb_struct,'Heterogen'))
    pdbout = [pdbout;writeHET(pdb_struct.Heterogen)];
end
%write HETNAM
if(isfield(pdb_struct,'HeterogenName'))
    pdbout = [pdbout;writeHETNAM_SYN(pdb_struct.HeterogenName,'HETNAM')];
end
%write HETSYN
if(isfield(pdb_struct,'HeterogenSynonym'))
    pdbout = [pdbout;writeHETNAM_SYN(pdb_struct.HeterogenSynonym,'HETSYN')];
end
%write FORMUL
if(isfield(pdb_struct,'Formula'))
    pdbout = [pdbout;writeFormula(pdb_struct.Formula)];
end
%write HELIX
if(isfield(pdb_struct,'Helix'))
    pdbout = [pdbout;writeHelix(pdb_struct.Helix)];
end
%write SHEET
if(isfield(pdb_struct,'Sheet'))
    pdbout = [pdbout;writeSheet(pdb_struct.Sheet)];
end
%write TURN
if(isfield(pdb_struct,'Turn'))
    pdbout = [pdbout;writeTurn(pdb_struct.Turn)];
end
%write SSBOND
if(isfield(pdb_struct,'SSBond'))
    pdbout = [pdbout;writeSSbond(pdb_struct.SSBond)];
end
%write LINK
if(isfield(pdb_struct,'Link'))
    pdbout = [pdbout;writeLINK(pdb_struct.Link)];
end
%write HYDBND
if(isfield(pdb_struct,'HydrogenBond'))
    pdbout = [pdbout;writeHYDBND(pdb_struct.HydrogenBond)];
end
%write SLTBRG
if(isfield(pdb_struct,'SaltBridge'))
    pdbout = [pdbout;writeSaltBridge(pdb_struct.SaltBridge)];
end
%write CISPEP
if(isfield(pdb_struct,'CISPeptides'))
    pdbout = [pdbout;writeCISPEP(pdb_struct.CISPeptides)];
end
%write SITE
if(isfield(pdb_struct,'Site'))
    pdbout = [pdbout;writeSite(pdb_struct.Site.SiteDetail)];
end
%write CRYST1
if(isfield(pdb_struct,'Cryst1'))
    pdbout = [pdbout;writeCRYST1(pdb_struct.Cryst1)];
end
%write ORIGX
if(isfield(pdb_struct,'OriginX'))
    pdbout = [pdbout;writeOriginX(pdb_struct.OriginX)];
end
%write SCALE
if(isfield(pdb_struct,'Scale'))
    pdbout = [pdbout;writeScale(pdb_struct.Scale)];
end
%write MTRIX
if(isfield(pdb_struct,'Matrix'))
    pdbout = [pdbout;writeMTRIX(pdb_struct.Matrix)];
end
%write TVECT
if(isfield(pdb_struct,'TranslationVector'))
    pdbout = [pdbout;writeTVECT(pdb_struct.TranslationVector)];
end

%write MODEL and ATOM records
if(isfield(pdb_struct,'Model'))
    if numel(pdb_struct.Model) == 1 && ~isfield(pdb_struct.Model, 'MDLSerNo')
        pdbout = [pdbout;writeModel(pdb_struct.Model)];
    else
        modelLine = blanks(80);
        endModelLine = blanks(80);
        for i = 1:numel(pdb_struct.Model)
            if isfield(pdb_struct.Model, 'MDLSerNo');
                modelSerNo = pdb_struct.Model(i).MDLSerNo;
            else
                modelSerNo = i;
            end
            out = writeModel(pdb_struct.Model(i));
            modelLine(1:14) = ['MODEL     ',...
                sprintf('%4d',modelSerNo)];
            endModelLine(1:6) = 'ENDMDL';
            pdbout = [pdbout; modelLine;out;endModelLine]; %#ok
        end
    end    
end

%write CONECT
if(isfield(pdb_struct,'Connectivity'))
    pdbout = [pdbout;writeConnect(pdb_struct.Connectivity)];
end

%write MASTER and END
if(isfield(pdb_struct,'Master'))
    pdbout = [pdbout;writeMaster(pdb_struct.Master)];
end

%----------------------------------------------------------%
function z = writeSEQRES(m_seqres)
% Number of chains
no_of_chains = length(m_seqres);
out = sum(arrayfun(@(x) ceil(length(x.ResidueNames)/52), m_seqres));
z = repmat(blanks(80), out, 1);
lc = 1;
y = [];
for i=1:no_of_chains
    % Parse residues so there are at most 11 residues on a line
    full_seq = m_seqres(i).ResidueNames;
    len_of_chain = length(full_seq);
    no_of_lines = floor(len_of_chain/52);
    r = mod(len_of_chain,52);
    if(no_of_lines)
        seqs = reshape(full_seq(1:(len_of_chain-r)),52,no_of_lines)';
        seqs = seqs(:,1:51);
        y = repmat(blanks(80),no_of_lines,1);
        y(:,1:70)= [reshape(sprintf('SEQRES  %2d ',1:no_of_lines),11,no_of_lines)' ...
            repmat(sprintf('%c %4d  ',...
            m_seqres(i).ChainID,...
            m_seqres(i).NumOfResidues),...
            no_of_lines,1) seqs];
    end
    if(r)
        % Write remaining part of the sequence
        x = [sprintf('SEQRES  %2d %c %4d  ',no_of_lines+1,m_seqres(i).ChainID,...
            m_seqres(i).NumOfResidues) full_seq(len_of_chain-r+1:len_of_chain)];
        z(lc:lc+no_of_lines,:)=[y;[x blanks(51-r+10)]];
        lc = lc + no_of_lines + 1;
        y=[];
    end
end

%-----------------------------------------------------------------%
function z = writeCISPEP(mRec)
noOfLines = length(mRec);
z = repmat(blanks(80),noOfLines,1);
for i=1:noOfLines
    z(i,1:59) = ['CISPEP ' ...
        sprintf('%3d %3s %1s %4d%1s   %3s %1s %4d%1s       %3d       %6.2f',...
        mRec(i).serNum,...
        mRec(i).ResName1,...
        mRec(i).chainID1,...
        mRec(i).seqNum1,...
        mRec(i).icode1,...
        mRec(i).ResName2,...
        mRec(i).chainID2,...
        mRec(i).seqNum2,...
        mRec(i).icode2,...
        mRec(i).modNum,...
        mRec(i).measure)];
end

%----------------------------------------------------------%
function z = writeCRYST1(mRec)
z = blanks(80);
z(1:70) = sprintf('CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %11s%4d',...
    mRec.a,...
    mRec.b,...
    mRec.c,...
    mRec.alpha,...
    mRec.beta,...
    mRec.gamma,...
    mRec.sGroup,...
    mRec.z);

%---------------------------------------------------------------------%
function z = writeCaveat(mRec)
%for CAVEAT
[noOfLines lineLen] = size(mRec.comment);
z = repmat(blanks(80),noOfLines,1);
z(:,1:6) = repmat('CAVEAT',noOfLines,1);
%first line
z(1,12:lineLen+19) = [mRec.idCode '    ' mRec.comment(1,:)];
if(noOfLines > 1)
    z(2:end,9:10) = reshape(sprintf('%2d',2:noOfLines),2,noOfLines-1)';
    z(2:end,12:lineLen+19) = [ repmat(mRec.idCode,noOfLines-1,1) '    '...
     mRec.comment(2:noOfLines,:)];
end

%---------------------------------------------------------------------%
function z=writeConnect(mRec)
n = length(mRec);
z = repmat(blanks(80),n,1);
for i =1:n
    y = blanks(80);
    temp = ['CONECT' sprintf('%5d',mRec(i).AtomSerNo) sprintf('%5d', ...
        mRec(i).BondAtomList)];
    tLen = length(temp);
    y(1,1:tLen) = temp;
    noH = length(mRec(i).HydrogenAtomList);
    noS = length(mRec(i).SaltBridgeAtom);
    x = blanks(30);
    for j = 1:noH
        idx = 1;
        x(idx:idx+4) = sprintf('%5d',mRec(i).HydrogenAtomList(j));
        idx = idx + 5;
        if(~mod(j,2))
            idx = idx + 5; %#ok<NASGU>
        end
    end

    for j = 1:noS
        switch j
            case 1
                x(11:15) = sprintf('%5d',mRec(i).SaltBridgeAtom(j));
            case 2
                x(26:30) = sprintf('%5d',mRec(i).SaltBridgeAtom(j));
        end
    end
    y(32:61) = x;
    z(i,1:length(y)) = y;
end


function z = writeDBREF(mRec)
noOfLines = length(mRec);
z = repmat(blanks(80),noOfLines,1);
for i=1:noOfLines
    z(i,1:68) = ['DBREF  ' ...
        sprintf('%4s %1s %4d%1s %4d%1s %6s %8s %12s %5d%1s %5d%1s',...
        mRec(i).idCode,...
        mRec(i).chainID,...
        mRec(i).seqBegin,...
        mRec(i).insertBegin,...
        mRec(i).seqEnd,...
        mRec(i).insertEnd,...
        mRec(i).database,...
        mRec(i).dbAccession,...
        mRec(i).dbIdCode,...
        mRec(i).dbseqBegin,...
        mRec(i).idbnsBeg,...
        mRec(i).dbseqEnd,...
        mRec(i).dbinsEnd)];
end

%----------------------------------------------------------------------%
function z = writeFootnoteN(mRec,remarkNo)
%for FTNOTE
[m n] = size(mRec);
z = repmat(blanks(80),m,1);
z(:,1:11+n) = [repmat(sprintf('FTNOTE %3d ',remarkNo),m,1) mRec];

%--------------------------------------------------------------------------
function z = writeFormula(mRec)
% FORMULA
out = lineStats(mRec,'ChemForm');
noOfFormula = size(out,1);
z = repmat(blanks(80),sum(out(:,1)),1);
lc = 1;
for i=1:noOfFormula
    noOfLines = out(i,1);
    lineLen = out(i,2);
    y = repmat(blanks(80),noOfLines,1);
    y(:,1:16) = repmat(['FORMUL  ',...
        sprintf('%2d  %3s ',...
        mRec(i).CompNo,...
        mRec(i).hetID)],...
        noOfLines,1);
    %first line
    %check if hetID is 'HOH' or 'UNX'
    if(strcmp(mRec(i).hetID,'HOH') || strcmp(mRec(i).hetID,'UNX'))
        y(1,19:18+lineLen) = mRec(i).ChemForm(1,:);
        if noOfLines > 1
            y(2:noOfLines,17:18+lineLen) = [reshape(sprintf('%2d',2:noOfLines), 2,noOfLines-1)',...
                        mRec(i).ChemForm(2:noOfLines,:)];
        end
    else
        %shifted one char to take care of the missing asterisk '*'
        y(1,20:19+lineLen) = mRec(i).ChemForm(1,:);
        if noOfLines > 1
            y(2:noOfLines,17:19+lineLen) = [reshape(sprintf('%2d',2:noOfLines),2,noOfLines-1)',...
                repmat(' ',noOfLines-1,1),...
                mRec(i).ChemForm(2:noOfLines,:)];
        end
    end
    z(lc:lc+noOfLines-1,:) = y;
    lc = lc + noOfLines;
end

%------------------------------------------------------------------%
function z = writeHET(mRec)
noOfHet = length(mRec);
z = repmat(blanks(80),noOfHet,1);
for i=1:noOfHet
    z(i,1:30+length(mRec(i).text)) = ['HET    ' ...
        sprintf('%3s  %1s%4d%1s  %5d     %s',...
        mRec(i).hetID,...
        mRec(i).ChainID,...
        mRec(i).seqNum,...
        mRec(i).iCode,...
        mRec(i).numHetAtoms,...
        mRec(i).text)];
end

%--------------------------------------------------------------------%
function z = writeHETNAM_SYN(mRec,recName)
%HETNAM, HETSYN

if(strcmp(recName,'HETNAM'))
	out = lineStats(mRec,'ChemName');
    noOfHET = size(out,1);
    z = repmat(blanks(80),sum(out(:,1)),1);
    lc = 1;
    for i=1:noOfHET
        noOfLines = out(i,1);
        lineLen = out(i,2);
        y = repmat(blanks(80),noOfLines,1);
        y(:,1:6) = repmat(recName,noOfLines,1);

        %first line
        y(1,12:15+lineLen) = [mRec(i).hetID ' ' mRec(i).ChemName(1,:)];
        if noOfLines > 1
            y(2:noOfLines,9:15+lineLen) = [reshape(sprintf('%2d',2:noOfLines),2,noOfLines-1)' ...
                repmat(sprintf(' %3s ',mRec(i).hetID),noOfLines-1,1),...
                mRec(i).ChemName(2:noOfLines,:)];
        end
        z(lc:lc+noOfLines-1,:) = y;
    end
elseif(strcmp(recName,'HETSYN'))
    out = lineStats(mRec,'hetSynonyms');
    noOfHET = size(out,1);
    z = repmat(blanks(80),sum(out(:,1)),1);
    lc = 1;
    for i=1:noOfHET
        noOfLines = out(i,1);
        lineLen = out(i,2);
        y = repmat(blanks(80),noOfLines,1);
        y(:,1:6) = repmat(recName,noOfLines,1);
        % First line
        y(1,12:15+lineLen) = [mRec(i).hetID ' ' mRec(i).hetSynonyms(1,:)];
        if noOfLines > 1
            y(2:noOfLines,9:15+lineLen) = [reshape(sprintf('%2d',2:noOfLines),2,noOfLines-1)' ...
                repmat(sprintf(' %3s ',mRec(i).hetID),noOfLines-1,1),...
                mRec(i).hetSynonyms(2:noOfLines,:)];
        end
        z(lc:lc+noOfLines-1,:) = y;
    end
end

function z = writeHYDBND(mRec)
noOfLines = length(mRec);
z = repmat(blanks(80),noOfLines,1);
for i=1:noOfLines
    z(i,1:72) = ['HYDBND      ' ...
        sprintf('%4s%1s%3s %1s%5d%1s %4s%1s %1s%5d%1s %4s%1s%3s %1s%5d%1s%6s %6s',...
        mRec(i).AtomName1,...
        mRec(i).altLoc1,...
        mRec(i).resName1,...
        mRec(i).Chain1,...
        mRec(i).resSeq1,...
        mRec(i).ICode1,...
        mRec(i).nameH,...
        mRec(i).altLocH,...
        mRec(i).ChainH,...
        mRec(i).resSeqH,...
        mRec(i).iCodeH,...
        mRec(i).name2,...
        mRec(i).altLoc2,...
        mRec(i).resName2,...
        mRec(i).chainID2,...
        mRec(i).resSeq2,...
        mRec(i).iCode2,...
        mRec(i).sym1,...
        mRec(i).sym2)];
end

%-----------------------------------------------------%
function z = writeHeader(mRec)
z = blanks(80);
z(1:66) = sprintf('HEADER    %-40s%9s   %4s',...
    mRec.classification,...
    mRec.depDate,...
    mRec.idCode);

%-------------------------------------------------------%
function z = writeHelix(mRec)
noOfLines = length(mRec);
z = repmat(blanks(80),noOfLines,1);
for i=1:noOfLines
    z(i,1:76) = ['HELIX  ' ...
        sprintf('%3d %3s %3s %1s %4d%1s %3s %1s %4d%1s%2d%30s %5d',...
        mRec(i).serNum,...
        mRec(i).helixID,...
        mRec(i).initResName,...
        mRec(i).initChainID,...
        mRec(i).initSeqNum,...
        mRec(i).initICode,...
        mRec(i).endResName,...
        mRec(i).endChainID,...
        mRec(i).endSeqNum,...
        mRec(i).endICode,...
        mRec(i).helixClass,...
        mRec(i).comment,...
        mRec(i).length)];
end

%--------------------------------------------------------%
function z = writeJournal(mJournal)
noOfJournalEntries = size(mJournal,1);
out = sum(arrayfun(@(x)  size(x.Author,1) + size(x.Title,1) + ...
            size(x.Editor,1) + size(x.Reference,1) + size(x.Publisher,1)+1, ...
            mJournal));
z = repmat(blanks(80),out,1);
lc = 1;
for i=1:noOfJournalEntries
    [m n] = size(mJournal(i).Author);
    y = repmat(blanks(80),m,1);
    y(:,1:n+3+16) = [repmat('JRNL        AUTH',m,1),...
        reshape(sprintf('%2d ',1:m),3,m)'  mJournal(i).Author];
    y(1,17:18) = '  ';%remove continuation number for the first author line
     z(lc:lc+m-1,:) = y;
    lc = lc + m;
    %Optional
    if(~isempty(mJournal(i).Title))
        [m n] = size(mJournal(i).Title);
        y = repmat(blanks(80),m,1);
        y(:,1:n+3+16) = [repmat('JRNL        TITL',m,1),...
            reshape(sprintf('%2d ',1:m),3,m)'  mJournal(i).Title];
        y(1,17:18) = '  ';%remove continuation number for the first Title line
        z(lc:lc+m-1,:) = y;
        lc = lc + m;
    end
    %Optional
    if(~isempty(mJournal(i).Editor))
        [m n] = size(mJournal(i).Editor);
        y = repmat(blanks(80),m,1);
        y(:,1:n+3+16) = [repmat('JRNL        EDIT',m,1),...
            reshape(sprintf('%2d ',1:m),3,m)'  mJournal(i).Editor];
        y(1,17:18) = '  ';%remove continuation number for the first Edit line
        z(lc:lc+m-1,:) = y;
        lc = lc + m;
    end
    [m n] = size(mJournal(i).Reference);
    y = repmat(blanks(80),m,1);
    y(:,1:n+3+16) = [repmat('JRNL        REF ',m,1),...
        reshape(sprintf('%2d ',1:m),3,m)'  mJournal(i).Reference];
    y(1,17:18) = '  ';
    z(lc:lc+m-1,:) = y;
    lc = lc + m;
    %Optional
    if(~isempty(mJournal(i).Publisher))
        [m n] = size(mJournal(i).Publisher);
        y = repmat(blanks(80),m,1);
        y(:,1:n+3+16) = [repmat('JRNL        PUBL',m,1),...
            reshape(sprintf('%2d ',1:m),3,m)'  mJournal(i).Publisher];
        y(1,17:18) = '  ';%remove continuation number for the first Edit line
        z(lc:lc+m-1,:) = y;
        lc = lc + m;
    end
    y = blanks(80);
    y(1,1:16) = 'JRNL        REFN';
    y(1,20:70) = mJournal(i).CitationReference;
    z(lc,:) = y;
    lc = lc + 1;
end


function z = writeLINK(mRec)
noOfLines = length(mRec);
z = repmat(blanks(80),noOfLines,1);
for i=1:noOfLines
    z(i,1:72) = ['LINK        ' ...
        sprintf('%4s%1s%3s %1s%4d%1s               %4s%1s%3s %1s%4d%1s  %6s %6s',...
        mRec(i).remove1,...
        mRec(i).altLoc1,...
        mRec(i).resName1,...
        mRec(i).chainID1,...
        mRec(i).resSeq1,...
        mRec(i).iCode1,...
        mRec(i).AtomName2,...
        mRec(i).altLoc2,...
        mRec(i).resName2,...
        mRec(i).chainID2,...
        mRec(i).resSeq2,...
        mRec(i).iCode2,...
        mRec(i).sym1,...
        mRec(i).sym2)];
end

%-----------------------------------------------------------%
function z = writeMODRES(mRec)
noOfLines = length(mRec);
z = repmat(blanks(80),noOfLines,1);
for i=1:noOfLines
    n = length(mRec(i).comment);
    z(i,1:29+n) = ['MODRES ' ...
        sprintf('%4s %3s %1s %4d%1s %3s  %s',...
        mRec(i).idCode,...
        mRec(i).resName,...
        mRec(i).chainID,...
        mRec(i).seqNum,...
        mRec(i).iCode,...
        mRec(i).stdRes,...
        mRec(i).comment)];
end

%----------------------------------------------------------%
function z= writeMTRIX(mRec)
n = length(mRec);
z = repmat(blanks(80),n,1);
for i =1:n
    z(i,1:55) = ['MTRIX' num2str(i) ' ',...
        sprintf('%3d%10.6f%10.6f%10.6f     %10.5f',...
        mRec(i).SerNo,...
        mRec(i).Mn1,...
        mRec(i).Mn2,...
        mRec(i).Mn3,...
        mRec(i).Vn)];
    if(mRec(i).iGiven == 1)
        z(i,60) = num2str(1);
    end
end


function z = writeMaster(mRec)
z = repmat(blanks(80),2,1);
z(1,1:70) = sprintf('MASTER    %5d    0%5d%5d%5d%5d%5d%5d%5d%5d%5d%5d',...
    mRec.numREMARK,...
    mRec.numHET,...
    mRec.numHelix,...
    mRec.numSheet,...
    mRec.numTurn,...
    mRec.numSite,...
    mRec.numXform,...
    mRec.numCoord,...
    mRec.numTer,...
    mRec.numConect,...
    mRec.numSeq);
z(2,1:3) = 'END';

%-----------------------------------------------------%
function z = writeOBS_SPR(mRec,recName)
%for OBSLTE and SPRSDE

switch recName
    case 'OBSLTE'
        [noOfLines lineLen] = size(mRec.rIdCode);
        z = repmat(blanks(80),noOfLines,1);
        z(:,1:6) = repmat(recName,noOfLines,1);
        %first line
        z(1,12:lineLen+31) = [mRec.repDate ' ' mRec.idCode '      ' mRec.rIdCode(1,:)];
        if(noOfLines > 1)
            z(2:end,9:10) = reshape(sprintf('%2d',2:noOfLines),2,noOfLines-1)';
            z(2:end,12:lineLen+31) = [ repmat([mRec.repDate ' ' mRec.idCode],noOfLines-1,1),...
                '      ' mRec.rIdCode(2:noOfLines,:)];
        end

    case 'SPRSDE'
        [noOfLines lineLen] = size(mRec.sIdCode);
        z = repmat(blanks(80),noOfLines,1);
        z(:,1:6) = repmat(recName,noOfLines,1);
        %first line
        z(1,12:lineLen+31) = [mRec.Supersededdate ' ' mRec.idCode '      ' mRec.sIdCode(1,:)];
        if(noOfLines > 1)
            z(2:end,9:10) = reshape(sprintf('%2d',2:noOfLines),2,noOfLines-1)';
            z(2:end,12:lineLen+31) = [ repmat([mRec.Supersededdate ' ' mRec.idCode],noOfLines-1,1),...
                '      ' mRec.sIdCode(2:noOfLines,:)];
        end
end

function z= writeOriginX(mRec)
n = length(mRec);
z = repmat(blanks(80),n,1);
for i =1:n
    z(i,1:55) = ['ORIGX' num2str(i) '    ',...
        sprintf('%10.6f%10.6f%10.6f     %10.5f',...
        mRec(i).On1,...
        mRec(i).On2,...
        mRec(i).On3,...
        mRec(i).Tn)];
end

function z = writeREVDAT(mRec)
%REVDAT
out = lineStats(mRec,'record');
noOfRev = size(out,1);
z = repmat(blanks(80),sum(out(:,1)),1);
lc = 1;
for i=1:noOfRev
    noOfLines = out(i,1);
    lineLen = out(i,2);
    if(noOfLines < 1)
        noOfLines = 1;
    end
    y = repmat(blanks(80),noOfLines,1);
    y(:,1:10) = repmat(['REVDAT ' sprintf('%3d',mRec(i).modNum)],noOfLines,1);

    %first line
    if(lineLen)
        y(1,14:39+lineLen) = sprintf('%9s %5s   %d       %s',...
            mRec(i).modDate,...
            mRec(i).modId,...
            mRec(i).modType(1),...
            mRec(i).record(1,:));
    else
        y(1,14:39) = sprintf('%9s %5s   %d       ',...
            mRec(i).modDate,...
            mRec(i).modId,...
            mRec(i).modType(1));
    end
    if noOfLines > 1
        y(2:noOfLines,11:12) = reshape(sprintf('%2d',2:noOfLines),2,noOfLines-1)';
        y(2:noOfLines,32:39+lineLen) = [num2str(mRec(i).modType(2:end)),...
            repmat(blanks(7),noOfLines-1,1),...
            mRec(i).record(2:end,:)];
    end
    z(lc:lc+noOfLines-1,:) = y;
    lc = lc + noOfLines;
end

function z = writeRemark1(mRec)
%REMARK 1
noOfRemark1 = length(mRec);
out = sum(arrayfun(@(x)  size(x.Author,1) + size(x.Title,1) + ...
            size(x.Editor,1) + size(x.Reference,1) + size(x.Publisher,1) + ...
            size(x.CitationReference,1), mRec));
z = repmat(blanks(80),out + noOfRemark1 + 1,1);
z(1,1:10) = 'REMARK   1';
lc = 2;
for i=1:noOfRemark1
    y = blanks(80);
    y(1:70) = ['REMARK   1 REFERENCE ' sprintf('%-49d',i)];
    z(lc,:) = y;
    lc = lc + 1;
    %Optional AUTHOR
    if(~isempty(mRec(i).Author))
        [m n] = size(mRec(i).Author);
        y = repmat(blanks(80),m,1);
        y(:,1:n+19) = [repmat('REMARK   1  AUTH',m,1),...
            reshape(sprintf('%2d ',1:m),3,m)',...
            mRec(i).Author];
        y(1,17:18) = '  ';%remove continuation number for the first author line
        z(lc:lc+m-1,:) = y;
        lc = lc + m;
    end
    %Optional TITLE
    if(~isempty(mRec(i).Title))
        [m n] = size(mRec(i).Title);
        y = repmat(blanks(80),m,1);
        y(:,1:n+19) = [repmat('REMARK   1  TITL',m,1),...
            reshape(sprintf('%2d ',1:m),3,m)'...
            mRec(i).Title];
        y(1,17:18) = '  ';%remove continuation number for the first Title line
        z(lc:lc+m-1,:) = y;
        lc = lc + m;
    end
    %Optional EDIT
    if(~isempty(mRec(i).Editor))
        [m n] = size(mRec(i).Editor);
        y = repmat(blanks(80),m,1);
        y(:,1:n+19) = [repmat('REMARK   1  EDIT',m,1)...
            reshape(sprintf('%2d ',1:m),3,m)'...
            mRec(i).Editor];
        y(1,17:18) = '  ';%remove continuation number for the first Edit line
        z(lc:lc+m-1,:) = y;
        lc = lc + m;
    end
    [m n] = size(mRec(i).Reference);
    y = repmat(blanks(80),m,1);
    y(:,1:n+19) = [repmat('REMARK   1  REF ',m,1)...
        reshape(sprintf('%2d ',1:m),3,m)'...
        mRec(i).Reference];
    y(1,17:18) = '  ';
    z(lc:lc+m-1,:) = y;
    lc = lc + m;
    %Optional PUBL
    if(~isempty(mRec(i).Publisher))
        [m n] = size(mRec(i).Publisher);
        y = repmat(blanks(80),m,1);
        y(:,1:n+19) = [repmat('REMARK   1  PUBL',m,1)...
            reshape(sprintf('%2d ',1:m),3,m)'...
            mRec(i).Publisher];
        y(1,17:18) = '  ';%remove continuation number for the first Edit line
         z(lc:lc+m-1,:) = y;
        lc = lc + m;
    end
    if(~isempty(mRec(i).CitationReference))
        [m n] = size(mRec(i).CitationReference);
        y = repmat(blanks(80),m,1);
        y(:,1:n+19) = [repmat('REMARK   1  REFN',m,1)...
            reshape(sprintf('%2d ',1:m),3,m)'...
            mRec(i).CitationReference];
        y(1,17:18) = '  ';
        z(lc:lc+m-1,:) = y;
        lc = lc + m;
    end
end

%---------------------------------------------------------%
function z = writeRemark2(mRec)
%for REMARK 2
z = blanks(80);
z(1:10) = 'REMARK   2';
if(isfield(mRec,'Resolution'))
    y = blanks(80);
    y(1:38) = ['REMARK   2 RESOLUTION.' sprintf('%5.2f ANGSTROMS.',mRec.Resolution)];
    z = [z;y];
end
[m n] = size(mRec.Detail);
y = repmat(blanks(80),m,1);
y(:,1:11+n) = [repmat('REMARK   2 ',m,1) mRec.Detail];
z=[z;y];

%----------------------------------------------------------%
function z = writeRemark3(mRec)
%for REMARK 3
z = blanks(80);
z(1:10) = 'REMARK   3';
[m n] = size(mRec.Refinement);
y = repmat(blanks(80),m,1);
y(:,1:11+n) = [repmat('REMARK   3 ',m,1) mRec.Refinement];
z=[z;y];
%----------------------------------------------------------%
function z = writeRemarkN(mRec,remarkNo)
%for REMARK 4-999
[m n] = size(mRec);
z = repmat(blanks(80),m,1);
z(:,1:11+n) = [repmat(sprintf('REMARK %3d ',remarkNo),m,1) mRec];
%----------------------------------------------------------%
function z = writeSSbond(mRec)
noOfLines = length(mRec);
z = repmat(blanks(80),noOfLines,1);
for i=1:noOfLines
    z(i,1:72) = ['SSBOND ' ...
        sprintf('%3d %3s %1s %4d%1s   %3s %1s %4d%1s                       %6s %6s',...
        mRec(i).serNum,...
        mRec(i).resName1,...
        mRec(i).chainID1,...
        mRec(i).seqNum1,...
        mRec(i).icode1,...
        mRec(i).resName2,...
        mRec(i).chainID2,...
        mRec(i).seqNum2,...
        mRec(i).icode2,...
        mRec(i).sym1,...
        mRec(i).sym2)];
    z(i,:) = strrep(z(i,:),'NaN','   ');
end

%----------------------------------------------------------%
function z = writeSaltBridge(mRec)
noOfLines = length(mRec);
z = repmat(blanks(80),noOfLines,1);
for i=1:noOfLines
    z(i,1:72) = ['SLTBRG      ' ...
        sprintf('%4s%1s%3s %1s%4d%1s               %4s%1s%3s %1s%4d%1s  %6s %6s',...
        mRec(i).AtomName1,...
        mRec(i).altLoc1,...
        mRec(i).resName1,...
        mRec(i).chainID1,...
        mRec(i).resSeq1,...
        mRec(i).iCode1,...
        mRec(i).AtomName2,...
        mRec(i).altLoc2,...
        mRec(i).resName2,...
        mRec(i).chainID2,...
        mRec(i).resSeq2,...
        mRec(i).iCode2,...
        mRec(i).sym1,...
        mRec(i).sym2)];
    z(i,:) = strrep(z(i,:),'NaN','   ');
end

%---------------------------------------------%
function z= writeScale(mRec)
n = length(mRec);
z = repmat(blanks(80),n,1);
for i =1:n
    z(i,1:55) = ['SCALE' num2str(i) '    '...
        sprintf('%10.6f%10.6f%10.6f     %10.5f',...
        mRec(i).Sn1,...
        mRec(i).Sn2,...
        mRec(i).Sn3,...
        mRec(i).Un)];
end


function z = writeSeqAdv(mSeqAdv)
noOfSeqAdv = length(mSeqAdv);
z = repmat(blanks(80),noOfSeqAdv,1);
for i=1:noOfSeqAdv
    z(i,1:49+length(mSeqAdv(i).conflict)) = ['SEQADV ' ...
        sprintf('%4s %3s %1s %4d%1s %4s %9s %3s %5d %s',...
        mSeqAdv(i).idCode,...
        mSeqAdv(i).resName,...
        mSeqAdv(i).chainID,...
        mSeqAdv(i).seqNum,...
        mSeqAdv(i).iCode,...
        mSeqAdv(i).database,...
        mSeqAdv(i).dbIdCode,...
        mSeqAdv(i).dbRes,...
        mSeqAdv(i).dbSeq,...
        mSeqAdv(i).conflict)];
end


function z = writeSheet(mRec)
noOfLines = length(mRec);
z = repmat(blanks(80),noOfLines,1);
for i=1:noOfLines
    z(i,1:70) = ['SHEET  ' ...
        sprintf('%3d %3s%2d %3s %1s%4d%1s %3s %1s%4d%1s%2d %4s%3s %1s%4d%1s %4s%3s %1s%4d%1s',...
        mRec(i).strand,mRec(i).sheetID,...
        mRec(i).numStrands,...
        mRec(i).initResName,...
        mRec(i).initChainID,...
        mRec(i).initSeqNum,...
        mRec(i).initICode,...
        mRec(i).endResName,...
        mRec(i).endChainID,...
        mRec(i).endSeqNum,...
        mRec(i).endICode,...
        mRec(i).sense,...
        mRec(i).curAtom,...
        mRec(i).curResName,...
        mRec(i).curChainId,...
        mRec(i).curResSeq,...
        mRec(i).curICode,...
        mRec(i).prevAtom,...
        mRec(i).prevResName,...
        mRec(i).prevChainId,...
        mRec(i).prevResSeq,...
        mRec(i).prevICode)];
    z(i,:) = strrep(z(i,:),'NaN','   ');
end

function z = writeSingleRecMultipleLine(mRec,recordName)
% for AUTHOR, COMPND, EXPDAT, KEYWDS, SOURCE, TITLE
[noOfLines lineLen] = size(mRec);
z = repmat(blanks(80),noOfLines,1);
z(:,1:length(recordName)) = repmat(recordName,noOfLines,1);
%first line
z(1,11:10+lineLen) = mRec(1,:);
if(noOfLines > 1)
    z(2:end,9:10) = reshape(sprintf('%2d',2:noOfLines),2,noOfLines-1)';
    %shift lines 2:n one character to match format in pdb files
    z(2:end,12:11+lineLen) = mRec(2:noOfLines,:);
end

function z = writeSite(mRec)
% SITE
out = lineStats(mRec,'SeqNo');
noOfSite = size(out,1);
z = repmat(blanks(80),sum(out(:,1)),1);
lc = 1;

for i=1:noOfSite
    noOfLines = out(i,:);
    y = repmat(blanks(80),noOfLines,1);
    y(:,1:17) = repmat(['SITE       ' sprintf('%3s %2d',mRec(i).SiteName,mRec(i).NumRes)],...
        noOfLines,1);
    for j = 1:noOfLines
        y(j,8:10) = sprintf('%3d',mRec(i).SeqNo(j));
    end
    %ResDet
    n = mRec(i).NumRes;
    k = 1;
    idx = 19;
    for j=1:n
        y(k,idx:idx+9) = sprintf('%3s %1s%4d%1s',...
            mRec(i).ResDet(j).ResName,...
            mRec(i).ResDet(j).ChainID,...
            mRec(i).ResDet(j).ResSeqNo,...
            mRec(i).ResDet(j).InsCode);
        idx = idx + 11;
        if(~mod(j,4))
            k = k + 1;
            idx = 19;
        end
    end
    z(lc:lc+noOfLines-1,:) = y;
    lc = lc + noOfLines;
end

%------------------------------------------------------------------------%
function z = writeTVECT(mRec)
noOfLines = length(mRec);
z = repmat(blanks(80),noOfLines,1);
for i=1:noOfLines
    z(i,1:70) = ['TVECT  ' ...
        sprintf('%3d%10.5f%10.5f%10.5f%30s',...
        mRec(i).SerNo,...
        mRec(i).t1,...
        mRec(i).t2,...
        mRec(i).t3,...
        mRec(i).text)];
    z(i,:) = strrep(z(i,:),'NaN','   ');
end

%---------------------------------------------------------------%
function z = writeTurn(mRec)
noOfLines = length(mRec);
z = repmat(blanks(80),noOfLines,1);
for i=1:noOfLines
    z(i,1:70) = ['TURN   ' ...
        sprintf('%3d %3s %3s %1s%4d%1s %3s %1s%4d%1s    %30s',...
        mRec(i).seq,mRec(i).turnId,...
        mRec(i).initResName,...
        mRec(i).initChainId,...
        mRec(i).initSeqNum,...
        mRec(i).initICode,...
        mRec(i).endResName,...
        mRec(i).endChainId,...
        mRec(i).endSeqNum,...
        mRec(i).endICode,...
        mRec(i).comment)];
    z(i,:) = strrep(z(i,:),'NaN','   ');
end
%--------------------------------------------------------------------------
function out = lineStats(mRec,fieldName)
recCount = length(mRec);
out = zeros(recCount,2);
[out(:,1) out(:,2)] = arrayfun(@(x) size(x.(fieldName)),mRec);


%----------------------------------------------------------%
function out = writeModel(model_struct)
% Returns coordinate records in a Model field
AtomTypeNames = {'Atom',...
                 'HeterogenAtom',...
                 'AtomSD',...
                 'AnisotropicTemp',...
                 'AnisotropicTempSD',...
                 'Terminal'};
             
ntypes = numel(AtomTypeNames);
type_counts = zeros(ntypes, 1);

atomSerNo = [];
for i = 1:ntypes
    if isfield(model_struct, AtomTypeNames{i})
        type_counts(i) = numel(model_struct.(AtomTypeNames{i}));
        atomSerNo = [atomSerNo; getAtomSerNo( model_struct.(AtomTypeNames{i}))]; %#ok
    end
end

total_counts = sum(type_counts);
atomRecords = repmat(blanks(80),total_counts,1);

% ATOM
start_count = 1;
end_count = type_counts(1);
count = 1;
for i = start_count:end_count
    x = model_struct.Atom(count);
    atomRecords(i,:) = ['ATOM  ',...
        sprintf('%5d %2s%1s%1s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s%2s',...
        x.AtomSerNo,...
        x.AtomNameStruct.chemSymbol,...
        x.AtomNameStruct.remoteInd,...
        x.AtomNameStruct.branch, ...
        x.altLoc,x.resName,...
        x.chainID,...
        x.resSeq,...
        x.iCode,...
        x.X,...
        x.Y,...
        x.Z,...
        x.occupancy,...
        x.tempFactor,...
        x.segID,...
        x.element,...
        x.charge)];
    count = count +1;
end

% HETATM
start_count = end_count + 1;
end_count = end_count + type_counts(2);
count = 1;
for i = start_count:end_count
    x = model_struct.HeterogenAtom(count);
    atomRecords(i,:) = ['HETATM',...
        sprintf('%5d %2s%1s%1s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s%2s',...
        x.AtomSerNo,...
        x.AtomNameStruct.chemSymbol,...
        x.AtomNameStruct.remoteInd,...
        x.AtomNameStruct.branch,...
        x.altLoc,...
        x.resName,...
        x.chainID,...
        x.resSeq,...
        x.iCode,...
        x.X,...
        x.Y,...
        x.Z,...
        x.occupancy,...
        x.tempFactor,...
        x.segID,...
        x.element,...
        x.charge)];
    count = count +1;
end

% SIGATM
start_count = end_count + 1;
end_count = end_count + type_counts(3);
count = 1;
for i = start_count:end_count
    x = model_struct.AtomSD(count);
    atomRecords(i,:) = ['SIGATM',...
        sprintf('%5d %2s%1s%1s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s%2s',...
        x.AtomSerNo,...
        x.AtomNameStruct.chemSymbol,...
        x.AtomNameStruct.remoteInd,...
        x.AtomNameStruct.branch,...
        x.altLoc,...
        x.resName,...
        x.chainID,...
        x.resSeq,...
        x.iCode,...
        x.sigX,...
        x.sigY,...
        x.sigZ,...
        x.sigOcc,...
        x.sigTemp,...
        x.segID,...
        x.element,...
        x.charge)];
    count = count +1;
end

% ANISOU
start_count = end_count + 1;
end_count = end_count + type_counts(4);
count = 1;
for i = start_count:end_count
x = model_struct.AnisotropicTemp(count);
           atomRecords(i,:) = ['ANISOU',...
                sprintf('%5d %2s%1s%1s%1s%3s %1s%4d%1s %7d%7d%7d%7d%7d%7d  %4s%2s%2s',...
                x.AtomSerNo,...
                x.AtomNameStruct.chemSymbol,...
                x.AtomNameStruct.remoteInd,...
                x.AtomNameStruct.branch,...
                x.altLoc,...
                x.resName,...
                x.chainID,...
                x.resSeq,...
                x.iCode,...
                x.U00,...
                x.U11,...
                x.U22,...
                x.U01,...
                x.U02,...
                x.U12,...
                x.segID,...
                x.element,...
                x.charge)];
  count = count +1;
end

% SIGUIJ
start_count = end_count + 1;
end_count = end_count + type_counts(5);
count = 1;
for i = start_count:end_count
    x = model_struct.AnisotropicTempSD(count);
    atomRecords(i,:) = ['SIGUIJ',...
        sprintf('%5d %2s%1s%1s%1s%3s %1s%4d%1s %7d%7d%7d%7d%7d%7d  %4s%2s%2s',...
        x.AtomSerNo,...
        x.AtomNameStruct.chemSymbol,...
        x.AtomNameStruct.remoteInd,...
        x.AtomNameStruct.branch,...
        x.altLoc,...
        x.resName,...
        x.chainID,...
        x.resSeq,...
        x.iCode,...
        x.SIG11,...
        x.SIG22,...
        x.SIG33,...
        x.SIG12,...
        x.SIG13,...
        x.SIG23,...
        x.segID,...
        x.element,...
        x.charge)];
    count = count +1;
end

% TER
start_count = end_count + 1;
end_count = end_count + type_counts(6);
count = 1;
for i = start_count:end_count
    x = model_struct.Terminal(count);
    atomRecords(i,1:27) = ['TER   ',...
        sprintf('%5d      %3s %1s%4d%1s',...
        x.SerialNo,...
        x.resName,...
        x.chainID,...
        x.resSeq,...
        x.iCode)];
    count = count +1;
end



% Sort the output according to the atomSerNo
[atomNO_sort, sort_idx] = sort(atomSerNo);
out = atomRecords(sort_idx, :);

%---------------------------------------------------------%
function atomSerNo = getAtomSerNo(atom_struct)
% Returns AtomSerNo for records in atom_struct
atomSerNo = [];
N = numel(atom_struct);

serNoFieldName = 'AtomSerNo';
if isfield(atom_struct, 'SerialNo')
    serNoFieldName = 'SerialNo';
end

if N ~= 0;
    atomSerNo = zeros(N, 1);
    for i = 1:N
        atomSerNo(i) = atom_struct(i).(serNoFieldName);
    end
end

%---------------------------------------------------------%
function isoldstruct = isOldpdb_st(oldpdb_st)
% Check if the input structure is an old PDB structure

isoldstruct = false;
% Field names for the coordinate information.
cofields = {'Atom',...
            'AtomSD',...
            'AnisotropicTemp',...
            'AnisotropicTempSD',...
            'Terminal',...
            'HeterogenAtom'};
        
realfields = isfield(oldpdb_st, cofields); 
if any(realfields)
    isoldstruct = true;
end
