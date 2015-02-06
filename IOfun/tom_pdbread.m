function PDB_struct = tom_pdbread(filename)
%tom_pdbread reads a Protein Data Bank file into a MATLAB structure.
%
%  Adapted by SN for the TOM-Suite from the original file PDBRead by S.A.
%  Nikumbh.
%  PDB_struct = PDBRead(filename) reads the file corresponding to filename and stores the information
%  contained in this file in the PDB_struct.
%  e.g PDB_struct = PDBRead('PDBSilk.txt')
%
%  The file that is being read by PDBRead should be compatible to the PDB file format described at the 
%  RCSB (Research Collaboratory for Structural Bioinformatics) web site. For more information about this 
%  format, please visit the following URL:
%  http://www.rcsb.org/pdb/docs/format/pdbguide2.2/guide2.2_frame.html
%  
%  The file 'PDBSilk.txt' mentioned in the sample code above was obtained from the following web site :
%  http://cmm.info.nih.gov/modeling/pdb_at_a_glance.html
%  
%  The current version of this function can recognize the standard record types mentioned at the above URL.
%  The information stored in the output structure can be extracted by using the . (dot) operator similar to 
%  a normal structure in MATLAB. 
%  e.g PDB_struct.HEADER will give the information about the HEADER record type.
%  
%  Each record type has format fields. For more information about the fields,please visit the above mentioned URL.
%  Information for each field for a record type can be accessed by using the . (dot) operator.
%  e.g. PDB_struct.HEADER.depDate gives the deposition date for the particular molecule structure.
%  The field names used for a particular record type are according to the information contained  at the 
%  RCSB web site.
%
%  Author : Sachin A. Nikumbh 
%  Email address : snikumbh@mathworks.com
%  For comments regarding this function please contact the author at above email address.
%  Adapted for TOM and Matlab 6.5 by SN
%    Copyright (c) 2004
%    TOM toolbox for Electron Tomography
%    Max-Planck-Institute for Biochemistry
%    Dept. Molecular Structural Biology
%    82152 Martinsried, Germany
%    http://www.biochem.mpg.de/tom


fid=fopen(filename,'r');

if fid == -1,
    error('Unable to open specified file');
else
    PDB_struct = [];
    PDB_struct = Initialize_struct(PDB_struct); % If the file is opened successfully, initialize the PDB structure
    
    % SingContd is the array for the FLAGS for the Single Continued type records
    % SingContd(1)=OBSLTE ,SingContd(2)=TITLE, SingContd(3)=CAVEAT, SingContd(4)=COMPND 
    % SingContd(5)=SOURCE, SingContd(6)=KEYWDS, SingContd(7)=EXPDTA, SingContd(8)=AUTHOR
    % SingContd(9)=SPRSDE  
    % These flags will be set when the first line containing these flags appears 
    SingContd = zeros(1,9);
    
    % Following are some of the initializations necessary for various record types mentioned in the comments
    JRNL_Entry = 1; % This is the count for the number of the JRNL records in the PDB file
    JRNLSubRecords = zeros(1,6); % The array of flags for the subrecords. JRNLSubRecords(1)=AUTH, JRNLSubRecords(2)=TITL, JRNLSubRecords(3)=EDIT
                                 %                                        JRNLSubRecords(4)=REF, JRNLSubRecords(5)=PUBL , JRNLSubRecords(6)=REFN  
    OldNoOfKeywds = 0; % KEYDS
    OldNoOfEntries = 0; % OBSLTE  
    OldNoOfAuthors = 0; % AUTHOR
    NumOfREVDAT = 0; % REVDAT
    NumOfDBREF = 0; % DBREF
    NumOfSEQADV = 0; % SEQADV
    NumOfMODRES = 0; % MODRES
    NumOfHET = 0; % HET
    NumOfHETNAM = 0; % HETNAM
    NumOfHELIX = 0; % HELIX
    NumOfSHEET = 0; % SHEET
    NumOfTURN = 0; %TURN
    NumOfSSBOND = 0; % SSBOND
    NumOfLINK = 0; % LINK
    NumOfHYDBND = 0; % HYDBND
    NumOfSLTBRG =0; % SLTBRG
    NumOfCISPEP = 0; % CISPEP
    NumOfTVECT = 0; % TVECT
    NumOfCONECT = 0; % CONECT
    NumOfSIGATM = 0; % SIGATM
    NumOfANISOU = 0; % ANISOU
    NumOfSIGUIJ = 0; % SIGUIJ
    NumOfATOM = 0; % ATOM
    NumOfHETSYN = 0; % HETSYN
    NumOfFORMUL = 0; % FORMUL
    NumOfHETNAM = 0; % HETNAM
    NumOfHETATM = 0; % HETATM
    NumOfMODEL = 0; % MODEL
    NumOfENDMDL = 0; % ENDMDL
    NumOfTER = 0; % TER
    JRNLAuth = 0;
    
    MODELFlag = 0; % This is the flag to mark the beginning of the MODEL record. 
                   % This will be set when MODEL record appears and will be reset to 0 when ENDMDL record appears
    
    CurRes = 0; % changed form '' to 0 FF 
    PrevRes = 0; % changed form '' to 0 FF 
    
    CurSITEName = '';
    PrevSITEName = '';
    
    CurHetIDHETSYN = '';
    PrevHetIDHETSYN = '';
    
    CurHetIDFORMUL = '';
    PrevHetIDFORMUL = '';
    
    CurHetIDHETNAM = '';
    PrevHetIDHETNAM = '';
    
    CurREMARK = 0;
    PrevREMARK = 0;
    
    TmpStruct = struct('ResName',{},'ChainID',{''},'ResSeqNo',{},'InsCode',{''}); % This is for SITE record
end

disp('Extracting the information............');

while 1
    
    tline = fgetl(fid);
    
    if ~ischar(tline)
        break; % For end of file recognition
    end
    
    if size(tline)>0 % Omit the empty lines to avoid error of invalid matrix index.
             
            sz = size(tline);
            tline = [tline blanks(80-sz(2))]; % RCSB web site format requires each line to have 80 characters. This avoids exceeding the matrix dimension for lines with  
                                              % less than 80 characters.  
            Record_name = upper(tline(1:6));
            Record_name = deblank(Record_name); % Assuming that the record name will be left alligned (as mentioned in the RCSB file format doc,remove trailing blanks
    
            if strncmp(Record_name,'ORIGX',5) | strncmp(Record_name,'SCALE',5) | strncmp(Record_name,'MTRIX',5) % This is done to take care of ORIGX1,ORIGX2,ORIGX3 
                                                                                                                % and similarly for SCALE and MTRIX
                Record_name = Record_name(1:5);                                                                 
            end
       
            switch Record_name
        
                %Single/Mandatory
                case 'HEADER' 
                    PDB_struct.HEADER = struct('name',{'HEADER'},'classification',{deblank(tline(11:50))},'depDate',{tline(51:59)},'idCode',{tline(63:66)});
                   
                %Single Continued/Optional : mandatory in withdrawn entries    
                case 'OBSLTE'
            
                    if ~SingContd(1)
                        SingContd(1) = 1; % set the flag
                        PDB_struct.OBSLTE = struct('name',{'OBSLTE'},'repDate',{tline(12:20)},'idCode',{tline(22:25)},'rIdCode',{removeblanks(tline(32:70))});
                    else
                        PDB_struct.OBSLTE.rIdCode = strvcat(PDB_struct.OBSLTE.rIdCode,removeblanks(tline(32:70)));
                    end
        
                %Single Continued/Mandatory    
                case 'TITLE'
                    PDB_struct.TITLE.title = strvcat(PDB_struct.TITLE.title,removeblanks(tline(11:70)));
                    
                %Single Continued/Optional               
                case 'CAVEAT'
                    if ~SingContd(3)
                        SingContd(3) = 1;
                        PDB_struct.CAVEAT =  struct('name',{'CAVEAT'},'idCode',{tline(12:15)},'comment',{removeblanks(tline(20:70))})
                    else
                        PDB_struct.CAVEAT.comment = strvcat(PDB_struct.CAVEAT.comment,removeblanks(tline(20:70)));
                    end
                        
                %Single Continued/Mandatory       
                case 'COMPND'
                    PDB_struct.COMPND.comp_description = strvcat(PDB_struct.COMPND.comp_description,removeblanks(tline(11:70)));
            
                %Single Continued/Mandatory     
                case 'SOURCE'
                    PDB_struct.SOURCE.src_description = strvcat(PDB_struct.SOURCE.src_description,removeblanks(tline(11:70)));
                           
                %Single Continued/Mandatory     
                case 'KEYWDS'
                    PDB_struct.KEYWDS.KeywdsList = strvcat(PDB_struct.KEYWDS.KeywdsList,removeblanks(tline(11:70)));
                           
                %Single Continued/Mandatory     
                case 'EXPDTA'
                    PDB_struct.EXPDTA.technique = strvcat(PDB_struct.EXPDTA.technique,removeblanks(tline(11:70)));
                               
                %Single Continued/Mandatory                   
                case 'AUTHOR'
                    PDB_struct.AUTHOR.AuthorsList = strvcat(PDB_struct.AUTHOR.AuthorsList,removeblanks(tline(11:70)));
                
                %Multiple/Mandatory    
                case 'REVDAT'
                    NumOfREVDAT = NumOfREVDAT+1;
                    PDB_struct.REVDAT(NumOfREVDAT) = struct('name',{'REVDAT'},'modNum',{str2num(tline(8:10))},'modDate',{tline(14:22)},'modId',{tline(24:28)},...
                                                            'modType',{str2num(tline(32))},'record',{removeblanks(tline(40:66))});
            
                %Single Continued/Optional      
                case 'SPRSDE'
                    
                    if ~SingContd(9)
                        OldEntries = 0;
                        SingContd(9) = 1;
                        PDB_struct.SPRSDE = struct('name',{'SPRSDE'},'sprsdedate',{tline(12:20)},'idCode',{tline(22:25)},'sIdCode',{removeblanks(tline(32:70))});
                    else
                        PDB_struct.SPRSDE.sIdCode = strvcat(PDB_struct.SPRSDE.sIdCode,removeblanks(tline(32:70)));
                    end
                
                %Other/Optional : This record has following sub-records: AUTH, TITL,EDIT,REF,PUBL,REFN        
                case 'JRNL'
                    
                    SubRecord = tline(13:16);
                    SubRecord = deblank(SubRecord); % Remove the trailing blanks. Needed for REF
                    
                    switch SubRecord
                        
                        case 'AUTH'
                            
                            if ~JRNLAuth
                                JRNLAuth = 1;
                                PDB_struct.JRNL.Entry(JRNL_Entry) = struct('AUTH',{''},'TITL',{''},'EDIT',{''},'REF',{''},'PUBL',{''},'REFN',{''});
                                PDB_struct.JRNL.Entry(JRNL_Entry).AUTH = strvcat(PDB_struct.JRNL.Entry(JRNL_Entry).AUTH,removeblanks(tline(20:70)));
                            else
                                PDB_struct.JRNL.Entry(JRNL_Entry).AUTH = strvcat(PDB_struct.JRNL.Entry(JRNL_Entry).AUTH,removeblanks(tline(20:70)));                                
                            end
                                                                               
                        case 'TITL'
                            PDB_struct.JRNL.Entry(JRNL_Entry).TITL = strvcat(PDB_struct.JRNL.Entry(JRNL_Entry).TITL,removeblanks(tline(20:70)));
                                                       
                        case 'EDIT'
                            PDB_struct.JRNL.Entry(JRNL_Entry).EDIT = strvcat(PDB_struct.JRNL.Entry(JRNL_Entry).EDIT,removeblanks(tline(20:70)));
                                                        
                        case 'REF'
                            PDB_struct.JRNL.Entry(JRNL_Entry).REF = strvcat(PDB_struct.JRNL.Entry(JRNL_Entry).REF,removeblanks(tline(20:70)));
                                                        
                        case 'PUBL'
                            PDB_struct.JRNL.Entry(JRNL_Entry).PUBL = strvcat(PDB_struct.JRNL.Entry(JRNL_Entry).PUBL,removeblanks(tline(20:70)));
                                                        
                        case 'REFN'
                            PDB_struct.JRNL.Entry(JRNL_Entry).REFN = strvcat(PDB_struct.JRNL.Entry(JRNL_Entry).REFN,removeblanks(tline(20:70)));
                            JRNL_Entry = JRNL_Entry+1; % REFN is the last subrecord and it is a single line record
                            JRNLAuth = 0;
                                                                            
                        otherwise
                            %disp('Invalid subrecord type');
                    end
                    
                    PDB_struct.JRNL.NoOfJRNLS = JRNL_Entry-1;
                
                % Some of the REMARK records are mandatory and some are optional    
                case 'REMARK'
                    
                    RemarkNo = str2num(tline(7:10));
                    
                    %Other/Optional
                    if RemarkNo == 1
                        
                        if strcmp(tline(12:20),'REFERENCE')
                            PDB_struct.REMARK1.NoOfJRNLS = str2num(tline(22:70));
                            PDB_struct.REMARK1.JRNLEntry(PDB_struct.REMARK1.NoOfJRNLS) = struct('AUTH',{''},'TITL',{''},'EDIT',{''},'REF',{''},'PUBL',{''},'REFN',{''});
                        else
                            SubRecord = tline(13:16);
                            SubRecord = deblank(SubRecord); % Remove the trailing blanks. Needed for REF
                    
                            switch SubRecord
                        
                                case 'AUTH'
                                    PDB_struct.REMARK1.JRNLEntry(PDB_struct.REMARK1.NoOfJRNLS).AUTH = strvcat(PDB_struct.REMARK1.JRNLEntry(PDB_struct.REMARK1.NoOfJRNLS).AUTH,...
                                                                                                                                           removeblanks(tline(20:70)));
                                                                               
                                case 'TITL'
                                    PDB_struct.REMARK1.JRNLEntry(PDB_struct.REMARK1.NoOfJRNLS).TITL = strvcat(PDB_struct.REMARK1.JRNLEntry(PDB_struct.REMARK1.NoOfJRNLS).TITL,...
                                                                                                                                           removeblanks(tline(20:70)));
                                                       
                                case 'EDIT'
                                    PDB_struct.REMARK1.JRNLEntry(PDB_struct.REMARK1.NoOfJRNLS).EDIT = strvcat(PDB_struct.REMARK1.JRNLEntry(PDB_struct.REMARK1.NoOfJRNLS).EDIT,...
                                                                                                                                           removeblanks(tline(20:70)));
                                                        
                                case 'REF'
                                    PDB_struct.REMARK1.JRNLEntry(PDB_struct.REMARK1.NoOfJRNLS).REF = strvcat(PDB_struct.REMARK1.JRNLEntry(PDB_struct.REMARK1.NoOfJRNLS).REF,...
                                                                                                                                          removeblanks(tline(20:70)));
                                                        
                                case 'PUBL'
                                    PDB_struct.REMARK1.JRNLEntry(PDB_struct.REMARK1.NoOfJRNLS).PUBL = strvcat(PDB_struct.REMARK1.JRNLEntry(PDB_struct.REMARK1.NoOfJRNLS).PUBL,...
                                                                                                                                           removeblanks(tline(20:70)));
                                                        
                                case 'REFN'
                                    PDB_struct.REMARK1.JRNLEntry(PDB_struct.REMARK1.NoOfJRNLS).REFN = strvcat(PDB_struct.REMARK1.JRNLEntry(PDB_struct.REMARK1.NoOfJRNLS).REFN,...
                                                                                                                                           removeblanks(tline(20:70)));
                                                                               
                                otherwise
                                    %disp('Invalid subrecord type');
                            end
                                               
                        end
                    
                    % Other/Mandatory:  Following code assumes a single occurance of the REMARK2 record type   
                    elseif RemarkNo == 2
                        
                        if strcmp(tline(12:22),'RESOLUTION.')
                            
                            if strcmp(tline(28:37),'ANGSTROMS.')
                                PDB_struct.REMARK2.Resolution = str2num(tline(23:27));
                            else
                                PDB_struct.REMARK2.Detail = strvcat(PDB_struct.REMARK2.Detail, removeblanks(tline(12:70)));
                            end
                            
                        else
                            PDB_struct.REMARK2.Detail = strvcat(PDB_struct.REMARK2.Detail, removeblanks(tline(12:70)));
                        end
                                                                        
                    % Other/Mandatory:  Following code assumes a single occurance of the REMARK2 record type        
                    elseif RemarkNo == 3
                                         
                        PDB_struct.REMARK3.Refinement = strvcat(PDB_struct.REMARK3.Refinement,removeblanks(tline(12:70)));
                            
                    else
                        
                        CurREMARK = RemarkNo;
                        tmpRemark = sprintf('%d',CurREMARK);
                        temp_struct = ['PDB_struct.' 'REMARK' tmpRemark '.comment'];
                        
                        if CurREMARK ~= PrevREMARK
                            temp_str = [temp_struct '= tline(12:70);'];
                            PrevREMARK = CurREMARK;
                        else
                            temp_str = [temp_struct '= strvcat(' temp_struct ',tline(12:70));'];
                        end
                        
                        eval(temp_str);
                        
                    end
                                   
                %Multiple/Optional     
                case 'DBREF'
                    NumOfDBREF = NumOfDBREF+1;
                    PDB_struct.DBREF(NumOfDBREF) = struct('name',{'DBREF'},'idCode',{tline(8:11)},'chainID',{tline(13)},'seqBegin',{str2num(tline(15:18))},...
                                                           'insertBegin',{tline(19)},'seqEnd',{str2num(tline(21:24))},'insertEnd',{tline(25)},'database',{tline(27:32)},...
                                                           'dbAccession',{tline(34:41)},'dbIdCode',{tline(43:54)},'dbseqBegin',{str2num(tline(56:60))},'idbnsBeg',{tline(61)},...
                                                            'dbseqEnd',{str2num(tline(63:67))},'dbinsEnd',{tline(68)});
                %Multiple/Optional     
                case 'SEQADV'
                    NumOfSEQADV = NumOfSEQADV+1;
                    PDB_struct.SEQADV(NumOfSEQADV) = struct('name',{'SEQADV'},'idCode',{tline(8:11)},'resName',{tline(13:15)},'chainID',{tline(17)},'seqNum',{str2num(tline(19:22))},...
                                                            'iCode',{tline(23)},'database',{tline(25:28)},'dbIdCode',{tline(30:38)},'dbRes',{tline(40:42)},'dbSeq',{str2num(tline(44:48))},...
                                                            'conflict',{removeblanks(tline(50:70))});
                %Multiple/Optional
                case 'SEQRES'
                    
                    if isspace(tline(12))
                        CurRes = sprintf('%d',str2num(tline(14:17)));
                    else
                        CurRes = tline(12);                        
                    end
                    
                    if CurRes ~= PrevRes % 
                        PDB_struct.SEQRES.NoOfResChain = PDB_struct.SEQRES.NoOfResChain + 1;
                        PDB_struct.SEQRES.ResChainDetail(PDB_struct.SEQRES.NoOfResChain).NoOfResidue = str2num(tline(14:17));
                        PDB_struct.SEQRES.ResChainDetail(PDB_struct.SEQRES.NoOfResChain).ChainID = tline(12);
                        PDB_struct.SEQRES.ResChainDetail(PDB_struct.SEQRES.NoOfResChain).AminoAcids = removeblanks(GetAminoAcids(tline(20:70)));
                        PrevRes = CurRes;
                    else
if ~isequal(PDB_struct.SEQRES.NoOfResChain,0)
PDB_struct.SEQRES.ResChainDetail(PDB_struct.SEQRES.NoOfResChain).AminoAcids = strcat(PDB_struct.SEQRES.ResChainDetail(PDB_struct.SEQRES.NoOfResChain).AminoAcids,GetAminoAcids(tline(20:70)));
                    end
                    end
                
                %Multiple/Optional    
                case 'MODRES'
                    NumOfMODRES = NumOfMODRES+1;
                    PDB_struct.MODRES(NumOfMODRES) = struct('name',{'MODRES'},'idCode',{tline(8:11)},'resName',{tline(13:15)},'chainID',{tline(17)},'seqNum',{str2num(tline(19:22))},...
                                                            'iCode',{tline(23)},'stdRes',{tline(25:27)},'comment',{removeblanks(tline(30:70))});            
                                                        
                %Multiple/Optional                                        
                case 'HET'
                    NumOfHET = NumOfHET+1;
                    PDB_struct.HET(NumOfHET) = struct('name',{'HET'},'hetID',{tline(8:10)},'ChainID',{tline(13)},'seqNum',{str2num(tline(14:17))},'iCode',{tline(18)},...
                                                      'numHetAtoms',{str2num(tline(21:25))},'text',{removeblanks(tline(31:70))});
                    
                %Multiple Continued/Optional                                   
                case 'HETNAM'  
                    CurHetIDHETNAM = tline(12:14);
                    
                    if ~strcmp(CurHetIDHETNAM,PrevHetIDHETNAM)
                        NumOfHETNAM = NumOfHETNAM + 1;
                        PDB_struct.HETNAM(NumOfHETNAM).hetID = CurHetIDHETNAM;
                        PDB_struct.HETNAM(NumOfHETNAM).ChemName = removeblanks(tline(16:70));
                        PrevHetIDHETNAM = CurHetIDHETNAM;
                    else
                        PDB_struct.HETNAM(NumOfHETNAM).ChemName = strvcat(PDB_struct.HETNAM(NumOfHETNAM).ChemName,removeblanks(tline(16:70)));
                    end
                        
                %Multiple/Optional    
                case 'HETSYN'
                    CurHetIDHETSYN = tline(12:14);
                    
                    if ~strcmp(CurHetIDHETSYN,PrevHetIDHETSYN)
                        NumOfHETSYN = NumOfHETSYN+1;
                        PDB_struct.HETSYN(NumOfHETSYN).hetID = CurHetIDHETSYN;
                        PDB_struct.HETSYN(NumOfHETSYN).hetSynonyms = {removeblanks(tline(16:70))};
                        PrevHetIDHETSYN = CurHetIDHETSYN;
                    else
                        PDB_struct.HETSYN(NumOfHETSYN).hetSynonyms = strvcat(PDB_struct.HETSYN(NumOfHETSYN).hetSynonyms,removeblanks(tline(16:70)));
                    end
                        
                %Multiple Continued/Optional           
                case 'FORMUL' 
                    CurHetIDFORMUL = tline(13:15);
                    
                    if ~strcmp(CurHetIDFORMUL,PrevHetIDFORMUL)
                        NumOfFORMUL = NumOfFORMUL+1;
                        PDB_struct.FORMUL(NumOfFORMUL).CompNo = str2num(tline(9:10));
                        PDB_struct.FORMUL(NumOfFORMUL).hetID = tline(13:15);
                        PDB_struct.FORMUL(NumOfFORMUL).ChemForm = removeblanks(tline(19:70));
                        PrevHetIDFORMUL = CurHetIDFORMUL;
                    else
                        PDB_struct.FORMUL(NumOfFORMUL).ChemForm = strvcat(PDB_struct.FORMUL(NumOfFORMUL).ChemForm,removeblanks(tline(19:70)));                    
                    end
                            
                %Multiple/Optional     
                case 'HELIX'
                    NumOfHELIX = NumOfHELIX+1;
                    PDB_struct.HELIX(NumOfHELIX) = struct('name',{'HELIX'},'serNum',{str2num(tline(8:10))},'helixID',{tline(12:14)},'initResName',{tline(16:18)},...
                                                          'initChainID',{tline(20)},'initSeqNum',{str2num(tline(22:25))},'initICode',{tline(26)},'endResName',{tline(28:30)},...
                                                          'endChainID',{tline(32)},'endSeqNum',{str2num(tline(34:37))},'endICode',{tline(38)},'helixClass',{str2num(tline(39:40))},...
                                                          'comment',{tline(41:70)},'length',{str2num(tline(72:76))});                    
  
                %Multiple/Optional                                  
                case 'SHEET'
                    NumOfSHEET = NumOfSHEET+1;
                    PDB_struct.SHEET(NumOfSHEET) = struct('name',{'SHEET'},'strand',{str2num(tline(8:10))},'sheetID',{tline(12:14)},'numStrands',{str2num(tline(15:16))},...
                                                          'initResName',{tline(18:20)},'initChainID',{tline(22)},'initSeqNum',{str2num(tline(23:26))},'initICode',{tline(27)},...
                                                          'endResName',{tline(29:31)},'endChainID',{tline(33)},'endSeqNum',{str2num(tline(34:37))},'endICode',{tline(38)},...
                                                          'sense',{str2num(tline(39:40))},'curAtom',{tline(42:45)},'curResName',{tline(46:48)},'curChainId',{tline(50)},...
                                                          'curResSeq',{str2num(tline(51:54))},'curICode',{tline(55)},'prevAtom',{tline(57:60)},'prevResName',{tline(61:63)},...
                                                          'prevChainId',{tline(65)},'prevResSeq',{str2num(tline(66:69))},'prevICode',{tline(70)});
                    
            
                %Multiple/Optional                                       
                case 'TURN'
                    NumOfTURN = NumOfTURN+1;
                    PDB_struct.TURN(NumOfTURN) = struct('name',{'TURN'},'seq',{str2num(tline(8:10))},'turnId',{tline(12:14)},'initResName',{tline(16:18)},'initSeqNum',{str2num(tline(21:24))},...
                                                        'initICode',{tline(25)},'endResName',{tline(27:29)},'endSeqNum',{str2num(tline(32:35))},'endICode',{tline(36)},'comment',{tline(41:70)});
                                                    
                %Multiple/Optional 
                case 'SSBOND'
                    NumOfSSBOND = NumOfSSBOND+1;
                    PDB_struct.SSBOND(NumOfSSBOND) = struct('name',{'SSBOND'},'serNum',{str2num(tline(8:10))},'resName1',{tline(12:14)},'chainID1',{tline(16)},'seqNum1',{str2num(tline(18:21))},...
                                                            'icode1',{tline(22)},'resName2',{tline(26:28)},'chainID2',{tline(30)},'seqNum2',{str2num(tline(32:35))},...
                                                            'icode2',{tline(36)},'sym1',{tline(60:65)},'sym2',{tline(67:72)});
            
                %Multiple/Optional                                         
                case 'LINK'
                    NumOfLINK = NumOfLINK+1;
                    PDB_struct.LINK(NumOfLINK) = struct('name',{'LINK'},'AtomName1',{tline(13:16)},'altLoc1',{tline(17)},'resName1',{tline(18:20)},'chainID1',{tline(22)},...
                                                        'resSeq1',{str2num(tline(23:26))},'iCode1',{tline(27)},'AtomName2',{tline(43:46)},'altLoc2',{tline(47)},'resName2',{tline(48:50)},...
                                                        'chainID2',{tline(52)},'resSeq2',{str2num(tline(53:56))},'iCode2',{tline(57)},'sym1',{tline(60:65)},'sym2',{tline(67:72)});
                    
            
                %Multiple/Optional                                    
                case 'HYDBND'
                    NumOfHYDBND = NumOfHYDBND+1;
                    PDB_struct.HYDBND(NumOfHYDBND) = struct('name',{'HYDBND'},'AtomName1',{tline(13:16)},'altLoc1',{tline(17)},'resName1',{tline(18:20)},'Chain1',{tline(22)},...
                                                            'resSeq1',{str2num(tline(23:27))},'ICode1',{tline(28)},'nameH',{tline(30:33)},'altLocH',{tline(34)},'ChainH',{tline(36)},...
                                                            'resSeqH',{str2num(tline(37:41))},'iCodeH',{tline(42)},'name2',{tline(44:47)},'altLoc2',{tline(48)},'resName2',{tline(49:51)},...
                                                            'chainID2',{tline(53)},'resSeq2',{str2num(tline(54:58))},'iCode2',{tline(59)},'sym1',{tline(60:65)},'sym2',{tline(67:72)});
                    
                %Multiple/Optional  
                case 'SLTBRG'
                    NumOfSLTBRG = NumOfSLTBRG+1;
                    PDB_struct.SLTBRG(NumOfSLTBRG) = struct('name',{'SLTBRG'},'AtomName1',{tline(13:16)},'altLoc1',{tline(17)},'resName1',{tline(18:20)},'chainID1',{tline(22)},...
                                                            'resSeq1',{str2num(tline(23:26))},'iCode1',{tline(27)},'AtomName2',{tline(43:46)},'altLoc2',{tline(47)},'resName2',{tline(48:50)},...
                                                            'chainID2',{tline(52)},'resSeq2',{str2num(tline(53:56))},'iCode2',{tline(57)},'sym1',{tline(60:65)},'sym2',{tline(67:72)});
            
                %Multiple/Optional                                        
                case 'CISPEP'
                    NumOfCISPEP = NumOfCISPEP+1;
                    PDB_struct.CISPEP(NumOfCISPEP) = struct('name',{'CISPEP'},'serNum',{str2num(tline(8:10))},'ResName1',{tline(12:14)},'chainID1',{tline(16)},'seqNum1',{str2num(tline(18:21))},...
                                                            'icode1',{tline(22)},'ResName2',{tline(26:28)},'chainID2',{tline(30)},'seqNum2',{tline(32:35)},'icode2',{tline(36)},...
                                                            'modNum',{str2num(tline(44:46))},'measure',{str2num(tline(54:59))});
                    
            
                %Multiple/Optional                                        
                case 'SITE'
                    CurSITEName = tline(12:14);
                    
                    if ~strcmp(CurSITEName,PrevSITEName)
                        ResNos = 0;
                        PDB_struct.SITE.NoofSite = PDB_struct.SITE.NoofSite+1;
                        PDB_struct.SITE.SITEDetail(PDB_struct.SITE.NoofSite).SiteName = removeblanks(tline(12:14)); 
                        PDB_struct.SITE.SITEDetail(PDB_struct.SITE.NoofSite).NoOfRes = str2num(tline(16:17));
                        [PDB_struct.SITE.SITEDetail(PDB_struct.SITE.NoofSite).ResDet ResNos] = GetResidueStruct(TmpStruct,tline(19:61),ResNos);
                        PrevSITEName = CurSITEName;
                    else
                        [PDB_struct.SITE.SITEDetail(PDB_struct.SITE.NoofSite).ResDet ResNos] = GetResidueStruct(PDB_struct.SITE.SITEDetail(PDB_struct.SITE.NoofSite).ResDet,...
                                                                                                                   tline(19:61),ResNos);
                    end
                                           
                %Single/Mandatory    
                case 'CRYST1' %Fields in this record: Header(record name), a,b,c (all 3 in Angstrom),alpha,beta,gamma(all 3 in degrees),sGroup,z
                    PDB_struct.CRYST1=struct('name',{'CRYST1'},'a',{str2num(tline(7:15))},'b',{str2num(tline(16:24))},'c',{str2num(tline(25:33))},...
                                             'alpha',{str2num(tline(34:40))},'beta',{str2num(tline(41:47))},'gamma',{str2num(tline(48:54))},...
                                             'sGroup',{tline(56:66)},'z',{str2num(tline(67:70))});
            
                %Single/Mandatory    
                case 'ORIGX' %Fields in this record: Header(record name), O[n][1](O11), O[n][2](O12), O[n][3](O13), T[n](T1)
                               
                    ORIG_name = ['ORIGX',tline(6)];
            
                    PDB_struct.ORIGX(str2num(tline(6))).name = ORIG_name;
                    PDB_struct.ORIGX(str2num(tline(6))).On1 = str2num(tline(11:20));            
                    PDB_struct.ORIGX(str2num(tline(6))).On2 = str2num(tline(21:30)); 
                    PDB_struct.ORIGX(str2num(tline(6))).On3 = str2num(tline(31:40)); 
                    PDB_struct.ORIGX(str2num(tline(6))).Tn = str2num(tline(46:55)); 
            
                %Single/Mandatory       
                case 'SCALE' 
                               
                    SCALE_name = ['SCALE',tline(6)];
            
                    PDB_struct.SCALE(str2num(tline(6))).name = SCALE_name;
                    PDB_struct.SCALE(str2num(tline(6))).Sn1 = str2num(tline(11:20)); 
                    PDB_struct.SCALE(str2num(tline(6))).Sn2 = str2num(tline(21:30)); 
                    PDB_struct.SCALE(str2num(tline(6))).Sn3 = str2num(tline(31:40)); 
                    PDB_struct.SCALE(str2num(tline(6))).Un = str2num(tline(46:55)); 
        
                %Single/Optional: Mandatory if the complete unit must be generated from the given coordinates using non-crystallographic symmetry  
                case 'MTRIX'
                                
                    MTRIX_name = ['MTRIX',tline(6)];
            
                    PDB_struct.MTRIX(str2num(tline(6))).name = MTRIX_name;
                    PDB_struct.MTRIX(str2num(tline(6))).SerNo = str2num(tline(8:10));
                    PDB_struct.MTRIX(str2num(tline(6))).Mn1 = str2num(tline(11:20)); 
                    PDB_struct.MTRIX(str2num(tline(6))).Mn2 = str2num(tline(21:30)); 
                    PDB_struct.MTRIX(str2num(tline(6))).Mn3 = str2num(tline(31:40)); 
                    PDB_struct.MTRIX(str2num(tline(6))).Vn = str2num(tline(46:55)); 
                    PDB_struct.MTRIX(str2num(tline(6))).iGiven = str2num(tline(60));
                   
                %Multiple/Optional    
                case 'TVECT'
                    NumOfTVECT = NumOfTVECT+1;
                    PDB_struct.TVECT(NumOfTVECT) = struct('name',{'TVECT'},'SerNo',{str2num(tline(8:10))},'t1',{str2num(tline(11:20))},'t2',{str2num(tline(21:30))},...
                                                          't3',{str2num(tline(31:40))},'text',{tline(41:70)});
                    
            
                % Group/Optional                                      
                case 'MODEL'
                    NumOfMODEL = NumOfMODEL+1;
                    MODELFlag = 1;
                    PDB_struct.MODEL(NumOfMODEL) = struct('name',{'MODEL'},'MDLSerNo',{str2num(tline(11:14))});
                    
                %Multiple/Optional     
                case 'ATOM'
                    NumOfATOM = NumOfATOM+1;
                    PDB_struct.ATOM(NumOfATOM) = struct('name',{'ATOM'},'AtomSerNo',{str2num(tline(7:11))},'AtomName',{tline(13:16)},'altLoc',{tline(17)},'resName',{tline(18:20)},...
                                                            'chainID',{tline(22)},'resSeq',{str2num(tline(23:26))},'iCode',{tline(27)},'X',{str2num(tline(31:38))},...
                                                            'Y',{str2num(tline(39:46))},'Z',{str2num(tline(47:54))},'occupancy',{str2num(tline(55:60))},'tempFactor',{str2num(tline(61:66))},...
                                                            'segID',{tline(73:76)},'element',{tline(77:78)},'charge',{tline(79:80)});
                    
            
                %Multiple/Optional    
                case 'SIGATM'
                    NumOfSIGATM = NumOfSIGATM+1;
                    PDB_struct.SIGATM(NumOfSIGATM) = struct('name',{'SIGATM'},'AtomSerNo',{str2num(tline(7:11))},'AtomName',{tline(13:16)},'altLoc',{tline(17)},'resName',{tline(18:20)},...
                                                            'chainID',{tline(22)},'resSeq',{str2num(tline(23:26))},'iCode',{tline(27)},'sigX',{str2num(tline(31:38))},...
                                                            'sigY',{str2num(tline(39:46))},'sigZ',{str2num(tline(47:54))},'sigOcc',{str2num(tline(55:60))},'sigTemp',{str2num(tline(61:66))},...
                                                            'segID',{tline(73:76)},'element',{tline(77:78)},'charge',{tline(79:80)});
            
                %Multiple/Optional                                         
                case 'ANISOU'
                    NumOfANISOU = NumOfANISOU+1;
                    PDB_struct.ANISOU(NumOfANISOU) = struct('name',{'ANISOU'},'AtomSerNo',{str2num(tline(7:11))},'AtomName',{tline(13:16)},'altLoc',{tline(17)},'resName',{tline(18:20)},...
                                                            'chainID',{tline(22)},'resSeq',{str2num(tline(23:26))},'iCode',{tline(27)},'U00',{str2num(tline(29:35))},'U11',{str2num(tline(36:42))},...
                                                            'U22',{str2num(tline(43:49))},'U01',{str2num(tline(50:56))},'U02',{str2num(tline(57:63))},'U12',{str2num(tline(64:70))},...
                                                            'segID',{tline(73:76)},'element',{tline(77:78)},'charge',{tline(79:80)});
                    
            
                %Multiple/Optional                                         
                case 'SIGUIJ'
                    NumOfSIGUIJ = NumOfSIGUIJ+1;
                    PDB_struct.SIGUIJ(NumOfSIGUIJ) = struct('name',{'SIGUIJ'},'AtomSerNo',{str2num(tline(7:11))},'AtomName',{tline(13:16)},'altLoc',{tline(17)},'resName',{tline(18:20)},...
                                                            'chainID',{tline(22)},'resSeq',{str2num(tline(23:26))},'iCode',{tline(27)},'SIG11',{str2num(tline(29:35))},'SIG22',{str2num(tline(36:42))},...
                                                            'SIG33',{str2num(tline(43:49))},'SIG12',{str2num(tline(50:56))},'SIG13',{str2num(tline(57:63))},'SIG23',{str2num(tline(64:70))},...
                                                            'segID',{tline(73:76)},'element',{tline(77:78)},'charge',{tline(79:80)});
            
                % Group/Optional                                        
                case 'TER'
                    NumOfTER = NumOfTER + 1;
                    PDB_struct.TER(NumOfTER) = struct('name',{'TER'},'SerialNo',{str2num(tline(7:11))},'resName',{tline(18:20)},'chainID',{tline(22)},'resSeq',{str2num(tline(23:26))},...
                                                      'iCode',{tline(27)});
                    
            
                %Multiple Continued/Optional     
                case 'HETATM'
                    NumOfHETATM = NumOfHETATM+1;
                    PDB_struct.HETATM(NumOfHETATM) = struct('name',{'HETATM'},'AtomSerNo',{str2num(tline(7:11))},'AtomName',{tline(13:16)},'altLoc',{tline(17)},'resName',{tline(18:20)},...
                                                            'chainID',{tline(22)},'resSeq',{str2num(tline(23:26))},'iCode',{tline(27)},'X',{str2num(tline(31:38))},...
                                                            'Y',{str2num(tline(39:46))},'Z',{str2num(tline(47:54))},'occupancy',{str2num(tline(55:60))},'tempFactor',{str2num(tline(61:66))},...
                                                            'segID',{tline(73:76)},'element',{tline(77:78)},'charge',{tline(79:80)});
                           
                % Group/Optional                                        
                case 'ENDMDL'
                    
                    MODELFlag = 0; % reset the MODEL flag
                    NumOfENDMDL = NumOfENDMDL + 1;
                    PDB_struct.ENDMDL(NumOfENDMDL) = struct('name',{'ENDMDL'},'RelMODELNo',{NumOfMODEL});
                                      
                %Multiple/Optional    
                case 'CONECT'
                    
                    NumOfCONECT = NumOfCONECT+1;
                    temp_a = str2num(tline(7:11));
                    temp_b = GetAtomList(tline(12:31));
                    temp_c = GetAtomList([tline(32:41) char(32) tline(47:56)]);
                    temp_d = GetAtomList([tline(42:46) char(32) tline(57:61)]);
                    PDB_struct.CONECT(NumOfCONECT) = struct('name',{'CONECT'},'AtomSerNo',{temp_a},'BondAtomList',{temp_b},'HydAtomList',{temp_c},...
                                                            'SaltBdgAtom',{temp_d});
                    
                    
                %Single/Mandatory   
                case 'MASTER'
                    PDB_struct.MASTER = struct('name',{'MASTER'},'numREMARK',{str2num(tline(11:15))},'numHET',{str2num(tline(21:25))},'numHelix',{str2num(tline(26:30))},...
                                               'numSheet',{str2num(tline(31:35))},'numTurn',{str2num(tline(36:40))},'numSite',{str2num(tline(41:45))},'numXform',{str2num(tline(46:50))},...
                                               'numCoord',{str2num(tline(51:55))},'numTer',{str2num(tline(56:60))},'numConect',{str2num(tline(61:65))},'numSeq',{str2num(tline(66:70))});
               
                %Single/Mandatory       
                case 'END'
                    PDB_struct.END.name = 'END';
                     
                otherwise
                    %disp('The file contains invalid record type');
                    
            end % for the SWITCH statement 
     end % for the IF statement checking the empty string
end % for the WHILE loop


%Initialize all the components of the structure
function PDB_struct = Initialize_struct(PDB_struct)

PDB_struct.TITLE = struct('name',{'TITLE'},'title',{''});
PDB_struct.COMPND = struct('name',{'COMPND'},'comp_description',{''});
PDB_struct.SOURCE = struct('name',{'SOURCE'},'src_description',{''});
PDB_struct.KEYWDS = struct('name',{'KEYWDS'},'KeywdsList',{''});;
PDB_struct.EXPDTA = struct('name',{'EXPDTA'},'technique',{''});
PDB_struct.AUTHOR = struct('name',{'AUTHOR'},'AuthorsList',{''});

PDB_struct.JRNL.Entry = struct('AUTH',{''},'TITL',{''},'EDIT',{''},'REF',{''},'PUBL',{''},'REFN',{''});
PDB_struct.JRNL.name = 'JRNL';
PDB_struct.JRNL.NoOfJRNLS = 0;

PDB_struct.REMARK1.name = 'REMARK1';
PDB_struct.REMARK1.NoOfJRNLS = 0;

PDB_struct.REMARK2 = struct('name',{'REMARK2'},'Detail',{''},'Resolution',{0});
PDB_struct.REMARK3 = struct('name',{'REMARK3'},'Refinement',{''});

PDB_struct.SEQRES.ResChainDetail = struct('NoOfResidue',{},'ChainID',{''},'AminoAcids',{''});
PDB_struct.SEQRES.NoOfResChain = 0;
PDB_struct.SEQRES.name = 'SEQRES';

PDB_struct.HETNAM = struct('name',{'HETNAM'},'hetID',{''},'ChemName',{''});
PDB_struct.HETSYN = struct('name',{'HETSYN'},'hetID',{''},'hetSynonyms',{''});

PDB_struct.FORMUL = struct('name',{'FORMUL'},'CompNo',{0},'hetID',{''},'ChemForm',{''});

ResDetail = struct('ResName',{},'ChainID',{''},'ResSeqNo',{},'InsCode',{''});
PDB_struct.SITE.SITEDetail = struct('SiteName',{''},'NoOfRes',{0},'ResDet',{ResDetail});
PDB_struct.SITE.name = 'SITE';
PDB_struct.SITE.NoofSite = 0;

PDB_struct.CRYST1 = [];
PDB_struct.ORIGX = [];
PDB_struct.SCALE = [];
PDB_struct.MTRIX = [];
PDB_struct.MASTER = [];
PDB_struct.END = [];

% REMOVEBLANKS removes both leading and trailing blanks. This function is written by Steve Simon
function out = removeblanks(in)
[r,c] = find( (in~=0) & ~isspace(in) );
if isempty(c),
    out = in([]);
else
    out = in(:,min(c):max(c));
end

function OutList = GetAtomList(InString)
OutList = [];
str = removeblanks(InString);
while size(str)>0
    [token,rem] = strtok(str);
    token = removeblanks(token);
    rem = removeblanks(rem);
    OutList = [OutList str2num(token)];
    str = rem;
end

function OutAcid = GetAminoAcids(InAcid)
OutAcid = strrep(InAcid,'ALA','A');
OutAcid = strrep(OutAcid,'ARG','R');
OutAcid = strrep(OutAcid,'ASN','N');
OutAcid = strrep(OutAcid,'ASP','D');
OutAcid = strrep(OutAcid,'ASX','B');
OutAcid = strrep(OutAcid,'CYS','C');
OutAcid = strrep(OutAcid,'GLN','Q');
OutAcid = strrep(OutAcid,'GLU','E');
OutAcid = strrep(OutAcid,'GLX','Z');
OutAcid = strrep(OutAcid,'GLY','G');
OutAcid = strrep(OutAcid,'HIS','H');
OutAcid = strrep(OutAcid,'ILE','I');
OutAcid = strrep(OutAcid,'LEU','L');
OutAcid = strrep(OutAcid,'LYS','K');
OutAcid = strrep(OutAcid,'MET','M');
OutAcid = strrep(OutAcid,'PHE','F');
OutAcid = strrep(OutAcid,'PRO','P');
OutAcid = strrep(OutAcid,'SER','S');
OutAcid = strrep(OutAcid,'THR','T');
OutAcid = strrep(OutAcid,'TRP','W');
OutAcid = strrep(OutAcid,'TYR','Y');
OutAcid = strrep(OutAcid,'VAL','V');
OutAcid = strrep(OutAcid,'UNK',' ');

OutAcid = OutAcid(~isspace(OutAcid));


function [OutStruct,OutNum] = GetResidueStruct(TmpStruct,InString,InNum)

a=1; b=10;
Count = 0;
sz = size(InString);

while b <= sz(2)
    test_str = removeblanks(InString(a:b));
    InNum = InNum + 1;   
    while size(test_str)>0
        [token,rem] = strtok(test_str);
        token = removeblanks(token);
        rem = removeblanks(rem);
        Count = Count + 1;
       
        if Count==1
            TmpStruct(InNum).ResName = {token};
        elseif Count==2
            TmpStruct(InNum).ChainID = {token};    
        elseif Count==3
            TmpStruct(InNum).ResSeqNo = {str2num(token)};
        else 
            TmpStruct(InNum).InsCode = {token};
        end
                
        test_str = rem;
    end   
    a=a+11;
    b=b+11;
    Count = 0;
end

OutNum = InNum;
OutStruct = TmpStruct;