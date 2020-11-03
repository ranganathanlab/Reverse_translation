lib=fastaread('Inputs/localkeep_New_Proteins.fasta');
localgene = load('Inputs/local_gene.mat').gene;
hits=blastreadlocal('Inputs/blast_local2.txt', 8);

for i=1:numel(lib)
seq=lib(i).Sequence; 
seq(seq=='-')=[]; 
lib(i).protein=seq; 
lib(i).len=numel(seq);
end
% now, retrieve the sequences from ncbi % load('rid.mat');

keep=zeros(size(lib)); 
for i=1:numel(hits)
    percent=hits(i).Hits(1).HSPs(1).Identities(1).Percent; 
    if percent==100
        keep(i)=1;
    end
    [mper, ind]=max(percent); 
    if mper ==100
        hits(i).Hits(1:ind-1)=[]; keep(i)=1;
        disp([num2str(i) ' Corrected']);
    else
        disp([num2str(i) ' Fail']);
    end
end

per=zeros(1,1); 
%for j=1:1
j=1;
per(j)=hits(i).Hits(j).HSPs(1).Identities(1).Percent;
    for i=1:numel(lib) 
        region=hits(i).Hits(1).HSPs(1).SubjectIndices; 
        if region(2) < region(1) % reverse frame
            region=sort(region);
            sequence{i}=seqrcomplement(localgene(i,region(1):region(2)));
        else
            sequence{i}=localgene(i,region(1):region(2));
        end
        
        if i==1
            sequence{i} = 'aatttcatttacaaggcaaaagcactgtacccatacgacgctgatgacgatgatgcttacgaaatctcatttgaacaaaatgaaatcctacaagtctctgacattgaaggcagatggtggaaggcaagaagggcaaacggtgaaacaggtattattccaagcaattatgttcaactaatcgatggt';
        end
        
        prot{i}=nt2aa(sequence{i},'AlternativeStartCodons','False','ACGTOnly','False');
        [a1,a2,a3]=nwalign(lib(i).Sequence, prot{i});
        lib(i).a2=a2;
        %align(i)=a2(:,:);
        match(i)=sum(a2(2,:)=='|')-lib(i).len; 
        disp(i);
    end
% turns out, some 16 hits are not good matches. i'll remove them from lib.
%lib(match<0)=[]; 
%hits(match<0)=[]; 
 %genbank(match<0)=[]; 
%prot(match<0)=[]; 
%sequence(match<0)=[]; 
for i=1:numel(lib)
    lib(i).dna=sequence{i};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%save('lib_local.mat','lib')
% correct degenerate nucleotides in the list
%lib(32).dna(45)='G';
%lib(81).dna(42)='A';
lib1=lib;
save('Outputs/lib_local.mat','lib')

sc_codon_usage = load('s_cerevisiae_codon_usage.mat').sc_codon_usage;
codon_freq = sc_codon_usage.codon_frequencies;
aa=sc_codon_usage.aa;
codon_table = sc_codon_usage.aa_codon;
resites={'GGATCC','GAATTC','AGCGCT','GCAGTG','CACTGC','AAAAA','GGGGG','CCCC C','TTTTT'}; 

%for i=1:numel(lib) 
%    lib(i).oligo_assemble=resite_remove(lib(i).dna,resites); 
%    if mod(i,1000)==0
%        disp(num2str(i));
%    end
%end

function seq=resite_remove(dna, resites)
% SEQ = RESITE_remove(DNA, RESITES) accepts a DNA sequence that is in
% frame, and finds restriction sites found in the cell array RESITES. If
% none are found, the dna sequence is returned unchanged to the variable
% SEQ. If RESITES are found, then synonymous substitutions are made to dna % to remove RESITES, and this edited dna sequence is returned to SEQ
% subramaniansk@gmail.com
% codons are picked from yeast codon usage. 
sc_codon_usage=load('s_cerevisiae_codon_usage.mat').sc_codon_usage;
codon_freq = sc_codon_usage.codon_frequencies;
aa=sc_codon_usage.aa;
codon_table = sc_codon_usage.aa_codon;
seq=upper(dna);
resites=upper(resites);
pro=nt2aa(seq,'AlternativeStartCodons','False');
re_flag =0;
while re_flag < 1
for i=1:numel(resites) loc=strfind(seq, resites{i}); if numel(loc)>=1
amino acid
for j=1:numel(loc)
pos=ceil(loc(j)/3);
% pick pos or pos+1 as the codon to change if (pos < numel(pro)-1) && (rand()>0.6)
pos=pos+1;
end
cod=codon_table{aa==pro(pos)};% get all codons for the
% pick a codon randomly, based on its frequency.
cod_freq=codon_freq(aa==pro(pos),:); cod_freq(cod_freq==0)=[];
cod_freq=cumsum(cod_freq);cod_freq=cod_freq/cod_freq(end); pick_cod=cod{find(rand() < cod_freq,1,'first')};
% switch in the new codon
seq(pos*3-2:pos*3)=pick_cod;
end
end
end
% we've gone through the sequence once and swiched codons where ever % we saw RE sites. this might not have removed all RE sites because % (a) the random switch might put back the WT codon i.e. create no
% change
% (b) the switch might itself create a new RE site.
% So, check for the presense of re sites again and set re_flag. re_flag =10;
for i=1:numel(resites)
if strfind(seq,resites{i}) re_flag=0;
end
end
end