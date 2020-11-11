lib=fastaread('Outputs/localkeep_New_Proteins.fasta');
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

