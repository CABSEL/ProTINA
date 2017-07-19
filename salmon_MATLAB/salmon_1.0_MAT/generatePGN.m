function [pgn,ppn2] = generatePGN(GList,tftg,ppi,ppi_thre,tftg_thre,ptf_thre)

%-------------------------------------------------------------------------
%DESCRIPTION:
%           Generate protein-gene regulatory network from gene regulatory 
%           network and protein protein interaction network
%
%INPUT ARGUMENTS:
%glist           gene ID file. The tfnet and glist file should use the same 
%                system for gene names. The genes in glist should have the 
%                same order as the expression data.
%pdi             gene regulatory network. The first column is the gene names 
%                of transcription factors and the second column is gene being 
%                regulated; the third column is optional, if present, it is 
%                the score for each regulation.
%ppi             protein protein interaction network. The first two columns 
%                are the interacting proteins; the third column is optional, 
%                if present, it is the score for each interaction.
%ppi_thre        threshold for protein-protein interaction. This only has 
%                effect if ppi has a third score column. Interactions with 
%                score lower than the threshold will be ignored.
%tftg_thre       threshold for transcription factor-gene regulation. This 
%                only has effect if pdi has a third score column. Interactions 
%                with score lower than the threshold will be ignored.
%ptf_thre        threshold for protein-transcription factor interaction. This 
%                only has effect if ppi has a third score column. Interactions 
%                with score lower than the threshold will be ignored. 
%
%OUTPUT ARGUMENTS:
%pgn             a sparse matrix representing protein gene regulatory network.
%--------------------------------------------------------------------------

    if nargin<4 || isempty(ppi_thre)
        ppi_thre = 0;
    end
    
    if nargin<5 || isempty(tftg_thre)
        tftg_thre = 0;
    end
    
    if nargin<6 || isempty(ptf_thre)
        ptf_thre = 0;
    end
    
    n = length(GList);
    
    %%%%%%%%%%%%%%%%%%%%%%%% GRN CONSTRUCTION %%%%%%%%%%%%%%%%%%%%%%%%
    %%% Check whether tftg has 3rd column (= score)
    if size(tftg,2)<3
        tftg = [tftg,num2cell(ones(size(tftg,1),1))];
        tftg_thre = 0;
    end
    
    %%%% Trim TF-gene interactions (high score & existing in GList)
    tftg = tftg(cell2mat(tftg(:,3)) >= tftg_thre,:);
    mem1 = ismember(upper(tftg(:,1)),upper(GList));
    mem2 = ismember(upper(tftg(:,2)),upper(GList));
    tftg = tftg(mem1 & mem2,:);

    %%%% GRN construction
    grn = zeros(n);
    tf = unique(tftg(:,1));

    for j=1:length(tf)
        [~,tgi]=ismember(upper(tftg(strcmpi(tf(j),tftg(:,1)),2)),upper(GList));
        grn(strcmpi(tf(j),GList),tgi) = 1;
    end

    grn = sparse(grn);

    %%%% summary %%%%
    tfn = length(tf);
    dgn = nnz(sum(grn~=0,1));
    edgen = nnz(grn);
    fprintf('TF-Gene network (%d TFs, %d genes, %d interactions) has been generated.\n',tfn,dgn,edgen);
    
    %%%%%%%%%%%%%%%%%%% PPI NETWORK CONSTRUCTION %%%%%%%%%%%%%%%%%%%%
    %%% Check whether tftg has 3rd column (= score)
    if size(ppi,2)<3
       ppi = [ppi,num2cell(ones(size(ppi,1),1))];
       ppi_thre = 0;
    end

    %%%% Trim PPIs (high score & existing in GList)
    mem1 = ismember(upper(ppi(:,1)),upper(GList));
    mem2 = ismember(upper(ppi(:,2)),upper(GList));
    ppi = ppi(mem1 & mem2,:);

    % index of ppi
    [~,ic1] = ismember(upper(ppi(:,1)),upper(GList));
    [~,ic2] = ismember(upper(ppi(:,2)),upper(GList));
    ppi_idx = [ic1,ic2,cell2mat(ppi(:,3:end))];

    % Remove duplicated interactions with lower scores
    ppi_idx(:,1:2) = sort(ppi_idx(:,1:2),2);
    ppi_idx = flipud(sortrows(ppi_idx));
    [~,uni] = unique(ppi_idx(:,1:2),'rows','stable');
    ppi_idx = flipud(ppi_idx(uni,:));

    % Remove self-loops
    selfloop = ppi_idx(:,1)-ppi_idx(:,2)==0;
    ppi_idx(selfloop,:) = [];


    %%%% Prot1 - TF network %%%%
    tf_idx = find(sum(grn~=0,2));
    e1 = find(ismember(ppi_idx(:,1),tf_idx) | ismember(ppi_idx(:,2),tf_idx));
    ppi_ly1 = ppi_idx(e1,:);
    proti1 = setdiff(unique(ppi_ly1(:,1:2)),tf_idx);

    % set the order (prot1-TF)
    x = find(ismember(ppi_ly1(:,1),tf_idx));
    ppi_ly1(x,[1,2]) = ppi_ly1(x,[2,1]);

    % If there is a TF1-TF2 interaction, copy as TF2-TF1 as well.
    copyi = ismember(ppi_ly1(:,1),tf_idx) & ismember(ppi_ly1(:,2),tf_idx);
    ppi_ly1 = [ppi_ly1; ppi_ly1(copyi,[2,1,3:end])];

    ppn1 = zeros(n);
    tfi2 = unique(ppi_ly1(:,2));
    for j=1:length(tfi2)
        ind = find(ppi_ly1(:,2)==tfi2(j));
        ppn1(ppi_ly1(ind,1),tfi2(j)) = ppi_ly1(ind,3);
    end
    
    ppn1(ppn1<ptf_thre)=0;
    ppn1 = sparse(ppn1);

    %%%% summary %%%%
    prot1n = nnz(sum(ppn1~=0,2));
    tfn = length(tfi2);
    edgen = nnz(ppn1);
    fprintf('Protein-TF network (%d proteins, %d TFs, %d interactions) has been generated.\n',prot1n,tfn,edgen);
    

    %%%% Prot2 - Prot1 network %%%%
    e2 = find(ismember(ppi_idx(:,1),proti1) | ismember(ppi_idx(:,2),proti1));
    e2 = setdiff(e2,e1);
    ppi_ly2 = ppi_idx(e2,:);

    % set the order (prot2-prot1)
    y = find(ismember(ppi_ly2(:,1),proti1));
    ppi_ly2(y,[1,2]) = ppi_ly2(y,[2,1]);

    % If there is Prot1a-prot1b interactions, copy as Prot1b-Prot1a as well.
    copyi = ismember(ppi_ly2(:,1),proti1) & ismember(ppi_ly2(:,2),proti1);
    ppi_ly2 = [ppi_ly2; ppi_ly2(copyi,[2,1,3])];

    ppn2 = zeros(n);
    proti1b = unique(ppi_ly2(:,2));
    for j=1:length(proti1b)
        ind = find(ppi_ly2(:,2)==proti1b(j));
        ppn2(ppi_ly2(ind,1),proti1b(j)) = ppi_ly2(ind,3);
    end
    
    ppn2(ppn2<ppi_thre)=0;
    ppn2 = sparse(ppn2);

    %%%% summary %%%%
    prot2n = nnz(sum(ppn2~=0,2));
    prot1n = nnz(sum(ppn2~=0,1));
    edgen = nnz(ppn2);
    fprintf('Protein-Protein network (%d upstream proteins, %d proteins, %d interactions) has been generated.\n',prot2n,prot1n,edgen);

    %%
    %%%%%%%%%%%%%%%%%%% PROTEIN-GENE NETWORK CONSTRUCTION %%%%%%%%%%%%%%%%%%%%
    pgn = (ppn2+eye(n))*ppn1;
    pgn = (pgn + eye(n))*grn;

    %%%% summary %%%%
    protn = nnz(sum(pgn~=0,2));
    dgn = nnz(sum(pgn~=0,1));
    edgen = nnz(pgn);
    fprintf('Protein-Gene network (%d proteins and TFs, %d downstream genes, %d interactions) has been generated.\n',protn,dgn,edgen);

    pgn = sparse(pgn);    
end