function slopeMat = generateSlope(lfc,tp,group)
%--------------------------------------------------------------------------
%DESCRIPTION:
%            Construct a time slope matrix from time series log2fc data 
%
%INPUT ARGUMENTS:
%lfc             log2fc expression data as matrix or dataframe, each row 
%                represents a gene and each column represents a sample
%tp              timepoint of each sample.
%group           group index of each sample. The same sample/experiment at 
%                different time points should have the same group index. 
%                Repitition experiment should have different group index. 
%                Together with time point, this information is needed to 
%                calculate time slope matrix. 
%
%OUTPUT ARGUMENTS:
%slopeMat        time slope matrix
%--------------------------------------------------------------------------
if nargin~=3
    error('The number of input arguments is incorrect.');
end

if length(group)~=size(lfc,2) || length(tp)~=size(lfc,2)
    error('Dimension of group or timepoint does not match the number of samples in the expression data.');
end

[group,grp_ord] = sort(group);
tp = tp(grp_ord);
lfc = lfc(:,grp_ord);

grp_unique = unique(group);

if length(grp_unique)==length(group)
    error('Slope matrix cannot be calculated because each group has only one sample.');
end

slopeMat = [];

for i=1:length(grp_unique)
    grp_index = find(group==grp_unique(i));
    if length(unique(tp(grp_index)))~=1
        if length(unique(tp(grp_index)))~=length(tp(grp_index))
           error('Replicated samples of some time points in one group. Samples in each group shoud have distinct time points or the same time point.');
        end
        [tp_sorted,tp_order] = sort(tp(grp_index));
        grp_index_order = grp_index(tp_order);
        tp(grp_index) = tp_sorted;
        lfc(:,grp_index) = lfc(:,grp_index_order);
        
        current_tp = tp(grp_index);
        current_lfc = lfc(:,grp_index);
        
        l = length(current_tp);
        if l>2
           for j=1:l
               if j==1
                  slop_sp = findiff(current_tp(1:3),current_lfc(:,1:3),1);
               elseif j==l
                  slop_sp = findiff(current_tp((l-2):l),current_lfc(:,(l-2):l),3);
               else
                  slop_sp = findiff(current_tp((j-1):(j+1)),current_lfc(:,(j-1):(j+1)),2);
               end
               slopeMat = horzcat(slopeMat, slop_sp);
           end
        else
           slop_sp = (current_lfc(:,2)-current_lfc(:,1))/(current_tp(2)-current_tp(1));
           slopeMat = horzcat(slopeMat, slop_sp);
        end
    end
end

end