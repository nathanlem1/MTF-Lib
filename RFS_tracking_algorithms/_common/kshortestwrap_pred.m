function [paths,costs]= kshortestwrap_pred(rs,k)

if k==0
    paths= [];
    costs= [];
    return;
end
  
  ns= length(rs); 
  [ds,is]= sort(rs(:),1,'descend'); %sort weights in decreasing order

  CM= 0*ones(ns,ns); %cost matrix for paths !REMEMBER ZERO COST DENOTES NO ARC CONNECTION I.E. INF COST BECAUSE OF SPARSE MATRIX INPUT FOR KSHORTEST PATH
  
  for i=1:ns
     CM(1:i-1,i)= ds(i); %only allow jumps to higher numbered nodes, inf costs (equiv to zero cost for sparse representation) on lower diag prohibit reverse jumps (hence cycles)
  end
   
  CMPad= 0*ones(ns+2,ns+2); %extra 2 states for start and finish points !REMEMBER ZERO COST DENOTES NO ARC CONNECTION I.E. INF COST BECAUSE OF SPARSE MATRIX INPUT FOR KSHORTEST PATH
  CMPad(1,2:end-1)= ds'; %must enter into one of original nodes OR
  CMPad(1,end)= eps(0); %exit immediately indicating no node selection (all target die) (eps used to denote zero cost or free jumo, as zero is used for inf in sparse format)
  CMPad(2:end-1,end)= eps(0); %must exit at last node at no cost (eps used to denote zero cost or free jump, as zero is used for inf in sparse format)
  CMPad(2:end-1,2:end-1)= CM; %cost for original nodes
  
  [paths,costs]= kShortestPath_any(CMPad,1,ns+2,k); %do k-shortest path
  
  for p=1:length(paths)
      if isequal(paths{p},[1 ns+2])
          paths{p}= [];  %0 indicates no nodes selected
      else
          paths{p}= paths{p}(2:end-1)-1; %strip dummy entry and finish nodes
          paths{p}= is(paths{p}); %convert index back to unsorted input
      end
  end
 
 