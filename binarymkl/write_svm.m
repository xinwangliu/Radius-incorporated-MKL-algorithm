function write_svm(fname, M, targetmat, issymm, subsample_pos, subsample_neg, fact)

    p = size(M,3);
    m = size(targetmat,1);
    n = size(targetmat,2);
    fid = fopen (fname, 'w');
    assert (fid~=-1);
    totalpos = sum(sum(targetmat==1))-n*issymm;
    totalneg = sum(sum(targetmat==-1));
    if(nargin < 6)
       fact = 1;
    end

    if (subsample_pos ==0 && subsample_neg ==0)
       if (totalpos > totalneg)
	  subsample_neg = 1*fact;
	  subsample_pos = totalneg/totalpos*fact;
       else
	  subsample_pos = 1*fact;
	  subsample_neg = totalpos/totalneg*fact;
       end
    end

    pos = 0; neg = 0;
    for i=1:m
       if (issymm) j_start = i+1;
       else  j_start = 1; end
       for j=j_start:n
	  if ((targetmat(i,j)==1 && subsample_pos > rand) || (targetmat(i,j)==-1 && subsample_neg > rand))
	     str = sprintf ('%d ', targetmat(i,j));
	     for k=1:p
		str = sprintf ('%s%d:%f ', str, k, M(i,j,k));
	     end
	     str = sprintf ('%s%d:%f ', str, p+1, 0.1); %const 1 for bias term
	     fprintf (fid, '%s\n', str);

	     if (targetmat(i,j)==1) pos = pos+1; end
	     if (targetmat(i,j)==-1) neg = neg+1; end
	  end
       end
       fprintf ('%d   \r', i); 
    end
    fclose (fid);

    fprintf ('written: %d +ve, %d -ve K-examples\n', pos, neg);



