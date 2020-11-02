function w = read_model (fname, dim)

    fid = fopen (fname, 'r');
    assert (fid~=-1);

    index = [];
    val = [];
    val = fscanf (fid, '%d:%f');
    w = [];
    for i=1:2:length(val)
       w(round(val(i))) = val(i+1);
    end
    if (length(w) < dim)
       l = length(w);
       w(l+1:dim) = 0;
    end
    w = w(1:dim);
    fclose (fid);


