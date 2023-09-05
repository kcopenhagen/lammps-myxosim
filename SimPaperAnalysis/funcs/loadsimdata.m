function bds = loadsimdata(fname)
    A = readmatrix(fname,"NumHeaderLines",9,"FileType","text");
    
    id = num2cell(A(:,1));
    mol = num2cell(A(:,2));
    xs = num2cell(A(:,3));
    ys = num2cell(A(:,4));
    zs = num2cell(A(:,5));
    mux = num2cell(A(:,6));
    muy = num2cell(A(:,7));
    muz = num2cell(A(:,8));
    vx = num2cell(A(:,9));
    vy = num2cell(A(:,10));
    vz = num2cell(A(:,11));
    af = num2cell(A(:,12));
    bds = struct('id',id,'mol',mol,'xs',xs,'ys',ys,'zs',zs, ...
        'mux',mux,'muy',muy,'muz',muz,'vx',vx,'vy',vy,'vz',vz,'af',af);
    
end