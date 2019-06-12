% converts F matrix output from operator inference into bilinear H form
% probably could be more efficient but I don't know how much I'm going to
% be using this code so I won't optimize right now (18 April 2018)
function H = F2H(F)

[n,~] = size(F);

[ii,jj,vv] = find(F);
iH = [];
jH = [];
vH = [];

bb = 0;
for i = 1:n     % loop thru chunks of F 
    cc = bb+n+1-i;
    sel = jj>bb & jj<=cc;   % find indices corresponding to said chunk
    itemp = ii(sel);
    jtemp = jj(sel);
    vtemp = vv(sel);
    for j = 1:length(jtemp)
        sj = jtemp(j)-bb;
        if sj==1 % then it's a diagonal term and doesn't need to be split
            iH = [iH; itemp(j)];
            jH = [jH; (i-1)*n+i+sj-1];
            vH = [vH; vtemp(j)];
        else % it's a cross term and needs to be split
            iH = [iH; itemp(j); itemp(j)];
            jH = [jH; (i-1)*n+i+sj-1; (i+sj-2)*n+i];
            vH = [vH; vtemp(j)/2; vtemp(j)/2];
        end
    end
    bb = cc;
end
H = sparse(iH,jH,vH,n,n^2);

% halfF = F/2;
% H = zeros(n,n^2);
% for i = 1:n
%     for j = i:n
%         istart = sum(n:-1:n-i+2);
%         H(:,(i-1)*n+j) = H(:,(i-1)*n+j) + halfF(:,istart+j-i+1);
%         H(:,(j-1)*n+i) = H(:,(j-1)*n+i) + halfF(:,istart+j-i+1);
%     end
% end