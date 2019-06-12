function X2 = get_x_sq(X)
w = size(X,2);
Xi = cell(w,1);
for j = 1:w
    Xi{j} = repmat(X(:,j),1,w-j+1).*X(:,j:w);
end
X2 = cat(2,Xi{:});
