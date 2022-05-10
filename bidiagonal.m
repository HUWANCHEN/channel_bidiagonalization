function [Q,B,P]=bidiagonal(D)
%************************
% D = QBP';
%************************
[m,n]=size(D);
P = eye(n);
Q = eye(m);
if m<n
    for j = m:-1:1
        %先更新行 使当前行的最后一个元素为非零 右乘酉矩阵 更新行
        colInd = n+j-m; %即将要更新的列的index
        x = D(j,1:colInd)';
        [v,b,ephi] = house(flipud(x));
        vv = flip(v); %满足要求的householder向量
        D(1:j,1:colInd) = ephi'*(D(1:j,1:colInd)-b*D(1:j,1:colInd)*vv*vv');
        P(:,1:colInd) = ephi'*(P(:,1:colInd)-b*P(:,1:colInd)*vv*vv');
        %再更新列 使当前列的最后一个元素为非零 左乘酉矩阵 更新列
        if j>1 %如果更新到了第一行 则不需要更新列
            y = D(1:j-1,colInd);
            [v,b,ephi] = house(flipud(y));
            vv = flip(v);
            D(1:j-1,1:colInd) = ephi*(D(1:j-1,1:colInd)-b*vv*vv'*D(1:j-1,1:colInd));
            Q(1:j-1,:) = ephi*(Q(1:j-1,:)-b*vv*vv'*Q(1:j-1,:));
        end
    end
else
    for j = 1:n
        [v,b,ephi] = house(D(j:m,j));
        D(j:m,j:n) = ephi*(D(j:m,j:n) - b*v*v'*D(j:m,j:n));
        Q(j:m,:) = ephi*(Q(j:m,:) - b*v*v'*Q(j:m,:));
        if j< n
            [v,b,ephi] = house(D(j,j+1:n)');
            D(j:m,j+1:n) = ephi'*(D(j:m,j+1:n) - b*D(j:m,j+1:n)*v*v');
            P(:,j+1:n) = ephi'*(P(:,j+1:n) - b*P(:,j+1:n)*v*v');
        end
    end
    
end
Q = Q';
B=D;
end