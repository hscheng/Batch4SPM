% 
% y a m*n variable 
% x cov variable
function resid = multicol_resid(y,x)
    col_num = size(y,2);
    resid = zeros(size(y));
    for ii = 1:col_num
        [~,~,resid_temp] = regress(y(:,ii),x);
        resid(:,ii) = resid_temp;
    end
end