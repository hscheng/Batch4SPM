function [resid,z_resid] = new_resid(y,x)
    resid = zeros(size(y));
    [dim_x,dim_y] = size(y);
    for ii = 1:dim_y
        paras = glmfit(x,y(:,ii));
        t1 = repmat(paras(2:end)',dim_x,1).*x;
        t2 = (sum(t1'))';
        t3 = t2+paras(1);
        resid(:,ii) = y(:,ii)-t3;
    end
	z_resid = (resid - repmat(mean(resid),dim_x,1))./repmat(std(resid),dim_x,1);
end