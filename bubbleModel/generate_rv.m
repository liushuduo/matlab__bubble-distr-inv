function X=generate_rv(N,x,p)

% Generates and rv according to the pdf described by x and p

c=cumtrapz(x,p);    % computes the cdf

c=c/c(end);         % Normalises, so that p does not have to be a normalised pdf, i.e. it can be improper

Y=rand(N,1);
X=Y;                % Also set X equal to zero

for k=1:N
    ind=find(c<Y(k),1,'last');
    %disp([c(ind) c(ind+1) Y(k)])
    if (ind==N)     % Covering the case where the values is at the end of the array
        X(k)=x(ind);
        continue
    end
    if (c(ind)~=c(ind+1))
        X(k)=(Y(k)-c(ind))*(x(ind+1)-x(ind))/(c(ind+1)-c(ind))+x(ind);  % Linear interpolation
    else            % If c deos not change, usually means c==1
        X(k)=x(ind);
    end
end
