function [ vec_norm,mean_theta,std_theta ] = vecnorm_calc(EndptList,i2)

%%% Calculate vectors from centroid to endpoints, sum them together and get
%%% vecotr norm to understand directionality


vec = zeros(length(EndptList),2);

% calculate vector norms from centroid to endpoint
for i = 1:size(EndptList,1)
    coords = EndptList(i,:);
    vec(i,:) = coords-i2;
end

vec_norm = norm(sum(vec));

theta_rads=[];
% Loop through all unique row pairs
for i = 1:size(vec, 1) - 1
    for j = i + 1:size(vec, 1)
        u = vec(i, :);  % First row of the pair
        v = vec(j, :);  % Second row of the pair
        CosTheta = max(min(dot(u,v)/(norm(u)*norm(v)),1),-1);;  % get cos_theta
        theta_rads = [theta_rads,acos(CosTheta)];
    end
end

mean_theta = mean(theta_rads);
std_theta = std(theta_rads);
