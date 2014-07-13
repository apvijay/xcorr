function [z z_val z_idx]  = xcorrn(x,y)
% Calculates the n-dimensional cross correlation between two matrices.
% Works only up to 3d now. Matrices are non-negative.
% 
% Inputs:
% x, y = non-negative n-dim matrices
% 
% Outputs:
% z = correlation matrix
% z_val = max value of z
% z_idx = index where z_val occurs
%

siz_x = size(x);
siz_y = size(y);

if sum(abs(siz_x - siz_y)) ~= 0
    fprintf('Input dimensions are not equal\n');
end

% Make the sizes 1 for 2nd and 3rd dimensions if needed
if numel(siz_x) == 1
    siz_x = [siz_x 1 1];
%     siz_y = [siz_y 1 1];
end

if numel(siz_x) == 2
    siz_x = [siz_x 1];
%     siz_y = [siz_y 1];
end

% For every shift, find element-by-element by multiplication
z = zeros(2*siz_x - 1);
for i1 = -siz_x(1):siz_x(1)
    for i2 = -siz_x(2):siz_x(2)
        for i3 = -siz_x(3):siz_x(3)
            z(i1+siz_x(1)+1,i2+siz_x(2)+1,i3+siz_x(3)+1) = ...
                sum(reshape(x .* shift_nd(y,i1,i2,i3),[],1));
        end
    end
end

% Normalise the correlation matrix
z = z / sqrt(sum(x(:).^2) * sum(y(:).^2));
[z_val z_idx] = max(z(:)); % max value and index (in 1d) 
[z_idx1 z_idx2 z_idx3] = ind2sub(size(z), z_idx); % index in 3d
z_idx = [z_idx1 z_idx2 z_idx3]';
end

function x_s = shift_nd(x,i1,i2,i3)
% Shifting a matrix in each dimension
% Works only for 2d, 3d matrices

siz_x = size(x);
siz_1 = siz_x(1);

if numel(siz_x) > 1
    siz_2 = siz_x(2);
else
    siz_2 = 1;
end
if numel(siz_x) > 2
    siz_3 = siz_x(3);
else
    siz_3 = 1;
end

x_s = x;

if siz_1 > 1
    if i1 > 0 % horizontal right shift
        x_s = [zeros(siz_1,i1,siz_3) x_s(:,1:siz_2-i1,:)];
    elseif i1 < 0 % horizontal left shift
        i1 = -i1;
        x_s = [x_s(:,i1+1:siz_2,:) zeros(siz_1,i1,siz_3)];
    else % no shift
%         x_s = x_s;
    end
end

if siz_2 > 1
    if i2 > 0 % vertical bottom shift
        x_s = [zeros(i2,siz_2,siz_3);x_s(1:siz_1-i2,:,:)];
    elseif i2 < 0 % vetical top shift
        i2 = -i2;
        x_s = [x_s(i2+1:siz_1,:,:);zeros(i2,siz_2,siz_3)];
    else % no shift
%         x_s = x_s;
    end
end

if siz_3 > 1
    if i3 > 0 % go-in shift
        x_s(:,:,i3+1:end) = x_s(:,:,1:siz_3-i3);
        x_s(:,:,1:i3) = zeros(siz_1,siz_2,i3);
    elseif i3 < 0 % come-out shift
        i3 = -i3;
        x_s(:,:,1:siz_3-i3) = x_s(:,:,i3+1:siz_3);
        x_s(:,:,siz_3-i3+1:end) = zeros(siz_1,siz_2,i3);
    else % no shift
%         x_s = x_s;
    end
end
end