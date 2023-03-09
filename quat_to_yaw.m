% input is quaternion q=[w,x,y,z], suppport vector
function yaw=quat_to_yaw(w,x,y,z)
dcm10 = 2 * (x .* y + w .* z);
dcm00 = w.*w + x.*x - y.*y - z.*z;
yaw = atan2(dcm10, dcm00);
end