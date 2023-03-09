% input is quaternion q=[w,x,y,z], suppport vector
function roll=quat_to_roll(w,x,y,z)
dcm21 = 2 * (w .* x + y .* z);
dcm22 = w.*w - x.*x - y.*y + z.*z;
roll = atan2(dcm21, dcm22);
end