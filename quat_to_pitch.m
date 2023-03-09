% input is quaternion q=[w,x,y,z], suppport vector
function pitch=quat_to_pitch(w,x,y,z)
dcm20 = 2 * (x .* z - w .* y);
pitch = asin(-dcm20);
end