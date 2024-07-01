function [t_new,xyz_new]=add_virtual_fifo_topic_data(cur_dataset)
t = cur_dataset.timestamp_sample;
dt = cur_dataset.dt;
samples = cur_dataset.samples;
scale = cur_dataset.scale;
total_samples = 0;

N=length(t);
for i=1:N
    total_samples=total_samples + samples(i);

end
t_new = zeros(total_samples,1);
xyz_new = zeros(total_samples, 3);
sample = 0;


X=cur_dataset{:,6:25};
Y=cur_dataset{:,38:57};
Z=cur_dataset{:,70:89};
for i=1:N
    start=samples(i)*(i-1)+1;
    last=samples(i)*i;
    % x=cur_dataset{i,6:6+samples(i)-1}*scale(i);
    x=X(i,1:samples(i))*scale(i);
    xyz_new(start:last,1)=x';

    % y=cur_dataset{i,38:38+samples(i)-1}*scale(i);
    y=Y(i,1:samples(i))*scale(i);
    xyz_new(start:last,2)=y';

    % z=cur_dataset{i,70:70+samples(i)-1}*scale(i);
    z=Z(i,1:samples(i))*scale(i);
    xyz_new(start:last,3)=z';

    for s= 1:samples(i)
        t_new(sample+s) = t(i)-(samples(i)-s)*dt(i);
    end
    sample=sample + samples(i);

end

% t_new=t;
% xyz_new=0;
end