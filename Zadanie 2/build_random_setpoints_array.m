function setPoints = build_random_setpoints_array(workpoint, len, interval, min, max)
    setPoints = workpoint.y*ones(1, len);
    for i = interval:interval:len
        j = randi([1, 2], 1);
        setPoints(j, i:end) = min(j) + (max(j) - min(j))*rand();
    end
end

