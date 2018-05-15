clear;


sensing_data = csvread('day_data.csv');

keySet = {};
valueSet = {};

map = containers.Map;
r=2;
cluster_data_count=20;

vals = (sensing_data(r+1, 1:cluster_data_count)*10)-(sensing_data(r, 1:cluster_data_count)*10);
ids = 1:20;

function [map] = compressionLZW (ids, vals)

for i=1:1:cluster_data_count
    if vals(i) ~= 0
        if map.isKey(vals(i)) ~= 1 % 없을 경우
            keySet{end+1}=vals(i);
            valueSet = [valueSet sprintf('%d',ids(i))];
            map = containers.Map(keySet, valueSet);
        else
            map(vals(i)) = sprintf('%s,%d',map(vals(i)), ids(i));
            keySet = map.keys;
            valueSet = map.values;
            map = containers.Map(map.keys, valueSet);
        end
    end
end

keySet = map.keys;
for i=1:1:length(map)
    if i ~= 1
        keySet{i} = keySet{i} - keySet{1};
    end
end

map = containers.Map(keySet, valueSet);

end

