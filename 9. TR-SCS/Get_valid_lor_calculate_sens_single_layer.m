function valid_lor = Get_valid_lor_calculate_sens_single_layer(ring_number, crystal_per_ring, min_crystal_diffenrence, max_ring_difference, num_parts, current_part)

% Attention:  num_parts must be exact devision by total_detectors !!!!!
% Can be calculated by total_detectors = sum(1 : det_per_ring) - det_per_ring;

% crystal_pair = formDetectorIndices(ring_number * crystal_per_ring);
crystal_pair = formDetectorIndices_partial(ring_number * crystal_per_ring, num_parts, current_part);
crystal_pair = single(crystal_pair - 1);

data_length = length(crystal_pair);

ring1 = double(zeros(data_length, 1));
ring2 = double(zeros(data_length, 1));
crystal1 = double(zeros(data_length, 1));
crystal2 = double(zeros(data_length, 1));
phi = double(zeros(data_length, 1));
u = double(zeros(data_length, 1));

for i = 1 : data_length
    ring1(i) = floor(crystal_pair(i, 1) / crystal_per_ring);
    ring2(i) = floor(crystal_pair(i, 2) / crystal_per_ring);
    crystal1(i) = mod(crystal_pair(i, 1), crystal_per_ring);
    crystal2(i) = mod(crystal_pair(i, 2), crystal_per_ring);
    phi(i) = floor((mod((crystal1(i) + crystal2(i) + floor(crystal_per_ring / 2)), crystal_per_ring)) / 2);
    if (((crystal1(i) + crystal2(i)) < (floor(3 * crystal_per_ring / 2 ))) && ((crystal1(i) + crystal2(i)) >= floor(crystal_per_ring / 2 )))
        u(i) = abs(crystal1(i) - crystal2(i)) -  floor(crystal_per_ring / 2) + floor((crystal_per_ring - 1) / 2);
    else
        u(i) = -abs(crystal1(i) - crystal2(i)) +  floor(crystal_per_ring / 2) + floor((crystal_per_ring - 1) / 2);
    end
end

% crystal_difference = abs(crystal2 - crystal1);
ring_difference = abs(ring2 - ring1);
% valid_lor_index = (crystal_difference < (crystal_per_ring - min_crystal_diffenrence)) & (crystal_difference > min_crystal_diffenrence) & (ring_difference < max_ring_difference);

valid_lor_index = (u <= (crystal_per_ring - min_crystal_diffenrence)) & (u >= min_crystal_diffenrence) & (ring_difference <= max_ring_difference);
valid_lor = crystal_pair(valid_lor_index, :) + 1;



end



function L_part = formDetectorIndices_partial(det_per_ring, num_parts, current_part)

% Attention:  num_parts must be exact devision by total_detectors !!!!!
% Can be calculated by total_detectors = sum(1 : det_per_ring) - det_per_ring;

    total_detectors = sum(1 : det_per_ring) - det_per_ring;
    detectors_per_part = ceil(total_detectors / num_parts);
    start_idx = (current_part - 1) * detectors_per_part + 1;
    end_idx = min(current_part * detectors_per_part, total_detectors);

    L_part = zeros(end_idx - start_idx + 1, 2, 'uint32');

    % Preallocate the output array
    % current_indices = zeros(det_per_ring * (det_per_ring - 1) / 2, 2, 'uint32');
    % current_indices = sparse(det_per_ring * (det_per_ring - 1) / 2, 2);
    current_indices = zeros(floor(detectors_per_part * 1.5), 2, 'uint32');

    idx = zeros(det_per_ring, 1);
    idx(1) = 1;
    for i = 1 : det_per_ring
        idx(i + 1) = idx(i) + det_per_ring - i;
    end
    clear i

    lower_threshold = find(idx <= start_idx);
    upper_threshold = find(idx >= end_idx);
    start_crystal = lower_threshold(end);
    end_crystal = upper_threshold(1);

    % Generate all indices at once using vectorized operations
    current_idx = 1;
    for i = start_crystal : end_crystal
        current_indices(current_idx : current_idx + det_per_ring - i - 1, :) = [repelem(i, det_per_ring - i)', (i + 1 : det_per_ring)'];
        current_idx = current_idx + det_per_ring - i;
    end

    % Copy the relevant indices to the output array
    L_part(1 : end_idx - start_idx + 1, :) = current_indices((start_idx - idx(start_crystal) + 1) : (start_idx - idx(start_crystal) + detectors_per_part), :);

end


function L = formDetectorIndices(det_per_ring)
%FORMDETECTORINDICES Outputs the detector index for the specified detector
%numbers
% This is a helper function
L = zeros(sum(1 : det_per_ring) - det_per_ring, 2, 'int32');
jh = int32(1);
for kk = int32(1) : (det_per_ring)
    if exist('OCTAVE_VERSION','builtin') == 0 && exist('repelem', 'builtin') == 0
        L(jh:(jh + (det_per_ring) - kk),:) = [repeat_elem((kk), det_per_ring-(kk-1)), ((kk):det_per_ring)'];
    elseif exist('OCTAVE_VERSION','builtin') == 5
        L(jh:(jh + (det_per_ring) - kk),:) = [repelem((kk), det_per_ring-(kk-1)), ((kk):det_per_ring)'];
    else
        L(jh:(jh + (det_per_ring) - kk -1),:) = [repelem((kk), det_per_ring-(kk))', ((kk + 1) : det_per_ring)'];
    end
    jh = jh + (det_per_ring) -kk;
end
L(L(:,1) == 0,:) = [];
end

