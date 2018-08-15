function amp = ampdist(spk_times)

counter = 0;

while ~isempty(spk_times)
    counter = counter + 1;
    spk = spk_times(1);
    eq_inds = (spk_times == spk);
    amp(counter) = sum(eq_inds);
    spk_times(eq_inds) = [];
end

% disp('finish')