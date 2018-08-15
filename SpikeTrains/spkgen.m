% =========================================================================
% Spike generator (Humphries)
% =========================================================================

% This function has changed to consider the correlation deviation and also
% the frequency modulation function.

function Spk = spkgen(t, N, r, alpha)
% N = 100;      % number of spkie trains
% r = 50        % firing rate (Hz)
% alpha = 0.5;  % correlation

t = t./1000; % to use the correct units!
dt = t(2)-t(1);
T = zeros(1, N);
p = r * N * dt;
U = rand(size(t));
Spk = U<=p;
% Spk = zeros(length(t),1);
% for i = 1:length(t)
%     p=r(i)*dt;  % Application of frequency modulation function
%     U = rand;
%     R = (U <= p); % the reference spike train
%     
%     for j = 1:N
%         if rand <= alpha
%             T(j) = R;
%         else
%             T(j) = (rand <= p);
%         end
%     end
%     Spk(i) = sum(T);
% end
