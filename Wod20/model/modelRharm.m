function [FId,Ld,Lh,Kd] = modelRharm(t,om) %Nrom sigZf,Lhm,Kd,T,
%% Projekt modelu - oblicz FId
% za³. funkcji harmonicznej
% Nrom = [1 2 5 20 36 42 134 236 500 600];
k=0; Ld=length(t); Lh=length(om); 
for (nh = 1:Lh)
    k = k + 1; FId(:, k) = sin(om(nh) * t'); k = k + 1; FId(:, k) = cos(om(nh) * t');
end

Kd = k;

end