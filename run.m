d = nmpc();

% figure(1)
% for i = d.p.nodenames
%     plot(squeeze(sum(d.s.dij(i, :, :), 2)))
%     hold on
% end
% hold off
% fig = gcf;
% title('Evolution of dij')
% xlabel('Time (min)')
% ylabel('Customers (PAX)')
% 
% figure(2)
% veh = randsample(d.p.veh, 5);
% for i = 1:length(veh)
%     plot(d.s.qk(veh(i), :))
%     hold on
% end
% hold off
% fig = gcf;
% title('Battery of 5 random vehicles')
% xlabel('Time (min)')
% ylabel('Battery (%)')