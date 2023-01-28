function plotWF(rvals, FGvals)
% Plot the wavefunction components as a function of radial coordinates
plot(rvals, FGvals(1, :), 'r')
hold on
plot(rvals, FGvals(2,:), 'b')
ylabel("Wavefunction")
xlabel("r")
legend("f", "g")
end