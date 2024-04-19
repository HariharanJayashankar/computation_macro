%%%%%%%%%%%%%%%%%%%%%%%%%%
% some graphs for testing
%%%%%%%%%%%%%%%%%%%%%%%%%%
[pmesh, amesh] = meshgrid(p.pgrid, p.agrid);
figure
surf(pmesh, amesh, V','EdgeColor','none')
xlabel("Price")
ylabel('$a_{it}$', 'Interpreter', 'latex')
title("Value Func")
saveas(gcf, "figures/valuefunc.png")

figure
surf(pmesh, amesh, polp','EdgeColor','none')
xlabel("Price")
ylabel('$a_{it}$', 'Interpreter', 'latex')
title("Policy for Price - Conditional on Changing")
saveas(gcf, "figures/polp.png")

pflexprice = flexsoln(w, p);
colormap(uint8([250, 192, 213;255, 255, 255]))
imagesc(p.agrid, log(p.pgrid), pollamb)
set(gca,'YDir','normal')
hold on;
plot(p.agrid, log(pflexprice), 'r')
ylabel("$\log(p_{it})$", 'Interpreter', 'latex')
xlabel('$a_{it}$', 'Interpreter', 'latex')
title("Policy for Price Change")
colorbar;
saveas(gcf, "figures/pollamb.png")

figure
plot(p.pgrid, pdist);
hold on
plot(p.pgrid, pdistalt, 'o');
legend("omega", "omegahat")
title("Marginal distribution of prices")
saveas(gcf, "figures/pdist.png")

figure
colormap(flipud(winter))
imagesc(p.agrid, log(p.pgrid), omega)
set(gca,'YDir','normal')
ylabel("$\log(p_{it})$", 'Interpreter', 'latex')
xlabel('$a_{it}$', 'Interpreter', 'latex')
title("p, a Distribution")
colorbar;
saveas(gcf, "figures/omega.png")

