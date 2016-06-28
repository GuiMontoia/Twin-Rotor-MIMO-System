incremento_u_11 = 0.1;
incremento_u_22 = 0.1;
u_11 = double(u1_ini);
u_22 = double(u2_ini);
cond_ini=double([ x1 x2 x3 x4 x5 x6 x7]');
% cond_ini=double([0.698 0 0.785 0 2.23  -4.95 -0.445]');
% cond_ini = x;
t_sim = 300;
%%
% set_param('modelo_completo','MaxConsecutiveZCsMsg','none');
sim modelo_completo;

subplot(211)
plot(salida_no_lineal1.Time,salida_no_lineal1.Data, ...
    salida_lineal1.Time,salida_lineal1.Data+cond_ini(1),'-.');
title('Cabeceo')
legend('Modelo no lineal','Modelo linealizado')
xlabel('tiempo (s)')
ylabel('angulo (rad)')

% % plot(x11.Time,x11.Data);
% % plot(x12.Time,x12.Data+cond_ini(1));
subplot(212)
plot(salida_no_lineal3.Time,salida_no_lineal3.Data, ...
    salida_lineal3.Time,salida_lineal3.Data+cond_ini(3),'-.');
% % ylim([0.5,0.8])
legend('Modelo no lineal','Modelo linealizado')
title('Guiñada')
ylabel('angulo (rad)')
xlabel('tiempo (s)')

% plot(x11.Time,x11.Data)
% ylim([0.5,0.8])
% plot(x32.Time,x32.Data+cond_ini(3))
