

time_steps = 100;


set_val = 13*ones(time_steps,1);
act_val = zeros(time_steps,1);
act_val(1) = 25;


Kp = 0.1;
Ki = 0.2;
Kd = 0;

integ = zeros(size(act_val));
delta = integ;
f1 = figure;

for i=2:time_steps
    
    err(i) = set_val(i)- act_val(i-1);
    
    integ(i) = integ(i-1) + Ki * err(i);
    
    
    
    delta(i) = Kp * err(i) + Ki * integ(i);
    act_val(i) = act_val(i-1) + delta(i);
    

 
    
    
end

    figure(f1);
    clf
    subplot(2,1,1)
    plot(act_val,'black','DisplayName','actual');
    hold on
    plot(set_val,'blue:','DisplayName','set point');
    subplot(2,1,2);
    plot(err,'g','DisplayName','err'); hold on
    plot(integ,'b','DisplayName','integrator');
    plot(delta,'y','DisplayName','delta applied');
    drawnow
    
    legend