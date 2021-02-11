% Function to calculate the Analytical Solution
function [analy_temp] = analytical(time,JM)

syms x
str = input('Enter IC fucntion for analytical analysis in x: ','s') ;
f = function_handle.empty; 
f = eval(['@(x)', str]);

xl = linspace(0,1,JM);
analy_temp = zeros(1,JM);

for i = 1:1:JM
        for n = 1:2:50
            if (n == 1)
                BnFunc = (f * sin(n*pi*x));
                Bn = eval(2*int(BnFunc, 0, 1));
                analy_temp_sum = Bn * sin(n*pi*xl(i)) * exp((-(n*pi)^2) * time);
            else
                 BnFunc = (-(x^2) + x) * sin(n*pi*x);
                 Bn = eval(2*int(BnFunc, 0, 1));
                 analy_temp_sum =  analy_temp_sum + (Bn * sin(n*pi*xl(i)) * exp((-(n*pi)^2) * time));
            end
        end
        analy_temp(i) = analy_temp_sum;
end

end