function d = ISSF_CBF_disturbance(x)
            d = 10 * [2 * sin(x(1)); 3 * sin(0.5*x(2))]+x(1)+x(2);
end