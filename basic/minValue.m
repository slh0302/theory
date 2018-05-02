function alph = minValue(f, point, dk)
    numOfvar = length(point);
    var_x = num2cell(sym('x',[1, numOfvar]));
    syms alp;
    f_alph = subs(f,  var_x,  num2cell(point + alp * dk));
    g_alph = jacobian(f_alph);
    alph = double(solve(g_alph, alp));
end