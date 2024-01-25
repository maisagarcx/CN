function A = find_A_from_L_U(L,U,pivot)
    PT = find_transp(pivot);
    A=PT*L*U;
end
