function max_n = maxMag(vetor, n)
    for i=1:1:n
        vetor(i) = round((abs(vetor(i))),15);
    end
    max_n = max(abs(vetor));
end
