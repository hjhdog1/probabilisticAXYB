function W = q2W(q)
    W = [q(1) -q(2) -q(3) -q(4);
        q(2) q(1) q(4) -q(3);
        q(3) -q(4) q(1) q(2);
        q(4) q(3) -q(2) q(1)];

end

