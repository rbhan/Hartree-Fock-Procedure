function hartree_fock

% 1. Get the nuclear repulsion energy
    Enuc = 8.90770810;
    
% 2. Use Psi to compute the one-electron integrals
    S = xlsread('han_hf_overlap.xlsx'); % overlap integrals
    T = xlsread('han_hf_kinetic.xlsx'); % kinetic energy integrals
    V = xlsread('han_hf_potential.xlsx'); % potential energy integrals
    h = T + V; % core Hamiltonian (1-e integrals)
    nbf = size(S,1); % each row is an atomic orbital
  % rewrite 2 electron integrals as 4 index tensor
    two_e = xlsread('han_hf_2eints.xlsx'); % two electron integrals
    Int = zeros(nbf, nbf, nbf, nbf); % 4 index tensor
    for a = 1:length(two_e)
        p = two_e(a,1)+1;
        q = two_e(a,2)+1;
        r = two_e(a,3)+1;
        s = two_e(a,4)+1;
        term = two_e(a,5);
    %[pq|rs] =[qp|rs] =[pq|sr] =[qp|sr] =[rs|pq] =[sr|pq] =[rs|qp] =[sr|qp]
        Int(p,q,r,s) = term;
        Int(q,p,r,s) = term;
        Int(p,q,s,r) = term;
        Int(q,p,s,r) = term;
        Int(r,s,p,q) = term;
        Int(s,r,p,q) = term;
        Int(r,s,q,p) = term;
        Int(s,r,q,p) = term;        
    end
    
% 3. Construct the orthogonalizing matrix invsqS
    % [V,D] = eig(A); diagonal matrix D = eVals and matrix V = right eVecs, s.t. A*V = V*D.
    [U, Lambda] = eig(S); % diagonalize S matrix
    % sort eVals and corresponding eVecs in descending order
    U = fliplr(U);
    Lambda = flipud(fliplr(Lambda));
        %S1 = U*Lambda*U';
    invsqLambda = Lambda;
    for j = 1:length(Lambda)
        invsqLambda(j,j) = Lambda(j,j)^(-0.5);
    end
    invsqS = U*invsqLambda*U';
    
% 4. Construct an initial (guess) density matrix
    D0 = zeros(nbf, nbf);
    F0_ortho = invsqS'*h*invsqS; % core Fock matrix in orthogonalized basis
    [C0_ortho, eps0] = eig(F0_ortho); % diagonalize orthogonalized Fock
    % sort eVals/eVecs in ascending order
    [~,perm0] = sort(diag(eps0));
    eps0 = eps0(perm0,perm0);
    C0_ortho = C0_ortho(:,perm0);
    C0 = invsqS*C0_ortho; % initial SCF eVec matrix in original basis
    for m = 1:nbf % rows
        for n = 1:nbf
            for i = 1:5
                D0(m,n) = D0(m,n) + C0(m,i)*C0(n,i);
            end
        end
    end
    E0 = Enuc; % set initial energy
    E = 0;
    dE = abs(E - E0);
    del = 10^-6;
    
iter = 0;
while dE > del;
X = sprintf('%%%%%%%% This is iteration number %d %%%%%',iter);
disp(X)
% 5. Perform the SCF iterations
    F = zeros(nbf, nbf);
    for p = 1:nbf
        for q = 1:nbf
            sumAO = 0;
            for r = 1:nbf
                for s = 1:nbf
                    sumAO = sumAO + D0(r,s)*(2*Int(p,q,r,s)-Int(p,r,q,s));
                end
            end
            F(p,q) = h(p,q) + sumAO; % forming the new Fock matrix
        end
    end
%%%%%%%%%%% generating new values
    % calculate old E0
    sumE = 0;
    for a = 1:nbf
        for b = 1:nbf
            sumE = sumE + D0(a,b)*(h(a,b) + F(a,b));
        end
    end
    E = sumE + Enuc; % calculating the electronic energy
    % generate new F, C, D
    F_ortho = invsqS'*F*invsqS; % Fock matrix transformed to orthonormal basis
    [C_ortho, eps] = eig(F_ortho); % diagonalize new orthonormal Fock
    % sort eVals/eVecs in ascending order
    [~,perm] = sort(diag(eps));
    eps = eps(perm,perm);
    C_ortho = C_ortho(:,perm);
    C = invsqS*C_ortho; % new SCF eVec matrix in orthonormal basis
    D = zeros(nbf, nbf); % construct new density matrix
    for m = 1:nbf % rows
        for n = 1:nbf
            for i = 1:5
                D(m,n) = D(m,n) + C(m,i)*C(n,i);
            end
        end
    end
%%%%%%%%%%% reset all initial guesses to most recent calculated vals
    D0 = D;
    E0 = E; % set old E to new E
    % calculate new E
    sumE = 0;
    for a = 1:nbf
        for b = 1:nbf
            sumE = sumE + D(a,b)*(h(a,b) + F(a,b));
        end
    end
    E = sumE + Enuc % calculating the electronic energy
    dE = abs(E - E0)
%%%%%%%%%%% count number of iterations
    iter = iter + 1;
end

X1 = sprintf('%%%%%%%% Converged in %d iterations %%%%%',iter);
disp(X1);
X2 = sprintf('%%%%%%%% Optimized H-F energy: %d %%%%%',E);
disp(X2);

end
