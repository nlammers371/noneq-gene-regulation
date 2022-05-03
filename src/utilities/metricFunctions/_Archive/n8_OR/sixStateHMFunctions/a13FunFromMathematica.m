function a23Sym = a23FunFromMathematica(RSym)

nStates = size(RSym,1);
syms c positive
syms cw positive
syms b positive

for i = 1:nStates
    for j = 1:nStates
        if i~=j 
            string = char(RSym(i,j));
            if length(string)>1
                syms(string,'positive')
            end
        end
    end
end

a13Sym = 