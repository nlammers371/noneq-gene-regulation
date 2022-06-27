function rate = sixStateProductionRate_v2(RSym)

nStates = size(RSym,1);
syms c positive
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

rate = ka.*(cw.*kb.*(ka.*wab+ki.*wib).*(ki.*(ki+(ka+(c+cw).*kb).*wba).* ...
  wib+(ki.*(ku+ka.*wab)+(ka+cw.*kb).*ku.*wba+ka.*(ka+(c+cw).*kb).* ...
  wab.*wba).*wua)+a.^2.*ku.*wua.*(c.^2.*kb.^2.*wba.*(ka.*wab+ki.* ...
  wib)+(ka+ki).*ku.*(ki.*wib+(ku+ka.*wab).*wua)+c.*kb.*(ka.^2.*wab.* ...
  wba+ki.^2.*wib+ka.*ki.*(wab+wba.*wib)+ki.*ku.*(wba.*wib+wua)+ka.* ...
  ku.*(wba+wab.*wua)))+a.*(c.^2.*kb.^2.*wba.*(ka.*wab+ki.*wib).*( ...
  ki.*wib+ka.*wab.*wua)+c.*kb.*(ki.*wib.*(ka.*ki.*wab+ka.*(ku+(ka+ ...
  cw.*kb).*wab).*wba+ki.*(ki+(ka+cw.*kb+ku).*wba).*wib)+(ka.*wab.*( ...
  ka.*ki.*wab+(ka+2.*cw.*kb).*ku.*wba+ka.*(ka+cw.*kb).*wab.*wba)+ ...
  ki.*(ka.*ku.*wab+ki.*(ku+ka.*wab)+2.*cw.*kb.*ku.*wba+ka.*(ka+cw.* ...
  kb+ku).*wab.*wba).*wib).*wua+ka.*ku.*wab.*(ki+ka.*wab).*wua.^2)+ ...
  ku.*(ki.*wib+(ku+ka.*wab).*wua).*((ka+ki).*(ki.*wib+ka.*wab.*wua)+ ...
  cw.*kb.*(ka.*wba+ki.*wba.*wib+ki.*wua+ka.*wab.*wua)))).^(-1).*( ...
  cw.*kb.*wab.*(ki.*(ki+(ka+(c+cw).*kb).*wba).*wib+(ki.*(ku+ka.*wab) ...
  +(ka+cw.*kb).*ku.*wba+ka.*(ka+(c+cw).*kb).*wab.*wba).*wua)+a.^2.* ...
  ku.*wua.*(c.*kb.*ku.*wba+c.^2.*kb.^2.*wab.*wba+ki.*ku.*wib+ku.*( ...
  ku+ka.*wab).*wua+c.*kb.*wab.*(ki+ka.*wba+ku.*wua))+a.*(c.^2.* ...
  kb.^2.*wab.*wba.*(ki.*wib+ka.*wab.*wua)+ku.*(ki.*wib+(ku+ka.*wab) ...
  .*wua).*(ki.*wib+ka.*wab.*wua+cw.*kb.*(wba+wab.*wua))+c.*kb.*( ...
  ki.^2.*wab.*wib+ki.*(ku+ka.*wab+cw.*kb.*wab).*wba.*wib+ki.*wab.*( ...
  ka.*wab+ku.*wib).*wua+wab.*wua.*((ka+2.*cw.*kb).*ku.*wba+ka.*(ka+ ...
  cw.*kb).*wab.*wba+ka.*ku.*wab.*wua))));