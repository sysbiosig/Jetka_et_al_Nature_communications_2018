function result = Two_JAKSTATt_extr4_small_ver2_stimulus(t,i)
    stim={};
    stim{1}=1*(t>0&t<5);
    stim{2}=1*(t>0&t<5);

    result=stim{i};
end

