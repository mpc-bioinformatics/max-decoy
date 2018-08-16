peps = [pep1, pep2, pep3, pep4, pep5]
mcs [0,1,2,3]

012
013
023
123

def glue()
    new_peps = []
    pep = ""
    for i in 0..mcs.length
        pep += peps[mcs[i]]
        if (i + 1) > mcs.length or mcs[i+1] != mcs[i] + 1   # if (i + 1) is not an index of msc and if the next entry in msc is not 
            new_peps << pep
            pep = ""
        end

    end
end