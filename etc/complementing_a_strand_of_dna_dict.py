import sys
#loop
with open(sys.argv[1], 'r') as f:
    def complementary_dna_dict(dna):
        #dictionary_preemble
        complementStrand = ""
        
        complementNuc = {"G":"C", "C":"G", "T":"A", "A":"T"}

        contents=f.readlines()
        for nuc in reversed(contents):
            complementStrand += complementNuc[nuc]
        #output
        print(complementStrand)
f.close()

