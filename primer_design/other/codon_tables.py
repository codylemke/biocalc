class CodonTable:
    def __init__(self, table_code):
        codon_tables = {
            1: ["FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                "---M------**--*----M---------------M----------------------------", "Standard"],
            2: ["FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
                "----------**--------------------MMMM----------**---M------------", "Vertebrate Mitochondrial"],
            3: ["FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                "----------**----------------------MM----------------------------", "Yeast Mitochondrial"],
            4: ["FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                "--MM------**-------M------------MMMM---------------M------------",
                "Mold Mitochondrial; Protozoan Mitochondrial; Coelenterate Mitochondrial; Mycoplasma; Spiroplasma"],
            5: ["FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",
                "---M------**--------------------MMMM---------------M------------", "Invertebrate Mitochondrial"],
            6: ["FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                "--------------*--------------------M----------------------------",
                "Ciliate, Dasycladacean and Hexamita Nuclear"],
            9: ["FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
                "----------**-----------------------M---------------M------------",
                "Echinoderm and Flatworm Mitochondrial"],
            10: ["FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                 "----------**-----------------------M----------------------------", "Euplotid Nuclear"],
            11: ["FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                 "---M------**--*----M------------MMMM---------------M------------",
                 "Bacterial, Archaeal and Plant Plastid"],
            12: ["FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                 "----------**--*----M---------------M----------------------------", "Alternative Yeast Nuclear"],
            13: ["FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG",
                 "---M------**----------------------MM---------------M------------", "Ascidian Mitochondrial"],
            14: ["FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
                 "-----------*-----------------------M----------------------------",
                 "Alternative Flatworm Mitochondrial"],
            16: ["FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                 "----------*---*--------------------M----------------------------", "Chlorophycean Mitochondrial"],
            21: ["FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
                 "----------**-----------------------M---------------M------------", "Trematode Mitochondrial"],
            22: ["FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                 "------*---*---*--------------------M----------------------------",
                 "Scenedesmus obliquus Mitochondrial"],
            23: ["FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                 "--*-------**--*-----------------M--M---------------M------------", "Thraustochytrium Mitochondrial"],
            24: ["FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",
                 "---M------**-------M---------------M---------------M------------", "Pterobranchia Mitochondrial"],
            25: ["FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                 "---M------**-----------------------M---------------M------------",
                 "Candidate Division SR1 and Gracilibacteria"],
            26: ["FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                 "----------**--*----M---------------M----------------------------", "Pachysolen tannophilus Nuclear"],
            27: ["FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                 "--------------*--------------------M----------------------------", "Karyorelict Nuclear"],
            28: ["FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                 "----------**--*--------------------M----------------------------", "Condylostoma Nuclear"],
            29: ["FFLLSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                 "--------------*--------------------M----------------------------", "Mesodinium Nuclear"],
            30: ["FFLLSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                 "--------------*--------------------M----------------------------", "Peritrich Nuclear"],
            31: ["FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                 "----------**-----------------------M----------------------------", "Blastocrithidia Nuclear"]
        }
        bases = "TCAG"
        codons = [a + b + c for a in bases for b in bases for c in bases]
        amino_acids = codon_tables[table_code][0]
        start_stop = codon_tables[table_code][1]
        self.name = codon_tables[table_code][2]

        self._table = dict(zip(codons, amino_acids))
        self._start_stop = dict(zip(codons, start_stop))
        self._aa_to_codons = dict()  # keys are aa chars, values are lists of codon sequences
        for (codon, aa) in self._table.items():
            if aa not in self._aa_to_codons:
                self._aa_to_codons[aa] = list()
            self._aa_to_codons[aa].append(codon)

    # def load_codon_usage

    def special(self, key):
        key = str(key).upper()
        return self._start_stop[key]

    def aa_to_codons(self, aa):
        return self._aa_to_codons[aa]

    def __getitem__(self, key):  # key is a codon, returns an amino acid.
        key = str(key).upper()
        return self._table[key]


eukaryote = {"ATT": "I",
             "ATC": "I",
             "ATA": "I",
             "CTT": "L",
             "CTC": "L",
             "CTA": "L",
             "CTG": "L",
             "TTA": "L",
             "TTG": "L",
             "GTT": "V",
             "GTC": "V",
             "GTA": "V",
             "GTG": "V",
             "TTT": "F",
             "TTC": "F",
             "ATG": "M",
             "TGT": "C",
             "TGC": "C",
             "GCT": "A",
             "GCC": "A",
             "GCA": "A",
             "GCG": "A",
             "GGT": "G",
             "GGC": "G",
             "GGA": "G",
             "GGG": "G",
             "CCT": "P",
             "CCC": "P",
             "CCA": "P",
             "CCG": "P",
             "ACT": "T",
             "ACC": "T",
             "ACA": "T",
             "ACG": "T",
             "TCT": "S",
             "TCC": "S",
             "TCA": "S",
             "TCG": "S",
             "AGT": "S",
             "AGC": "S",
             "TAT": "Y",
             "TAC": "Y",
             "TGG": "W",
             "CAA": "Q",
             "CAG": "Q",
             "AAT": "N",
             "AAC": "N",
             "CAT": "H",
             "CAC": "H",
             "GAA": "E",
             "GAG": "E",
             "GAT": "D",
             "GAC": "D",
             "AAA": "K",
             "AAG": "K",
             "CGT": "R",
             "CGC": "R",
             "CGA": "R",
             "CGG": "R",
             "AGA": "R",
             "AGG": "R",
             "TAA": "Stop",
             "TAG": "Stop",
             "TGA": "Stop"}
