"""Title

Description

Usage Example:
    example
"""
# Local Modules

# Global Constants
adapters = {
    5: 'GGTTAACCGCAATGAAGACTG',
    3: 'GTGTCTTCTAACGCCAATTGG'}

overhangs = {
    1: {5: 'CCTC', 3: 'CATA'},
    2: {5: 'CATA', 3: 'AAAA'},
    3: {5: 'AAAA', 3: 'AGGA'},
    4: {5: 'AGGA', 3: 'AGTG'},
    5: {5: 'AGTG', 3: 'CAGC'},
    6: {5: 'CAGC', 3: 'TGAA'},
    7: {5: 'TGAA', 3: 'ATTA'},
    8: {5: 'ATTA', 3: 'AATC'},
    9: {5: 'AATC', 3: 'CCAG'}}


# FUNCTIONS -------------------------------------------------------------------
def append_gge_adapters(sequence, organism, module):
    """Appends the appropriate adapters for ordering"""
    # For M1, simply input only the full sequence of the homology arm
    # For M2, simply input only the full sequence of the promoter
    # For M3, input only the full sequence desired to be between the end of
    # the promoter and the start of the kozak sequence for eukaryotes or the
    # end of the promoter and the end of the RBS for eukaryotes.
    # For M4, simply input the full sequence of the n-tag from the start codon
    # to the end of the tag
    # For M5, simply input the full sequence from the start codon to the stop
    # codon
    # For M6, simply input the full sequence from the beginning of the c-tag
    # to the stop codon
    # For M7, simply input the full sequence from the end of the stop codon to
    # the beginning of the terminator
    # For M8, simply input the full sequence from the beginning to the end of
    # the terminator
    # For M9, simply input the full sequence of the homology arm
    def determine_left_overhang():
        """Returns the sequence of the portion of the overhang that should be
        appended to the sequence"""
        if sequence[:4] == overhangs[module][5]:
            overhang = ''
        elif sequence[:3] == overhangs[module][5][1:]:
            overhang = overhangs[module][5][:1]
        elif sequence[:2] == overhangs[module][5][2:]:
            overhang = overhangs[module][5][:2]
        elif sequence[0] == overhangs[module][5][3:]:
            overhang = overhangs[module][5][:3]
        else:
            overhang = overhangs[module][5]
        return overhang

    def determine_right_overhang():
        """Returns the sequence of the portion of the overhang that should be
        appended to the sequence"""
        if sequence[-4:] == overhangs[module][3]:
            overhang = ''
        elif sequence[-3:] == overhangs[module][3][:3]:
            overhang = overhangs[module][3][-1]
        elif sequence[-2:] == overhangs[module][3][:2]:
            overhang = overhangs[module][3][-2:]
        elif sequence[-1] == overhangs[module][3][:1]:
            overhang = overhangs[module][3][-3:]
        else:
            overhang = overhangs[module][3]
        return overhang

    left_overhang = determine_left_overhang()
    right_overhang = determine_right_overhang()
    if module == 1: # Homology Arm 1
        final_sequence = adapters[5]+left_overhang+sequence+right_overhang+adapters[3]
    if module == 2: # Promoter
        final_sequence = adapters[5]+left_overhang+sequence+right_overhang+adapters[3]
    if module == 3: # 5'UTR
        if organism == 'e_coli':
            for index in range(0, len(sequence)):
                if sequence[-(4+index):-index] == 'AGGA':
                    sequence = sequence[:-(4+index)]
            final_sequence = adapters[5]+left_overhang+sequence+overhangs[3][3]+adapters[3]
        elif organism == 'yeast' or organism == 'mammalian':
            final_sequence = adapters[5]+left_overhang+sequence+right_overhang+adapters[3]
        else:
            raise ValueError
    if module == 4: # TIS/n-tag
        if sequence[:3] == 'ATG':
            sequence = sequence[3:]
        if organism == 'e_coli':
            tis = 'AGGAGAGCAGCTATG'
        elif organism == 'e_coli_enhanced':
            tis = 'AGGAGAGCAGCTATGCAGCTT'
        elif organism == 'yeast':
            tis = 'AGGAAAAAAATGTCT'
        elif organism == 'mammalian':
            tis = 'AGGAGCCACCATGGGC'
        else:
            raise ValueError
        final_sequence = adapters[5]+tis+sequence+'AGTAGTG'+adapters[3]
    if module == 5: # Gene
        if sequence[:3] == 'ATG':
            sequence = sequence[3:]
        if sequence[-3:] == 'TAA' or sequence[-3:] == 'TGA' or sequence[-3:] == 'TAG':
            sequence = sequence[:-3]
        final_sequence = adapters[5]+'AGTGGT'+sequence+'GGTAGCAGC'+adapters[3]
    if module == 6: # c-tag
        if sequence[-3:] != 'TAA':
            if sequence[-3:] == 'TGA' or sequence[-3:] == 'TAG':
                sequence = sequence[:-3]+'TAA'
            else:
                sequence = sequence+'TAA'
        final_sequence = adapters[5]+overhangs[6][5]+sequence+overhangs[6][3]+adapters[3]
    if module == 7: # 3'UTR
        final_sequence = adapters[5]+left_overhang+sequence+right_overhang+adapters[3]
    if module == 8: # Terminator
        final_sequence = adapters[5]+left_overhang+sequence+right_overhang+adapters[3]
    if module == 9: # Homology Arm 2
        final_sequence = adapters[5]+left_overhang+sequence+right_overhang+adapters[3]
    return final_sequence
