"""Takes the sequences from a fasta file and appends GGE adapters to them for
Module 5

Description

Usage Example:
    example
"""
# Standard Library
import pathlib
import re
# Local Modules
from models import DNA
# Global Constants
ROOT_DIR = pathlib.Path(__file__).resolve().parent.parent


# FUNCTIONS -------------------------------------------------------------------
def parse_fasta(fasta_file):
    """Creates Sequence objects for each sequence in a fasta file"""
    with open(fasta_file.open(mode='r')) as file:
        fasta_contents = file.read()
    fasta_entries = fasta_contents.split('>')
    for entry in fasta_entries:
        lines = entry.split('\n')
        header = lines[0]
        sequence_list = list()
        for line in lines[1:]:
            sequence_list.append(line)
        sequence = ''.join(sequence_list)
        header_fields = header.split('|')
        standard_headers = {
            'lcl': 'local',
            'bbs': 'GenInfo backbone seqid',
            'bbm': 'GenInfo backbone moltype',
            'gim': 'GenInfo import ID',
            'gb': 'GenBank',
            'emb': 'EMBL',
            'pir': 'PIR',
            'sp': 'SWISS-PROT',
            'pat': 'patent',
            'pgp': 'pre-grant patent',
            'ref': 'RefSeq',
            'gnl': 'general database reference',
            'gi': 'GenInfo integrated database',
            'dbj': 'DDBJ',
            'prf': 'PRF',
            'pdb': 'PDB',
            'tpg': 'third-party GenBank',
            'tpe': 'third-party EMBL',
            'tpd': 'third-party DDBJ',
            'tr': 'TrEMBL',
        }
        if len(header_fields) > 1:
            if header_fields[0] in standard_headers:
                sequence_type = standard_headers[header_fields[0]]
                accession = header_fields[1]
                if sequence_type in ['SWISS-PROT', 'TrEMBL']:
                    name = header_fields[2].split()[0]
                    organism = re.search(r'(?<=OS\=)(.+)(?=\sOX\=)', header_fields[2])
                    taxid = re.search(r'(?<=OX\=)(.+)(?=\sGN\=)', header_fields[2])

    return

# SCRIPT ----------------------------------------------------------------------
if __name__ == '__main__':
    pass