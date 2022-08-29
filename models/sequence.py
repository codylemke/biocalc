"""Contains the Sequence object

Description...

Usage Example:
    ...example
"""
# Standard Libraries
import pathlib
import re
# Global Constants
ROOT_DIR = pathlib.Path(__file__).resolve().parent.parent


# CLASSES ---------------------------------------------------------------------
class Sequence:
    """Representation of DNA"""
    def __init__(self, sequence, name='generic_sequence', sequence_type=None,
    source_organism=None, sequence_class= None, description=None,
    accession=None, code_number=None, author=None, comments=None,
    references=None, files=None):
        """Constructor Function"""
        self.name = name
        self.sequence = sequence
        self.sequence_type = sequence_type
        self.source_organism = source_organism
        self.sequence_class = sequence_class
        self.description = description
        self.accession = accession
        self.code_number = code_number
        self.author = author
        self.comments = comments
        self.references = references
        self.files = files
        return
    
    @property
    def sequence(self):
        """Sequence validation"""
        return self._sequence

    @sequence.setter
    def sequence(self, seq):
        from .protein import Protein
        from .nucleotide import Nucleotide
        seq = (seq.replace(' ', '').replace('\n', '').replace('-', '')
               .replace('_', '').strip().upper())
        if not all(char in Nucleotide.standard_bases
                   or char in Nucleotide.ambiguous_bases
                   for char in seq):
            if all(char in Protein.standard_residues
                   or char in Protein.ambiguous_residues
                   for char in seq):
                raise ValueError(
                    'The sequence entered appears to be a protein sequence.')
            else:
                raise ValueError(
                    'The sequence entered does not appear to be biological.')
        if self.nucleotide_type == 'generic':
            if 'U' in seq and 'T' not in seq:
                self.nucleotide_type = 'RNA'
            elif 'T' in seq and 'U' not in seq:
                self.nucleotide_type = 'DNA'
            else:
                raise AttributeError(textwrap.dedent("""\
                    The `nucleotide_type` could not be determined.
                    Please manually specify the `nucleotide_type` as a keyword argument.
                    eg. DNA or RNA"""))
        elif self.nucleotide_type == 'DNA':
            if 'U' in seq and 'T' not in seq:
                raise AttributeError(
                    'The `nucleotide_type` attribute was '
                    'specified as DNA but appears to be RNA')
        elif self.nucleotide_type == 'RNA':
            if 'T' in seq and 'U' not in seq:
                raise AttributeError(
                    'The `nucleotide_type` attribute was '
                    'specified as RNA but appears to be DNA')
        self._sequence = seq
        return

    @classmethod
    def parse_fasta(cls, fasta_file):
        """Creates sequence objects for each sequence in a fasta file"""
        sequence_objects = list()
        with fasta_file.open(mode='r') as file:
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
                    origin = standard_headers[header_fields[0]]
                    accession = header_fields[1]
                    if origin in ['SWISS-PROT', 'TrEMBL']:
                        name = header_fields[2].split()[0]
                        organism = re.search(r'(?<=OS\=)(.+)(?=\sOX\=)', header_fields[2])
                        taxid = re.search(r'(?<=OX\=)(.+)(?=\sGN\=)', header_fields[2])
                        sequence_object = cls(name=name, sequence=sequence, origin=origin,
                        organism=organism, taxid=taxid)
                    else:
                        sequence_object = cls(name=accession, sequence=sequence, accession=accession,)
            else:
                sequence_object = cls(name=header, sequence=sequence)
            sequence_objects.append(sequence_object)
        return sequence_objects