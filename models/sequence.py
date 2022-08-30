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
    """Representation of a biological sequence"""
    def __init__(self, sequence, name='generic_sequence', sequence_type=None,
                 source=None, accession=None, accession_db=None,
                 accession_index=None):
        """Constructor Function"""
        self.name = name
        self.sequence_type = sequence_type
        self.sequence = sequence
        self.source = source
        self.accession = accession
        self.accession_db = accession_db
        self.accession_index= accession_index
        return
    
    @property
    def sequence(self):
        """Sequence validation"""
        return self._sequence

    @sequence.setter
    def sequence(self, seq):
        from .nucleotide import Nucleotide
        from .protein import Protein
        seq = (seq.replace(' ', '').replace('\n', '').replace('-', '')
               .replace('_', '').strip().upper())
        if self.sequence_type == 'nucleotide':
            if not all(char in Nucleotide.standard_bases
                       or char in Nucleotide.ambiguous_bases
                       or char in Nucleotide.nonstandard_bases for char in seq):
                raise ValueError('Sequence contains non-nucleotide characters')
        elif self.sequence_type == 'protein':
            if not all(char in Protein.standard_residues
                       or char in Protein.ambiguous_residues
                       or char in Protein.nonstandard_residues for char in seq):
                raise ValueError('Sequence contains non-amino acid characters')
        elif self.sequence_type is None:
            if all(char in Nucleotide.standard_bases for char in seq):
                self.sequence_type = 'nucleotide'
            elif all(char in Protein.standard_residues for char in seq):
                self.sequence_type = 'protein'
            else:
                raise ValueError('The sequence type could not be determined.')
        self._sequence = seq
        return

    @property
    def length(self):
        """Return the number of base pairs in the sequence"""
        char_count = len(self.sequence)
        return char_count

    @classmethod
    def parse_fasta(cls, fasta_file, sequence_type):
        """Creates sequence objects for each sequence in a fasta file"""
        
        def main():
            sequence_objects = list()
            with fasta_file.open(mode='r') as file:
                fasta_contents = file.read()
            fasta_entries = fasta_contents.split('>')
            for entry in fasta_entries:
                lines = entry.split('\n')
                header = lines[0]
                sequence = ''.join(lines[1:]).strip()
                header_fields = header.split('|')
                if len(header_fields) > 1:
                    header_dict = parse_fasta_header(header_fields)
            return
            
        def parse_fasta_header(header_fields):
            """Parses the fasta file header for conventional annotations"""
            standard_headers = {
            # abbreviation: [database [conventional annotations]]
                'lcl': ['local', ['accession']],
                'bbs': ['GenInfo backbone seqid', ['accession']],
                'bbm': ['GenInfo backbone moltype', ['accession']],
                'gim': ['GenInfo import ID', ['accession']],
                'gb': ['GenBank', ['accession', 'locus']],
                'emb': ['EMBL', ['accession', 'locus']],
                'pir': ['PIR', ['accession', 'locus']],
                'sp': ['SWISS-PROT', ['accession', 'locus']],
                'pat': ['patent', ['country', 'patent', 'sequence_number']],
                'pgp': ['pre-grant patent', ['country', 'application_number', 'sequence_number']],
                'ref': ['RefSeq', ['accession', 'name']],
                'gnl': ['general database reference', ['database', 'accession']],
                'gi': ['GenInfo integrated database', ['accession']],
                'dbj': ['DDBJ', ['accession', 'locus']],
                'prf': ['PRF', ['acccession', 'name']],
                'pdb': ['PDB', ['accession', 'chain']],
                'tpg': ['third-party GenBank', ['accession', 'name']],
                'tpe': ['third-party EMBL', ['acccession', 'name']],
                'tpd': ['third-party DDBJ', ['accession', 'name']],
                'tr': ['TrEMBL', ['accession', 'name']]}
            header_dict = dict()
            try:
                fasta_format = standard_headers[header_fields[0]]
            except KeyError as err:
                raise KeyError('Header type is unrecognized.') from err
            else:
                database = fasta_format[0]
                if 'accession' in fasta_format[1]:
                    header_dict['accession'] = header_fields[fasta_format[1].index('accession')]
                if 'locus' in fasta_format[1]:
                    header_dict['locus'] = header_fields[fasta_format[1].index('locus')]
                if 'name' in fasta_format[1]:
                    header_dict['name'] = header_fields[fasta_format[1].index('name')]
                if 'country' in fasta_format[1]:
                    header_dict['country'] = header_fields[fasta_format[1].index('country')]
                if 'patent' in fasta_format[1]:
                    header_dict['patent'] = header_fields[fasta_format[1].index('patent')]
                if 'sequence_number' in fasta_format[1]:
                    header_dict['sequence_number'] = header_fields[fasta_format[1].index('sequence_number')]
                if 'chain' in fasta_format[1]:
                    header_dict['chain'] = header_fields[fasta_format[1].index('chain')]
                if database in ['SWISS-PROT', 'TrEMBL']:
                    details = header_dict['name']
                    header_dict['name'] = details.split()[0]
                    header_dict['protein_name'] = re.search(r'(?<=\W*)(.+)(?=\sOS\=)', details)
                    header_dict['organism_name'] = re.search(r'(?<=OS\=)(.+)(?=\sOX\=)', details)
                    if 'GN=' in details:
                        header_dict['organism_taxid'] = re.search(r'(?<=OX\=)(.+)(?=\sGN\=)', details)
                        header_dict['gene_name'] = re.search(r'(?<=GN\=)(.+)(?=\sPE\=)', details)
                    else:
                        header_dict['organism_taxid'] = re.search(r'(?<=OX\=)(.+)(?=\sPE\=)', details)
                    header_dict['protein_existence'] = re.search(r'(?<=PE\=)(.+)(?=\sSV\=)', details)
                    header_dict['sequence_version'] = re.search(r'(?<=SV\=)(.+)$', details)
            return header_dict
            
        return main()
