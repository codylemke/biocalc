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
    
    # def __new__(cls, *args, **kwargs):
    #     """Returns children object instead of Sequence object depending on
    #       the nature of the sequence"""
    #     if kwargs.get('sequence_type') == 'protein':
    #         from .protein import Protein
    #         obj = Protein()
    #     elif kwargs.get('sequence_type') == 'nucleotide':
    #         from .nucleotide import Nucleotide
    #         obj = Nucleotide()
    #     else:
    #         obj = super().__new__(cls)
    #         return obj
    
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
            fasta_entries = fasta_contents.split('>')[1:]
            for entry in fasta_entries:
                lines = entry.split('\n')
                header = lines[0]
                sequence = ''.join(lines[1:]).strip()
                # header_fields = header.split('|')
                # if len(header_fields) > 1:
                #     header_dict = parse_fasta_header(header_fields)
                #     name = header_dict['name']
                # else:
                name = header
                if sequence_type == 'nucleotide':
                    sequence_object = cls(name=name, sequence=sequence)
                elif sequence_type == 'protein':
                    sequence_object = cls(name=name, sequence=sequence)
                elif sequence_type is None:
                    sequence_object = cls(name=name, sequence=sequence)
                else:
                    raise Exception
                sequence_objects.append(sequence_object)
            return sequence_objects

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
                'sp': ['SWISS-PROT', ['accession', 'name']],
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
                'tr': ['TrEMBL', ['accession', 'name']],
                'ENA': ['EBI', ['accession', 'name']]}
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

    @staticmethod
    def create_fasta(sequence_objects):
        """Creates a fasta file including all input DNA sequences"""
        fasta_list = list()
        for sequence_object in sequence_objects:
            fasta_list.append(f'>{sequence_object.name}\n{sequence_object.sequence}\n')
        fasta_string = ''.join(fasta_list).strip()
        return fasta_string

    @staticmethod
    def create_csv(sequence_objects):
        """Creates a fasta file including all input DNA sequences"""
        csv_list = list()
        for sequence_object in sequence_objects:
            csv_list.append(f'{sequence_object.name},{sequence_object.sequence}\n')
        csv_string = ''.join(csv_list).strip()
        return csv_string

    @classmethod
    def fasta_to_csv(cls, fasta_file, sequence_type, output_name):
        """Creates a csv file for sequences from a fasta file"""
        sequences = cls.parse_fasta(fasta_file, sequence_type)
        csv_string = cls.create_csv(sequences)
        output_file = ROOT_DIR / 'scripts' / 'outputs' / f'{output_name}.csv'
        output_file.touch()
        with output_file.open(mode='w') as file:
            file.write(csv_string)
        return

    @staticmethod
    def align_sequences(sequence_objects, algorithm):
        """Generates a multiple sequence alignment using mafft and returns
        the location of the output file"""
        import subprocess
        if all (sequence.sequence_type == 'nucleotide' for sequence in sequence_objects):
            alignment_type = 'nucleotide'
        elif all (sequence.sequence_type == 'protein' for sequence in sequence_objects):
            alignmetn_type = 'protein'
        else:
            raise ValueError('Not all input sequences are of the same type')
            
        fasta_string = Sequence.create_fasta(sequence_objects)
        sequences_path = ROOT_DIR / 'static' / 'temp' / 'sequences.fasta'
        with open(sequences_path, mode='w', encoding='utf-8') as sequences_fasta:
            sequences_fasta.write(fasta_string)
        alignment_path = ROOT_DIR / 'static' / 'temp' / 'alignment.fasta'
        alignment_path.touch()
        subprocess.run(f'mafft --auto --maxiterate 100 --thread 4 {str(sequences_path)} > {str(alignment_path)}', check=True, shell=True)
        with open(alignment_path, mode='r', encoding='utf-8') as alignment_fasta:
            alignment_string = alignment_fasta.read()
        # sequences_path.unlink()
        # alignment_path.unlink()
        return alignment_string

    @staticmethod
    def clustalo_align(sequence_objects):
        """Placeholder"""
        return

    @staticmethod
    def muscle_align(sequence_objects):
        """Placeholder"""
        return

    @staticmethod
    def build_fasttree(fasta_alignment):
        """Generates a phylogenetic tree using fasttree from a multiple
        sequence alignment in fasta format"""
        return

    @staticmethod
    def fasta_to_tree(fasta_file):
        """Returns a phylogenetic tree from a fasta file containing multiple
        sequences."""
        return