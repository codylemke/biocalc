import dnacauldron as dc
from Bio import SeqIO
import os.path
import argparse
import sys


class DNAAssembly:

    @classmethod
    def read_assembly_file(cls, infile_handle):
        """
          reads a tsv-formatted intermediate file
          returns a list of DNAAssemblyIntermediate objects
        """
        out = list()

        for (i, line) in enumerate(infile_handle):
            line = line.strip()
            if line != "":
                parts = line.split("\t")  # TODO: might have the first field be the well number
                construct_name = parts[0]
                frags = parts[1:]
                n = cls(construct_name, frags)
                out.append(n)
        return out

    def __init__(self, name, parts):
        self.name = name
        self.parts = parts

    def gg_assemble(self, sequence_repository, enzyme):
        assembly = dc.Type2sRestrictionAssembly(name=self.name, parts=self.parts, enzyme=enzyme, expected_constructs=1,
                                                expect_no_unused_parts=True)
        simulation = assembly.simulate(sequence_repository=sequence_repository)
        if len(simulation.construct_records) == 1:
            return simulation.construct_records
        else:
            print("reaction with %d possible products: %s" % (len(simulation.construct_records), self.name),
                  file=sys.stderr)
            return simulation.construct_records
            # return None


def main(input_handle, output_handle, params):
    parts_library = dict()

    assembly_file = DNAAssembly.read_assembly_file(input_handle)

    for fn in params.input_sequences:
        extension = os.path.splitext(fn)[1][1:]
        for record in SeqIO.parse(fn, extension):
            name = record.name
            assert name not in parts_library
            parts_library[name] = record

    sequence_repository = dc.SequenceRepository({'parts': parts_library})

    out_constructs = list()
    for assembly in assembly_file:
        prods = assembly.gg_assemble(sequence_repository, params.enzyme)
        for prod in prods:
            # if prod is not None:
            out_constructs.append(prod)

    SeqIO.write(out_constructs, output_handle, "gb")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", default=None, help="Assembly file describing which parts should go together.")
    parser.add_argument("-o", default=None, help="genbank file to save the assembled product sequences to.")
    parser.add_argument("input_sequences", default=None, nargs='+',
                        help="genbank file of input sequences.") # TODO: allow fasta files or platemaps also.
    parser.add_argument("--enzyme", type=str, default=None, help="Name of TypeIIS enzyme used for "
                                                                 "the golden gate assembly.")

    args = parser.parse_args()

    if args.i is not None:
        input_handle = open(args.i, "r")
    else:
        input_handle = sys.stdin

    if args.o is not None:
        output_handle = open(args.o, "w")
    else:
        output_handle = sys.stdout

    main(input_handle, output_handle, args)

    if args.i is not None:
        input_handle.close()

    if args.o is not None:
        output_handle.close()
