#! /usr/bin/python
"""
TODO: this is kind of broken, in terms of being broken for a specific application.
"""
import sys
import argparse
from Bio import SeqIO
from Bio.SeqUtils import MeltingTemp


def main(inpt, output, params):
    sides = set()
    if (params.direction == "both") or (params.direction == "f"):
        sides.add("forward")
    if (params.direction == "both") or (params.direction == "r"):
        sides.add("reverse")

    for seq in SeqIO.parse(inpt, "fasta"):
        max_len = params.max_length
        if max_len > len(seq):
            max_len = len(seq)
        for side in sides:
            # for side in ("reverse",):
            template = seq.seq
            if side == "reverse":
                template = template.reverse_complement()
            best_primer = dict()  # seq, tm
            for l in range(params.min_length, max_len + 1):
                sub = template[0:l]
                tm = MeltingTemp.Tm_NN(sub, Na=params.Na, K=params.K, Tris=params.Tris, Mg=params.Mg, dNTPs=params.dNTPs
                                       , saltcorr=6, dnac1=params.primer_conc / 2, dnac2=params.primer_conc / 2)
                if (sub[-1] in "GC") or (not params.gc_clamp):
                    if (tm >= params.minTm) and (tm <= params.maxTm):
                        if "seq" not in best_primer:
                            best_primer['seq'] = sub
                            best_primer['tm'] = tm
                        elif abs(params.targetTm - tm) < abs(params.targetTm - best_primer['tm']):
                            best_primer['seq'] = sub
                            best_primer['tm'] = tm
            if len(best_primer) == 0:
                print("could not find %s primer for seq %s" % (side, seq.id), file=sys.stderr)
            else:
                # print(">%s_tm_%0.1f_len_%d\t%s\t%s" \
                # % (side[0],best_primer['tm'],len(best_primer['seq']),seq.id,best_primer['seq']), file=output)
                if params.output_format == "fasta":
                    out_linestart = ">"
                    out_separator = "\n"
                else:
                    out_linestart = ""
                    out_separator = "\t"
                print("%s%s_tm_%0.1f_len_%d %s%s%s" % (
                  out_linestart, side[0], best_primer['tm'], len(best_primer['seq']), seq.id, out_separator,
                  params.overhang_5 + best_primer['seq'] + params.overhang_3), file=output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", default=None)
    parser.add_argument("-o", default=None)
    parser.add_argument("--minTm", type=float, default=54,
                        help="Will report an error if a primer can't be found with a Tm of at least this.")
    parser.add_argument("--maxTm", type=float, default=70,
                        help="Will report an error if a primer can't be found with a Tm of at most this.")
    parser.add_argument("--targetTm", type=float, default=65, help="Will try to find a primer with a Tm clost to this.")
    parser.add_argument("--min_length", type=int, default=12,
                        help="minimum acceptable length for the core primer (without the overhangs attached)")
    parser.add_argument("--max_length", type=int, default=25,
                        help="maximum acceptable length for the core primer (without the overhangs attached)")
    parser.add_argument("--direction", type=str, default="both", choices={"both", "f", "r"},
                        help="which direction to design primers for")
    parser.add_argument("--output_format", type=str, default="table", choices={"table", "fasta"},
                        help="format to write output file in")
    parser.add_argument('--gc_clamp', action="store_true", default=False)

    parser.add_argument("--Na", type=float, default=50, help="Sodium ion concentration [50]")
    parser.add_argument("--K", type=float, default=0, help="Potassium ion concentration [0]")
    parser.add_argument("--Tris", type=float, default=0, help="Tris ion concentration [0]")
    parser.add_argument("--Mg", type=float, default=100, help="Magnesium ion concentration [100]")
    parser.add_argument("--dNTPs", type=float, default=50, help="dNTP ion concentration [50]")
    parser.add_argument("--primer_conc", type=float, default=250, help="primer concentration [250]")

    parser.add_argument("--overhang_5", type=str, default="", help="5' overhang to attach to all primers")
    parser.add_argument("--overhang_3", type=str, default="", help="3' overhang to attach to all primers")

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
