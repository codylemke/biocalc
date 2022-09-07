import sys
from collections import OrderedDict


def parse_diamond(filename,
                  header=None):
    if header is None:
        header = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart",
                  "send", "evalue", "bitscore"]
    out = list()
    # TODO:need to replace this with an iterator so that the entire file is not stored in memory
    with open(filename) as infile:
        for line in infile:
            line = line.strip()
            rec = OrderedDict()
            parts = line.split("\t", len(header) - 1)
            for i, k in enumerate(header):
                rec[k] = parts[i]
            out.append(rec)
    return out


def parse_diamond_iterator(filename,
                           header=None):
    if header is None:
        header = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend",
                  "sstart", "send", "evalue", "bitscore"]
    # out = list()
    with open(filename) as infile:
        for line in infile:
            line = line.strip()
            rec = OrderedDict()
            parts = line.split("\t", len(header) - 1)
            for i, k in enumerate(header):
                rec[k] = parts[i]
            yield rec


def parse_fasta(filename):
    """
        input: the name of a fasta file
        output: a dict where names are sequence id's (id line before first space), and values are the sequence in fasta
                format (">id length\nsequence")
    """
    prev_len = 0
    prev_name = None
    prev_seq = ""
    out = dict()
    with open(filename) as input_handle:
        for line in input_handle:
            line = line.strip()

            if line[0] == ">":
                parts = line.split(None, 1)
                name = parts[0][1:]
                if prev_name is not None:
                    out[prev_name] = ">" + prev_name + " " + str(prev_len) + " \n" + prev_seq
                prev_len = 0
                prev_name = name
                prev_seq = ""
            else:
                prev_len += len(line)
                prev_seq += line
        if prev_name is not None:
            out[prev_name] = ">" + prev_name + " " + str(prev_len) + " \n" + prev_seq
    return out


def parse_clstr(filename):
    cluster_to_members = dict()  # probably could be a set of cluster names, not a dict of sets.
    member_to_cluster = dict()
    with open(filename) as infile:
        members = set()
        captain = None
        for line in infile:
            line = line.strip()
            if line[0] == ">":
                if (captain is not None) and (len(members) > 0):
                    cluster_to_members[captain] = members
                    for member in members:
                        member_to_cluster[member] = captain
                    members = set()
                    captain = None
                elif len(members) > 0:
                    print("Error, cannot find ", file=sys.stderr)
            else:
                parts = line.split(None, 3)
                name = parts[2][1:-3]
                members.add(name)
                if parts[3] == "*":
                    captain = name
        if (captain is not None) and (len(members) > 0):
            cluster_to_members[captain] = members
            for member in members:
                member_to_cluster[member] = captain
    return cluster_to_members, member_to_cluster


def _open_if_is_name(filename_or_handle):
    out = filename_or_handle
    input_type = "handle"
    try:
        out = open(filename_or_handle, "r")
        input_type = "name"
    except TypeError:
        pass
    except Exception as e:
        raise e

    return out, input_type


def parse_simple_list(filename):
    """
        splits a text file on newlines, returns a set with the stripped lines as entries
    """
    out = set()
    (infile, input_type) = _open_if_is_name(filename)

    for line in infile:
        out.add(line.strip())

    if input_type == "name":
        infile.close()
    return out


class SimpleParser:
    def __init__(self, filename, strip=True):
        """
            Generator to parse a file line-by-line, stripping lines and returning any that aren't just whitespace.
        """
        (self.infile, self.input_type) = _open_if_is_name(filename)
        self.strip = strip

    def __iter__(self):
        return self

    def __del__(self):
        if self.input_type == "name":
            self.infile.close()

    def __next__(self):
        line = self.infile.readline()
        while line != '':
            if self.strip:
                line = line.strip()
            if len(line) == 0:
                line = self.infile.readline()
            else:
                return line

        if self.input_type == "name":
            self.infile.close()
        raise StopIteration
