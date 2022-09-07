from Bio import SeqIO
import pandas as pd

input = open(r"..\Experiments\211019 BM16\input\BM16_4.3_pETDuet translation.fasta", "r")
seqs = list(SeqIO.parse(input, "fasta"))

seq = seqs[0].seq
seq_str = str(seq)
# print(seqs)

all_mutations = []
ids = []
for idx, char in enumerate(seq_str.split('*')[0]):
    if not char == '*':
        mutation_code = char + str(idx+1) + 'ARNDCEQGHILKMFPSTWYV'
        all_mutations.append(mutation_code)
        ids.append('TEST_' + str(idx+1).zfill(4))
    else:
        print('stop codon')
        break

df_out = pd.DataFrame({"Construct Name": ids, "Target Mutation": all_mutations})
exp_dir = r"..\Experiments\\211019 BM16\\input\\"
out_file_name = "TEST_input_ALL_mutations.xlsx"
df_out.to_excel(exp_dir + out_file_name, index=False)
