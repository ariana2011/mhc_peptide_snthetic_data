# %%
import random
import pandas as pd
import copy
import argparse
random.seed("1024")




parser = argparse.ArgumentParser(description='Give the minimum and maxim lenghth of peptide and MHC moelcule')

parser.add_argument('--sample_num',type=int,default = 100 ,required=True, help = "number of true samples to generate")
parser.add_argument('--decoys',type=int,default = 99 ,required=False, help = "number of decoys per sample")
parser.add_argument('--peptide_min',type=int,default = 7 ,required=False, help = "peptide maximum length")
parser.add_argument('--peptide_max',type=int,default = 15 ,required=False,  help = "peptide maximum length")
parser.add_argument('--mhc_min',type=int,default = 15 ,required=False,  help = "MHC minimum length")
parser.add_argument('--mhc_max',type=int,default = 35 ,required=False,  help = "MHC maximum length")
args = parser.parse_args()


args = parser.parse_args()








# %%
# print(len("YYAGYREKYRQTDVNKLYLRYDSYTWAEWAYEWY"))


# %%
amino_acids = "G A L M F W K Q E S P V I C Y H R N D T"
amino_acid_letter = amino_acids.split()
# print(len(amino_acid_letter)-1)
# amino_acid_letter



def random_letter_maker(N):
    return(''.join(random.choice(amino_acid_letter) for _ in range(N)))
# %%


amino_acids

# %%


def peptide_decoy_maker(a_list_peptide, count=99):
    """
    args:
        a_list(list): a list of splitted peptides and amino acids, index 1,3,5 are the spot of binding
        count(int): number of decoys per sample

    returns:
        temp_set(set): a set of decoys
        labels(list): labels of the decoys ("False")
    """
    temp_set = set()

    # print(d)
    labels = ["False" for i in range(count)]

    while len(temp_set) < count:
        number_of_letters = random.choice([1, 2])
        a_list = copy.deepcopy(a_list_peptide)

        if number_of_letters == 1:
            n = random.choice([1, 3, 5])
            d = random.choice(range(len(amino_acid_letter)-1))

            if a_list[n][0] != amino_acids[d]:  # not to choose the same letter
                replaced_text = a_list[n].replace(
                    a_list[n][0], amino_acid_letter[d])
                a_list[n] = replaced_text
                seq = "".join(a_list)
                seq = seq.replace(" ", "")
                temp_set.add(seq)

        elif number_of_letters == 2:
            n = random.choices([1, 3, 5], k=2)
            d_1 = random.choice(range(len(amino_acid_letter)-1))
            d_2 = random.choice(range(len(amino_acid_letter)-1))

            # not to choose the same letter
            if a_list[n[0]][0] != amino_acids[d_1] and a_list[n[1]][0] != amino_acids[d_2]:
                replaced_text_1 = a_list[n[0]].replace(
                    a_list[n[0]][0], amino_acid_letter[d_1])
                replaced_text_2 = a_list[n[1]].replace(
                    a_list[n[1]][0], amino_acid_letter[d_2])
                a_list[n[0]] = replaced_text_1
                a_list[n[1]] = replaced_text_2
                seq = "".join(a_list)
                seq = seq.replace(" ", "")
                temp_set.add(seq)

    return list(temp_set), labels


# peptide_decoy_maker(['', 'Y', '', 'A', '', 'Q', 'VNWWPS'], count=100)

# %%


def peptide_mhc_sequence(peptide_part_len, binding_pairs, index_selection):
    
    first_part_peptide = random_letter_maker(peptide_part_len[0])
    second_part_peptide = random_letter_maker(peptide_part_len[1])
    third_part_peptide = random_letter_maker(peptide_part_len[2])
    # fourth_part_peptide = random_letter_maker(peptide_part_len[3])

    first_amino_peptide = binding_pairs[0][index_selection]
    second_amino_peptide = binding_pairs[1][index_selection]
    third_amino_peptide = binding_pairs[2][index_selection]

    sequences_list_peptide = [first_part_peptide, first_amino_peptide, second_part_peptide,
                              second_amino_peptide, third_part_peptide, third_amino_peptide]

    peptide_str = "".join(sequences_list_peptide)

    return peptide_str, sequences_list_peptide


# %%
def peptide_maker(amino_pairs=[("A", "G"), ("Y", "E"), ("Q", "N")], peptide_part_len=(2, 4, 2, 2),
                  mhc_part_len=(4, 3, 3, 6), n_flank_len=15, c_flank_len=15, repeat=False, random_choice=True,
                  decoy_count=5):
    """
    This function makes peptied and its corresponding protein.

    :param amino_pairs: a list of tuples of pairing amino acids e.g. [("A","G"),("Y","E"),("Q","N")]
    :param peptide_part_len: a tuple of peptide length around each amino acid first_seq|peptide1|second_seq|peptide2|third_seq|peptide3|foruth_seq|foruth_seq
    :param mhc_part_len: a tuple of peptide length around each amino acid first_seq|peptide1|second_seq|peptide2|third_seq|peptide3|foruth_seq|foruth_seq
    :param repeat: if repate is true, two mathing seq would occure in mhc (NOT IN THIS VERSION)
    :random_choice: if True, the order of binding letters tuples will be random 

    :returns: amindo acid sequence, n_flanc, c_flank, mhc sequence
    :raises ValueError: raises an exception
    """

    if random_choice:
        binding_pairs = random.sample(amino_pairs, k=len(
            amino_pairs))  # choosing 3 pairs randomly
    else:
        binding_pairs = amino_pairs

    # checking the length of petide
    peptide_lenghth = sum(peptide_part_len) + len(binding_pairs)
    if peptide_lenghth < args.peptide_min and peptide_lenghth > args.peptide_max:
        raise ValueError(
            'sum of peptide length should be 7<= peptide length <=15')

    mhc_lenghth = sum(mhc_part_len) + len(binding_pairs)
    if mhc_lenghth < args.mhc_min and mhc_lenghth > args.mhc_max:
        raise ValueError(
            'sum of peptide length should be 7<= peptide length <=15')

    # if not repeat:
    #     if len(mhc_part_len)!= 4:
    #         raise ValueError('mhc_part_len tuple should have 4 values')

    index_selection = random.choice([0, 1])
    pair_selection_mhc = 0 if index_selection == 1 else 1

    random_bind_seq = random.choice([0, 1])

    #n_flank and c_flank
    n_flanc = random_letter_maker(n_flank_len)
    c_flanc = random_letter_maker(c_flank_len)

    # peptide sequence
    
    peptidseq_str, peptidseq_list = peptide_mhc_sequence(
        peptide_part_len, binding_pairs, index_selection)
    Mhcseq_str, Mhcseq_list = peptide_mhc_sequence(
        mhc_part_len, binding_pairs, pair_selection_mhc)

    peptide_or_mhc_choice = random.choice(["peptidseq_list", "Mhcseq_list"])

    peptide_or_mhc_decoy = peptidseq_list if peptide_or_mhc_choice == "peptidseq_list" else Mhcseq_list
    peptide_or_mhc_other = Mhcseq_str if peptide_or_mhc_choice == "peptidseq_list" else peptidseq_str

    peptide_or_mhc_other = [peptide_or_mhc_other for i in range(decoy_count)]
    decoys, deocys_labels = peptide_decoy_maker(
        peptide_or_mhc_decoy, decoy_count)
    # print(peptide_or_mhc_other)
    # print(decoys, deocys_labels)

    if peptide_or_mhc_choice == "peptidseq_list":
        peptide_side_decoy = decoys
        unchanged_mhc = peptide_or_mhc_other
        return (peptidseq_str, Mhcseq_str, "True"), (peptide_side_decoy, unchanged_mhc, deocys_labels)
    else:
        mhc_side_decoy = decoys
        unchanged_peptide = peptide_or_mhc_other
        return (peptidseq_str, Mhcseq_str, "True"), (unchanged_peptide, mhc_side_decoy, deocys_labels)


# %%
def multiple_seq_maker(pairs=[("A", "G"), ("Y", "E"), ("Q", "N"), ("P", "T"), ("M", "K"), ("F", "W")],
                       len_binding_pairs=3, number_of_decoys=99):

    numbers = [1, 2, 3, 4, 5, 6, 7, 8, 9]



    while True:

        peptide_nonbinding_lenghts = random.choices(
            numbers, k=len_binding_pairs)
        Mhc_nonbinding_length = random.choices(numbers, k=len_binding_pairs)

        peptide_lenghth = sum(peptide_nonbinding_lenghts) + len_binding_pairs
        mhc_lenghth = sum(Mhc_nonbinding_length) + len_binding_pairs

        if (peptide_lenghth > args.peptide_min and peptide_lenghth < args.peptide_max) and (mhc_lenghth > args.mhc_min and mhc_lenghth < args.mhc_max):
            break

    seqs, deqs = peptide_maker(amino_pairs=pairs, peptide_part_len=peptide_nonbinding_lenghts,
                               mhc_part_len=Mhc_nonbinding_length, n_flank_len=15, c_flank_len=15, repeat=False, random_choice=True,
                               decoy_count=number_of_decoys)

    return seqs, deqs


# %%
# multiple_seq_maker()


#%%
def main(true_counts, decoy_counts):
    true_seqs = []
    decoys = []
    for i in range(true_counts):
        seqs, deqs = multiple_seq_maker(number_of_decoys=decoy_counts)
        true_seqs.append(seqs)
        decoys.append(deqs)
    peptide_decs = []
    mhc_decs = []
    labels_decs = []
    for sublist in decoys:
        for j in range(len(sublist[0])):
            peptide_decs.append(sublist[0][j])
            mhc_decs.append(sublist[1][j])
            labels_decs.append(sublist[2][j])

    df_decoys = pd.DataFrame(list(zip(peptide_decs, mhc_decs, labels_decs)), columns=[
                             "peptide", "mhc", "label"])
    df_not_decoys = pd.DataFrame(
        true_seqs, columns=["peptide", "mhc", "label"])

    df_decoys.to_csv("decoys_dataset.csv")
    df_not_decoys.to_csv("true_dataset.csv")


# %%
if __name__ == "__main__":
    true_counts = args.sample_num
    decoy_counts = args.decoys  # per true_counts

    main(true_counts, decoy_counts)

