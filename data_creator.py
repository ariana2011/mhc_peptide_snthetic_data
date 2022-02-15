# %%
import random
from unicodedata import decimal
import pandas as pd
import copy
import argparse
import logging
logging.basicConfig(level=logging.DEBUG)
logging.debug('This will get logged')
logging.disable()
# random.seed("1024")


parser = argparse.ArgumentParser(
    description='Give the minimum and maxim lenghth of peptide and MHC moelcule')

parser.add_argument('--sample_num', type=int, default=100,
                    required=True, help="number of true samples to generate")
parser.add_argument('--decoys', type=int, default=99,
                    required=False, help="number of decoys per sample")
parser.add_argument('--peptide_min', type=int, default=7,
                    required=False, help="peptide maximum length")
parser.add_argument('--peptide_max', type=int, default=15,
                    required=False,  help="peptide maximum length")
parser.add_argument('--mhc_min', type=int, default=15,
                    required=False,  help="MHC minimum length")
parser.add_argument('--mhc_max', type=int, default=35,
                    required=False,  help="MHC maximum length")
args = parser.parse_args()


args = parser.parse_args()


# %%
amino_acids = "G A L M F W K Q E S P V I C Y H R N D T"
amino_acid_letter = amino_acids.split()
# logging.disable(logging.INFO)


class PeptideMhcPair():

    def __init__(self, pairs):

        self.pairs = pairs

    @staticmethod
    def random_letter_maker(N):
        return(''.join(random.choice(amino_acid_letter) for _ in range(N)))

    def _peptide_mhc_sequence(self, peptide_part_len, binding_pairs, index_selection):
        """
        A function utilized for creating peptids OR Mhc sequence
        args:
            peptide_part_len("tupl or list"): a tuble or list indiciating the number of random letter between the binding pairs e.g.: (2,4,3,4) ==> aa b aaaa b aaa b aaaa
            binding_pairs(list of tuples): a list of the tuples of binding amino acid pair letter e.g. [("A", "G"), ("Y", "E"), ("Q", "N")]
            index_selection(int): a radnom number among 0,1 to select the direction of binding from amino acid pairs e.g. ("A", "G") ==> 0 =="A" , 1 == "G"
        """

        first_part_peptide = PeptideMhcPair.random_letter_maker(
            peptide_part_len[0])
        second_part_peptide = PeptideMhcPair.random_letter_maker(
            peptide_part_len[1])
        third_part_peptide = PeptideMhcPair.random_letter_maker(
            peptide_part_len[2])
        fourth_part_peptide = PeptideMhcPair.random_letter_maker(
            peptide_part_len[3])

        first_amino_peptide = binding_pairs[0][index_selection]
        second_amino_peptide = binding_pairs[1][index_selection]
        third_amino_peptide = binding_pairs[2][index_selection]

        sequences_list_peptide = [first_part_peptide, first_amino_peptide, second_part_peptide,
                                  second_amino_peptide, third_part_peptide, third_amino_peptide, fourth_part_peptide]

        peptide_str = "".join(sequences_list_peptide)

        return peptide_str, sequences_list_peptide



    def flank_maker(self, len_flank):
        return PeptideMhcPair.random_letter_maker(len_flank)
    
    

    def peptide_maker(self, amino_pairs=[("A", "G"), ("Y", "E"), ("Q", "N")], peptide_part_len=(2, 4, 2, 2),
                      mhc_part_len=(4, 3, 3, 20), n_flank_len=15, c_flank_len=15, random_choice=True):
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


        index_selection = random.choice([0, 1])
        pair_selection_mhc = 0 if index_selection == 1 else 1

        random_bind_seq = random.choice([0, 1])

        #n_flank and c_flank
        n_flanc = self.flank_maker(n_flank_len)
        c_flanc = self.flank_maker(c_flank_len)

        peptidseq_str, peptidseq_list = self._peptide_mhc_sequence(
            peptide_part_len, binding_pairs, index_selection)
        Mhcseq_str, Mhcseq_list = self._peptide_mhc_sequence(
            mhc_part_len, binding_pairs, pair_selection_mhc)

        return peptidseq_str, peptidseq_list, Mhcseq_str, Mhcseq_list

    def _peptide_decoy_maker(self, a_list_peptide):
        """
        args:
            a_list(list): a list of splitted peptides and amino acids, index 1,3,5 are the spots of binding
            count(int): number of decoys per sample

        returns:
            seq(str): altered version of the originial sequence
        """
        # print(a_list_peptide)
        number_of_letters = random.choice([1, 2])
        a_list = copy.deepcopy(a_list_peptide)

        
        while True:
        
            if number_of_letters == 1:
                n = random.choice([1, 3, 5])
                d = random.choice(range(len(amino_acid_letter)))

                # if a_list[n][0] != amino_acid_letter[d]:  # not to choose the same letter: Does not work (15/2)
                replaced_text = a_list[n].replace(
                a_list[n][0], amino_acid_letter[d])
                a_list[n] = replaced_text

            elif number_of_letters == 2:
                n = random.choices([1, 3, 5], k=2)
                d_1 = random.choice(range(len(amino_acid_letter)))
                d_2 = random.choice(range(len(amino_acid_letter)))

                # not to choose the same letter
                # if a_list[n[0]][0] != amino_acids[d_1] and a_list[n[1]][0] != amino_acids[d_2]:
                replaced_text_1 = a_list[n[0]].replace(
                    a_list[n[0]][0], amino_acid_letter[d_1])
                replaced_text_2 = a_list[n[1]].replace(
                    a_list[n[1]][0], amino_acid_letter[d_2])
                a_list[n[0]] = replaced_text_1
                a_list[n[1]] = replaced_text_2

            if "".join(a_list_peptide) != "".join(a_list):
                break
                # print("".join(a_list_peptide) , "HHHHHHHHHHHHHHi","".join(a_list))
        seq = "".join(a_list)
        seq = seq.replace(" ", "")
        return seq

    def decoy_output(self, peptidseq_str, peptidseq_list, Mhcseq_str, Mhcseq_list, decoy_count=10):
        """

        return:
            temp_set(list): a listed set of (peptide_decoy, mhc_decoy,"False")
        """

        temp_set = set()

        while len(temp_set) < decoy_count:
            peptide_or_mhc_choice = random.choice(["peptidseq_list", "Mhcseq_list", "both"])

            if peptide_or_mhc_choice == "peptidseq_list":
                logging.info("making peptide decoy")
                peptide_decoy = self._peptide_decoy_maker(peptidseq_list)
                mhc_decoy = Mhcseq_str

            elif peptide_or_mhc_choice == "Mhcseq_list":
                logging.info("making MHC decoy")
                mhc_decoy = self._peptide_decoy_maker(Mhcseq_list)
                peptide_decoy = peptidseq_str

            elif peptide_or_mhc_choice == "both":
                logging.info("making both peptide and MHC decoy")
                peptide_decoy = self._peptide_decoy_maker(peptidseq_list)
                mhc_decoy = self._peptide_decoy_maker(Mhcseq_list)
                
            if peptide_decoy is not None and mhc_decoy is not None:
                temp_set.add((peptide_decoy, mhc_decoy, 'False'))

        return list(temp_set)


# test = PeptideMhcPair(peptide_min=7, peptide_max=15, mhc_min=10, mhc_max=35, pairs=[
#                       ("A", "G"), ("Y", "E"), ("Q", "N"), ("P", "T"), ("M", "K"), ("F", "W")])
# peptidseq_str, peptidseq_list, Mhcseq_str, Mhcseq_list = test.peptide_maker()
# a = test._peptide_decoy_maker(['DM', 'A', 'TLDI', 'Q', 'HI', 'Y', 'QK'])

# c = test.decoy_output(peptidseq_str, peptidseq_list,
#                       Mhcseq_str, Mhcseq_list, decoy_count=10)
# print(len(c))
# print("aaaa",c)


class dataWrapper():
    def __init__(self, peptide_min, peptide_max, mhc_min, mhc_max, sample_count=1, decoy_count= 9, pairs = 6, binding_pos_count = 3):
        
        

        self.sample_count = sample_count
        self.decoy_count = decoy_count
        self.pairs =  self.random_pairs(pairs) if type(pairs) is int else pairs

        self.peptide_min = peptide_min
        self.peptide_max = peptide_max
        self.mhc_min = mhc_min
        self.mhc_max = mhc_max
        self.binding_pos_count = binding_pos_count
        
        
    def _min_max_checker(self, parts_len, min_len, max_len, mhc_or_peptide = None):
        
        """
        a funciton to check the pre_specified values for min and max length of peptide or mhc

        """

        peptide_lenghth = sum(parts_len) + self.binding_pos_count
        

        if peptide_lenghth < min_len or peptide_lenghth > max_len:
            raise ValueError(
                f'the sequence of {mhc_or_peptide} should have {min_len} minimum and {max_len} maximum length, {mhc_or_peptide} length adds up to {peptide_lenghth}')
            
    def random_pairs(self, length = 6):
        """
        length(int): number of random pairs for binding 
        """
        return list({tuple(random.choices(amino_acid_letter, k=2)) for letters in range(length)})

    def random_indexes_letter(self, length=4):
        """
        args:
            length(int): number of spots between binding indeces. for 3 binding sites we have 4 spots
        """
        numbers = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        return random.choices(numbers, k=length)

    def multiple_seq_maker(self, pairs=None,number_of_decoys = 99,  peptide_min=7, peptide_max=15,
                           mhc_min=10, mhc_max=35, peptide_part_len= None, mhc_part_len=None):

        if pairs is None:
            pairs = self.pairs
              
        peptid_mhc_obj = PeptideMhcPair( pairs = pairs)
        peptidseq_str, peptidseq_list, Mhcseq_str, Mhcseq_list = peptid_mhc_obj.peptide_maker(amino_pairs= pairs, peptide_part_len = peptide_part_len,
                                                                                             mhc_part_len = mhc_part_len )
        decoys = peptid_mhc_obj.decoy_output(
            peptidseq_str, peptidseq_list, Mhcseq_str, Mhcseq_list, decoy_count = self.decoy_count)

        return(peptidseq_str, Mhcseq_str,"True"),decoys


    def print_to_files(self, true_counts, peptide_indeces, mhc_indeces, tail_peptide_min, tail_peptide_max, tail_mhc_min, tail_mhc_max):
        true_seqs = []
        decoys = []

        for _ in range(true_counts):
            
            peptide_indeces = [peptide_indeces[0],peptide_indeces[1], peptide_indeces[2], random.randint(tail_peptide_min, tail_peptide_max)]
            mhc_indeces = [mhc_indeces[0],mhc_indeces[1], mhc_indeces[2], random.randint(tail_mhc_min, tail_mhc_max)]
            
            
            self._min_max_checker(peptide_indeces[:3]+[tail_peptide_max], self.peptide_min, self.peptide_max, mhc_or_peptide = "peptide")
            self._min_max_checker(mhc_indeces[:3]+[tail_mhc_max], self.mhc_min, self.mhc_max, mhc_or_peptide= "Mhc")
            
            seqs, deqs = self.multiple_seq_maker( number_of_decoys = self.decoy_count, peptide_part_len = peptide_indeces, mhc_part_len = mhc_indeces)
            true_seqs.append(seqs)
            decoys.append(deqs)
            
        peptide_decs = []
        mhc_decs = []
        labels_decs = []
        for sublist in decoys:
            for j in range(len(sublist[0])):
                peptide_decs.append(sublist[j][0])
                mhc_decs.append(sublist[j][1])
                labels_decs.append(sublist[j][2])

        df_decoys = pd.DataFrame(list(zip(peptide_decs, mhc_decs, labels_decs)), columns=[
                                "peptide", "mhc", "label"])
        df_not_decoys = pd.DataFrame(
            true_seqs, columns=["peptide", "mhc", "label"])
        
        print(len(df_not_decoys['peptide']))
        print(len(df_not_decoys['peptide'].unique()))
        print(len(df_not_decoys['mhc']))
        print(len(df_not_decoys['mhc'].unique()))
        
        # for i,j in zip(df_not_decoys['peptide'].values, df_not_decoys['mhc'].values):
        #     # print(i,j)
        #     for b,d in zip(df_decoys['peptide'].values, df_decoys['mhc'].values):
        #         if i == b and j == d:
    
        #             print(i, b)
        #             print(j, d)

        df_decoys.to_csv("decoys_dataset.csv")
        df_not_decoys.to_csv("true_dataset.csv")




def main():
    sample_num = args.sample_num
    decoys_ratio = args.decoys
    min_peptide_len = args.peptide_min
    max_peptide_len = args.peptide_max
    min_mhc_len = args.mhc_min
    max_mhc_len = args.mhc_max
    
    
    tail_peptide_min = 1
    tail_peptide_max = 3
    tail_mhc_min = 5
    tail_mhc_max = 20
    peptide_indeces = [2,4,3]  #for manual indexing
    mhc_indeces =  [7,4,3]    #for manual indexing
    # mhc_indeces = [7,5,3, 30]
    # pairs = [('K', 'V'), ('K', 'L'), ('C', 'M'), ('I', 'T'), ('P', 'W'), ('A', 'M'), ('E', 'D'), ('L', 'Q'), 
    #          ('S', 'N'), ('A', 'L'), ('H', 'T'), ('L', 'F'), ('F', 'M'), ('G', 'T'), ('F', 'G'), ('R', 'N'), ('D', 'M'), 
    #          ('E', 'S'), ('F', 'V'), ('I', 'F'), ('Y', 'V'), ('G', 'Q'), ('W', 'W'), ('I', 'H'), ('E', 'I'), ('F', 'P'),
    #          ('Q', 'H'), ('M', 'C'), ('C', 'A'), ('E', 'P'), ('K', 'D'), ('A', 'A'), ('D', 'E'), ('V', 'L'), ('R', 'V'),
    #          ('G', 'M'), ('E', 'Y'), ('G', 'I'), ('G', 'R'), ('L', 'E')]
    
    # pairs = pairs[:30]
    pairs = 40
    test2 = dataWrapper(min_peptide_len, max_peptide_len, min_mhc_len, max_mhc_len, sample_count = sample_num, decoy_count = decoys_ratio, pairs = pairs)
    print(test2.pairs)
    # s2 = test2.random_pairs(6)
    # a, b =test2.multiple_seq_maker(pairs=None,peptide_part_len= [2,5,3, random.randint(1, 4)],mhc_part_len = [7,5,3, random.randint(9, 15)])
    test2.print_to_files(sample_num, peptide_indeces, mhc_indeces, tail_peptide_min,tail_peptide_max,tail_mhc_min, tail_mhc_max)

main()

