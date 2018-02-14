# For second calculation
class MultiRaceProbs(object):
    """Constructor
           Intialize an instance of ` MultiRaceProbs` with a a total probability 
           of the comb of races.,best haplotype and best probability.
           """
    def __init__(self, prob_comb, best_hapl, best_prob):
        self.prob_comb = prob_comb
        self.best_hapl = best_hapl
        self.best_prob = best_prob

    def update_prob(self, new_prob):
        self.prob_comb += new_prob

    def set_best_hapl(self, best_hapl, best_prob):
        # Change the best haplotype to a new
        self.best_hapl = best_hapl
        self.best_prob = best_prob

    def check_max(self, prob):
        # Check if the new prob is better
        return self.best_prob < prob

    def get_best_hapl(self):
        return self.best_hapl

    def get_prob_comb(self):
        return self.prob_comb


class GetBestOption(object):
    def __init__(self, file_name):
        """Constructor
        Intialize an instance of `GetBestOption` with a file name
        """
        self.data = file_name

    # Get the best option
    def get_best_options(self):

        # Open the files
        file = open(self.data, 'r')
        fout1_name = str(self.data)+'_best_options_a'
        fout1 = open(fout1_name, 'w')
        fout2_name =  str(self.data)+'_best_options_b'
        fout2 = open(fout2_name, 'w')
        rows = file.readlines()
        pre_patient = ""
        data_patient=[]

        # Go over the data
        for i in range(len(rows)):
            # If this is empty row
            print(rows[i])
            if rows[i] == '\n':
                continue
            #Rows[i]=ID, Haplotype1, Probability 1, Race 1, Haplotype2, Probability 2, Race 2
            # Create the format
            rows[i] = rows[i].strip('\n')
            data_list = rows[i].split(',')

            #data_list[0]= id
            # If this is a new patient
            if pre_patient != data_list[0] and pre_patient != "":

                # Find the options
                option_a = ','.join(self.get_race_hap_a(data_patient))
                option_b = ','.join(self.get_race_hap_b(data_patient))

                # Write the files
                fout2.write((option_b)+'\n')
                fout1.write((option_a)+'\n')

                # Init to the next patient
                data_patient = []

            data_patient.append(data_list)

            # data_list[0]= id
            pre_patient = data_list[0]

        ### To the last patient
        # Find the options
        option_a = ','.join(self.get_race_hap_a(data_patient))
        option_b = ','.join(self.get_race_hap_b(data_patient))

        # Write the files
        fout2.write((option_b) + '\n')
        fout1.write((option_a) + '\n')

        file.close()
        fout1.close()
        fout2.close()

    # First
    # Input:list of all the lists of ID, Haplotype1, Probability 1, Race 1, Haplotype2, Probability 2, Race 2
    # # of one patient
    def get_race_hap_a(self, probs):
        # Init
        max = 0
        best = None
        # Go over the rows
        for i in range(len(probs)):

            # The data
            data = probs[i]

            # data[2] = prob1 and data[5]=prob2
            prob1 = float(data[2])
            prob2 = float(data[5])

            # Calculate the prob of the row
            current = prob1 * prob2

            # If we found new mex
            if current > max:
                # Init the max
                max = current
                best = i
        # Return the best row
        return (probs[best])

    # Input:list of all the lists of ID, Haplotype1, Probability 1, Race 1, Haplotype2, Probability 2, Race 2
    # # of one patient
    def get_race_hap_b(self, probs):
        # Init
        dict_pops = {}
		
		# Go over all the pairs of hapl and their data
        for i in range(len(probs)):
            # The data
            data = probs[i]

            # Key to the dict of all the pops
            key = ""
            # data[2] = prob1 and data[5]=prob2
            prob1 = float(data[2])
            prob2 = float(data[5])

            # Calculate the prob of the row
            current_prob = prob1 * prob2

            # Set the races
            race1 = data[3]
            race2 = data[6]

            # Find the race combination for the dict
            # The format of a  key is race1race2
            comb_data = self.find_pops_comb(race1, race2,
                                            dict_pops)

			# If this is the first time we see this comb
            if not comb_data:
                # Set a new key
                key = race1+race2,
                # Init the data to be the best we found and as total prob.
                dict_pops[key] = MultiRaceProbs(current_prob, data, current_prob)
            else:
				# Check if we found a better row for this population comb
                if comb_data.check_max(current_prob):
				#Change the data to better one we found
                  comb_data.set_best_hapl(data, current_prob)
				# Add the prob to the probs we had
                comb_data.update_prob(current_prob)
        return self.find_max(dict_pops)

    # We want to find if the combination is in out dict
    # We dont know the order- so we try the right one.
    # We assume each pair of pops appear in one way only
    def find_pops_comb(self,key1,key2,dict_pops):
        # Init
        key=""
        comb_data= None
        # Find the current comb
        if (key1 + key2) in dict_pops.keys():
            key = key1 + key2
        elif (key2 + key1) in dict_pops.keys():
            key = key2 + key1

        # If it's not a first time we meet this comb
        if key != "":
            comb_data = dict_pops[key]

        return comb_data


    def find_max(self, dict_pops):
        best_peir = []
        best_prob = 0
        # Go over all the dict
        for key in dict_pops:
            if dict_pops[key].get_prob_comb() > best_prob:
                # Set the best current prob and the pair.
                best_prob = dict_pops[key].get_prob_comb()
                best_peir = dict_pops[key].get_best_hapl()
        return best_peir
