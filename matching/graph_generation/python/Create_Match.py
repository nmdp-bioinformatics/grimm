
'''
This code  generate for each donor all the 9/10,8/10,2/10,1/2. If you don't want
to compute the 8/10 options- set flag_eight to False in
"from_ten_to_nine_options" function.
'''

# This library contains functions of csv files,allows us to read and write
import csv

# This library contains functions to go effectively over a list
import itertools

# Save a dictionary into a pickle file.
import pickle

# This function find if the file exist
import os.path

import datetime


# Virable to count the nodes
COUNTIDNODE = 1


def create_dictionary():
    # Init the path of the wanted file
    path = 'mytest.txt_out_best_options_a'

    # Init
    flag = False
    global COUNTIDNODE
    DICT = {}

    # The min freq we want to consider
    min_freq = 0.00000001
    # Open the data
    data = open(path)
    reader = csv.reader(data)

    # If this is no the first time we create the dictionary
    # Else- the dict is empty and the COUNTIDNODE=0
    if os.path.exists('save_dict_graph.p'):
        # Load this dict
        DICT = pickle.load(open("save_dict_graph.p", "rb"))
        # Init the ID of the nodes
        COUNTIDNODE = max(DICT.values())+1

    # IF THIS IS THE FIRST TIME WE RUN THE PROGRAM
    else:
        # Init the dict
        DICT = {}

        # Init the ID of the nodes
        COUNTIDNODE = 1

    nodes_list = []
    edges_list = []
    donor_2_ten_list = []

    # Go over all the csv file
    # The format of each row is :ID, Haplotype1, Probability 1, Race 1, Haplotype2, Probability 2, Race 2
    for row in reader:
        # If it is a header-continue
        if not flag:
            flag = True
            continue
        # If this is an empty row
        if row == [] :
            continue
        # Comp the freq of the full gen
         # row[3] = freq1 and row[6] =freq2
        freq= float(row[3])*float(row[6])

        # Skip the line if the freq no good
        if freq<min_freq:
            continue

        # Add the data as the dna of 10 and the id to the dictionary if this is new
        if row[0] not in DICT.keys():
            DICT.update({row[0]: COUNTIDNODE})

            # Add the node
            nodes_list.append([COUNTIDNODE, row[0], 'Donor'])

            # Precede the node id!!!
            COUNTIDNODE += 1

        # Create the sorted data out of the sting
        # Convert to array
        # row[1] = hapl1 and row[4] = hapl2
        list_row1 = row[1].split('~')
        list_row2 = row[4].split('~')

        # Sort the two list into one list
        # Example: ['A*01:01', 'A*03:01', 'B*07:02', 'B*08:01', 'C*07:01', 'C*07:02', 'DQB1*02:01', 'DQB1*06:02', 'DRB1*03:01', 'DRB1*15:01']
        list_row_sorted = (sorted(list_row1+list_row2))

        # Check and correct null alleles
        list_row_sorted_fix = null_alleles_fix(list_row_sorted)

        # Convert to string
        sorted_string = ','.join(str(e) for e in list_row_sorted_fix)

        # Check if the data in the dictionary
        if sorted_string not in DICT.keys():
            DICT.update({sorted_string: COUNTIDNODE})

            # Add the data to the list of nodes
            nodes_list.append([COUNTIDNODE,(sorted_string) ,"MUUG"])

            # Precede the node id!!! MUST DO IT AT EVERY ADDING TO THE NODE!!!
            COUNTIDNODE += 1

            # Create all the 9-10 links and add it to the dictio   nary
            # In this function the countNode may grow!!!!!!!
            # If you don't want the nine-eight options- send false
            from_ten_to_nine_options(list_row_sorted, sorted_string,
                                     DICT, nodes_list, edges_list, True)

            # Create all the links with allels and add it to the dictionary
            # In this function the countNode may grow!!!!!!!
            from_ten_to_allels(list_row_sorted, sorted_string,  DICT,
                               nodes_list, edges_list)

        # Add to the list of donor to gen , contains a frequency
        donor_2_ten_list.append([DICT[sorted_string], DICT[row[0]],
                                 "DONOR_TO_MUUG", freq])

    # Close the file
    data.close()

    write_results(nodes_list, edges_list, donor_2_ten_list)
    # Save the dict to next usage
    pickle.dump(DICT, open("save_dict_graph.p", "wb"))
    print(COUNTIDNODE)


def write_results(nodes_list, edges_list,donor_2_ten_list):
    time_str = str(datetime.datetime.now().strftime("%m_%d_%y_%H_%M"))
    with open('Nodes_Data_'+time_str+'.csv', 'w') as newFile:
        # Create the fields names in the new file
        fieldnames = [':nodeId:ID(node)', 'name', 'type:LABEL']
        writer = csv.DictWriter(newFile, fieldnames=fieldnames)
        writer.writeheader()
        writer = csv.writer(newFile)
        writer.writerows(nodes_list)
    # Open a file to write the links between patient and G-10
    with open('Edges_Data_Donor2MUUG_'+time_str+'.csv', 'w') as newFile:
        # Create the fields names in the new file
        fieldnames = [':END_ID(node)', ':START_ID(node)', ':TYPE', 'frequency:FLOAT']
        writer = csv.DictWriter(newFile, fieldnames=fieldnames)
        writer.writeheader()
        writer = csv.writer(newFile)
        writer.writerows(donor_2_ten_list)

    # Open a file to write the links between all others
    with open('Edges_Data_'+time_str+'.csv', 'w') as newFile:
        # Create the fields names in the new file
        fieldnames = [':END_ID(node)', ':START_ID(node)', ':TYPE']
        writer = csv.DictWriter(newFile, fieldnames=fieldnames)
        writer.writeheader()
        writer = csv.writer(newFile)
        writer.writerows(edges_list)

# Check and correct the null alleles
def null_alleles_fix(data):
    # Go over the list
    for i in range(0, 10):
        # If the last char is "N" - the null allele
        string_allele = str(data[i])
        if string_allele[-1:] == "N":
            # If this is a even place-take the next
            if i % 2 == 0:
                data[i] = data[i+1]
            else:
                # Take the previous
                data[i] = data[i-1]
    return data


# Check all the ten to nine options
# flag_eight-If you don't want the nine-eight options- set it as False
def from_ten_to_nine_options(data, string_ten, DICT, nodes_list, edges_list,flag_eight):
    global COUNTIDNODE
    string_data = ""

    # Go over the array
    for i in itertools.combinations(data, 9):
        string_data = ','.join(i)

        if string_data not in DICT.keys():
            # Add to the dict
            DICT[string_data] = COUNTIDNODE

            # Add to nodes list
            nodes_list.append([COUNTIDNODE, string_data, "MUUG_1"])

            # Precede the node id!!!!!
            COUNTIDNODE += 1

            # Find all nine-to eight -only if flag_eight is true
            if flag_eight:
                from_nine_to_eight(i, string_data,  DICT, nodes_list, edges_list)

        edges_list.append([DICT[string_data], DICT[string_ten], "MUUG_TO_MUUG_1"])
    return


# This function is find all the 8 to 10 options and add it!
def from_nine_to_eight(data, string_nine, DICT,nodes_list, edges_list ):
    global COUNTIDNODE
    string_data = ""

    # Go over the array
    for i in itertools.combinations(data, 8):
        string_data = ','.join(i)

        if string_data not in DICT.keys():
            # Add to the dict
            DICT[string_data] = COUNTIDNODE

            # Add to node list
            nodes_list.append([COUNTIDNODE, string_data, "MUUG_2"])
            # Precede the node id!!!!!
            COUNTIDNODE += 1

        edges_list.append([DICT[string_data], DICT[string_nine], "MUUG_1_TO_MUUG_2"])


# Find all the allels
def from_ten_to_allels(data, string_ten, DICT, nodes_list, edges_list):
    # Init
    global COUNTIDNODE
    string_data = ""
    arr = []
    place = 0
    end = place + 2

    # While it is not the end of the string
    while end <= 10:
        # Take the allele
        arr = list(itertools.islice(data, place, end))

        # Create the string
        string_data = ','.join(arr)

        # If the allele not in the dict
        if string_data not in DICT.keys():
            DICT[string_data] = COUNTIDNODE

            nodes_list.append([COUNTIDNODE, string_data, "SLUG"])
            # Precede the node id!!!!!
            COUNTIDNODE += 1
            from_two_to_one(arr, string_data, DICT,nodes_list, edges_list)


        edges_list.append([DICT[string_data], DICT[string_ten], "MUUG_TO_SLUG"])

        # Precede the interval
        place += 2
        end += 2


def from_two_to_one(data, allels, DICT,nodes_list, edges_list):
    # Init
    global COUNTIDNODE
    string_data = ""
    flag = 0

    # Go over the allele
    for i in range(2):
        # Create string
        string_data = str(data[i])

        # If it is in the dict
        if string_data not in DICT.keys():
            DICT[string_data] = COUNTIDNODE

            nodes_list.append([COUNTIDNODE, string_data, "Allele"])

            # Precede the node id!!!!!
            COUNTIDNODE += 1

        edges_list.append([DICT[string_data], DICT[allels], "SLUG_2_ALLELE"])

        # If this is the same allele
        if data[0] == data[1]:
            break


#Main function
def main():
    # Run the function

    create_dictionary()
    return


if __name__ == '__main__':
    main()
