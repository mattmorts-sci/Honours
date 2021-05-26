# Function for assinging unique ids to phytozome blast output
# Generates four new files: '{filename}_details_mod.txt' the details file with the unique id added and the PACid modified to remove 'PAC:'
#                           '{filename}_new.fasta' the fasta file with unique id added to start of the header
#                           '{filename}_blinded.fasta' the fast file with the header replaced with only the unique id
#                           '{filename}_index.txt' an index file with the unique id, PACid and header

import pandas as pd
from log_func import log
from append_file import append_file
from species_fix import sp_fix


def phyto(file_list):  # Creates the function for parsing phytozome output files
    '''
    Takes in the fasta and details file names without the file extentions
    Three args, fasta file name, details file name, SAL variant
    All files must be .txt and do not provide the file extentions
    '''

# variables for list of files to append
    new_filelist = []
    index_filelist = []
    blinded_filelist = []
    mod_filelist = []

# Takes a list of files to apply the function to and iterates the rest of the function through
# each set of files in flist
    with open(f'{file_list}.txt') as flist:
        for i in flist:
            ls = ()
            ls = i.split()
            fastaFile = ls[0]
            detailsFile = ls[1]
            variant = ls[2]

# Opens the unid file (file with list of unique ids), assign each id from the file as an element
# in a list called unique_id
            with open('unique_id.txt', 'r') as file:
                unique_id = file.read().split(' ')

            unique_id.pop()  # removes a ' ' element from the end of the file
            unid_init_len = len(unique_id)  # Checks initial length of the unique_id file

# Open the fasta file, assign unique id by popping from list to each sequence header, saves as a
# new file
            ofile = open(f'{fastaFile}_new.txt', 'w')  # Same fasta headers but prefixed with
            # unique id
            ofile_2 = open(f'{fastaFile}_index.txt', 'w')  # Index file with the fields below
            ofile_2.write('PACid' + '\t' + 'Unid' + '\t' + 'SAL_variant' + '\t' + 'Datasource' +
                          '\t' + 'Header' + '\n')
            ofile_3 = open(f'{fastaFile}_blinded.txt', 'w')  # Fasta file with headers only as
            # unique id

# Adds the unique id to the start of the header in the 'new' fasta file
            with open(f'{fastaFile}.txt', 'r') as myfasta:

                count = 0

                for line in myfasta:
                    if line.startswith('>'):
                        if count == 0:
                            ofile.write('>' + unique_id.pop() + ' ' + line.lstrip('>'))
                        else:
                            ofile.write('\n>' + unique_id.pop() + ' ' + line.lstrip('>'))
                        count += 1
                    else:
                        if line.endswith('*\n'):
                            ofile.write(line.strip('*\n'))
                        else:
                            ofile.write(line.strip('\n'))
                        count += 1

                log(f'{fastaFile}_new.txt was created from {fastaFile}.txt and {count} lines written')

            ofile.close()

# Writes the index file
            with open(f'{fastaFile}_new.txt', 'r') as myNewfasta:

                count = 1

                for line in myNewfasta:
                    if line.startswith('>'):
                        i = line.split()
                        pacid = i[-1][5:-1]
                        unid = i[0][1:]
                        header = i[1:-2]
                        hd = ''
                        for i in header:
                            hd = hd + i + ' '
                        ofile_2.write(str(pacid) + '\t' + str(unid) + '\t' + variant + '\t' +
                                      'Phytozome' + '\t' + str(hd) + '\n')
                        count += 1
                    else:
                        continue

                log(f'{fastaFile}_index.txt was created from {fastaFile}_new.txt and {count} '
                    + 'lines were written')

# Creates a new fasta file with header = unid only
                myNewfasta.seek(0)

                count = 0

                for line in myNewfasta:
                    if line.startswith('>'):
                        ofile_3.write(line[:6] + '\n')
                    else:
                        ofile_3.write(line)
                    count += 1

                log(f'{fastaFile}_blinded.txt was created from {fastaFile}_new.txt and {count}'
                    + f' lines were written')

            ofile_2.close()
            ofile_3.close()

# Overwrites the unique_id.txt file minus the unique values that were used
            unid_final_len = len(unique_id)

            log(f'{unid_init_len - unid_final_len} unique ids were assigned from unique_id.txt, ' +
                f'file has {unid_final_len} ids remaining')

            ofile = open('unique_id.txt', 'w')

            for i in unique_id:
                ofile.write(str(i) + ' ')

            ofile.close()

# Open details file, in pandas, repalces the 'PAC:xxxxxxxx' with 'xxxxxxxx' and add value (unid)
# from dict as a new column based off key (pac id), save as new file

            details = pd.read_csv(f'{detailsFile}.txt', sep='\t')

            # Removes'PAC:' from column values
            details['PACid'] = details['PACid'].str.slice_replace(0, 4)

            index_file = pd.read_csv(f'{fastaFile}_index.txt', sep='\t')

            index_file = index_file[['Unid', 'PACid', 'SAL_variant']]  # Selects these three columns

            # Coverts this column type to str to match other dataframe
            index_file['PACid'] = index_file['PACid'].apply(str)

            # Merges both dataframes adding in the unique id column to the details file
            details = pd.merge(index_file, details, how='inner', on='PACid')

            details.to_csv(f'{detailsFile}_mod.txt', sep='\t', index=False)

            log(f'{len(details)} sequence details written to {detailsFile}_mod.txt')

# Adds each of the files created in the function above to gloabl lists, these lists are then used
# to append each file to 'phytozome all' files
            new_filelist.append(f'{fastaFile}_new')
            index_filelist.append(f'{fastaFile}_index')
            blinded_filelist.append(f'{fastaFile}_blinded')
            mod_filelist.append(f'{detailsFile}_mod')

    append_file('phytozome_new_all', new_filelist)
    log(f'{len(new_filelist)} files were appended to phytozome_new_all.txt')
    append_file('phytozome_index_all', index_filelist)
    log(f'{len(index_filelist)} files were appended to phytozome_index_all.txt')
    append_file('phytozome_blinded_all', blinded_filelist)
    log(f'{len(blinded_filelist)} files were appended to phytozome_blinded_all.txt')
    append_file('phytozome_details_mod_all', mod_filelist)
    log(f'{len(mod_filelist)} files were appended to phytozome_details_mod_all.txt')

# Runs the species fix function (phytozome only)
    sp_fix()

# Adds the each sequence to the details file against the corresponding row based off the unique id

    detailsCor = pd.read_csv('phytozome_index_all.txt', sep='\t')
    detailsCor['Unid'] = detailsCor['Unid'].apply(str)

    k = ''
    v = ''
    seq_dict = {}

    # Creates a dictionary where the unique id is added as the key and the sequence as the value
    with open('phytozome_blinded_all.txt') as seqs:
        for line in seqs:
            if line.startswith('>'):
                k = line.lstrip('>').strip('\n')
            else:
                v = line.strip('\n')
            seq_dict.update({k: v})

    seqDF = pd.DataFrame.from_dict(seq_dict, orient='index')
    seqDF.reset_index(level=0, inplace=True)
    seqDF.columns = ['Unid', 'Sequence']
    indexSeq = detailsCor.merge(seqDF, on=['Unid'], how='inner').fillna('NaN')

    indexSeq.to_csv('phytozome_master.txt', sep='\t', index=False)
    log(
        f'phytozome_master.txt was created from phytozome_index_all.txt and {len(indexSeq)} were added')
