import pandas as pd
from log_func import log
from append_file import append_file


def onekp(file_list):  # Creates the function for parsing phytozome output files
    '''
    Takes in the fasta and details file names without the file extentions
    Three args, fasta file name, details file name, SAL variant
    All files must be .txt and do not provide the file extentions

    More detailed documentation is provided in the phyto_func version.
    '''

    new_filelist = []
    index_filelist = []
    blinded_filelist = []
    mod_filelist = []

# Iterates the rest of the function through each set of files in flist
    with open(f'{file_list}.txt') as flist:
        for i in flist:
            ls = ()
            ls = i.split()
            fastaFile = ls[0]
            detailsFile = ls[1]
            variant = ls[2]

# Open unid file, assign each number as element in a list
            with open('unique_id.txt', 'r') as file:
                unique_id = file.read().split(' ')

            unique_id.pop()  # removes a ' ' element from the end of the file
            unid_init_len = len(unique_id)  # Checks initial length of the unique_id file

# Open the fasta file, assign unid by popping from list to each sequence header, save as new file
            ofile = open(f'{fastaFile}_new.txt', 'w')
            ofile_2 = open(f'{fastaFile}_index.txt', 'w')
            ofile_3 = open(f'{fastaFile}_blinded.txt', 'w')

            ofile_2.write('Onekp_index_id' + '\t' + 'Scaffold' + '\t' + 'Unid' + '\t' +
                          'SAL_variant' + '\t' + 'Datasource' + '\t' + 'Header' + '\t' +
                          'Subject Seq-id' + '\n')

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
                        i = line.split('_')
                        key = line.split(' ')[1]
                        temp = i[0].split('|')
                        unid = temp[0][1:-4]
                        onekp_index_id = temp[2]
                        scaffold = i[2].split(' ')[0]
                        line = line.strip()
                        ofile_2.write(str(onekp_index_id) + '\t' + str(scaffold) + '\t' + str(unid)
                                      + '\t' + variant + '\t' + 'Onekp' + '\t' + str(line) + '\t' +
                                      str(key) + '\n')
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

            index_file = pd.read_csv(f'{fastaFile}_index.txt', sep='\t')

# Selects these four columns
            index_file = index_file[['Unid', 'Onekp_index_id',
                                     'Scaffold', 'SAL_variant', 'Subject Seq-id']]

# Merges both dataframes adding in the unnique id column to the details file
            details = pd.merge(index_file, details, how='inner', on='Subject Seq-id')
            del details['Subject Seq-id']

            details.to_csv(f'{detailsFile}_mod.txt', sep='\t', index=False)

            log(f'{len(details)} sequence details written to {detailsFile}_mod.txt')

            new_filelist.append(f'{fastaFile}_new')
            index_filelist.append(f'{fastaFile}_index')
            blinded_filelist.append(f'{fastaFile}_blinded')
            mod_filelist.append(f'{detailsFile}_mod')

    append_file('onekp_new_all', new_filelist)
    log(f'{len(new_filelist)} files were appended to onekp_new_all.txt')
    append_file('onekp_index_all', index_filelist)
    log(f'{len(index_filelist)} files were appended to onekp_index_all.txt')
    append_file('onekp_blinded_all', blinded_filelist)
    log(f'{len(blinded_filelist)} files were appended to onekp_blinded_all.txt')
    append_file('onekp_details_mod_all', mod_filelist)
    log(f'{len(mod_filelist)} files were appended to onekp_details_mod_all.txt')

# Adds the sequence to the details file

    detailsCor = pd.read_csv('onekp_index_all.txt', sep='\t')
    detailsCor['Unid'] = detailsCor['Unid'].apply(str)

    k = ''
    v = ''
    seq_dict = {}
    with open('onekp_blinded_all.txt') as seqs:
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

    indexSeq.to_csv('onekp_master.txt', sep='\t', index=False)
    log(f'onekp_master.txt was created from onekp_index_all.txt and {len(indexSeq)} were added')
