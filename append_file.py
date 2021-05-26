def append_file(name, lst):
    '''
    This function takes a list of files and appends them to a new txt file.
    The function takes two arguments, the name (without the file extension) of the new file
    and a list of files without the file extension
    '''
    from log_func import log

    with open(f'{name}.txt', 'a') as op_file:

        for i in lst:
            count = 0
            with open(f'{i}.txt', 'r') as temp_file:
                for line in temp_file:
                    op_file.write(line)
                    count += 1
                log(f'{count} lines were added to {name}.txt')

    op_file.close()
