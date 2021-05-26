# Unique id generation. This generates 89,999 unique ids,
# shuffles them and allows them to be popped off for assignment

# This has already been run, DO NOT RE-RUN it or it will overwrite the file

# randomly shuffle a sequence
from random import shuffle
from random import seed

# seed to provide reproducibility
seed(1)

# prepare a sequence
ids = [i for i in range(10000, 99999)]

# randomly shuffle the sequence
shuffle(ids)

ofile = open('unique_id.txt', 'w')

for i in ids:
    ofile.write(str(i) + ' ')

ofile.close()
