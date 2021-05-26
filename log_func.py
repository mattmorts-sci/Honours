def log(log_data):
    '''
    A simple logging function. Whatever string argument is passed to the function it will append
    to the 'Honours_log.txt' file with a datetime stamp
    '''
    from datetime import datetime
    date_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')

    with open('Honours_log.txt', 'a') as log:
        log.write(date_time + '\t' + log_data + '\n')
