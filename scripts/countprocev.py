import sys 

# calcola il numero di eventi processati
# dal log file di crab che si ottiene col
# seguente comando
# grep Events [crab_directory]/res/*.stdout > evtnum.log

in_file = open("evtnum.log","r")
i = 0
processed_events = 0
while 1:
    text = in_file.readline()
    if len(text) == 0:
        break
#    print text
    text.rstrip()
    subpieces = text.split('=')
    num_of_pieces = len(subpieces)
#    print num_of_pieces
    my_piece = subpieces[1]
#    print my_piece
    carnum = my_piece.find("passed")
#    print carnum
    temp_processed_events = my_piece[0:carnum]
#    print temp_processed_events
#    processed_events.rstrip
#    print processed_events
    processed_events += int(temp_processed_events)
    #    files.append(text)
    i = i + 1
#    print processed_events
in_file.close()
print i
print processed_events

