import sys
import numpy as np
import re

def main(database_file=None, query_file=None):
    """
    usage:

        python ./BLAST_Algorithm.py <databases_file> <query_sequence_file>
    """   

    # read databases and their names in lists
    databases = []
    database_names = []
    database_count = 0
    database = ""
    file1 = open(database_file,'r')
    for line in file1:
        if (">" in line):
            database_count = database_count + 1
            database_names.append(line)
            if (database_count != 1):
                databases.append(database)
                database = ""
        else:
            database = database + line.split('\n')[0]
    

    # read query sequence into string
    query_sequence = ""
    file2 = open(query_file,'r')
    for line in file2:
        if (">"  not in line):
            query_sequence = query_sequence + line.split('\n')[0]
    
    # stores the ktuples found in the query sequence and their locations 
    query_ktuples = []
    query_ktuples_locations = []
    # using a k size of 28
    k = 28
    q_length = len(query_sequence)
    for q in range(0, q_length - (k - 1)):
        ktuple = query_sequence[q:q+k-1]
        if (ktuple not in query_ktuples):
            loc_list = [m.start() for m in re.finditer(ktuple, query_sequence)]
            loc = ""
            for l in loc_list:
                if (loc != ""):
                    loc = loc + ";" + str(l)
                else:
                    loc = loc + str(l)
            query_ktuples.append(ktuple)
            query_ktuples_locations.append(loc)

    # stores the ktuples occurring in the databases and their locations
    database_ktuples = []
    db_ktuples_locations = []
    for db in databases:
        db_length = len(db)
        db_num = databases.index(db)
        for s in range(0, db_length - (k - 1)):
            ktuple = db[s:s+k-1]
            if ktuple in query_ktuples:
                loc = str(db_num) + "," + str(s)
                if (ktuple in database_ktuples):
                    x = database_ktuples.index(ktuple)
                    db_ktuples_locations[x] = db_ktuples_locations[x] + ";" + loc
                if (ktuple not in database_ktuples):
                    database_ktuples.append(ktuple)
                    db_ktuples_locations.append(loc)

    q_matches = []
    db_matches = []
    # for all the kmers found in the query sequence:
    for t in range(0, len(query_ktuples)):
        kmer = query_ktuples[t]
        # if the kmer is also found in the database:
        if kmer in database_ktuples:
            # the kmer positions in the query and database locations list are found and the locations of the matching kmers are stored 
            # in the q_matches and db_matches lists
            q_match_loc = query_ktuples_locations[t]
            q_matches.append(q_match_loc)
            i = database_ktuples.index(kmer)
            db_match_loc = db_ktuples_locations[i] 
            db_matches.append(db_match_loc)


    # creates empty lists to store positions that have already been aligned so the same regions are not aligned again
    already_aligned_q_start = []
    already_aligned_q_end = []
    already_aligned_db_seq = []
    already_aligned_db_seq_start = []
    already_aligned_db_seq_end = []
    # iterates through all the kmers in the query sequence that matched the database
    for pos in range(0, len(q_matches)):
        q_locations = q_matches[pos]
        db_locations = db_matches[pos]
        db_loc_list = db_locations.split(";")
        q_loc_list = q_locations.split(";")
        # iterates through all the positions the kmer was found in the the query sequence
        for q_loc in q_loc_list:
            q_loc = int(q_loc)
            # iterates through all the positions the kmer was found in the databases
            for db_loc in db_loc_list:
                best_q_seq = ""
                best_db_seq = ""
                best_score = 0
                db_num = int(db_loc.split(',')[0])
                db_start = int(db_loc.split(',')[1])
                # only aligns the two sequences if they have not already been aligned to each other
                part_of_earlier_alignment = "False"
                if (len(already_aligned_q_start) != 0):
                    for at in range(0, len(already_aligned_q_start)):
                        if (q_loc >= already_aligned_q_start[at]) and (already_aligned_q_end[at] > q_loc) and (db_num == already_aligned_db_seq[at]) and (db_start >= already_aligned_db_seq_start[at]) and (already_aligned_db_seq_end[at] > db_start):
                            part_of_earlier_alignment = "True"      
                if (part_of_earlier_alignment == "False"):
                    db = databases[db_num]
                    db_name = database_names[db_num]
                    # seeds the search where the database and query match
                    q_seq = query_sequence[q_loc:q_loc+k]
                    db_seq = db[db_start:db_start+k]
                    db_length = len(db)
                    # aligns the database and query and gets the alignment score
                    align_score = align_seqs(db_seq, q_seq)
                    best_score = align_score
                    ex_r = 0
                    end_met = "False"
                    # extends search to the right until the score is less than the original score
                    while (align_score >= best_score) and (end_met == "False"):
                        ex_r = ex_r + 1 
                        q_seq = query_sequence[q_loc:q_loc+k+ex_r]
                        db_seq = db[db_start:db_start+k+ex_r]
                        # aligns the extended sequences
                        align_score = align_seqs(db_seq, q_seq)
                        # if the alignment score is the best score so far 
                        if align_score >= best_score:
                            # the alignment score is the new high score
                            best_score = align_score
                            # store the query and database sequences to generate the consensus later
                            best_q_seq = q_seq
                            best_db_seq = db_seq 
                        if (len(q_seq) >= (q_length - q_loc)) or (len(q_seq) >= (db_length - db_start)):
                            end_met = "True"
                    if (best_score > align_score):
                        ex_r = ex_r - 1
                    ex_l = 0
                    align_score = best_score
                    # extends the search to the left until the score is less than the original score
                    if (q_loc != 0) and (db_start != 0):
                        end_met = "False"
                        while (align_score >= best_score) and (end_met == "False"):
                            ex_l = ex_l + 1
                            q_seq = query_sequence[q_loc-ex_l:q_loc+k+ex_r]
                            db_seq = db[db_start-ex_l:db_start+k+ex_r]
                            # aligns the extends sequences
                            align_score = align_seqs(db_seq, q_seq)
                            # if the alignment score is the best score so far 
                            if align_score >= best_score:
                                best_score = align_score
                                best_q_seq = q_seq
                                best_db_seq = db_seq
                            if (len(q_seq) >= (q_length - q_loc + ex_l)) and (len(q_seq) >= (db_length - db_start + ex_l)):
                                end_met = "True"
                        if (best_score > align_score):
                            ex_l = ex_l - 1
                    # appending final positions to lists for already aligned positions
                    q_final_start = q_loc - ex_l
                    q_final_end = q_loc + k + ex_r
                    db_final_start = db_start - ex_l
                    db_final_end = db_start + k + ex_r
                    already_aligned_q_start.append(q_final_start)
                    already_aligned_q_end.append(q_final_end)
                    already_aligned_db_seq.append(db_num)
                    already_aligned_db_seq_start.append(db_final_start)
                    already_aligned_db_seq_end.append(db_final_end)
                    # get a consensus sequence for the best alignment
                    print("Match at position ", q_loc, " in the query sequence and position ", db_start, "in database sequence ", db_name)
                    get_consensus(best_db_seq, best_q_seq)          

def align_seqs(db_seq, q_seq):
    match = 1
    mismatch = -2
    open_gap = -2
    ex_gap = -2
    # find the lengths of the sequences and uses them to intialize the score and traceback matrices
    matrix_length = len(db_seq) + 1
    matrix_height = len(q_seq) + 1
    score_matrix = np.empty((matrix_height, matrix_length))
    traceback_matrix = np.chararray((matrix_height, matrix_length))
    # sets the score in the left-top most cell of the matrix to 0 and the traceback to "s" for start
    score_matrix[0][0] = 0
    traceback_matrix[0][0] = "s"
    # uses the gap penalties to fill the top row of the matrices
    for x in range(1, matrix_length):
        score_matrix[0][x] = 0
        traceback_matrix[0][x] = "h"      
    # uses the gap penalties to fill the leftmost column of the matrices
    for y in range(1, matrix_height):
        score_matrix[y][0] = 0
        traceback_matrix[y][0] = "v"       
    # for each row in the matrix
    for i in range(1, matrix_height): 
        # for each cell in the row
        for j in range(1, matrix_length):
            vertical = 0
            horizontal = 0
            # if the cell above the current cell is being traced back vertically, then the vertical score for this cell is found by extending the gap
            if (traceback_matrix[i - 1][j] == b'v'):
                vertical = score_matrix[i - 1][j] + ex_gap
            # if the cell above the current cell is not being traced back vertically, then the vertical score for this cell is found by opening a gap
            else:
                vertical = score_matrix[i - 1][j] + open_gap
            if (0 > vertical):
                vertical = 0
            # if the cell left of the current cell is being traced back horizontally, then the horizontal score for this cell is found by 
            # extending the gap
            if (traceback_matrix[i][j - 1] == b'h'):
                horizontal = score_matrix[i][j - 1] + ex_gap
            # if the cell left of the current cell is not being traced back horizontally, then the horizontal score for this cell is found by 
            # opening a gap
            else:
                horizontal = score_matrix[i][j - 1] + open_gap
            if (0 > horizontal):
                horizontal = 0
            diagonal = 0
            # if the sequences match at this position then the diagonal score for this cell is found by adding the match score to the score 
            # of the diagonal cell
            if (db_seq[j - 1] == q_seq[i - 1]):
                diagonal = score_matrix[i - 1][j - 1] + match
            # if the sequences do not match at this position then the diagonal score for this cell is found by adding the mismatch score to the score 
            # of the diagonal cell
            elif (db_seq[j - 1] != db_seq[i - 1]):
                diagonal = score_matrix[i - 1][j - 1] + mismatch
            if (0 > diagonal):
                diagonal = 0
            # finds which score is the highest and then enters the highest score and the corresponding trace to the matrices
            max = vertical
            trace = "v"
            if (horizontal > max):
                max = horizontal
                trace = "h"
            if (diagonal >= max):
                max = diagonal
                trace = "d"
            score_matrix[i][j] = max
            traceback_matrix[i][j] = trace

    # the score for the alignment is the highest score in the matrix
    high_score = 0
    for o in range(0, matrix_height): 
        for p in range(0, matrix_length):
            if (score_matrix[o][p] >= high_score):
                high_score = score_matrix[o][p]
    align_score = high_score
    return(align_score)

def get_consensus(best_db_seq, best_q_seq):
    # aligning the sequences again
    match = 1
    mismatch = -2
    open_gap = -2
    ex_gap = -2
    # find the lengths of the sequences and uses them to intialize the score and traceback matrices
    matrix_length = len(best_db_seq) + 1
    matrix_height = len(best_q_seq) + 1
    score_matrix = np.empty((matrix_height, matrix_length))
    traceback_matrix = np.chararray((matrix_height, matrix_length))
    # sets the score in the left-top most cell of the matrix to 0 and the traceback to "s" for start
    score_matrix[0][0] = 0
    traceback_matrix[0][0] = "s"
    # uses the gap penalties to fill the top row of the matrices
    for x in range(0, matrix_length):
        score_matrix[0][x] = 0
        traceback_matrix[0][x] = "h"        
    # uses the gap penalties to fill the leftmost column of the matrices
    for y in range(0, matrix_height):
        score_matrix[y][0] = 0
        traceback_matrix[y][0] = "v"       
    # for each row in the matrix
    for i in range(1, matrix_height): 
        # for each cell in the row
        for j in range(1, matrix_length):
            vertical = 0
            horizontal = 0
            # if the cell above the current cell is being traced back vertically, then the vertical score for this cell is found by extending the gap
            if (traceback_matrix[i - 1][j] == b'v'):
                vertical = score_matrix[i - 1][j] + ex_gap
            # if the cell above the current cell is not being traced back vertically, then the vertical score for this cell is found by opening a gap
            else:
                vertical = score_matrix[i - 1][j] + open_gap
            if (0 > vertical):
                vertical = 0
            # if the cell left of the current cell is being traced back horizontally, then the horizontal score for this cell is found by 
            # extending the gap
            if (traceback_matrix[i][j - 1] == b'h'):
                horizontal = score_matrix[i][j - 1] + ex_gap
            # if the cell left of the current cell is not being traced back horizontally, then the horizontal score for this cell is found by 
            # opening a gap
            else:
                horizontal = score_matrix[i][j - 1] + open_gap
            if (0 > horizontal):
                horizontal = 0
            diagonal = 0
            # if the sequences match at this position then the diagonal score for this cell is found by adding the match score to the score 
            # of the diagonal cell
            if (best_db_seq[j - 1] == best_q_seq[i - 1]):
                diagonal = score_matrix[i - 1][j - 1] + match
            # if the sequences do not match at this position then the diagonal score for this cell is found by adding the mismatch score to the score 
            # of the diagonal cell
            elif (best_db_seq[j - 1] != best_q_seq[i - 1]):
                diagonal = score_matrix[i - 1][j - 1] + mismatch
            if (0 > diagonal):
                diagonal = 0
            # finds which score is the highest and then enters the highest score and the corresponding trace to the matrices
            max = vertical
            trace = "v"
            if (horizontal > max):
                max = horizontal
                trace = "h"
            if (diagonal >= max):
                max = diagonal
                trace = "d"
            score_matrix[i][j] = max
            traceback_matrix[i][j] = trace
    # finds the last position of the local alignment
    high_score = 0
    for o in range(0, matrix_height): 
        for p in range(0, matrix_length):
            if (score_matrix[o][p] >= high_score):
                high_score = score_matrix[o][p]
    r = o
    c = p
	# intializes empty strings to store the aligned sequences in
    best_db_seq_aligned = ""
    best_q_seq_aligned = ""
    while (score_matrix[r][c] > 0):
		# if the traceback for the current cell is vertical:
        if (traceback_matrix[r][c] == b'v'):
			# a gap is added to the beginning of the aligned best_db_seq
            best_db_seq_aligned = "-" + best_db_seq_aligned
			# the nucleotide from best_q_seq is added to the beginning of the aligned best_q_seq
            best_q_seq_aligned = best_q_seq[r - 1] + best_q_seq_aligned
			# the current row in the traceback is moved up 1
            r = r - 1
		# if the traceback for the current cell is horizontal:
        elif (traceback_matrix[r][c] == b'h'):
			# the nucleotide from best_db_seq is added to the beginnng of the aligned best_db_seq
            best_db_seq_aligned = best_db_seq[c - 1] + best_db_seq_aligned
			# a gap is added to the beginning of the aligned best_q_seq
            best_q_seq_aligned = "-" + best_q_seq_aligned
			# the current column in the traceback is moved left by 1
            c = c - 1
		# if the traceback for the current cell is diagonal:
        elif (traceback_matrix[r][c] == b'd'):
			# the nucleotide from best_db_seq is added to the beginnng of the aligned best_db_seq
            best_db_seq_aligned = best_db_seq[c - 1] + best_db_seq_aligned
			# the nucleotide from best_q_seq is added to the beginning of the aligned best_q_seq
            best_q_seq_aligned = best_q_seq[r - 1] + best_q_seq_aligned
			# the current row in the traceback is moved up by 1
            r = r - 1
			# the current column in the traceback is moved left by 1
            c = c - 1
    consensus_seq = ""
    con_seq_length = len(best_db_seq_aligned)
	# at each position in the alignment:
    for z in range(0, con_seq_length):
		# if the sequences match then that nucleotide is used in the consensus sequence
        if (best_db_seq_aligned[z] == best_q_seq_aligned[z]):
            consensus_seq = consensus_seq + best_db_seq_aligned[z]
		# if there is a gap in one of the sequences then the nucleotide from the other sequence is used in the consensus
        elif (best_db_seq_aligned[z] == '-'):
            consensus_seq = consensus_seq + best_q_seq_aligned[z]
        elif (best_q_seq_aligned[z] == '-'):
            consensus_seq = consensus_seq + best_db_seq_aligned[z]
        # if the nucleotides mismatch then the IUPAC symbol for that pairwise combination is found and used in the consensus sequence
        else:
            nuc_pos = ""
            if (best_q_seq_aligned[z] < best_db_seq_aligned[z]):
                nuc_pos = nuc_pos + best_q_seq_aligned[z] + best_db_seq_aligned[z]
            else:
                nuc_pos = nuc_pos + best_db_seq_aligned[z] + best_q_seq_aligned[z]
            if (nuc_pos == "AG"):
                consensus_seq = consensus_seq + "R"
            elif (nuc_pos == "CT"):
                consensus_seq = consensus_seq + "Y"
            elif (nuc_pos == "AC"):
                consensus_seq = consensus_seq + "M"
            elif (nuc_pos == "GT"):
                consensus_seq = consensus_seq + "K"
            elif (nuc_pos == "CG"):
                consensus_seq = consensus_seq + "S"
            elif (nuc_pos == "AT"):
                consensus_seq = consensus_seq + "W"

	# prints the alignment to the terminal 
    print("Consensus Seq:\t", consensus_seq)
    print("Database Seq: \t", best_db_seq_aligned)  
    print("Query Sequence:\t", best_q_seq_aligned) 
    print("\n")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(main.__doc__)
        exit()

    main(
    database_file = sys.argv[1],
    query_file = sys.argv[2],
    )