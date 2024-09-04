#include <string>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <bits/stdc++.h>

/*
usage:
	g++ UPGMA_Tree_Generation.cpp	
	./a.out
*/

using namespace std;
// creates matrices to store all of the alignments and their scores
string aligned_sequences[100][4];
int aligned_scores[100];
int align_count = 0;

// function to find the consensus sequence
string get_consensus(string seq1,string seq2)
{
	// creating vectors that contain all of the possible pairwise combinations resulting in the IUPAC symobls
	vector<std::string> B_combos{ "BC", "BG", "BT", "GY", "CK", "ST", "BY", "BK", "BS", "KS", "KY", "SY" };
	vector<std::string> D_combos{ "AD", "DG", "DT", "RT", "AK", "GW", "DR", "DK", "DW", "KR", "RW", "KW" };
	vector<std::string> H_combos{ "AH", "CH", "HT", "AY", "MT", "CW", "HY", "HM", "HW", "MY", "WY", "MW" };
	vector<std::string> V_combos{ "AV", "CV", "GV", "CR", "GM", "AS", "RV", "MV", "SV", "MR", "RS", "MS" };
	vector<std::string> N_combos{ "BD", "BH", "BV", "DH", "DV", "HV", "AB", "CD", "GH", "TV", "RY", "MK", "SW", "BR", "BM", "BW", "DY", "DM",
		"DS", "HR", "HK", "HS", "VY", "KV", "VW", "AN", "CN", "GN", "NT", "NR", "NY", "MN", "KN", "NS", "NW", "BN", "DN", "HN", "VN", "NN" };

	string consensus_seq = "";
	string seq1_aligned;
	string seq2_aligned;

	// finding the alignment in the matrix of aligned sequences
	for (int x = 0; x < align_count; x++)
	{
		if ((aligned_sequences[x][0] == seq1) && (aligned_sequences[x][1] == seq2))
		{
			seq1_aligned = aligned_sequences[x][2];
			seq2_aligned = aligned_sequences[x][3];
		}
	}
	int seq_length = seq1_aligned.length();
	// at each position in the alignment:
	for (int y = 0; y < seq_length; y++)
	{
		// if the sequences match then that nucleotide is used in the consensus sequence
		if (seq1_aligned[y] == seq2_aligned[y])
		{
			consensus_seq = consensus_seq + seq1_aligned[y];
		}
		// if there is a gap in one of the sequences then the nucleotide from the other sequence is used in the consensus
		else if (seq1_aligned[y] == '-')
		{
			consensus_seq = consensus_seq + seq2_aligned[y];
		}
		else if (seq2_aligned[y] == '-')
		{
			consensus_seq = consensus_seq + seq1_aligned[y];
		}
		// if the nucleotides mismatch then the IUPAC symbol for that pairwise combination is found and used in the consensus sequence
		else
		{
			string nuc_pos = "";
			if (seq2_aligned[y] < seq1_aligned[y])
			{
				nuc_pos = nuc_pos + seq2_aligned[y] + seq1_aligned[y];
			}
			else
			{
				nuc_pos = nuc_pos + seq1_aligned[y] + seq2_aligned[y];
			}
			if ((nuc_pos == "AG") || (nuc_pos == "AR") || (nuc_pos == "GR"))
			{
				consensus_seq = consensus_seq + "R";
			}
			else if ((nuc_pos == "CT") || (nuc_pos == "CY") || (nuc_pos == "TY"))
			{
				consensus_seq = consensus_seq + "Y";
			}
			else if ((nuc_pos == "AC") || (nuc_pos == "AM") || (nuc_pos == "CM"))
			{
				consensus_seq = consensus_seq + "M";
			}
			else if ((nuc_pos == "GT") || (nuc_pos == "GK") || (nuc_pos == "KT"))
			{
				consensus_seq = consensus_seq + "K";
			}
			else if ((nuc_pos == "CG") || (nuc_pos == "CS") || (nuc_pos == "GS"))
			{
				consensus_seq = consensus_seq + "S";
			}
			else if ((nuc_pos == "AT") || (nuc_pos == "AW") || (nuc_pos == "TW"))
			{
				consensus_seq = consensus_seq + "W";
			}
			else if (std::find(std::begin(B_combos), std::end(B_combos), nuc_pos) != std::end(B_combos))
			{
				consensus_seq = consensus_seq + "B";
			}
			else if (std::find(std::begin(D_combos), std::end(D_combos), nuc_pos) != std::end(D_combos))
			{
				consensus_seq = consensus_seq + "D";
			}
			else if (std::find(std::begin(H_combos), std::end(H_combos), nuc_pos) != std::end(H_combos))
			{
				consensus_seq = consensus_seq + "H";
			}
			else if (std::find(std::begin(V_combos), std::end(V_combos), nuc_pos) != std::end(V_combos))
			{
				consensus_seq = consensus_seq + "V";
			}
			else if (std::find(std::begin(N_combos), std::end(N_combos), nuc_pos) != std::end(N_combos))
			{
				consensus_seq = consensus_seq + "N";
			}
		}
	}
	return consensus_seq;
}

int nucleotide_align(string seq1, string seq2, int match, int mismatch, int open_gap, int ex_gap)
{
	int align_score;
	string already_aligned;
	int align_found_pos = 0;
	// finds whether the alignment has already been calculated
	for (int a = 0; a < align_count; a++)
	{
		if ((aligned_sequences[a][0] == seq1) && (aligned_sequences[a][1] == seq2))
		{
			already_aligned = "yes";
			align_found_pos = a;
		}
	}
	// if the alignment has already been performed the score is found in the matrix and returned to the main function
	if (already_aligned == "yes")
	{
		align_score = aligned_scores[align_found_pos];
		return align_score;
	}
	// the alignment is performed if it hasn't been already
	else
	{
		// vector containing all of the possible mismatched nucleotide pairwise combinations
		vector<std::string> mismatched_nucs{ "AC", "AG", "AT", "CG", "CT", "GT", "CR", "TR", "AY", "GY", "MT", "GM",
		"AK", "CK", "AS", "ST", "GW", "CW", "AB", "CD", "GH", "TV", "RY", "KM", "SW" };

		// find the lengths of the sequences and uses them to intialize the score and traceback matrices
		int matrix_length = seq1.length() + 1;
		int matrix_height = seq2.length() + 1;
		int score_matrix[matrix_height][matrix_length];
		string traceback_matrix[matrix_height][matrix_length];
		// sets the score in the left-top most cell of the matrix to 0 and the traceback to "s" for start
		score_matrix[0][0] = 0;
		traceback_matrix[0][0] = "s";
		// uses the gap penalties to fill the top row of the matrices
		for (int x = 1; x < matrix_length; x++)
		{
			score_matrix[0][x] = open_gap + ex_gap * x;
			traceback_matrix[0][x] = "h";
		}
		// uses the gap penalties to fill the leftmost column of the matrices
		for (int y = 1; y < matrix_height; y++)
		{
			score_matrix[y][0] = open_gap + ex_gap * y;
			traceback_matrix[y][0] = "v";
		}
		// for each row in the matrix
		for (int i = 1; i < matrix_height; i++)
		{
			// for each cell in the row
			for (int j = 1; j < matrix_length; j++)
			{
				int vertical;
				int horizontal;
				// if the cell above the current cell is being traced back vertically, then the vertical score for this cell is found by extending the gap
				if (traceback_matrix[i - 1][j] == "v")
				{
					vertical = score_matrix[i - 1][j] + ex_gap;
				}
				// if the cell above the current cell is not being traced back vertically, then the vertical score for this cell is found by opening a gap
				else
				{
					vertical = score_matrix[i - 1][j] + open_gap;
				}
				// if the cell left of the current cell is being traced back horizontally, then the horizontal score for this cell is found by 
				// extending the gap
				if (traceback_matrix[i][j - 1] == "h")
				{
					horizontal = score_matrix[i][j - 1] + ex_gap;
				}
				// if the cell left of the current cell is being traced back horizontally, then the horizontal score for this cell is found by 
				// extending the gap
				else
				{
					horizontal = score_matrix[i][j - 1] + open_gap;
				}
				int diagonal;
				string nuc_pos = "";
				// determining whether the nucleotides at the positions are matched or mismatched to then implement the match or mismatch score 
				// to calculate the diagonal score
				if (seq2[i - 1] < seq1[j - 1])
				{
					nuc_pos = nuc_pos + seq2[i - 1] + seq1[j - 1];
				}
				else
				{
					nuc_pos = nuc_pos + seq1[j - 1] + seq2[i - 1];
				}
				if (std::find(std::begin(mismatched_nucs), std::end(mismatched_nucs), nuc_pos) != std::end(mismatched_nucs))
				{
					diagonal = score_matrix[i - 1][j - 1] + mismatch;
				}
				else
				{
					diagonal = score_matrix[i - 1][j - 1] + match;
				}
				// finds which score is the highest and then enters the highest score and the corresponding trace to the matrices
				int max = vertical;
				string trace = "v";
				if (horizontal > max)
				{
					max = horizontal;
					trace = "h";
				}
				if (diagonal >= max)
				{
					max = diagonal;
					trace = "d";
				}
				score_matrix[i][j] = max;
				traceback_matrix[i][j] = trace;
			}
		}
		// the score for the alignment is in the bottom right-most cell of the score matrix
		align_score = score_matrix[matrix_height - 1][matrix_length - 1];
		// assign the current row and column positions in the traceback to r and c
		int r = matrix_height - 1;
		int c = matrix_length - 1;
		// intializes empty strings to store the aligned sequences in
		string seq1_aligned = "";
		string seq2_aligned = "";
		do {
			// if the traceback for the current cell is vertical:
			if (traceback_matrix[r][c] == "v")
			{
				// a gap is added to the beginning of the aligned seq1
				seq1_aligned = "-" + seq1_aligned;
				// the nucleotide from seq2 is added to the beginning of the aligned seq2
				seq2_aligned = seq2[r - 1] + seq2_aligned;
				// the current row in the traceback is moved up 1
				r = r - 1;
			}
			// if the traceback for the current cell is horizontal:
			else if (traceback_matrix[r][c] == "h")
			{
				// the nucleotide from seq1 is added to the beginnng of the aligned seq1
				seq1_aligned = seq1[c - 1] + seq1_aligned;
				// a gap is added to the beginning of the aligned seq2
				seq2_aligned = "-" + seq2_aligned;
				// the current column in the traceback is moved left by 1
				c = c - 1;
			}
			// if the traceback for the current cell is diagonal:
			else if (traceback_matrix[r][c] == "d")
			{
				// the nucleotide from seq1 is added to the beginnng of the aligned seq1
				seq1_aligned = seq1[c - 1] + seq1_aligned;
				// the nucleotide from seq2 is added to the beginning of the aligned seq2
				seq2_aligned = seq2[r - 1] + seq2_aligned;
				// the current row in the traceback is moved up 1
				r = r - 1;
				// the current column in the traceback is moved left by 1
				c = c - 1;
			}
			// repeat until the top-left corner of the matrix is reached
		} while ((r > 0) && (c > 0));

		// the unaligned, aligned sequences, and the alignment score are added to the alignment matrices
		// the align score is returned to the main function
		string alignment = seq1_aligned + "\n" + seq2_aligned + "\n";
		aligned_sequences[align_count][0] = seq1;
		aligned_sequences[align_count][1] = seq2;
		aligned_sequences[align_count][2] = seq1_aligned;
		aligned_sequences[align_count][3] = seq2_aligned;
		aligned_scores[align_count] = align_score;
		align_count++;
		return align_score;
	}
}

int main()
{
	// prompts user to enter the filename 
	string filename;
	cout << "Enter fasta file:\n";
	cin >> filename;
	ifstream fin;
	fin.open(filename);

	string sequences[20];
	int seq_count = 0;

	string seq_names[20];

	// read in the sequences and their names from the file into matrices
	while (getline(fin, filename))
	{
		if (filename.find(">") != std::string::npos) {
			seq_count++;
			int split = filename.find_first_of("\r");
			seq_names[seq_count - 1] = filename.substr(0, split);
		}
		else {
			int split = filename.find_first_of("\r");
			sequences[seq_count-1] += filename.substr(0, split);
		}
	}

	fin.close();

	int match;
	int mismatch;
	int open_gap;
	int ex_gap;

	// prompt user to enter the match, mismatch, opening gap, and extending gap scores
	cout << "Enter the match score:\n";
	cin >> match;

	cout << "Enter the mismatch penalty:\n";
	cin >> mismatch;

	cout << "Enter the opening gap penalty:\n";
	cin >> open_gap;

	cout << "Enter the extending gap penalty:\n";
	cin >> ex_gap;

	// creating a matrix that contains all the sequences in the pool
	int total_seqs = seq_count;
	string pool[total_seqs];
	for (int a = 0; a < total_seqs; a++)
	{
		pool[a] = sequences[a];
	}

	do {
		int high_score = -5000;
		string best_seq1;
		string best_seq2;
		// for each possible combination of sequences in the pool:
		for (int i = 0; i < seq_count - 1; i++)
		{
			for (int j = i + 1; j < seq_count; j++)
			{
				int score;
				string seq1 = pool[i];
				string seq2 = pool[j];
				// the alignment is found
				score = nucleotide_align(seq1, seq2, match, mismatch, open_gap, ex_gap);
				// if the score for the alignment is greater than the current high score then the alignment becomes the newest highest scoring
				// alignment
				if (score >= high_score)
				{
					high_score = score;
					best_seq1 = seq1;
					best_seq2 = seq2;
				}
			}
		}
		int m = 1;
		// removing the sequences that had the highest alignment score from the pool
		for (int c = 0; c < seq_count; c++)
		{
			if ((pool[c] == best_seq1) || (pool[c] == best_seq2))
			{
				for (int k = c; k < seq_count - m; k++)
				{
					pool[k] = pool[k + 1];
				}
				pool[seq_count - m] = "";
				m++;
				if ((pool[c] == best_seq1) || (pool[c] == best_seq2))
				{
					for (int k = c; k < seq_count - m; k++)
					{
						pool[k] = pool[k + 1];
					}
					pool[seq_count - m] = "";
					m++;
				}
			}

		}
		seq_count--;
		// adding the consensus sequence generated from the best alignment to the pool
		pool[seq_count - 1] = get_consensus(best_seq1, best_seq2);
		// repeat until left with only the one consensus sequence in the pool
	} while (seq_count > 1);

	string MSA[total_seqs];

	// stores the sequences aligned to the consensus sequence in a matrix
	for (int z = 0; z < total_seqs; z++)
	{
		int a_score = nucleotide_align(pool[0], sequences[z], match, mismatch, open_gap, ex_gap);
		MSA[z] = aligned_sequences[align_count - 1][3];
	}
	// creates matrices to store tree clusters and their distances
	string tree_calcs[total_seqs - 1];
	float tree_distances[total_seqs - 1];
	int count = 0;

	// creates a distance matrix
	float distances_new[total_seqs][total_seqs];
	string clusters_new[total_seqs];
	// fills the distance matrix using Kimura
	for (int m = 0; m < total_seqs; m++)
	{
		clusters_new[m] = seq_names[m];
		for (int n = 0; n < total_seqs; n++)
		{
			int transversions = 0;
			int transitions = 0;
			if (n >= m)
			{
				distances_new[m][n] = 0;
			}
			else
			{
				string seq1 = MSA[m];
				string seq2 = MSA[n];
				int seq_length = seq1.length();
				for (int o = 0; o < seq_length; o++)
				{
					string nuc_pos = "";
					if (seq2[o] < seq1[o])
					{
						nuc_pos = nuc_pos + seq2[o] + seq1[o];
					}
					else
					{
						nuc_pos = nuc_pos + seq1[o] + seq2[o];
					}
					if ((nuc_pos == "AG") || (nuc_pos == "CT"))
					{
						transversions++;
					}
					else if ((nuc_pos == "AC") || (nuc_pos == "AT") || (nuc_pos == "CG") || (nuc_pos == "GT"))
					{
						transitions++;
					}
				}
				float ts = transitions * 1.0 / seq_length;
				float tv = transversions * 1.0 / seq_length;
				distances_new[m][n] = 0.5 * log(1 / (1 - (2 * ts) - tv)) + (0.25 * log(1 / (1 - (2 * tv))));
			}
		}
	}

	float distances_old[total_seqs][total_seqs];
	string clusters_old[total_seqs];

	// finds smallest distance, creates the new cluster, calculates the new distance matrix, repeats until all clusters have been joined together
	for (int s = total_seqs; s > 2; s--)
	{
		int pos1;
		int pos2;
		float min = 50000;
		for (int t = 0; t < s; t++)
		{
			clusters_old[t] = clusters_new[t];
			for (int u = 0; u < s; u++)
			{
				distances_old[t][u] = distances_new[t][u];
				if ((distances_old[t][u] < min) && (u != t) && (distances_old[t][u]))
				{
					min = distances_old[t][u];
					pos1 = u;
					pos2 = t;
				}
			}
		}
		int cluster_pos;
		int emp_pos;
		if (pos1 < pos2)
		{
			cluster_pos = pos1;
			emp_pos = pos2;
		}
		else if (pos2 < pos1)
		{
			cluster_pos = pos2;
			emp_pos = pos1;
		}
		string clus2 = "no";
		string form_cluster = "";
		for (int f = 0; f < s - 1; f++)
		{
			if (f == cluster_pos)
			{
				clusters_new[f] = clusters_old[cluster_pos] + "\t" + clusters_old[emp_pos];
				form_cluster = clusters_old[cluster_pos] + "\t" + clusters_old[emp_pos];
				tree_calcs[count] = clusters_new[f];
				tree_distances[count] = min;
				count++;
			}
			else if (f == emp_pos)
			{
				clus2 = "yes";
				if (f != s - 1)
				{
					clusters_new[f] = clusters_old[f + 1];
				}
			}
			else if (clus2 == "no")
			{
				clusters_new[f] = clusters_old[f];
			}
			else if (clus2 == "yes")
			{
				clusters_new[f] = clusters_old[f + 1];
			}

		}
		for (int d = 0; d < s - 1; d++)
		{
			for (int e = 0; e < s - 1; e++)
			{
				if (e >= d)
				{
					distances_new[d][e] = 0;
				}
				else if ((form_cluster == tree_calcs[count - 1]) || (form_cluster == tree_calcs[count - 1]))
				{
					int new_calc;
					float new_dist = 0;
					for (int b = 0; b < s; b++)
					{
						if (clusters_old[b] == clusters_new[d])
						{
							new_calc = b;
						}
					}
					if (new_calc > cluster_pos)
					{
						new_dist = new_dist + distances_old[new_calc][cluster_pos];
					}
					else if (new_calc < cluster_pos)
					{
						new_dist = new_dist + distances_old[cluster_pos][new_calc];
					}
					if (new_calc > emp_pos)
					{
						new_dist = new_dist + distances_old[new_calc][emp_pos];
					}
					else if (new_calc < emp_pos)
					{
						new_dist = new_dist + distances_old[emp_pos][new_calc];
					}
					distances_new[d][e] = new_dist / 2;
				}
				else
				{
					int old_d;
					int old_e;
					for (int j = 0; j < s; j++)
					{
						if (clusters_old[j] == clusters_new[d])
						{
							old_d = j;
						}
						if (clusters_old[j] == clusters_new[e])
						{
							old_e = j;
						}
					}
					distances_new[d][e] = distances_old[old_d][old_e];
				}
			}
		}
	}
	tree_calcs[total_seqs - 2] = clusters_new[0] + "\t" + clusters_new[1];
	tree_distances[total_seqs - 2] = distances_new[1][0];

	// initializing a matrix to create the tree in
	int rows = total_seqs + ((total_seqs - 1) * 3);
	int columns = total_seqs * 3 + 1;
	string tree[rows][columns];
	for (int r1 = 0; r1 < rows; r1++)
	{
		for (int c1 = 0; c1 < columns; c1++)
		{
			tree[r1][c1] = "     ";
		}
	}
	// adding the first two branches to the tree
	int name_length = tree_calcs[0].length();
	int split = tree_calcs[0].find_first_of("\t");
	string branch1 = tree_calcs[0].substr(0, split);
	string branch2 = tree_calcs[0].substr(split + 1, name_length - 1);
	tree[0][columns - 1] = branch1;
	tree[4][columns - 1] = branch2;
	for (int a = columns - 2; a >= columns - 4; a--)
	{
		tree[0][a] = "-----";
		tree[4][a] = "-----";
	}
	string dis = to_string(tree_distances[0]);
	int round_dis = dis.find_first_of(".");
	dis = dis.substr(0, round_dis + 3);
	tree[1][columns - 4] = "|";
	tree[2][columns - 4] = "|" + dis;
	tree[3][columns - 4] = "|";
	int left_most = columns - 4;
	int down_most = 4;
	int mid_point[total_seqs - 1];
	mid_point[0] = 2;
	// adding the rest of the branches on to the first two
	for (int branch = 1; branch < total_seqs - 1; branch++)
	{
		string cluster_name = tree_calcs[branch];
		string last_cluster = tree_calcs[branch - 1];
		name_length = cluster_name.length();
		// if the next cluster being added to the tree is contains the last cluster
		if (cluster_name.find(last_cluster) != std::string::npos)
		{
			int last_midpoint = mid_point[branch - 1];
			int clus_loc = cluster_name.find(last_cluster);
			if (clus_loc == 0)
			{
				int start_new = last_cluster.length() + 1;
				branch1 = cluster_name.substr(start_new, name_length - 1);
			}
			else
			{
				branch1 = cluster_name.substr(0, clus_loc - 1);
			}
			int old_branch;
			string connect_old = "no";
			for (int check = 0; check < branch; check++)
			{
				if (branch1 == tree_calcs[check])
				{
					old_branch = check;
					connect_old = "yes";
				}
			}
			// adding a branch to the last cluster on the tree
			if (connect_old == "no")
			{
				tree[down_most + 4][columns - 1] = branch1;
				for (int b = columns - 2; b >= left_most - 3; b--)
				{
					tree[down_most + 4][b] = "-----";
					if (b < left_most)
					{
						tree[last_midpoint][b] = "-----";
					}
				}
				for (int gh = last_midpoint + 1; gh <= down_most + 3; gh++)
				{
					tree[gh][left_most - 3] = "|    ";
				}
				mid_point[branch] = (mid_point[branch - 1] + down_most + 4) / 2;
				dis = to_string(tree_distances[branch]);
				round_dis = dis.find_last_of(".");
				dis = dis.substr(0, round_dis + 3);
				int new_midpoint = mid_point[branch];
				tree[new_midpoint][left_most - 3] = "|" + dis;
				left_most = left_most - 3;
				down_most = down_most + 4;
			}
			// connecting two preexisting clusters on the tree to each other
			else if (connect_old == "yes")
			{
				int old_midpoint = mid_point[old_branch];
				for (int c = left_most - 1; c >= left_most - 3; c--)
				{
					tree[last_midpoint][c] = "-----";
				}
				int old_line = left_most - 3;
				string check_line = tree[old_midpoint][old_line];
				while (check_line.find("|") == std::string::npos)
				{
					tree[old_midpoint][old_line] = "-----";
				}
				for (int de = old_midpoint + 1; de <= last_midpoint - 1; de++)
				{
					tree[de][left_most - 3] = "|    ";
				}
				int new_midpoint = (old_midpoint + last_midpoint) / 2;
				mid_point[branch] = (old_midpoint + last_midpoint) / 2;
				dis = to_string(tree_distances[branch]);
				round_dis = dis.find_last_of(".");
				dis = dis.substr(0, round_dis + 3);
				tree[new_midpoint][left_most - 3] = "|" + dis;
				left_most = left_most - 3;
			}
		}
		else
		{
			// check whether the cluster contains any of the clusters present on the tree so far
			int old_branch;
			string connect_old = "no";
			for (int past = 0; past < branch; past++)
			{
				if (cluster_name.find(tree_calcs[past]) != std::string::npos)
				{
					old_branch = past;
					connect_old = "yes";
				}
			}
			// adding a new cluster not yet connected to the other on the tree
			if (connect_old == "no")
			{
				int split = cluster_name.find_first_of("\t");
				string branch1 = cluster_name.substr(0, split);
				string branch2 = cluster_name.substr(split + 1, name_length - 1);
				tree[down_most + 4][columns - 1] = branch1;
				tree[down_most + 8][columns - 1] = branch2;
				for (int a = columns - 2; a >= columns - 4; a--)
				{
					tree[down_most + 4][a] = "-----";
					tree[down_most + 8][a] = "-----";
				}
				dis = to_string(tree_distances[branch]);
				round_dis = dis.find_first_of(".");
				dis = dis.substr(0, round_dis + 3);
				tree[down_most + 5][columns - 4] = "|";
				tree[down_most + 6][columns - 4] = "|" + dis;
				tree[down_most + 7][columns - 4] = "|";
				mid_point[branch] = down_most + 6;
				down_most = down_most + 8;
			}
			// adds new branch to cluster that is not the previous cluster added by shifting all the other clusters down
			else if (connect_old == "yes")
			{
				int found_last;
				for (int fl = 0; fl <= down_most; fl++)
				{
					if (cluster_name.find(tree[fl][columns - 1]) != std::string::npos)
					{
						found_last = fl;
					}
				}
				for (int r3 = down_most; r3 > found_last; r3--)
				{
					for (int c3 = left_most; c3 < columns; c3++)
					{
						tree[r3 + 4][c3] = tree[r3][c3];
					}
				}
				int clus_loc = cluster_name.find(tree_calcs[old_branch]);
				if (clus_loc == 0)
				{
					int start_new = last_cluster.length() + 1;
					branch1 = cluster_name.substr(start_new, name_length - 1);
				}
				else
				{
					branch1 = cluster_name.substr(0, clus_loc - 1);
				}
				tree[found_last + 4][columns - 1] = branch1;
				int old_midpoint = mid_point[old_branch];
				for (int b = columns - 2; b >= left_most - 3; b--)
				{
					tree[found_last + 4][b] = "-----";
					if (b < left_most)
					{
						tree[old_midpoint][b] = "-----";
					}
				}
				for (int gh = old_midpoint + 1; gh <= found_last + 3; gh++)
				{
					tree[gh][left_most - 3] = "|    ";
				}
				int new_midpoint = (old_midpoint + found_last + 4) / 2;
				mid_point[branch] = new_midpoint;
				dis = to_string(tree_distances[branch]);
				round_dis = dis.find_last_of(".");
				dis = dis.substr(0, round_dis + 3);
				tree[new_midpoint][left_most - 3] = "|" + dis;
				left_most = left_most - 3;
				down_most = down_most + 4;
			}
		}

	}
	// prints tree to terminal and text file named "UPGMA_tree.txt"
	ofstream fout;
	string outfilename = "UPGMA_tree.txt";
	fout.open(outfilename);
	for (int r2 = 0; r2 < rows; r2++)
	{
		for (int c2 = 0; c2 < columns; c2++)
		{
			fout << tree[r2][c2];
			cout << tree[r2][c2];
		}
		fout << "\n";
		cout << "\n";
	}
	fout.close();


}