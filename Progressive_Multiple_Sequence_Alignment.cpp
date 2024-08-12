#include <string>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <bits/stdc++.h>

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

	// print the alignment to the terminal and to a text file named "MSA_output.txt"
	ofstream fout;
	string outfilename = "MSA_output.txt";
	fout.open(outfilename);

	for (int z = 0; z < total_seqs; z++)
	{
		int a_score = nucleotide_align(pool[0], sequences[z], match, mismatch, open_gap, ex_gap);
		if (z == 0)
		{ 
			fout << "Consensus Sequence\t" <<  aligned_sequences[align_count - 1][2]  << "\n";
			cout << "\nConsensus Sequence\t" << aligned_sequences[align_count - 1][2] << "\n";
		}
		fout << seq_names[z] << "\t" << aligned_sequences[align_count - 1][3]  << "\n";
		cout << seq_names[z] << "\t" << aligned_sequences[align_count - 1][3] << "\n";
	}
	fout.close();
}