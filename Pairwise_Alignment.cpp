#include <string>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

void nucleotide_align(string seq1, string seq2)
{
	int match;
	int mismatch;
	int gap;

	cout << "Enter the match score:\n";
	cin >> match;

	cout << "Enter the mismatch penalty:\n";
	cin >> mismatch;

	cout << "Enter the gap penalty:\n";
	cin >> gap;

	cout << "seq1: " << seq1 << "\n";
	cout << "seq2: " << seq2 << "\n";
	cout << "match: " << match << " mismatch: " << mismatch << " gap: " << gap << "\n";

	int matrix_length = seq1.length() + 1;
	// cout << "matrix length: " << matrix_length << "\n";
	int matrix_height = seq2.length() + 1;
	// cout << "matrix height: " << matrix_height << "\n";
	int score_matrix[matrix_height][matrix_length];
	string traceback_matrix[matrix_height][matrix_length];
	score_matrix[0][0] = 0;
	traceback_matrix[0][0] = "s";
	for (int x = 1; x < matrix_length; x++)
	{
		score_matrix[0][x] = gap * x;
		traceback_matrix[0][x] = "h";
	}
	for (int y = 0; y < matrix_height; y++)
	{
		score_matrix[y][0] = gap * y;
		traceback_matrix[y][0] = "v";
	}
	for (int i = 1; i < matrix_height; i++)
	{
		for (int j = 1; j < matrix_length; j++)
		{
			int vertical = score_matrix[i - 1][j] + gap;
			int horizontal = score_matrix[i][j - 1] + gap;
			int diagonal;
			if (seq1[j - 1] == seq2[i - 1])
			{
				diagonal = score_matrix[i - 1][j - 1] + match;
			}
			else if (seq1[j - 1] != seq2[i - 1])
			{
				diagonal = score_matrix[i - 1][j - 1] + mismatch;
			}
			int max = vertical;
			string trace = "v";
			if (horizontal > max)
			{
				max = horizontal;
				trace = "h";
			}
			else if (diagonal >= max)
			{
				max = diagonal;
				trace = "d";
			}
			score_matrix[i][j] = max;
			traceback_matrix[i][j] = trace;
		}
	}

	//for (int e = 0; e < matrix_height; e++)
	//{
	//	for (int t = 0; t < matrix_length; t++)
	//	{
	//		cout << traceback_matrix[e][t] << "\t";
	//	}
	//	cout << "\n";
	//}

	int align_score = score_matrix[matrix_height - 1][matrix_length - 1];
	int r = matrix_height - 1;
	int c = matrix_height - 1;
	string seq_alignment = traceback(seq1, seq2, matrix_height, matrix_length, score_matrix, traceback_matrix, r, c);

}

void protein_align(string seq1, string seq2)
{

	int gap;
	cout << "Enter the gap penalty:\n";
	cin >> gap;

	string pam_file = "pam2.txt";
	string pam100[21][21];
	ifstream pin;
	pin.open(pam_file);
	for (int i = 0; i < 21; i++)
	{
		for (int j = 0; j < 21; j++)
		{
			if ((i == 0) && (j == 0)) { pam100[i][j] = "X"; }
			else {
				pin >> pam100[i][j];
			}

		}
	}
	pin.close();

	int matrix_length = seq1.length() + 1;
	cout << "matrix length: " << matrix_length << "\n";
	int matrix_height = seq2.length() + 1;
	cout << "matrix height: " << matrix_height << "\n";
	int score_matrix[matrix_height][matrix_length];
	string traceback_matrix[matrix_height][matrix_length];
	score_matrix[0][0] = 0;
	traceback_matrix[0][0] = "s";
	for (int x = 1; x < matrix_length; x++)
	{
		score_matrix[0][x] = gap * x;
		traceback_matrix[0][x] = "h";
	}
	for (int y = 0; y < matrix_height; y++)
	{
		score_matrix[y][0] = gap * y;
		traceback_matrix[y][0] = "v";
	}
	for (int i = 1; i < matrix_height; i++)
	{
		for (int j = 1; j < matrix_length; j++)
		{
			int vertical = score_matrix[i - 1][j] + gap;
			int horizontal = score_matrix[i][j - 1] + gap;
			// char pos1 = seq1[j - 1];
			string pos1(1, seq1[j - 1]);
			// char pos2 = seq2[i - 1];
			string pos2(1, seq2[i - 1]);
			int matrix_c;
			int matrix_r;
			// cout << pos1 << "\t" << pos2 << "\n";
			for (int l = 0; l < 21; l++)
			{
				if (pam100[0][l] == pos1)
				{
					matrix_c = l;
				}
			}
			for (int p = 0; p < 21; p++)
			{
				if (pam100[p][0] == pos2)
				{
					matrix_r = p;
				}
			}
			int diagonal = score_matrix[i - 1][j - 1] + stoi(pam100[matrix_r][matrix_c]);

			int max = vertical;
			string trace = "v";
			if (horizontal > max) 
			{
				max = horizontal;
				trace = "h";
			}
			else if (diagonal >= max)
			{
				max = diagonal;
				trace = "d";
			}
			score_matrix[i][j] = max;
			traceback_matrix[i][j] = trace;
		}
	}
	int end_value = score_matrix[matrix_height - 1][matrix_length - 1];
	int r = matrix_height - 1;
	int c = matrix_height - 1;
	string protein_alignment = traceback(seq1, seq2, matrix_height, matrix_length, score_matrix, traceback_matrix, r, c);
	
}


string traceback(string seq1, string seq2, int matrix_height, int matrix_length, int score_matrix[matrix_height][matrix_length], string traceback_matrix[matrix_height][matrix_length], int r, int c)
{
	string seq1_aligned = "";
	string seq2_aligned = "";
	do {

		if (traceback_matrix[r][c] == "v")
		{
			seq1_aligned = "-" + seq1_aligned;
			seq2_aligned = seq2[r - 1] + seq2_aligned;
			r = r - 1;
		}
		else if (traceback_matrix[r][c] == "h")
		{
			seq1_aligned = seq1[c - 1] + seq1_aligned;
			seq2_aligned = "-" + seq2_aligned;
			c = c - 1;
		}
		else if (traceback_matrix[r][c] == "d")
		{
			seq1_aligned = seq1[c - 1] + seq1_aligned;
			seq2_aligned = seq2[r - 1] + seq2_aligned;
			r = r - 1;
			c = c - 1;
		}
	} while ((r > 0) && (c > 0));

	string alignment =  seq1_aligned + "\n" + seq2_aligned + "\n";
	cout << alignment;
	return alignment; 
}

int main()
{

	string filename;
	cout << "Enter fasta file:\n";
	cin >> filename;
	ifstream fin;
	fin.open(filename);

	string seq1;
	string seq2;
	int seq_count = 0;

	while (getline(fin, filename))
	{
		if (filename.find(">") != std::string::npos) {
			seq_count++;
		}
		else {
			if (seq_count == 1) {
				int split = filename.find_first_of("\r");
				seq1 = seq1 + filename.substr(0, split);
			}
			else if (seq_count == 2) {
				int split = filename.find_first_of("\r");
				seq2 = seq2 + filename.substr(0, split);
			}
		}
	}
	fin.close();

	string aligntype;
	cout << "Enter whether nucleotide or protein:\n";
	cin >> aligntype;
	if (aligntype == "protein")
	{
		protein_align(seq1, seq2);
	}
	else if (aligntype == "nucleotide")
	{
		nucleotide_align(seq1, seq2);
	}

}