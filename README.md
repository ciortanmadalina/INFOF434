## INFOF434 exercises
http://homepages.ulb.ac.be/~dgonze/TEACHING/tp_needle.html

### align.py
Aim of the work

The aim of this work is to implement the algorithm of Needleman & Wunsh to align two protein sequences. The program should take as input two protein sequences (ex: sequences.txt), a substitution matrix (ex: pam250.tab), a gap penalty (assumed constant, e.g. -6), and no end_gap penalty (i.e. score=0). The program should return the score and the corresponding alignment.

The sequences and parameters proposed here are the following (but feel free to test other sequences):
       HEAGAWGHEE		(constant) gap penalty = -6
       PAWHEAE			end_gap penalty = 0

Implementation

The main steps of your program are:
Reading the input sequences file and storing the sequence in vectors.
Reading the scoring matrix file and storing it in a matrix.
Initializing the scores matrix.
Computing iteratively the scores.
Backtracking and deducing the alignment.
Display the score and the alignment on the screen.

Check and run

What score do you obtain for the sequences, substitution matrix, and gap penalty given above? Which alignment do you obtain?

You should obtain a score = 23 and the following alignment:
       HEAGAWGHEE-
       ---PAW-HEAE
Run the progam for other sequences, scoring matrices and gap penalties.

Statistics

How would you evaluate the significance of the alignment score obtained? Propose and implement a method to associate a p-value to the obtained alignment score.

To go further...

In this first step, we have considered only linear gap penalty. The case of affine gap penalty (in which different penalties are given to the opening and to the extension of gaps) is left to the more experienced programmers.