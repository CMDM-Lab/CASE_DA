# CASE_DA
Codes and data for the journal article "Quantum-Inspired Approach to Natural Product Structure Elucidation"
CASE: Computer-Assisted Structure Elucidation, DA: Digital Annealer

# example.py
The example code to run DA to find solutions using Ophiopogon japonicus scaffold 1 and configuration 1 shown in the paper

# data format of .cIdx file
Only the data used in in this study is explained.

s0000009929(scaffold ID)	8	8
4362.736772999938	67.04219899999998
4(number of configurations)
4(number of sidechains for this configuration)	0(sidechain at position 0)	1(position 1)	20(position 20)	2(position 2)  5	ALLDB0000009095	ALLDB0000055239	ALLDB0000142322	ALLDB0000160953	ALLDB0000186912	
5	0	1	20	2	5  1	ALLDB0000054755	
5	0	1	20	2	16  1	ALLDB0000159623	
3	0	1	20  1	ALLDB0000172974	
21(positions in this scaffold)
0(position 0)	2(possible sidechains)
[*]C([H])([H])[H]	(SMILES of this sidechain) 6 1(valence)	1.3362791E-5(probability)	15.023475000000001(molecular weight) 0.0
[*]C([H])=O	2 1	4.4542635E-6	29.002740000000003 0.0
1 2
...
