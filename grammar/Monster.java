package grammar;

import biology.Base;
import math.Probability;

import java.io.IOException;
import java.nio.file.FileSystems;
import java.nio.file.Files;
import java.util.List;

public class Monster{
	public static SCFG Gs;

	public static Nonterminal S;
	public static Nonterminal T;
	public static Nonterminal C;
	public static Nonterminal A;
	public static Nonterminal P;
	public static Nonterminal L;
	public static Nonterminal F;
	public static Nonterminal H;
	public static Nonterminal G;
	public static Nonterminal B;
	public static Nonterminal M;
	public static Nonterminal O;
	public static Nonterminal N;
	public static Nonterminal U;
	public static Nonterminal Z;

	public static Terminal dot;
	public static Terminal open;
	public static Terminal close;

	public static int min_HEL;
	public static int min_HL;
	public static int min_PS;

	public static double[] em1;   //emission probabilities for bases
	public static double[][] em2; //emission probabilities for base pairs

	public static void init(String conffile){
		double values[] = {};
		try{
			List<String> filelist = Files.readAllLines(FileSystems.getDefault().getPath(conffile));
			values = new double[51];
			int i = 0;
			for(String s : filelist){
				if(s.length() == 0 || s.startsWith("#")) continue;
				values[i] = Double.parseDouble(s);
				i++;
			}
		}catch(IOException e){
			System.out.println("Error reading file");
			System.exit(0);
		}

		int iter = 0;
		min_HL = (int) values[iter++]; // common choices: 1 or 3
		min_HEL = (int) values[iter++]; // common choices: 1 or 2
		min_PS = 2 * min_HEL + min_HL;

		em1 = new double[Base.count];
		em1[Base.adenine.getId()] = Probability.fromNumber(values[iter++]);
		em1[Base.cytosine.getId()] = Probability.fromNumber(values[iter++]);
		em1[Base.guanine.getId()] = Probability.fromNumber(values[iter++]);
		em1[Base.uracil.getId()] = Probability.fromNumber(values[iter++]);

		em2 = new double[Base.count][Base.count];
		em2[Base.adenine.getId()][Base.adenine.getId()] = Probability.fromNumber(values[iter++]);
		em2[Base.adenine.getId()][Base.cytosine.getId()] = Probability.fromNumber(values[iter++]);
		em2[Base.adenine.getId()][Base.guanine.getId()] = Probability.fromNumber(values[iter++]);
		em2[Base.adenine.getId()][Base.uracil.getId()] = Probability.fromNumber(values[iter++]);
		em2[Base.cytosine.getId()][Base.adenine.getId()] = Probability.fromNumber(values[iter++]);
		em2[Base.cytosine.getId()][Base.cytosine.getId()] = Probability.fromNumber(values[iter++]);
		em2[Base.cytosine.getId()][Base.guanine.getId()] = Probability.fromNumber(values[iter++]);
		em2[Base.cytosine.getId()][Base.uracil.getId()] = Probability.fromNumber(values[iter++]);
		em2[Base.guanine.getId()][Base.adenine.getId()] = Probability.fromNumber(values[iter++]);
		em2[Base.guanine.getId()][Base.cytosine.getId()] = Probability.fromNumber(values[iter++]);
		em2[Base.guanine.getId()][Base.guanine.getId()] = Probability.fromNumber(values[iter++]);
		em2[Base.guanine.getId()][Base.uracil.getId()] = Probability.fromNumber(values[iter++]);
		em2[Base.uracil.getId()][Base.adenine.getId()] = Probability.fromNumber(values[iter++]);
		em2[Base.uracil.getId()][Base.cytosine.getId()] = Probability.fromNumber(values[iter++]);
		em2[Base.uracil.getId()][Base.guanine.getId()] = Probability.fromNumber(values[iter++]);
		em2[Base.uracil.getId()][Base.uracil.getId()] = Probability.fromNumber(values[iter++]);


		S = new Nonterminal('S');
		T = new Nonterminal('T');
		C = new Nonterminal('C');
		A = new Nonterminal('A');
		P = new Nonterminal('P');
		L = new Nonterminal('L');
		F = new Nonterminal('F');
		H = new Nonterminal('H');
		G = new Nonterminal('G');
		B = new Nonterminal('B');
		M = new Nonterminal('M');
		O = new Nonterminal('O');
		N = new Nonterminal('N');
		U = new Nonterminal('U');
		Z = new Nonterminal('Z');

		dot = new Terminal('.');
		open = new Terminal('(');
		close = new Terminal(')');


		// Note that the order of these rules is fixed by constraints set by the used Earley-Parser. DO NOT CHANGE!
		Rule[] rules = new Rule[] {
			new Rule(Z, new GrammarSymbol[] {dot}, Probability.fromNumber(values[iter++])),          // 0
			min_HEL == 1 ? new Rule(A, new GrammarSymbol[] {open, L, close}, Probability.fromNumber(values[iter++])) : new Rule(A, new GrammarSymbol[] {open, open, L, close, close}, Probability.fromNumber(values[iter++])),
			new Rule(P, new GrammarSymbol[] {open, L, close}, Probability.fromNumber(values[iter++])),
			new Rule(C, new GrammarSymbol[]{Z, C},    Probability.fromNumber(values[iter++])),
			new Rule(C, new GrammarSymbol[]{Z},       Probability.fromNumber(values[iter++])),
			new Rule(H, new GrammarSymbol[]{Z, H},    Probability.fromNumber(values[iter++])),  // 5
			new Rule(H, new GrammarSymbol[]{Z},       Probability.fromNumber(values[iter++])),
			new Rule(B, new GrammarSymbol[]{Z, B},    Probability.fromNumber(values[iter++])),
			new Rule(B, new GrammarSymbol[]{Z},       Probability.fromNumber(values[iter++])),
			new Rule(U, new GrammarSymbol[]{Z, U},    Probability.fromNumber(values[iter++])),
			new Rule(U, new GrammarSymbol[]{},        Probability.fromNumber(values[iter++])),  // 10
			new Rule(T, new GrammarSymbol[]{C},       Probability.fromNumber(values[iter++])),
			new Rule(T, new GrammarSymbol[]{A},       Probability.fromNumber(values[iter++])),
			new Rule(T, new GrammarSymbol[]{C, A},    Probability.fromNumber(values[iter++])),
			new Rule(T, new GrammarSymbol[]{A, T},    Probability.fromNumber(values[iter++])),
			new Rule(T, new GrammarSymbol[]{C, A, T}, Probability.fromNumber(values[iter++])),  // 15
			min_HL == 1 ? new Rule(F, new GrammarSymbol[]{H},       Probability.fromNumber(values[iter++])): new Rule(F, new GrammarSymbol[]{Z, Z, H}, Probability.fromNumber(values[iter++])),
			new Rule(G, new GrammarSymbol[]{B, A},    Probability.fromNumber(values[iter++])),
			new Rule(G, new GrammarSymbol[]{A, B},    Probability.fromNumber(values[iter++])),
			new Rule(G, new GrammarSymbol[]{B, A, B}, Probability.fromNumber(values[iter++])),
			new Rule(M, new GrammarSymbol[]{U, A, O}, Probability.fromNumber(values[iter++])),       // 20
			new Rule(O, new GrammarSymbol[]{U, A, N}, Probability.fromNumber(values[iter++])),
			new Rule(N, new GrammarSymbol[]{U, A, N}, Probability.fromNumber(values[iter++])),
			new Rule(N, new GrammarSymbol[]{U},       Probability.fromNumber(values[iter++])),
			new Rule(L, new GrammarSymbol[]{F},       Probability.fromNumber(values[iter++])),
			new Rule(L, new GrammarSymbol[]{P},       Probability.fromNumber(values[iter++])),  // 25
			new Rule(L, new GrammarSymbol[]{G},       Probability.fromNumber(values[iter++])),
			new Rule(L, new GrammarSymbol[]{M},       Probability.fromNumber(values[iter++])),
			new Rule(S, new GrammarSymbol[]{T},       Probability.fromNumber(values[iter++]))
		};

		Gs = new SCFG(15, new Terminal[] {dot, open, close}, rules);
	}
}
