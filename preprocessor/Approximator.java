package preprocessor;

import biology.Base;
import grammar.SCFG;
import grammar.Monster;
import grammar.Rule;
import grammar.GrammarSymbol;
import grammar.Terminal;
import grammar.Nonterminal;
import math.Probability;

import java.util.Arrays;
import java.util.Iterator;


public class Approximator{
	private Base[] rnasequence;
	private SCFG grammar;
	private double[] rf1;         //rf[u] -> log(relative frequency of base u)
	private double[][][] rf2;     //rf[dist, p1, p2] -> log(relative frequency of base pair p1, p2 at distance dist)
	private double ema1;          //sequence-dependent averaged emission probabilities for bases
	private double[] ema2;        //sequence-dependent averaged emission probabilities for base pairs
	private double inside[][][];  //inside[dist, p, q] -> inside item of rule p, scanned by q steps with distance dist
	private double outside[][][]; //outside[dist, p, q] -> outside item of rule p, scanned by q steps with distance dist
	private double alpha[][];     // alpha[X][dist] -> nonterminal X on distance dist
	private double beta[][];      // beta[X][dist] -> nonterminal X on distance dist
	private double prod_em[][];   //prod_ema[i][j] -> product of emission probabilities for r_i, ..., r_j
	/*
	Notice:
	Arrays depending on a 'dist' parameter have length+2 or length+3 many places, that
	is from 0 to n+1 resp. from 0 to n+2. While distances 0 and n+1 and n+2 do not actualy make sense (in most cases), coding turns out to be
	easier that way. Of course the additional sentinel places are filled with 'impossible'.
	 */

	public Approximator(Base[] rna, SCFG g){
		rnasequence = rna;
		grammar = g;
	}

	private void compute_relative_frequencies(){
		// cf. p.223
		int length = rnasequence.length;
		rf1 = new double[Base.count]; // assuming java initializes that with 0.0
		rf2 = new double[length+2][Base.count][Base.count]; // assuming java initializes that with 0.0

		int sum = 0;

		for(Base b : rnasequence){
			rf1[b.getId()] += 1;
			sum += 1;
		}
		assert sum == length;
		for(int i=0; i<Base.count; i++){
			rf1[i] = Probability.fromRatio(rf1[i], sum);
		}


		for(int dist=1; dist<=length; dist++){
			sum = 0;
			for(int i=0; i<length; i++){
				int j = dist + i - 1;
				// considering positions i and j with distance dist
				if(i + 1 + Monster.min_HL <= j && j < length){
					Base b1 = rnasequence[i];
					Base b2 = rnasequence[j];
					//if(Base.validpair(b1, b2)){
						// not valid, not counted
						rf2[dist][b1.getId()][b2.getId()] += 1;
						sum += 1;
					//}
				}
			}

			for(int p1=0; p1<Base.count; p1++){
				for(int p2=0; p2<Base.count; p2++){
					if(rf2[dist][p1][p2] == 0){
						// not counted, not valid
						rf2[dist][p1][p2] = Probability.impossible;
					}else{
						rf2[dist][p1][p2] = Probability.fromRatio(rf2[dist][p1][p2], sum);
					}
				}
			}
		}
	}

	private void compute_emissions(){
		// cf. p.224
		int length = rnasequence.length;

		ema2 = new double[length+2];

		//ema1 (dist = 1)
		double result = Probability.impossible;
		for(int u=0; u<Base.count; u++){
			double summand = Probability.multiply(Monster.em1[u], rf1[u]);
			result = Probability.add(result, summand);
		}
		ema1 = result;


		//ema2 (0 <= dist <= n+1)
		ema2[0] = Probability.impossible;
		ema2[length+1] = Probability.impossible;
		for(int dist=1; dist<=length; dist++){
			result = Probability.impossible;

			for(int p1=0; p1<Base.count; p1++){
				for(int p2=0; p2<Base.count; p2++){
					double summand = Probability.multiply(Monster.em2[p1][p2], rf2[dist][p1][p2]);
					result = Probability.add(result, summand);
				}
			}

			ema2[dist] = result;
		}

		//cf. p.235
		prod_em = new double[length][length];
		for(int i=0; i<length; i++){
			// lower matrix half needs to be filled
			Arrays.fill(prod_em[i], Probability.impossible);
		}
		for(int start=0; start<length; start++){
			double current = Probability.sure;
			for(int end=start; end<length; end++){
				double em = Monster.em1[rnasequence[end].getId()];
				current = Probability.multiply(current, em);
				prod_em[start][end] = current;
			}
		}
	}

	private void compute_inside_items(){
		// cf. p.226
		int length = rnasequence.length;
		int rc = grammar.getRuleCount();

		inside = new double[length+2][rc][];
		for(int i=0; i<=length+1; i++){
			for(int j=0; j<rc; j++){
				int size = grammar.getRule(j).length();
				inside[i][j] = new double[size+1];
				Arrays.fill(inside[i][j], Probability.impossible);
			}
		}

		// Predict
		for(int p=0; p<rc; p++){
			Rule r = grammar.getRule(p);
			if(r.length() == 0){
				inside[1][p][0] = r.getProbability();
			}else{
				inside[1][p][0] = Probability.sure;
			}
		}

		for(int dist=1; dist<=length+1; dist++){
			for(int p=0; p<rc; p++){
				Rule rule = grammar.getRule(p);
				int card = rule.length();
				for(int q=1; q<=card; q++){ // q=0 has been done before - and the loop body only makes sense for q >= 1
					GrammarSymbol lastsymbol = rule.rightside_char(q-1);
					if(lastsymbol instanceof Terminal){
						//Scan
						if(lastsymbol == Monster.dot){
							inside[dist][p][q] = Probability.multiply(ema1, inside[dist-1][p][q-1]);
						}else if(lastsymbol == Monster.open){
							inside[dist][p][q] = inside[dist-1][p][q-1];
						}else{ // Monster.close
							inside[dist][p][q] = Probability.multiply(ema2[dist-1], inside[dist-1][p][q-1]);
						}

						if(q == card){
							inside[dist][p][q] = Probability.multiply(inside[dist][p][q], rule.getProbability());
						}
					}else{ // Nonterminal
						//Complete
						Nonterminal y = (Nonterminal) lastsymbol;
						double result = Probability.impossible;
						for(int k=1; k<=dist; k++){
							double factor = inside[k][p][q-1];

							double innersum = Probability.impossible;
							Iterator<Integer> it = grammar.getRules(y);
							while(it.hasNext()){
								int index = it.next();
								Rule yrule = grammar.getRule(index);
								double innersummand = inside[dist-k+1][index][yrule.length()];
								innersum = Probability.add(innersum, innersummand);
							}

							double summand = Probability.multiply(factor, innersum);
							result = Probability.add(result, summand);
						}
						inside[dist][p][q] = result;

						if(q == card){
							inside[dist][p][q] = Probability.multiply(inside[dist][p][q], rule.getProbability());
						}
					}
				}
			}
		}
	}

	private void compute_outside_items(){
		// cf. p.229
		int length = rnasequence.length;
		int rc = grammar.getRuleCount();

		outside = new double[length+3][rc][];
		for(int i=0; i<length+3; i++){
			for(int j=0; j<rc; j++){
				int size = grammar.getRule(j).length();
				outside[i][j] = new double[size+1];
				Arrays.fill(outside[i][j], Probability.impossible);
			}
		}

		// S->T is the only rule of S
		Iterator<Integer> it = grammar.getRules(Monster.S);
		outside[length+1][it.next()][1] = Probability.sure;

		for(int dist=length+1; dist>=1; dist--){
			for(int p=rc-1; p>=0; p--){
				Rule rule = grammar.getRule(p);
				int card = rule.length();
				for(int q=card; q>0; q--){  //We "do nothing" for q = 0 anyways
					GrammarSymbol lastsymbol = rule.rightside_char(q-1);
					if(lastsymbol instanceof Terminal){
						//Scan reverse
						if(lastsymbol == Monster.dot){
							outside[dist][p][q-1] = Probability.multiply(ema1, outside[dist+1][p][q]);
						}else if(lastsymbol == Monster.open){
							outside[dist][p][q-1] = outside[dist+1][p][q];
						}else{  // Monster.close
							outside[dist][p][q-1] = Probability.multiply(ema2[dist], outside[dist+1][p][q]);
						}

						if(q == card){
							outside[dist][p][q-1] = Probability.multiply(outside[dist][p][q-1], rule.getProbability());
						}
					}else{ //Nonterminal
						//Complete reverse
						double fact = Probability.sure;
						if(q == card) fact = rule.getProbability();

						for(int k=1; k<=dist; k++){
							double sum = Probability.impossible;
							Iterator<Integer> yit = grammar.getRules((Nonterminal) lastsymbol);
							while(yit.hasNext()){
								int yrule = yit.next();
								int yrulelength = grammar.getRule(yrule).length();
								sum = Probability.add(sum, inside[dist-k+1][yrule][yrulelength]);
							}
							double additive = Probability.multiply(Probability.multiply(outside[dist][p][q], sum), fact);
							outside[k][p][q-1] = Probability.add(outside[k][p][q-1], additive);

							yit = grammar.getRules((Nonterminal) lastsymbol);
							while(yit.hasNext()){
								int yrule = yit.next();
								int yrulelength = grammar.getRule(yrule).length();

								additive = Probability.multiply(Probability.multiply(outside[dist][p][q], inside[k][p][q-1]), fact);
								outside[dist-k+1][yrule][yrulelength] = Probability.add(outside[dist-k+1][yrule][yrulelength], additive);
							}
						}
					}
				}
			}
		}
	}

	private void compute_alpha(){
		int length = rnasequence.length;
		alpha = new double[grammar.getNonterminalCount()][length+1];
		Iterator<Nonterminal> ntit = grammar.getNonterminals();
		while(ntit.hasNext()){
			Nonterminal x = ntit.next();
			for(int dist=0; dist<=length; dist++){
				double result = Probability.impossible;

				Iterator<Integer> rit = grammar.getRules(x);
				while(rit.hasNext()){
					int rindex = rit.next();
					double summand = inside[dist+1][rindex][grammar.getRule(rindex).length()];
					result = Probability.add(result, summand);
				}
				alpha[x.getNumber()][dist] = result;
			}
		}
	}

	private void compute_beta(){
		int length = rnasequence.length;
		beta = new double[grammar.getNonterminalCount()][length+1];
		Iterator<Nonterminal> ntit = grammar.getNonterminals();
		while(ntit.hasNext()){
			Nonterminal x = ntit.next();
			for(int dist=0; dist<=length; dist++){
				double result = Probability.impossible;

				Iterator<Integer> rit = grammar.getRules(x);
				while(rit.hasNext()){
					int rindex = rit.next();
					double cmp = outside[dist+1][rindex][grammar.getRule(rindex).length()];
					result = Probability.max(result, cmp);
				}
				beta[x.getNumber()][dist] = result;
			}
		}
	}

	public void compute_all(){
		compute_relative_frequencies();
		compute_emissions();
		compute_inside_items();
		compute_outside_items();
		compute_alpha();
		compute_beta();
	}

	public double alpha_hat(Nonterminal x, int i, int j){
		//cf. p.235
		assert x != Monster.S && x != Monster.L && x!= Monster.H && x != Monster.Z;
		double result = alpha[x.getNumber()][j-i+1];
		double factor = Probability.sure;

		if(x == Monster.C || x == Monster.F || x == Monster.B || x == Monster.U){
			double num = prod_em[i][j];
			double denom = Probability.power(ema1, j - i + 1);
			factor = Probability.puredivide(num, denom);
		}else if(x == Monster.A && j>=Monster.min_HEL-1 && i <= rnasequence.length-Monster.min_HEL && j-i+1 >= 2*Monster.min_HEL){
			double num = Probability.sure;
			double denom = Probability.sure;
			for(int k=0; k<=Monster.min_HEL-1; k++){
				num = Probability.multiply(num, Monster.em2[rnasequence[i+k].getId()][rnasequence[j-k].getId()]);
				denom = Probability.multiply(denom, ema2[(j-k)-(i+k)+1]);
			}
			factor = Probability.puredivide(num, denom);
		}else if(x == Monster.P){
			double num = Monster.em2[rnasequence[i].getId()][rnasequence[j].getId()];
			double denom = ema2[j-i+1];
			factor = Probability.puredivide(num, denom);
		}

		return Probability.multiply(result, factor);
	}

	public double beta_hat(Nonterminal x, int i, int j){
		//cf. p.235
		assert x != Monster.S && x != Monster.C && x != Monster.P && x != Monster.F && x != Monster.H && x != Monster.B && x != Monster.U && x != Monster.Z;
		double result = beta[x.getNumber()][j-i+1];
		double factor = Probability.sure;

		if((x == Monster.L || x == Monster.G || x == Monster.M) && i >= Monster.min_HEL && j < rnasequence.length-Monster.min_HEL){
			double num = Probability.sure;
			double denom = Probability.sure;
			for(int k=0; k<=Monster.min_HEL-1; k++){
				num = Probability.multiply(num, Monster.em2[rnasequence[i - Monster.min_HEL + k].getId()][rnasequence[j + Monster.min_HEL - k].getId()]);
				denom = Probability.multiply(denom, ema2[(j + Monster.min_HEL - k) - (i - Monster.min_HEL + k) + 1]);
			}
			factor = Probability.puredivide(num, denom);
		}

		return Probability.multiply(result, factor);
	}
}
