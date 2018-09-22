package preprocessor;

import biology.Base;
import grammar.Nonterminal;
import grammar.SCFG;


public class StochasticModel{
	private Approximator approx;

	public StochasticModel(Base[] rna, SCFG g){
		approx = new Approximator(rna, g);
	}

	public void compute(){
		approx.compute_all();
	}

	public double inside(Nonterminal nt, int start, int end){
		return approx.alpha_hat(nt, start, end);
	}

	public double outside(Nonterminal nt, int start, int end){
		return approx.beta_hat(nt, start, end);
	}
}
