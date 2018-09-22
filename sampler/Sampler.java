package sampler;

import biology.Base;
import biology.Secondary;
import biology.SecondarySequence;
import grammar.Monster;
import grammar.Nonterminal;
import grammar.Rule;
import math.Probability;
import preprocessor.StochasticModel;
import supplement.PairDict;


import java.util.Arrays;

public class Sampler{
	private boolean verbose;
	private Base[] sequence;
	private Secondary[] folded;
	private int[] brackets;
	private double sampleprobability;
	private boolean stembuilding;
	private boolean exteriorloopfinished;
	private StochasticModel model;
	private SetSampler sam;
	private final int max_hairpin = 30;
	private final int max_bulge = 30;
	private final int max_strand = 30;
	private int[][][] firstPossMLStarts; //firstPossMLStarts[k-1][start][end] = min start position for k-th multiloop element ending at end and starting at least at start
	private boolean[][] furtherPSpossible; //furtherPSpossible[start][end] = true iff further paired substructure possible between start and end
	private PairDict pseudohelices;
	private int lasthelix_start;
	private int lasthelix_end;
	private int lasthelix_length;

	public Sampler(Base[] sequence, StochasticModel model, long seed, boolean verbose){
		this.verbose = verbose;
		this.sequence = sequence;
		this.model = model;
		sam = new SetSampler(seed);
		folded = new Secondary[sequence.length];
		Arrays.fill(folded, Secondary.undecided);
		brackets = new int[sequence.length];
		Arrays.fill(brackets, -1);
		sampleprobability = Probability.sure;
		exteriorloopfinished = false;
		stembuilding = false;
		pseudohelices = new PairDict(sequence.length);
	}

	public void prepare(){
		// cf. p.248
		firstPossMLStarts = new int[3][sequence.length][sequence.length];
		for(int[][] i : firstPossMLStarts){
			for(int[] j : i){
				Arrays.fill(j, sequence.length); // "infinity"
			}
		}

		for(int end=sequence.length-1; end>=0; end--){
			int min1 = sequence.length;
			int min2 = sequence.length;
			int min3 = sequence.length;
			for(int minallowedstart=end; minallowedstart>=0; minallowedstart--){
				double inside, outside;

				// k = 1
				if(minallowedstart <= end-2*Monster.min_PS+1){
					inside = model.inside(Monster.M, minallowedstart, end);
					outside = model.outside(Monster.M, minallowedstart, end);
					if(Probability.isPossible(inside) && Probability.isPossible(outside)){
						min1 = minallowedstart;
					}
					firstPossMLStarts[0][minallowedstart][end] = min1;
				}

				// k = 2
				if(minallowedstart <= end-Monster.min_PS+1){
					inside = model.inside(Monster.O, minallowedstart, end);
					outside = model.outside(Monster.O, minallowedstart, end);
					if(Probability.isPossible(inside) && Probability.isPossible(outside)){
						min2 = minallowedstart;
					}
					firstPossMLStarts[1][minallowedstart][end] = min2;
				}

				// k = 3
				// if(minallowedstart <= end+1) <- always true
				inside = model.inside(Monster.N, minallowedstart, end);
				outside = model.outside(Monster.N, minallowedstart, end);
				if(Probability.isPossible(inside) && Probability.isPossible(outside)){
					min3 = minallowedstart;
				}
				firstPossMLStarts[2][minallowedstart][end] = min3;
			}
		}


		furtherPSpossible = new boolean[sequence.length][sequence.length]; // assuming java initializes that with false
		// "simple" dynamic programming - that's what she said
		for(int d=Monster.min_PS-1; d<sequence.length; d++){
			for(int startpos=0; startpos<sequence.length-d; startpos++){
				int endpos = startpos + d;
				double in = model.inside(Monster.A, startpos, endpos);
				double out = model.outside(Monster.A, startpos, endpos);
				boolean hit = Probability.isPossible(Probability.multiply(in, out));
				if(d == Monster.min_PS-1){
					furtherPSpossible[startpos][endpos] = hit;
				}else{
					furtherPSpossible[startpos][endpos] = hit || furtherPSpossible[startpos+1][endpos] || furtherPSpossible[startpos][endpos-1];
				}
			}
		}
	}

	public void fold_everything(){
		// S -> T
		accountRule(28);
		int unpstart = 0;
		int unpend = sequence.length-1;
		do{
			foldPairedSubstructure(unpstart, unpend, false);

			for(unpstart=0; unpstart<sequence.length; unpstart++){
				if(folded[unpstart] == Secondary.undecided) break;
			}
			for(unpend=unpstart; unpend<sequence.length; unpend++){
				if(folded[unpend] != Secondary.undecided) break;
			}
			unpend--;
		}while(unpstart != sequence.length);

		if(!exteriorloopfinished){
			// we may have "forgotten" to close the loop
			// T -> C
			accountRule(11);
			exteriorloopfinished = true;
		}

		// multiply emission probabilities
		for(int i=0; i<sequence.length; i++){
			if(folded[i] == Secondary.unpaired){
				sampleprobability = Probability.multiply(sampleprobability, Monster.em1[sequence[i].getId()]);
			}else if(folded[i] == Secondary.open){
				int j = brackets[i];
				sampleprobability = Probability.multiply(sampleprobability, Monster.em2[sequence[i].getId()][sequence[j].getId()]);
			}
		}
	}

	private SampleStatus foldPairedSubstructure(int start, int end, boolean in_ml){
		SampleStatus status = foldHairpinLoop(start, end, in_ml);
		boolean completed = status.isCompleted();
		while(!completed){
			status = expandStructure(start, end, status);
			completed = status.isCompleted();
		}
		return status;
	}

	private SampleStatus foldHairpinLoop(int start, int end, boolean in_ml){
		sam.restart();

		//pcHL, cf. p.253
		for(int i = start + Monster.min_HEL; i <= end - Monster.min_HEL; i++){
			int j = i + Monster.min_HL - 1;
			if(!checkUnfolded(i - Monster.min_HEL, j + Monster.min_HEL)) continue;
			for(; j <= i + max_hairpin - 1 && j <= end - Monster.min_HEL; j++){
				if(folded[j + Monster.min_HEL] != Secondary.undecided) break;
				double p1 = model.outside(Monster.L, i, j);
				double p2 = model.inside(Monster.F, i, j);
				double p3 = Monster.Gs.getRule(24).getProbability();
				double prob = Probability.multiply(p1, Probability.multiply(p2, p3));
				if(Probability.isPossible(prob)){
					//System.out.println("HL possible at " + i + ", " + j + " with prob " + prob);
					sam.add(new Structure(i, j, prob, Structure.hairpinloop, Structure.hl_no));
				}
			}
		}

		if(sam.size() == 0){
			debug_s("No hairpin possible on " + start + ", " + end);
			return completeStructure(start, end, in_ml);
		}else{
			debug_s(sam.size()+" hairpins possible on "+start+", "+end);
		}

		Structure loop = sam.uniformsample();
		int loopstart = loop.getStart();
		int loopend = loop.getEnd();
		debug_s("drawn HL ("+loopstart+", "+loopend+") -> " + loop.getProbability());

		for(int i=loopstart; i<=loopend; i++){
			folded[i] = Secondary.unpaired;
		}
		lasthelix_length = 0;
		lasthelix_start = loop.getStart();
		lasthelix_end = loop.getEnd();

		// L -> F
		// (F -> Z^(min_HL-1)H)
		// (loopend-loopstart+1)-(min_HL-1)-1 times H -> ZH
		// H -> Z
		accountRule(24);
		accountRule(16);
		accountRule(5, (loopend-loopstart+1)-(Monster.min_HL-1)-1);
		accountRule(6);

		return new SampleStatus(loop.getStart(), loop.getEnd(), false);
	}

	private SampleStatus expandStructure(int start, int end, SampleStatus status){
		int ii = status.getII(), jj = status.getJJ();

		//left(right)Bound is the first folded position left (right) of ii (jj). Possibly -1 (sequence.length).
		int leftBound = ii-1, rightBound = jj+1;
		while(leftBound>=0 && folded[leftBound] == Secondary.undecided){
			leftBound--;
		}
		while(rightBound<sequence.length && folded[rightBound] == Secondary.undecided){
			rightBound++;
		}
		// this should not occur according to the author - but it does. Rarely, but it does.
		//assert leftBound < 0 || folded[leftBound] == Secondary.open || folded[leftBound] == Secondary.close;

		//left(right)Bound_p is the first paired position left (right) of ii (jj). Possibly -1 (sequence.length).
		int leftBound_p = ii-1, rightBound_p = jj+1;
		while(leftBound_p>=0 && folded[leftBound_p] != Secondary.open && folded[leftBound_p] != Secondary.close){
			leftBound_p--;
		}
		while(rightBound_p<sequence.length && folded[rightBound_p] != Secondary.open && folded[rightBound_p] != Secondary.close){
			rightBound_p++;
		}

		int lowerbound;
		if(lasthelix_length == 0 && lasthelix_start == ii){
			lowerbound = ii-1;
		}else{
			lowerbound = leftBound;
		}

		int upperbound;
		if((lasthelix_length == 0 && lasthelix_end == jj) || pseudohelices.existsEnd(jj)){
			upperbound = jj+1;
		}else{
			upperbound = rightBound;
		}

		// ml_lowerupperbound is the endpoint of the last helix that surrounds [ii, jj]
		// ml_upperupperbound is the endpoint of the first helix that surrounds [ii, jj]
		// Note: The endpoint of a helix is always outside (lucky coincidence!)
		// TODO: kind of makes no sense. lower ist always greater than upper, is it?
		int ml_lowerupperbound = -1;
		int ml_upperupperbound = -1;
		for(int i=jj+1; i<sequence.length; i++){
			if((folded[i] == Secondary.close && brackets[i] < ii) || pseudohelices.existsEnd(i)){
				ml_lowerupperbound = i;
			}
			// in the meantime Abs[] has been applied to the pseudohelices
			if((folded[i] == Secondary.close && brackets[i] < ii) || (pseudohelices.existsEnd(i) && pseudohelices.minStart(i) < ii)){
				if(ml_upperupperbound == -1) ml_upperupperbound = i;
			}
		}
		ml_lowerupperbound = Math.max(ml_lowerupperbound+1, Math.max(Monster.min_HEL+2*Monster.min_PS, jj+1));
		if(ml_upperupperbound == -1) ml_upperupperbound = sequence.length;
		ml_upperupperbound = Math.min(ml_upperupperbound+1, sequence.length-Monster.min_HEL);
		if(ml_lowerupperbound > ml_upperupperbound) ml_lowerupperbound = ml_upperupperbound;
		ml_lowerupperbound = Math.min(ml_lowerupperbound, end+1);
		ml_upperupperbound = Math.min(ml_upperupperbound, end+1);

		int ml_lowerbound = lowerbound;
		if(ml_lowerbound > -1){
			int substructstart = brackets[ml_lowerbound];
			while(substructstart > 0 && folded[substructstart - 1] == Secondary.unpaired){
				//TODO: we might need to actually store the unpregs start points and stop at the first one. or maybe not. whatever. It's probably fine.
				substructstart--;
			}
			if(substructstart < 0 || !pseudohelices.existsStart(substructstart)) ml_lowerbound += Monster.min_HEL;
		}else{
			ml_lowerbound += Monster.min_HEL;
		}
		ml_lowerbound = Math.max(ml_lowerbound, start-1);

		int el_lowerbound = lowerbound;
		if(el_lowerbound >= 0 && pseudohelices.existsEnd(el_lowerbound)){
			el_lowerbound += Monster.min_HEL;
		}
		if(el_lowerbound < -1) el_lowerbound = -1;

		int el_upperbound = upperbound;
		if(el_upperbound < sequence.length && pseudohelices.existsStart(el_upperbound)){
			el_upperbound -= Monster.min_HEL;
		}
		if(el_upperbound > end+1) el_upperbound = end + 1;

		boolean possiblechoiceforEL = arepossiblechoicesforEL(ii, jj, end);

		//precpairing is true iff a paired substructure ending at ii-1 can be folded on [start, ii-1] or already exists
		boolean precpairing = false;
		if(el_lowerbound == ii-1) precpairing = true;
		//trying to fold something on (k, ii-1)
		for(int k=el_lowerbound+1; !precpairing && k<=(ii-1)-(Monster.min_PS-1); k++){
			double prob;
			if((ii-1 - k)+1 < Monster.min_HL+2*Monster.min_HEL) continue;

			prob = Probability.multiply(model.inside(Monster.A, k, ii-1), model.outside(Monster.A, k, ii-1));
			if(Probability.isPossible(prob)){
				precpairing = true;
			}
		}


		sam.restart();
		pcSP(start, end, ii, jj);
		if(0 <= lasthelix_length && lasthelix_length < Monster.min_HEL){
			// Then we only want basepairs - to close the structure
			stembuilding = true;
			if(sam.size() == 0){
				// Uh oh - we don't have any. Yet.
				int missingbps = Monster.min_HEL - lasthelix_length;
				sam.add(new Structure(ii-missingbps, jj+missingbps, Probability.sure, Structure.stackedpair, Structure.bp_special));
			}
		}else{
			pcBI(start, end, ii, jj, leftBound, rightBound);
			pcML(ii, jj, ml_lowerbound, ml_lowerupperbound, ml_upperupperbound);
			if(possiblechoiceforEL){
				pcEL(ii, jj, el_lowerbound, el_upperbound, precpairing);
			}
		}

		if(sam.size() == 0){
			debug_s("No expanding possible on ("+start+", "+end+") for ("+ii+", "+jj+")");
			if(stembuilding){
				// stem is finished. We only accounted rules L -> P -> (L), so now we must take min_HEL of
				// them back because the closing pair is always built by A -> (^min_HEL  L  )^min_HEL
				revertRule(25, Monster.min_HEL);
				revertRule(2, Monster.min_HEL);
				accountRule(1);
				stembuilding = false;
			}
			return new SampleStatus(ii, jj, true);
		}

		Structure choice = sam.uniformsample();
		if(choice.getType() != Structure.stackedpair){
			if(stembuilding){
				// stem is finished. We only accounted rules L -> P -> (L), so now we must take min_HEL of
				// them back because the closing pair is always built by A -> (^min_HEL  L  )^min_HEL
				revertRule(25, Monster.min_HEL);
				revertRule(2, Monster.min_HEL);
				accountRule(1);
			}
			stembuilding = false;
		}
		switch(choice.getType()){
		case Structure.stackedpair:
			debug_s("drawn SP ("+choice.getStart()+", "+choice.getEnd()+") -> " + choice.getProbability());
			return constructBasepairs(choice);
		case Structure.bulgeinterior:
			debug_s("drawn BI ("+choice.getStart()+", "+choice.getEnd()+") -> " + choice.getProbability());
			return constructBulgeInterior(choice, status);
		case Structure.multiloop:
			debug_s("drawn ML ("+choice.getStart()+", "+choice.getEnd()+") -> " + choice.getProbability());
			return constructMultiloop(start, end, choice, status);
		default:  // Structure.exteriorloop
			debug_s("drawn EL ("+choice.getStart()+", "+choice.getEnd()+") -> " + choice.getProbability());
			return constructExterior(choice, status);
		}
	}

	private void pcSP(int start, int end, int ii, int jj){
		// cf. p.253
		int i, j;

		//pcSP_A
		i = ii - Monster.min_HEL;
		j = jj + Monster.min_HEL;
		if(start <= i && i <= j && j <= end && checkUnfolded(i, i+Monster.min_HEL-1) && checkUnfolded(j-Monster.min_HEL+1, j)){
			double p1 = model.inside(Monster.A, i, j);
			double p2 = model.outside(Monster.A, i, j);
			double prob = Probability.multiply(p1, p2);
			if(Probability.isPossible(prob)){
				sam.add(new Structure(i, j, prob, Structure.stackedpair, Structure.bp_a));
			}
		}

		//pcSP_P
		i = ii - 1;
		j = jj + 1;
		if(start <= i && i <= j && j <= end && folded[i] == Secondary.undecided && folded[j] == Secondary.undecided){
			double p1 = model.inside(Monster.P, i, j);
			double p2 = model.outside(Monster.L, i, j);
			double p3 = Monster.Gs.getRule(25).getProbability();
			double prob = Probability.multiply(Probability.multiply(p1, p2), p3);
			if(Probability.isPossible(prob)){
				sam.add(new Structure(i, j, prob, Structure.stackedpair, Structure.bp_p));
			}
		}
	}

	private void pcBI(int start, int end, int ii, int jj, int leftBound, int rightBound){
		int i, j;

		//pcBI_(G->BA)
		j = jj;
		if(j <= end-Monster.min_HEL && j+Monster.min_HEL < rightBound){
			for(int k=1; k<=max_bulge; k++){
				i = ii - k;
				if(i < start+Monster.min_HEL || i-Monster.min_HEL <= leftBound) break;
				double p1 = model.outside(Monster.G, i, j);
				double p2 = model.inside(Monster.B, i, ii - 1);
				double p3 = model.inside(Monster.A, ii, j);
				double p4 = Monster.Gs.getRule(17).getProbability();
				double prob = Probability.multiply(p1, Probability.multiply(p2, Probability.multiply(p3, p4)));
				if(Probability.isPossible(prob)){
					sam.add(new Structure(i, j, prob, Structure.bulgeinterior, Structure.bi_ba));
				}
			}
		}

		//pcBI_(G->AB)
		i = ii;
		if(i >= start+Monster.min_HEL && i-Monster.min_HEL > leftBound){
			for(int k=1; k<=max_bulge; k++){
				j = jj + k;
				if(j > end-Monster.min_HEL || j+Monster.min_HEL >= rightBound) break;
				double p1 = model.outside(Monster.G, i, j);
				double p2 = model.inside(Monster.A, i, jj);
				double p3 = model.inside(Monster.B, jj+1, j);
				double p4 = Monster.Gs.getRule(18).getProbability();
				double prob = Probability.multiply(p1, Probability.multiply(p2, Probability.multiply(p3, p4)));
				if(Probability.isPossible(prob)){
					sam.add(new Structure(i, j, prob, Structure.bulgeinterior, Structure.bi_ab));
				}
			}
		}

		//pcBI_(G->BAB)
		for(int k1=1; k1<=max_bulge; k1++){
			i = ii - k1;
			if(i < start+Monster.min_HEL || i-Monster.min_HEL <= leftBound) break;
			for(int k2=1; k2<=max_bulge; k2++){
				j = jj + k2;
				if(j > end-Monster.min_HEL || j+Monster.min_HEL >= rightBound) break;
				double p1 = model.outside(Monster.G, i, j);
				double p2 = model.inside(Monster.B, i, ii-1);
				double p3 = model.inside(Monster.A, ii, jj);
				double p4 = model.inside(Monster.B, jj+1, j);
				double p5 = Monster.Gs.getRule(19).getProbability();
				double prob = Probability.multiply(p1, Probability.multiply(p2, Probability.multiply(p3, Probability.multiply(p4, p5))));
				if(Probability.isPossible(prob)){
					sam.add(new Structure(i, j, prob, Structure.bulgeinterior, Structure.bi_bab));
				}
			}
		}
	}

	private void pcML(int ii, int jj, int ml_lowerbound, int ml_lowerupperbound, int ml_upperupperbound){
		int i, j;
		//pcML
		// usedpositions is the set of *all* used positions
		// freepositions consists of jj, lowerupperbound-1, upperupperbound-1 and *all* undecided positions from lowerupperbound-1 to upperupperbound-1
		// freepositions is sorted, if that is important

		//k = 1
		for(j=Math.max(jj+Monster.min_PS, ml_lowerupperbound-1); j<=ml_upperupperbound-1; j++){
			if(folded[j] != Secondary.undecided) continue;
			int iend = Math.max(ii-max_strand, firstPossMLStarts[0][Math.min(ml_lowerbound+1, sequence.length-1)][j]);
			for(i=ii; i>=iend; i--){
				if(ispossiblechoiceforML(ml_lowerbound+1, jj, i, j, 1)){
					double prob = ioRevMk(i, j, ii, jj, 1);
					if(Probability.isPossible(prob)){
						sam.add(new Structure(i, j, prob, Structure.multiloop, Structure.ml_1));
						//debug("added ML1 "+i+", "+j+" -> "+prob);
					}
				}
			}
		}

		//k = 2
		for(j=Math.max(jj, ml_lowerupperbound-1); j<=ml_upperupperbound-1; j++){
			if(folded[j] != Secondary.undecided) continue;
			int iend = Math.max(ii-max_strand, firstPossMLStarts[1][Math.min(ml_lowerbound+1 + Monster.min_PS, sequence.length-1)][j]);
			for(i=ii; i>=iend; i--){
				if(ispossiblechoiceforML(ml_lowerbound+1, jj, i, j, 2)){
					double prob = ioRevMk(i, j, ii, jj, 2);
					if(Probability.isPossible(prob)){
						sam.add(new Structure(i, j, prob, Structure.multiloop, Structure.ml_2));
						//debug("added ML2 "+i+", "+j+" -> "+prob);
					}
				}
			}
		}

		//k = 3
		for(j=Math.max(jj, ml_lowerupperbound-1); j<=ml_upperupperbound-1; j++){
			if(folded[j] != Secondary.undecided) continue;
			int iend = Math.max(ii-max_strand, firstPossMLStarts[2][Math.min(ml_lowerbound+1 + 2*Monster.min_PS, sequence.length-1)][j]);
			for(i=ii; i>=iend; i--){
				if(ispossiblechoiceforML(ml_lowerbound+1, jj, i, j, 3)){
					double prob = ioRevMk(i, j, ii, jj, 3);
					if(Probability.isPossible(prob)){
						sam.add(new Structure(i, j, prob, Structure.multiloop, Structure.ml_3));
						//debug("added ML3 "+i+", "+j+" -> "+prob);
					}
				}
			}
		}
	}

	private void pcEL(int ii, int jj, int el_lowerbound, int el_upperbound, boolean precpairing){
		//pcE1 (T->C)
		// can not occur

		//pcE2 (T->A)
		if(ii >= el_lowerbound+1 && jj == sequence.length-1 && precpairing){
			double p1 = model.outside(Monster.T, ii, jj);
			double p2 = model.inside(Monster.A, ii, jj);
			double p3 = Monster.Gs.getRule(12).getProbability();
			double prob = Probability.multiply(p1, Probability.multiply(p2, p3));
			if(jj-ii+1 >= Monster.min_PS && Probability.isPossible(prob)){
				sam.add(new Structure(ii, jj, prob, Structure.exteriorloop, Structure.el_a));
			}
		}

		//pcE3 (T->CA)
		if(jj == sequence.length-1){
			for(int i=el_lowerbound+1; i<=ii-1; i++){
				double p1 = model.outside(Monster.T, i, jj);
				double p2 = model.inside(Monster.C, i, ii-1);
				double p3 = model.inside(Monster.A, ii, jj);
				double p4 = Monster.Gs.getRule(13).getProbability();
				double prob = Probability.multiply(p1, Probability.multiply(p2, Probability.multiply(p3, p4)));
				if(jj-i+1 >= Monster.min_PS+1 && Probability.isPossible(prob)){
					sam.add(new Structure(i, jj, prob, Structure.exteriorloop, Structure.el_ca));
				}
			}
		}

		//pcE4 (T->AT)
		if(ii >= el_lowerbound+1 && jj <= Math.min(el_upperbound-1, sequence.length-2) && precpairing){
			double p1 = model.outside(Monster.T, ii, sequence.length-1);
			double p2 = model.inside(Monster.A, ii, jj);
			double p3 = model.inside(Monster.T, jj+1, sequence.length-1);
			double p4 = Monster.Gs.getRule(14).getProbability();
			double prob = Probability.multiply(p1, Probability.multiply(p2, Probability.multiply(p3, p4)));
			if(sequence.length-1-ii+1 >= Monster.min_PS+1 && Probability.isPossible(prob)){
				sam.add(new Structure(ii, sequence.length-1, prob, Structure.exteriorloop, Structure.el_at));
			}
		}

		//pcE5 (T->CAT)
		if(jj <= Math.min(el_upperbound-1, sequence.length-2)){
			for(int i=el_lowerbound+1; i<=ii-1; i++){
				double p1 = model.outside(Monster.T, i, sequence.length-1);
				double p2 = model.inside(Monster.C, i, ii-1);
				double p3 = model.inside(Monster.A, ii, jj);
				double p4 = model.inside(Monster.T, jj+1, sequence.length-1);
				double p5 = Monster.Gs.getRule(15).getProbability();
				double prob = Probability.multiply(p1, Probability.multiply(p2, Probability.multiply(p3, Probability.multiply(p4, p5))));
				if(sequence.length-1-i+1 >= Monster.min_PS+2 && Probability.isPossible(prob)){
					sam.add(new Structure(i, sequence.length-1, prob, Structure.exteriorloop, Structure.el_cat));
				}
			}
		}
	}

	private SampleStatus completeStructure(int start, int end, boolean in_ml){
		debug_s("Completing " + start + ", " + end + (in_ml ? " inside ML" : ""));

		for(int k=start; k<=end; k++){
			if(folded[k] == Secondary.undecided){
				folded[k] = Secondary.unpaired;
				if(in_ml){
					// U -> ZU
					accountRule(9);
					// We already do each U -> eps in <constructMultiloop>
				}else{
					if((k < end && folded[k+1] == Secondary.undecided) || (k < sequence.length-1 && folded[k+1] == Secondary.unpaired) || (start > 0 && folded[start-1] == Secondary.unpaired)){
						// C -> ZC
						accountRule(3);
					}else{
						// C -> Z
						accountRule(4);
					}
				}
			}
		}
		return new SampleStatus(-1, sequence.length, true);
	}

	private SampleStatus constructBasepairs(Structure struct){
		// cf. p.251
		int i = struct.getStart();
		int j = struct.getEnd();

		assert lasthelix_length != -1;
		switch(struct.getSubtype()){
		case Structure.bp_a:
			for(int k=0; k<Monster.min_HEL; k++){
				folded[i+k] = Secondary.open;
				folded[j-k] = Secondary.close;
				brackets[i+k] = j-k;
				brackets[j-k] = i+k;
			}
			lasthelix_start = i;
			lasthelix_end = j;
			lasthelix_length += Monster.min_HEL;

			// We always do as if the helix was extended by L -> P -> (L) and when it is
			// finished, we take min_HEL of it back and replace it by one A -> (^min_HEL  L  )^min_HEL

			// L -> P
			// P -> (L)
			accountRule(25, Monster.min_HEL);
			accountRule(2, Monster.min_HEL);
			break;
		case Structure.bp_p:
			folded[i] = Secondary.open;
			folded[j] = Secondary.close;
			brackets[i] = j;
			brackets[j] = i;
			lasthelix_start = i;
			lasthelix_end = j;
			lasthelix_length += 1;

			// L -> P
			// P -> (L)
			accountRule(25);
			accountRule(2);
			break;
		case Structure.bp_special:
			//delete last helix. We replace it with this one.
			for(int k=0; k<lasthelix_length; k++){
				folded[lasthelix_start+k] = Secondary.undecided;
				folded[lasthelix_end-k] = Secondary.undecided;
				brackets[lasthelix_start+k] = -1;
				brackets[lasthelix_end-k] = -1;
			}
			// we sampled with the wrong rule - we revert that
			// L -> P
			revertRule(25, lasthelix_length);
			// P -> (L)
			revertRule(2, lasthelix_length);

			for(int k=0; k<Monster.min_HEL; k++){
				folded[i+k] = Secondary.open;
				folded[j-k] = Secondary.close;
				brackets[i+k] = j-k;
				brackets[j-k] = i+k;
			}
			lasthelix_start = i;
			lasthelix_end = j;
			lasthelix_length = Monster.min_HEL;

			// see bp_a case above
			accountRule(25, Monster.min_HEL);
			break;
		}

		return new SampleStatus(i, j, false);
	}

	private SampleStatus constructBulgeInterior(Structure struct, SampleStatus status){
		// cf. p.251
		int i = struct.getStart(), j = struct.getEnd();
		int ii = status.getII(), jj = status.getJJ();

		switch(struct.getSubtype()){
		case Structure.bi_ba:
			for(int k=0; k<(ii-1)-i+1; k++){
				folded[i+k] = Secondary.unpaired;
			}
			// L -> G
			// G -> BA
			// ii-i-1 times B -> ZB
			// B -> Z
			accountRule(26);
			accountRule(17);
			accountRule(7, ii-i-1);
			accountRule(8);
			break;
		case Structure.bi_ab:
			for(int k=0; k<j-(jj+1)+1; k++){
				folded[jj+1+k] = Secondary.unpaired;
			}
			// L -> G
			// G -> AB
			// j-jj-1 times B -> ZB
			// B -> Z
			accountRule(26);
			accountRule(18);
			accountRule(7, j-jj-1);
			accountRule(8);

			break;
		case Structure.bi_bab:
			for(int k=0; k<(ii-1)-i+1; k++){
				folded[i+k] = Secondary.unpaired;
			}
			for(int k=0; k<j-(jj+1)+1; k++){
				folded[jj+1+k] = Secondary.unpaired;
			}
			// L -> G
			// G -> BAB
			// (ii-i-1) + (j-jj-1) times B -> ZB
			// 2 times B -> Z
			accountRule(26);
			accountRule(19);
			accountRule(7, (ii-i-1) + (j-jj-1));
			accountRule(8, 2);

			break;
		}
		lasthelix_start = i;
		lasthelix_end = j;
		lasthelix_length = 0;

		return new SampleStatus(i, j, false);
	}

	private SampleStatus constructExterior(Structure struct, SampleStatus status){
		// cf. p.251
		int i = struct.getStart(), j = struct.getEnd();
		int ii = status.getII();

		int subtype = struct.getSubtype();

		if(subtype == Structure.el_c){
			accountRule(11);
			exteriorloopfinished = true;
		}else if(subtype == Structure.el_ca){
			accountRule(13);
			exteriorloopfinished = true;
		}else if(subtype == Structure.el_cat){
			accountRule(15);
		}else if(subtype == Structure.el_a){
			accountRule(12);
			exteriorloopfinished = true;
		}else if(subtype == Structure.el_at){
			accountRule(14);
		}

		if(subtype == Structure.el_c || subtype == Structure.el_ca || subtype == Structure.el_cat){
			for(int k=0; k<(ii-1)-i+1; k++){
				folded[i+k] = Secondary.unpaired;
			}
			// ii-i-1 times C -> ZC
			// C -> Z
			accountRule(3, ii-i-1);
			accountRule(4);
		}
		// "Nothing to do" for subtype a

		return new SampleStatus(i, j, true);
	}

	private SampleStatus constructMultiloop(int start, int end, Structure struct, SampleStatus status){
		// cf. p.252
		int structstart = struct.getStart();
		int ii = status.getII(), jj;
		for(int k=0; k<(ii-1)-structstart+1; k++){
			folded[structstart+k] = Secondary.unpaired;
		}
		// ii-structstart times U -> ZU
		accountRule(9, ii-structstart);
		// see below for U -> eps

		ii = struct.getStart();
		jj = struct.getEnd();
		if(pseudohelices.exists(ii, jj)){
			return new SampleStatus(ii, jj, true);
		}
		pseudohelices.add(ii, jj);
		lasthelix_start = ii;
		lasthelix_end = jj;
		lasthelix_length = -1;

		int mlstart, mlend;
		if(struct.getSubtype() == Structure.ml_1){
			mlstart = ii;
			mlend = jj;
		}else{
			mlend = jj;
			mlstart = ii-1;
			while(mlstart >= 0 && folded[mlstart] != Secondary.open && folded[mlstart] != Secondary.close){
				mlstart--;
			}
			mlstart += 1;
			mlstart = firstPossMLStarts[0][mlstart + Monster.min_HEL][mlend];
		}
		if(mlend == end){
			int k;
			if(struct.getSubtype() == Structure.ml_1){
				k = 1;
			}else if(struct.getSubtype() == Structure.ml_2){
				k = 2;
			}else{  // ml_3
				k = 3;
			}
			return new SampleStatus(ii, jj, true, k);
		}

		int numberofsubstructs = 1;
		while(! (ii == -1 && jj == sequence.length)){
			debug_s("Starting recursive folding of " + mlstart + ", "+ mlend);
			SampleStatus res = foldPairedSubstructure(mlstart, mlend, true);
			numberofsubstructs += 1;
			debug_s("Finished recursive folding of " + mlstart + ", " + mlend);
			// assert res.getK() != -1;
			ii = res.getII();
			jj = res.getJJ();
			int minsubstructstart = pseudohelices.minStart(mlend);
			if(jj == mlend && res.getK() == 1){
				if(ii == minsubstructstart){
					mlstart = ii;
				}
			}
			if(jj == sequence.length-1){
				mlstart = lasthelix_end + 1 + Monster.min_HEL;
				if(mlstart >= sequence.length) mlstart = sequence.length - 1;
				mlstart = Math.min(firstPossMLStarts[0][mlstart][mlend], minsubstructstart);
			}
		}
		numberofsubstructs -= 1;
		ii = mlstart;
		jj = mlend;
		lasthelix_start = ii;
		lasthelix_end = jj;
		lasthelix_length = 0;

		// L -> M
		// M -> UAO
		// O -> UAN
		// (numberofsubstructs-2) times N -> UAN
		// N -> U
		// (numberofsubstructs+1) times U -> eps
		accountRule(27);
		accountRule(20);
		accountRule(21);
		accountRule(22, numberofsubstructs-2);
		accountRule(23);
		accountRule(10, numberofsubstructs+1);

		return new SampleStatus(ii, jj, false);
	}

	private boolean checkUnfolded(int start, int end){
		for(int k=start; k<=end; k++){
			if(k < 0 || k >= folded.length) continue;
			if(folded[k] != Secondary.undecided) return false;
		}
		return true;
	}

	private double ioRevMk(int i, int j, int ii, int jj, int k){
		int threshold, rule;
		Nonterminal first, last;
		if(k == 1){
			threshold = 2 * Monster.min_PS;
			rule = 20;
			first = Monster.M;
			last = Monster.O;
		}else if(k == 2){
			threshold = Monster.min_PS;
			rule = 21;
			first = Monster.O;
			last = Monster.N;
		}else if(k == 3){
			threshold = Monster.min_PS;
			rule = 22;
			first = Monster.N;
			last = Monster.N;
		}else{
			assert false;
			return 0;
		}

		if(j-i+1 < threshold) return Probability.impossible;

		double p1 = model.outside(first, i, j);
		double p2 = model.inside(Monster.U, i, ii-1);
		double p3 = model.inside(Monster.A, ii, jj);
		double p4 = model.inside(last, jj+1, j);
		double p5 = Monster.Gs.getRule(rule).getProbability();

		return Probability.multiply(p1, p2, p3, p4, p5);
	}

	private double ioRevML(int start, int end, int k){
		Nonterminal nt = null;
		if(k == 1){
			nt = Monster.M;
		}else if(k == 2){
			nt = Monster.O;
		}else if(k == 3){
			nt = Monster.N;
		}
		double p1 = model.outside(nt, start, end);
		double p2 = model.inside(nt, start, end);
		return Probability.multiply(p1, p2);
	}

	private boolean ispossiblechoiceforML(int start, int jj, int i, int j, int k){
		if(i<start + (k-1)*Monster.min_PS || (j-i+1) < (2-(k-1)) * Monster.min_PS || !Probability.isPossible(ioRevML(i, j, k))) return false;
		int mlstart = firstPossMLStarts[0][start][j];
		if(j-mlstart+1 < 2*Monster.min_PS) return false;
		//test foundation
		for(int l=1; l<=Monster.min_HEL; l++){
			if(mlstart-l < 0 || j+l >= sequence.length) return false;
			if(folded[mlstart-l] != Secondary.undecided || folded[j+l] != Secondary.undecided) return false;
		}

		if(k == 1){
			if(furtherPSpossible[jj+1][j]){
				return true;
			}
		}else if(k >= 2){
			if(furtherPSpossible[mlstart][i-1]){
				return true;
			}
		}
		return false;
	}

	private boolean arepossiblechoicesforEL(int ii, int jj, int end){
		for(int k=ii-1; k>=0; k--){
			if((folded[k] == Secondary.open && brackets[k] > jj) || (pseudohelices.existsStart(k) && pseudohelices.maxEnd(k) > jj)){
				return false;
			}
		}
		int mlstart = jj + 1 + Monster.min_HEL;
		if(mlstart >= sequence.length) mlstart = sequence.length - 1;
		mlstart = firstPossMLStarts[0][mlstart][end];
		int currstart = pseudohelices.minStart(end);
		if(currstart == sequence.length) return true;
		if (mlstart+Monster.min_PS-1 < currstart && furtherPSpossible[mlstart][currstart-1]) return true;
		return false;
	}

	private boolean validChoicekx(int k, Nonterminal X, int start, int jj, int i, int j){
		int fp1 = firstPossMLStarts[k-1][Math.min(start+(k-1)*Monster.min_PS, sequence.length-1)][j];
		int fp2 = firstPossMLStarts[0][start][j];
		return i >= start + (k-1) * Monster.min_PS &&
			(j-i+1) >= (2 - (k-1)) * Monster.min_PS &&
			Probability.isPossible(Probability.multiply(model.outside(X, i, j), model.inside(X, i, j))) &&
			i >= fp1 &&
			j - fp2 + 1 >= 2 * Monster.min_PS &&
			(
				(k == 1 && jj+1<sequence.length && furtherPSpossible[jj+1][j]) || (k >= 2 && furtherPSpossible[fp2][i-1])
			) &&
			checkUnfolded(fp2-Monster.min_HEL, fp2-1) &&
			checkUnfolded(j+1, j+Monster.min_HEL);
	}

	private void accountRule(int index){
		accountRule(index, 1);
	}

	private void revertRule(int index){
		revertRule(index, 1);
	}

	private void accountRule(int index, int number){
		Rule r = Monster.Gs.getRule(index);
		String d = r.toString();
		if(number != 1) d = number + " times " + d;
		debug_r(d);

		double p;
		if(number == 1){
			p = r.getProbability();
		}else{
			p = Probability.power(r.getProbability(), number);
		}
		sampleprobability = Probability.multiply(sampleprobability, p);
	}

	private void revertRule(int index, int number){
		Rule r = Monster.Gs.getRule(index);
		String d = r.toString();
		if(number > 1) d = number + " times " + d;
		debug_r("reverting "+d);

		double p;
		if(number == 1){
			p = r.getProbability();
		}else{
			p = Probability.power(r.getProbability(), number);
		}
		sampleprobability = Probability.divide(sampleprobability, p);
	}

	private void debug_r(String s){
		if(verbose){
			System.out.print("[R] ");
			System.out.println(s);
		}
	}

	private void debug_s(String s){
		if(verbose){
			System.out.print("[S] ");
			System.out.println(s);
		}
	}

	public SecondarySequence getFoldedSequence(){
		return new SecondarySequence(folded, brackets);
	}

	public double getProbability(){
		return sampleprobability;
	}

	public String toString(){
		String result = "";
		for(Secondary c : folded){
			result += c.toString();
		}
		return result;
	}
}
