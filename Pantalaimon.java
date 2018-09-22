import biology.Base;
import biology.RNAParser;
import biology.SecondarySequence;
import grammar.Monster;
import math.Probability;
import math.Statistics;
import preprocessor.StochasticModel;
import sampler.Sampler;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Random;


public class Pantalaimon{

	public static void main(String[] argv){
		if(argv.length != 2){
			usage();
			System.exit(0);
		}

		String conf = argv[0];
		if(conf.equals("hl1hel1")){
			Monster.init(Monster.class.getResource("monster_hl1_hel1").getPath());
		}else if(conf.equals("hl3hel2")){
			Monster.init(Monster.class.getResource("monster_hl3_hel2").getPath());
		}else{
			Monster.init(conf);
		}

		Base[] sequence = RNAParser.parse_rna(argv[1]);
		SecondarySequence solution = null;

		System.out.print("Preprocessing sequence ... ");
		StochasticModel model = new StochasticModel(sequence, Monster.Gs);
		model.compute();
		System.out.println("Done!");
		System.out.println("Given primary structure: "+RNAParser.emit(sequence));

		Random random = new Random();
		BufferedReader in = new BufferedReader(new InputStreamReader(System.in));
		String input;

		System.out.println("Ready to play! Enter command 'help' to get a helpful message");

		while(true){
			System.out.print("Enter command: > ");
			try{
				input = in.readLine();
			}catch(IOException e){
				input = "";
			}
			if(input.equals("exit")){
				break;
			}else if(input.equals("help")){
				System.out.println("List of commands:");
				System.out.println("    exit             exit the program");
				System.out.println("                     fold a secondary structure");
				System.out.println("    <seed>           fold the (unique) sequence that corresponds to the given seed");
				System.out.println("    batch <number>   fold number many secondary structures");
				System.out.println("    mp <number>      fold number many secondary structures and display the one with highest probability");
				System.out.println("    solution <seq>   supply a (suspected) solution, so foldings can be compared to it");
				System.out.println("    stat <number>    fold number many secondary structures and display the two with highest sensitivity and ppv. Only possible if a solution was supplied");
			}else if(input.length() == 0){
				long seed = random.nextLong();
				System.out.println("Using 'random' seed " + seed);
				Sampler s = new Sampler(sequence, model, seed, true);
				s.prepare();
				s.fold_everything();

				SecondarySequence predicted = s.getFoldedSequence();

				System.out.println(s.toString());
				System.out.println(RNAParser.emit(sequence));
				System.out.println("seed: "+seed+", probability: "+s.getProbability());
				if(solution != null){
					Statistics stat = predicted.compare(solution);
					System.out.println("sensitivity: "+stat.sensitivity()+", PPV: "+stat.ppv());
					System.out.println(solution.toString());
				}
			}else if(input.startsWith("batch ")){
				int number;
				try{
					number = Integer.parseInt(input.substring(6));
				}catch(NumberFormatException e){
					System.out.println("No valid number detected");
					number = 0;
				}

				for(int iter=0; iter<number; iter++){
					long seed = Math.abs(random.nextLong());

					Sampler s = new Sampler(sequence, model, seed, false);
					s.prepare();
					s.fold_everything();

					SecondarySequence predicted = s.getFoldedSequence();

					String score = " (seed: "+seed+", probability: "+s.getProbability();
					if(solution != null){
						Statistics stat = predicted.compare(solution);
						score += ", sensitivity: "+stat.sensitivity()+", PPV: "+stat.ppv();
					}
					System.out.println(s.toString() + score + ")");
				}
			}else if(input.startsWith("mp ")){
				int number;
				try{
					number = Integer.parseInt(input.substring(3));
				}catch(NumberFormatException e){
					System.out.println("No valid number detected");
					number = 0;
				}

				SecondarySequence best_seq = null;
				double best_prob = Probability.impossible;
				long best_seed = 0;

				for(int iter = 0; iter < number; iter++){
					long seed = Math.abs(random.nextLong());

					Sampler s = new Sampler(sequence, model, seed, false);
					s.prepare();
					s.fold_everything();

					SecondarySequence predicted = s.getFoldedSequence();
					double newprob = s.getProbability();
					if(Probability.greater(best_prob, newprob)){
						best_prob = newprob;
						best_seed = seed;
						best_seq = predicted;
					}
				}
				if(best_seq != null){
					System.out.println("");
					System.out.println("MP structure:");
					System.out.println(best_seq.toString());
					System.out.println(RNAParser.emit(sequence));
					System.out.println("seed: " + best_seed + ", probability: " + best_prob);
					if(solution != null){
						Statistics stat = best_seq.compare(solution);
						System.out.println("sensitivity: "+stat.sensitivity()+", PPV: "+stat.ppv());
						System.out.println(solution.toString());
					}
				}
			}else if(input.startsWith("solution ")){
				if(input.length() - 9 != sequence.length){
					System.out.println("no possible sequence found. Check the length of your string");
				}else{
					String folding = input.substring(9);
					solution = new SecondarySequence(folding);
					System.out.println("assuming solution: "+solution.toString());
				}
			}else if(input.startsWith("stat ")){
				int number;
				try{
					number = Integer.parseInt(input.substring(5));
				}catch(NumberFormatException e){
					System.out.println("no valid number detected");
					number = 0;
				}
				if(solution == null){
					System.out.println("a solution must be given to use this command");
				}else{
					Statistics best_sens_stat = null;
					Statistics best_ppv_stat = null;
					SecondarySequence best_sens_seq = null;
					SecondarySequence best_ppv_seq = null;
					long best_sens_seed = -1;
					long best_ppv_seed = -1;
					double best_ppv_prob = -1;
					double best_sens_prob = -1;

					for(int iter=0; iter<number; iter++){
						long seed = Math.abs(random.nextLong());

						Sampler s = new Sampler(sequence, model, seed, false);
						s.prepare();
						s.fold_everything();

						SecondarySequence predicted = s.getFoldedSequence();
						Statistics stat = predicted.compare(solution);

						if(best_sens_stat == null || stat.sensitivity() > best_sens_stat.sensitivity()){
							best_sens_seq = predicted;
							best_sens_seed = seed;
							best_sens_stat = stat;
							best_sens_prob = s.getProbability();
						}

						if(best_ppv_stat == null || stat.ppv() > best_ppv_stat.ppv()){
							best_ppv_seq = predicted;
							best_ppv_seed = seed;
							best_ppv_stat = stat;
							best_ppv_prob = s.getProbability();
						}
					}
					if(best_sens_seq != null){
						System.out.println();
						System.out.println("SENS structure:");
						System.out.println(best_sens_seq.toString());
						System.out.println("seed: " + best_sens_seed + ", probability: " + best_sens_prob);
						System.out.println("sensitivity: " + best_sens_stat.sensitivity() + ", PPV: " + best_sens_stat.ppv());

						System.out.println("PPV structure:");
						System.out.println(best_ppv_seq.toString());
						System.out.println("seed: " + best_ppv_seed + ", probability: " + best_ppv_prob);
						System.out.println("sensitivity: " + best_ppv_stat.sensitivity() + ", PPV: " + best_ppv_stat.ppv());

						System.out.println();
						System.out.println(RNAParser.emit(sequence));
						System.out.println(solution.toString());
					}
				}
			}else{
				long seed;
				try{
					seed = Long.parseLong(input);
					System.out.println("Using seed " + seed);
					Sampler s = new Sampler(sequence, model, seed, true);
					s.prepare();
					s.fold_everything();

					SecondarySequence predicted = s.getFoldedSequence();

					System.out.println(s.toString());
					System.out.println(RNAParser.emit(sequence));
					System.out.println("seed: "+seed+", probability: "+s.getProbability());

					if(solution != null){
						Statistics stat = predicted.compare(solution);
						System.out.println("sensitivity: "+stat.sensitivity()+", PPV: "+stat.ppv());
						System.out.println(solution.toString());
					}
				}catch(NumberFormatException e){
					System.out.println("No valid command. Enter 'help' for a helpful message");
				}
			}

		}
		System.out.println("Game over");
	}

	private static void usage(){
		System.out.println("usage: java Pantalaimon  grammar_parameters  sequence");
		System.out.println("  where");
		System.out.println("    grammar_parameters is a file containing transition and emission probabilities. Use 'hl1hel1' or 'hl3hel2' to use one of the shipped configurations");
		System.out.println("    sequence is the primary structure, i.e. a word consisting of A, C, G, U (case insensitive)");
	}
}
