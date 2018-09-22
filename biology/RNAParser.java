package biology;

public class RNAParser{
	public static Base[] parse_rna(String rna){
		int l = rna.length();
		Base[] result = new Base[l];

		for(int i=0; i<l; i++){
			char c = rna.charAt(i);
			if(c == 'A' || c == 'a'){
				result[i] = Base.adenine;
			}else if(c == 'C' || c == 'c'){
				result[i] = Base.cytosine;
			}else if(c == 'G' || c == 'g'){
				result[i] = Base.guanine;
			}else if(c == 'U' || c == 'u'){
				result[i] = Base.uracil;
			}else{
				result[i] = null;
			}
		}

		return result;
	}

	public static String emit(Base[] sequence){
		String result = "";
		for(Base b : sequence){
			result += b.toString();
		}
		return result;
	}
}
