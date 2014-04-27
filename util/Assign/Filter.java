import java.util.*;
import java.io.*;

public class Filter {
	public static void main(String args[]) throws IOException {
		BufferedReader in = new BufferedReader(new FileReader(args[0]));
		BufferedReader filterIn = new BufferedReader(new FileReader(args[1]));
		
		ArrayList<Integer> filter = new ArrayList<Integer>();
		String line;
		
		while((line = filterIn.readLine()) != null) {
			if(!line.trim().equals("")) {
				filter.add(Integer.parseInt(line.split(" ")[0]));
			}
		}
		
		filterIn.close();
		
		int total = Integer.parseInt(in.readLine());
		System.out.println(total - filter.size());
		
		for(int i = 0; i < total; i++) {
			line = in.readLine();
			if(!filter.contains(i)) System.out.println(line);
		}
	}
}
