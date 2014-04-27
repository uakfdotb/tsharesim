import java.io.*;
import java.util.*;

public class Order {
	public static void main(String args[]) throws IOException {
		//first, read the vertices
		HashMap<Integer, String> map = new HashMap<Integer, String>();
		BufferedReader in = new BufferedReader(new InputStreamReader(System.in));
		
		int n = Integer.parseInt(in.readLine());
		
		String line = "";
		while((line = in.readLine()) != null) {
			String[] parts = line.split(": ");
			int id = Integer.parseInt(parts[0]);
			map.put(id, parts[1].replace(",", " "));
			
			if(id < 0 || id >= n) {
				System.err.println("warning: bad id: " + id);
			}
		}
		
		System.out.println(n);
		
		for(int i = 0; i < n; i++) {
			if(map.containsKey(i)) {
				System.out.println(map.get(i));
			} else {
				System.err.println("warning: nonexistent at " + i);
			}
		}
	}
}
