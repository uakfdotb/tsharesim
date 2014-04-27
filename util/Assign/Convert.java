//converts from edges.out to graph.gra
import java.io.*;
import java.util.*;

public class Convert {
	public static void main(String args[]) throws IOException {
		BufferedReader in = new BufferedReader(new InputStreamReader(System.in));
		int num_edges = Integer.parseInt(in.readLine());
		
		HashMap<Integer, ArrayList<Distance>> map = new HashMap<Integer, ArrayList<Distance>>();
		
		String line;
		while((line = in.readLine()) != null) {
			String[] parts = line.split(" ");
			int a = Integer.parseInt(parts[0]);
			int b = Integer.parseInt(parts[1]);
			double distance = Double.parseDouble(parts[2]);
			
			if(!map.containsKey(a)) {
				map.put(a, new ArrayList<Distance>());
			}
			
			if(!map.containsKey(b)) {
				map.put(b, new ArrayList<Distance>());
			}
			
			map.get(a).add(new Distance(b, distance));
			map.get(b).add(new Distance(a, distance));
		}
		
		int num_vertices = num_edges;
		
		for(int i = 0; i < num_edges; i++) {
			if(!map.containsKey(i)) {
				System.err.println("detected " + i + " vertices total");
				num_vertices = i;
				break;
			}
		}
		
		System.out.println(num_vertices);
		
		for(int i = 0; i < num_vertices; i++) {
			System.out.print(i + ":");
			
			Collections.sort(map.get(i));
			for(int j = 0; j < map.get(i).size(); j++) {
				System.out.print(" " + map.get(i).get(j).id + " " + map.get(i).get(j).d);
			}
			
			System.out.println();
		}
	}
}

class Distance implements Comparable {
	int id;
	double d;
	
	public Distance(int id, double d) {
		this.id = id;
		this.d = d;
	}
	
	public int compareTo(Object o) {
		return id - ((Distance) o).id;
	}
}
