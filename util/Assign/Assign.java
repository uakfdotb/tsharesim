import java.io.*;
import java.util.ArrayList;

public class Assign {
	public static void main(String args[]) throws IOException {
		//first, read the vertices
		ArrayList<Vertex> vertices = new ArrayList<Vertex>();
		BufferedReader in = new BufferedReader(new FileReader(args[0])); //args[0] is the vertices.out file
		String line = "";
		
		int counter = 0;
		
		in.readLine(); //first line is the number of vertices
		while((line = in.readLine()) != null) {
			//String[] parts = line.split(": ");
			
			//if(parts.length == 2) {
				//int id = Integer.parseInt(parts[0]);
				
				//parts = parts[1].split(",");
				String[] parts = line.split(" ");
				
				if(parts.length == 2) {
					double x = Double.parseDouble(parts[0]);
					double y = Double.parseDouble(parts[1]);
					
					vertices.add(new Vertex(counter, x, y));
					counter++;
				} else {
					System.out.println("Error: skipping invalid vertex: " + line);
				}
			//} else {
			//	System.out.println("Error: skipping invalid vertex: " + line);
			//}
		}
		
		//now assign the customers
		in = new BufferedReader(new InputStreamReader(System.in));
		
		while((line = in.readLine()) != null) {
			String[] parts = line.split("\t");
			
			if(parts.length == 5) {
				long time = Long.parseLong(parts[0]);
				double sx = Double.parseDouble(parts[1]);
				double sy = Double.parseDouble(parts[2]);
				double ex = Double.parseDouble(parts[3]);
				double ey = Double.parseDouble(parts[4]);
				
				int vertexStart = find(vertices, sx, sy);
				int vertexDest = find(vertices, ex, ey);
				System.out.println(time + "\t" + vertexStart + "\t" + vertexDest);
			} else {
				System.out.println("Error: line doesn't have five parts: " + line);
			}
		}
	}
	
	public static int find(ArrayList<Vertex> vertices, double x, double y) {
		//find the closest vertex to (x, y)
		Vertex minVertex = null;
		double minDistance = Double.MAX_VALUE;
		
		for(Vertex v : vertices) {
			double d = distance(v.x, v.y, x, y);
			
			if(d < minDistance) {
				minVertex = v;
				minDistance = d;
			}
		}
		
		return minVertex.id;
	}
	
	public static double distance(double sx, double sy, double ex, double ey) {
		double dx = ex - sx;
		double dy = ey - sy;
		return Math.sqrt(dx * dx + dy * dy);
	}
}

class Vertex {
	int id;
	double x;
	double y;
	
	public Vertex(int id, double x, double y) {
		this.id = id;
		this.x = x;
		this.y = y;
	}
}
