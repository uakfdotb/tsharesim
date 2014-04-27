import java.util.*;
import java.io.*;

public class Complete {
	public static Vertex[] vertices;

	public static void main(String args[]) throws IOException {
		//construct the graph
		BufferedReader in = new BufferedReader(new InputStreamReader(System.in));
		
		int numVertices = Integer.parseInt(in.readLine());
		vertices = new Vertex[numVertices];
		
		for(int i = 0; i < numVertices; i++) {
			vertices[i] = new Vertex();
			String[] parts = in.readLine().split(":")[1].split(" ");
			
			for(int j = 1; j < parts.length - 1; j += 2) {
				vertices[i].neighbors.add(Integer.parseInt(parts[j]));
				vertices[i].distances.add(Double.parseDouble(parts[j + 1]));
			}
		}
		
		//done constructing graph
		
		graToEdges();
	}
	
	public static void graToEdges() {
		//convert to edge list
		int totalEdges = 0; //this counts each edge twice (undirected graph)
		
		for(int i = 0; i < vertices.length; i++) {
			totalEdges += vertices[i].neighbors.size();
		}
		
		System.out.println(totalEdges / 2);
		
		for(int i = 0; i < vertices.length; i++) {
			for(int j = 0; j < vertices[i].neighbors.size(); j++) {
				//don't duplicate edges
				if(i < vertices[i].neighbors.get(j)) System.out.println(i + " " + vertices[i].neighbors.get(j) + " " + vertices[i].distances.get(j));
			}
		}
	}
	
	public static void detectEdgeProblems() {
		for(int i = 0; i < vertices.length; i++) {
			for(int j = 0; j < vertices[i].neighbors.size(); j++) {
				//search for corresponding edges
				if(!vertices[vertices[i].neighbors.get(j)].neighbors.contains(i)) System.out.println("corresponding " + i + "," + j);
				
				//search for loop edges
				if(vertices[i].neighbors.get(j) == i) System.out.println("loop " + i);
			}
		}
	}
	
	public static void removeSeparated() {
		//flood fill at vertices[0]
		boolean[] visited = floodfill(0);
		
		//now visited shows false for unvisited vertices; remove these
		// we have to update vertex indexes because we are updating them
		// first create a map to new indexes
		
		int[] idmap = new int[vertices.length];
		int idcounter = 0;
		
		for(int i = 0; i < vertices.length; i++) {
			if(visited[i]) {
				idmap[i] = idcounter;
				idcounter++;
			} else {
				System.err.println(i + " is unvisited");
				idmap[i] = -1;
			}
		}
		
		//now write new vertices to stdout using idmap
		// idcounter is now the total number of new vertices (max index + 1)
		System.out.println(idcounter);
		
		for(int i = 0; i < vertices.length; i++) {
			if(idmap[i] >= 0) { //equivalent to if(visited[i])
				System.out.print(idmap[i] + ":");
				Vertex vertex = vertices[i];
			
				for(int j = 0; j < vertex.neighbors.size(); j++) {
					//while we're at it, remove loops
					if(i == vertex.neighbors.get(j)) continue;
					
					//we don't have to check idmap for -1
					// because if neighbors of a visited vertex are unvisited vertices
					// then those unvisited vertices should have been visited
					System.out.print(" " + idmap[vertex.neighbors.get(j)] + " " + vertex.distances.get(j));
				}
			
				System.out.println();
			}
		}
	}
	
	public static boolean[] floodfill(int x) {
		boolean[] visited = new boolean[vertices.length];
		for(int i = 0; i < visited.length; i++) {
			visited[i] = false;
		}
		
		visited[x] = true;
		floodfill(x, visited);
		return visited;
	}
	
	public static void floodfill(int x, boolean[] visited) {
		Vertex vertex = vertices[x];
		
		for(int i = 0; i < vertex.neighbors.size(); i++) {
			int id = vertex.neighbors.get(i);
			
			if(!visited[id]) {
				visited[id] = true;
				floodfill(id, visited);
			}
		}
	}
}

class Vertex {
	ArrayList<Integer> neighbors;
	ArrayList<Double> distances;

	public Vertex() {
		neighbors = new ArrayList<Integer>();
		distances = new ArrayList<Double>();
	}
}
