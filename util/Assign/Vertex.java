import java.io.*;
import java.util.ArrayList;

public class Vertex {
	//some identifier for this vertex
	// usually index in some array that this vertex is found in
	int id;

	//coordinate information on this vertex (GPS system)
	double x;
	double y;
	
	//list of neighbors that this vertex is attached to
	ArrayList<Neighbor> neighbors;

	public Vertex() {
		id = -1;
		x = Double.NaN;
		y = Double.NaN;
		neighbors = new ArrayList<Neighbor>();
	}

	public Vertex(double x, double y) {
		id = -1;
		this.x = x;
		this.y = y;
		neighbors = new ArrayList<Neighbor>();
	}

	public Vertex(int id) {
		this.id = id;
		x = Double.NaN;
		y = Double.NaN;
		neighbors = new ArrayList<Neighbor>();
	}

	public Vertex(int id, double x, double y) {
		this.id = id;
		this.x = x;
		this.y = y;
		neighbors = new ArrayList<Neighbor>();
	}
	
	public ArrayList<Neighbor> getNeighbors() {
		return neighbors;
	}
	
	//assumes distance is Euclidean distance from here to parameter vertex
	public void addNeighobr(Vertex vertex) {
		neighbors.add(new Neighbor(vertex, euclideanDistanceTo(vertex)));
	}
	
	public void addNeighbor(Vertex vertex, double distance) {
		neighbors.add(new Neighbor(vertex, distance));
	}

	public double euclideanDistanceTo(Vertex vertex) {
		//calculate Euclidean distance from here to parameter vertex
		return euclideanDistance(x, y, vertex.x, vertex.y);
	}
	
	public static double euclideanDistance(double sx, double sy, double ex, double ey) {
		double dx = sx - ex;
		double dy = sy - ey;
		return Math.sqrt(dx * dx + dy * dy);
	}
	
	//loads vertices from a .gra file
	//format is:
	// first line contains number of vertices
	// the nth line (excluding first line) is for the nth vertex, at index/id n - 1
	//  line contains "id: <neighbor1 id> <neighbor1 distance> <neighbor2 id> <neighbor2 distance> ..."
	//gra format holds no information about location of vertex itself
	//also note that distances are not necessarily Euclidean
	// and that the graph may be directed
	public static Vertex[] loadFromGra(String filename) throws IOException {
		BufferedReader in = new BufferedReader(new FileReader(filename));
		
		int numVertices = Integer.parseInt(in.readLine());
		Vertex[] vertices = new Vertex[numVertices];
		
		//allocate memory so that objects are properly linked via neighbors
		for(int i = 0; i < numVertices; i++) {
			vertices[i] = new Vertex(i);
		}
		
		//now read in adjacency list
		for(int i = 0; i < numVertices; i++) {
			String[] parts = in.readLine().split(":")[1].split(" ");
			
			for(int j = 1; j < parts.length - 1; j += 2) {
				int neighborId = Integer.parseInt(parts[j]);
				double distance = Double.parseDouble(parts[j + 1]);
				
				vertices[i].addNeighbor(vertices[neighborId], distance);
			}
		}
		
		in.close();
		return vertices;
	}
	
	//loads vertices from vertex list and edge list files
	//vertex list contains n, the number of vertices, followed by n lines containing position of a vertex
	// each veretx line is "x y", where x and y are doubles
	//edge list contains m, the number of edges, followed by m lines connecting a vertex i with a vertex j with distance d
	// each edge line is "i j d", where i and j are integers and d is a double
	//if parameter undirected is true, i should either always be < or always be > than j
	public static Vertex[] loadFromLists(String vertexFilename, String edgeFilename, boolean undirected) throws IOException {
		//first create vertices
		BufferedReader in = new BufferedReader(new FileReader(vertexFilename));
		
		int numVertices = Integer.parseInt(in.readLine());
		Vertex[] vertices = new Vertex[numVertices];
		
		for(int i = 0; i < numVertices; i++) {
			String[] position = in.readLine().split(" ");
			double x = Double.parseDouble(position[0]);
			double y = Double.parseDouble(position[1]);
			
			vertices[i] = new Vertex(i, x, y);
		}
		
		in.close();
		
		//now create edges
		in = new BufferedReader(new FileReader(edgeFilename));
		int numEdges = Integer.parseInt(in.readLine());
		
		for(int i = 0; i < numEdges; i++) {
			String[] line = in.readLine().split(" ");
			int a = Integer.parseInt(line[0]);
			int b = Integer.parseInt(line[1]);
			double d = Double.parseDouble(line[2]);
			
			vertices[a].addNeighbor(vertices[b], d);
			
			if(undirected) {
				vertices[b].addNeighbor(vertices[a], d);
			}
		}
		
		in.close();
		return vertices;
	}
	
	public static void saveAsGra(Vertex[] vertices, String filename) throws IOException {
		PrintWriter out = new PrintWriter(new FileWriter(filename), true);
		
		out.println(vertices.length);
		
		for(Vertex vertex : vertices) {
			out.print(vertex.id + ":");
			
			for(Neighbor neighbor : vertex.getNeighbors()) {
				out.printf(" %d %.7f", neighbor.vertex.id, neighbor.distance);
			}
			
			out.println();
		}
		
		out.close();
	}
	
	public static void saveAsLists(Vertex[] vertices, String vertexFilename, String edgeFilename) throws IOException {
		//write vertices first
		PrintWriter out = new PrintWriter(new FileWriter(vertexFilename), true);
		
		out.println(vertices.length);
		
		for(Vertex vertex : vertices) {
			out.println(vertex.x + "," + vertex.y);
		}
		
		out.close();
		
		//now write edges
		out = new PrintWriter(new FileWriter(edgeFilename), true);
		
		//first find total number of edges
		int totalEdges = 0;
		for(Vertex vertex : vertices) {
			totalEdges += vertex.getNeighbors().size();
		}
		
		out.println(totalEdges);
		
		//now print edges
		for(Vertex vertex : vertices) {
			for(Neighbor neighbor : vertex.getNeighbors()) {
				out.printf("%d %d %.7f", vertex.id, neighbor.vertex.id, neighbor.distance);
				out.println();
			}
		}
		
		out.close();
	}
}

class Neighbor {
	Vertex vertex; //the vertex neighbor
	double distance; //distance to neighbor (normally Euclidean distance)
	
	public Neighbor(Vertex vertex, double distance) {
		this.vertex = vertex;
		this.distance = distance;
	}
}
