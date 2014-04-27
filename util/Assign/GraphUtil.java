public class GraphUtil {
	public static double DELTA = 0.00000001;

	public static void main(String args[]) throws java.io.IOException {
		if(args.length < 3) {
			System.out.println("usage: GraphUtil vertices-file edges-file action");
			System.out.println("\taction may be one of: bound, grid, grid-est, togra");
			System.out.println("\tactions may require more parameters");
			return;
		}
		
		String vertexFilename = args[0];
		String edgeFilename = args[1];
		String action = args[2].toLowerCase(); //case insensitive
		
		//first, construct the graph
		Vertex[] vertices = Vertex.loadFromLists(vertexFilename, edgeFilename, true); //true for undirected graph
		
		if(action.equals("bound") || action.equals("grid") || action.equals("grid-est")) {
			//find the bounding rectangle for the graph
			// if we want to grid, we have to find bound first
			double minX = Double.MAX_VALUE;
			double maxX = Double.MIN_VALUE;
			double minY = Double.MAX_VALUE;
			double maxY = Double.MIN_VALUE;
			
			for(Vertex vertex : vertices) {
				minX = Math.min(minX, vertex.x);
				maxX = Math.max(maxX, vertex.x);
				minY = Math.min(minY, vertex.y);
				maxY = Math.max(maxY, vertex.y);
			}
			
			System.out.printf("bounding rectangle: (%.7f, %.7f) to (%.7f, %.7f)", minX, minY, maxX, maxY);
			System.out.println();
			
			if(action.equals("grid") || action.equals("grid-est")) {
				//split graph into a grid of user-set rows and columns
				// grid will calculate the total length of edges within each cell of the grid
				// grid-est will estimate this by assuming entire edge of a vertex inside a cell is within the cell
				boolean estimate = action.equals("grid-est");
			
				if(args.length < 5) {
					System.out.println("usage: GraphUtil vertices-file edges-file grid rows(x) cols(y)");
					return;
				}
				
				double xLength = maxX - minX;
				double yLength = maxY - minY;
			
				int rows = Integer.parseInt(args[3]);
				int cols = Integer.parseInt(args[4]);
				
				double gridXLength = xLength / rows;
				double gridYLength = yLength / cols;
			
				double[][] gridDistance = new double[rows][cols];
			
				//default all distances to 0 (unnecessary but for cleanness)
				for(int i = 0; i < rows; i++) {
					for(int j = 0; j < cols; j++) {
						gridDistance[i][j] = 0;
					}
				}
			
				//go through each edge
				for(Vertex vertex : vertices) {
					for(Neighbor neighbor : vertex.getNeighbors()) {
						//undirected graph, so don't do everything twice
						if(vertex.id >= neighbor.vertex.id) continue;
						
						//normalize x and y; simple translation doesn't affect lengths
						double sx = vertex.x - minX;
						double sy = vertex.y - minY;
						double ex = neighbor.vertex.x - minX;
						double ey = neighbor.vertex.y - minY;
							
						if(estimate) {
							//estimation simply uses the beginning vertices' grid
							int gridX = (int) (sx / gridXLength);
							int gridY = (int) (sy / gridYLength);
							gridDistance[gridX][gridY] += Vertex.euclideanDistance(vertex.x, vertex.y, neighbor.vertex.x, neighbor.vertex.y);
						} else {
							//calculate slope of line from s to e
							// also do horizontal and vertical checks
							boolean horizontal = sy == ey;
							boolean vertical = sx == ex;
						
							double slope = 0;
							if(!vertical) slope = (ey - sy) / (ex - sx);
							else if(ey > sy) slope = Double.POSITIVE_INFINITY;
							else slope = Double.NEGATIVE_INFINITY;
					
							//see whether we go left/right and up/down from s to e
							int xSign = 1;
							int ySign = 1;
					
							if(ex < sx) xSign = -1;
							if(ey < sy) ySign = -1;
					
							//iterate from s until we reach e and have recorded all lengths
							int counter = 0;
							while(true) {
								//find next position of s
								// this may be either a result of a collision with row or collision with col
								// or if we have reached e
							
								//first assume collisions and find next row / col info that we'll hit
								double nextRowX = 0, nextRowY = 0, nextColX = 0, nextColY = 0;
								double nx, ny; //final next values
							
								if(!vertical) {
									nextRowX = ((int) ((sx + xSign * gridXLength) / gridXLength)) * gridXLength;
									
									if(nextRowX - sx < DELTA) {
										nextRowX += xSign * gridXLength;
									}
									
									nextRowY = sy + (nextRowX - sx) * slope;
								}
							
								if(!horizontal) {
									nextColY = ((int) ((sy + ySign * gridYLength) / gridYLength)) * gridYLength;
									
									if(nextColY - sy < DELTA) {
										nextColY += ySign * gridYLength;
									}
									
									nextColX = sx + (nextColY - sy) / slope;
								}
							
								//determine which row / col info is correct
								if(vertical) {
									nx = nextColX;
									ny = nextColY;
								} else if(horizontal) {
									nx = nextRowX;
									ny = nextRowY;
								} else {
									//this is the only case where we have to really choose
									// choice is which is closer to our current position
									if(Math.abs(nextColX - sx) < Math.abs(nextRowX - sx)) {
										nx = nextColX;
										ny = nextColY;
									} else {
										nx = nextRowX;
										ny = nextRowY;
									}
								}
							
								//now make sure that e is not closer to s than n
								// because it may be horizontal/vertical, easiest way is to use Euclidean (instead of absolute value)
								if(Vertex.euclideanDistance(sx, sy, ex, ey) <= Vertex.euclideanDistance(sx, sy, nx, ny)) {
									nx = ex;
									ny = ey;
								}
								
								//also check if n and e are very close
								else if(Math.abs(nx - ex) < DELTA && Math.abs(ny - ey) < DELTA) {
									nx = ex;
									ny = ey;
								}
							
								//now, n is guaranteed to be in same grid as s
								// so add length to that grid
								// we find grid by using midpoint because s and n might both be on the edge
								double mx = (nx + sx) / 2;
								double my = (ny + sy) / 2;
								int gridX = (int) (mx / gridXLength);
								int gridY = (int) (my / gridYLength);
								
								gridDistance[gridX][gridY] += Vertex.euclideanDistance(sx, sy, nx, ny);
							
								//are we done?
								if(nx == ex && ny == ey) break;
							
								//for next iteration, reset s
								sx = nx;
								sy = ny;
							}
						}
					}
				}
				
				double total = 0;
				//now print grid distances
				for(int i = 0; i < rows; i++) {
					for(int j = 0; j < cols; j++) {
						System.out.printf("(%.7f, %.7f) to (%.7f, %.7f) -> (%d, %d) to (%d, %d): %.7f", minX + i * gridXLength, minY + j * gridYLength, minX + (i + 1) * gridXLength, minY + (j + 1) * gridYLength, i, j, i + 1, j + 1, gridDistance[i][j]);
						System.out.println();
						
						total += gridDistance[i][j];
					}
				}
				
				System.out.println("total distance: " + total);
			}
		} else if(action.equals("togra")) {
			//convert to gra format (see Vertex.java)
			
			if(args.length < 4) {
				System.out.println("usage: GraphUtil vertices-file edges-file togra output-file");
				return;
			}
			
			Vertex.saveAsGra(vertices, args[3]);
		}
	}
}
