import java.io.*;

public class FileParser {
	public static void main(String args[]) throws IOException {
		File targetDir = new File(args[0]);
		File[] files = targetDir.listFiles();
		
		for(File file : files) {
			if(file.isDirectory()) continue;
			
			//determine the capacity of this run
			String[] parameterSplit = file.getName().split("_");
			int capacity = Integer.parseInt(parameterSplit[1]);
			
			BufferedReader in = new BufferedReader(new FileReader(file));
			PrintWriter out = new PrintWriter(new FileWriter(file.getAbsolutePath() + ".csv"));
		
			String line;
			int counter = 0;
			int totalCounter = 0;
			
			double[] column_total = new double[4 + capacity];
			int[] column_position = new int[4 + capacity];
			column_position[0] = 5 + capacity;
			column_position[1] = 4 + capacity;
			column_position[2] = 3 + capacity;
			column_position[3] = 2 + capacity;

			for(int i = 0; i < capacity; i++) {
				column_position[i + 4] = capacity - i;
			}
			
			double[] column_tmp = new double[4 + capacity];
		
			while((line = in.readLine()) != null) {
				if(!line.startsWith("it")) {
					continue;
				}
				
				String[] parts = line.split(" ");
				
				if(parts.length < column_position[0]) {
					System.out.println("error: " + line);
					continue;
				}
				
				for(int i = 0; i < column_position.length; i++) {
					String target = parts[parts.length - column_position[i]];
					
					if(target.isEmpty() || target.equals("nan")) {
						column_tmp[i] = 0;
					} else {
						try {
							column_tmp[i] = Double.parseDouble(target);
						} catch(NumberFormatException e) {
							column_tmp[i] = 0;
							System.out.println(file.getAbsolutePath() + ": " + target);
							e.printStackTrace();
						}
					}
					
					column_total[i] += column_tmp[i];
				}
			
				if(counter % 100 == 99) {
					print(column_tmp, out);
				}
			
				counter++;
				totalCounter++;
			}
			
			double[] column_average = new double[column_position.length];
			
			for(int i = 0; i < column_position.length; i++) {
				column_average[i] = column_total[i] / totalCounter;
			}
			
			print(column_total, out);
			print(column_average, out);
			
			out.close();
			in.close();
		}
	}

	public static void print(double[] array, PrintWriter out) {
		//print first n-1 columns
		for(int i = 0; i < array.length - 1; i++) {
			out.print(array[i] + ",");
		}
		
		//print last column
		out.println(array[array.length - 1]);
	}
}
