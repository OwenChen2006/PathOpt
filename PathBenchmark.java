import java.io.*;
import java.util.*;

/**
 * Simple benchmark to test some paths
 * in paths from us cities that are close
 * since slowRoute for such cities might be fast enough
 */

public class PathBenchmark {

     public Map<String,Point> readData(String filename) throws IOException{
        Scanner s = new Scanner(new File(filename));
        Map<String,Point> map = new HashMap<>();
        while (s.hasNextLine()) {
            String line = s.nextLine();
            String[] data = line.split(",");
            String name = data[0] + " " + data[1];
            map.put(name, new Point(Double.parseDouble(data[2]),
                                      Double.parseDouble(data[3])));
        }
        s.close();
        return map;
    }
    
    public void benchmark() throws FileNotFoundException, IOException{
        String[][] dataPairs = {
            new String[]{"White Plains NY", "Boston MA"},
            new String[]{"Wilmington NC", "Cherokee NC"},
            new String[]{"Phoenix AZ", "Las Vegas NV"},
            new String[]{"Missoula MT", "Boise ID"}
        };
        String usFile = "data/usa.graph";
        String usCities = "data/uscities.csv";
        GraphProcessor gp = new GraphProcessor();
        gp.initialize(new FileInputStream(usFile));
        Map<String,Point> map = readData(usCities);

        for(String[] pairs : dataPairs) {
            Point start = gp.nearestPoint(map.get(pairs[0]));
            Point end = gp.nearestPoint(map.get(pairs[1]));

            double stime = System.nanoTime();
            List<Point> route = gp.route(start,end);
            double distance = gp.routeDistance(route);
            double etime = System.nanoTime();
            double time = (etime-stime)/1e9;
            System.out.printf("%1.3f secs, %1.3f at %d\t %s to %s\n",
                              time,distance,route.size(),pairs[0],pairs[1]);
        }
    }
    public static void main(String[] args)  {
        PathBenchmark bench = new PathBenchmark();
        try {
            bench.benchmark();
        } catch (IOException e) {
            System.out.println("error reading data");
        }
    }
}
