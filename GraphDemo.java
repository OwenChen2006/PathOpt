/**
 * Demonstrates the calculation of shortest paths in the US Highway
 * network, showing the functionality of GraphProcessor and using
 * Visualize
 * @author v 0.9 Brandon Fain
 * @author Owen Astrachan
 */
import java.util.*;
import java.io.*;

public class GraphDemo {
    
    /**
     * Keys in the map are locations like "Durham NC" or "Portland OR" or "Portland ME",
     * the corresponding value is the Point for that location label
     */
    private Map<String, Point> myMap;
    GraphProcessor myGP;
    Visualize myViz;

    public GraphDemo(GraphProcessor gp, Visualize viz){
        myMap = new HashMap<>();
        myGP = gp;
        myViz = viz;
    }

    /**
     * File should be in uscities.csv format: town,state,latitude,longitude
     * i.e., Provo,UT,40.2457,-111.6457
     * @param filename is the name of a properly formatted file
     * @throws IOException if file reading failes
     * @return name of all locations as space delimited strings: "Provo UT"
     */
    public ArrayList<String> readData(String filename) throws IOException{
        Scanner s = new Scanner(new File(filename));
        ArrayList<String> list = new ArrayList<>();
        while (s.hasNextLine()) {
            String line = s.nextLine();
            String[] data = line.split(",");
            String name = data[0] + " " + data[1];
            list.add(name);
            myMap.put(name, new Point(Double.parseDouble(data[2]),
                                      Double.parseDouble(data[3])));
        }
        s.close();
        return list;
    }

    public void showRoute(Point start, Point end){
        double stime = System.nanoTime();
        List<Point> path = myGP.route(start, end);
        double etime = System.nanoTime();
        double dtime = (etime-stime)/1e9;

        double dist = myGP.routeDistance(path);
        System.out.printf("start: %s, end: %s\n",
                          start,end);
        System.out.printf("shortest path has %d points\n",path.size());
        System.out.printf("shortest path is %2.3f in length\n",dist);
        System.out.printf("pq time = %1.3f\n",dtime);

        myViz.drawRoute(path);
    }

    /**
     * Let user choose start and end cities, then
     * find and visualize route between these locations 
     * @param list is choices users will choose from
     */
    public void userInteract(ArrayList<String> list) {
        
        String start = "Miami FL";
        String end = "Seattle WA";

        start = CitySelector.showDialog(list, "Choose Start City");
        end = CitySelector.showDialog(list, "Choose End City");
        System.out.printf("route from %s to %s\n",start,end);

        StdDraw.setTitle(String.format("201 Route from %s to %s",start,end));
      
        // find points in graph near starting point
        Point nearStart = myGP.nearestPoint(myMap.get(start));
        Point nearEnd = myGP.nearestPoint(myMap.get(end));
        showRoute(nearStart,nearEnd);       
    }
   
    public static void main(String[] args) throws IOException {
        String usaCityFile = "data/uscities.csv";

        String[] durhamData = {"images/durham.png", 
                               "data/durham.vis",
                               "data/durham.graph"};
        String[] usaData = {"images/usa.png",
                            "data/usa.vis",
                            //"data/usa-lower48.graph"};
                            "data/usa.graph"};

        String[] simpleData = {"images/simple.png",
                               "data/simple.vis",
                               "data/simple.graph"};

        // useThisData can point to durham, simple, or usa, modify to test
        String[] useThisData = usaData;
        GraphProcessor gp = new GraphProcessor();
        gp.initialize(new FileInputStream(useThisData[2]));
        Visualize viz = new Visualize(useThisData[1],useThisData[0]);
        GraphDemo gd = new GraphDemo(gp,viz);
        ArrayList<String> list = gd.readData(usaCityFile);
        Collections.sort(list);
        gd.userInteract(list);
    }
}