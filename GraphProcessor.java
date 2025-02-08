import java.io.*;
import java.util.*;
import java.security.InvalidAlgorithmParameterException;

/**
 * Models a weighted graph of latitude-longitude points
 * and supports various distance and routing operations.
 * To do: Add your name(s) as additional authors
 * @author Brandon Fain
 * @author Owen Astrachan modified in Fall 2023
 *
 */
public class GraphProcessor {

    // include instance variables here
    int myVertexCount;
    int myEdgeCount;
    Map<Point,List<Point>> myGraph;

    public GraphProcessor(){
        myGraph = new HashMap<>();
    }

    /**
     * Creates and initializes a graph from a source data
     * file in the .graph format. Should be called
     * before any other methods work.
     * @param file a FileInputStream of the .graph file
     * @throws IOException if file not found or error reading
     */

     public void initialize(FileInputStream file) throws IOException {
        // TODO implement by reading info and creating graph
        Scanner reader = new Scanner(file);
        ArrayList<Point> pointList = new ArrayList<>();
        try {
            int numVertices = reader.nextInt();
            int numEdges = reader.nextInt();
            reader.nextLine();
            
            for (int i = 0; i < numVertices; i++) {

                String[] vertexInfo = reader.nextLine().split(" ");
                double myLatitude = Double.parseDouble(vertexInfo[1]);
                double myLongitude = Double.parseDouble(vertexInfo[2]);

                Point myPoint = new Point(myLatitude, myLongitude);
                pointList.add(myPoint);

                myGraph.putIfAbsent(myPoint, new ArrayList<Point>());
            }

            for (int i = 0; i < numEdges; i++) {
                
                String[] edgeInfo = reader.nextLine().split(" ");
                int indexU = Integer.parseInt(edgeInfo[0]);
                int indexV = Integer.parseInt(edgeInfo[1]);

                Point u = pointList.get(indexU);
                Point v = pointList.get(indexV);
                myGraph.get(u).add(v);
                myGraph.get(v).add(u);
            
            }
        } catch (IllegalStateException | NoSuchElementException | NumberFormatException e){
            throw new IOException("Could not read .graph file");
        }
    }


    /**
     * Searches for the point in the graph that is closest in
     * straight-line distance to the parameter point p
     * @param p is a point, not necessarily in the graph
     * @return The closest point in the graph to p
     */
    public Point nearestPoint(Point p) {
        Point close = null;
        double dist = Double.MAX_VALUE;
 
        for (Point p1 : myGraph.keySet()) {
         double realDistance = p.distance(p1);
         if (realDistance < dist) {
             dist = realDistance;
             close = p1;
         }
        }
 
         return close;
     }
 


    /**
     * Calculates the total distance along the route, summing
     * the distance between the first and the second Points, 
     * the second and the third, ..., the second to last and
     * the last. Distance returned in miles.
     * @param start Beginning point. May or may not be in the graph.
     * @param end Destination point May or may not be in the graph.
     * @return The distance to get from start to end
     */
    public double routeDistance(List<Point> route) {
        double d = 0.0;
        for (int i = route.size() -1; i > 0; i--) {
            d += route.get(i).distance(route.get(i - 1));
        }
        return d;
    }

    

    /**
     * Checks if input points are part of a connected component
     * in the graph, that is, can one get from one to the other
     * only traversing edges in the graph
     * @param p1 one point
     * @param p2 another point
     * @return true if and only if p2 is reachable from p1 (and vice versa)
     */
    public boolean connected(Point p1, Point p2) {
        HashSet<Point> vis = new HashSet<Point>();
        Stack<Point> stack = new Stack<>();
        boolean result = false;
        stack.push(p1);
        
        while (stack.size() != 0) {
            Point lookingAt = stack.pop();

            if (lookingAt.equals(p2)) {
                result = true;
                return result;
            }
            vis.add(lookingAt);
            for (Point surroundings: myGraph.get(lookingAt)) {
                if (!vis.contains(surroundings)) {
                    stack.push(surroundings);
                }
            }
        }
        return result;
    }



    /**
     * Return a list of all points that lead to end
     * @param predMap contains predecessor info such that
     * predMap.get(end) is vertex before end in route
     * assumption: if predMap.get(v) == null, then v is start of route
     * @param end last value in route
     * @return list of points that lead to end
     */
    private List<Point> findPath(Map<Point, Point> predMap, Point end) {
        LinkedList<Point> path = new LinkedList<>();
        Point curr = end;
    
        while (curr != null) {
            path.addFirst(curr);
            curr = predMap.get(curr);
        }
    
        return path;
    }
    
   
    private Point minPoint(Map<Point,Double> map, Set<Point> visited){
        Point minP = null;
        double minD = Double.MAX_VALUE;
        for(Point p : map.keySet()){
            if (! visited.contains(p) && map.get(p) < minD){
                minD = map.get(p);
                minP = p;
            }
        }
        return minP;
    }


    /**
     * Slow implementation of shortest path
     * @param start
     * @param end
     * @return list of points from start to end
     */
    public List<Point> slowRoute(Point start, Point end){

        if (! myGraph.keySet().contains(start) || ! connected(start, end)){
            throw new IllegalArgumentException("not connected");
        }

        HashMap<Point,Point> predMap = new HashMap<>();
        predMap.put(start,null);
        HashSet<Point> visited = new HashSet<>();
        HashMap<Point,Double> distanceMap = new HashMap<>();
        distanceMap.put(start,0.0);

        while (distanceMap.size() < myGraph.size()){
            Point current = minPoint(distanceMap,visited);
            visited.add(current);
            
            if (current.equals(end)){
                break;
            }
            double d = distanceMap.get(current);
            for(Point p : myGraph.get(current)){
                double weight = p.distance(current);
                double newDistance = d + weight;

                if (!visited.contains(p) &&  newDistance < distanceMap.getOrDefault(p,Double.MAX_VALUE)){
                    distanceMap.put(p,newDistance);
                    predMap.put(p,current);
                }
            }
        }
        List<Point> list =  findPath(predMap,end);
        return list;
    }


    /**
     * Returns the shortest path, traversing the graph, that begins at start
     * and terminates at end, including start and end as the first and last
     * points in the returned list. If there is no such route, either because
     * start is not connected to end or because start equals end, throws an
     * exception.
     * @param start Beginning point.
     * @param end Destination point.
     * @return The shortest path [start, ..., end].
     * @throws IllegalArgumentException if there is no such route, 
     * either because start is not connected to end or because start equals end.
     */
    
    
     public List<Point> route(Point start, Point end) throws IllegalArgumentException {
      
        if (!myGraph.containsKey(start)) {
            throw new IllegalArgumentException("Start point not found in the graph.");
        }
        if (!myGraph.containsKey(end)) {
            throw new IllegalArgumentException("End point not found in the graph.");
        }
        if (start.equals(end)) {
            throw new IllegalArgumentException("Start and end points must be different.");
        }
        if (!connected(start, end)) {
            throw new IllegalArgumentException("No path exists between the start and end points.");
        }
    
        HashMap<Point, Double> distanceMap = new HashMap<>();
        HashMap<Point, Point> predecessorMap = new HashMap<>();
    
      
        PriorityQueue<Point> priorityQueue = new PriorityQueue<>(
            (p1, p2) -> Double.compare(distanceMap.getOrDefault(p1, Double.MAX_VALUE),
                                       distanceMap.getOrDefault(p2, Double.MAX_VALUE))
        );
    
   
        distanceMap.put(start, 0.0);
        priorityQueue.add(start);
    

        while (!priorityQueue.isEmpty()) {
            Point current = priorityQueue.remove();
    
          
            for (Point neighbor : myGraph.get(current)) {
                double newDistance = distanceMap.getOrDefault(current, Double.MAX_VALUE) + current.distance(neighbor);
    
                // Update distance and predecessor if a shorter path is found
                if (newDistance < distanceMap.getOrDefault(neighbor, Double.MAX_VALUE)) {
                    distanceMap.put(neighbor, newDistance);
                    priorityQueue.add(neighbor);
                    predecessorMap.put(neighbor, current);
                }
            }
        }
    
      
        List<Point> path = new ArrayList<>();
        path.add(end);
        Point current = end;
    
        while (!current.equals(start)) {
            current = predecessorMap.get(current);
            if (current == null) { 
                throw new IllegalArgumentException("No path exists between the start and end points.");
            }
            path.add(current);
        }
    

        Collections.reverse(path);
        return path;
    }
    
    
    
   
    
    
    
    public static void main(String[] args) throws FileNotFoundException, IOException {
        String name = "data/usa.graph";
        name = "data/simple.graph";
        GraphProcessor gp = new GraphProcessor();
        gp.initialize(new FileInputStream(name));
        System.out.println("running GraphProcessor");
    }


    
}
